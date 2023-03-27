#include "WireCellGen/WireGenerator.h"
#include "WireCellUtil/Intersection.h"
#include "WireCellUtil/NamedFactory.h"

#include <cmath>
#include <iostream>  // debugging
#include <vector>

WIRECELL_FACTORY(WireGenerator, WireCell::WireGenerator, WireCell::IWireGenerator)

using namespace WireCell;
using namespace std;

namespace WireCell {
    class GenWire : public IWire {
        WirePlaneId m_wpid;
        int m_index;
        Ray m_ray;

       public:
        GenWire(WirePlaneId wpid, int index, const Ray& ray)
          : m_wpid(wpid)
          , m_index(index)
          , m_ray(ray)
        {
        }
        virtual ~GenWire() {}

        int ident() const
        {
            int iplane = m_wpid.index();
            ++iplane;
            return iplane * 100000 + m_index;
        }

        WirePlaneId planeid() const { return m_wpid; }

        int index() const { return m_index; }

        int channel() const { return ident(); }

        WireCell::Ray ray() const { return m_ray; }

        /// fixme: supply
        int segment() const { return 0; }
        int face() const { return 0; }
        int apa() const { return 0; }

        void set_index(int ind) { m_index = ind; }
        void set_planeid(WirePlaneId wpid) { m_wpid = wpid; }
    };
}  // namespace WireCell

static GenWire* make_wire(int index, const Point& point, const Point& proto, const Ray& bounds)
{
    const Point pt1 = point;
    const Point pt2 = pt1 + proto;
    const Ray wireray(pt1, pt2);

    Ray hits;
    int hitmask = box_intersection(bounds, pt1, ray_unit(wireray), hits);
    if (3 != hitmask) {
        return 0;
    }
    // ray should point generally towards +Y
    if (hits.first.y() > hits.second.y()) {
        hits = Ray(hits.second, hits.first);
    }
    return new GenWire(WirePlaneId(kUnknownLayer), index, hits);
}

struct SortByIndex {
    inline bool operator()(const GenWire* lhs, const GenWire* rhs) { return lhs->index() < rhs->index(); }
};

static void make_one_plane(IWire::vector& returned_wires, WirePlaneId wpid, const Ray& bounds, const Ray& step)
{
    const Vector xaxis(1, 0, 0);
    const Point starting_point = step.first;
    const Vector pitch = step.second - starting_point;
    const Vector proto = pitch.cross(xaxis).norm();

    std::vector<GenWire*> these_wires;

    int pos_index = 0;
    Point offset = starting_point;
    while (true) {  // go in positive pitch direction
        GenWire* wire = make_wire(pos_index, offset, proto, bounds);
        if (!wire) {
            break;
        }
        these_wires.push_back(wire);
        offset = wire->center() + pitch;
        pos_index += 1;
    }

    int neg_index = -1;  // now go in negative pitch direction
    const Vector neg_pitch = -1.0 * pitch;
    offset = these_wires[0]->center() + neg_pitch;  // start one below first upward going one
    while (true) {                                  // go in negative pitch direction
        GenWire* wire = make_wire(neg_index, offset, proto, bounds);
        if (!wire) {
            break;
        }
        these_wires.push_back(wire);
        offset = wire->center() + neg_pitch;
        neg_index -= 1;
    }

    // order by index
    std::sort(these_wires.begin(), these_wires.end(), SortByIndex());

    // load in to store and fix up index and plane
    for (size_t ind = 0; ind < these_wires.size(); ++ind) {
        GenWire* pwire = these_wires[ind];
        pwire->set_index(ind);
        pwire->set_planeid(wpid);
        returned_wires.push_back(IWire::pointer(pwire));
    }

    // std::cerr << "Made "<<store.size()<<" wires for plane " << plane << std::endl;
    // std::cerr << "step = " << step << std::endl;
    // std::cerr << "bounds = " << bounds << std::endl;
}

bool WireGenerator::operator()(const input_pointer& wp, output_pointer& wires)
{
    if (!wp) {
        wires = nullptr;
        return true;
    }
    IWire::vector* pwires = new IWire::vector;
    make_one_plane(*pwires, WirePlaneId(kUlayer), wp->bounds(), wp->pitchU());
    make_one_plane(*pwires, WirePlaneId(kVlayer), wp->bounds(), wp->pitchV());
    make_one_plane(*pwires, WirePlaneId(kWlayer), wp->bounds(), wp->pitchW());
    wires = output_pointer(pwires);
    return true;
}

WireGenerator::WireGenerator() {}

WireGenerator::~WireGenerator() {}
