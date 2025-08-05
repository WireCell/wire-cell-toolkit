/** An envelope modeled with polygons stacked along a given axis.

    The envelope is described a list of polygonal slabs.

    Each slab is defined by:

    - two coordinate values along the axis giving the span of the slab.

    - a set of 2D points (the polygon) in the plane orthogonal to the axis.

    The 2D points are ordered according to the right hand rule.  If axis is the
    X axis then the points are ordered (y,z).  If axis is the Y axis (z,x).  If
    axis is the Z axis then (x,y).

    Note, slabs are checked for point containment in the order that they are
    provided.  A point is contained in the envelope as soon as the first slab is
    found to contain the point.


    This corresponds to "ToyFiducial" from WC prototype.
*/

#include "WireCellIface/IFiducial.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/BoundingBox.h"

#include <vector>
// #include <iostream>             // only debug

// Implementation is totally local to this comp. unit so no need for namespacing.
class PolyFiducial;

WIRECELL_FACTORY(PolyFiducial, PolyFiducial,
                 WireCell::IFiducial, WireCell::IConfigurable)


using namespace WireCell;

// "ray casting algorithm" aka "pnpoly"
static bool is_inside(double x, double y, const std::vector<double>& cx, const std::vector<double>& cy) {
    bool inside = false;
    const size_t n = cx.size();
    for (size_t i=0, j=n-1; i<n; j=i++) {

        if (((cy[i] > y) != (cy[j] > y)) &&
            // Calculate x-coordinate of intersection with the ray
            (x < (cx[j] - cx[i]) * (y - cy[i]) / (cy[j] - cy[i]) + cx[i])) {
            inside = !inside;
        }
    }
    return inside;
}        


struct PolySlab {
    BoundingBox bb;             // bounding box of 3D points of the slab
    std::vector<double> a,b;    // 2D points.

    bool inside(const Point& p, int axis=0) const {
        if (! bb.inside(p)) return false;

        const int i = axis;   // along axis
        const int j = (i+1)%3;  // along transverse coord a.
        const int k = (i+2)%3;  // along transverse coord b.

        return is_inside(p[j], p[k], a, b);
    }

};

class PolyFiducial : public IFiducial, public IConfigurable {

    int m_axis{0};
    BoundingBox m_bb;
    std::vector<PolySlab> m_slabs;

public:

    PolyFiducial() {}
    virtual ~PolyFiducial() {}

    virtual Configuration default_configuration() const {
        Configuration cfg;
        // The axis of the slabs. x=0, y=1, z=2
        cfg["axis"] = 0;        

        // Each slab is an object like:
        // {min:0, max:1, corners:[ [a1,b1], [a2,b2], ...]}
        // min/max give coordinate values along the axis.
        // corners is array of pairs of coordinates in the transverse plane.
        cfg["slabs"] = Json::arrayValue;
        return cfg;
    }

    virtual void configure(const Configuration& cfg) {
        m_axis = get(cfg, "axis", m_axis);

        const int i = m_axis;   // along axis
        const int j = (i+1)%3;  // along transverse coord a.
        const int k = (i+2)%3;  // along transverse coord b.

        m_slabs.clear();
        for (const auto& jslab : cfg["slabs"]) {
            PolySlab slab;

            Point p1,p2;
            p1[i] = get<double>(jslab, "min",0);
            p2[i] = get<double>(jslab, "max",0);
            if (p1[i] == p2[i]) {
                raise<ValueError>("poly slab has no width");
            }

            const auto& jcorners = jslab["corners"];
            if (jcorners.size() < 3) {
                raise<ValueError>("poly slab must have at least 3 corners");
            }

            for (const auto& jab : jcorners) {
                p1[j] = p2[j] = jab[0].asDouble();
                p1[k] = p2[k] = jab[1].asDouble();

                slab.bb(p1);
                slab.bb(p2);
                m_bb(slab.bb.bounds());
                slab.a.push_back(p1[j]);
                slab.b.push_back(p1[k]);               
            }
            m_slabs.push_back(slab);
        }
    }


    // IFiducial
    virtual bool contained(const Point& point) const {
        if (! m_bb.inside(point)) return false;

        for (const auto& slab : m_slabs) {
            if (slab.inside(point)) return true;
        }
        return false;
    }
};
