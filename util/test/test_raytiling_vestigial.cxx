/* https://github.com/WireCell/wire-cell-toolkit/issues/196
 */

#include "WireCellUtil/RayClustering.h"
#include "WireCellUtil/RayTiling.h"
#include <WireCellUtil/DetectorWires.h>
#include <WireCellUtil/Testing.h>
#include <WireCellUtil/RaySvg.h>
#include <WireCellUtil/Exceptions.h>
#include <vector>
#include <iostream>
#include <fstream>

using namespace WireCell;

int main(int argc, char* argv[])
{
    const double nudge=1e-3;

    WireSchema::StoreDB storedb;
    std::vector<size_t> nwires = {2400, 2400, 3456};
    std::vector<size_t> wire_offset = {0, 2400, 4800};
    std::vector<RaySvg::wire_vector_t> plane_wires;

    for (size_t ind=0; ind<3; ++ind) {
        auto& plane = get_append(storedb, ind, 0, 0, 0);
        size_t ini = storedb.wires.size();
        size_t got = DetectorWires::plane(storedb, plane, DetectorWires::uboone[ind]);
        AssertMsg(got == nwires[ind], "wrong number of wires in plane");

        auto beg = storedb.wires.begin()+ini;
        plane_wires.push_back(RaySvg::wire_vector_t(beg, beg+got));
    }
    DetectorWires::flat_channels(storedb);

    WireSchema::Store store(std::make_shared<WireSchema::StoreDB>(storedb));

    std::vector<std::vector<WireSchema::Wire>> wireobjs(3);
    std::vector<Vector> pitches(3);
    for (auto plane : store.planes()) {
        wireobjs[plane.ident] = store.wires(plane);
        pitches[plane.ident] = store.mean_pitch(plane);
    }

    std::string fname = argv[0];
    fname += ".json.bz2";
    WireSchema::dump(fname.c_str(), store);

    auto raypairs = WireSchema::ray_pairs(store, store.faces()[0]);
    // const int nlayers = raypairs.size();

    RayGrid::Coordinates coords(raypairs);
    RaySvg::Geom geom{coords, plane_wires};

    const double hit=1.0;
    std::vector<std::pair<int,int>> tries = {
        {1,1727},               // single blob
        {3,1727},               // 3 vestigial
    };

    auto span_size = [](const auto& s) -> int {
        int siz = s.bounds.second-s.bounds.first;
        std::cerr << "span_size="<<siz<<" [" << s.bounds.first << "," << s.bounds.second << "]\n";
        return siz;
    };

    auto do_blobs = [&](const RayGrid::activities_t& activities, const std::string& fname) -> RayGrid::blobs_t {
        auto blobs = RayGrid::make_blobs(coords, activities, nudge);

        std::cerr << fname << "\n";
        std::cerr << "do_blobs: nblobs=" << blobs.size() << " nactivites=" << activities.size() << "\n";

        auto top = RaySvg::svg_full(geom, activities, blobs);
        std::ofstream ofstr(fname);
        ofstr << svggpp::dumps(top);
        return blobs;
    };        

    auto do_one = [&](int span, int iwire) -> RayGrid::blobs_t {
        RayGrid::activities_t activities;
        activities.emplace_back(0, 1, hit);
        activities.emplace_back(1, 1, hit);
        activities.emplace_back(2, span+1, hit, nwires[0]/2);
        activities.emplace_back(3, span+2, hit, nwires[1]/2);
        activities.emplace_back(4, span, hit, iwire);

        std::stringstream ss;
        ss << argv[0] << "-one_s" << span << "-w" << iwire << ".svg";
        auto blobs = do_blobs(activities, ss.str());
        std::cerr << "do_one: nblobs=" << blobs.size() << " nactivites=" << activities.size() << "\n";
        return blobs;
    };

    auto do_three = [&](const std::vector<int> spans, const std::vector<int>& iwires) -> RayGrid::blobs_t {
        RayGrid::activities_t activities;
        activities.emplace_back(0, 1, hit);
        activities.emplace_back(1, 1, hit);
        std::stringstream ss;
        ss << argv[0] << "-three";
        for (int iplane=0; iplane<3; ++iplane) {
            int span = spans[iplane];
            int iwire = iwires[iplane];
            activities.emplace_back(iplane+2, span, hit, iwire);
            ss << "_p" << iplane << "-s" << span << "-w" << iwire;
        }
        ss << ".svg";
        auto blobs = do_blobs(activities, ss.str());
        std::cerr << "do_three: nblobs=" << blobs.size() << " nactivites=" << activities.size() << "\n";
        return blobs;
    };

    auto do_shot = [&](const std::vector<Point>& points, const std::string& name) -> RayGrid::blobs_t {
        std::cerr << "do_shot: npoints="<<points.size() << " " << name << "\n";

        RayGrid::activities_t activities;
        activities.emplace_back(0, 1, 1);
        activities.emplace_back(1, 1, 1);
        
        for (size_t ilayer=2; ilayer<5; ++ilayer) {
            const auto& pdir = coords.pitch_dirs()[ilayer];
            const auto& cen = coords.centers()[ilayer];

            std::vector<RayGrid::Activity::value_t> measure;

            for (const auto& point : points) {
                const auto vrel = point - cen;
                double pit = vrel.dot(pdir);
                size_t pind = coords.pitch_index(pit, ilayer);
                if (measure.size() <= pind) {
                    measure.resize(pind+1, 0);
                }
                measure[pind] += 1;
            }
            activities.push_back(RayGrid::Activity(ilayer, {measure.begin(), measure.end()}));
        }
        std::stringstream ss;
        ss << argv[0] << "-shot-" << name << ".svg";
        auto blobs = do_blobs(activities, ss.str());
        AssertMsg(blobs.size() > 0, "No shot blobs");
        return blobs;
    };

    // this is a dorky function
    auto spread = [](const Point& point, double size, int num=10) -> std::vector<Point> {
        const double spacing = size/num;

        std::vector<Point> ret;
        for (int ix=0; ix<num; ++ix) {
            for (int iy=0; iy<num; ++iy) {
                for (int iz=0; iz<num; ++iz) {
                    ret.emplace_back(point.x() + (num/2.0 - ix)*spacing,
                                     point.y() + (num/2.0 - iy)*spacing,
                                     point.z() + (num/2.0 - iz)*spacing);
                }
            }
        }
        return ret;
    };

    // vestigial wires should produce 3 corners and 2 strips per plane
    {                          
        auto blobs = do_one(3, 1727);

        if (blobs.size() != 1) {
            THROW(ValueError() << errmsg{"Unexpected number of blobs"});
        }
        const auto& blob = blobs[0];
        if (blob.corners().size() != 3) {
            THROW(ValueError() << errmsg{"Unexpected number of blob corners"});
        }
        const auto& strips = blob.strips();
        if (strips.size() != 5) {
            THROW(ValueError() << errmsg{"Unexpected number of blob strips"});
        }
        for (const auto& strip : strips) {
            int size = span_size(strip);
            std::cerr << "test one: " << strip.layer << " " << size << "\n";
            if (strip.layer <=1 ) {
                AssertMsg(size == 1, "Active area strip wrong size");
            }
            else {
                AssertMsg(size == 2, "Active wire strip wrong size");
            }
        }
    }

    { // should be same
        auto blobs = do_three({2,2,2}, {1200, 1200, 1728});
        AssertMsg(blobs.size(), "No blobs found for do_three same as ab ove");
    }

    // single-corner blob

    {
        auto blobs = do_one(1, 1727);
        if (blobs.size()) {
            THROW(ValueError() << errmsg{"No blobs expected with single-corner blob"});
        }
    }
    {
        auto blobs = do_one(3, 1725);
        if (blobs.size()) {
            THROW(ValueError() << errmsg{"No blobs expected with single-corner blob"});
        }
    }

    // some specifics
    { 
        auto blobs = do_three({6,6,5}, {338, 338, 8});
        AssertMsg(blobs.size(), "No blobs found for three specific");
    }
    // some discrepencies with WCP found by haiwang
    {
        auto blobs = do_three({6,6,5}, {338, 338, 8});
        AssertMsg(blobs.size(), "No blobs found for three specific");
    }

    // some points
    {
        Point origin = coords.centers()[4] + 0.1*3456*3*units::mm*coords.pitch_dirs()[4];
        auto points = spread(origin, 10*units::mm);
        AssertMsg(points.size(), "Got no points");
        auto blobs = do_shot(points, "left");
        AssertMsg(blobs.size(), "Got no blobs");
    }
    {
        Point origin = coords.centers()[4] + 0.5*3456*3*units::mm*coords.pitch_dirs()[4];
        auto points = spread(origin, 10*units::mm);
        AssertMsg(points.size(), "Got no points");
        auto blobs = do_shot(points, "middle");
        AssertMsg(blobs.size(), "Got no blobs");
    }
    {
        Point origin = coords.centers()[4] + 0.9*3456*3*units::mm*coords.pitch_dirs()[4];
        auto points = spread(origin, 10*units::mm);
        AssertMsg(points.size(), "Got no points");
        auto blobs = do_shot(points, "right");
        AssertMsg(blobs.size(), "Got no blobs");
    }
    

    return 0;
}
