#include "WireCellUtil/RayClustering.h"
#include "WireCellUtil/RayTiling.h"
#include <WireCellUtil/DetectorWires.h>
#include <WireCellUtil/Testing.h>
#include <WireCellUtil/RaySvg.h>
#include <vector>
#include <iostream>
#include <fstream>

using namespace WireCell;

int main(int argc, char* argv[])
{
    WireSchema::StoreDB storedb;
    std::vector<size_t> nwires = {2400, 2400, 3456};
    std::vector<size_t> wire_offset = {0, 2400, 4800};
    std::vector<RaySvg::wire_vector_t> plane_wires;

    for (size_t ind=0; ind<3; ++ind) {
        auto& plane = get_append(storedb, ind, 0, 0, 0);
        size_t ini = storedb.wires.size();
        size_t got = DetectorWires::plane(storedb, plane, DetectorWires::uboone[ind]);
        std::cerr << ind << ": " << nwires[ind] << "\n";
        Assert(got == nwires[ind]);
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

    RayGrid::Coordinates coords(raypairs);
    RaySvg::Geom geom{coords, plane_wires};

    const double hit=1.0;
    for (size_t span=1; span<10; span +=2 ) 
    {
        RayGrid::activities_t activities;
        activities.emplace_back(0, 1, hit);
        activities.emplace_back(1, 1, hit);
        activities.emplace_back(2, span+1, hit, nwires[0]/2);
        activities.emplace_back(3, span+2, hit, nwires[1]/2);
        activities.emplace_back(4, span, hit, nwires[2]/2);
        // we overwrite this last in each step of the following scan

        for (size_t iwire=0; iwire<nwires[2]-span; ++iwire) {
            activities[4] = RayGrid::Activity(4, span, hit, iwire);
            auto blobs = RayGrid::make_blobs(coords, activities);

            if (blobs.empty()) {
                continue;
            }
            std::cerr << iwire << ": " << blobs.size() << "\n";

            std::stringstream ss;
            ss << argv[0] << "-span" << span << ".svg";

            auto top = RaySvg::svg_full(geom, activities, blobs);
            std::ofstream ofstr(ss.str());
            ofstr << svggpp::dumps(top);

        }

    }
    return 0;
}

