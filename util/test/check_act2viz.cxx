#include "WireCellUtil/String.h"
#include "WireCellUtil/Testing.h"
#include <WireCellUtil/DetectorWires.h>
#include <WireCellUtil/RaySvg.h>

#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>

namespace po = boost::program_options;
using namespace WireCell;
using namespace svggpp;
using namespace WireCell::String;

template<typename D>
void usage(const D& d, const std::string& msg = "")
{
    std::cerr << d << std::endl;
    std::cerr << msg << std::endl;
    exit (-1);
}

WireSchema::Store load_microboone_wires(const std::string& fname, int cor=-1)
{
    if (cor < 0) {
        cor = (int)WireSchema::Correction::pitch;
    }
    return WireSchema::load(fname.c_str(), (WireSchema::Correction) cor);
}

WireSchema::Store generate_microboone_wires()
{
    //
    // MicroBooNE'ish
    //

    WireSchema::StoreDB storedb;
    std::vector<RaySvg::wire_vector_t> plane_wires;

    for (size_t ind=0; ind<3; ++ind) {
        auto& plane = get_append(storedb, ind, 0, 0, 0);
        DetectorWires::plane(storedb, plane, DetectorWires::uboone[ind]);

        // size_t ini = storedb.wires.size();
        // size_t got = DetectorWires::plane(storedb, plane, DetectorWires::uboone[ind]);
        // AssertMsg(got == nwires[ind], "wrong number of wires in plane");

        // // std::cerr << "add wires: " << ini << " + " << got << "\n";
        // auto beg = storedb.wires.begin()+ini;
        // plane_wires.push_back(RaySvg::wire_vector_t(beg, beg+got));
    }
    DetectorWires::flat_channels(storedb);
    return WireSchema::Store(std::make_shared<WireSchema::StoreDB>(storedb));
}

std::vector<RaySvg::wire_vector_t> get_plane_wires(WireSchema::Store& store)
{
    std::vector<RaySvg::wire_vector_t> ret;
    for (const auto& plane : store.planes()) {
        auto wires = store.wires(plane);
        ret.push_back(RaySvg::wire_vector_t(wires.begin(), wires.end()));
    }
    return ret;
}

int main(int argc, char* argv[])
{

    po::options_description desc("Read in activities, run tiling, generate SVG viz and dump blobs\n\tcheck_act2viz [options] activities.txt\n\nOptions");
    desc.add_options()("help,h", "produce help message")
        ("nudge,n", po::value<double>(),
         "the 'nudge' value used in tiling, default=1e-3")

        ("wires,w", po::value<std::string>(),
         "use the given microboone wires file instead of generating microboone'sh wires")

        ("correction,c", po::value<int>(),
         "correction level for wires if read from file, 1: load, 2: +order, 3: +pitch. default=3 (all)")

        ("dump,d", po::value<std::string>(),
         "dump blobs to given text file")

        ("output,o", po::value<std::string>(),
         "output vis to given SVG file")

        ("activities,a", po::value<std::string>(),
         "provide an activities file")

        ;

    po::positional_options_description pos;
    pos.add("activities", -1);

    po::variables_map opts;
    po::store(po::command_line_parser(argc, argv).
              options(desc).positional(pos).run(), opts);
    po::notify(opts);

    if (opts.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    double nudge=1e-3;
    if (opts.count("nudge")) {
        nudge = opts["nudge"].as<double>();
    }

    std::string wires_file{""};
    if (opts.count("wires")) {
        wires_file = opts["wires"].as<std::string>();
    }
    int wires_correction = -1;
    if (opts.count("correction")) {
        wires_correction = opts["correction"].as<int>();
    }

    std::string output = "/dev/stdout";
    if (opts.count("output")) {
        output = opts["output"].as<std::string>();
    }

    std::string dump{""};
    if (opts.count("dump")) {
        dump = opts["dump"].as<std::string>();
    }

    std::string activities_filename;
    if (opts.count("activities")) {
        activities_filename = opts["activities"].as<std::string>();
    }
    else usage(desc, "no activities file given");

    std::ifstream fstr(activities_filename);
    if (!fstr) usage(desc, "failed to open activities file");


    WireSchema::Store store;
    if (wires_file.empty()) {
        store = generate_microboone_wires();
    }
    else {
        store = load_microboone_wires(wires_file, wires_correction);
    }

    auto plane_wires = get_plane_wires(store);
    const std::vector<size_t> nwires = {2400, 2400, 3456};
    for (size_t ind=0; ind<3; ++ind) {
        AssertMsg(plane_wires[ind].size() == nwires[ind], "wrong number of wires");
    }

    auto raypairs = WireSchema::ray_pairs(store, store.faces()[0]);

    RaySvg::Geom geom{RayGrid::Coordinates(raypairs), plane_wires};

    std::vector<std::vector<RayGrid::Activity::value_t>> measures(5);
    measures[0].push_back(1);
    measures[1].push_back(1);

    const std::vector<int> nrays = {2, 2, 2400, 2400, 3456};

    // Parse file into mesures/activities
    while (fstr) {
        std::string line;
        std::getline(fstr, line);
        auto chunks = split(line, " ");

        if (chunks.empty()) continue;
        if (!endswith(chunks[0], "strip")) continue;

        int level = atoi(chunks[1].substr(1).c_str());
        if (level < 0 or level >= 5) {
            std::cerr << "bogus level! line:\n" << line;
            return -1;
        }
        if (level < 2) { continue; }

        auto pp = split(chunks[2].substr(6,chunks.size() - 8), ",");
        int pmin = std::max(0,              atoi(pp[0].c_str()));
        int pmax = std::min(nrays[level]-1, atoi(pp[1].c_str()));

        auto& m = measures[level];
        if (m.size() < (size_t)pmax) {
            m.resize(pmax, 0);
        }
        // A strip's upper most bound is not "in" the strip activities
        for (int pind=pmin; pind<pmax; ++pind) {
            m[pind] += 1;
        }
        // std::cerr << "level=" << level << " pind=["<<pmin<<","<<pmax<<"]\n";
    }
    RayGrid::activities_t activities;
    for (int layer=0; layer<5; ++layer) {
        const auto& m = measures[layer];
        activities.push_back(RayGrid::Activity(layer, {m.begin(), m.end()}));
    }

    auto blobs = RayGrid::make_blobs(geom.coords, activities, nudge);
    std::cerr << "nactivities=" << activities.size()
              << " nblobs=" << blobs.size() << "\n";
    if (! dump.empty()) {
        std::ofstream out(dump);
        for (const auto& blob : blobs) {
            out << blob.as_string();
        }
    }

    auto top = RaySvg::svg_full(geom, activities, blobs);

    std::ofstream ofstr(output);
    ofstr << svggpp::dumps(top);

    return 0;
}
