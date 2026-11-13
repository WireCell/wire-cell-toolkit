// Tests for YZMap (lookup, edge clamp, grid checks) and DepoTransformYZ
// (per-plane YZ charge scale inside DepoTransform).
//
// The DepoTransformYZ cases use the MicroBooNE wires and field response from
// wire-cell-data (WIRECELL_PATH, detectors.jsonnet), under their own instance names.

#include "WireCellGen/DepoTransform.h"
#include "WireCellGen/DepoTransformYZ.h"
#include "WireCellGen/YZMap.h"

#include "WireCellAux/SimpleDepo.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellAux/Testing.h"

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IFrame.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <map>
#include <string>
#include <vector>

using namespace WireCell;

namespace {

    // Write a map {key: {plane: grid[z][y]}} to a file in tdir and return its path.
    std::string write_map(const Persist::TempDir& tdir, const std::string& name,
                          const std::map<std::string, std::map<int, std::vector<std::vector<double>>>>& m)
    {
        Json::Value jtop;
        for (const auto& [key, planes] : m) {
            for (const auto& [plane, grid] : planes) {
                Json::Value jgrid(Json::arrayValue);
                for (const auto& row : grid) {
                    Json::Value jrow(Json::arrayValue);
                    for (double v : row) jrow.append(v);
                    jgrid.append(jrow);
                }
                jtop[key][std::to_string(plane)] = jgrid;
            }
        }
        const std::string path = (tdir.path / name).string();
        Persist::dump(path, jtop);
        return path;
    }

    Configuration yzmap_cfg(const std::string& filename, int nbinsz, int nbinsy,
                            double bin_width, double bin_height, double zoffset, double yoffset)
    {
        Configuration cfg;
        cfg["filename"] = filename;
        cfg["nbinsz"] = nbinsz;
        cfg["nbinsy"] = nbinsy;
        cfg["bin_width"] = bin_width;
        cfg["bin_height"] = bin_height;
        cfg["zoffset"] = zoffset;
        cfg["yoffset"] = yoffset;
        return cfg;
    }

    // Configure a YZMap through the factory so components can find it by type:name.
    void make_yzmap(const std::string& name, const Configuration& cfg)
    {
        auto icfg = Factory::lookup<IConfigurable>("YZMap", name);
        icfg->configure(cfg);
    }

    const double tick = 0.5 * units::us;
    const int nticks = 2000;

    // A MicroBooNE anode with its own wire schema, field response and plane
    // impact responses, under names no other test uses.  (The aux Testing
    // anodes helper reconfigures the shared "AnodePlane:<ident>" and the default
    // "WireSchemaFile" and "FieldResponse" instances, which other tests rely on.)
    struct Detector {
        IAnodePlane::pointer anode;
        std::string anode_tn;
        std::vector<std::string> pirs;
    };

    const Detector& detector()
    {
        static Detector det;
        if (det.anode) return det;

        // FftwDFT (Aux), the Gen components and FieldResponse (SigProc) are
        // found through their plugins, loaded at run time as Testing::anodes() does.
        Testing::load_plugins({"WireCellAux", "WireCellGen", "WireCellSigProc"});
        Testing::get_dft();  // creates the FftwDFT instance the responses look up

        auto kd = Testing::detectors()["uboone"];
        REQUIRE(!kd.isNull());
        {
            auto icfg = Factory::lookup<IConfigurable>("WireSchemaFile", "yzwires");
            auto cfg = icfg->default_configuration();
            cfg["filename"] = kd["wires"];
            icfg->configure(cfg);
        }
        {
            auto icfg = Factory::lookup<IConfigurable>("FieldResponse", "yzfr");
            auto cfg = icfg->default_configuration();
            cfg["filename"] = kd["fields"].isArray() ? kd["fields"][0] : kd["fields"];
            icfg->configure(cfg);
        }
        {
            auto icfg = Factory::lookup<IConfigurable>("AnodePlane", "yzanode");
            auto cfg = icfg->default_configuration();
            cfg["ident"] = 0;
            cfg["wire_schema"] = "WireSchemaFile:yzwires";
            cfg["faces"][0]["anode"] = 10 * units::cm - 6 * units::mm;
            cfg["faces"][0]["response"] = 10 * units::cm - 6 * units::mm;
            cfg["faces"][0]["cathode"] = 2.5604 * units::m;
            icfg->configure(cfg);
        }
        det.anode_tn = "AnodePlane:yzanode";
        det.anode = Factory::find_tn<IAnodePlane>(det.anode_tn);

        for (int plane = 0; plane < 3; ++plane) {
            const std::string name = "yzpir" + std::to_string(plane);
            auto icfg = Factory::lookup<IConfigurable>("PlaneImpactResponse", name);
            auto cfg = icfg->default_configuration();
            cfg["field_response"] = "FieldResponse:yzfr";
            cfg["plane"] = plane;
            cfg["nticks"] = nticks;
            cfg["tick"] = tick;
            icfg->configure(cfg);
            det.pirs.push_back("PlaneImpactResponse:" + name);
        }
        return det;
    }

    Configuration transform_cfg(const Configuration& base, const std::string& anode_tn,
                                const std::vector<std::string>& pirs)
    {
        Configuration cfg = base;
        cfg["anode"] = anode_tn;
        cfg["pirs"] = Json::arrayValue;
        for (const auto& p : pirs) cfg["pirs"].append(p);
        cfg["fluctuate"] = false;
        cfg["tick"] = tick;
        cfg["start_time"] = 0.0;
        cfg["readout_time"] = nticks * tick;
        cfg["drift_speed"] = 1.6 * units::mm / units::us;
        return cfg;
    }

    // Sum of |charge| per plane index of one frame.
    std::map<int, double> plane_sums(const IFrame::pointer& frame, const IAnodePlane::pointer& anode)
    {
        std::map<int, double> sums;
        for (const auto& trace : *frame->traces()) {
            const int plane = anode->resolve(trace->channel()).index();
            for (float q : trace->charge()) sums[plane] += std::abs(q);
        }
        return sums;
    }

    IFrame::pointer run(IDepoFramer& framer, const IDepo::vector& depos)
    {
        IDepoSet::pointer in = std::make_shared<Aux::SimpleDepoSet>(0, depos);
        IFrame::pointer out;
        REQUIRE(framer(in, out));
        REQUIRE(out);
        return out;
    }

    IDepo::pointer make_depo(double y, double z)
    {
        // x is inside the sensitive volume; DepoTransform does not drift.
        return std::make_shared<Aux::SimpleDepo>(200 * units::us, Point(1 * units::m, y, z), 10000.0,
                                                 nullptr, 1 * units::mm, 1 * units::mm);
    }
}  // namespace

TEST_CASE("yzmap lookup and edge clamp")
{
    Persist::TempDir tdir("yzmap-%%%%-%%%%", true);
    // grid[z][y]: 2 z bins x 2 y bins of 1 cm starting at (y, z) = (0, 0)
    auto fname = write_map(tdir, "map.json", {{"A", {{0, {{1, 2}, {3, 4}}}}}});
    Gen::YZMap yzmap;
    yzmap.configure(yzmap_cfg(fname, 2, 2, 1 * units::cm, 1 * units::cm, 0, 0));

    const double cm = units::cm;
    CHECK(yzmap.value("A", 0, 0.5 * cm, 0.5 * cm) == 1);
    CHECK(yzmap.value("A", 0, 1.5 * cm, 0.5 * cm) == 2);
    CHECK(yzmap.value("A", 0, 0.5 * cm, 1.5 * cm) == 3);
    CHECK(yzmap.value("A", 0, 1.5 * cm, 1.5 * cm) == 4);

    // Outside the map: the nearest edge bin.
    CHECK(yzmap.value("A", 0, -5 * cm, -5 * cm) == 1);
    CHECK(yzmap.value("A", 0, 100 * cm, 0.5 * cm) == 2);
    CHECK(yzmap.value("A", 0, 0.5 * cm, 100 * cm) == 3);
    CHECK(yzmap.value("A", 0, 100 * cm, 100 * cm) == 4);
    // Exactly at the upper edge (bin index == nbins).
    CHECK(yzmap.value("A", 0, 2 * cm, 2 * cm) == 4);

    // Unknown anode or plane.
    CHECK(yzmap.value("B", 0, 0.5 * cm, 0.5 * cm) == 0);
    CHECK(yzmap.value("A", 1, 0.5 * cm, 0.5 * cm) == 0);
}

TEST_CASE("yzmap grid size must match the configured binning")
{
    Persist::TempDir tdir("yzmap-%%%%-%%%%", true);
    auto fname = write_map(tdir, "map.json", {{"A", {{0, {{1, 2}, {3, 4}}}}}});
    {
        Gen::YZMap yzmap;
        CHECK_THROWS_AS(yzmap.configure(yzmap_cfg(fname, 2, 3, 1, 1, 0, 0)), ValueError);
    }
    {
        Gen::YZMap yzmap;
        CHECK_THROWS_AS(yzmap.configure(yzmap_cfg(fname, 3, 2, 1, 1, 0, 0)), ValueError);
    }
    {
        Gen::YZMap yzmap;
        CHECK_THROWS_AS(yzmap.configure(yzmap_cfg(fname, 2, 2, 0, 1, 0, 0)), ValueError);
    }
}

TEST_CASE("depotransformyz scales each plane by its own map value")
{
    const auto& det = detector();
    auto anode = det.anode;
    const auto& anode_tn = det.anode_tn;
    const auto& pirs = det.pirs;

    Persist::TempDir tdir("yzmap-%%%%-%%%%", true);
    // Two z bins of 520 cm (MicroBooNE wires span z 0-1037 cm), one y bin.
    const std::map<int, std::vector<double>> val = {{0, {0.5, 1.5}}, {1, {1.0, 0.8}}, {2, {2.0, 1.25}}};
    std::map<int, std::vector<std::vector<double>>> grids;
    for (const auto& [plane, v] : val) grids[plane] = {{v[0]}, {v[1]}};
    make_yzmap("yztest", yzmap_cfg(write_map(tdir, "map.json", {{anode_tn, grids}}),
                                   2, 1, 520 * units::cm, 400 * units::cm, 0, 200 * units::cm));

    Gen::DepoTransform plain;
    plain.configure(transform_cfg(plain.default_configuration(), anode_tn, pirs));

    Gen::DepoTransformYZ scaled;
    auto cfg = transform_cfg(scaled.default_configuration(), anode_tn, pirs);
    cfg["yzmap"] = "YZMap:yztest";
    scaled.configure(cfg);

    for (int zbin = 0; zbin < 2; ++zbin) {
        const double z = (zbin == 0 ? 200 : 800) * units::cm;
        IDepo::vector depos{make_depo(10 * units::cm, z)};
        auto s0 = plane_sums(run(plain, depos), anode);
        auto s1 = plane_sums(run(scaled, depos), anode);
        for (int plane = 0; plane < 3; ++plane) {
            CAPTURE(zbin);
            CAPTURE(plane);
            REQUIRE(s0[plane] > 0);
            CHECK(s1[plane] / s0[plane] == doctest::Approx(val.at(plane)[zbin]).epsilon(1e-4));
        }
    }

    // A depo in each z bin at once: each plane sums both, scaled separately.
    IDepo::vector depos{make_depo(-50 * units::cm, 200 * units::cm), make_depo(50 * units::cm, 800 * units::cm)};
    auto a0 = plane_sums(run(plain, {depos[0]}), anode);
    auto b0 = plane_sums(run(plain, {depos[1]}), anode);
    auto s1 = plane_sums(run(scaled, depos), anode);
    for (int plane = 0; plane < 3; ++plane) {
        CAPTURE(plane);
        const double expect = val.at(plane)[0] * a0[plane] + val.at(plane)[1] * b0[plane];
        // Induction-plane |q| sums are not additive where the two signals overlap;
        // the depos are 6 m apart in z and 1 m in y, so they do not.
        CHECK(s1[plane] == doctest::Approx(expect).epsilon(1e-4));
    }
}

TEST_CASE("depotransformyz with a unit map equals DepoTransform")
{
    const auto& det = detector();
    auto anode = det.anode;
    const auto& anode_tn = det.anode_tn;
    const auto& pirs = det.pirs;

    Persist::TempDir tdir("yzmap-%%%%-%%%%", true);
    std::map<int, std::vector<std::vector<double>>> grids;
    for (int plane = 0; plane < 3; ++plane) grids[plane] = {{1.0}};
    make_yzmap("yzunit", yzmap_cfg(write_map(tdir, "map.json", {{anode_tn, grids}}),
                                   1, 1, 1100 * units::cm, 400 * units::cm, 0, 200 * units::cm));

    Gen::DepoTransform plain;
    plain.configure(transform_cfg(plain.default_configuration(), anode_tn, pirs));
    Gen::DepoTransformYZ scaled;
    auto cfg = transform_cfg(scaled.default_configuration(), anode_tn, pirs);
    cfg["yzmap"] = "YZMap:yzunit";
    scaled.configure(cfg);

    IDepo::vector depos{make_depo(0, 300 * units::cm)};
    auto f0 = run(plain, depos);
    auto f1 = run(scaled, depos);
    REQUIRE(f0->traces()->size() == f1->traces()->size());
    for (size_t ind = 0; ind < f0->traces()->size(); ++ind) {
        const auto& t0 = f0->traces()->at(ind);
        const auto& t1 = f1->traces()->at(ind);
        CHECK(t0->channel() == t1->channel());
        CHECK(t0->tbin() == t1->tbin());
        CHECK(t0->charge() == t1->charge());
    }
}

TEST_CASE("depotransformyz requires a map entry for every plane")
{
    const auto& det = detector();
    auto anode = det.anode;
    const auto& anode_tn = det.anode_tn;
    const auto& pirs = det.pirs;

    Persist::TempDir tdir("yzmap-%%%%-%%%%", true);
    // Planes 0 and 1 only.
    std::map<int, std::vector<std::vector<double>>> grids = {{0, {{1.0}}}, {1, {{1.0}}}};
    make_yzmap("yzpartial", yzmap_cfg(write_map(tdir, "map.json", {{anode_tn, grids}}),
                                      1, 1, 1100 * units::cm, 400 * units::cm, 0, 200 * units::cm));

    Gen::DepoTransformYZ scaled;
    auto cfg = transform_cfg(scaled.default_configuration(), anode_tn, pirs);

    cfg["yzmap"] = "";
    CHECK_THROWS_AS(scaled.configure(cfg), ValueError);

    cfg["yzmap"] = "YZMap:yzpartial";
    CHECK_THROWS_AS(scaled.configure(cfg), ValueError);  // no plane 2

    cfg["yzmap_key"] = "AnodePlane:999";
    CHECK_THROWS_AS(scaled.configure(cfg), ValueError);  // no such key
}
