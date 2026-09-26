// doc pdvd/113 -- ImproveCluster_2's retile_mode knob.
//
// The Steiner stage's retiler fabricates activity (dead channels and good charge
// within 20 cm, a +-3-wire x +-3-slice disc painted along two Dijkstra paths)
// before re-tiling the cluster.  retile_mode removes that fabrication step by
// step: "no_paint", "footprint", "none" (no tiling: each original blob is
// re-sampled in place).  These cases pin:
//   1. the C++ default is "full", the historical retile, so a config without the
//      key is byte-identical;
//   2. an unknown mode fails on its own message (validated before the base
//      configure, which needs live detector volumes);
//   3. "none" rebuilds each blob from its wire bounds (ResampleLive helpers) and
//      samples it with the PRODUCTION PDHD/PDVD charge_stepped setting,
//      disable_mix_dead_cell = false: the rebuilt blob samples exactly like the
//      tiled blob (doctest_resample_live.cxx pins only the `true` setting);
//   4. why the retile's forced cells turn into points: under
//      disable_mix_dead_cell = false a non-stepped wire at charge EXACTLY 0 (the
//      retile's sentinel) is kept, while a live wire below the 4000 cut is not.

#include "WireCellClus/ResampleLive.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IBlobSampler.h"
#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/RayHelpers.h"
#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/doctest.h"

#include <map>
#include <string>
#include <vector>

using namespace WireCell;
namespace RL = WireCell::Clus::ResampleLive;

namespace {

    // PDHD apa1 as its production clustering config builds it (the fixture of
    // doctest_resample_live.cxx, under this file's own component names).
    IAnodePlane::pointer d113_anode()
    {
        static IAnodePlane::pointer anode;
        if (anode) return anode;
        PluginManager& pm = PluginManager::instance();
        pm.add("WireCellAux");
        pm.add("WireCellGen");
        pm.add("WireCellClus");
        {
            auto icfg = Factory::lookup<IConfigurable>("WireSchemaFile", "d113_pdhd_wires");
            auto cfg = icfg->default_configuration();
            cfg["filename"] = "protodunehd-wires-larsoft-v1.json.bz2";
            icfg->configure(cfg);
        }
        auto icfg = Factory::lookup<IConfigurable>("AnodePlane", "d113_pdhd_apa1");
        auto cfg = icfg->default_configuration();
        cfg["ident"] = 1;
        cfg["nimpacts"] = 10;
        cfg["wire_schema"] = "WireSchemaFile:d113_pdhd_wires";
        cfg["faces"][0] = Json::nullValue;
        cfg["faces"][1]["anode"] = 3520.945;
        cfg["faces"][1]["response"] = 3430.465;
        cfg["faces"][1]["cathode"] = 1.5875;
        icfg->configure(cfg);
        anode = Factory::find<IAnodePlane>("AnodePlane", "d113_pdhd_apa1");
        return anode;
    }

    IAnodeFace::pointer d113_face(const IAnodePlane::pointer& anode)
    {
        IAnodeFace::pointer face;
        for (const auto& f : anode->faces()) {
            if (f) face = f;
        }
        return face;
    }

    // The PR job's retile sampler as pdhd/protodunevd clus.jsonnet bs_live_face
    // builds it for strategy_name 'charge_stepped'.
    IBlobSampler::pointer d113_sampler(const std::string& name)
    {
        auto icfg = Factory::lookup<IConfigurable>("BlobSampler", name);
        auto cfg = icfg->default_configuration();
        Json::Value one(Json::objectValue);
        one["name"] = "charge_stepped";
        one["disable_mix_dead_cell"] = false;
        cfg["strategy"] = Json::Value(Json::arrayValue);
        cfg["strategy"].append(one);
        cfg["extra"] = Json::Value(Json::arrayValue);
        for (const char* e : {".*wire_index", ".*charge_val", ".*charge_unc", "wpid"}) {
            cfg["extra"].append(e);
        }
        icfg->configure(cfg);
        return Factory::find<IBlobSampler>("BlobSampler", name);
    }

    void d113_identical(const PointCloud::Dataset& a, const PointCloud::Dataset& b)
    {
        REQUIRE(a.keys() == b.keys());
        REQUIRE(a.size_major() == b.size_major());
        for (const auto& key : a.keys()) {
            auto aa = a.get(key);
            auto bb = b.get(key);
            REQUIRE(aa->dtype() == bb->dtype());
            const auto ab = aa->bytes();
            const auto bbs = bb->bytes();
            CHECK_MESSAGE(std::vector<std::byte>(ab.begin(), ab.end()) ==
                              std::vector<std::byte>(bbs.begin(), bbs.end()),
                          "array " << key << " differs");
        }
    }

    // The largest blob tiled from a dense 4 cm x 4 cm patch of apa1.
    RayGrid::Blob d113_tiled_blob(const IAnodeFace::pointer& face)
    {
        const auto& coords = face->raygrid();
        std::vector<Point> pts;
        for (double y = 2840; y <= 2880; y += 1.0) {
            for (double z = 1132; z <= 1172; z += 1.0) {
                pts.emplace_back(3520.945 * units::mm, y * units::mm, z * units::mm);
            }
        }
        auto measures = RayGrid::make_measures(coords, pts);
        auto activities = RayGrid::make_activities(coords, measures);
        auto blobs = RayGrid::make_blobs(coords, activities);
        REQUIRE(!blobs.empty());
        auto width = [](const RayGrid::Strip& s) { return s.bounds.second - s.bounds.first; };
        size_t best = 0;
        for (size_t i = 1; i < blobs.size(); ++i) {
            const auto& a = blobs[i].strips();
            const auto& b = blobs[best].strips();
            if (width(a[2]) * width(a[4]) > width(b[2]) * width(b[4])) best = i;
        }
        return blobs[best];
    }

    RL::plane_bounds_t d113_bounds(const RayGrid::Blob& shape)
    {
        const auto& s = shape.strips();
        REQUIRE(s.size() == 5);
        return {std::make_pair((int) s[2].bounds.first, (int) s[2].bounds.second),
                std::make_pair((int) s[3].bounds.first, (int) s[3].bounds.second),
                std::make_pair((int) s[4].bounds.first, (int) s[4].bounds.second)};
    }
}

TEST_CASE("pdvd doc113: retile_mode defaults to the full retile")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("ImproveCluster_2", "d113_retile_mode_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE_MESSAGE(cfg.isMember("retile_mode"), "missing knob: retile_mode");
    CHECK(cfg["retile_mode"].asString() == "full");
}

TEST_CASE("pdvd doc113: an unknown retile_mode is refused by name")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("ImproveCluster_2", "d113_retile_mode_bad");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    cfg["retile_mode"] = "nopaint";   // a plausible typo of "no_paint"
    std::string msg;
    try {
        icfg->configure(cfg);
    }
    catch (const ValueError& e) {
        msg = errstr(e);
    }
    CHECK_MESSAGE(msg.find("unknown retile_mode 'nopaint'") != std::string::npos, "got: " << msg);
}

TEST_CASE("pdvd doc113: retile_mode none -- a rebuilt blob samples like the tiled blob under production charge_stepped")
{
    auto anode = d113_anode();
    REQUIRE(anode);
    auto face = d113_face(anode);
    REQUIRE(face);
    const auto planes = face->planes();

    const RayGrid::Blob tiled = d113_tiled_blob(face);
    const RL::plane_bounds_t b = d113_bounds(tiled);
    const RayGrid::Blob rebuilt = RL::shape_from_bounds(face->raygrid(), b);
    REQUIRE(d113_bounds(rebuilt) == b);

    // Charges above and below the 4000 cut, plus one dead (0, 1e12) wire per
    // plane inside the bounds, so all three admission branches are exercised.
    ISlice::map_t activity;
    for (int p = 0; p < 3; ++p) {
        const auto& wires = planes[p]->wires();
        for (int wi = b[p].first - 2; wi < b[p].second + 2; ++wi) {
            auto ich = anode->channel(wires[wi]->channel());
            REQUIRE(ich);
            if (wi == b[p].first + 1) {
                activity[ich] = ISlice::value_t(RL::dead_charge, RL::dead_error);
                continue;
            }
            activity[ich] = ISlice::value_t((wi % 3 == 0) ? 2500.0 : 6000.0 + wi, 40.0);
        }
    }
    auto slice = std::make_shared<Aux::SimpleSlice>(nullptr, 7, 7 * 2 * units::ms, 2 * units::ms, activity);
    IBlob::pointer iblob_tiled = std::make_shared<Aux::SimpleBlob>(0, 1.0f, 0.0f, tiled, slice, face);
    IBlob::pointer iblob_rebuilt = std::make_shared<Aux::SimpleBlob>(0, 0.0f, 0.0f, rebuilt, slice, face);

    auto bs = d113_sampler("d113_cs_prod");
    auto [ds_t, aux_t] = bs->sample_blob(iblob_tiled, 0);
    auto [ds_r, aux_r] = bs->sample_blob(iblob_rebuilt, 0);
    MESSAGE("charge_stepped (disable_mix_dead_cell=false): " << ds_t.size_major()
            << " points tiled, " << ds_r.size_major() << " rebuilt");
    REQUIRE(ds_t.size_major() > 0);
    d113_identical(ds_t, ds_r);
    d113_identical(aux_t, aux_r);
}

TEST_CASE("pdvd doc113: under disable_mix_dead_cell=false a zero-charge wire is admitted, a weak live wire is not")
{
    auto anode = d113_anode();
    auto face = d113_face(anode);
    REQUIRE(face);
    const auto planes = face->planes();
    const RayGrid::Blob shape = d113_tiled_blob(face);
    const RL::plane_bounds_t b = d113_bounds(shape);

    // In the widest plane (charge_stepped's "max" view) the second strip wire is
    // not a stepped ("must") wire.  Give it charge exactly 0 in one slice and a
    // live 2000 in the other; every other wire is live at 6000.
    int pmax = 0;
    for (int p = 1; p < 3; ++p) {
        if (b[p].second - b[p].first > b[pmax].second - b[pmax].first) pmax = p;
    }
    REQUIRE(b[pmax].second - b[pmax].first >= 4);
    auto build = [&](double q_probe) {
        ISlice::map_t activity;
        for (int p = 0; p < 3; ++p) {
            const auto& wires = planes[p]->wires();
            for (int wi = b[p].first; wi < b[p].second; ++wi) {
                auto ich = anode->channel(wires[wi]->channel());
                REQUIRE(ich);
                const bool probe = (p == pmax && wi == b[p].first + 1);
                activity[ich] = ISlice::value_t(probe ? q_probe : 6000.0, 40.0);
            }
        }
        auto slice = std::make_shared<Aux::SimpleSlice>(nullptr, 0, 0.0, 2 * units::ms, activity);
        return (IBlob::pointer) std::make_shared<Aux::SimpleBlob>(0, 0.0f, 0.0f, shape, slice, face);
    };
    auto bs = d113_sampler("d113_cs_zero");
    auto [ds_zero, aux_zero] = bs->sample_blob(build(0.0), 0);
    auto [ds_weak, aux_weak] = bs->sample_blob(build(2000.0), 0);
    MESSAGE("probe wire at 0: " << ds_zero.size_major() << " points; at 2000: " << ds_weak.size_major());
    REQUIRE(ds_weak.size_major() > 0);
    CHECK(ds_zero.size_major() > ds_weak.size_major());
}
