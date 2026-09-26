// doc sbnd_xin/pr/149 round 2 -- ClusteringResampleLive's pure helpers
// (clus/inc/WireCellClus/ResampleLive.h).
//
// The visitor rebuilds, from a saved point-cloud tree, the blob shape and slice
// activity a BlobSampler reads.  The runtime identity gate (a 'stepped' resample
// reproduces the clustering job's saved cloud bit for bit) proves that on real
// events; these cases pin the two rebuilds in isolation so a regression names its
// piece:
//   1. the strip layout (two dummy layers, then U/V/W at 2/3/4);
//   2. a blob rebuilt from a tiled blob's wire bounds samples EXACTLY like the
//      tiled blob, for stepped and for charge_stepped with disable_mix_dead_cell
//      true (the prototype PR job's setting) -- every array, bit for bit;
//   3. the activity precedence: live beats dead, dead-only is (0, 1e12), and a
//      wire neither live nor dead stays absent (the negative control shows why:
//      writing absent wires as (0, 1e12) flips charge_stepped's dead-plane test);
//   4. the scalar refresh touches only the four point-dependent arrays.

#include "WireCellClus/ResampleLive.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IBlobSampler.h"
#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleSlice.h"
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

    // PDHD apa1 as its production clustering config builds it (the same fixture
    // as doctest_blob_sampler_wrapped_channel.cxx, under this file's own names).
    IAnodePlane::pointer pr149r2_anode()
    {
        static IAnodePlane::pointer anode;
        if (anode) return anode;
        PluginManager& pm = PluginManager::instance();
        pm.add("WireCellAux");
        pm.add("WireCellGen");
        pm.add("WireCellClus");
        {
            auto icfg = Factory::lookup<IConfigurable>("WireSchemaFile", "pr149r2_pdhd_wires");
            auto cfg = icfg->default_configuration();
            cfg["filename"] = "protodunehd-wires-larsoft-v1.json.bz2";
            icfg->configure(cfg);
        }
        auto icfg = Factory::lookup<IConfigurable>("AnodePlane", "pr149r2_pdhd_apa1");
        auto cfg = icfg->default_configuration();
        cfg["ident"] = 1;
        cfg["nimpacts"] = 10;
        cfg["wire_schema"] = "WireSchemaFile:pr149r2_pdhd_wires";
        cfg["faces"][0] = Json::nullValue;
        cfg["faces"][1]["anode"] = 3520.945;
        cfg["faces"][1]["response"] = 3430.465;
        cfg["faces"][1]["cathode"] = 1.5875;
        icfg->configure(cfg);
        anode = Factory::find<IAnodePlane>("AnodePlane", "pr149r2_pdhd_apa1");
        return anode;
    }

    // The LAST non-null face, as doctest_blob_sampler_wrapped_channel.cxx picks it:
    // the null-configured faces[0] is not reliably a null pointer.
    IAnodeFace::pointer live_face(const IAnodePlane::pointer& anode)
    {
        IAnodeFace::pointer face;
        for (const auto& f : anode->faces()) {
            if (f) face = f;
        }
        return face;
    }

    // The clustering job's extra list (sbnd/clus.jsonnet bs_live_face).
    IBlobSampler::pointer make_sampler(const std::string& name, const std::string& strategy)
    {
        auto icfg = Factory::lookup<IConfigurable>("BlobSampler", name);
        auto cfg = icfg->default_configuration();
        Json::Value one(Json::objectValue);
        one["name"] = strategy;
        if (strategy == "charge_stepped") {
            one["disable_mix_dead_cell"] = true;   // the prototype PR job's resample
        }
        cfg["strategy"] = Json::Value(Json::arrayValue);
        cfg["strategy"].append(one);
        cfg["extra"] = Json::Value(Json::arrayValue);
        for (const char* e : {".*wire_index", ".*charge_val", ".*charge_unc", "wpid"}) {
            cfg["extra"].append(e);
        }
        icfg->configure(cfg);
        return Factory::find<IBlobSampler>("BlobSampler", name);
    }

    // Every array of `a` equals the same-named array of `b`, element for element.
    void require_identical(const PointCloud::Dataset& a, const PointCloud::Dataset& b)
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
    RayGrid::Blob tiled_blob(const IAnodeFace::pointer& face)
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

    RL::plane_bounds_t bounds_of(const RayGrid::Blob& shape)
    {
        const auto& s = shape.strips();
        REQUIRE(s.size() == 5);
        return {std::make_pair((int) s[2].bounds.first, (int) s[2].bounds.second),
                std::make_pair((int) s[3].bounds.first, (int) s[3].bounds.second),
                std::make_pair((int) s[4].bounds.first, (int) s[4].bounds.second)};
    }
}

TEST_CASE("pr149 r2: strips_from_bounds lays out two dummy layers then U/V/W")
{
    const RL::plane_bounds_t b = {std::make_pair(10, 14), std::make_pair(20, 21), std::make_pair(30, 39)};
    const auto strips = RL::strips_from_bounds(b);
    REQUIRE(strips.size() == 5);
    for (int l = 0; l < 5; ++l) {
        CHECK(strips[l].layer == l);
    }
    CHECK(strips[0].bounds == std::make_pair(0, 1));
    CHECK(strips[1].bounds == std::make_pair(0, 1));
    CHECK(strips[2].bounds == std::make_pair(10, 14));
    CHECK(strips[3].bounds == std::make_pair(20, 21));
    CHECK(strips[4].bounds == std::make_pair(30, 39));
}

TEST_CASE("pr149 r2: a blob rebuilt from its wire bounds samples exactly like the tiled blob")
{
    auto anode = pr149r2_anode();
    REQUIRE(anode);
    auto face = live_face(anode);
    REQUIRE(face);
    const auto planes = face->planes();

    const RayGrid::Blob tiled = tiled_blob(face);
    const RL::plane_bounds_t b = bounds_of(tiled);
    const RayGrid::Blob rebuilt = RL::shape_from_bounds(face->raygrid(), b);
    REQUIRE(bounds_of(rebuilt) == b);

    // A live slice: a mix of charges above and below charge_stepped's 4000 cut,
    // so its all-wires branch keeps some non-stepped wires and drops others.
    ISlice::map_t activity;
    for (int p = 0; p < 3; ++p) {
        const auto& wires = planes[p]->wires();
        for (int wi = b[p].first - 2; wi < b[p].second + 2; ++wi) {
            auto ich = anode->channel(wires[wi]->channel());
            REQUIRE(ich);
            activity[ich] = ISlice::value_t((wi % 3 == 0) ? 2500.0 : 6000.0 + wi, 40.0);
        }
    }
    auto slice = std::make_shared<Aux::SimpleSlice>(nullptr, 7, 7 * 2 * units::ms, 2 * units::ms, activity);
    IBlob::pointer iblob_tiled = std::make_shared<Aux::SimpleBlob>(0, 1.0f, 0.0f, tiled, slice, face);
    IBlob::pointer iblob_rebuilt = std::make_shared<Aux::SimpleBlob>(0, 0.0f, 0.0f, rebuilt, slice, face);

    for (const std::string strategy : {"stepped", "charge_stepped"}) {
        auto bs = make_sampler("pr149r2_" + strategy, strategy);
        auto [ds_t, aux_t] = bs->sample_blob(iblob_tiled, 0);
        auto [ds_r, aux_r] = bs->sample_blob(iblob_rebuilt, 0);
        MESSAGE(strategy << ": " << ds_t.size_major() << " points from the tiled blob, "
                << ds_r.size_major() << " from the rebuilt one");
        REQUIRE(ds_t.size_major() > 0);
        require_identical(ds_t, ds_r);
        require_identical(aux_t, aux_r);
    }
}

TEST_CASE("pr149 r2: compose_activity -- live beats dead, dead is (0, 1e12), absent stays absent")
{
    // wires 0..6: 1 and 2 live, 2 and 3 in the dead-wind registry, the rest neither
    std::map<int, std::pair<double, double>> live_row = {{1, {5000.0, 50.0}}, {2, {3000.0, 30.0}}};
    auto live = [&](int w, double& q, double& e) {
        auto it = live_row.find(w);
        if (it == live_row.end()) return false;
        q = it->second.first;
        e = it->second.second;
        return true;
    };
    auto dead = [](int w) { return w == 2 || w == 3; };

    const auto out = RL::compose_activity(0, 6, live, dead);
    REQUIRE(out.size() == 3);
    CHECK(out[0].wire == 1);
    CHECK(out[0].charge == 5000.0);
    CHECK(out[0].error == 50.0);
    CHECK(out[1].wire == 2);   // live AND dead-listed: live wins
    CHECK(out[1].charge == 3000.0);
    CHECK(out[1].error == 30.0);
    CHECK(out[2].wire == 3);   // dead only
    CHECK(out[2].charge == RL::dead_charge);
    CHECK(out[2].error == RL::dead_error);
    CHECK(RL::compose_activity(5, 4, live, dead).empty());   // empty range
}

TEST_CASE("pr149 r2: an absent wire must not be written dead (charge_stepped's dead-plane test)")
{
    auto anode = pr149r2_anode();
    auto face = live_face(anode);
    REQUIRE(face);
    const auto planes = face->planes();
    const RayGrid::Blob shape = tiled_blob(face);
    const RL::plane_bounds_t b = bounds_of(shape);

    // The widest plane is charge_stepped's "max" view.  Every wire of the blob
    // is live at 6000 except, in that plane, the FIRST strip wire, which is
    // absent (a live channel with no signal), and the SECOND, which is live at
    // 2000 -- below the 4000 cut and not a stepped ("must") wire.
    int pmax = 0;
    for (int p = 1; p < 3; ++p) {
        if (b[p].second - b[p].first > b[pmax].second - b[pmax].first) pmax = p;
    }
    REQUIRE(b[pmax].second - b[pmax].first >= 4);
    auto build = [&](bool absent_as_dead) {
        ISlice::map_t activity;
        for (int p = 0; p < 3; ++p) {
            const auto& wires = planes[p]->wires();
            for (int wi = b[p].first; wi < b[p].second; ++wi) {
                auto ich = anode->channel(wires[wi]->channel());
                REQUIRE(ich);
                double q = 6000.0;
                if (p == pmax && wi == b[p].first) {
                    if (absent_as_dead) activity[ich] = ISlice::value_t(RL::dead_charge, RL::dead_error);
                    continue;
                }
                if (p == pmax && wi == b[p].first + 1) q = 2000.0;
                activity[ich] = ISlice::value_t(q, 40.0);
            }
        }
        auto slice = std::make_shared<Aux::SimpleSlice>(nullptr, 0, 0.0, 2 * units::ms, activity);
        return (IBlob::pointer) std::make_shared<Aux::SimpleBlob>(0, 0.0f, 0.0f, shape, slice, face);
    };
    auto bs = make_sampler("pr149r2_cs_absent", "charge_stepped");
    auto [ds_absent, aux_absent] = bs->sample_blob(build(false), 0);
    auto [ds_dead, aux_dead] = bs->sample_blob(build(true), 0);
    MESSAGE("charge_stepped points: absent wire left out " << ds_absent.size_major()
            << ", absent wire written dead " << ds_dead.size_major());
    REQUIRE(ds_absent.size_major() > 0);
    // Written dead, the boundary wire marks the max plane bad and drops its
    // threshold to 0, so the 2000 wire is kept; left absent, it is cut.  The
    // same slice read two ways samples differently -- which is why the visitor
    // never writes an absent wire as dead.
    CHECK(ds_dead.size_major() > ds_absent.size_major());
}

TEST_CASE("pr149 r2: refresh_scalar replaces only the point-dependent arrays")
{
    PointCloud::Dataset scalar;
    scalar.add("charge", PointCloud::Array({12.5}));
    scalar.add("center_x", PointCloud::Array({1.0}));
    scalar.add("center_y", PointCloud::Array({2.0}));
    scalar.add("center_z", PointCloud::Array({3.0}));
    scalar.add("npoints", PointCloud::Array({4}));
    scalar.add("max_wire_interval", PointCloud::Array({3}));

    PointCloud::Dataset fresh;
    fresh.add("charge", PointCloud::Array({0.0}));
    fresh.add("center_x", PointCloud::Array({10.0}));
    fresh.add("center_y", PointCloud::Array({20.0}));
    fresh.add("center_z", PointCloud::Array({30.0}));
    fresh.add("npoints", PointCloud::Array({40}));
    fresh.add("max_wire_interval", PointCloud::Array({1}));

    RL::refresh_scalar(scalar, fresh);
    CHECK(scalar.get("charge")->elements<double>()[0] == 12.5);   // blob property kept
    CHECK(scalar.get("max_wire_interval")->elements<int>()[0] == 3);
    CHECK(scalar.get("center_x")->elements<double>()[0] == 10.0);
    CHECK(scalar.get("center_y")->elements<double>()[0] == 20.0);
    CHECK(scalar.get("center_z")->elements<double>()[0] == 30.0);
    CHECK(scalar.get("npoints")->elements<int>()[0] == 40);
    CHECK(scalar.keys() == fresh.keys());
}
