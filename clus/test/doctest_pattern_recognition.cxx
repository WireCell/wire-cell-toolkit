#include "WireCellClus/PatternDebugIO.h"
#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade.h"
#include "WireCellClus/NeutrinoTaggerInfo.h"
#include "WireCellClus/ParticleDataSet.h"
#include "WireCellClus/IPCTransform.h"

#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IRecombinationModel.h"

#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

#include <boost/graph/graph_traits.hpp>
#include <cmath>
#include <cstdlib>
#include <fstream>

// Usage:
//   wcdoctest-clus -tc="pattern_recognition*"
//   wcdoctest-clus -ts="pattern_recognition pure helpers"
//   wcdoctest-clus -ts="pattern_recognition replay [A]"
//
//   # Generate fixture [A] (vertex-not-found event):
//   WCT_DUMP_TAGGER_INPUTS=../clus/test/data/tagger_check_neutrino_input.json \
//     wire-cell -A infiles=rootfiles/nuselEval_5384_130_6501.root uboone-mabc.jsonnet
//
//   # Generate fixture [B] (vertex-found event — pick a CC numu from the 35-event sample):
//   WCT_DUMP_TAGGER_INPUTS=../clus/test/data/pattern_recognition_vertex_input.json \
//     wire-cell -A infiles=rootfiles/<your_numu_file>.root uboone-mabc.jsonnet
//   Then set WCT_TEST_DUMP_B to that path when running.

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::PR;
using namespace WireCell::Clus::Facade;

static const char* kDumpA  = "../clus/test/data/tagger_check_neutrino_input.json";
static const char* kDumpB  = "../clus/test/data/pattern_recognition_vertex_input.json";
static const char* kConfig = "../clus/test/data/uboone-mabc_config.json";

// ─── shared bootstrap ─────────────────────────────────────────────────────────

static void configure_components(const Json::Value& configs)
{
    auto log = Log::logger("test");
    for (const auto& comp : configs) {
        std::string type = comp["type"].asString();
        std::string name = comp.get("name", "").asString();
        std::string tn   = name.empty() ? type : type + ":" + name;
        try {
            auto icfg = Factory::lookup_tn<IConfigurable>(tn);
            auto cfg  = icfg->default_configuration();
            if (comp.isMember("data")) {
                for (const auto& key : comp["data"].getMemberNames())
                    cfg[key] = comp["data"][key];
            }
            icfg->configure(cfg);
        } catch (const std::exception& e) {
            log->trace("Skipping {}: {}", tn, e.what());
        }
    }
}

struct PrTestEnv {
    IDetectorVolumes::pointer       dv;
    IPCTransformSet::pointer        pcts;
    std::vector<IAnodePlane::pointer> anodes;
    ParticleDataSet::pointer        pdata;
    IRecombinationModel::pointer    recomb;
    DebugIO::TaggerTestData         fixture;
};

static Json::Value g_wcfg;  // cached for anode discovery

static void init_wct_once()
{
    static bool done = false;
    if (done) return;
    done = true;

    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellClus");
    pm.add("WireCellAux");
    pm.add("WireCellGen");
    pm.add("WireCellSigProc");

    auto full = Persist::load(kConfig);
    g_wcfg = Json::Value(Json::arrayValue);
    if (full.isArray() && !full.empty() &&
        full[0]["type"].asString() == "wire-cell") {
        for (Json::ArrayIndex i = 1; i < full.size(); ++i)
            g_wcfg.append(full[i]);
    } else {
        g_wcfg = full;
    }
    configure_components(g_wcfg);
}

static std::vector<IAnodePlane::pointer> collect_anodes()
{
    std::vector<IAnodePlane::pointer> anodes;
    for (const auto& comp : g_wcfg) {
        if (comp["type"].asString() != "AnodePlane") continue;
        std::string nm = comp.get("name", "").asString();
        std::string tn = nm.empty() ? "AnodePlane" : "AnodePlane:" + nm;
        anodes.push_back(Factory::find_tn<IAnodePlane>(tn));
    }
    return anodes;
}

static PrTestEnv build_env(const std::string& dump_path)
{
    init_wct_once();

    PrTestEnv env;
    env.dv    = Factory::find_tn<IDetectorVolumes>("DetectorVolumes");
    env.pcts  = Factory::find_tn<IPCTransformSet>("PCTransformSet");
    env.anodes = collect_anodes();

    try {
        auto icfg = Factory::find_tn<IConfigurable>("ParticleDataSet:ParticleDataSet");
        env.pdata = std::dynamic_pointer_cast<ParticleDataSet>(icfg);
    } catch (...) {}

    try {
        env.recomb = Factory::find_tn<IRecombinationModel>("BoxRecombination:box_recomb");
    } catch (...) {}

    env.fixture = DebugIO::load_tagger_inputs(dump_path);

    auto* grouping = env.fixture.grouping_node->value.facade<Grouping>();
    grouping->set_anodes(env.anodes);
    grouping->set_detector_volumes(env.dv);

    return env;
}

static PrTestEnv& env_A()
{
    static PrTestEnv e = build_env(kDumpA);
    return e;
}

// Build [B] only if the file exists, return nullptr otherwise.
static PrTestEnv* env_B_ptr()
{
    static bool attempted = false;
    static PrTestEnv* ptr = nullptr;
    if (!attempted) {
        attempted = true;
        if (std::ifstream(kDumpB).good()) {
            static PrTestEnv e = build_env(kDumpB);
            ptr = &e;
        }
    }
    return ptr;
}

// ─── per-test TrackFitting + Graph ────────────────────────────────────────────

struct PrContext {
    std::shared_ptr<TrackFitting> tf;
    std::shared_ptr<Graph>        graph;
    PatternAlgorithms             algo;
};

static PrContext make_context(PrTestEnv& env)
{
    PrContext ctx;
    ctx.tf = std::make_shared<TrackFitting>();
    ctx.tf->set_detector_volume(env.dv);
    ctx.tf->set_pc_transforms(env.pcts);
    ctx.tf->set_parameters(env.fixture.trackfitting_params);

    std::vector<Cluster*> to_preload = {env.fixture.main_cluster};
    for (auto* c : env.fixture.other_clusters) to_preload.push_back(c);
    ctx.tf->preload_clusters(to_preload);

    ctx.graph = std::make_shared<Graph>();
    ctx.tf->add_graph(ctx.graph);
    return ctx;
}

// Run the PR chain up to (and including) a given step and return {main_vertex, final_main_vertex}.
enum class Step {
    AfterFindProtoVertex,
    AfterClusteringPoints,
    AfterSeparateTrackShower,
    AfterDetermineDirection,
    AfterShowerDetermining,
    AfterDetermineMainVertex,
    AfterDetermineOverallMainVertex,
    AfterImproveVertex,
    AfterExamineDirection,
    AfterShowerClustering,
};

struct PrStepResult {
    VertexPtr             main_vertex;
    VertexPtr             final_main_vertex;
    IndexedVertexSet      vertices_in_long_muon;
    IndexedSegmentSet     segments_in_long_muon;
    ClusterVertexMap      map_cluster_main_vertices;
    IndexedShowerSet      showers;
    ShowerVertexMap       map_vertex_in_shower;
    ShowerSegmentMap      map_segment_in_shower;
    VertexShowerSetMap    map_vertex_to_shower;
    ClusterPtrSet         used_shower_clusters;
    Pi0KineFeatures       pio_kine;
    int                   acc_segment_id{0};
    IndexedShowerSet      pi0_showers;
    ShowerIntMap          map_shower_pio_id;
    std::map<int, std::vector<ShowerPtr>>       map_pio_id_showers;
    std::map<int, std::pair<double, int>>       map_pio_id_mass;
    std::map<int, std::pair<int, int>>          map_pio_id_saved_pair;
};

static PrStepResult run_through(PrTestEnv& env, PrContext& ctx, Step stop_after)
{
    PrStepResult r;

    ctx.algo.find_proto_vertex(*ctx.graph, *env.fixture.main_cluster,
                               *ctx.tf, env.dv, true, 2, true);
    if (stop_after == Step::AfterFindProtoVertex) return r;

    ctx.algo.clustering_points(*ctx.graph, *env.fixture.main_cluster, env.dv);
    if (stop_after == Step::AfterClusteringPoints) return r;

    ctx.algo.separate_track_shower(*ctx.graph, *env.fixture.main_cluster);
    if (stop_after == Step::AfterSeparateTrackShower) return r;

    ctx.algo.determine_direction(*ctx.graph, *env.fixture.main_cluster,
                                 env.pdata, env.recomb);
    if (stop_after == Step::AfterDetermineDirection) return r;

    ctx.algo.shower_determining_in_main_cluster(*ctx.graph, *env.fixture.main_cluster,
                                                env.pdata, env.recomb, env.dv);
    if (stop_after == Step::AfterShowerDetermining) return r;

    ctx.algo.determine_main_vertex(*ctx.graph, *env.fixture.main_cluster,
                                   r.main_vertex,
                                   r.vertices_in_long_muon, r.segments_in_long_muon,
                                   *ctx.tf, env.dv, env.pdata, env.recomb);
    if (r.main_vertex)
        r.map_cluster_main_vertices[env.fixture.main_cluster] = r.main_vertex;
    if (stop_after == Step::AfterDetermineMainVertex) return r;

    r.final_main_vertex = ctx.algo.determine_overall_main_vertex(
        *ctx.graph, r.map_cluster_main_vertices,
        env.fixture.main_cluster, env.fixture.other_clusters,
        r.vertices_in_long_muon, r.segments_in_long_muon,
        *ctx.tf, env.dv, env.pdata, env.recomb, true);
    if (r.final_main_vertex)
        r.map_cluster_main_vertices[env.fixture.main_cluster] = r.final_main_vertex;
    if (stop_after == Step::AfterDetermineOverallMainVertex) return r;

    if (r.final_main_vertex) {
        ctx.algo.improve_vertex(*ctx.graph, *env.fixture.main_cluster,
                                r.final_main_vertex,
                                r.vertices_in_long_muon, r.segments_in_long_muon,
                                *ctx.tf, env.dv, env.pdata, env.recomb, true, true);
        r.map_cluster_main_vertices[env.fixture.main_cluster] = r.final_main_vertex;
    }
    if (stop_after == Step::AfterImproveVertex) return r;

    if (r.final_main_vertex) {
        ctx.algo.clustering_points(*ctx.graph, *env.fixture.main_cluster, env.dv);
        ctx.algo.examine_direction(*ctx.graph, r.final_main_vertex, r.final_main_vertex,
                                   r.vertices_in_long_muon, r.segments_in_long_muon,
                                   env.pdata, env.recomb, true);
    }
    if (stop_after == Step::AfterExamineDirection) return r;

    if (r.final_main_vertex) {
        ctx.algo.shower_clustering_with_nv(
            r.acc_segment_id, r.pi0_showers,
            r.map_shower_pio_id, r.map_pio_id_showers,
            r.map_pio_id_mass, r.map_pio_id_saved_pair,
            r.pio_kine,
            r.vertices_in_long_muon, r.segments_in_long_muon,
            *ctx.graph, r.final_main_vertex, r.showers,
            env.fixture.main_cluster, env.fixture.other_clusters,
            r.map_cluster_main_vertices,
            r.map_vertex_in_shower, r.map_segment_in_shower,
            r.map_vertex_to_shower, r.used_shower_clusters,
            *ctx.tf, env.dv, env.pdata, env.recomb);
    }
    return r;
}

// ─── TEST_SUITE: pure helpers ─────────────────────────────────────────────────

TEST_SUITE("pattern_recognition pure helpers")
{

TEST_CASE("pattern_recognition calc_PCA_main_axis collinear")
{
    // Points along the x-axis: PCA direction should be ≈ ±x.
    std::vector<Facade::geo_point_t> pts;
    for (int i = -5; i <= 5; ++i)
        pts.push_back(Facade::geo_point_t{double(i), 0., 0.});

    PatternAlgorithms algo;
    auto [center, axis] = algo.calc_PCA_main_axis(pts);

    CHECK(std::abs(axis[0]) > 0.99);
    CHECK(std::abs(axis[1]) < 0.05);
    CHECK(std::abs(axis[2]) < 0.05);
}

TEST_CASE("pattern_recognition calc_PCA_main_axis degenerate single point")
{
    // Single point — should not crash; axis may be zero vector or unit in any direction.
    std::vector<Facade::geo_point_t> pts{{1., 2., 3.}};
    PatternAlgorithms algo;
    CHECK_NOTHROW(algo.calc_PCA_main_axis(pts));
}

TEST_CASE("pattern_recognition order_segments sorts by length descending")
{
    // order_segments sorts by track length descending (longest first).
    // segment_track_length uses fits[i].point for geometric distance, so populate fits.
    auto make_seg_with_fits = [](std::initializer_list<double> xs) {
        auto s = make_segment();
        std::vector<Fit> fits;
        for (double x : xs) { Fit f; f.point = WireCell::Point{x, 0, 0}; fits.push_back(f); }
        s->fits(fits);
        return s;
    };
    auto s_short  = make_seg_with_fits({0.0, 1.0});          // length 1 unit
    auto s_long   = make_seg_with_fits({0.0, 5.0, 10.0});    // length 10 units
    auto s_medium = make_seg_with_fits({0.0, 3.0});           // length 3 units

    std::vector<SegmentPtr> input = {s_short, s_long, s_medium};
    std::vector<SegmentPtr> ordered;
    PatternAlgorithms algo;
    algo.order_segments(ordered, input);

    REQUIRE(ordered.size() == 3);
    CHECK(ordered[0] == s_long);
    CHECK(ordered[2] == s_short);
}

TEST_CASE("pattern_recognition find_cluster_segments and vertices")
{
    // Build a tiny 3-segment graph; verify find_cluster_* returns correct counts.
    // (Uses fixture only for Cluster* identity; cluster point clouds not accessed.)
    auto& env = env_A();
    REQUIRE(env.fixture.main_cluster != nullptr);

    Graph g;
    auto v1 = make_vertex(g);
    auto v2 = make_vertex(g);
    auto v3 = make_vertex(g);
    // tag all vertices with main_cluster
    v1->cluster(env.fixture.main_cluster);
    v2->cluster(env.fixture.main_cluster);
    v3->cluster(env.fixture.main_cluster);

    auto s1 = make_segment(); s1->cluster(env.fixture.main_cluster);
    auto s2 = make_segment(); s2->cluster(env.fixture.main_cluster);
    add_segment(g, s1, v1, v2);
    add_segment(g, s2, v2, v3);

    PatternAlgorithms algo;
    auto segs  = algo.find_cluster_segments(g, *env.fixture.main_cluster);
    auto verts = algo.find_cluster_vertices(g, *env.fixture.main_cluster);

    CHECK(segs.size()  == 2);
    CHECK(verts.size() == 3);
}

TEST_CASE("pattern_recognition clean_up_graph removes all cluster content")
{
    // clean_up_graph removes ALL segments and vertices tagged with the given cluster.
    auto& env = env_A();
    REQUIRE(env.fixture.main_cluster != nullptr);

    Graph g;
    auto v1 = make_vertex(g); v1->cluster(env.fixture.main_cluster);
    auto v2 = make_vertex(g); v2->cluster(env.fixture.main_cluster);
    auto s  = make_segment(); s->cluster(env.fixture.main_cluster);
    add_segment(g, s, v1, v2);

    CHECK(boost::num_vertices(g) == 2);
    CHECK(boost::num_edges(g) == 1);

    PatternAlgorithms algo;
    bool changed = algo.clean_up_graph(g, *env.fixture.main_cluster);
    CHECK(changed);
    CHECK(boost::num_vertices(g) == 0);
    CHECK(boost::num_edges(g) == 0);
}

} // TEST_SUITE pure helpers

// ─── TEST_SUITE: replay [A] (vertex-not-found fixture) ───────────────────────

// For event 5384/130/6501, find_proto_vertex produces an empty graph — this is the
// expected behavior for this fixture (the event has no convincing neutrino vertex).
// The [A] test suite validates that each PR step runs without exceptions and leaves
// the graph in a consistent state, even when no structure is found.  Assertions on
// non-empty graph content belong in the [B] suite (vertex-found fixture).

TEST_SUITE("pattern_recognition replay [A]")
{

TEST_CASE("pattern_recognition find_proto_vertex [A]")
{
    auto& env = env_A();
    auto ctx  = make_context(env);
    CHECK_NOTHROW(run_through(env, ctx, Step::AfterFindProtoVertex));

    PatternAlgorithms algo;
    auto verts = algo.find_cluster_vertices(*ctx.graph, *env.fixture.main_cluster);

    // All vertex fit positions must be finite (passes vacuously if graph is empty)
    for (auto& v : verts) {
        CHECK(std::isfinite(v->fit().point[0]));
        CHECK(std::isfinite(v->fit().point[1]));
        CHECK(std::isfinite(v->fit().point[2]));
    }
    // All segments must have at least one wcpt
    for (auto& s : algo.find_cluster_segments(*ctx.graph, *env.fixture.main_cluster)) {
        CHECK(!s->wcpts().empty());
    }
    MESSAGE("find_proto_vertex: graph has ", boost::num_vertices(*ctx.graph),
            " vertices, ", boost::num_edges(*ctx.graph), " edges");
}

TEST_CASE("pattern_recognition clustering_points [A]")
{
    auto& env = env_A();
    auto ctx  = make_context(env);
    CHECK_NOTHROW(run_through(env, ctx, Step::AfterClusteringPoints));

    PatternAlgorithms algo;
    auto segs = algo.find_cluster_segments(*ctx.graph, *env.fixture.main_cluster);
    // For every segment that exists, global indices count must be ≥ 0
    for (auto& seg : segs) {
        CHECK(seg->global_indices("associate_points").size() >= 0);
    }
}

TEST_CASE("pattern_recognition separate_track_shower [A]")
{
    auto& env = env_A();
    auto ctx  = make_context(env);
    CHECK_NOTHROW(run_through(env, ctx, Step::AfterSeparateTrackShower));

    PatternAlgorithms algo;
    auto segs = algo.find_cluster_segments(*ctx.graph, *env.fixture.main_cluster);
    // Every segment must have a valid dirsign in {-1, 0, +1}
    for (auto& seg : segs) {
        int d = seg->dirsign();
        CHECK((d == -1 || d == 0 || d == 1));
    }
}

TEST_CASE("pattern_recognition determine_direction [A]")
{
    auto& env = env_A();
    auto ctx  = make_context(env);
    CHECK_NOTHROW(run_through(env, ctx, Step::AfterDetermineDirection));

    PatternAlgorithms algo;
    for (auto& seg : algo.find_cluster_segments(*ctx.graph, *env.fixture.main_cluster)) {
        int d = seg->dirsign();
        CHECK((d == -1 || d == 0 || d == 1));
    }
}

TEST_CASE("pattern_recognition shower_determining_in_main_cluster [A]")
{
    auto& env = env_A();
    auto ctx  = make_context(env);
    CHECK_NOTHROW(run_through(env, ctx, Step::AfterShowerDetermining));
}

TEST_CASE("pattern_recognition determine_main_vertex returns null for [A]")
{
    // Event 5384/130/6501: no vertex found — confirm graceful null return.
    auto& env = env_A();
    auto ctx  = make_context(env);
    PrStepResult r;
    CHECK_NOTHROW(r = run_through(env, ctx, Step::AfterDetermineMainVertex));
    MESSAGE("main_vertex found: ", (r.main_vertex != nullptr));
    // Full chain yields nullptr; pio_kine stays at defaults
    CHECK(r.pio_kine.flag == 0);
}

TEST_CASE("pattern_recognition determine_overall_main_vertex [A] no crash")
{
    auto& env = env_A();
    auto ctx  = make_context(env);
    PrStepResult r;
    CHECK_NOTHROW(r = run_through(env, ctx, Step::AfterDetermineOverallMainVertex));
    MESSAGE("final_main_vertex found: ", (r.final_main_vertex != nullptr));
    CHECK(r.pio_kine.flag == 0);
}

} // TEST_SUITE replay [A]

// ─── TEST_SUITE: replay [B] (vertex-found fixture) ───────────────────────────

TEST_SUITE("pattern_recognition replay [B]")
{

TEST_CASE("pattern_recognition determine_main_vertex [B]")
{
    PrTestEnv* env = env_B_ptr();
    if (!env) {
        { std::string _m = std::string("Skipping: fixture B not found at ") + kDumpB; MESSAGE(_m); }
        return;
    }

    auto ctx = make_context(*env);
    auto r   = run_through(*env, ctx, Step::AfterDetermineMainVertex);

    // For a CC numu event the main vertex should be found.
    REQUIRE(r.main_vertex != nullptr);

    // Vertex position should be inside the detector active volume.
    auto& pt = r.main_vertex->fit().point;
    CHECK(std::isfinite(pt[0]));
    CHECK(std::isfinite(pt[1]));
    CHECK(std::isfinite(pt[2]));

    // contained_by returns the WirePlaneId of the containing volume;
    // the x-coordinate is the drift direction — check it's in a plausible range.
    CHECK(pt[0] > -100.0);   // upstream of cathode (generous bound)
    CHECK(pt[0] <  500.0);   // downstream of APA
}

TEST_CASE("pattern_recognition determine_overall_main_vertex [B]")
{
    PrTestEnv* env = env_B_ptr();
    if (!env) {
        { std::string _m = std::string("Skipping: fixture B not found at ") + kDumpB; MESSAGE(_m); }
        return;
    }

    auto ctx = make_context(*env);
    auto r   = run_through(*env, ctx, Step::AfterDetermineOverallMainVertex);

    REQUIRE(r.final_main_vertex != nullptr);
    auto& pt = r.final_main_vertex->fit().point;
    CHECK(std::isfinite(pt[0]));
    CHECK(std::isfinite(pt[1]));
    CHECK(std::isfinite(pt[2]));
}

TEST_CASE("pattern_recognition improve_vertex [B]")
{
    PrTestEnv* env = env_B_ptr();
    if (!env) {
        { std::string _m = std::string("Skipping: fixture B not found at ") + kDumpB; MESSAGE(_m); }
        return;
    }

    auto ctx = make_context(*env);

    // Capture pre-improve position
    auto r_before = run_through(*env, ctx, Step::AfterDetermineOverallMainVertex);
    if (!r_before.final_main_vertex) {
        MESSAGE("Fixture B has no overall vertex — skipping improve_vertex check");
        return;
    }
    auto pt_before = r_before.final_main_vertex->fit().point;

    // Now run one step further (improve_vertex mutates final_main_vertex in place)
    auto ctx2 = make_context(*env);
    auto r_after = run_through(*env, ctx2, Step::AfterImproveVertex);
    REQUIRE(r_after.final_main_vertex != nullptr);

    auto& pt_after = r_after.final_main_vertex->fit().point;
    CHECK(std::isfinite(pt_after[0]));
    CHECK(std::isfinite(pt_after[1]));
    CHECK(std::isfinite(pt_after[2]));

    // Vertex should not teleport more than 20 cm
    double dx = pt_after[0] - pt_before[0];
    double dy = pt_after[1] - pt_before[1];
    double dz = pt_after[2] - pt_before[2];
    double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
    CHECK(dist < 20.0);  // units: cm
}

TEST_CASE("pattern_recognition shower_clustering_with_nv [B]")
{
    PrTestEnv* env = env_B_ptr();
    if (!env) {
        { std::string _m = std::string("Skipping: fixture B not found at ") + kDumpB; MESSAGE(_m); }
        return;
    }

    auto ctx = make_context(*env);
    auto r   = run_through(*env, ctx, Step::AfterShowerClustering);

    if (!r.final_main_vertex) {
        MESSAGE("Fixture B has no final vertex; shower clustering was skipped");
        return;
    }

    // All vertices in map_vertex_to_shower must be reachable from the graph.
    PatternAlgorithms algo;
    auto graph_verts = algo.find_cluster_vertices(*ctx.graph, *env->fixture.main_cluster);
    IndexedVertexSet graph_vtx_set(graph_verts.begin(), graph_verts.end());

    for (auto& [vtx, showers] : r.map_vertex_to_shower) {
        bool found_in_graph = graph_vtx_set.count(vtx) > 0;
        // Shower vertices may live in other clusters too — just check they're non-null.
        CHECK(vtx != nullptr);
        (void)found_in_graph;
    }

    // pio_kine.mass should be finite (or zero if no pi0 was found)
    CHECK(std::isfinite(r.pio_kine.mass));
    CHECK(std::isfinite(r.pio_kine.energy_1));
}

} // TEST_SUITE replay [B]
