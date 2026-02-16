#include "WireCellClus/PatternDebugIO.h"
#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade.h"

#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellClus/IPCTransform.h"

#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

#include <cstdlib>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::PR;

// These are commands to save data and run tests ... 
// #Generate dump data - Run your full pipeline with:
// WCT_DUMP_INIT_FIRST_SEGMENT=./tmp/test_data.json wire-cell -l stderr -A kind=both  -A beezip=mabc_0.zip -A initial_index="0" -A initial_runNo="5384" -A initial_subRunNo="130" -A initial_eventNo="6501" -A infiles=rootfiles/nuselEval_5384_130_6501.root uboone-mabc.jsonnet
// # Run geometry-only test (no DV/PCTS needed):
// WCT_TEST_DUMP=./tmp/test_data.json wcdoctest-clus -tc="init_first_segment geometry only"
// # generate test configuration
// wcsonnet -o ./tmp/test_config.json uboone-mabc.jsonnet 
// # Run full end-to-end test (needs WCT config for DV/PCTS):
// WCT_TEST_DUMP=./tmp/test_data.json WCT_TEST_CONFIG=./tmp/test_config.json wcdoctest-clus -tc="init_first_segment end-to-end"




/// Helper: configure WCT components from a JSON config array.
/// The config should be an array of objects, each with "type", optional "name",
/// and "data" keys (standard WCT configuration format).
/// Components that can't be found are silently skipped, so you can pass the
/// full evaluated jsonnet config without filtering.
static void configure_components(const Json::Value& configs)
{
    auto log = Log::logger("test");
    for (const auto& comp : configs) {
        std::string type = comp["type"].asString();
        std::string name = comp.get("name", "").asString();
        std::string tn = type;
        if (!name.empty()) tn += ":" + name;

        try {
            auto icfg = Factory::lookup_tn<IConfigurable>(tn);
            auto cfg = icfg->default_configuration();

            // Overlay saved data onto defaults
            if (comp.isMember("data")) {
                for (const auto& key : comp["data"].getMemberNames()) {
                    cfg[key] = comp["data"][key];
                }
            }
            icfg->configure(cfg);
            log->debug("Configured {}", tn);
        } catch (const std::exception& e) {
            log->debug("Skipping {}: {}", tn, e.what());
        }
    }
}

TEST_CASE("init_first_segment geometry only")
{
    // This test validates the geometry/selection logic without needing
    // DetectorVolumes or TrackFitting.
    //
    // Usage:
    //   - Without env var: uses default data file in clus/test/data/
    //   - With env var (for development):
    //     WCT_TEST_DUMP=/path/to/dump.json wcdoctest-clus -tc="init_first_segment geometry only"

    // Default path to test data (relative to build directory)
    const char* default_dump = "../clus/test/data/init_first_segment_input.json";

    // Use env var if set, otherwise fall back to default
    const char* dump_path = std::getenv("WCT_TEST_DUMP");
    if (!dump_path) dump_path = default_dump;

    MESSAGE("Using dump file: ", dump_path);

    auto data = DebugIO::load_init_first_segment_inputs(dump_path);

    REQUIRE(data.cluster != nullptr);
    REQUIRE(data.main_cluster != nullptr);

    // Verify steiner_pc was loaded
    REQUIRE(data.cluster->has_pc("steiner_pc"));
    const auto& spc = data.cluster->get_pc("steiner_pc");
    CHECK(spc.size_major() > 0);
    MESSAGE("Steiner PC has ", spc.size_major(), " points");

    // Verify steiner_graph was loaded
    REQUIRE(data.cluster->has_graph("steiner_graph"));
    const auto& graph = data.cluster->get_graph("steiner_graph");
    CHECK(boost::num_vertices(graph) > 0);
    CHECK(boost::num_edges(graph) > 0);
    MESSAGE("Steiner graph: ", boost::num_vertices(graph), " vertices, ",
            boost::num_edges(graph), " edges");

    // Use precomputed boundary indices (computing them requires the anode)
    auto boundary_indices = data.boundary_steiner_indices;
    REQUIRE(boundary_indices.first != boundary_indices.second);
    MESSAGE("Boundary indices: ", boundary_indices.first, ", ", boundary_indices.second);

    // Test rough path computation
    const auto& scope = data.cluster->get_default_scope();
    MESSAGE("Scope name: ", scope.name, " with coords: ", scope.coords.size(), ", array size:", scope.coords[0].size());
    const auto& x_coords = spc.get(scope.coords.at(0))->elements<double>();
    const auto& y_coords = spc.get(scope.coords.at(1))->elements<double>();
    const auto& z_coords = spc.get(scope.coords.at(2))->elements<double>();

    Facade::geo_point_t first_pt(
        x_coords[boundary_indices.first],
        y_coords[boundary_indices.first],
        z_coords[boundary_indices.first]);
    Facade::geo_point_t second_pt(
        x_coords[boundary_indices.second],
        y_coords[boundary_indices.second],
        z_coords[boundary_indices.second]);

    MESSAGE("Boundary point 1: (", first_pt.x(), ", ", first_pt.y(), ", ", first_pt.z(), ")");
    MESSAGE("Boundary point 2: (", second_pt.x(), ", ", second_pt.y(), ", ", second_pt.z(), ")");

    // Test KNN on main_cluster steiner PC
    if (data.main_cluster->has_pc("steiner_pc")) {
        auto knn1 = data.main_cluster->kd_steiner_knn(1, first_pt, "steiner_pc");
        auto knn2 = data.main_cluster->kd_steiner_knn(1, second_pt, "steiner_pc");
        CHECK(!knn1.empty());
        CHECK(!knn2.empty());
        MESSAGE("KNN distances: ", std::sqrt(knn1[0].second), ", ", std::sqrt(knn2[0].second));
    }

    // Test rough path
    PatternAlgorithms pattern_algos;
    auto path_points = pattern_algos.do_rough_path(*data.cluster, first_pt, second_pt);
    CHECK(path_points.size() > 1);
    MESSAGE("Rough path has ", path_points.size(), " points");
}


TEST_CASE("init_first_segment end-to-end")
{
    // Full end-to-end test that calls init_first_segment with real data.
    // Requires DetectorVolumes and PCTransforms to be configured.
    //
    // Usage:
    //   - Without env vars: uses default data files in clus/test/data/
    //   - With env vars (for development):
    //     WCT_TEST_DUMP=/path/to/dump.json WCT_TEST_CONFIG=/path/to/config.json \
    //       wcdoctest-clus -tc="init_first_segment end-to-end"

    // Default paths to test data (relative to build directory)
    const char* default_dump = "../clus/test/data/init_first_segment_input.json";
    const char* default_config = "../clus/test/data/uboone-mabc_config.json";

    // Use env vars if set, otherwise fall back to defaults
    const char* dump_path = std::getenv("WCT_TEST_DUMP");
    const char* config_path = std::getenv("WCT_TEST_CONFIG");
    if (!dump_path) dump_path = default_dump;
    if (!config_path) config_path = default_config;

    MESSAGE("Using dump file: ", dump_path);
    MESSAGE("Using config file: ", config_path);

    // Load plugins
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellClus");
    pm.add("WireCellAux");
    pm.add("WireCellGen");
    pm.add("WireCellSigProc");

    // Load and apply WCT configuration
    auto wcfg_full = Persist::load(config_path);

    // Extract component configs from wire-cell JSON format
    // Format: [{"type":"wire-cell", "data":{...}}, <component configs>...]
    // The first item is metadata, rest are component configs
    Json::Value wcfg(Json::arrayValue);
    if (wcfg_full.isArray() && !wcfg_full.empty() &&
        wcfg_full[0]["type"].asString() == "wire-cell") {
        // Skip first item (metadata), rest are component configs
        for (Json::ArrayIndex i = 1; i < wcfg_full.size(); ++i) {
            wcfg.append(wcfg_full[i]);
        }
    } else {
        // Assume it's already just the configs array
        wcfg = wcfg_full;
    }

    configure_components(wcfg);

    // Get configured components
    auto dv = Factory::find_tn<IDetectorVolumes>("DetectorVolumes");
    auto pcts = Factory::find_tn<IPCTransformSet>("PCTransformSet");
    REQUIRE(dv);
    REQUIRE(pcts);

    // Collect anodes from the config (same way MultiAlgBlobClustering does)
    std::vector<IAnodePlane::pointer> anodes;
    for (const auto& comp : wcfg) {
        if (comp["type"].asString() == "AnodePlane") {
            std::string tn = "AnodePlane";
            std::string name = comp.get("name", "").asString();
            if (!name.empty()) tn += ":" + name;
            anodes.push_back(Factory::find_tn<IAnodePlane>(tn));
        }
    }
    REQUIRE(!anodes.empty());
    MESSAGE("Found ", anodes.size(), " AnodePlane(s)");

    // Load cluster data
    auto data = DebugIO::load_init_first_segment_inputs(dump_path);
    REQUIRE(data.cluster != nullptr);

    // Set anodes and detector volumes on the grouping (required by Blob/Cluster methods)
    auto* grouping = data.grouping_node->value.facade<Facade::Grouping>();
    grouping->set_anodes(anodes);
    grouping->set_detector_volumes(dv);

    // Create and configure TrackFitting
    auto tf = std::make_shared<TrackFitting>();
    tf->set_detector_volume(dv);
    tf->set_pc_transforms(pcts);
    tf->set_parameters(data.trackfitting_params);

    // Run init_first_segment
    PR::Graph pr_graph;
    PatternAlgorithms pattern_algos;
    auto seg = pattern_algos.init_first_segment(
        pr_graph, *data.cluster, data.main_cluster,
        *tf, dv, data.flag_back_search);

    // Validate results
    if (seg) {
        MESSAGE("init_first_segment succeeded");
        CHECK(boost::num_vertices(pr_graph) >= 2);

        // Check vertex positions
        auto [vbegin, vend] = boost::vertices(pr_graph);
        for (auto vit = vbegin; vit != vend; ++vit) {
            auto vtx = pr_graph[*vit].vertex;
            if (vtx) {
                MESSAGE("Vertex fit point: (",
                    vtx->fit().point.x(), ", ",
                    vtx->fit().point.y(), ", ",
                    vtx->fit().point.z(), ")");
                MESSAGE("  dQ=", vtx->fit().dQ, " dx=", vtx->fit().dx);
            }
        }
    } else {
        WARN("init_first_segment returned nullptr (path too short or fitting failed)");
    }
}
