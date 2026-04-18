#include "WireCellClus/PatternDebugIO.h"
#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade.h"
#include "WireCellClus/NeutrinoTaggerInfo.h"
#include "WireCellClus/ParticleDataSet.h"

#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IRecombinationModel.h"
#include "WireCellClus/IPCTransform.h"

#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

#include <cstdlib>
#include <cmath>

// Usage:
//   # Generate fixture (run from qlport/ directory):
//   WCT_DUMP_TAGGER_INPUTS=./tmp/tagger_check_neutrino_input.json \
//     wire-cell -A infiles=rootfiles/nuselEval_5384_130_6501.root uboone-mabc.jsonnet
//
//   # Generate WCT config:
//   wcsonnet -o ./tmp/uboone-mabc_config.json uboone-mabc.jsonnet
//
//   # Run test (default fixtures from clus/test/data/):
//   wcdoctest-clus -tc="tagger_check_neutrino end-to-end"
//
//   # Run with custom fixtures:
//   WCT_TEST_DUMP=/path/to/tagger_check_neutrino_input.json \
//   WCT_TEST_CONFIG=/path/to/uboone-mabc_config.json \
//     wcdoctest-clus -tc="tagger_check_neutrino end-to-end"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::PR;
using namespace WireCell::Clus::Facade;

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
            if (comp.isMember("data")) {
                for (const auto& key : comp["data"].getMemberNames()) {
                    cfg[key] = comp["data"][key];
                }
            }
            icfg->configure(cfg);
            log->trace("Configured {}", tn);
        } catch (const std::exception& e) {
            log->trace("Skipping {}: {}", tn, e.what());
        }
    }
}


TEST_CASE("tagger_check_neutrino end-to-end")
{
    // Default paths (relative to build directory, matching init_first_segment precedent).
    const char* default_dump   = "../clus/test/data/tagger_check_neutrino_input.json";
    const char* default_config = "../clus/test/data/uboone-mabc_config.json";

    const char* dump_path   = std::getenv("WCT_TEST_DUMP");
    const char* config_path = std::getenv("WCT_TEST_CONFIG");
    if (!dump_path)   dump_path   = default_dump;
    if (!config_path) config_path = default_config;

    std::string s_dump_path   = dump_path;
    std::string s_config_path = config_path;
    MESSAGE("Dump file:   ", s_dump_path);
    MESSAGE("Config file: ", s_config_path);

    // Load plugins
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellClus");
    pm.add("WireCellAux");
    pm.add("WireCellGen");
    pm.add("WireCellSigProc");

    // Load and apply WCT configuration
    auto wcfg_full = Persist::load(config_path);
    Json::Value wcfg(Json::arrayValue);
    if (wcfg_full.isArray() && !wcfg_full.empty() &&
        wcfg_full[0]["type"].asString() == "wire-cell") {
        for (Json::ArrayIndex i = 1; i < wcfg_full.size(); ++i) {
            wcfg.append(wcfg_full[i]);
        }
    } else {
        wcfg = wcfg_full;
    }
    configure_components(wcfg);

    // Retrieve configured geometry/physics components
    auto dv   = Factory::find_tn<IDetectorVolumes>("DetectorVolumes");
    auto pcts = Factory::find_tn<IPCTransformSet>("PCTransformSet");
    REQUIRE(dv);
    REQUIRE(pcts);

    // Collect anodes
    std::vector<IAnodePlane::pointer> anodes;
    for (const auto& comp : wcfg) {
        if (comp["type"].asString() == "AnodePlane") {
            std::string tn = "AnodePlane";
            std::string nm = comp.get("name", "").asString();
            if (!nm.empty()) tn += ":" + nm;
            anodes.push_back(Factory::find_tn<IAnodePlane>(tn));
        }
    }
    REQUIRE(!anodes.empty());

    // Get particle dataset and recombination model
    ParticleDataSet::pointer pdata;
    IRecombinationModel::pointer recomb;
    try {
        auto icfg = Factory::find_tn<IConfigurable>("ParticleDataSet:ParticleDataSet");
        pdata = std::dynamic_pointer_cast<ParticleDataSet>(icfg);
    } catch (...) {}
    REQUIRE(pdata);

    try {
        recomb = Factory::find_tn<IRecombinationModel>("BoxRecombination:box_recomb");
    } catch (...) {}
    REQUIRE(recomb);

    // Load cluster fixture
    auto data = DebugIO::load_tagger_inputs(dump_path);
    REQUIRE(data.cluster != nullptr);
    REQUIRE(data.main_cluster != nullptr);

    // Set geometry on grouping
    auto* grouping = data.grouping_node->value.facade<Grouping>();
    grouping->set_anodes(anodes);
    grouping->set_detector_volumes(dv);

    // Create and configure TrackFitting
    auto tf = std::make_shared<TrackFitting>();
    tf->set_detector_volume(dv);
    tf->set_pc_transforms(pcts);
    tf->set_parameters(data.trackfitting_params);

    // Preload charge data (re-derived from blob 3D PCs + geometry)
    {
        std::vector<Cluster*> to_preload;
        to_preload.push_back(data.main_cluster);
        for (auto* c : data.other_clusters) to_preload.push_back(c);
        tf->preload_clusters(to_preload);
    }

    // --- Pattern recognition (mirrors TaggerCheckNeutrino::visit()) ---
    auto pr_graph = std::make_shared<Graph>();
    tf->add_graph(pr_graph);

    PatternAlgorithms pattern_algos;

    IndexedVertexSet   vertices_in_long_muon;
    IndexedSegmentSet  segments_in_long_muon;
    ClusterVertexMap   map_cluster_main_vertices;
    VertexPtr          main_vertex = nullptr;

    // Initial PR on main cluster
    pattern_algos.find_proto_vertex(*pr_graph, *data.main_cluster, *tf, dv, true, 2, true);
    pattern_algos.clustering_points(*pr_graph, *data.main_cluster, dv);
    pattern_algos.separate_track_shower(*pr_graph, *data.main_cluster);
    pattern_algos.determine_direction(*pr_graph, *data.main_cluster, pdata, recomb);
    pattern_algos.shower_determining_in_main_cluster(*pr_graph, *data.main_cluster, pdata, recomb, dv);
    pattern_algos.determine_main_vertex(*pr_graph, *data.main_cluster,
                                        main_vertex, vertices_in_long_muon,
                                        segments_in_long_muon, *tf, dv, pdata, recomb);
    if (main_vertex) {
        map_cluster_main_vertices[data.main_cluster] = main_vertex;
        main_vertex = nullptr;
    }

    // Determine overall neutrino vertex (traditional path; no DL in test config)
    VertexPtr final_main_vertex = pattern_algos.determine_overall_main_vertex(
        *pr_graph, map_cluster_main_vertices, data.main_cluster, data.other_clusters,
        vertices_in_long_muon, segments_in_long_muon, *tf, dv, pdata, recomb, true);
    if (final_main_vertex) {
        map_cluster_main_vertices[data.main_cluster] = final_main_vertex;
    }

    MESSAGE("final_main_vertex found: ", (final_main_vertex != nullptr));

    // Showers needed by taggers (may be empty if vertex not found; taggers handle that)
    IndexedShowerSet showers;
    ShowerVertexMap  map_vertex_in_shower;
    ShowerSegmentMap map_segment_in_shower;
    VertexShowerSetMap map_vertex_to_shower;
    ClusterPtrSet    used_shower_clusters;
    int              acc_segment_id = 0;
    IndexedShowerSet pi0_showers;
    ShowerIntMap     map_shower_pio_id;
    std::map<int, std::vector<ShowerPtr>> map_pio_id_showers;
    std::map<int, std::pair<double, int>> map_pio_id_mass;
    std::map<int, std::pair<int, int>>    map_pio_id_saved_pair;
    Pi0KineFeatures  pio_kine{};

    if (final_main_vertex) {
        pattern_algos.improve_vertex(*pr_graph, *data.main_cluster, final_main_vertex,
                                     vertices_in_long_muon, segments_in_long_muon,
                                     *tf, dv, pdata, recomb, true, true);
        map_cluster_main_vertices[data.main_cluster] = final_main_vertex;

        pattern_algos.clustering_points(*pr_graph, *data.main_cluster, dv);
        pattern_algos.examine_direction(*pr_graph, final_main_vertex, final_main_vertex,
                                        vertices_in_long_muon, segments_in_long_muon,
                                        pdata, recomb, true);

        pattern_algos.shower_clustering_with_nv(
            acc_segment_id, pi0_showers,
            map_shower_pio_id, map_pio_id_showers,
            map_pio_id_mass, map_pio_id_saved_pair,
            pio_kine,
            vertices_in_long_muon, segments_in_long_muon,
            *pr_graph, final_main_vertex, showers,
            data.main_cluster, data.other_clusters,
            map_cluster_main_vertices,
            map_vertex_in_shower, map_segment_in_shower,
            map_vertex_to_shower, used_shower_clusters,
            *tf, dv, pdata, recomb);
    }

    // --- Run taggers ---
    TaggerInfo tagger_info;
    pattern_algos.init_tagger_info(tagger_info);

    std::vector<Cluster*> all_clusters;
    all_clusters.push_back(data.main_cluster);
    for (auto* c : data.other_clusters) all_clusters.push_back(c);

    if (final_main_vertex) {
        pattern_algos.cosmic_tagger(*pr_graph, final_main_vertex,
                                    showers, map_segment_in_shower, map_vertex_to_shower,
                                    segments_in_long_muon,
                                    data.main_cluster, all_clusters, dv,
                                    tagger_info);

        pattern_algos.numu_tagger(*pr_graph, final_main_vertex,
                                  showers, segments_in_long_muon,
                                  data.main_cluster, tagger_info);
    }

    // --- Assertions ---

    // If a vertex was found, cosmic tagger must have run and filled its fields.
    if (final_main_vertex) {
        CHECK(tagger_info.cosmic_filled == doctest::Approx(1.0f));
        CHECK(std::isfinite(tagger_info.cosmic_flag));
        CHECK(std::isfinite(tagger_info.numu_cc_flag));
    }

    // Top-level flags must be finite regardless of whether a vertex was found.
    CHECK(std::isfinite(tagger_info.cosmic_flag));
    CHECK(std::isfinite(tagger_info.gap_flag));
    CHECK(std::isfinite(tagger_info.mip_flag));
    CHECK(std::isfinite(tagger_info.mip_quality_flag));

    // SSM sentinels should still hold the default -999 (SSM tagger not called).
    CHECK(tagger_info.ssm_Nsm == doctest::Approx(-999.0f));

    MESSAGE("cosmic_filled = ", tagger_info.cosmic_filled);
    MESSAGE("cosmic_flag   = ", tagger_info.cosmic_flag);
    MESSAGE("numu_cc_flag  = ", tagger_info.numu_cc_flag);
}
