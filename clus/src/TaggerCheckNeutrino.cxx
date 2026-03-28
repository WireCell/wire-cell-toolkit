#include "WireCellClus/TaggerCheckNeutrino.h"
#include "WireCellClus/NeutrinoPatternBase.h" // pattern recognition ...
#include "WireCellClus/PatternDebugIO.h"      // debug dump/load

#include "WireCellUtil/Persist.h"

#include <cstdlib>

class TaggerCheckNeutrino;
WIRECELL_FACTORY(TaggerCheckNeutrino, TaggerCheckNeutrino,
                 WireCell::INamed, WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::Clus::PR;


struct edge_base_t {
    typedef boost::edge_property_tag kind;
};

void TaggerCheckNeutrino::configure(const WireCell::Configuration& config)
{
    m_grouping_name = get(config, "grouping_name", m_grouping_name);
    m_trackfitting_config_file = get(config, "trackfitting_config_file", m_trackfitting_config_file);
    m_perf = get(config, "perf", m_perf);
    auto dl_weights_raw = get(config, "dl_weights", m_dl_weights);
    if (!dl_weights_raw.empty()) {
        m_dl_weights = Persist::resolve(dl_weights_raw);
        if (m_dl_weights.empty()) {
            SPDLOG_LOGGER_WARN(log, "TaggerCheckNeutrino: dl_weights path not found: {}", dl_weights_raw);
        }
    }
    m_dl_vtx_cut = get(config, "dl_vtx_cut", m_dl_vtx_cut);
    m_dQdx_scale  = get(config, "dQdx_scale",  m_dQdx_scale);
    m_dQdx_offset = get(config, "dQdx_offset", m_dQdx_offset);

    if (!m_trackfitting_config_file.empty()) {
        load_trackfitting_config(m_trackfitting_config_file);
    }

    NeedDV::configure(config);
    NeedPCTS::configure(config);
    NeedRecombModel::configure(config);
    NeedParticleData::configure(config);
}

Configuration TaggerCheckNeutrino::default_configuration() const
{
    Configuration cfg;
    cfg["grouping"] = m_grouping_name;
    cfg["detector_volumes"] = "DetectorVolumes";
    cfg["pc_transforms"] = "PCTransformSet";
    cfg["recombination_model"] = "BoxRecombination";  
    cfg["particle_dataset"] = "ParticleDataSet"; 

    cfg["trackfitting_config_file"] = "";
    cfg["perf"] = m_perf;
    cfg["dl_weights"] = "";      // empty = DL vertex disabled
    cfg["dl_vtx_cut"] = 20.0;   // mm (= 2 cm)
    cfg["dQdx_scale"]  = 0.1;   // dQ scale factor for SCN network input
    cfg["dQdx_offset"] = -1000.0; // dQ offset for SCN network input

    return cfg;
}

void TaggerCheckNeutrino::visit(Ensemble& ensemble) const
{
    // Configure the track fitter with detector volume
    m_track_fitter->set_detector_volume(m_dv);
    m_track_fitter->set_pc_transforms(m_pcts); 

    // Get the specified grouping (default: "live")
    auto groupings = ensemble.with_name(m_grouping_name);
    if (groupings.empty()) {
        return;
    }
    
    auto& grouping = *groupings.at(0);
    
    // Find clusters that have the main_cluster flag (set by clustering_recovering_bundle)
    Cluster* main_cluster = nullptr;
    std::vector<Cluster*> other_clusters;  // beam_flash clusters that are not the main cluster

    int nclusters = grouping.nchildren();
    int n_main_clusters = 0;
    int n_in_beam_clusters = 0;
    for (auto* cluster : grouping.children()) {
        if (cluster->get_flag(Flags::main_cluster)) {
            main_cluster = cluster;
            n_main_clusters ++;
        }
        if (cluster->get_flag(Flags::beam_flash)) n_in_beam_clusters++;
    }
    for (auto* cluster : grouping.children()) {
        if (cluster != main_cluster && cluster->get_flag(Flags::beam_flash)) {
            other_clusters.push_back(cluster);
        }
    }

    SPDLOG_LOGGER_DEBUG(log, "Found {} clusters, {} main clusters, {} in-beam clusters, {} of blobs in main cluster id {}", nclusters, n_main_clusters, n_in_beam_clusters, main_cluster->nchildren(), main_cluster->get_cluster_id());


    // // Debug dump (only when env var is set)
    // if (main_cluster) {
    //     if (const char* dump_path = std::getenv("WCT_DUMP_INIT_FIRST_SEGMENT")) {
    //         DebugIO::dump_init_first_segment_inputs(
    //             dump_path, *main_cluster, main_cluster, true, *m_track_fitter);
    //     }
    // }

    SPDLOG_LOGGER_DEBUG(log, "Number of Main Clusters: {}", n_main_clusters);

    std::set<VertexPtr> vertices_in_long_muon;
    std::set<SegmentPtr> segments_in_long_muon;
    VertexPtr main_vertex = nullptr;
    std::map<Cluster*, VertexPtr> map_cluster_main_vertices;

    // Pre-load charge data for all beam-flash clusters once so that
    // do_multi_tracking calls throughout pattern recognition can use
    // flag_force_load_data=false and avoid redundant prepare_data() calls.
    {
        std::vector<WireCell::Clus::Facade::Cluster*> clusters_to_preload;
        clusters_to_preload.push_back(main_cluster);
        for (auto* c : other_clusters) clusters_to_preload.push_back(c);
        m_track_fitter->preload_clusters(clusters_to_preload);
    }

    // Create PRGraph and first segment
    auto pr_graph = std::make_shared<WireCell::Clus::PR::Graph>();
    m_track_fitter->add_graph(pr_graph);

    WireCell::Clus::PR::PatternAlgorithms pattern_algos;
    pattern_algos.m_perf = m_perf;
    m_track_fitter->set_perf(m_perf);

    {
        // initial pattern recognitions
        pattern_algos.find_proto_vertex(*pr_graph, *main_cluster, *m_track_fitter, m_dv, true, 2, true);

        // shower related operations
        pattern_algos.clustering_points(*pr_graph, *main_cluster, m_dv);
        pattern_algos.separate_track_shower(*pr_graph, *main_cluster);
        
        // direction determination
        pattern_algos.determine_direction(*pr_graph, *main_cluster, particle_data(), m_recomb_model);

        // shower clustering
        pattern_algos.shower_determining_in_main_cluster(*pr_graph, *main_cluster, particle_data(), m_recomb_model, m_dv);

        // main vertex determination
        pattern_algos.determine_main_vertex(*pr_graph, *main_cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, *m_track_fitter, m_dv, particle_data(), m_recomb_model);

        if (main_vertex !=nullptr){
            map_cluster_main_vertices[main_cluster] = main_vertex;
            main_vertex = nullptr;
        }
    }


    // Loop over other (non-main) beam-flash clusters
    if (!other_clusters.empty()) {
        for (auto* cluster : other_clusters) {
            if (cluster->get_length() > 6 * units::cm) {
                // std::cout << "Long Cluster " << cluster->get_cluster_id() << " " << cluster->nchildren() << std::endl;
                // Long cluster: break tracks and do 2 rounds of other-track finding
                pattern_algos.find_proto_vertex(*pr_graph, *cluster, *m_track_fitter, m_dv, true, 2, false);
                pattern_algos.clustering_points(*pr_graph, *cluster, m_dv);
                pattern_algos.separate_track_shower(*pr_graph, *cluster);
                pattern_algos.determine_direction(*pr_graph, *cluster, particle_data(), m_recomb_model);
                pattern_algos.shower_determining_in_main_cluster(*pr_graph, *cluster, particle_data(), m_recomb_model, m_dv);
                pattern_algos.determine_main_vertex(*pr_graph, *cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, *m_track_fitter, m_dv, particle_data(), m_recomb_model);
                if (main_vertex != nullptr) {
                    map_cluster_main_vertices[cluster] = main_vertex;
                    main_vertex = nullptr;
                }
            } else {
                // Short cluster: no track breaking, 1 round; fall back to init_point_segment if needed
                if (!pattern_algos.find_proto_vertex(*pr_graph, *cluster, *m_track_fitter, m_dv, false, 1, false)) {
                    // std::cout << "Point Cluster " << cluster->get_cluster_id() << " " << cluster->nchildren() <<std::endl;
                    pattern_algos.init_point_segment(*pr_graph, *cluster, *m_track_fitter, m_dv);
                }
                pattern_algos.clustering_points(*pr_graph, *cluster, m_dv);
                pattern_algos.separate_track_shower(*pr_graph, *cluster);
                pattern_algos.determine_direction(*pr_graph, *cluster, particle_data(), m_recomb_model);
                pattern_algos.shower_determining_in_main_cluster(*pr_graph, *cluster, particle_data(), m_recomb_model, m_dv);
                pattern_algos.determine_main_vertex(*pr_graph, *cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, *m_track_fitter, m_dv, particle_data(), m_recomb_model);
                if (main_vertex != nullptr) {
                    map_cluster_main_vertices[cluster] = main_vertex;
                    main_vertex = nullptr;
                }
            }
        }

        // Deghost across all beam-flash clusters (main + others)
        std::vector<Cluster*> all_clusters;
        all_clusters.push_back(main_cluster);
        all_clusters.insert(all_clusters.end(), other_clusters.begin(), other_clusters.end());
        
        pattern_algos.deghosting(*pr_graph, map_cluster_main_vertices, all_clusters, *m_track_fitter, m_dv);
    }
     
    // Determine the overall neutrino vertex.
    // If DL weights are configured, try DL first (matches prototype flag_dl_vtx logic).
    // Fall back to traditional algorithm if DL is disabled or does not change the vertex.
    // DL path updates map_cluster_main_vertices[main_cluster] directly (by-ref parameter).
    // Traditional path returns the chosen vertex; capture it and sync to the map.
    VertexPtr final_main_vertex = nullptr;
    bool flag_dl_changed = false;
    if (!m_dl_weights.empty()) {
        flag_dl_changed = pattern_algos.determine_overall_main_vertex_DL(
            *pr_graph, map_cluster_main_vertices, main_cluster, other_clusters,
            vertices_in_long_muon, segments_in_long_muon,
            *m_track_fitter, m_dv, particle_data(), m_recomb_model,
            m_dl_weights, m_dl_vtx_cut, m_dQdx_scale, m_dQdx_offset);
    }
    if (!flag_dl_changed) {
        final_main_vertex = pattern_algos.determine_overall_main_vertex(
            *pr_graph, map_cluster_main_vertices, main_cluster, other_clusters,
            vertices_in_long_muon, segments_in_long_muon,
            *m_track_fitter, m_dv, particle_data(), m_recomb_model, true);
        if (final_main_vertex) {
            map_cluster_main_vertices[main_cluster] = final_main_vertex;
        }
    }

    // Retrieve the chosen neutrino vertex regardless of which path ran
    {
        auto it = map_cluster_main_vertices.find(main_cluster);
        if (it != map_cluster_main_vertices.end()) {
            final_main_vertex = it->second;
        }
    }

    

    // // Post-vertex refinement (matches prototype block after determine_overall_main_vertex):
    // //   1. Minuit-based vertex position fit
    // //   2. Re-cluster EM shower points with refined vertex
    // //   3. Re-examine track directions (flag_final=true)
    // //   4. Re-separate tracks and showers
    // std::size_t n_main_cluster_vertices = 0;
    // std::size_t n_main_cluster_segments = 0;
    // std::size_t n_main_cluster_fit_points = 0;
    // for (const auto& nd : PR::graph_nodes(*pr_graph)) {
    //     const auto& vtx = (*pr_graph)[nd].vertex;
    //     if (vtx && vtx->cluster() == main_cluster) {
    //         ++n_main_cluster_vertices;
    //     }
    // }
    // for (const auto& ed : PR::ordered_edges(*pr_graph)) {
    //     const auto& seg = (*pr_graph)[ed].segment;
    //     if (seg && seg->cluster() == main_cluster) {
    //         ++n_main_cluster_segments;
    //         n_main_cluster_fit_points += seg->fits().size();
    //     }
    // }
    // SPDLOG_LOGGER_DEBUG(log,
    //                     "Debug Cluster {} has vertices={} segments={} fit_points={} in PR graph",
    //                     main_cluster->get_cluster_id(), n_main_cluster_vertices,
    //                     n_main_cluster_segments, n_main_cluster_fit_points);

    int acc_segment_id = 0;
    std::set<ShowerPtr> pi0_showers;
    std::map<ShowerPtr, int> map_shower_pio_id;
    std::map<int, std::vector<ShowerPtr>> map_pio_id_showers;
    std::map<int, std::pair<double, int>> map_pio_id_mass;
    std::map<int, std::pair<int, int>> map_pio_id_saved_pair;
    std::map<VertexPtr, ShowerPtr> map_vertex_in_shower;
    std::map<SegmentPtr, ShowerPtr> map_segment_in_shower;
    std::map<VertexPtr, std::set<ShowerPtr>> map_vertex_to_shower;
    std::set<Facade::Cluster*> used_shower_clusters;
    std::set<ShowerPtr> showers;

    if (final_main_vertex) {
   

        pattern_algos.improve_vertex(*pr_graph, *main_cluster, final_main_vertex,
                                     vertices_in_long_muon, segments_in_long_muon,
                                     *m_track_fitter, m_dv, particle_data(), m_recomb_model,
                                     true, true);
        // improve_vertex may update final_main_vertex pointer; sync back to map
        map_cluster_main_vertices[main_cluster] = final_main_vertex;

        pattern_algos.clustering_points(*pr_graph, *main_cluster, m_dv);

        pattern_algos.examine_direction(*pr_graph, final_main_vertex, final_main_vertex,
                                        vertices_in_long_muon, segments_in_long_muon,
                                        particle_data(), m_recomb_model, true);

        SPDLOG_LOGGER_DEBUG(log, "Overall main vertex cluster={}", main_cluster->get_cluster_id());
        // pattern_algos.print_segs_info(*pr_graph, *main_cluster, final_main_vertex);

        pattern_algos.separate_track_shower(*pr_graph, *main_cluster);

        pattern_algos.shower_clustering_with_nv(acc_segment_id, pi0_showers,
                                                map_shower_pio_id, map_pio_id_showers,
                                                map_pio_id_mass, map_pio_id_saved_pair,
                                                vertices_in_long_muon, segments_in_long_muon,
                                                *pr_graph, final_main_vertex, showers,
                                                main_cluster, other_clusters,
                                                map_cluster_main_vertices,
                                                map_vertex_in_shower, map_segment_in_shower,
                                                map_vertex_to_shower, used_shower_clusters,
                                                *m_track_fitter, m_dv, particle_data(),
                                                m_recomb_model);

        
    }

    // n_main_cluster_vertices = 0;
    // n_main_cluster_segments = 0;
    // n_main_cluster_fit_points = 0;
    // for (const auto& nd : PR::graph_nodes(*pr_graph)) {
    //     const auto& vtx = (*pr_graph)[nd].vertex;
    //     if (vtx && vtx->cluster() == main_cluster) {
    //         ++n_main_cluster_vertices;
    //     }
    // }
    // for (const auto& ed : PR::ordered_edges(*pr_graph)) {
    //     const auto& seg = (*pr_graph)[ed].segment;
    //     if (seg && seg->cluster() == main_cluster) {
    //         ++n_main_cluster_segments;
    //         n_main_cluster_fit_points += seg->fits().size();
    //     }
    // }
    // SPDLOG_LOGGER_DEBUG(log,
    //                     "Debug Cluster {} has vertices={} segments={} fit_points={} in PR graph",
    //                     main_cluster->get_cluster_id(), n_main_cluster_vertices,
    //                     n_main_cluster_segments, n_main_cluster_fit_points);


    // Mark each cluster's main vertex so bee output can identify it
    // for (auto& [cluster, vtx] : map_cluster_main_vertices) {
    //     if (vtx) {
    //         vtx->set_flags(PR::VertexFlags::kNeutrinoVertex);

    //         const auto& wcpt = vtx->wcpt().point;
    //         if (vtx->fit().valid()) {
    //             const auto& fitpt = vtx->fit().point;
    //             SPDLOG_LOGGER_DEBUG(log,
    //                                 "Cluster {} neutrino vertex wcpt=({:.2f}, {:.2f}, {:.2f}) cm fit=({:.2f}, {:.2f}, {:.2f}) cm",
    //                                 cluster->get_cluster_id(),
    //                                 wcpt.x() / units::cm, wcpt.y() / units::cm, wcpt.z() / units::cm,
    //                                 fitpt.x() / units::cm, fitpt.y() / units::cm, fitpt.z() / units::cm);
    //         }
    //         else {
    //             SPDLOG_LOGGER_DEBUG(log,
    //                                 "Cluster {} neutrino vertex wcpt=({:.2f}, {:.2f}, {:.2f}) cm fit=invalid",
    //                                 cluster->get_cluster_id(),
    //                                 wcpt.x() / units::cm, wcpt.y() / units::cm, wcpt.z() / units::cm);
    //         }
    //     }
    // }

    // Store TrackFitting in the grouping for later access by bee output and tracking sink
    grouping.set_track_fitting(m_track_fitter);
}

void TaggerCheckNeutrino::load_trackfitting_config(const std::string& config_file)
{
    try {
        // Load JSON file
        std::ifstream file(config_file);
        if (!file.is_open()) {
            std::cerr << "TaggerCheckNeutrino: Cannot open config file: " << config_file << std::endl;
            return;
        }
        
        Json::Value root;
        Json::CharReaderBuilder builder;
        std::string errs;
        
        if (!Json::parseFromStream(builder, file, &root, &errs)) {
            std::cerr << "TaggerCheckNeutrino: Failed to parse JSON: " << errs << std::endl;
            return;
        }
        
        // Apply each parameter from the JSON file
        for (const auto& param_name : root.getMemberNames()) {
            if (param_name.substr(0, 1) == "_") continue;  // Skip comments
            
            try {
                double value = root[param_name].asDouble();
                m_track_fitter->set_parameter(param_name, value);
                // SPDLOG_LOGGER_DEBUG(log, "Set {} = {}", param_name, value);
            } catch (const std::exception& e) {
                std::cerr << "TaggerCheckNeutrino: Failed to set parameter " << param_name 
                        << ": " << e.what() << std::endl;
            }
        }
        
        SPDLOG_LOGGER_DEBUG(log, "Successfully loaded TrackFitting configuration");
        
    } catch (const std::exception& e) {
        std::cerr << "TaggerCheckNeutrino: Exception loading config: " << e.what() << std::endl;
        std::cerr << "TaggerCheckNeutrino: Using default TrackFitting parameters" << std::endl;
    }
}

