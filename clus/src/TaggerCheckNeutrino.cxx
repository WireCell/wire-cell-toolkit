#include "WireCellClus/TaggerCheckNeutrino.h"
#include "WireCellClus/NeutrinoPatternBase.h" // pattern recognition ...
#include "WireCellClus/PatternDebugIO.h"      // debug dump/load

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

        // // shower clustering
        // pattern_algos.shower_determining_in_main_cluster(*pr_graph, *main_cluster, particle_data(), m_recomb_model, m_dv);

        // // main vertex determination
        // pattern_algos.determine_main_vertex(*pr_graph, *main_cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, *m_track_fitter, m_dv, particle_data(), m_recomb_model, true);

        // if (main_vertex !=0){
        //     map_cluster_main_vertices[main_cluster] = main_vertex;
        //     main_vertex = 0;
        // }
    }


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

