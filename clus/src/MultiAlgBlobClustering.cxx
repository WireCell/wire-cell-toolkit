#include "WireCellClus/MultiAlgBlobClustering.h"
#include "WireCellClus/Facade_Summary.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRShower.h"
#include "WireCellClus/TrackFitting.h"


#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMdataset.h"
#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/ExecMon.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/GraphTools.h"

#include <map>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdio>

WIRECELL_FACTORY(MultiAlgBlobClustering, WireCell::Clus::MultiAlgBlobClustering, WireCell::INamed,
                 WireCell::ITensorSetFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;
using WireCell::GraphTools::mir;

MultiAlgBlobClustering::MultiAlgBlobClustering()
  : Aux::Logger("MultiAlgBlobClustering", "clus")
//   , m_bee_dead("channel-deadarea", 1*units::mm, 3) // tolerance, minpts
{
}


static
std::string format_path(
    std::string path,
    const std::string& name,
    int ident,
    const std::map<std::string, std::string> subpaths)
{
    auto it = subpaths.find(name);
    if (it == subpaths.end()) {
        path += "/" + name;
    }
    else {
        path += it->second;
    }
    if (path.find("%") == std::string::npos) {
        return path;
    }
    return String::format(path, ident);
}

std::string MultiAlgBlobClustering::inpath(const std::string& name, int ident)
{
    return format_path(m_inpath, name, ident, m_insubpaths);
}
std::string MultiAlgBlobClustering::outpath(const std::string& name, int ident)
{
    return format_path(m_outpath, name, ident, m_outsubpaths);
}


void MultiAlgBlobClustering::configure(const WireCell::Configuration& cfg)
{
    m_groupings = convert(cfg["groupings"], m_groupings);

    m_inpath = get(cfg, "inpath", m_inpath);
    m_outpath = get(cfg, "outpath", m_outpath);

    for (const auto& jsp : cfg["insubpaths"]) {
        m_insubpaths[jsp["name"].asString()] = jsp["subpath"].asString();
    }
    for (const auto& jsp : cfg["outsubpaths"]) {
        m_outsubpaths[jsp["name"].asString()] = jsp["subpath"].asString();
    }

    {
        auto jcid = cfg["cluster_id_order"];
        if (jcid.isString()) {
            m_clusters_id_order = jcid.asString();
        }
    }

    if (cfg.isMember("bee_dir")) {
        SPDLOG_LOGGER_DEBUG(log, "the 'bee_dir' option is no longer supported, instead use 'bee_zip' to name a .zip file");
    }
    std::string bee_zip = get<std::string>(cfg, "bee_zip", "mabc.zip");
    // Add new configuration option for initial index
    m_initial_index = get<int>(cfg, "initial_index", m_initial_index);

    //std::cout << "Xin: " << m_initial_index << " " << bee_zip << std::endl;
    m_sink.reset(bee_zip, m_initial_index);  // Use the new reset with initial index

    // Configure RSE numbers
    if (cfg.isMember("use_config_rse")) {
        m_use_config_rse = get(cfg, "use_config_rse", false);
        if (m_use_config_rse) {
            // Only read RSE if we're using configured values
            m_runNo = get(cfg, "runNo", m_runNo);
            m_subRunNo = get(cfg, "subRunNo", m_subRunNo);
            m_eventNo = get(cfg, "eventNo", m_eventNo);
            
             // Set RSE in sink during configuration
            m_sink.set_rse(m_runNo, m_subRunNo, m_eventNo);
        }
    }

    m_grouping2file_prefix = get(cfg, "grouping2file_prefix", m_grouping2file_prefix);

    m_save_deadarea = get(cfg, "save_deadarea", m_save_deadarea);

    m_dead_live_overlap_offset = get(cfg, "dead_live_overlap_offset", m_dead_live_overlap_offset);

    for (auto jtn : cfg["pipeline"]) {
        std::string tn = jtn.asString();
        SPDLOG_LOGGER_DEBUG(log, "configuring clustering method: {}", tn);
        auto imeth = Factory::find_tn<IEnsembleVisitor>(tn);
        m_pipeline.emplace_back(EnsembleVisitor{tn, imeth});
    }

    m_perf = get(cfg, "perf", m_perf);

    for (const auto& aname : cfg["anodes"]) {
        auto anode = Factory::find_tn<IAnodePlane>(aname.asString());
        m_anodes.push_back(anode);
    }

    m_dv = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());

    m_dump_json = get<bool>(cfg, "dump_json", false);

    // Configure bee points sets
    if (cfg.isMember("bee_points_sets")) {
        auto bee_points_sets = cfg["bee_points_sets"];
        for (const auto& bps : bee_points_sets) {
            BeePointsConfig bpc;
            bpc.name = get<std::string>(bps, "name", "");
            bpc.detector = get<std::string>(bps, "detector", "uboone");
            bpc.algorithm = get<std::string>(bps, "algorithm", bpc.name);
            bpc.pcname = get<std::string>(bps, "pcname", "3d");
            bpc.grouping = get<std::string>(bps, "grouping", "live");
            bpc.visitor = get<std::string>(bps, "visitor", "");
            bpc.filter = get<int>(bps, "filter", 1); // 1 for on, 0 for off, -1 for inverse filter
            
            // Get coordinates
            if (bps.isMember("coords")) {
                for (const auto& coord : bps["coords"]) {
                    bpc.coords.push_back(coord.asString());
                }
            } else {
                // Default coordinates
                bpc.coords = {"x", "y", "z"};
            }
            
            bpc.individual = get<bool>(bps, "individual", false);
            bpc.dQdx_scale = get<double>(bps, "dQdx_scale", 1.0);
            bpc.dQdx_offset = get<double>(bps, "dQdx_offset", 0.0);
            bpc.use_associate_points = get<bool>(bps, "use_associate_points", false);
            bpc.use_graph_vertices = get<bool>(bps, "use_graph_vertices", false);

            m_bee_points_configs.push_back(bpc);
            
            
            // If individual, also initialize bee points for each APA and face
            if (bpc.individual) {
                for (const auto& anode : m_anodes) {
                    int apa = anode->ident();
                    // Initialize the outer map if it doesn't exist
                    if (m_bee_points[bpc.name].by_apa_face.find(apa) == 
                        m_bee_points[bpc.name].by_apa_face.end()) {
                        m_bee_points[bpc.name].by_apa_face[apa] = std::map<int, Bee::Points>();
                    }
                    
                    // Initialize bee points for each face
                    for (size_t face_index = 0; face_index < anode->faces().size(); ++face_index) {
                        int face = anode->faces()[face_index]->which();
                        std::string algo_name = String::format("%s-apa%d-face%d", bpc.algorithm.c_str(), apa,  face);
                        // std::cout << "Test: Individual: " << algo_name << std::endl;
                        m_bee_points[bpc.name].by_apa_face[apa][face] =  Bee::Points(bpc.detector, algo_name);
                    }
                }
            }else{
                m_bee_points[bpc.name].global.detector(bpc.detector);
                m_bee_points[bpc.name].global.algorithm(String::format("%s-global", bpc.name));
                // std::cout << "Test: Global: " << m_bee_points[bpc.name].global.algorithm() << std::endl;
            }
            
            SPDLOG_LOGGER_DEBUG(log, "Configured bee points set: {}, algorithm: {}, individual: {}", 
                        bpc.name, bpc.algorithm, bpc.individual ? "true" : "false");
        }
    } 

    // Configure particle-flow bee output (mc tree)
    if (cfg.isMember("bee_pf")) {
        for (const auto& pf : cfg["bee_pf"]) {
            BeePFConfig pfc;
            pfc.name     = get<std::string>(pf, "name",     "mc");
            pfc.visitor  = get<std::string>(pf, "visitor",  "");
            pfc.grouping = get<std::string>(pf, "grouping", "live");
            m_bee_pf_configs.push_back(pfc);
            m_bee_pf_trees[pfc.name] = Bee::ParticleTree(pfc.name);
            SPDLOG_LOGGER_DEBUG(log, "Configured bee_pf: name={} visitor={}", pfc.name, pfc.visitor);
        }
    }

    // Initialize patches for each APA and face
    if (m_save_deadarea) {
        for (const auto& anode : m_anodes) {
            int apa = anode->ident();
            
            // Initialize the outer map if it doesn't exist
            if (m_bee_dead_patches.find(apa) == 
                m_bee_dead_patches.end()) {
                m_bee_dead_patches[apa] = std::map<int, Bee::Patches>();
            }
            
            // Initialize patches for each face
            for (size_t face_index = 0; face_index < anode->faces().size(); ++face_index) {
                int face = anode->faces()[face_index]->which();
                std::string name = String::format("channel-deadarea-apa%d-face%d", apa, face);
                m_bee_dead_patches[apa].insert({face,Bee::Patches(name, 1*units::mm, 3)}); // Same parameters as the global one
            }
        }
    }
}

WireCell::Configuration MultiAlgBlobClustering::default_configuration() const
{
    Configuration cfg;

    assign(cfg["groupings"], m_groupings);

    cfg["inpath"] = m_inpath;
    cfg["outpath"] = m_outpath;

    // repeat defaults as literals just incase some "clever" person tries to
    // call this method AFTER configure() as that method mutates m_inlive, etc.
    cfg["inlive"] = "/live";
    cfg["outlive"] = "/live";
    cfg["indead"] = "/dead";
    cfg["outdead"] = "/dead";

    // cfg["bee_dir"] = m_bee_dir;
    cfg["bee_zip"] = "mabc.zip";
    cfg["save_deadarea"] = m_save_deadarea;

    // Add the new parameter to default configuration
    cfg["initial_index"] = m_initial_index;

    cfg["dead_live_overlap_offset"] = m_dead_live_overlap_offset;

    cfg["use_config_rse"] = false;  // By default, don't use configured RSE
    cfg["runNo"] = m_runNo;
    cfg["subRunNo"] = m_subRunNo;
    cfg["eventNo"] = m_eventNo;

    return cfg;
}

void MultiAlgBlobClustering::finalize()
{
    flush();
    m_sink.close();
}

static void reset_bee(int ident, WireCell::Bee::Points& bpts)
{
    int run=0, evt=0;
    if (ident > 0) {
        run = (ident >> 16) & 0x7fff;
        evt = (ident) & 0xffff;
    }
    bpts.reset(evt, 0, run);
}

void MultiAlgBlobClustering::flush(WireCell::Bee::Points& bpts, int ident)
{
    if (bpts.empty()) return;

    m_sink.write(bpts);
    reset_bee(ident, bpts);
}

void MultiAlgBlobClustering::flush(int ident)
{
    // flush(m_bee_img, ident);
    // flush(m_bee_ld,  ident);
     // Flush all bee points sets

     for (auto& [name, apa_bpts] : m_bee_points) {
         // C++17 can not use structured bindings in lambda capture list.
         const std::string the_name = name;

        // Find the configuration for this name to check if it's individual
        auto it = std::find_if(m_bee_points_configs.begin(), m_bee_points_configs.end(),
                              [&the_name](const BeePointsConfig& cfg) { return cfg.name == the_name; });
        
        bool individual = (it != m_bee_points_configs.end()) ? it->individual : false;
        
        if (individual) {
            // Write individual bee points
            for (auto& [anode_id, face_map] : apa_bpts.by_apa_face) {
                for (auto& [face, bpts] : face_map) {
                    if (!bpts.empty()) {
                        m_sink.write(bpts);
                        // Clear after writing
                        int run = 0, evt = 0;
                        if (ident > 0) {
                            run = (ident >> 16) & 0x7fff;
                            evt = (ident) & 0xffff;
                        }
                        bpts.reset(evt, 0, run);
                    }
                }
            }
        } else {
            // Write global bee points
            if (!apa_bpts.global.empty()) {
                m_sink.write(apa_bpts.global);
                // Clear after writing
                int run = 0, evt = 0;
                if (ident > 0) {
                    run = (ident >> 16) & 0x7fff;
                    evt = (ident) & 0xffff;
                }
                apa_bpts.global.reset(evt, 0, run);
            }
        }
    }


    // if (m_save_deadarea && m_bee_dead.size()) {
    //     m_bee_dead.flush();
    //     m_sink.write(m_bee_dead);
    //     m_bee_dead.clear();
    // }
    if (m_save_deadarea) {

        // Flush individual patches
        for (auto& [apa, face_map] : m_bee_dead_patches) {
            for (auto& [face, patches] : face_map) {
                if (patches.size()) {
                    patches.flush();
                    m_sink.write(patches);
                    patches.clear();
                }
            }
        }
    }

    // Flush particle-flow mc trees
    for (auto& [name, tree] : m_bee_pf_trees) {
        if (!tree.empty()) {
            m_sink.write(tree);
            tree.reset();
        }
    }

    m_last_ident = ident;
}



// Helper function remains the same as in the previous response

void MultiAlgBlobClustering::fill_bee_points(const std::string& name, const Grouping& grouping)
{
    // std::cout << "Test: " << name << " " << grouping.wpids().size() << std::endl;
    
    if (m_bee_points.find(name) == m_bee_points.end()) {
        SPDLOG_LOGGER_WARN(log, "Bee points set '{}' not found, skipping", name);
        return;
    }
    
    auto& apa_bpts = m_bee_points[name];
    
    // Find the configuration for this name
    auto it = std::find_if(m_bee_points_configs.begin(), m_bee_points_configs.end(),
                          [&name](const BeePointsConfig& cfg) { return cfg.name == name; });
    
    if (it == m_bee_points_configs.end()) {
        SPDLOG_LOGGER_WARN(log, "Configuration for bee points set '{}' not found, skipping", name);
        return;
    }
    
    const auto& config = *it;
    
    // Reset RSE values for all points objects
    if (m_use_config_rse) {
        apa_bpts.global.rse(m_runNo, m_subRunNo, m_eventNo);
        for (auto& [apa, face_map] : apa_bpts.by_apa_face) {
            for (auto& [face, bpts] : face_map) {
                bpts.rse(m_runNo, m_subRunNo, m_eventNo);
            }
        }
    } else {
        // Use the default approach with ident
        int run = 0, evt = 0;
        if (m_last_ident > 0) {
            run = (m_last_ident >> 16) & 0x7fff;
            evt = (m_last_ident) & 0xffff;
        }
        apa_bpts.global.reset(evt, 0, run);
        for (auto& [anode_id, face_map] : apa_bpts.by_apa_face) {
            for (auto& [face, bpts] : face_map) {
                bpts.reset(evt, 0, run);
            }
        }
    }
    
    auto wpids = grouping.wpids();



    if (config.individual){ // fill in the individual APA
        for (auto wpid: wpids) {
            int apa = wpid.apa();
            int face = wpid.face();
            auto it = apa_bpts.by_apa_face.find(apa);
            if (it != apa_bpts.by_apa_face.end()) {
                auto it2 = it->second.find(face);
                if (it2 != it->second.end()) {
                    for (const auto* cluster : grouping.children()) {
                        fill_bee_points_from_cluster(it2->second, *cluster, config.pcname, config.coords, config.filter);
                    }
                }
            }
        }
    }else{ // fill in the global
        // std::cout << "Test: " << name << " " << grouping.wpids().size() << " " << grouping.nchildren() << std::endl;

        for (const auto* cluster : grouping.children()) {
            fill_bee_points_from_cluster(apa_bpts.global, *cluster, config.pcname, config.coords, config.filter);
        }
    }
}

// Fill bee points from PRGraph track trajectories
void MultiAlgBlobClustering::fill_bee_points_from_pr_graph(const std::string& name, const Grouping& grouping)
{
    if (m_bee_points.find(name) == m_bee_points.end()) {
        SPDLOG_LOGGER_WARN(log, "Bee points set '{}' not found for PR graph, skipping", name);
        return;
    }

    auto& apa_bpts = m_bee_points[name];

    // Find the configuration for this name
    auto it = std::find_if(m_bee_points_configs.begin(), m_bee_points_configs.end(),
                          [&name](const BeePointsConfig& cfg) { return cfg.name == name; });

    if (it == m_bee_points_configs.end()) {
        SPDLOG_LOGGER_WARN(log, "Configuration for bee points set '{}' not found, skipping", name);
        return;
    }

    const auto& config = *it;

    // Reset RSE values for all points objects
    if (m_use_config_rse) {
        apa_bpts.global.rse(m_runNo, m_subRunNo, m_eventNo);
        for (auto& [apa, face_map] : apa_bpts.by_apa_face) {
            for (auto& [face, bpts] : face_map) {
                bpts.rse(m_runNo, m_subRunNo, m_eventNo);
            }
        }
    } else {
        // Use the default approach with ident
        int run = 0, evt = 0;
        if (m_last_ident > 0) {
            run = (m_last_ident >> 16) & 0x7fff;
            evt = (m_last_ident) & 0xffff;
        }
        apa_bpts.global.reset(evt, 0, run);
        for (auto& [anode_id, face_map] : apa_bpts.by_apa_face) {
            for (auto& [face, bpts] : face_map) {
                bpts.reset(evt, 0, run);
            }
        }
    }

    // Get the PRGraph from the grouping
    auto pr_graph = grouping.get_pr_graph();
    if (!pr_graph) {
        SPDLOG_LOGGER_WARN(log, "No PR graph found in grouping for bee points set '{}'", name);
        return;
    }

    SPDLOG_LOGGER_DEBUG(log, "Filling bee points '{}' from PR graph with {} vertices and {} edges",
               name, boost::num_vertices(*pr_graph), boost::num_edges(*pr_graph));

    // Iterate through all segments (edges) in the graph
    int segment_count = 0;
    for (auto edge_desc : mir(boost::edges(*pr_graph))) {
        const auto& edge_bundle = (*pr_graph)[edge_desc];
        auto segment = edge_bundle.segment;

        if (!segment) continue;

        // Encode ID as cluster_id * 1000 + segment graph index for global uniqueness
        const int cluster_id = segment->cluster() ? segment->cluster()->get_cluster_id() : 0;
        const int encoded_id = cluster_id * 1000 + static_cast<int>(segment->get_graph_index());

        if (config.use_associate_points) {
            // --- shower_track mode: use associated points, charge from shower flags ---
            const bool is_shower =
                segment->flags_any(PR::SegmentFlags::kShowerTrajectory) ||
                segment->flags_any(PR::SegmentFlags::kShowerTopology) ||
                (segment->has_particle_info() && std::abs(segment->particle_info()->pdg()) == 11);
            const double charge = is_shower ? 15000.0 : 0.0;

            auto dpc = segment->dpcloud("associate_points");
            if (!dpc) {
                segment_count++;
                continue;
            }
            for (const auto& dp : dpc->get_points()) {
                WireCell::Point point(dp.x, dp.y, dp.z);
                apa_bpts.global.append(point, charge, cluster_id, encoded_id);
            }
        } else {
            // --- default mode: use fitted points with dQdx scale/offset ---
            const auto& fits = segment->fits();

            if (fits.empty()) {
                segment_count++;
                continue;
            }

            SPDLOG_LOGGER_DEBUG(log, "Segment {} has {} fitted points", encoded_id, fits.size());

            for (const auto& fit : fits) {
                if (!fit.valid()) continue;

                const auto& point = fit.point;
                double charge = fit.dQ;
                charge = charge * config.dQdx_scale + config.dQdx_offset;
                if (charge < 0) charge = 0;

                if (config.individual) {
                    if (fit.paf.first >= 0 && fit.paf.second >= 0) {
                        int apa = fit.paf.first;
                        int face = fit.paf.second;
                        auto it_apa = apa_bpts.by_apa_face.find(apa);
                        if (it_apa != apa_bpts.by_apa_face.end()) {
                            auto it_face = it_apa->second.find(face);
                            if (it_face != it_apa->second.end()) {
                                it_face->second.append(point, charge, cluster_id, encoded_id);
                            }
                        }
                    }
                } else {
                    apa_bpts.global.append(point, charge, cluster_id, encoded_id);
                }
            }
        }

        segment_count++;
    }

    SPDLOG_LOGGER_DEBUG(log, "Filled bee points '{}' from {} segments", name, segment_count);
}


void MultiAlgBlobClustering::fill_bee_vertices_from_pr_graph(const std::string& name, const Facade::Grouping& grouping)
{
    if (m_bee_points.find(name) == m_bee_points.end()) {
        SPDLOG_LOGGER_WARN(log, "Bee points set '{}' not found for graph vertices, skipping", name);
        return;
    }

    auto& apa_bpts = m_bee_points[name];

    // Reset RSE
    if (m_use_config_rse) {
        apa_bpts.global.rse(m_runNo, m_subRunNo, m_eventNo);
    } else {
        int run = 0, evt = 0;
        if (m_last_ident > 0) {
            run = (m_last_ident >> 16) & 0x7fff;
            evt = (m_last_ident) & 0xffff;
        }
        apa_bpts.global.reset(evt, 0, run);
    }

    auto pr_graph = grouping.get_pr_graph();
    if (!pr_graph) {
        SPDLOG_LOGGER_WARN(log, "No PR graph found in grouping for vertices bee set '{}'", name);
        return;
    }

    int vertex_count = 0;
    for (auto node_desc : PR::ordered_nodes(*pr_graph)) {
        const auto& node_bundle = (*pr_graph)[node_desc];
        auto vertex = node_bundle.vertex;
        if (!vertex) { ++vertex_count; continue; }

        // Encode ID as cluster_id * 1000 + vertex graph index for global uniqueness
        const int cluster_id = vertex->cluster() ? vertex->cluster()->get_cluster_id() : 0;
        const int encoded_id = cluster_id * 1000 + static_cast<int>(vertex->get_graph_index());

        const WireCell::Point& point = vertex->fit().point;
        const double charge = vertex->flags_any(PR::VertexFlags::kNeutrinoVertex) ? 15000.0 : 0.0;
        apa_bpts.global.append(point, charge, cluster_id, encoded_id);
        ++vertex_count;
    }

    SPDLOG_LOGGER_DEBUG(log, "Filled bee vertices '{}' from {} vertices", name, vertex_count);
}


// Helper: map PDG code to a short display name for the Bee particle tree.
static std::string pf_pdg_to_name(int pdg)
{
    switch (pdg) {
        case  11: return "e-";
        case -11: return "e+";
        case  13: return "mu-";
        case -13: return "mu+";
        case  22: return "gamma";
        case  211: return "pi+";
        case -211: return "pi-";
        case 2212: return "p";
        case 2112: return "n";
        case  321: return "K+";
        case -321: return "K-";
        default:   return "particle";
    }
}


// Hierarchical particle-flow tree matching the prototype "mc" Bee JSON format.
//
// Output is a bare JSON array of jsTree nodes, each node:
//   { "id":N, "text":"name  KE MeV",
//     "data":{"start":[x,y,z],"end":[x,y,z]},
//     "children":[...] }
// Leaf nodes (no children) additionally carry "icon":"jstree-file".
//
// Algorithm (mirrors prototype NeutrinoID::fill_particle_tree):
//   1. BFS from main_vertex through non-shower track segments, establishing
//      parent-child relationships among segments and recording which track
//      segment "arrived at" each vertex (vtx_incoming_seg map).
//   2. Disconnected track segments (not reachable from main_vertex) are
//      collected as additional root-level nodes so nothing is lost.
//   3. Each shower is attached under its parent track segment according to
//      start_connection_type:
//        type 1 (direct)   – nested directly as a leaf child.
//        type 2/3 (indirect/gap) – an intermediate pseudo-gamma leaf is
//                                  inserted first, then the shower under it.
//   4. Node IDs follow the prototype convention: cluster_id*1000 + seg_id.
void MultiAlgBlobClustering::fill_bee_pf_tree(const BeePFConfig& cfg,
                                               const Facade::Grouping& grouping,
                                               bool flag_print)
{
    flag_print = true;
    
    auto map_it = m_bee_pf_trees.find(cfg.name);
    if (map_it == m_bee_pf_trees.end()) {
        SPDLOG_LOGGER_WARN(log, "bee_pf tree storage '{}' not found", cfg.name);
        return;
    }
    auto& tree = map_it->second;

    auto pr_graph = grouping.get_pr_graph();
    if (!pr_graph) return;

    auto tf = grouping.get_track_fitting();
    if (!tf) return;

    auto main_vertex = tf->get_main_vertex();
    if (!main_vertex) {
        SPDLOG_LOGGER_DEBUG(log, "fill_bee_pf_tree '{}': no main vertex, skipping", cfg.name);
        return;
    }
    const auto* main_cluster = main_vertex->cluster();

    const auto& showers            = tf->get_showers();
    const auto& pi0_showers        = tf->get_pi0_showers();
    const auto& map_shower_pio_id  = tf->get_map_shower_pio_id();
    const auto& map_pio_id_mass    = tf->get_map_pio_id_mass();
    PR::IndexedSegmentSet conn4_skip_segs;

    // --- Vertex → node-descriptor map ---
    std::map<PR::VertexPtr, PR::node_descriptor, PR::VertexIndexCmp> vtx_to_nd;
    for (auto nd : PR::ordered_nodes(*pr_graph)) {
        if (auto vtx = (*pr_graph)[nd].vertex) vtx_to_nd[vtx] = nd;
    }
    if (!vtx_to_nd.count(main_vertex)) return;

    // --- Segment → shower map; collect all shower segments ---
    std::map<PR::SegmentPtr, PR::ShowerPtr, PR::SegmentIndexCmp> seg_to_shower;
    PR::IndexedSegmentSet shower_segs;
    for (const auto& shower : showers) {
        PR::IndexedVertexSet sv; PR::IndexedSegmentSet ss;
        shower->fill_sets(sv, ss, /*flag_exclude_start_segment=*/false);
        for (const auto& seg : ss) { seg_to_shower[seg] = shower; shower_segs.insert(seg); }

        auto [_, conn_type] = shower->get_start_vertex_and_type();
        if (conn_type == 4) {
            for (const auto& seg : ss) conn4_skip_segs.insert(seg);
            if (auto start_seg = shower->start_segment()) conn4_skip_segs.insert(start_seg);
        }
    }

    // --- BFS from main_vertex through track-only (non-shower) segments ---
    //   seg_parent[S] = nullptr  →  S is a root (direct daughter of neutrino vtx)
    //   seg_parent[S] = P        →  S is a child of P
    //   seg_endpoints[S] = {near_vtx, far_vtx}  (near = toward main_vertex)
    //   vtx_incoming_seg[V] = S  →  S is the segment that first reached V from main_vertex
    std::map<PR::SegmentPtr, PR::SegmentPtr,  PR::SegmentIndexCmp> seg_parent;
    std::map<PR::SegmentPtr, std::vector<PR::SegmentPtr>, PR::SegmentIndexCmp> seg_children;
    std::map<PR::SegmentPtr, std::pair<PR::VertexPtr,PR::VertexPtr>, PR::SegmentIndexCmp> seg_endpoints;
    std::map<PR::VertexPtr,  PR::SegmentPtr,  PR::VertexIndexCmp>  vtx_incoming_seg;

    PR::IndexedVertexSet visited_vtxs;
    PR::IndexedSegmentSet used_segs = shower_segs;   // pre-mark showers as visited

    visited_vtxs.insert(main_vertex);
    std::vector<std::pair<PR::VertexPtr, PR::SegmentPtr>> bfs_cur;

    // Seed BFS: all non-shower edges adjacent to main_vertex
    for (auto [eit, end] = boost::out_edges(vtx_to_nd.at(main_vertex), *pr_graph);
         eit != end; ++eit) {
        auto seg = (*pr_graph)[*eit].segment;
        if (!seg || used_segs.count(seg) || conn4_skip_segs.count(seg)) continue;
        auto far = PR::find_other_vertex(*pr_graph, seg, main_vertex);
        if (!far) continue;
        used_segs.insert(seg);
        seg_parent[seg]    = nullptr;
        seg_endpoints[seg] = {main_vertex, far};
        vtx_incoming_seg[far] = seg;
        bfs_cur.push_back({far, seg});
    }

    while (!bfs_cur.empty()) {
        std::vector<std::pair<PR::VertexPtr, PR::SegmentPtr>> bfs_next;
        for (auto& [cur_vtx, inc_seg] : bfs_cur) {
            if (visited_vtxs.count(cur_vtx)) continue;
            visited_vtxs.insert(cur_vtx);
            auto nd_it = vtx_to_nd.find(cur_vtx);
            if (nd_it == vtx_to_nd.end()) continue;
            for (auto [eit, end] = boost::out_edges(nd_it->second, *pr_graph);
                 eit != end; ++eit) {
                auto seg = (*pr_graph)[*eit].segment;
                if (!seg || used_segs.count(seg) || conn4_skip_segs.count(seg)) continue;
                auto far = PR::find_other_vertex(*pr_graph, seg, cur_vtx);
                if (!far) continue;
                used_segs.insert(seg);
                seg_parent[seg]    = inc_seg;
                seg_endpoints[seg] = {cur_vtx, far};
                seg_children[inc_seg].push_back(seg);
                if (!vtx_incoming_seg.count(far)) vtx_incoming_seg[far] = seg;
                bfs_next.push_back({far, seg});
            }
        }
        bfs_cur = std::move(bfs_next);
    }

    // // Log disconnected non-shower track segments (not added to particle flow).
    // for (auto edge_desc : mir(boost::edges(*pr_graph))) {
    //     auto seg = (*pr_graph)[edge_desc].segment;
    //     if (!seg || seg_to_shower.count(seg) || seg_parent.count(seg)) continue;
    //     auto [va, vb] = PR::find_vertices(*pr_graph, seg);
    //     const int cluster_id = seg->cluster() ? seg->cluster()->get_cluster_id() : -1;
    //     const int graph_idx  = static_cast<int>(seg->get_graph_index());
    //     const auto& fits = seg->fits();
    //     const double length  = fits.empty() ? -1.0
    //         : PR::walk_length(fits.begin(), fits.end(),
    //                           [](const PR::Fit& f) -> WireCell::Point { return f.point; }) / units::cm;
    //     std::string pi_name = "?";
    //     if (seg->has_particle_info()) {
    //         pi_name = seg->particle_info()->name();
    //     }
    //     std::cout << "[fill_bee_pf_tree] DISCONNECTED track seg"
    //               << "  cluster=" << cluster_id
    //               << "  graph_idx=" << graph_idx
    //               << "  encoded_id=" << (cluster_id >= 0 ? cluster_id * 1000 + graph_idx : graph_idx)
    //               << "  length_cm=" << std::fixed << std::setprecision(2) << length
    //               << "  particle=" << pi_name
    //               << "  has_va=" << (va ? 1 : 0) << " " << (va ? va->fit().point : WireCell::Point(0,0,0)) << " " <<  va->wcpt().point << " " << seg->fits().size()  << " " << seg->fits()[0].point << " " << seg->fits()[1].point
    //               << "  has_vb=" << (vb ? 1 : 0) << " " << (vb ? vb->fit().point : WireCell::Point(0,0,0)) << " " << vb->wcpt().point 
    //               << "\n";
    // }

    // --- Extend vtx_incoming_seg through shower vertex sets (mirrors prototype) ---
    // The prototype guarantees a shower's start_vtx is always picked from
    // (main-cluster vertices ∪ existing shower vertices) during clustering, so it is
    // always reachable at fill time.  We replicate that guarantee here by propagating
    // vtx_incoming_seg into every vertex belonging to an already-resolved shower, then
    // repeating to a fixed point so that showers nested inside other showers resolve too.
    //
    // Two flavors of "resolved":
    //   a) start_vtx == main_vertex (or in root_reachable_vtxs)
    //      → shower hangs from root; add its vertices to root_reachable_vtxs
    //   b) start_vtx in vtx_incoming_seg
    //      → shower hangs from a track segment; extend vtx_incoming_seg with its vertices
    PR::IndexedVertexSet root_reachable_vtxs;
    std::map<PR::VertexPtr, PR::ShowerPtr, PR::VertexIndexCmp> vtx_to_parent_shower;
    {
        bool any_added = true;
        while (any_added) {
            any_added = false;
            for (const auto& shower : showers) {
                auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
                if (conn_type == 4) continue;
                if (!start_vtx) continue;

                const bool at_main  = (start_vtx == main_vertex);
                const bool at_root  = root_reachable_vtxs.count(start_vtx) > 0;
                const bool at_track = vtx_incoming_seg.count(start_vtx) > 0;
                if (!at_main && !at_root && !at_track) continue;

                PR::SegmentPtr parent_seg = at_track ? vtx_incoming_seg.at(start_vtx) : nullptr;
                PR::ShowerPtr parent_shower = (at_main || at_root) ? shower : nullptr;

                PR::IndexedVertexSet sv; PR::IndexedSegmentSet ss;
                shower->fill_sets(sv, ss, /*flag_exclude_start_segment=*/false);
                for (const auto& vtx : sv) {
                    if (vtx == main_vertex) continue;
                    if (at_main || at_root) {
                        // hangs from root → add to root_reachable_vtxs
                        if (!root_reachable_vtxs.count(vtx)) {
                            root_reachable_vtxs.insert(vtx);
                            vtx_to_parent_shower[vtx] = parent_shower;
                            any_added = true;
                        }
                    } else {
                        // hangs from a track segment → extend vtx_incoming_seg
                        if (!vtx_incoming_seg.count(vtx) && !root_reachable_vtxs.count(vtx)) {
                            vtx_incoming_seg[vtx] = parent_seg;
                            any_added = true;
                        }
                    }
                }
            }
        }
    }

    // --- Attach each shower to its parent (track segment, other shower, or root) ---
    // type 1 (direct):    nested directly under parent as shower leaf
    // type 2/3 (indirect): a pseudo-gamma node is inserted between parent and shower
    // shower_parent_vtx[shower] = the connection vertex used by that shower
    using ShowerSegMap = std::map<PR::SegmentPtr,
                                  std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>>,
                                  PR::SegmentIndexCmp>;
    using ShowerShowerMap = std::map<PR::ShowerPtr,
                                     std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>>,
                                     PR::ShowerIndexCmp>;
    ShowerSegMap seg_direct_showers;    // seg → [(shower, conn_vtx)]
    ShowerSegMap seg_indirect_showers;  // seg → [(shower, conn_vtx)]
    ShowerShowerMap shower_direct_showers;    // shower → [(child_shower, conn_vtx)]
    ShowerShowerMap shower_indirect_showers;  // shower → [(child_shower, conn_vtx)]
    std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>> root_direct_showers;
    std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>> root_indirect_showers;

    for (const auto& shower : showers) {
        auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();

        // type 4 = "not clearly connected"; skip entirely (prototype behaviour)
        if (conn_type == 4) {
            if (flag_print) {
                auto start_seg = shower->start_segment();
                const auto* cl = start_seg ? start_seg->cluster() : nullptr;
                std::cout << "[fill_bee_pf_tree] SKIP shower (conn_type=4)"
                          << "  pdg=" << shower->get_particle_type()
                          << "  ke=" << shower->get_kine_best() / units::MeV << " MeV"
                          << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                          << "\n";
            }
            continue;
        }

        bool direct = (conn_type == 1);

        if (!start_vtx || start_vtx == main_vertex) {
            if (flag_print) {
                auto start_seg = shower->start_segment();
                const auto* cl = start_seg ? start_seg->cluster() : nullptr;
                WireCell::Point sp = shower->get_start_point();
                std::cout << "[fill_bee_pf_tree] ROOT shower"
                          << "  conn_type=" << conn_type
                          << "  reason=" << (!start_vtx ? "null_start_vtx" : "start_vtx==main_vertex")
                          << "  pdg=" << shower->get_particle_type()
                          << "  ke=" << shower->get_kine_best() / units::MeV << " MeV"
                          << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                          << "  start=(" << sp.x()/units::cm << "," << sp.y()/units::cm << "," << sp.z()/units::cm << ") cm"
                          << "\n";
            }
            auto& vec = direct ? root_direct_showers : root_indirect_showers;
            vec.push_back({shower, main_vertex});
        } else {
            auto it = vtx_incoming_seg.find(start_vtx);
            if (it != vtx_incoming_seg.end()) {
                if (flag_print) {
                    auto start_seg = shower->start_segment();
                    const auto* cl = start_seg ? start_seg->cluster() : nullptr;
                    WireCell::Point sp = shower->get_start_point();
                    const int parent_seg_id = (it->second->id() >= 0)
                                                ? it->second->id()
                                                : static_cast<int>(it->second->get_graph_index());
                    std::cout << "[fill_bee_pf_tree] SEGMENT-attached shower"
                              << "  conn_type=" << conn_type
                              << "  parent_seg=" << parent_seg_id
                              << "  pdg=" << shower->get_particle_type()
                              << "  ke=" << shower->get_kine_best() / units::MeV << " MeV"
                              << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                              << "  start=(" << sp.x()/units::cm << "," << sp.y()/units::cm << "," << sp.z()/units::cm << ") cm"
                              << "\n";
                }
                auto& mp = direct ? seg_direct_showers : seg_indirect_showers;
                mp[it->second].push_back({shower, start_vtx});
            } else {
                // start_vtx not reachable from main_vertex → check for parent shower
                if (root_reachable_vtxs.count(start_vtx)) {
                    // start_vtx is inside a root-level shower → attach to that parent shower
                    auto parent_shower_it = vtx_to_parent_shower.find(start_vtx);
                    if (parent_shower_it != vtx_to_parent_shower.end()) {
                        PR::ShowerPtr parent_shower = parent_shower_it->second;
                        if (flag_print) {
                            auto start_seg = shower->start_segment();
                            const auto* cl = start_seg ? start_seg->cluster() : nullptr;
                            WireCell::Point sp = shower->get_start_point();
                            std::cout << "[fill_bee_pf_tree] SHOWER-attached shower (via parent shower vtx)"
                                      << "  conn_type=" << conn_type
                                      << "  parent_shower_pdg=" << parent_shower->get_particle_type()
                                      << "  pdg=" << shower->get_particle_type()
                                      << "  ke=" << shower->get_kine_best() / units::MeV << " MeV"
                                      << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                                      << "  start=(" << sp.x()/units::cm << "," << sp.y()/units::cm << "," << sp.z()/units::cm << ") cm"
                                      << "\n";
                        }
                        auto& mp = direct ? shower_direct_showers : shower_indirect_showers;
                        mp[parent_shower].push_back({shower, start_vtx});
                    } else {
                        if (flag_print) {
                            auto start_seg = shower->start_segment();
                            const auto* cl = start_seg ? start_seg->cluster() : nullptr;
                            WireCell::Point sp = shower->get_start_point();
                            std::cout << "[fill_bee_pf_tree] ROOT shower (via root-reachable shower vtx, no parent found)"
                                      << "  conn_type=" << conn_type
                                      << "  pdg=" << shower->get_particle_type()
                                      << "  ke=" << shower->get_kine_best() / units::MeV << " MeV"
                                      << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                                      << "  start=(" << sp.x()/units::cm << "," << sp.y()/units::cm << "," << sp.z()/units::cm << ") cm"
                                      << "\n";
                        }
                        auto& vec = direct ? root_direct_showers : root_indirect_showers;
                        vec.push_back({shower, main_vertex});
                    }
                } else {
                    // start_vtx truly isolated from main_vertex → fallback to root
                    if (flag_print) {
                        auto start_seg = shower->start_segment();
                        const auto* cl = start_seg ? start_seg->cluster() : nullptr;
                        WireCell::Point sp = shower->get_start_point();
                        std::cout << "[fill_bee_pf_tree] ROOT shower (fallback: start_vtx not in BFS tree)"
                                  << "  conn_type=" << conn_type
                                  << "  pdg=" << shower->get_particle_type()
                                  << "  ke=" << shower->get_kine_best() / units::MeV << " MeV"
                                  << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                                  << "  start=(" << sp.x()/units::cm << "," << sp.y()/units::cm << "," << sp.z()/units::cm << ") cm"
                                  << "\n";
                    }
                    auto& vec = direct ? root_direct_showers : root_indirect_showers;
                    vec.push_back({shower, start_vtx});
                }
            }
        }
    }

    // --- Helpers ---
    auto get_vtx_pt = [](PR::VertexPtr v) -> WireCell::Point {
        return v->fit().valid() ? v->fit().point : v->wcpt().point;
    };

    // ID following prototype convention: cluster_id * 1000 + seg_id
    auto seg_display_id = [](PR::SegmentPtr seg) -> int {
        int sid = seg->id();
        if (sid < 0) sid = static_cast<int>(seg->get_graph_index());
        const auto* cl = seg->cluster();
        return cl ? cl->get_cluster_id() * 1000 + sid : sid;
    };

    int next_id = 1;  // fallback counter for nodes without a natural ID

    auto make_node = [&](int id,
                         const std::string& text,
                         const WireCell::Point& start,
                         const WireCell::Point& end) -> Configuration {
        Configuration node;
        node["id"]   = id;
        node["text"] = text;
        Configuration dj;
        dj["start"][0] = start.x() / units::cm;
        dj["start"][1] = start.y() / units::cm;
        dj["start"][2] = start.z() / units::cm;
        dj["end"][0]   = end.x() / units::cm;
        dj["end"][1]   = end.y() / units::cm;
        dj["end"][2]   = end.z() / units::cm;
        node["data"] = dj;
        node["children"] = Json::arrayValue;
        return node;
    };

    auto format_mev = [](double energy) -> std::string {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.2f", energy / units::MeV);
        return std::string(buf);
    };

    // Forward declare as std::function to handle mutual recursion with make_shower_leaf
    std::function<void(Configuration&, const std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>>&, 
                       const std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>>&, PR::VertexPtr)> append_showers;

    auto make_shower_leaf = [&](PR::ShowerPtr shower) -> Configuration {
        const int pdg = shower->get_particle_type();
        const std::string ke = format_mev(shower->get_kine_best());
        auto [svtx, sconn] = shower->get_start_vertex_and_type();
        auto start_seg = shower->start_segment();
        int id = start_seg ? seg_display_id(start_seg) : (next_id++);
        auto node = make_node(id,
                              pf_pdg_to_name(pdg) + "  " + ke + " MeV",
                              shower->get_start_point(), shower->get_end_point());
        if (flag_print) {
            const auto* cl = start_seg ? start_seg->cluster() : nullptr;
            std::cout << "[fill_bee_pf_tree] ADD shower-leaf"
                      << "  id=" << id
                      << "  pdg=" << pdg
                      << "  conn_type=" << sconn
                      << "  ke=" << ke << " MeV"
                      << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                      << "  has_start_vtx=" << (svtx ? 1 : 0)
                      << "\n";
        }
        
        // Append any showers attached to this shower
        static const std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>> empty_showers;
        auto di_it  = shower_direct_showers.find(shower);
        auto ind_it = shower_indirect_showers.find(shower);
        append_showers(node["children"],
                       di_it  != shower_direct_showers.end()   ? di_it->second  : empty_showers,
                       ind_it != shower_indirect_showers.end() ? ind_it->second : empty_showers,
                       svtx);
        
        if (node["children"].empty()) node["icon"] = "jstree-file";
        return node;
    };

    // Helper: insert one gamma pseudo-node + shower-leaf into a parent children array.
    auto append_gamma_shower = [&](Configuration& parent_children, PR::ShowerPtr sh, PR::VertexPtr conn_vtx) {
        PR::VertexPtr cv = conn_vtx;
        WireCell::Point gstart = cv ? get_vtx_pt(cv) : sh->get_start_point();
        WireCell::Point gend   = sh->get_start_point();
        const std::string gamma_ke = format_mev(sh->get_kine_best());
        const int gamma_id = next_id++;
        auto gamma = make_node(gamma_id, "gamma  " + gamma_ke + " MeV", gstart, gend);
        if (flag_print) {
            std::cout << "[fill_bee_pf_tree] ADD pseudo-gamma"
                      << "  id=" << gamma_id
                      << "  ke=" << gamma_ke << " MeV"
                      << "  child_pdg=" << sh->get_particle_type()
                      << "\n";
        }
        gamma["children"].append(make_shower_leaf(sh));
        if (gamma["children"].empty()) gamma["icon"] = "jstree-file";
        parent_children.append(gamma);
    };

    // Append all showers (direct + indirect via pseudo-gamma) into a children array,
    // given the connection vertex for the indirect case.
    // Pi0 showers are grouped by pi0_id and rendered as: pi0 node → gamma → shower_leaf.
    append_showers = [&](Configuration& children,
                               const std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>>& direct,
                               const std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>>& indirect,
                               PR::VertexPtr fallback_conn_vtx) {
        // --- Non-pi0 direct showers ---
        for (auto& [sh, _] : direct) {
            if (pi0_showers.count(sh)) continue;
            children.append(make_shower_leaf(sh));
        }
        // --- Non-pi0 indirect showers (pseudo-gamma) ---
        for (auto& [sh, conn_vtx] : indirect) {
            if (pi0_showers.count(sh)) continue;
            PR::VertexPtr cv = conn_vtx ? conn_vtx : fallback_conn_vtx;
            append_gamma_shower(children, sh, cv);
        }

        // --- Pi0 showers: group by pi0_id, emit one pi0 node per group ---
        // Collect from both direct and indirect lists.
        std::map<int, std::vector<std::pair<PR::ShowerPtr, PR::VertexPtr>>> pi0_groups;
        for (auto& [sh, cv] : direct)   if (pi0_showers.count(sh)) pi0_groups[map_shower_pio_id.at(sh)].push_back({sh, cv});
        for (auto& [sh, cv] : indirect) if (pi0_showers.count(sh)) pi0_groups[map_shower_pio_id.at(sh)].push_back({sh, cv});

        for (auto& [pi0_id, group] : pi0_groups) {
            auto mass_it = map_pio_id_mass.find(pi0_id);
            const double pi0_ke = (mass_it != map_pio_id_mass.end()) ? mass_it->second.first : 0.0;
            // Pi0 sits at the connection vertex (point particle: start == end)
            PR::VertexPtr conn_vtx = group[0].second ? group[0].second : fallback_conn_vtx;
            WireCell::Point pi0_pt = conn_vtx ? get_vtx_pt(conn_vtx) : WireCell::Point(0,0,0);
            const int pi0_node_id = next_id++;
            auto pi0_node = make_node(pi0_node_id, "pi0  " + format_mev(pi0_ke) + " MeV", pi0_pt, pi0_pt);
            if (flag_print) {
                std::cout << "[fill_bee_pf_tree] ADD pi0-node"
                          << "  id=" << pi0_node_id
                          << "  pi0_id=" << pi0_id
                          << "  ke=" << format_mev(pi0_ke) << " MeV"
                          << "  n_showers=" << group.size()
                          << "\n";
            }
            for (auto& [sh, cv] : group) {
                PR::VertexPtr gcv = cv ? cv : fallback_conn_vtx;
                append_gamma_shower(pi0_node["children"], sh, gcv);
            }
            if (pi0_node["children"].empty()) pi0_node["icon"] = "jstree-file";
            children.append(pi0_node);
        }
    };

    // Build the full JSON subtree for a track segment (recursive).
    std::function<Configuration(PR::SegmentPtr)> build_seg_node =
        [&](PR::SegmentPtr seg) -> Configuration {

        std::string pname = "particle";
        std::string ke_str = "0.00";
        if (seg->has_particle_info()) {
            auto pi = seg->particle_info();
            pname  = pi->name();
            ke_str = format_mev(pi->kinetic_energy());
        }

        auto ep_it = seg_endpoints.find(seg);
        WireCell::Point start_pt, end_pt;
        PR::VertexPtr far_vtx = nullptr;
        if (ep_it != seg_endpoints.end()) {
            start_pt = ep_it->second.first  ? get_vtx_pt(ep_it->second.first)  : WireCell::Point(0,0,0);
            end_pt   = ep_it->second.second ? get_vtx_pt(ep_it->second.second) : start_pt;
            far_vtx  = ep_it->second.second;
        }

        if (flag_print) {
            const auto* seg_cluster = seg->cluster();
            const bool in_main_cluster = (main_cluster && seg_cluster == main_cluster);
            const bool is_shower_seg = (seg_to_shower.count(seg) > 0);
            auto pit = seg_parent.find(seg);
            auto parent_seg_dbg = (pit != seg_parent.end()) ? pit->second : nullptr;
            const std::string parent_text = parent_seg_dbg ? std::to_string(seg_display_id(parent_seg_dbg)) : "ROOT";
            std::cout << "[fill_bee_pf_tree] ADD track-node"
                      << "  seg=" << seg_display_id(seg)
                      << "  parent=" << parent_text
                      << "  name=" << pname
                      << "  ke=" << ke_str << " MeV"
                      << "  cluster=" << (seg_cluster ? std::to_string(seg_cluster->get_cluster_id()) : "?")
                      << "  in_main_cluster=" << (in_main_cluster ? 1 : 0)
                      << "  is_shower_seg=" << (is_shower_seg ? 1 : 0)
                      << "  start=(" << start_pt.x()/units::cm << "," << start_pt.y()/units::cm << "," << start_pt.z()/units::cm << ") cm"
                      << "  end=(" << end_pt.x()/units::cm << "," << end_pt.y()/units::cm << "," << end_pt.z()/units::cm << ") cm"
                      << "\n";
        }

        auto node = make_node(seg_display_id(seg),
                              pname + "  " + ke_str + " MeV",
                              start_pt, end_pt);

        // Attach showers connected at the far-end vertex of this segment
        static const std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>> empty_showers;
        auto di_it  = seg_direct_showers.find(seg);
        auto ind_it = seg_indirect_showers.find(seg);
        append_showers(node["children"],
                       di_it  != seg_direct_showers.end()   ? di_it->second  : empty_showers,
                       ind_it != seg_indirect_showers.end() ? ind_it->second : empty_showers,
                       far_vtx);

        // Recurse into child track segments
        auto ch_it = seg_children.find(seg);
        if (ch_it != seg_children.end()) {
            for (const auto& child : ch_it->second) {
                if (conn4_skip_segs.count(child)) continue;
                node["children"].append(build_seg_node(child));
            }
        }

        if (node["children"].empty()) node["icon"] = "jstree-file";
        return node;
    };

    // --- Assemble top-level particles array ---
    Configuration particles = Json::arrayValue;

    // Showers attached directly to the neutrino vertex
    append_showers(particles, root_direct_showers, root_indirect_showers, main_vertex);

    // Root track segments (direct daughters of the neutrino vertex).
    // Disconnected segments (orphaned fragments unreachable from main_vertex)
    // are now skipped entirely to avoid adding zero-energy orphaned particles.
    for (auto& [seg, parent] : seg_parent) {
        if (parent != nullptr) continue;   // skip non-roots
        if (conn4_skip_segs.count(seg)) continue;
        particles.append(build_seg_node(seg));
    }

    tree.set_particles(particles);
    SPDLOG_LOGGER_DEBUG(log, "fill_bee_pf_tree '{}': {} top-level particles",
                        cfg.name, particles.size());
}


// Helper function to fill bee points from a single cluster
void MultiAlgBlobClustering::fill_bee_points_from_cluster(
    Bee::Points& bpts, const Cluster& cluster, 
    const std::string& pcname, const std::vector<std::string>& coords, int filter)
{
    int clid = cluster.get_cluster_id(); //bpts.back_cluster_id() + 1;

    // std::cout << "Test: " << bpts.size() << " " << bpts.back_cluster_id() << " " <<  clid << std::endl;

    if (pcname == "steiner_pc"){
        // Export Steiner points ... 
        // std::cout << "Exporting Steiner points for cluster ID: " << clid << " " << cluster.nchildren() << std::endl;

        auto& steiner_pc = cluster.get_pc(pcname);
        if (steiner_pc.empty()) {
            return;
        }
        // Get coordinate arrays from the point cloud
        const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
        const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>(); 
        const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
        const auto& flag_steiner_terminal = steiner_pc.get("flag_steiner_terminal")->elements<int>();

        // std::cout << "Steiner Test: " << x_coords.size() << " " << y_coords.size() << " " << z_coords.size() << std::endl;

         for (size_t i = 0; i < x_coords.size(); ++i) {
            // Create point from steiner point cloud
            Point vtx(x_coords[i], y_coords[i], z_coords[i]);

            // Get the point index from the default scope
            auto point_index = cluster.get_closest_point_index(vtx);
            
            auto charge_result = cluster.calc_charge_wcp(point_index, 4000, true);
            double point_charge = charge_result.second; // Extract the charge value from the pair

            if (flag_steiner_terminal[i]) {
                bpts.append(Point(x_coords[i], y_coords[i], z_coords[i]), point_charge, clid, 1);  // terminals  ... 
            }else{
                bpts.append(Point(x_coords[i], y_coords[i], z_coords[i]), point_charge, clid, 0); // non-terminals ...
            }
         }


    }else{
        // Get the scope
        Scope scope = {pcname, coords};
        
        auto filter_scope = cluster.get_scope_filter(scope);

        // std::cout << "Test: " << cluster.get_cluster_id() << " " << clid << " " << scope << " " << filter_scope << std::endl;

        bool use_scope = true;
        if (filter == 1) {
            use_scope = filter_scope;
        }
        else if (filter == 0) {
            use_scope = true; // ignore filter_scope, always true
        }
        else if (filter == -1) {
            use_scope = !filter_scope;
        }

        if (use_scope) {
            // Access the points through the cluster's scoped view
            const WireCell::PointCloud::Tree::ScopedView<double>& sv = cluster.sv<double>(scope);
            const auto& spcs = sv.pcs();

            // Dead-channel threshold: uncertainty > 1e10 flags a dead wire
            // (matches PointTreeBuilding m_dead_threshold and Facade_Cluster::is_wire_dead).
            const double dead_threshold = 1e10;

            // For each scoped point cloud (one per blob), compute per-point wire charge.
            // Each point gets q = mean(Q_U, Q_V, Q_W) over the non-dead planes for the
            // specific wires that define this 3D point — matching the prototype formula.
            // The global cluster point index is obtained via a KD-tree nearest-neighbour
            // lookup (exact match: every scoped-view point is also in the cluster PC).
            for (size_t spc_idx = 0; spc_idx < spcs.size(); ++spc_idx) {
                const auto& spc = spcs[spc_idx];
                auto x = spc.get().get(coords[0])->elements<double>();
                auto y = spc.get().get(coords[1])->elements<double>();
                auto z = spc.get().get(coords[2])->elements<double>();

                const size_t size = x.size();
                for (size_t ind = 0; ind < size; ++ind) {
                    // Resolve global point index via spatial lookup.
                    // charge_value() caches all per-plane charge vectors on first call,
                    // so subsequent lookups are O(1) vector reads.
                    const WireCell::Point pt(x[ind], y[ind], z[ind]);
                    const size_t pt_idx = cluster.get_closest_point_index(pt);

                    // Per-plane charge mean (prototype formula), excluding dead planes.
                    double sum = 0.0;
                    int nplanes = 0;
                    for (int plane = 0; plane < 3; ++plane) {
                        if (!cluster.is_wire_dead(pt_idx, plane, dead_threshold)) {
                            sum += cluster.charge_value(pt_idx, plane);
                            ++nplanes;
                        }
                    }
                    const double point_charge = (nplanes > 0) ? sum / nplanes : 0.0;

                    bpts.append(pt, point_charge, clid, clid);
                }
            }

        }
    }

}




void MultiAlgBlobClustering::fill_bee_patches_from_grouping(
    const WireCell::Clus::Facade::Grouping& grouping)
{
    // auto wpids = grouping.wpids();

    // For each cluster in the grouping
    for (const auto* cluster : grouping.children()) {
        // Get the wpids to determine which APA and face this cluster belongs to

        fill_bee_patches_from_cluster(*cluster);

        
        // if (!wpids.empty()) {
        //     // Store patches by APA and face
        //     for (auto wpid : wpids) {
        //         int apa = wpid.apa();
        //         int face = wpid.face();
        //        
        //     }
        // } 
    }
}


// Helper function to fill patches from a single cluster
void MultiAlgBlobClustering::fill_bee_patches_from_cluster(
    const WireCell::Clus::Facade::Cluster& cluster)
{
    int first_slice = -1;
    
    // Get the underlying node that contains this cluster
    const auto* cluster_node = cluster.node();
    if (!cluster_node) {
        SPDLOG_LOGGER_WARN(log, "Cannot access node for cluster");
        return;
    }
    
    // Iterate through child nodes (blobs)
    for (const auto* bnode : cluster_node->children()) {
        auto wpid = bnode->value.facade<Blob>()->wpid();
        int apa = wpid.apa();
        int face = wpid.face();


        auto it_apa = m_bee_dead_patches.find(apa);
        if (it_apa != m_bee_dead_patches.end()) {
            auto it_face = it_apa->second.find(face);
            if (it_face != it_apa->second.end()) {
                auto & patches = it_face->second;

                // Access the local point clouds in the node
                const auto& lpcs = bnode->value.local_pcs();
                
                // Get the scalar PC to find the slice index
                if (lpcs.find("scalar") == lpcs.end()) {
                    continue;  // Skip if no scalar PC
                }
                const auto& pc_scalar = lpcs.at("scalar");
                
                // Get slice_index_min
                if (!pc_scalar.get("slice_index_min")) {
                    continue;  // Skip if no slice_index_min
                }
                int slice_index_min = pc_scalar.get("slice_index_min")->elements<int>()[0];
                
                // Set first_slice if not already set
                if (first_slice < 0) {
                    first_slice = slice_index_min;
                }
                
                // Skip blobs not on the first slice
                if (slice_index_min != first_slice) continue;
                
                // Access the corner point cloud
                if (lpcs.find("corner") == lpcs.end()) {
                    continue;  // Skip if no corner PC
                }
                const auto& pc_corner = lpcs.at("corner");
                
                // Get y and z coordinates
                if (!pc_corner.get("y") || !pc_corner.get("z")) {
                    continue;  // Skip if missing y or z
                }
                const auto& y = pc_corner.get("y")->elements<double>();
                const auto& z = pc_corner.get("z")->elements<double>();
                
                // Add to patches
                patches.append(y.begin(), y.end(), z.begin(), z.end());
            }
        }        
    }
}


struct Perf {
    bool enable;
    Log::logptr_t log;
    ExecMon em;

    Perf(bool e, Log::logptr_t l, const std::string& t = "starting MultiAlgBlobClustering")
      : enable(e)
      , log(l)
      , em(t)
    {
    }

    ~Perf()
    {
        if (!enable) return;
        SPDLOG_LOGGER_DEBUG(log, "MultiAlgBlobClustering performance summary:\n{}", em.summary());
    }

    void operator()(const std::string& ctx)
    {
        if (!enable) return;
        em(ctx);
    }

    void dump(const std::string& ctx, const Ensemble& ensemble, bool shallow = true, bool mon = true)
    {
        if (!enable) return;
        if (mon) (*this)(ctx);

        SPDLOG_LOGGER_DEBUG(log, "{} ensemble with {} groupings:", ctx, ensemble.nchildren());

        for (const auto* grouping : ensemble.children()) {

            {
                auto name = grouping->get_name();
                size_t npoints_total = 0;
                size_t nzero = 0;
                size_t count = 0;
                for (const auto* cluster : grouping->children()) {
                    int n = cluster->npoints();
                    if (n == 0) {
                        ++nzero;
                    }
                    npoints_total += n;
                    // SPDLOG_LOGGER_DEBUG(log, "loaded cluster {} with {} points out of {}", count, n, npoints_total);
                    ++count;
                    // std::cout << "Xin: " << name << " loaded cluster " << count << " with " << n << "points and " << cluster->nchildren() << "blobs" << std::endl;
                }

                

                SPDLOG_LOGGER_DEBUG(log, "\tgrouping \"{}\": {}, {} points and {} clusters with no points",
                           name, *grouping, npoints_total, nzero);
                (void)count;
            }

            if (shallow) continue;

            auto children = grouping->children();  // copy
            sort_clusters(children);
            size_t count = 0;
            for (const auto* cluster : children) {
                bool sane = cluster->sanity(log);
                SPDLOG_LOGGER_DEBUG(log, "\t\tcluster {} {} sane:{}", count++, *cluster, sane);
            }
        }
    }
};


Grouping& MultiAlgBlobClustering::load_grouping(
    Ensemble& ensemble,
    const std::string& name,
    const std::string& path,
    const ITensorSet::pointer ints)
{
    const auto& tens = *ints->tensors();
    try {
        ensemble.add_grouping_node(name, as_pctree(tens, path));
    }
    catch (WireCell::KeyError& err) {
        SPDLOG_LOGGER_WARN(log, "No pc-tree at tensor datapath {}, making empty", path);
        ensemble.make_grouping(name);
    }
        
    Grouping* grouping = ensemble.with_name(name).at(0);
    if (!grouping) {
        raise<KeyError>("failed to make grouping node %s at %s", name, path);
    }

    grouping->enumerate_idents();
    grouping->set_anodes(m_anodes);
    grouping->set_detector_volumes(m_dv);
    return *grouping;
}

bool MultiAlgBlobClustering::operator()(const input_pointer& ints, output_pointer& outts)
{
    outts = nullptr;
    if (!ints) {
        flush();
        SPDLOG_LOGGER_DEBUG(log, "EOS at call {}", m_count++);
        return true;
    }

    Perf perf{m_perf, log};

    const int ident = ints->ident();
    SPDLOG_LOGGER_DEBUG(log, "loading tensor set ident={} (last={})", ident, m_last_ident);
    if (m_last_ident < 0) {     // first time.
        if (m_use_config_rse) {
            // Set RSE in the sink
            m_sink.set_rse(m_runNo, m_subRunNo, m_eventNo);
        }
        // Use default behavior
        // reset_bee(ident, m_bee_img);
        // reset_bee(ident, m_bee_ld);
        m_last_ident = ident;
    }
    else if (m_last_ident != ident) {
        flush(ident);
        if (m_use_config_rse) {
            // Update event number for next event
            m_eventNo++;
            // Update RSE in sink
            m_sink.set_rse(m_runNo, m_subRunNo, m_eventNo);
        }
    }
    // else do nothing when ident is unchanged.


    Points::node_t root;
    Ensemble& ensemble = *root.value.facade<Ensemble>();

    for (const auto& gname : m_groupings) {
        const auto datapath = inpath(gname, ident);
        load_grouping(ensemble, gname, datapath, ints);
        perf.dump("loaded " + gname, ensemble);
    }    

    if (m_save_deadarea) {
        auto gs = ensemble.with_name("live");
        if (gs.size()) {
            // Fill patches from the dead grouping
            fill_bee_patches_from_grouping(*gs[0]);
            perf("dump dead regions to bee");
        }
    }

    perf.dump("pre clustering", ensemble);

    for (const auto& config : m_bee_points_configs) {
        if (config.name != "img") {
            continue;
        }
        auto gs = ensemble.with_name("live");
        if (gs.empty()) {
            continue;
        }
        fill_bee_points(config.name, *gs[0]);
    }

    perf.dump("start clustering", ensemble);

    // THE MAIN LOOP
    for (const auto& cmeth : m_pipeline) {
        cmeth.meth->visit(ensemble);
        perf.dump(cmeth.name, ensemble);

        for (auto* grouping : ensemble.children()) {
            grouping->enumerate_idents(m_clusters_id_order);
        }

        // Dump bee points right after specific visitor runs
        for (const auto& config : m_bee_points_configs) {
            if (config.name == "img") continue;
            if (config.visitor.empty() || config.visitor != cmeth.name) continue;

            auto gs = ensemble.with_name(config.grouping);
            if (gs.empty()) {
                continue;
            }

            // Check if this visitor produced a PRGraph that we should save
            auto pr_graph = gs[0]->get_pr_graph();

            // std::cout << "Test: Visitor: " << cmeth.name << " Grouping: " << config.grouping << " " << pr_graph << std::endl;

            if (pr_graph) {
                if (config.use_graph_vertices) {
                    fill_bee_vertices_from_pr_graph(config.name, *gs[0]);
                } else {
                    // Fill bee points from PRGraph (for track trajectories)
                    fill_bee_points_from_pr_graph(config.name, *gs[0]);
                }
                // std::cout << "Filled bee points from PR graph for visitor: " << cmeth.name << " grouping: " << config.grouping << std::endl;
            } else {
                // Fill bee points from clusters normally
                fill_bee_points(config.name, *gs[0]);
                // std::cout << "Filled bee points from clusters for visitor: " << cmeth.name << " grouping: " << config.grouping << std::endl;
            }
        }

        // Particle-flow dump triggered by the same visitor
        for (const auto& pf_cfg : m_bee_pf_configs) {
            if (pf_cfg.visitor.empty() || pf_cfg.visitor != cmeth.name) continue;
            auto pf_gs = ensemble.with_name(pf_cfg.grouping);
            if (pf_gs.empty()) continue;
            const auto& pf_grouping = *pf_gs[0];
            auto tf = pf_grouping.get_track_fitting();
            if (!tf) continue;
            fill_bee_pf_tree(pf_cfg, pf_grouping);
        }
    }

    //
    // At this point, the ensemble may have more or fewer groupings just "live"
    // and "dead" including no groupings at all.  But for now, we assume the
    // original "live" and "dead" still exist and with their original facades.
    // Famous last words....
    //
    

    // Fill all configured bee points sets (except those with visitor-specific handling)
    for (const auto& config : m_bee_points_configs) {
        if(config.name == "img") continue;

        // Skip configs with visitor specified - they were already handled in the visitor loop
        if (!config.visitor.empty()) continue;

        auto gs = ensemble.with_name(config.grouping);
        if (gs.empty()) {
            continue;
        }
        fill_bee_points(config.name, *gs[0]);

    }
    perf("dump live clusters to bee");

    if (m_grouping2file_prefix.size()) {
        std::string fname = String::format("%s-%d.npz", m_grouping2file_prefix, m_count);
        auto live = ensemble.with_name("live");
        grouping2file(*live[0], fname);
    }
    auto grouping_names = ensemble.names();

    if (m_dump_json) {
        for (const auto& name : grouping_names) {
            auto gs = ensemble.with_name(name);
            Persist::dump(String::format("%s-summary-%d.json", name, ident),
                          json_summary(*gs[0]), true);
        }
    }

    SPDLOG_LOGGER_DEBUG(log, "Produce pctrees with {} groupings", grouping_names.size());
    
    ITensor::vector outtens;
    for (const auto& name : grouping_names) {

        // This next bit may look a little weird and it is so some explanation
        // is warranted.  Originally, we had disembodied "root" grouping nodes,
        // live and dead.  To clean up the clustering api we added the
        // "ensemble" as root node with children consisting of grouping nodes.
        // At the time of writing, the as_tensors() does not like serializing
        // non-root nodes I do not want to debug right now.  And, I do not want
        // the "ensemble" concept to leak out from the MABC+clustering context.
        // So, I remove each grouping child node from the ensemble prior to
        // serializing.  The remove gives an auto_ptr so the node is destructed
        // as this loop progresses.
        auto gs = ensemble.with_name(name);
        auto& grouping = *gs[0];
        auto node = ensemble.remove_child(grouping);
        auto tens = as_tensors(*node, outpath(name, ident));
        outtens.insert(outtens.end(), tens.begin(), tens.end());
        SPDLOG_LOGGER_DEBUG(log, "Produce {} tensors for grouping {}", tens.size(), name);
    }
    outts = as_tensorset(outtens, ident);

    perf("done");

    return true;
}
