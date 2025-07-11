#include "CreateSteinerGraph.h"
#include "SteinerGrapher.h"

#include "WireCellClus/Graphs.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellClus/ClusteringFuncs.h"


WIRECELL_FACTORY(CreateSteinerGraph, WireCell::Clus::Steiner::CreateSteinerGraph,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)


using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

Steiner::CreateSteinerGraph::CreateSteinerGraph()
    : Aux::Logger("CreateSteinerGraph", "clus")
{
}
Steiner::CreateSteinerGraph::~CreateSteinerGraph()
{
}


void Steiner::CreateSteinerGraph::configure(const WireCell::Configuration& cfg)
{
    m_grouping_name = get(cfg, "grouping", m_grouping_name);
    m_graph_name = get(cfg, "graph", m_graph_name);
    m_replace = get(cfg, "replace", m_replace);

    NeedDV::configure(cfg);
    NeedPCTS::configure(cfg);

    m_grapher_config.dv = m_dv;
    m_grapher_config.pcts = m_pcts;
    const std::string retiler_tn = get<std::string>(cfg, "retiler", "RetileCluster");
    m_grapher_config.retile = Factory::find_tn<IPCTreeMutate>(retiler_tn);
}

Configuration Steiner::CreateSteinerGraph::default_configuration() const
{
    Configuration cfg;
    // Build the Steiner graph for clusters in this grouping.
    cfg["grouping"] = m_grouping_name;
    // Name of the resulting graph on the cluster
    cfg["graph"] = m_graph_name;
    // If true, replace any pre-existing graph with that name, else do
    // nothing if one already exists.
    cfg["replace"] = m_replace;

    return cfg;
}


void Steiner::CreateSteinerGraph::visit(Ensemble& ensemble) const
{
    auto& grouping = *ensemble.with_name(m_grouping_name).at(0);
    
    // Container to hold clusters after the initial filter
    std::vector<Cluster*> filtered_clusters;

    for (auto* cluster : grouping.children()) {
        if (cluster->get_flag(Flags::beam_flash)){
            filtered_clusters.push_back(cluster);


        }
    }

   
    for (auto* cluster : filtered_clusters) {
        // separate the clusters into separated pieces ...
        auto cc =cluster->get_pcarray("isolated", "perblob");
        // convert span to vector
        std::vector<int> cc_vec(cc.begin(), cc.end());
        // std::cout << "Xin: " << cluster->ident() << " " << cc_vec.size() << std::endl;
        std::cout << "Xin: " << cluster->get_flash().time()/units::us << " " << cluster->nchildren() << " " << cluster->npoints() <<  " " << std::endl;
        
        if (cc_vec.size() < 2) continue;
        auto splits = grouping.separate(cluster, cc_vec);

         // Apply the scope filter settings to all new clusters
        for (auto& [id, new_cluster] : splits) {
            // Store the split/group ID as a scalar
            new_cluster->set_scalar<int>("split_id", id);
            // Optionally also store the original parent's ident
            new_cluster->set_scalar<int>("parent_ident", cluster->ident());
            std::cout << "Xin1: " << new_cluster->get_flash().time()/units::us << " " << new_cluster->nchildren() << " " << new_cluster->npoints() <<  " " << std::endl;
        }

        std::cout << "Xin1: " << cluster->get_flash().time()/units::us << " " << cluster->nchildren() << " " << cluster->npoints() <<  " " << std::endl;


        Steiner::Grapher sg(*cluster, m_grapher_config, log);
        bool already = cluster->has_graph(m_graph_name);
        if (already || m_replace) {
            auto cell_points_map = sg.form_cell_points_map();
            auto& graph = sg.get_graph("basic_pid");

            std::cout << "Xin2: " << cell_points_map.size() << " Graph vertices: " << boost::num_vertices(graph) << ", edges: " << boost::num_edges(graph) << std::endl;
            
            // sg.establish_same_blob_steiner_edges("basic_pid", true);

            // std::cout << "Xin2: " << cell_points_map.size() << " Graph vertices: " << boost::num_vertices(graph) << ", edges: " << boost::num_edges(graph) << std::endl;
            // sg.remove_same_blob_steiner_edges("basic_pid");
            // std::cout << "Xin2: " << cell_points_map.size() << " Graph vertices: " << boost::num_vertices(graph) << ", edges: " << boost::num_edges(graph) << std::endl;

            // Test Steiner_Graph ...
            

            auto pair_points = cluster->get_two_boundary_wcps();
            std::cout << "Xin3: " << pair_points.first.x() << " " 
                    << pair_points.first.y() << " " 
                    << pair_points.first.z() << " | "
                    << pair_points.second.x() << " " 
                    << pair_points.second.y() << " " 
                    << pair_points.second.z() << std::endl;
            auto first_index  =   cluster->get_closest_point_index(pair_points.first);
            auto second_index =   cluster->get_closest_point_index(pair_points.second);
            std::cout << "Xin3: " << first_index << " " << second_index << std::endl;
            std::vector<size_t> path_point_indices = cluster->graph_algorithms("basic_pid").shortest_path(first_index, second_index);
            std::cout << "Xin3: " << path_point_indices.size() << std::endl;
            // for (const auto& idx : path_point_indices) {
            //     auto point = cluster->point3d(idx);
            //     std::cout << "Xin4: " << point.x() << " " << point.y() << " " << point.z() << std::endl;
            // }

            sg.create_steiner_tree(cluster, path_point_indices, "basic_pid", "steiner_graph", true, "steiner_pc");
            const auto& steiner_point_cloud = sg.get_point_cloud("steiner_pc");
            const auto& steiner_graph = sg.get_graph("steiner_graph");
            auto& flag_terminals = sg.get_flag_steiner_terminal();
            size_t num_true_terminals = std::count(flag_terminals.begin(), flag_terminals.end(), true);
            const auto& new_to_old = sg.get_new_to_old_mapping();

            std::cout << "Xin2: " << cell_points_map.size() << " Graph vertices: " << boost::num_vertices(steiner_graph) << ", edges: " << boost::num_edges(steiner_graph) << " " << steiner_point_cloud.size_major()  <<  " " << num_true_terminals << std::endl;
            
            auto pair_idx = cluster->get_two_boundary_steiner_graph_idx("steiner_graph", "steiner_pc", false);
            std::cout << "Xin3: " << pair_idx.first << " " << pair_idx.second << std::endl;

            auto kd_results = cluster->kd_steiner_knn(1, pair_points.first);
            auto kd_points = cluster->kd_steiner_points(kd_results);
            std::cout << "Xin4: " << kd_points.size() << " " << *kd_points.begin() << std::endl;

            // auto edge_weight_map = get(boost::edge_weight, steiner_graph);
            // for (auto edge_it = boost::edges(steiner_graph); edge_it.first != edge_it.second; ++edge_it.first) {
            //     auto edge = *edge_it.first;
            //     auto src = boost::source(edge, steiner_graph);
            //     auto tgt = boost::target(edge, steiner_graph);

            //     // Use the new-to-old mapping to get original vertex IDs
                
            //     auto src_id = new_to_old.at(src);
            //     auto tgt_id = new_to_old.at(tgt);

            //     // Get the edge weight using the proper accessor
            //     auto weight = edge_weight_map[edge];

            //     std::cout << "Edge from vertex " << cluster->point3d(src_id) << " to " << cluster->point3d(tgt_id) << " with weight " << weight << " " << flag_terminals[src] << " " << flag_terminals[tgt] << std::endl;
            // }

            // for (const auto& [cell, points] : cell_points_map) {
            //     // std::cout << "Xin2 Cell: " << cell->slice_index_min() << " " << cell->u_wire_index_min() << " " << cell->v_wire_index_min() << " " << cell->w_wire_index_min() << " has " << points.size() << " points." << std::endl;
            //     // for (const auto& point : points) {
            //         // auto info = cluster->calc_charge_wcp(point);
            //         // std::cout << "Xin2 Point: " << point << " " << info.first << " " << info.second << std::endl;
            //     // }
            //     std::vector<const Blob*> single_blob = {cell};
            //     auto blob_peaks = sg.find_peak_point_indices(single_blob, "basic_pid", true);
            //     std::cout << "Xin2: " << cell->slice_index_min() << " " << cell->u_wire_index_min() << " " << cell->v_wire_index_min() << " " << cell->w_wire_index_min()  << "  " << points.size() << "  " << blob_peaks.size() <<std::endl;

            // }
            // auto steiner_terminals = sg.find_steiner_terminals("basic_pid");
            // std::cout << "Xin3: " << steiner_terminals.size() << std::endl;
            // auto extrem_points = cluster->get_extreme_wcps();
            // std::cout << "Xin4: " << extrem_points.size() << std::endl;
            // for (const auto& pts : extrem_points) {
            //     for (const auto& pt : pts) {
            //         std::cout << "Extreme point: ("
            //                   << pt.x() << ", "
            //                   << pt.y() << ", "
            //                   << pt.z() << ")" << std::endl;
            //     }
            // }

            // auto extreme_boundary_points = cluster->get_two_boundary_wcps(1);
            // std::cout << "Xin3: " << extreme_boundary_points.first.x() << " " 
            //           << extreme_boundary_points.first.y() << " " 
            //           << extreme_boundary_points.first.z() << " | "
            //           << extreme_boundary_points.second.x() << " " 
            //           << extreme_boundary_points.second.y() << " " 
            //           << extreme_boundary_points.second.z() << std::endl;

            // for (const auto& [cell, points] : cell_points_map) {
            //     auto total_charge = cell->estimate_total_charge();
            //     auto min_charge = cell->estimate_minimum_charge();

            //     std::cout << "Xin2 Cell: " << cell->slice_index_min() << " " << cell->u_wire_index_min() << " " << cell->v_wire_index_min() << " " << cell->w_wire_index_min() 
            //               << " has " << points.size() << " points, total charge: " << total_charge 
            //               << ", min charge: " << min_charge << " " << cell->get_wire_charge(0, cell->u_wire_index_min()) << " " << cell->get_wire_charge_error(0, cell->u_wire_index_min()) 
            //               << std::endl;
            //     // cell->check_dead_wire_consistency();
            // }

            //  std::cout << grouping.get_dead_winds1(0,0,0).size() << " " << grouping.get_dead_winds(0,0,0).size() << " " << grouping.get_dead_winds1(0,0,1).size() << " " << grouping.get_dead_winds(0,0,1).size() << " " << grouping.get_dead_winds1(0,0,2).size() << " " << grouping.get_dead_winds(0,0,2).size() << std::endl;


            //    auto gr = sg.create_steiner_graph();
           // Do we do any tests, eg on num_vertices()?
           //cluster->give_graph(m_graph_name, std::move(gr));
        }

            
            
    }
}
