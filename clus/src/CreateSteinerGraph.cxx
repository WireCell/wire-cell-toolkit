#include "CreateSteinerGraph.h"
#include "SteinerGrapher.h"

#include "WireCellClus/Graphs.h"
#include <chrono>
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

    m_perf = get(cfg, "perf", m_perf);
    m_grapher_config.dv = m_dv;
    m_grapher_config.pcts = m_pcts;
    m_grapher_config.perf = m_perf; // propagate perf flag into Grapher instances
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
    // If true, print per-step timing to stdout.
    cfg["perf"] = m_perf;

    return cfg;
}


void Steiner::CreateSteinerGraph::visit(Ensemble& ensemble) const
{
    using Clock = std::chrono::steady_clock;
    using MS = std::chrono::duration<double, std::milli>;
    auto t_visit_start = Clock::now();
    auto t0 = Clock::now();

    auto& grouping = *ensemble.with_name(m_grouping_name).at(0);
    
    // Container to hold clusters after the initial filter
    std::vector<Cluster*> filtered_clusters;

    Cluster* main_cluster = nullptr;

    for (auto* cluster : grouping.children()) {
        // check scope 
        auto& default_scope = cluster->get_default_scope();
        auto& raw_scope = cluster->get_raw_scope();
        // if scope is not raw, apply filter ...
        if (default_scope.hash()!=raw_scope.hash() && (!cluster->get_scope_filter(default_scope)) ) continue;

        if (cluster->get_flag(Flags::beam_flash)){
            filtered_clusters.push_back(cluster);
            if (cluster->get_flag(Flags::main_cluster)) {
                main_cluster = cluster;
            }
        }
    }

    SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph: {} clusters with beam_flash flag. {}", filtered_clusters.size(), main_cluster->ident());
    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: filter clusters took {} ms", MS(Clock::now() - t0).count());

    if (main_cluster != nullptr){
        // // test steiner ...
        // {
        //     Steiner::Grapher sg(*main_cluster, m_grapher_config, log);
        //     auto& graph = sg.get_graph("basic_pid"); // ensure basic_pid graph exists
        //     std::cout << "CreateSteinerGraph: " << boost::num_vertices(graph) << " vertices " << boost::num_edges(graph) << " edges." << std::endl;
        //     std::vector<size_t> path_point_indices ;
        //     sg.create_steiner_tree(main_cluster, path_point_indices, "basic_pid", "steiner_graph", false, "steiner_pc");
        //     const auto& steiner_point_cloud = sg.get_point_cloud("steiner_pc");
        //     const auto& steiner_graph = sg.get_graph("steiner_graph");
        //     auto& flag_terminals = sg.get_flag_steiner_terminal();
        //     size_t num_true_terminals = std::count(flag_terminals.begin(), flag_terminals.end(), true);

        //     std::cout << "CreateSteinerGraph: " << "steiner_graph with " 
        //               << boost::num_vertices(steiner_graph) << " vertices and "
        //               << boost::num_edges(steiner_graph) << " edges." << " " << flag_terminals.size() << " " << num_true_terminals << std::endl;
        // }


        if (m_grapher_config.retile ) {
            // Call the mutate function with the appropriate configuration, create new cluster
            t0 = Clock::now();
            auto new_node = m_grapher_config.retile->mutate(*main_cluster->node());
            auto new_cluster_1 = new_node->value.facade<Cluster>();
            auto& new_cluster = grouping.make_child();
            new_cluster.take_children(*new_cluster_1);  // Move all blobs from improved cluster
            new_cluster.from(*main_cluster);
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: retile->mutate took {} ms", MS(Clock::now() - t0).count());

            // create the new graph
            t0 = Clock::now();
            new_cluster.find_graph("ctpc_ref_pid", *main_cluster, m_dv, m_pcts);
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: find_graph(ctpc_ref_pid) took {} ms", MS(Clock::now() - t0).count());


            t0 = Clock::now();
            Steiner::Grapher sg(new_cluster, m_grapher_config, log);
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: Grapher construction took {} ms", MS(Clock::now() - t0).count());

            auto& graph = sg.get_graph("ctpc_ref_pid");
            SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph: ctpc_ref_pid with {} vertices and {} edges.", boost::num_vertices(graph), boost::num_edges(graph));

            t0 = Clock::now();
            sg.establish_same_blob_steiner_edges("ctpc_ref_pid", false);
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: establish_same_blob_steiner_edges took {} ms", MS(Clock::now() - t0).count());
            SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph: ctpc_ref_pid with {} vertices and {} edges.", boost::num_vertices(graph), boost::num_edges(graph));

            t0 = Clock::now();
            auto pair_points = new_cluster.get_two_boundary_wcps();
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: get_two_boundary_wcps took {} ms", MS(Clock::now() - t0).count());
            SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph: {} {} {} | {} {} {}", pair_points.first.x(), pair_points.first.y(), pair_points.first.z(), pair_points.second.x(), pair_points.second.y(), pair_points.second.z());

            t0 = Clock::now();
            auto first_index  =   new_cluster.get_closest_point_index(pair_points.first);
            auto second_index =   new_cluster.get_closest_point_index(pair_points.second);
            std::vector<size_t> path_point_indices = new_cluster.graph_algorithms("ctpc_ref_pid").shortest_path(first_index, second_index);
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: shortest_path took {} ms", MS(Clock::now() - t0).count());
            SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph: {} {} # of points along path: {}", first_index, second_index, path_point_indices.size());
            
            t0 = Clock::now();
            sg.remove_same_blob_steiner_edges("ctpc_ref_pid");
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: remove_same_blob_steiner_edges took {} ms", MS(Clock::now() - t0).count());
            SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph: ctpc_ref_pid with {} vertices and {} edges.", boost::num_vertices(graph), boost::num_edges(graph));

            // path_point_indices belong to new_cluster, on which sg is based
            // main_cluster is a reference to filter points ...
            t0 = Clock::now();
            sg.create_steiner_tree(main_cluster, path_point_indices, "ctpc_ref_pid", "steiner_graph", false, "steiner_pc");
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: create_steiner_tree took {} ms", MS(Clock::now() - t0).count());
            const auto& steiner_point_cloud = sg.get_point_cloud("steiner_pc");
            const auto& steiner_graph = sg.get_graph("steiner_graph");
            auto& flag_terminals = sg.get_flag_steiner_terminal();
            size_t num_true_terminals = std::count(flag_terminals.begin(), flag_terminals.end(), true);

            SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph: steiner_graph with {} vertices and {} edges. {} {} {}", boost::num_vertices(steiner_graph), boost::num_edges(steiner_graph), steiner_point_cloud.size(), flag_terminals.size(), num_true_terminals);

            // pass the new_cluster's steiner_graph and stener_pc to the main cluster
            t0 = Clock::now();
            Steiner::Grapher main_sg(*main_cluster, m_grapher_config, log);
            main_sg.transfer_pc(sg, "steiner_pc", "steiner_pc");
            main_sg.transfer_graph(sg, "steiner_graph", "steiner_graph");
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: transfer_pc/graph took {} ms", MS(Clock::now() - t0).count());

            // test ... 
            auto pair_idx = main_cluster->get_two_boundary_steiner_graph_idx("steiner_graph", "steiner_pc", false);
            // std::cout << "Xin3: " << pair_idx.first << " " << pair_idx.second << " " << pair_points.first << std::endl;
            auto kd_results = main_cluster->kd_steiner_knn(1, pair_points.first);
            auto kd_points = main_cluster->kd_steiner_points(kd_results);
            // std::cout << "Xin4: " << kd_points.size() << " " << (*kd_points.begin()).first << " " << (*kd_points.begin()).second.first << " " << (*kd_points.begin()).second.second << std::endl;

            // delete new cluster from grouping after usage ...
            t0 = Clock::now();
            auto* new_cluster_ptr = &new_cluster;
            grouping.destroy_child(new_cluster_ptr, true);
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: destroy_child took {} ms", MS(Clock::now() - t0).count());
        }
    }

    // Create Steiner graphs for associated (non-main) beam_flash clusters.
    // In the prototype, Improve_PR3DCluster_2 (retiling) is called inside create_steiner_graph
    // for every cluster, so we replicate the full main-cluster pipeline here.
    if (m_grapher_config.retile) {
        for (auto* cluster : filtered_clusters) {
            if (cluster == main_cluster) continue;

            t0 = Clock::now();
            auto new_node = m_grapher_config.retile->mutate(*cluster->node());
            auto new_cluster_1 = new_node->value.facade<Cluster>();
            auto& new_cluster = grouping.make_child();
            new_cluster.take_children(*new_cluster_1);
            new_cluster.from(*cluster);
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [assoc {}] retile->mutate took {} ms", cluster->get_cluster_id(), MS(Clock::now() - t0).count());

            t0 = Clock::now();
            new_cluster.find_graph("ctpc_ref_pid", *cluster, m_dv, m_pcts);
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [assoc {}] find_graph(ctpc_ref_pid) took {} ms", cluster->get_cluster_id(), MS(Clock::now() - t0).count());

            t0 = Clock::now();
            Steiner::Grapher sg(new_cluster, m_grapher_config, log);
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [assoc {}] Grapher construction took {} ms", cluster->get_cluster_id(), MS(Clock::now() - t0).count());

            auto& graph = sg.get_graph("ctpc_ref_pid");
            SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph [assoc {}]: ctpc_ref_pid with {} vertices and {} edges.", cluster->get_cluster_id(), boost::num_vertices(graph), boost::num_edges(graph));

            t0 = Clock::now();
            sg.establish_same_blob_steiner_edges("ctpc_ref_pid", false);
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [assoc {}] establish_same_blob_steiner_edges took {} ms", cluster->get_cluster_id(), MS(Clock::now() - t0).count());

            t0 = Clock::now();
            auto pair_points = new_cluster.get_two_boundary_wcps();
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [assoc {}] get_two_boundary_wcps took {} ms", cluster->get_cluster_id(), MS(Clock::now() - t0).count());

            t0 = Clock::now();
            auto first_index  = new_cluster.get_closest_point_index(pair_points.first);
            auto second_index = new_cluster.get_closest_point_index(pair_points.second);
            std::vector<size_t> path_point_indices = new_cluster.graph_algorithms("ctpc_ref_pid").shortest_path(first_index, second_index);
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [assoc {}] shortest_path took {} ms", cluster->get_cluster_id(), MS(Clock::now() - t0).count());

            t0 = Clock::now();
            sg.remove_same_blob_steiner_edges("ctpc_ref_pid");
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [assoc {}] remove_same_blob_steiner_edges took {} ms", cluster->get_cluster_id(), MS(Clock::now() - t0).count());

            t0 = Clock::now();
            sg.create_steiner_tree(cluster, path_point_indices, "ctpc_ref_pid", "steiner_graph", false, "steiner_pc");
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [assoc {}] create_steiner_tree took {} ms", cluster->get_cluster_id(), MS(Clock::now() - t0).count());

            const auto& steiner_point_cloud = sg.get_point_cloud("steiner_pc");
            const auto& steiner_graph = sg.get_graph("steiner_graph");
            auto& flag_terminals = sg.get_flag_steiner_terminal();
            size_t num_true_terminals = std::count(flag_terminals.begin(), flag_terminals.end(), true);
            SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph [assoc {}]: steiner_graph with {} vertices and {} edges. {} {} {}",
                                cluster->get_cluster_id(), boost::num_vertices(steiner_graph), boost::num_edges(steiner_graph),
                                steiner_point_cloud.size(), flag_terminals.size(), num_true_terminals);

            t0 = Clock::now();
            Steiner::Grapher cluster_sg(*cluster, m_grapher_config, log);
            cluster_sg.transfer_pc(sg, "steiner_pc", "steiner_pc");
            cluster_sg.transfer_graph(sg, "steiner_graph", "steiner_graph");
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [assoc {}] transfer_pc/graph took {} ms", cluster->get_cluster_id(), MS(Clock::now() - t0).count());

            t0 = Clock::now();
            auto* new_cluster_ptr = &new_cluster;
            grouping.destroy_child(new_cluster_ptr, true);
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [assoc {}] destroy_child took {} ms", cluster->get_cluster_id(), MS(Clock::now() - t0).count());
        }
    }

    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: visit() TOTAL took {} ms", MS(Clock::now() - t_visit_start).count());

    // for (auto* cluster : filtered_clusters) {
        // bool already = cluster->has_graph(m_graph_name);
        // if (already || m_replace) {
     
        //     const auto& new_to_old = sg.get_new_to_old_mapping();
        //     std::cout << "Xin2: " << cell_points_map.size() << " Graph vertices: " << boost::num_vertices(steiner_graph) << ", edges: " << boost::num_edges(steiner_graph) << " " << steiner_point_cloud.size_major()  <<  " " << num_true_terminals << std::endl;


        //     // auto edge_weight_map = get(boost::edge_weight, steiner_graph);
        //     // for (auto edge_it = boost::edges(steiner_graph); edge_it.first != edge_it.second; ++edge_it.first) {
        //     //     auto edge = *edge_it.first;
        //     //     auto src = boost::source(edge, steiner_graph);
        //     //     auto tgt = boost::target(edge, steiner_graph);
        //     //     // Use the new-to-old mapping to get original vertex IDs
        //     //     auto src_id = new_to_old.at(src);
        //     //     auto tgt_id = new_to_old.at(tgt);

        //     //     // Get the edge weight using the proper accessor
        //     //     auto weight = edge_weight_map[edge];

        //     //     std::cout << "Edge from vertex " << cluster->point3d(src_id) << " to " << cluster->point3d(tgt_id) << " with weight " << weight << " " << flag_terminals[src] << " " << flag_terminals[tgt] << std::endl;
        //     // }

        //     // for (const auto& [cell, points] : cell_points_map) {
        //     //     // std::cout << "Xin2 Cell: " << cell->slice_index_min() << " " << cell->u_wire_index_min() << " " << cell->v_wire_index_min() << " " << cell->w_wire_index_min() << " has " << points.size() << " points." << std::endl;
        //     //     // for (const auto& point : points) {
        //     //         // auto info = cluster->calc_charge_wcp(point);
        //     //         // std::cout << "Xin2 Point: " << point << " " << info.first << " " << info.second << std::endl;
        //     //     // }
        //     //     std::vector<const Blob*> single_blob = {cell};
        //     //     auto blob_peaks = sg.find_peak_point_indices(single_blob, "basic_pid", true);
        //     //     std::cout << "Xin2: " << cell->slice_index_min() << " " << cell->u_wire_index_min() << " " << cell->v_wire_index_min() << " " << cell->w_wire_index_min()  << "  " << points.size() << "  " << blob_peaks.size() <<std::endl;

        //     // }
        //     // auto steiner_terminals = sg.find_steiner_terminals("basic_pid");
        //     // std::cout << "Xin3: " << steiner_terminals.size() << std::endl;
        //     // auto extrem_points = cluster->get_extreme_wcps();
        //     // std::cout << "Xin4: " << extrem_points.size() << std::endl;
        //     // for (const auto& pts : extrem_points) {
        //     //     for (const auto& pt : pts) {
        //     //         std::cout << "Extreme point: ("
        //     //                   << pt.x() << ", "
        //     //                   << pt.y() << ", "
        //     //                   << pt.z() << ")" << std::endl;
        //     //     }
        //     // }

        //     // auto extreme_boundary_points = cluster->get_two_boundary_wcps(1);
        //     // std::cout << "Xin3: " << extreme_boundary_points.first.x() << " " 
        //     //           << extreme_boundary_points.first.y() << " " 
        //     //           << extreme_boundary_points.first.z() << " | "
        //     //           << extreme_boundary_points.second.x() << " " 
        //     //           << extreme_boundary_points.second.y() << " " 
        //     //           << extreme_boundary_points.second.z() << std::endl;

        //     // for (const auto& [cell, points] : cell_points_map) {
        //     //     auto total_charge = cell->estimate_total_charge();
        //     //     auto min_charge = cell->estimate_minimum_charge();

        //     //     std::cout << "Xin2 Cell: " << cell->slice_index_min() << " " << cell->u_wire_index_min() << " " << cell->v_wire_index_min() << " " << cell->w_wire_index_min() 
        //     //               << " has " << points.size() << " points, total charge: " << total_charge 
        //     //               << ", min charge: " << min_charge << " " << cell->get_wire_charge(0, cell->u_wire_index_min()) << " " << cell->get_wire_charge_error(0, cell->u_wire_index_min()) 
        //     //               << std::endl;
        //     //     // cell->check_dead_wire_consistency();
        //     // }

        //     //  std::cout << grouping.get_dead_winds1(0,0,0).size() << " " << grouping.get_dead_winds(0,0,0).size() << " " << grouping.get_dead_winds1(0,0,1).size() << " " << grouping.get_dead_winds(0,0,1).size() << " " << grouping.get_dead_winds1(0,0,2).size() << " " << grouping.get_dead_winds(0,0,2).size() << std::endl;

    // }
}
