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

    SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph: {} clusters with beam_flash flag. main={}", filtered_clusters.size(), main_cluster ? main_cluster->ident() : -1);
    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: filter clusters took {} ms", MS(Clock::now() - t0).count());

    // Helper that runs the full retile→graph→Steiner pipeline for one cluster.
    // `src` is the cluster to retile; `ref` is the reference cluster used for
    // find_graph, create_steiner_tree, and the final graph transfer.
    // For the main cluster src==ref; for associated clusters src==ref as well
    // (retiled copy is a fresh child, but the graph is seeded from the same
    // cluster and transferred back to it).
    // `is_main` triggers the extra kd_steiner_knn probe done only for the main
    // cluster.  Returns true if a steiner_graph was successfully produced.
    auto process_cluster_steiner = [&](Cluster* src, bool is_main) -> bool {
        const std::string tag = is_main
            ? std::string("main ") + std::to_string(src->ident())
            : std::string("assoc ") + std::to_string(src->get_cluster_id());

        t0 = Clock::now();
        auto new_node = m_grapher_config.retile->mutate(*src->node());
        auto new_cluster_1 = new_node->value.facade<Cluster>();
        auto& new_cluster = grouping.make_child();
        new_cluster.take_children(*new_cluster_1);
        new_cluster.from(*src);
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [{}] retile->mutate took {} ms", tag, MS(Clock::now() - t0).count());

        t0 = Clock::now();
        new_cluster.find_graph("ctpc_ref_pid", *src, m_dv, m_pcts);
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [{}] find_graph(ctpc_ref_pid) took {} ms", tag, MS(Clock::now() - t0).count());

        t0 = Clock::now();
        Steiner::Grapher sg(new_cluster, m_grapher_config, log);
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [{}] Grapher construction took {} ms", tag, MS(Clock::now() - t0).count());

        auto& graph = sg.get_graph("ctpc_ref_pid");
        SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph [{}]: ctpc_ref_pid with {} vertices and {} edges.", tag, boost::num_vertices(graph), boost::num_edges(graph));

        t0 = Clock::now();
        sg.establish_same_blob_steiner_edges("ctpc_ref_pid", false);
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [{}] establish_same_blob_steiner_edges took {} ms", tag, MS(Clock::now() - t0).count());
        SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph [{}]: ctpc_ref_pid with {} vertices and {} edges.", tag, boost::num_vertices(graph), boost::num_edges(graph));

        t0 = Clock::now();
        auto pair_points = new_cluster.get_two_boundary_wcps();
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [{}] get_two_boundary_wcps took {} ms", tag, MS(Clock::now() - t0).count());
        SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph [{}]: {} {} {} | {} {} {}", tag,
            pair_points.first.x(), pair_points.first.y(), pair_points.first.z(),
            pair_points.second.x(), pair_points.second.y(), pair_points.second.z());

        t0 = Clock::now();
        auto first_index  = new_cluster.get_closest_point_index(pair_points.first);
        auto second_index = new_cluster.get_closest_point_index(pair_points.second);
        std::vector<size_t> path_point_indices = new_cluster.graph_algorithms("ctpc_ref_pid").shortest_path(first_index, second_index);
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [{}] shortest_path took {} ms", tag, MS(Clock::now() - t0).count());
        SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph [{}]: {} {} # of points along path: {}", tag, first_index, second_index, path_point_indices.size());

        t0 = Clock::now();
        sg.remove_same_blob_steiner_edges("ctpc_ref_pid");
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [{}] remove_same_blob_steiner_edges took {} ms", tag, MS(Clock::now() - t0).count());
        SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph [{}]: ctpc_ref_pid with {} vertices and {} edges.", tag, boost::num_vertices(graph), boost::num_edges(graph));

        t0 = Clock::now();
        sg.create_steiner_tree(src, path_point_indices, "ctpc_ref_pid", "steiner_graph", false, "steiner_pc");
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [{}] create_steiner_tree took {} ms", tag, MS(Clock::now() - t0).count());

        auto* new_cluster_ptr = &new_cluster;
        if (!new_cluster.has_graph("steiner_graph")) {
            SPDLOG_LOGGER_WARN(log, "CreateSteinerGraph: create_steiner_tree produced no steiner_graph for {}, skipping transfer", tag);
            // Fall through: for the main cluster, still process associated clusters.
            grouping.destroy_child(new_cluster_ptr, true);
            return false;
        }

        const auto& steiner_point_cloud = sg.get_point_cloud("steiner_pc");
        const auto& steiner_graph       = sg.get_graph("steiner_graph");
        auto& flag_terminals = sg.get_flag_steiner_terminal();
        size_t num_true_terminals = std::count(flag_terminals.begin(), flag_terminals.end(), true);
        SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph [{}]: steiner_graph with {} vertices and {} edges. {} {} {}", tag,
            boost::num_vertices(steiner_graph), boost::num_edges(steiner_graph),
            steiner_point_cloud.size(), flag_terminals.size(), num_true_terminals);

        t0 = Clock::now();
        Steiner::Grapher ref_sg(*src, m_grapher_config, log);
        ref_sg.transfer_pc(sg, "steiner_pc", "steiner_pc");
        ref_sg.transfer_graph(sg, "steiner_graph", "steiner_graph");
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [{}] transfer_pc/graph took {} ms", tag, MS(Clock::now() - t0).count());

        if (is_main) {
            // Extra probe done only for the main cluster.
            (void)src->get_two_boundary_steiner_graph_idx("steiner_graph", "steiner_pc", false);
            auto kd_results = src->kd_steiner_knn(1, pair_points.first);
            (void)src->kd_steiner_points(kd_results);
        }

        t0 = Clock::now();
        grouping.destroy_child(new_cluster_ptr, true);
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: [{}] destroy_child took {} ms", tag, MS(Clock::now() - t0).count());
        return true;
    };

    if (m_grapher_config.retile) {
        if (main_cluster != nullptr) {
            process_cluster_steiner(main_cluster, /*is_main=*/true);
        }

        // Associated (non-main) beam_flash clusters.
        for (auto* cluster : filtered_clusters) {
            if (cluster == main_cluster) continue;
            process_cluster_steiner(cluster, /*is_main=*/false);
        }
    }

    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "CreateSteinerGraph timing: visit() TOTAL took {} ms", MS(Clock::now() - t_visit_start).count());

}
