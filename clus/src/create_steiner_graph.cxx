#include "make_graphs.h"

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Graphs.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Grouping.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PointTree.h"


class CreateSteinerGraph;
WIRECELL_FACTORY(CreateSteinerGraph, CreateSteinerGraph,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using WireCell::Clus::Facade::Ensemble;
using WireCell::Clus::Facade::Cluster;
using WireCell::Clus::Facade::Grouping;
using graph_type = WireCell::Clus::Graphs::Weighted::Graph;
using namespace WireCell::PointCloud::Tree; // Points, Scope, etc

class CreateSteinerGraph : public IConfigurable, public Clus::IEnsembleVisitor {
    std::string m_grouping_name{"live"};
    std::string m_graph_name{"steiner"};
    bool m_replace{true};

public:
    virtual void configure(const WireCell::Configuration& cfg) {
        m_grouping_name = get(cfg, "grouping", m_grouping_name);
        m_graph_name = get(cfg, "graph", m_graph_name);
        m_replace = get(cfg, "replace", m_replace);
    }
    virtual Configuration default_configuration() const {
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

    virtual void visit(Ensemble& ensemble) const;

private:

};


// Get a COPY of a graph held by a cluster
static graph_type get_graph(const Cluster& cluster, const std::string& flavor)
{
    const auto& graph = cluster.get_graph(flavor);
    graph_type copy;
    boost::copy_graph(graph, copy);
    return copy;
}

static void establish_same_blob_steiner_edges(graph_type& new_graph, const Cluster& cluster /*,...*/)
{
}

static void improve_cluster(Cluster& new_cluster, const Cluster& orig_cluster /*,...*/)
{
}

static void improve_cluster_1(Cluster& new_cluster, const Cluster& orig_cluster /*,...*/)
{
}

// Reimplement Improve_PR3DCluster_2
static void improve_cluster_2(Cluster& new_cluster, const Cluster& orig_cluster
                              /* , int nrebin, int frame_length, double unit_dis,...*/)
{
    auto graph = get_graph(orig_cluster, "ctpc");
    establish_same_blob_steiner_edges(graph, new_cluster /*,...*/);

    // QUESTION: in WCP we have:
    // pr3dcluster->dijkstra_shortest_paths(wcps1.first);
    // pr3dcluster->cal_shortest_path(wcps1.second);
    // pr3dcluster->remove_same_mcell_steiner_edges();
    // pr3dcluster->Del_graph();
    //
    // Why?  Here, the local "graph" holds the steiner edges and the
    // dijk. info is made on demand via GA so I think we do not need all this.
    // But, maybe we need this on the new cluster:?
    new_cluster.give_graph("ctpc", std::move(graph));

    improve_cluster_1(new_cluster, orig_cluster); // ....

    improve_cluster(new_cluster, orig_cluster);
}



static void calc_sampling_points(Cluster& new_cluster /*, ...*/)
{
}

// Fixme: I don't want to return a copy, but for now, that's the mock up.
static PointCloud::Dataset get_point_cloud(const Cluster& cluster)
{
    return PointCloud::Dataset {};
}

static graph_type create_steiner_tree(Cluster& new_cluster, PointCloud::Dataset& point_cloud_steiner /*,...*/)
{
    return graph_type{};
}

static graph_type create_steiner_graph(const Cluster& orig_cluster /*,...*/)
{
    Points::node_t new_cluster_node;
    Cluster& new_cluster = *new_cluster_node.value.facade<Cluster>();

    // Replaces Improve_PR3DCluster_2
    improve_cluster_2(new_cluster, orig_cluster/*, ...*/);

    // Replaces WCP function of same name
    calc_sampling_points(new_cluster /*, ...*/);

    // This replaces Create_point_cloud()
    auto new_pc = get_point_cloud(new_cluster); // which PC is this?

    // This replaces Create_graph(ctpc, pc)
    auto new_graph = WireCell::Clus::Graphs::make_graph_basic(new_cluster /*,new_pc*/); // fixme; this is almost certainly not the graph to make.  which one (relaxed? ctpc?), how to provide args?

    // Replaces as establish_same_mcell_steiner_edges()
    establish_same_blob_steiner_edges(new_graph, new_cluster /*,...*/);

    /*
     * Note: no need (?) to "prime" dijkstra as they are calculated when needed.
     * Also, no need to "remove_same_mcell_steiner_edges()" as the new_* objects
     * are all temporaries.
     */

    // Replaces WCP's Create_steiner_tree
    // steiner tree with some basic cuts ...
    PointCloud::Dataset point_cloud_steiner;
    // fixme: needs some more input
    auto graph_steiner = create_steiner_tree(new_cluster, point_cloud_steiner /*,...*/);
    auto orig_graph = get_graph(orig_cluster, "...");
    establish_same_blob_steiner_edges(orig_graph, orig_cluster  /*,...*/);

    PointCloud::Dataset point_cloud_steiner_terminal;
    // transfer u, v, w points from point_cloud_"steiner.

    return graph_steiner;
}

void CreateSteinerGraph::visit(Ensemble& ensemble) const
{
    auto& grouping = *ensemble.with_name(m_grouping_name).at(0);
    for (auto* cluster : grouping.children()) {
        auto& graph = cluster->get_graph(m_graph_name); // creates empty if not existing
        if (boost::num_vertices(graph) == 0 || m_replace) {
            graph = create_steiner_graph(*cluster);
        }
    }
}
