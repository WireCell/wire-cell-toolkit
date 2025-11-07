#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Grouping.h"

#include "connect_graphs.h"
#include "make_graphs.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::Clus::Graphs;


Weighted::Graph WireCell::Clus::Graphs::make_graph_closely(
    const Cluster& cluster) 
{
    Weighted::Graph graph(cluster.npoints());
    connect_graph_closely(cluster, graph);
    return graph;
}

Weighted::Graph WireCell::Clus::Graphs::make_graph_closely_pid(
    const Cluster& cluster) 
{
    Weighted::Graph graph(cluster.npoints());
    connect_graph_closely_pid(cluster, graph);
    return graph;
}

Weighted::Graph WireCell::Clus::Graphs::make_graph_basic(
    const Cluster& cluster) 
{
    auto graph = make_graph_closely(cluster);
    connect_graph(cluster, graph);
    return graph;
}

Weighted::Graph WireCell::Clus::Graphs::make_graph_basic_pid(
    const Cluster& cluster,
    const Cluster& ref_cluster) 
{
    auto graph = make_graph_closely_pid(cluster);
    
    connect_graph_with_reference(cluster, ref_cluster, graph);
    
    return graph;
}

Weighted::Graph WireCell::Clus::Graphs::make_graph_ctpc(
    const Cluster& cluster,
    IDetectorVolumes::pointer dv, 
    IPCTransformSet::pointer pcts) 
{
    auto graph = make_graph_closely(cluster);
    connect_graph_ctpc(cluster, dv, pcts, graph);
    connect_graph(cluster, graph);
    return graph;
}

Weighted::Graph WireCell::Clus::Graphs::make_graph_ctpc_pid(
        const Cluster& cluster,
        const Cluster& ref_cluster,
        IDetectorVolumes::pointer dv, 
        IPCTransformSet::pointer pcts)
{
    // Start with close connections
    auto graph = make_graph_closely_pid(cluster);
    
    // Add CTPC connections with reference filtering
    connect_graph_ctpc_with_reference(cluster, ref_cluster, dv, pcts, graph);
    connect_graph_with_reference(cluster, ref_cluster, graph);

    return graph;
}


Weighted::Graph WireCell::Clus::Graphs::make_graph_relaxed(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv, 
    IPCTransformSet::pointer pcts) 
{
    auto graph = make_graph_closely(cluster);
    connect_graph_relaxed(cluster, dv, pcts, graph);
    return graph;
}
