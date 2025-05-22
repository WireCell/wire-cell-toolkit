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


Weighted::GraphPtr WireCell::Clus::Graphs::make_graph_closely(
    const Cluster& cluster) 
{
    auto graph = std::make_unique<Weighted::Graph>(cluster.npoints());
    connect_graph_closely(cluster, *graph);
    return graph;
}

Weighted::GraphPtr WireCell::Clus::Graphs::make_graph_basic(
    const Cluster& cluster) 
{
    auto graph = make_graph_closely(cluster);
    connect_graph(cluster, *graph);
    return graph;
}

Weighted::GraphPtr WireCell::Clus::Graphs::make_graph_ctpc(
    const Cluster& cluster,
    IDetectorVolumes::pointer dv, 
    IPCTransformSet::pointer pcts) 
{
    auto graph = make_graph_closely(cluster);
    connect_graph_ctpc(cluster, dv, pcts, *graph);
    connect_graph(cluster, *graph);
    return graph;
}

Weighted::GraphPtr WireCell::Clus::Graphs::make_graph_overclustering_protection(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv, 
    IPCTransformSet::pointer pcts) 
{
    auto graph = make_graph_closely(cluster);
    connect_graph_overclustering_protection(cluster, dv, pcts, *graph);
    return graph;
}
