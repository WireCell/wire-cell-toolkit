#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Grouping.h"

#include "connect_graphs.h"
#include "make_graphs.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;


Graph::Ident::graph_ptr WireCell::Clus::make_graph_closely(
    const Cluster& cluster) 
{
    auto gptr = std::make_unique<Graph::Ident::graph_type>(cluster.npoints());
    connect_graph_closely(cluster, *gptr);
    return gptr;
}

Graph::Ident::graph_ptr WireCell::Clus::make_graph_basic(
    const Cluster& cluster) 
{
    auto gptr = make_graph_closely(cluster);
    connect_graph(cluster, *gptr);
    return gptr;
}

Graph::Ident::graph_ptr WireCell::Clus::make_graph_ctpc(
    const Cluster& cluster,
    IDetectorVolumes::pointer dv, 
    IPCTransformSet::pointer pcts) 
{
    auto gptr = make_graph_closely(cluster);
    connect_graph_ctpc(cluster, dv, pcts, *gptr);
    connect_graph(cluster, *gptr);
    return gptr;
}

Graph::Ident::graph_ptr WireCell::Clus::make_graph_overclustering_protection(
    const Cluster& cluster,
    IDetectorVolumes::pointer dv, 
    IPCTransformSet::pointer pcts) 
{
    auto gptr = make_graph_closely(cluster);
    connect_graph_overclustering_protection(cluster, dv, pcts, *gptr);
    return gptr;
}

