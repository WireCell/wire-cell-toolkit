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

Weighted::Graph WireCell::Clus::Graphs::make_graph_relaxed_fast(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv,
    IPCTransformSet::pointer pcts)
{
    auto graph = make_graph_closely(cluster);
    static const RelaxedFastCfg fast_cfg;   // defaults, connect_graphs.h
    connect_graph_relaxed(cluster, dv, pcts, graph, &fast_cfg);
    return graph;
}

Weighted::Graph WireCell::Clus::Graphs::make_graph_relaxed_pid(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv,
    IPCTransformSet::pointer pcts)
{
    auto graph = make_graph_closely_pid(cluster);
    connect_graph_relaxed_pid(cluster, dv, pcts, graph);
    return graph;
}

Weighted::Graph WireCell::Clus::Graphs::make_graph_relaxed_strict(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv,
    IPCTransformSet::pointer pcts)
{
    auto graph = make_graph_closely(cluster);
    connect_graph_relaxed_strict(cluster, dv, pcts, graph);
    return graph;
}

Weighted::Graph WireCell::Clus::Graphs::make_graph_relaxed_strict_img(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv,
    IPCTransformSet::pointer pcts)
{
    auto graph = make_graph_closely(cluster);
    connect_graph_relaxed_strict(cluster, dv, pcts, graph, /*image_check=*/true);
    return graph;
}

Weighted::Graph WireCell::Clus::Graphs::make_graph_relaxed_strict_img_2d(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv,
    IPCTransformSet::pointer pcts)
{
    auto graph = make_graph_closely(cluster);
    connect_graph_relaxed_strict(cluster, dv, pcts, graph, /*image_check=*/true, /*two_d_check=*/true);
    return graph;
}

Weighted::Graph WireCell::Clus::Graphs::make_graph_relaxed_strict_img_2d_wfloor(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv,
    IPCTransformSet::pointer pcts)
{
    auto graph = make_graph_closely(cluster);
    connect_graph_relaxed_strict(cluster, dv, pcts, graph, /*image_check=*/true, /*two_d_check=*/true,
                                  /*floor_w_override=*/true);
    return graph;
}

Weighted::Graph WireCell::Clus::Graphs::make_graph_relaxed_strict_img_2d_rescue(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv,
    IPCTransformSet::pointer pcts)
{
    auto graph = make_graph_closely(cluster);
    connect_graph_relaxed_strict(cluster, dv, pcts, graph, /*image_check=*/true, /*two_d_check=*/true,
                                  /*floor_w_override=*/true, /*two_d_rescue=*/true);
    return graph;
}

Weighted::Graph WireCell::Clus::Graphs::make_graph_relaxed_strict_img_2d_rescue_long(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv,
    IPCTransformSet::pointer pcts)
{
    auto graph = make_graph_closely(cluster);
    connect_graph_relaxed_strict(cluster, dv, pcts, graph, /*image_check=*/true, /*two_d_check=*/true,
                                  /*floor_w_override=*/true, /*two_d_rescue=*/true,
                                  /*long_check=*/true, /*long_min_planes=*/1);
    return graph;
}

Weighted::Graph WireCell::Clus::Graphs::make_graph_relaxed_strict_img_2d_rescue_long2(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv,
    IPCTransformSet::pointer pcts)
{
    auto graph = make_graph_closely(cluster);
    connect_graph_relaxed_strict(cluster, dv, pcts, graph, /*image_check=*/true, /*two_d_check=*/true,
                                  /*floor_w_override=*/true, /*two_d_rescue=*/true,
                                  /*long_check=*/true, /*long_min_planes=*/2);
    return graph;
}

Weighted::Graph WireCell::Clus::Graphs::make_graph_relaxed_strict_img_2d_rescue_long_wtrack(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv,
    IPCTransformSet::pointer pcts)
{
    auto graph = make_graph_closely(cluster);
    connect_graph_relaxed_strict(cluster, dv, pcts, graph, /*image_check=*/true, /*two_d_check=*/true,
                                  /*floor_w_override=*/true, /*two_d_rescue=*/true,
                                  /*long_check=*/true, /*long_min_planes=*/1,
                                  /*w_track_excuse=*/true);
    return graph;
}
