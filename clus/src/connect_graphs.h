#ifndef WIRECELLCLUS_PRIVATE_CONNECT_GRAPHS
#define WIRECELLCLUS_PRIVATE_CONNECT_GRAPHS

#include "WireCellClus/Graphs.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellClus/IPCTransform.h"
#include "WireCellClus/Facade_Cluster.h"

namespace WireCell::Clus::Graphs {


    // See make_graphs.h for bundling of construction and connecting.

    void connect_graph(const Facade::Cluster& cluster, 
                       Weighted::Graph& graph);

    void connect_graph_ctpc(const Facade::Cluster& cluster,
                            IDetectorVolumes::pointer dv,
                            Clus::IPCTransformSet::pointer pcts,
                            Weighted::Graph& graph);

    void connect_graph_closely(const Facade::Cluster& cluster,
                               Weighted::Graph& graph);


    // ne' overclustering protection
    void connect_graph_relaxed(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv, 
        IPCTransformSet::pointer pcts,
        Weighted::Graph& graph);


}
    
#endif
