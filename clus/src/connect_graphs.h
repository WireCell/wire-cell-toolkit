#ifndef WIRECELLCLUS_PRIVATE_CONNECT_GRAPHS
#define WIRECELLCLUS_PRIVATE_CONNECT_GRAPHS

#include "WireCellClus/Graph.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellClus/IPCTransform.h"


namespace WireCell::Clus {
    namespace Facade {
        class Cluster;
    }


    // See make_graphs.h for bundling of construction and connecting.
    void connect_graph(const Facade::Cluster& cluster, 
                       Graph::Ident::graph_type& graph);

    void connect_graph_ctpc(const Facade::Cluster& cluster,
                            IDetectorVolumes::pointer dv,
                            Clus::IPCTransformSet::pointer pcts,
                            Graph::Ident::graph_type& graph);

    void connect_graph_closely(const Facade::Cluster& cluster,
                               Graph::Ident::graph_type& graph);

    void connect_graph_overclustering_protection(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv, 
        IPCTransformSet::pointer pcts,
        Graph::Ident::graph_type& graph);


}

#endif
