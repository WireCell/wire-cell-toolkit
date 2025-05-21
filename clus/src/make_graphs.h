#ifndef WIRECELLCLUS_PRIVATE_MAKE_GRAPHS
#define WIRECELLCLUS_PRIVATE_MAKE_GRAPHS

#include "WireCellClus/Graph.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellClus/IPCTransform.h"


namespace WireCell::Clus {
    namespace Facade {
        class Cluster;
    }

    // factory functions wrapping up construction and various connect_graph*
    // functions.

    // just closely connected.
    Graph::Ident::graph_ptr make_graph_closely(
        const Facade::Cluster& cluster);

    // closely + basic connection
    Graph::Ident::graph_ptr make_graph_basic(
        const Facade::Cluster& cluster);

    // closely + ctpc connection
    Graph::Ident::graph_ptr make_graph_ctpc(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv, 
        IPCTransformSet::pointer pcts);

    // closely + overclustering protection connection
    Graph::Ident::graph_ptr make_graph_overclustering_protection(
        const Facade::Cluster& clus,
        IDetectorVolumes::pointer dv, 
        IPCTransformSet::pointer pcts);


}

#endif
