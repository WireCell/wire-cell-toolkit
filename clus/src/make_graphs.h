#ifndef WIRECELLCLUS_PRIVATE_MAKE_GRAPHS
#define WIRECELLCLUS_PRIVATE_MAKE_GRAPHS

#include "WireCellClus/Graphs.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellClus/IPCTransform.h"
#include "WireCellClus/Facade_Cluster.h"


namespace WireCell::Clus::Graphs {

    // factory functions wrapping up construction and various connect_graph*
    // functions.

    // just closely connected.
    Weighted::Graph make_graph_closely(
        const Facade::Cluster& cluster);

    // closely + basic connection
    Weighted::Graph make_graph_basic(
        const Facade::Cluster& cluster);

    // closely + ctpc connection
    Weighted::Graph make_graph_ctpc(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv, 
        IPCTransformSet::pointer pcts);

    // closely + relaxed (overclustering protection)
    Weighted::Graph make_graph_relaxed(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv, 
        IPCTransformSet::pointer pcts);

}

#endif
