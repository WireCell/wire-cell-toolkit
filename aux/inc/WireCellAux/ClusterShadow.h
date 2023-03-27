/** Cluster Shadow

    See aux/docs/cluster-shadow.org for details.
 */
#ifndef WIRECELL_AUX_CLUSTERSHADOW
#define WIRECELL_AUX_CLUSTERSHADOW

#include "WireCellIface/ICluster.h"
#include "WireCellUtil/Rectangles.h"
#include "WireCellAux/BlobShadow.h"

namespace WireCell::Aux::ClusterShadow {

    struct Graph {
        // The shadow type code w=wire, c=channel copied from bsgraph.
        char stype;
    };

    // The "coverage" is a Rectangles to describe each tomographic
    // view of the cluster.  An interval on the Y-axis spans either
    // channel or wire indices (integral number).  An interval on the
    // X-axis spans a time slice (floating point number).  In
    // principle the time slice spans need not be uniform.  The value
    // of the Rectangles map are blob vertex descriptors into the
    // original ICluster graph.  The coverages of a pair of clusters
    // in shadow may be compared with set-theoretic operations
    // (difference, intersection, etc).
    using coverage_t = Rectangles<cluster_vertex_t, double, int>;

    // A node summarizes a geometric cluster.
    struct Node {
        // The region of each view covered by the cluster, indexed by
        // wpid.index().
        std::vector<coverage_t> coverage;
 
        // Number of blobs in the cluster.
        int nblobs;

        // The sum of (charge) value and its uncertainty over the blobs
        using value_t = WireCell::Measurement::float32;
        value_t value;

    };
    
    // An edge summarizes the collection of blob shadows in a view and
    // between two clusters.
    struct Edge {
        // The wire plane in which the shadow exists.  The
        // wpid.face() value is only relevant for a wire shadow.
        WirePlaneId wpid{0};
    };

    using graph_t = boost::adjacency_list<boost::setS, boost::vecS,
                                          boost::undirectedS,
                                          Node, Edge, Graph>;
    using vdesc_t = boost::graph_traits<graph_t>::vertex_descriptor;
    using edesc_t = boost::graph_traits<graph_t>::edge_descriptor;


    // This records a map from a blob vertex descriptor of a ICluster
    // graph to the cluster shadow graph vertex (ie, a cluster).  
    using blob_cluster_map_t = std::unordered_map<cluster_vertex_t, vdesc_t>;

    // Return a cluster shadow graph generated from an ICluster graph
    // and a blob shadow graph.  The "clusters" will be filled to
    // associate each b-type cluster graph vertex with a vertex in the
    // returned cluster shadow graph.
    graph_t shadow(const cluster_graph_t& cgraph,
                   const BlobShadow::graph_t& bsgraph,
                   blob_cluster_map_t & clusters);

    // As above but no cluster map.
    graph_t shadow(const cluster_graph_t& cgraph,
                   const BlobShadow::graph_t& bsgraph);

}


#endif
