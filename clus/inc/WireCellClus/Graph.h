#ifndef WIRECELLCLUS_GRAPH
#define WIRECELLCLUS_GRAPH

#include "WireCellUtil/Graph.h"

// Graph types used in clus.
namespace WireCell::Clus::Graph {

    /// An "ident" graph type has a vertex property of integer type.  This
    /// integer is intended to represent the "identity" of some externally
    /// provided data, entity or instance.  Typically, the integer is used as an
    /// index in an array of values.  The integer vertex property of an "ident"
    /// graph is (of course) distinct from the vertex descriptor (also an
    /// integer).  The vertex "ident" is user-provided and in general not the
    /// same value of the vertex "descriptor".  The descriptor value is
    /// determined internally by boost::graph.  A single precision floating
    /// point "weight" edge property is provided.
    namespace Ident {


        struct VertexProp {
            int ident;
        };
        using EdgeProp = boost::property<boost::edge_weight_t, float>; 

        using graph_type = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                                                  VertexProp, EdgeProp>;
        using vertex_descriptor = boost::graph_traits<graph_type>::vertex_descriptor;
        using edge_descriptor = boost::graph_traits<graph_type>::edge_descriptor;
        
        /// Return the ident values in vertex order.
        std::vector<int> idents(const graph_type& g);

        /// Return the edge weights in edge order.
        std::vector<float> weights(const graph_type& g);

        /// Return the edges as matched pairs of vertex indices.  The first
        /// vector in the pair holds the "tail" vertices, the second holds the
        /// "head" vertices.
        std::pair<std::vector<int>, std::vector<int>>
        edges(const graph_type& g);

    }

    // A weighted-edge graph
    namespace Weighted {
        using EdgeProp = boost::property<boost::edge_weight_t, double>;
        using graph_type = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                                                 boost::no_property, EdgeProp>;
        using vertex_descriptor = boost::graph_traits<graph_type>::vertex_descriptor;
        using edge_descriptor = boost::graph_traits<graph_type>::edge_descriptor;
    }        


}

#endif

