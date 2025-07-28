#include "WireCellClus/PRGraph.h"

namespace WireCell::Clus::PR {

    GraphProperty& get_graph_bundle(graph_type& sg) {
        // Note, this may look totally sketchy because it is!  There is
        // apparently no other way to get the underlying graph bundle than to go
        // behind the subgraph's back and access its underlying (nonsubgraph)
        // graph.  Boost's subgraph.hpp file has this curious code:
        //
        // public: // Probably shouldn't be public....
        //     Graph m_graph;
        //
        // It is reasonable to expect one day that public is changed to private
        // (or protected).  When that day comes, this function will break.  I
        // can only hope that such a change also brings a way to get the graph
        // bundle property or the underlying graph.  This function will then
        // need a compile-time switch to pick between old and new access
        // methods.
        return sg.m_graph[boost::graph_bundle];
    }
    const GraphProperty& get_graph_bundle(const graph_type& sg) {
        return sg.m_graph[boost::graph_bundle];
    }
}

    
