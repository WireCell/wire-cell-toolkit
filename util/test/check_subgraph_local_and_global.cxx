#include "WireCellUtil/Subgraph.h"

#include <iostream>
#include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/subgraph.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/properties.hpp>

// Bundled properties and graph types
struct VertexProperties { char name; };
struct EdgeProperties { double weight; };
using BaseGraph = boost::adjacency_list<
    boost::setS, boost::set, boost::undirectedS,
    boost::property<boost::vertex_index_t, int, VertexProperties>,
    boost::property<boost::edge_index_t, int, EdgeProperties>
>;
using SGraph = boost::subgraph<BaseGraph>;
using VertexDescriptor = typename boost::graph_traits<SGraph>::vertex_descriptor;
using EdgeDescriptor = typename boost::graph_traits<SGraph>::edge_descriptor;

int main() {
    // S0: The root graph
    SGraph S0(4); // Call to constructor for initial vertices
    S0[(VertexDescriptor)0].name = 'a';
    S0[(VertexDescriptor)1].name = 'b';
    S0[(VertexDescriptor)2].name = 'c';
    S0[(VertexDescriptor)3].name = 'd';

    // Free function calls for add_edge
    boost::add_edge(0, 1, S0);
    boost::add_edge(1, 2, S0);
    boost::add_edge(2, 3, S0);
    
    S0[boost::edge(0, 1, S0).first].weight = 10.0;
    S0[boost::edge(1, 2, S0).first].weight = 20.0;
    S0[boost::edge(2, 3, S0).first].weight = 30.0;

    std::cout << "S0 (Root) Edges:" << std::endl;
    for (auto e : boost::make_iterator_range(boost::edges(S0))) {
        // When using 'S0' (the root subgraph), its local_to_global and global_to_local map to itself.
        // So, S0.local_to_global(desc) just returns desc.
        std::cout << "  Global (" << S0.local_to_global(boost::source(e, S0))
                  << ", " << S0.local_to_global(boost::target(e, S0))
                  << ") | " << S0[boost::source(e, S0)].name << " -- " << S0[boost::target(e, S0)].name
                  << " [Weight: " << S0[e].weight << "]" << std::endl;
    }
    std::cout << std::endl;

    // S1: Child of S0, with vertices 'b' and 'c'
    SGraph& S1 = S0.create_subgraph(); // create_subgraph() is a member function of boost::subgraph
    // Free function calls for add_vertex (adding existing vertices to a subgraph)
    boost::add_vertex( (VertexDescriptor)1, S1); // Add global '1' ('b') from S0 to S1. Local descriptor is '0'.
    boost::add_vertex( (VertexDescriptor)2, S1); // Add global '2' ('c') from S0 to S1. Local descriptor is '1'.
    
    std::cout << "S1 (Child of S0) Edges:" << std::endl;
    for (auto e : boost::make_iterator_range(boost::edges(S1))) {
        std::cout << "  Local (" << boost::source(e, S1) << ", " << boost::target(e, S1)
                  << ") | Parent (" << S1.local_to_global(boost::source(e, S1))
                  << ", " << S1.local_to_global(boost::target(e, S1))
                  << ") | " << S1[boost::source(e, S1)].name << " -- " << S1[boost::target(e, S1)].name
                  << " [Weight: " << S1[e].weight << "]" << std::endl;
    }
    std::cout << std::endl;

    // S2: Child of S1, with vertex 'c'
    SGraph& S2 = S1.create_subgraph(); // create_subgraph() is a member function of boost::subgraph
    // To add 'c' from the S1 context, we need to get its local descriptor in S1.
    VertexDescriptor c_desc_in_S1 = S1.global_to_local((VertexDescriptor)2); // global_to_local is a member function
    // Free function call for add_vertex (adding existing vertex to a subgraph)
    boost::add_vertex(c_desc_in_S1, S2); // Add the vertex from S1's perspective.
    
    std::cout << "S2 (Child of S1) Edges:" << std::endl;
    if (boost::num_edges(S2) == 0) {
        std::cout << "  (no edges - 'c' is only vertex, no self-loops)" << std::endl;
    } else {
        for (auto e : boost::make_iterator_range(boost::edges(S2))) {
            std::cout << "  Local (" << boost::source(e, S2) << ", " << boost::target(e, S2)
                      << ") | Parent (" << S2.local_to_global(boost::source(e, S2))
                      << ", " << S2.local_to_global(boost::target(e, S2))
                      << ") | " << S2[boost::source(e, S2)].name << " -- " << S2[boost::target(e, S2)].name
                      << " [Weight: " << S2[e].weight << "]" << std::endl;
        }
    }
    std::cout << "S2 has " << boost::num_vertices(S2) << " vertex(es)." << std::endl;
    std::cout << "S2 local to S1: " << S2.local_to_global((VertexDescriptor)0) << std::endl;
    // To get to S0 (root) from S2, we must go S2 -> S1 -> S0
    std::cout << "S2 local to S0: " << S1.local_to_global(S2.local_to_global((VertexDescriptor)0)) << std::endl;
    
    return 0;
}
