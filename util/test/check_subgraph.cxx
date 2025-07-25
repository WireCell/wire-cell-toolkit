#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/graph_utility.hpp>

// Define the bundled properties
struct VertexProperties {
    char name;
};

struct EdgeProperties {
    double weight;
};

using BaseGraph = boost::adjacency_list<
    boost::vecS,                // vertex
    boost::vecS,                // edge
    boost::undirectedS,
    boost::property<boost::vertex_index_t, int, VertexProperties>,
    boost::property<boost::edge_index_t, int, EdgeProperties>
    >;

using MainGraph = boost::subgraph<BaseGraph>;
using SubGraph = MainGraph;

using VertexDescriptor = typename boost::graph_traits<MainGraph>::vertex_descriptor;
using EdgeDescriptor = typename boost::graph_traits<MainGraph>::edge_descriptor;


// Helper function to print edges with properties and descriptor info
template <typename G>
void print_graph_edges_with_properties(const G& g, const std::string& title) {
    std::cout << "--- " << title << " ---" << std::endl;
    if (boost::num_edges(g) == 0) {
        std::cout << "  (no edges)" << std::endl;
    } else {

        for (auto e : boost::make_iterator_range(boost::edges(g))) {
            VertexDescriptor u = boost::source(e, g);
            VertexDescriptor v = boost::target(e, g);
            char u_name = g[u].name;
            char v_name = g[v].name;
            double weight = g[e].weight;
            
            // Demonstrate getting the global descriptor from a local one
            if (!g.is_root()) {
                VertexDescriptor u_global = g.local_to_global(u);
                VertexDescriptor v_global = g.local_to_global(v);
                std::cout << "  Local (" << u << ", " << v << ") -> Global (" << u_global << ", " << v_global << ") | "
                          << u_name << " -- " << v_name << " [Weight: " << weight << "]" << std::endl;
            } else {
                std::cout << "  (" << u << ", " << v << ") | "
                          << u_name << " -- " << v_name << " [Weight: " << weight << "]" << std::endl;
            }
        }
    }
    std::cout << std::endl;
}

// External function to perform the "break edge" operation
template <typename G>
void break_edge(G& g, typename boost::graph_traits<G>::edge_descriptor e) {

    VertexDescriptor u = boost::source(e, g);
    VertexDescriptor v = boost::target(e, g);
    double old_weight = g[e].weight;

    boost::remove_edge(e, g);

    VertexDescriptor x = boost::add_vertex(g);
    g[x].name = 'x';

    auto e1 = boost::add_edge(u, x, g).first;
    auto e2 = boost::add_edge(x, v, g).first;
    g[e1].weight = old_weight / 2.0;
    g[e2].weight = old_weight / 2.0;
}

int main() {
    MainGraph main_graph(4);
    main_graph[(VertexDescriptor)0].name = 'a';
    main_graph[(VertexDescriptor)1].name = 'b';
    main_graph[(VertexDescriptor)2].name = 'c';
    main_graph[(VertexDescriptor)3].name = 'd';

    auto e_ab = boost::add_edge(0, 1, main_graph).first;
    auto e_bc = boost::add_edge(1, 2,  main_graph).first;
    auto e_cd = boost::add_edge(2, 3, main_graph).first;
    main_graph[e_ab].weight = 10.0;
    main_graph[e_bc].weight = 20.0;
    main_graph[e_cd].weight = 30.0;


    print_graph_edges_with_properties(main_graph, "Initial Main Graph Edges (Global Descriptors)");

    // Create a subgraph for the "b -- c -- d" path
    SubGraph& sub_graph = main_graph.create_subgraph();
    
    // Add global vertex descriptors to the subgraph
    boost::add_vertex(1, sub_graph); // Global 'b' becomes local '0'
    boost::add_vertex(2, sub_graph); // Global 'c' becomes local '1'
    boost::add_vertex(3, sub_graph); // Global 'd' becomes local '2'

    print_graph_edges_with_properties(sub_graph, "Initial Subgraph Edges (Local Descriptors)");

    // Let's get the local descriptor for edge 'b -- c'
    // First, find the global descriptor for the edge from the main graph.
    // The main_graph is a subgraph, so we use global_to_local to convert
    // the global index to a local descriptor for this specific root subgraph.
    // In this simple case, global descriptor IS the local descriptor for the root.
    auto global_b = 1;
    auto global_c = 2;
    auto local_e_bc = boost::edge(global_b, global_c, sub_graph).first;

    // Now, we break the edge from within the context of the subgraph.
    // The `break_edge` function will modify the underlying main_graph.
    std::cout << "Performing 'break edge' on 'b -- c' from the subgraph's perspective..." << std::endl;
    break_edge(sub_graph, local_e_bc);

    // Demonstrate the new state of both the main graph and the subgraph
    print_graph_edges_with_properties(main_graph, "Main Graph after 'break edge'");
    print_graph_edges_with_properties(sub_graph, "Subgraph after 'break edge'");

    return 0;
}
