#include "WireCellUtil/Logging.h"


// Boost Graph Library includes
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp> // For print_graph
#include <boost/graph/filtered_graph.hpp> // For filtered_graph

#include <iostream>
#include <set>
#include <functional> // For std::function


using spdlog::debug;

// --- Your Boost Graph Related Types ---
using edge_weight_type = double;
using Graph = boost::adjacency_list<
    boost::vecS,          // vertices
    boost::vecS,          // edges
    boost::undirectedS,   // edge direction (none)
    boost::property<boost::vertex_index_t, size_t>,
    boost::property<boost::edge_weight_t, edge_weight_type>
    >;
using graph_type = Graph;
using vertex_type = boost::graph_traits<Graph>::vertex_descriptor;
using edge_type = boost::graph_traits<Graph>::edge_descriptor;

// A set of unique vertices or edges;
using vertex_set = std::set<vertex_type>;
using edge_set = std::set<edge_type>;

// Filtered graphs and their predicates.
using vertex_predicate = std::function<bool(vertex_type)>;
using edge_predicate = std::function<bool(edge_type)>;
using filtered_graph_type = boost::filtered_graph<graph_type, edge_predicate, vertex_predicate>;

// --- Your reduce_edges Function ---
filtered_graph_type reduce_edges(const graph_type& graph, const edge_set& edges, bool accept)
{
    // The 'filter' lambda captures 'edges' and 'accept' by reference.
    // This lambda will be called for each edge when iterating the filtered graph.
    auto filter = [&](edge_type edge) {
        // This debug statement will only print if the predicate is actually invoked.
        debug("  [Predicate Call] Checking edge: {} -- {}",
              boost::source(edge, graph), boost::target(edge, graph));
        return accept == (edges.count(edge) > 0);
    };

    edge_predicate epred = filter;
    vertex_predicate vpred = boost::keep_all(); // Keep all vertices regardless of edges

    // Construct the filtered graph view.
    return filtered_graph_type(graph, epred, vpred);
}

int main()
{
    // Create a sample graph
    graph_type g(5); // 5 vertices (0, 1, 2, 3, 4)

    // Add some edges
    edge_type e01, e02, e12, e13, e24;
    bool b;

    boost::tie(e01, b) = add_edge(0, 1, 1.0, g);
    boost::tie(e02, b) = add_edge(0, 2, 2.0, g);
    boost::tie(e12, b) = add_edge(1, 2, 3.0, g);
    boost::tie(e13, b) = add_edge(1, 3, 4.0, g);
    boost::tie(e24, b) = add_edge(2, 4, 5.0, g);

    std::cout << "Original Graph Edges:" << std::endl;
    // Iterate and print edges of the original graph to confirm setup
    auto edge_it_orig = boost::edges(g);
    for (; edge_it_orig.first != edge_it_orig.second; ++edge_it_orig.first) {
        std::cout << "  " << boost::source(*edge_it_orig.first, g) << " -- "
                  << boost::target(*edge_it_orig.first, g) << std::endl;
    }
    std::cout << "Number of edges in original graph: " << boost::num_edges(g) << std::endl << std::endl;


    // --- Test Case 1: Accept only a specific set of edges ---
    edge_set edges_to_retain;
    edges_to_retain.insert(e01);
    edges_to_retain.insert(e13);

    std::cout << "Reducing edges (accepting e01, e13):" << std::endl;
    filtered_graph_type fg_accept = reduce_edges(g, edges_to_retain, true);

    // ************* IMPORTANT: THIS IS WHERE THE PREDICATE IS INVOKED *************
    // Iterate over the edges of the filtered graph to trigger predicate calls
    std::cout << "Edges in filtered graph (accepting):" << std::endl;
    auto edge_it_filtered_accept = boost::edges(fg_accept);
    for (; edge_it_filtered_accept.first != edge_it_filtered_accept.second; ++edge_it_filtered_accept.first) {
        std::cout << "  " << boost::source(*edge_it_filtered_accept.first, fg_accept) << " -- "
                  << boost::target(*edge_it_filtered_accept.first, fg_accept) << std::endl;
    }
    std::cout << "Number of edges in filtered graph (accepting): " << boost::num_edges(fg_accept) << std::endl << std::endl;


    // --- Test Case 2: Exclude a specific set of edges ---
    edge_set edges_to_exclude;
    edges_to_exclude.insert(e01); // Exclude 0-1
    edges_to_exclude.insert(e24); // Exclude 2-4

    std::cout << "Reducing edges (excluding e01, e24):" << std::endl;
    filtered_graph_type fg_exclude = reduce_edges(g, edges_to_exclude, false);

    // ************* IMPORTANT: THIS IS WHERE THE PREDICATE IS INVOKED *************
    // Iterate over the edges of the filtered graph to trigger predicate calls
    std::cout << "Edges in filtered graph (excluding):" << std::endl;
    auto edge_it_filtered_exclude = boost::edges(fg_exclude);
    for (; edge_it_filtered_exclude.first != edge_it_filtered_exclude.second; ++edge_it_filtered_exclude.first) {
        std::cout << "  " << boost::source(*edge_it_filtered_exclude.first, fg_exclude) << " -- "
                  << boost::target(*edge_it_filtered_exclude.first, fg_exclude) << std::endl;
    }
    std::cout << "Number of edges in filtered graph (excluding): " << boost::num_edges(fg_exclude) << std::endl << std::endl;

    return 0;
}

