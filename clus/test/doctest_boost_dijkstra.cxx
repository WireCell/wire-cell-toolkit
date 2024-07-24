#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/graphviz.hpp> // fails to compile
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graphviz.hpp> // OK
#include <vector>
#include <iostream>

// Define graph and property types
struct TestEdgeProp {
    int dist;
};
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::no_property, TestEdgeProp> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;

TEST_CASE("standalone dijkstra") {
    Graph graph;
    // Add vertices and edges here
    // ...
    std::vector<Vertex> parents(boost::num_vertices(graph));
    std::vector<int> distances(boost::num_vertices(graph));
    Vertex v0 = 0;  // Starting vertex
    boost::dijkstra_shortest_paths(
        graph, v0,
        boost::weight_map(boost::get(&TestEdgeProp::dist, graph))
            .predecessor_map(&parents[0])
            .distance_map(&distances[0])
    );
    // Process results here
    // ...
}