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

    // Add vertices
    Vertex v0 = boost::add_vertex(graph);
    Vertex v1 = boost::add_vertex(graph);
    Vertex v2 = boost::add_vertex(graph);
    Vertex v3 = boost::add_vertex(graph);
    Vertex v4 = boost::add_vertex(graph);

    // Add edges
    Edge e0 = boost::add_edge(v0, v1, TestEdgeProp{1}, graph).first;
    Edge e1 = boost::add_edge(v1, v4, TestEdgeProp{1}, graph).first;
    Edge e2 = boost::add_edge(v0, v2, TestEdgeProp{1}, graph).first;
    Edge e3 = boost::add_edge(v2, v4, TestEdgeProp{2}, graph).first;

    std::vector<Vertex> parents(boost::num_vertices(graph));
    std::vector<int> distances(boost::num_vertices(graph));
    boost::dijkstra_shortest_paths(
        graph, v0,
        boost::weight_map(boost::get(&TestEdgeProp::dist, graph))
            .predecessor_map(&parents[0])
            .distance_map(&distances[0])
    );
    std::cout << "parents: ";
    for (size_t i = 0; i < parents.size(); i++) {
        std::cout << i << "->" << parents[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "distances: ";
    for (size_t i = 0; i < distances.size(); i++) {
        std::cout << i << "->" << distances[i] << " ";
    }
    std::cout << std::endl;
}
