#include <iostream>
#include "WireCellUtil/Graph.h"



// Define the vertex properties
struct VertexProperties {
    int id;
    std::string name;
};

// Define the graph type
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexProperties> graph_t;

// Function to print the graph
template <typename Graph>
void printGraph(const Graph& graph)
{
    using vdesc_t = typename Graph::vertex_descriptor;
    // Iterate over the vertices of the graph
    auto vertexRange = boost::vertices(graph);
    for (auto it = vertexRange.first; it != vertexRange.second; ++it) {
        // Get the current vertex
        vdesc_t vertex = *it;

        // Print the vertex ID
        std::cout << "Vertex " << vertex << " " << graph[vertex].name << ": ";

        // Iterate over the adjacent vertices of the current vertex
        auto adjacentRange = boost::adjacent_vertices(vertex, graph);
        for (auto adjIt = adjacentRange.first; adjIt != adjacentRange.second; ++adjIt) {
            // Get the adjacent vertex
            vdesc_t adjVertex = *adjIt;

            // Print the adjacent vertex ID
            std::cout << adjVertex << " ";
        }

        std::cout << std::endl;
    }
}


int main()
{
    // Create a graph
    graph_t graph;

    // Add vertices
    graph_t::vertex_descriptor v1 = boost::add_vertex(VertexProperties{1, "V1"}, graph);
    graph_t::vertex_descriptor v2 = boost::add_vertex(VertexProperties{2, "V2"}, graph);
    graph_t::vertex_descriptor v3 = boost::add_vertex(VertexProperties{3, "V3"}, graph);

    // Add edges
    boost::add_edge(v1, v2, graph);
    boost::add_edge(v1, v3, graph);

    // Print the graph
    std::cout << "graph\n";
    printGraph(graph);

    using VFiltered =
        typename boost::filtered_graph<graph_t, boost::keep_all, std::function<bool(graph_t::vertex_descriptor)> >;
    VFiltered fg1(graph, {}, [&](auto v) { return v % 2 == 0;; });
    std::cout << "fg1\n";
    printGraph(fg1);

    graph_t copiedGraph;
    using desc_map_t = std::unordered_map<graph_t::vertex_descriptor, graph_t::vertex_descriptor>;
    using pm_desc_map_t = boost::associative_property_map<desc_map_t>;
    desc_map_t o2c;
    boost::copy_graph(fg1, copiedGraph, boost::orig_to_copy(pm_desc_map_t(o2c)));
    std::cout << "copiedGraph\n";
    printGraph(copiedGraph);
    for (const auto& [o, c] : o2c) {
        std::cout << o << " -> " << c << std::endl;
    }

    return 0;
}