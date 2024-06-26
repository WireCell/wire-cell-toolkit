#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <iostream>

// Define the graph type
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

// Define the filter for vertices
struct VertexFilter {
    bool operator()(const Graph::vertex_descriptor& v) const {
        return v % 2 == 0; // Keep only even vertices
    }
};

// Define the filter for edges
struct EdgeFilter {
    EdgeFilter(const Graph& g) : graph(g) {}

    bool operator()(const Graph::edge_descriptor& e) const {
        Graph::vertex_descriptor src = boost::source(e, graph);
        // Graph::vertex_descriptor tgt = boost::target(e, graph);

        return src % 2 == 0; // Keep only edges connecting even src vertices
    }

    const Graph& graph;
};

template<class Gr>
void print(const Gr& graph) {
    // Print the filtered graph
    for (auto ed : boost::make_iterator_range(boost::edges(graph))) {
        std::cout << ed << " ";
    }
    for (auto vt : boost::make_iterator_range(boost::vertices(graph))) {
        std::cout << vt << " ";
    }
    std::cout << std::endl;
}

int main() {
    // Create a sample graph
    Graph g(6);
    boost::add_edge(0, 1, g);
    boost::add_edge(0, 2, g);
    boost::add_edge(1, 2, g);

    // Create the filtered graph
    VertexFilter vertex_filter;
    EdgeFilter edge_filter(g);
    boost::filtered_graph<Graph, boost::keep_all, VertexFilter> fg_v(g, {}, vertex_filter);
    boost::filtered_graph<Graph, EdgeFilter, boost::keep_all> fg_e(g, edge_filter, {});
    
    std::cout << "input graph:" << std::endl;
    print(g);
    std::cout << "filtered graph vertex:" << std::endl;
    print(fg_v);
    std::cout << "filtered graph edge:" << std::endl;
    print(fg_e);

    return 0;
}