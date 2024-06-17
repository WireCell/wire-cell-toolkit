#include <iostream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>

using namespace boost;

int main() {
    // Define the graph
    typedef adjacency_list<vecS, vecS, directedS> Graph;
    Graph g(5);

    // Add some edges to the graph
    add_edge(0, 1, g);
    add_edge(0, 2, g);
    add_edge(1, 3, g);
    add_edge(2, 3, g);
    add_edge(3, 4, g);

    // Define the visitor object
    class bfs_visitor : public default_bfs_visitor {
    public:
        void discover_vertex(int v, const Graph& g) const {
            std::cout << v << " ";
        }
    } vis;

    // Perform the breadth-first search starting from vertex 0
    breadth_first_search(g, 1, visitor(vis));
    // output: 1, 3, 4

    return 0;
}