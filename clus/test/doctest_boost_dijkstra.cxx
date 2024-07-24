
// #include <boost/graph/graphviz.hpp>

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <vector>
#include <iostream>

// Define graph and property types

namespace TESTFacade {
    using namespace boost;
    struct VertexProp {
        int index;
    };
    struct EdgeProp {
        float dist;
    };
    typedef adjacency_list<vecS, vecS, undirectedS, VertexProp, EdgeProp> MCUGraph;
    typedef graph_traits<MCUGraph>::vertex_descriptor vertex_descriptor;
    typedef graph_traits<MCUGraph>::edge_descriptor edge_descriptor;

    namespace Cluster {
        std::vector<vertex_descriptor> parents;
        std::vector<int> distances;
    }  // namespace Cluster
}  // namespace Facade

using namespace TESTFacade;
using namespace TESTFacade::Cluster;
TEST_CASE("standalone dijkstra")
{
    MCUGraph graph(5);
    // Add vertices and edges here
    // ...
    parents.resize(boost::num_vertices(graph));
    distances.resize(boost::num_vertices(graph));
    vertex_descriptor v0 = vertex((size_t) 1, graph);  // Starting vertex
    const auto& param =
    boost::weight_map(boost::get(&EdgeProp::dist, graph))
    .predecessor_map(&parents[0])
    .distance_map(&distances[0]);
    std::cout << "param's C++ type: " << typeid(param).name() << std::endl;
    boost::dijkstra_shortest_paths(graph, v0, param);
}