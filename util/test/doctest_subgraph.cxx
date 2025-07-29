#include "WireCellUtil/Subgraph.h"
#include "WireCellUtil/doctest.h"

#include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/subgraph.hpp>

#include <iostream>

struct VertexBundle {
    std::string name;
};

struct EdgeBundle {
    double weight;
};

struct GraphBundle {
    std::string description;
};

template<typename VertexCollection, typename EdgeCollection>
void do_tests() {

    using BaseGraph = boost::adjacency_list<
        VertexCollection, EdgeCollection, boost::directedS,
        boost::property<boost::vertex_index_t, int, VertexBundle>,
        boost::property<boost::edge_index_t, int, EdgeBundle>,
        GraphBundle
    >;
    using Graph = boost::subgraph<BaseGraph>;

    Graph g;

    auto& gb = g.graph_bundle();
    gb.description = "Main Graph";
    std::cout << "Main graph description: " << gb.description << std::endl;
    REQUIRE(gb.description == "Main Graph");

    // Add vertices and edges with bundled properties to the main graph
    auto v1 = boost::add_vertex(g);
    g[v1].name = "Node A";
    auto v2 = boost::add_vertex(g);
    g[v2].name = "Node B";
    auto e1 = boost::add_edge(v1, v2, g).first;
    g[e1].weight = 1.5;

    std::cout << "Vertex 1 name: " << g[v1].name << std::endl;
    std::cout << "Edge 1 weight: " << g[e1].weight << std::endl;

    // Create a subgraph
    auto& sg = g.create_subgraph();
    // the subgraph is a distinct graph and thus has its own graph bundle
    REQUIRE(sg.graph_bundle().description == "");

    // Add vertices from main graph to subgraph (they are conceptually shared)
    auto sv1 = add_vertex(v1, sg);
    sg[sv1].name = "Subnode A";
    auto sv2 = add_vertex(v2, sg);

    // Properties are all the same so child subgraph will reset values on parent subgraph.
    REQUIRE(g[v1].name == "Subnode A");
    // Child sees same property as parent has.
    REQUIRE(sg[sv2].name == "Node B");

    auto [se1, have_edge] = sg.find_edge(e1);
    REQUIRE(have_edge);

    REQUIRE(sg[se1].weight == 1.5);

}

TEST_CASE("subgraph vv") {
    do_tests<boost::vecS, boost::vecS>();
}
TEST_CASE("subgraph sv") {
    do_tests<boost::setS, boost::vecS>();
}
TEST_CASE("subgraph vs") {
    do_tests<boost::vecS, boost::setS>();
}
TEST_CASE("subgraph s") {
    do_tests<boost::setS, boost::setS>();
}
