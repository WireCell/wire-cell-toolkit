#include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/subgraph.hpp>

#include "WireCellUtil/Subgraph.h"
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

using BaseGraph = boost::adjacency_list<
        boost::setS, boost::setS, boost::directedS,
        boost::property<boost::vertex_index_t, int, VertexBundle>,
        boost::property<boost::edge_index_t, int, EdgeBundle>,
        GraphBundle
    >;
using Graph = boost::subgraph<BaseGraph>;

static
GraphBundle& get_graph_bundle(Graph& g)
{
    // subgraph appears to firewall graph level property bundle.  Standard
    // accesses do not work.

    // Accessing graph-level bundled property of the main graph
    // g[boost::graph_bundle].description = "Main Graph";
    // boost::get_property(boost::graph_bundle, g).description = "Main Graph";
    // GraphBundle& gb = g[boost::graph_bundle_type<BaseaGraph>()];
    // GraphBundle& gb = g.graph_bundle(); // no member named ....
    
    // Finally something that "works":
    return g.m_graph[boost::graph_bundle];

    // But, obviously this is VERY SKETCHY for future versions.  As of boost
    // 1.85, we see this in the subgraph.hpp source:
    //
    // public: // Probably shouldn't be public....
    //     Graph m_graph;
    //
    // At this point, I think we can only hope that if m_graph is made private
    // in the future there is a replacement.
}



int main() {
    Graph g;

    auto& gb = get_graph_bundle(g);
    gb.description = "Main Graph";
    std::cout << "Main graph description: " << gb.description << std::endl;

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

    auto& sgb = get_graph_bundle(sg);
    sgb.description = "Child subgraph 1";
    std::cout << "Subgraph description: " << sgb.description << std::endl;

    std::cout << "Subgraph description: " << sg.graph_bundle().description << std::endl;

    // Add vertices from main graph to subgraph (they are conceptually shared)
    auto sv1 = add_vertex(v1, sg);
    sg[sv1].name = "Subnode A";
    auto sv2 = add_vertex(v2, sg);
    auto se1 = add_edge(sv1, sv2, sg).first;
    sg[se1].weight = 2.5;

    // Accessing bundled properties of vertices and edges within the subgraph
    // These are the same objects as in the main graph
    std::cout << "Subgraph Vertex 1 name: " << sg[sv1].name << std::endl;
    std::cout << "Subgraph Edge 1 weight: " << sg[se1].weight << std::endl;
    std::cout << "Graph Vertex 1 name: " << g[v1].name << std::endl;
    std::cout << "Graph Edge 1 weight: " << g[e1].weight << std::endl;


    return 0;
}
