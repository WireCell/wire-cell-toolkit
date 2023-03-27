#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/String.h"

#include <iostream>
#include <fstream>

using namespace WireCell::GraphTools;
using namespace WireCell::String;

struct Node {
    char code;
    int id;
};

struct Edge {
    size_t count;
};

using ugraph_t = boost::adjacency_list<boost::setS, boost::vecS,
                                       boost::undirectedS, Node, Edge>;
using uvdesc_t = boost::graph_traits<ugraph_t>::vertex_descriptor;
using uedesc_t = boost::graph_traits<ugraph_t>::edge_descriptor;


int main()
{
    ugraph_t gr;
    auto s0 = boost::add_vertex({'s', 0}, gr);
    auto s1 = boost::add_vertex({'s', 1}, gr);
    auto b0 = boost::add_vertex({'b', 0}, gr);
    auto w0 = boost::add_vertex({'w', 0}, gr);
    auto b1 = boost::add_vertex({'b', 1}, gr);
    auto b2 = boost::add_vertex({'b', 2}, gr);
    auto w1 = boost::add_vertex({'w', 1}, gr);
    auto w2 = boost::add_vertex({'w', 2}, gr);
    auto w3 = boost::add_vertex({'w', 3}, gr);

    boost::add_edge(b0, s0, gr);
    boost::add_edge(b0, w0, gr);
    boost::add_edge(b0, w2, gr);
    boost::add_edge(s0, b1, gr);
    boost::add_edge(b1, w1, gr);
    boost::add_edge(b0, w3, gr);
    boost::add_edge(b1, w3, gr);
    boost::add_edge(s1, b2, gr);
    boost::add_edge(b2, w0, gr);
    boost::add_edge(b2, w1, gr);

    std::cerr << dotify(gr) << std::endl;


    return 0;
}
