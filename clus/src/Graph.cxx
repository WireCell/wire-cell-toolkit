#include "WireCellClus/Graph.h"      // boost graph
#include "WireCellUtil/GraphTools.h" // mir/vertex_range

using namespace WireCell::Clus::Graph;
using namespace WireCell::GraphTools;
        
std::vector<int> Ident::idents(const graph_type& g)
{
    std::vector<int> ret;
    ret.reserve(boost::num_vertices(g));
    for (const auto& vtx : vertex_range(g)) {
        ret.push_back(g[vtx].ident);
    }
    return ret;
}

std::vector<float> Ident::weights(const graph_type& g)
{
    std::vector<float> ret;
    ret.reserve(boost::num_edges(g));
    auto weight_map = boost::get(boost::edge_weight, g);
    for (const auto& edge : mir(boost::edges(g))) {
        float weight = weight_map[edge];
        ret.push_back(weight);
    }
    return ret;
}

std::pair<std::vector<int>, std::vector<int>>
Ident::edges(const graph_type& g)
{
    std::pair<std::vector<int>, std::vector<int>> ret;
    const size_t nedges = boost::num_edges(g);
    ret.first.reserve(nedges);
    ret.second.reserve(nedges);
    for (const auto& edge : mir(boost::edges(g))) {
        ret.first.push_back(boost::source(edge, g));
        ret.second.push_back(boost::target(edge, g));
    }
    return ret;
}


