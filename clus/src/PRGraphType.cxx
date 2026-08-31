#include "WireCellClus/PRGraphType.h"
#include "WireCellUtil/GraphTools.h"

namespace WireCell::Clus::PR {

    node_vector graph_nodes(Graph& g) {
        return node_vector(boost::vertices(g).first, boost::vertices(g).second);
    }

    node_vector ordered_nodes(Graph& g) {
        auto nodes = graph_nodes(g);
        std::sort(nodes.begin(), nodes.end(), [&g](const node_descriptor& a, const node_descriptor& b) {
            return g[a].index < g[b].index;
        });
        return nodes;
    }

    std::vector<edge_descriptor> ordered_edges(Graph& g) {
        auto [ebegin, eend] = boost::edges(g);
        std::vector<edge_descriptor> edges(ebegin, eend);
        std::sort(edges.begin(), edges.end(), [&g](const edge_descriptor& a, const edge_descriptor& b) {
            return g[a].index < g[b].index;
        });
        return edges;
    }

}
