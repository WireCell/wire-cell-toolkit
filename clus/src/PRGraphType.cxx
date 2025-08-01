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

}
