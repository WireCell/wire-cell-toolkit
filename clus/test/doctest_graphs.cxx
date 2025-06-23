#include "WireCellClus/Graphs.h"
#include "WireCellUtil/GraphTools.h"

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

using namespace WireCell::Clus::Graphs::Weighted;
using WireCell::GraphTools::edge_range;
using spdlog::debug;

// static
// filtered_graph_type reduce(const graph_type& graph, const vertex_set& vertices, bool accept) 
// {
//     auto filter = [&](vertex_type vtx) {
//         return accept == (vertices.find(vtx) != vertices.end());
//     };
//     return filtered_graph_type(graph, boost::keep_all(), filter);
// }

static
filtered_graph_type reduce_edges(const graph_type& graph, const edge_set& edges, bool accept=true);
static
filtered_graph_type reduce_edges(const graph_type& graph, const edge_set& edges, bool accept) 
{
    auto filter = [&](edge_type edge) {
        debug("check edge: {} -- {}", boost::source(edge, graph), boost::target(edge, graph));
        return accept == (edges.count(edge) > 0);
    };
    edge_predicate epred = filter;
    vertex_predicate vpred = boost::keep_all();
    return filtered_graph_type(graph, epred, vpred);
}

// static
// filtered_graph_type weight_threshold(const graph_type& graph, double threshold, bool accept) 
// {
//     auto weight_map = get(boost::edge_weight, graph);
//     auto filter = [&](edge_type edge) {
//         return accept == (get(weight_map, edge) >= threshold);
//     };
//     return filtered_graph_type(graph, filter, boost::keep_all());
// }


TEST_CASE("clus graphs")
{
    graph_type graph(3);
    REQUIRE(boost::num_vertices(graph) == 3);

    auto [e0,ok0] = boost::add_edge(0, 1, graph);
    REQUIRE(ok0);
    auto [e1,ok1] = boost::add_edge(0, 2, graph);
    REQUIRE(ok1);

    REQUIRE (boost::num_edges(graph) == 2);
    for (const auto& edge : edge_range(graph)) {
        debug("made edge: {} -- {}", boost::source(edge, graph), boost::target(edge, graph));
    }

    edge_set to_filter = {e0};

    auto filtered_with_e0 = reduce_edges(graph, to_filter);
    // CHECK (boost::num_edges(filtered_with_e0) == 1);
    // Gemini lied to me, num_edges returns original number, not filtered!
    CHECK (boost::num_edges(filtered_with_e0) == 2);
    for (const auto& edge : edge_range(filtered_with_e0)) {
        debug("have edge: {} -- {}", boost::source(edge, filtered_with_e0), boost::target(edge, filtered_with_e0));
    }
    CHECK (boost::num_edges(filtered_with_e0) == 2);
    auto filtered_without_e0 = reduce_edges(graph, to_filter, false);
    CHECK (boost::num_edges(filtered_without_e0) == 2);
}
