#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/GraphTools.h"

#include "WireCellAux/ClusterHelpers.h"

#include "WireCellUtil/Graph.h"

#include <unordered_map>

using namespace WireCell;
using spdlog::debug;

// Define a struct for the vertex properties
struct VertexProperties {
    size_t ident;
};
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, VertexProperties> graph_t;


namespace {
    // function make SimpleCluster with two connected blobs
    graph_t graph1()
    {
        // make a graph
        graph_t graph;
        boost::add_vertex({1}, graph);
        boost::add_vertex({2}, graph);
        // add edges
        boost::add_edge(0, 1, graph);
        return graph;
    }
    graph_t graph2()
    {
        // make a graph
        graph_t graph;
        boost::add_vertex({10}, graph);
        boost::add_vertex({20}, graph);
        // add edges
        boost::add_edge(0, 1, graph);
        return graph;
    }

    // print the graph into a string
    std::string print_graph(const graph_t& graph)
    {
        std::stringstream ss;
        for (auto const& vtx : GraphTools::mir(boost::vertices(graph))) {
            const auto& node = graph[vtx];
            ss << node.ident << " ";
        }
        ss << "\n";
        for (auto const& edge : GraphTools::mir(boost::edges(graph))) {
            const auto& src = boost::source(edge, graph);
            const auto& tgt = boost::target(edge, graph);
            ss << src << " " << tgt << "\n";
        }
        return ss.str();
    }
}


TEST_CASE("test merge cluster")
{
    spdlog::set_level(spdlog::level::debug); // Set global log level to debug
    const auto g1 = graph1();
    const auto g2 = graph2();
    const auto merged = Aux::merge_graphs<graph_t>({g2, g1});
    debug("g1 {}", print_graph(g1));
    debug("g2 {}", print_graph(g2));
    debug("merged {}", print_graph(merged));
}