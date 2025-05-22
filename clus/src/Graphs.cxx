#include "WireCellClus/Graphs.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Graphs;

Weighted::ShortestPaths::ShortestPaths(size_t source, const std::vector<size_t> predecessors)
    : m_source(source)
    , m_predecessors(predecessors)
{
}

const std::vector<size_t>&
Weighted::ShortestPaths::path(size_t destination) const
{
    auto& path = m_paths[destination]; // construct and hold
    if (path.size()) {
        return path;
    }    

    path.push_back(destination);
    size_t prev = destination;
    for (size_t vertex = m_predecessors[destination]; vertex != m_source; vertex = m_predecessors[vertex]) 
    {
        path.push_back(vertex);
        if (vertex == prev) {
            break;
        }
        prev = vertex;
    }
    path.push_back(m_source);
    std::reverse(path.begin(), path.end());

    return path;
}

Weighted::ShortestPathsGraph::ShortestPathsGraph(GraphPtr&& graph) : m_graph(std::move(graph)) {}

const Weighted::ShortestPaths&
Weighted::ShortestPathsGraph::paths(size_t source) const
{
    auto it = m_sps.find(source);
    if (it != m_sps.end()) {
        return it->second;
    }

    const size_t nvtx = boost::num_vertices(*m_graph);
    std::vector<size_t> predecessors(nvtx); 
    std::vector<double> distances(nvtx); // ignore

    const auto& param = weight_map(get(boost::edge_weight, *m_graph))
				   .predecessor_map(&predecessors[0])
				   .distance_map(&distances[0]);
    boost::dijkstra_shortest_paths(*m_graph, source, param);

    auto got = m_sps.emplace(source, Weighted::ShortestPaths(source, predecessors));
    return got.first->second;
}

const std::vector<size_t>&
Weighted::ShortestPathsGraph::path(size_t source, size_t destination) const
{
    return paths(source).path(destination);
}

