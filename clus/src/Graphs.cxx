#include "WireCellClus/Graphs.h"
#include "PAAL.h"


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

Weighted::Voronoi Weighted::voronoi(const Weighted::graph_type& graph,
                                    const std::vector<Weighted::vertex_type>& terminals)
{
    Voronoi result;
    const size_t npoints = boost::num_vertices(graph);
    auto index = get(boost::vertex_index, graph);

    result.terminal.resize(npoints); // nearest_terminal
    auto nearest_terminal_map = boost::make_iterator_property_map(result.terminal.begin(), index);
    for (auto terminal : terminals) {
        nearest_terminal_map[terminal] = terminal;
    }

    auto edge_weight = get(boost::edge_weight, graph);

    result.distance.resize(npoints);
    auto distance_map = boost::make_iterator_property_map(result.distance.begin(), index);

    result.last_edge.resize(npoints);
    auto last_edge = boost::make_iterator_property_map(result.last_edge.begin(), index);

    boost::dijkstra_shortest_paths(
        graph, terminals.begin(), terminals.end(),
        boost::dummy_property_map(),
        distance_map,
        edge_weight,
        index,
        PAAL::less(),
        boost::closed_plus<edge_weight_type>(),
        std::numeric_limits<edge_weight_type>::max(), 0,
        boost::make_dijkstra_visitor(
            PAAL::make_nearest_recorder(
                nearest_terminal_map, last_edge, boost::on_edge_relaxed{})));
    return result;
}

Weighted::GraphAlgorithms::GraphAlgorithms(const Graph& graph, size_t max_cache_size) 
    : m_graph(graph), m_max_cache_size(max_cache_size) 
{
    if (m_max_cache_size == 0) {
        m_max_cache_size = 1; // Ensure at least 1 entry can be cached
    }
}

void Weighted::GraphAlgorithms::update_cache_access(size_t source) const
{
    auto it = m_sps.find(source);
    if (it != m_sps.end()) {
        // Move to front of access order list (most recently used)
        m_access_order.erase(it->second.first);
        m_access_order.push_front(source);
        it->second.first = m_access_order.begin();
    }
}

void Weighted::GraphAlgorithms::evict_oldest_if_needed() const
{
    while (m_sps.size() >= m_max_cache_size) {
        // Remove least recently used (back of list)
        size_t oldest = m_access_order.back();
        m_access_order.pop_back();
        m_sps.erase(oldest);
    }
}

const Weighted::ShortestPaths&
Weighted::GraphAlgorithms::shortest_paths(size_t source) const
{
    auto it = m_sps.find(source);
    if (it != m_sps.end()) {
        // Cache hit - update access order
        update_cache_access(source);
        return it->second.second;
    }

    // Cache miss - need to evict if cache is full
    evict_oldest_if_needed();

    // Calculate shortest paths using Dijkstra
    const size_t nvtx = boost::num_vertices(m_graph);
    std::vector<size_t> predecessors(nvtx); 
    std::vector<Weighted::dijkstra_distance_type> distances(nvtx);

    const auto& param = weight_map(get(boost::edge_weight, m_graph))
                       .predecessor_map(&predecessors[0])
                       .distance_map(&distances[0]);
    boost::dijkstra_shortest_paths(m_graph, source, param);

    // Add to front of access order list
    m_access_order.push_front(source);
    
    // Insert into cache with iterator to list position
    auto result = m_sps.emplace(source, 
        std::make_pair(m_access_order.begin(), 
                      Weighted::ShortestPaths(source, predecessors)));

    return result.first->second.second;
}

const std::vector<size_t>&
Weighted::GraphAlgorithms::shortest_path(size_t source, size_t destination) const
{
    return shortest_paths(source).path(destination);
}

const std::vector<size_t>&
Weighted::GraphAlgorithms::connected_components() const
{
    if (m_cc.empty()) {
        m_cc.resize(boost::num_vertices(m_graph));
        boost::connected_components(m_graph, &m_cc[0]);
    }
    return m_cc;
}


void Weighted::GraphAlgorithms::clear_cache() const
{
    m_sps.clear();
    m_access_order.clear();
}


Weighted::filtered_graph_type Weighted::GraphAlgorithms::reduce(const vertex_set& vertices, bool accept) const
{
    auto filter = [&](vertex_type vtx) {
        return accept == (vertices.find(vtx) != vertices.end());
    };
    return Weighted::filtered_graph_type(m_graph, boost::keep_all(), filter);
}

Weighted::filtered_graph_type Weighted::GraphAlgorithms::reduce(const edge_set& edges, bool accept) const
{
    auto filter = [&](edge_type edge) {
        return accept == (edges.find(edge) != edges.end());
    };
    return Weighted::filtered_graph_type(m_graph, filter, boost::keep_all());
}
Weighted::filtered_graph_type Weighted::GraphAlgorithms::weight_threshold(double threshold, bool accept) const
{
    auto weight_map = get(boost::edge_weight, m_graph);
    auto filter = [&](edge_type edge) {
        return accept == (get(weight_map, edge) >= threshold);
    };
    return Weighted::filtered_graph_type(m_graph, filter, boost::keep_all());
}
