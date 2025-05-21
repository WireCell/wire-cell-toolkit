#include "WireCellClus/Graph.h"      // boost graph
#include "WireCellUtil/GraphTools.h" // Mir/vertex_range

using namespace WireCell;
using namespace WireCell::Clus::Graph;
using namespace WireCell::GraphTools;
        
// Return value in map at key "name" or "name"
static std::string getname(const name_map_t& p, const std::string& name)
{
    auto it = p.find(name);
    if (it == p.end()) {
        return name;
    }
    return it->second;
}

bool Ident::load(const std::string& name, const store_t& store, Ident::graph_type& g,
                 const name_map_t& props)
{
    auto vit = store.find(vertices_name(name));
    if (vit == store.end()) {
        return false;
    }
    auto eit = store.find(edges_name(name));
    if (eit == store.end()) {
        return false;
    }
    
    g.clear();

    const auto& idents = vit->second.get(getname(props, "ident"))->elements<ident_type>();
    for (const auto& ident : idents) {
        boost::add_vertex(Ident::VertexProp{ident}, g);
    }
    
    const auto& weights = eit->second.get(getname(props, "weight"))->elements<weight_type>();
    const auto& tails = eit->second.get(getname(props, "tail"))->elements<vertex_descriptor>();
    const auto& heads = eit->second.get(getname(props, "head"))->elements<vertex_descriptor>();
    const size_t ne = weights.size();
    for(size_t ind=0; ind<ne; ++ind) {
        boost::add_edge(tails[ind], heads[ind], weights[ind], g);
    }
    return true;
}


bool Ident::save(const std::string& name, const Ident::graph_type& g, store_t& store,
                 const name_map_t& props)
{
    using Ident::ident_type;
    using Ident::weight_type;
    using Ident::vertex_descriptor;

    // Get the vertex/edge datasets, creating if missing.
    auto& vds = store[vertices_name(name)];
    auto& eds = store[edges_name(name)];

    size_t nv = boost::num_vertices(g);
    size_t ne = boost::num_edges(g);

    // It seems there is no way to directly iterate over values of a specific
    // property if the vertex/edge sets and so we bounce through temporary
    // vectors.

    auto idents = vds.allocate<ident_type>(getname(props, "ident"), {nv})->elements<ident_type>();

    auto weights = eds.allocate<weight_type>(getname(props, "weight"), {ne})->elements<weight_type>();

    auto tails = eds.allocate<vertex_descriptor>(getname(props, "tail"), {ne})->elements<vertex_descriptor>();
    auto heads = eds.allocate<vertex_descriptor>(getname(props, "head"), {ne})->elements<vertex_descriptor>();
    
    size_t ind=0;
    for (const auto& vtx : vertex_range(g)) {
        idents[ind++] = g[vtx].ident;
    }

    ind=0;
    for (const auto& edge : edge_range(g)) {
        weights[ind] = boost::get(boost::edge_weight_t(), g, edge);
        tails[ind] = boost::source(edge, g);
        heads[ind] = boost::target(edge, g);
        ++ind;
    }

    return true;
}



std::vector<Ident::ident_type> Ident::idents(const graph_type& g)
{
    std::vector<Ident::ident_type> ret;
    ret.reserve(boost::num_vertices(g));
    for (const auto& vtx : vertex_range(g)) {
        ret.push_back(g[vtx].ident);
    }
    return ret;
}

std::vector<Ident::weight_type> Ident::weights(const graph_type& g)
{
    std::vector<Ident::weight_type> ret;
    ret.reserve(boost::num_edges(g));
    auto weight_map = boost::get(boost::edge_weight, g);
    for (const auto& edge : mir(boost::edges(g))) {
        float weight = weight_map[edge];
        ret.push_back(weight);
    }
    return ret;
}

std::pair<std::vector<Ident::vertex_descriptor>, std::vector<Ident::vertex_descriptor>>
Ident::edges(const graph_type& g)
{
    std::pair<std::vector<Ident::vertex_descriptor>, std::vector<Ident::vertex_descriptor>> ret;
    const size_t nedges = boost::num_edges(g);
    ret.first.reserve(nedges);
    ret.second.reserve(nedges);
    for (const auto& edge : mir(boost::edges(g))) {
        ret.first.push_back(boost::source(edge, g));
        ret.second.push_back(boost::target(edge, g));
    }
    return ret;
}


DijkstraShortestPaths Ident::dijkstra_shortest_paths(const Ident::graph_type& graph, size_t point_index)
{
    const size_t nvtx = boost::num_vertices(graph);

    DijkstraShortestPaths ret{
        boost::vertex(point_index, graph),
        std::vector<vertex_descriptor>(nvtx),
        std::vector<int>(nvtex)
    };
    const auto& param = boost::weight_map(boost::get(boost::edge_weight, graph))
        .predecessor_map(&ret.parents[0])
        .distance_map(&ret.distances[0]);

    boost::dijkstra_shortest_paths(graph, ret.vertex, param);
    
    return ret;
}

