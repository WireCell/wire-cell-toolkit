#include "WireCellClus/NeutrinoPatternBase.h"

using namespace WireCell::Clus::PR;

std::set<VertexPtr> PatternAlgorithms::find_vertices(Graph& graph, const Facade::Cluster& cluster)
{
    std::set<VertexPtr> result;
    
    // Iterate through all vertices in the graph
    auto [vbegin, vend] = boost::vertices(graph);
    for (auto vit = vbegin; vit != vend; ++vit) {
        VertexPtr vtx = graph[*vit].vertex;
        
        // Check if this vertex belongs to the specified cluster
        if (vtx && vtx->cluster() && vtx->cluster() == &cluster) {
            result.insert(vtx);
        }
    }
    
    return result;
}

std::set<SegmentPtr> PatternAlgorithms::find_segments(Graph& graph, const Facade::Cluster& cluster)
{
    std::set<SegmentPtr> result;
    
    // Iterate through all edges (segments) in the graph
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr seg = graph[*eit].segment;
        
        // Check if this segment belongs to the specified cluster
        if (seg && seg->cluster() && seg->cluster() == &cluster) {
            result.insert(seg);
        }
    }
    
    return result;
}

bool PatternAlgorithms::clean_up_graph(Graph& graph, const Facade::Cluster& cluster)
{
    bool modified = false;
    
    // First, find and remove all segments associated with this cluster
    std::set<SegmentPtr> segments_to_remove = find_segments(graph, cluster);
    for (auto seg : segments_to_remove) {
        if (remove_segment(graph, seg)) {
            modified = true;
        }
    }
    
    // Then, find and remove all vertices associated with this cluster
    // Note: vertices that are still connected to other segments won't be removed
    // until their segments are removed first
    std::set<VertexPtr> vertices_to_remove = find_vertices(graph, cluster);
    for (auto vtx : vertices_to_remove) {
        if (remove_vertex(graph, vtx)) {
            modified = true;
        }
    }
    
    return modified;
}