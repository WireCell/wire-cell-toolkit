#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

void PatternAlgorithms::clustering_points(Graph& graph, Facade::Cluster& cluster, const IDetectorVolumes::pointer& dv, const std::string& cloud_name, double search_range, double scaling_2d){
    // Collect all segments that belong to this cluster
    std::set<SegmentPtr> segments;
    
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr seg = graph[*eit].segment;
        if (seg && seg->cluster() == &cluster) {
            segments.insert(seg);
        }
    }
    
    // Run clustering on the collected segments
    if (!segments.empty()) {
        clustering_points_segments(segments, dv, cloud_name, search_range, scaling_2d);
    }
}

void PatternAlgorithms::separate_track_shower(Graph&graph, Facade::Cluster& cluster) {
    // Iterate through all edges (segments) in the graph
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr seg = graph[*eit].segment;
        
        // Skip if segment is null or doesn't belong to this cluster
        if (!seg || seg->cluster() != &cluster) continue;
        
        // First check if segment is a shower topology
        segment_is_shower_topology(seg);
        
        // If not shower topology, check if it's a shower trajectory
        if (!seg->flags_any(SegmentFlags::kShowerTopology)) {
            segment_is_shower_trajectory(seg);
        }
    }
}
