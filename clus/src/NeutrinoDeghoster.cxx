#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

namespace {
    // Helper function to sort clusters by total length in descending order
    bool sortbysec(const std::pair<Facade::Cluster*, double>& a,
                   const std::pair<Facade::Cluster*, double>& b) {
        return (a.second > b.second);
    }
}

void PatternAlgorithms::order_clusters(Graph& graph, std::vector<Facade::Cluster*>& ordered_clusters, std::map<Facade::Cluster*, std::vector<SegmentPtr> >& map_cluster_to_segments, std::map<Facade::Cluster*, double>& map_cluster_total_length){
    // Clear output containers
    map_cluster_to_segments.clear();
    map_cluster_total_length.clear();
    ordered_clusters.clear();
    
    // Iterate through all segments in the graph
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr seg = graph[*eit].segment;
        
        if (!seg || !seg->cluster()) continue;
        
        // Get the segment's cluster
        Facade::Cluster* cluster = seg->cluster();
        
        // Calculate segment length
        double length = segment_track_length(seg);
        
        // Check if this is the first segment for this cluster
        if (map_cluster_total_length.find(cluster) == map_cluster_total_length.end()) {
            // First segment for this cluster - initialize
            std::vector<SegmentPtr> segments;
            segments.push_back(seg);
            map_cluster_to_segments[cluster] = segments;
            map_cluster_total_length[cluster] = length;
        } else {
            // Add to existing cluster
            map_cluster_to_segments[cluster].push_back(seg);
            map_cluster_total_length[cluster] += length;
        }
    }
    
    // Create a vector of pairs (cluster, total_length) for sorting
    std::vector<std::pair<Facade::Cluster*, double>> temp_pair_vec;
    for (auto it = map_cluster_total_length.begin(); it != map_cluster_total_length.end(); ++it) {
        temp_pair_vec.push_back(std::make_pair(it->first, it->second));
    }
    
    // Sort clusters by total length in descending order
    std::sort(temp_pair_vec.begin(), temp_pair_vec.end(), sortbysec);
    
    // Fill ordered_clusters with sorted results
    for (auto it = temp_pair_vec.begin(); it != temp_pair_vec.end(); ++it) {
        ordered_clusters.push_back(it->first);
    }
}
