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

void PatternAlgorithms::determine_direction(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model) {
    // Iterate through all edges (segments) in the graph
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr seg = graph[*eit].segment;
        
        // Skip if segment is null or doesn't belong to this cluster
        if (!seg || seg->cluster() != &cluster) continue;
        
        // Get the two vertices of this segment
        auto [start_v, end_v] = find_vertices(graph, seg);
        if (!start_v || !end_v) {
            std::cout << "Error in finding vertices for a segment" << std::endl;
            continue;
        }
        
        // Check if vertices match the segment endpoints (start_v should be at front, end_v at back)
        const auto& wcpts = seg->wcpts();
        if (wcpts.size() < 2) continue;
        
        auto front_pt = wcpts.front().point;
        auto back_pt = wcpts.back().point;
        
        // Determine which vertex is start and which is end based on point positions
        double dis_sv_front = ray_length(Ray{start_v->wcpt().point, front_pt});
        double dis_sv_back = ray_length(Ray{start_v->wcpt().point, back_pt});
        
        if (dis_sv_front > dis_sv_back) {
            std::swap(start_v, end_v);
        }
        
        // Count number of segments connected to each vertex
        int start_n = 0, end_n = 0;
        if (start_v->descriptor_valid()) {
            start_n = boost::degree(start_v->get_descriptor(), graph);
        }
        if (end_v->descriptor_valid()) {
            end_n = boost::degree(end_v->get_descriptor(), graph);
        }
        
        bool flag_print = false;
        // if (seg->cluster() == main_cluster) flag_print = true;
        
        if (seg->flags_any(SegmentFlags::kShowerTrajectory)) {
            // Trajectory shower
            segment_determine_shower_direction_trajectory(seg, start_n, end_n, particle_data, recomb_model, 43000/units::cm, flag_print);
        } else if (seg->flags_any(SegmentFlags::kShowerTopology)) {
            // Topology shower
            segment_determine_shower_direction(seg, particle_data, recomb_model);
        } else {
            // Track
            segment_determine_dir_track(seg, start_n, end_n, particle_data, recomb_model, 43000/units::cm, flag_print);
        }
    }
}

