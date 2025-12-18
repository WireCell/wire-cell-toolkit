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

std::pair<int, double> PatternAlgorithms::calculate_num_daughter_showers(Graph& graph, VertexPtr vertex, SegmentPtr segment, bool flag_count_shower) {
    int number_showers = 0;
    double acc_length = 0;
    
    std::set<VertexPtr> used_vertices;
    std::set<SegmentPtr> used_segments;
    
    std::vector<std::pair<VertexPtr, SegmentPtr>> segments_to_be_examined;
    segments_to_be_examined.push_back(std::make_pair(vertex, segment));
    used_vertices.insert(vertex);
    
    while(segments_to_be_examined.size() > 0) {
        std::vector<std::pair<VertexPtr, SegmentPtr>> temp_segments;
        for (auto it = segments_to_be_examined.begin(); it != segments_to_be_examined.end(); it++) {
            VertexPtr prev_vtx = it->first;
            SegmentPtr current_sg = it->second;
            
            if (used_segments.find(current_sg) != used_segments.end()) continue; // looked at it before
            
            // Check if segment is a shower (has kShowerTrajectory or kShowerTopology flags)
            bool is_shower = current_sg->flags_any(SegmentFlags::kShowerTrajectory) || 
                           current_sg->flags_any(SegmentFlags::kShowerTopology);
            
            if (is_shower || (!flag_count_shower)) {
                number_showers++;
                acc_length += segment_track_length(current_sg);
            }
            used_segments.insert(current_sg);
            
            VertexPtr curr_vertex = find_other_vertex(graph, current_sg, prev_vtx);
            if (used_vertices.find(curr_vertex) != used_vertices.end()) continue;
            
            // Get all segments connected to curr_vertex
            if (curr_vertex && curr_vertex->descriptor_valid()) {
                auto vd = curr_vertex->get_descriptor();
                auto edge_range = boost::out_edges(vd, graph);
                for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
                    SegmentPtr seg = graph[*eit].segment;
                    if (seg) {
                        temp_segments.push_back(std::make_pair(curr_vertex, seg));
                    }
                }
            }
            used_vertices.insert(curr_vertex);
        }
        segments_to_be_examined = temp_segments;
    }
    
    return std::make_pair(number_showers, acc_length);
} 

void PatternAlgorithms::examine_good_tracks(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data) {
    // Iterate through all edges (segments) in the graph
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr sg = graph[*eit].segment;
        
        // Skip if segment is null or doesn't belong to this cluster
        if (!sg || sg->cluster() != &cluster) continue;
        
        // Skip if segment is a shower
        if (sg->flags_any(SegmentFlags::kShowerTrajectory) || sg->flags_any(SegmentFlags::kShowerTopology)) continue;
        
        // Skip if no direction or weak direction
        if (sg->dirsign() == 0 || sg->dir_weak()) continue;
        
        // Get the two vertices of this segment
        auto [vertex1, vertex2] = find_vertices(graph, sg);
        if (!vertex1 || !vertex2) continue;
        
        // Determine start and end vertices based on segment direction
        VertexPtr start_vertex = nullptr, end_vertex = nullptr;
        const auto& wcpts = sg->wcpts();
        if (wcpts.size() < 2) continue;
        
        auto front_pt = wcpts.front().point;
        auto back_pt = wcpts.back().point;
        
        if (sg->dirsign() == 1) {
            // Direction is forward (from front to back)
            double dis1_front = ray_length(Ray{vertex1->wcpt().point, front_pt});
            double dis1_back = ray_length(Ray{vertex1->wcpt().point, back_pt});
            if (dis1_front < dis1_back) {
                start_vertex = vertex1;
                end_vertex = vertex2;
            } else {
                start_vertex = vertex2;
                end_vertex = vertex1;
            }
        } else if (sg->dirsign() == -1) {
            // Direction is backward (from back to front)
            double dis1_front = ray_length(Ray{vertex1->wcpt().point, front_pt});
            double dis1_back = ray_length(Ray{vertex1->wcpt().point, back_pt});
            if (dis1_front < dis1_back) {
                start_vertex = vertex2;
                end_vertex = vertex1;
            } else {
                start_vertex = vertex1;
                end_vertex = vertex2;
            }
        }
        
        if (!start_vertex || !end_vertex) continue;
        
        // Calculate number of daughter showers
        auto result_pair = calculate_num_daughter_showers(graph, start_vertex, sg);
        int num_daughter_showers = result_pair.first;
        double length_daughter_showers = result_pair.second;
        
        // Calculate maximum angle between this segment and others at end_vertex
        double max_angle = 0;
        WireCell::Point end_pt = end_vertex->fit().valid() ? end_vertex->fit().point : end_vertex->wcpt().point;
        WireCell::Vector dir1 = segment_cal_dir_3vector(sg, end_pt, 15*units::cm);
        WireCell::Vector drift_dir(1, 0, 0);
        double min_para_angle = 1e9;
        
        // Get all segments connected to end_vertex
        if (end_vertex->descriptor_valid()) {
            auto vd = end_vertex->get_descriptor();
            auto edge_range = boost::out_edges(vd, graph);
            for (auto eit2 = edge_range.first; eit2 != edge_range.second; ++eit2) {
                SegmentPtr sg1 = graph[*eit2].segment;
                if (!sg1 || sg1 == sg) continue;
                
                WireCell::Vector dir2 = segment_cal_dir_3vector(sg1, end_pt, 15*units::cm);
                double angle = std::acos(std::min(1.0, std::max(-1.0, dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude())))) / 3.1415926 * 180.0;
                if (angle > max_angle) max_angle = angle;
                
                angle = std::fabs(std::acos(std::min(1.0, std::max(-1.0, drift_dir.dot(dir2) / (drift_dir.magnitude() * dir2.magnitude())))) / 3.1415926 * 180.0 - 90.0);
                if (angle < min_para_angle) min_para_angle = angle;
            }
        }
        
        // Check if this track should be reclassified as an electron shower
        double drift_angle = std::fabs(std::acos(std::min(1.0, std::max(-1.0, drift_dir.dot(dir1) / (drift_dir.magnitude() * dir1.magnitude())))) / 3.1415926 * 180.0 - 90.0);
        double length = segment_track_length(sg);
        
        if ((num_daughter_showers >= 4 || (length_daughter_showers > 50*units::cm && num_daughter_showers >= 2)) &&
            (max_angle > 155 || (drift_angle < 15 && min_para_angle < 15 && min_para_angle + drift_angle < 25)) &&
            length < 15*units::cm) {
            
            // Reclassify as electron (PDG 11)
            auto pinfo = std::make_shared<Aux::ParticleInfo>(
                11,                                              // electron PDG
                particle_data->get_particle_mass(11),           // electron mass
                particle_data->pdg_to_name(11),                 // "e-"
                WireCell::D4Vector<double>(0, 0, 0, 0)         // zero 4-momentum (will be recalculated)
            );
            sg->particle_info(pinfo);
            
            // Reset direction and mark as weak
            sg->dirsign(0);
            sg->dir_weak(true);
        }
        
        // Debug output (commented out)
        // std::cout << sg->get_id() << " " << sg->particle_type() << " " << num_daughter_showers << " " 
        //           << length/units::cm << " " << max_angle << " " << min_para_angle << " " << drift_angle << std::endl;
    }
}

void PatternAlgorithms::fix_maps_multiple_tracks_in(Graph& graph, Facade::Cluster& cluster){
    // Iterate through all vertices in the graph
    auto [vbegin, vend] = boost::vertices(graph);
    for (auto vit = vbegin; vit != vend; ++vit) {
        VertexPtr vtx = graph[*vit].vertex;
        
        // Skip if vertex is null or doesn't belong to this cluster
        if (!vtx || !vtx->cluster() || vtx->cluster() != &cluster) continue;
        
        // Check how many segments are connected to this vertex
        if (!vtx->descriptor_valid()) continue;
        auto vd = vtx->get_descriptor();
        if (boost::degree(vd, graph) <= 1) continue;
        
        int n_in = 0;
        int n_in_shower = 0;
        std::vector<SegmentPtr> in_tracks;
        
        // Get vertex position
        WireCell::Point vtx_point = vtx->wcpt().point;
        
        // Iterate through all segments connected to this vertex
        auto edge_range = boost::out_edges(vd, graph);
        for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
            SegmentPtr sg = graph[*eit].segment;
            if (!sg) continue;
            
            // Determine if this vertex is at the front or back of the segment
            const auto& wcpts = sg->wcpts();
            if (wcpts.size() < 2) continue;
            
            auto front_pt = wcpts.front().point;
            auto back_pt = wcpts.back().point;
            
            double dis_front = ray_length(Ray{vtx_point, front_pt});
            double dis_back = ray_length(Ray{vtx_point, back_pt});
            
            bool flag_start = (dis_front < dis_back); // vertex is at the front of segment
            
            // Check if this segment is pointing "in" to the vertex
            // "in" means: (at front and direction is -1) OR (at back and direction is 1)
            if ((flag_start && sg->dirsign() == -1) || (!flag_start && sg->dirsign() == 1)) {
                n_in++;
                
                // Check if it's a shower
                if (sg->flags_any(SegmentFlags::kShowerTrajectory) || sg->flags_any(SegmentFlags::kShowerTopology)) {
                    n_in_shower++;
                } else {
                    in_tracks.push_back(sg);
                }
            }
        }
        
        // If there are multiple incoming tracks (not all showers), reset their directions
        if (n_in > 1 && n_in != n_in_shower) {
            for (auto it1 = in_tracks.begin(); it1 != in_tracks.end(); it1++) {
                (*it1)->dirsign(0);
                (*it1)->dir_weak(true);
            }
        }
    }
}

void PatternAlgorithms::fix_maps_shower_in_track_out(Graph& graph, Facade::Cluster& cluster){
    // Iterate through all vertices in the graph
    auto [vbegin, vend] = boost::vertices(graph);
    for (auto vit = vbegin; vit != vend; ++vit) {
        VertexPtr vtx = graph[*vit].vertex;
        
        // Skip if vertex is null or doesn't belong to this cluster
        if (!vtx || !vtx->cluster() || vtx->cluster() != &cluster) continue;
        
        // Check how many segments are connected to this vertex
        if (!vtx->descriptor_valid()) continue;
        auto vd = vtx->get_descriptor();
        if (boost::degree(vd, graph) <= 1) continue;
        
        std::vector<SegmentPtr> in_showers;
        bool flag_turn_shower_dir = false;
        
        // Get vertex position
        WireCell::Point vtx_point = vtx->wcpt().point;
        
        // Iterate through all segments connected to this vertex
        auto edge_range = boost::out_edges(vd, graph);
        for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
            SegmentPtr sg = graph[*eit].segment;
            if (!sg) continue;
            
            // Determine if this vertex is at the front or back of the segment
            const auto& wcpts = sg->wcpts();
            if (wcpts.size() < 2) continue;
            
            auto front_pt = wcpts.front().point;
            auto back_pt = wcpts.back().point;
            
            double dis_front = ray_length(Ray{vtx_point, front_pt});
            double dis_back = ray_length(Ray{vtx_point, back_pt});
            
            bool flag_start = (dis_front < dis_back); // vertex is at the front of segment
            
            // Check if segment is a shower
            bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) || 
                           sg->flags_any(SegmentFlags::kShowerTopology);
            
            // Check if this is an "incoming" segment (pointing into vertex)
            if ((flag_start && sg->dirsign() == -1) || (!flag_start && sg->dirsign() == 1)) {
                if (is_shower) {
                    in_showers.push_back(sg);
                }
            }
            // Check if this is an "outgoing" segment (pointing away from vertex)
            else if ((flag_start && sg->dirsign() == 1) || (!flag_start && sg->dirsign() == -1)) {
                // If it's an outgoing non-shower track with strong direction
                if (!is_shower && !sg->dir_weak()) {
                    flag_turn_shower_dir = true;
                }
            }
        }
        
        // If there's a strong outgoing track and incoming showers, flip shower directions
        if (flag_turn_shower_dir) {
            for (auto it1 = in_showers.begin(); it1 != in_showers.end(); it1++) {
                (*it1)->dirsign((*it1)->dirsign() * (-1));
                (*it1)->dir_weak(true);
            }
        }
    }
}

void PatternAlgorithms::improve_maps_one_in(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_strong_check){
    bool flag_update = true;
    std::set<VertexPtr> used_vertices;
    std::set<SegmentPtr> used_segments;
    
    while(flag_update) {
        flag_update = false;
        
        // Iterate through all vertices in the graph
        auto [vbegin, vend] = boost::vertices(graph);
        for (auto vit = vbegin; vit != vend; ++vit) {
            VertexPtr vtx = graph[*vit].vertex;
            
            // Skip if vertex is null or doesn't belong to this cluster
            if (!vtx || !vtx->cluster() || vtx->cluster() != &cluster) continue;
            
            // Check how many segments are connected to this vertex
            if (!vtx->descriptor_valid()) continue;
            auto vd = vtx->get_descriptor();
            if (boost::degree(vd, graph) <= 1) continue;
            
            // Skip if already processed
            if (used_vertices.find(vtx) != used_vertices.end()) continue;
            
            int n_in = 0;
            std::map<SegmentPtr, bool> map_sg_dir; // segment -> flag_start
            
            // Get vertex position
            WireCell::Point vtx_point = vtx->wcpt().point;
            
            // Iterate through all segments connected to this vertex
            auto edge_range = boost::out_edges(vd, graph);
            for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
                SegmentPtr sg = graph[*eit].segment;
                if (!sg) continue;
                
                // Skip if segment already processed
                if (used_segments.find(sg) != used_segments.end()) continue;
                
                // Determine if this vertex is at the front or back of the segment
                const auto& wcpts = sg->wcpts();
                if (wcpts.size() < 2) continue;
                
                auto front_pt = wcpts.front().point;
                auto back_pt = wcpts.back().point;
                
                double dis_front = ray_length(Ray{vtx_point, front_pt});
                double dis_back = ray_length(Ray{vtx_point, back_pt});
                
                bool flag_start = (dis_front < dis_back); // vertex is at the front of segment
                
                // Check if this is an "incoming" segment (pointing into vertex)
                if ((flag_start && sg->dirsign() == -1) || (!flag_start && sg->dirsign() == 1)) {
                    if (flag_strong_check) {
                        // Only count if direction is strong
                        if (!sg->dir_weak()) n_in++;
                    } else {
                        n_in++;
                    }
                }
                
                // Collect segments with no or weak direction
                if (sg->dirsign() == 0 || sg->dir_weak()) {
                    map_sg_dir[sg] = flag_start;
                }
            }
            
            // If no segments to change direction, mark vertex as used
            if (map_sg_dir.size() == 0) {
                used_vertices.insert(vtx);
            }
            
            // If there are incoming segments, set all weak/no-direction segments to point out
            if (n_in > 0) {
                for (auto it1 = map_sg_dir.begin(); it1 != map_sg_dir.end(); it1++) {
                    SegmentPtr sg = it1->first;
                    bool flag_start = it1->second;
                    
                    // Set direction to point away from vertex
                    if (flag_start) {
                        sg->dirsign(1);  // at front, point forward
                    } else {
                        sg->dirsign(-1); // at back, point backward
                    }
                    
                    // Recalculate 4-momentum if particle info exists
                    if (sg->has_particle_info()) {
                        int pdg_code = sg->particle_info()->pdg();
                        auto four_momentum = segment_cal_4mom(sg, pdg_code, particle_data, recomb_model);
                        
                        // Update particle info with new 4-momentum
                        auto pinfo = std::make_shared<Aux::ParticleInfo>(
                            pdg_code,
                            particle_data->get_particle_mass(pdg_code),
                            particle_data->pdg_to_name(pdg_code),
                            four_momentum
                        );
                        sg->particle_info(pinfo);
                    }
                    
                    sg->dir_weak(true);
                    used_segments.insert(sg);
                    flag_update = true;
                }
                used_vertices.insert(vtx);
            }
        }
    }
}

void PatternAlgorithms::improve_maps_shower_in_track_out(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_strong_check){
    bool flag_update = true;
    std::set<VertexPtr> used_vertices;
    std::set<SegmentPtr> used_segments;
    
    while(flag_update) {
        flag_update = false;
        
        // Iterate through all vertices in the graph
        auto [vbegin, vend] = boost::vertices(graph);
        for (auto vit = vbegin; vit != vend; ++vit) {
            VertexPtr vtx = graph[*vit].vertex;
            
            // Skip if vertex is null or doesn't belong to this cluster
            if (!vtx || !vtx->cluster() || vtx->cluster() != &cluster) continue;
            
            // Check how many segments are connected to this vertex
            if (!vtx->descriptor_valid()) continue;
            auto vd = vtx->get_descriptor();
            if (boost::degree(vd, graph) <= 1) continue;
            
            // Skip if already processed
            if (used_vertices.find(vtx) != used_vertices.end()) continue;
            
            // int n_in = 0;
            int n_in_shower = 0;
            std::vector<SegmentPtr> out_tracks;
            std::map<SegmentPtr, bool> map_no_dir_segments; // segment -> flag_start
            
            // Get vertex position
            WireCell::Point vtx_point = vtx->wcpt().point;
            
            // Iterate through all segments connected to this vertex
            auto edge_range = boost::out_edges(vd, graph);
            for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
                SegmentPtr sg = graph[*eit].segment;
                if (!sg) continue;
                
                // Determine if this vertex is at the front or back of the segment
                const auto& wcpts = sg->wcpts();
                if (wcpts.size() < 2) continue;
                
                auto front_pt = wcpts.front().point;
                auto back_pt = wcpts.back().point;
                
                double dis_front = ray_length(Ray{vtx_point, front_pt});
                double dis_back = ray_length(Ray{vtx_point, back_pt});
                
                bool flag_start = (dis_front < dis_back); // vertex is at the front of segment
                
                bool is_shower = sg->flags_any(SegmentFlags::kShowerTrajectory) || 
                               sg->flags_any(SegmentFlags::kShowerTopology);
                
                // Check if this is an "incoming" segment (pointing into vertex)
                if ((flag_start && sg->dirsign() == -1) || (!flag_start && sg->dirsign() == 1)) {
                    // n_in++;
                    if (is_shower) {
                        n_in_shower++;
                    }
                }
                // Check if this is an "outgoing" segment (pointing away from vertex)
                else if ((flag_start && sg->dirsign() == 1) || (!flag_start && sg->dirsign() == -1)) {
                    if (!is_shower) {
                        // Check if it's weak or has no particle type
                        bool no_particle_type = !sg->has_particle_info() || sg->particle_info()->pdg() == 0;
                        if (sg->dir_weak() || (no_particle_type && !flag_strong_check)) {
                            out_tracks.push_back(sg);
                        }
                    }
                }
                // Segment with no direction
                else if (sg->dirsign() == 0) {
                    map_no_dir_segments[sg] = flag_start;
                }
            }
            
            // If there are incoming showers and outgoing tracks or no-direction segments
            if (n_in_shower > 0 && (out_tracks.size() > 0 || map_no_dir_segments.size() > 0)) {
                // Reclassify outgoing tracks as electrons
                for (auto it1 = out_tracks.begin(); it1 != out_tracks.end(); it1++) {
                    SegmentPtr sg1 = *it1;
                    
                    // Set as electron (PDG 11)
                    int pdg_code = 11;
                    auto four_momentum = WireCell::D4Vector<double>(0, 0, 0, 0);
                    
                    // Recalculate 4-momentum if segment has valid energy
                    if (sg1->has_particle_info() && sg1->particle_info()->energy() > 0) {
                        four_momentum = segment_cal_4mom(sg1, pdg_code, particle_data, recomb_model);
                    }
                    
                    auto pinfo = std::make_shared<Aux::ParticleInfo>(
                        pdg_code,
                        particle_data->get_particle_mass(pdg_code),
                        particle_data->pdg_to_name(pdg_code),
                        four_momentum
                    );
                    sg1->particle_info(pinfo);
                    sg1->dirsign(0);
                    
                    flag_update = true;
                }
                
                // Process no-direction segments
                for (auto it1 = map_no_dir_segments.begin(); it1 != map_no_dir_segments.end(); it1++) {
                    SegmentPtr sg1 = it1->first;
                    if (used_segments.find(sg1) != used_segments.end()) continue;
                    
                    // If it's not already a shower, set as electron
                    bool is_shower = sg1->flags_any(SegmentFlags::kShowerTrajectory) || 
                                   sg1->flags_any(SegmentFlags::kShowerTopology);
                    
                    if (!is_shower) {
                        int pdg_code = 11;
                        auto four_momentum = WireCell::D4Vector<double>(0, 0, 0, 0);
                        
                        // Recalculate 4-momentum if segment has valid energy
                        if (sg1->has_particle_info() && sg1->particle_info()->energy() > 0) {
                            four_momentum = segment_cal_4mom(sg1, pdg_code, particle_data, recomb_model);
                        }
                        
                        auto pinfo = std::make_shared<Aux::ParticleInfo>(
                            pdg_code,
                            particle_data->get_particle_mass(pdg_code),
                            particle_data->pdg_to_name(pdg_code),
                            four_momentum
                        );
                        sg1->particle_info(pinfo);
                    }
                    
                    sg1->dir_weak(true);
                    used_segments.insert(sg1);
                    flag_update = true;
                }
            }
            
            used_vertices.insert(vtx);
        }
    }
}

