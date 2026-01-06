#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRShowerFunctions.h"
#include <Eigen/Dense>

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

void PatternAlgorithms::update_shower_maps(std::set<ShowerPtr>& showers,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters){
    // Clear all maps
    map_vertex_to_shower.clear();
    map_vertex_in_shower.clear();
    map_segment_in_shower.clear();
    used_shower_clusters.clear();
    
    // Iterate through all showers
    for (auto shower : showers) {
        // Map start vertex to shower
        auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
        if (start_vtx) {
            map_vertex_to_shower[start_vtx].insert(shower);
        }
        
        // Fill maps using TrajectoryView - iterate through all vertices and segments in the shower
        TrajectoryView& traj = shower->fill_maps();
        
        // Fill map_vertex_in_shower with all vertices in this shower
        for (auto vdesc : traj.nodes()) {
            auto vtx = traj.view_graph()[vdesc].vertex;
            if (vtx) {
                map_vertex_in_shower[vtx] = shower;
            }
        }
        
        // Fill map_segment_in_shower with all segments in this shower
        for (auto edesc : traj.edges()) {
            auto seg = traj.view_graph()[edesc].segment;
            if (seg) {
                map_segment_in_shower[seg] = shower;
            }
        }
    }
    
    // Collect all cluster IDs from segments in the map
    for (auto it = map_segment_in_shower.begin(); it != map_segment_in_shower.end(); it++) {
        auto seg = it->first;
        if (seg && seg->cluster()) {
            used_shower_clusters.insert(seg->cluster());
        }
    }
}

void PatternAlgorithms::shower_clustering_with_nv_in_main_cluster(Graph& graph, VertexPtr main_vertex, std::set<ShowerPtr>& showers,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters, std::set<VertexPtr>& vertices_in_long_muon, std::set<SegmentPtr>& segments_in_long_muon){
    if (!main_vertex) return;
    
    // Build map_vertex_segments from graph
    std::map<VertexPtr, std::set<SegmentPtr>> map_vertex_segments;
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr seg = graph[*eit].segment;
        if (!seg) continue;
        
        auto source_vdesc = boost::source(*eit, graph);
        auto target_vdesc = boost::target(*eit, graph);
        VertexPtr v1 = graph[source_vdesc].vertex;
        VertexPtr v2 = graph[target_vdesc].vertex;
        
        if (v1) map_vertex_segments[v1].insert(seg);
        if (v2) map_vertex_segments[v2].insert(seg);
    }
    
    // Step 1: Collect segments from showers starting at main_vertex
    std::set<SegmentPtr> used_segments;
    for (auto shower : showers) {
        auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
        if (start_vtx == main_vertex) {
            // Collect all segments from this shower
            for (auto edesc : shower->edges()) {
                auto seg = shower->view_graph()[edesc].segment;
                if (seg) {
                    used_segments.insert(seg);
                }
            }
            
            // If there's only 1 segment and it's short, clear used_segments
            if (used_segments.size() == 1 && segment_track_length(*used_segments.begin()) < 8 * units::cm) {
                used_segments.clear();
            }
        }
    }
    
    // Step 2: If used_segments is empty, try to create new showers
    if (used_segments.empty()) {
        std::set<ShowerPtr> del_showers;
        std::map<ShowerPtr, SegmentPtr> map_shower_max_sg;
        ShowerPtr max_shower = nullptr;
        double max_length = 0;
        
        // Get segments connected to main_vertex
        auto& main_vtx_segments = map_vertex_segments[main_vertex];
        
        for (auto sg : main_vtx_segments) {
            // Skip if segment is in a shower
            if (map_segment_in_shower.find(sg) != map_segment_in_shower.end()) continue;
            
            // Skip segments in long muon
            if (segments_in_long_muon.find(sg) != segments_in_long_muon.end()) continue;
            
            // Get segment properties
            double medium_dQ_dx = segment_median_dQ_dx(sg);
            double medium_dQ_dx_1 = medium_dQ_dx / (43e3 / units::cm);
            
            // Get particle type if available
            int particle_type = 0;
            if (sg->has_particle_info() && sg->particle_info()) {
                particle_type = sg->particle_info()->pdg();
            }
            
            // Skip segments with certain particle types and high dQ/dx
            if ((particle_type == 11) ||
                (particle_type == 2212 && (medium_dQ_dx_1 > 1.45 || medium_dQ_dx_1 > 2.7)) ||
                (particle_type == 211 && medium_dQ_dx_1 > 2.0)) {
                continue;
            }
            
            // Create a new shower
            ShowerPtr shower = std::make_shared<Shower>(graph);
            shower->set_start_vertex(main_vertex, 1);
            shower->set_start_segment(sg);
            shower->complete_structure_with_start_segment(used_segments);
            
            // Analyze shower structure
            int n_tracks = 0;
            int n_showers = 0;
            double total_length = 0;
            double max_seg_length = 0;
            SegmentPtr max_sg = nullptr;
            bool flag_good_track = false;
            
            // Iterate through shower segments
            for (auto edesc : shower->edges()) {
                auto sg1 = shower->view_graph()[edesc].segment;
                if (!sg1) continue;
                
                double length = segment_track_length(sg1);
                double medium_dQ_dx_sg = segment_median_dQ_dx(sg1);
                double medium_dQ_dx_norm = medium_dQ_dx_sg / (43e3 / units::cm);
                
                // Check if segment is a shower
                bool is_shower = sg1->flags_any(SegmentFlags::kShowerTrajectory) || 
                                sg1->flags_any(SegmentFlags::kShowerTopology);
                if (is_shower) n_showers++;
                
                n_tracks++;
                total_length += length;
                
                if (max_seg_length < length) {
                    max_seg_length = length;
                    max_sg = sg1;
                }
                
                // Check for good track properties
                if (!sg1->dir_weak() && (length > 3.6 * units::cm || (length > 2.4 * units::cm && medium_dQ_dx_norm > 2.5))) {
                    // Find end vertex of this segment
                    auto seg_edesc = sg1->get_descriptor();
                    auto source_vdesc = boost::source(seg_edesc, graph);
                    auto target_vdesc = boost::target(seg_edesc, graph);
                    VertexPtr v1 = graph[source_vdesc].vertex;
                    VertexPtr v2 = graph[target_vdesc].vertex;
                    
                    // Determine end vertex based on dirsign
                    VertexPtr end_vertex = nullptr;
                    if (sg1->dirsign() == 1) {
                        // Check which vertex is at the end
                        if (!sg1->fits().empty() && v1 && v2) {
                            auto& fits = sg1->fits();
                            double dist1 = (fits.back().point - (v1->fit().valid() ? v1->fit().point : v1->wcpt().point)).magnitude();
                            double dist2 = (fits.back().point - (v2->fit().valid() ? v2->fit().point : v2->wcpt().point)).magnitude();
                            end_vertex = (dist1 < dist2) ? v1 : v2;
                        }
                    } else if (sg1->dirsign() == -1) {
                        if (!sg1->fits().empty() && v1 && v2) {
                            auto& fits = sg1->fits();
                            double dist1 = (fits.front().point - (v1->fit().valid() ? v1->fit().point : v1->wcpt().point)).magnitude();
                            double dist2 = (fits.front().point - (v2->fit().valid() ? v2->fit().point : v2->wcpt().point)).magnitude();
                            end_vertex = (dist1 < dist2) ? v1 : v2;
                        }
                    }
                    
                    if (end_vertex && map_vertex_segments.find(end_vertex) != map_vertex_segments.end()) {
                        if (map_vertex_segments[end_vertex].size() > 1) {
                            bool flag_non_ele = false;
                            for (auto sg2 : map_vertex_segments[end_vertex]) {
                                if (sg2 == sg1) continue;
                                bool is_shower2 = sg2->flags_any(SegmentFlags::kShowerTrajectory) || 
                                                 sg2->flags_any(SegmentFlags::kShowerTopology);
                                if (!is_shower2) flag_non_ele = true;
                            }
                            
                            if (!flag_non_ele && map_vertex_segments[end_vertex].size() <= 3) {
                                flag_good_track = true;
                            }
                        } else {
                            flag_good_track = true;
                        }
                    }
                }
            }

            (void) n_showers; // n_showers is not used further
            
            // Count vertex types in shower
            int n_multi_vtx = 0;
            int n_two_vtx = 0;
            std::map<VertexPtr, int> vtx_segment_count;
            
            for (auto edesc : shower->edges()) {
                auto seg = shower->view_graph()[edesc].segment;
                if (!seg) continue;
                
                auto seg_edesc = seg->get_descriptor();
                auto source_vdesc = boost::source(seg_edesc, graph);
                auto target_vdesc = boost::target(seg_edesc, graph);
                VertexPtr v1 = graph[source_vdesc].vertex;
                VertexPtr v2 = graph[target_vdesc].vertex;
                
                if (v1) vtx_segment_count[v1]++;
                if (v2) vtx_segment_count[v2]++;
            }
            
            for (auto& [vtx, count] : vtx_segment_count) {
                if (count == 2) n_two_vtx++;
                else if (count > 2) n_multi_vtx++;
            }
            
            // Apply selection criteria
            if (!flag_good_track && n_multi_vtx > 0 && max_seg_length < 65 * units::cm &&
                ((total_length < n_tracks * 27 * units::cm && total_length < 85 * units::cm) ||
                 (total_length < n_tracks * 18 * units::cm && total_length < 95 * units::cm)) &&
                n_two_vtx < 3) {
                
                map_shower_max_sg[shower] = max_sg;
                if (shower->get_total_length() > max_length) {
                    max_shower = shower;
                    max_length = shower->get_total_length();
                }
            }
        }
        
        // Process selected showers
        for (auto& [shower, max_sg] : map_shower_max_sg) {
            if (shower == max_shower) {
                // Convert to EM shower (particle type 11 = electron)
                if (shower->start_segment() && shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info()) {
                    shower->start_segment()->particle_info()->set_pdg(11);
                }
                
                // Set avoid muon check flag on max segment
                if (max_sg) {
                    max_sg->set_flags(SegmentFlags::kAvoidMuonCheck);
                }
                
                // Add shower to collection
                showers.insert(shower);
                
                // Check for showers to delete that conflict with new shower
                std::set<VertexPtr> shower_vertices;
                for (auto vdesc : shower->nodes()) {
                    auto vtx = shower->view_graph()[vdesc].vertex;
                    if (vtx) shower_vertices.insert(vtx);
                }
                
                for (auto shower1 : showers) {
                    if (shower == shower1) continue;
                    auto [start_vtx1, conn_type1] = shower1->get_start_vertex_and_type();
                    if (conn_type1 == 1 && start_vtx1 && shower_vertices.find(start_vtx1) != shower_vertices.end()) {
                        del_showers.insert(shower1);
                    }
                }
            }
        }
        
        // Delete conflicting showers
        for (auto shower1 : del_showers) {
            auto [start_vtx1, conn_type1] = shower1->get_start_vertex_and_type();
            if (start_vtx1 != main_vertex) {
                showers.erase(shower1);
            }
        }
        
        // Update shower maps
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters);
    }
}

void PatternAlgorithms::shower_clustering_connecting_to_main_vertex(Graph& graph, VertexPtr main_vertex, std::set<ShowerPtr>& showers,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters){
    if (!main_vertex) return;
    
    // Build map_vertex_segments from graph
    std::map<VertexPtr, std::set<SegmentPtr>> map_vertex_segments;
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr seg = graph[*eit].segment;
        if (!seg) continue;
        
        auto source_vdesc = boost::source(*eit, graph);
        auto target_vdesc = boost::target(*eit, graph);
        VertexPtr v1 = graph[source_vdesc].vertex;
        VertexPtr v2 = graph[target_vdesc].vertex;
        
        if (v1) map_vertex_segments[v1].insert(seg);
        if (v2) map_vertex_segments[v2].insert(seg);
    }
    
    // Step 1: Collect segments from showers starting at main_vertex
    std::set<SegmentPtr> used_segments;
    for (auto shower : showers) {
        auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
        if (start_vtx == main_vertex) {
            // Collect all segments from this shower
            for (auto edesc : shower->edges()) {
                auto seg = shower->view_graph()[edesc].segment;
                if (seg) {
                    used_segments.insert(seg);
                }
            }
            
            // If there's only 1 segment and it's short, clear used_segments
            if (used_segments.size() == 1 && segment_track_length(*used_segments.begin()) < 8 * units::cm) {
                used_segments.clear();
            }
        }
    }
    
    // Step 2: If used_segments is empty, try to create new showers
    if (used_segments.empty()) {
        std::set<ShowerPtr> del_showers;
        std::map<ShowerPtr, SegmentPtr> map_shower_max_sg;
        ShowerPtr max_shower = nullptr;
        double max_length = 0;
        
        // Get segments connected to main_vertex
        auto& main_vtx_segments = map_vertex_segments[main_vertex];
        
        for (auto sg : main_vtx_segments) {
            // Calculate number of daughter showers for this segment
            auto pair_result = calculate_num_daughter_showers(graph, main_vertex, sg, true);
            
            // Skip if segment is in a shower
            if (map_segment_in_shower.find(sg) != map_segment_in_shower.end()) continue;
            
            // Get segment properties
            double medium_dQ_dx = segment_median_dQ_dx(sg);
            double medium_dQ_dx_1 = medium_dQ_dx / (43e3 / units::cm);
            
            // Get particle type if available
            int particle_type = 0;
            if (sg->has_particle_info() && sg->particle_info()) {
                particle_type = sg->particle_info()->pdg();
            }
            
            // Skip segments with certain particle types and high dQ/dx
            if ((particle_type == 11) ||
                (particle_type == 2212 && ((medium_dQ_dx_1 > 1.45 && pair_result.first <= 3) || medium_dQ_dx_1 > 2.7)) ||
                (particle_type == 211 && medium_dQ_dx_1 > 2.0)) {
                continue;
            }
            
            // Create a new shower
            ShowerPtr shower = std::make_shared<Shower>(graph);
            shower->set_start_vertex(main_vertex, 1);
            shower->set_start_segment(sg);
            shower->complete_structure_with_start_segment(used_segments);
            
            // Analyze shower structure
            int n_tracks = 0;
            int n_showers = 0;
            double total_length = 0;
            double max_seg_length = 0;
            SegmentPtr max_sg = nullptr;
            bool flag_good_track = false;
            
            // Iterate through shower segments
            for (auto edesc : shower->edges()) {
                auto sg1 = shower->view_graph()[edesc].segment;
                if (!sg1) continue;
                
                double length = segment_track_length(sg1);
                double medium_dQ_dx_sg = segment_median_dQ_dx(sg1);
                double medium_dQ_dx_norm = medium_dQ_dx_sg / (43e3 / units::cm);
                
                // Check if segment is a shower
                bool is_shower = sg1->flags_any(SegmentFlags::kShowerTrajectory) || 
                                sg1->flags_any(SegmentFlags::kShowerTopology);
                if (is_shower) n_showers++;
                
                n_tracks++;
                total_length += length;
                
                if (max_seg_length < length) {
                    max_seg_length = length;
                    max_sg = sg1;
                }
                
                // Check for good track properties
                if (!sg1->dir_weak() && (length > 3.6 * units::cm || (length > 2.4 * units::cm && medium_dQ_dx_norm > 2.5))) {
                    // Find end vertex of this segment
                    auto seg_edesc = sg1->get_descriptor();
                    auto source_vdesc = boost::source(seg_edesc, graph);
                    auto target_vdesc = boost::target(seg_edesc, graph);
                    VertexPtr v1 = graph[source_vdesc].vertex;
                    VertexPtr v2 = graph[target_vdesc].vertex;
                    
                    // Determine end vertex based on dirsign
                    // For dirsign == 1, direction is from front to back of fits
                    // For dirsign == -1, direction is from back to front of fits
                    VertexPtr end_vertex = nullptr;
                    if (sg1->dirsign() == 1) {
                        // End is at back of fits
                        if (!sg1->fits().empty() && v1 && v2) {
                            auto& fits = sg1->fits();
                            double dist1 = (fits.back().point - (v1->fit().valid() ? v1->fit().point : v1->wcpt().point)).magnitude();
                            double dist2 = (fits.back().point - (v2->fit().valid() ? v2->fit().point : v2->wcpt().point)).magnitude();
                            end_vertex = (dist1 < dist2) ? v1 : v2;
                        }
                    } else if (sg1->dirsign() == -1) {
                        // End is at front of fits (reversed direction)
                        if (!sg1->fits().empty() && v1 && v2) {
                            auto& fits = sg1->fits();
                            double dist1 = (fits.front().point - (v1->fit().valid() ? v1->fit().point : v1->wcpt().point)).magnitude();
                            double dist2 = (fits.front().point - (v2->fit().valid() ? v2->fit().point : v2->wcpt().point)).magnitude();
                            end_vertex = (dist1 < dist2) ? v1 : v2;
                        }
                    }
                    
                    if (end_vertex && map_vertex_segments.find(end_vertex) != map_vertex_segments.end()) {
                        if (map_vertex_segments[end_vertex].size() > 1) {
                            bool flag_non_ele = false;
                            for (auto sg2 : map_vertex_segments[end_vertex]) {
                                if (sg2 == sg1) continue;
                                bool is_shower2 = sg2->flags_any(SegmentFlags::kShowerTrajectory) || 
                                                 sg2->flags_any(SegmentFlags::kShowerTopology);
                                if (!is_shower2) flag_non_ele = true;
                            }
                            
                            if (!flag_non_ele && map_vertex_segments[end_vertex].size() <= 3) {
                                flag_good_track = true;
                            }
                        } else {
                            flag_good_track = true;
                        }
                    }
                }
            }
            
            (void) n_showers; // n_showers is not used further
            
            // Count vertex types in shower
            int n_multi_vtx = 0;
            int n_two_vtx = 0;
            std::map<VertexPtr, int> vtx_segment_count;
            
            for (auto edesc : shower->edges()) {
                auto seg = shower->view_graph()[edesc].segment;
                if (!seg) continue;
                
                auto seg_edesc = seg->get_descriptor();
                auto source_vdesc = boost::source(seg_edesc, graph);
                auto target_vdesc = boost::target(seg_edesc, graph);
                VertexPtr v1 = graph[source_vdesc].vertex;
                VertexPtr v2 = graph[target_vdesc].vertex;
                
                if (v1) vtx_segment_count[v1]++;
                if (v2) vtx_segment_count[v2]++;
            }
            
            for (auto& [vtx, count] : vtx_segment_count) {
                if (count == 2) n_two_vtx++;
                else if (count > 2) n_multi_vtx++;
            }
            
            // Apply selection criteria
            if (!flag_good_track && n_multi_vtx > 0 && max_seg_length < 65 * units::cm &&
                ((total_length < n_tracks * 27 * units::cm && total_length < 85 * units::cm) ||
                 (total_length < n_tracks * 18 * units::cm && total_length < 95 * units::cm)) &&
                n_two_vtx < 3) {
                
                map_shower_max_sg[shower] = max_sg;
                if (shower->get_total_length() > max_length) {
                    max_shower = shower;
                    max_length = shower->get_total_length();
                }
            }
        }
        
        // Process selected showers
        for (auto& [shower, max_sg] : map_shower_max_sg) {
            if (shower == max_shower) {
                // Convert to EM shower (particle type 11 = electron)
                std::cout << "Convert EM shower " << shower->start_segment()->id() << std::endl;
                if (shower->start_segment() && shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info()) {
                    shower->start_segment()->particle_info()->set_pdg(11);
                }
                
                // Set avoid muon check flag on max segment
                if (max_sg) {
                    max_sg->set_flags(SegmentFlags::kAvoidMuonCheck);
                }
                
                // Add shower to collection
                showers.insert(shower);
                
                // Check for showers to delete that conflict with new shower
                // Collect all vertices in the new shower
                std::set<VertexPtr> shower_vertices;
                for (auto vdesc : shower->nodes()) {
                    auto vtx = shower->view_graph()[vdesc].vertex;
                    if (vtx) shower_vertices.insert(vtx);
                }
                
                for (auto shower1 : showers) {
                    if (shower == shower1) continue;
                    auto [start_vtx1, conn_type1] = shower1->get_start_vertex_and_type();
                    if (conn_type1 == 1 && start_vtx1 && shower_vertices.find(start_vtx1) != shower_vertices.end()) {
                        del_showers.insert(shower1);
                    }
                }
            }
        }
        
        // Delete conflicting showers
        for (auto shower1 : del_showers) {
            auto [start_vtx1, conn_type1] = shower1->get_start_vertex_and_type();
            if (start_vtx1 != main_vertex) {
                showers.erase(shower1);
            }
        }
        
        // Update shower maps
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters);
    }
}
