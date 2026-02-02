#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRShowerFunctions.h"
#include <Eigen/Dense>

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

namespace {
    struct cluster_point_info {
        Facade::Cluster* cluster;
        double min_angle;
        double min_dis;
        VertexPtr min_vertex;
        WireCell::Point min_point;
    };
    
    bool sortbydis(const cluster_point_info &a, const cluster_point_info &b) {
        return (a.min_dis < b.min_dis);
    }
}

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

void PatternAlgorithms::shower_clustering_with_nv_from_main_cluster(Graph& graph, VertexPtr main_vertex, Facade::Cluster* main_cluster, std::set<ShowerPtr>& showers,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters){
    if (!main_vertex || !main_cluster) return;
    
    // Build map_segment_vertices from graph
    std::map<SegmentPtr, std::set<VertexPtr>> map_segment_vertices;
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr seg = graph[*eit].segment;
        if (!seg) continue;
        
        auto source_vdesc = boost::source(*eit, graph);
        auto target_vdesc = boost::target(*eit, graph);
        VertexPtr v1 = graph[source_vdesc].vertex;
        VertexPtr v2 = graph[target_vdesc].vertex;
        
        if (v1) map_segment_vertices[seg].insert(v1);
        if (v2) map_segment_vertices[seg].insert(v2);
    }
    
    std::map<ShowerPtr, WireCell::Vector> map_shower_dir;
    std::map<ShowerPtr, double> map_shower_angle_offset;
    WireCell::Vector drift_dir(1, 0, 0);
    SegmentPtr max_length_segment = nullptr;
    
    // Step 1: Find the maximum length segment in main_cluster
    {
        double max_length = 0;
        for (auto& [seg, vertices] : map_segment_vertices) {
            if (seg->cluster() != main_cluster) continue;
            double length = segment_track_length(seg);
            if (length > max_length && length > 6 * units::cm) {
                max_length = length;
                max_length_segment = seg;
            }
        }
    }
    
    // Step 2: Build map_shower_dir for showers in main_cluster
    for (auto& [seg, vertices] : map_segment_vertices) {
        if (seg->cluster() != main_cluster) continue;
        if (map_segment_in_shower.find(seg) == map_segment_in_shower.end()) continue;
        
        ShowerPtr shower = map_segment_in_shower[seg];
        
        // Skip long muons
        int particle_type = 0;
        if (shower->start_segment() && shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info()) {
            particle_type = shower->start_segment()->particle_info()->pdg();
        }
        if (std::abs(particle_type) == 13) continue;
        
        double total_length = shower->get_total_length();
        
        if (seg == shower->start_segment()) {
            auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
            WireCell::Point start_point = start_vtx ? (start_vtx->fit().valid() ? start_vtx->fit().point : start_vtx->wcpt().point) : WireCell::Point(0, 0, 0);
            
            bool is_shower_topology = seg->flags_any(SegmentFlags::kShowerTopology);
            double medium_dQ_dx = segment_median_dQ_dx(seg);
            double seg_length = segment_track_length(seg);
            
            if (is_shower_topology || shower->get_num_segments() > 2 || 
                (medium_dQ_dx > 43e3 / units::cm * 1.5 && seg_length > 0)) {
                if (seg_length > 10 * units::cm) {
                    WireCell::Vector dir_shower = segment_cal_dir_3vector(seg, start_point, 15 * units::cm);
                    map_shower_dir[shower] = dir_shower;
                } else {
                    WireCell::Vector dir_shower = shower_cal_dir_3vector(*shower, start_point, 15 * units::cm);
                    map_shower_dir[shower] = dir_shower;
                }
            } else if (shower->get_num_segments() <= 2) {
                if (total_length > 30 * units::cm) {
                    if (seg_length > 10 * units::cm) {
                        WireCell::Vector dir_shower = segment_cal_dir_3vector(seg, start_point, 15 * units::cm);
                        map_shower_dir[shower] = dir_shower;
                    } else {
                        WireCell::Vector dir_shower = shower_cal_dir_3vector(*shower, start_point, 15 * units::cm);
                        map_shower_dir[shower] = dir_shower;
                    }
                }
            }
            
            // Very large shower
            if (total_length > 100 * units::cm) {
                WireCell::Vector dir_shower = shower_cal_dir_3vector(*shower, start_point, 60 * units::cm);
                map_shower_dir[shower] = dir_shower;
            }
            
            // Max length segment or segment connected to main_vertex
            if (seg == max_length_segment && map_shower_dir.find(shower) == map_shower_dir.end()) {
                WireCell::Vector dir_shower = shower_cal_dir_3vector(*shower, start_point, 15 * units::cm);
                map_shower_dir[shower] = dir_shower;
            } else if (map_shower_dir.find(shower) == map_shower_dir.end() && 
                      map_segment_vertices[seg].find(main_vertex) != map_segment_vertices[seg].end()) {
                if (seg_length > 5 * units::cm) {
                    WireCell::Vector dir_shower = shower_cal_dir_3vector(*shower, start_point, 15 * units::cm);
                    map_shower_dir[shower] = dir_shower;
                }
            }
            
            // Check if parallel to drift direction
            if (map_shower_dir.find(shower) != map_shower_dir.end()) {
                map_shower_angle_offset[shower] = 0;
                double angle_to_drift = std::abs(map_shower_dir[shower].dot(drift_dir) / (map_shower_dir[shower].magnitude() * drift_dir.magnitude()));
                angle_to_drift = std::acos(std::clamp(angle_to_drift, -1.0, 1.0)) / M_PI * 180.0;
                if (std::abs(angle_to_drift - 90) < 5) {
                    map_shower_dir[shower] = shower_cal_dir_3vector(*shower, start_point, 50 * units::cm);
                    map_shower_angle_offset[shower] = 5;
                }
            }
        }
    }
    
    // Step 3: If no shower directions found, try to add segments based on closest distance
    if (map_shower_dir.empty()) {
        std::map<ShowerPtr, double> map_shower_length;
        for (auto shower : showers) {
            map_shower_length[shower] = shower->get_total_length();
        }
        
        bool flag_continue = true;
        while (flag_continue) {
            flag_continue = false;
            for (auto& [seg1, vertices] : map_segment_vertices) {
                if (seg1->cluster() == main_cluster) continue;
                if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) continue;
                
                double min_dis = 1e9;
                ShowerPtr min_shower = nullptr;
                
                for (auto shower : showers) {
                    int particle_type = 0;
                    if (shower->start_segment() && shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info()) {
                        particle_type = shower->start_segment()->particle_info()->pdg();
                    }
                    if (particle_type == 13) continue;
                    
                    if (segment_track_length(seg1) > 0.75 * map_shower_length[shower]) continue;
                    
                    double dis = shower_get_closest_dis(*shower, seg1);
                    if (dis < min_dis) {
                        min_dis = dis;
                        min_shower = shower;
                    }
                }
                
                if (min_shower && min_dis < 3.5 * units::cm) {
                    min_shower->add_segment(seg1);
                    map_shower_length[min_shower] = min_shower->get_total_length();
                    flag_continue = true;
                }
            }
            update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters);
        }
    }
    
    if (map_shower_dir.empty()) return;
    
    // Step 4: Examine other segments and add to showers based on angle and distance
    for (auto& [seg1, vertices] : map_segment_vertices) {
        if (seg1->cluster() == main_cluster) continue;
        if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) continue;
        
        double min_dis = 1e9;
        ShowerPtr min_shower = nullptr;
        
        for (auto& [shower, dir] : map_shower_dir) {
            auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
            if (!start_vtx) continue;
            
            WireCell::Point start_point = start_vtx->fit().valid() ? start_vtx->fit().point : start_vtx->wcpt().point;
            
            // Get closest point on segment to shower start vertex
            auto [dist, closest_pt] = segment_get_closest_point(seg1, start_point);
            
            // Vector from shower start to closest point
            WireCell::Vector v1(closest_pt.x() - start_point.x(), 
                               closest_pt.y() - start_point.y(), 
                               closest_pt.z() - start_point.z());
            
            double angle = std::acos(std::clamp(dir.dot(v1) / (dir.magnitude() * v1.magnitude()), -1.0, 1.0));
            angle = angle / M_PI * 180.0;
            
            double angle_offset = 0;
            if (map_shower_angle_offset.find(shower) != map_shower_angle_offset.end()) {
                angle_offset = map_shower_angle_offset[shower];
            }
            
            // Check angle and distance criteria
            if ((angle < 25.0 + angle_offset && dist < 80 * units::cm) ||
                (angle < 12.5 + angle_offset * 8 / 5 && dist < 130 * units::cm) ||
                (angle < 5 + angle_offset * 2 && dist < 200 * units::cm)) {
                
                double dis = std::pow(dist * std::cos(angle * M_PI / 180.0), 2) / std::pow(40 * units::cm, 2) + 
                            std::pow(dist * std::sin(angle * M_PI / 180.0), 2) / std::pow(5 * units::cm, 2);
                
                if (dis < min_dis) {
                    min_dis = dis;
                    min_shower = shower;
                }
            }
        }
        
        if (min_shower) {
            min_shower->add_segment(seg1);
        }
    }
    
    update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters);
}

void PatternAlgorithms::shower_clustering_with_nv_from_vertices(Graph& graph, VertexPtr main_vertex, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, std::set<ShowerPtr>& showers,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters, std::set<VertexPtr>& vertices_in_long_muon, std::set<SegmentPtr>& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
    if (!main_vertex || !main_cluster) return;
    
    // Build map_cluster_segments and map_segment_cluster
    std::map<Facade::Cluster*, std::vector<SegmentPtr>> map_cluster_segments;
    std::map<SegmentPtr, Facade::Cluster*> map_segment_cluster;
    std::map<SegmentPtr, std::set<VertexPtr>> map_segment_vertices;
    
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr seg = graph[*eit].segment;
        if (!seg || !seg->cluster()) continue;
        
        map_cluster_segments[seg->cluster()].push_back(seg);
        map_segment_cluster[seg] = seg->cluster();
        
        auto source_vdesc = boost::source(*eit, graph);
        auto target_vdesc = boost::target(*eit, graph);
        VertexPtr v1 = graph[source_vdesc].vertex;
        VertexPtr v2 = graph[target_vdesc].vertex;
        
        if (v1) map_segment_vertices[seg].insert(v1);
        if (v2) map_segment_vertices[seg].insert(v2);
    }
    
    // Step 1: Build map_cluster_center_point
    std::map<Facade::Cluster*, std::pair<WireCell::Point, double>> map_cluster_center_point;
    
    for (auto cluster : other_clusters) {
        auto it1 = map_cluster_segments.find(cluster);
        if (it1 == map_cluster_segments.end()) continue;
        
        double acc_length = 0;
        double acc_length1 = 0;
        WireCell::Point p(0, 0, 0);
        int np = 0;
        
        for (auto seg : it1->second) {
            if (map_segment_in_shower.find(seg) != map_segment_in_shower.end()) continue;
            
            bool is_shower = seg->flags_any(SegmentFlags::kShowerTrajectory) || seg->flags_any(SegmentFlags::kShowerTopology);
            int particle_type = 0;
            if (seg->has_particle_info() && seg->particle_info()) {
                particle_type = seg->particle_info()->pdg();
            }
            
            if (is_shower || particle_type == 0 || 
                ((std::abs(particle_type) == 13 || std::abs(particle_type) == 211) && seg->dir_weak())) {
                double length = segment_track_length(seg);
                acc_length += length;
                
                const auto& fits = seg->fits();
                for (const auto& fit : fits) {
                    p.set(p.x() + fit.point.x(), p.y() + fit.point.y(), p.z() + fit.point.z());
                    np++;
                }
            }
            
            if (particle_type != 11 && !is_shower) {
                acc_length1 += segment_track_length(seg);
            }
        }
        
        if ((acc_length > 1.0 * units::cm && acc_length >= acc_length1) || acc_length > 10 * units::cm) {
            if (np > 0) {
                p.set(p.x() / np, p.y() / np, p.z() / np);
            }
            map_cluster_center_point[cluster] = std::make_pair(p, acc_length);
        }
    }
    
    // Step 2: List main cluster vertices
    std::vector<VertexPtr> main_cluster_vertices;
    auto [vbegin, vend] = boost::vertices(graph);
    for (auto vit = vbegin; vit != vend; ++vit) {
        VertexPtr vtx = graph[*vit].vertex;
        if (!vtx || !vtx->cluster() || vtx->cluster() != main_cluster) continue;
        
        if (vtx != main_vertex) {
            if (vertices_in_long_muon.find(vtx) != vertices_in_long_muon.end()) continue;
            if (map_vertex_in_shower.find(vtx) != map_vertex_in_shower.end()) continue;
        }
        main_cluster_vertices.push_back(vtx);
    }
    
    // Step 3: Analyze each cluster against main cluster vertices
    std::map<Facade::Cluster*, cluster_point_info> map_cluster_pi;
    
    // Get wpid_params for point cloud creation
    auto* grouping = main_cluster->grouping();
    if (!grouping) return;
    const auto& wpids = grouping->wpids();
    std::map<WirePlaneId, std::tuple<Facade::geo_point_t, double, double, double>> wpid_params;
    std::map<WirePlaneId, std::pair<Facade::geo_point_t, double>> wpid_U_dir, wpid_V_dir, wpid_W_dir;
    std::set<int> apas;
    Facade::compute_wireplane_params(wpids, dv, wpid_params, wpid_U_dir, wpid_V_dir, wpid_W_dir, apas);
    
    for (auto& [cluster, center_pair] : map_cluster_center_point) {
        WireCell::Point center_p = center_pair.first;
        
        // Create point cloud from cluster segments
        std::vector<std::pair<Facade::geo_point_t, WirePlaneId>> point_plane_pairs;
        for (auto seg : map_cluster_segments[cluster]) {
            for (const auto& fit : seg->fits()) {
                WirePlaneId wpid = dv->contained_by(fit.point);
                point_plane_pairs.emplace_back(fit.point, wpid);
            }
        }
        auto dpc_points = Facade::make_points_direct(cluster, dv, wpid_params, point_plane_pairs, false);
        auto pcloud = std::make_shared<Facade::DynamicPointCloud>(wpid_params);
        pcloud->add_points(dpc_points);
        
        cluster_point_info min_pi;
        min_pi.cluster = cluster;
        min_pi.min_angle = 90;
        min_pi.min_dis = 1e9;
        min_pi.min_vertex = nullptr;
        
        cluster_point_info main_pi;
        main_pi.cluster = cluster;
        main_pi.min_vertex = main_vertex;
        
        for (auto vtx : main_cluster_vertices) {
            WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
            
            // Get closest point using KD-tree
            auto& kd3d = pcloud->kd3d();
            std::vector<double> query = {vtx_pt.x(), vtx_pt.y(), vtx_pt.z()};
            auto results = kd3d.knn(1, query);
            
            if (results.empty()) continue;
            
            const size_t idx = results[0].first;
            const double dis = std::sqrt(results[0].second);  // KD-tree returns squared distance
            const auto& closest_pt_data = pcloud->get_points()[idx];
            WireCell::Point closest_pt(closest_pt_data.x, closest_pt_data.y, closest_pt_data.z);
            
            WireCell::Vector v1(closest_pt.x() - vtx_pt.x(),
                               closest_pt.y() - vtx_pt.y(),
                               closest_pt.z() - vtx_pt.z());
            WireCell::Vector v2(center_p.x() - closest_pt.x(),
                               center_p.y() - closest_pt.y(),
                               center_p.z() - closest_pt.z());
            
            WireCell::Point near_center = pcloud->get_center_point_radius(closest_pt, 2 * units::cm);
            WireCell::Vector v3(near_center.x() - closest_pt.x(),
                               near_center.y() - closest_pt.y(),
                               near_center.z() - closest_pt.z());
            
            double angle = std::acos(std::clamp(v1.dot(v2) / (v1.magnitude() * v2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
            double angle3 = std::acos(std::clamp(v1.dot(v3) / (v1.magnitude() * v3.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
            
            if (angle < 30 || (dis < 5 * units::cm && angle < 45)) {
                angle = std::min(angle, angle3);
            }
            
            if (angle < 7.5) {
                if (dis * std::sin(angle / 180.0 * M_PI) < min_pi.min_dis * std::sin(min_pi.min_angle / 180.0 * M_PI) && angle < 90) {
                    min_pi.min_angle = angle;
                    min_pi.min_dis = dis;
                    min_pi.min_vertex = vtx;
                    min_pi.min_point = closest_pt;
                }
            } else {
                if (angle < min_pi.min_angle) {
                    min_pi.min_angle = angle;
                    min_pi.min_dis = dis;
                    min_pi.min_vertex = vtx;
                    min_pi.min_point = closest_pt;
                }
            }
            
            if (vtx == main_vertex) {
                main_pi.min_angle = angle;
                main_pi.min_dis = dis;
                main_pi.min_point = closest_pt;
            }
        }
        
        if (!min_pi.min_vertex) {
            min_pi.min_angle = main_pi.min_angle;
            min_pi.min_vertex = main_vertex;
            min_pi.min_point = main_pi.min_point;
            min_pi.min_dis = main_pi.min_dis;
        }
        
        WireCell::Point main_vtx_pt = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
        WireCell::Point min_vtx_pt = min_pi.min_vertex->fit().valid() ? min_pi.min_vertex->fit().point : min_pi.min_vertex->wcpt().point;
        double vtx_dis = (main_vtx_pt - min_vtx_pt).magnitude();
        
        if (main_pi.min_angle < min_pi.min_angle + 3 && main_pi.min_dis < min_pi.min_dis * 1.2 &&
            (min_pi.min_angle > 0.9 * main_pi.min_angle || vtx_dis < 1.5 * units::cm)) {
            map_cluster_pi[cluster] = main_pi;
        } else {
            map_cluster_pi[cluster] = min_pi;
        }
    }
    
    // Step 4: Sort by distance
    std::vector<cluster_point_info> vec_pi;
    for (auto& [cluster, pi] : map_cluster_pi) {
        vec_pi.push_back(pi);
    }
    std::sort(vec_pi.begin(), vec_pi.end(), sortbydis);
    
    std::map<Facade::Cluster*, VertexPtr> map_cluster_associated_vertex;
    for (const auto& pi : vec_pi) {
        if (pi.min_angle < 10) {
            map_cluster_associated_vertex[pi.cluster] = pi.min_vertex;
        }
    }
    
    // Step 5: Process each cluster
    for (const auto& pi : vec_pi) {
        Facade::Cluster* cluster = pi.cluster;
        VertexPtr vertex = pi.min_vertex;
        WireCell::Point point = pi.min_point;
        SegmentPtr sg1 = nullptr;
        double angle = pi.min_angle;
        
        if (angle > 50 && pi.min_dis > 6 * units::cm) continue;
        if (angle > 60) continue;
        
        // Find segment at the point
        for (auto seg : map_cluster_segments[cluster]) {
            if (map_segment_in_shower.find(seg) != map_segment_in_shower.end()) continue;
            auto [dis, closest_pt] = segment_get_closest_point(seg, point);
            if (dis < 0.01 * units::cm) {
                sg1 = seg;
                break;
            }
        }
        
        if (!sg1) continue;
        
        // Create new shower
        ShowerPtr shower = std::make_shared<Shower>(graph);
        shower->set_start_vertex(vertex, 2);
        showers.insert(shower);
        
        // Check if point is at segment endpoint
        const auto& fits = sg1->fits();
        if (!fits.empty()) {
            double dist_front = (fits.front().point - point).magnitude();
            double dist_back = (fits.back().point - point).magnitude();
            
            if (dist_front < 0.01 * units::cm || dist_back < 0.01 * units::cm) {
                shower->set_start_segment(sg1, true);
            } else {
                // Break segment at point
                auto [success, seg_pair, new_vtx] = break_segment(graph, sg1, point, particle_data, recomb_model, dv);
                
                if (!success || !new_vtx) {
                    shower->set_start_segment(sg1);
                } else {
                    // Determine which segment to use based on direction to vertex
                    WireCell::Point vtx_pt = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
                    WireCell::Vector v3(point.x() - vtx_pt.x(), point.y() - vtx_pt.y(), point.z() - vtx_pt.z());
                    
                    WireCell::Vector v1 = segment_cal_dir_3vector(seg_pair.first, point, 5 * units::cm);
                    WireCell::Vector v2 = segment_cal_dir_3vector(seg_pair.second, point, 5 * units::cm);
                    
                    double angle1 = std::acos(std::clamp(v1.dot(v3) / (v1.magnitude() * v3.magnitude()), -1.0, 1.0));
                    double angle2 = std::acos(std::clamp(v2.dot(v3) / (v2.magnitude() * v3.magnitude()), -1.0, 1.0));
                    
                    if (angle1 < angle2) {
                        shower->set_start_segment(seg_pair.first);
                    } else {
                        shower->set_start_segment(seg_pair.second);
                    }
                }
            }
        } else {
            shower->set_start_segment(sg1);
        }
        
        // Set direction based on vertex proximity
        auto start_seg = shower->start_segment();
        if (start_seg) {
            const auto& seg_fits = start_seg->fits();
            if (!seg_fits.empty()) {
                WireCell::Point vtx_pt = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
                double dis1 = (vtx_pt - seg_fits.front().point).magnitude();
                double dis2 = (vtx_pt - seg_fits.back().point).magnitude();
                
                if (dis1 < dis2) {
                    start_seg->dirsign(1);
                } else {
                    start_seg->dirsign(-1);
                }
            }
            
            // Set particle type to electron if needed
            int pdg = 0;
            if (start_seg->has_particle_info() && start_seg->particle_info()) {
                pdg = start_seg->particle_info()->pdg();
            }
            if (pdg == 0 || std::abs(pdg) == 13) {
                auto four_momentum = segment_cal_4mom(start_seg, 11, particle_data, recomb_model);
                
                // Create ParticleInfo for electron
                auto pinfo = std::make_shared<Aux::ParticleInfo>(
                    11,                                          // electron PDG
                    particle_data->get_particle_mass(11),       // electron mass
                    particle_data->pdg_to_name(11),             // "electron"
                    four_momentum                                // 4-momentum
                );
                
                // Store particle info in start_segment
                start_seg->particle_info(pinfo);
            }
        }
        
        // Complete shower structure
        std::set<SegmentPtr> used_segments;
        shower->complete_structure_with_start_segment(used_segments);
        
        // Calculate shower direction
        auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
        WireCell::Point start_pt = start_vtx ? (start_vtx->fit().valid() ? start_vtx->fit().point : start_vtx->wcpt().point) : point;
        
        WireCell::Vector dir_shower = segment_cal_dir_3vector(shower->start_segment(), point, 15 * units::cm);
        WireCell::Vector dir_main(point.x() - start_pt.x(), point.y() - start_pt.y(), point.z() - start_pt.z());
        
        if (std::acos(std::clamp(dir_shower.dot(dir_main) / (dir_shower.magnitude() * dir_main.magnitude()), -1.0, 1.0)) / M_PI * 180.0 > 30) {
            auto [_, test_p] = shower_get_closest_point(*shower, start_pt);
            dir_shower = shower_cal_dir_3vector(*shower, test_p, 30 * units::cm);
        }
        if (dir_shower.magnitude() < 0.001) dir_shower = dir_main;
        
        // Add segments from other clusters
        for (auto& [seg1, vertices] : map_segment_vertices) {
            if (seg1->cluster() == main_cluster) continue;
            if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) continue;
            if (seg1->cluster() == shower->start_segment()->cluster()) continue;
            
            auto it1 = map_cluster_associated_vertex.find(map_segment_cluster[seg1]);
            
            auto [pair_dis, pair_point] = segment_get_closest_point(seg1, start_pt);
            WireCell::Vector v1(pair_point.x() - start_pt.x(), pair_point.y() - start_pt.y(), pair_point.z() - start_pt.z());
            WireCell::Vector v2(pair_point.x() - point.x(), pair_point.y() - point.y(), pair_point.z() - point.z());
            
            double angle_v1 = std::acos(std::clamp(dir_shower.dot(v1) / (dir_shower.magnitude() * v1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
            double angle_v2 = std::acos(std::clamp(dir_shower.dot(v2) / (dir_shower.magnitude() * v2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
            
            const auto& seg1_fits = seg1->fits();
            double tmp_shower_dis = !seg1_fits.empty() ? (shower->start_segment()->fits().front().point - seg1_fits.front().point).magnitude() : 1e9;
            double close_shower_dis = shower_get_closest_dis(*shower, seg1);
            
            if (angle_v2 > 30) continue;
            
            if ((angle_v1 < 25 && (pair_dis < 80 * units::cm || close_shower_dis < 25 * units::cm)) ||
                (angle_v2 < 25 && (tmp_shower_dis < 40 * units::cm || close_shower_dis < 25 * units::cm)) ||
                (angle_v1 < 12.5 && (pair_dis < 120 * units::cm || close_shower_dis < 40 * units::cm)) ||
                (angle_v2 < 12.5 && (tmp_shower_dis < 80 * units::cm || close_shower_dis < 40 * units::cm))) {
                
                if (it1 != map_cluster_associated_vertex.end() && seg1->cluster() != shower->start_segment()->cluster()) {
                    if (it1->second != vertex) {
                        double dis1 = shower_get_dis(*shower, seg1);
                        if (dis1 > 25 * units::cm && dis1 > pair_dis * 0.4) {
                            continue;
                        }
                    }
                }
                shower->add_segment(seg1);
            }
        }
        
        // Update particle type
        shower->update_particle_type(particle_data, recomb_model);
        
        bool tmp_flag = (shower->start_vertex() == main_vertex);
        std::cout << "Separated shower: " << shower->start_segment()->cluster()->get_cluster_id() * 1000 + shower->start_segment()->id() 
                  << " " << (shower->start_segment()->has_particle_info() && shower->start_segment()->particle_info() ? shower->start_segment()->particle_info()->pdg() : 0)
                  << " " << shower->get_num_segments() << " " << tmp_flag << " " << pi.min_dis / units::cm << std::endl;
        
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters);
        
        // Iteratively add segments based on distance
        {
            std::map<ShowerPtr, double> map_shower_length;
            std::map<ShowerPtr, WireCell::Vector> map_shower_dir;
            
            for (auto shower1 : showers) {
                map_shower_length[shower1] = shower1->get_total_length();
                auto [start_vtx1, _] = shower1->get_start_vertex_and_type();
                WireCell::Point start_pt1 = start_vtx1 ? (start_vtx1->fit().valid() ? start_vtx1->fit().point : start_vtx1->wcpt().point) : WireCell::Point(0, 0, 0);
                auto [__, test_p] = shower_get_closest_point(*shower1, start_pt1);
                map_shower_dir[shower1] = shower_cal_dir_3vector(*shower1, test_p, 30 * units::cm);
            }
            
            bool flag_continue = true;
            while (flag_continue) {
                flag_continue = false;
                for (auto& [seg1, vertices] : map_segment_vertices) {
                    if (seg1->cluster() == main_cluster) continue;
                    if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) continue;
                    
                    double min_dis = 1e9;
                    ShowerPtr min_shower = nullptr;
                    
                    for (auto shower1 : showers) {
                        if (segment_track_length(seg1) > 0.75 * map_shower_length[shower1]) continue;
                        
                        auto [start_vtx1, _] = shower1->get_start_vertex_and_type();
                        WireCell::Point start_point1 = start_vtx1 ? (start_vtx1->fit().valid() ? start_vtx1->fit().point : start_vtx1->wcpt().point) : WireCell::Point(0, 0, 0);
                        auto [__, test_p] = segment_get_closest_point(seg1, start_point1);
                        auto [___, test_p1] = shower_get_closest_point(*shower1, start_point1);
                        
                        WireCell::Vector tmp_dir(test_p.x() - test_p1.x(), test_p.y() - test_p1.y(), test_p.z() - test_p1.z());
                        double angle = std::acos(std::clamp(tmp_dir.dot(map_shower_dir[shower1]) / (tmp_dir.magnitude() * map_shower_dir[shower1].magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        
                        double dis = shower_get_closest_dis(*shower1, seg1);
                        if (dis < min_dis && angle < 45) {
                            min_dis = dis;
                            min_shower = shower1;
                        }
                    }
                    
                    if (min_shower && min_dis < 3.5 * units::cm) {
                        min_shower->add_segment(seg1);
                        map_shower_length[min_shower] = min_shower->get_total_length();
                        flag_continue = true;
                    }
                }
                update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters);
            }
        }
    }
    
    std::cout << "With separated-cluster shower: " << showers.size() << std::endl;
}

void PatternAlgorithms::examine_merge_showers(std::set<ShowerPtr>& showers, VertexPtr main_vertex,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters, std::set<VertexPtr>& vertices_in_long_muon, std::set<SegmentPtr>& segments_in_long_muon, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){

    if (!main_vertex || map_vertex_to_shower.find(main_vertex) == map_vertex_to_shower.end()) {
        return;
    }
    
    // Get all showers starting at main vertex
    auto& main_vertex_showers = map_vertex_to_shower[main_vertex];
    
    // Track which showers have been merged
    std::set<ShowerPtr> showers_to_remove;
    bool flag_update = false;
    
    // Iterate through all showers from main vertex
    for (auto shower1 : main_vertex_showers) {
        // Skip if already processed
        if (showers_to_remove.find(shower1) != showers_to_remove.end()) continue;
        
        // Skip muons (particle type 13)
        if (shower1->get_particle_type() == 13) continue;
        
        // Check start vertex type
        auto [start_vtx1, start_type1] = shower1->get_start_vertex_and_type();
        if (start_type1 != 1) continue;
        
        // Calculate direction for shower1 (100 cm distance cut)
        WireCell::Point start_point1 = shower1->get_start_point();
        WireCell::Vector dir1 = shower_cal_dir_3vector(*shower1, start_point1, 100 * units::cm);
        
        // Look for candidate showers to merge
        for (auto shower2 : main_vertex_showers) {
            // Skip if same shower
            if (shower1 == shower2) continue;
            
            // Skip if already processed
            if (showers_to_remove.find(shower2) != showers_to_remove.end()) continue;
            
            // Skip muons
            if (shower2->get_particle_type() == 13) continue;
            
            // Check start vertex type
            auto [start_vtx2, start_type2] = shower2->get_start_vertex_and_type();
            if (start_type2 != 2) continue;
            
            // Calculate direction for shower2
            WireCell::Point start_point2 = shower2->get_start_point();
            WireCell::Vector dir2 = shower_cal_dir_3vector(*shower2, start_point2, 100 * units::cm);
            
            // Calculate angle between directions
            double cos_angle = dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude());
            double angle = std::acos(std::max(-1.0, std::min(1.0, cos_angle)));
            double angle_deg = angle * 180.0 / M_PI;
            
            // Merge if angle is less than 10 degrees
            if (angle_deg < 10.0) {
                // Merge shower2 into shower1
                shower1->add_shower(*shower2);
                
                // Update particle type and kinematics
                shower1->update_particle_type(particle_data, recomb_model);
                shower1->calculate_kinematics(particle_data, recomb_model);
                
                // Recalculate kinetic charge
                double kine_charge = cal_kine_charge(shower1, graph, track_fitter, dv);
                shower1->set_kine_charge(kine_charge);
                shower1->set_flag_kinematics(true);
                
                // Mark shower2 for removal
                showers_to_remove.insert(shower2);
                flag_update = true;
            }
        }
    }
    
    // Remove merged showers from the main set
    for (auto shower : showers_to_remove) {
        showers.erase(shower);
    }
    
    // Update shower maps if any merges occurred
    if (flag_update) {
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, 
                          map_vertex_to_shower, used_shower_clusters);
    }
}


void PatternAlgorithms::shower_clustering_in_other_clusters(Graph& graph, VertexPtr main_vertex, std::set<ShowerPtr>& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, std::map<Facade::Cluster*, VertexPtr> map_cluster_main_vertices,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_save){
    
    if (!main_vertex || !main_cluster) return;
    
    // Build map_vertex_segments and map_segment_vertices
    std::map<VertexPtr, std::set<SegmentPtr>> map_vertex_segments;
    std::map<SegmentPtr, std::set<VertexPtr>> map_segment_vertices;
    
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr seg = graph[*eit].segment;
        if (!seg) continue;
        
        auto source_vdesc = boost::source(*eit, graph);
        auto target_vdesc = boost::target(*eit, graph);
        VertexPtr v1 = graph[source_vdesc].vertex;
        VertexPtr v2 = graph[target_vdesc].vertex;
        
        if (v1) {
            map_vertex_segments[v1].insert(seg);
            map_segment_vertices[seg].insert(v1);
        }
        if (v2) {
            map_vertex_segments[v2].insert(seg);
            map_segment_vertices[seg].insert(v2);
        }
    }
    
    // Build map_cluster_length
    std::map<Facade::Cluster*, double> map_cluster_length;
    for (auto cluster : other_clusters) {
        map_cluster_length[cluster] = cluster->get_length();
    }
    
    // Collect vertices in main cluster as well as existing showers
    std::vector<VertexPtr> vertices;
    for (auto& [vtx, seg_set] : map_vertex_segments) {
        if ((vtx->cluster() && vtx->cluster()->get_cluster_id() == main_cluster->get_cluster_id()) ||
            map_vertex_in_shower.find(vtx) != map_vertex_in_shower.end()) {
            vertices.push_back(vtx);
        }
    }
    
    // Process clusters in map_cluster_main_vertices
    for (auto& [cluster, vertex] : map_cluster_main_vertices) {
        if (used_shower_clusters.find(cluster) != used_shower_clusters.end()) continue;
        if (map_cluster_length[cluster] < 4 * units::cm) continue;
        
        double min_dis = 1e9;
        VertexPtr min_vertex = nullptr;
        double main_dis = 1e9;
        
        for (auto vtx : vertices) {
            WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
            
            // Get closest distance using cluster's get_closest_dis method
            double dis = cluster->get_closest_dis(vtx_pt);
            
            if (dis < min_dis) {
                min_dis = dis;
                min_vertex = vtx;
            }
            
            if (vtx == main_vertex) {
                main_dis = dis;
            }
        }
        
        if (min_dis > 0.8 * main_dis) {
            min_dis = main_dis;
            min_vertex = main_vertex;
        }
        
        // Find a shower segment starting at vertex
        SegmentPtr sg = nullptr;
        if (map_vertex_segments.find(vertex) != map_vertex_segments.end()) {
            for (auto seg : map_vertex_segments[vertex]) {
                if (seg->flags_any(SegmentFlags::kShowerTrajectory) || 
                    seg->flags_any(SegmentFlags::kShowerTopology)) {
                    sg = seg;
                    break;
                }
            }
        }
        
        int connection_type = 3;
        
        if (sg) {
            // Create new shower
            ShowerPtr shower = std::make_shared<Shower>(graph);
            shower->set_start_vertex(min_vertex, connection_type);
            shower->set_start_segment(sg);
            
            // Set direction based on distance to vertex
            WireCell::Point vertex_pt = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
            const auto& fits = sg->fits();
            if (!fits.empty()) {
                double dis1 = (vertex_pt - fits.front().point).magnitude();
                double dis2 = (vertex_pt - fits.back().point).magnitude();
                
                if (dis1 < dis2) {
                    sg->dirsign(1);
                } else {
                    sg->dirsign(-1);
                }
            }
            
            // Complete shower structure
            std::set<SegmentPtr> used_segments;
            shower->complete_structure_with_start_segment(used_segments);
            
            // Calculate shower direction
            WireCell::Vector dir_shower = shower_cal_dir_3vector(*shower, vertex_pt, 15 * units::cm);
            
            // Cluster with the rest - add segments based on angle and distance
            for (auto& [seg1, vertices_set] : map_segment_vertices) {
                if (seg1->cluster() == main_cluster) continue;
                if (map_segment_in_shower.find(seg1) != map_segment_in_shower.end()) continue;
                if (seg1->cluster() == shower->start_segment()->cluster()) continue;
                
                // Find the closest point
                auto [pair_dis, pair_point] = segment_get_closest_point(seg1, vertex_pt);
                WireCell::Vector v1(pair_point.x() - vertex_pt.x(), 
                                   pair_point.y() - vertex_pt.y(), 
                                   pair_point.z() - vertex_pt.z());
                
                double angle = std::acos(std::clamp(dir_shower.dot(v1) / (dir_shower.magnitude() * v1.magnitude()), -1.0, 1.0));
                angle = angle / M_PI * 180.0;
                
                if ((angle < 25 && pair_dis < 80 * units::cm) ||
                    (angle < 12.5 && pair_dis < 120 * units::cm)) {
                    shower->add_segment(seg1);
                }
            }
            
            // Update particle type
            shower->update_particle_type(particle_data, recomb_model);
            
            // Check with other showers and merge if needed
            std::vector<ShowerPtr> showers_to_be_removed;
            for (auto shower1 : showers) {
                WireCell::Point start_pt1 = shower1->get_start_point();
                WireCell::Vector dir_shower1 = shower_cal_dir_3vector(*shower1, start_pt1, 15 * units::cm);
                WireCell::Vector dir2(start_pt1.x() - vertex_pt.x(),
                                     start_pt1.y() - vertex_pt.y(),
                                     start_pt1.z() - vertex_pt.z());
                
                double angle = std::acos(std::clamp(dir_shower.dot(dir_shower1) / (dir_shower.magnitude() * dir_shower1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                double angle1 = std::acos(std::clamp(dir_shower.dot(dir2) / (dir_shower.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                
                if ((angle < 25 && angle1 < 15 && dir2.magnitude() < 80 * units::cm) ||
                    (angle < 12.5 && angle1 < 7.5 && dir2.magnitude() < 120 * units::cm)) {
                    shower->add_shower(*shower1);
                    showers_to_be_removed.push_back(shower1);
                }
            }
            
            for (auto shower_to_remove : showers_to_be_removed) {
                showers.erase(shower_to_remove);
            }
            
            showers.insert(shower);
        }
    }
    
    update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, 
                      map_vertex_to_shower, used_shower_clusters);
    
    // Process remaining other_clusters not in map_cluster_main_vertices
    for (auto cluster : other_clusters) {
        if (used_shower_clusters.find(cluster) != used_shower_clusters.end()) continue;
        
        double min_dis = 1e9;
        VertexPtr min_vertex = nullptr;
        double main_dis = 1e9;
        
        for (auto vtx : vertices) {
            WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
            
            // Get closest distance using cluster's get_closest_dis method
            double dis = cluster->get_closest_dis(vtx_pt);
            
            if (dis < min_dis) {
                min_dis = dis;
                min_vertex = vtx;
            }
            
            if (vtx == main_vertex) {
                main_dis = dis;
            }
        }
        
        if (min_dis > 0.8 * main_dis) {
            min_dis = main_dis;
            min_vertex = main_vertex;
        }
        
        // Find a segment from this cluster
        SegmentPtr sg = nullptr;
        for (auto& [seg, vertices_set] : map_segment_vertices) {
            if (seg->cluster() != cluster) continue;
            sg = seg;
            break;
        }
        
        int connection_type = 3;
        if (min_dis > 80 * units::cm) {
            connection_type = 4;
        }

        if (!flag_save) connection_type = 4;
        
        if (sg) {
            // Create new shower
            ShowerPtr shower = std::make_shared<Shower>(graph);
            shower->set_start_vertex(min_vertex, connection_type);
            shower->set_start_segment(sg);
            
            // Set direction if not already set
            if (sg->dirsign() == 0) {
                // Get start and end vertices
                auto seg_edesc = sg->get_descriptor();
                auto source_vdesc = boost::source(seg_edesc, graph);
                auto target_vdesc = boost::target(seg_edesc, graph);
                VertexPtr v1 = graph[source_vdesc].vertex;
                VertexPtr v2 = graph[target_vdesc].vertex;
                
                if (v1 && v2) {
                    if (map_vertex_segments[v1].size() == 1 && map_vertex_segments[v2].size() > 1) {
                        sg->dirsign(1);
                    } else if (map_vertex_segments[v1].size() > 1 && map_vertex_segments[v2].size() == 1) {
                        sg->dirsign(-1);
                    } else {
                        // Examine vertices based on distance to main_vertex
                        WireCell::Point main_vtx_pt = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
                        const auto& fits = sg->fits();
                        if (!fits.empty()) {
                            double dis1 = (main_vtx_pt - fits.front().point).magnitude();
                            double dis2 = (main_vtx_pt - fits.back().point).magnitude();
                            
                            if (dis1 < dis2) {
                                sg->dirsign(1);
                            } else {
                                sg->dirsign(-1);
                            }
                        }
                    }
                }
            }
            
            // Set particle type to electron if needed
            int particle_type = 0;
            if (sg->has_particle_info() && sg->particle_info()) {
                particle_type = sg->particle_info()->pdg();
            }
            
            if (particle_type == 0 || 
                (std::abs(particle_type) == 13 && segment_track_length(sg) < 40 * units::cm && sg->dir_weak())) {
                auto four_momentum = segment_cal_4mom(sg, 11, particle_data, recomb_model);
                
                // Create ParticleInfo for electron
                auto pinfo = std::make_shared<Aux::ParticleInfo>(
                    11,                                          // electron PDG
                    particle_data->get_particle_mass(11),       // electron mass
                    particle_data->pdg_to_name(11),             // "electron"
                    four_momentum                                // 4-momentum
                );
                
                sg->particle_info(pinfo);
            }
            
            // Complete shower structure
            std::set<SegmentPtr> used_segments;
            shower->complete_structure_with_start_segment(used_segments);
            showers.insert(shower);
        }
    }
    
    update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, 
                      map_vertex_to_shower, used_shower_clusters);
}


void PatternAlgorithms::examine_shower_1(Graph& graph, VertexPtr main_vertex, std::set<ShowerPtr>& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, std::map<Facade::Cluster*, VertexPtr> map_cluster_main_vertices,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){

    if (!main_vertex) return;
    
    // Build map_vertex_segments
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
    
    // Check if there is already a large EM shower connecting to main_vertex
    bool flag_skip = false;
    auto it = map_vertex_to_shower.find(main_vertex);
    if (it != map_vertex_to_shower.end()) {
        for (auto shower : it->second) {
            if (shower->start_segment() && shower->start_segment()->has_particle_info() && 
                shower->start_segment()->particle_info()->pdg() != 11) continue;
            
            auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
            if (conn_type != 1) continue;
            
            double energy = 0;
            if (shower->get_kine_best() != 0) {
                energy = shower->get_kine_best();
            } else {
                energy = shower->get_kine_charge();
            }
            if (energy > 80 * units::MeV) flag_skip = true;
        }
    }
    
    bool flag_added = false;
    
    if (!flag_skip) {
        std::set<ShowerPtr> used_showers;
        std::map<SegmentPtr, std::set<ShowerPtr>> map_segment_showers;
        std::map<SegmentPtr, ShowerPtr> map_segment_new_shower;
        std::set<SegmentPtr> used_segments;
        std::set<ShowerPtr> del_showers;
        
        // Loop over segments at main_vertex
        if (map_vertex_segments.find(main_vertex) != map_vertex_segments.end()) {
            WireCell::Point main_vtx_pt = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
            
            for (auto sg : map_vertex_segments[main_vertex]) {
                // Skip strong direction segments or certain particle types
                if (!sg->dir_weak()) continue;
                
                int particle_type = 0;
                if (sg->has_particle_info() && sg->particle_info()) {
                    particle_type = sg->particle_info()->pdg();
                }
                
                double medium_dQ_dx = segment_median_dQ_dx(sg);
                if (particle_type == 2212 && medium_dQ_dx / (43e3 / units::cm) > 1.6) continue;
                if (particle_type == 11) continue;
                
                // Form a new shower
                ShowerPtr shower1 = std::make_shared<Shower>(graph);
                shower1->set_start_vertex(main_vertex, 1);
                shower1->set_start_segment(sg);
                shower1->complete_structure_with_start_segment(used_segments);
                
                WireCell::Vector dir1 = segment_cal_dir_3vector(sg, main_vtx_pt, 15 * units::cm);
                
                // Check against existing showers
                for (auto shower : showers) {
                    double energy = shower->get_kine_charge();
                    double min_dis = shower_get_closest_dis(*shower, sg);
                    
                    auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
                    
                    if (conn_type > 2 && min_dis > 3 * units::cm) continue;
                    if (shower->start_segment() && shower->start_segment()->has_particle_info() &&
                        shower->start_segment()->particle_info()->pdg() != 11) continue;
                    if (start_vtx == main_vertex && conn_type == 1) continue;
                    
                    // Get shower1's vertices
                    TrajectoryView& traj1 = shower1->fill_maps();
                    std::set<VertexPtr> shower1_vertices;
                    for (auto vdesc : traj1.nodes()) {
                        auto vtx = traj1.view_graph()[vdesc].vertex;
                        if (vtx) shower1_vertices.insert(vtx);
                    }
                    
                    if (conn_type == 1 && energy > 80 * units::MeV) {
                        if (shower1_vertices.find(start_vtx) == shower1_vertices.end()) continue;
                        else {
                            WireCell::Vector dir3 = shower_cal_dir_3vector(*shower, shower->get_start_point(), 15 * units::cm);
                            double angle = std::acos(std::clamp(dir3.dot(dir1) / (dir3.magnitude() * dir1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                            if (angle > 30) continue;
                        }
                    } else {
                        WireCell::Vector dir3 = shower_cal_dir_3vector(*shower, shower->get_start_point(), 15 * units::cm);
                        
                        // Find closest vertex in shower1 to shower start point
                        WireCell::Point min_point;
                        double min_vtx_dis = 1e9;
                        for (auto vtx3 : shower1_vertices) {
                            WireCell::Point vtx3_pt = vtx3->fit().valid() ? vtx3->fit().point : vtx3->wcpt().point;
                            double dis = (vtx3_pt - shower->get_start_point()).magnitude();
                            if (dis < min_vtx_dis) {
                                min_vtx_dis = dis;
                                min_point = vtx3_pt;
                            }
                        }
                        
                        WireCell::Vector dir4(shower->get_start_point().x() - min_point.x(),
                                             shower->get_start_point().y() - min_point.y(),
                                             shower->get_start_point().z() - min_point.z());
                        
                        double angle3 = std::acos(std::clamp(dir3.dot(dir1) / (dir3.magnitude() * dir1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        double angle4 = std::acos(std::clamp(dir4.dot(dir1) / (dir4.magnitude() * dir1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        double tmp_angle = std::min(angle3, angle4);
                        
                        if (energy > 25 * units::MeV && tmp_angle > 40) continue;
                    }
                    
                    if (used_showers.find(shower) != used_showers.end()) continue;
                    
                    auto [shower_dis, shower_point] = shower_get_closest_point(*shower, main_vtx_pt);
                    WireCell::Vector dir2(shower_point.x() - main_vtx_pt.x(),
                                         shower_point.y() - main_vtx_pt.y(),
                                         shower_point.z() - main_vtx_pt.z());
                    
                    double angle = std::acos(std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    
                    if ((angle < 15 && min_dis < 36 * units::cm) ||
                        (angle < 10 && min_dis < 46 * units::cm) ||
                        (angle < 7.5)) {
                        map_segment_showers[sg].insert(shower);
                        used_showers.insert(shower);
                    }
                }
                
                if (map_segment_showers.find(sg) != map_segment_showers.end()) {
                    map_segment_new_shower[sg] = shower1;
                } else {
                    // Shower1 will be deleted by shared_ptr
                }
            }
        }
        
        // Process segments with associated showers
        for (auto& [sg, associated_showers] : map_segment_showers) {
            ShowerPtr shower1 = map_segment_new_shower[sg];
            int num_showers = associated_showers.size();
            
            double max_energy = 0;
            double total_energy = 0;
            for (auto shower : associated_showers) {
                double energy = shower->get_kine_charge();
                if (energy > max_energy) max_energy = energy;
                total_energy += energy;
            }
            
            // Analyze shower1 structure
            TrajectoryView& traj = shower1->fill_maps();
            std::set<VertexPtr> shower1_vertices;
            for (auto vdesc : traj.nodes()) {
                auto vtx = traj.view_graph()[vdesc].vertex;
                if (vtx) shower1_vertices.insert(vtx);
            }
            
            double max_length = 0;
            SegmentPtr max_sg = nullptr;
            int n_tracks = 0;
            int n_showers = 0;
            double total_length = 0;
            bool flag_good_track = false;
            
            for (auto edesc : traj.edges()) {
                auto sg1 = traj.view_graph()[edesc].segment;
                if (!sg1) continue;
                
                double length = segment_track_length(sg1);
                double medium_dQ_dx = segment_median_dQ_dx(sg1) / (43e3 / units::cm);
                
                if (!sg1->dir_weak()) {
                    // Find end vertex
                    auto seg_edesc = sg1->get_descriptor();
                    auto source_vdesc = boost::source(seg_edesc, graph);
                    auto target_vdesc = boost::target(seg_edesc, graph);
                    VertexPtr v1 = graph[source_vdesc].vertex;
                    VertexPtr v2 = graph[target_vdesc].vertex;
                    
                    VertexPtr end_vertex = nullptr;
                    if (sg1->dirsign() == 1) {
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
                                bool is_shower = sg2->flags_any(SegmentFlags::kShowerTrajectory) || 
                                                sg2->flags_any(SegmentFlags::kShowerTopology);
                                if (!is_shower) flag_non_ele = true;
                            }
                            if (!flag_non_ele && map_vertex_segments[end_vertex].size() <= 3) {
                                flag_good_track = true;
                            }
                        } else {
                            flag_good_track = true;
                        }
                    }
                }
                
                bool is_shower = sg1->flags_any(SegmentFlags::kShowerTrajectory) || 
                                sg1->flags_any(SegmentFlags::kShowerTopology);
                if (is_shower) n_showers++;
                
                n_tracks++;
                total_length += length;
                if (max_length < length) {
                    max_length = length;
                    max_sg = sg1;
                }
                (void)medium_dQ_dx; // To avoid unused variable warning
            }
            (void)flag_good_track;
            (void)n_showers;           

            // Check if should skip
            bool flag_skip_segment = false;
            for (auto shower_check : showers) {
                auto [start_vtx_check, conn_type_check] = shower_check->get_start_vertex_and_type();
                if (conn_type_check == 1 && associated_showers.find(shower_check) == associated_showers.end()) {
                    if (shower1_vertices.find(start_vtx_check) != shower1_vertices.end() && 
                        start_vtx_check != main_vertex && 
                        shower_check->get_kine_charge() > 60 * units::MeV) {
                        flag_skip_segment = true;
                    }
                }
            }
            
            // Decide whether to create this shower
            if (total_length < 70 * units::cm && 
                ((n_tracks == 1 && total_length < 60 * units::cm) ||
                 (n_tracks == 1 && total_length < 65 * units::cm && num_showers > 3 && total_energy > 150 * units::MeV) ||
                 total_length < n_tracks * 36 * units::cm) &&
                (total_energy > 50 * units::MeV || total_energy / units::MeV > total_length / units::cm * 0.75) &&
                !flag_skip_segment) {
                
                // Set particle type to electron
                if (shower1->start_segment() && shower1->start_segment()->has_particle_info() && 
                    shower1->start_segment()->particle_info()) {
                    shower1->start_segment()->particle_info()->set_pdg(11);
                }
                shower1->start_segment()->set_flags(SegmentFlags::kAvoidMuonCheck);
                shower1->update_particle_type(particle_data, recomb_model);
                
                // Merge associated showers
                for (auto shower : associated_showers) {
                    del_showers.insert(shower);
                    shower1->add_shower(*shower);
                }
                
                shower1->calculate_kinematics(particle_data, recomb_model);
                double kine_charge = cal_kine_charge(shower1, graph, track_fitter, dv);
                shower1->set_kine_charge(kine_charge);
                shower1->set_flag_kinematics(true);
                
                showers.insert(shower1);
                std::cout << "Create a new low-energy shower: " << kine_charge / units::MeV << " MeV" << std::endl;
                flag_added = true;
            }
        }
        
        // Remove deleted showers
        for (auto shower : del_showers) {
            showers.erase(shower);
        }
        
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, 
                          map_vertex_to_shower, used_shower_clusters);
    }
    
    // Second part: merge existing showers if nothing was added
    if (!flag_added) {
        std::map<ShowerPtr, std::set<ShowerPtr>> map_shower_showers;
        WireCell::Point main_vtx_pt = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
        
        auto it = map_vertex_to_shower.find(main_vertex);
        if (it != map_vertex_to_shower.end()) {
            for (auto shower : it->second) {
                if (shower->start_segment() && shower->start_segment()->has_particle_info() &&
                    shower->start_segment()->particle_info()->pdg() != 11) continue;
                
                auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
                if (conn_type != 1) continue;
                
                WireCell::Vector dir1 = shower_cal_dir_3vector(*shower, main_vtx_pt, 15 * units::cm);
                
                for (auto shower1 : showers) {
                    double energy = shower1->get_kine_charge();
                    double min_dis = shower_get_closest_dis(*shower1, shower->start_segment());
                    
                    auto [start_vtx1, conn_type1] = shower1->get_start_vertex_and_type();
                    
                    if (conn_type1 > 2 && min_dis > 3 * units::cm) continue;
                    if (shower1->start_segment() && shower1->start_segment()->has_particle_info() &&
                        shower1->start_segment()->particle_info()->pdg() != 11) continue;
                    if (start_vtx1 == main_vertex && conn_type1 == 1) continue;
                    
                    // Get shower's vertices
                    TrajectoryView& traj_shower = shower->fill_maps();
                    std::set<VertexPtr> shower_vertices;
                    for (auto vdesc : traj_shower.nodes()) {
                        auto vtx = traj_shower.view_graph()[vdesc].vertex;
                        if (vtx) shower_vertices.insert(vtx);
                    }
                    
                    if (conn_type1 == 1 && shower_vertices.find(start_vtx1) == shower_vertices.end()) continue;
                    
                    if (shower1->get_total_length() < 3 * units::cm) continue;
                    
                    auto [shower_dis, shower_point] = shower_get_closest_point(*shower1, main_vtx_pt);
                    WireCell::Vector dir2(shower_point.x() - main_vtx_pt.x(),
                                         shower_point.y() - main_vtx_pt.y(),
                                         shower_point.z() - main_vtx_pt.z());
                    
                    double angle = std::acos(std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    
                    WireCell::Vector dir3 = shower_cal_dir_3vector(*shower1, shower1->get_start_point(), 30 * units::cm);
                    double angle1 = std::acos(std::clamp(dir2.dot(dir3) / (dir2.magnitude() * dir3.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    
                    if (angle < 15 && angle1 < 15 && min_dis < 28 * units::cm) {
                        map_shower_showers[shower].insert(shower1);
                    }
                    (void)energy; // To avoid unused variable warning
                }
            }
        }
        
        // Find the shower combination with maximum energy
        std::set<ShowerPtr> del_showers;
        ShowerPtr max_shower = nullptr;
        double max_energy = 0;
        
        for (auto& [shower, associated_showers] : map_shower_showers) {
            double acc_energy = 0;
            if (shower->get_kine_best() != 0) {
                acc_energy += shower->get_kine_best();
            } else {
                acc_energy += shower->get_kine_charge();
            }
            
            for (auto shower1 : associated_showers) {
                if (shower1->get_kine_best() != 0) {
                    acc_energy += shower1->get_kine_best();
                } else {
                    acc_energy += shower1->get_kine_charge();
                }
            }
            
            if (acc_energy > max_energy) {
                max_energy = acc_energy;
                max_shower = shower;
            }
        }
        
        if (max_shower) {
            for (auto shower1 : map_shower_showers[max_shower]) {
                max_shower->add_shower(*shower1);
                del_showers.insert(shower1);
            }
            
            max_shower->calculate_kinematics(particle_data, recomb_model);
            max_shower->start_segment()->set_flags(SegmentFlags::kAvoidMuonCheck);
            double kine_charge = cal_kine_charge(max_shower, graph, track_fitter, dv);
            max_shower->set_kine_charge(kine_charge);
            max_shower->set_flag_kinematics(true);
        }
        
        // Remove deleted showers
        for (auto shower : del_showers) {
            showers.erase(shower);
        }
        
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, 
                          map_vertex_to_shower, used_shower_clusters);
    }
}


void PatternAlgorithms::examine_showers(Graph& graph, VertexPtr main_vertex, std::set<ShowerPtr>& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, std::map<Facade::Cluster*, VertexPtr> map_cluster_main_vertices,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){
    
    if (!main_vertex) return;
    
    // Build map_vertex_segments
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
    
    std::map<SegmentPtr, ShowerPtr> map_merge_seg_shower;
    WireCell::Vector drift_dir(1, 0, 0);
    std::set<ShowerPtr> del_showers;
    
    WireCell::Point main_vtx_pt = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
    
    // Loop over segments at main_vertex
    if (map_vertex_segments.find(main_vertex) != map_vertex_segments.end()) {
        for (auto sg : map_vertex_segments[main_vertex]) {
            // Skip if already in shower and not a muon
            if (map_segment_in_shower.find(sg) != map_segment_in_shower.end()) {
                if (sg->has_particle_info() && sg->particle_info()->pdg() != 13) continue;
            }
            
            double sg_length = segment_track_length(sg);
            
            // Skip long segments
            if ((sg_length > 45 * units::cm && !sg->dir_weak()) || sg_length > 55 * units::cm) continue;
            
            // Find other vertex
            auto seg_edesc = sg->get_descriptor();
            auto source_vdesc = boost::source(seg_edesc, graph);
            auto target_vdesc = boost::target(seg_edesc, graph);
            VertexPtr v1 = graph[source_vdesc].vertex;
            VertexPtr v2 = graph[target_vdesc].vertex;
            VertexPtr vtx = (v1 == main_vertex) ? v2 : v1;
            if (!vtx) continue;
            
            auto daughter_result = calculate_num_daughter_showers(graph, main_vertex, sg, false);
            double daughter_length = daughter_result.second;
            
            bool flag_checked = false;
            
            // Case I: Check showers at the other vertex
            if (map_vertex_to_shower.find(vtx) != map_vertex_to_shower.end()) {
                WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
                WireCell::Vector dir1 = segment_cal_dir_3vector(sg, vtx_pt, 15 * units::cm);
                WireCell::Vector dir1_1 = segment_cal_dir_3vector(sg, main_vtx_pt, 15 * units::cm);
                bool flag_tmp_connected = false;
                
                // First pass: check for high-energy directly connected showers
                for (auto shower : map_vertex_to_shower[vtx]) {
                    if (shower->start_segment() && shower->start_segment()->has_particle_info() &&
                        shower->start_segment()->particle_info()->pdg() != 11) continue;
                    
                    auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
                    double Eshower = (shower->get_kine_best() != 0) ? shower->get_kine_best() : shower->get_kine_charge();
                    
                    if (conn_type == 1 && Eshower > 60 * units::MeV) flag_tmp_connected = true;
                }
                
                // Second pass: check merging conditions
                for (auto shower : map_vertex_to_shower[vtx]) {
                    if (shower->start_segment() && shower->start_segment()->has_particle_info() &&
                        shower->start_segment()->particle_info()->pdg() != 11) continue;
                    
                    auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
                    WireCell::Vector dir2 = shower_cal_dir_3vector(*shower, shower->get_start_point(), 100 * units::cm);
                    
                    double Eshower = (shower->get_kine_best() != 0) ? shower->get_kine_best() : shower->get_kine_charge();
                    
                    // Check shower composition
                    TrajectoryView& traj = shower->fill_maps();
                    double tmp_total_length = 0;
                    double tmp_track_length = 0;
                    
                    for (auto edesc : traj.edges()) {
                        auto sg1 = traj.view_graph()[edesc].segment;
                        if (!sg1) continue;
                        if (sg1->cluster() != shower->start_segment()->cluster()) continue;
                        
                        double length = segment_track_length(sg1);
                        tmp_total_length += length;
                        if (!sg1->dir_weak()) tmp_track_length += length;
                    }
                    
                    if (tmp_track_length > 3 * units::cm && tmp_track_length > 0.25 * tmp_total_length) continue;
                    
                    if (conn_type == 1 && Eshower > 100 * units::MeV) flag_checked = true;
                    
                    double angle1 = std::acos(std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    double angle2 = std::acos(std::clamp(dir1_1.dot(dir2) / (dir1_1.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    double tmp_angle = std::min(180.0 - angle1, angle2);
                    
                    auto [closest_dis, closest_pt] = segment_get_closest_point(sg, shower->get_start_point());
                    
                    // Special checks for connection type 2
                    if (conn_type == 2) {
                        if ((!sg->dir_weak() && sg_length > 3 * units::cm) || flag_tmp_connected) {
                            if (closest_dis < 8 * units::cm && Eshower > 75 * units::MeV && tmp_angle < 6) {
                                // Allow this condition
                            } else {
                                continue;
                            }
                        }
                        if (closest_dis > 20 * units::cm && Eshower < 150 * units::MeV && tmp_angle > 2.5) continue;
                    }
                    
                    // Main merging conditions
                    if ((Eshower > 800 * units::MeV && tmp_angle < 30) ||
                        (Eshower > 150 * units::MeV && tmp_angle < 10) ||
                        (Eshower > 150 * units::MeV && tmp_angle < 18 && conn_type == 1 && sg->dir_weak()) ||
                        (Eshower > 100 * units::MeV && tmp_angle < 10 && sg_length < 25 * units::cm) ||
                        (Eshower > 250 * units::MeV && tmp_angle < 15) ||
                        (Eshower > 360 * units::MeV && tmp_angle < 25) ||
                        (Eshower > 100 * units::MeV && Eshower <= 150 * units::MeV && tmp_angle < 15 && 
                         sg_length < 25 * units::cm && flag_checked) ||
                        (Eshower > 60 * units::MeV && conn_type == 2 && sg->dir_weak() &&
                         ((tmp_angle < 15 && closest_dis < 18 * units::cm) || 
                          (tmp_angle < 17.5 && closest_dis < 6 * units::cm)) &&
                         sg_length < 15 * units::cm) ||
                        (Eshower > 60 * units::MeV && conn_type == 2 && tmp_angle < 7.5 && 
                         closest_dis < 8 * units::cm && sg_length < 20 * units::cm)) {
                        map_merge_seg_shower[sg] = shower;
                        continue;
                    }
                }
            }
            
            if (map_merge_seg_shower.find(sg) != map_merge_seg_shower.end()) continue;
            if (flag_checked) continue;
            if (!sg->dir_weak() && sg_length > 6 * units::cm && daughter_length < 40 * units::cm) continue;
            
            // Case II: Check showers at main_vertex (not directly connected)
            if (map_vertex_to_shower.find(main_vertex) != map_vertex_to_shower.end()) {
                WireCell::Vector dir1 = segment_cal_dir_3vector(sg, main_vtx_pt, 15 * units::cm);
                
                for (auto shower : map_vertex_to_shower[main_vertex]) {
                    if (shower->start_segment() && shower->start_segment()->has_particle_info() &&
                        shower->start_segment()->particle_info()->pdg() != 11) continue;
                    
                    auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
                    if (conn_type == 1) continue; // Skip directly connected
                    
                    WireCell::Vector dir2(shower->get_start_point().x() - main_vtx_pt.x(),
                                         shower->get_start_point().y() - main_vtx_pt.y(),
                                         shower->get_start_point().z() - main_vtx_pt.z());
                    WireCell::Vector dir3 = shower_cal_dir_3vector(*shower, shower->get_start_point(), 100 * units::cm);
                    
                    auto [min_dis, closest_pt] = segment_get_closest_point(sg, shower->get_start_point());
                    
                    double angle_dir2 = std::acos(std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    double angle_dir3 = std::acos(std::clamp(dir1.dot(dir3) / (dir1.magnitude() * dir3.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    double angle_drift1 = std::acos(std::clamp(dir1.dot(drift_dir) / (dir1.magnitude() * drift_dir.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    double angle_drift3 = std::acos(std::clamp(dir3.dot(drift_dir) / (dir3.magnitude() * drift_dir.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    
                    if (((shower->get_kine_charge() > 80 * units::MeV && angle_dir2 < 10) ||
                         (shower->get_kine_charge() > 50 * units::MeV && angle_dir2 < 3) ||
                         (shower->get_kine_charge() > 80 * units::MeV && angle_dir3 < 6 && angle_dir2 < 17.5) ||
                         (shower->get_kine_charge() > 80 * units::MeV && angle_dir3 < 6 && 
                          std::fabs(90 - angle_drift1) < 10 && std::fabs(90 - angle_drift3) < 10 && angle_dir2 < 30)) &&
                        (sg_length > 5 * units::cm || (sg_length > 3 * units::cm && min_dis < 2.0 * units::cm))) {
                        map_merge_seg_shower[sg] = shower;
                        continue;
                    }
                }
            }
            
            // Case III: Check other showers
            WireCell::Vector dir1 = segment_cal_dir_3vector(sg, main_vtx_pt, 15 * units::cm);
            
            for (auto shower : showers) {
                auto [start_vtx_shower, conn_type_shower] = shower->get_start_vertex_and_type();
                if (start_vtx_shower == vtx || start_vtx_shower == main_vertex) continue;
                
                if (shower->start_segment() && shower->start_segment()->has_particle_info() &&
                    shower->start_segment()->particle_info()->pdg() != 11) continue;
                
                if (conn_type_shower <= 2) {
                    WireCell::Vector dir2(shower->get_start_point().x() - main_vtx_pt.x(),
                                         shower->get_start_point().y() - main_vtx_pt.y(),
                                         shower->get_start_point().z() - main_vtx_pt.z());
                    WireCell::Vector dir3 = shower_cal_dir_3vector(*shower, shower->get_start_point(), 100 * units::cm);
                    
                    double angle_dir2 = std::acos(std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    double angle_dir3 = std::acos(std::clamp(dir2.dot(dir3) / (dir2.magnitude() * dir3.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    
                    if ((((shower->get_kine_charge() > 80 * units::MeV && angle_dir2 < 15) ||
                          (shower->get_kine_charge() > 50 * units::MeV && angle_dir2 < 5)) &&
                         angle_dir3 < 15 && (angle_dir2 + angle_dir3) < 25 && sg_length > 5 * units::cm)) {
                        map_merge_seg_shower[sg] = shower;
                        continue;
                    }
                }
            }
        }
    }
    
    // Mark showers for deletion if their start segment should be merged
    for (auto shower : showers) {
        SegmentPtr sg1 = shower->start_segment();
        if (sg1 && map_merge_seg_shower.find(sg1) != map_merge_seg_shower.end()) {
            sg1->set_flags(SegmentFlags::kAvoidMuonCheck);
            del_showers.insert(shower);
        }
    }
    
    // Delete marked showers first
    if (del_showers.size() != 0) {
        for (auto shower1 : del_showers) {
            showers.erase(shower1);
        }
        del_showers.clear();
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, 
                          map_vertex_to_shower, used_shower_clusters);
    }
    
    // Perform the merging
    std::set<ShowerPtr> updated_showers;
    for (auto& [sg, shower] : map_merge_seg_shower) {
        std::cout << "EM shower modification: " << shower->start_segment()->id() << " -> " << sg->id() << std::endl;
        updated_showers.insert(shower);
        
        auto [pair_vertex, pair_conn_type] = shower->get_start_vertex_and_type();
        
        if (pair_conn_type != 1) {
            shower->add_segment(sg);
            shower->set_start_vertex(main_vertex, 1);
            shower->set_start_segment(sg);
            shower->set_start_point(main_vtx_pt);
            std::set<SegmentPtr> tmp_used_segments;
            shower->complete_structure_with_start_segment(tmp_used_segments);
            if (segment_track_length(sg) > 44 * units::cm || sg->dir_weak()) {
                sg->set_flags(SegmentFlags::kAvoidMuonCheck);
            }
        } else {
            shower->add_segment(sg);
            shower->set_start_vertex(main_vertex, 1);
            shower->set_start_segment(sg);
            shower->set_start_point(main_vtx_pt);
            std::set<SegmentPtr> tmp_used_segments;
            shower->complete_structure_with_start_segment(tmp_used_segments);
            if (shower->get_num_main_segments() >= 3) {
                sg->set_flags(SegmentFlags::kAvoidMuonCheck);
            }
        }
        
        // Add other showers that connect to this shower
        TrajectoryView& traj = shower->fill_maps();
        std::set<VertexPtr> shower_vertices;
        for (auto vdesc : traj.nodes()) {
            auto vtx = traj.view_graph()[vdesc].vertex;
            if (vtx) shower_vertices.insert(vtx);
        }
        
        for (auto shower1 : showers) {
            if (shower == shower1) continue;
            auto [start_vtx1, conn_type1] = shower1->get_start_vertex_and_type();
            if (conn_type1 == 1 && start_vtx1 != main_vertex) {
                if (shower_vertices.find(start_vtx1) != shower_vertices.end()) {
                    shower->add_shower(*shower1);
                    del_showers.insert(shower1);
                }
            }
        }
        
        // Update particle type and kinematics
        if (sg->has_particle_info() && sg->particle_info()) {
            sg->particle_info()->set_pdg(11);
        }
        shower->update_particle_type(particle_data, recomb_model);
        shower->calculate_kinematics(particle_data, recomb_model);
        double kine_charge = cal_kine_charge(shower, graph, track_fitter, dv);
        shower->set_kine_charge(kine_charge);
        shower->set_flag_kinematics(true);
    }
    
    // Delete merged showers (not at main_vertex)
    for (auto shower1 : del_showers) {
        auto [start_vtx1, conn_type1] = shower1->get_start_vertex_and_type();
        if (start_vtx1 != main_vertex) {
            showers.erase(shower1);
        }
    }
    del_showers.clear();
    
    // Check other showers and merge with updated showers
    for (auto shower : updated_showers) {
        WireCell::Point shower_start = shower->get_start_point();
        WireCell::Vector dir1 = shower_cal_dir_3vector(*shower, shower_start, 25 * units::cm);
        
        for (auto shower1 : showers) {
            if (updated_showers.find(shower1) != updated_showers.end()) continue;
            if (shower1->start_segment() && shower1->start_segment()->has_particle_info() &&
                shower1->start_segment()->particle_info()->pdg() != 11) continue;
            
            auto [start_vtx1, conn_type1] = shower1->get_start_vertex_and_type();
            if (conn_type1 == 2) {
                if (del_showers.find(shower1) != del_showers.end()) continue;
                
                auto [shower_vtx, shower_conn] = shower->get_start_vertex_and_type();
                WireCell::Point shower_vtx_pt = shower_vtx->fit().valid() ? shower_vtx->fit().point : shower_vtx->wcpt().point;
                
                WireCell::Vector dir2(shower1->get_start_point().x() - shower_vtx_pt.x(),
                                     shower1->get_start_point().y() - shower_vtx_pt.y(),
                                     shower1->get_start_point().z() - shower_vtx_pt.z());
                WireCell::Vector dir3 = shower_cal_dir_3vector(*shower1, shower1->get_start_point(), 25 * units::cm);
                
                double angle_dir2 = std::acos(std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                double angle_dir3 = std::acos(std::clamp(dir1.dot(dir3) / (dir1.magnitude() * dir3.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                
                if (angle_dir2 < 10 && angle_dir3 < 20) {
                    shower->add_shower(*shower1);
                    shower->calculate_kinematics(particle_data, recomb_model);
                    double kine_charge = cal_kine_charge(shower, graph, track_fitter, dv);
                    shower->set_kine_charge(kine_charge);
                    shower->set_flag_kinematics(true);
                    del_showers.insert(shower1);
                }
            }
        }
    }
    
    // Final deletion
    for (auto shower1 : del_showers) {
        showers.erase(shower1);
    }
    
    if (map_merge_seg_shower.size() > 0) {
        update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower, 
                          map_vertex_to_shower, used_shower_clusters);
    }
    
    // Finally call examine_shower_1
    examine_shower_1(graph, main_vertex, showers, main_cluster, other_clusters, map_cluster_main_vertices,
                    map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower, used_shower_clusters,
                    track_fitter, dv, particle_data, recomb_model);
}


void PatternAlgorithms::id_pi0_with_vertex(int acc_segment_id, std::set<ShowerPtr>& pi0_showers, std::map<ShowerPtr, int>& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass,  std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Graph& graph, VertexPtr main_vertex, std::set<ShowerPtr>& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, std::map<Facade::Cluster*, VertexPtr> map_cluster_main_vertices,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){

    if (!main_vertex) return;
    
    // Build map_vertex_segments
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
    
    // Figure out all disconnected showers
    std::set<ShowerPtr> disconnected_showers;
    std::map<ShowerPtr, WireCell::Vector> map_shower_dir;
    
    for (auto& [vtx, shower_set] : map_vertex_to_shower) {
        for (auto shower : shower_set) {
            auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
            
            if (conn_type == 2 && std::abs(shower->get_particle_type()) != 13) {
                disconnected_showers.insert(shower);
                map_shower_dir[shower] = shower_cal_dir_3vector(*shower, shower->get_start_point(), 15 * units::cm);
            } else if (conn_type == 1) {
                WireCell::Point start_vtx_pt = start_vtx->fit().valid() ? start_vtx->fit().point : start_vtx->wcpt().point;
                map_shower_dir[shower] = shower_cal_dir_3vector(*shower, start_vtx_pt, 15 * units::cm);
            }
        }
    }
    
    // Build candidate vertices
    std::set<VertexPtr> candidate_vertices;
    candidate_vertices.insert(main_vertex);
    
    for (auto& [vtx, shower_set] : map_vertex_to_shower) {
        bool flag_add = true;
        auto it_in_shower = map_vertex_in_shower.find(vtx);
        if (it_in_shower != map_vertex_in_shower.end()) {
            flag_add = false;
            auto [start_vtx, conn_type] = it_in_shower->second->get_start_vertex_and_type();
            if (vtx == start_vtx) flag_add = true;
        }
        if (flag_add) candidate_vertices.insert(vtx);
    }
    
    // Map shower pairs to masses and vertices
    std::map<std::pair<ShowerPtr, ShowerPtr>, std::vector<std::pair<double, VertexPtr>>> map_shower_pair_mass_vertex;
    
    for (auto vtx : candidate_vertices) {
        std::vector<ShowerPtr> tmp_showers;
        std::map<ShowerPtr, WireCell::Vector> local_dirs;
        
        WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
        
        // Add directly connected showers (type 1, not muon)
        auto it2 = map_vertex_to_shower.find(vtx);
        if (it2 != map_vertex_to_shower.end()) {
            for (auto shower : it2->second) {
                auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
                if (conn_type == 1 && std::abs(shower->get_particle_type()) != 13) {
                    tmp_showers.push_back(shower);
                    local_dirs[shower] = shower->get_init_dir();
                }
            }
        }
        
        // Add disconnected showers within angle
        for (auto shower : disconnected_showers) {
            WireCell::Vector dir1 = map_shower_dir[shower];
            WireCell::Vector dir2(shower->get_start_point().x() - vtx_pt.x(),
                                 shower->get_start_point().y() - vtx_pt.y(),
                                 shower->get_start_point().z() - vtx_pt.z());
            
            auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
            
            if (start_vtx == vtx) {
                tmp_showers.push_back(shower);
                local_dirs[shower] = shower->get_init_dir();
            } else {
                double angle = std::acos(std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                if (angle < 30) {
                    tmp_showers.push_back(shower);
                    local_dirs[shower] = dir2;
                }
            }
        }
        
        // Calculate pi0 masses for all pairs
        if (tmp_showers.size() > 1) {
            for (size_t i = 0; i < tmp_showers.size(); i++) {
                ShowerPtr shower_1 = tmp_showers[i];
                WireCell::Vector dir1 = local_dirs[shower_1];
                
                for (size_t j = i + 1; j < tmp_showers.size(); j++) {
                    ShowerPtr shower_2 = tmp_showers[j];
                    WireCell::Vector dir2 = local_dirs[shower_2];
                    
                    double angle = std::acos(std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0));
                    double mass_pio = std::sqrt(4 * shower_1->get_kine_charge() * shower_2->get_kine_charge() * 
                                               std::pow(std::sin(angle / 2.0), 2));
                    
                    auto [start_vtx_1, conn_type_1] = shower_1->get_start_vertex_and_type();
                    auto [start_vtx_2, conn_type_2] = shower_2->get_start_vertex_and_type();
                    
                    if (conn_type_1 == 1 && conn_type_2 == 1) continue;
                    
                    map_shower_pair_mass_vertex[std::make_pair(shower_1, shower_2)].push_back(std::make_pair(mass_pio, vtx));
                }
            }
        }
    }
    
    // Store pi0 kinematic variables (not implemented as member variables yet, just processing)
    WireCell::Point main_vtx_pt = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
    
    // Find pi0 iteratively
    // static int acc_segment_id = 100000;  // Static counter for pi0 IDs
    
    while (map_shower_pair_mass_vertex.size() > 0) {
        double mass_diff = 1e9;
        double mass_save = 0;
        ShowerPtr shower_1 = nullptr;
        ShowerPtr shower_2 = nullptr;
        double mass_offset = 10 * units::MeV;
        VertexPtr vtx = nullptr;
        double mass_penalty = 0;
        double tmp_mass_penalty = 0;
        
        // Find best pi0 candidate
        for (auto& [shower_pair, mass_vtx_vec] : map_shower_pair_mass_vertex) {
            for (auto& [mass, candidate_vtx] : mass_vtx_vec) {
                auto [start_vtx_1, conn_type_1] = shower_pair.first->get_start_vertex_and_type();
                auto [start_vtx_2, conn_type_2] = shower_pair.second->get_start_vertex_and_type();
                
                if (conn_type_1 == 2 && conn_type_2 == 2) {
                    tmp_mass_penalty = 6 * units::MeV;
                } else {
                    tmp_mass_penalty = 0;
                }
                
                if (mass - 135 * units::MeV + mass_offset < 35 * units::MeV && 
                    mass - 135 * units::MeV + mass_offset > -25 * units::MeV) {
                    if (std::abs(mass - 135 * units::MeV + mass_offset) - tmp_mass_penalty < 
                        std::abs(mass_diff) - mass_penalty) {
                        mass_diff = mass - 135 * units::MeV + mass_offset;
                        mass_penalty = tmp_mass_penalty;
                        mass_save = mass;
                        shower_1 = shower_pair.first;
                        shower_2 = shower_pair.second;
                        vtx = candidate_vtx;
                    }
                }
            }
        }
        
        // If found a good pi0, mark it
        if (mass_diff < 35 * units::MeV && mass_diff > -25 * units::MeV) {
            pi0_showers.insert(shower_1);
            pi0_showers.insert(shower_2);
            
            int pio_id = acc_segment_id;
            acc_segment_id++;
            
            map_shower_pio_id[shower_1] = pio_id;
            map_shower_pio_id[shower_2] = pio_id;
            map_pio_id_mass[pio_id] = std::make_pair(mass_save, 1);
            map_pio_id_showers[pio_id].push_back(shower_1);
            map_pio_id_showers[pio_id].push_back(shower_2);
            
            // Update shower start vertices if needed
            auto [start_vtx_1, conn_type_1] = shower_1->get_start_vertex_and_type();
            if (start_vtx_1 != vtx) {
                shower_1->set_start_vertex(vtx, 2);
                shower_1->calculate_kinematics(particle_data, recomb_model);
            }
            
            auto [start_vtx_2, conn_type_2] = shower_2->get_start_vertex_and_type();
            if (start_vtx_2 != vtx) {
                shower_2->set_start_vertex(vtx, 2);
                shower_2->calculate_kinematics(particle_data, recomb_model);
            }
            
            std::cout << "Pi0 found with mass: " << mass_save / units::MeV << " MeV with " 
                      << shower_1->get_kine_charge() / units::MeV << " MeV + " 
                      << shower_2->get_kine_charge() / units::MeV << " MeV" << std::endl;
        } else {
            break;
        }
        
        // Remove pairs involving these showers
        std::vector<std::pair<ShowerPtr, ShowerPtr>> to_be_removed;
        for (auto& [shower_pair, mass_vtx_vec] : map_shower_pair_mass_vertex) {
            if (shower_pair.first == shower_1 || shower_pair.first == shower_2 ||
                shower_pair.second == shower_1 || shower_pair.second == shower_2) {
                to_be_removed.push_back(shower_pair);
            }
        }
        for (auto& pair : to_be_removed) {
            map_shower_pair_mass_vertex.erase(pair);
        }
    }
    
    // Find pi0 vertices and change incoming muons to pions
    std::set<VertexPtr> pi0_vertices;
    for (auto shower : pi0_showers) {
        auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
        pi0_vertices.insert(start_vtx);
    }
    
    for (auto vtx : pi0_vertices) {
        if (map_vertex_segments.find(vtx) == map_vertex_segments.end()) continue;
        
        WireCell::Point vtx_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
        
        for (auto sg : map_vertex_segments[vtx]) {
            // Determine if segment starts or ends at this vertex
            bool flag_start = false;
            if (!sg->fits().empty()) {
                double dist_front = (sg->fits().front().point - vtx_pt).magnitude();
                double dist_back = (sg->fits().back().point - vtx_pt).magnitude();
                flag_start = (dist_front < dist_back);
            }
            
            int dirsign_val = sg->dirsign();
            
            // Check if segment is coming in (opposite direction)
            bool is_incoming = (flag_start && dirsign_val == -1) || (!flag_start && dirsign_val == 1);
            
            if (is_incoming && sg->has_particle_info()) {
                int pdg = sg->particle_info()->pdg();
                if (std::abs(pdg) == 13 || pdg == 0) {
                    // Change muon to pion
                    sg->particle_info()->set_pdg(211);
                    sg->particle_info()->set_mass(139.57 * units::MeV);  // Pion mass
                }
            }
        }
    }
}


void PatternAlgorithms::id_pi0_without_vertex(int acc_segment_id, std::set<ShowerPtr>& pi0_showers, std::map<ShowerPtr, int>& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass,  std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Graph& graph, VertexPtr main_vertex, std::set<ShowerPtr>& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, std::map<Facade::Cluster*, VertexPtr> map_cluster_main_vertices,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){

    if (!main_vertex) return;
    
    // Build map_vertex_segments and segments_in_long_muon
    std::map<VertexPtr, std::set<SegmentPtr>> map_vertex_segments;
    std::set<SegmentPtr> segments_in_long_muon;  // Placeholder - would need proper implementation
    
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
    
    // Check main vertex conditions
    if (map_vertex_segments[main_vertex].size() > 2) return;
    
    if (map_vertex_segments[main_vertex].size() > 0) {
        auto first_seg = *map_vertex_segments[main_vertex].begin();
        auto last_seg = *map_vertex_segments[main_vertex].rbegin();
        
        if ((map_segment_in_shower.find(first_seg) == map_segment_in_shower.end() &&
             map_segment_in_shower.find(last_seg) == map_segment_in_shower.end()) ||
            segments_in_long_muon.find(first_seg) != segments_in_long_muon.end() ||
            segments_in_long_muon.find(last_seg) != segments_in_long_muon.end()) {
            return;
        }
    }
    
    // Build good_showers set
    std::set<ShowerPtr> good_showers;
    {
        auto it = map_vertex_to_shower.find(main_vertex);
        if (it != map_vertex_to_shower.end()) {
            for (auto shower : it->second) {
                if (pi0_showers.find(shower) != pi0_showers.end()) return;
                
                auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
                if (conn_type == 1) {
                    good_showers.insert(shower);
                }
            }
        }
        
        if (good_showers.size() > 1) {
            ShowerPtr max_shower = nullptr;
            double max_energy = 0;
            for (auto shower : good_showers) {
                double energy = shower->get_kine_charge();
                if (energy > max_energy) {
                    max_energy = energy;
                    max_shower = shower;
                }
            }
            good_showers.clear();
            if (max_shower) good_showers.insert(max_shower);
        }
    }
    
    // Check if we have exactly 2 segments at main vertex
    if (map_vertex_segments[main_vertex].size() == 2) {
        bool flag_return = true;
        int num_showers = 0;
        
        for (auto sg : map_vertex_segments[main_vertex]) {
            if (map_segment_in_shower.find(sg) == map_segment_in_shower.end()) {
                double sg_length = segment_track_length(sg);
                if (sg_length < 1.2 * units::cm && (sg->dirsign() == 0 || sg->dir_weak())) {
                    flag_return = false;
                }
            } else {
                num_showers++;
            }
        }
        
        if (flag_return && num_showers == map_vertex_segments[main_vertex].size()) {
            flag_return = false;
        }
        if (flag_return) return;
    }
    
    // Build map of showers to rays (lines)
    std::map<ShowerPtr, WireCell::Ray> map_shower_ray;
    WireCell::Point main_vtx_pt = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
    
    // Add showers from main vertex
    auto it_main = map_vertex_to_shower.find(main_vertex);
    if (it_main != map_vertex_to_shower.end()) {
        for (auto shower : it_main->second) {
            if (shower->get_particle_type() == 13) continue;
            if (shower->get_total_length() < 3 * units::cm) continue;
            if (pi0_showers.find(shower) != pi0_showers.end()) continue;
            
            WireCell::Point test_p = shower->get_start_point();
            WireCell::Vector dir = shower_cal_dir_3vector(*shower, test_p, 15 * units::cm);
            WireCell::Point p2(test_p.x() + dir.x(), test_p.y() + dir.y(), test_p.z() + dir.z());
            map_shower_ray[shower] = WireCell::Ray(test_p, p2);
        }
    }
    
    // Add showers from other vertices
    for (auto& [vtx, shower_set] : map_vertex_to_shower) {
        if (vtx == main_vertex) continue;
        
        for (auto shower : shower_set) {
            if (shower->get_particle_type() == 13) continue;
            if (shower->get_total_length() < 3 * units::cm) continue;
            if (pi0_showers.find(shower) != pi0_showers.end()) continue;
            
            auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
            if (conn_type != 3) continue;
            
            if (!shower->start_segment()->flags_any(SegmentFlags::kShowerTrajectory) &&
                !shower->start_segment()->flags_any(SegmentFlags::kShowerTopology)) continue;
            
            auto [closest_dis, test_p] = shower_get_closest_point(*shower, main_vtx_pt);
            WireCell::Vector dir = shower_cal_dir_3vector(*shower, test_p, 15 * units::cm);
            WireCell::Point p2(test_p.x() + dir.x(), test_p.y() + dir.y(), test_p.z() + dir.z());
            map_shower_ray[shower] = WireCell::Ray(test_p, p2);
        }
    }
    
    if (map_shower_ray.size() > 1) {
        // Calculate pi0 masses for shower pairs
        std::map<std::pair<ShowerPtr, ShowerPtr>, std::pair<double, WireCell::Point>> map_shower_pair_mass_point;
        
        for (auto it = map_shower_ray.begin(); it != map_shower_ray.end(); it++) {
            ShowerPtr shower_1 = it->first;
            WireCell::Ray ray1 = it->second;
            double length_1 = shower_1->get_total_length();
            
            for (auto it1 = it; it1 != map_shower_ray.end(); it1++) {
                ShowerPtr shower_2 = it1->first;
                if (shower_1 == shower_2) continue;
                
                WireCell::Ray ray2 = it1->second;
                if (ray1.first == ray2.first) continue;
                
                double length_2 = shower_2->get_total_length();
                
                WireCell::Vector dir1 = ray_vector(ray1);
                WireCell::Vector dir2 = ray_vector(ray2);
                double angle_between = std::acos(std::clamp(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()), -1.0, 1.0));
                if (angle_between == 0) continue;
                
                auto [p1_closest, p2_closest] = ray_closest_points(ray1, ray2);
                WireCell::Point center((p1_closest.x() + p2_closest.x()) / 2.0,
                                      (p1_closest.y() + p2_closest.y()) / 2.0,
                                      (p1_closest.z() + p2_closest.z()) / 2.0);
                
                if (length_1 > 15 * units::cm && length_2 > 15 * units::cm) {
                    WireCell::Vector dir_to_shower1(ray1.first.x() - center.x(),
                                                   ray1.first.y() - center.y(),
                                                   ray1.first.z() - center.z());
                    WireCell::Vector dir_to_shower2(ray2.first.x() - center.x(),
                                                   ray2.first.y() - center.y(),
                                                   ray2.first.z() - center.z());
                    
                    if (dir_to_shower1.magnitude() < 3 * units::cm) dir_to_shower1 = dir1;
                    if (dir_to_shower2.magnitude() < 3 * units::cm) dir_to_shower2 = dir2;
                    
                    double angle1 = std::acos(std::clamp(dir_to_shower1.dot(dir1) / (dir_to_shower1.magnitude() * dir1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    double angle2 = std::acos(std::clamp(dir_to_shower2.dot(dir2) / (dir_to_shower2.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                    
                    if (angle1 > 25 || angle2 > 25) continue;
                    
                    double angle = std::acos(std::clamp(dir_to_shower1.dot(dir_to_shower2) / (dir_to_shower1.magnitude() * dir_to_shower2.magnitude()), -1.0, 1.0));
                    double mass_pio = std::sqrt(4 * shower_1->get_kine_charge() * shower_2->get_kine_charge() * 
                                               std::pow(std::sin(angle / 2.0), 2));
                    map_shower_pair_mass_point[std::make_pair(shower_1, shower_2)] = std::make_pair(mass_pio, center);
                    
                } else if (length_1 > 15 * units::cm || length_2 > 15 * units::cm) {
                    WireCell::Vector dir_to_c1, dir_to_c2;
                    
                    if (length_1 > length_2) {
                        center = WireCell::Point((p1_closest.x() + p2_closest.x()) / 2.0,
                                                (p1_closest.y() + p2_closest.y()) / 2.0,
                                                (p1_closest.z() + p2_closest.z()) / 2.0);
                        
                        auto [dis2, test_p] = shower_get_closest_point(*shower_2, center);
                        WireCell::Vector dir3 = shower_cal_dir_3vector(*shower_2, test_p, 15 * units::cm);
                        WireCell::Point p3(test_p.x() + dir3.x(), test_p.y() + dir3.y(), test_p.z() + dir3.z());
                        WireCell::Ray ray3(test_p, p3);
                        
                        auto [new_p1, new_p2] = ray_closest_points(ray1, ray3);
                        center = new_p1;
                        
                        dir_to_c1 = WireCell::Vector(ray1.first.x() - center.x(),
                                                     ray1.first.y() - center.y(),
                                                     ray1.first.z() - center.z());
                        dir_to_c2 = WireCell::Vector(test_p.x() - center.x(),
                                                     test_p.y() - center.y(),
                                                     test_p.z() - center.z());
                        
                        double angle1 = std::acos(std::clamp(dir_to_c1.dot(dir1) / (dir_to_c1.magnitude() * dir1.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        double angle2 = std::acos(std::clamp(dir_to_c2.dot(dir3) / (dir_to_c2.magnitude() * dir3.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        if (angle1 > 25 || angle2 > 25) continue;
                        
                    } else {
                        center = WireCell::Point((p1_closest.x() + p2_closest.x()) / 2.0,
                                                (p1_closest.y() + p2_closest.y()) / 2.0,
                                                (p1_closest.z() + p2_closest.z()) / 2.0);
                        
                        auto [dis1, test_p] = shower_get_closest_point(*shower_1, center);
                        WireCell::Vector dir3 = shower_cal_dir_3vector(*shower_1, test_p, 15 * units::cm);
                        WireCell::Point p3(test_p.x() + dir3.x(), test_p.y() + dir3.y(), test_p.z() + dir3.z());
                        WireCell::Ray ray3(test_p, p3);
                        
                        auto [new_p1, new_p2] = ray_closest_points(ray3, ray2);
                        center = new_p2;
                        
                        dir_to_c2 = WireCell::Vector(ray2.first.x() - center.x(),
                                                     ray2.first.y() - center.y(),
                                                     ray2.first.z() - center.z());
                        dir_to_c1 = WireCell::Vector(test_p.x() - center.x(),
                                                     test_p.y() - center.y(),
                                                     test_p.z() - center.z());
                        
                        double angle1 = std::acos(std::clamp(dir_to_c1.dot(dir3) / (dir_to_c1.magnitude() * dir3.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        double angle2 = std::acos(std::clamp(dir_to_c2.dot(dir2) / (dir_to_c2.magnitude() * dir2.magnitude()), -1.0, 1.0)) / M_PI * 180.0;
                        if (angle1 > 25 || angle2 > 25) continue;
                    }
                    
                    double angle = std::acos(std::clamp(dir_to_c1.dot(dir_to_c2) / (dir_to_c1.magnitude() * dir_to_c2.magnitude()), -1.0, 1.0));
                    double mass_pio = std::sqrt(4 * shower_1->get_kine_charge() * shower_2->get_kine_charge() * 
                                               std::pow(std::sin(angle / 2.0), 2));
                    map_shower_pair_mass_point[std::make_pair(shower_1, shower_2)] = std::make_pair(mass_pio, center);
                    
                } else {
                    break;
                }
            }
        }
        
        // Find best pi0 candidate
        double mass_diff = 1e9;
        double mass_save = 0;
        ShowerPtr shower_1 = nullptr;
        ShowerPtr shower_2 = nullptr;
        double mass_offset = 10 * units::MeV;
        WireCell::Point vtx_point;
        
        for (auto& [shower_pair, mass_point] : map_shower_pair_mass_point) {
            if (std::abs(mass_point.first - 135 * units::MeV + mass_offset) < mass_diff) {
                if (good_showers.find(shower_pair.first) == good_showers.end() && 
                    good_showers.find(shower_pair.second) == good_showers.end()) continue;
                
                shower_1 = shower_pair.first;
                shower_2 = shower_pair.second;
                mass_diff = std::abs(mass_point.first - 135 * units::MeV + mass_offset);
                mass_save = mass_point.first;
                vtx_point = mass_point.second;
            }
        }
        
        // If found good pi0, update everything
        if (mass_diff < 60 * units::MeV && shower_1 && shower_2) {
            pi0_showers.insert(shower_1);
            pi0_showers.insert(shower_2);
            
            int pio_id = acc_segment_id;
            acc_segment_id++;
            
            map_shower_pio_id[shower_1] = pio_id;
            map_shower_pio_id[shower_2] = pio_id;
            map_pio_id_mass[pio_id] = std::make_pair(mass_save, 2);
            map_pio_id_showers[pio_id].push_back(shower_1);
            map_pio_id_showers[pio_id].push_back(shower_2);
            
            // Update main vertex position (hack) - set to reconstructed pi0 decay point
            main_vertex->fit().point = vtx_point;
            main_vertex->fit().dQ = 0;
            
            // Add other segments from main_vertex to showers
            auto [start_vtx_1, conn_type_1] = shower_1->get_start_vertex_and_type();
            if (start_vtx_1 == main_vertex && conn_type_1 == 1) {
                for (auto sg : map_vertex_segments[main_vertex]) {
                    if (sg == shower_1->start_segment()) continue;
                    shower_1->add_segment(sg);
                }
            }
            
            auto [start_vtx_2, conn_type_2] = shower_2->get_start_vertex_and_type();
            if (start_vtx_2 == main_vertex && conn_type_2 == 1) {
                for (auto sg : map_vertex_segments[main_vertex]) {
                    if (sg == shower_2->start_segment()) continue;
                    shower_2->add_segment(sg);
                }
            }
            
            shower_1->set_start_vertex(main_vertex, 2);
            shower_1->calculate_kinematics(particle_data, recomb_model);
            
            shower_2->set_start_vertex(main_vertex, 2);
            shower_2->calculate_kinematics(particle_data, recomb_model);
            
            update_shower_maps(showers, map_vertex_in_shower, map_segment_in_shower,
                              map_vertex_to_shower, used_shower_clusters);
            
            std::cout << "Pi0 (displaced vertex) found with mass: " << mass_save / units::MeV 
                      << " MeV with " << shower_1->get_kine_charge() / units::MeV << " MeV + " 
                      << shower_2->get_kine_charge() / units::MeV << " MeV" << std::endl;
        }
    }
}


void PatternAlgorithms::shower_clustering_with_nv(int acc_segment_id, std::set<ShowerPtr>& pi0_showers, std::map<ShowerPtr, int>& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass,  std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, std::set<VertexPtr>& vertices_in_long_muon, std::set<SegmentPtr>& segments_in_long_muon, Graph& graph, VertexPtr main_vertex, std::set<ShowerPtr>& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, std::map<Facade::Cluster*, VertexPtr> map_cluster_main_vertices,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model){

    
    
    // Connect to the main cluster
    shower_clustering_with_nv_in_main_cluster(graph, main_vertex, showers, 
                                              map_vertex_in_shower, map_segment_in_shower, 
                                              map_vertex_to_shower, used_shower_clusters,
                                              vertices_in_long_muon, segments_in_long_muon);
    
    // Examine things connecting to the main vertex
    shower_clustering_connecting_to_main_vertex(graph, main_vertex, showers,
                                                map_vertex_in_shower, map_segment_in_shower,
                                                map_vertex_to_shower, used_shower_clusters);
    
    // Shower clustering from main cluster
    shower_clustering_with_nv_from_main_cluster(graph, main_vertex, main_cluster, showers,
                                                map_vertex_in_shower, map_segment_in_shower,
                                                map_vertex_to_shower, used_shower_clusters);
    
    // Shower clustering from vertices
    shower_clustering_with_nv_from_vertices(graph, main_vertex, main_cluster, other_clusters, showers,
                                           map_vertex_in_shower, map_segment_in_shower,
                                           map_vertex_to_shower, used_shower_clusters,
                                           vertices_in_long_muon, segments_in_long_muon,
                                           track_fitter, dv, particle_data, recomb_model);
    
    // Calculate shower kinematics
    calculate_shower_kinematics(showers, vertices_in_long_muon, segments_in_long_muon,
                                graph, track_fitter, dv, particle_data, recomb_model);
    
    // Examine and merge showers
    examine_merge_showers(showers, main_vertex, map_vertex_in_shower, map_segment_in_shower,
                         map_vertex_to_shower, used_shower_clusters,
                         vertices_in_long_muon, segments_in_long_muon,
                         graph, track_fitter, dv, particle_data, recomb_model);
    
    // Check remaining clusters
    shower_clustering_in_other_clusters(graph, main_vertex, showers, main_cluster, other_clusters,
                                       map_cluster_main_vertices, map_vertex_in_shower,
                                       map_segment_in_shower, map_vertex_to_shower,
                                       used_shower_clusters, track_fitter, dv,
                                       particle_data, recomb_model, true);
    
    // Calculate shower kinematics again
    calculate_shower_kinematics(showers, vertices_in_long_muon, segments_in_long_muon,
                                graph, track_fitter, dv, particle_data, recomb_model);
    
    // Examine shower trunk and add to shower
    examine_showers(graph, main_vertex, showers, main_cluster, other_clusters,
                   map_cluster_main_vertices, map_vertex_in_shower, map_segment_in_shower,
                   map_vertex_to_shower, used_shower_clusters,
                   track_fitter, dv, particle_data, recomb_model);
    
    // Identify pi0 with vertex
    id_pi0_with_vertex(acc_segment_id, pi0_showers, map_shower_pio_id, map_pio_id_showers, map_pio_id_mass,
                      map_pio_id_saved_pair, graph, main_vertex, showers, main_cluster,
                      other_clusters, map_cluster_main_vertices, map_vertex_in_shower,
                      map_segment_in_shower, map_vertex_to_shower, used_shower_clusters,
                      track_fitter, dv, particle_data, recomb_model);
    
    // Identify pi0 without vertex (displaced vertex)
    id_pi0_without_vertex(acc_segment_id, pi0_showers, map_shower_pio_id, map_pio_id_showers,
                         map_pio_id_mass, map_pio_id_saved_pair, graph, main_vertex, showers,
                         main_cluster, other_clusters, map_cluster_main_vertices,
                         map_vertex_in_shower, map_segment_in_shower, map_vertex_to_shower,
                         used_shower_clusters, track_fitter, dv, particle_data, recomb_model);
}