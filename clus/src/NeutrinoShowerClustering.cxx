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