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

void PatternAlgorithms::deghost_clusters(Graph& graph, std::vector<Facade::Cluster*>& all_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    // Order clusters by total segment length
    std::map<Facade::Cluster*, std::vector<SegmentPtr>> map_cluster_to_segments;
    std::map<Facade::Cluster*, double> map_cluster_total_length;
    std::vector<Facade::Cluster*> ordered_clusters;
    order_clusters(graph, ordered_clusters, map_cluster_to_segments, map_cluster_total_length);
    
    if (ordered_clusters.empty()) return;
    
    // Get first cluster's grouping to access wpids and dead channel info
    auto* first_grouping = ordered_clusters[0]->grouping();
    if (!first_grouping) return;
    
    // Build wpid_params
    const auto& wpids = first_grouping->wpids();
    std::map<WirePlaneId, std::tuple<Facade::geo_point_t, double, double, double>> wpid_params;
    std::map<WirePlaneId, std::pair<Facade::geo_point_t, double>> wpid_U_dir, wpid_V_dir, wpid_W_dir;
    std::set<int> apas;
    Facade::compute_wireplane_params(wpids, dv, wpid_params, wpid_U_dir, wpid_V_dir, wpid_W_dir, apas);
    
    // Get dead channel maps
    std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_u_index; 
    std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_v_index; 
    std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_w_index;
    
    for (const auto& wpid : wpids) {
        int apa = wpid.apa();
        int face = wpid.face();
        af_dead_u_index[apa][face] = first_grouping->get_dead_winds(apa, face, 0);
        af_dead_v_index[apa][face] = first_grouping->get_dead_winds(apa, face, 1);
        af_dead_w_index[apa][face] = first_grouping->get_dead_winds(apa, face, 2);
    }
    
    // Create global point clouds
    auto global_point_cloud = std::make_shared<Facade::DynamicPointCloud>(wpid_params);
    auto global_steiner_point_cloud = std::make_shared<Facade::DynamicPointCloud>(wpid_params);
    auto global_skeleton_cloud = std::make_shared<Facade::DynamicPointCloud>(wpid_params);
    
    // Add points from clusters not in ordered list
    for (auto cluster : all_clusters) {
        if (map_cluster_total_length.find(cluster) == map_cluster_total_length.end()) {
            global_point_cloud->add_points(Facade::make_points_cluster(cluster, wpid_params, true));
        }
    }
    
    std::vector<Facade::Cluster*> to_be_removed_clusters;
    
    // Process ordered clusters
    for (size_t i = 0; i < ordered_clusters.size(); i++) {
        if (i == 0) {
            // First cluster: add all its points
            global_point_cloud->add_points(Facade::make_points_cluster(ordered_clusters[i], wpid_params, true));
            global_steiner_point_cloud->add_points(Facade::make_points_cluster_steiner(ordered_clusters[i], wpid_params, true));
            
            // Add skeleton points from segments
            auto it = map_cluster_to_segments.find(ordered_clusters[i]);
            if (it != map_cluster_to_segments.end()) {
                for (auto seg : it->second) {
                    // Create point-plane pairs from segment fits
                    std::vector<std::pair<Facade::geo_point_t, WirePlaneId>> point_plane_pairs;
                    for (const auto& fit : seg->fits()) {
                        WirePlaneId wpid = dv->contained_by(fit.point);
                        point_plane_pairs.emplace_back(fit.point, wpid);
                    }
                    global_skeleton_cloud->add_points(Facade::make_points_direct(ordered_clusters[i], dv, wpid_params, point_plane_pairs, true));
                }
            }
        } else {
            // Process subsequent clusters
            Facade::Cluster* cluster = ordered_clusters[i];
            int num_dead[3] = {0, 0, 0};
            int num_unique[3] = {0, 0, 0};
            int num_total_points = 0;
            
            double dis_cut = 1.2 * units::cm;
            
            auto it = map_cluster_to_segments.find(cluster);
            if (it != map_cluster_to_segments.end()) {
                for (auto seg : it->second) {
                    for (const auto& fit : seg->fits()) {
                        Facade::geo_point_t test_point = fit.point;
                        num_total_points++;
                        
                        WirePlaneId test_wpid = dv->contained_by(test_point);
                        int apa = test_wpid.apa();
                        int face = test_wpid.face();
                        
                        // Get point in raw coordinates for dead channel check
                        auto transform = track_fitter.get_pc_transforms()->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
                        double cluster_t0 = cluster->get_cluster_t0();
                        Facade::geo_point_t p_raw = transform->backward(test_point, cluster_t0, face, apa);
                        
                        // const auto& winds = cluster->wire_indices();
                        
                        // U plane
                        bool flag_dead = false;
                        if (af_dead_u_index.find(apa) != af_dead_u_index.end() && 
                            af_dead_u_index[apa].find(face) != af_dead_u_index[apa].end()) {
                            flag_dead = cluster->grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0);
                        }
                        
                        if (!flag_dead) {
                            auto results = global_point_cloud->get_closest_2d_point_info(test_point, 0, face, apa);
                            if (std::get<0>(results) <= dis_cut / 2.) {
                                // Overlap with global point cloud
                            } else {
                                results = global_steiner_point_cloud->get_closest_2d_point_info(test_point, 0, face, apa);
                                if (std::get<0>(results) <= dis_cut * 2. / 3.) {
                                    // Overlap with steiner cloud
                                } else {
                                    results = global_skeleton_cloud->get_closest_2d_point_info(test_point, 0, face, apa);
                                    if (std::get<0>(results) <= dis_cut * 6. / 4.) {
                                        // Overlap with skeleton cloud
                                    } else {
                                        num_unique[0]++;
                                    }
                                }
                            }
                        } else {
                            num_dead[0]++;
                        }
                        
                        // V plane
                        flag_dead = false;
                        if (af_dead_v_index.find(apa) != af_dead_v_index.end() && 
                            af_dead_v_index[apa].find(face) != af_dead_v_index[apa].end()) {
                            flag_dead = cluster->grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1);
                        }
                        
                        if (!flag_dead) {
                            auto results = global_point_cloud->get_closest_2d_point_info(test_point, 1, face, apa);
                            if (std::get<0>(results) <= dis_cut / 2.) {
                            } else {
                                results = global_steiner_point_cloud->get_closest_2d_point_info(test_point, 1, face, apa);
                                if (std::get<0>(results) <= dis_cut * 2. / 3.) {
                                } else {
                                    results = global_skeleton_cloud->get_closest_2d_point_info(test_point, 1, face, apa);
                                    if (std::get<0>(results) <= dis_cut * 6. / 4.) {
                                    } else {
                                        num_unique[1]++;
                                    }
                                }
                            }
                        } else {
                            num_dead[1]++;
                        }
                        
                        // W plane
                        flag_dead = false;
                        if (af_dead_w_index.find(apa) != af_dead_w_index.end() && 
                            af_dead_w_index[apa].find(face) != af_dead_w_index[apa].end()) {
                            flag_dead = cluster->grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2);
                        }
                        
                        if (!flag_dead) {
                            auto results = global_point_cloud->get_closest_2d_point_info(test_point, 2, face, apa);
                            if (std::get<0>(results) <= dis_cut / 2.) {
                            } else {
                                results = global_steiner_point_cloud->get_closest_2d_point_info(test_point, 2, face, apa);
                                if (std::get<0>(results) <= dis_cut * 2. / 3.) {
                                } else {
                                    results = global_skeleton_cloud->get_closest_2d_point_info(test_point, 2, face, apa);
                                    if (std::get<0>(results) <= dis_cut * 6. / 4.) {
                                    } else {
                                        num_unique[2]++;
                                    }
                                }
                            }
                        } else {
                            num_dead[2]++;
                        }
                    }
                }
            }
            
            // Calculate percentages
            bool flag_add = true;
            if (num_total_points > 0) {
                double unique_percent_u = num_unique[0] * 1.0 / num_total_points;
                double unique_percent_v = num_unique[1] * 1.0 / num_total_points;
                double unique_percent_w = num_unique[2] * 1.0 / num_total_points;
                
                double dead_percent_u = num_dead[0] * 1.0 / num_total_points;
                double dead_percent_v = num_dead[1] * 1.0 / num_total_points;
                double dead_percent_w = num_dead[2] * 1.0 / num_total_points;
                
                double max_unique_percent = std::max({unique_percent_u, unique_percent_v, unique_percent_w});
                double min_unique_percent = std::min({unique_percent_u, unique_percent_v, unique_percent_w});
                double ave_unique_percent = (unique_percent_u + unique_percent_v + unique_percent_w) / 3.;
                double max_dead_percent = std::max({dead_percent_u, dead_percent_v, dead_percent_w});
                
                // Apply ghosting criteria
                if ((max_dead_percent >= 0.8 && max_unique_percent <= 0.35 && ave_unique_percent <= 0.16 && min_unique_percent <= 0.08) ||
                    (max_unique_percent <= 0.1 && ave_unique_percent <= 0.05 && min_unique_percent <= 0.025) ||
                    (max_dead_percent < 0.8 && max_dead_percent >= 0.7 && max_unique_percent <= 0.2 && ave_unique_percent <= 0.1 && min_unique_percent <= 0.05)) {
                    flag_add = false;
                }
                
                // Additional check for one dead plane
                if ((num_dead[0] == num_total_points || num_dead[1] == num_total_points || num_dead[2] == num_total_points) &&
                    ((num_unique[0] == 0 && num_unique[1] == 0) || (num_unique[0] == 0 && num_unique[2] == 0) || (num_unique[2] == 0 && num_unique[1] == 0)) &&
                    flag_add && max_unique_percent < 0.75) {
                    flag_add = false;
                }
            }
            
            if (flag_add) {
                // Add to global clouds
                global_point_cloud->add_points(Facade::make_points_cluster(cluster, wpid_params, true));
                global_steiner_point_cloud->add_points(Facade::make_points_cluster_steiner(cluster, wpid_params, true));
                
                auto it = map_cluster_to_segments.find(cluster);
                if (it != map_cluster_to_segments.end()) {
                    for (auto seg : it->second) {
                        std::vector<std::pair<Facade::geo_point_t, WirePlaneId>> point_plane_pairs;
                        for (const auto& fit : seg->fits()) {
                            WirePlaneId wpid = dv->contained_by(fit.point);
                            point_plane_pairs.emplace_back(fit.point, wpid);
                        }
                        global_skeleton_cloud->add_points(Facade::make_points_direct(cluster, dv, wpid_params, point_plane_pairs, true));
                    }
                }
            } else {
                to_be_removed_clusters.push_back(cluster);
            }
        }
    }
    
    // Remove segments from ghosted clusters
    for (auto cluster : to_be_removed_clusters) {
        auto it = map_cluster_to_segments.find(cluster);
        if (it != map_cluster_to_segments.end()) {
            for (auto seg : it->second) {
                remove_segment(graph, seg);
            }
        }
    }
    
    // Clean up orphaned vertices
    std::vector<VertexPtr> tmp_vertices;
    auto [vbegin, vend] = boost::vertices(graph);
    for (auto vit = vbegin; vit != vend; ++vit) {
        VertexPtr vtx = graph[*vit].vertex;
        if (vtx && boost::out_degree(*vit, graph) == 0) {
            tmp_vertices.push_back(vtx);
        }
    }
    
    for (auto vtx : tmp_vertices) {
        remove_vertex(graph, vtx);
    }
}

