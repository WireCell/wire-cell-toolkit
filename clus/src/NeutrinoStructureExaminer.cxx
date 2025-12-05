#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

void PatternAlgorithms::examine_structure(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    // Change 2 to 1 (merge two segments into one straight segment)
    if (examine_structure_2(graph, cluster, track_fitter, dv)) {
        track_fitter.do_multi_tracking(true, true, true);
    }
    
    // Straighten 1 (replace curved segments with straight lines)
    if (examine_structure_1(graph, cluster, track_fitter, dv)) {
        track_fitter.do_multi_tracking(true, true, true);
    }
}

bool PatternAlgorithms::examine_structure_1(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    // Look at each segment, if the more straight one is better, change it
    bool flag_update = false;
    
    // Get transform and grouping from track_fitter
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster.get_scope_transform(cluster.get_default_scope()));
    double cluster_t0 = cluster.get_cluster_t0();
    auto grouping = cluster.grouping();
    
    if (!transform || !grouping) {
        return false;
    }
    
    // Iterate through all edges (segments) in the graph
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr sg = graph[*eit].segment;
        
        // Skip if segment doesn't belong to this cluster
        if (!sg || sg->cluster() != &cluster) continue;
        
        // Get segment properties
        double length = segment_track_length(sg);
        double medium_dQ_dx = segment_median_dQ_dx(sg) / (43e3/units::cm);
        
        // Check if segment is short enough and has reasonable dQ/dx
        if (length < 5*units::cm || (length < 8*units::cm && medium_dQ_dx > 1.5)) {
            // Get the two vertices of this segment
            auto [vtx1, vtx2] = find_vertices(graph, sg);
            if (!vtx1 || !vtx2) continue;
            
            const auto& wcpts = sg->wcpts();
            if (wcpts.size() < 2) continue;
            
            // Get start and end points
            Facade::geo_point_t start_p = wcpts.front().point;
            Facade::geo_point_t end_p = wcpts.back().point;
            
            // Check the track by testing points along a straight line
            double step_size = 0.6 * units::cm;
            double distance = std::sqrt(std::pow(start_p.x() - end_p.x(), 2) + 
                                       std::pow(start_p.y() - end_p.y(), 2) + 
                                       std::pow(start_p.z() - end_p.z(), 2));
            int ncount = std::round(distance / step_size);
            
            std::vector<Facade::geo_point_t> new_pts;
            bool flag_replace = true;
            int n_bad = 0;
            
            // Test points along the straight line
            for (int i = 1; i < ncount; i++) {
                Facade::geo_point_t test_p(
                    start_p.x() + (end_p.x() - start_p.x()) / ncount * i,
                    start_p.y() + (end_p.y() - start_p.y()) / ncount * i,
                    start_p.z() + (end_p.z() - start_p.z()) / ncount * i
                );
                new_pts.push_back(test_p);
                
                // Check if this point is good
                auto test_wpid = dv->contained_by(test_p);
                if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                    auto temp_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                    if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*units::cm, 0, 0)) {
                        n_bad++;
                    }
                }
                
                if (n_bad > 1) {
                    flag_replace = false;
                    break;
                }
            }
            
            // If the straight line is better, replace the segment path
            if (flag_replace) {
                // Get steiner point cloud
                const auto& steiner_pc = cluster.get_pc("steiner_pc");
                const auto& coords = cluster.get_default_scope().coords;
                const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
                const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
                const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
                
                // Build new WCPoint list with start point
                std::vector<WCPoint> new_wcpts;
                WCPoint start_wcp = wcpts.front();
                WCPoint end_wcp = wcpts.back();
                
                new_wcpts.push_back(start_wcp);
                
                // Distance threshold for considering points as "same" (0.01 cm)
                const double distance_threshold = 0.01 * units::cm;
                
                // Add intermediate points from steiner cloud
                for (size_t i = 0; i < new_pts.size(); i++) {
                    auto knn_results = cluster.kd_steiner_knn(1, new_pts[i], "steiner_pc");
                    if (!knn_results.empty()) {
                        size_t idx = knn_results[0].first;
                        WCPoint wcp;
                        wcp.point = Facade::geo_point_t(x_coords[idx], y_coords[idx], z_coords[idx]);
                        
                        // Only add if different from last point (using distance comparison)
                        double dist_to_last = ray_length(Ray{wcp.point, new_wcpts.back().point});
                        if (dist_to_last > distance_threshold) {
                            new_wcpts.push_back(wcp);
                        }
                        
                        // Stop if we reached the end point (using distance comparison)
                        double dist_to_end = ray_length(Ray{wcp.point, end_wcp.point});
                        if (dist_to_end < distance_threshold) break;
                    }
                }
                
                // Add end point if not already added (using distance comparison)
                double dist_last_to_end = ray_length(Ray{new_wcpts.back().point, end_wcp.point});
                if (dist_last_to_end > distance_threshold) {
                    new_wcpts.push_back(end_wcp);
                }
                
                // Update the segment with new points
                sg->wcpts(new_wcpts);
                
                flag_update = true;
                std::cout << "Cluster: " << cluster.ident() << " replace Track Content with Straight Line for segment with length " << length/units::cm << " cm" << std::endl;
            }
        }
    }
    
    return flag_update;
}


bool PatternAlgorithms::examine_structure_2(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    bool flag_update = false;
    
    // Get transform and grouping from track_fitter
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster.get_scope_transform(cluster.get_default_scope()));
    double cluster_t0 = cluster.get_cluster_t0();
    auto grouping = cluster.grouping();
    
    if (!transform || !grouping) {
        return false;
    }
    
    bool flag_continue = true;
    while (flag_continue) {
        flag_continue = false;
        
        // Iterate through all vertices in the graph
        auto [vbegin, vend] = boost::vertices(graph);
        for (auto vit = vbegin; vit != vend; ++vit) {
            VertexPtr vtx = graph[*vit].vertex;
            
            // Skip if vertex doesn't belong to this cluster
            if (!vtx || vtx->cluster() != &cluster) continue;
            
            // Check if this vertex has exactly 2 connected segments
            auto vd = vtx->get_descriptor();
            if (boost::degree(vd, graph) != 2) continue;
            
            // Get the two segments connected to this vertex
            auto edge_range = boost::out_edges(vd, graph);
            auto eit = edge_range.first;
            SegmentPtr sg1 = graph[*eit].segment;
            ++eit;
            SegmentPtr sg2 = graph[*eit].segment;
            
            if (!sg1 || !sg2) continue;
            
            // double length1 = segment_track_length(sg1);
            // double length2 = segment_track_length(sg2);
            
            // Get the other vertices (endpoints)
            VertexPtr vtx1 = find_other_vertex(graph, sg1, vtx);
            VertexPtr vtx2 = find_other_vertex(graph, sg2, vtx);
            
            if (!vtx1 || !vtx2) continue;
            
            // Use fitted points if available, otherwise use wcpt points
            Facade::geo_point_t start_p = vtx1->fit().valid() ? vtx1->fit().point : vtx1->wcpt().point;
            Facade::geo_point_t end_p = vtx2->fit().valid() ? vtx2->fit().point : vtx2->wcpt().point;
            
            // Test points along a straight line between the two endpoints
            double step_size = 0.6 * units::cm;
            double distance = std::sqrt(std::pow(start_p.x() - end_p.x(), 2) + 
                                       std::pow(start_p.y() - end_p.y(), 2) + 
                                       std::pow(start_p.z() - end_p.z(), 2));
            int ncount = std::round(distance / step_size);
            
            std::vector<Facade::geo_point_t> new_pts;
            bool flag_replace = true;
            int n_bad = 0;
            
            // Test points along the straight line
            for (int i = 1; i < ncount; i++) {
                Facade::geo_point_t test_p(
                    start_p.x() + (end_p.x() - start_p.x()) / ncount * i,
                    start_p.y() + (end_p.y() - start_p.y()) / ncount * i,
                    start_p.z() + (end_p.z() - start_p.z()) / ncount * i
                );
                new_pts.push_back(test_p);
                
                // Check if this point is good
                auto test_wpid = dv->contained_by(test_p);
                if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                    auto temp_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                    if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*units::cm, 0, 0)) {
                        n_bad++;
                    }
                }
                
                if (n_bad > 1) {
                    flag_replace = false;
                    break;
                }
            }
            
            // If the straight line is better, merge the two segments
            if (flag_replace) {
                std::cout << "Cluster: " << cluster.ident() << " Merge two segments with a straight one, vtx at (" 
                          << vtx->wcpt().point.x()/units::cm << ", "
                          << vtx->wcpt().point.y()/units::cm << ", "
                          << vtx->wcpt().point.z()/units::cm << ") cm" << std::endl;
                
                // Check if the two endpoint vertices are at the same position
                double dist_vtx1_vtx2 = ray_length(Ray{vtx1->wcpt().point, vtx2->wcpt().point});
                const double distance_threshold = 0.01 * units::cm;
                
                if (dist_vtx1_vtx2 < distance_threshold) {
                    // The two endpoint vertices are the same, merge vtx2 into vtx1
                    merge_vertex_into_another(graph, vtx2, vtx1, dv);
                } else {
                    // Create a new segment with straight line path
                    // Get steiner point cloud
                    const auto& steiner_pc = cluster.get_pc("steiner_pc");
                    const auto& coords = cluster.get_default_scope().coords;
                    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
                    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
                    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
                    
                    // Build new WCPoint list with start point
                    std::vector<WCPoint> new_wcpts;
                    WCPoint start_wcp = vtx1->wcpt();
                    WCPoint end_wcp = vtx2->wcpt();
                    
                    new_wcpts.push_back(start_wcp);
                    
                    // Add intermediate points from steiner cloud
                    for (size_t i = 0; i < new_pts.size(); i++) {
                        auto knn_results = cluster.kd_steiner_knn(1, new_pts[i], "steiner_pc");
                        if (!knn_results.empty()) {
                            size_t idx = knn_results[0].first;
                            WCPoint wcp;
                            wcp.point = Facade::geo_point_t(x_coords[idx], y_coords[idx], z_coords[idx]);
                            
                            // Only add if different from last point (using distance comparison)
                            double dist_to_last = ray_length(Ray{wcp.point, new_wcpts.back().point});
                            if (dist_to_last > distance_threshold) {
                                new_wcpts.push_back(wcp);
                            }
                        }
                    }
                    
                    // Add end point if not already added (using distance comparison)
                    double dist_last_to_end = ray_length(Ray{new_wcpts.back().point, end_wcp.point});
                    if (dist_last_to_end > distance_threshold) {
                        new_wcpts.push_back(end_wcp);
                    }
                    
                    // Create new segment with the straight line path
                    auto new_seg = make_segment();
                    new_seg->wcpts(new_wcpts).cluster(&cluster).dirsign(0);
                    
                    // Add the new segment to the graph
                    add_segment(graph, new_seg, vtx1, vtx2);
                }
                
                // Delete the old segments and middle vertex
                remove_segment(graph, sg1);
                remove_segment(graph, sg2);
                remove_vertex(graph, vtx);
                
                flag_update = true;
                flag_continue = true;
                break;
            }
        }
    }
    
    return flag_update;
}

bool PatternAlgorithms::examine_structure_3(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    bool flag_update = false;
    bool flag_continue = true;
    
    while (flag_continue) {
        flag_continue = false;
        
        // Iterate through all vertices in the graph
        auto [vbegin, vend] = boost::vertices(graph);
        for (auto vit = vbegin; vit != vend; ++vit) {
            VertexPtr vtx = graph[*vit].vertex;
            
            // Skip if vertex doesn't belong to this cluster
            if (!vtx || vtx->cluster() != &cluster) continue;
            
            // Check if this vertex has exactly 2 connected segments
            auto vd = vtx->get_descriptor();
            if (boost::degree(vd, graph) != 2) continue;
            
            // Get the two segments connected to this vertex
            auto edge_range = boost::out_edges(vd, graph);
            auto eit = edge_range.first;
            SegmentPtr sg1 = graph[*eit].segment;
            ++eit;
            SegmentPtr sg2 = graph[*eit].segment;
            
            if (!sg1 || !sg2) continue;
            
            // Get the other vertices (endpoints)
            VertexPtr vtx1 = find_other_vertex(graph, sg1, vtx);
            VertexPtr vtx2 = find_other_vertex(graph, sg2, vtx);
            
            if (!vtx1 || !vtx2) continue;
            
            // Get the vertex position (use fit if available, otherwise wcpt)
            WireCell::Point vtx_point = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
            
            // Calculate direction vectors at different distances
            WireCell::Vector dir1 = segment_cal_dir_3vector(sg1, vtx_point, 10*units::cm);
            WireCell::Vector dir2 = segment_cal_dir_3vector(sg2, vtx_point, 10*units::cm);
            
            WireCell::Vector dir3 = segment_cal_dir_3vector(sg1, vtx_point, 3*units::cm);
            WireCell::Vector dir4 = segment_cal_dir_3vector(sg2, vtx_point, 3*units::cm);
            
            // Calculate angles (180 - angle between directions)
            double angle_10cm = (3.1415926 - std::acos(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()))) / 3.1415926 * 180.0;
            double angle_3cm = (3.1415926 - std::acos(dir3.dot(dir4) / (dir3.magnitude() * dir4.magnitude()))) / 3.1415926 * 180.0;
            
            // Check if segments are nearly collinear (small angles)
            if (angle_10cm < 18 && angle_3cm < 27) {
                std::cout << "Cluster: " << cluster.ident() << " Merge two segments into one according to angle " 
                          << angle_10cm << "° (10cm) and " << angle_3cm << "° (3cm)" << std::endl;
                
                // Merge the two segments by combining their WCPoint lists
                const auto& wcpts1 = sg1->wcpts();
                const auto& wcpts2 = sg2->wcpts();
                
                std::vector<WCPoint> merged_wcpts;
                const double distance_threshold = 0.01 * units::cm;
                
                // Determine how to merge based on which endpoints connect
                double dist_front1_front2 = ray_length(Ray{wcpts1.front().point, wcpts2.front().point});
                double dist_front1_back2 = ray_length(Ray{wcpts1.front().point, wcpts2.back().point});
                double dist_back1_front2 = ray_length(Ray{wcpts1.back().point, wcpts2.front().point});
                double dist_back1_back2 = ray_length(Ray{wcpts1.back().point, wcpts2.back().point});
                
                if (dist_front1_front2 < distance_threshold) {
                    // front1 connects to front2: reverse wcpts2, skip duplicate, add wcpts1
                    for (auto it = wcpts2.rbegin(); it != wcpts2.rend(); ++it) {
                        merged_wcpts.push_back(*it);
                    }
                    // Skip first point of wcpts1 as it's the same as last added
                    for (size_t i = 1; i < wcpts1.size(); ++i) {
                        merged_wcpts.push_back(wcpts1[i]);
                    }
                } else if (dist_front1_back2 < distance_threshold) {
                    // front1 connects to back2: add wcpts2, skip duplicate, add wcpts1
                    for (const auto& wcp : wcpts2) {
                        merged_wcpts.push_back(wcp);
                    }
                    // Skip first point of wcpts1 as it's the same as last added
                    for (size_t i = 1; i < wcpts1.size(); ++i) {
                        merged_wcpts.push_back(wcpts1[i]);
                    }
                } else if (dist_back1_front2 < distance_threshold) {
                    // back1 connects to front2: add wcpts1, skip duplicate, add wcpts2
                    for (const auto& wcp : wcpts1) {
                        merged_wcpts.push_back(wcp);
                    }
                    // Skip first point of wcpts2 as it's the same as last added
                    for (size_t i = 1; i < wcpts2.size(); ++i) {
                        merged_wcpts.push_back(wcpts2[i]);
                    }
                } else if (dist_back1_back2 < distance_threshold) {
                    // back1 connects to back2: add wcpts1, skip duplicate, reverse wcpts2
                    for (const auto& wcp : wcpts1) {
                        merged_wcpts.push_back(wcp);
                    }
                    // Skip last point of wcpts2 (reverse order) as it's the same as last added
                    for (auto it = wcpts2.rbegin() + 1; it != wcpts2.rend(); ++it) {
                        merged_wcpts.push_back(*it);
                    }
                }
                
                // Create new segment with merged points
                auto new_seg = make_segment();
                new_seg->wcpts(merged_wcpts).cluster(&cluster).dirsign(0);
                
                // Add the new segment to the graph
                add_segment(graph, new_seg, vtx1, vtx2);
                
                // Delete the old segments and middle vertex
                remove_segment(graph, sg1);
                remove_segment(graph, sg2);
                remove_vertex(graph, vtx);
                
                flag_update = true;
                flag_continue = true;
                break;
            }
        }
    }
    
    return flag_update;
}

bool PatternAlgorithms::examine_structure_4(VertexPtr vertex, bool flag_final_vertex, Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    bool flag_update = false;
    
    // Get transform and grouping from track_fitter
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster.get_scope_transform(cluster.get_default_scope()));
    double cluster_t0 = cluster.get_cluster_t0();
    auto grouping = cluster.grouping();
    
    if (!transform || !grouping) {
        return false;
    }
    
    // Check if vertex belongs to this cluster
    if (!vertex || vertex->cluster() != &cluster) return false;
    
    // Check vertex degree
    auto vd = vertex->get_descriptor();
    int degree = boost::degree(vd, graph);
    if (degree < 2 && !flag_final_vertex) return false;
    
    // Get steiner point cloud and flag terminals
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
    const auto& flag_terminals = steiner_pc.get("flag_steiner_terminal")->elements<int>();
    
    // Get vertex position (use fit if available, otherwise wcpt)
    WireCell::Point vtx_point = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
    
    // Find candidate wcps within 6 cm radius
    auto candidate_results = cluster.kd_steiner_radius(6*units::cm, vtx_point, "steiner_pc");
    
    double max_dis = 0;
    WireCell::Point max_wcp_point;
    bool found_candidate = false;
    
    // Loop over candidate points
    for (const auto& [idx, dist_sq] : candidate_results) {
        // Check if this is a steiner terminal
        if (!flag_terminals[idx]) continue;
        
        WireCell::Point test_p(x_coords[idx], y_coords[idx], z_coords[idx]);
        double dis = std::sqrt(std::pow(test_p.x() - vtx_point.x(), 2) + 
                              std::pow(test_p.y() - vtx_point.y(), 2) + 
                              std::pow(test_p.z() - vtx_point.z(), 2));
        
        // Find minimum distances to all segments
        double min_dis = 1e9;
        double min_dis_u = 1e9;
        double min_dis_v = 1e9;
        double min_dis_w = 1e9;
        
        // Check against all segments in the graph
        auto [ebegin, eend] = boost::edges(graph);
        for (auto eit = ebegin; eit != eend; ++eit) {
            SegmentPtr sg = graph[*eit].segment;
            if (!sg) continue;
            
            // Get 3D closest point distance
            auto [dist_3d, closest_pt] = segment_get_closest_point(sg, test_p, "fit");
            if (dist_3d < min_dis) min_dis = dist_3d;
            
            // Get closest wcpt distance
            const auto& wcpts = sg->wcpts();
            for (const auto& wcp : wcpts) {
                double tmp_dis = std::sqrt(std::pow(wcp.point.x() - test_p.x(), 2) +
                                          std::pow(wcp.point.y() - test_p.y(), 2) +
                                          std::pow(wcp.point.z() - test_p.z(), 2));
                if (tmp_dis < min_dis) min_dis = tmp_dis;
            }
            
            // Get 2D distances - need to determine apa and face
            auto test_wpid = dv->contained_by(test_p);
            if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                auto [dist_u, dist_v, dist_w] = segment_get_closest_2d_distances(sg, test_p, test_wpid.apa(), test_wpid.face(), "fit");
                if (dist_u < min_dis_u) min_dis_u = dist_u;
                if (dist_v < min_dis_v) min_dis_v = dist_v;
                if (dist_w < min_dis_w) min_dis_w = dist_w;
            }
        }
        
        // Check distance criteria
        if (min_dis > 0.9*units::cm && 
            min_dis_u + min_dis_v + min_dis_w > 1.8*units::cm &&
            ((min_dis_u > 0.8*units::cm && min_dis_v > 0.8*units::cm) ||
             (min_dis_u > 0.8*units::cm && min_dis_w > 0.8*units::cm) ||
             (min_dis_v > 0.8*units::cm && min_dis_w > 0.8*units::cm))) {
            
            // Test points along straight line for good points
            double step_size = 0.3 * units::cm;
            WireCell::Point start_p = vtx_point;
            WireCell::Point end_p = test_p;
            double distance = std::sqrt(std::pow(start_p.x() - end_p.x(), 2) + 
                                       std::pow(start_p.y() - end_p.y(), 2) + 
                                       std::pow(start_p.z() - end_p.z(), 2));
            int ncount = std::round(distance / step_size);
            
            bool flag_pass = true;
            int n_bad = 0;
            
            for (int j = 1; j < ncount; j++) {
                WireCell::Point test_p1(
                    start_p.x() + (end_p.x() - start_p.x()) / ncount * j,
                    start_p.y() + (end_p.y() - start_p.y()) / ncount * j,
                    start_p.z() + (end_p.z() - start_p.z()) / ncount * j
                );
                
                auto test_wpid = dv->contained_by(test_p1);
                if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                    auto temp_p_raw = transform->backward(test_p1, cluster_t0, test_wpid.face(), test_wpid.apa());
                    if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.3*units::cm, 0, 0)) {
                        n_bad++;
                    }
                }
                
                if (n_bad > 0) {
                    flag_pass = false;
                    break;
                }
            }
            
            if (flag_pass) {
                if (max_dis < dis) {
                    max_dis = dis;
                    max_wcp_point = test_p;
                    max_wcp_idx = idx;
                    found_candidate = true;
                }
            }
        }
    }
    
    // If we found a good candidate, create a new segment
    if (found_candidate && max_dis > 1.6*units::cm) {
        // Create new vertex at the terminal point
        VertexPtr v1 = make_vertex(graph);
        WCPoint new_wcp;
        new_wcp.point = max_wcp_point;
        v1->wcpt(new_wcp);
        v1->cluster(&cluster);
        
        // Build wcpoint list for the new segment
        std::vector<WCPoint> wcp_list;
        wcp_list.push_back(vertex->wcpt());
        
        // Add intermediate points with 1 cm steps
        double step_size = 1.0 * units::cm;
        WireCell::Point start_p = vtx_point;
        WireCell::Point end_p = max_wcp_point;
        double distance = std::sqrt(std::pow(start_p.x() - end_p.x(), 2) + 
                                   std::pow(start_p.y() - end_p.y(), 2) + 
                                   std::pow(start_p.z() - end_p.z(), 2));
        int ncount = std::round(distance / step_size);
        
        const double distance_threshold = 0.01 * units::cm;
        
        for (int j = 1; j < ncount; j++) {
            WireCell::Point tmp_p(
                start_p.x() + (end_p.x() - start_p.x()) / ncount * j,
                start_p.y() + (end_p.y() - start_p.y()) / ncount * j,
                start_p.z() + (end_p.z() - start_p.z()) / ncount * j
            );
            
            auto knn_results = cluster.kd_steiner_knn(1, tmp_p, "steiner_pc");
            if (!knn_results.empty()) {
                size_t idx = knn_results[0].first;
                WCPoint wcp;
                wcp.point = WireCell::Point(x_coords[idx], y_coords[idx], z_coords[idx]);
                
                // Only add if different from last point
                double dist_to_last = ray_length(Ray{wcp.point, wcp_list.back().point});
                if (dist_to_last > distance_threshold) {
                    wcp_list.push_back(wcp);
                }
            }
        }
        
        // Add end point if not already added
        WCPoint end_wcp;
        end_wcp.point = max_wcp_point;
        double dist_last_to_end = ray_length(Ray{wcp_list.back().point, end_wcp.point});
        if (dist_last_to_end > distance_threshold) {
            wcp_list.push_back(end_wcp);
        }
        
        std::cout << "Cluster: " << cluster.ident() << " Add a track to the main vertex, " 
                  << wcp_list.size() << " points" << std::endl;
        
        // Create new segment
        auto sg1 = make_segment();
        sg1->wcpts(wcp_list).cluster(&cluster).dirsign(0);
        
        // Add segment to graph
        add_segment(graph, sg1, v1, vertex);
        
        flag_update = true;
    }
    
    return flag_update;
}