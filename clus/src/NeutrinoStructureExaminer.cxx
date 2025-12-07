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
                    // max_wcp_idx = idx;
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


bool PatternAlgorithms::crawl_segment(Graph& graph, Facade::Cluster& cluster, SegmentPtr seg, VertexPtr vertex, TrackFitting& track_fitter, IDetectorVolumes::pointer dv ){
    bool flag = false;
    
    // Validate that segment, vertex, and cluster all match
    if (!seg || !vertex || seg->cluster() != &cluster || vertex->cluster() != &cluster) {
        return flag;
    }
    
    // Step 1: Find points at ~3cm distance from vertex on other connected segments
    std::map<SegmentPtr, Facade::geo_point_t> map_segment_point;
    
    if (!vertex->descriptor_valid()) return flag;
    auto vd = vertex->get_descriptor();
    auto edge_range = boost::out_edges(vd, graph);
    
    for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
        SegmentPtr sg = graph[*eit].segment;
        if (!sg || sg == seg) continue;
        
        const auto& fits = sg->fits();
        if (fits.empty()) continue;
        
        Facade::geo_point_t min_point = fits.front().point;
        double min_dis = 1e9;
        
        for (size_t i = 0; i < fits.size(); i++) {
            double dis = std::fabs(ray_length(Ray{fits[i].point, vertex->fit().point}) - 3.0 * units::cm);
            if (dis < min_dis) {
                min_dis = dis;
                min_point = fits[i].point;
            }
        }
        map_segment_point[sg] = min_point;
    }
    
    // Step 2: Determine which end of seg connects to vertex
    const auto& seg_wcpts = seg->wcpts();
    const auto& seg_fits = seg->fits();
    if (seg_wcpts.size() < 2) return flag;
    
    bool flag_start = false;
    double dis_front = ray_length(Ray{vertex->wcpt().point, seg_wcpts.front().point});
    double dis_back = ray_length(Ray{vertex->wcpt().point, seg_wcpts.back().point});
    
    if (dis_front < dis_back) {
        flag_start = true;
    }
    
    // Step 3: Build list of points to test (from vertex end, excluding endpoints)
    std::vector<Facade::geo_point_t> pts_to_be_tested;
    
    if (flag_start) {
        for (size_t i = 1; i + 1 < seg_fits.size(); i++) {
            pts_to_be_tested.push_back(seg_fits[i].point);
        }
    } else {
        for (int i = int(seg_fits.size()) - 2; i > 0; i--) {
            pts_to_be_tested.push_back(seg_fits[i].point);
        }
    }
    
    if (pts_to_be_tested.empty()) return flag;
    
    // Step 4: Test points for good connectivity
    double step_size = 0.3 * units::cm;
    int max_bin = -1;
    
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster.get_scope_transform(cluster.get_default_scope()));
    double cluster_t0 = cluster.get_cluster_t0();
    
    for (size_t i = 0; i < pts_to_be_tested.size(); i++) {
        int n_bad = 0;
        Facade::geo_point_t end_p = pts_to_be_tested[i];
        
        for (auto it = map_segment_point.begin(); it != map_segment_point.end(); it++) {
            Facade::geo_point_t start_p = it->second;
            double distance = ray_length(Ray{start_p, end_p});
            int ncount = std::round(distance / step_size);
            
            for (int j = 1; j < ncount; j++) {
                Facade::geo_point_t test_p(
                    start_p.x() + (end_p.x() - start_p.x()) / ncount * j,
                    start_p.y() + (end_p.y() - start_p.y()) / ncount * j,
                    start_p.z() + (end_p.z() - start_p.z()) / ncount * j
                );
                
                // Check if point is good
                auto test_wpid = dv->contained_by(test_p);
                if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                    auto temp_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                    if (!cluster.grouping()->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2 * units::cm, 0, 0)) {
                        n_bad++;
                    }
                }
            }
        }
        
        if (n_bad == 0) max_bin = i;
    }
    
    // Step 5: Update segment and vertex if good point found
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
    
    while (max_bin >= 0) {
        // Find closest steiner point to the test point
        auto vtx_new_knn = cluster.kd_steiner_knn(1, pts_to_be_tested[max_bin], "steiner_pc");
        if (vtx_new_knn.empty()) {
            max_bin--;
            continue;
        }
        
        size_t new_idx = vtx_new_knn[0].first;
        Facade::geo_point_t vtx_new_point(x_coords[new_idx], y_coords[new_idx], z_coords[new_idx]);
        
        // Check if new point is valid (not the same as current endpoints)
        if (flag_start && ray_length(Ray{vtx_new_point, seg_wcpts.back().point}) < 0.01 * units::cm) {
            max_bin--;
            continue;
        }
        if (!flag_start && ray_length(Ray{vtx_new_point, seg_wcpts.front().point}) < 0.01 * units::cm) {
            max_bin--;
            continue;
        }
        if (ray_length(Ray{vtx_new_point, vertex->wcpt().point}) < 0.01 * units::cm) {
            break;
        }
        
        // Update current segment
        std::vector<Facade::geo_point_t> new_path;
        if (flag_start) {
            // Keep points from back to new vertex
            double dis_limit = ray_length(Ray{vtx_new_point, seg_wcpts.back().point});
            for (int idx = seg_wcpts.size() - 1; idx >= 0; idx--) {
                double dis = ray_length(Ray{seg_wcpts[idx].point, seg_wcpts.back().point});
                if (dis < dis_limit) {
                    new_path.push_back(seg_wcpts[idx].point);
                }
            }
            std::reverse(new_path.begin(), new_path.end());
            if (new_path.size() > 1 && ray_length(Ray{new_path.front(), vtx_new_point}) < 0.01 * units::cm) {
                new_path.erase(new_path.begin());
            }
            new_path.insert(new_path.begin(), vtx_new_point);
        } else {
            // Keep points from front to new vertex
            double dis_limit = ray_length(Ray{vtx_new_point, seg_wcpts.front().point});
            for (size_t idx = 0; idx < seg_wcpts.size(); idx++) {
                double dis = ray_length(Ray{seg_wcpts[idx].point, seg_wcpts.front().point});
                if (dis < dis_limit) {
                    new_path.push_back(seg_wcpts[idx].point);
                }
            }
            if (new_path.size() > 1 && ray_length(Ray{new_path.back(), vtx_new_point}) < 0.01 * units::cm) {
                new_path.pop_back();
            }
            new_path.push_back(vtx_new_point);
        }
        
        // Replace segment with updated path
        auto other_vertex = find_other_vertex(graph, seg, vertex);
        if (!other_vertex) break;
        
        remove_segment(graph, seg);
        auto new_seg = create_segment_for_cluster(cluster, dv, new_path, 0);
        if (!new_seg) break;
        
        add_segment(graph, new_seg, flag_start ? other_vertex : vertex, 
                                     flag_start ? vertex : other_vertex);
        seg = new_seg;
        
        // Update other connected segments
        for (auto it = map_segment_point.begin(); it != map_segment_point.end(); it++) {
            SegmentPtr other_sg = it->first;
            const auto& other_wcpts = other_sg->wcpts();
            if (other_wcpts.empty()) continue;
            
            bool flag_front = (ray_length(Ray{other_wcpts.front().point, vertex->wcpt().point}) < 
                              ray_length(Ray{other_wcpts.back().point, vertex->wcpt().point}));
            Facade::geo_point_t min_p = it->second;
            
            // Find closest point in other segment to min_p
            size_t min_idx = 0;
            double min_dis = 1e9;
            for (size_t j = 0; j < other_wcpts.size(); j++) {
                double dis = ray_length(Ray{min_p, other_wcpts[j].point});
                if (dis < min_dis) {
                    min_dis = dis;
                    min_idx = j;
                }
            }
            
            // Build new path from vtx_new_point to min point
            Facade::geo_point_t min_wcpt_point = other_wcpts[min_idx].point;
            auto path_points = do_rough_path(cluster, vtx_new_point, min_wcpt_point);
            
            // Combine with rest of segment
            std::vector<Facade::geo_point_t> combined_path;
            if (flag_front) {
                combined_path = path_points;
                for (size_t j = min_idx + 1; j < other_wcpts.size(); j++) {
                    combined_path.push_back(other_wcpts[j].point);
                }
            } else {
                for (size_t j = 0; j < min_idx; j++) {
                    combined_path.push_back(other_wcpts[j].point);
                }
                std::reverse(path_points.begin(), path_points.end());
                combined_path.insert(combined_path.end(), path_points.begin(), path_points.end());
            }
            
            if (combined_path.size() <= 1) continue;
            
            // Replace other segment
            auto other_v2 = find_other_vertex(graph, other_sg, vertex);
            if (!other_v2) continue;
            
            remove_segment(graph, other_sg);
            auto new_other_seg = create_segment_for_cluster(cluster, dv, combined_path, 0);
            if (!new_other_seg) continue;
            
            add_segment(graph, new_other_seg, flag_front ? vertex : other_v2,
                                              flag_front ? other_v2 : vertex);
        }
        
        // Update vertex position
        vertex->wcpt().point = vtx_new_point;
        if (vertex->fit().valid()) {
            vertex->fit().point = vtx_new_point;
        }
        
        // Perform multi-tracking
        track_fitter.do_multi_tracking(true, true, true);
        
        flag = true;
        break;
    }
    
    return flag;
}

void PatternAlgorithms::examine_segment(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    // Step 1: Examine short segments with multiple connections at both ends
    auto [ebegin, eend] = boost::edges(graph);
    std::vector<SegmentPtr> segments_to_examine;
    
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr sg = graph[*eit].segment;
        if (!sg || sg->cluster() != &cluster) continue;
        
        double length = segment_track_length(sg);
        if (length > 4 * units::cm) continue;
        
        auto [v1, v2] = find_vertices(graph, sg);
        if (!v1 || !v2) continue;
        
        // Check both vertices have at least 2 connections
        auto v1d = v1->get_descriptor();
        auto v2d = v2->get_descriptor();
        if (boost::degree(v1d, graph) < 2 || boost::degree(v2d, graph) < 2) continue;
        
        segments_to_examine.push_back(sg);
    }
    
    // Examine each short segment for potential crawling
    for (auto sg : segments_to_examine) {
        auto [v1, v2] = find_vertices(graph, sg);
        if (!v1 || !v2) continue;
        
        std::vector<VertexPtr> cand_vertices = {v1, v2};
        
        for (size_t i = 0; i < 2; i++) {
            VertexPtr vtx = cand_vertices[i];
            if (!vtx->descriptor_valid()) continue;
            
            // Calculate direction of current segment at this vertex
            Facade::geo_point_t vtx_point = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
            auto dir1 = segment_cal_dir_3vector(sg, vtx_point, 2 * units::cm);
            
            double max_angle = 0;
            double min_angle = 180;
            
            // Check angles with other connected segments
            auto vd = vtx->get_descriptor();
            auto edge_range = boost::out_edges(vd, graph);
            
            for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
                SegmentPtr sg1 = graph[*eit].segment;
                if (!sg1 || sg1 == sg) continue;
                
                auto dir3 = segment_cal_dir_3vector(sg1, vtx_point, 2 * units::cm);
                
                // Calculate angle between directions
                double dot_product = dir1.dot(dir3);
                double mag1 = dir1.magnitude();
                double mag3 = dir3.magnitude();
                
                if (mag1 > 0 && mag3 > 0) {
                    double cos_angle = dot_product / (mag1 * mag3);
                    // Clamp to [-1, 1] to handle numerical errors
                    if (cos_angle > 1.0) cos_angle = 1.0;
                    if (cos_angle < -1.0) cos_angle = -1.0;
                    double angle = std::acos(cos_angle) / 3.1415926 * 180.0;
                    
                    if (angle > max_angle) max_angle = angle;
                    if (angle < min_angle) min_angle = angle;
                }
            }
            
            // If angles indicate a sharp turn, try to crawl the segment
            if (max_angle > 150 && min_angle > 105) {
                crawl_segment(graph, cluster, sg, vtx, track_fitter, dv);
            }
        }
    }
    
    // Step 2: Merge vertices at the same position
    std::set<VertexPtr> all_vertices;
    auto [vbegin, vend] = boost::vertices(graph);
    for (auto vit = vbegin; vit != vend; ++vit) {
        VertexPtr vtx = graph[*vit].vertex;
        if (vtx && vtx->cluster() == &cluster) {
            all_vertices.insert(vtx);
        }
    }
    
    bool flag_merge = true;
    while (flag_merge) {
        flag_merge = false;
        VertexPtr vtx1 = nullptr;
        VertexPtr vtx2 = nullptr;
        
        for (auto it1 = all_vertices.begin(); it1 != all_vertices.end(); ++it1) {
            vtx1 = *it1;
            for (auto it2 = it1; it2 != all_vertices.end(); ++it2) {
                vtx2 = *it2;
                
                // Check if two different vertices are at the same position
                if (vtx1 != vtx2 && 
                    ray_length(Ray{vtx1->wcpt().point, vtx2->wcpt().point}) < 0.01 * units::cm) {
                    
                    // Find segments to remove (connected to both vertices)
                    std::vector<SegmentPtr> to_be_removed_segments;
                    
                    if (vtx2->descriptor_valid()) {
                        auto v2d = vtx2->get_descriptor();
                        auto edge_range = boost::out_edges(v2d, graph);
                        
                        for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
                            SegmentPtr sg = graph[*eit].segment;
                            if (!sg) continue;
                            
                            // Check if this segment is also connected to vtx1
                            auto [seg_v1, seg_v2] = find_vertices(graph, sg);
                            if ((seg_v1 == vtx1 || seg_v2 == vtx1) && 
                                (seg_v1 == vtx2 || seg_v2 == vtx2)) {
                                to_be_removed_segments.push_back(sg);
                            }
                        }
                    }
                    
                    // Merge vtx2 into vtx1
                    if (merge_vertex_into_another(graph, vtx2, vtx1, dv)) {
                        // Remove duplicate segments
                        for (auto sg : to_be_removed_segments) {
                            remove_segment(graph, sg);
                        }
                        
                        all_vertices.erase(vtx2);
                        flag_merge = true;
                        break;
                    }
                }
            }
            if (flag_merge) break;
        }
    }
    
    // Step 3: Remove duplicate segments (same endpoints)
    std::set<SegmentPtr> segments_to_be_removed;
    
    for (auto it = all_vertices.begin(); it != all_vertices.end(); ++it) {
        VertexPtr vtx = *it;
        if (!vtx->descriptor_valid()) continue;
        
        auto vd = vtx->get_descriptor();
        auto edge_range = boost::out_edges(vd, graph);
        
        std::vector<SegmentPtr> tmp_segments;
        for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
            SegmentPtr sg = graph[*eit].segment;
            if (sg) tmp_segments.push_back(sg);
        }
        
        // Compare all pairs of segments
        for (size_t i = 0; i < tmp_segments.size(); i++) {
            auto [v1_i, v2_i] = find_vertices(graph, tmp_segments[i]);
            if (!v1_i || !v2_i) continue;
            
            for (size_t j = i + 1; j < tmp_segments.size(); j++) {
                auto [v1_j, v2_j] = find_vertices(graph, tmp_segments[j]);
                if (!v1_j || !v2_j) continue;
                
                // Check if segments have the same endpoints
                bool same_endpoints = false;
                
                double dis_v1i_v1j = ray_length(Ray{v1_i->wcpt().point, v1_j->wcpt().point});
                double dis_v1i_v2j = ray_length(Ray{v1_i->wcpt().point, v2_j->wcpt().point});
                double dis_v2i_v1j = ray_length(Ray{v2_i->wcpt().point, v1_j->wcpt().point});
                double dis_v2i_v2j = ray_length(Ray{v2_i->wcpt().point, v2_j->wcpt().point});
                
                if ((dis_v1i_v1j < 0.01 * units::cm && dis_v2i_v2j < 0.01 * units::cm) ||
                    (dis_v1i_v2j < 0.01 * units::cm && dis_v2i_v1j < 0.01 * units::cm)) {
                    same_endpoints = true;
                }
                
                if (same_endpoints) {
                    segments_to_be_removed.insert(tmp_segments[j]);
                }
            }
        }
    }
    
    // Remove duplicate segments
    for (auto sg : segments_to_be_removed) {
        remove_segment(graph, sg);
    }
    
    // Step 4: Remove isolated vertices (no connections)
    std::set<VertexPtr> vertices_to_be_removed;
    auto [vbegin2, vend2] = boost::vertices(graph);
    
    for (auto vit = vbegin2; vit != vend2; ++vit) {
        VertexPtr vtx = graph[*vit].vertex;
        if (vtx && vtx->cluster() == &cluster) {
            if (vtx->descriptor_valid()) {
                auto vd = vtx->get_descriptor();
                if (boost::degree(vd, graph) == 0) {
                    vertices_to_be_removed.insert(vtx);
                }
            }
        }
    }
    
    for (auto vtx : vertices_to_be_removed) {
        remove_vertex(graph, vtx);
    }
}


bool PatternAlgorithms::examine_vertices_1p(Graph&graph, VertexPtr v1, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    if (!v1 || !v2 || !v1->cluster() || v1->cluster() != v2->cluster()) {
        return false;
    }
    
    auto& cluster = *v1->cluster();
    
    // Find the segment between v1 and v2
    SegmentPtr sg = find_segment(graph, v1, v2);
    if (!sg) return false;
    
    // Check that v1 has exactly 2 connections
    if (!v1->descriptor_valid()) return false;
    auto v1d = v1->get_descriptor();
    if (boost::degree(v1d, graph) != 2) return false;
    
    // Get vertex positions
    Facade::geo_point_t v1_p = v1->fit().valid() ? v1->fit().point : v1->wcpt().point;
    Facade::geo_point_t v2_p = v2->fit().valid() ? v2->fit().point : v2->wcpt().point;
    
    // Get wpid for coordinate conversion
    auto v1_wpid = dv->contained_by(v1_p);
    auto v2_wpid = dv->contained_by(v2_p);
    if (v1_wpid.face() == -1 || v1_wpid.apa() == -1 || 
        v2_wpid.face() == -1 || v2_wpid.apa() == -1) {
        return false;
    }
    
    int v1_apa = v1_wpid.apa();
    int v1_face = v1_wpid.face();
   
    int v2_apa = v2_wpid.apa();
    int v2_face = v2_wpid.face();
    
    // Get time normalization factor from first blob
    auto first_blob = cluster.children()[0];
    int ntime_ticks = first_blob->slice_index_max() - first_blob->slice_index_min();
    if (ntime_ticks <= 0) ntime_ticks = 1;  // avoid division by zero
    
    // Convert vertices to U/V/W/T coordinates
    auto [v1_t_raw, v1_u_ch] = cluster.grouping()->convert_3Dpoint_time_ch(v1_p, v1_apa, v1_face, 0);
    auto [v1_t_u, v1_v_ch] = cluster.grouping()->convert_3Dpoint_time_ch(v1_p, v1_apa, v1_face, 1);
    auto [v1_t_v, v1_w_ch] = cluster.grouping()->convert_3Dpoint_time_ch(v1_p, v1_apa, v1_face, 2);
    
    auto [v2_t_raw, v2_u_ch] = cluster.grouping()->convert_3Dpoint_time_ch(v2_p, v2_apa, v2_face, 0);
    auto [v2_t_u, v2_v_ch] = cluster.grouping()->convert_3Dpoint_time_ch(v2_p, v2_apa, v2_face, 1);
    auto [v2_t_v, v2_w_ch] = cluster.grouping()->convert_3Dpoint_time_ch(v2_p, v2_apa, v2_face, 2);
    
    // Normalize time by ntime_ticks to account for convention difference
    double v1_t = double(v1_t_raw) / ntime_ticks;
    double v2_t = double(v2_t_raw) / ntime_ticks;
    
    double v1_u = v1_u_ch;
    double v1_v = v1_v_ch;
    double v1_w = v1_w_ch;
    
    double v2_u = v2_u_ch;
    double v2_v = v2_v_ch;
    double v2_w = v2_w_ch;
    
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster.get_scope_transform(cluster.get_default_scope()));
    double cluster_t0 = cluster.get_cluster_t0();
    
    int ncount_close = 0;
    int ncount_dead = 0;
    int ncount_line = 0;
    
    // Check each plane (U, V, W)
    for (int pind = 0; pind < 3; pind++) {
        double v1_wire, v2_wire;
        if (pind == 0) { v1_wire = v1_u; v2_wire = v2_u; }
        else if (pind == 1) { v1_wire = v1_v; v2_wire = v2_v; }
        else { v1_wire = v1_w; v2_wire = v2_w; }
        
        // Check if vertices are close in this 2D projection
        double dist_2d = std::sqrt(std::pow(v1_wire - v2_wire, 2) + std::pow(v1_t - v2_t, 2));
        
        if (dist_2d < 2.5) {
            ncount_close++;
        } else {
            // Check if segment points are all in dead region
            bool flag_dead = true;
            const auto& seg_fits = sg->fits();

            for (size_t i = 0; i < seg_fits.size(); i++) {
                auto test_wpid = dv->contained_by(seg_fits[i].point);
                auto p_raw = transform->backward(seg_fits[i].point, cluster_t0, test_wpid.face(), test_wpid.apa());
                if (!cluster.grouping()->get_closest_dead_chs(p_raw, 1, test_wpid.apa(), test_wpid.face(), pind)) {
                    flag_dead = false;
                    break;
                }
            }
            
            if (flag_dead) {
                ncount_dead++;
            } else {
                // Check if the third view forms a line
                // Find the other segment connected to v1
                SegmentPtr sg1 = nullptr;
                auto edge_range = boost::out_edges(v1d, graph);
                for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
                    SegmentPtr temp_sg = graph[*eit].segment;
                    if (temp_sg && temp_sg != sg) {
                        sg1 = temp_sg;
                        break;
                    }
                }
                
                if (!sg1) continue;
                
                const auto& pts_2 = sg1->fits();
                
                // Direction vector from v1 to v2 in this 2D projection
                Facade::geo_vector_t v1_2d(v2_wire - v1_wire, v2_t - v1_t, 0);
                Facade::geo_vector_t v2_2d(0, 0, 0);
                double min_dis = 1e9;
                Facade::geo_point_t start_p = v2_p;
                Facade::geo_point_t end_p;
                
                // Find point on sg1 that's approximately 9 units away from v1
                for (size_t i = 0; i < pts_2.size(); i++) {
                    auto test_wpid = dv->contained_by(pts_2[i].point);
                    auto [p_t_raw, p_wire_ch] = cluster.grouping()->convert_3Dpoint_time_ch(pts_2[i].point, test_wpid.apa(), test_wpid.face(), pind);
                    double p_t = double(p_t_raw) / ntime_ticks;
                    double p_wire = p_wire_ch;
                    
                    Facade::geo_vector_t v3(p_wire - v1_wire, p_t - v1_t, 0);
                    double dis = std::fabs(v3.magnitude() - 9.0);
                    if (dis < min_dis) {
                        min_dis = dis;
                        v2_2d = v3;
                        end_p = pts_2[i].point;
                    }
                }
                
                // Check angle between v1_2d and v2_2d
                double mag1 = v1_2d.magnitude();
                double mag2 = v2_2d.magnitude();
                double angle = 180.0;
                
                if (mag1 > 0 && mag2 > 0) {
                    double cos_angle = v1_2d.dot(v2_2d) / (mag1 * mag2);
                    if (cos_angle > 1.0) cos_angle = 1.0;
                    if (cos_angle < -1.0) cos_angle = -1.0;
                    angle = 180.0 - std::acos(cos_angle) / 3.1415926 * 180.0;
                }
                
                if (angle < 30.0 || (mag1 < 8.0 && angle < 35.0)) {
                    ncount_line++;
                } else {
                    // Check if path from v2 to end_p is good
                    double step_size = 0.6 * units::cm;
                    double path_length = ray_length(Ray{start_p, end_p});
                    int ncount = std::round(path_length / step_size);
                    int n_bad = 0;
                    
                    for (int i = 1; i < ncount; i++) {
                        Facade::geo_point_t test_p(
                            start_p.x() + (end_p.x() - start_p.x()) / ncount * i,
                            start_p.y() + (end_p.y() - start_p.y()) / ncount * i,
                            start_p.z() + (end_p.z() - start_p.z()) / ncount * i
                        );
                        
                        auto test_wpid = dv->contained_by(test_p);
                        if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                            auto temp_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            if (!cluster.grouping()->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2 * units::cm, 0, 0)) {
                                n_bad++;
                            }
                        }
                    }
                    
                    if (n_bad <= 1) {
                        ncount_line++;
                    }
                }
            }
        }
    }
    
    // Decision logic
    if (ncount_close >= 2 ||
        (ncount_close == 1 && ncount_dead == 1 && ncount_line >= 1) ||
        (ncount_close == 1 && ncount_dead == 2) ||
        (ncount_close == 1 && ncount_line >= 2) ||
        ncount_line >= 3) {
        return true;
    }
    
    return false;
}

bool PatternAlgorithms::examine_vertices_1(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    bool flag_continue = false;
    
    VertexPtr v1 = nullptr;
    VertexPtr v2 = nullptr;
    VertexPtr v3 = nullptr;
    
    // Iterate through all vertices in the graph
    auto [vbegin, vend] = boost::vertices(graph);
    for (auto vit = vbegin; vit != vend; ++vit) {
        VertexPtr vtx = graph[*vit].vertex;
        if (!vtx || vtx->cluster() != &cluster) continue;
        
        // Check if vertex has exactly 2 connections (potential candidate)
        if (boost::degree(*vit, graph) != 2) continue;
        
        // Get the two connected segments
        std::vector<SegmentPtr> connected_segments;
        auto [ebegin, eend] = boost::out_edges(*vit, graph);
        for (auto eit = ebegin; eit != eend; ++eit) {
            SegmentPtr seg = graph[*eit].segment;
            if (seg) connected_segments.push_back(seg);
        }
        
        if (connected_segments.size() != 2) continue;
        
        // Check each segment
        for (auto seg : connected_segments) {
            // Only consider short segments (<4cm)
            if (segment_track_length(seg) > 4.0 * units::cm) continue;
            
            // Find the other vertex of this segment
            VertexPtr vtx1 = find_other_vertex(graph, seg, vtx);
            if (!vtx1) continue;
            
            // The other vertex must have at least 2 connections
            if (!vtx1->descriptor_valid()) continue;
            auto vd1 = vtx1->get_descriptor();
            if (boost::degree(vd1, graph) < 2) continue;
            
            // Check if these two vertices represent the same physical point
            if (examine_vertices_1p(graph, vtx, vtx1, track_fitter, dv)) {
                v1 = vtx;
                v2 = vtx1;
                
                // Find v3: the vertex on the other side of v1
                for (auto seg2 : connected_segments) {
                    if (seg2 == seg) continue;
                    VertexPtr vtx_other = find_other_vertex(graph, seg2, v1);
                    if (vtx_other) {
                        v3 = vtx_other;
                        break;
                    }
                }
                
                flag_continue = true;
                break;
            }
        }
        
        if (flag_continue) break;
    }
    
    // Merge vertices if found
    if (v1 && v2 && v3) {
        // Find the segments to be removed
        SegmentPtr sg = find_segment(graph, v1, v2);  // segment between v1 and v2
        SegmentPtr sg1 = find_segment(graph, v1, v3); // segment between v1 and v3
        
        if (!sg || !sg1) {
            return false;
        }
        
        // Create new segment from v3 to v2 using Steiner graph shortest path
        auto path_points = do_rough_path(cluster, v3->wcpt().point, v2->wcpt().point);
        
        if (path_points.size() < 2) {
            return false;
        }
        
        // Create the new segment
        auto sg2 = create_segment_for_cluster(cluster, dv, path_points, 0);
        if (!sg2) {
            return false;
        }
        
        std::cout << "Cluster: " << cluster.ident() << " Merge Vertices Type I " 
                  << "combining two segments into new segment" << std::endl;
        
        // Add new segment to graph
        add_segment(graph, sg2, v2, v3);
        
        // Remove old segments and vertex
        remove_segment(graph, sg);
        remove_segment(graph, sg1);
        remove_vertex(graph, v1);
        
        // Update tracking
        track_fitter.do_multi_tracking(true, true, true);
    }
    
    return flag_continue;
}

bool PatternAlgorithms::examine_vertices_2(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    bool flag_continue = false;
    
    VertexPtr v1 = nullptr;
    VertexPtr v2 = nullptr;
    SegmentPtr sg = nullptr;
    
    // Iterate through all edges (segments) in the graph
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr segment = graph[*eit].segment;
        if (!segment || segment->cluster() != &cluster) continue;
        
        // Get the two vertices of this segment
        auto vertices = find_vertices(graph, segment);
        VertexPtr vtx1 = vertices.first;
        VertexPtr vtx2 = vertices.second;
        
        if (!vtx1 || !vtx2) continue;
        
        // Get positions (prefer fit point over wcpt)
        Facade::geo_point_t p1 = vtx1->fit().valid() ? vtx1->fit().point : vtx1->wcpt().point;
        Facade::geo_point_t p2 = vtx2->fit().valid() ? vtx2->fit().point : vtx2->wcpt().point;
        
        // Calculate distance between vertices
        double dis = ray_length(Ray{p1, p2});
        
        // Check if vertices are very close (<0.45cm) or moderately close (<1.5cm with both having degree 2)
        if (dis < 0.45 * units::cm) {
            v1 = vtx1;
            v2 = vtx2;
            sg = segment;
            flag_continue = true;
            break;
        } else if (dis < 1.5 * units::cm) {
            // Check if both vertices have exactly 2 connections
            if (!vtx1->descriptor_valid() || !vtx2->descriptor_valid()) continue;
            auto vd1 = vtx1->get_descriptor();
            auto vd2 = vtx2->get_descriptor();
            
            if (boost::degree(vd1, graph) == 2 && boost::degree(vd2, graph) == 2) {
                v1 = vtx1;
                v2 = vtx2;
                sg = segment;
                flag_continue = true;
                break;
            }
        }
    }
    
    // Merge vertices if found
    if (v1 && v2 && sg) {
        std::cout << "Cluster: " << cluster.ident() << " Merge Vertices Type II" << std::endl;
        
        // Remove the segment between v1 and v2
        remove_segment(graph, sg);
        
        // Collect all segments connected to v2 (excluding the one we just removed)
        std::vector<SegmentPtr> v2_segments;
        if (v2->descriptor_valid()) {
            auto vd2 = v2->get_descriptor();
            auto [ebegin2, eend2] = boost::out_edges(vd2, graph);
            for (auto eit2 = ebegin2; eit2 != eend2; ++eit2) {
                SegmentPtr seg2 = graph[*eit2].segment;
                if (seg2) {
                    v2_segments.push_back(seg2);
                }
            }
        }
        
        // For each segment connected to v2, create a new segment from v3 to v1
        for (auto old_seg : v2_segments) {
            // Find the other vertex (v3) connected to v2 through this segment
            VertexPtr v3 = find_other_vertex(graph, old_seg, v2);
            if (!v3) continue;
            
            // Create new segment from v3 to v1 using Steiner graph shortest path
            auto path_points = do_rough_path(cluster, v3->wcpt().point, v1->wcpt().point);
            
            if (path_points.size() < 2) continue;
            
            // Create the new segment
            auto new_seg = create_segment_for_cluster(cluster, dv, path_points, 0);
            if (!new_seg) continue;
            
            // Add new segment to graph
            add_segment(graph, new_seg, v3, v1);
            
            // Remove old segment
            remove_segment(graph, old_seg);
        }
        
        // Remove v2 vertex
        remove_vertex(graph, v2);
        
        // Clean up isolated vertices (vertices with no connections)
        std::vector<VertexPtr> isolated_vertices;
        auto [vbegin, vend] = boost::vertices(graph);
        for (auto vit = vbegin; vit != vend; ++vit) {
            if (boost::degree(*vit, graph) == 0) {
                VertexPtr vtx = graph[*vit].vertex;
                if (vtx) {
                    isolated_vertices.push_back(vtx);
                }
            }
        }
        
        for (auto vtx : isolated_vertices) {
            remove_vertex(graph, vtx);
        }
        
        // Update tracking
        track_fitter.do_multi_tracking(true, true, true);
    }
    
    return flag_continue;
}


bool PatternAlgorithms::examine_vertices_4p(Graph&graph, VertexPtr v1, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    // Find the segment between v1 and v2
    SegmentPtr sg1 = find_segment(graph, v1, v2);
    if (!sg1) return true;
    
    bool flag = true;
    
    // Get cluster information
    auto cluster = v1->cluster();
    if (!cluster) return true;
    
    // Get transform and grouping
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster->get_scope_transform(cluster->get_default_scope()));
    double cluster_t0 = cluster->get_cluster_t0();
    auto grouping = cluster->grouping();
    
    if (!transform || !grouping) return true;
    
    // Get v1 and v2 positions
    Facade::geo_point_t v1_point = v1->fit().valid() ? v1->fit().point : v1->wcpt().point;
    Facade::geo_point_t v2_point = v2->fit().valid() ? v2->fit().point : v2->wcpt().point;
    
    // Check segments of v1 with respect to v2
    if (!v1->descriptor_valid()) return true;
    auto vd1 = v1->get_descriptor();
    
    auto [ebegin, eend] = boost::out_edges(vd1, graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr sg = graph[*eit].segment;
        if (!sg || sg == sg1) continue;
        
        // Get segment points
        const auto& pts = sg->wcpts();
        if (pts.empty()) continue;
        
        // Find point on segment approximately 3cm from v1
        Facade::geo_point_t min_point = pts.front().point;
        double min_dis = 1e9;
        
        for (size_t i = 0; i < pts.size(); i++) {
            double dis = std::fabs(ray_length(Ray{pts[i].point, v1_point}) - 3.0 * units::cm);
            if (dis < min_dis) {
                min_dis = dis;
                min_point = pts[i].point;
            }
        }
        
        // Test connectivity from min_point to v2
        double step_size = 0.3 * units::cm;
        Facade::geo_point_t start_p = min_point;
        Facade::geo_point_t end_p = v2_point;
        
        double distance = ray_length(Ray{start_p, end_p});
        int ncount = std::round(distance / step_size);
        int n_bad = 0;
        
        for (int i = 1; i < ncount; i++) {
            Facade::geo_point_t test_p(
                start_p.x() + (end_p.x() - start_p.x()) / ncount * i,
                start_p.y() + (end_p.y() - start_p.y()) / ncount * i,
                start_p.z() + (end_p.z() - start_p.z()) / ncount * i
            );
            
            // Check if test point is in good region
            auto test_wpid = dv->contained_by(test_p);
            if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                auto temp_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2 * units::cm, 0, 0)) {
                    n_bad++;
                }
            }
        }
        
        // If any bad points found, return false
        if (n_bad != 0) {
            flag = false;
            break;
        }
    }
    
    return flag;
}

bool PatternAlgorithms::examine_vertices_4(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    bool flag_continue = false;
    
    // Drift direction (X direction)
    Facade::geo_vector_t drift_dir_abs(1, 0, 0);
    
    // Get steiner point cloud for later use
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
    
    // Iterate through all segments in the graph
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr sg = graph[*eit].segment;
        if (!sg || sg->cluster() != &cluster) continue;
        
        const auto& pts = sg->wcpts();
        if (pts.size() < 2) continue;
        
        // Get vertices
        auto pair_vertices = find_vertices(graph, sg);
        VertexPtr v1 = pair_vertices.first;
        VertexPtr v2 = pair_vertices.second;
        if (!v1 || !v2) continue;
        
        // Calculate segment direction
        Facade::geo_vector_t tmp_dir(
            pts.front().point.x() - pts.back().point.x(),
            pts.front().point.y() - pts.back().point.y(),
            pts.front().point.z() - pts.back().point.z()
        );
        
        // Calculate direct length between endpoints
        double direct_length = ray_length(Ray{pts.front().point, pts.back().point});
        // double track_length = segment_track_length(sg);
        
        // Calculate angle with drift direction (in degrees)
        double tmp_dir_mag = tmp_dir.magnitude();
        double angle = 90.0; // default perpendicular
        if (tmp_dir_mag > 0) {
            double cos_angle = drift_dir_abs.dot(tmp_dir) / tmp_dir_mag;
            // Clamp to [-1, 1] to avoid numerical issues with acos
            cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
            angle = std::acos(cos_angle) * 180.0 / M_PI;
        }
        
        // Check conditions: short segment OR perpendicular to drift
        if (direct_length < 2.0 * units::cm || 
            (tmp_dir.magnitude() < 3.5 * units::cm && std::fabs(angle - 90.0) < 10)) {
            
            // Check v1 first
            if (!v1->descriptor_valid() || !v2->descriptor_valid()) continue;
            auto vd1 = v1->get_descriptor();
            auto vd2 = v2->get_descriptor();
            
            if (boost::degree(vd1, graph) >= 2 && examine_vertices_4p(graph, v1, v2, track_fitter, dv)) {
                // Merge v1's segments to v2
                
                // Get v2 position
                Facade::geo_point_t vtx_new_point = v2->wcpt().point;
                
                // Collect segments connected to v1 (except sg)
                std::vector<SegmentPtr> v1_segments;
                auto [e1begin, e1end] = boost::out_edges(vd1, graph);
                for (auto e1it = e1begin; e1it != e1end; ++e1it) {
                    SegmentPtr sg1 = graph[*e1it].segment;
                    if (sg1 && sg1 != sg) {
                        v1_segments.push_back(sg1);
                    }
                }
                
                // Process each segment connected to v1
                for (auto sg1 : v1_segments) {
                    const auto& vec_wcps = sg1->wcpts();
                    if (vec_wcps.empty()) continue;
                    
                    // Determine which end connects to v1
                    bool flag_front = (ray_length(Ray{vec_wcps.front().point, v1->wcpt().point}) <
                                      ray_length(Ray{vec_wcps.back().point, v1->wcpt().point}));
                    
                    // Get v1 position
                    Facade::geo_point_t v1_point = v1->fit().valid() ? v1->fit().point : v1->wcpt().point;
                    
                    // Find point ~3cm from v1 on this segment
                    WCPoint min_wcp = vec_wcps.front();
                    double min_dis = 1e9;
                    
                    // Calculate max distance to determine dis_cut
                    double max_dis = std::max(
                        ray_length(Ray{vec_wcps.front().point, v1_point}),
                        ray_length(Ray{vec_wcps.back().point, v1_point})
                    );
                    double dis_cut = 0;
                    double default_dis_cut = 2.5 * units::cm;
                    if (max_dis > 2 * default_dis_cut) dis_cut = default_dis_cut;
                    
                    for (size_t j = 0; j < vec_wcps.size(); j++) {
                        double dis1 = ray_length(Ray{vec_wcps[j].point, v1_point});
                        double dis = std::fabs(dis1 - 3.0 * units::cm);
                        if (dis < min_dis && dis1 > dis_cut) {
                            min_wcp = vec_wcps[j];
                            min_dis = dis;
                        }
                    }
                    
                    // Build new path from v2 to min_wcp using Steiner graph
                    std::list<WCPoint> new_list;
                    new_list.push_back(v2->wcpt());
                    
                    // Add intermediate points
                    double dis_step = 2.0 * units::cm;
                    int ncount = std::round(ray_length(Ray{vtx_new_point, min_wcp.point}) / dis_step);
                    if (ncount < 2) ncount = 2;
                    
                    for (int qx = 1; qx < ncount; qx++) {
                        Facade::geo_point_t tmp_p(
                            vtx_new_point.x() + (min_wcp.point.x() - vtx_new_point.x()) / ncount * qx,
                            vtx_new_point.y() + (min_wcp.point.y() - vtx_new_point.y()) / ncount * qx,
                            vtx_new_point.z() + (min_wcp.point.z() - vtx_new_point.z()) / ncount * qx
                        );
                        
                        // Find closest steiner point
                        auto knn_results = cluster.kd_steiner_knn(1, tmp_p, "steiner_pc");
                        if (!knn_results.empty()) {
                            size_t idx = knn_results[0].first;
                            WCPoint tmp_wcp;
                            tmp_wcp.point = Facade::geo_point_t(x_coords[idx], y_coords[idx], z_coords[idx]);
                            
                            double dist_to_steiner = ray_length(Ray{tmp_wcp.point, tmp_p});
                            if (dist_to_steiner > 0.3 * units::cm) continue;
                            
                            // Check if not duplicate
                            bool is_duplicate = (ray_length(Ray{tmp_wcp.point, new_list.back().point}) < 0.01 * units::cm);
                            bool is_min_wcp = (ray_length(Ray{tmp_wcp.point, min_wcp.point}) < 0.01 * units::cm);
                            
                            if (!is_duplicate && !is_min_wcp) {
                                new_list.push_back(tmp_wcp);
                            }
                        }
                    }
                    new_list.push_back(min_wcp);
                    
                    // Combine with rest of segment
                    std::list<WCPoint> old_list(vec_wcps.begin(), vec_wcps.end());
                    
                    if (flag_front) {
                        // Remove points up to min_wcp from front
                        while (!old_list.empty() && 
                               ray_length(Ray{old_list.front().point, min_wcp.point}) > 0.01 * units::cm) {
                            old_list.pop_front();
                        }
                        if (!old_list.empty()) old_list.pop_front();
                        
                        // Prepend new path
                        for (auto it = new_list.rbegin(); it != new_list.rend(); ++it) {
                            old_list.push_front(*it);
                        }
                    } else {
                        // Remove points up to min_wcp from back
                        while (!old_list.empty() && 
                               ray_length(Ray{old_list.back().point, min_wcp.point}) > 0.01 * units::cm) {
                            old_list.pop_back();
                        }
                        if (!old_list.empty()) old_list.pop_back();
                        
                        // Append new path
                        for (auto it = new_list.rbegin(); it != new_list.rend(); ++it) {
                            old_list.push_back(*it);
                        }
                    }
                    
                    // Update segment with new points
                    std::vector<WCPoint> new_wcpts(old_list.begin(), old_list.end());
                    sg1->wcpts(new_wcpts);
                    
                    // Find other vertex and update connection
                    VertexPtr v3 = find_other_vertex(graph, sg1, v1);
                    if (v3) {
                        remove_segment(graph, sg1);
                        add_segment(graph, sg1, v2, v3);
                    }
                }
                
                // Remove v1 and sg
                remove_vertex(graph, v1);
                remove_segment(graph, sg);
                
                flag_continue = true;
                std::cout << "Cluster: " << cluster.ident() << " Merge Vertices Type III" << std::endl;
                track_fitter.do_multi_tracking(true, true, true);
                break;
                
            } else if (boost::degree(vd2, graph) >= 2 && examine_vertices_4p(graph, v2, v1, track_fitter, dv)) {
                // Merge v2's segments to v1 (symmetric case)
                
                // Get v1 position
                Facade::geo_point_t vtx_new_point = v1->wcpt().point;
                
                // Collect segments connected to v2 (except sg)
                std::vector<SegmentPtr> v2_segments;
                auto [e2begin, e2end] = boost::out_edges(vd2, graph);
                for (auto e2it = e2begin; e2it != e2end; ++e2it) {
                    SegmentPtr sg1 = graph[*e2it].segment;
                    if (sg1 && sg1 != sg) {
                        v2_segments.push_back(sg1);
                    }
                }
                
                // Process each segment connected to v2
                for (auto sg1 : v2_segments) {
                    const auto& vec_wcps = sg1->wcpts();
                    if (vec_wcps.empty()) continue;
                    
                    // Determine which end connects to v2
                    bool flag_front = (ray_length(Ray{vec_wcps.front().point, v2->wcpt().point}) <
                                      ray_length(Ray{vec_wcps.back().point, v2->wcpt().point}));
                    
                    // Get v2 position
                    Facade::geo_point_t v2_point = v2->fit().valid() ? v2->fit().point : v2->wcpt().point;
                    
                    // Find point ~3cm from v2 on this segment
                    WCPoint min_wcp = vec_wcps.front();
                    double min_dis = 1e9;
                    
                    // Calculate max distance to determine dis_cut
                    double max_dis = std::max(
                        ray_length(Ray{vec_wcps.front().point, v2_point}),
                        ray_length(Ray{vec_wcps.back().point, v2_point})
                    );
                    double dis_cut = 0;
                    double default_dis_cut = 2.5 * units::cm;
                    if (max_dis > 2 * default_dis_cut) dis_cut = default_dis_cut;
                    
                    for (size_t j = 0; j < vec_wcps.size(); j++) {
                        double dis1 = ray_length(Ray{vec_wcps[j].point, v2_point});
                        double dis = std::fabs(dis1 - 3.0 * units::cm);
                        if (dis < min_dis && dis1 > dis_cut) {
                            min_wcp = vec_wcps[j];
                            min_dis = dis;
                        }
                    }
                    
                    // Build new path from v1 to min_wcp using Steiner graph
                    std::list<WCPoint> new_list;
                    new_list.push_back(v1->wcpt());
                    
                    // Add intermediate points
                    double dis_step = 2.0 * units::cm;
                    int ncount = std::round(ray_length(Ray{vtx_new_point, min_wcp.point}) / dis_step);
                    if (ncount < 2) ncount = 2;
                    
                    for (int qx = 1; qx < ncount; qx++) {
                        Facade::geo_point_t tmp_p(
                            vtx_new_point.x() + (min_wcp.point.x() - vtx_new_point.x()) / ncount * qx,
                            vtx_new_point.y() + (min_wcp.point.y() - vtx_new_point.y()) / ncount * qx,
                            vtx_new_point.z() + (min_wcp.point.z() - vtx_new_point.z()) / ncount * qx
                        );
                        
                        // Find closest steiner point
                        auto knn_results = cluster.kd_steiner_knn(1, tmp_p, "steiner_pc");
                        if (!knn_results.empty()) {
                            size_t idx = knn_results[0].first;
                            WCPoint tmp_wcp;
                            tmp_wcp.point = Facade::geo_point_t(x_coords[idx], y_coords[idx], z_coords[idx]);
                            
                            double dist_to_steiner = ray_length(Ray{tmp_wcp.point, tmp_p});
                            if (dist_to_steiner > 0.3 * units::cm) continue;
                            
                            // Check if not duplicate
                            bool is_duplicate = (ray_length(Ray{tmp_wcp.point, new_list.back().point}) < 0.01 * units::cm);
                            bool is_min_wcp = (ray_length(Ray{tmp_wcp.point, min_wcp.point}) < 0.01 * units::cm);
                            
                            if (!is_duplicate && !is_min_wcp) {
                                new_list.push_back(tmp_wcp);
                            }
                        }
                    }
                    new_list.push_back(min_wcp);
                    
                    // Combine with rest of segment
                    std::list<WCPoint> old_list(vec_wcps.begin(), vec_wcps.end());
                    
                    if (flag_front) {
                        // Remove points up to min_wcp from front
                        while (!old_list.empty() && 
                               ray_length(Ray{old_list.front().point, min_wcp.point}) > 0.01 * units::cm) {
                            old_list.pop_front();
                        }
                        if (!old_list.empty()) old_list.pop_front();
                        
                        // Prepend new path
                        for (auto it = new_list.rbegin(); it != new_list.rend(); ++it) {
                            old_list.push_front(*it);
                        }
                    } else {
                        // Remove points up to min_wcp from back
                        while (!old_list.empty() && 
                               ray_length(Ray{old_list.back().point, min_wcp.point}) > 0.01 * units::cm) {
                            old_list.pop_back();
                        }
                        if (!old_list.empty()) old_list.pop_back();
                        
                        // Append new path
                        for (auto it = new_list.rbegin(); it != new_list.rend(); ++it) {
                            old_list.push_back(*it);
                        }
                    }
                    
                    // Update segment with new points
                    std::vector<WCPoint> new_wcpts(old_list.begin(), old_list.end());
                    sg1->wcpts(new_wcpts);
                    
                    // Find other vertex and update connection
                    VertexPtr v3 = find_other_vertex(graph, sg1, v2);
                    if (v3) {
                        remove_segment(graph, sg1);
                        add_segment(graph, sg1, v1, v3);
                    }
                }
                
                // Remove v2 and sg
                remove_vertex(graph, v2);
                remove_segment(graph, sg);
                
                flag_continue = true;
                std::cout << "Cluster: " << cluster.ident() << " Merge Vertices Type III" << std::endl;
                track_fitter.do_multi_tracking(true, true, true);
                break;
            }
        }
        
        if (flag_continue) break;
    }
    
    return flag_continue;
}

void PatternAlgorithms::examine_vertices(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    bool flag_continue = true;
    
    while (flag_continue) {
        flag_continue = false;
        
        // Examine and clean up segment topology
        examine_segment(graph, cluster, track_fitter, dv);
        
        // Merge vertex if the kink is not at right location (Type I)
        flag_continue = flag_continue || examine_vertices_1(graph, cluster, track_fitter, dv);
        
        // Count vertices in the graph
        size_t num_vertices = boost::num_vertices(graph);
        
        // Merge vertices if they are too close (Type II) - only if more than 2 vertices
        if (num_vertices > 2) {
            flag_continue = flag_continue || examine_vertices_2(graph, cluster, track_fitter, dv);
        }
        
        // Merge vertices if they are reasonably close (Type III/IV)
        flag_continue = flag_continue || examine_vertices_4(graph, cluster, track_fitter, dv);
    }
}












