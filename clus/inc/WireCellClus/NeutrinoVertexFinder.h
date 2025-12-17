#include "WireCellClus/NeutrinoPatternBase.h"

bool WireCell::Clus::PR::PatternAlgorithms::search_for_vertex_activities(Graph& graph, VertexPtr vertex, std::set<SegmentPtr>& segments_set, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, double search_range){
    // Get steiner point cloud and terminal flags
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
    const auto& flag_steiner_terminal = steiner_pc.get("flag_steiner_terminal")->elements<int>();
    
    // Get transform and grouping for point validation
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster.get_scope_transform(cluster.get_default_scope()));
    double cluster_t0 = cluster.get_cluster_t0();
    auto grouping = cluster.grouping();
    
    if (!transform || !grouping) return false;
    
    // Get vertex position
    WireCell::Point vtx_point = vertex->fit().valid() ? vertex->fit().point : vertex->wcpt().point;
    
    // Collect directions from existing segments
    std::vector<WireCell::Vector> saved_dirs;
    for (auto seg : segments_set) {
        WireCell::Vector dir = vertex_segment_get_dir(vertex, seg, graph, 5*units::cm);
        if (dir.magnitude() != 0) {
            saved_dirs.push_back(dir);
        }
    }
    
    // Get candidate points within search range
    auto candidate_results = cluster.kd_steiner_radius(search_range, vtx_point, "steiner_pc");
    
    double max_dis = 0;
    size_t max_idx = 0;
    bool found = false;
    
    // First round: look for points with good angular separation and charge
    for (const auto& [idx, dist_sq] : candidate_results) {
        if (!flag_steiner_terminal[idx]) continue;
        
        // Skip if this is the vertex itself
        if (ray_length(Ray{WireCell::Point(x_coords[idx], y_coords[idx], z_coords[idx]), vertex->wcpt().point}) < 0.01*units::cm) continue;
        
        double dis = std::sqrt(dist_sq);
        WireCell::Point test_p(x_coords[idx], y_coords[idx], z_coords[idx]);
        
        // Find minimum distance to all segments
        double min_dis = 1e9;
        double min_dis_u = 1e9;
        double min_dis_v = 1e9;
        double min_dis_w = 1e9;
        
        auto test_wpid = dv->contained_by(test_p);
        if (test_wpid.face() == -1 || test_wpid.apa() == -1) continue;
        
        for (auto it = boost::edges(graph).first; it != boost::edges(graph).second; ++it) {
            SegmentPtr sg = graph[*it].segment;
            if (!sg || sg->cluster() != &cluster) continue;
            
            auto [dist_3d, closest_pt] = segment_get_closest_point(sg, test_p, "fit");
            if (dist_3d < min_dis) min_dis = dist_3d;
            
            auto [dist_u, dist_v, dist_w] = segment_get_closest_2d_distances(sg, test_p, test_wpid.apa(), test_wpid.face(), "fit");
            if (dist_u < min_dis_u) min_dis_u = dist_u;
            if (dist_v < min_dis_v) min_dis_v = dist_v;
            if (dist_w < min_dis_w) min_dis_w = dist_w;
        }
        
        if (min_dis > 0.6*units::cm && min_dis_u + min_dis_v + min_dis_w > 1.2*units::cm) {
            WireCell::Vector dir(test_p.x() - vtx_point.x(), test_p.y() - vtx_point.y(), test_p.z() - vtx_point.z());
            double sum_angle = 0;
            double min_angle = 1e9;
            
            for (size_t j = 0; j < saved_dirs.size(); j++) {
                double angle = std::acos(dir.dot(saved_dirs[j]) / (dir.magnitude() * saved_dirs[j].magnitude())) / 3.1415926 * 180.0;
                sum_angle += angle;
                if (angle < min_angle) min_angle = angle;
            }
            
            // Get average charge
            double sum_charge = 0;
            int ncount = 0;
            auto test_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
            
            for (int plane = 0; plane < 3; plane++) {
                if (!grouping->get_closest_dead_chs(test_p_raw, plane, test_wpid.apa(), test_wpid.face())) {
                    sum_charge += grouping->get_ave_charge(test_p_raw, 0.3*units::cm, plane, test_wpid.apa(), test_wpid.face());
                    ncount++;
                }
            }
            if (ncount != 0) sum_charge /= ncount;
            
            if ((sum_angle) * (sum_charge + 1e-9) > max_dis) {
                max_dis = (sum_angle) * (sum_charge + 1e-9);
                max_idx = idx;
                found = true;
            }
        }
    }
    
    // Second round: if nothing found, use relaxed criteria
    if (max_dis == 0) {
        for (const auto& [idx, dist_sq] : candidate_results) {
            if (!flag_steiner_terminal[idx]) continue;
            
            // Skip if this is the vertex itself
            if (ray_length(Ray{WireCell::Point(x_coords[idx], y_coords[idx], z_coords[idx]), vertex->wcpt().point}) < 0.01*units::cm) continue;
            
            double dis = std::sqrt(dist_sq);
            WireCell::Point test_p(x_coords[idx], y_coords[idx], z_coords[idx]);
            
            // Find minimum distance to all segments
            double min_dis = 1e9;
            double min_dis_u = 1e9;
            double min_dis_v = 1e9;
            double min_dis_w = 1e9;
            
            auto test_wpid = dv->contained_by(test_p);
            if (test_wpid.face() == -1 || test_wpid.apa() == -1) continue;
            
            for (auto it = boost::edges(graph).first; it != boost::edges(graph).second; ++it) {
                SegmentPtr sg = graph[*it].segment;
                if (!sg || sg->cluster() != &cluster) continue;
                
                auto [dist_3d, closest_pt] = segment_get_closest_point(sg, test_p, "fit");
                if (dist_3d < min_dis) min_dis = dist_3d;
                
                auto [dist_u, dist_v, dist_w] = segment_get_closest_2d_distances(sg, test_p, test_wpid.apa(), test_wpid.face(), "fit");
                if (dist_u < min_dis_u) min_dis_u = dist_u;
                if (dist_v < min_dis_v) min_dis_v = dist_v;
                if (dist_w < min_dis_w) min_dis_w = dist_w;
            }
            
            if (min_dis > 0.36*units::cm && min_dis_u + min_dis_v + min_dis_w > 0.8*units::cm) {
                // Get average charge
                double sum_charge = 0;
                int ncount = 0;
                auto test_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                
                for (int plane = 0; plane < 3; plane++) {
                    if (!grouping->get_closest_dead_chs(test_p_raw, plane, test_wpid.apa(), test_wpid.face())) {
                        sum_charge += grouping->get_ave_charge(test_p_raw, 0.3*units::cm, plane, test_wpid.apa(), test_wpid.face());
                        ncount++;
                    }
                }
                if (ncount != 0) sum_charge /= ncount;
                
                if (min_dis + (min_dis_u + min_dis_v + min_dis_w) / std::sqrt(3.0) > max_dis && sum_charge > 20000) {
                    max_dis = min_dis + (min_dis_u + min_dis_v + min_dis_w) / std::sqrt(3.0);
                    max_idx = idx;
                    found = true;
                }
            }
        }
    }
    
    // If a good candidate was found, create new vertex and segment
    if (found && max_dis != 0) {
        WireCell::Point max_point(x_coords[max_idx], y_coords[max_idx], z_coords[max_idx]);
        
        // Create new vertex at the found point
        auto v1 = make_vertex(graph);
        WCPoint new_wcp;
        new_wcp.point = max_point;
        v1->wcpt(new_wcp).cluster(&cluster);
        
        // Build path from vertex to new vertex using steiner point cloud
        WireCell::Point mid_p(
            (vertex->wcpt().point.x() + max_point.x()) / 2.0,
            (vertex->wcpt().point.y() + max_point.y()) / 2.0,
            (vertex->wcpt().point.z() + max_point.z()) / 2.0
        );
        
        auto [mid_idx, mid_pt] = cluster.get_closest_wcpoint(mid_p);
        
        std::list<WireCell::Point> wcp_list;
        wcp_list.push_back(vertex->wcpt().point);
        
        if (ray_length(Ray{mid_pt, wcp_list.back()}) > 0.01*units::cm) {
            wcp_list.push_back(mid_pt);
        }
        
        if (ray_length(Ray{max_point, wcp_list.back()}) > 0.01*units::cm) {
            wcp_list.push_back(max_point);
        }
        
        if (wcp_list.size() > 1) {
            std::cout << "Cluster: " << cluster.ident() << " Vertex Activity Found at " << mid_p << std::endl;
            
            // Convert to vector for segment creation
            std::vector<WireCell::Point> path_points(wcp_list.begin(), wcp_list.end());
            
            // Create new segment
            auto sg1 = create_segment_for_cluster(cluster, dv, path_points, 0);
            if (sg1) {
                add_segment(graph, sg1, v1, vertex);
                return true;
            }
        }
    }
    
    return false;
}

