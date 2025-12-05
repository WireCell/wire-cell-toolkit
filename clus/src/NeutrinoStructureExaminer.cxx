#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

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

