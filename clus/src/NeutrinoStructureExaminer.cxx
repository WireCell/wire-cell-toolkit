#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"

#include "WireCellAux/Logger.h"

#include <algorithm>

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

void PatternAlgorithms::examine_structure(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    // Change 2 to 1 (merge two segments into one straight segment)
    if (examine_structure_2(graph, cluster, track_fitter, dv)) {
        track_fitter.do_multi_tracking(true, true, false, false, false, &cluster);
    }

    // Straighten 1 (replace curved segments with straight lines)
    if (examine_structure_1(graph, cluster, track_fitter, dv)) {
        track_fitter.do_multi_tracking(true, true, false, false, false, &cluster);
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
            
            // Get the original WCPoints (needed for endpoint metadata later)
            const auto& wcpts = sg->wcpts();
            if (wcpts.size() < 2) continue;
            
            // Always use wcpts endpoints for the straight-line scan direction.
            // Using fits here caused a mismatch: fits can span a larger range than
            // wcpts, so the scan would go beyond wcpts_back, the stop condition
            // would never fire, and wcpts_back would be force-appended as a reverse
            // kink at the end of the new path.
            Facade::geo_point_t start_p = wcpts.front().point;
            Facade::geo_point_t end_p   = wcpts.back().point;
            
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
                
                // Check if this point is good.  Points outside any TPC (face==-1)
                // are treated as bad — prototype's is_good_point returns false there.
                auto test_wpid = dv->contained_by(test_p);
                if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                    auto temp_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                    if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*units::cm, 0, 0)) {
                        n_bad++;
                    }
                } else {
                    n_bad++;  // outside all TPCs → bad sample (B.3)
                }

                if (n_bad > 1) {
                    flag_replace = false;
                    break;
                }
            }

            // If the straight line is better, replace the segment path
            if (flag_replace) {
                if (!cluster.has_pc("steiner_pc")) continue;
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
                std::vector<Facade::geo_point_t> main_pts_1;
                for (const auto& wcp : new_wcpts) main_pts_1.push_back(wcp.point);
                create_segment_point_cloud(sg, main_pts_1, dv, "main");

                flag_update = true;
                // std::cout << "Cluster: " << cluster.ident() << " replace Track Content with Straight Line for segment with length " << length/units::cm << " cm" << std::endl;
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

        // Iterate in insertion order for deterministic results
        for (const auto& vd_cur : ordered_nodes(graph)) {
            VertexPtr vtx = graph[vd_cur].vertex;

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
            
            // Always use wcpts endpoints for the straight-line scan direction.
            // Using fit().point here caused a mismatch: fit positions can differ from
            // wcpt positions, so the scan would miss the wcpt anchor, the 0.01 cm
            // stop/append condition would not fire correctly, and end_wcp would be
            // force-appended as a kink at the end of the new path.
            Facade::geo_point_t start_p = vtx1->wcpt().point;
            Facade::geo_point_t end_p   = vtx2->wcpt().point;
            
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
                
                // Check if this point is good.  Points outside any TPC (face==-1)
                // are treated as bad — prototype's is_good_point returns false there.
                auto test_wpid = dv->contained_by(test_p);
                if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                    auto temp_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                    if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*units::cm, 0, 0)) {
                        n_bad++;
                    }
                } else {
                    n_bad++;  // outside all TPCs → bad sample (B.3)
                }

                if (n_bad > 1) {
                    flag_replace = false;
                    break;
                }
            }

            // If the straight line is better, merge the two segments
            if (flag_replace) {
                // std::cout << "Cluster: " << cluster.ident() << " Merge two segments with a straight one, vtx at (" 
                //           << vtx->wcpt().point.x()/units::cm << ", "
                //           << vtx->wcpt().point.y()/units::cm << ", "
                //           << vtx->wcpt().point.z()/units::cm << ") cm" << std::endl;
                
                // Check if the two endpoint vertices are at the same position
                double dist_vtx1_vtx2 = ray_length(Ray{vtx1->wcpt().point, vtx2->wcpt().point});
                const double distance_threshold = 0.01 * units::cm;
                
                if (dist_vtx1_vtx2 < distance_threshold) {
                    // vtx1 and vtx2 are co-located: no new segment is created.
                    // merge_vertex_into_another(vtx2, vtx1) is deferred until AFTER
                    // sg1/sg2/vtx are removed below — calling it now would crash
                    // because merge_vertex_into_another removes and re-adds sg2 as a
                    // vtx1-vtx edge, but that slot is already occupied by sg1; with
                    // setS edges add_segment aliases sg2's descriptor to sg1's, so the
                    // subsequent remove_segment(sg2) destroys the shared edge. (B.7)
                } else {
                    // Create a new segment with straight line path
                    if (!cluster.has_pc("steiner_pc")) continue;
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
                    std::vector<Facade::geo_point_t> main_pts_2;
                    for (const auto& wcp : new_wcpts) main_pts_2.push_back(wcp.point);
                    create_segment_point_cloud(new_seg, main_pts_2, dv, "main");

                    // Add the new segment to the graph
                    add_segment(graph, new_seg, vtx1, vtx2);
                }
                
                // Delete the old segments and middle vertex
                remove_segment(graph, sg1);
                remove_segment(graph, sg2);
                remove_vertex(graph, vtx);

                // B.7: Co-located case — now that sg1/sg2/vtx are gone, safely
                // merge vtx2 into vtx1 to preserve vtx2's other segments
                // (matching prototype behavior of re-parenting vtx2's connections).
                // The aliasing risk is gone: sg2 (vtx2→vtx) and sg1 (vtx1→vtx)
                // are already removed, so no parallel-edge conflict can occur.
                if (dist_vtx1_vtx2 < distance_threshold && vtx2->descriptor_valid()) {
                    if (boost::degree(vtx2->get_descriptor(), graph) > 0) {
                        merge_vertex_into_another(graph, vtx2, vtx1, dv);
                    } else {
                        remove_vertex(graph, vtx2);
                    }
                }

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

        // Iterate in insertion order for deterministic results
        for (const auto& vd_cur : ordered_nodes(graph)) {
            VertexPtr vtx = graph[vd_cur].vertex;

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
            
            // Calculate direction vectors at 10cm
            WireCell::Vector dir1 = segment_cal_dir_3vector(sg1, vtx_point, 10*units::cm);
            WireCell::Vector dir2 = segment_cal_dir_3vector(sg2, vtx_point, 10*units::cm);
            
            if (dir1.magnitude() == 0 || dir2.magnitude() == 0) continue;

            // Calculate 10cm angle (180 - angle between directions)
            double angle_10cm = (3.1415926 - std::acos(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude()))) / 3.1415926 * 180.0;
            
            // Only compute 3cm vectors and angle if 10cm angle passes
            if (angle_10cm < 18) {
                WireCell::Vector dir3 = segment_cal_dir_3vector(sg1, vtx_point, 3*units::cm);
                WireCell::Vector dir4 = segment_cal_dir_3vector(sg2, vtx_point, 3*units::cm);
                
                if (dir3.magnitude() == 0 || dir4.magnitude() == 0) continue;
                
                double angle_3cm = (3.1415926 - std::acos(dir3.dot(dir4) / (dir3.magnitude() * dir4.magnitude()))) / 3.1415926 * 180.0;
                
                // Check if segments are nearly collinear (small 3cm angle)
                if (angle_3cm < 27) {
                    SPDLOG_LOGGER_DEBUG(s_log, "examine_structure: Cluster: {} Merge two segments into one according to angle {}° (10cm) and {}° (3cm)", cluster.ident(), angle_10cm, angle_3cm);
                    
                    // Merge the two segments by combining their WCPoint lists.
                    // Use graph topology to determine orientation — i.e. which wcpt
                    // endpoint of each segment is at the shared vertex (vtx) — rather
                    // than a geometric distance threshold.  A distance threshold (the
                    // prior approach, 0.01 cm) can silently fail when ES1/ES2 kNN
                    // projection shifts wcpts by up to ~0.15 cm. (E.2)
                    const auto& wcpts1 = sg1->wcpts();
                    const auto& wcpts2 = sg2->wcpts();

                    // find_vertices returns (v_front, v_back) where v_front is the
                    // vertex nearest to wcpts.front(). Compare by pointer identity.
                    auto [sg1_front_vtx, sg1_back_vtx] = find_vertices(graph, sg1);
                    auto [sg2_front_vtx, sg2_back_vtx] = find_vertices(graph, sg2);
                    bool sg1_vtx_front = (sg1_front_vtx == vtx);  // vtx at wcpts1.front()
                    bool sg2_vtx_front = (sg2_front_vtx == vtx);  // vtx at wcpts2.front()

                    std::vector<WCPoint> merged_wcpts;
                    merged_wcpts.reserve(wcpts1.size() + wcpts2.size());

                    if (sg1_vtx_front && sg2_vtx_front) {
                        // front1 ↔ vtx ↔ front2: reverse wcpts2, then wcpts1 from idx 1
                        for (auto it = wcpts2.rbegin(); it != wcpts2.rend(); ++it)
                            merged_wcpts.push_back(*it);
                        for (size_t i = 1; i < wcpts1.size(); ++i)
                            merged_wcpts.push_back(wcpts1[i]);
                    } else if (sg1_vtx_front && !sg2_vtx_front) {
                        // front1 ↔ vtx ↔ back2: wcpts2 forward, then wcpts1 from idx 1
                        for (const auto& wcp : wcpts2)
                            merged_wcpts.push_back(wcp);
                        for (size_t i = 1; i < wcpts1.size(); ++i)
                            merged_wcpts.push_back(wcpts1[i]);
                    } else if (!sg1_vtx_front && sg2_vtx_front) {
                        // back1 ↔ vtx ↔ front2: wcpts1 forward, then wcpts2 from idx 1
                        for (const auto& wcp : wcpts1)
                            merged_wcpts.push_back(wcp);
                        for (size_t i = 1; i < wcpts2.size(); ++i)
                            merged_wcpts.push_back(wcpts2[i]);
                    } else {
                        // back1 ↔ vtx ↔ back2: wcpts1 forward, then reverse wcpts2 from rbegin+1
                        for (const auto& wcp : wcpts1)
                            merged_wcpts.push_back(wcp);
                        for (auto it = wcpts2.rbegin() + 1; it != wcpts2.rend(); ++it)
                            merged_wcpts.push_back(*it);
                    }

                    if (merged_wcpts.empty()) continue;  // degenerate: both wcpts lists empty

                    // Create new segment with merged points
                    auto new_seg = make_segment();
                    new_seg->wcpts(merged_wcpts).cluster(&cluster).dirsign(0);
                    std::vector<Facade::geo_point_t> main_pts_3;
                    for (const auto& wcp : merged_wcpts) main_pts_3.push_back(wcp.point);
                    create_segment_point_cloud(new_seg, main_pts_3, dv, "main");

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

    if (!cluster.has_pc("steiner_pc")) return false;

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
        // Use sqrt of pre-computed squared distance from kd_steiner_radius (E.4)
        double dis = std::sqrt(dist_sq);

        // B.2/B.4: resolve (apa,face) for test_p once, outside the segment loop.
        // Terminals outside every TPC are rejected immediately — without this guard,
        // min_dis_{u,v,w} stay at 1e9 and the 2D-distance criterion trivially passes.
        auto test_wpid = dv->contained_by(test_p);
        if (test_wpid.face() == -1 || test_wpid.apa() == -1) continue;

        // Find minimum distances to all segments.
        // The graph is per-grouping (one graph for all clusters in the event),
        // so this iterates segments from ALL clusters — matching the prototype's
        // unconstrained map_segment_vertices scan. (B.6)
        double min_dis = 1e9;
        double min_dis_u = 1e9;
        double min_dis_v = 1e9;
        double min_dis_w = 1e9;

        // B.5: use ordered_edges for deterministic tie-breaking in max_dis selection
        for (const auto& ed : ordered_edges(graph)) {
            SegmentPtr sg = graph[ed].segment;
            if (!sg) continue;

            // Get closest point distance from both fitted points and wcpts
            auto [dist_3d, closest_pt] = segment_get_closest_point(sg, test_p, "fit");
            if (dist_3d < min_dis) min_dis = dist_3d;

            // Also check wcpts for closest distance
            auto [dist_wcpt, closest_wcp] = segment_get_closest_point(sg, test_p, "main");
            if (dist_wcpt < min_dis) min_dis = dist_wcpt;

            // 2D distances — only meaningful within the same APA/face as test_p
            auto [dist_u, dist_v, dist_w] = segment_get_closest_2d_distances(sg, test_p, test_wpid.apa(), test_wpid.face(), "fit");
            if (dist_u < min_dis_u) min_dis_u = dist_u;
            if (dist_v < min_dis_v) min_dis_v = dist_v;
            if (dist_w < min_dis_w) min_dis_w = dist_w;
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

                // B.3: outside-TPC steps count as bad (prototype's is_good_point
                // returns false for points with no wires/charge)
                auto test_wpid1 = dv->contained_by(test_p1);
                if (test_wpid1.face() != -1 && test_wpid1.apa() != -1) {
                    auto temp_p_raw = transform->backward(test_p1, cluster_t0, test_wpid1.face(), test_wpid1.apa());
                    if (!grouping->is_good_point(temp_p_raw, test_wpid1.apa(), test_wpid1.face(), 0.3*units::cm, 0, 0)) {
                        n_bad++;
                    }
                } else {
                    n_bad++;  // outside all TPCs → bad sample (B.3)
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
        
        SPDLOG_LOGGER_DEBUG(s_log, "examine_structure: Cluster: {} Add a track to the main vertex, {} points", cluster.ident(), wcp_list.size());
        
        // Create new segment
        auto sg1 = make_segment();
        sg1->wcpts(wcp_list).cluster(&cluster).dirsign(0);
        std::vector<Facade::geo_point_t> main_pts_4;
        for (const auto& wcp : wcp_list) main_pts_4.push_back(wcp.point);
        create_segment_point_cloud(sg1, main_pts_4, dv, "main");

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
    
    const auto& seg_front_point = seg_wcpts.front().point;
    const auto& seg_back_point  = seg_wcpts.back().point;
    bool flag_start = false;
    double dis_front = ray_length(Ray{vertex->wcpt().point, seg_front_point});
    double dis_back = ray_length(Ray{vertex->wcpt().point, seg_back_point});
    
    if (dis_front < dis_back) {
        flag_start = true;
    }
    
    // Step 3: Build list of points to test (from vertex end, excluding endpoints)
    std::vector<Facade::geo_point_t> pts_to_be_tested;
    pts_to_be_tested.reserve(seg_fits.size());
    
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
                        if (n_bad > 0) break;
                    }
                }
            }
            if (n_bad > 0) break;
        }
        
        if (n_bad == 0) max_bin = i;
    }
    
    // Step 5: Update segment and vertex if good point found
    if (!cluster.has_pc("steiner_pc")) return false;
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
        if (flag_start && ray_length(Ray{vtx_new_point, seg_back_point}) < 0.01 * units::cm) {
            max_bin--;
            continue;
        }
        if (!flag_start && ray_length(Ray{vtx_new_point, seg_front_point}) < 0.01 * units::cm) {
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
            double dis_limit = ray_length(Ray{vtx_new_point, seg_back_point});
            std::vector<Facade::geo_point_t> tmp_path;
            for (int idx = seg_wcpts.size() - 1; idx >= 0; idx--) {
                double dis = ray_length(Ray{seg_wcpts[idx].point, seg_back_point});
                if (dis < dis_limit) {
                    tmp_path.push_back(seg_wcpts[idx].point);
                }
            }
            if (tmp_path.size() > 1 && ray_length(Ray{tmp_path.back(), vtx_new_point}) < 0.01 * units::cm) {
                tmp_path.pop_back();
            }
            tmp_path.push_back(vtx_new_point);
            std::reverse(tmp_path.begin(), tmp_path.end());
            new_path = std::move(tmp_path);
        } else {
            // Keep points from front to new vertex
            double dis_limit = ray_length(Ray{vtx_new_point, seg_front_point});
            for (size_t idx = 0; idx < seg_wcpts.size(); idx++) {
                double dis = ray_length(Ray{seg_wcpts[idx].point, seg_front_point});
                if (dis < dis_limit) {
                    new_path.push_back(seg_wcpts[idx].point);
                }
            }
            if (new_path.size() > 1 ) {
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
        track_fitter.do_multi_tracking(true, true, false, false, false, &cluster);
        
        flag = true;
        break;
    }
    
    return flag;
}

void PatternAlgorithms::examine_segment(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    // Step 1: Examine short segments with multiple connections at both ends
    std::vector<SegmentPtr> segments_to_examine;

    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;
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
    // Use a vector in insertion order for deterministic iteration
    std::vector<VertexPtr> all_vertices;
    for (const auto& vd : ordered_nodes(graph)) {
        VertexPtr vtx = graph[vd].vertex;
        if (vtx && vtx->cluster() == &cluster) {
            all_vertices.push_back(vtx);
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
                        
                        all_vertices.erase(std::remove(all_vertices.begin(), all_vertices.end(), vtx2), all_vertices.end());
                        flag_merge = true;
                        break;
                    }
                }
            }
            if (flag_merge) break;
        }
    }
    
    // Step 3: Remove duplicate segments (same endpoints)
    std::vector<SegmentPtr> segments_to_be_removed;

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
        
        // Build canonical key set for O(E) duplicate detection using vertex pointers
        using SegmentKey = std::string;
        std::unordered_set<SegmentKey> seen_keys;
        
        for (size_t i = 0; i < tmp_segments.size(); i++) {
            auto [v1_i, v2_i] = find_vertices(graph, tmp_segments[i]);
            if (!v1_i || !v2_i) continue;
            
            // Build canonical key: smaller pointer address first
            SegmentKey key = std::to_string(std::min(uintptr_t(v1_i.get()), uintptr_t(v2_i.get())))
                           + ":"
                           + std::to_string(std::max(uintptr_t(v1_i.get()), uintptr_t(v2_i.get())));
            
            if (seen_keys.count(key)) {
                // Duplicate: a segment with the same endpoint pair was already seen
                if (std::find(segments_to_be_removed.begin(), segments_to_be_removed.end(), tmp_segments[i]) == segments_to_be_removed.end()) {
                    segments_to_be_removed.push_back(tmp_segments[i]);
                }
            } else {
                seen_keys.insert(key);
            }
        }
    }
    
    // Remove duplicate segments
    for (auto sg : segments_to_be_removed) {
        remove_segment(graph, sg);
    }
    
    // Step 4: Remove isolated vertices (no connections)
    std::vector<VertexPtr> vertices_to_be_removed;
    for (const auto& vd : ordered_nodes(graph)) {
        VertexPtr vtx = graph[vd].vertex;
        if (vtx && vtx->cluster() == &cluster) {
            if (vtx->descriptor_valid() && boost::degree(vd, graph) == 0) {
                vertices_to_be_removed.push_back(vtx);
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
    
    // Time normalization factor.
    // convert_3Dpoint_time_ch returns tind = round(time/tick) in raw ADC tick units,
    // and Blob::slice_index is also in raw ADC tick units (see Facade_Blob.h comment
    // "unit: tick").  For a typical blob that spans exactly one readout bin:
    //   ntime_ticks = nrebin  (e.g. 4 for SBND/MicroBooNE)
    // This makes tind/ntime_ticks equivalent to the prototype's
    //   offset_t + slope_xt * x  (time-slice index, in units comparable to wire pitch)
    // which is what the 2D-distance threshold of 2.5 is calibrated to.
    // Assumption: children()[0] spans exactly one readout bin.  This holds for all
    // standard wire-cell chunked blobs.  If the first blob unusually spans multiple
    // readout bins (ntime_ticks > nrebin), the threshold comparison becomes too loose
    // by that factor — safe but imprecise.
    auto first_blob = cluster.children()[0];
    int ntime_ticks = first_blob->slice_index_max() - first_blob->slice_index_min();
    if (ntime_ticks <= 0) ntime_ticks = 1;  // guard against degenerate blob

    // Convert vertices to U/V/W/T coordinates.
    // convert_3Dpoint_time_ch returns (tind, wind): tind is view-independent
    // (drift2time depends only on point[0] = x, not on the plane index), so tind
    // is identical for all three pind values for the same point.  We read it once
    // (pind=0) and discard the time component from the pind=1,2 calls.
    auto [v1_t_raw, v1_u_ch] = cluster.grouping()->convert_3Dpoint_time_ch(v1_p, v1_apa, v1_face, 0);
    auto [v1_t1,    v1_v_ch] = cluster.grouping()->convert_3Dpoint_time_ch(v1_p, v1_apa, v1_face, 1);
    auto [v1_t2,    v1_w_ch] = cluster.grouping()->convert_3Dpoint_time_ch(v1_p, v1_apa, v1_face, 2);
    (void)v1_t1; (void)v1_t2;  // tind is view-independent; only wire channels differ

    auto [v2_t_raw, v2_u_ch] = cluster.grouping()->convert_3Dpoint_time_ch(v2_p, v2_apa, v2_face, 0);
    auto [v2_t1,    v2_v_ch] = cluster.grouping()->convert_3Dpoint_time_ch(v2_p, v2_apa, v2_face, 1);
    auto [v2_t2,    v2_w_ch] = cluster.grouping()->convert_3Dpoint_time_ch(v2_p, v2_apa, v2_face, 2);
    (void)v2_t1; (void)v2_t2;

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
                
    const auto& seg_fits = sg->fits();

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
                
                // Find point on sg1 that's approximately 9 units away from v1.
                // p_t_raw from convert_3Dpoint_time_ch is view-independent (tind depends
                // only on x / drift_speed), so passing pind here is harmless — p_t_raw
                // equals what pind=0 would return.  p_wire_ch is the plane-correct channel.
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

bool PatternAlgorithms::examine_vertices_1(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex){
    bool flag_continue = false;
    
    VertexPtr v1 = nullptr;
    VertexPtr v2 = nullptr;
    VertexPtr v3 = nullptr;
    
    // Iterate in insertion order for deterministic results
    for (const auto& vd_cur : ordered_nodes(graph)) {
        VertexPtr vtx = graph[vd_cur].vertex;
        if (!vtx || vtx->cluster() != &cluster) continue;

        // Check if vertex has exactly 2 connections (potential candidate)
        if (boost::degree(vd_cur, graph) != 2) continue;

        // Get the two connected segments and cache their neighbor vertices,
        // avoiding redundant find_other_vertex calls later.
        std::vector<std::pair<SegmentPtr, VertexPtr>> seg_vtx_pairs;
        auto [ebegin, eend] = boost::out_edges(vd_cur, graph);
        for (auto eit = ebegin; eit != eend; ++eit) {
            SegmentPtr seg = graph[*eit].segment;
            if (seg) seg_vtx_pairs.emplace_back(seg, graph[boost::target(*eit, graph)].vertex);
        }
        
        if (seg_vtx_pairs.size() != 2) continue;
        
        // Check each segment
        for (auto& [seg, vtx1] : seg_vtx_pairs) {
            // Only consider short segments (<4cm)
            if (segment_track_length(seg) > 4.0 * units::cm) continue;
            
            // vtx1 is the cached neighbor vertex — no graph traversal needed
            if (!vtx1) continue;
            
            // The other vertex must have at least 2 connections
            if (!vtx1->descriptor_valid()) continue;
            auto vd1 = vtx1->get_descriptor();
            if (boost::degree(vd1, graph) < 2) continue;
            
            // Check if these two vertices represent the same physical point
            if (examine_vertices_1p(graph, vtx, vtx1, track_fitter, dv)) {
                v1 = vtx;
                v2 = vtx1;
                
                // Find v3: neighbor already cached in the other pair entry
                for (auto& [seg2, vtx_other] : seg_vtx_pairs) {
                    if (seg2 == seg) continue;
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
        if (v1 == main_vertex) {
          return false;
        }
        // Find the segments to be removed
        SegmentPtr sg = find_segment(graph, v1, v2);  // segment between v1 and v2
        SegmentPtr sg1 = find_segment(graph, v1, v3); // segment between v1 and v3
        
        if (!sg || !sg1) {
            return false;
        }
        
        // Create new segment from v3 to v2 using Steiner graph shortest path
        auto path_points = do_rough_path(cluster, v3->wcpt().point, v2->wcpt().point);
        
        if (path_points.size() < 2) {
            return flag_continue;
        }
        
        // Create the new segment
        auto sg2 = create_segment_for_cluster(cluster, dv, path_points, 0);
        if (!sg2) {
            return flag_continue;
        }
        
        // std::cout << "Cluster: " << cluster.ident() << " Merge Vertices Type I " 
        //           << "combining two segments into new segment" << std::endl;
        
        // Add new segment to graph
        add_segment(graph, sg2, v2, v3);
        
        // Remove old segments and vertex
        remove_segment(graph, sg);
        remove_segment(graph, sg1);
        remove_vertex(graph, v1);
        
        // Update tracking
        track_fitter.do_multi_tracking(true, true, false, false, false, &cluster);
    }
    
    return flag_continue;
}

bool PatternAlgorithms::examine_vertices_2(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex){
    bool flag_continue = false;
    
    VertexPtr v1 = nullptr;
    VertexPtr v2 = nullptr;
    SegmentPtr sg = nullptr;
    
    // Iterate in insertion order for deterministic results
    for (const auto& ed_cur : ordered_edges(graph)) {
        SegmentPtr segment = graph[ed_cur].segment;
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
        if (v2!= main_vertex){
            // std::cout << "Cluster: " << cluster.ident() << " Merge Vertices Type II" << std::endl;
            
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
            
            // Clean up isolated vertices (vertices with no connections) that belong
            // to this cluster.  Must filter by cluster: the graph may contain vertices
            // from other clusters, and we must not touch them here.
            std::vector<VertexPtr> isolated_vertices;
            for (const auto& vd : ordered_nodes(graph)) {
                if (boost::degree(vd, graph) == 0) {
                    VertexPtr vtx = graph[vd].vertex;
                    if (vtx && vtx->cluster() == &cluster) {
                        isolated_vertices.push_back(vtx);
                    }
                }
            }

            for (auto vtx : isolated_vertices) {
                remove_vertex(graph, vtx);
            }
            
            // Update tracking
            track_fitter.do_multi_tracking(true, true, false, false, false, &cluster);
        }else{
            flag_continue = false;
        }
    }
    
    return flag_continue;
}


bool PatternAlgorithms::examine_vertices_4p(Graph&graph, VertexPtr v1, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    // Find the segment between v1 and v2.
    // sg1 may be nullptr if the direct edge was already removed by an earlier step in
    // the same examine_vertices_4 pass.  When sg1 is nullptr the loop below skips
    // nothing (sg != nullptr is true for all surviving segments), so all of v1's
    // remaining segments are tested against v2 — this is conservative but correct.
    SegmentPtr sg1 = find_segment(graph, v1, v2);
    
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
        
        // Get segment points (fit points, matching prototype's get_point_vec())
        const auto& pts = sg->fits();
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

bool PatternAlgorithms::examine_vertices_4(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex){
    bool flag_continue = false;
    
    // Drift direction (X direction)
    Facade::geo_vector_t drift_dir_abs(1, 0, 0);
    
    if (!cluster.has_pc("steiner_pc")) return false;

    // Get steiner point cloud for later use
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
    
    // Snapshot all segments before iterating — the loop body calls remove_segment/add_segment
    // which would invalidate boost::edges iterators (even with listS storage the removed
    // edge's list node is freed, making the live iterator dangle on the next ++eit).
    std::vector<SegmentPtr> all_segments;
    for (auto [eit, eend] = boost::edges(graph); eit != eend; ++eit) {
        SegmentPtr s = graph[*eit].segment;
        if (s && s->cluster() == &cluster) all_segments.push_back(s);
    }

    for (auto& sg : all_segments) {
        if (!sg->descriptor_valid()) continue; // may have been removed in an earlier iteration
        
        const auto& pts = sg->fits();
        if (pts.size() < 2) continue;

        // Get vertices
        auto pair_vertices = find_vertices(graph, sg);
        VertexPtr v1 = pair_vertices.first;
        VertexPtr v2 = pair_vertices.second;
        if (!v1 || !v2) continue;

        // Calculate segment direction (using fit points, matching prototype's get_point_vec())
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
            
            if (boost::degree(vd1, graph) >= 2 && examine_vertices_4p(graph, v1, v2, track_fitter, dv) && v1 != main_vertex) {
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
                    
                    // Combine with rest of segment.
                    // The prototype trims old_list by WCPoint::index equality, which
                    // assumes min_wcp.index appears exactly once in vec_wcps.  Here we
                    // use a 0.01 cm position tolerance instead, which is safer when
                    // WCPoint::index is not reliably populated.  If min_wcp were somehow
                    // absent from old_list (should not happen for well-formed Steiner data),
                    // the while-loop would drain old_list to empty — the empty() guard
                    // ensures the subsequent pop_front/pop_back is a no-op.
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
                    std::vector<Facade::geo_point_t> main_pts_ev1;
                    for (const auto& wcp : new_wcpts) main_pts_ev1.push_back(wcp.point);
                    create_segment_point_cloud(sg1, main_pts_ev1, dv, "main");

                    // Find other vertex and update connection
                    VertexPtr v3 = find_other_vertex(graph, sg1, v1);
                    if (v3) {
                        remove_segment(graph, sg1);
                        add_segment(graph, sg1, v2, v3);
                    }
                }
                
                // Remove sg first, then v1 (remove_vertex implicitly removes adjacent edges,
                // so sg must be removed before v1 to avoid double-free on the edge descriptor)
                remove_segment(graph, sg);
                remove_vertex(graph, v1);
                
                flag_continue = true;
                // std::cout << "Cluster: " << cluster.ident() << " Merge Vertices Type III" << std::endl;
                track_fitter.do_multi_tracking(true, true, false, false, false, &cluster);
                break;
                
            } else if (boost::degree(vd2, graph) >= 2 && examine_vertices_4p(graph, v2, v1, track_fitter, dv) && v2 != main_vertex) {
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
                    
                    // Combine with rest of segment.
                    // The prototype trims old_list by WCPoint::index equality, which
                    // assumes min_wcp.index appears exactly once in vec_wcps.  Here we
                    // use a 0.01 cm position tolerance instead, which is safer when
                    // WCPoint::index is not reliably populated.  If min_wcp were somehow
                    // absent from old_list (should not happen for well-formed Steiner data),
                    // the while-loop would drain old_list to empty — the empty() guard
                    // ensures the subsequent pop_front/pop_back is a no-op.
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
                    std::vector<Facade::geo_point_t> main_pts_ev2;
                    for (const auto& wcp : new_wcpts) main_pts_ev2.push_back(wcp.point);
                    create_segment_point_cloud(sg1, main_pts_ev2, dv, "main");

                    // Find other vertex and update connection
                    VertexPtr v3 = find_other_vertex(graph, sg1, v2);
                    if (v3) {
                        remove_segment(graph, sg1);
                        add_segment(graph, sg1, v1, v3);
                    }
                }
                
                // Remove sg first, then v2 (remove_vertex implicitly removes adjacent edges,
                // so sg must be removed before v2 to avoid double-free on the edge descriptor)
                remove_segment(graph, sg);
                remove_vertex(graph, v2);
                
                flag_continue = true;
                // std::cout << "Cluster: " << cluster.ident() << " Merge Vertices Type III" << std::endl;
                track_fitter.do_multi_tracking(true, true, false, false, false, &cluster);
                break;
            }
        }
        
        if (flag_continue) break;
    }
    
    return flag_continue;
}

void PatternAlgorithms::examine_vertices(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex){
    bool flag_continue = true;
    
    while (flag_continue) {
        flag_continue = false;
        
        // Examine and clean up segment topology
        examine_segment(graph, cluster, track_fitter, dv);
        
        // Merge vertex if the kink is not at right location (Type I)
        flag_continue = flag_continue || examine_vertices_1(graph, cluster, track_fitter, dv, main_vertex);
        
        // Count only vertices belonging to this cluster (matches prototype's
        // find_vertices(cluster).size() check — graph may hold other clusters' vertices).
        size_t num_cluster_vertices = 0;
        for (const auto& vd : ordered_nodes(graph)) {
            VertexPtr vtx = graph[vd].vertex;
            if (vtx && vtx->cluster() == &cluster) num_cluster_vertices++;
        }

        // Merge vertices if they are too close (Type II) - only if cluster has more than 2 vertices
        if (num_cluster_vertices > 2) {
            flag_continue = flag_continue || examine_vertices_2(graph, cluster, track_fitter, dv, main_vertex);
        }
        
        // Merge vertices if they are reasonably close (Type III/IV)
        flag_continue = flag_continue || examine_vertices_4(graph, cluster, track_fitter, dv, main_vertex);
    }
}

void PatternAlgorithms::examine_partial_identical_segments(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    bool flag_continue = true;
    
    
    
    while (flag_continue) {
        flag_continue = false;
        
        // Iterate in insertion order for deterministic results
        for (const auto& vd_cur : ordered_nodes(graph)) {
            VertexPtr vtx = graph[vd_cur].vertex;
            if (!vtx || vtx->cluster() != &cluster) continue;

            // Only process vertices with more than 2 connections
            size_t degree = boost::degree(vd_cur, graph);
            if (degree <= 2) continue;

            // Find pair of segments with maximum overlap distance
            SegmentPtr max_sg1 = nullptr;
            SegmentPtr max_sg2 = nullptr;
            double max_dis = 0;
            Facade::geo_point_t max_point;

            // Collect out-edges in insertion order for deterministic pair selection
            std::vector<edge_descriptor> out_edges_sorted;
            {
                auto [oe_begin, oe_end] = boost::out_edges(vd_cur, graph);
                out_edges_sorted.assign(oe_begin, oe_end);
                std::sort(out_edges_sorted.begin(), out_edges_sorted.end(),
                          [&graph](const edge_descriptor& a, const edge_descriptor& b) {
                              return graph[a].index < graph[b].index;
                          });
            }
            for (auto eit1 = out_edges_sorted.begin(); eit1 != out_edges_sorted.end(); ++eit1) {
                SegmentPtr sg1 = graph[*eit1].segment;
                if (!sg1) continue;

                // Use fit points, matching prototype's get_point_vec()
                const auto& pts_1 = sg1->fits();
                if (pts_1.empty()) continue;

                // Order points from vertex outward
                std::vector<Facade::geo_point_t> test_pts;
                bool front_is_vtx = (ray_length(Ray{pts_1.front().point, vtx->wcpt().point}) <
                                    ray_length(Ray{pts_1.back().point, vtx->wcpt().point}));

                if (front_is_vtx) {
                    for (const auto& pt : pts_1) {
                        test_pts.push_back(pt.point);
                    }
                } else {
                    for (auto it = pts_1.rbegin(); it != pts_1.rend(); ++it) {
                        test_pts.push_back(it->point);
                    }
                }

                if (test_pts.empty()) continue;

                // Compare with other segments
                for (auto eit2 = std::next(eit1); eit2 != out_edges_sorted.end(); ++eit2) {
                    SegmentPtr sg2 = graph[*eit2].segment;
                    if (!sg2) continue;

                    // Check overlap along sg1
                    for (size_t k = 0; k < test_pts.size(); k++) {
                        auto [closest_dis, closest_pt] = segment_get_closest_point(sg2, test_pts[k], "fit");

                        if (closest_dis < 0.3 * units::cm) {
                            // Get position for vertex
                            Facade::geo_point_t vtx_point = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
                            double dis = ray_length(Ray{test_pts[k], vtx_point});

                            if (dis > max_dis) {
                                max_dis = dis;
                                max_point = test_pts[k];
                                max_sg1 = sg1;
                                max_sg2 = sg2;
                            }
                        } else {
                            break; // No longer overlapping
                        }
                    }
                }
            }
            
            // If significant overlap found (>5cm), split the vertex
            if (max_dis > 5.0 * units::cm) {
                // Find closest existing vertex to max_point
                double min_dis = 1e9;
                VertexPtr min_vertex = nullptr;
                
                auto [vbegin2, vend2] = boost::vertices(graph);
                for (auto vit2 = vbegin2; vit2 != vend2; ++vit2) {
                    VertexPtr vtx1 = graph[*vit2].vertex;
                    if (!vtx1 || vtx1->cluster() != &cluster) continue;
                    
                    Facade::geo_point_t vtx1_point = vtx1->fit().valid() ? vtx1->fit().point : vtx1->wcpt().point;
                    double dis = ray_length(Ray{max_point, vtx1_point});
                    
                    if (dis < min_dis) {
                        min_dis = dis;
                        min_vertex = vtx1;
                    }
                }
                
                if (min_dis < 0.3 * units::cm) {
                    // Merge to existing vertex
                    
                    // Check if there's already a segment between min_vertex and vtx
                    SegmentPtr good_segment = find_segment(graph, min_vertex, vtx);
                    
                    // If no existing segment, create one
                    if (!good_segment) {
                        auto path_points = do_rough_path(cluster, min_vertex->wcpt().point, vtx->wcpt().point);
                        if (path_points.size() >= 2) {
                            auto sg3 = create_segment_for_cluster(cluster, dv, path_points, 0);
                            if (sg3) {
                                add_segment(graph, sg3, min_vertex, vtx);
                            }
                        }
                    }
                    
                    // Reconnect max_sg1 to min_vertex
                    if (max_sg1 && max_sg1 != good_segment) {
                        VertexPtr tmp_vtx = find_other_vertex(graph, max_sg1, vtx);
                        if (tmp_vtx && tmp_vtx != min_vertex) {
                            auto path_points = do_rough_path(cluster, min_vertex->wcpt().point, tmp_vtx->wcpt().point);
                            if (path_points.size() >= 2) {
                                std::vector<WCPoint> new_wcpts;
                                for (const auto& p : path_points) {
                                    WCPoint wcp;
                                    wcp.point = p;
                                    new_wcpts.push_back(wcp);
                                }
                                max_sg1->wcpts(new_wcpts);
                                std::vector<Facade::geo_point_t> main_pts_sg1a;
                                for (const auto& wcp : new_wcpts) main_pts_sg1a.push_back(wcp.point);
                                create_segment_point_cloud(max_sg1, main_pts_sg1a, dv, "main");

                                remove_segment(graph, max_sg1);
                                add_segment(graph, max_sg1, min_vertex, tmp_vtx);
                            }
                        } else if (tmp_vtx == min_vertex) {
                            remove_segment(graph, max_sg1);
                        }
                    }
                    
                    // Reconnect max_sg2 to min_vertex
                    if (max_sg2 && max_sg2 != good_segment) {
                        VertexPtr tmp_vtx = find_other_vertex(graph, max_sg2, vtx);
                        if (tmp_vtx && tmp_vtx != min_vertex) {
                            auto path_points = do_rough_path(cluster, min_vertex->wcpt().point, tmp_vtx->wcpt().point);
                            if (path_points.size() >= 2) {
                                std::vector<WCPoint> new_wcpts;
                                for (const auto& p : path_points) {
                                    WCPoint wcp;
                                    wcp.point = p;
                                    new_wcpts.push_back(wcp);
                                }
                                max_sg2->wcpts(new_wcpts);
                                std::vector<Facade::geo_point_t> main_pts_sg2a;
                                for (const auto& wcp : new_wcpts) main_pts_sg2a.push_back(wcp.point);
                                create_segment_point_cloud(max_sg2, main_pts_sg2a, dv, "main");

                                remove_segment(graph, max_sg2);
                                add_segment(graph, max_sg2, min_vertex, tmp_vtx);
                            }
                        } else if (tmp_vtx == min_vertex) {
                            remove_segment(graph, max_sg2);
                        }
                    }
                    
                    track_fitter.do_multi_tracking(true, true, false, false, false, &cluster);
                    
                } else {
                    // Create new vertex at split point
                    if (!cluster.has_pc("steiner_pc")) continue;
                    // Get steiner point cloud
                    const auto& steiner_pc = cluster.get_pc("steiner_pc");
                    const auto& coords = cluster.get_default_scope().coords;
                    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
                    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
                    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();

                    // Find closest steiner point
                    auto knn_results = cluster.kd_steiner_knn(1, max_point, "steiner_pc");
                    if (!knn_results.empty()) {
                        size_t idx = knn_results[0].first;
                        Facade::geo_point_t vtx_new_point(x_coords[idx], y_coords[idx], z_coords[idx]);
                        
                        // Create new vertex
                        auto vtx2 = make_vertex(graph);
                        WCPoint new_wcp;
                        new_wcp.point = vtx_new_point;
                        vtx2->wcpt(new_wcp).cluster(&cluster);
                        
                        // Create segment between new vertex and original vertex
                        auto path_points = do_rough_path(cluster, vtx_new_point, vtx->wcpt().point);
                        if (path_points.size() >= 2) {
                            auto sg3 = create_segment_for_cluster(cluster, dv, path_points, 0);
                            if (sg3) {
                                add_segment(graph, sg3, vtx2, vtx);
                            }
                        }
                        
                        // Reconnect max_sg1 to new vertex
                        if (max_sg1) {
                            VertexPtr tmp_vtx = find_other_vertex(graph, max_sg1, vtx);
                            if (tmp_vtx) {
                                auto path_points1 = do_rough_path(cluster, vtx_new_point, tmp_vtx->wcpt().point);
                                if (path_points1.size() >= 2) {
                                    std::vector<WCPoint> new_wcpts;
                                    for (const auto& p : path_points1) {
                                        WCPoint wcp;
                                        wcp.point = p;
                                        new_wcpts.push_back(wcp);
                                    }
                                    max_sg1->wcpts(new_wcpts);
                                    std::vector<Facade::geo_point_t> main_pts_sg1b;
                                    for (const auto& wcp : new_wcpts) main_pts_sg1b.push_back(wcp.point);
                                    create_segment_point_cloud(max_sg1, main_pts_sg1b, dv, "main");

                                    remove_segment(graph, max_sg1);
                                    add_segment(graph, max_sg1, vtx2, tmp_vtx);
                                }
                            }
                        }
                        
                        // Reconnect max_sg2 to new vertex
                        if (max_sg2) {
                            VertexPtr tmp_vtx = find_other_vertex(graph, max_sg2, vtx);
                            if (tmp_vtx) {
                                auto path_points2 = do_rough_path(cluster, vtx_new_point, tmp_vtx->wcpt().point);
                                if (path_points2.size() >= 2) {
                                    std::vector<WCPoint> new_wcpts;
                                    for (const auto& p : path_points2) {
                                        WCPoint wcp;
                                        wcp.point = p;
                                        new_wcpts.push_back(wcp);
                                    }
                                    max_sg2->wcpts(new_wcpts);
                                    std::vector<Facade::geo_point_t> main_pts_sg2b;
                                    for (const auto& wcp : new_wcpts) main_pts_sg2b.push_back(wcp.point);
                                    create_segment_point_cloud(max_sg2, main_pts_sg2b, dv, "main");

                                    remove_segment(graph, max_sg2);
                                    add_segment(graph, max_sg2, vtx2, tmp_vtx);
                                }
                            }
                        }
                        
                        track_fitter.do_multi_tracking(true, true, false, false, false, &cluster);
                    }
                }
                
                flag_continue = true;
            }
            
            if (flag_continue) break;
        }
    }
}

Facade::geo_point_t PatternAlgorithms::get_local_extension(Facade::Cluster& cluster, const Facade::geo_point_t& wcp){
    // Determine which point cloud to use
    std::string pc_name = "steiner_pc";
    
    // Get local direction using Hough transform
    Facade::geo_vector_t dir1 = cluster.vhough_transform(wcp, 10.0 * units::cm);
    dir1 = dir1 * (-1.0);  // Reverse direction
    
    // Drift direction
    Facade::geo_vector_t drift_dir(1, 0, 0);
    
    // Calculate angle with drift direction (in degrees)
    double dir1_mag = dir1.magnitude();
    if (dir1_mag == 0) return wcp;
    
    double cos_angle = drift_dir.dot(dir1) / dir1_mag;
    cos_angle = std::max(-1.0, std::min(1.0, cos_angle));  // Clamp to [-1, 1]
    double angle = std::acos(cos_angle) * 180.0 / M_PI;
    
    // If angle is close to perpendicular to drift (90° ± 7.5°), return original point
    if (std::fabs(angle - 90.0) < 7.5) {
        return wcp;
    }
    
    // Get nearby points within 10 cm radius
    auto kd_results = cluster.kd_steiner_radius(10.0 * units::cm, wcp, pc_name);
    
    if (kd_results.empty()) {
        return wcp;
    }
    
    // Get point coordinates
    const auto& pc = cluster.get_pc(pc_name);
    const auto& coords = cluster.get_default_scope().coords;
    const auto& x_coords = pc.get(coords.at(0))->elements<double>();
    const auto& y_coords = pc.get(coords.at(1))->elements<double>();
    const auto& z_coords = pc.get(coords.at(2))->elements<double>();
    
    // Find point with maximum projection along dir1
    double max_val = 0;
    geo_point_t result = wcp;
    
    for (const auto& [idx, dist] : kd_results) {
        Facade::geo_point_t pt(x_coords[idx], y_coords[idx], z_coords[idx]);
        
        // Calculate projection along dir1
        double val = dir1.x() * (pt.x() - wcp.x()) + 
                     dir1.y() * (pt.y() - wcp.y()) + 
                     dir1.z() * (pt.z() - wcp.z());
        
        if (val > max_val) {
            max_val = val;
            result = pt;
        }
    }
    
    return result;
}

void PatternAlgorithms::examine_vertices_3(Graph& graph, Facade::Cluster& main_cluster, std::pair<VertexPtr, VertexPtr> main_cluster_initial_pair_vertices, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    // Examine main_cluster_initial_pair_vertices to see if they need extension
    std::vector<VertexPtr> temp_vertices;
    if (main_cluster_initial_pair_vertices.first) temp_vertices.push_back(main_cluster_initial_pair_vertices.first);
    if (main_cluster_initial_pair_vertices.second) temp_vertices.push_back(main_cluster_initial_pair_vertices.second);
    
    bool flag_refit = false;
    
    for (size_t i = 0; i < temp_vertices.size(); i++) {
        VertexPtr vtx = temp_vertices[i];
        
        // Check if vertex is still in graph
        if (!vtx->descriptor_valid()) continue;
        auto vd = vtx->get_descriptor();
        
        // Only process vertices with exactly one connected segment
        if (boost::degree(vd, graph) != 1) continue;
        
        // Get the single connected segment
        SegmentPtr sg = nullptr;
        auto [ebegin, eend] = boost::out_edges(vd, graph);
        for (auto eit = ebegin; eit != eend; ++eit) {
            sg = graph[*eit].segment;
            break;
        }
        if (!sg) continue;
        
        const auto& wcps = sg->wcpts();
        if (wcps.empty()) continue;
        
        // Determine which end of segment connects to vertex
        bool flag_start = (ray_length(Ray{wcps.front().point, vtx->wcpt().point}) <
                          ray_length(Ray{wcps.back().point, vtx->wcpt().point}));
        
        // Get the other end of the segment
        WCPoint wcp2 = flag_start ? wcps.back() : wcps.front();
        
        // Try to extend the vertex using get_local_extension
        auto wcp1_point = get_local_extension(main_cluster, vtx->wcpt().point);
        
        // Check if extension found a different point
        bool same_as_vtx = (ray_length(Ray{wcp1_point, vtx->wcpt().point}) < 0.01 * units::cm);
        bool same_as_wcp2 = (ray_length(Ray{wcp1_point, wcp2.point}) < 0.01 * units::cm);
        
        if (same_as_vtx || same_as_wcp2) continue;
        
        // Create new path from extended point to other end
        std::vector<Facade::geo_point_t> path_points;
        if (flag_start) {
            path_points = do_rough_path(main_cluster, wcp1_point, wcp2.point);
        } else {
            path_points = do_rough_path(main_cluster, wcp2.point, wcp1_point);
        }
        
        // Only update if new path is not too much longer than original
        if (path_points.size() > 0 && path_points.size() < wcps.size() * 2) {
            // Update vertex position
            vtx->wcpt().point = wcp1_point;
            
            // Update segment path
            std::vector<WCPoint> new_wcpts;
            for (const auto& p : path_points) {
                WCPoint wcp;
                wcp.point = p;
                new_wcpts.push_back(wcp);
            }
            sg->wcpts(new_wcpts);
            std::vector<Facade::geo_point_t> main_pts_ev3;
            for (const auto& wcp : new_wcpts) main_pts_ev3.push_back(wcp.point);
            create_segment_point_cloud(sg, main_pts_ev3, dv, "main");

            flag_refit = true;
        }
    }
    
    if (flag_refit) {
        track_fitter.do_multi_tracking(true, true, false, false, false, &main_cluster);
    }

    // Find and remove redundant short segments
    // Use vector in insertion order for deterministic iteration
    std::vector<SegmentPtr> segments_to_be_removed;

    for (const auto& ed : ordered_edges(graph)) {
        SegmentPtr sg = graph[ed].segment;
        if (!sg || sg->cluster() != &main_cluster) continue;
        
        auto pair_vertices = find_vertices(graph, sg);
        if (!pair_vertices.first || !pair_vertices.second) continue;
        
        // Check if either vertex has only one connection
        auto vd1 = pair_vertices.first->get_descriptor();
        auto vd2 = pair_vertices.second->get_descriptor();
        bool v1_single = (boost::degree(vd1, graph) == 1);
        bool v2_single = (boost::degree(vd2, graph) == 1);
        
        if (!v1_single && !v2_single) continue;
        
        // Check if segment is short
        double direct_length = segment_track_direct_length(sg);
        if (direct_length >= 5.0 * units::cm) continue;
        
        // Check if all points on this segment are close to other segments in 2D
        // (use fit points, matching prototype's get_point_vec()).
        //
        // Efficiency: two early-exit shortcuts keep the inner work small:
        //   (a) Per-point inner loop: once all three per-view minima are already
        //       below 0.6 cm no remaining segment can make the point "unique" —
        //       break immediately.
        //   (b) Per-segment outer loop: as soon as one unique point is found the
        //       segment will NOT be removed — no need to test remaining points.
        // Pre-snapshot the other-segment list once (outside the point loop) so
        // we avoid re-iterating boost::edges and re-filtering on every point.
        std::vector<SegmentPtr> other_segs;
        for (auto [e2b, e2e] = boost::edges(graph); e2b != e2e; ++e2b) {
            SegmentPtr sg1 = graph[*e2b].segment;
            if (sg1 && sg1 != sg) other_segs.push_back(sg1);
        }

        const auto& pts = sg->fits();
        int num_unique = 0;

        for (size_t i = 0; i < pts.size(); i++) {
            // Get APA and face for this point
            auto wpid = dv->contained_by(pts[i].point);
            if (wpid.apa() == -1 || wpid.face() == -1) continue;

            double min_u = 1e9;
            double min_v = 1e9;
            double min_w = 1e9;

            for (SegmentPtr sg1 : other_segs) {
                auto [dist_u, dist_v, dist_w] = segment_get_closest_2d_distances(sg1, pts[i].point, wpid.apa(), wpid.face(), "fit");

                if (dist_u < min_u) min_u = dist_u;
                if (dist_v < min_v) min_v = dist_v;
                if (dist_w < min_w) min_w = dist_w;

                // (a) Point already fully covered in all views — no need to
                //     check remaining segments.
                if (min_u <= 0.6 * units::cm &&
                    min_v <= 0.6 * units::cm &&
                    min_w <= 0.6 * units::cm) break;
            }

            // If point is far from all other segments in any view, it's unique
            if (min_u > 0.6 * units::cm || min_v > 0.6 * units::cm || min_w > 0.6 * units::cm) {
                num_unique++;
                break;  // (b) Segment has a unique point — it will not be removed.
            }
        }
        
        // If no unique points, mark for removal
        if (num_unique == 0) {
            segments_to_be_removed.push_back(sg);
        }
    }

    // Collect vertices that will need examination after removal
    // Use vector in insertion order for deterministic call order into examine_structure_4
    std::vector<VertexPtr> can_vertices;
    
    for (auto sg : segments_to_be_removed) {
        auto pair_vertices = find_vertices(graph, sg);
        
        if (pair_vertices.first && pair_vertices.first->descriptor_valid()) {
            auto vd1 = pair_vertices.first->get_descriptor();
            if (boost::degree(vd1, graph) > 1 &&
                std::find(can_vertices.begin(), can_vertices.end(), pair_vertices.first) == can_vertices.end()) {
                can_vertices.push_back(pair_vertices.first);
            }
        }

        if (pair_vertices.second && pair_vertices.second->descriptor_valid()) {
            auto vd2 = pair_vertices.second->get_descriptor();
            if (boost::degree(vd2, graph) > 1 &&
                std::find(can_vertices.begin(), can_vertices.end(), pair_vertices.second) == can_vertices.end()) {
                can_vertices.push_back(pair_vertices.second);
            }
        }
        
        remove_segment(graph, sg);
    }
    
    // Remove isolated vertices that belong to this cluster.
    // Must filter by cluster: the graph may hold vertices from other clusters
    // that happen to be degree-0 (pre-existing orphans); do not touch them here.
    bool flag_cont = true;
    while (flag_cont) {
        flag_cont = false;
        for (const auto& vd : ordered_nodes(graph)) {
            if (boost::degree(vd, graph) == 0) {
                VertexPtr vtx = graph[vd].vertex;
                if (vtx && vtx->cluster() == &main_cluster) {
                    remove_vertex(graph, vtx);
                    flag_cont = true;
                    break;
                }
            }
        }
    }
    
    // Examine remaining vertices with examine_structure_4
    for (auto vtx : can_vertices) {
        if (vtx->descriptor_valid()) {
            examine_structure_4(vtx, false, graph, main_cluster, track_fitter, dv);
        }
    }
    
    // Refit if segments were removed
    if (segments_to_be_removed.size() > 0) {
        track_fitter.do_multi_tracking(true, true, false, false, false, &main_cluster);
    }
}


        
bool PatternAlgorithms::examine_structure_final_1(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv) {
    // Merge two segments if a direct connection is better
    bool flag_update = false;
    bool flag_continue = true;
    
    // Get transform and grouping from track_fitter
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster.get_scope_transform(cluster.get_default_scope()));
    double cluster_t0 = cluster.get_cluster_t0();
    auto grouping = cluster.grouping();
    
    if (!transform || !grouping) {
        return false;
    }
    
    while (flag_continue) {
        flag_continue = false;
        
        // Iterate in insertion order for deterministic results
        for (const auto& vd_cur : ordered_nodes(graph)) {
            VertexPtr vtx = graph[vd_cur].vertex;

            // Skip if vertex doesn't belong to this cluster
            if (!vtx || vtx->cluster() != &cluster) continue;

            // Skip the main vertex
            if (vtx == main_vertex) continue;

            // Only consider vertices with exactly 2 connections
            auto vd = vtx->get_descriptor();
            if (boost::degree(vd, graph) != 2) continue;
            
            // Get the two segments connected to this vertex
            auto edge_range = boost::out_edges(vd, graph);
            auto eit = edge_range.first;
            SegmentPtr sg1 = graph[*eit].segment;
            ++eit;
            SegmentPtr sg2 = graph[*eit].segment;
            
            if (!sg1 || !sg2) continue;
            
            // Check if segments have identical endpoints (same start and end points)
            const auto& wcpts1 = sg1->wcpts();
            const auto& wcpts2 = sg2->wcpts();
            
            if (wcpts1.size() < 2 || wcpts2.size() < 2) continue;
            
            // Check if segments are identical (same endpoints)
            double dist_front_front = ray_length(Ray{wcpts1.front().point, wcpts2.front().point});
            double dist_back_back = ray_length(Ray{wcpts1.back().point, wcpts2.back().point});
            double dist_front_back = ray_length(Ray{wcpts1.front().point, wcpts2.back().point});
            double dist_back_front = ray_length(Ray{wcpts1.back().point, wcpts2.front().point});
            
            if ((dist_front_front < 0.1*units::cm && dist_back_back < 0.1*units::cm) ||
                (dist_front_back < 0.1*units::cm && dist_back_front < 0.1*units::cm)) {
                // Segments are identical, delete one
                s_log->debug("examine_structure_final_1: cluster {} removing duplicate segment at vtx ({:.2f},{:.2f},{:.2f})",
                    cluster.ident(), vtx->wcpt().point.x()/units::cm, vtx->wcpt().point.y()/units::cm, vtx->wcpt().point.z()/units::cm);
                remove_segment(graph, sg2);
                flag_update = true;
                flag_continue = true;
                break;
            }
            
            // Get segment lengths
            // double length1 = segment_track_length(sg1);
            // double length2 = segment_track_length(sg2);
            
            // Get the other vertices
            VertexPtr vtx1 = find_other_vertex(graph, sg1, vtx);
            VertexPtr vtx2 = find_other_vertex(graph, sg2, vtx);
            
            if (!vtx1 || !vtx2) continue;
            
            // Get start and end points (use fit if available, otherwise wcpt)
            Facade::geo_point_t start_p = vtx1->fit().valid() ? vtx1->fit().point : vtx1->wcpt().point;
            Facade::geo_point_t end_p = vtx2->fit().valid() ? vtx2->fit().point : vtx2->wcpt().point;
            
            // Check the straight line path
            double step_size = 0.6 * units::cm;
            double distance = ray_length(Ray{start_p, end_p});
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
            
            // If the straight line is better, replace the two segments with one new segment
            if (flag_replace) {
                // Use helper function to merge the two segments
                if (merge_two_segments_into_one(graph, sg1, vtx, sg2, dv)) {
                    s_log->debug("examine_structure_final_1: cluster {} merged two segments through vtx ({:.2f},{:.2f},{:.2f})",
                        cluster.ident(), vtx->wcpt().point.x()/units::cm, vtx->wcpt().point.y()/units::cm, vtx->wcpt().point.z()/units::cm);
                    flag_update = true;
                    flag_continue = true;
                    break;
                }
            }
        }
    } // while continue
    
    if (flag_update) {
        track_fitter.do_multi_tracking(true, true, false, false, false, &cluster);
    }
    
    return flag_update;
}

bool PatternAlgorithms::examine_structure_final_1p(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    bool flag_update = false;
    
    // Check if main_vertex has exactly 2 connected segments
    if (!main_vertex->descriptor_valid()) return flag_update;
    auto vd = main_vertex->get_descriptor();
    if (boost::degree(vd, graph) != 2) return flag_update;
    
    // Get the two segments connected to main_vertex
    auto edge_range = boost::out_edges(vd, graph);
    auto eit = edge_range.first;
    SegmentPtr sg1 = graph[*eit].segment;
    ++eit;
    SegmentPtr sg2 = graph[*eit].segment;
    
    if (!sg1 || !sg2) return flag_update;
    
    // Get main vertex position
    WireCell::Point main_vtx_point = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
    
    // Calculate direction vectors for both segments
    WireCell::Vector dir1 = segment_cal_dir_3vector(sg1, main_vtx_point, 15*units::cm);
    WireCell::Vector dir2 = segment_cal_dir_3vector(sg2, main_vtx_point, 15*units::cm);
    
    // Calculate angle between directions (in degrees)
    // Prototype: dir1.Angle(dir2)*180/pi — angle between the two direction vectors
    // Both dirs point away from main_vertex; collinear segments give ~180° (antiparallel)
    double angle = std::acos(dir1.dot(dir2) / (dir1.magnitude() * dir2.magnitude())) / 3.1415926 * 180.0;
    
    // Get segment lengths
    double length1 = segment_track_length(sg1);
    double length2 = segment_track_length(sg2);
    
    // Only proceed if segments are nearly collinear (angle > 175 degrees)
    if (angle > 175) {
        // Get transform and grouping for point validation
        const auto transform = track_fitter.get_pc_transforms()->pc_transform(
            cluster.get_scope_transform(cluster.get_default_scope()));
        // double cluster_t0 = cluster.get_cluster_t0();
        auto grouping = cluster.grouping();
        
        if (!transform || !grouping) return flag_update;
        
        if (length1 < 6*units::cm && length1 < length2) {
            // sg1 is short - merge it into sg2
            VertexPtr vtx = find_other_vertex(graph, sg1, main_vertex);
            if (!vtx) return flag_update;
            
            const auto& vec_wcps = sg2->wcpts();
            const auto& vec_wcps1 = sg1->wcpts();
            
            if (vec_wcps.empty() || vec_wcps1.empty()) return flag_update;
            
            // Determine which end of sg2 connects to main_vertex
            // Use wcpt().point (steiner node) not fit().point — segment wcpts are anchored at steiner nodes
            WireCell::Point main_wcpt_point = main_vertex->wcpt().point;
            bool flag_front = (ray_length(Ray{vec_wcps.front().point, main_wcpt_point}) < 0.01*units::cm);

            // Determine which end of sg1 connects to main_vertex
            bool flag_front1 = (ray_length(Ray{vec_wcps1.front().point, main_wcpt_point}) < 0.01*units::cm);

            // Create a list to merge the wcpts
            std::list<WCPoint> old_list;
            std::copy(vec_wcps.begin(), vec_wcps.end(), std::back_inserter(old_list));

            // Merge sg1 points into sg2 based on orientation
            if (flag_front && flag_front1) {
                // Both connect at front - add sg1 in order to front of old_list
                for (auto it1 = vec_wcps1.begin(); it1 != vec_wcps1.end(); it1++) {
                    if (ray_length(Ray{(*it1).point, old_list.front().point}) > 0.01*units::cm) {
                        old_list.push_front(*it1);
                    }
                }
            } else if (flag_front && (!flag_front1)) {
                // sg2 front connects, sg1 back connects - add sg1 in reverse to front of old_list
                for (auto it1 = vec_wcps1.rbegin(); it1 != vec_wcps1.rend(); it1++) {
                    if (ray_length(Ray{(*it1).point, old_list.front().point}) > 0.01*units::cm) {
                        old_list.push_front(*it1);
                    }
                }
            } else if ((!flag_front) && flag_front1) {
                // sg2 back connects, sg1 front connects - add sg1 in order to back of old_list
                for (auto it1 = vec_wcps1.begin(); it1 != vec_wcps1.end(); it1++) {
                    if (ray_length(Ray{(*it1).point, old_list.back().point}) > 0.01*units::cm) {
                        old_list.push_back(*it1);
                    }
                }
            } else if ((!flag_front) && (!flag_front1)) {
                // Both connect at back - add sg1 in reverse to back of old_list
                for (auto it1 = vec_wcps1.rbegin(); it1 != vec_wcps1.rend(); it1++) {
                    if (ray_length(Ray{(*it1).point, old_list.back().point}) > 0.01*units::cm) {
                        old_list.push_back(*it1);
                    }
                }
            }
            
            // Update sg2's wcpts with merged list
            std::vector<WCPoint> new_wcpts;
            new_wcpts.reserve(old_list.size());
            std::copy(std::begin(old_list), std::end(old_list), std::back_inserter(new_wcpts));
            sg2->wcpts(new_wcpts);
            std::vector<Facade::geo_point_t> main_pts_sg2m;
            for (const auto& wcp : new_wcpts) main_pts_sg2m.push_back(wcp.point);
            create_segment_point_cloud(sg2, main_pts_sg2m, dv, "main");

            // Update main_vertex to vtx's position
            WCPoint vtx_wcp = vtx->wcpt();
            main_vertex->wcpt(vtx_wcp);
            if (vtx->fit().valid()) {
                main_vertex->fit(vtx->fit());
            }
            
            // Reconnect all segments from vtx to main_vertex (except sg1)
            std::vector<SegmentPtr> vtx_segments;
            if (vtx->descriptor_valid()) {
                auto vtx_vd = vtx->get_descriptor();
                auto [vtx_ebegin, vtx_eend] = boost::out_edges(vtx_vd, graph);
                for (auto vtx_eit = vtx_ebegin; vtx_eit != vtx_eend; ++vtx_eit) {
                    SegmentPtr seg = graph[*vtx_eit].segment;
                    if (seg && seg != sg1) {
                        vtx_segments.push_back(seg);
                    }
                }
            }
            
            for (auto seg : vtx_segments) {
                VertexPtr other_vtx = find_other_vertex(graph, seg, vtx);
                if (other_vtx && other_vtx != main_vertex) {
                    remove_segment(graph, seg);
                    add_segment(graph, seg, main_vertex, other_vtx);
                }
            }
            
            // Delete sg1 and vtx
            remove_segment(graph, sg1);
            remove_vertex(graph, vtx);
            s_log->debug("examine_structure_final_1p: cluster {} merged short sg1 ({:.2f} cm) into sg2 at main_vtx ({:.2f},{:.2f},{:.2f})",
                cluster.ident(), length1/units::cm,
                main_vtx_point.x()/units::cm, main_vtx_point.y()/units::cm, main_vtx_point.z()/units::cm);
            flag_update = true;

        } else if (length2 < 6*units::cm && length2 < length1) {
            // sg2 is short - merge it into sg1
            VertexPtr vtx = find_other_vertex(graph, sg2, main_vertex);
            if (!vtx) return flag_update;
            
            const auto& vec_wcps = sg1->wcpts();
            const auto& vec_wcps1 = sg2->wcpts();
            
            if (vec_wcps.empty() || vec_wcps1.empty()) return flag_update;
            
            // Determine which end of sg1 connects to main_vertex
            // Use wcpt().point (steiner node) not fit().point — segment wcpts are anchored at steiner nodes
            WireCell::Point main_wcpt_point = main_vertex->wcpt().point;
            bool flag_front = (ray_length(Ray{vec_wcps.front().point, main_wcpt_point}) < 0.01*units::cm);

            // Determine which end of sg2 connects to main_vertex
            bool flag_front1 = (ray_length(Ray{vec_wcps1.front().point, main_wcpt_point}) < 0.01*units::cm);

            // Create a list to merge the wcpts
            std::list<WCPoint> old_list;
            std::copy(vec_wcps.begin(), vec_wcps.end(), std::back_inserter(old_list));

            // Merge sg2 points into sg1 based on orientation
            if (flag_front && flag_front1) {
                // Both connect at front - add sg2 in order to front of old_list
                for (auto it1 = vec_wcps1.begin(); it1 != vec_wcps1.end(); it1++) {
                    if (ray_length(Ray{(*it1).point, old_list.front().point}) > 0.01*units::cm) {
                        old_list.push_front(*it1);
                    }
                }
            } else if (flag_front && (!flag_front1)) {
                // sg1 front connects, sg2 back connects - add sg2 in reverse to front of old_list
                for (auto it1 = vec_wcps1.rbegin(); it1 != vec_wcps1.rend(); it1++) {
                    if (ray_length(Ray{(*it1).point, old_list.front().point}) > 0.01*units::cm) {
                        old_list.push_front(*it1);
                    }
                }
            } else if ((!flag_front) && flag_front1) {
                // sg1 back connects, sg2 front connects - add sg2 in order to back of old_list
                for (auto it1 = vec_wcps1.begin(); it1 != vec_wcps1.end(); it1++) {
                    if (ray_length(Ray{(*it1).point, old_list.back().point}) > 0.01*units::cm) {
                        old_list.push_back(*it1);
                    }
                }
            } else if ((!flag_front) && (!flag_front1)) {
                // Both connect at back - add sg2 in reverse to back of old_list
                for (auto it1 = vec_wcps1.rbegin(); it1 != vec_wcps1.rend(); it1++) {
                    if (ray_length(Ray{(*it1).point, old_list.back().point}) > 0.01*units::cm) {
                        old_list.push_back(*it1);
                    }
                }
            }
            
            // Update sg1's wcpts with merged list
            std::vector<WCPoint> new_wcpts;
            new_wcpts.reserve(old_list.size());
            std::copy(std::begin(old_list), std::end(old_list), std::back_inserter(new_wcpts));
            sg1->wcpts(new_wcpts);
            std::vector<Facade::geo_point_t> main_pts_sg1m;
            for (const auto& wcp : new_wcpts) main_pts_sg1m.push_back(wcp.point);
            create_segment_point_cloud(sg1, main_pts_sg1m, dv, "main");

            // Update main_vertex to vtx's position
            WCPoint vtx_wcp = vtx->wcpt();
            main_vertex->wcpt(vtx_wcp);
            if (vtx->fit().valid()) {
                main_vertex->fit(vtx->fit());
            }
            
            // Reconnect all segments from vtx to main_vertex (except sg2)
            std::vector<SegmentPtr> vtx_segments;
            if (vtx->descriptor_valid()) {
                auto vtx_vd = vtx->get_descriptor();
                auto [vtx_ebegin, vtx_eend] = boost::out_edges(vtx_vd, graph);
                for (auto vtx_eit = vtx_ebegin; vtx_eit != vtx_eend; ++vtx_eit) {
                    SegmentPtr seg = graph[*vtx_eit].segment;
                    if (seg && seg != sg2) {
                        vtx_segments.push_back(seg);
                    }
                }
            }
            
            for (auto seg : vtx_segments) {
                VertexPtr other_vtx = find_other_vertex(graph, seg, vtx);
                if (other_vtx && other_vtx != main_vertex) {
                    remove_segment(graph, seg);
                    add_segment(graph, seg, main_vertex, other_vtx);
                }
            }
            
            // Delete sg2 and vtx
            remove_segment(graph, sg2);
            remove_vertex(graph, vtx);
            s_log->debug("examine_structure_final_1p: cluster {} merged short sg2 ({:.2f} cm) into sg1 at main_vtx ({:.2f},{:.2f},{:.2f})",
                cluster.ident(), length2/units::cm,
                main_vtx_point.x()/units::cm, main_vtx_point.y()/units::cm, main_vtx_point.z()/units::cm);
            flag_update = true;
        }

        // If we updated, redo multi-tracking
        if (flag_update) {
            track_fitter.do_multi_tracking(true, true, false, false, false, &cluster);
        }
    }

    return flag_update;
}

bool PatternAlgorithms::examine_structure_final_2(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv) {
    bool flag_updated = false;
    
    if (!main_vertex || !main_vertex->descriptor_valid()) return flag_updated;
    
    // Get transform and grouping for point validation
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster.get_scope_transform(cluster.get_default_scope()));
    double cluster_t0 = cluster.get_cluster_t0();
    auto grouping = cluster.grouping();
    
    if (!transform || !grouping) return flag_updated;
    
    // Continue looping until no more updates
    bool flag_continue = true;
    while (flag_continue) {
        flag_continue = false;
        bool flag_update = false;
        
        auto main_vd = main_vertex->get_descriptor();
        
        // Loop over all segments connected to main_vertex
        auto [ebegin, eend] = boost::out_edges(main_vd, graph);
        for (auto eit = ebegin; eit != eend; ++eit) {
            SegmentPtr sg = graph[*eit].segment;
            if (!sg) continue;
            
            // Find the other vertex of this segment
            VertexPtr vtx1 = find_other_vertex(graph, sg, main_vertex);
            if (!vtx1 || !vtx1->descriptor_valid()) continue;
            
            // Skip if either vertex has only 1 connection
            auto vtx1_vd = vtx1->get_descriptor();
            if (boost::degree(vtx1_vd, graph) == 1 || boost::degree(main_vd, graph) == 1) continue;
            
            // Check distance between vertices
            WireCell::Point main_vtx_point = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
            WireCell::Point vtx1_point = vtx1->fit().valid() ? vtx1->fit().point : vtx1->wcpt().point;
            
            double dis = ray_length(Ray{main_vtx_point, vtx1_point});
            
            if (dis < 2.0*units::cm) {
                // Check if vtx1 can be merged into main_vertex
                flag_update = true;
                
                // Check all segments connected to vtx1 (except sg)
                auto [vtx1_ebegin, vtx1_eend] = boost::out_edges(vtx1_vd, graph);
                for (auto vtx1_eit = vtx1_ebegin; vtx1_eit != vtx1_eend; ++vtx1_eit) {
                    SegmentPtr sg1 = graph[*vtx1_eit].segment;
                    if (!sg1 || sg1 == sg) continue;
                    
                    const auto& sg1_wcpts = sg1->wcpts();
                    if (sg1_wcpts.empty()) continue;

                    // Determine which end of sg connects to vtx1
                    // Use vtx1->wcpt().point (steiner node), not vtx1_point (fit), to match sg wcpts
                    const auto& sg_wcpts = sg->wcpts();
                    if (sg_wcpts.empty()) continue;

                    WireCell::Point vtx1_wcpt_point = vtx1->wcpt().point;
                    bool flag_start = (ray_length(Ray{sg_wcpts.front().point, vtx1_wcpt_point}) < 0.01*units::cm);

                    // Find point at ~3cm from vtx1
                    WireCell::Point min_point = sg1_wcpts.front().point;
                    double min_dis = 1e9;
                    int min_index = 0;

                    for (size_t i = 0; i < sg1_wcpts.size(); i++) {
                        double dis = std::fabs(ray_length(Ray{sg1_wcpts.at(i).point, vtx1_point}) - 3*units::cm);
                        if (dis < min_dis) {
                            min_dis = dis;
                            min_point = sg1_wcpts.at(i).point;
                            min_index = i;
                        }
                    }

                    // Check connectivity from min_point to vtx1
                    bool flag_connect = true;
                    if (flag_start) {
                        for (int i = min_index; i >= 0; i--) {
                            auto test_wpid = dv->contained_by(sg1_wcpts.at(i).point);
                            if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                                auto temp_p_raw = transform->backward(sg1_wcpts.at(i).point, cluster_t0, test_wpid.face(), test_wpid.apa());
                                if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*units::cm, 0, 0)) {
                                    flag_connect = false;
                                    break;
                                }
                            }
                        }
                    } else {
                        for (size_t i = min_index; i < sg1_wcpts.size(); i++) {
                            auto test_wpid = dv->contained_by(sg1_wcpts.at(i).point);
                            if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                                auto temp_p_raw = transform->backward(sg1_wcpts.at(i).point, cluster_t0, test_wpid.face(), test_wpid.apa());
                                if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*units::cm, 0, 0)) {
                                    flag_connect = false;
                                    break;
                                }
                            }
                        }
                    }
                    
                    // Check path from min_point to main_vertex
                    if (flag_connect) {
                        double step_size = 0.3 * units::cm;
                        WireCell::Point start_p = min_point;
                        WireCell::Point end_p = main_vtx_point;
                        int ncount = std::round(ray_length(Ray{start_p, end_p}) / step_size);
                        int n_bad = 0;
                        
                        for (int i = 1; i < ncount; i++) {
                            WireCell::Point test_p(
                                start_p.x() + (end_p.x() - start_p.x()) / ncount * i,
                                start_p.y() + (end_p.y() - start_p.y()) / ncount * i,
                                start_p.z() + (end_p.z() - start_p.z()) / ncount * i
                            );
                            
                            auto test_wpid = dv->contained_by(test_p);
                            if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                                auto temp_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                                if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*units::cm, 0, 0)) {
                                    n_bad++;
                                }
                            }
                        }
                        if (n_bad > 0) flag_update = false;
                    }
                }
                
                // Check if sg is solid in all three views (if vtx1 has only 2 connections)
                if ((!flag_update) && boost::degree(vtx1_vd, graph) == 2) {
                    const auto& sg_wcpts2 = sg->wcpts();
                    for (size_t i = 0; i < sg_wcpts2.size(); i++) {
                        WireCell::Point test_p = sg_wcpts2.at(i).point;

                        auto test_wpid = dv->contained_by(test_p);
                        if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                            auto temp_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*units::cm, 0, 0)) {
                                flag_update = true;
                            }
                        }

                        // Check midpoint
                        if (i + 1 != sg_wcpts2.size()) {
                            WireCell::Point mid_p(
                                test_p.x() + (sg_wcpts2.at(i+1).point.x() - test_p.x()) / 2.,
                                test_p.y() + (sg_wcpts2.at(i+1).point.y() - test_p.y()) / 2.,
                                test_p.z() + (sg_wcpts2.at(i+1).point.z() - test_p.z()) / 2.
                            );
                            
                            auto mid_wpid = dv->contained_by(mid_p);
                            if (mid_wpid.face() != -1 && mid_wpid.apa() != -1) {
                                auto mid_p_raw = transform->backward(mid_p, cluster_t0, mid_wpid.face(), mid_wpid.apa());
                                if (!grouping->is_good_point(mid_p_raw, mid_wpid.apa(), mid_wpid.face(), 0.2*units::cm, 0, 0)) {
                                    flag_update = true;
                                }
                            }
                        }
                    }
                }
                
                // Perform the merge
                if (flag_update) {
                    s_log->debug("examine_structure_final_2: cluster {} merging vtx ({:.2f},{:.2f},{:.2f}) into main_vtx ({:.2f},{:.2f},{:.2f}) dis={:.2f} cm",
                        cluster.ident(),
                        vtx1_point.x()/units::cm, vtx1_point.y()/units::cm, vtx1_point.z()/units::cm,
                        main_vtx_point.x()/units::cm, main_vtx_point.y()/units::cm, main_vtx_point.z()/units::cm,
                        dis/units::cm);

                    // Use helper function to merge vtx1 into main_vertex
                    merge_vertex_into_another(graph, vtx1, main_vertex, dv);
                    
                    break;
                }
            }
        }
        
        // If updated, redo tracking and continue loop
        if (flag_update) {
            flag_continue = true;
            flag_updated = true;
            track_fitter.do_multi_tracking(true, true, false, false, false, &cluster);
        }
    }
    
    return flag_updated;
}

bool PatternAlgorithms::examine_structure_final_3(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv) {
    bool flag_updated = false;
    
    if (!main_vertex || !main_vertex->descriptor_valid()) return flag_updated;
    
    // Get transform and grouping for point validation
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster.get_scope_transform(cluster.get_default_scope()));
    double cluster_t0 = cluster.get_cluster_t0();
    auto grouping = cluster.grouping();
    
    if (!transform || !grouping) return flag_updated;
    
    // Continue looping until no more updates
    bool flag_continue = true;
    while (flag_continue) {
        flag_continue = false;
        bool flag_update = false;
        
        auto main_vd = main_vertex->get_descriptor();
        
        // Loop over all segments connected to main_vertex
        auto [ebegin, eend] = boost::out_edges(main_vd, graph);
        for (auto eit = ebegin; eit != eend; ++eit) {
            SegmentPtr sg = graph[*eit].segment;
            if (!sg) continue;
            
            // Find the other vertex of this segment
            VertexPtr vtx1 = find_other_vertex(graph, sg, main_vertex);
            if (!vtx1 || !vtx1->descriptor_valid()) continue;
            
            // Skip if vtx1 has only 1 connection
            auto vtx1_vd = vtx1->get_descriptor();
            if (boost::degree(vtx1_vd, graph) == 1) continue;
            
            // Check distance between vertices
            WireCell::Point main_vtx_point = main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;
            WireCell::Point vtx1_point = vtx1->fit().valid() ? vtx1->fit().point : vtx1->wcpt().point;
            
            double dis = ray_length(Ray{main_vtx_point, vtx1_point});
            
            if (dis < 2.5*units::cm) {
                // Check if main_vertex can be merged into vtx1
                flag_update = true;
                
                // Check all segments connected to main_vertex (except sg)
                auto [main_ebegin, main_eend] = boost::out_edges(main_vd, graph);
                for (auto main_eit = main_ebegin; main_eit != main_eend; ++main_eit) {
                    SegmentPtr sg1 = graph[*main_eit].segment;
                    if (!sg1 || sg1 == sg) continue;
                    
                    const auto& sg1_wcpts = sg1->wcpts();
                    if (sg1_wcpts.empty()) continue;

                    // Determine which end of sg connects to main_vertex
                    const auto& sg_wcpts = sg->wcpts();
                    if (sg_wcpts.empty()) continue;

                    WireCell::Point main_wcpt_point = main_vertex->wcpt().point;
                    bool flag_start = (ray_length(Ray{sg_wcpts.front().point, main_wcpt_point}) < 0.01*units::cm);

                    // Find point at ~3cm from main_vertex
                    WireCell::Point min_point = sg1_wcpts.front().point;
                    double min_dis = 1e9;
                    int min_index = 0;

                    for (size_t i = 0; i < sg1_wcpts.size(); i++) {
                        double dis = std::fabs(ray_length(Ray{sg1_wcpts.at(i).point, main_vtx_point}) - 3*units::cm);
                        if (dis < min_dis) {
                            min_dis = dis;
                            min_point = sg1_wcpts.at(i).point;
                            min_index = i;
                        }
                    }

                    // Check connectivity from min_point to main_vertex
                    bool flag_connect = true;
                    if (flag_start) {
                        for (int i = min_index; i >= 0; i--) {
                            auto test_wpid = dv->contained_by(sg1_wcpts.at(i).point);
                            if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                                auto temp_p_raw = transform->backward(sg1_wcpts.at(i).point, cluster_t0, test_wpid.face(), test_wpid.apa());
                                if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*units::cm, 0, 0)) {
                                    flag_connect = false;
                                    break;
                                }
                            }
                        }
                    } else {
                        for (size_t i = min_index; i < sg1_wcpts.size(); i++) {
                            auto test_wpid = dv->contained_by(sg1_wcpts.at(i).point);
                            if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                                auto temp_p_raw = transform->backward(sg1_wcpts.at(i).point, cluster_t0, test_wpid.face(), test_wpid.apa());
                                if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*units::cm, 0, 0)) {
                                    flag_connect = false;
                                    break;
                                }
                            }
                        }
                    }
                    
                    // Check path from min_point to vtx1
                    if (flag_connect) {
                        double step_size = 0.3 * units::cm;
                        WireCell::Point start_p = min_point;
                        WireCell::Point end_p = vtx1_point;
                        int ncount = std::round(ray_length(Ray{start_p, end_p}) / step_size);
                        int n_bad = 0;
                        
                        for (int i = 1; i < ncount; i++) {
                            WireCell::Point test_p(
                                start_p.x() + (end_p.x() - start_p.x()) / ncount * i,
                                start_p.y() + (end_p.y() - start_p.y()) / ncount * i,
                                start_p.z() + (end_p.z() - start_p.z()) / ncount * i
                            );
                            
                            auto test_wpid = dv->contained_by(test_p);
                            if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                                auto temp_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                                if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.3*units::cm, 0, 0)) {
                                    n_bad++;
                                }
                            }
                        }
                        if (n_bad > 0) flag_update = false;
                    }
                }
                
                // Perform the merge
                if (flag_update) {
                    s_log->debug("examine_structure_final_3: cluster {} merging main_vtx ({:.2f},{:.2f},{:.2f}) into vtx ({:.2f},{:.2f},{:.2f}) dis={:.2f} cm",
                        cluster.ident(),
                        main_vtx_point.x()/units::cm, main_vtx_point.y()/units::cm, main_vtx_point.z()/units::cm,
                        vtx1_point.x()/units::cm, vtx1_point.y()/units::cm, vtx1_point.z()/units::cm,
                        dis/units::cm);
                    
                    // Collect segments to update
                    std::vector<SegmentPtr> segments_to_update;
                    auto [main_ebegin2, main_eend2] = boost::out_edges(main_vd, graph);
                    for (auto main_eit = main_ebegin2; main_eit != main_eend2; ++main_eit) {
                        SegmentPtr sg1 = graph[*main_eit].segment;
                        if (sg1 && sg1 != sg) {
                            segments_to_update.push_back(sg1);
                        }
                    }
                    
                    // Process each segment connected to main_vertex (except sg)
                    for (auto sg1 : segments_to_update) {
                        WCPoint vtx_new_wcp = vtx1->wcpt();
                        std::vector<WCPoint> vec_wcps = sg1->wcpts();
                        
                        if (vec_wcps.empty()) continue;
                        
                        // Determine orientation: compare against steiner node position, not fit position
                        WireCell::Point main_wcpt_point2 = main_vertex->wcpt().point;
                        bool flag_front = (ray_length(Ray{vec_wcps.front().point, main_wcpt_point2}) < 0.01*units::cm);
                        
                        // Find point at ~3cm from main_vertex
                        WCPoint min_wcp = vec_wcps.front();
                        double min_dis = 1e9;
                        
                        for (size_t j = 0; j < vec_wcps.size(); j++) {
                            double dis1 = ray_length(Ray{vec_wcps.at(j).point, main_vtx_point});
                            double dis = std::fabs(dis1 - 3.0*units::cm);
                            if (dis < min_dis) {
                                min_wcp = vec_wcps.at(j);
                                min_dis = dis;
                            }
                        }
                        
                        // Build shortest path from vtx1 to min_wcp
                        std::list<WCPoint> new_list;
                        new_list.push_back(vtx_new_wcp);
                        
                        // Add intermediate points using steiner point cloud
                        {
                            double dis_step = 1.0*units::cm;
                            int ncount = std::round(ray_length(Ray{vtx_new_wcp.point, min_wcp.point}) / dis_step);
                            if (ncount < 2) ncount = 2;
                            
                            for (int qx = 1; qx < ncount; qx++) {
                                WireCell::Point tmp_p(
                                    vtx_new_wcp.point.x() + (min_wcp.point.x() - vtx_new_wcp.point.x()) / ncount * qx,
                                    vtx_new_wcp.point.y() + (min_wcp.point.y() - vtx_new_wcp.point.y()) / ncount * qx,
                                    vtx_new_wcp.point.z() + (min_wcp.point.z() - vtx_new_wcp.point.z()) / ncount * qx
                                );
                                
                                auto [tmp_idx, tmp_wcp_pt] = cluster.get_closest_wcpoint(tmp_p);
                                WCPoint tmp_wcp;
                                tmp_wcp.point = tmp_wcp_pt;
                                
                                // Check distance
                                if (ray_length(Ray{tmp_wcp.point, tmp_p}) > 0.3*units::cm) continue;
                                
                                // Check if different from last point and min_wcp
                                if (ray_length(Ray{tmp_wcp.point, new_list.back().point}) > 0.01*units::cm && 
                                    ray_length(Ray{tmp_wcp.point, min_wcp.point}) > 0.01*units::cm) {
                                    new_list.push_back(tmp_wcp);
                                }
                            }
                        }
                        new_list.push_back(min_wcp);
                        
                        // Merge with existing wcpts
                        std::list<WCPoint> old_list;
                        std::copy(vec_wcps.begin(), vec_wcps.end(), std::back_inserter(old_list));
                        
                        if (flag_front) {
                            // Remove points up to min_wcp from front
                            while (old_list.size() > 0 && ray_length(Ray{old_list.front().point, min_wcp.point}) > 0.01*units::cm) {
                                old_list.pop_front();
                            }
                            if (old_list.size() > 0) old_list.pop_front();
                            
                            // Add new_list to front in reverse
                            for (auto it = new_list.rbegin(); it != new_list.rend(); it++) {
                                old_list.push_front(*it);
                            }
                        } else {
                            // Remove points up to min_wcp from back
                            while (old_list.size() > 0 && ray_length(Ray{old_list.back().point, min_wcp.point}) > 0.01*units::cm) {
                                old_list.pop_back();
                            }
                            if (old_list.size() > 0) old_list.pop_back();
                            
                            // Add new_list to back in reverse
                            for (auto it = new_list.rbegin(); it != new_list.rend(); it++) {
                                old_list.push_back(*it);
                            }
                        }
                        
                        // Update segment wcpts
                        std::vector<WCPoint> new_wcpts;
                        new_wcpts.reserve(old_list.size());
                        std::copy(std::begin(old_list), std::end(old_list), std::back_inserter(new_wcpts));
                        sg1->wcpts(new_wcpts);
                        std::vector<Facade::geo_point_t> main_pts_sg1f;
                        for (const auto& wcp : new_wcpts) main_pts_sg1f.push_back(wcp.point);
                        create_segment_point_cloud(sg1, main_pts_sg1f, dv, "main");
                    }

                    // Update main_vertex to vtx1's position
                    main_vertex->wcpt(vtx1->wcpt());
                    if (vtx1->fit().valid()) {
                        main_vertex->fit(vtx1->fit());
                    }
                    
                    // Reconnect segments from vtx1 to main_vertex (except sg)
                    std::vector<SegmentPtr> vtx1_segments;
                    auto [vtx1_ebegin2, vtx1_eend2] = boost::out_edges(vtx1_vd, graph);
                    for (auto vtx1_eit = vtx1_ebegin2; vtx1_eit != vtx1_eend2; ++vtx1_eit) {
                        SegmentPtr sg1 = graph[*vtx1_eit].segment;
                        if (sg1 && sg1 != sg) {
                            vtx1_segments.push_back(sg1);
                        }
                    }
                    
                    for (auto sg1 : vtx1_segments) {
                        VertexPtr tt_vtx = find_other_vertex(graph, sg1, vtx1);
                        if (tt_vtx && tt_vtx != main_vertex) {
                            remove_segment(graph, sg1);
                            add_segment(graph, sg1, main_vertex, tt_vtx);
                        } else {
                            // Self-loop case
                            remove_segment(graph, sg1);
                        }
                    }
                    
                    // Delete sg first, then vtx1.
                    // IMPORTANT: boost::remove_vertex automatically removes all
                    // incident edges, so removing vtx1 first would silently free
                    // the sg edge descriptor, leaving sg->descriptor_valid()==true
                    // while the underlying edge is gone.  The subsequent
                    // remove_segment(sg) would then call boost::remove_edge on a
                    // freed descriptor → double-free crash.
                    remove_segment(graph, sg);
                    remove_vertex(graph, vtx1);
                    
                    break;
                }
            }
        }
        
        // If updated, redo tracking and continue loop
        if (flag_update) {
            flag_continue = true;
            flag_updated = true;
            track_fitter.do_multi_tracking(true, true, false, false, false, &cluster);
        }
    }
    
    return flag_updated;
}
   

bool PatternAlgorithms::examine_structure_final(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv) {
    examine_structure_final_1(graph, main_vertex, cluster, track_fitter, dv);
    examine_structure_final_1p(graph, main_vertex, cluster, track_fitter, dv);
    examine_structure_final_2(graph, main_vertex, cluster, track_fitter, dv);
    examine_structure_final_3(graph, main_vertex, cluster, track_fitter, dv);
    return true;
}

