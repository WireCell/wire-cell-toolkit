#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

std::set<VertexPtr> PatternAlgorithms::find_cluster_vertices(Graph& graph, const Facade::Cluster& cluster)
{
    std::set<VertexPtr> result;
    
    // Iterate through all vertices in the graph
    auto [vbegin, vend] = boost::vertices(graph);
    for (auto vit = vbegin; vit != vend; ++vit) {
        VertexPtr vtx = graph[*vit].vertex;
        
        // Check if this vertex belongs to the specified cluster
        if (vtx && vtx->cluster() && vtx->cluster() == &cluster) {
            result.insert(vtx);
        }
    }
    
    return result;
}

std::set<SegmentPtr> PatternAlgorithms::find_cluster_segments(Graph& graph, const Facade::Cluster& cluster)
{
    std::set<SegmentPtr> result;
    
    // Iterate through all edges (segments) in the graph
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        SegmentPtr seg = graph[*eit].segment;
        
        // Check if this segment belongs to the specified cluster
        if (seg && seg->cluster() && seg->cluster() == &cluster) {
            result.insert(seg);
        }
    }
    
    return result;
}

bool PatternAlgorithms::clean_up_graph(Graph& graph, const Facade::Cluster& cluster)
{
    bool modified = false;
    
    // First, find and remove all segments associated with this cluster
    std::set<SegmentPtr> segments_to_remove = find_cluster_segments(graph, cluster);
    for (auto seg : segments_to_remove) {
        if (remove_segment(graph, seg)) {
            modified = true;
        }
    }
    
    // Then, find and remove all vertices associated with this cluster
    // Note: vertices that are still connected to other segments won't be removed
    // until their segments are removed first
    std::set<VertexPtr> vertices_to_remove = find_cluster_vertices(graph, cluster);
    for (auto vtx : vertices_to_remove) {
        if (remove_vertex(graph, vtx)) {
            modified = true;
        }
    }
    
    return modified;
}

std::vector<Facade::geo_point_t> PatternAlgorithms::do_rough_path(const Facade::Cluster& cluster,Facade::geo_point_t& first_point, Facade::geo_point_t& last_point){  
        // Find closest indices in the steiner point cloud
        auto first_knn_results = cluster.kd_steiner_knn(1, first_point, "steiner_pc");
        auto last_knn_results = cluster.kd_steiner_knn(1, last_point, "steiner_pc");
        
        auto first_index = first_knn_results[0].first;  // Get the index from the first result
        auto last_index = last_knn_results[0].first;   // Get the index from the first result
 
        // 4. Use Steiner graph to find the shortest path
        const std::vector<size_t>& path_indices = 
            cluster.graph_algorithms("steiner_graph").shortest_path(first_index, last_index);
            
        std::vector<Facade::geo_point_t> path_points;
        const auto& steiner_pc = cluster.get_pc("steiner_pc");
        const auto& coords = cluster.get_default_scope().coords;
        const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
        const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
        const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();

        for (size_t idx : path_indices) {
            path_points.emplace_back(x_coords[idx], y_coords[idx], z_coords[idx]);
        }
        return path_points;
}

SegmentPtr PatternAlgorithms::create_segment_for_cluster(WireCell::Clus::Facade::Cluster& cluster, IDetectorVolumes::pointer dv, const std::vector<Facade::geo_point_t>& path_points, int dir){
     // Step 3: Prepare segment data
    std::vector<PR::WCPoint> wcpoints;
    // const auto transform = m_pcts->pc_transform(cluster.get_scope_transform(cluster.get_default_scope()));
    // Step 4: Create segment connecting the vertices
    auto segment = PR::make_segment();
    
    // create and associate Dynamic Point Cloud
    for (const auto& point : path_points) {
        PR::WCPoint wcp;
        wcp.point = point; 
        wcpoints.push_back(wcp);
    }

    // Step 5: Configure the segment
    segment->wcpts(wcpoints).cluster(&cluster).dirsign(dir); // direction: +1, 0, or -1
            
    // auto& wcpts = segment->wcpts();
    // for (size_t i=0;i!=path_points.size(); i++){
    //     std::cout << "A: " << i << " " << path_points.at(i) << " " << wcpts.at(i).point << std::endl;
    // }
    create_segment_point_cloud(segment, path_points, dv, "main");

    return segment;
}

SegmentPtr PatternAlgorithms::create_segment_from_vertices(Graph& graph, Facade::Cluster& cluster, VertexPtr v1, VertexPtr v2, IDetectorVolumes::pointer dv){
     // Create Segment using the vertices to derive a path 
    auto path_points = do_rough_path(cluster, v1->wcpt().point, v2->wcpt().point);
    
    // Check if path has enough points (similar to WCPPID check)
    if (path_points.size() <= 1) {
        return nullptr;
    }
    
    auto seg = create_segment_for_cluster(cluster, dv, path_points);
    WireCell::Clus::PR::add_segment(graph, seg, v1, v2);
    return seg;
}



SegmentPtr PatternAlgorithms::init_first_segment(Graph& graph, Facade::Cluster& cluster, Facade::Cluster* main_cluster,TrackFitting& track_fitter, IDetectorVolumes::pointer dv, bool flag_back_search)
{
    // Get two boundary points from the cluster
    auto boundary_indices = cluster.get_two_boundary_steiner_graph_idx("steiner_graph", "steiner_pc");

    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();

    // Add the two boundary points as additional extreme point groups
    Facade::geo_point_t boundary_point_first(x_coords[boundary_indices.first], 
                                y_coords[boundary_indices.first], 
                                z_coords[boundary_indices.first]);
    Facade::geo_point_t boundary_point_second(x_coords[boundary_indices.second], 
                                y_coords[boundary_indices.second], 
                                z_coords[boundary_indices.second]);
    Facade::geo_point_t first_pt = boundary_point_first;
    Facade::geo_point_t second_pt = boundary_point_second;
    
    // Determine the starting point based on whether this is the main cluster or not
    if (cluster.get_flag(Facade::Flags::main_cluster)) {
        // Main cluster: start from downstream (or upstream if flag_back_search)
        if (flag_back_search) {
            // Start from high z (upstream/backward)
            if (first_pt.z() < second_pt.z()) {
                std::swap(first_pt, second_pt);
            }
        } else {
            // Start from low z (downstream/forward)
            if (first_pt.z() > second_pt.z()) {
                std::swap(first_pt, second_pt);
            }
        }
    } else if (main_cluster) {
        // Non-main cluster: start from the point closest to main cluster
        // Find closest distances to main cluster's Steiner point cloud
        auto knn1 = main_cluster->kd_steiner_knn(1, first_pt, "steiner_pc");
        auto knn2 = main_cluster->kd_steiner_knn(1, second_pt, "steiner_pc");
        
        if (!knn1.empty() && !knn2.empty()) {
            double dis1 = std::sqrt(knn1[0].second);
            double dis2 = std::sqrt(knn2[0].second);
            
            // Start from the point closer to main cluster
            if (dis2 < dis1) {
                std::swap(first_pt, second_pt);
            }
        }
    }
    
    // Create vertices for the endpoints
    VertexPtr v1 = make_vertex(graph);
    v1->wcpt().point = first_pt;
    v1->cluster(&cluster);
    VertexPtr v2 = make_vertex(graph);
    v2->wcpt().point = second_pt;
    v2->cluster(&cluster);

    auto seg = create_segment_from_vertices(graph, cluster, v1, v2, dv);
    if (!seg) {
        remove_vertex(graph, v1);
        remove_vertex(graph, v2);
        return nullptr;
    }

    // // Create Segment using the vertices to derive a path 
    // auto path_points = do_rough_path(cluster, first_pt, second_pt);
    // // Check if path has enough points (similar to WCPPID check)
    // if (path_points.size() <= 1) {     
    // }
    // auto seg = create_segment_for_cluster(cluster, dv, path_points);
    // WireCell::Clus::PR::add_segment(graph, seg, v1, v2);

    // perform fitting ...
    track_fitter.add_segment(seg);
    track_fitter.do_single_tracking(seg, true, true);
    const auto& fine_path = track_fitter.get_fine_tracking_path();
    const auto& dQ_vec = track_fitter.get_dQ();
    const auto& dx_vec = track_fitter.get_dx();
    const auto& pu_vec = track_fitter.get_pu();
    const auto& pv_vec = track_fitter.get_pv();
    const auto& pw_vec = track_fitter.get_pw();
    const auto& pt_vec = track_fitter.get_pt();
    const auto& chi2_vec = track_fitter.get_reduced_chi2();
    
    if (fine_path.size()>1) {
        v1->fit().point = fine_path.front().first;
        if (!dQ_vec.empty()) v1->fit().dQ = dQ_vec.front();
        if (!dx_vec.empty()) v1->fit().dx = dx_vec.front();
        if (!pu_vec.empty()) v1->fit().pu = pu_vec.front();
        if (!pv_vec.empty()) v1->fit().pv = pv_vec.front();
        if (!pw_vec.empty()) v1->fit().pw = pw_vec.front();
        if (!pt_vec.empty()) v1->fit().pt = pt_vec.front();
        if (!chi2_vec.empty()) v1->fit().reduced_chi2 = chi2_vec.front();
        
        v2->fit().point = fine_path.back().first;
        if (!dQ_vec.empty()) v2->fit().dQ = dQ_vec.back();
        if (!dx_vec.empty()) v2->fit().dx = dx_vec.back();
        if (!pu_vec.empty()) v2->fit().pu = pu_vec.back();
        if (!pv_vec.empty()) v2->fit().pv = pv_vec.back();
        if (!pw_vec.empty()) v2->fit().pw = pw_vec.back();
        if (!pt_vec.empty()) v2->fit().pt = pt_vec.back();
        if (!chi2_vec.empty()) v2->fit().reduced_chi2 = chi2_vec.back();

        // Set fit information for segment
        std::vector<Fit> fits;
        for (size_t i = 0; i < fine_path.size(); ++i) {
            Fit fit;
            fit.point = fine_path[i].first;
            if (i < dQ_vec.size()) fit.dQ = dQ_vec[i];
            if (i < dx_vec.size()) fit.dx = dx_vec[i];
            if (i < pu_vec.size()) fit.pu = pu_vec[i];
            if (i < pv_vec.size()) fit.pv = pv_vec[i];
            if (i < pw_vec.size()) fit.pw = pw_vec[i];
            if (i < pt_vec.size()) fit.pt = pt_vec[i];
            if (i < chi2_vec.size()) fit.reduced_chi2 = chi2_vec[i];
            fits.push_back(fit);
        }
        seg->fits(fits);
    }else{
        // Tracking failed, clean up
        remove_segment(graph, seg);
        remove_vertex(graph, v1);
        remove_vertex(graph, v2);
        return nullptr;
    }
    
    return seg;
}


std::pair<Facade::geo_point_t,  size_t> PatternAlgorithms::proto_extend_point(const Facade::Cluster& cluster, Facade::geo_point_t& p, Facade::geo_vector_t& dir, Facade::geo_vector_t& dir_other, bool flag_continue){
    const double step_dis = 1.0 * units::cm;
    
    // Get steiner point cloud data
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& steiner_x = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& steiner_y = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& steiner_z = steiner_pc.get(coords.at(2))->elements<double>();
    
    // Find closest point in steiner point cloud
    auto curr_knn_results = cluster.kd_steiner_knn(1, p, "steiner_pc");
    if (curr_knn_results.empty()) {
        return std::make_pair(p, -1); // No steiner point found ...
    }
    
    size_t curr_index = curr_knn_results[0].first;
    Facade::geo_point_t curr_wcp(steiner_x[curr_index], steiner_y[curr_index], steiner_z[curr_index]);
    Facade::geo_point_t next_wcp = curr_wcp;
    
    // Save starting position and direction
    Facade::geo_point_t saved_start_wcp = curr_wcp;
    Facade::geo_vector_t saved_dir = dir;
    
    // Forward search
    while(flag_continue){
        flag_continue = false;
        
        for (int i = 0; i != 3; i++){
            Facade::geo_point_t test_p(
                curr_wcp.x() + dir.x() * step_dis * (i + 1),
                curr_wcp.y() + dir.y() * step_dis * (i + 1),
                curr_wcp.z() + dir.z() * step_dis * (i + 1)
            );
            
            // Try steiner point cloud first
            auto next_knn_steiner = cluster.kd_steiner_knn(1, test_p, "steiner_pc");
            if (!next_knn_steiner.empty()) {
                size_t next_index = next_knn_steiner[0].first;
                next_wcp = Facade::geo_point_t(steiner_x[next_index], steiner_y[next_index], steiner_z[next_index]);
                Facade::geo_vector_t dir2(
                    next_wcp.x() - curr_wcp.x(), 
                    next_wcp.y() - curr_wcp.y(), 
                    next_wcp.z() - curr_wcp.z()
                );
                
                double mag2 = dir2.magnitude();
                if (mag2 != 0) {
                    double angle = std::acos(dir2.dot(dir) / mag2) / 3.1415926 * 180.0;
                    if (angle < 25.0) {
                        flag_continue = true;
                        curr_wcp = next_wcp;
                        curr_index = next_index;
                        dir = dir2 + dir * 5.0 * units::cm; // momentum trick
                        dir = dir / dir.magnitude();
                        break;
                    }
                }
            }
            
            // Try regular point cloud
            auto closest_result = cluster.get_closest_wcpoint(test_p);
            // size_t regular_index = closest_result.first;
            next_wcp = closest_result.second;
            
            Facade::geo_vector_t dir1(
                next_wcp.x() - curr_wcp.x(), 
                next_wcp.y() - curr_wcp.y(), 
                next_wcp.z() - curr_wcp.z()
            );
            
            double mag1 = dir1.magnitude();
            if (mag1 != 0) {
                double angle = std::acos(dir1.dot(dir) / mag1) / 3.1415926 * 180.0;
                if (angle < 17.5) {
                    flag_continue = true;
                    curr_wcp = next_wcp;
                    // For regular point cloud, we need to find it in steiner cloud again
                    auto updated_knn = cluster.kd_steiner_knn(1, curr_wcp, "steiner_pc");
                    if (!updated_knn.empty()) {
                        curr_index = updated_knn[0].first;
                    }
                    dir = dir1 + dir * 5.0 * units::cm; // momentum trick
                    dir = dir / dir.magnitude();
                    break;
                }
            }
        }
    }
    
    // Ensure we return the steiner point cloud position
    Facade::geo_point_t test_p(curr_wcp.x(), curr_wcp.y(), curr_wcp.z());
    auto final_knn = cluster.kd_steiner_knn(1, test_p, "steiner_pc");
    if (!final_knn.empty()) {
        curr_index = final_knn[0].first;
        curr_wcp = Facade::geo_point_t(steiner_x[curr_index], steiner_y[curr_index], steiner_z[curr_index]);
    }
    
    // Return: point, (cloud_type=2 for steiner, point_index)
    return std::make_pair(curr_wcp,  curr_index);
}

bool PatternAlgorithms::proto_break_tracks(const Facade::Cluster& cluster, const Facade::geo_point_t& first_wcp, const Facade::geo_point_t& curr_wcp, const Facade::geo_point_t& last_wcp, std::list<Facade::geo_point_t>& wcps_list1, std::list<Facade::geo_point_t>& wcps_list2, bool flag_pass_check){
    
    // Calculate distances
    double dis1 = std::sqrt(std::pow(curr_wcp.x() - first_wcp.x(), 2) + 
                            std::pow(curr_wcp.y() - first_wcp.y(), 2) + 
                            std::pow(curr_wcp.z() - first_wcp.z(), 2));
    double dis2 = std::sqrt(std::pow(curr_wcp.x() - last_wcp.x(), 2) + 
                            std::pow(curr_wcp.y() - last_wcp.y(), 2) + 
                            std::pow(curr_wcp.z() - last_wcp.z(), 2));
    
    // Check if distances are sufficient or if we should pass the check
    if ((dis1 > 1.0 * units::cm && dis2 > 1.0 * units::cm) || flag_pass_check) {
        // Find shortest path from first_wcp to curr_wcp using steiner graph
        Facade::geo_point_t first_point = first_wcp;
        Facade::geo_point_t curr_point = curr_wcp;
        auto path1 = do_rough_path(cluster, first_point, curr_point);
        // Convert vector to list
        wcps_list1.clear();
        for (const auto& pt : path1) {
            wcps_list1.push_back(pt);
        }
        
        // Find shortest path from curr_wcp to last_wcp using steiner graph
        Facade::geo_point_t curr_point2 = curr_wcp;
        Facade::geo_point_t last_point = last_wcp;
        auto path2 = do_rough_path(cluster, curr_point2, last_point);
        // Convert vector to list
        wcps_list2.clear();
        for (const auto& pt : path2) {
            wcps_list2.push_back(pt);
        }
        
        // Remove overlapping points at the junction
        // Count how many points overlap from the end of list1 and beginning of list2
        int count = 0;
        if (!wcps_list1.empty() && !wcps_list2.empty()) {
            // Compare points from the end of list1 with the beginning of list2
            // Use reverse iterator for list1 and forward iterator for list2
            auto it1 = wcps_list1.rbegin();  // Start from the back of list1
            auto it2 = wcps_list2.begin();   // Start from the front of list2
            
            while (it1 != wcps_list1.rend() && it2 != wcps_list2.end()) {
                // Check if points are the same (within tolerance)
                double dx = it1->x() - it2->x();
                double dy = it1->y() - it2->y();
                double dz = it1->z() - it2->z();
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (dist < 0.01 * units::cm) {  // same point
                    count++;
                    ++it1;
                    ++it2;
                } else {
                    break;  // no more overlapping points
                }
            }
            
            // Remove overlapping points (keep one copy at the junction)
            for (int i = 0; i < count; i++) {
                if (i + 1 != count) {  // Keep the last overlapping point
                    if (!wcps_list1.empty()) wcps_list1.pop_back();
                    if (!wcps_list2.empty()) wcps_list2.pop_front();
                }
            }
        }
        
        // Check if we have valid paths
        if (wcps_list1.size() <= 1 || wcps_list2.size() <= 1) {
            return false;
        }
        
        return true;
    } else {
        return false;
    }
}

bool PatternAlgorithms::replace_segment_and_vertex(Graph& graph, SegmentPtr& seg, VertexPtr old_vertex, VertexPtr new_vertex, IDetectorVolumes::pointer dv){
    // Get the cluster from the old segment
    auto cluster = seg->cluster();
    if (!cluster) {
        return false;
    }
    
    // Get the other vertex connected to this segment (the one we'll keep)
    VertexPtr other_vertex = find_other_vertex(graph, seg, old_vertex);
    if (!other_vertex) {
        return false;
    }
    
    // Create new segment with the path points
    SegmentPtr new_seg = create_segment_from_vertices(graph, *cluster, other_vertex, new_vertex, dv);
    if (!new_seg) {
        return false;
    }
    
    // Remove the old segment (this will disconnect it from the graph)
    remove_segment(graph, seg);
    
    // Remove the old vertex if it no longer has any connected segments
    if (old_vertex->descriptor_valid()) {
        auto vd = old_vertex->get_descriptor();
        if (boost::degree(vd, graph) == 0) {
            remove_vertex(graph, old_vertex);
        }
    }
        
    // Update the output parameter
    seg = new_seg;
    
    return true;
}



bool PatternAlgorithms::replace_segment_and_vertex(Graph& graph, SegmentPtr& seg, VertexPtr& vtx, std::list<Facade::geo_point_t>& path_point_list, Facade::geo_point_t& break_point, IDetectorVolumes::pointer dv) {
    // // Check that the vertex is only connected to one segment
    // if (!vtx->descriptor_valid()) {
    //     return false;
    // }
    // auto vd = vtx->get_descriptor();
    // if (boost::degree(vd, graph) != 1) {
    //     return false;  // Vertex is connected to more than one segment, cannot replace
    // }
    
    // Get the cluster from the old segment
    auto cluster = seg->cluster();
    if (!cluster) {
        return false;
    }
    
    // Get the other vertex connected to this segment (the one we'll keep)
    VertexPtr other_vertex = find_other_vertex(graph, seg, vtx);
    if (!other_vertex) {
        return false;
    }
    
    // Create new vertex at the break point
    VertexPtr new_vtx = make_vertex(graph);
    new_vtx->wcpt().point = break_point;
    new_vtx->cluster(cluster);
    
    // Convert list to vector for create_segment_for_cluster
    std::vector<Facade::geo_point_t> path_points;
    for (const auto& pt : path_point_list) {
        path_points.push_back(pt);
    }
    
    // Check if path has enough points
    if (path_points.size() <= 1) {
        remove_vertex(graph, new_vtx);
        return false;
    }
    
    // Create new segment with the path points
    SegmentPtr new_seg = create_segment_for_cluster(*cluster, dv, path_points, seg->dirsign());
    if (!new_seg) {
        remove_vertex(graph, new_vtx);
        return false;
    }
    
    // Remove the old segment (this will disconnect it from the graph)
    remove_segment(graph, seg);

    // Remove the old vertex, if it no longer has any connected segments  
    auto vd = vtx->get_descriptor();
    if (boost::degree(vd, graph) == 0) remove_vertex(graph, vtx);
    
    // Add the new segment connecting other_vertex and new_vtx
    add_segment(graph, new_seg, other_vertex, new_vtx);
    
    // Update the output parameters
    seg = new_seg;
    vtx = new_vtx;
    
    return true;
}

 bool PatternAlgorithms::break_segment_into_two(Graph& graph, VertexPtr vtx1, SegmentPtr seg, VertexPtr vtx2, std::list<Facade::geo_point_t>& path_point_list1, Facade::geo_point_t& break_point, std::list<Facade::geo_point_t>& path_point_list2, IDetectorVolumes::pointer dv){
    // Get the cluster from the old segment
    auto cluster = seg->cluster();
    if (!cluster) {
        return false;
    }
    
    // Verify that vtx1 and vtx2 are the endpoints of seg
    auto [v1, v2] = find_vertices(graph, seg);
    if ((v1 != vtx1 || v2 != vtx2) && (v1 != vtx2 || v2 != vtx1)) {
        return false;  // The provided vertices don't match the segment endpoints
    }
    
    // Create new vertex at the break point
    VertexPtr new_vtx = make_vertex(graph);
    new_vtx->wcpt().point = break_point;
    new_vtx->cluster(cluster);
    
    // Convert lists to vectors for create_segment_for_cluster
    std::vector<Facade::geo_point_t> path_points1;
    for (const auto& pt : path_point_list1) {
        path_points1.push_back(pt);
    }
    
    std::vector<Facade::geo_point_t> path_points2;
    for (const auto& pt : path_point_list2) {
        path_points2.push_back(pt);
    }
    
    // Check if paths have enough points
    if (path_points1.size() <= 1 || path_points2.size() <= 1) {
        remove_vertex(graph, new_vtx);
        return false;
    }
    
    // Create first new segment with path_points1
    SegmentPtr new_seg1 = create_segment_for_cluster(*cluster, dv, path_points1, seg->dirsign());
    if (!new_seg1) {
        remove_vertex(graph, new_vtx);
        return false;
    }
    
    // Create second new segment with path_points2
    SegmentPtr new_seg2 = create_segment_for_cluster(*cluster, dv, path_points2, seg->dirsign());
    if (!new_seg2) {
        remove_vertex(graph, new_vtx);
        return false;
    }
    
    // Remove the old segment
    remove_segment(graph, seg);
    
    // Add the first new segment connecting vtx1 and new_vtx
    add_segment(graph, new_seg1, vtx1, new_vtx);
    
    // Add the second new segment connecting new_vtx and vtx2
    add_segment(graph, new_seg2, new_vtx, vtx2);
    
    return true;
 }

 bool PatternAlgorithms::break_segments(Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, std::vector<SegmentPtr>& remaining_segments, float dis_cut) {
    bool flag_modified = false;
    int count = 0;
    std::set<size_t> saved_break_wcp_indices;
    
    while(!remaining_segments.empty() && count < 2) {
        SegmentPtr curr_sg = remaining_segments.back();
        auto cluster = curr_sg->cluster();
        remaining_segments.pop_back();
        
        // Get the two vertices of this segment
        auto [start_v, end_v] = find_vertices(graph, curr_sg);
        if (!start_v || !end_v) {
            continue;
        }
        
        // Check if vertices match the segment endpoints
        const auto& wcpts = curr_sg->wcpts();
        if (wcpts.size() < 2) continue;
        
        auto front_pt = wcpts.front().point;
        auto back_pt = wcpts.back().point;
        
        // Determine which vertex is start and which is end based on point positions
        double dis_sv_front = ray_length(Ray{start_v->wcpt().point, front_pt});
        double dis_sv_back = ray_length(Ray{start_v->wcpt().point, back_pt});
        
        if (dis_sv_front > dis_sv_back) {
            std::swap(start_v, end_v);
        }
        
        // Initialize the start test point
        Facade::geo_point_t break_wcp = start_v->wcpt().point;
        const auto& point_vec = curr_sg->wcpts();
        Facade::geo_point_t test_start_p = point_vec.front().point;
        
        if (dis_cut > 0) {
            for (size_t i = 0; i < point_vec.size(); ++i) {
                double dis = ray_length(Ray{point_vec[i].point, point_vec.front().point});
                if (dis > dis_cut) {
                    test_start_p = point_vec[i].point;
                    break;
                }
            }
        }
        
        // Search for kinks and extend the break point
        while(ray_length(Ray{start_v->wcpt().point, break_wcp}) <= 1.0 * units::cm &&
              ray_length(Ray{end_v->wcpt().point, break_wcp}) > 1.0 * units::cm) {
            
            auto kink_tuple = segment_search_kink(curr_sg, test_start_p, "fit");
            auto& [kink_point, dir1, dir2, flag_continue] = kink_tuple;
            
            if (dir1.magnitude() != 0) {
                // Find the extreme point
                Facade::geo_vector_t dir1_geo(dir1.x(), dir1.y(), dir1.z());
                Facade::geo_vector_t dir2_geo(dir2.x(), dir2.y(), dir2.z());
                Facade::geo_point_t kink_geo(kink_point.x(), kink_point.y(), kink_point.z());
                
                auto [break_pt, break_idx] = proto_extend_point(*cluster, kink_geo, dir1_geo, dir2_geo, flag_continue);
                break_wcp = break_pt;
                
                // Check if we've seen this break point before
                if (saved_break_wcp_indices.find(break_idx) != saved_break_wcp_indices.end()) {
                    test_start_p = kink_geo;
                    kink_tuple = segment_search_kink(curr_sg, test_start_p, "fit");
                    auto& [kink_point2, dir1_2, dir2_2, flag_continue2] = kink_tuple;
                    Facade::geo_vector_t dir1_geo2(dir1_2.x(), dir1_2.y(), dir1_2.z());
                    Facade::geo_vector_t dir2_geo2(dir2_2.x(), dir2_2.y(), dir2_2.z());
                    Facade::geo_point_t kink_geo2(kink_point2.x(), kink_point2.y(), kink_point2.z());
                    auto [break_pt2, break_idx2] = proto_extend_point(*cluster, kink_geo2, dir1_geo2, dir2_geo2, flag_continue2);
                    break_wcp = break_pt2;
                    break_idx = break_idx2;
                } else {
                    saved_break_wcp_indices.insert(break_idx);
                }
                
                if (ray_length(Ray{start_v->wcpt().point, break_wcp}) <= 1.0 * units::cm &&
                    ray_length(Ray{end_v->wcpt().point, break_wcp}) > 1.0 * units::cm) {
                    test_start_p = kink_geo;
                }
            } else {
                break;
            }
        }
        
        // Check if we should break the segment
        if (ray_length(Ray{start_v->wcpt().point, break_wcp}) > 1.0 * units::cm) {
            std::list<Facade::geo_point_t> wcps_list1;
            std::list<Facade::geo_point_t> wcps_list2;
            
            bool flag_break;
            bool flag_pass_check = false;
            
            // Check if end vertex is close to break point and has only one connection
            if (ray_length(Ray{end_v->wcpt().point, break_wcp}) < 1.0 * units::cm) {
                auto vd = end_v->get_descriptor();
                if (boost::degree(vd, graph) == 1) {
                    flag_pass_check = true;
                }
            }
            
            flag_break = proto_break_tracks(*cluster, start_v->wcpt().point, break_wcp, 
                                           end_v->wcpt().point, wcps_list1, wcps_list2, flag_pass_check);
            
            if (flag_break) {
                // Check geometry constraints
                Facade::geo_vector_t tv1 = end_v->wcpt().point - start_v->wcpt().point;
                Facade::geo_vector_t tv2 = end_v->wcpt().point - break_wcp;
                
                double min_dis = 1e9;
                for (const auto& wcp : wcps_list1) {
                    double dis = ray_length(Ray{wcp, end_v->wcpt().point});
                    if (dis < min_dis) min_dis = dis;
                }
                
                double angle = std::acos(tv1.dot(tv2) / (tv1.magnitude() * tv2.magnitude())) / 3.1415926 * 180.0;
                
                // Check if we should replace end vertex instead of breaking
                if (min_dis / units::cm < 1.5 && angle > 120) {
                    auto vd = end_v->get_descriptor();
                    if (boost::degree(vd, graph) == 1) {
                        // Replace segment and end vertex
                        SegmentPtr new_seg = curr_sg;
                        VertexPtr new_vtx = end_v;
                        if (replace_segment_and_vertex(graph, new_seg, new_vtx, wcps_list1, break_wcp, dv)) {
                            flag_modified = true;
                            // Perform tracking
                            // track_fitter.add_graph(&graph); added already
                            track_fitter.do_multi_tracking(true, true, false);
                        }
                    }
                } else {
                    // Break segment into two
                    if (break_segment_into_two(graph, start_v, curr_sg, end_v, wcps_list1, break_wcp, wcps_list2, dv)) {
                        flag_modified = true;
                        // Perform tracking
                        // track_fitter.add_graph(&graph); added already
                        track_fitter.do_multi_tracking(true, true, false);
                        
                        // Find the new segment connecting to end_v and add it back to remaining_segments
                        auto [new_vtx, new_end_v] = find_vertices(graph, curr_sg);
                        // The new segment created from wcps_list2 connects new_vtx to end_v
                        // We need to find it
                        auto erange = boost::out_edges(end_v->get_descriptor(), graph);
                        for (auto eit = erange.first; eit != erange.second; ++eit) {
                            SegmentPtr seg = graph[*eit].segment;
                            if (seg && seg->cluster() == curr_sg->cluster()) {
                                remaining_segments.push_back(seg);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    
    return flag_modified;
}


bool PatternAlgorithms::merge_two_segments_into_one(Graph& graph, SegmentPtr& seg1, VertexPtr& vtx, SegmentPtr& seg2, IDetectorVolumes::pointer dv){
    // Get cluster from seg1 (should be same as seg2)
    auto cluster = seg1->cluster();
    if (!cluster || cluster != seg2->cluster()) {
        return false;
    }
    
    // Get the other vertices (not vtx) from seg1 and seg2
    VertexPtr vtx1 = find_other_vertex(graph, seg1, vtx);
    VertexPtr vtx2 = find_other_vertex(graph, seg2, vtx);
    
    if (!vtx1 || !vtx2) {
        return false;
    }
    
    // Create new segment from vtx1 to vtx2
    SegmentPtr new_seg = create_segment_from_vertices(graph, *cluster, vtx1, vtx2, dv);
    if (!new_seg) {
        return false;
    }
    
    // Delete old segments
    remove_segment(graph, seg1);
    remove_segment(graph, seg2);
    
    // Delete the middle vertex
    remove_vertex(graph, vtx);
    
    // Update output parameter
    seg1 = new_seg;
    
    return true;
}
