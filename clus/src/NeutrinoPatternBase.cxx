#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

std::set<VertexPtr> PatternAlgorithms::find_vertices(Graph& graph, const Facade::Cluster& cluster)
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

std::set<SegmentPtr> PatternAlgorithms::find_segments(Graph& graph, const Facade::Cluster& cluster)
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
    std::set<SegmentPtr> segments_to_remove = find_segments(graph, cluster);
    for (auto seg : segments_to_remove) {
        if (remove_segment(graph, seg)) {
            modified = true;
        }
    }
    
    // Then, find and remove all vertices associated with this cluster
    // Note: vertices that are still connected to other segments won't be removed
    // until their segments are removed first
    std::set<VertexPtr> vertices_to_remove = find_vertices(graph, cluster);
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
    
    auto seg = create_segment_for_cluster(cluster, dv, path_points, 1);
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
    // auto seg = create_segment_for_cluster(cluster, dv, path_points, 1);
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
