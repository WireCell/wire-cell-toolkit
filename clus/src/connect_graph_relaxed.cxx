
#include "WireCellClus/Graphs.h"
#include "WireCellClus/IPCTransform.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Grouping.h"

#include "connect_graphs.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Graphs;
using namespace WireCell::Clus::Facade;

void Graphs::connect_graph_relaxed(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv, 
    IPCTransformSet::pointer pcts,
    Weighted::Graph& graph)    
{
    const bool use_ctpc = true;
    const auto* grouping = cluster.grouping();

    // Get all the wire plane IDs from the grouping
    const auto& wpids = grouping->wpids();

    // Key: pair<APA, face>, Value: drift_dir, angle_u, angle_v, angle_w
    std::map<WirePlaneId , std::tuple<geo_point_t, double, double, double>> wpid_params;
    std::map<WirePlaneId, geo_point_t> wpid_U_dir;
    std::map<WirePlaneId, geo_point_t> wpid_V_dir;
    std::map<WirePlaneId, geo_point_t> wpid_W_dir;
    std::set<int> apas;
    for (const auto& wpid : wpids) {
        int apa = wpid.apa();
        int face = wpid.face();
        apas.insert(apa);

        // Create wpids for all three planes with this APA and face
        WirePlaneId wpid_u(kUlayer, face, apa);
        WirePlaneId wpid_v(kVlayer, face, apa);
        WirePlaneId wpid_w(kWlayer, face, apa);
     
        // Get drift direction based on face orientation
        int face_dirx = dv->face_dirx(wpid_u);
        geo_point_t drift_dir(face_dirx, 0, 0);
        
        // Get wire directions for all planes
        Vector wire_dir_u = dv->wire_direction(wpid_u);
        Vector wire_dir_v = dv->wire_direction(wpid_v);
        Vector wire_dir_w = dv->wire_direction(wpid_w);

        // Calculate angles
        double angle_u = std::atan2(wire_dir_u.z(), wire_dir_u.y());
        double angle_v = std::atan2(wire_dir_v.z(), wire_dir_v.y());
        double angle_w = std::atan2(wire_dir_w.z(), wire_dir_w.y());

        wpid_params[wpid] = std::make_tuple(drift_dir, angle_u, angle_v, angle_w);
        wpid_U_dir[wpid] = geo_point_t(0, cos(angle_u), sin(angle_u));
        wpid_V_dir[wpid] = geo_point_t(0, cos(angle_v), sin(angle_v));
        wpid_W_dir[wpid] = geo_point_t(0, cos(angle_w), sin(angle_w));
    }

    // this drift direction is only used to calculate isochronous case, so this is OK ...
    const geo_vector_t drift_dir_abs(1, 0, 0); 


    // Form connected components
    std::vector<int> component(num_vertices(graph));
    const size_t num = connected_components(graph, &component[0]);
    
    if (num <= 1) return;

    // Create point clouds using connected components
    std::vector<std::shared_ptr<Simple3DPointCloud>> pt_clouds;
    std::vector<std::vector<size_t>> pt_clouds_global_indices;
    
    // Create ordered components  
    std::vector<ComponentInfo> ordered_components;
    ordered_components.reserve(component.size());
    for (size_t i = 0; i < component.size(); ++i) {
        ordered_components.emplace_back(i);
    }
    
    // Assign vertices to components
    for (size_t i = 0; i < component.size(); ++i) {
        ordered_components[component[i]].add_vertex(i);
    }

    // Sort components by minimum vertex index
    std::sort(ordered_components.begin(), ordered_components.end(),
        [](const ComponentInfo& a, const ComponentInfo& b) {
            return a.min_vertex < b.min_vertex;
        });

    // Create point clouds for each component
    const auto& points = cluster.points();
    for (const auto& comp : ordered_components) {
        auto pt_cloud = std::make_shared<Simple3DPointCloud>();

        std::vector<size_t> global_indices;
        
        for (size_t vertex_idx : comp.vertex_indices) {
            pt_cloud->add({points[0][vertex_idx], points[1][vertex_idx], points[2][vertex_idx]});
            global_indices.push_back(vertex_idx);
        }
        pt_clouds.push_back(pt_cloud);
        pt_clouds_global_indices.push_back(global_indices);
    }

    // Initialize distance metrics 
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis(num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_mst(num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir1(num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir2(num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir_mst(num, std::vector<std::tuple<int, int, double>>(num));

    // Initialize all distances to inf
    for (size_t j = 0; j != num; j++) {
        for (size_t k = 0; k != num; k++) {
            index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_mst[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir_mst[j][k] = std::make_tuple(-1, -1, 1e9);
        }
    }

    // Calculate distances between components
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            // Get closest points between components
            index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));

            // Skip small clouds
            if ((num < 100 && pt_clouds.at(j)->get_num_points() > 100 && pt_clouds.at(k)->get_num_points() > 100 &&
                 (pt_clouds.at(j)->get_num_points() + pt_clouds.at(k)->get_num_points()) > 400) ||
                (pt_clouds.at(j)->get_num_points() > 500 && pt_clouds.at(k)->get_num_points() > 500)) {
                
                // Get closest points and calculate directions
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));

                geo_vector_t dir1 = cluster.vhough_transform(p1, 30 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), 
                                                pt_clouds_global_indices.at(j));
                geo_vector_t dir2 = cluster.vhough_transform(p2, 30 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k),
                                                pt_clouds_global_indices.at(k)); 
                dir1 = dir1 * -1;
                dir2 = dir2 * -1;

                std::pair<int, double> result1 = pt_clouds.at(k)->get_closest_point_along_vec(p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);

                if (result1.first >= 0) {
                    index_index_dis_dir1[j][k] = std::make_tuple(std::get<0>(index_index_dis[j][k]), 
                                                                result1.first, result1.second);
                }

                std::pair<int, double> result2 = pt_clouds.at(j)->get_closest_point_along_vec(p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm); 

                if (result2.first >= 0) {
                    index_index_dis_dir2[j][k] = std::make_tuple(result2.first,
                                                                std::get<1>(index_index_dis[j][k]), 
                                                                result2.second);
                }
            }
            // Now check the path 

            {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis[j][k])));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis/step_dis + 1;

                

                // Track different types of "bad" points
                int num_bad[4] = {0,0,0,0};   // more than one of three are bad
                int num_bad1[4] = {0,0,0,0};  // at least one of three are bad
                int num_bad2[3] = {0,0,0};    // number of dead channels

                // Check points along path
                for (int ii = 0; ii != num_steps; ii++) {
                    geo_point_t test_p(
                        p1.x() + (p2.x() - p1.x())/num_steps*(ii + 1),
                        p1.y() + (p2.y() - p1.y())/num_steps*(ii + 1),
                        p1.z() + (p2.z() - p1.z())/num_steps*(ii + 1)
                    );

                    // Test point quality using grouping parameters
                    std::vector<int> scores;
                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa()!=-1){
                            geo_point_t test_p_raw = test_p;
                            if (cluster.get_default_scope().hash() != cluster.get_raw_scope().hash()){
                                const auto transform = pcts->pc_transform(cluster.get_scope_transform());
                                double cluster_t0 = cluster.get_cluster_t0();
                                test_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            scores = grouping->test_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            
                            // Check overall quality
                            if (scores[0] + scores[3] + scores[1] + scores[4] + (scores[2]+scores[5])*2 < 3) {
                                num_bad[0]++;
                            }
                            if (scores[0]+scores[3]==0) num_bad[1]++;
                            if (scores[1]+scores[4]==0) num_bad[2]++;
                            if (scores[2]+scores[5]==0) num_bad[3]++;

                            if (scores[3]!=0) num_bad2[0]++;
                            if (scores[4]!=0) num_bad2[1]++;
                            if (scores[5]!=0) num_bad2[2]++;
                            
                            if (scores[0] + scores[3] + scores[1] + scores[4] + (scores[2]+scores[5]) < 3) {
                                num_bad1[0]++;
                            }
                            if (scores[0]+scores[3]==0) num_bad1[1]++;
                            if (scores[1]+scores[4]==0) num_bad1[2]++;
                            if (scores[2]+scores[5]==0) num_bad1[3]++;
                        }
                    }
                }

                auto test_wpid = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv);

                // Calculate angles between directions
                geo_vector_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
                geo_vector_t tempV5;

                double angle1 = tempV1.angle(wpid_U_dir.at(test_wpid)); 
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2)) * sin(angle1),
                        0);
                angle1 = tempV5.angle(drift_dir_abs);

                double angle2 = tempV1.angle(wpid_V_dir.at(test_wpid));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2)) * sin(angle2),
                        0);
                angle2 = tempV5.angle(drift_dir_abs);

                double angle1p = tempV1.angle(wpid_W_dir.at(test_wpid));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2)) * sin(angle1p),
                        0); 
                angle1p = tempV5.angle(drift_dir_abs);

                tempV5.set(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
                double angle3 = tempV5.angle(drift_dir_abs);

                bool flag_strong_check = true;

                // Define constants for readability
                constexpr double pi = 3.141592653589793;
                constexpr double perp_angle_tol = 10.0/180.0*pi;
                constexpr double wire_angle_tol = 12.5/180.0*pi;
                constexpr double perp_angle = pi/2.0;
                constexpr double invalid_dist = 1e9;

                if (fabs(angle3 - perp_angle) < perp_angle_tol) {
                    geo_vector_t tempV2 = cluster.vhough_transform(p1, 15*units::cm);
                    geo_vector_t tempV3 = cluster.vhough_transform(p2, 15*units::cm);
                    
                    if (fabs(tempV2.angle(drift_dir_abs) - perp_angle) < perp_angle_tol &&
                        fabs(tempV3.angle(drift_dir_abs) - perp_angle) < perp_angle_tol) {
                        flag_strong_check = false;
                    }
                }
                else if (angle1 < wire_angle_tol || angle2 < wire_angle_tol || angle1p < wire_angle_tol) {
                    flag_strong_check = false;
                }

                // Helper function to check if ratio exceeds threshold
                auto exceeds_ratio = [](int val, int steps, double ratio = 0.75) {
                    return val >= ratio * steps;
                };

                // Helper function to invalidate distance
                auto invalidate_distance = [&]() {
                    index_index_dis[j][k] = std::make_tuple(-1, -1, invalid_dist);
                };

                if (flag_strong_check) {
                    if (num_bad1[0] > 7 || (num_bad1[0] > 2 && exceeds_ratio(num_bad1[0], num_steps))) {
                        invalidate_distance();
                    }
                }
                else {
                    bool parallel_angles = (angle1 < wire_angle_tol && angle2 < wire_angle_tol) ||
                                        (angle1p < wire_angle_tol && angle1 < wire_angle_tol) ||
                                        (angle1p < wire_angle_tol && angle2 < wire_angle_tol);

                    if (parallel_angles) {
                        if (num_bad[0] > 7 || (num_bad[0] > 2 && exceeds_ratio(num_bad[0], num_steps))) {
                            invalidate_distance();
                        }
                    }
                    else if (angle1 < wire_angle_tol) {
                        int sum_bad = num_bad[2] + num_bad[3];
                        if (sum_bad > 9 || (sum_bad > 2 && exceeds_ratio(sum_bad, num_steps)) || num_bad[3] >= 3) {
                            invalidate_distance();
                        }
                    }
                    else if (angle2 < wire_angle_tol) {
                        int sum_bad = num_bad[1] + num_bad[3];
                        if (sum_bad > 9 || (sum_bad > 2 && exceeds_ratio(sum_bad, num_steps)) || num_bad[3] >= 3) {
                            invalidate_distance();
                        }
                    }
                    else if (angle1p < wire_angle_tol) {
                        int sum_bad = num_bad[2] + num_bad[1];
                        if (sum_bad > 9 || (sum_bad > 2 && exceeds_ratio(sum_bad, num_steps))) {
                            invalidate_distance();
                        }
                    }
                    else if (num_bad[0] > 7 || (num_bad[0] > 2 && exceeds_ratio(num_bad[0], num_steps))) {
                        invalidate_distance();
                    }
                }
            }

            // Now check path again ... 
            if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir1[j][k])); 
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k])));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir1[j][k])); 
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + 
                                pow(p1.y() - p2.y(), 2) + 
                                pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis/step_dis + 1;
                int num_bad = 0;
                int num_bad1 = 0;

                // Check intermediate points along path
                for (int ii = 0; ii != num_steps; ii++) {
                    geo_point_t test_p(
                        p1.x() + (p2.x() - p1.x())/num_steps*(ii + 1),
                        p1.y() + (p2.y() - p1.y())/num_steps*(ii + 1),
                        p1.z() + (p2.z() - p1.z())/num_steps*(ii + 1)
                    );

                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa()!=-1){
                            geo_point_t test_p_raw = test_p;
                            if (cluster.get_default_scope().hash() != cluster.get_raw_scope().hash()){
                                const auto transform = pcts->pc_transform(cluster.get_scope_transform());
                                double cluster_t0 = cluster.get_cluster_t0();
                                test_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }                            
                            const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            if (!good_point) {
                                num_bad++;
                            }
                            if (!grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face(), 0.6*units::cm, 1, 0)) {
                                num_bad1++;
                            }
                        }
                    }
                }
                
                auto test_wpid = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv);

                // Calculate angles
                geo_vector_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
                geo_vector_t tempV5;
                
                double angle1 = tempV1.angle(wpid_U_dir.at(test_wpid));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle1),
                        0);
                angle1 = tempV5.angle(drift_dir_abs);
                
                double angle2 = tempV1.angle(wpid_V_dir.at(test_wpid));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle2),
                        0);
                angle2 = tempV5.angle(drift_dir_abs);
                
                tempV5.set(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
                double angle3 = tempV5.angle(drift_dir_abs);
                
                double angle1p = tempV1.angle(wpid_W_dir.at(test_wpid));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle1p),
                        0);
                angle1p = tempV5.angle(drift_dir_abs);

                const double pi = 3.141592653589793;
                if (fabs(angle3 - pi/2) < 10.0/180.0*pi || 
                    angle1 < 12.5/180.0*pi ||
                    angle2 < 12.5/180.0*pi || 
                    angle1p < 7.5/180.0*pi) {
                    // Parallel or prolonged case
                    if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75*num_steps)) {
                        index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }
                else {
                    if (num_bad1 > 7 || (num_bad1 > 2 && num_bad1 >= 0.75*num_steps)) {
                        index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }
            }

            //Now check path again ... 
            // Now check the path...
            if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir2[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k])));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir2[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + 
                                pow(p1.y() - p2.y(), 2) + 
                                pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis/step_dis + 1;
                int num_bad = 0;
                int num_bad1 = 0;

                // Check points along path
                for (int ii = 0; ii != num_steps; ii++) {
                    geo_point_t test_p(
                        p1.x() + (p2.x() - p1.x())/num_steps*(ii + 1),
                        p1.y() + (p2.y() - p1.y())/num_steps*(ii + 1),
                        p1.z() + (p2.z() - p1.z())/num_steps*(ii + 1)
                    );

                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa()!=-1){
                            geo_point_t test_p_raw = test_p;
                            if (cluster.get_default_scope().hash() != cluster.get_raw_scope().hash()){
                                const auto transform = pcts->pc_transform(cluster.get_scope_transform());
                                double cluster_t0 = cluster.get_cluster_t0();
                                test_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            if (!good_point) {
                                num_bad++;
                            }
                            if (!grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face(), 0.6*units::cm, 1, 0)) {
                                num_bad1++;
                            }
                        }
                    }
                }

                auto test_wpid = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv);

                // Calculate angles between directions
                geo_vector_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
                geo_vector_t tempV5;

                double angle1 = tempV1.angle(wpid_U_dir.at(test_wpid));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle1),
                        0);
                angle1 = tempV5.angle(drift_dir_abs);

                double angle2 = tempV1.angle(wpid_V_dir.at(test_wpid));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle2),
                        0);
                angle2 = tempV5.angle(drift_dir_abs);

                tempV5.set(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
                double angle3 = tempV5.angle(drift_dir_abs);

                double angle1p = tempV1.angle(wpid_W_dir.at(test_wpid));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle1p),
                        0);
                angle1p = tempV5.angle(drift_dir_abs);

                const double pi = 3.141592653589793;
                bool is_parallel = fabs(angle3 - pi/2) < 10.0/180.0*pi || 
                                angle1 < 12.5/180.0*pi ||
                                angle2 < 12.5/180.0*pi || 
                                angle1p < 7.5/180.0*pi;

                if (is_parallel) {
                    // Parallel or prolonged case
                    if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75*num_steps)) {
                        index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }
                else {
                    if (num_bad1 > 7 || (num_bad1 > 2 && num_bad1 >= 0.75*num_steps)) {
                        index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }
            }
        }
    }

    // deal with MST of first type
    {
        Weighted::Graph temp_graph(num);
        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j;
                int index2 = k;
                if (std::get<0>(index_index_dis[j][k]) >= 0) {
                    if (!boost::edge(index1, index2, temp_graph).second) {
                        add_edge(index1, index2, std::get<2>(index_index_dis[j][k]), temp_graph);
                    }
                }
            }
        }

        // Process MST
        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_mst);
    }

    // MST of the direction ...
    {
        Weighted::Graph temp_graph(num);

        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j;
                int index2 = k;
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0 || std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    if (!boost::edge(index1, index2, temp_graph).second) {
                        // Add edge with minimum distance from both directions
                    add_edge(
                        index1, index2,
                        std::min(std::get<2>(index_index_dis_dir1[j][k]), std::get<2>(index_index_dis_dir2[j][k])),
                        temp_graph);
                    }
                }
            }
        }

        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_dir_mst);

    }

    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            if (std::get<2>(index_index_dis[j][k]) < 3 * units::cm) {
                index_index_dis_mst[j][k] = index_index_dis[j][k];
            }

            // establish the path ...
            if (std::get<0>(index_index_dis_mst[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_mst[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_mst[j][k]));
                float dis;
                if (std::get<2>(index_index_dis_mst[j][k]) > 5 * units::cm) {
                    dis = std::get<2>(index_index_dis_mst[j][k]);
                }
                else {
                    dis = std::get<2>(index_index_dis_mst[j][k]);
                }
                if (!boost::edge(gind1, gind2, graph).second) {
                    /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                }
            }

            if (std::get<0>(index_index_dis_dir_mst[j][k]) >= 0) {
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k]));
                    float dis;
                    if (std::get<2>(index_index_dis_dir1[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir1[j][k]) * 1.1;
                    }
                    else {
                        dis = std::get<2>(index_index_dis_dir1[j][k]);
                    }
                    if(!boost::edge(gind1, gind2, graph).second) {
                        /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                    }
                }
                if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k]));
                    float dis;
                    if (std::get<2>(index_index_dis_dir2[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir2[j][k]) * 1.1;
                    }
                    else {
                        dis = std::get<2>(index_index_dis_dir2[j][k]);
                    }
                    // }
                    if (!boost::edge(gind1, gind2, graph).second) {
                        /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                    }
                }
            }

        }  // k
    }  // j
    
}

 bool Graphs::check_connectivity(const Facade::Cluster& cluster,  IDetectorVolumes::pointer dv, IPCTransformSet::pointer pcts, std::tuple<int, int, double>& index_index_dis, std::shared_ptr<Facade::Simple3DPointCloud> pc1, std::vector<size_t> pc1_global_index, std::shared_ptr<Facade::Simple3DPointCloud> pc2, std::vector<size_t> pc2_global_index,
    double step_size, bool flag_strong_check){
    
    // Check if indices are valid
    if (std::get<0>(index_index_dis) == -1 || std::get<1>(index_index_dis) == -1) return false;
    
    // Get points from cluster
    // const auto& points = cluster.points();
    
    // Get the two points from the point cloud
    int idx1 = std::get<0>(index_index_dis);
    int idx2 = std::get<1>(index_index_dis);
    
    geo_point_t p1= pc1->point(idx1);
    geo_point_t p2= pc2->point(idx2);
    
    // Get wire plane IDs for the two points
    auto wpid_p1 = cluster.wire_plane_id(pc1_global_index.at(idx1));
    auto wpid_p2 = cluster.wire_plane_id(pc2_global_index.at(idx2));
    auto wpid_pc = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv);
    
    int apa1 = wpid_p1.apa();
    int face1 = wpid_p1.face();
    int apa2 = wpid_p2.apa();
    int face2 = wpid_p2.face();
    int apa3 = wpid_pc.apa();
    int face3 = wpid_pc.face();
    
    // Get grouping for CTPC checks
    const auto* grouping = cluster.grouping();
    
    // Calculate directions using VHoughTrans equivalent (vhough_transform)
    // Use point clouds pc1 and pc2 for local direction calculation
    geo_vector_t dir1 = cluster.vhough_transform(p1, 15*units::cm, Cluster::HoughParamSpace::theta_phi, pc1, pc1_global_index);
    dir1 = dir1 * -1;
    
    geo_vector_t dir2 = cluster.vhough_transform(p2, 15*units::cm, Cluster::HoughParamSpace::theta_phi, pc2, pc2_global_index);
    dir2 = dir2 * -1;
    
    geo_vector_t dir3(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
    
    // Check directions using the check_direction function
    std::vector<bool> flag_1 = check_direction(cluster, dir1, apa1, face1);
    std::vector<bool> flag_2 = check_direction(cluster, dir2, apa2, face2);
    
    // For dir3, use either apa/face from p1 or p2 (use p1 as reference)
    std::vector<bool> flag_3 = check_direction(cluster, dir3, apa3, face3);
    
    bool flag_prolonged_u = false;
    bool flag_prolonged_v = false;
    bool flag_prolonged_w = false;
    bool flag_parallel = false;
    
    // Check if prolonged along wire directions or parallel to drift
    if (flag_3.at(0) && (flag_1.at(0) || flag_2.at(0))) flag_prolonged_u = true;
    if (flag_3.at(1) && (flag_1.at(1) || flag_2.at(1))) flag_prolonged_v = true;
    if (flag_3.at(2) && (flag_1.at(2) || flag_2.at(2))) flag_prolonged_w = true;
    if (flag_3.at(3) && (flag_1.at(3) && flag_2.at(3))) flag_parallel = true;
    
    // Calculate distance and number of steps
    double dis = std::sqrt(std::pow(p1.x() - p2.x(), 2) + 
                          std::pow(p1.y() - p2.y(), 2) + 
                          std::pow(p1.z() - p2.z(), 2));
    int num_steps = std::round(dis / step_size);
    
    if (num_steps == 0) num_steps = 1;
    
    int num_bad[5] = {0, 0, 0, 0, 0};
    
    double radius_cut = 0.6 * units::cm;
    if (step_size < radius_cut) radius_cut = step_size;
    
    // Check points along the path
    for (int i = 0; i != num_steps; i++) {
        geo_point_t test_p(
            p1.x() + (p2.x() - p1.x()) / (num_steps + 1.0) * (i + 1),
            p1.y() + (p2.y() - p1.y()) / (num_steps + 1.0) * (i + 1),
            p1.z() + (p2.z() - p1.z()) / (num_steps + 1.0) * (i + 1)
        );
        
        // Get wire plane ID for test point
        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
        
        if (test_wpid.apa() == -1) continue;
        
        // Transform point if needed
        geo_point_t test_p_raw = test_p;
        if (cluster.get_default_scope().hash() != cluster.get_raw_scope().hash()) {
            const auto transform = pcts->pc_transform(cluster.get_scope_transform());
            double cluster_t0 = cluster.get_cluster_t0();
            test_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
        }
        
        // Test point quality with appropriate radius
        double test_radius;
        if (i == 0 || i + 1 == num_steps) {
            test_radius = dis / (num_steps + 1.0) * 0.98;
        } else {
            if (flag_strong_check) {
                test_radius = dis / (num_steps + 1.0);
            } else {
                test_radius = radius_cut;
            }
        }
        
        // Get detailed scores for this point
        std::vector<int> scores = grouping->test_good_point(test_p_raw, test_wpid.apa(), test_wpid.face(), test_radius);
        
        int num_bad_details = 0;
        
        // Check U plane (indices 0=live, 3=dead)
        if (scores.at(0) + scores.at(3) == 0) {
            if (!flag_prolonged_u) num_bad[0]++;
            num_bad_details++;
        }
        
        // Check V plane (indices 1=live, 4=dead)
        if (scores.at(1) + scores.at(4) == 0) {
            if (!flag_prolonged_v) num_bad[1]++;
            num_bad_details++;
        }
        
        // Check W plane (collection, indices 2=live, 5=dead)
        if (scores.at(2) + scores.at(5) == 0) {
            if (!flag_prolonged_w) num_bad[2]++;
            num_bad_details++;
        }
        
        // Count overall bad points
        if (flag_parallel) {
            // Parallel case: more than one plane bad
            if (num_bad_details > 1) num_bad[3]++;
        } else {
            // Non-parallel: any plane bad
            if (num_bad_details > 0) num_bad[3]++;
        }
    }
    
    // Strong check - very strict criteria
    if (flag_strong_check && ((num_bad[0] + num_bad[1] + num_bad[2]) > 0 || num_bad[3] >= 2)) {
        return false;
    }
    
    // Prolonged case - allow some bad points but not too many
    if (num_bad[0] <= 2 && num_bad[1] <= 2 && num_bad[2] <= 2 &&
        (num_bad[0] + num_bad[1] + num_bad[2] <= 3) && 
        num_bad[0] < 0.1 * num_steps && 
        num_bad[1] < 0.1 * num_steps && 
        num_bad[2] < 0.1 * num_steps &&
        (num_bad[0] + num_bad[1] + num_bad[2]) < 0.15 * num_steps) {
        
        // Special case: if prolonged in all three directions, check overall quality
        if (flag_prolonged_u && flag_prolonged_v && flag_prolonged_w) {
            if (num_bad[3] >= 0.6 * num_steps) return false;
        }
        
        return true;
    } 
    // Alternative case - overall good quality
    else if (num_bad[3] <= 2 && num_bad[3] < 0.1 * num_steps) {
        return true;
    }
    
    return false;
 }

void Graphs::connect_graph_relaxed_pid(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv, 
        IPCTransformSet::pointer pcts,
        Weighted::Graph& graph){
    
    // const auto* grouping = cluster.grouping();
    const geo_vector_t drift_dir_abs(1, 0, 0);
    
    // Form connected components
    std::vector<int> component(num_vertices(graph));
    const size_t num = connected_components(graph, &component[0]);
    
    if (num <= 1) return;
    
    // Create point clouds using connected components
    std::vector<std::shared_ptr<Simple3DPointCloud>> pt_clouds;
    std::vector<std::vector<size_t>> pt_clouds_global_indices;
    
    // Create ordered components
    std::vector<ComponentInfo> ordered_components;
    ordered_components.reserve(component.size());
    for (size_t i = 0; i < component.size(); ++i) {
        ordered_components.emplace_back(i);
    }
    
    // Assign vertices to components
    for (size_t i = 0; i < component.size(); ++i) {
        ordered_components[component[i]].add_vertex(i);
    }
    
    // Sort components by minimum vertex index
    std::sort(ordered_components.begin(), ordered_components.end(),
        [](const ComponentInfo& a, const ComponentInfo& b) {
            return a.min_vertex < b.min_vertex;
        });
    
    // Create point clouds for each component
    const auto& points = cluster.points();
    for (const auto& comp : ordered_components) {
        auto pt_cloud = std::make_shared<Simple3DPointCloud>();
        std::vector<size_t> global_indices;
        
        for (size_t vertex_idx : comp.vertex_indices) {
            pt_cloud->add({points[0][vertex_idx], points[1][vertex_idx], points[2][vertex_idx]});
            global_indices.push_back(vertex_idx);
        }
        pt_clouds.push_back(pt_cloud);
        pt_clouds_global_indices.push_back(global_indices);
    }
    
    // Initialize distance metrics
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis(num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir1(num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir2(num, std::vector<std::tuple<int, int, double>>(num));
    
    // Initialize all distances to inf
    for (size_t j = 0; j != num; j++) {
        for (size_t k = 0; k != num; k++) {
            index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
        }
    }
    
    // Calculate distances between components with connectivity checks
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            // Get closest points between components
            std::tuple<int, int, double> temp_index_index_dis = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));
            
            if (std::get<0>(temp_index_index_dis) != -1) {
                index_index_dis[j][k] = temp_index_index_dis;
                
                // Check connectivity
                bool flag = check_connectivity(cluster, dv, pcts, index_index_dis[j][k], 
                                              pt_clouds.at(j), pt_clouds_global_indices.at(j), pt_clouds.at(k), pt_clouds_global_indices.at(k));
                
                // Special case for very close distances with large point clouds
                if (std::get<2>(temp_index_index_dis) <= 0.9 * units::cm && 
                    pt_clouds.at(j)->get_num_points() > 200 && 
                    pt_clouds.at(k)->get_num_points() > 200) {
                    
                    if (!flag) {
                        // Try to find better connection points nearby
                        geo_point_t test_p1 = pt_clouds.at(k)->point(std::get<1>(temp_index_index_dis));
                        geo_point_t test_p2 = pt_clouds.at(j)->point(std::get<0>(temp_index_index_dis));
                        
                        auto temp_wcps1 = 
                            pt_clouds.at(j)->get_closest_wcpoints_radius(test_p1, 
                                std::get<2>(temp_index_index_dis) + 0.9 * units::cm);
                        auto temp_wcps2 = 
                            pt_clouds.at(k)->get_closest_wcpoints_radius(test_p2, 
                                std::get<2>(temp_index_index_dis) + 0.9 * units::cm);
                        
                        // Try different point combinations
                        for (size_t kk1 = 0; kk1 < temp_wcps1.size() && !flag; kk1++) {
                            for (size_t kk2 = 0; kk2 < temp_wcps2.size() && !flag; kk2++) {
                                double dis = std::sqrt(
                                    std::pow(temp_wcps1[kk1].second.x() - temp_wcps2[kk2].second.x(), 2) +
                                    std::pow(temp_wcps1[kk1].second.y() - temp_wcps2[kk2].second.y(), 2) +
                                    std::pow(temp_wcps1[kk1].second.z() - temp_wcps2[kk2].second.z(), 2));
                                
                                std::tuple<int, int, double> temp_tuple = 
                                    std::make_tuple(temp_wcps2[kk2].first, temp_wcps1[kk1].first, dis);
                                
                                if (check_connectivity(cluster, dv, pcts, temp_tuple, 
                                                      pt_clouds.at(j), pt_clouds_global_indices.at(j), pt_clouds.at(k), pt_clouds_global_indices.at(k), 
                                                      0.3 * units::cm, true)) {
                                    flag = true;
                                    index_index_dis[j][k] = temp_tuple;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if (!flag) {
                    index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
                }
                index_index_dis[k][j] = index_index_dis[j][k];
                
                // Calculate directional connections
                if (std::get<0>(temp_index_index_dis) != -1) {
                    geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(temp_index_index_dis));
                    geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(temp_index_index_dis));
                    
                    // Direction from p1
                    geo_vector_t dir1 = cluster.vhough_transform(p1, 30 * units::cm, 
                        Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                    dir1 = dir1 * -1;
                    
                    std::pair<int, double> result1 = pt_clouds.at(k)->get_closest_point_along_vec(
                        p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                    
                    // If no result and perpendicular to drift, try longer hough
                    if (result1.first < 0 && 
                        std::fabs(dir1.angle(drift_dir_abs) * 180.0 / M_PI - 90.0) < 10.0) {
                        if (std::fabs(dir1.angle(drift_dir_abs) * 180.0 / M_PI - 90.0) < 5.0)
                            dir1 = cluster.vhough_transform(p1, 80 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                        else if (std::fabs(dir1.angle(drift_dir_abs) * 180.0 / M_PI - 90.0) < 10.0)
                            dir1 = cluster.vhough_transform(p1, 50 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                        dir1 = dir1 * -1;
                        result1 = pt_clouds.at(k)->get_closest_point_along_vec(
                            p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                    }
                    
                    if (result1.first >= 0) {
                        index_index_dis_dir1[j][k] = std::make_tuple(
                            std::get<0>(index_index_dis[j][k]), result1.first, result1.second);
                        
                        if (!check_connectivity(cluster, dv, pcts, index_index_dis_dir1[j][k], 
                                               pt_clouds.at(j), pt_clouds_global_indices.at(j), pt_clouds.at(k), pt_clouds_global_indices.at(k))) {
                            index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                        }
                        index_index_dis_dir1[k][j] = index_index_dis_dir1[j][k];
                    }
                    
                    // Direction from p2
                    geo_vector_t dir2 = cluster.vhough_transform(p2, 30 * units::cm, 
                        Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                    dir2 = dir2 * -1;
                    
                    std::pair<int, double> result2 = pt_clouds.at(j)->get_closest_point_along_vec(
                        p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                    
                    if (result2.first < 0 && 
                        std::fabs(dir2.angle(drift_dir_abs) * 180.0 / M_PI - 90.0) < 10.0) {
                        if (std::fabs(dir2.angle(drift_dir_abs) * 180.0 / M_PI - 90.0) < 5.0)
                            dir2 = cluster.vhough_transform(p2, 80 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                        else if (std::fabs(dir2.angle(drift_dir_abs) * 180.0 / M_PI - 90.0) < 10.0)
                            dir2 = cluster.vhough_transform(p2, 50 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                        dir2 = dir2 * -1;
                        result2 = pt_clouds.at(j)->get_closest_point_along_vec(
                            p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                    }
                    
                    if (result2.first >= 0) {
                        index_index_dis_dir2[j][k] = std::make_tuple(
                            result2.first, std::get<1>(index_index_dis[j][k]), result2.second);
                        
                        if (!check_connectivity(cluster, dv, pcts, index_index_dis_dir2[j][k], 
                                               pt_clouds.at(j), pt_clouds_global_indices.at(j), pt_clouds.at(k), pt_clouds_global_indices.at(k))) {
                            index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                        }
                        index_index_dis_dir2[k][j] = index_index_dis_dir2[j][k];
                    }
                }
            }
        }
    }
    
    // Examine middle path for all three connection types
    double step_dis = 1.0 * units::cm;
    
    auto examine_middle_path = [&](std::vector<std::vector<std::tuple<int, int, double>>>& index_dis_array) {
        std::map<std::pair<int, int>, std::set<int>> map_add_connections;
        
        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                if (std::get<0>(index_dis_array[j][k]) >= 0) {
                    int idx1 = std::get<0>(index_dis_array[j][k]);
                    int idx2 = std::get<1>(index_dis_array[j][k]);
                    
                    geo_point_t wp1(points[0][pt_clouds_global_indices[j][idx1]], 
                                   points[1][pt_clouds_global_indices[j][idx1]], 
                                   points[2][pt_clouds_global_indices[j][idx1]]);
                    geo_point_t wp2(points[0][pt_clouds_global_indices[k][idx2]], 
                                   points[1][pt_clouds_global_indices[k][idx2]], 
                                   points[2][pt_clouds_global_indices[k][idx2]]);
                    
                    double length = std::sqrt(
                        std::pow(wp1.x() - wp2.x(), 2) + 
                        std::pow(wp1.y() - wp2.y(), 2) + 
                        std::pow(wp1.z() - wp2.z(), 2));
                    
                    if (length > 3 * units::cm) {
                        std::set<int> connections;
                        int ncount = std::round(length / step_dis);
                        
                        for (int qx = 1; qx < ncount; qx++) {
                            geo_point_t test_p(
                                wp1.x() + (wp2.x() - wp1.x()) * qx / ncount,
                                wp1.y() + (wp2.y() - wp1.y()) * qx / ncount,
                                wp1.z() + (wp2.z() - wp1.z()) * qx / ncount);
                            
                            for (size_t qx1 = 0; qx1 != num; qx1++) {
                                if (qx1 == j || qx1 == k) continue;
                                
                                // Skip small components, check only >= 50 points
                                if (pt_clouds.at(qx1)->get_closest_dis(test_p) < 0.6 * units::cm &&
                                    pt_clouds.at(qx1)->get_num_points() >= 50) {
                                    connections.insert(qx1);
                                }
                            }
                        }
                        
                        if (!connections.empty()) {
                            map_add_connections[std::make_pair(j, k)] = connections;
                        }
                    }
                }
            }
        }
        
        // Iteratively disconnect paths blocked by intermediate components
        bool flag_continue = true;
        while (flag_continue) {
            flag_continue = false;
            std::set<std::pair<int, int>> used_pairs;
            
            for (auto it = map_add_connections.begin(); it != map_add_connections.end(); it++) {
                int j = it->first.first;
                int k = it->first.second;
                bool flag_disconnect = true;
                
                for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++) {
                    int qx = *it1;
                    if ((std::get<0>(index_index_dis[j][qx]) != -1 || 
                         std::get<0>(index_index_dis_dir1[j][qx]) != -1 || 
                         std::get<0>(index_index_dis_dir2[j][qx]) != -1) &&
                        (std::get<0>(index_index_dis[k][qx]) != -1 || 
                         std::get<0>(index_index_dis_dir1[k][qx]) != -1 || 
                         std::get<0>(index_index_dis_dir2[k][qx]) != -1)) {
                        flag_disconnect = false;
                        break;
                    }
                }
                
                if (flag_disconnect) {
                    flag_continue = true;
                    index_dis_array[j][k] = std::make_tuple(-1, -1, 1e9);
                    index_dis_array[k][j] = index_dis_array[j][k];
                    used_pairs.insert(it->first);
                }
            }
            
            for (auto it = used_pairs.begin(); it != used_pairs.end(); it++) {
                map_add_connections.erase(*it);
            }
        }
    };
    
    // Examine all three connection types
    examine_middle_path(index_index_dis);
    examine_middle_path(index_index_dis_dir1);
    examine_middle_path(index_index_dis_dir2);
    
    // Final assembly: add edges to graph
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            // Add closest distance connections
            if (std::get<0>(index_index_dis[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis[j][k]));
                
                if (!boost::edge(gind1, gind2, graph).second) {
                    add_edge(gind1, gind2, std::get<2>(index_index_dis[j][k]), graph);
                }
            }
            
            // Add directional connection 1
            if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k]));
                
                float dis;
                if (std::get<2>(index_index_dis_dir1[j][k]) > 5 * units::cm) {
                    dis = std::get<2>(index_index_dis_dir1[j][k]) * 1.2;
                } else {
                    dis = std::get<2>(index_index_dis_dir1[j][k]);
                }
                
                if (!boost::edge(gind1, gind2, graph).second) {
                    add_edge(gind1, gind2, dis, graph);
                }
            }
            
            // Add directional connection 2
            if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k]));
                
                float dis;
                if (std::get<2>(index_index_dis_dir2[j][k]) > 5 * units::cm) {
                    dis = std::get<2>(index_index_dis_dir2[j][k]) * 1.2;
                } else {
                    dis = std::get<2>(index_index_dis_dir2[j][k]);
                }
                
                if (!boost::edge(gind1, gind2, graph).second) {
                    add_edge(gind1, gind2, dis, graph);
                }
            }
        }
    }
}