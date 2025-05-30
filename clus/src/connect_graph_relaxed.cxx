
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
                                double cluster_t0 = cluster.get_flash().time();
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
                                double cluster_t0 = cluster.get_flash().time();
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
                                double cluster_t0 = cluster.get_flash().time();
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
                    add_edge(index1, index2, std::get<2>(index_index_dis[j][k]), temp_graph);
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
                    add_edge(
                        index1, index2,
                        std::min(std::get<2>(index_index_dis_dir1[j][k]), std::get<2>(index_index_dis_dir2[j][k])),
                        temp_graph);
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
                /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
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
                    /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
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
                    /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                }
            }

        }  // k
    }  // j
    
}

