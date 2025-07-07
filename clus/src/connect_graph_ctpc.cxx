#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Grouping.h"

#include "connect_graphs.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

void Graphs::connect_graph_ctpc(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv,
    Clus::IPCTransformSet::pointer pcts,
    Weighted::Graph& graph)
{
    // This used to be the body of Cluster::Connect_graph(dv,pcts,use_ctpc).
    const bool use_ctpc=true;
    const auto* grouping = cluster.grouping();

    // now form the connected components
    std::vector<int> component(num_vertices(graph));
    const size_t num = connected_components(graph, &component[0]);

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

    if (num <= 1) return;

    std::vector<std::shared_ptr<Simple3DPointCloud>> pt_clouds;
    std::vector<std::vector<size_t>> pt_clouds_global_indices; // can use to access wpid ...

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

    /// DEBUGONLY:
    if (0) {
        for (size_t i = 0; i != num; i++) {
            std::cout << *pt_clouds.at(i) << std::endl;
            std::cout << "global indices: ";
            for (size_t j = 0; j != pt_clouds_global_indices.at(i).size(); j++) {
                std::cout << pt_clouds_global_indices.at(i).at(j) << " ";
            }
            std::cout << std::endl;
        }
    }

    // Initiate dist. metrics
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis(
        num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_mst(
        num, std::vector<std::tuple<int, int, double>>(num));

    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir1(
        num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir2(
        num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir_mst(
        num, std::vector<std::tuple<int, int, double>>(num));

    for (size_t j = 0; j != num; j++) {
        for (size_t k = 0; k != num; k++) {
            index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_mst[j][k] = std::make_tuple(-1, -1, 1e9);

            index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir_mst[j][k] = std::make_tuple(-1, -1, 1e9);
        }
    }

    // Calc. dis, dis_dir1, dis_dir2
    // check against the closest distance ...
    // no need to have MST ...
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));

            if ((num < 100 && pt_clouds.at(j)->get_num_points() > 100 && pt_clouds.at(k)->get_num_points() > 100 &&
                 (pt_clouds.at(j)->get_num_points() + pt_clouds.at(k)->get_num_points()) > 400) ||
                (pt_clouds.at(j)->get_num_points() > 500 && pt_clouds.at(k)->get_num_points() > 500)) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));

                geo_point_t dir1 = cluster.vhough_transform(p1, 30 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                geo_point_t dir2 = cluster.vhough_transform(p2, 30 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                dir1 = dir1 * -1;
                dir2 = dir2 * -1;

                std::pair<int, double> result1 = pt_clouds.at(k)->get_closest_point_along_vec(
                    p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);

                if (result1.first >= 0) {
                    index_index_dis_dir1[j][k] =
                        std::make_tuple(std::get<0>(index_index_dis[j][k]), result1.first, result1.second);
                }

                std::pair<int, double> result2 = pt_clouds.at(j)->get_closest_point_along_vec(
                    p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);

                if (result2.first >= 0) {
                    index_index_dis_dir2[j][k] =
                        std::make_tuple(result2.first, std::get<1>(index_index_dis[j][k]), result2.second);
                }
            }

            // Now check the path ...
            {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis[j][k])));

                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis / step_dis + 1;
                int num_bad = 0;
                geo_point_t test_p;
                for (int ii = 0; ii != num_steps; ii++) {
                    test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                               p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                               p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));

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
                            if (!good_point) num_bad++;
                        }
                    }
                }

                if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                    index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
                }
            }

            // Now check the path ...
            if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir1[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k])));

                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir1[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis / step_dis + 1;
                int num_bad = 0;
                geo_point_t test_p;
                for (int ii = 0; ii != num_steps; ii++) {
                    test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                               p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                               p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));
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
                            if (!good_point) num_bad++;
                        }
                    }
                }

                if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                    index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                }
            }

            // Now check the path ...
            if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir2[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k])));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir2[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis / step_dis + 1;
                int num_bad = 0;
                geo_point_t test_p;
                for (int ii = 0; ii != num_steps; ii++) {
                    test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                               p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                               p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));
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
                            if (!good_point) num_bad++;
                        }
                    }
                }

                if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                    index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
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
                    if (!boost::edge(gind1, gind2, graph).second) {
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
                    if (!boost::edge(gind1, gind2, graph).second) {
                        /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                    }
                }
            }

        }  // k
    }  // j
}




using namespace WireCell::Clus::Facade;

void Graphs::connect_graph_ctpc_with_reference(
    const Facade::Cluster& cluster,
    const Facade::Cluster& ref_cluster,
    IDetectorVolumes::pointer dv,
    Clus::IPCTransformSet::pointer pcts,
    Weighted::Graph& graph)
{
    // Enable CTPC functionality (combining ctpc logic with reference filtering)
    const bool use_ctpc = true;
    const auto* grouping = cluster.grouping();
    
    // Drift direction for directional analysis (equivalent to prototype's drift_dir)
    geo_vector_t drift_dir_abs(1, 0, 0);
    
    // Form connected components from existing graph
    std::vector<int> component(num_vertices(graph));
    const size_t num = connected_components(graph, &component[0]);

    // Create ordered components structure (same as baseline)
    std::vector<ComponentInfo> ordered_components;
    ordered_components.reserve(component.size());
    for (size_t i = 0; i < component.size(); ++i) {
        ordered_components.emplace_back(i);
    }

    // Assign vertices to components
    for (size_t i = 0; i < component.size(); ++i) {
        ordered_components[component[i]].add_vertex(i);
    }

    // Sort components by minimum vertex index for deterministic behavior
    std::sort(ordered_components.begin(), ordered_components.end(), 
        [](const ComponentInfo& a, const ComponentInfo& b) {
            return a.min_vertex < b.min_vertex;
        });

    if (num <= 1) return;

    // Initialize point cloud containers for each component
    std::vector<std::shared_ptr<Simple3DPointCloud>> pt_clouds;
    std::vector<std::vector<size_t>> pt_clouds_global_indices;
    
    for (const auto& comp : ordered_components) {
        auto pt_cloud = std::make_shared<Simple3DPointCloud>();
        pt_clouds.push_back(pt_cloud);
        pt_clouds_global_indices.push_back(std::vector<size_t>());
    }
    
    const auto& points = cluster.points();
    std::set<size_t> excluded_points;  // Track excluded points (prototype's excluded_points)
    
    // Check if reference cluster is valid and not empty
    bool use_reference_filtering = (ref_cluster.is_valid() && ref_cluster.npoints() > 0);
    
    // REFERENCE FILTERING PHASE - equivalent to prototype's filtering logic
    for (size_t i = 0; i < component.size(); ++i) {
        bool should_exclude = false;
        
        // Phase 1: Check point quality (equivalent to prototype's mcell->IsPointGood)
        if (!is_point_good(cluster, i, 2)) {
            should_exclude = true;
        } 
        // Phase 2: Reference cluster distance filtering (only if ref_cluster is not empty)
        else if (use_reference_filtering) {
            const auto& ref_kd = ref_cluster.kd3d();  
            double temp_min_dis = 0;
            geo_point_t temp_p(points[0][i], points[1][i], points[2][i]);
            std::vector<double> query_point = {temp_p.x(), temp_p.y(), temp_p.z()};
            auto knn_result = ref_kd.knn(1, query_point);
            
            if (!knn_result.empty()) {
                temp_min_dis = std::sqrt(knn_result[0].second);  // knn returns squared distance
            }
            
            // Key filtering criterion from prototype: >= 1.0 cm means exclude
            if (temp_min_dis >= 1.0 * units::cm) {
                should_exclude = true;
            }
        }
        // If ref_cluster is empty, only use point quality check
        
        if (should_exclude) {
            excluded_points.insert(i);
        } else {
            // Add to appropriate component cloud
            size_t comp_idx = component[i];
            pt_clouds.at(comp_idx)->add({points[0][i], points[1][i], points[2][i]});
            pt_clouds_global_indices.at(comp_idx).push_back(i);
        }
    }
    
    // Store excluded points in cluster cache (matches prototype's excluded_points)
    const_cast<Cluster&>(cluster).set_excluded_points(excluded_points);

    /// DEBUGONLY:
    if (0) {
        for (size_t i = 0; i != num; i++) {
            std::cout << *pt_clouds.at(i) << std::endl;
            std::cout << "global indices: ";
            for (size_t j = 0; j != pt_clouds_global_indices.at(i).size(); j++) {
                std::cout << pt_clouds_global_indices.at(i).at(j) << " ";
            }
            std::cout << std::endl;
        }
    }

    // Initialize distance metric containers (same structure as baseline)
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis(
        num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_mst(
        num, std::vector<std::tuple<int, int, double>>(num));

    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir1(
        num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir2(
        num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir_mst(
        num, std::vector<std::tuple<int, int, double>>(num));

    // Initialize all distances to invalid/infinite
    for (size_t j = 0; j != num; j++) {
        for (size_t k = 0; k != num; k++) {
            index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_mst[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir_mst[j][k] = std::make_tuple(-1, -1, 1e9);
        }
    }

    // DISTANCE CALCULATION AND CTPC PATH VALIDATION
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            // Find closest points between components
            index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));

            // Enhanced directional analysis for large components (from prototype logic)
            if ((num < 100 && pt_clouds.at(j)->get_num_points() > 100 && pt_clouds.at(k)->get_num_points() > 100 &&
                 (pt_clouds.at(j)->get_num_points() + pt_clouds.at(k)->get_num_points()) > 400) ||
                (pt_clouds.at(j)->get_num_points() > 500 && pt_clouds.at(k)->get_num_points() > 500)) {
                
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));
                
                // Use cluster's vhough_transform method for directional analysis
                geo_vector_t dir1 = cluster.vhough_transform(p1, 30 * units::cm, 
                    Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                geo_vector_t dir2 = cluster.vhough_transform(p2, 30 * units::cm, 
                    Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                dir1 = dir1 * -1;
                dir2 = dir2 * -1;

                // Directional search from p1 towards p2
                std::pair<int, double> result1 = pt_clouds.at(k)->get_closest_point_along_vec(
                    p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                
                // Enhanced drift direction analysis (from prototype)
                if (result1.first < 0) {
                    double angle_deg = dir1.angle(drift_dir_abs) * 180.0 / M_PI;
                    if (std::abs(angle_deg - 90.0) < 10.0) {
                        // Direction nearly perpendicular to drift - try longer hough transform
                        if (std::abs(angle_deg - 90.0) < 5.0) {
                            dir1 = cluster.vhough_transform(p1, 80 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                        } else if (std::abs(angle_deg - 90.0) < 10.0) {
                            dir1 = cluster.vhough_transform(p1, 50 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                        }
                        dir1 = dir1 * -1;
                        result1 = pt_clouds.at(k)->get_closest_point_along_vec(
                            p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                    }
                }

                if (result1.first >= 0) {
                    index_index_dis_dir1[j][k] = std::make_tuple(
                        std::get<0>(index_index_dis[j][k]), result1.first, result1.second);
                }

                // Directional search from p2 towards p1
                std::pair<int, double> result2 = pt_clouds.at(j)->get_closest_point_along_vec(
                    p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);

                // Enhanced drift direction analysis for dir2
                if (result2.first < 0) {
                    double angle_deg2 = dir2.angle(drift_dir_abs) * 180.0 / M_PI;
                    if (std::abs(angle_deg2 - 90.0) < 10.0) {
                        if (std::abs(angle_deg2 - 90.0) < 5.0) {
                            dir2 = cluster.vhough_transform(p2, 80 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                        } else if (std::abs(angle_deg2 - 90.0) < 10.0) {
                            dir2 = cluster.vhough_transform(p2, 50 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                        }
                        dir2 = dir2 * -1;
                        result2 = pt_clouds.at(j)->get_closest_point_along_vec(
                            p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                    }
                }

                if (result2.first >= 0) {
                    index_index_dis_dir2[j][k] = std::make_tuple(
                        result2.first, std::get<1>(index_index_dis[j][k]), result2.second);
                }
            }

            // CTPC PATH VALIDATION - Check basic distance path
            {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis[j][k])));

                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis / step_dis + 1;
                int num_bad = 0;
                geo_point_t test_p;
                
                for (int ii = 0; ii != num_steps; ii++) {
                    test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                               p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                               p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));

                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa() != -1) {
                            geo_point_t test_p_raw = test_p;
                            if (cluster.get_default_scope().hash() != cluster.get_raw_scope().hash()) {
                                const auto transform = pcts->pc_transform(cluster.get_scope_transform());
                                double cluster_t0 = cluster.get_flash().time();
                                test_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            if (!good_point) num_bad++;
                        }
                    }
                }

                if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                    index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
                }
            }

            // CTPC PATH VALIDATION - Check directional path 1
            if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir1[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k])));

                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir1[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis / step_dis + 1;
                int num_bad = 0;
                geo_point_t test_p;
                
                for (int ii = 0; ii != num_steps; ii++) {
                    test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                               p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                               p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));
                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa() != -1) {
                            geo_point_t test_p_raw = test_p;
                            if (cluster.get_default_scope().hash() != cluster.get_raw_scope().hash()) {
                                const auto transform = pcts->pc_transform(cluster.get_scope_transform());
                                double cluster_t0 = cluster.get_flash().time();
                                test_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            if (!good_point) num_bad++;
                        }
                    }
                }

                if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                    index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                }
            }

            // CTPC PATH VALIDATION - Check directional path 2
            if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir2[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k])));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir2[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis / step_dis + 1;
                int num_bad = 0;
                geo_point_t test_p;
                
                for (int ii = 0; ii != num_steps; ii++) {
                    test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                               p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                               p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));
                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa() != -1) {
                            geo_point_t test_p_raw = test_p;
                            if (cluster.get_default_scope().hash() != cluster.get_raw_scope().hash()) {
                                const auto transform = pcts->pc_transform(cluster.get_scope_transform());
                                double cluster_t0 = cluster.get_flash().time();
                                test_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            if (!good_point) num_bad++;
                        }
                    }
                }

                if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                    index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                }
            }
        }
    }

    // Build MST for basic distances
    {
        Weighted::Graph temp_graph(num);
        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j, index2 = k;
                if (std::get<0>(index_index_dis[j][k]) >= 0) {
                    if (!boost::edge(index1, index2, temp_graph).second) {
                        add_edge(index1, index2, std::get<2>(index_index_dis[j][k]), temp_graph);
                    }
                }
            }
        }
        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_mst);
    }

    // Build MST for directional distances
    {
        Weighted::Graph temp_graph(num);
        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j, index2 = k;
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0 || 
                    std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    if (!boost::edge(index1, index2, temp_graph).second) {
                        add_edge(index1, index2,
                            std::min(std::get<2>(index_index_dis_dir1[j][k]),
                                   std::get<2>(index_index_dis_dir2[j][k])), temp_graph);
                    }
                }
            }
        }
        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_dir_mst);
    }

    // Final graph construction phase
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            // Add short distance connections directly to MST
            if (std::get<2>(index_index_dis[j][k]) < 3 * units::cm) {
                index_index_dis_mst[j][k] = index_index_dis[j][k];
            }

            // Add MST basic distance edges to graph
            if (std::get<0>(index_index_dis_mst[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_mst[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_mst[j][k]));

                float dis;
                if (std::get<2>(index_index_dis_mst[j][k]) > 5 * units::cm) {
                    dis = std::get<2>(index_index_dis_mst[j][k]);
                } else {
                    dis = std::get<2>(index_index_dis_mst[j][k]);
                }
                if (!boost::edge(gind1, gind2, graph).second) {
                    /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                }
            }

            // Add MST directional edges to graph (with penalty for longer distances)
            if (std::get<0>(index_index_dis_dir_mst[j][k]) >= 0) {
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k]));

                    float dis;
                    if (std::get<2>(index_index_dis_dir1[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir1[j][k]) * 1.1;  // Matches ctpc baseline penalty
                    } else {
                        dis = std::get<2>(index_index_dis_dir1[j][k]);
                    }
                    if (!boost::edge(gind1, gind2, graph).second) {
                        /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                    }
                }
                
                if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k]));

                    float dis;
                    if (std::get<2>(index_index_dis_dir2[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir2[j][k]) * 1.1;  // Matches ctpc baseline penalty
                    } else {
                        dis = std::get<2>(index_index_dis_dir2[j][k]);
                    }
                    if (!boost::edge(gind1, gind2, graph).second) {
                        /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                    }
                }
            }
        }
    }
}