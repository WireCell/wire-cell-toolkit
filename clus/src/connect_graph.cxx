#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Grouping.h"

#include "connect_graphs.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

void Graphs::connect_graph(const Cluster& cluster, Weighted::Graph& graph)
{
    // This used to be the body of Cluster::Connect_graph().

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
    std::vector<std::vector<size_t>> pt_clouds_global_indices;
    // use this to link the global index to the local index
    // Create point clouds using ordered components
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

    Weighted::Graph temp_graph(num);

    for (size_t j=0;j!=num;j++){
      for (size_t k=j+1;k!=num;k++){
            index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));

            int index1 = j;
            int index2 = k;

            if (!boost::edge(index1, index2, temp_graph).second) 
            /*auto edge =*/ add_edge(index1,index2, std::get<2>(index_index_dis[j][k]), temp_graph);
      }
    }

    process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_mst);

    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            if (std::get<2>(index_index_dis[j][k])<3*units::cm){
                index_index_dis_mst[j][k] = index_index_dis[j][k];
            }

            if (num < 100)
            if (pt_clouds.at(j)->get_num_points()>100 && pt_clouds.at(k)->get_num_points()>100 &&
                (pt_clouds.at(j)->get_num_points()+pt_clouds.at(k)->get_num_points()) > 400){
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
        }
    }

    // MST for the directionality ...
    {
        Weighted::Graph temp_graph(num);

        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j;
                int index2 = k;
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0 || std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    
                    if (!boost::edge(index1, index2, temp_graph).second) 
                    add_edge(
                        index1, index2,
                        std::min(std::get<2>(index_index_dis_dir1[j][k]), std::get<2>(index_index_dis_dir2[j][k])),
                        temp_graph);
                }
            }
        }

        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_dir_mst);
    
    }

    // now complete graph according to the direction
    // according to direction ...
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
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
                        dis = std::get<2>(index_index_dis_dir1[j][k]) * 1.2;
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
                        dis = std::get<2>(index_index_dis_dir2[j][k]) * 1.2;
                    }
                    else {
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


using namespace WireCell::Clus::Facade;


void Graphs::connect_graph_with_reference(
    const Facade::Cluster& cluster,
    const Facade::Cluster& ref_cluster,
    Weighted::Graph& graph)
{
    // Drift direction (used in prototype for angle checks)
    geo_vector_t drift_dir_abs(1, 0, 0);
    
    // now form the connected components
    std::vector<int> component(num_vertices(graph));
    const size_t num = connected_components(graph, &component[0]);

    // Create ordered components (same as baseline)
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
    std::vector<std::vector<size_t>> pt_clouds_global_indices;
    
    // Initialize pt_clouds for each component (same as baseline)
    for (size_t comp_idx = 0; comp_idx < ordered_components.size(); ++comp_idx) {
        auto pt_cloud = std::make_shared<Simple3DPointCloud>();
        pt_clouds.push_back(pt_cloud);
        pt_clouds_global_indices.push_back(std::vector<size_t>());
    }
    
    const auto& points = cluster.points();
    std::set<size_t> excluded_points;  // Track excluded points
    
    // Check if reference cluster is empty
    bool use_reference_filtering = (ref_cluster.is_valid() && ref_cluster.npoints() > 0);
    
    // Process each point with reference filtering (matches prototype exactly)
    for (size_t i = 0; i < component.size(); ++i) {
        bool should_exclude = false;
        
        // Check if point is good (equivalent to prototype's mcell->IsPointGood)
        if (!is_point_good(cluster, i, 2)) {
            should_exclude = true;
        } else if (use_reference_filtering) {
            // Only check distance to reference cluster if it's not empty
            const auto& ref_kd = ref_cluster.kd3d();  // Use reference cluster's KD-tree
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
        // If ref_cluster is empty, we skip the distance check and only use the point quality check
        
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
    // Note: When ref_cluster is empty, excluded_points will only contain points that fail the quality check
    const_cast<Cluster&>(cluster).set_excluded_points(excluded_points);


    // Initiate dist. metrics (same as baseline)
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

    // Distance calculation (same as baseline)
    Weighted::Graph temp_graph(num);

    for (size_t j=0;j!=num;j++){
      for (size_t k=j+1;k!=num;k++){
            index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));

            int index1 = j;
            int index2 = k;

            if (!boost::edge(index1, index2, temp_graph).second) 
            /*auto edge =*/ add_edge(index1,index2, std::get<2>(index_index_dis[j][k]), temp_graph);
      }
    }

    process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_mst);

    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            if (std::get<2>(index_index_dis[j][k])<3*units::cm){
                index_index_dis_mst[j][k] = index_index_dis[j][k];
            }

            if (num < 100)
            if (pt_clouds.at(j)->get_num_points()>100 && pt_clouds.at(k)->get_num_points()>100 &&
                (pt_clouds.at(j)->get_num_points()+pt_clouds.at(k)->get_num_points()) > 400){
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));
            
                // Use cluster's vhough_transform method with drift direction awareness
                geo_vector_t dir1 = cluster.vhough_transform(p1, 30 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                geo_vector_t dir2 = cluster.vhough_transform(p2, 30 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                dir1 = dir1 * -1;
                dir2 = dir2 * -1;

                std::pair<int, double> result1 = pt_clouds.at(k)->get_closest_point_along_vec(
                    p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                
                    // If no result and direction is nearly perpendicular to drift, try longer hough transform as in prototype
                    double angle_deg = dir1.angle(drift_dir_abs) * 180.0 / M_PI;
                    if (result1.first < 0 && std::abs(angle_deg - 90.0) < 10.0) {
                        if (std::abs(angle_deg - 90.0) < 5.0)
                            dir1 = cluster.vhough_transform(p1, 80 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                        else if (std::abs(angle_deg - 90.0) < 10.0)
                            dir1 = cluster.vhough_transform(p1, 50 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                        dir1 = dir1 * -1;
                        result1 = pt_clouds.at(k)->get_closest_point_along_vec(
                            p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                    }

                if (result1.first >= 0) {
                    index_index_dis_dir1[j][k] =
                        std::make_tuple(std::get<0>(index_index_dis[j][k]), result1.first, result1.second);
                }

                std::pair<int, double> result2 = pt_clouds.at(j)->get_closest_point_along_vec(
                    p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);

                // Additional drift direction check (from prototype, though isochronous search was commented out)
                // If no result and direction is nearly perpendicular to drift, try longer hough transform as in prototype
                double angle_deg2 = dir2.angle(drift_dir_abs) * 180.0 / M_PI;
                if (result2.first < 0 && std::abs(angle_deg2 - 90.0) < 10.0) {
                    if (std::abs(angle_deg2 - 90.0) < 5.0)
                        dir2 = cluster.vhough_transform(p2, 80 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                    else if (std::abs(angle_deg2 - 90.0) < 10.0)
                        dir2 = cluster.vhough_transform(p2, 50 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                    dir2 = dir2 * -1;
                    result2 = pt_clouds.at(j)->get_closest_point_along_vec(
                        p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                }

              

                if (result2.first >= 0) {
                    index_index_dis_dir2[j][k] =
                        std::make_tuple(result2.first, std::get<1>(index_index_dis[j][k]), result2.second);
                }
            }
        }
    }

    // MST for the directionality ... (same as baseline)
    {
        Weighted::Graph temp_graph(num);

        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j;
                int index2 = k;
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0 || std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    
                    if (!boost::edge(index1, index2, temp_graph).second) 
                    add_edge(
                        index1, index2,
                        std::min(std::get<2>(index_index_dis_dir1[j][k]), std::get<2>(index_index_dis_dir2[j][k])),
                        temp_graph);
                }
            }
        }

        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_dir_mst);
    
    }

    // now complete graph according to the direction (same as baseline)
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
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
                        dis = std::get<2>(index_index_dis_dir1[j][k]) * 1.2;
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
                        dis = std::get<2>(index_index_dis_dir2[j][k]) * 1.2;
                    }
                    else {
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

// Helper function equivalent to prototype's mcell->IsPointGood
bool Graphs::is_point_good(const Cluster& cluster, size_t point_index, int ncut) {
    double charge_u = cluster.charge_value(point_index, 0);
    double charge_v = cluster.charge_value(point_index, 1); 
    double charge_w = cluster.charge_value(point_index, 2);
    
    int ncount = 0;
    if (charge_u > 10) ncount++;
    if (charge_v > 10) ncount++;
    if (charge_w > 10) ncount++;
    
    return ncount >= ncut;
}

std::vector<bool> check_direction(const Facade::Cluster& cluster, Facade::geo_vector_t& v1, int apa, int face, double angle_cut_1, double angle_cut_2){
    // Get grouping to access wire geometry
    auto grouping = cluster.grouping();
    if (!grouping) {
        // Return all false if no grouping available
        return std::vector<bool>(4, false);
    }
    
    // Get wire angles from grouping for this APA and face
    auto [angle_u, angle_v, angle_w] = grouping->wire_angles(apa, face);
    
    // Get drift direction from grouping
    int drift_dirx = grouping->get_drift_dir().at(apa).at(face);
    Facade::geo_vector_t drift_dir_abs(fabs(drift_dirx), 0, 0);
    
    // Construct wire direction vectors
    // U wire: angle_u from Y axis in YZ plane
    Facade::geo_vector_t U_dir(0, std::cos(angle_u), std::sin(angle_u));
    // V wire: angle_v from Y axis in YZ plane  
    Facade::geo_vector_t V_dir(0, std::cos(angle_v), std::sin(angle_v));
    // W wire: angle_w from Y axis in YZ plane
    Facade::geo_vector_t W_dir(0, std::cos(angle_w), std::sin(angle_w));
    
    // Project v1 onto YZ plane
    Facade::geo_vector_t tempV1(0, v1.y(), v1.z());
    Facade::geo_vector_t tempV5;
    
    // Prolonged U - project onto plane perpendicular to U wire direction
    double angle1 = tempV1.angle(U_dir);
    tempV5 = Facade::geo_vector_t(
        std::fabs(v1.x()),
        std::sqrt(v1.y()*v1.y() + v1.z()*v1.z()) * std::sin(angle1),
        0
    );
    angle1 = tempV5.angle(drift_dir_abs);
    
    // Prolonged V - project onto plane perpendicular to V wire direction
    double angle2 = tempV1.angle(V_dir);
    tempV5 = Facade::geo_vector_t(
        std::fabs(v1.x()),
        std::sqrt(v1.y()*v1.y() + v1.z()*v1.z()) * std::sin(angle2),
        0
    );
    angle2 = tempV5.angle(drift_dir_abs);
    
    // Prolonged W - project onto plane perpendicular to W wire direction
    double angle3 = tempV1.angle(W_dir);
    tempV5 = Facade::geo_vector_t(
        std::fabs(v1.x()),
        std::sqrt(v1.y()*v1.y() + v1.z()*v1.z()) * std::sin(angle3),
        0
    );
    angle3 = tempV5.angle(drift_dir_abs);
    
    // Parallel - angle with respect to drift direction
    double angle4 = v1.angle(drift_dir);
    
    std::vector<bool> results(4, false);
    
    // Check if prolonged along U wire (< 12.5 degrees)
    if (angle1 < angle_cut_1 / 180.0 * M_PI) results.at(0) = true;
    // Check if prolonged along V wire (< 12.5 degrees)
    if (angle2 < angle_cut_1 / 180.0 * M_PI) results.at(1) = true;
    // Check if prolonged along W wire (< 12.5 degrees)
    if (angle3 < angle_cut_1 / 180.0 * M_PI) results.at(2) = true;
    // Check if perpendicular to drift (within 10 degrees of 90 degrees)
    if (std::fabs(angle4 - M_PI/2.0) < angle_cut_2 / 180.0 * M_PI) results.at(3) = true;
    
    return results;
}
