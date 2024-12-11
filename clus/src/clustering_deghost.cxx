#include <WireCellClus/ClusteringFuncs.h>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;

// The original developers do not care.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

// #define __DEBUG__
#ifdef __DEBUG__
#define LogDebug(x) std::cout << "[yuhw]: " << __LINE__ << " : " << x << std::endl
#else
#define LogDebug(x)
#endif


void WireCell::PointCloud::Facade::clustering_deghost(Grouping& live_grouping,
                                  const bool use_ctpc, double length_cut)
{
    std::map<int, std::pair<double, double>>& dead_u_index = live_grouping.get_dead_winds(0, 0);
    std::map<int, std::pair<double, double>>& dead_v_index = live_grouping.get_dead_winds(0, 1);
    std::map<int, std::pair<double, double>>& dead_w_index = live_grouping.get_dead_winds(0, 2);

    std::vector<Cluster *> live_clusters = live_grouping.children();  // copy
    // sort the clusters by length using a lambda function
    std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster *cluster1, const Cluster *cluster2) {
        return cluster1->get_length() > cluster2->get_length();
    });

    const auto &tp = live_grouping.get_params();
    // this is for 4 time slices
    // double time_slice_width = tp.nticks_live_slice * tp.tick_drift;

    // Create two point clouds ...
    // One for the points ... point --> index --> cluster (vector) ...
    // The other for the skeleton of each track ...  point --> index --> cluster (vector)
    // Both cloud needs to be dynamic, keep adding things into it as we improve the knowledge
    auto global_point_cloud = std::make_shared<DynamicPointCloud>(tp.angle_u, tp.angle_v, tp.angle_w);
    auto global_skeleton_cloud = std::make_shared<DynamicPointCloud>(tp.angle_u, tp.angle_v, tp.angle_w);

    std::vector<Cluster *> to_be_removed_clusters;
    // std::set<std::pair<Cluster *, Cluster *>> to_be_merged_pairs;
    cluster_connectivity_graph_t  g;
    std::unordered_map<int, int> ilive2desc;  // added live index to graph descriptor
    std::map<const Cluster*, int> map_cluster_index;
    for (const Cluster* live : live_grouping.children()) {
        size_t ilive = map_cluster_index.size();
        map_cluster_index[live] = ilive;
        ilive2desc[ilive] = boost::add_vertex(ilive, g);
    }

    for (size_t i = 0; i != live_clusters.size(); i++) {
        if (i == 0) {
            // fill anyway ...
            // live_clusters.at(i)->Create_point_cloud();
            global_point_cloud->add_points(live_clusters.at(i), 0);
            if (live_clusters.at(i)->get_length() >
                30 * units::cm) {  // should be the default for most of them ...
                live_clusters.at(i)->construct_skeleton(use_ctpc);
                global_skeleton_cloud->add_points(live_clusters.at(i), 1);
            }
            else {
                global_skeleton_cloud->add_points(live_clusters.at(i), 0);
            }
        }
        else {
            // start the process to add things in and perform deghosting ...
            Cluster *cluster = live_clusters.at(i);
            const auto& winds = cluster->wire_indices();

            if (length_cut == 0 || live_clusters.at(i)->get_length() < length_cut) {
                // cluster->Create_point_cloud();
                // WCP::WCPointCloud<double> &cloud = cluster->get_point_cloud()->get_cloud();
                // int num_total_points = cloud.pts.size();  // total number of points
                const size_t num_total_points = cluster->npoints();  // total number of points
                size_t num_dead[3] = {0, 0, 0};              // dead wires in each view
                size_t num_unique[3] = {0, 0, 0};            // points that are unique (not agree with any other clusters)
                std::map<const Cluster *, int> map_cluster_num[3];

                double dis_cut = 1.2 * units::cm;

                for (size_t j = 0; j != num_total_points; j++) {
                    // geo_point_t test_point(cluster->point3d(j).x(), cloud.pts.at(j).y, cloud.pts.at(j).z);
                    geo_point_t test_point = cluster->point3d(j);

                    bool flag_dead = false;

                    #ifdef __DEBUG__
                    if (num_total_points == 134) {
                        std::cout << "point: " << j << " " << test_point << " " << winds[0][j] << " " << winds[1][j] << " " << winds[2][j] << std::endl;
                    }
                    #endif

                    

                    if (dead_u_index.find(winds[0][j]) != dead_u_index.end()) {
                        if (cluster->point3d(j).x() >= dead_u_index[winds[0][j]].first &&
                            cluster->point3d(j).x() <= dead_u_index[winds[0][j]].second) {
                            flag_dead = true;
                        }
                    }

                    

                    if (!flag_dead) {
                        std::tuple<double, const Cluster *, size_t> results =
                            global_point_cloud->get_closest_2d_point_info(test_point, 0);

                        // if (cluster->nchildren()==801 && j==0)  std::cout  << j << " AU " << test_point << " " << std::get<0>(results) << " " << std::get<1>(results)->get_length()/units::cm << std::endl;

                        if (std::get<0>(results) <= dis_cut / 3.) {
                            if (map_cluster_num[0].find(std::get<1>(results)) == map_cluster_num[0].end()) {
                                map_cluster_num[0][std::get<1>(results)] = 1;
                            }
                            else {
                                map_cluster_num[0][std::get<1>(results)]++;
                            }
                        }
                        else {
                            results = global_skeleton_cloud->get_closest_2d_point_info(test_point, 0);

                            // if (cluster->nchildren()==801 && j==0)  std::cout  << j << " BU " << test_point << " " << std::get<0>(results) << " " << std::get<1>(results)->get_length()/units::cm << std::endl;

                            if (std::get<0>(results) <= dis_cut * 2.0) {
                                if (map_cluster_num[0].find(std::get<1>(results)) == map_cluster_num[0].end()) {
                                    map_cluster_num[0][std::get<1>(results)] = 1;
                                }
                                else {
                                    map_cluster_num[0][std::get<1>(results)]++;
                                }
                            }
                            else {
                                num_unique[0]++;
                            }
                        }
                    }
                    else {
                        num_dead[0]++;
                    }

                    flag_dead = false;

                    if (dead_v_index.find(winds[1][j]) != dead_v_index.end()) {
                        #ifdef __DEBUG__
                        if (num_total_points == 134) {
                            std::cout << "dead_v_index: " << winds[1][j] << " " << dead_v_index[winds[1][j]].first << " " << dead_v_index[winds[1][j]].second << std::endl;
                        }
                        #endif
                        if (cluster->point3d(j).x() >= dead_v_index[winds[1][j]].first &&
                            cluster->point3d(j).x() <= dead_v_index[winds[1][j]].second) {
                            flag_dead = true;
                        }
                    }

                    if (!flag_dead) {
                        #ifdef __DEBUG__
                        if (num_total_points == 134) {
                            for (size_t i = 0; i != global_point_cloud->get_num_points(); i++) {
                                const auto p3d = global_point_cloud->point3d(i);
                                const auto p2d0 = global_point_cloud->point2d(0, i);
                                const auto p2d1 = global_point_cloud->point2d(1, i);
                                LogDebug("global_point_cloud: " << i << " 3d " << p3d << " 2dp0 " << p2d0[0] << " " << p2d0[1] << " 2dp1 " << p2d1[0] << " " << p2d1[1]);
                            }
                        }
                        #endif
                        std::tuple<double, const Cluster *, size_t> results =
                            global_point_cloud->get_closest_2d_point_info(test_point, 1);

                        // if (cluster->nchildren()==801 && j==0)  std::cout  << j << " AV " << test_point << " " << std::get<0>(results) << " " << std::get<1>(results)->get_length()/units::cm << std::endl;

                        #ifdef __DEBUG__
                        if (num_total_points == 134) {
                            const auto c = std::get<1>(results);
                            const auto p = global_point_cloud->point3d(std::get<2>(results));
                            LogDebug("results: cluster " << c->npoints() << " point " << p << " dist " << std::get<0>(results)/units::mm << " dist_cut: " << dis_cut);
                        }
                        #endif
                        if (std::get<0>(results) <= dis_cut / 3.) {
                            if (map_cluster_num[1].find(std::get<1>(results)) == map_cluster_num[1].end()) {
                                map_cluster_num[1][std::get<1>(results)] = 1;
                            }
                            else {
                                map_cluster_num[1][std::get<1>(results)]++;
                            }
                        }
                        else {
                            results = global_skeleton_cloud->get_closest_2d_point_info(test_point, 1);

                            // if (cluster->nchildren()==801 && j==0)  std::cout  << j << " BV " << test_point << " " << std::get<0>(results) << " " << std::get<1>(results)->get_length()/units::cm << std::endl;

                            if (std::get<0>(results) <= dis_cut * 2.0) {
                                if (map_cluster_num[1].find(std::get<1>(results)) == map_cluster_num[1].end()) {
                                    map_cluster_num[1][std::get<1>(results)] = 1;
                                }
                                else {
                                    map_cluster_num[1][std::get<1>(results)]++;
                                }
                            }
                            else {
                                num_unique[1]++;
                            }
                        }
                    }
                    else {
                        num_dead[1]++;
                    }

                    flag_dead = false;

                    if (dead_w_index.find(winds[2][j]) != dead_w_index.end()) {
                        if (cluster->point3d(j).x() >= dead_w_index[winds[2][j]].first &&
                            cluster->point3d(j).x() <= dead_w_index[winds[2][j]].second) {
                            flag_dead = true;
                        }
                    }

                    if (!flag_dead) {
                        std::tuple<double, const Cluster *, size_t> results =
                            global_point_cloud->get_closest_2d_point_info(test_point, 2);

                        // if (cluster->nchildren()==801 && j==0)  std::cout  << j << " AW " << test_point << " " << std::get<0>(results) << " " << std::get<1>(results)->get_length()/units::cm << std::endl;

                        if (std::get<0>(results) <= dis_cut / 3.) {
                            if (map_cluster_num[2].find(std::get<1>(results)) == map_cluster_num[2].end()) {
                                map_cluster_num[2][std::get<1>(results)] = 1;
                            }
                            else {
                                map_cluster_num[2][std::get<1>(results)]++;
                            }
                        }
                        else {
                            results = global_skeleton_cloud->get_closest_2d_point_info(test_point, 2);

                            // if (cluster->nchildren()==801 && j==0)  std::cout  << j << " BW " << test_point <<  " " <<std::get<0>(results) << " " << std::get<1>(results)->get_length()/units::cm << std::endl;

                            if (std::get<0>(results) <= dis_cut * 2.0) {
                                if (map_cluster_num[2].find(std::get<1>(results)) == map_cluster_num[2].end()) {
                                    map_cluster_num[2][std::get<1>(results)] = 1;
                                }
                                else {
                                    map_cluster_num[2][std::get<1>(results)]++;
                                }
                            }
                            else {
                                num_unique[2]++;
                            }
                        }
                    }
                    else {
                        num_dead[2]++;
                    }

                    // if (cluster->nchildren()==801 && j==0) std::cout << j << " " << test_point << " " << cluster->point3d(j).x() << " " << flag_dead << " " << num_unique[0] << " " << num_dead[0] << " " << num_unique[1] << " " << num_dead[1] << " " << num_unique[2] << " " << num_dead[2] << std::endl;

                    if ((num_unique[1] + num_unique[0] + num_unique[2]) >= 0.24 * num_total_points &&
                        (num_unique[1] + num_unique[0] + num_unique[2]) > 25)
                        break;
                }
                LogDebug("num_total_points = " << num_total_points);
                LogDebug("num_unique[0] = " << num_unique[0] << ", num_unique[1] = " << num_unique[1] << ", num_unique[2] = " << num_unique[2]);
                LogDebug("num_dead[0] = " << num_dead[0] << ", num_dead[1] = " << num_dead[1] << ", num_dead[2] = " << num_dead[2]);

                bool flag_save = false;

                // if (cluster->nchildren()==801) std::cout << cluster->get_length()/units::cm << " " << num_total_points << " " << num_unique[0] << " " << num_dead[0] << " " << num_unique[1] << " " << num_dead[1] << " " << num_unique[2] << " " << num_dead[2] << " " << std::endl;

                if (((num_unique[0] <= 0.1 * (num_total_points - num_dead[0]) ||
                      num_unique[0] <= 0.1 * num_total_points && num_unique[0] <= 8) &&
                         (num_unique[1] <= 0.1 * (num_total_points - num_dead[1]) ||
                          num_unique[1] <= 0.1 * num_total_points && num_unique[1] <= 8) &&
                         (num_unique[2] <= 0.1 * (num_total_points - num_dead[2]) ||
                          num_unique[2] <= 0.1 * num_total_points && num_unique[2] <= 8) &&
                         ((num_unique[0] + num_unique[1] + num_unique[2]) <=
                              0.05 * (num_total_points - num_dead[0] + num_total_points - num_dead[1] +
                                      num_total_points - num_dead[2]) ||
                          (num_unique[0] + num_unique[1] + num_unique[2]) < 0.15 * num_total_points &&
                              (num_unique[0] + num_unique[1] + num_unique[2]) <= 8) ||
                     num_unique[0] == 0 && num_unique[1] == 0 && num_unique[2] < 0.24 * num_total_points ||
                     num_unique[0] == 0 && num_unique[2] == 0 && num_unique[1] < 0.24 * num_total_points ||
                     num_unique[2] == 0 && num_unique[1] == 0 && num_unique[0] < 0.24 * num_total_points ||
                     num_unique[0] == 0 && (num_unique[1] + num_unique[2]) < 0.12 * num_total_points * 2 ||
                     num_unique[1] == 0 && (num_unique[0] + num_unique[2]) < 0.12 * num_total_points * 2 ||
                     num_unique[2] == 0 && (num_unique[1] + num_unique[0]) < 0.12 * num_total_points * 2 ||
                     (num_unique[1] + num_unique[0] + num_unique[2]) < 0.24 * num_total_points &&
                         (num_unique[1] < 0.02 * num_total_points || num_unique[0] < 0.02 * num_total_points ||
                          num_unique[2] < 0.02 * num_total_points)) &&
                    (num_unique[0] + num_unique[1] + num_unique[2]) <= 500) {
                    flag_save = false;
                    LogDebug("pass the first cut " << num_total_points);

                    // now try to compare
                    // find the maximal for each map
                    const Cluster *max_cluster_u = 0, *max_cluster_v = 0, *max_cluster_w = 0;
                    int max_value_u = 0, max_value_v = 0, max_value_w = 0;
                    for (auto it = map_cluster_num[0].begin(); it != map_cluster_num[0].end(); it++) {
                        if (it->second > max_value_u) {
                            max_value_u = it->second;
                            max_cluster_u = it->first;
                        }
                    }
                    for (auto it = map_cluster_num[1].begin(); it != map_cluster_num[1].end(); it++) {
                        if (it->second > max_value_v) {
                            max_value_v = it->second;
                            max_cluster_v = it->first;
                        }
                    }
                    for (auto it = map_cluster_num[2].begin(); it != map_cluster_num[2].end(); it++) {
                        if (it->second > max_value_w) {
                            max_value_w = it->second;
                            max_cluster_w = it->first;
                        }
                    }
                    bool flag_remove = true;
                    LogDebug("max_value_u: " << max_value_u << ", max_value_v: " << max_value_v << ", max_value_w: " << max_value_w);

                    if (max_cluster_u == max_cluster_v && max_value_u > 0.8 * (num_total_points - num_dead[0]) &&
                        max_value_v > 0.8 * (num_total_points - num_dead[1])) {
                        if (map_cluster_num[2].find(max_cluster_u) != map_cluster_num[2].end()) {
                            if (map_cluster_num[2][max_cluster_u] > 0.65 * (num_total_points - num_dead[2])) {
                                // ToyPointCloud *cloud1 = cluster->get_point_cloud();
                                // ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
                                // std::tuple<int, int, double> temp_results = cloud1->get_closest_points(cloud2);
                                std::tuple<int, int, double> temp_results = cluster->get_closest_points(*max_cluster_u);
                                if (std::get<2>(temp_results) < 20 * units::cm) {
                                    // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster_u));
                                    boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                    ilive2desc[map_cluster_index[max_cluster_u]], g);
                                    flag_remove = false;
                                }
                            }
                        }
                        else {
                            if (num_total_points == num_dead[2] && max_cluster_u != 0) {
                                // ToyPointCloud *cloud1 = cluster->get_point_cloud();
                                // ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
                                // std::tuple<int, int, double> temp_results = cloud1->get_closest_points(cloud2);
                                std::tuple<int, int, double> temp_results = cluster->get_closest_points(*max_cluster_u);
                                if (std::get<2>(temp_results) < 20 * units::cm) {
                                    // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster_u));
                                    boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                    ilive2desc[map_cluster_index[max_cluster_u]], g);
                                    flag_remove = false;
                                }
                            }
                        }
                    }
                    else if (max_cluster_u == max_cluster_w && max_value_u > 0.8 * (num_total_points - num_dead[0]) &&
                             max_value_w > 0.8 * (num_total_points - num_dead[2])) {
                        if (map_cluster_num[1].find(max_cluster_u) != map_cluster_num[1].end()) {
                            if (map_cluster_num[1][max_cluster_u] > 0.65 * (num_total_points - num_dead[1])) {
                                // ToyPointCloud *cloud1 = cluster->get_point_cloud();
                                // ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
                                // std::tuple<int, int, double> temp_results = cloud1->get_closest_points(cloud2);
                                std::tuple<int, int, double> temp_results = cluster->get_closest_points(*max_cluster_u);
                                if (std::get<2>(temp_results) < 20 * units::cm) {
                                    // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster_u));
                                    boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                    ilive2desc[map_cluster_index[max_cluster_u]], g);
                                    flag_remove = false;
                                }
                            }
                        }
                        else {
                            if (num_total_points == num_dead[1] && max_cluster_u != 0) {
                                // ToyPointCloud *cloud1 = cluster->get_point_cloud();
                                // ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
                                // std::tuple<int, int, double> temp_results = cloud1->get_closest_points(cloud2);
                                std::tuple<int, int, double> temp_results = cluster->get_closest_points(*max_cluster_u);
                                if (std::get<2>(temp_results) < 20 * units::cm) {
                                    // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster_u));
                                    boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                    ilive2desc[map_cluster_index[max_cluster_u]], g);
                                    flag_remove = false;
                                }
                            }
                        }
                    }
                    else if (max_cluster_w == max_cluster_v && max_value_w > 0.8 * (num_total_points - num_dead[2]) &&
                             max_value_v > 0.8 * (num_total_points - num_dead[1])) {
                        if (map_cluster_num[0].find(max_cluster_w) != map_cluster_num[0].end()) {
                            if (map_cluster_num[0][max_cluster_w] > 0.65 * (num_total_points - num_dead[0])) {
                                // ToyPointCloud *cloud1 = cluster->get_point_cloud();
                                // ToyPointCloud *cloud2 = max_cluster_w->get_point_cloud();
                                // std::tuple<int, int, double> temp_results = cloud1->get_closest_points(cloud2);
                                std::tuple<int, int, double> temp_results = cluster->get_closest_points(*max_cluster_w);
                                if (std::get<2>(temp_results) < 20 * units::cm) {
                                    // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster_w));
                                    boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                    ilive2desc[map_cluster_index[max_cluster_w]], g);
                                    flag_remove = false;
                                }
                            }
                        }
                        else {
                            if (num_total_points == num_dead[0] && max_cluster_w != 0) {
                                // ToyPointCloud *cloud1 = cluster->get_point_cloud();
                                // ToyPointCloud *cloud2 = max_cluster_w->get_point_cloud();
                                // std::tuple<int, int, double> temp_results = cloud1->get_closest_points(cloud2);
                                std::tuple<int, int, double> temp_results = cluster->get_closest_points(*max_cluster_w);
                                if (std::get<2>(temp_results) < 20 * units::cm) {
                                    // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster_w));
                                    boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                    ilive2desc[map_cluster_index[max_cluster_w]], g);
                                    flag_remove = false;
                                }
                            }
                        }
                    }
                    if (flag_remove) {
                        // if (cluster->get_length() > 100*units::cm) std::cout << "Remove cluster 1: " << cluster->nchildren() << " " << cluster->get_length()/units::cm << std::endl;
                        to_be_removed_clusters.push_back(cluster);
                    }
                }
                else {

                    flag_save = true;
                    if ((num_unique[0] + num_unique[1] + num_unique[2]) /
                                (num_total_points - num_dead[0] + num_total_points - num_dead[1] + num_total_points -
                                 num_dead[2] + 1e-9) <
                            0.15 &&
                        cluster->get_length() < 25 * units::cm) {
                        const Cluster *max_cluster_u = 0, *max_cluster_v = 0, *max_cluster_w = 0;
                        int max_value_u = 0, max_value_v = 0, max_value_w = 0;
                        for (auto it = map_cluster_num[0].begin(); it != map_cluster_num[0].end(); it++) {
                            if (it->second > max_value_u) {
                                max_value_u = it->second;
                                max_cluster_u = it->first;
                            }
                        }
                        for (auto it = map_cluster_num[1].begin(); it != map_cluster_num[1].end(); it++) {
                            if (it->second > max_value_v) {
                                max_value_v = it->second;
                                max_cluster_v = it->first;
                            }
                        }
                        for (auto it = map_cluster_num[2].begin(); it != map_cluster_num[2].end(); it++) {
                            if (it->second > max_value_w) {
                                max_value_w = it->second;
                                max_cluster_w = it->first;
                            }
                        }

                        /* std::cout << max_cluster_u << " " << max_value_u/(num_total_points-num_dead[0]+1e-9) << " "
                         */
                        /* 	  << max_cluster_v << " " << max_value_v/(num_total_points-num_dead[1]+1e-9) << " " */
                        /* 	  << max_cluster_w << " " << max_value_w/(num_total_points-num_dead[2]+1e-9) <<
                         * std::endl; */

                        if ((max_cluster_u == max_cluster_v && max_cluster_v == max_cluster_w) ||
                            (max_cluster_u == max_cluster_v && max_cluster_w == 0) ||
                            (max_cluster_w == max_cluster_v && max_cluster_u == 0) ||
                            (max_cluster_u == max_cluster_w && max_cluster_v == 0)) {
                            //  std::cout << cluster->get_cluster_id() << " " << (num_unique[0]+num_unique[1] +
                            //  num_unique[2])/(num_total_points - num_dead[0] + num_total_points - num_dead[1] +
                            //  num_total_points - num_dead[2]+1e-9) << " " <<
                            //  (max_value_u+max_value_v+max_value_w)/(num_total_points  + num_total_points  +
                            //  num_total_points +1e-9) << std::endl;

                            if ((max_value_u + max_value_v + max_value_w) /
                                    (num_total_points + num_total_points + num_total_points + 1e-9) >
                                0.25) {
                                flag_save = false;
                                if (max_cluster_u != 0) {
                                    // ToyPointCloud *cloud1 = cluster->get_point_cloud();
                                    // ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
                                    // std::tuple<int, int, double> temp_results = cloud1->get_closest_points(cloud2);
                                    std::tuple<int, int, double> temp_results = cluster->get_closest_points(*max_cluster_u);
                                    if (std::get<2>(temp_results) < 20 * units::cm) {
                                        // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster_u));
                                        boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                        ilive2desc[map_cluster_index[max_cluster_u]], g);
                                    }
                                    else {
                                        // if (cluster->get_length() > 100*units::cm) std::cout << "Remove cluster 2: " << cluster->nchildren() << " " << cluster->get_length()/units::cm << std::endl;
                                        to_be_removed_clusters.push_back(cluster);
                                    }
                                }
                                else if (max_cluster_v != 0) {
                                    // ToyPointCloud *cloud1 = cluster->get_point_cloud();
                                    // ToyPointCloud *cloud2 = max_cluster_v->get_point_cloud();
                                    // std::tuple<int, int, double> temp_results = cloud1->get_closest_points(cloud2);
                                    std::tuple<int, int, double> temp_results = cluster->get_closest_points(*max_cluster_v);
                                    if (std::get<2>(temp_results) < 20 * units::cm) {
                                        // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster_v));
                                        boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                        ilive2desc[map_cluster_index[max_cluster_v]], g);
                                    }
                                    else {
                                        // if (cluster->get_length() > 100*units::cm) std::cout << "Remove cluster 3: " << cluster->nchildren() << " " << cluster->get_length()/units::cm << std::endl;
                                        to_be_removed_clusters.push_back(cluster);
                                    }
                                }
                                else if (max_cluster_w != 0) {
                                    // ToyPointCloud *cloud1 = cluster->get_point_cloud();
                                    // ToyPointCloud *cloud2 = max_cluster_w->get_point_cloud();
                                    // std::tuple<int, int, double> temp_results = cloud1->get_closest_points(cloud2);
                                    std::tuple<int, int, double> temp_results = cluster->get_closest_points(*max_cluster_w);
                                    if (std::get<2>(temp_results) < 20 * units::cm) {
                                        // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster_w));
                                        boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                        ilive2desc[map_cluster_index[max_cluster_w]], g);
                                    }
                                    else {
                                        // if (cluster->get_length() > 100*units::cm) std::cout << "Remove cluster 4: " << cluster->nchildren() << " " << cluster->get_length()/units::cm << std::endl;
                                        to_be_removed_clusters.push_back(cluster);
                                    }
                                }
                            }
                        }
                        else if (max_cluster_u == max_cluster_v && max_cluster_u != 0) {
                            if ((max_value_u + max_value_v + map_cluster_num[2][max_cluster_u]) /
                                    (num_total_points + num_total_points + num_total_points + 1e-9) >
                                0.25) {
                                flag_save = false;
                                // ToyPointCloud *cloud1 = cluster->get_point_cloud();
                                // ToyPointCloud *cloud2 = max_cluster_u->get_point_cloud();
                                // std::tuple<int, int, double> temp_results = cloud1->get_closest_points(cloud2);
                                std::tuple<int, int, double> temp_results = cluster->get_closest_points(*max_cluster_u);
                                if (std::get<2>(temp_results) < 20 * units::cm) {
                                    // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster_u));
                                    boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                    ilive2desc[map_cluster_index[max_cluster_u]], g);
                                }
                                else {
                                    // if (cluster->get_length() > 100*units::cm) std::cout << "Remove cluster 5: " << cluster->nchildren() << " " << cluster->get_length()/units::cm << std::endl;
                                    to_be_removed_clusters.push_back(cluster);
                                }
                            }
                        }
                        else if (max_cluster_v == max_cluster_w && max_cluster_v != 0) {
                            if ((map_cluster_num[0][max_cluster_v] + max_value_v + max_value_w) /
                                    (num_total_points + num_total_points + num_total_points + 1e-9) >
                                0.25) {
                                flag_save = false;
                                // ToyPointCloud *cloud1 = cluster->get_point_cloud();
                                // ToyPointCloud *cloud2 = max_cluster_v->get_point_cloud();
                                // std::tuple<int, int, double> temp_results = cloud1->get_closest_points(cloud2);
                                std::tuple<int, int, double> temp_results = cluster->get_closest_points(*max_cluster_v);
                                if (std::get<2>(temp_results) < 20 * units::cm) {
                                    // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster_v));
                                    boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                    ilive2desc[map_cluster_index[max_cluster_v]], g);
                                }
                                else {
                                    // if (cluster->get_length() > 100*units::cm) std::cout << "Remove cluster 6: " << cluster->nchildren() << " " << cluster->get_length()/units::cm << std::endl;
                                    to_be_removed_clusters.push_back(cluster);
                                }
                            }
                        }
                        else if (max_cluster_u == max_cluster_w && max_cluster_w != 0) {
                            if ((max_value_u + map_cluster_num[1][max_cluster_w] + max_value_w) /
                                    (num_total_points + num_total_points + num_total_points + 1e-9) >
                                0.25) {
                                flag_save = false;
                                // ToyPointCloud *cloud1 = cluster->get_point_cloud();
                                // ToyPointCloud *cloud2 = max_cluster_w->get_point_cloud();
                                // std::tuple<int, int, double> temp_results = cloud1->get_closest_points(cloud2);
                                std::tuple<int, int, double> temp_results = cluster->get_closest_points(*max_cluster_w);
                                if (std::get<2>(temp_results) < 20 * units::cm) {
                                    // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster_w));
                                    boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                    ilive2desc[map_cluster_index[max_cluster_w]], g);
                                }
                                else {
                                    // if (cluster->get_length() > 100*units::cm) std::cout << "Remove cluster 7: " << cluster->nchildren() << " " << cluster->get_length()/units::cm << std::endl;
                                    to_be_removed_clusters.push_back(cluster);
                                }
                            }
                        }
                    }
                    // two cases, merge clusters or remove clusters
                }

                if (flag_save) {
                    // live_clusters.at(i)->Create_point_cloud();
                    global_point_cloud->add_points(live_clusters.at(i), 0);
                    if (live_clusters.at(i)->get_length() > 30 * units::cm) {
                        live_clusters.at(i)->construct_skeleton(use_ctpc);
                        global_skeleton_cloud->add_points(live_clusters.at(i), 1);
                    }
                }
            }
            else {
                // live_clusters.at(i)->Create_point_cloud();
                global_point_cloud->add_points(live_clusters.at(i), 0);
                if (live_clusters.at(i)->get_length() > 30 * units::cm) {
                    live_clusters.at(i)->construct_skeleton(use_ctpc);
                    global_skeleton_cloud->add_points(live_clusters.at(i), 1);
                }
            }
        }
        LogDebug("Cluster " << i << " " << live_clusters.at(i)->nchildren() << " " << live_clusters.at(i)->npoints());
        LogDebug("global_point_cloud: " << global_point_cloud->get_num_points() << " global_skeleton_cloud: " << global_skeleton_cloud->get_num_points());
    }

    // merge clusters
    cluster_set_t new_clusters;
    merge_clusters(g, live_grouping, new_clusters);

    // remove clusters
    LogDebug("to_be_removed_clusters.size() = " << to_be_removed_clusters.size());
    for (auto live : to_be_removed_clusters) {
        // std::cout << "Remove cluster " << live->nchildren() << " " << live->get_length()/units::cm << std::endl;
        live_grouping.remove_child(*live);
    } 

    // // merge clusters
    // std::vector<std::set<Cluster *>> merge_clusters;
    // for (auto it = to_be_merged_pairs.begin(); it != to_be_merged_pairs.end(); it++) {
    //     Cluster *cluster1 = (*it).first;
    //     Cluster *cluster2 = (*it).second;
    //     //  std::cout << cluster1 << " " << cluster2 << " " << cluster1->get_cluster_id() << " " <<
    //     //  cluster2->get_cluster_id() << std::endl;

    //     bool flag_new = true;
    //     std::vector<std::set<Cluster *>> temp_set;
    //     for (auto it1 = merge_clusters.begin(); it1 != merge_clusters.end(); it1++) {
    //         std::set<Cluster *> &clusters = (*it1);
    //         if (clusters.find(cluster1) != clusters.end() || clusters.find(cluster2) != clusters.end()) {
    //             clusters.insert(cluster1);
    //             clusters.insert(cluster2);
    //             flag_new = false;
    //             temp_set.push_back(clusters);
    //             // break;
    //         }
    //     }
    //     if (flag_new) {
    //         std::set<Cluster *> clusters;
    //         clusters.insert(cluster1);
    //         clusters.insert(cluster2);
    //         merge_clusters.push_back(clusters);
    //     }
    //     if (temp_set.size() > 1) {
    //         // merge them further ...
    //         std::set<Cluster *> clusters;
    //         for (size_t i = 0; i != temp_set.size(); i++) {
    //             for (auto it1 = temp_set.at(i).begin(); it1 != temp_set.at(i).end(); it1++) {
    //                 clusters.insert(*it1);
    //             }
    //             merge_clusters.erase(find(merge_clusters.begin(), merge_clusters.end(), temp_set.at(i)));
    //         }
    //         merge_clusters.push_back(clusters);
    //     }
    // }

    // // merge clusters into new clusters, delete old clusters
    // for (auto it = merge_clusters.begin(); it != merge_clusters.end(); it++) {
    //     std::set<Cluster *> &clusters = (*it);
    //     Cluster *ncluster = new Cluster((*clusters.begin())->get_cluster_id());
    //     live_clusters.push_back(ncluster);
    //     for (auto it1 = clusters.begin(); it1 != clusters.end(); it1++) {
    //         Cluster *ocluster = *(it1);
    //         // std::cout << ocluster->get_cluster_id() << " ";
    //         SMGCSelection &mcells = ocluster->get_mcells();
    //         for (auto it2 = mcells.begin(); it2 != mcells.end(); it2++) {
    //             SlimMergeGeomCell *mcell = (*it2);
    //             // std::cout << ocluster->get_cluster_id() << " " << mcell << std::endl;
    //             int time_slice = mcell->GetTimeSlice();
    //             ncluster->AddCell(mcell, time_slice);
    //         }
    //         live_clusters.erase(find(live_clusters.begin(), live_clusters.end(), ocluster));
    //         cluster_length_map.erase(ocluster);
    //         delete ocluster;
    //     }
    //     std::vector<int> range_v1 = ncluster->get_uvwt_range();
    //     double length_1 = sqrt(2. / 3. *
    //                                (pow(pitch_u * range_v1.at(0), 2) + pow(pitch_v * range_v1.at(1), 2) +
    //                                 pow(pitch_w * range_v1.at(2), 2)) +
    //                            pow(time_slice_width * range_v1.at(3), 2));
    //     cluster_length_map[ncluster] = length_1;
    //     // std::cout << std::endl;
    // }

    // // delete clusters ...
    // for (auto it = to_be_removed_clusters.begin(); it != to_be_removed_clusters.end(); it++) {
    //     Cluster *ocluster = *it;
    //     live_clusters.erase(find(live_clusters.begin(), live_clusters.end(), ocluster));
    //     cluster_length_map.erase(ocluster);
    //     delete ocluster;
    // } 
}
