#include <WireCellClus/ClusteringFuncs.h>

// The original developers do not care.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"


using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;


// #define __DEBUG__
#ifdef __DEBUG__
#define LogDebug(x) std::cout << "[yuhw]: " << __LINE__ << " : " << x << std::endl
#else
#define LogDebug(x)
#endif

void WireCell::PointCloud::Facade::clustering_connect1(Grouping& live_grouping)
{
    const auto &tp = live_grouping.get_params();
    auto global_point_cloud = std::make_shared<DynamicPointCloud>(tp.angle_u, tp.angle_v, tp.angle_w);
    for (const Cluster *cluster : live_grouping.children()) {
        global_point_cloud->add_points(cluster, 0);
    }
    std::map<int, std::pair<double, double>>& dead_u_index = live_grouping.get_dead_winds(0, 0);
    std::map<int, std::pair<double, double>>& dead_v_index = live_grouping.get_dead_winds(0, 1);
    std::map<int, std::pair<double, double>>& dead_w_index = live_grouping.get_dead_winds(0, 2);

    LogDebug("global_point_cloud.get_num_points() " << global_point_cloud->get_num_points());
    LogDebug("dead_u_index.size() " << dead_u_index.size() << " dead_v_index.size() " << dead_v_index.size() << " dead_w_index.size() " << dead_w_index.size());
    // sort the clusters length ...
    std::vector<Cluster *> live_clusters = live_grouping.children();  // copy
    std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster *cluster1, const Cluster *cluster2) {
        return cluster1->get_length() > cluster2->get_length();
    });

    // double time_slice_width = tp.nticks_live_slice * tp.tick_drift;

    auto global_skeleton_cloud = std::make_shared<DynamicPointCloud>(tp.angle_u, tp.angle_v, tp.angle_w);

    double extending_dis = 50 * units::cm;
    double angle = 7.5;
    double loose_dis_cut = 7.5 * units::cm;

    // std::set<std::pair<Cluster *, Cluster *>> to_be_merged_pairs;
    cluster_connectivity_graph_t  g;
    std::unordered_map<int, int> ilive2desc;  // added live index to graph descriptor
    std::map<const Cluster*, int> map_cluster_index;
    for (const Cluster* live : live_grouping.children()) {
        size_t ilive = map_cluster_index.size();
        map_cluster_index[live] = ilive;
        ilive2desc[ilive] = boost::add_vertex(ilive, g);
    }

    // WCP::WCPointCloud<double> &global_cloud = global_skeleton_cloud->get_cloud();
    // std::vector<int> global_cloud_indices_u;
    // std::vector<int> global_cloud_indices_v;
    // std::vector<int> global_cloud_indices_w;

    std::map<const Cluster *, geo_point_t> map_cluster_dir1;
    std::map<const Cluster *, geo_point_t> map_cluster_dir2;

    geo_point_t drift_dir(1, 0, 0);
    const auto [angle_u,angle_v,angle_w] = live_grouping.wire_angles();
    geo_point_t U_dir(0,cos(angle_u),sin(angle_u));
    geo_point_t V_dir(0,cos(angle_v),sin(angle_v));
    geo_point_t W_dir(0,cos(angle_w),sin(angle_w));
    // geo_point_t U_dir(0, cos(60. / 180. * 3.1415926), sin(60. / 180. * 3.1415926));
    // geo_point_t V_dir(0, cos(60. / 180. * 3.1415926), -sin(60. / 180. * 3.1415926));
    // geo_point_t W_dir(0, 1, 0);

    for (size_t i = 0; i != live_clusters.size(); i++) {
        Cluster *cluster = live_clusters.at(i);
        assert (cluster->npoints() > 0); // preempt segfault in get_two_extreme_points()

        // if (cluster->get_length()/units::cm>5){
        //     std::cout << "Connect 0: " << cluster->get_length()/units::cm << " " << cluster->get_center() << std::endl;
        // }

        LogDebug("#b " << cluster->nchildren() << " length " << cluster->get_length());

        #ifdef __DEBUG__
        if (cluster->nchildren() == 84) break;
        #endif
        // cluster->Create_point_cloud();

        std::pair<geo_point_t, geo_point_t> extreme_points = cluster->get_two_extreme_points();
        LogDebug("#b " << cluster->nchildren() << " extreme_points " << extreme_points.first << " " << extreme_points.second);
        geo_point_t main_dir(extreme_points.second.x() - extreme_points.first.x(),
                          extreme_points.second.y() - extreme_points.first.y(),
                          extreme_points.second.z() - extreme_points.first.z());
        geo_point_t dir1, dir2;
        bool flag_para_1 = false;
        bool flag_prol_1 = false;
        bool flag_para_2 = false;
        bool flag_prol_2 = false;

        if (main_dir.magnitude() > 10 * units::cm &&
            fabs(main_dir.angle(drift_dir) - 3.1415926 / 2.) < 5 * 3.1415926 / 180.) {
            dir1 = main_dir;
            dir1 = dir1* -1;
            dir2 = main_dir;
        }
        else if (cluster->get_length() > 25 * units::cm) {
            dir1 = cluster->vhough_transform(extreme_points.first, 80 * units::cm);
            if (dir1.magnitude() != 0) dir1 = dir1.norm();
            if (fabs(dir1.angle(drift_dir) - 3.1415926 / 2.) < 5 * 3.1415926 / 180.) {
                dir1.set(dir1.x(), (extreme_points.second.y() - extreme_points.first.y()) / main_dir.magnitude(),
                            (extreme_points.second.z() - extreme_points.first.z()) / main_dir.magnitude());
                dir1 = dir1 * -1;
            }
            dir2 = cluster->vhough_transform(extreme_points.second, 80 * units::cm);
            if (dir2.magnitude() != 0) dir2 = dir2.norm();
            if (fabs(dir2.angle(drift_dir) - 3.1415926 / 2.) < 5 * 3.1415926 / 180.) {
                dir2.set(dir2.x(), (extreme_points.second.y() - extreme_points.first.y()) / main_dir.magnitude(),
                            (extreme_points.second.z() - extreme_points.first.z()) / main_dir.magnitude());
            }
            if (dir1.dot(main_dir) > 0) dir1 *= -1;
            if (dir2.dot(dir1) > 0) dir2 *= -1;
        }
        else {
            dir1 = global_point_cloud->vhough_transform(extreme_points.first, extending_dis);
            dir2 = global_point_cloud->vhough_transform(extreme_points.second, extending_dis);
            if (dir1.dot(main_dir) > 0) dir1 *= -1;
            if (dir2.dot(dir1) > 0) dir2 *= -1;
        }
        LogDebug("#b " << cluster->nchildren() << " dir1 " << dir1 << " dir2 " << dir2);

        bool flag_add_dir1 = true;
        bool flag_add_dir2 = true;
        map_cluster_dir1[cluster] = dir1;
        map_cluster_dir2[cluster] = dir2;

        if (fabs(dir1.angle(drift_dir) - 3.1415926 / 2.) < 7.5 * 3.1415926 / 180.) {
            flag_para_1 = true;
        }
        else {
            geo_point_t tempV1(0, dir1.y(), dir1.z());
            geo_point_t tempV5;
            double angle1 = tempV1.angle(U_dir);
            tempV5.set(fabs(dir1.x()), sqrt(pow(dir1.y(), 2) + pow(dir1.z(), 2)) * sin(angle1), 0);
            angle1 = tempV5.angle(drift_dir);

            if (angle1 < 7.5 / 180. * 3.1415926) {
                flag_prol_1 = true;
            }
            else {
                angle1 = tempV1.angle(V_dir);
                tempV5.set(fabs(dir1.x()), sqrt(pow(dir1.y(), 2) + pow(dir1.z(), 2)) * sin(angle1), 0);
                angle1 = tempV5.angle(drift_dir);

                if (angle1 < 7.5 / 180. * 3.1415926) {
                    flag_prol_1 = true;
                }
                else {
                    angle1 = tempV1.angle(W_dir);
                    tempV5.set(fabs(dir1.x()), sqrt(pow(dir1.y(), 2) + pow(dir1.z(), 2)) * sin(angle1), 0);
                    angle1 = tempV5.angle(drift_dir);

                    if (angle1 < 7.5 / 180. * 3.1415926) {
                        flag_prol_1 = true;
                    }
                }
            }
        }

        if (fabs(dir2.angle(drift_dir) - 3.1415926 / 2.) < 7.5 * 3.1415926 / 180.) {
            flag_para_2 = true;
        }
        else {
            geo_point_t tempV2(0, dir2.y(), dir2.z());
            geo_point_t tempV6;
            double angle2 = tempV2.angle(U_dir);
            tempV6.set(fabs(dir2.x()), sqrt(pow(dir2.y(), 2) + pow(dir2.z(), 2)) * sin(angle2), 0);
            angle2 = tempV6.angle(drift_dir);
            if (angle2 < 7.5 / 180. * 3.1415926) {
                flag_prol_2 = true;
            }
            else {
                angle2 = tempV2.angle(V_dir);
                tempV6.set(fabs(dir2.x()), sqrt(pow(dir2.y(), 2) + pow(dir2.z(), 2)) * sin(angle2), 0);
                angle2 = tempV6.angle(drift_dir);
                if (angle2 < 7.5 / 180. * 3.1415926) {
                    flag_prol_2 = true;
                }
                else {
                    angle2 = tempV2.angle(W_dir);
                    tempV6.set(fabs(dir2.x()), sqrt(pow(dir2.y(), 2) + pow(dir2.z(), 2)) * sin(angle2), 0);
                    angle2 = tempV6.angle(drift_dir);
                    if (angle2 < 7.5 / 180. * 3.1415926) {
                        flag_prol_2 = true;
                    }
                }
            }
        }

        if ((flag_para_1 || flag_prol_1) && cluster->get_length() < 15 * units::cm) {
            flag_add_dir1 = true;
        }
        else if (cluster->get_length() >= 15 * units::cm) {
            flag_add_dir1 = true;
        }
        else {
            flag_add_dir1 = false;
        }

        if ((flag_para_2 || flag_prol_2) && cluster->get_length() < 15 * units::cm) {
            flag_add_dir2 = true;
        }
        else if (cluster->get_length() >= 15 * units::cm) {
            flag_add_dir2 = true;
        }
        else {
            flag_add_dir2 = false;
        }

        if (fabs(dir1.angle(drift_dir) - 3.1415926 / 2.) < 5 * 3.1415926 / 180.) {
            flag_para_1 = true;
        }
        else {
            flag_para_1 = false;
        }
        if (fabs(dir2.angle(drift_dir) - 3.1415926 / 2.) < 5 * 3.1415926 / 180.) {
            flag_para_2 = true;
        }
        else {
            flag_para_2 = false;
        }
        
        LogDebug("#b " << cluster->nchildren() << " flag_para_1 " << flag_para_1 << " flag_prol_1 " << flag_prol_1 << " flag_para_2 " << flag_para_2 << " flag_prol_2 " << flag_prol_2);

        

        if (i == 0) {
            if (flag_para_1 || flag_prol_1) {
                global_skeleton_cloud->add_points(cluster, extreme_points.first, dir1, extending_dis * 3, 1.2 * units::cm,
                                                angle);
                dir1 *= -1;
                global_skeleton_cloud->add_points(cluster, extreme_points.first, dir1, extending_dis * 3, 1.2 * units::cm,
                                                angle);
            }
            else {
                global_skeleton_cloud->add_points(cluster, extreme_points.first, dir1, extending_dis, 1.2 * units::cm,
                                                angle);
                dir1 *= -1;
                global_skeleton_cloud->add_points(cluster, extreme_points.first, dir1, extending_dis, 1.2 * units::cm,
                                                angle);
            }

            if (flag_para_2 || flag_prol_2) {
                global_skeleton_cloud->add_points(cluster, extreme_points.second, dir2, extending_dis * 3.0,
                                                1.2 * units::cm, angle);
                dir2 *= -1;
                global_skeleton_cloud->add_points(cluster, extreme_points.second, dir2, extending_dis * 3.0,
                                                1.2 * units::cm, angle);
            }
            else {
                global_skeleton_cloud->add_points(cluster, extreme_points.second, dir2, extending_dis, 1.2 * units::cm,
                                                angle);
                dir2 *= -1;
                global_skeleton_cloud->add_points(cluster, extreme_points.second, dir2, extending_dis, 1.2 * units::cm,
                                                angle);
            }
        }
        else {
            if (cluster->get_length() < 100 * units::cm ||
                fabs(dir2.angle(drift_dir) - 3.1415926 / 2.) < 5 * 3.1415926 / 180. &&
                    fabs(dir1.angle(drift_dir) - 3.1415926 / 2.) < 5 * 3.1415926 / 180. &&
                    cluster->get_length() < 200 * units::cm) {
                // WCP::WCPointCloud<double> &cloud = cluster->get_point_cloud()->get_cloud();
                LogDebug("#b " << cluster->nchildren() << " gsc " << global_skeleton_cloud->get_num_points());
                int num_total_points = cluster->npoints();
                const auto& winds = cluster->wire_indices();
                int num_unique[3] = {0, 0, 0};            // points that are unique (not agree with any other clusters)
                std::map<const Cluster *, int> map_cluster_num[3];
                for (int j = 0; j != num_total_points; j++) {
                    geo_point_t test_point(cluster->point3d(j).x(), cluster->point3d(j).y(), cluster->point3d(j).z());

                    bool flag_dead = false;
                    if (dead_u_index.find(winds[0][j]) != dead_u_index.end()) {
                        if (cluster->point3d(j).x() >= dead_u_index[winds[0][j]].first &&
                            cluster->point3d(j).x() <= dead_u_index[winds[0][j]].second) {
                            flag_dead = true;
                        }
                    }

                    if (!flag_dead) {
                        std::vector<std::tuple<double, const Cluster *, size_t>> results =
                            global_skeleton_cloud->get_2d_points_info(test_point, loose_dis_cut, 0);
                        LogDebug("#b " << cluster->nchildren() << " test_point " << test_point << " loose_dis_cut " << loose_dis_cut << " results.size() " << results.size());
                        bool flag_unique = true;
                        if (results.size() > 0) {
                            std::set<const Cluster *> temp_clusters;
                            for (size_t k = 0; k != results.size(); k++) {
                                // LogDebug("#b " << cluster->nchildren() << " results.at(k) " << std::get<0>(results.at(k)) << " " << global_skeleton_cloud->dist_cut(0,std::get<2>(results.at(k))));
                                if (std::get<0>(results.at(k)) <
                                    global_skeleton_cloud->dist_cut(0,std::get<2>(results.at(k)))) {
                                    flag_unique = false;
                                    temp_clusters.insert(std::get<1>(results.at(k)));
                                }
                            }
                            for (auto it = temp_clusters.begin(); it != temp_clusters.end(); it++) {
                                if (map_cluster_num[0].find(*it) == map_cluster_num[0].end()) {
                                    map_cluster_num[0][*it] = 1;
                                }
                                else {
                                    map_cluster_num[0][*it]++;
                                }
                            }
                        }
                        if (flag_unique) num_unique[0]++;
                    }
                    else {
                        std::vector<std::tuple<double, const Cluster *, size_t>> results =
                            global_skeleton_cloud->get_2d_points_info(test_point, loose_dis_cut, 0);
                        bool flag_unique = true;
                        if (results.size() > 0) {
                            std::set<const Cluster *> temp_clusters;
                            for (size_t k = 0; k != results.size(); k++) {
                                if (std::get<0>(results.at(k)) < loose_dis_cut / 3. * 2.) {
                                    flag_unique = false;
                                    temp_clusters.insert(std::get<1>(results.at(k)));
                                }
                            }
                            for (auto it = temp_clusters.begin(); it != temp_clusters.end(); it++) {
                                if (map_cluster_num[0].find(*it) == map_cluster_num[0].end()) {
                                    map_cluster_num[0][*it] = 1;
                                }
                                else {
                                    map_cluster_num[0][*it]++;
                                }
                            }
                        }
                        if (flag_unique) num_unique[0]++;
                    }

                    flag_dead = false;
                    if (dead_v_index.find(winds[1][j]) != dead_v_index.end()) {
                        if (cluster->point3d(j).x() >= dead_v_index[winds[1][j]].first &&
                            cluster->point3d(j).x() <= dead_v_index[winds[1][j]].second) {
                            flag_dead = true;
                        }
                    }

                    if (!flag_dead) {
                        std::vector<std::tuple<double, const Cluster *, size_t>> results =
                            global_skeleton_cloud->get_2d_points_info(test_point, loose_dis_cut, 1);
                        bool flag_unique = true;
                        if (results.size() > 0) {
                            std::set<const Cluster *> temp_clusters;
                            for (size_t k = 0; k != results.size(); k++) {
                                if (std::get<0>(results.at(k)) <
                                    global_skeleton_cloud->dist_cut(1,std::get<2>(results.at(k)))) {
                                    flag_unique = false;
                                    temp_clusters.insert(std::get<1>(results.at(k)));
                                }
                            }
                            for (auto it = temp_clusters.begin(); it != temp_clusters.end(); it++) {
                                if (map_cluster_num[1].find(*it) == map_cluster_num[1].end()) {
                                    map_cluster_num[1][*it] = 1;
                                }
                                else {
                                    map_cluster_num[1][*it]++;
                                }
                            }
                        }
                        if (flag_unique) num_unique[1]++;
                    }
                    else {
                        std::vector<std::tuple<double, const Cluster *, size_t>> results =
                            global_skeleton_cloud->get_2d_points_info(test_point, loose_dis_cut, 1);
                        bool flag_unique = true;
                        if (results.size() > 0) {
                            std::set<const Cluster *> temp_clusters;
                            for (size_t k = 0; k != results.size(); k++) {
                                if (std::get<0>(results.at(k)) < loose_dis_cut / 3. * 2.) {
                                    flag_unique = false;
                                    temp_clusters.insert(std::get<1>(results.at(k)));
                                }
                            }
                            for (auto it = temp_clusters.begin(); it != temp_clusters.end(); it++) {
                                if (map_cluster_num[1].find(*it) == map_cluster_num[1].end()) {
                                    map_cluster_num[1][*it] = 1;
                                }
                                else {
                                    map_cluster_num[1][*it]++;
                                }
                            }
                        }
                        if (flag_unique) num_unique[1]++;
                    }

                    flag_dead = false;
                    if (dead_w_index.find(winds[2][j]) != dead_w_index.end()) {
                        if (cluster->point3d(j).x() >= dead_w_index[winds[2][j]].first &&
                            cluster->point3d(j).x() <= dead_w_index[winds[2][j]].second) {
                            flag_dead = true;
                        }
                    }

                    if (!flag_dead) {
                        std::vector<std::tuple<double, const Cluster *, size_t>> results =
                            global_skeleton_cloud->get_2d_points_info(test_point, loose_dis_cut, 2);
                        bool flag_unique = true;
                        if (results.size() > 0) {
                            std::set<const Cluster *> temp_clusters;
                            for (size_t k = 0; k != results.size(); k++) {
                                if (std::get<0>(results.at(k)) <
                                    global_skeleton_cloud->dist_cut(2,std::get<2>(results.at(k)))) {
                                    flag_unique = false;
                                    temp_clusters.insert(std::get<1>(results.at(k)));
                                }
                            }
                            for (auto it = temp_clusters.begin(); it != temp_clusters.end(); it++) {
                                if (map_cluster_num[2].find(*it) == map_cluster_num[2].end()) {
                                    map_cluster_num[2][*it] = 1;
                                }
                                else {
                                    map_cluster_num[2][*it]++;
                                }
                            }
                        }
                        if (flag_unique) num_unique[2]++;
                    }
                    else {
                        std::vector<std::tuple<double, const Cluster *, size_t>> results =
                            global_skeleton_cloud->get_2d_points_info(test_point, loose_dis_cut, 2);
                        bool flag_unique = true;
                        if (results.size() > 0) {
                            std::set<const Cluster *> temp_clusters;
                            for (size_t k = 0; k != results.size(); k++) {
                                if (std::get<0>(results.at(k)) < loose_dis_cut / 3. * 2.) {
                                    flag_unique = false;
                                    temp_clusters.insert(std::get<1>(results.at(k)));
                                }
                            }
                            for (auto it = temp_clusters.begin(); it != temp_clusters.end(); it++) {
                                if (map_cluster_num[2].find(*it) == map_cluster_num[2].end()) {
                                    map_cluster_num[2][*it] = 1;
                                }
                                else {
                                    map_cluster_num[2][*it]++;
                                }
                            }
                        }
                        if (flag_unique) num_unique[2]++;
                    }
                } // loop over points
                LogDebug("num_unique " << num_unique[0] << " " << num_unique[1] << " " << num_unique[2]);
                LogDebug("map_cluster_num " << map_cluster_num[0].size() << " " << map_cluster_num[1].size() << " " << map_cluster_num[2].size());
                if (cluster->nchildren() == 84) {
                    for (auto it = map_cluster_num[0].begin(); it != map_cluster_num[0].end(); it++) {
                        LogDebug("map_cluster_num[0] " << it->first->nchildren() << " " << it->second);
                    }
                }

                bool flag_merge = false;

                {
                    const Cluster *max_cluster_u = 0, *max_cluster_v = 0, *max_cluster_w = 0;
                    int max_value_u[3] = {0, 0, 0};
                    int max_value_v[3] = {0, 0, 0};
                    int max_value_w[3] = {0, 0, 0};

                    int max_value[3] = {0, 0, 0};
                    const Cluster *max_cluster = 0;

                    for (auto it = map_cluster_num[0].begin(); it != map_cluster_num[0].end(); it++) {
                        if (it->second > max_value_u[0]) {
                            max_value_u[0] = it->second;
                            max_cluster_u = it->first;

                            if (map_cluster_num[1].find(max_cluster_u) != map_cluster_num[1].end()) {
                                max_value_u[1] = map_cluster_num[1][max_cluster_u];
                            }
                            else {
                                max_value_u[1] = 0;
                            }

                            if (map_cluster_num[2].find(max_cluster_u) != map_cluster_num[2].end()) {
                                max_value_u[2] = map_cluster_num[2][max_cluster_u];
                            }
                            else {
                                max_value_u[2] = 0;
                            }
                        }
                    }
                    for (auto it = map_cluster_num[1].begin(); it != map_cluster_num[1].end(); it++) {
                        if (it->second > max_value_v[1]) {
                            max_value_v[1] = it->second;
                            max_cluster_v = it->first;

                            if (map_cluster_num[0].find(max_cluster_v) != map_cluster_num[0].end()) {
                                max_value_v[0] = map_cluster_num[0][max_cluster_v];
                            }
                            else {
                                max_value_v[0] = 0;
                            }
                            if (map_cluster_num[2].find(max_cluster_v) != map_cluster_num[2].end()) {
                                max_value_v[2] = map_cluster_num[2][max_cluster_v];
                            }
                            else {
                                max_value_v[2] = 0;
                            }
                        }
                    }
                    for (auto it = map_cluster_num[2].begin(); it != map_cluster_num[2].end(); it++) {
                        if (it->second > max_value_w[2]) {
                            max_value_w[2] = it->second;
                            max_cluster_w = it->first;

                            if (map_cluster_num[1].find(max_cluster_w) != map_cluster_num[1].end()) {
                                max_value_w[1] = map_cluster_num[1][max_cluster_w];
                            }
                            else {
                                max_value_w[1] = 0;
                            }
                            if (map_cluster_num[0].find(max_cluster_w) != map_cluster_num[0].end()) {
                                max_value_w[0] = map_cluster_num[0][max_cluster_w];
                            }
                            else {
                                max_value_w[0] = 0;
                            }
                        }
                    }

                    if ((max_value_u[0] > 0.33 * num_total_points || max_value_u[0] > 100) &&
                        (max_value_u[1] > 0.33 * num_total_points || max_value_u[1] > 100) &&
                        (max_value_u[2] > 0.33 * num_total_points || max_value_u[2] > 100)) {
                        if (max_value_u[0] + max_value_u[1] + max_value_u[2] >
                            max_value[0] + max_value[1] + max_value[2]) {
                            max_value[0] = max_value_u[0];
                            max_value[1] = max_value_u[1];
                            max_value[2] = max_value_u[2];
                            max_cluster = max_cluster_u;
                        }
                    }
                    if ((max_value_v[0] > 0.33 * num_total_points || max_value_v[0] > 100) &&
                        (max_value_v[1] > 0.33 * num_total_points || max_value_v[1] > 100) &&
                        (max_value_v[2] > 0.33 * num_total_points || max_value_v[2] > 100)) {
                        if (max_value_v[0] + max_value_v[1] + max_value_v[2] >
                            max_value[0] + max_value[1] + max_value[2]) {
                            max_value[0] = max_value_v[0];
                            max_value[1] = max_value_v[1];
                            max_value[2] = max_value_v[2];
                            max_cluster = max_cluster_v;
                        }
                    }
                    if ((max_value_w[0] > 0.33 * num_total_points || max_value_w[0] > 100) &&
                        (max_value_w[1] > 0.33 * num_total_points || max_value_w[1] > 100) &&
                        (max_value_w[2] > 0.33 * num_total_points || max_value_w[2] > 100)) {
                        if (max_value_w[0] + max_value_w[1] + max_value_w[2] >
                            max_value[0] + max_value[1] + max_value[2]) {
                            max_value[0] = max_value_w[0];
                            max_value[1] = max_value_w[1];
                            max_value[2] = max_value_w[2];
                            max_cluster = max_cluster_w;
                        }
                    }

                    // if (max_cluster != 0)
                    // if (fabs(cluster->get_length()/units::cm - 50) < 5 ){
                    //     std::cout << "Check: " << cluster->get_length()/units::cm << " " << max_cluster->get_length()/units::cm << " " << cluster->get_center() << " " << max_cluster->get_center() << " " << max_value[0] << " " << max_value[1] << " " << max_value[2] << " " << num_total_points << " " << num_unique[0] << " " << num_unique[1] << " " << num_unique[2] << " "  << std::endl;
                    // }

                    // if overlap a lot merge
                    if ((max_value[0] + max_value[1] + max_value[2]) >
                            0.75 * (num_total_points + num_total_points + num_total_points) &&
                        ((num_unique[1] + num_unique[0] + num_unique[2]) < 0.24 * num_total_points ||
                         ((num_unique[1] + num_unique[0] + num_unique[2]) < 0.45 * num_total_points &&
                          (num_unique[1] + num_unique[0] + num_unique[2]) < 25))) {

                        

                        if (fabs(dir1.angle(map_cluster_dir1[max_cluster]) - 3.1415926 / 2.) >= 70 * 3.1415926 / 180. ||
                            fabs(dir1.angle(map_cluster_dir2[max_cluster]) - 3.1415926 / 2.) >= 70 * 3.1415926 / 180. ||
                            fabs(dir2.angle(map_cluster_dir1[max_cluster]) - 3.1415926 / 2.) >= 70 * 3.1415926 / 180. ||
                            fabs(dir2.angle(map_cluster_dir2[max_cluster]) - 3.1415926 / 2.) >= 70 * 3.1415926 / 180.) {
                            flag_merge = true;
                            // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster));
                            boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                            ilive2desc[map_cluster_index[max_cluster]], g);
                            // std::cout << "Connect 1 1: " << cluster->get_length()/units::cm << " " << max_cluster->get_length()/units::cm << " " << cluster->get_center() << " " << max_cluster->get_center() << std::endl;
                            // curr_cluster = max_cluster;
                        }

                        LogDebug("max_cluster_u #b " << max_cluster_u->nchildren() << " max_cluster_v #b " << max_cluster_v->nchildren() << " max_cluster_w #b " << max_cluster_w->nchildren());
                        LogDebug("max_cluster #b " << max_cluster->nchildren() << " map_cluster_dir1 " << map_cluster_dir1[max_cluster] << " map_cluster_dir2 " << map_cluster_dir2[max_cluster]);
                        if (fabs(dir1.angle(map_cluster_dir1[max_cluster]) - 3.1415926 / 2.) < 75 * 3.1415926 / 180. &&
                            fabs(dir1.angle(map_cluster_dir2[max_cluster]) - 3.1415926 / 2.) < 75 * 3.1415926 / 180.) {
                            flag_add_dir1 = false;
                        }
                        else {
                            flag_add_dir1 = true;
                        }
                        if (fabs(dir2.angle(map_cluster_dir1[max_cluster]) - 3.1415926 / 2.) < 75 * 3.1415926 / 180. &&
                            fabs(dir2.angle(map_cluster_dir2[max_cluster]) - 3.1415926 / 2.) < 75 * 3.1415926 / 180.) {
                            flag_add_dir2 = false;
                        }
                        else {
                            flag_add_dir2 = true;
                        }
                    }

                    if ((max_value[0] + max_value[1] + max_value[2]) > 300 && !flag_merge) {
                        if (cluster->get_length() > 25 * units::cm ||
                            max_cluster->get_length() > 25 * units::cm) {
                            // if overlap significant, compare the PCA
                            // cluster->Calc_PCA();
                            geo_point_t p1_c = cluster->get_center();
                            geo_point_t p1_dir(cluster->get_pca_axis(0).x(), cluster->get_pca_axis(0).y(),
                                            cluster->get_pca_axis(0).z());
                            // max_cluster->Calc_PCA();
                            geo_point_t p2_c = max_cluster->get_center();
                            geo_point_t p2_dir(max_cluster->get_pca_axis(0).x(), max_cluster->get_pca_axis(0).y(),
                                            max_cluster->get_pca_axis(0).z());

                            double angle_diff = p1_dir.angle(p2_dir) / 3.1415926 * 180.;
                            double angle1_drift = p1_dir.angle(drift_dir) / 3.1415926 * 180.;
                            double angle2_drift = p2_dir.angle(drift_dir) / 3.1415926 * 180.;
                            Ray l1(p1_c, p1_c+p1_dir);
                            Ray l2(p2_c, p2_c+p2_dir);
                            // double dis = l1.closest_dis(l2);
                            double dis = ray_length(ray_pitch(l1, l2));
                            double dis1 =
                                sqrt(pow(p1_c.x() - p2_c.x(), 2) + pow(p1_c.y() - p2_c.y(), 2) + pow(p1_c.z() - p2_c.z(), 2));

                            // if (fabs(cluster->get_length()/units::cm - 50) < 5 ){
                            //     std::cout << "Check 1: " << cluster->get_length()/units::cm << " " << max_cluster->get_length()/units::cm << " " << angle_diff << " " << dis/units::cm << " " << dis1/units::cm << " " << angle1_drift << " " << angle2_drift << std::endl;
                            // }

                            if ((angle_diff < 5 || angle_diff > 175 ||
                                 fabs(angle1_drift - 90) < 5 && fabs(angle2_drift - 90) < 5 &&
                                     fabs(angle1_drift - 90) + fabs(angle2_drift - 90) < 6 &&
                                     (angle_diff < 30 || angle_diff > 150)) &&
                                    dis < 1.5 * units::cm ||
                                (angle_diff < 10 || angle_diff > 170) && dis < 0.9 * units::cm &&
                                    dis1 > (cluster->get_length() + max_cluster->get_length()) / 3.) {
                                // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster));
                                boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                ilive2desc[map_cluster_index[max_cluster]], g);

                                // std::cout << "Connect 1 2: " << cluster->get_length()/units::cm << " " << max_cluster->get_length()/units::cm << " " << cluster->get_center() << " " << max_cluster->get_center() << std::endl;
                                // curr_cluster = max_cluster;
                                flag_merge = true;
                            }
                            else if (((angle_diff < 5 || angle_diff > 175) && dis < 2.5 * units::cm ||
                                      (angle_diff < 10 || angle_diff > 170) && dis < 1.2 * units::cm) &&
                                     dis1 > (cluster->get_length() + max_cluster->get_length()) / 3.) {
                                // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster));
                                boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                ilive2desc[map_cluster_index[max_cluster]], g);

                                // std::cout << "Connect 1 3: " << cluster->get_length()/units::cm << " " << max_cluster->get_length()/units::cm << " " << cluster->get_center() << " " << max_cluster->get_center() << std::endl;
                                flag_merge = true;
                            }

                            if ((fabs(dir2.angle(drift_dir) - 3.1415926 / 2.) < 5 * 3.1415926 / 180. &&
                                 fabs(dir1.angle(drift_dir) - 3.1415926 / 2.) < 5 * 3.1415926 / 180.) &&
                                (max_value[0] + max_value[1] + max_value[2]) >
                                    0.7 * (num_total_points + num_total_points + num_total_points)) {
                                // to_be_merged_pairs.insert(std::make_pair(cluster, max_cluster));
                                boost::add_edge(ilive2desc[map_cluster_index[cluster]],
                                                ilive2desc[map_cluster_index[max_cluster]], g);
                                // std::cout << "Connect 1 4: " << cluster->get_length()/units::cm << " " << max_cluster->get_length()/units::cm << " " << cluster->get_center() << " " << max_cluster->get_center() << std::endl;
                                flag_merge = true;
                            }
                        }
                    }
                }
            }  // length cut ...

            LogDebug("#b " << cluster->nchildren() << " flag_add_dir1 " << flag_add_dir1 << " flag_add_dir2 " << flag_add_dir2);

            // if (cluster->get_length()/units::cm>5){
            //     std::cout << "Connect 0-1: " << cluster->get_length()/units::cm << " " << cluster->get_center() << " " << flag_add_dir1 << " " << flag_add_dir2 << " " << flag_para_1 << " " << flag_prol_1 << " " << flag_para_2 << " " << flag_prol_2 << " " << extreme_points.first << " " << extreme_points.second << " " << dir1 << " " << dir2 << " " << extending_dis << " " << angle << std::endl;
            // }

            if (flag_add_dir1) {
                // add extension points in ...
                if (flag_para_1 || flag_prol_1) {
                    global_skeleton_cloud->add_points(cluster, extreme_points.first, dir1, extending_dis * 3,
                                                    1.2 * units::cm, angle);
                    dir1 *= -1;
                    global_skeleton_cloud->add_points(cluster, extreme_points.first, dir1, extending_dis * 3,
                                                    1.2 * units::cm, angle);
                }
                else {
                    global_skeleton_cloud->add_points(cluster, extreme_points.first, dir1, extending_dis, 1.2 * units::cm,
                                                    angle);
                    dir1 *= -1;
                    global_skeleton_cloud->add_points(cluster, extreme_points.first, dir1, extending_dis, 1.2 * units::cm,
                                                    angle);
                }
            }

            if (flag_add_dir2) {
                if (flag_para_2 || flag_prol_2) {
                    global_skeleton_cloud->add_points(cluster, extreme_points.second, dir2, extending_dis * 3.0,
                                                    1.2 * units::cm, angle);
                    dir2 *= -1;
                    global_skeleton_cloud->add_points(cluster, extreme_points.second, dir2, extending_dis * 3.0,
                                                    1.2 * units::cm, angle);
                }
                else {
                    global_skeleton_cloud->add_points(cluster, extreme_points.second, dir2, extending_dis,
                                                    1.2 * units::cm, angle);
                    dir2 *= -1;
                    global_skeleton_cloud->add_points(cluster, extreme_points.second, dir2, extending_dis,
                                                    1.2 * units::cm, angle);
                }
            }
        }  // not the first cluster ...
    }  // loop over clusters ...

    LogDebug("#edges " << boost::num_edges(g) << " #vertices " << boost::num_vertices(g));

    // merge clusters

    /**
     * round1
    */
    // to_be_merged_pairs -> merge_clusters
    // std::vector<std::set<Cluster *>> merge_clusters;
    // for (auto it = to_be_merged_pairs.begin(); it != to_be_merged_pairs.end(); it++) {
    //     Cluster *cluster1 = (*it).first;
    //     Cluster *cluster2 = (*it).second;
    //     //  LogDebug(cluster1 << " " << cluster2 << " " << cluster1->get_cluster_id() << " " <<
    //     //  cluster2->get_cluster_id());

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

    // // merge_clusters -> new_clusters
    // WCP::ClusterSelection new_clusters;

    // // merge clusters into new clusters, delete old clusters
    // for (auto it = merge_clusters.begin(); it != merge_clusters.end(); it++) {
    //     std::set<Cluster *> &clusters = (*it);
    //     Cluster *ncluster = new Cluster((*clusters.begin())->get_cluster_id());
    //     live_clusters.push_back(ncluster);

    //     new_clusters.push_back(ncluster);

    //     for (auto it1 = clusters.begin(); it1 != clusters.end(); it1++) {
    //         Cluster *ocluster = *(it1);
    //         // LogDebug(ocluster->get_cluster_id() << " ";
    //         SMGCSelection &mcells = ocluster->get_mcells();
    //         for (auto it2 = mcells.begin(); it2 != mcells.end(); it2++) {
    //             SlimMergeGeomCell *mcell = (*it2);
    //             // LogDebug(ocluster->get_cluster_id() << " " << mcell);
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
    //     // LogDebug(std::endl;
    // }
    /**
     * end of round1
    */

    cluster_set_t new_clusters;
    merge_clusters(g, live_grouping, new_clusters);
    live_clusters.clear();
    live_clusters = live_grouping.children();  // copy
    std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster *cluster1, const Cluster *cluster2) {
        return cluster1->get_length() > cluster2->get_length();
    });
    cluster_connectivity_graph_t  g2;
    ilive2desc.clear();  // added live index to graph descriptor
    map_cluster_index.clear();
    for (const Cluster* live : live_grouping.children()) {
        size_t ilive = map_cluster_index.size();
        map_cluster_index[live] = ilive;
        ilive2desc[ilive] = boost::add_vertex(ilive, g2);
    }

    // to_be_merged_pairs.clear(); // clear it for other usage ...
    for (auto it = new_clusters.begin(); it != new_clusters.end(); it++) {
        const Cluster *cluster_1 = (*it);
        // cluster_1->Calc_PCA();
        geo_point_t p1_c = cluster_1->get_center();
        geo_point_t p1_dir(cluster_1->get_pca_axis(0).x(), cluster_1->get_pca_axis(0).y(), cluster_1->get_pca_axis(0).z());
        Ray l1(p1_c, p1_c+p1_dir);
        for (auto it1 = live_clusters.begin(); it1 != live_clusters.end(); it1++) {
            Cluster *cluster_2 = (*it1);
            if (cluster_2->get_length() < 3 * units::cm) continue;
            if (cluster_2 == cluster_1) continue;

            if (cluster_1->get_length() > 25 * units::cm || cluster_2->get_length() > 25 * units::cm ||
                (cluster_1->get_length() + cluster_2->get_length()) > 30 * units::cm) {
                // cluster_2->Calc_PCA();
                geo_point_t p2_c = cluster_2->get_center();
                geo_point_t p2_dir(cluster_2->get_pca_axis(0).x(), cluster_2->get_pca_axis(0).y(),
                                cluster_2->get_pca_axis(0).z());

                geo_point_t cc_dir(p2_c.x() - p1_c.x(), p2_c.y() - p1_c.y(), p2_c.z() - p1_c.z());

                double angle_diff = fabs(p1_dir.angle(p2_dir) - 3.1415926 / 2.) / 3.1415926 * 180.;
                double angle_diff1 = fabs(cc_dir.angle(p1_dir) - 3.1415926 / 2.) / 3.1415926 * 180;
                double angle_diff2 = fabs(cc_dir.angle(p2_dir) - 3.1415926 / 2.) / 3.1415926 * 180;

                Ray l2(p2_c, p2_c+p2_dir);
                // double dis = l1.closest_dis(l2);
                double dis = ray_length(ray_pitch(l1, l2));

                double dis1 = sqrt(pow(p1_c.x() - p2_c.x(), 2) + pow(p1_c.y() - p2_c.y(), 2) + pow(p1_c.z() - p2_c.z(), 2));

                if (p1_dir.magnitude() != 0) p1_dir = p1_dir.norm();
                if (p2_dir.magnitude() != 0) p2_dir = p2_dir.norm();

                // bool flag_merge = false;

                // if (cluster_1->get_length()>300*units::cm) std::cout << "Check 2: " << cluster_1->get_length()/units::cm << " " << cluster_2->get_length()/units::cm << " " << angle_diff << " " << dis/units::cm << " " << dis1/units::cm << " " << angle_diff1 << " " << angle_diff2 << std::endl;

                if (((angle_diff > 85) && (angle_diff1 > 90 - 1.5 * (90 - angle_diff)) &&
                         (angle_diff2 > 90 - 1.5 * (90 - angle_diff)) && dis < 2.5 * units::cm ||
                     (angle_diff > 80) && angle_diff1 > 80 && angle_diff2 > 80 && dis < 1.2 * units::cm) &&
                    dis1 > (cluster_2->get_length() + cluster_1->get_length()) / 3.) {
                    // to_be_merged_pairs.insert(std::make_pair(cluster_1, cluster_2));
                    boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
                                    ilive2desc[map_cluster_index[cluster_2]], g2);
                    // std::cout << "Connect 2: " << cluster_1->get_length()/units::cm << " " << cluster_2->get_length()/units::cm << " " << cluster_1->get_center() << " " << cluster_2->get_center() << std::endl;
                    // flag_merge = true;
                }
                else if ((angle_diff > 87) && (angle_diff1 > 90 - 1.5 * (90 - angle_diff)) &&
                         (angle_diff2 > 90 - 1.5 * (90 - angle_diff)) && dis < 4.0 * units::cm &&
                         dis1 > (cluster_2->get_length() + cluster_1->get_length()) / 2. &&
                         cluster_2->get_length() > 15 * units::cm &&
                         cluster_1->get_length() > 15 * units::cm &&
                         cluster_2->get_length() + cluster_1->get_length() > 45 * units::cm) {
                    // to_be_merged_pairs.insert(std::make_pair(cluster_1, cluster_2));
                    boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
                                    ilive2desc[map_cluster_index[cluster_2]], g2);
                    // std::cout << "Connect 2: " << cluster_1->get_length()/units::cm << " " << cluster_2->get_length()/units::cm << " " << cluster_1->get_center() << " " << cluster_2->get_center() << std::endl;
                    // flag_merge = true;
                }
            }
        }
    }

    new_clusters.clear();
    merge_clusters(g2, live_grouping, new_clusters);
    /**
     * round2
    */
    // merge_clusters.clear(); // clear it for other usage ...
    // for (auto it = to_be_merged_pairs.begin(); it != to_be_merged_pairs.end(); it++) {
    //     Cluster *cluster1 = (*it).first;
    //     Cluster *cluster2 = (*it).second;
    //     //  LogDebug(cluster1 << " " << cluster2 << " " << cluster1->get_cluster_id() << " " <<
    //     //  cluster2->get_cluster_id());

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

    // for (auto it = merge_clusters.begin(); it != merge_clusters.end(); it++) {
    //     std::set<Cluster *> &clusters = (*it);
    //     Cluster *ncluster = new Cluster((*clusters.begin())->get_cluster_id());
    //     live_clusters.push_back(ncluster);
    //     for (auto it1 = clusters.begin(); it1 != clusters.end(); it1++) {
    //         Cluster *ocluster = *(it1);
    //         // LogDebug(ocluster->get_cluster_id() << " ";
    //         SMGCSelection &mcells = ocluster->get_mcells();
    //         for (auto it2 = mcells.begin(); it2 != mcells.end(); it2++) {
    //             SlimMergeGeomCell *mcell = (*it2);
    //             // LogDebug(ocluster->get_cluster_id() << " " << mcell);
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
    //     // LogDebug(std::endl;
    // }
    /**
     * end of round2
    */
}
