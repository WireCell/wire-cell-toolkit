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


void WireCell::PointCloud::Facade::clustering_neutrino(Grouping &live_grouping, int num_try)
{
    std::vector<Cluster *> live_clusters = live_grouping.children();  // copy
    // sort the clusters by length using a lambda function
    std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster *cluster1, const Cluster *cluster2) {
        return cluster1->get_length() < cluster2->get_length();
    });

    const auto &mp = live_grouping.get_params();
    // this is for 4 time slices
    double time_slice_width = mp.nticks_live_slice * mp.tick_drift;

    geo_point_t drift_dir(1, 0, 0);
    geo_point_t vertical_dir(0, 1, 0);
    geo_point_t beam_dir(0, 0, 1);

    // find all the clusters that are inside the box ...
    std::vector<Cluster *> contained_clusters;
    std::vector<Cluster *> candidate_clusters;

    for (size_t i = 0; i != live_clusters.size(); i++) {
        Cluster *cluster = live_clusters.at(i);
        // cluster->Create_point_cloud();
        std::pair<geo_point_t, geo_point_t> hl_wcps = cluster->get_highest_lowest_points();
        std::pair<geo_point_t, geo_point_t> fb_wcps = cluster->get_front_back_points();
        std::pair<geo_point_t, geo_point_t> el_wcps = cluster->get_earliest_latest_points();

        // std::cout << cluster->get_cluster_id()  << " " << cluster->get_length() /units::cm << " " <<
        // el_wcps.first.x()/units::cm << " " << el_wcps.second.x()/units::cm << std::endl;

        // if (el_wcps.first.x() < -1 * units::cm || el_wcps.second.x() > 257 * units::cm ||
        if (el_wcps.first.x() < mp.FV_xmin - mp.FV_xmin_margin || el_wcps.second.x() > mp.FV_xmax + mp.FV_xmax_margin || cluster->get_length() < 6.0 * units::cm)
            continue;

        bool flag_fy = false;
        bool flag_by = false;
        bool flag_fx = false;
        bool flag_bx = false;
        bool flag_fz = false;
        bool flag_bz = false;

        std::vector<geo_point_t> saved_wcps;
        if (hl_wcps.first.y() > mp.FV_ymax) {
            saved_wcps.push_back(hl_wcps.first);
            flag_fy = true;
        }

        if (hl_wcps.second.y() < mp.FV_ymin) {
            saved_wcps.push_back(hl_wcps.second);
            flag_by = true;
        }

        if (fb_wcps.first.z() > mp.FV_zmax) {
            bool flag_save = true;
            for (size_t j = 0; j != saved_wcps.size(); j++) {
                double dis = sqrt(pow(saved_wcps.at(j).x() - fb_wcps.first.x(), 2) +
                                  pow(saved_wcps.at(j).y() - fb_wcps.first.y(), 2) +
                                  pow(saved_wcps.at(j).z() - fb_wcps.first.z(), 2));
                if (dis < 15 * units::cm) {
                    flag_save = false;
                    break;
                }
            }
            if (flag_save) {
                saved_wcps.push_back(fb_wcps.first);
                flag_bz = true;
            }
        }

        if (fb_wcps.second.z() < mp.FV_zmin) {
            bool flag_save = true;
            for (size_t j = 0; j != saved_wcps.size(); j++) {
                double dis = sqrt(pow(saved_wcps.at(j).x() - fb_wcps.second.x(), 2) +
                                  pow(saved_wcps.at(j).y() - fb_wcps.second.y(), 2) +
                                  pow(saved_wcps.at(j).z() - fb_wcps.second.z(), 2));
                if (dis < 15 * units::cm) {
                    flag_save = false;
                    break;
                }
            }
            if (flag_save) {
                saved_wcps.push_back(fb_wcps.second);
                flag_fz = true;
            }
        }

        if (el_wcps.first.x() < mp.FV_xmin) {
            bool flag_save = true;
            for (size_t j = 0; j != saved_wcps.size(); j++) {
                double dis = sqrt(pow(saved_wcps.at(j).x() - el_wcps.first.x(), 2) +
                                  pow(saved_wcps.at(j).y() - el_wcps.first.y(), 2) +
                                  pow(saved_wcps.at(j).z() - el_wcps.first.z(), 2));
                if (dis < 15 * units::cm) {
                    flag_save = false;
                    break;
                }
            }
            if (flag_save) {
                saved_wcps.push_back(el_wcps.first);
                flag_fx = true;
            }
        }

        if (el_wcps.second.x() > mp.FV_xmax) {
            bool flag_save = true;
            for (size_t j = 0; j != saved_wcps.size(); j++) {
                double dis = sqrt(pow(saved_wcps.at(j).x() - el_wcps.second.x(), 2) +
                                  pow(saved_wcps.at(j).y() - el_wcps.second.y(), 2) +
                                  pow(saved_wcps.at(j).z() - el_wcps.second.z(), 2));
                if (dis < 15 * units::cm) {
                    flag_save = false;
                    break;
                }
            }
            if (flag_save) {
                saved_wcps.push_back(el_wcps.second);
                flag_bx = true;
            }
        }
        if (saved_wcps.size() <= 1) {
            candidate_clusters.push_back(cluster);
            contained_clusters.push_back(cluster);
        }

        if (saved_wcps.size() >= 2 && (flag_fx && flag_bx || flag_fy && flag_by || flag_fz && flag_bz)) {
        }
        else {
            contained_clusters.push_back(cluster);
        }
    }

    /// TODO: replace with graph? edges between closest clusters, edges are weighted by distance
    std::map<Cluster *, std::pair<Cluster *, double>> cluster_close_cluster_map;
    // calculate the closest distance??? ...
    for (size_t i = 0; i != live_clusters.size(); i++) {
        Cluster *cluster1 = live_clusters.at(i);
        // ToyPointCloud *cloud1 = cluster1->get_point_cloud();
        for (size_t j = i + 1; j != live_clusters.size(); j++) {
            Cluster *cluster2 = live_clusters.at(j);
            // ToyPointCloud *cloud2 = cluster2->get_point_cloud();

            // std::tuple<int, int, double> results = cloud2->get_closest_points(cloud1);
            std::tuple<int, int, double> results = cluster2->get_closest_points(*cluster1);
            double dis = std::get<2>(results);

            if (cluster_close_cluster_map.find(cluster1) == cluster_close_cluster_map.end()) {
                cluster_close_cluster_map[cluster1] = std::make_pair(cluster2, dis);
            }
            else {
                if (dis < cluster_close_cluster_map[cluster1].second)
                    cluster_close_cluster_map[cluster1] = std::make_pair(cluster2, dis);
            }

            if (cluster_close_cluster_map.find(cluster2) == cluster_close_cluster_map.end()) {
                cluster_close_cluster_map[cluster2] = std::make_pair(cluster1, dis);
            }
            else {
                if (dis < cluster_close_cluster_map[cluster2].second)
                    cluster_close_cluster_map[cluster2] = std::make_pair(cluster1, dis);
            }
        }
    }

    //  std::cout << contained_clusters.size() << " " << candidate_clusters.size() << std::endl;

    std::set<std::pair<Cluster *, Cluster *>> to_be_merged_pairs;

    std::set<Cluster *> used_clusters;

    std::map<Cluster *, std::shared_ptr<Simple3DPointCloud>> cluster_cloud_map;

    std::map<Cluster *, geo_point_t> cluster_dir1_map;
    std::map<Cluster *, geo_point_t> cluster_dir2_map;

    // for (auto it = candidate_clusters.begin(); it != candidate_clusters.end(); it++) {
    //     Cluster *cluster1 = (*it);
    //     cluster1->Calc_PCA();
    // }

    // ignore very small ones?
    // two short ones, NC pi0 case
    // one short one and one big one, CC pi0
    for (auto it = candidate_clusters.begin(); it != candidate_clusters.end(); it++) {
        Cluster *cluster1 = (*it);
        // cluster1->Create_point_cloud();
        // ToyPointCloud *cloud1 = cluster1->get_point_cloud();
        for (auto it1 = contained_clusters.begin(); it1 != contained_clusters.end(); it1++) {
            Cluster *cluster2 = (*it1);
            // can not be the same
            if (cluster2 == cluster1) continue;
            if (cluster2->get_length() > 150 * units::cm) {
                geo_point_t dir1(cluster2->get_pca_axis(0).x(), cluster2->get_pca_axis(0).y(),
                                 cluster2->get_pca_axis(0).z());
                if (fabs(dir1.angle(vertical_dir) - 3.1415926 / 2.) / 3.1415926 * 180. > 80) continue;
            }
            // cluster2->Create_point_cloud();
            // ToyPointCloud *cloud2 = cluster2->get_point_cloud();

            // std::tuple<int, int, double> results = cloud2->get_closest_points(cloud1);
            std::tuple<int, int, double> results = cluster2->get_closest_points(*cluster1);
            double dis = std::get<2>(results);
            if ((dis > 80 * units::cm || cluster1->get_length() > 80 * units::cm && dis > 10 * units::cm)) continue;

            if (cluster_cloud_map.find(cluster1) == cluster_cloud_map.end()) {
                // cluster1->Calc_PCA();
                geo_point_t center = cluster1->get_center();
                geo_point_t main_dir(cluster1->get_pca_axis(0).x(), cluster1->get_pca_axis(0).y(),
                                     cluster1->get_pca_axis(0).z());
                main_dir = main_dir.norm();

                // ToyPointCloud *cloud1_ext = new ToyPointCloud(angle_u, angle_v, angle_w);
                auto cloud1_ext = std::make_shared<Simple3DPointCloud>();
                cluster_cloud_map[cluster1] = cloud1_ext;
                // WCP::PointVector pts;
                std::vector<geo_point_t> pts;
                std::pair<geo_point_t, geo_point_t> extreme_pts = cluster1->get_two_extreme_points();

                if (cluster1->get_length() > 25 * units::cm) {
                    extreme_pts.first = cluster1->calc_ave_pos(extreme_pts.first, 5 * units::cm);
                    extreme_pts.second = cluster1->calc_ave_pos(extreme_pts.second, 5 * units::cm);
                }

                geo_point_t dir1 = cluster1->vhough_transform(extreme_pts.first, 30 * units::cm);
                geo_point_t dir2 = cluster1->vhough_transform(extreme_pts.second, 30 * units::cm);

                cluster_dir1_map[cluster1] = dir1;
                cluster_dir2_map[cluster1] = dir2;

                bool flag_enable_temp = false;
                std::pair<geo_point_t, geo_point_t> temp_extreme_pts;
                geo_point_t temp_dir1;
                geo_point_t temp_dir2;
                int num_clusters = 0;

                if (cluster1->nnearby(extreme_pts.first, 15 * units::cm) <= 75 && cluster1->npoints() > 75 ||
                    cluster1->nnearby(extreme_pts.second, 15 * units::cm) <= 75 && cluster1->npoints() > 75 ||
                    cluster1->get_pca_value(1) > 0.022 * cluster1->get_pca_value(0) &&
                        cluster1->get_length() > 45 * units::cm) {
                    // std::vector<Cluster *> sep_clusters = Separate_2(cluster1, 2.5 * units::cm);
                    const double orig_cluster_length = cluster1->get_length();
                    std::cout  << "[neutrino] cluster1->npoints() " << cluster1->npoints() << " " << cluster1->point(0) << std::endl;
                    const auto b2id = Separate_2(cluster1, 2.5 * units::cm);
                    // false: do not remove the cluster1
                    auto sep_clusters = live_grouping.separate(cluster1, b2id, false);
                    assert(cluster1 != nullptr);
                    Cluster *largest_cluster = 0;
                    int max_num_points = 0;
                    for (auto [id, sep_cluster] : sep_clusters) {
                        if (sep_cluster->npoints() > max_num_points) {
                            max_num_points = sep_cluster->npoints();
                            largest_cluster = sep_cluster;
                        }
                    }

                    temp_extreme_pts = largest_cluster->get_two_extreme_points();
                    center = largest_cluster->get_center();
                    main_dir.set(largest_cluster->get_pca_axis(0).x(), largest_cluster->get_pca_axis(0).y(),
                                 largest_cluster->get_pca_axis(0).z());
                    num_clusters = sep_clusters.size();
                    // largest_cluster->Create_point_cloud();

                    if (orig_cluster_length > 25 * units::cm) {
                        temp_extreme_pts.first = largest_cluster->calc_ave_pos(temp_extreme_pts.first, 5 * units::cm);
                        temp_extreme_pts.second = largest_cluster->calc_ave_pos(temp_extreme_pts.second, 5 * units::cm);
                    }

                    flag_enable_temp = true;
                    temp_dir1 = largest_cluster->vhough_transform(temp_extreme_pts.first, 30 * units::cm);
                    temp_dir2 = largest_cluster->vhough_transform(temp_extreme_pts.second, 30 * units::cm);

                    // for (size_t j = 0; j != sep_clusters.size(); j++) {
                    //     delete sep_clusters.at(j);
                    // }
                    // merge back ...
                    // cluster1 = &(live_grouping.make_child());
                    for (size_t j = 0; j != sep_clusters.size(); j++) {
                        cluster1->take_children(*sep_clusters.at(j), true);
                        live_grouping.destroy_child(sep_clusters.at(j));
                        assert(sep_clusters.at(j) == nullptr);
                    }
                    std::cout  << "[neutrino] cluster1->npoints() " << cluster1->npoints() << " " << cluster1->point(0) << std::endl;
                }

                dir1 *= -1;
                dir2 *= -1;
                dir1 = dir1.norm();
                dir2 = dir2.norm();
                if (flag_enable_temp) {
                    temp_dir1 *= -1;
                    temp_dir2 *= -1;
                    temp_dir1 = temp_dir1.norm();
                    temp_dir2 = temp_dir2.norm();
                }

                bool flag_add1 = true;
                if (cluster1->nnearby(extreme_pts.first, 15 * units::cm) <= 75 &&
                        cluster1->get_length() > 60 * units::cm ||
                    flag_enable_temp && num_clusters >= 4 &&
                        cluster1->get_pca_value(1) > 0.022 * cluster1->get_pca_value(0))
                    flag_add1 = false;

                bool flag_add2 = true;
                if (cluster1->nnearby(extreme_pts.second, 15 * units::cm) <= 75 &&
                        cluster1->get_length() > 60 * units::cm ||
                    flag_enable_temp && num_clusters >= 4 &&
                        cluster1->get_pca_value(1) > 0.022 * cluster1->get_pca_value(0))
                    flag_add2 = false;

                for (size_t j = 0; j != 150; j++) {
                    if (flag_add1) {
                        geo_point_t pt1(extreme_pts.first.x() + dir1.x() * (j + 1) * 0.5 * units::cm,
                                        extreme_pts.first.y() + dir1.y() * (j + 1) * 0.5 * units::cm,
                                        extreme_pts.first.z() + dir1.z() * (j + 1) * 0.5 * units::cm);
                        pts.push_back(pt1);
                    }
                    if (flag_add2) {
                        geo_point_t pt2(extreme_pts.second.x() + dir2.x() * (j + 1) * 0.5 * units::cm,
                                        extreme_pts.second.y() + dir2.y() * (j + 1) * 0.5 * units::cm,
                                        extreme_pts.second.z() + dir2.z() * (j + 1) * 0.5 * units::cm);
                        pts.push_back(pt2);
                    }
                    if (flag_enable_temp) {
                        geo_point_t pt1(temp_extreme_pts.first.x() + temp_dir1.x() * (j + 1) * 0.5 * units::cm,
                                        temp_extreme_pts.first.y() + temp_dir1.y() * (j + 1) * 0.5 * units::cm,
                                        temp_extreme_pts.first.z() + temp_dir1.z() * (j + 1) * 0.5 * units::cm);
                        pts.push_back(pt1);
                        geo_point_t pt2(temp_extreme_pts.second.x() + temp_dir2.x() * (j + 1) * 0.5 * units::cm,
                                        temp_extreme_pts.second.y() + temp_dir2.y() * (j + 1) * 0.5 * units::cm,
                                        temp_extreme_pts.second.z() + temp_dir2.z() * (j + 1) * 0.5 * units::cm);
                        pts.push_back(pt2);
                    }

                    if ((!flag_add1) && (!flag_add2) && (!flag_enable_temp)) {
                        pts.push_back(extreme_pts.first);
                        pts.push_back(extreme_pts.second);
                    }

                    if (cluster1->get_length() < 60 * units::cm) {
                        geo_point_t temp1(extreme_pts.first.x() - center.x(), extreme_pts.first.y() - center.y(),
                                          extreme_pts.first.z() - center.z());
                        double length1 = temp1.dot(main_dir);
                        if (length1 > 0) {
                            geo_point_t pt3(center.x() + main_dir.x() * (length1 + (j + 1) * 0.5 * units::cm),
                                            center.y() + main_dir.y() * (length1 + (j + 1) * 0.5 * units::cm),
                                            center.z() + main_dir.z() * (length1 + (j + 1) * 0.5 * units::cm));
                            pts.push_back(pt3);
                        }
                        else {
                            geo_point_t pt3(center.x() + main_dir.x() * (length1 - (j + 1) * 0.5 * units::cm),
                                            center.y() + main_dir.y() * (length1 - (j + 1) * 0.5 * units::cm),
                                            center.z() + main_dir.z() * (length1 - (j + 1) * 0.5 * units::cm));
                            pts.push_back(pt3);
                        }

                        geo_point_t temp2(extreme_pts.second.x() - center.x(), extreme_pts.second.y() - center.y(),
                                          extreme_pts.second.z() - center.z());
                        double length2 = temp2.dot(main_dir);

                        if (length2 > 0) {
                            geo_point_t pt4(center.x() + main_dir.x() * (length2 + (j + 1) * 0.5 * units::cm),
                                            center.y() + main_dir.y() * (length2 + (j + 1) * 0.5 * units::cm),
                                            center.z() + main_dir.z() * (length2 + (j + 1) * 0.5 * units::cm));
                            pts.push_back(pt4);
                        }
                        else {
                            geo_point_t pt4(center.x() + main_dir.x() * (length2 - (j + 1) * 0.5 * units::cm),
                                            center.y() + main_dir.y() * (length2 - (j + 1) * 0.5 * units::cm),
                                            center.z() + main_dir.z() * (length2 - (j + 1) * 0.5 * units::cm));
                            pts.push_back(pt4);
                        }
                    }
                }
                // cloud1_ext->AddPoints(pts);
                for (size_t j = 0; j != pts.size(); j++) {
                    cloud1_ext->add({pts.at(j).x(), pts.at(j).y(), pts.at(j).z()});
                }
                // cloud1_ext->build_kdtree_index();
            }

            if (cluster_cloud_map.find(cluster2) == cluster_cloud_map.end()) {
                // cluster2->Calc_PCA();
                geo_point_t center = cluster2->get_center();
                geo_point_t main_dir(cluster2->get_pca_axis(0).x(), cluster2->get_pca_axis(0).y(),
                                     cluster2->get_pca_axis(0).z());
                main_dir = main_dir.norm();

                // ToyPointCloud *cloud2_ext = new ToyPointCloud(angle_u, angle_v, angle_w);
                auto cloud2_ext = std::make_shared<Simple3DPointCloud>();
                cluster_cloud_map[cluster2] = cloud2_ext;
                // WCP::PointVector pts;
                std::vector<geo_point_t> pts;
                std::pair<geo_point_t, geo_point_t> extreme_pts = cluster2->get_two_extreme_points();

                if (cluster2->get_length() > 25 * units::cm) {
                    extreme_pts.first = cluster2->calc_ave_pos(extreme_pts.first, 5 * units::cm);
                    extreme_pts.second = cluster2->calc_ave_pos(extreme_pts.second, 5 * units::cm);
                }

                geo_point_t dir1 = cluster2->vhough_transform(extreme_pts.first, 30 * units::cm);
                geo_point_t dir2 = cluster2->vhough_transform(extreme_pts.second, 30 * units::cm);

                cluster_dir1_map[cluster2] = dir1;
                cluster_dir2_map[cluster2] = dir2;

                bool flag_enable_temp = false;
                std::pair<geo_point_t, geo_point_t> temp_extreme_pts;
                geo_point_t temp_dir1;
                geo_point_t temp_dir2;
                int num_clusters = 0;

                if (cluster2->nnearby(extreme_pts.first, 15 * units::cm) <= 75 && cluster2->npoints() > 75 ||
                    cluster2->nnearby(extreme_pts.second, 15 * units::cm) <= 75 && cluster2->npoints() > 75 ||
                    cluster2->get_pca_value(1) > 0.022 * cluster2->get_pca_value(0) &&
                        cluster2->get_length() > 45 * units::cm) {
                    // std::vector<Cluster *> sep_clusters = Separate_2(cluster2, 2.5 * units::cm);
                    std::cout  << "[neutrino] cluster2->npoints() " << cluster2->npoints() << " " << cluster2->point(0) << std::endl;
                    const double orig_cluster_length = cluster2->get_length();
                    const auto b2id = Separate_2(cluster2, 2.5 * units::cm);
                    auto sep_clusters = live_grouping.separate(cluster2, b2id, false);
                    assert(cluster2 != nullptr);
                    Cluster *largest_cluster = 0;
                    int max_num_points = 0;
                    for (auto [id, sep_cluster] : sep_clusters) {
                        if (sep_cluster->npoints() > max_num_points) {
                            max_num_points = sep_cluster->npoints();
                            largest_cluster = sep_cluster;
                        }
                    }
                    temp_extreme_pts = largest_cluster->get_two_extreme_points();
                    center = largest_cluster->get_center();
                    main_dir.set(largest_cluster->get_pca_axis(0).x(), largest_cluster->get_pca_axis(0).y(),
                                 largest_cluster->get_pca_axis(0).z());
                    num_clusters = sep_clusters.size();

                    // largest_cluster->Create_point_cloud();

                    if (orig_cluster_length > 25 * units::cm) {
                        temp_extreme_pts.first = largest_cluster->calc_ave_pos(temp_extreme_pts.first, 5 * units::cm);
                        temp_extreme_pts.second = largest_cluster->calc_ave_pos(temp_extreme_pts.second, 5 * units::cm);
                    }

                    flag_enable_temp = true;
                    temp_dir1 = largest_cluster->vhough_transform(temp_extreme_pts.first, 30 * units::cm);
                    temp_dir2 = largest_cluster->vhough_transform(temp_extreme_pts.second, 30 * units::cm);

                    // for (size_t j = 0; j != sep_clusters.size(); j++) {
                    //     delete sep_clusters.at(j);
                    // }
                    // merge back ...
                    // cluster2 = &(live_grouping.make_child());
                    for (size_t j = 0; j != sep_clusters.size(); j++) {
                        cluster2->take_children(*sep_clusters.at(j), true);
                        live_grouping.destroy_child(sep_clusters.at(j));
                        assert(sep_clusters.at(j) == nullptr);
                    }
                    std::cout  << "[neutrino] cluster2->npoints() " << cluster2->npoints() << " " << cluster2->point(0) << std::endl;
                }

                dir1 *= -1;
                dir2 *= -1;
                dir1 = dir1.norm();
                dir2 = dir2.norm();
                if (flag_enable_temp) {
                    temp_dir1 *= -1;
                    temp_dir2 *= -1;
                    temp_dir1 = temp_dir1.norm();
                    temp_dir2 = temp_dir2.norm();
                }

                bool flag_add1 = true;
                if (cluster2->nnearby(extreme_pts.first, 15 * units::cm) <= 75 &&
                        cluster2->get_length() > 60 * units::cm ||
                    flag_enable_temp && num_clusters >= 4 &&
                        cluster2->get_pca_value(1) > 0.022 * cluster2->get_pca_value(0))
                    flag_add1 = false;
                bool flag_add2 = true;
                if (cluster2->nnearby(extreme_pts.second, 15 * units::cm) <= 75 &&
                        cluster2->get_length() > 60 * units::cm ||
                    flag_enable_temp && num_clusters >= 4 &&
                        cluster2->get_pca_value(1) > 0.022 * cluster2->get_pca_value(0))
                    flag_add2 = false;

                // std::cout << flag_add1 << " " << flag_add2 << " " << dir1.x() << " " << dir1.y() << " " <<
                // dir1.z()
                // << " " << dir2.x() << " " << dir2.y() << " " << dir2.z() << std::endl;

                for (size_t j = 0; j != 150; j++) {
                    if (flag_add1) {
                        geo_point_t pt1(extreme_pts.first.x() + dir1.x() * (j + 1) * 0.5 * units::cm,
                                        extreme_pts.first.y() + dir1.y() * (j + 1) * 0.5 * units::cm,
                                        extreme_pts.first.z() + dir1.z() * (j + 1) * 0.5 * units::cm);
                        pts.push_back(pt1);
                    }
                    if (flag_add2) {
                        geo_point_t pt2(extreme_pts.second.x() + dir2.x() * (j + 1) * 0.5 * units::cm,
                                        extreme_pts.second.y() + dir2.y() * (j + 1) * 0.5 * units::cm,
                                        extreme_pts.second.z() + dir2.z() * (j + 1) * 0.5 * units::cm);
                        pts.push_back(pt2);
                    }

                    if (flag_enable_temp) {
                        geo_point_t pt1(temp_extreme_pts.first.x() + temp_dir1.x() * (j + 1) * 0.5 * units::cm,
                                        temp_extreme_pts.first.y() + temp_dir1.y() * (j + 1) * 0.5 * units::cm,
                                        temp_extreme_pts.first.z() + temp_dir1.z() * (j + 1) * 0.5 * units::cm);
                        pts.push_back(pt1);
                        geo_point_t pt2(temp_extreme_pts.second.x() + temp_dir2.x() * (j + 1) * 0.5 * units::cm,
                                        temp_extreme_pts.second.y() + temp_dir2.y() * (j + 1) * 0.5 * units::cm,
                                        temp_extreme_pts.second.z() + temp_dir2.z() * (j + 1) * 0.5 * units::cm);
                        pts.push_back(pt2);
                    }

                    if ((!flag_add1) && (!flag_add2)) {
                        pts.push_back(extreme_pts.first);
                        pts.push_back(extreme_pts.second);
                    }

                    if (cluster2->get_length() < 60 * units::cm) {
                        geo_point_t temp1(extreme_pts.first.x() - center.x(), extreme_pts.first.y() - center.y(),
                                          extreme_pts.first.z() - center.z());
                        double length1 = temp1.dot(main_dir);
                        if (length1 > 0) {
                            geo_point_t pt3(center.x() + main_dir.x() * (length1 + (j + 1) * 0.5 * units::cm),
                                            center.y() + main_dir.y() * (length1 + (j + 1) * 0.5 * units::cm),
                                            center.z() + main_dir.z() * (length1 + (j + 1) * 0.5 * units::cm));
                            pts.push_back(pt3);
                        }
                        else {
                            geo_point_t pt3(center.x() + main_dir.x() * (length1 - (j + 1) * 0.5 * units::cm),
                                            center.y() + main_dir.y() * (length1 - (j + 1) * 0.5 * units::cm),
                                            center.z() + main_dir.z() * (length1 - (j + 1) * 0.5 * units::cm));
                            pts.push_back(pt3);
                        }

                        geo_point_t temp2(extreme_pts.second.x() - center.x(), extreme_pts.second.y() - center.y(),
                                          extreme_pts.second.z() - center.z());
                        double length2 = temp2.dot(main_dir);

                        if (length2 > 0) {
                            geo_point_t pt4(center.x() + main_dir.x() * (length2 + (j + 1) * 0.5 * units::cm),
                                            center.y() + main_dir.y() * (length2 + (j + 1) * 0.5 * units::cm),
                                            center.z() + main_dir.z() * (length2 + (j + 1) * 0.5 * units::cm));
                            pts.push_back(pt4);
                        }
                        else {
                            geo_point_t pt4(center.x() + main_dir.x() * (length2 - (j + 1) * 0.5 * units::cm),
                                            center.y() + main_dir.y() * (length2 - (j + 1) * 0.5 * units::cm),
                                            center.z() + main_dir.z() * (length2 - (j + 1) * 0.5 * units::cm));
                            pts.push_back(pt4);
                        }
                    }
                }
                // cloud2_ext->AddPoints(pts);
                for (size_t j = 0; j != pts.size(); j++) {
                    // cloud2_ext->add(pts.at(j));
                    cloud2_ext->add({pts.at(j).x(), pts.at(j).y(), pts.at(j).z()});
                }
                // cloud2_ext->build_kdtree_index();
            }

            // ToyPointCloud *cloud1_ext = cluster_cloud_map[cluster1];
            // ToyPointCloud *cloud2_ext = cluster_cloud_map[cluster2];
            auto cloud1_ext = cluster_cloud_map[cluster1];
            auto cloud2_ext = cluster_cloud_map[cluster2];

            int merge_type = 0;
            bool flag_merge = false;
            {
                std::tuple<int, int, double> results_1 = cloud1_ext->get_closest_points(*cluster2);
                geo_point_t test_pt = cloud1_ext->point(std::get<0>(results_1));
                geo_point_t test_pt1 = cluster2->point3d(std::get<1>(results_1));
                double dis1 = std::get<2>(results_1);
                double dis2 = cluster1->get_closest_dis(test_pt);

                if (dis1 < std::min(std::max(4.5 * units::cm, dis2 * sin(15 / 180. * 3.1415926)), 12 * units::cm) &&
                        (cluster2->get_length() > 25 * units::cm || cluster1->get_length() <= cluster2->get_length()) ||
                    dis1 < std::min(std::max(2.5 * units::cm, dis2 * sin(10 / 180. * 3.1415926)), 10 * units::cm) ||
                    dis1 < std::min(std::max(4.5 * units::cm, dis2 * sin(25 / 180. * 3.1415926)), 12 * units::cm) &&
                        dis < 30 * units::cm && dis2 < 30 * units::cm && cluster1->get_length() > 15 * units::cm &&
                        cluster2->get_length() > 15 * units::cm ||
                    cluster1->get_length() > 45 * units::cm && dis1 < 16 * units::cm &&
                        fabs(test_pt.x() - test_pt1.x()) < 3.2 * units::cm &&
                        (fabs(drift_dir.angle(cluster_dir1_map[cluster1]) - 3.1415926 / 2.) / 3.1415926 * 180. < 5 ||
                         fabs(drift_dir.angle(cluster_dir2_map[cluster1]) - 3.1415926 / 2.) / 3.1415926 * 180. < 5)) {
                    // std::cout << test_pt1.x()/units::cm << " " << test_pt1.y()/units::cm << " " <<
                    // test_pt1.z()/units::cm
                    // << std::endl;

                    {  // std::tuple<int,int,double> results =  cloud2->get_closest_points(cloud1);
                        // geo_point_t test_pt2(cloud2->get_cloud().pts.at(std::get<0>(results)).x(),
                        //                cloud2->get_cloud().pts.at(std::get<0>(results)).y(),
                        //                cloud2->get_cloud().pts.at(std::get<0>(results)).z());
                        // geo_point_t test_pt3(cloud1->get_cloud().pts.at(std::get<1>(results)).x(),
                        //                cloud1->get_cloud().pts.at(std::get<1>(results)).y(),
                        //                cloud1->get_cloud().pts.at(std::get<1>(results)).z());
                        geo_point_t test_pt2 = cluster2->point3d(std::get<0>(results));
                        geo_point_t test_pt3 = cluster1->point3d(std::get<1>(results));
                        test_pt3 = cluster1->calc_ave_pos(test_pt3, 5 * units::cm);
                        geo_point_t temp_dir(test_pt2.x() - test_pt3.x(), test_pt2.y() - test_pt3.y(),
                                             test_pt2.z() - test_pt3.z());
                        double angle_diff1 =
                            fabs(temp_dir.angle(cluster_dir1_map[cluster1]) - 3.1415926 / 2.) / 3.1415926 * 180.;
                        double angle_diff2 =
                            fabs(temp_dir.angle(cluster_dir2_map[cluster1]) - 3.1415926 / 2.) / 3.1415926 * 180.;
                        if ((angle_diff1 > 65 || angle_diff2 > 65) &&
                            (dis * sin((90 - angle_diff1) / 180. * 3.1415926) < 4.5 * units::cm ||
                             dis * sin((90 - angle_diff2) / 180. * 3.1415926) < 4.5 * units::cm)) {
                            if (!cluster2->judge_vertex(test_pt1)) {
                                test_pt1 = test_pt2;
                            }
                        }
                    }

                    if (cluster1->get_length() > 25 * units::cm &&
                        cluster1->get_pca_value(1) < 0.0015 * cluster1->get_pca_value(0)) {
                        flag_merge = false;

                        if (dis < 0.5 * units::cm && dis1 < 1.5 * units::cm && dis2 < 1.5 * units::cm)
                            flag_merge = cluster2->judge_vertex(test_pt1, 0.5, 0.6);
                    }
                    else {
                        if (cluster2->get_length() < 30 * units::cm) {
                            flag_merge = true;
                            if (cluster1->get_length() > 15 * units::cm &&
                                cluster1->get_pca_value(1) < 0.012 * cluster1->get_pca_value(0)) {
                                if (dis1 > std::max(2.5 * units::cm, dis2 * sin(7.5 / 180. * 3.1415926)))
                                    flag_merge = false;
                            }

                            if (cluster1->get_length() > 150 * units::cm) {
                                flag_merge = false;
                            }
                        }
                        else if (JudgeSeparateDec_1(cluster2, drift_dir, cluster2->get_length(), time_slice_width)) {
                            if (dis2 < 5 * units::cm) {
                                flag_merge = cluster2->judge_vertex(test_pt1, 2. / 3.);
                            }
                            else if (dis < 0.5 * units::cm) {
                                flag_merge = cluster2->judge_vertex(test_pt1, 0.5, 0.6);
                            }
                            else {
                                flag_merge = cluster2->judge_vertex(test_pt1);
                            }

                            if (cluster1->get_length() > 15 * units::cm &&
                                cluster1->get_pca_value(1) < 0.012 * cluster1->get_pca_value(0)) {
                                if (dis1 > std::max(2.5 * units::cm, dis2 * sin(7.5 / 180. * 3.1415926)))
                                    flag_merge = false;
                            }
                        }
                        else {
                            if (dis2 < 5 * units::cm) {
                                flag_merge = cluster2->judge_vertex(test_pt1, 2. / 3.);
                            }
                            else if (dis < 0.5 * units::cm) {
                                flag_merge = cluster2->judge_vertex(test_pt1, 0.5, 0.6);
                            }
                            else {
                                flag_merge = cluster2->judge_vertex(test_pt1);
                            }

                            if (cluster1->get_length() > 15 * units::cm &&
                                cluster1->get_pca_value(1) < 0.012 * cluster1->get_pca_value(0)) {
                                if (dis1 > std::max(3.5 * units::cm, dis2 * sin(7.5 / 180. * 3.1415926)))
                                    flag_merge = false;
                            }

                            if (flag_merge && cluster2->get_length() > 200 * units::cm && dis2 < 12 * units::cm &&
                                cluster2->get_pca_value(1) < 0.0015 * cluster2->get_pca_value(0)) {
                                geo_point_t cluster2_dir(cluster2->get_pca_axis(0).x(), cluster2->get_pca_axis(0).y(),
                                                         cluster2->get_pca_axis(0).z());
                                if (fabs(cluster2_dir.angle(vertical_dir) / 3.1415926 * 180. - 3.1415926 / 2.) /
                                            3.1415926 * 180. >
                                        45 &&
                                    fabs(cluster2_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 20)
                                    flag_merge = false;
                            }
                        }
                        merge_type = 1;
                    }
                    //

                    if (cluster_close_cluster_map[cluster1].second < 1.2 * units::cm &&
                        cluster_close_cluster_map[cluster1].first != cluster2 &&
                        cluster_close_cluster_map[cluster1].first->get_length() > 60 * units::cm &&
                        cluster1->get_pca_value(1) > 0.012 * cluster1->get_pca_value(0) && dis1 > 0.6 * units::cm) {
                        flag_merge = false;
                    }

                    if (test_pt1.y() > 112 * units::cm && dis < 5 * units::cm && dis1 < 3.0 * units::cm &&
                        cluster1->get_length() > 60 * units::cm && cluster2->get_length() > 60 * units::cm)
                        flag_merge = false;

                    if (flag_merge && cluster1->get_length() > 150 * units::cm &&
                        cluster2->get_length() > 150 * units::cm &&
                        (cluster1->get_pca_value(1) < 0.03 * cluster1->get_pca_value(0) ||
                         cluster2->get_pca_value(1) < 0.03 * cluster2->get_pca_value(0))) {
                        // protect against two long tracks ...
                        // cluster1->Calc_PCA();
                        // cluster2->Calc_PCA();
                        geo_point_t temp_dir1(cluster1->get_pca_axis(0).x(), cluster1->get_pca_axis(0).y(),
                                              cluster1->get_pca_axis(0).z());
                        geo_point_t temp_dir2(cluster2->get_pca_axis(0).x(), cluster2->get_pca_axis(0).y(),
                                              cluster2->get_pca_axis(0).z());
                        if (fabs(temp_dir1.angle(temp_dir2) - 3.1415926 / 2.) < 60 / 180. * 3.1415926)
                            flag_merge = false;
                    }
                }

                if (dis < 1.8 * units::cm && cluster1->get_length() < 75 * units::cm &&
                    cluster2->get_length() < 75 * units::cm &&
                    (cluster1->get_length() + cluster2->get_length()) < 120 * units::cm) {
                    flag_merge = true;
                    merge_type = 2;
                }
            }

            if (!flag_merge) {
                std::tuple<int, int, double> results_2 = cloud1_ext->get_closest_points(*cloud2_ext);
                // geo_point_t test_pt(cloud1_ext->get_cloud().pts.at(std::get<0>(results_2)).x(),
                //                     cloud1_ext->get_cloud().pts.at(std::get<0>(results_2)).y(),
                //                     cloud1_ext->get_cloud().pts.at(std::get<0>(results_2)).z());
                geo_point_t test_pt = cloud1_ext->point(std::get<0>(results_2));
                // geo_point_t test_pt1(cloud2_ext->get_cloud().pts.at(std::get<1>(results_2)).x(),
                //                      cloud2_ext->get_cloud().pts.at(std::get<1>(results_2)).y(),
                //                      cloud2_ext->get_cloud().pts.at(std::get<1>(results_2)).z());
                geo_point_t test_pt1 = cloud2_ext->point(std::get<1>(results_2));
                double dis1 = std::get<2>(results_2);
                // double dis2 = cloud1->get_closest_dis(test_pt);
                double dis2 = cluster1->get_closest_dis(test_pt);
                // double dis3 = cloud2->get_closest_dis(test_pt1);
                double dis3 = cluster2->get_closest_dis(test_pt1);
                if (dis1 < std::min(std::max(4.5 * units::cm, (dis2 + dis3) / 2. * sin(15 / 180. * 3.1415926)),
                                    12 * units::cm) &&
                    dis1 < std::min(std::max(4.5 * units::cm, (dis3 + dis2) / 2. * sin(15 / 180. * 3.1415926)),
                                    12 * units::cm) &&
                    dis2 + dis3 < 72 * units::cm && cluster2->get_length() < 60 * units::cm &&
                    cluster1->get_length() < 60 * units::cm) {
                    flag_merge = true;
                    merge_type = 3;
                }
                else if (dis2 + dis3 < 90 * units::cm && dis1 < 2.7 * units::cm && dis < 20 * units::cm &&
                         cluster2->get_length() > 30 * units::cm && cluster1->get_length() > 30 * units::cm) {
                    // cluster1->Calc_PCA();
                    // cluster2->Calc_PCA();
                    if (cluster1->get_pca_value(1) > 0.0015 * cluster1->get_pca_value(0) &&
                        cluster2->get_pca_value(1) > 0.0015 * cluster2->get_pca_value(0)) {
                        flag_merge = true;
                        merge_type = 3;
                    }
                }

                if (test_pt1.y() > 112 * units::cm && dis2 < 5 * units::cm && dis1 < 3.0 * units::cm &&
                    dis3 < 5 * units::cm && cluster1->get_length() > 60 * units::cm &&
                    cluster2->get_length() > 60 * units::cm)
                    flag_merge = false;
            }

            if (flag_merge) {
                bool flag_proceed = true;
                if (merge_type == 1) {
                    if (used_clusters.find(cluster1) != used_clusters.end()) flag_proceed = false;
                }
                else if (merge_type == 3) {
                    if (used_clusters.find(cluster1) != used_clusters.end() &&
                        used_clusters.find(cluster2) != used_clusters.end())
                        flag_proceed = false;
                }
                else if (merge_type == 2) {
                    if (used_clusters.find(cluster1) != used_clusters.end() ||
                        used_clusters.find(cluster2) != used_clusters.end())
                        flag_proceed = false;
                }

                if (flag_proceed) {
                    to_be_merged_pairs.insert(std::make_pair(cluster1, cluster2));
                    if (merge_type == 1) {
                        used_clusters.insert(cluster1);
                    }
                    else if (merge_type == 3) {
                        used_clusters.insert(cluster1);
                        used_clusters.insert(cluster2);
                    }
                }
            }
        }
    }

    // for (auto it = cluster_cloud_map.begin(); it != cluster_cloud_map.end(); it++) {
    //     delete it->second;
    // }

    // prepare a graph ...
    typedef cluster_connectivity_graph_t Graph;
    Graph g;
    std::unordered_map<int, int> ilive2desc;  // added live index to graph descriptor
    std::map<const Cluster*, int> map_cluster_index;
    for (const Cluster* live : live_grouping.children()) {
        size_t ilive = map_cluster_index.size();
        map_cluster_index[live] = ilive;
        ilive2desc[ilive] = boost::add_vertex(ilive, g);
    }
    for (auto [cluster1, cluster2] : to_be_merged_pairs) {
        boost::add_edge(ilive2desc[map_cluster_index[cluster1]],
                        ilive2desc[map_cluster_index[cluster2]], g);
    }
    cluster_set_t new_clusters;
    merge_clusters(g, live_grouping, new_clusters);

    // merge clusters
    // std::vector<std::set<Cluster *>> merge_clusters;
    // for (auto it = to_be_merged_pairs.begin(); it != to_be_merged_pairs.end(); it++) {
    //     Cluster *cluster1 = (*it).first;
    //     Cluster *cluster2 = (*it).second;
    //     bool flag_new = true;
    //     std::vector<std::set<Cluster *>> temp_set;
    //     for (auto it1 = merge_clusters.begin(); it1 != merge_clusters.end(); it1++) {
    //         std::set<Cluster *> &clusters = (*it1);
    //         if (clusters.find(cluster1) != clusters.end() || clusters.find(cluster2) != clusters.end()) {
    //             clusters.insert(cluster1);
    //             clusters.insert(cluster2);
    //             flag_new = false;
    //             temp_set.push_back(clusters);
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

    // merge clusters into new clusters, delete old clusters
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
    // }
}
