#include <WireCellClus/ClusteringFuncs.h>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;


map_cluster_cluster_vec WireCell::PointCloud::Facade::clustering_isolated(Grouping& live_grouping)
{
    std::vector<Cluster *> live_clusters = live_grouping.children();  // copy
    const auto &mp = live_grouping.get_params();
    // this is for 4 time slices
    double time_slice_width = mp.nticks_live_slice * mp.tick_drift;
    geo_point_t drift_dir(1, 0, 0);

    int range_cut = 150;
    int length_cut = 20 * units::cm;

    std::vector<Cluster*> big_clusters;
    std::vector<Cluster*> small_clusters;

    for (size_t i = 0; i != live_clusters.size(); i++) {
        std::tuple<int, int, int, int> ranges_tuple = live_clusters.at(i)->get_uvwt_range();
        std::vector<int> ranges = {std::get<0>(ranges_tuple), std::get<1>(ranges_tuple), std::get<2>(ranges_tuple), std::get<3>(ranges_tuple)};
        int max = 0;
        for (int j = 0; j != 4; j++) {
            if (ranges.at(j) > max) max = ranges.at(j);
        }
        if (max < range_cut && live_clusters.at(i)->get_length() < length_cut) {
            small_clusters.push_back(live_clusters.at(i));
        }
        else {
            if (live_clusters.at(i)->get_length() < 60 * units::cm) {
                if (JudgeSeparateDec_1(live_clusters.at(i), drift_dir, live_clusters.at(i)->get_length(),
                                       time_slice_width)) {
                    std::vector<Cluster *> sep_clusters = Separate_2(live_clusters.at(i), 2.5 * units::cm);
                    int max = 0;
                    double max_length = 0;
                    for (auto it = sep_clusters.begin(); it != sep_clusters.end(); it++) {
                        // std::tuple<int, int, int, int> ranges = (*it)->get_uvwt_range();
                        std::tuple<int, int, int, int> ranges_tuple = live_clusters.at(i)->get_uvwt_range();
                        std::vector<int> ranges = {std::get<0>(ranges_tuple), std::get<1>(ranges_tuple), std::get<2>(ranges_tuple), std::get<3>(ranges_tuple)};
                        double length_1 = sqrt(2. / 3. *
                                                   (pow(mp.pitch_u * ranges.at(0), 2) + pow(mp.pitch_v * ranges.at(1), 2) +
                                                    pow(mp.pitch_w * ranges.at(2), 2)) +
                                               pow(time_slice_width * ranges.at(3), 2));
                        for (int j = 0; j != 4; j++) {
                            if (ranges.at(j) > max) {
                                max = ranges.at(j);
                                max_length = length_1;
                            }
                        }
                        if (max >= range_cut || max_length >= length_cut) break;
                    }

                    for (size_t j = 0; j != sep_clusters.size(); j++) {
                        delete sep_clusters.at(j);
                    }

                    if (max < range_cut && max_length < length_cut) {
                        small_clusters.push_back(live_clusters.at(i));
                    }
                    else {
                        big_clusters.push_back(live_clusters.at(i));
                    }
                }
                else {
                    big_clusters.push_back(live_clusters.at(i));
                }
            }
            else {
                big_clusters.push_back(live_clusters.at(i));
            }
        }
    }

    std::set<std::pair<Cluster *, Cluster *>> to_be_merged_pairs;

    // clustering small with big ones ...
    double small_big_dis_cut = 80 * units::cm;
    std::set<Cluster *> used_small_clusters;
    for (auto it = small_clusters.begin(); it != small_clusters.end(); it++) {
        Cluster *curr_cluster = (*it);
        // curr_cluster->Create_point_cloud();
        // ToyPointCloud *cloud1 = curr_cluster->get_point_cloud();
        double min_dis = 1e9;
        Cluster *min_dis_cluster = 0;

        for (auto it1 = big_clusters.begin(); it1 != big_clusters.end(); it1++) {
            Cluster *big_cluster = (*it1);
            // big_cluster->Create_point_cloud();
            // ToyPointCloud *cloud2 = big_cluster->get_point_cloud();
            // std::tuple<int, int, double> results = cloud2->get_closest_points(cloud1);
            std::tuple<int, int, double> results = big_cluster->get_closest_points(*curr_cluster);
            double dis = std::get<2>(results);
            if (dis < min_dis) {
                min_dis = dis;
                min_dis_cluster = big_cluster;
            }
        }
        if (min_dis < small_big_dis_cut) {
            to_be_merged_pairs.insert(std::make_pair(min_dis_cluster, curr_cluster));
            used_small_clusters.insert(curr_cluster);
        }
    }

    // small distance ...
    double small_small_dis_cut = 5 * units::cm;
    for (size_t i = 0; i != small_clusters.size(); i++) {
        Cluster *cluster1 = small_clusters.at(i);
        // ToyPointCloud *cloud1 = cluster1->get_point_cloud();
        for (size_t j = i + 1; j != small_clusters.size(); j++) {
            Cluster *cluster2 = small_clusters.at(j);
            // ToyPointCloud *cloud2 = cluster2->get_point_cloud();
            // std::tuple<int, int, double> results = cloud2->get_closest_points(cloud1);
            std::tuple<int, int, double> results = cluster2->get_closest_points(*cluster1);
            double dis = std::get<2>(results);
            if (dis < small_small_dis_cut) {
                if (used_small_clusters.find(cluster1) != used_small_clusters.end() &&
                        used_small_clusters.find(cluster2) == used_small_clusters.end() ||
                    used_small_clusters.find(cluster2) != used_small_clusters.end() &&
                        used_small_clusters.find(cluster1) == used_small_clusters.end()) {
                    to_be_merged_pairs.insert(std::make_pair(cluster1, cluster2));
                    used_small_clusters.insert(cluster1);
                    used_small_clusters.insert(cluster2);
                }
            }
        }
    }

    std::vector<Cluster *> remaining_small_clusters;
    for (auto it = small_clusters.begin(); it != small_clusters.end(); it++) {
        Cluster *curr_cluster = (*it);
        if (used_small_clusters.find(curr_cluster) == used_small_clusters.end())
            remaining_small_clusters.push_back(curr_cluster);
    }

    // clustering small with small ones ...
    small_small_dis_cut = 50 * units::cm;
    for (size_t i = 0; i != remaining_small_clusters.size(); i++) {
        Cluster *cluster1 = remaining_small_clusters.at(i);
        // ToyPointCloud *cloud1 = cluster1->get_point_cloud();
        for (size_t j = i + 1; j != remaining_small_clusters.size(); j++) {
            Cluster *cluster2 = remaining_small_clusters.at(j);
            // ToyPointCloud *cloud2 = cluster2->get_point_cloud();
            // std::tuple<int, int, double> results = cloud2->get_closest_points(cloud1);
            std::tuple<int, int, double> results = cluster2->get_closest_points(*cluster1);
            double dis = std::get<2>(results);
            if (dis < small_small_dis_cut) {
                to_be_merged_pairs.insert(std::make_pair(cluster1, cluster2));
            }
        }
    }

    // clustering big ones ...
    // cloud1 is the longer one
    // used_big_clusters holds the shorter one
    std::set<Cluster *> used_big_clusters;
    double big_dis_cut = 3 * units::cm;
    double big_dis_range_cut = 16 * units::cm;
    for (size_t i = 0; i != big_clusters.size(); i++) {
        Cluster *cluster1 = big_clusters.at(i);
        // cluster1->Create_point_cloud();
        for (size_t j = i + 1; j != big_clusters.size(); j++) {
            Cluster *cluster2 = big_clusters.at(j);
            // cluster2->Create_point_cloud();
            // ToyPointCloud *cloud1, *cloud2;
            double pca_ratio, small_cluster_length;
            // if (cluster1->get_length() > cluster2->get_length()) {
            //     // cloud1 = cluster1->get_point_cloud();
            //     // cloud2 = cluster2->get_point_cloud();
            //     if (used_big_clusters.find(cluster2) != used_big_clusters.end()) continue;
            //     // cluster2->Calc_pca();
            //     pca_ratio = cluster2->get_pca_value(1) / cluster2->get_pca_value(0);
            //     small_cluster_length = cluster2->get_length();
            // }
            // else {
            //     // cloud1 = cluster2->get_point_cloud();
            //     // cloud2 = cluster1->get_point_cloud();
            //     if (used_big_clusters.find(cluster1) != used_big_clusters.end()) continue;
            //     // cluster1->Calc_pca();
            //     pca_ratio = cluster1->get_pca_value(1) / cluster1->get_pca_value(0);
            //     small_cluster_length = cluster1->get_length();
            // }
            // make sure cluster1 is the longer one
            if (!(cluster1->get_length() > cluster2->get_length())) std::swap(cluster1, cluster2);
            if (used_big_clusters.find(cluster2) != used_big_clusters.end()) continue;
            pca_ratio = cluster2->get_pca_value(1) / cluster2->get_pca_value(0);
            small_cluster_length = cluster2->get_length();

            // std::tuple<int, int, double> results = cloud2->get_closest_points(cloud1);
            std::tuple<int, int, double> results = cluster2->get_closest_points(*cluster1);
            double min_dis = std::get<2>(results);
            if (min_dis < big_dis_cut) {
                bool flag_merge = true;

                int num_outside_points = 0;
                /* int num_close_points = 0; */
                // const int N = cloud2->get_num_points();
                const int N = cluster2->npoints();
                // WCP::WCPointCloud<double> &cloud = cloud2->get_cloud();
                for (int k = 0; k != N; k++) {
                    // Point test_p1(cloud.pts[k].x, cloud.pts[k].y, cloud.pts[k].z);
                    geo_point_t test_p1 = cluster2->point(k);
                    // double close_dis = cloud1->get_closest_dis(test_p1);
                    double close_dis = cluster1->get_closest_dis(test_p1);
                    if (close_dis > big_dis_range_cut) {
                        num_outside_points++;
                        if (num_outside_points > 400) break;
                    }
                    /* if (close_dis <= big_dis_cut) */
                    /*   num_close_points ++;  */
                }

                // std::cout << num_outside_points << " " << N << " " << num_close_points << std::endl;

                if (num_outside_points > 0.125 * N || num_outside_points > 400) flag_merge = false;

                //	std::cout << pca_ratio << " " << small_cluster_length/units::cm << std::endl;
                if (flag_merge && small_cluster_length > 60 * units::cm && pca_ratio < 0.0015) flag_merge = false;

                if (flag_merge) {
                    to_be_merged_pairs.insert(std::make_pair(cluster1, cluster2));
                    if (cluster1->get_length() < cluster2->get_length()) {
                        used_big_clusters.insert(cluster1);
                    }
                    else {
                        used_big_clusters.insert(cluster2);
                    }
                }
            }
        }
    }

    // merge clusters
    std::vector<std::set<Cluster *>> merge_clusters;
    std::set<Cluster *> used_clusters;

    for (auto it = to_be_merged_pairs.begin(); it != to_be_merged_pairs.end(); it++) {
        Cluster *cluster1 = (*it).first;
        Cluster *cluster2 = (*it).second;
        /* std::cout << cluster1 << " " << cluster2 << " " << cluster1->get_cluster_id() << " " <<
         * cluster2->get_cluster_id() << " " << cluster1->get_length()/units::cm << " " <<
         * cluster2->get_length()/units::cm << std::endl; */

        used_clusters.insert(cluster1);
        used_clusters.insert(cluster2);

        bool flag_new = true;
        std::vector<std::set<Cluster *>> temp_set;
        for (auto it1 = merge_clusters.begin(); it1 != merge_clusters.end(); it1++) {
            std::set<Cluster *> &clusters = (*it1);
            if (clusters.find(cluster1) != clusters.end() || clusters.find(cluster2) != clusters.end()) {
                clusters.insert(cluster1);
                clusters.insert(cluster2);
                flag_new = false;
                temp_set.push_back(clusters);
                // break;
            }
        }
        if (flag_new) {
            std::set<Cluster *> clusters;
            clusters.insert(cluster1);
            clusters.insert(cluster2);
            merge_clusters.push_back(clusters);
        }
        if (temp_set.size() > 1) {
            // merge them further ...
            std::set<Cluster *> clusters;
            for (size_t i = 0; i != temp_set.size(); i++) {
                for (auto it1 = temp_set.at(i).begin(); it1 != temp_set.at(i).end(); it1++) {
                    clusters.insert(*it1);
                }
                merge_clusters.erase(find(merge_clusters.begin(), merge_clusters.end(), temp_set.at(i)));
            }
            merge_clusters.push_back(clusters);
        }
    }

    for (auto it = live_clusters.begin(); it != live_clusters.end(); it++) {
        Cluster *cluster = *it;
        if (used_clusters.find(cluster) == used_clusters.end()) {
            std::set<Cluster *> temp_clusters;
            temp_clusters.insert(cluster);
            merge_clusters.push_back(temp_clusters);
        }
    }

    // new stuff ...
    map_cluster_cluster_vec results;
    for (auto it = merge_clusters.begin(); it != merge_clusters.end(); it++) {
        std::set<Cluster *> &cluster_set = (*it);
        double max_length = 0;
        Cluster *max_cluster;
        for (auto it1 = cluster_set.begin(); it1 != cluster_set.end(); it1++) {
            Cluster *temp_cluster = (*it1);
            if (temp_cluster->get_length() > max_length) {
                max_length = temp_cluster->get_length();
                max_cluster = temp_cluster;
            }
        }
        std::vector<std::pair<Cluster *, double>> temp;
        results[max_cluster] = temp;
        // max_cluster->Create_point_cloud();
        // ToyPointCloud *cloud1 = max_cluster->get_point_cloud();

        for (auto it1 = cluster_set.begin(); it1 != cluster_set.end(); it1++) {
            Cluster *temp_cluster = (*it1);
            if (temp_cluster == max_cluster) continue;

            // temp_cluster->Create_point_cloud();
            // ToyPointCloud *cloud2 = temp_cluster->get_point_cloud();
            // std::tuple<int, int, double> Tresults = cloud2->get_closest_points(cloud1);
            std::tuple<int, int, double> Tresults = temp_cluster->get_closest_points(*max_cluster);
            double dis = std::get<2>(Tresults);
            results[max_cluster].push_back(std::make_pair(temp_cluster, dis));
        }
    }

    return results;
}