#include <WireCellClus/ClusteringFuncs.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;

static bool JudgeSeparateDec_2(const Cluster* cluster, const geo_point_t& drift_dir,
                               std::vector<geo_point_t>& boundary_points, std::vector<geo_point_t>& independent_points,
                               const double cluster_length);

void WireCell::PointCloud::Facade::clustering_separate(Grouping& live_grouping,
                                                       std::map<int, std::pair<double, double>>& dead_u_index,
                                                       std::map<int, std::pair<double, double>>& dead_v_index,
                                                       std::map<int, std::pair<double, double>>& dead_w_index)
{
    std::vector<Cluster*> live_clusters = live_grouping.children();  // copy
    // for (size_t ilive = 0; ilive < 3; ++ilive) {
    for (size_t ilive = 0; ilive < live_clusters.size(); ++ilive) {
        const auto& live = live_clusters.at(ilive);
        const size_t nblobs = live->children().size();
        // std::cout << "ilive: " << ilive
        // << " nblobs " << nblobs
        // << " npoints " << live->npoints()
        // << " counter: " << global_counter_get_closest_wcpoint << std::endl;
        /// TODO: below is just for debugging, need to wait for real impl.
        {
            std::unordered_set<size_t> nblobs_to_check = {612, 725, 2414};
            if (nblobs_to_check.find(nblobs) == nblobs_to_check.end()) continue;
            live->Create_graph();
            // std::vector<geo_point_t> boundary_points, independent_points;
            // bool sep_Dec2 = JudgeSeparateDec_2(live, {1, 0, 0}, boundary_points, independent_points, live->get_length());
            // std::cout << "sep_Dec2: " << sep_Dec2 << std::endl;
        }
    }
}

static bool JudgeSeparateDec_2(const Cluster* cluster, const geo_point_t& drift_dir,
                               std::vector<geo_point_t>& boundary_points, std::vector<geo_point_t>& independent_points,
                               const double cluster_length)
{
    boundary_points = cluster->get_hull();
    std::vector<geo_point_t> hy_points;
    std::vector<geo_point_t> ly_points;
    std::vector<geo_point_t> hz_points;
    std::vector<geo_point_t> lz_points;
    std::vector<geo_point_t> hx_points;
    std::vector<geo_point_t> lx_points;

    std::set<int> independent_surfaces;

    for (size_t j = 0; j != boundary_points.size(); j++) {
        if (j == 0) {
            hy_points.push_back(boundary_points.at(j));
            ly_points.push_back(boundary_points.at(j));
            hz_points.push_back(boundary_points.at(j));
            lz_points.push_back(boundary_points.at(j));
            hx_points.push_back(boundary_points.at(j));
            lx_points.push_back(boundary_points.at(j));
        }
        else {
            geo_point_t test_p(boundary_points.at(j).x(), boundary_points.at(j).y(), boundary_points.at(j).z());
            if (cluster->nnearby(test_p, 15 * units::cm) > 75) {
                if (boundary_points.at(j).y() > hy_points.at(0).y()) hy_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).y() < ly_points.at(0).y()) ly_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).x() > hx_points.at(0).x()) hx_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).x() < lx_points.at(0).x()) lx_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).z() > hz_points.at(0).z()) hz_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).z() < lz_points.at(0).z()) lz_points.at(0) = boundary_points.at(j);
            }
        }
    }

    bool flag_outx = false;
    if (hx_points.at(0).x() > 257 * units::cm || lx_points.at(0).x() < -1 * units::cm) flag_outx = true;

    if (hy_points.at(0).y() > 101.5 * units::cm) {
        for (size_t j = 0; j != boundary_points.size(); j++) {
            if (boundary_points.at(j).y() > 101.5 * units::cm) {
                bool flag_save = true;
                for (size_t k = 0; k != hy_points.size(); k++) {
                    double dis = sqrt(pow(hy_points.at(k).x() - boundary_points.at(j).x(), 2) +
                                      pow(hy_points.at(k).y() - boundary_points.at(j).y(), 2) +
                                      pow(hy_points.at(k).z() - boundary_points.at(j).z(), 2));
                    if (dis < 25 * units::cm) {
                        if (boundary_points.at(j).y() > hy_points.at(k).y()) hy_points.at(k) = boundary_points.at(j);
                        flag_save = false;
                    }
                }
                if (flag_save) hy_points.push_back(boundary_points.at(j));
            }
        }
    }

    if (ly_points.at(0).y() < -99.5 * units::cm) {
        for (size_t j = 0; j != boundary_points.size(); j++) {
            if (boundary_points.at(j).y() < -99.5 * units::cm) {
                bool flag_save = true;
                for (size_t k = 0; k != ly_points.size(); k++) {
                    double dis = sqrt(pow(ly_points.at(k).x() - boundary_points.at(j).x(), 2) +
                                      pow(ly_points.at(k).y() - boundary_points.at(j).y(), 2) +
                                      pow(ly_points.at(k).z() - boundary_points.at(j).z(), 2));
                    if (dis < 25 * units::cm) {
                        if (boundary_points.at(j).y() < ly_points.at(k).y()) ly_points.at(k) = boundary_points.at(j);
                        flag_save = false;
                    }
                }
                if (flag_save) ly_points.push_back(boundary_points.at(j));
            }
        }
    }
    if (hz_points.at(0).z() > 1022 * units::cm) {
        for (size_t j = 0; j != boundary_points.size(); j++) {
            if (boundary_points.at(j).z() > 1022 * units::cm) {
                bool flag_save = true;
                for (size_t k = 0; k != hz_points.size(); k++) {
                    double dis = sqrt(pow(hz_points.at(k).x() - boundary_points.at(j).x(), 2) +
                                      pow(hz_points.at(k).y() - boundary_points.at(j).y(), 2) +
                                      pow(hz_points.at(k).z() - boundary_points.at(j).z(), 2));
                    if (dis < 25 * units::cm) {
                        if (boundary_points.at(j).z() > hz_points.at(k).z()) hz_points.at(k) = boundary_points.at(j);
                        flag_save = false;
                    }
                }
                if (flag_save) hz_points.push_back(boundary_points.at(j));
            }
        }
    }
    if (lz_points.at(0).z() < 15 * units::cm) {
        for (size_t j = 0; j != boundary_points.size(); j++) {
            if (boundary_points.at(j).z() < 15 * units::cm) {
                bool flag_save = true;
                for (size_t k = 0; k != lz_points.size(); k++) {
                    double dis = sqrt(pow(lz_points.at(k).x() - boundary_points.at(j).x(), 2) +
                                      pow(lz_points.at(k).y() - boundary_points.at(j).y(), 2) +
                                      pow(lz_points.at(k).z() - boundary_points.at(j).z(), 2));
                    if (dis < 25 * units::cm) {
                        if (boundary_points.at(j).z() < lz_points.at(k).z()) lz_points.at(k) = boundary_points.at(j);
                        flag_save = false;
                    }
                }
                if (flag_save) lz_points.push_back(boundary_points.at(j));
            }
        }
    }

    int num_outside_points = 0;
    int num_outx_points = 0;

    for (size_t j = 0; j != hy_points.size(); j++) {
        if (hy_points.at(j).x() >= 1 * units::cm && hy_points.at(j).x() <= 255 * units::cm &&
            hy_points.at(j).y() >= -99.5 * units::cm && hy_points.at(j).y() <= 101.5 * units::cm &&
            hy_points.at(j).z() >= 15 * units::cm && hy_points.at(j).z() <= 1022 * units::cm && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(hy_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(hy_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(hy_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(hy_points.at(j));
            if (hy_points.at(j).y() > 104 * units::cm) {
                independent_surfaces.insert(0);
            }
            else if (hy_points.at(j).y() < -99.5 * units::cm) {
                independent_surfaces.insert(1);
            }
            else if (hy_points.at(j).z() > 1025 * units::cm) {
                independent_surfaces.insert(2);
            }
            else if (hy_points.at(j).z() < 12 * units::cm) {
                independent_surfaces.insert(3);
            }
            else if (hy_points.at(j).x() > 255 * units::cm) {
                independent_surfaces.insert(4);
            }
            else if (hy_points.at(j).x() < 1 * units::cm) {
                independent_surfaces.insert(5);
            }

            if (hy_points.at(j).y() > 104 * units::cm || hy_points.at(j).y() < -99.5 * units::cm ||
                hy_points.at(j).z() < 12 * units::cm || hy_points.at(j).z() > 1025 * units::cm ||
                hy_points.at(j).x() < 1 * units::cm || hy_points.at(j).x() > 255 * units::cm)
                num_outside_points++;
            if (hy_points.at(j).x() < -1 * units::cm || hy_points.at(j).x() > 257 * units::cm) num_outx_points++;
        }
    }
    for (size_t j = 0; j != ly_points.size(); j++) {
        if (ly_points.at(j).x() >= 1 * units::cm && ly_points.at(j).x() <= 255 * units::cm &&
            ly_points.at(j).y() >= -99.5 * units::cm && ly_points.at(j).y() <= 101.5 * units::cm &&
            ly_points.at(j).z() >= 15 * units::cm && ly_points.at(j).z() <= 1022 * units::cm && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(ly_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(ly_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(ly_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(ly_points.at(j));

            if (ly_points.at(j).y() < -99.5 * units::cm) {
                independent_surfaces.insert(1);
            }
            else if (ly_points.at(j).y() > 104 * units::cm) {
                independent_surfaces.insert(0);
            }
            else if (ly_points.at(j).z() > 1025 * units::cm) {
                independent_surfaces.insert(2);
            }
            else if (ly_points.at(j).z() < 12 * units::cm) {
                independent_surfaces.insert(3);
            }
            else if (ly_points.at(j).x() > 255 * units::cm) {
                independent_surfaces.insert(4);
            }
            else if (ly_points.at(j).x() < 1 * units::cm) {
                independent_surfaces.insert(5);
            }

            if (ly_points.at(j).y() > 104 * units::cm || ly_points.at(j).y() < -99.5 * units::cm ||
                ly_points.at(j).z() < 12 * units::cm || ly_points.at(j).z() > 1025 * units::cm ||
                ly_points.at(j).x() < 1 * units::cm || ly_points.at(j).x() > 255 * units::cm)
                num_outside_points++;
            if (ly_points.at(j).x() < -1 * units::cm || ly_points.at(j).x() > 257 * units::cm) num_outx_points++;
        }
    }
    for (size_t j = 0; j != hz_points.size(); j++) {
        if (hz_points.at(j).x() >= 1 * units::cm && hz_points.at(j).x() <= 255 * units::cm &&
            hz_points.at(j).y() >= -99.5 * units::cm && hz_points.at(j).y() <= 101.5 * units::cm &&
            hz_points.at(j).z() >= 15 * units::cm && hz_points.at(j).z() <= 1022 * units::cm && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(hz_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(hz_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(hz_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(hz_points.at(j));

            if (hz_points.at(j).z() > 1025 * units::cm) {
                independent_surfaces.insert(2);
            }
            else if (hz_points.at(j).z() < 12 * units::cm) {
                independent_surfaces.insert(3);
            }
            else if (hz_points.at(j).y() > 104 * units::cm) {
                independent_surfaces.insert(0);
            }
            else if (hz_points.at(j).y() < -99.5 * units::cm) {
                independent_surfaces.insert(1);
            }
            else if (hz_points.at(j).x() > 255 * units::cm) {
                independent_surfaces.insert(4);
            }
            else if (hz_points.at(j).x() < 1 * units::cm) {
                independent_surfaces.insert(5);
            }

            if (hz_points.at(j).y() > 104 * units::cm || hz_points.at(j).y() < -99.5 * units::cm ||
                hz_points.at(j).z() < 12 * units::cm || hz_points.at(j).z() > 1025 * units::cm ||
                hz_points.at(j).x() < 1 * units::cm || hz_points.at(j).x() > 255 * units::cm)
                num_outside_points++;
            if (hz_points.at(j).x() < -1 * units::cm || hz_points.at(j).x() > 257 * units::cm) num_outx_points++;
        }
    }
    for (size_t j = 0; j != lz_points.size(); j++) {
        if (lz_points.at(j).x() >= 1 * units::cm && lz_points.at(j).x() <= 255 * units::cm &&
            lz_points.at(j).y() >= -99.5 * units::cm && lz_points.at(j).y() <= 101.5 * units::cm &&
            lz_points.at(j).z() >= 15 * units::cm && lz_points.at(j).z() <= 1022 * units::cm && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(lz_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(lz_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(lz_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(lz_points.at(j));

            if (lz_points.at(j).z() < 12 * units::cm) {
                independent_surfaces.insert(3);
            }
            else if (lz_points.at(j).z() > 1025 * units::cm) {
                independent_surfaces.insert(2);
            }
            else if (lz_points.at(j).y() > 104 * units::cm) {
                independent_surfaces.insert(0);
            }
            else if (lz_points.at(j).y() < -99.5 * units::cm) {
                independent_surfaces.insert(1);
            }
            else if (lz_points.at(j).x() > 255 * units::cm) {
                independent_surfaces.insert(4);
            }
            else if (lz_points.at(j).x() < 1 * units::cm) {
                independent_surfaces.insert(5);
            }

            if (lz_points.at(j).y() > 104 * units::cm || lz_points.at(j).y() < -99.5 * units::cm ||
                lz_points.at(j).z() < 12 * units::cm || lz_points.at(j).z() > 1025 * units::cm ||
                lz_points.at(j).x() < 1 * units::cm || lz_points.at(j).x() > 255 * units::cm)
                num_outside_points++;
            if (lz_points.at(j).x() < -1 * units::cm || lz_points.at(j).x() > 257 * units::cm) num_outx_points++;
        }
    }
    for (size_t j = 0; j != hx_points.size(); j++) {
        if (hx_points.at(j).x() >= 1 * units::cm && hx_points.at(j).x() <= 255 * units::cm &&
            hx_points.at(j).y() >= -99.5 * units::cm && hx_points.at(j).y() <= 101.5 * units::cm &&
            hx_points.at(j).z() >= 15 * units::cm && hx_points.at(j).z() <= 1022 * units::cm && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(hx_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(hx_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(hx_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(hx_points.at(j));

            if (hx_points.at(j).y() > 104 * units::cm || hx_points.at(j).y() < -99.5 * units::cm ||
                hx_points.at(j).z() < 12 * units::cm || hx_points.at(j).z() > 1025 * units::cm ||
                hx_points.at(j).x() < 1 * units::cm || hx_points.at(j).x() > 255 * units::cm)
                num_outside_points++;
            if (hx_points.at(j).x() < -1 * units::cm || hx_points.at(j).x() > 257 * units::cm) {
                num_outx_points++;
            }

            if (lx_points.at(j).x() > 255 * units::cm) {
                independent_surfaces.insert(4);
            }
            else if (lx_points.at(j).x() < 1 * units::cm) {
                independent_surfaces.insert(5);
            }
            else if (lx_points.at(j).y() > 104 * units::cm) {
                independent_surfaces.insert(0);
            }
            else if (lx_points.at(j).y() < -99.5 * units::cm) {
                independent_surfaces.insert(1);
            }
            else if (lx_points.at(j).z() > 1025 * units::cm) {
                independent_surfaces.insert(2);
            }
            else if (lx_points.at(j).z() < 12 * units::cm) {
                independent_surfaces.insert(3);
            }
        }
    }
    for (size_t j = 0; j != lx_points.size(); j++) {
        if (lx_points.at(j).x() >= 1 * units::cm && lx_points.at(j).x() <= 255 * units::cm &&
            lx_points.at(j).y() >= -99.5 * units::cm && lx_points.at(j).y() <= 101.5 * units::cm &&
            lx_points.at(j).z() >= 15 * units::cm && lx_points.at(j).z() <= 1022 * units::cm && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(lx_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(lx_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(lx_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(lx_points.at(j));

            if (lx_points.at(j).y() > 104 * units::cm || lx_points.at(j).y() < -99.5 * units::cm ||
                lx_points.at(j).z() < 12 * units::cm || lx_points.at(j).z() > 1025 * units::cm ||
                lx_points.at(j).x() < 1 * units::cm || lx_points.at(j).x() > 255 * units::cm)
                num_outside_points++;
            if (lx_points.at(j).x() < -1 * units::cm || lx_points.at(j).x() > 257 * units::cm) {
                num_outx_points++;
            }

            if (lx_points.at(j).x() < 1 * units::cm) {
                independent_surfaces.insert(5);
            }
            else if (lx_points.at(j).x() > 255 * units::cm) {
                independent_surfaces.insert(4);
            }
            else if (lx_points.at(j).y() > 104 * units::cm) {
                independent_surfaces.insert(0);
            }
            else if (lx_points.at(j).y() < -99.5 * units::cm) {
                independent_surfaces.insert(1);
            }
            else if (lx_points.at(j).z() > 1025 * units::cm) {
                independent_surfaces.insert(2);
            }
            else if (lx_points.at(j).z() < 12 * units::cm) {
                independent_surfaces.insert(3);
            }
        }
    }

    int num_far_points = 0;

    if (independent_points.size() == 2 && (independent_surfaces.size() > 1 || flag_outx)) {
        geo_vector_t dir_1(independent_points.at(1).x() - independent_points.at(0).x(),
                           independent_points.at(1).y() - independent_points.at(0).y(),
                           independent_points.at(1).z() - independent_points.at(0).z());
        dir_1 = dir_1.norm();
        for (size_t j = 0; j != boundary_points.size(); j++) {
            geo_vector_t dir_2(boundary_points.at(j).x() - independent_points.at(0).x(),
                               boundary_points.at(j).y() - independent_points.at(0).y(),
                               boundary_points.at(j).z() - independent_points.at(0).z());
            double angle_12 = dir_1.angle(dir_2);
            geo_vector_t dir_3 = dir_2 - dir_1 * dir_2.magnitude() * cos(angle_12);
            double angle_3 = dir_3.angle(drift_dir);
            // std::cout << dir_3.Mag()/units::cm << " " << fabs(angle_3-3.1415926/2.)/3.1415926*180. << " " <<
            // fabs(dir_3.X()/units::cm) << std::endl;
            if (fabs(angle_3 - 3.1415926 / 2.) / 3.1415926 * 180. < 7.5) {
                if (fabs(dir_3.x() / units::cm) > 14 * units::cm) num_far_points++;
                if (fabs(dir_1.angle(drift_dir) - 3.1415926 / 2.) / 3.1415926 * 180. > 15) {
                    if (dir_3.magnitude() > 20 * units::cm) num_far_points++;
                }
            }
            else {
                if (dir_3.magnitude() > 20 * units::cm) num_far_points++;
            }
        }

        // find the middle points and close distance ...
        geo_point_t middle_point((independent_points.at(1).x() + independent_points.at(0).x()) / 2.,
                                 (independent_points.at(1).y() + independent_points.at(0).y()) / 2.,
                                 (independent_points.at(1).z() + independent_points.at(0).z()) / 2.);
        const auto knn_res = cluster->kd_knn(1, middle_point);
        if (knn_res.size() != 1) {
            raise<ValueError>("knn_res.size() %d != 1", knn_res.size());
        }
        double middle_dis = knn_res[0].second;
        // std::cout << middle_dis/units::cm << " " << num_far_points << std::endl;
        if (middle_dis > 25 * units::cm) {
            num_far_points = 0;
        }
    }

    double max_x = -1e9, min_x = 1e9;
    double max_y = -1e9, min_y = 1e9;
    double max_z = -1e9, min_z = 1e9;
    for (auto it = independent_points.begin(); it != independent_points.end(); it++) {
        if ((*it).x() > max_x) max_x = (*it).x();
        if ((*it).x() < min_x) min_x = (*it).x();
        if ((*it).y() > max_y) max_y = (*it).y();
        if ((*it).y() < min_y) min_y = (*it).y();
        if ((*it).z() > max_z) max_z = (*it).z();
        if ((*it).z() < min_z) min_z = (*it).z();
        // std::cout << (*it).x()/units::cm << " " << (*it).y()/units::cm << " " << (*it).z()/units::cm << std::endl;
    }
    if (hx_points.size() > 0) {
        if (hx_points.at(0).x() > max_x + 10 * units::cm) max_x = hx_points.at(0).x();
    }
    if (lx_points.size() > 0) {
        if (lx_points.at(0).x() < min_x - 10 * units::cm) min_x = lx_points.at(0).x();
    }

    if (max_x - min_x < 2.5 * units::cm &&
        sqrt(pow(max_y - min_y, 2) + pow(max_z - min_z, 2) + pow(max_x - min_x, 2)) > 150 * units::cm) {
        independent_points.clear();
        return false;
    }
    if (max_x - min_x < 2.5 * units::cm && independent_points.size() == 2 && num_outx_points == 0) {
        independent_points.clear();
        return false;
    }

    if ((num_outside_points > 1 && independent_surfaces.size() > 1 ||
         num_outside_points > 2 && cluster_length > 250 * units::cm || num_outx_points > 0) &&
        (independent_points.size() > 2 || independent_points.size() == 2 && num_far_points > 0))
        return true;

    // about to return false ...
    independent_points.clear();

    for (size_t j = 0; j != hy_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(hy_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(hy_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(hy_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(hy_points.at(j));
    }

    for (size_t j = 0; j != ly_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(ly_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(ly_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(ly_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(ly_points.at(j));
    }

    for (size_t j = 0; j != hx_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(hx_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(hx_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(hx_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(hx_points.at(j));
    }

    for (size_t j = 0; j != lx_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(lx_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(lx_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(lx_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(lx_points.at(j));
    }

    for (size_t j = 0; j != hz_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(hz_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(hz_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(hz_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(hz_points.at(j));
    }

    for (size_t j = 0; j != lz_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(lz_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(lz_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(lz_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(lz_points.at(j));
    }

    return false;
}