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

void WireCell::PointCloud::Facade::clustering_ctpointcloud(Grouping& live_grouping){

    // test a few different functions and then print out ...
    std::cout << "Test CTPointCloud" << std::endl;

    // geo_point_t p(0*units::cm, 0*units::cm, 0*units::cm);
    // std::vector<Cluster *> live_clusters = live_grouping.children();
    // for (size_t i = 0; i != live_clusters.size(); i++) {
    //     auto results = live_clusters.at(i)->get_closest_wcpoint(p);
    //     if (live_clusters.at(i)->get_length()/units::cm >200)
    //     std::cout << "Test: " << live_clusters.at(i)->get_length()/units::cm << " " << results.second << std::endl;

    // }


    geo_point_t p(-1204.49*units::mm, -57.85*units::mm, 5635*units::mm);
    std::cout << "Test: " << p << std::endl;

    bool flag = live_grouping.is_good_point(p,0,0.6*units::cm, 1, 1);
    bool flag_wc = live_grouping.is_good_point_wc(p,0,0.6*units::cm, 1, 1);
    std::cout << "Test is_good_point: " << flag << " " << flag_wc << std::endl;

    auto closest_points_u = live_grouping.get_closest_points(p,0.6*units::cm, 0, 0);
    auto closest_points_v = live_grouping.get_closest_points(p,0.6*units::cm, 0, 1);
    auto closest_points_w = live_grouping.get_closest_points(p,0.6*units::cm, 0, 2);
    std::cout << "Test get_closest_points: " << closest_points_u.size() << " " << closest_points_v.size() << " " << closest_points_w.size() << std::endl;

    bool flag_dead_chs_u = live_grouping.get_closest_dead_chs(p,1,0,0);
    bool flag_dead_chs_v = live_grouping.get_closest_dead_chs(p,1,0,1);
    bool flag_dead_chs_w = live_grouping.get_closest_dead_chs(p,1,0,2);
    std::cout << "Test get_closest_dead_chs: " << flag_dead_chs_u << " " << flag_dead_chs_v << " " << flag_dead_chs_w << std::endl;

    auto time_ch_u = live_grouping.convert_3Dpoint_time_ch(p,0,0);
    auto time_ch_v = live_grouping.convert_3Dpoint_time_ch(p,0,1);
    auto time_ch_w = live_grouping.convert_3Dpoint_time_ch(p,0,2);
    std::cout << "Test convert_3Dpoint_time_ch: " << int(std::get<0>(time_ch_u)/4) << " " << std::get<1>(time_ch_u) <<  " " << std::get<1>(time_ch_v)+2400 << " " << std::get<1>(time_ch_w)+4800 << std::endl;

    std::cout << "Test Number of Points: " << live_grouping.get_num_points(0,0) << " " << live_grouping.get_num_points(0,1) << " " << live_grouping.get_num_points(0,2) << std::endl;

    auto num_planes = live_grouping.test_good_point(p,0,0.6*units::cm, 1);
    std::cout << "Test test_good_point: " << num_planes[0] << " " << num_planes[1] << " " << num_planes[2] << " " << num_planes[3] << " " << num_planes[4] << " " << num_planes[5] << std::endl;

}