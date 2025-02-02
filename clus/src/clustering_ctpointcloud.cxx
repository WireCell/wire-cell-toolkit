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

    std::cout << "Test Ave Charge: " << live_grouping.get_ave_3d_charge(p,1*units::cm,0) << " " << live_grouping.get_ave_charge(p,1*units::cm,0,0) << " " << live_grouping.get_ave_charge(p,1*units::cm,0,1) << " " << live_grouping.get_ave_charge(p,1*units::cm,0,2) << " " << std::endl; 

    auto point_u = live_grouping.convert_time_ch_2Dpoint(10*4, 10, 0, 0);
    auto point_v = live_grouping.convert_time_ch_2Dpoint(10*4, 10, 0, 1);
    auto point_w = live_grouping.convert_time_ch_2Dpoint(10*4, 10, 0, 2);

    std::cout << "Test 2D Conversion " << point_u.first << " " << point_u.second << " " << point_v.first << " " << point_v.second << " "  << point_w.first << " " << point_w.second << std::endl;

    auto dead_chs_u = live_grouping.get_overlap_dead_chs(10*4,1000*4,0,2400,0,0);
    auto dead_chs_v = live_grouping.get_overlap_dead_chs(10*4,1000*4,0,2400,0,1);
    auto dead_chs_w = live_grouping.get_overlap_dead_chs(10*4,1000*4,0,4800,0,2);
    std::cout << "Test Overlap dead chs: " << dead_chs_u.size() << " " << dead_chs_v.size() << " " << dead_chs_w.size() << std::endl;

    std::cout << "Test all dead chs: " << live_grouping.get_all_dead_chs(0,0).size() + live_grouping.get_all_dead_chs(0,1).size() + live_grouping.get_all_dead_chs(0,2).size() << std::endl;

    auto good_chs_u = live_grouping.get_overlap_good_ch_charge(10*4,1000*4,0,2400,0,0);
    auto good_chs_v = live_grouping.get_overlap_good_ch_charge(10*4,1000*4,0,2400,0,1);
    auto good_chs_w = live_grouping.get_overlap_good_ch_charge(10*4,1000*4,0,4800,0,2);

    std::cout << "Test all good chs: " << good_chs_u.size() << " " << good_chs_v.size() << " " << good_chs_w.size() << std::endl;

// run 5384 130 6501 results
//     Test CTPointCloud
// Test: (-1204.49 -57.85 5635)
// Test is_good_point: 1 1
// Test get_closest_points: 5 5 6
// Test get_closest_dead_chs: 0 0 0
// Test convert_3Dpoint_time_ch: 253 1294 3655 6678
// Test Number of Points: 42097 42335 42029
// Test test_good_point: 1 1 1 0 0 0
// Test Ave Charge: 3549.97 3189.11 4524.67 2936.12 
// Test 2D Conversion -1739.58 -984.392 -1739.58 -967.591 -1739.58 32.5
// Test Overlap dead chs: 33 40 22
// Test all dead chs: 1027
// Test all good chs: 11233 11063 12483

    std::cout << "Test new functions in Facade::Cluster " << std::endl;
    std::vector<Cluster *> live_clusters = live_grouping.children();
    for (size_t i = 0; i != live_clusters.size(); i++) {
        auto results = live_clusters.at(i)->get_closest_wcpoint(p);
        if (live_clusters.at(i)->get_length()/units::cm >239){
            std::cout << "Test: " << live_clusters.at(i)->get_length()/units::cm << " " << results.second << std::endl;
            auto points = live_clusters.at(i)->get_main_axis_points();
            std::cout << "Test: " << points.first << " " << points.second << std::endl;

            auto dir1 = live_clusters.at(i)->calc_dir(points.first, points.second, 10*units::cm);
            std::cout << "Test: " << dir1 << std::endl;

            std::vector<geo_point_t> points1;
            points1.push_back(points.first);
            points1.push_back(points.second);
            points1.push_back(p);
            live_clusters.at(i)->Calc_PCA(points1);
            std::cout << "Test: " << live_clusters.at(i)->get_center() << " " << live_clusters.at(i)->get_pca_axis(0) << " " << live_clusters.at(i)->get_pca_axis(1) << " " << live_clusters.at(i)->get_pca_axis(2) << std::endl;

            geo_point_t p2(0,0,0);
            auto dir2 = live_clusters.at(i)->calc_pca_dir(p2, points1);
            std::cout << "Test: " << dir2 << std::endl;

            auto p5 = live_clusters.at(i)->calc_ave_pos(points.first, 10);
            std::cout << "Test: " << p5 << std::endl;
        }
    }
}