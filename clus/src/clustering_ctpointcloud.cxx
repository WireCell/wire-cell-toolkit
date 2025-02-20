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

            // test shortest path 
            auto start_wcpoint_idx = live_clusters.at(i)->get_closest_point_index(points.first);
            auto end_wcpoint_idx = live_clusters.at(i)->get_closest_point_index(points.second);

            live_clusters.at(i)->dijkstra_shortest_paths(start_wcpoint_idx, true);
            live_clusters.at(i)->cal_shortest_path(end_wcpoint_idx);

            auto indices = live_clusters.at(i)->get_path_wcps();
            auto points2 = live_clusters.at(i)->indices_to_points(indices);
            std::cout << "Test shortest path: " << points2.size() << " " << points2.at(0) << " " << points2.at(points2.size()-1) << std::endl;

            std::vector<geo_point_t> points6;
            {geo_point_t temp_p(592.338, 1144.19, 1897); points6.push_back(temp_p);}
            {geo_point_t temp_p(592.338, 1135.53, 1900); points6.push_back(temp_p);}
            {geo_point_t temp_p(587.934, 1125.14, 1900); points6.push_back(temp_p);}
            {geo_point_t temp_p(583.53, 1119.08, 1898.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(581.328, 1113.88, 1898.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(581.328, 1106.96, 1898.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(576.924, 1104.36, 1900); points6.push_back(temp_p);}
            {geo_point_t temp_p(576.924, 1102.63, 1903); points6.push_back(temp_p);}
            {geo_point_t temp_p(576.924, 1099.16, 1903); points6.push_back(temp_p);}
            {geo_point_t temp_p(576.924, 1088.77, 1903); points6.push_back(temp_p);}
            {geo_point_t temp_p(576.924, 1087.04, 1906); points6.push_back(temp_p);}
            {geo_point_t temp_p(574.722, 1084.44, 1907.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(574.722, 1074.05, 1907.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(570.318, 1070.58, 1907.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(570.318, 1060.19, 1907.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(570.318, 1056.73, 1907.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(568.116, 1054.13, 1909); points6.push_back(temp_p);}
            {geo_point_t temp_p(568.116, 1043.74, 1909); points6.push_back(temp_p);}
            {geo_point_t temp_p(568.116, 1042, 1912); points6.push_back(temp_p);}
            {geo_point_t temp_p(563.712, 1034.21, 1910.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(563.712, 1023.82, 1910.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(563.712, 1016.89, 1910.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(559.308, 1011.69, 1913.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(559.308, 999.568, 1910.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(559.308, 990.042, 1912); points6.push_back(temp_p);}
            {geo_point_t temp_p(554.904, 980.516, 1913.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(552.702, 975.319, 1913.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(550.5, 973.587, 1913.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(550.5, 963.195, 1913.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(550.5, 952.803, 1913.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(550.5, 946.74, 1915); points6.push_back(temp_p);}
            {geo_point_t temp_p(548.298, 940.678, 1916.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(543.894, 933.75, 1916.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(543.894, 923.358, 1916.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(543.894, 919.894, 1916.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(541.692, 916.43, 1916.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(537.288, 906.037, 1916.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(537.288, 895.645, 1916.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(532.884, 893.913, 1919.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(532.884, 891.315, 1921); points6.push_back(temp_p);}
            {geo_point_t temp_p(532.884, 884.387, 1921); points6.push_back(temp_p);}
            {geo_point_t temp_p(532.884, 880.922, 1921); points6.push_back(temp_p);}
            {geo_point_t temp_p(530.682, 877.458, 1921); points6.push_back(temp_p);}
            {geo_point_t temp_p(530.682, 870.53, 1921); points6.push_back(temp_p);}
            {geo_point_t temp_p(528.48, 867.932, 1919.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(524.076, 861.004, 1922.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(524.076, 850.612, 1922.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(524.076, 840.219, 1922.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(524.076, 836.755, 1922.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(521.874, 829.827, 1922.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(519.672, 822.899, 1922.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(517.47, 812.507, 1922.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(517.47, 809.043, 1922.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(515.268, 798.65, 1922.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(513.066, 795.186, 1922.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(513.066, 788.258, 1922.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(510.864, 784.794, 1922.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(510.864, 775.267, 1924); points6.push_back(temp_p);}
            {geo_point_t temp_p(508.662, 772.67, 1925.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(506.46, 762.277, 1925.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(504.258, 751.885, 1925.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(499.854, 741.492, 1925.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(499.854, 731.1, 1925.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(499.854, 725.038, 1927); points6.push_back(temp_p);}
            {geo_point_t temp_p(495.45, 715.512, 1928.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(495.45, 708.583, 1928.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(493.248, 698.191, 1928.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(491.046, 687.799, 1928.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(488.844, 677.407, 1928.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(488.844, 673.942, 1928.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(484.44, 665.282, 1931.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(480.036, 661.818, 1931.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(480.036, 651.426, 1931.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(480.036, 641.034, 1931.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(480.036, 637.57, 1931.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(477.834, 630.641, 1931.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(477.834, 623.713, 1931.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(475.632, 613.321, 1931.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(475.632, 609.857, 1931.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(473.43, 598.599, 1930); points6.push_back(temp_p);}
            {geo_point_t temp_p(471.228, 585.609, 1934.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(466.824, 576.948, 1934.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(464.622, 571.752, 1934.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(464.622, 561.36, 1934.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(460.218, 554.431, 1937.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(460.218, 544.039, 1937.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(460.218, 533.646, 1937.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(460.218, 530.182, 1937.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(458.016, 519.79, 1937.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(455.814, 509.398, 1937.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(453.612, 499.005, 1937.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(451.41, 490.345, 1940.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(451.41, 486.881, 1940.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(449.208, 476.489, 1940.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(444.804, 470.427, 1942); points6.push_back(temp_p);}
            {geo_point_t temp_p(440.4, 464.365, 1943.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(438.198, 460.9, 1943.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(438.198, 450.508, 1943.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(433.794, 444.446, 1945); points6.push_back(temp_p);}
            {geo_point_t temp_p(433.794, 440.982, 1945); points6.push_back(temp_p);}
            {geo_point_t temp_p(431.592, 437.518, 1945); points6.push_back(temp_p);}
            {geo_point_t temp_p(431.592, 426.259, 1943.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(431.592, 415.867, 1943.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(431.592, 405.475, 1943.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(431.592, 402.011, 1943.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(429.39, 395.082, 1943.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(429.39, 389.02, 1945); points6.push_back(temp_p);}
            {geo_point_t temp_p(429.39, 382.092, 1945); points6.push_back(temp_p);}
            {geo_point_t temp_p(427.188, 379.494, 1946.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(422.784, 376.03, 1946.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(422.784, 365.638, 1946.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(422.784, 359.576, 1948); points6.push_back(temp_p);}
            {geo_point_t temp_p(418.38, 348.317, 1946.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(418.38, 344.853, 1946.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(416.178, 334.461, 1946.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(413.976, 327.532, 1946.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(411.774, 317.14, 1946.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(409.572, 306.748, 1946.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(407.37, 296.355, 1946.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(405.168, 289.427, 1946.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(400.764, 280.767, 1949.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(400.764, 270.375, 1949.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(400.764, 263.447, 1949.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(398.562, 256.518, 1949.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(396.36, 249.59, 1949.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(391.956, 239.198, 1949.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(387.552, 234.868, 1948); points6.push_back(temp_p);}
            {geo_point_t temp_p(387.552, 224.476, 1948); points6.push_back(temp_p);}
            {geo_point_t temp_p(387.552, 214.083, 1948); points6.push_back(temp_p);}
            {geo_point_t temp_p(387.552, 207.155, 1948); points6.push_back(temp_p);}
            {geo_point_t temp_p(385.35, 204.557, 1949.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(385.35, 201.093, 1949.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(383.148, 198.495, 1951); points6.push_back(temp_p);}
            {geo_point_t temp_p(380.946, 192.432, 1952.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(380.946, 185.504, 1952.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(378.744, 178.576, 1952.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(374.34, 168.184, 1952.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(374.34, 157.791, 1952.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(372.138, 150.863, 1952.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(369.936, 140.471, 1952.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(367.734, 135.275, 1955.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(367.734, 128.346, 1955.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(363.33, 124.884, 1955.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(363.33, 114.488, 1955.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(363.33, 107.564, 1955.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(361.128, 104.097, 1955.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(361.128, 98.0347, 1957); points6.push_back(temp_p);}
            {geo_point_t temp_p(358.926, 90.2434, 1955.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(356.724, 83.3143, 1955.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(354.522, 76.3853, 1955.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(352.32, 67.7261, 1958.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(350.118, 58.1973, 1960); points6.push_back(temp_p);}
            {geo_point_t temp_p(347.916, 46.9389, 1958.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(343.512, 43.4768, 1958.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(343.512, 34.8176, 1961.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(343.512, 27.8886, 1961.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(341.31, 24.4215, 1961.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(341.31, 20.9595, 1961.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(336.906, 14.0304, 1961.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(336.906, 3.63921, 1961.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(336.906, 0.172097, 1961.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(334.704, -5.01991, 1964.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(330.3, -15.4161, 1964.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(328.098, -20.612, 1967.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(328.098, -31.0043, 1967.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(328.098, -37.9326, 1967.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(323.694, -41.3967, 1967.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(323.694, -44.8608, 1967.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(323.694, -47.459, 1969); points6.push_back(temp_p);}
            {geo_point_t temp_p(319.29, -48.3252, 1970.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(319.29, -58.7164, 1970.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(317.088, -67.3774, 1970.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(317.088, -75.1712, 1969); points6.push_back(temp_p);}
            {geo_point_t temp_p(314.886, -80.3688, 1972); points6.push_back(temp_p);}
            {geo_point_t temp_p(312.684, -85.564, 1975); points6.push_back(temp_p);}
            {geo_point_t temp_p(310.482, -89.8941, 1973.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(310.482, -96.8223, 1973.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(308.28, -102.884, 1975); points6.push_back(temp_p);}
            {geo_point_t temp_p(308.28, -104.616, 1978); points6.push_back(temp_p);}
            {geo_point_t temp_p(306.078, -111.545, 1978); points6.push_back(temp_p);}
            {geo_point_t temp_p(301.674, -116.74, 1981); points6.push_back(temp_p);}
            {geo_point_t temp_p(301.674, -118.472, 1984); points6.push_back(temp_p);}
            {geo_point_t temp_p(301.674, -120.204, 1987); points6.push_back(temp_p);}
            {geo_point_t temp_p(301.674, -121.938, 1984); points6.push_back(temp_p);}
            {geo_point_t temp_p(301.674, -123.67, 1987); points6.push_back(temp_p);}
            {geo_point_t temp_p(299.472, -127.133, 1987); points6.push_back(temp_p);}
            {geo_point_t temp_p(299.472, -130.596, 1987); points6.push_back(temp_p);}
            {geo_point_t temp_p(299.472, -134.062, 1987); points6.push_back(temp_p);}
            {geo_point_t temp_p(299.472, -139.258, 1984); points6.push_back(temp_p);}
            {geo_point_t temp_p(297.27, -146.186, 1984); points6.push_back(temp_p);}
            {geo_point_t temp_p(297.27, -147.918, 1987); points6.push_back(temp_p);}
            {geo_point_t temp_p(297.27, -154.846, 1987); points6.push_back(temp_p);}
            {geo_point_t temp_p(297.27, -156.578, 1984); points6.push_back(temp_p);}
            {geo_point_t temp_p(292.866, -167.836, 1985.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(290.664, -176.497, 1985.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(290.664, -182.559, 1987); points6.push_back(temp_p);}
            {geo_point_t temp_p(286.26, -192.085, 1988.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(281.856, -195.549, 1988.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(277.452, -201.611, 1990); points6.push_back(temp_p);}
            {geo_point_t temp_p(277.452, -210.272, 1993); points6.push_back(temp_p);}
            {geo_point_t temp_p(273.048, -215.468, 1990); points6.push_back(temp_p);}
            {geo_point_t temp_p(273.048, -226.726, 1988.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(273.048, -237.119, 1988.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(273.048, -244.047, 1988.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(270.846, -249.243, 1991.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(268.644, -251.841, 1993); points6.push_back(temp_p);}
            {geo_point_t temp_p(264.24, -256.171, 1991.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(262.038, -257.903, 1991.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(262.038, -264.831, 1991.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(257.634, -269.161, 1993); points6.push_back(temp_p);}
            {geo_point_t temp_p(255.432, -272.626, 1993); points6.push_back(temp_p);}
            {geo_point_t temp_p(253.23, -281.286, 1996); points6.push_back(temp_p);}
            {geo_point_t temp_p(253.23, -283.018, 1999); points6.push_back(temp_p);}
            {geo_point_t temp_p(253.23, -289.946, 1999); points6.push_back(temp_p);}
            {geo_point_t temp_p(248.826, -299.472, 1997.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(248.826, -306.4, 1997.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(246.624, -307.267, 1999); points6.push_back(temp_p);}
            {geo_point_t temp_p(246.624, -310.731, 1999); points6.push_back(temp_p);}
            {geo_point_t temp_p(246.624, -318.525, 1997.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(244.422, -321.123, 1999); points6.push_back(temp_p);}
            {geo_point_t temp_p(244.422, -325.453, 1997.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(244.422, -332.381, 1997.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(244.422, -341.908, 1999); points6.push_back(temp_p);}
            {geo_point_t temp_p(240.018, -347.97, 2000.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(237.816, -354.032, 2002); points6.push_back(temp_p);}
            {geo_point_t temp_p(235.614, -358.362, 2000.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(231.21, -369.62, 1999); points6.push_back(temp_p);}
            {geo_point_t temp_p(231.21, -374.817, 2002); points6.push_back(temp_p);}
            {geo_point_t temp_p(229.008, -380.879, 2003.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(226.806, -390.405, 2005); points6.push_back(temp_p);}
            {geo_point_t temp_p(222.402, -396.467, 2003.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(222.402, -406.859, 2003.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(217.998, -412.056, 2003.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(217.998, -418.984, 2003.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(213.594, -422.448, 2003.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(209.19, -432.84, 2003.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(209.19, -439.769, 2003.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(209.19, -442.367, 2005); points6.push_back(temp_p);}
            {geo_point_t temp_p(209.19, -445.831, 2005); points6.push_back(temp_p);}
            {geo_point_t temp_p(206.988, -452.759, 2005); points6.push_back(temp_p);}
            {geo_point_t temp_p(204.786, -459.687, 2005); points6.push_back(temp_p);}
            {geo_point_t temp_p(200.382, -463.151, 2005); points6.push_back(temp_p);}
            {geo_point_t temp_p(198.18, -464.017, 2006.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(198.18, -464.883, 2008); points6.push_back(temp_p);}
            {geo_point_t temp_p(198.18, -470.945, 2006.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(195.978, -476.142, 2006.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(191.574, -488.266, 2003.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(189.372, -494.328, 2005); points6.push_back(temp_p);}
            {geo_point_t temp_p(184.968, -500.39, 2003.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(184.968, -501.256, 2005); points6.push_back(temp_p);}
            {geo_point_t temp_p(184.968, -502.988, 2008); points6.push_back(temp_p);}
            {geo_point_t temp_p(182.766, -509.916, 2008); points6.push_back(temp_p);}
            {geo_point_t temp_p(182.766, -514.247, 2006.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(180.564, -519.443, 2006.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(180.564, -526.371, 2006.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(176.16, -534.165, 2008); points6.push_back(temp_p);}
            {geo_point_t temp_p(171.756, -542.825, 2011); points6.push_back(temp_p);}
            {geo_point_t temp_p(171.756, -553.218, 2011); points6.push_back(temp_p);}
            {geo_point_t temp_p(167.352, -556.682, 2011); points6.push_back(temp_p);}
            {geo_point_t temp_p(167.352, -561.878, 2008); points6.push_back(temp_p);}
            {geo_point_t temp_p(165.15, -562.744, 2006.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(165.15, -573.137, 2006.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(165.15, -576.6, 2006.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(160.746, -586.993, 2009.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(158.544, -590.457, 2009.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(154.14, -593.921, 2009.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(151.938, -599.116, 2009.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(149.736, -600.849, 2009.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(149.736, -611.241, 2009.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(149.736, -621.634, 2009.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(145.332, -632.025, 2012.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(140.928, -637.222, 2012.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(138.726, -643.284, 2011); points6.push_back(temp_p);}
            {geo_point_t temp_p(134.322, -648.48, 2008); points6.push_back(temp_p);}
            {geo_point_t temp_p(134.322, -650.212, 2005); points6.push_back(temp_p);}
            {geo_point_t temp_p(129.918, -658.872, 2002); points6.push_back(temp_p);}
            {geo_point_t temp_p(125.514, -664.068, 2005); points6.push_back(temp_p);}
            {geo_point_t temp_p(123.312, -670.997, 2005); points6.push_back(temp_p);}
            {geo_point_t temp_p(123.312, -677.925, 2005); points6.push_back(temp_p);}
            {geo_point_t temp_p(118.908, -684.853, 2005); points6.push_back(temp_p);}
            {geo_point_t temp_p(118.908, -690.915, 2006.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(116.706, -700.442, 2008); points6.push_back(temp_p);}
            {geo_point_t temp_p(112.302, -706.504, 2009.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(110.1, -717.762, 2008); points6.push_back(temp_p);}
            {geo_point_t temp_p(107.898, -720.36, 2009.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(103.494, -724.69, 2008); points6.push_back(temp_p);}
            {geo_point_t temp_p(99.09, -731.619, 2008); points6.push_back(temp_p);}
            {geo_point_t temp_p(96.888, -738.547, 2008); points6.push_back(temp_p);}
            {geo_point_t temp_p(96.888, -748.939, 2008); points6.push_back(temp_p);}
            {geo_point_t temp_p(92.484, -756.733, 2012.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(88.08, -767.992, 2011); points6.push_back(temp_p);}
            {geo_point_t temp_p(88.08, -774.92, 2011); points6.push_back(temp_p);}
            {geo_point_t temp_p(83.676, -783.58, 2008); points6.push_back(temp_p);}
            {geo_point_t temp_p(83.676, -792.24, 2011); points6.push_back(temp_p);}
            {geo_point_t temp_p(81.474, -797.437, 2014); points6.push_back(temp_p);}
            {geo_point_t temp_p(81.474, -807.829, 2014); points6.push_back(temp_p);}
            {geo_point_t temp_p(79.272, -813.025, 2017); points6.push_back(temp_p);}
            {geo_point_t temp_p(77.07, -820.819, 2015.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(74.868, -823.417, 2017); points6.push_back(temp_p);}
            {geo_point_t temp_p(74.868, -825.15, 2020); points6.push_back(temp_p);}
            {geo_point_t temp_p(70.464, -831.211, 2021.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(70.464, -838.14, 2021.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(68.262, -839.872, 2024.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(68.262, -843.336, 2024.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(66.06, -850.264, 2024.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(61.656, -860.656, 2024.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(57.252, -867.585, 2027.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(52.848, -872.781, 2027.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(50.646, -877.977, 2027.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(50.646, -884.905, 2027.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(48.444, -890.101, 2030.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(44.04, -896.163, 2029); points6.push_back(temp_p);}
            {geo_point_t temp_p(44.04, -903.092, 2029); points6.push_back(temp_p);}
            {geo_point_t temp_p(39.636, -910.886, 2027.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(39.636, -917.814, 2027.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(35.232, -928.206, 2027.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(30.828, -934.269, 2029); points6.push_back(temp_p);}
            {geo_point_t temp_p(26.424, -939.465, 2026); points6.push_back(temp_p);}
            {geo_point_t temp_p(26.424, -948.125, 2029); points6.push_back(temp_p);}
            {geo_point_t temp_p(22.02, -951.589, 2029); points6.push_back(temp_p);}
            {geo_point_t temp_p(22.02, -961.981, 2029); points6.push_back(temp_p);}
            {geo_point_t temp_p(17.616, -965.445, 2029); points6.push_back(temp_p);}
            {geo_point_t temp_p(15.414, -968.91, 2029); points6.push_back(temp_p);}
            {geo_point_t temp_p(13.212, -973.24, 2030.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(13.212, -981.034, 2032); points6.push_back(temp_p);}
            {geo_point_t temp_p(13.212, -986.23, 2029); points6.push_back(temp_p);}
            {geo_point_t temp_p(8.808, -996.622, 2029); points6.push_back(temp_p);}
            {geo_point_t temp_p(4.404, -1000.95, 2033.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(4.404, -1007.88, 2033.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(2.202, -1018.27, 2033.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(0, -1020.87, 2035); points6.push_back(temp_p);}
            {geo_point_t temp_p(-4.404, -1022.6, 2038); points6.push_back(temp_p);}
            {geo_point_t temp_p(-8.808, -1030.4, 2036.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-8.808, -1037.33, 2036.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-13.212, -1040.79, 2036.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-13.212, -1047.72, 2036.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-17.616, -1050.32, 2038); points6.push_back(temp_p);}
            {geo_point_t temp_p(-22.02, -1055.51, 2041); points6.push_back(temp_p);}
            {geo_point_t temp_p(-22.02, -1063.31, 2039.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-22.02, -1073.7, 2039.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-26.424, -1084.09, 2039.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-28.626, -1091.02, 2039.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-33.03, -1094.48, 2039.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-33.03, -1104.88, 2039.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-37.434, -1115.27, 2039.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-41.838, -1117, 2039.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-44.04, -1122.2, 2039.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-44.04, -1129.12, 2039.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-46.242, -1136.05, 2039.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-50.646, -1142.98, 2042.5); points6.push_back(temp_p);}
            {geo_point_t temp_p(-50.646, -1153.37, 2042.5); points6.push_back(temp_p);}

            std::cout << "Test: " << points6.size() << " " << points6.at(0) << " " << points6.at(points6.size()-1) << std::endl;
            

            live_clusters.at(i)->organize_points_path_vec(points6,0.6*units::cm);
            std::cout << "Test: " << points6.size() << " " << points6.at(0) << " " << points6.at(points6.size()-1) << std::endl;

            live_clusters.at(i)->organize_path_points(points6,0.6*units::cm);
            std::cout << "Test: " << points6.size() << " " << points6.at(0) << " " << points6.at(points6.size()-1) << std::endl;

            live_clusters.at(i)->examine_graph(true);
        }
    }
}