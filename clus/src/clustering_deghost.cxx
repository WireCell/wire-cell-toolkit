#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"

class ClusteringDeghost;
WIRECELL_FACTORY(ClusteringDeghost, ClusteringDeghost,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)


using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;

// bool Cluster::construct_skeleton(IDetectorVolumes::pointer dv, IPCTransformSet::pointer pcts, const bool use_ctpc)
// {

static std::pair<size_t, size_t> skeleton_points_hilo(const Cluster& cluster)
{
    geo_point_t highest_wcp = cluster.point3d(0);
    geo_point_t lowest_wcp = cluster.point3d(0);
    size_t highest_index = 0;
    size_t lowest_index = 0;

    geo_point_t main_dir = cluster.get_pca().axis.at(0);
    main_dir = main_dir.norm();
    geo_point_t center = cluster.get_pca().center;
    geo_point_t temp_pt(highest_wcp.x() - center.x(), highest_wcp.y() - center.y(), highest_wcp.z() - center.z());
    double highest_value = temp_pt.dot(main_dir);
    double lowest_value = highest_value;

    for (int i = 1; i < cluster.npoints(); i++) {
        temp_pt.set(cluster.point3d(i).x() - center.x(),
                    cluster.point3d(i).y() - center.y(),
                    cluster.point3d(i).z() - center.z());
        double value = temp_pt.dot(main_dir);
        if (value > highest_value) {
            highest_value = value;
            highest_wcp = cluster.point3d(i);
            highest_index = i;
        }
        else if (value < lowest_value) {
            lowest_value = value;
            lowest_wcp = cluster.point3d(i);
            lowest_index = i;
        }
    }
    return std::make_pair(highest_index, lowest_index);
}

static std::vector<size_t> get_path_wcps(const Cluster& cluster, 
                                         IDetectorVolumes::pointer dv,
                                         IPCTransformSet::pointer pcts,
                                         bool use_ctpc)
{
    auto [hi, lo] = skeleton_points_hilo(cluster);
    if (use_ctpc) {
        return cluster.graph_algorithms("ctpc", dv, pcts).shortest_path(hi, lo);
    }
    else {
        return cluster.graph_algorithms().shortest_path(hi, lo);
    }
}



static void clustering_deghost(Grouping& live_grouping,
                               IDetectorVolumes::pointer dv,
                               IPCTransformSet::pointer pcts,
                               const Tree::Scope& scope,
                               const bool use_ctpc,
                               double length_cut = 0);

class ClusteringDeghost : public IConfigurable, public Clus::IEnsembleVisitor, private NeedDV, private NeedPCTS, private NeedScope {
public:
    ClusteringDeghost() {}
    virtual ~ClusteringDeghost() {}

    void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config);
        NeedScope::configure(config);
        
        use_ctpc_ = get(config, "use_ctpc", true);
        length_cut_ = get(config, "length_cut", 0);
    }
    virtual Configuration default_configuration() const {
        Configuration cfg;
        return cfg;
    }
    
    void visit(Ensemble& ensemble) const {
        auto& live = *ensemble.with_name("live").at(0);
        clustering_deghost(live, m_dv, m_pcts, m_scope,  use_ctpc_, length_cut_);
    }
    
private:
    double use_ctpc_{true};
    double length_cut_{0};
};


// The original developers do not care.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

// #define __DEBUG__
#ifdef __DEBUG__
#define LogDebug(x) std::cout << "[yuhw]: " << __LINE__ << " : " << x << std::endl
#else
#define LogDebug(x)
#endif

// This can handle entire APA (including all faces) data
static void clustering_deghost(
    Grouping& live_grouping,
    IDetectorVolumes::pointer dv,
    IPCTransformSet::pointer pcts,
    const Tree::Scope& scope,
    const bool use_ctpc, double length_cut)
{
    // Get all the wire plane IDs from the grouping
    const auto& wpids = live_grouping.dv_wpids();
    // Key: pair<APA, face>, Value: drift_dir, angle_u, angle_v, angle_w
    std::map<WirePlaneId , std::tuple<geo_point_t, double, double, double>> wpid_params;
    std::set<int> apas;

    std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_u_index; 
    std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_v_index; 
    std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_w_index; 

    // 
    //  NOTE, most of this can be replaced by a couple of function calls from DetUtils.h
    //

    for (const auto& wpid : wpids) {
        int apa = wpid.apa();
        int face = wpid.face();
        apas.insert(apa);

        // Create wpids for all three planes with this APA and face
        WirePlaneId wpid_u(kUlayer, face, apa);
        WirePlaneId wpid_v(kVlayer, face, apa);
        WirePlaneId wpid_w(kWlayer, face, apa);
     
        // Get drift direction based on face orientation
        int face_dirx = dv->face_dirx(wpid_u);
        geo_point_t drift_dir(face_dirx, 0, 0);
        
        // Get wire directions for all planes
        Vector wire_dir_u = dv->wire_direction(wpid_u);
        Vector wire_dir_v = dv->wire_direction(wpid_v);
        Vector wire_dir_w = dv->wire_direction(wpid_w);

        // Calculate angles
        double angle_u = std::atan2(wire_dir_u.z(), wire_dir_u.y());
        double angle_v = std::atan2(wire_dir_v.z(), wire_dir_v.y());
        double angle_w = std::atan2(wire_dir_w.z(), wire_dir_w.y());

        wpid_params[wpid] = std::make_tuple(drift_dir, angle_u, angle_v, angle_w);


        af_dead_u_index[apa][face] = live_grouping.get_dead_winds(apa, face, 0);
        af_dead_v_index[apa][face] = live_grouping.get_dead_winds(apa, face, 1);
        af_dead_w_index[apa][face] = live_grouping.get_dead_winds(apa, face, 2);
    }

    if (apas.size() > 1) {
        raise<ValueError>("apas.size() %d > 1", apas.size());
    }
  

    std::vector<Cluster *> live_clusters = live_grouping.children();  // copy

    // sort the clusters by length using a lambda function
    std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster *cluster1, const Cluster *cluster2) {
        return cluster1->get_length() > cluster2->get_length();
    });

    //
    // NOTE: these (and above code) can be replaced with a single call to
    // make_dynamicpointcloud(dv) from DetUtils.
    //

    // auto global_point_cloud_legacy = std::make_shared<DynamicPointCloudLegacy>(angle_u, angle_v, angle_w);
    auto global_point_cloud = std::make_shared<DynamicPointCloud>(wpid_params);
    // auto global_skeleton_cloud = std::make_shared<DynamicPointCloudLegacy>(angle_u, angle_v, angle_w);
    auto global_skeleton_cloud = std::make_shared<DynamicPointCloud>(wpid_params);
    // replace with the new DynamicPointCloud class

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

    bool flag_first = true;
    for (size_t i = 0; i != live_clusters.size(); i++) {
        if (live_clusters.at(i)->get_default_scope().hash() != scope.hash()) {
            live_clusters.at(i)->set_default_scope(scope);
            // std::cout << "Test: Set default scope: " << pc_name << " " << coords[0] << " " << coords[1] << " " << coords[2] << " " << cluster->get_default_scope().hash() << " " << scope.hash() << std::endl;
        }
        
        // if not within the scope filter, nor processing ...
        if (!live_clusters.at(i)->get_scope_filter(scope)) continue;
        
        if (flag_first) {
            // fill anyway ...
            // live_clusters.at(i)->Create_point_cloud();
            // global_point_cloud_legacy->add_points(live_clusters.at(i), 0);
            global_point_cloud->add_points(make_points_cluster(live_clusters.at(i), wpid_params, true));
            if (live_clusters.at(i)->get_length() >
                30 * units::cm) {  // should be the default for most of them ...
                const auto& path_wcps = get_path_wcps(*live_clusters.at(i), dv, pcts, use_ctpc);
                // global_skeleton_cloud->add_points(live_clusters.at(i), 1);
                global_skeleton_cloud->add_points(make_points_cluster_skeleton(live_clusters.at(i), dv, wpid_params, path_wcps, true));
            }
            else {
                // global_skeleton_cloud->add_points(live_clusters.at(i), 0);
                global_skeleton_cloud->add_points(make_points_cluster(live_clusters.at(i), wpid_params, true));
            }
            flag_first = false;
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
                    geo_point_t test_point = cluster->point3d(j);
                    auto test_wpid = cluster->wire_plane_id(j);
                    bool flag_dead = false;

                    #ifdef __DEBUG__
                    if (num_total_points == 134) {
                        std::cout << "point: " << j << " " << test_point << " " << winds[0][j] << " " << winds[1][j] << " " << winds[2][j] << std::endl;
                    }
                    #endif

                    
                    auto& dead_u_index = af_dead_u_index.at(test_wpid.apa()).at(test_wpid.face());
                    if (dead_u_index.find(winds[0][j]) != dead_u_index.end()) {
                        if (cluster->point3d_raw(j).x() >= dead_u_index[winds[0][j]].first &&
                            cluster->point3d_raw(j).x() <= dead_u_index[winds[0][j]].second) {
                            flag_dead = true;
                        }
                    }

                    

                    if (!flag_dead) {
                        std::tuple<double, const Cluster *, size_t> results =
                            global_point_cloud->get_closest_2d_point_info(test_point, 0, test_wpid.face(), test_wpid.apa()); 

                

                        if (std::get<0>(results) <= dis_cut / 3.) {
                            if (map_cluster_num[0].find(std::get<1>(results)) == map_cluster_num[0].end()) {
                                map_cluster_num[0][std::get<1>(results)] = 1;
                            }
                            else {
                                map_cluster_num[0][std::get<1>(results)]++;
                            }
                        }
                        else {
                            results = global_skeleton_cloud->get_closest_2d_point_info(test_point, 0, test_wpid.face(), test_wpid.apa()); 

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
                    auto& dead_v_index = af_dead_v_index.at(test_wpid.apa()).at(test_wpid.face());

                    if (dead_v_index.find(winds[1][j]) != dead_v_index.end()) {
                        #ifdef __DEBUG__
                        if (num_total_points == 134) {
                            std::cout << "dead_v_index: " << winds[1][j] << " " << dead_v_index[winds[1][j]].first << " " << dead_v_index[winds[1][j]].second << std::endl;
                        }
                        #endif
                        if (cluster->point3d_raw(j).x() >= dead_v_index[winds[1][j]].first &&
                            cluster->point3d_raw(j).x() <= dead_v_index[winds[1][j]].second) {
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
                            global_point_cloud->get_closest_2d_point_info(test_point, 1, test_wpid.face(), test_wpid.apa()); 

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
                            // results = global_skeleton_cloud->get_closest_2d_point_info(test_point, 1);
                            results = global_skeleton_cloud->get_closest_2d_point_info(test_point, 1, test_wpid.face(), test_wpid.apa()); 

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
                    auto& dead_w_index = af_dead_w_index.at(test_wpid.apa()).at(test_wpid.face());

                    if (dead_w_index.find(winds[2][j]) != dead_w_index.end()) {
                        if (cluster->point3d_raw(j).x() >= dead_w_index[winds[2][j]].first &&
                            cluster->point3d_raw(j).x() <= dead_w_index[winds[2][j]].second) {
                            flag_dead = true;
                        }
                    }

                    if (!flag_dead) {
                        std::tuple<double, const Cluster *, size_t> results =
                            global_point_cloud->get_closest_2d_point_info(test_point, 2, test_wpid.face(), test_wpid.apa()); 

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
                            // results = global_skeleton_cloud->get_closest_2d_point_info(test_point, 2);
                            results = global_skeleton_cloud->get_closest_2d_point_info(test_point, 2, test_wpid.face(), test_wpid.apa()); 

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
                // LogDebug("num_total_points = " << num_total_points);
                // LogDebug("num_unique[0] = " << num_unique[0] << ", num_unique[1] = " << num_unique[1] << ", num_unique[2] = " << num_unique[2]);
                // LogDebug("num_dead[0] = " << num_dead[0] << ", num_dead[1] = " << num_dead[1] << ", num_dead[2] = " << num_dead[2]);

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
                    // LogDebug("pass the first cut " << num_total_points);

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
                    // LogDebug("max_value_u: " << max_value_u << ", max_value_v: " << max_value_v << ", max_value_w: " << max_value_w);

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
                    // global_point_cloud_legacy->add_points(live_clusters.at(i), 0);
                    global_point_cloud->add_points(make_points_cluster(live_clusters.at(i), wpid_params, true));
                    if (live_clusters.at(i)->get_length() > 30 * units::cm) {
                        const auto& path_wcps = get_path_wcps(*live_clusters.at(i), dv, pcts, use_ctpc);
                        // global_skeleton_cloud->add_points(live_clusters.at(i), 1);
                        global_skeleton_cloud->add_points(make_points_cluster_skeleton(live_clusters.at(i), dv, wpid_params, path_wcps, true ));
                    }
                }
            }
            else {
                // live_clusters.at(i)->Create_point_cloud();
                // global_point_cloud_legacy->add_points(live_clusters.at(i), 0);
                global_point_cloud->add_points(make_points_cluster(live_clusters.at(i), wpid_params, true));
                if (live_clusters.at(i)->get_length() > 30 * units::cm) {
                    const auto& path_wcps = get_path_wcps(*live_clusters.at(i), dv, pcts, use_ctpc);
                    // global_skeleton_cloud->add_points(live_clusters.at(i), 1);
                    global_skeleton_cloud->add_points(make_points_cluster_skeleton(live_clusters.at(i), dv, wpid_params, path_wcps, true));
                }
            }
        }
        // LogDebug("Cluster " << i << " " << live_clusters.at(i)->n_blobs() << " " << live_clusters.at(i)->npoints());
        // LogDebug("global_point_cloud: " << global_point_cloud->get_num_points() << " global_skeleton_cloud: " << global_skeleton_cloud->get_num_points());
    }


    auto new_clusters = merge_clusters(g, live_grouping);

    // remove clusters
    LogDebug("to_be_removed_clusters.size() = " << to_be_removed_clusters.size());
    for (auto live : to_be_removed_clusters) {
        // std::cout << "Remove cluster " << live->nchildren() << " " << live->get_length()/units::cm << std::endl;
        live_grouping.destroy_child(live);
        assert(live == nullptr);
    } 







//      {
//    auto live_clusters = live_grouping.children(); // copy
//     // Process each cluster
//     for (size_t iclus = 0; iclus < live_clusters.size(); ++iclus) {
//         Cluster* cluster = live_clusters.at(iclus);
//         auto& scope = cluster->get_default_scope();
//         std::cout << "Test: " << iclus << " " << cluster->nchildren() << " " << scope.pcname << " " << scope.coords[0] << " " << scope.coords[1] << " " << scope.coords[2] << " " << cluster->get_scope_filter(scope)<< " " << cluster->get_pca().center << std::endl;
//     }
//   }

}
