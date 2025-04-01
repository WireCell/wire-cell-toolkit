#include <WireCellClus/ClusteringFuncs.h>
#include "WireCellUtil/ExecMon.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;

// This is a test function, not used in clustering
void WireCell::PointCloud::Facade::clustering_test(
    Grouping& live_grouping,
    const Grouping& dead_grouping,
    cluster_set_t& cluster_connected_dead,
    const IDetectorVolumes::pointer dv,
    const std::string& pc_name,                        // point cloud name
    const std::vector<std::string>& coords            // coordinate names
)
{
    using spdlog::debug;
  
    // form map from dead to set of live clusters ...
    std::map<const Cluster*, std::vector<const Cluster*>> dead_live_cluster_mapping;
    std::vector<const Cluster*> dead_cluster_order;
    std::map<const Cluster*, std::vector<std::vector<const Blob*>>> dead_live_mcells_mapping;

    std::vector<Cluster*> live_clusters = live_grouping.children(); // copy
    std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster *cluster1, const Cluster *cluster2) {
        return cluster1->get_length() > cluster2->get_length();
    });

    /// TEST: wpid
    for (const auto& wpid : live_grouping.wpids()) {
        SPDLOG_INFO("CTest live_grouping wpid {}", wpid.name());
    }
    for (size_t iclus = 0; iclus != live_clusters.size(); iclus++) {
        Cluster* cluster = live_clusters.at(iclus);
        const auto& wpids = cluster->wpids_blob();
        for (size_t i=0; i != wpids.size(); i++) {
            const auto& wpid = wpids.at(i);
            SPDLOG_INFO("CTest Cluster {} i {} name {}", iclus, i, wpid.name());
            break;
        }
        for (size_t iblob = 0; iblob != cluster->children().size(); iblob++) {
            const auto* blob = cluster->children().at(iblob);
            SPDLOG_INFO("CTest Cluster {} Blob {} blob->wpid().name() {}", iclus, iblob, blob->wpid().name());
            break;
        }
        break;
    }

    /// TEST: IDetectorVolumes
    {
        std::vector<WirePlaneLayer_t> layers = {kUlayer, kVlayer, kWlayer};
        for (const auto& gwpid : live_grouping.wpids()) {
            for (const auto& layer : layers) {
                WirePlaneId wpid(layer, gwpid.face(), gwpid.apa());
                int face_dirx = dv->face_dirx(wpid);
                Vector wire_direction = dv->wire_direction(wpid);
                double angle = std::atan2(wire_direction.z(), wire_direction.y());
                Vector pitch_vector = dv->pitch_vector(wpid);
                SPDLOG_INFO("CTest wpid.name {} face_dirx {} wire_direction {} angle rad:{}  deg:{} pitch_vector {}", wpid.name(), face_dirx, wire_direction, angle, angle*180/3.1415926, pitch_vector);
            }
        }
        // metadata
        {
            WirePlaneId all(0);
            SPDLOG_INFO("FV_xmax {}", dv->metadata(all)["FV_xmax"].asDouble());
            WirePlaneId a0f0pA(kAllLayers, 0, 0);
            SPDLOG_INFO("FV_xmax {}", dv->metadata(a0f0pA)["FV_xmax"].asDouble());
            Json::FastWriter fastWriter;
            SPDLOG_INFO("metadata(a0f0pA): {}", fastWriter.write(dv->metadata(a0f0pA)));
        }
    }

    /// TEST: points: wpid, merge 3d/2d
    {
        for (size_t iclus = 0; iclus != live_clusters.size(); iclus++) {
            Cluster* cluster = live_clusters.at(iclus);

            // TEST: kd2d with wpid
            {
                // expecting:
                // kd3d.ndim() 3 kd3d.npoints() 4248 (non-zero)
                // kd2dp0.ndim() 2 kd2dp0.npoints() 4248 (same as 3D)
                // kd2dp0_a0f1.ndim() 2 kd2dp0_a0f1.npoints() 0 (if we do not have a0f1)
                auto& kd3d = cluster->kd3d();
                SPDLOG_INFO("CTest Cluster {} kd3d.ndim() {} kd3d.npoints() {}", iclus, kd3d.ndim(), kd3d.npoints());
                auto& kd2dp0 = cluster->kd2d(0, 0, 0);
                SPDLOG_INFO("CTest Cluster {} kd2dp0.ndim() {} kd2dp0.npoints() {}", iclus, kd2dp0.ndim(), kd2dp0.npoints());
                auto& kd2dp0_a0f1 = cluster->kd2d(0, 0 , 0);
                SPDLOG_INFO("CTest Cluster {} kd2dp0_a0f1.ndim() {} kd2dp0_a0f1.npoints() {}", iclus, kd2dp0_a0f1.ndim(), kd2dp0_a0f1.npoints());
            }

            // TEST: flat_pc and points_property
            {
                auto& sv3d = cluster->sv3d();
                const auto fpc = sv3d.flat_pc("3d", {"uwire_index"});
                SPDLOG_INFO("CTest Cluster {} sv3d.keys().size() {} sv3d.size_major() {}",
                            iclus, fpc.keys().size(), fpc.size_major());
                const auto& uwire_index_fpc = fpc.get("uwire_index")->elements<int>();
                SPDLOG_INFO("CTest Cluster {} flat_pc uwire_index[0] {}", iclus, uwire_index_fpc[0]);
                const auto uwire_index_pp = cluster->points_property<int>("uwire_index");
                SPDLOG_INFO("CTest Cluster {} points_property uwire_index[0] {}", iclus, uwire_index_pp[0]);
            }

            // TEST: points_property
            {
                const auto x = cluster->points_property<double>("x");
                const auto y = cluster->points_property<double>("y");
                const auto z = cluster->points_property<double>("z");
                const auto wpid_ident = cluster->points_property<int>("wpid");
                for (size_t ipt=0; ipt!=x.size(); ipt++) {
                    WirePlaneId wpid(wpid_ident[ipt]);
                    SPDLOG_INFO("CTest Cluster {} wpid {} x {} y {} z {}", iclus, wpid.name(), x[ipt], y[ipt], z[ipt]);
                    break; // only one point
                }
            }

            break; // only one cluster
        }
    }

    /// TEST: DynamicPointCloud
    {
        // Get all the wire plane IDs from the grouping
        const auto& wpids = live_grouping.wpids();
        // Key: pair<APA, face>, Value: drift_dir, angle_u, angle_v, angle_w
        std::map<WirePlaneId , std::tuple<geo_point_t, double, double, double>> wpid_params;
        std::set<int> apas;
    
        std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_u_index; 
        std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_v_index; 
        std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_w_index; 
    
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

        auto [drift_dir, angle_u, angle_v, angle_w] = extract_geometry_params(live_grouping, dv);
        auto dpc = std::make_shared<DynamicPointCloud>(wpid_params);
        // auto dpcl = std::make_shared<DynamicPointCloudLegacy>(angle_u, angle_v, angle_w);
        double extending_dis = 50 * units::cm;
        double angle = 7.5;
        double loose_dis_cut = 7.5 * units::cm;
        geo_point_t dir1(1, 0, 0);
        for (size_t iclus = 0; iclus != live_clusters.size(); iclus++) {
            Cluster* cluster = live_clusters.at(iclus);
            const auto test_point = cluster->point3d(0);
            std::pair<geo_point_t, geo_point_t> extreme_points = cluster->get_two_extreme_points();
            // dpcl->add_points(cluster, extreme_points.first, dir1, extending_dis * 3, 1.2 * units::cm, angle);
            dpc->add_points(make_points_linear_extrapolation(cluster, extreme_points.first, dir1, extending_dis * 3, 1.2 * units::cm, angle, dv, wpid_params));
            // const auto results_legacy = dpcl->get_2d_points_info(test_point, loose_dis_cut, 0);
            const auto results = dpc->get_2d_points_info(test_point, loose_dis_cut, 0, 0, 0);
            // SPDLOG_INFO("CTest Cluster {} results_legacy.size() {} results.size() {}", iclus, results_legacy.size(), results.size());
        }
        for (size_t iclus = 0; iclus != live_clusters.size(); iclus++) {
            Cluster* cluster = live_clusters.at(iclus);
            std::pair<geo_point_t, geo_point_t> extreme_points = cluster->get_two_extreme_points();
            // dpcl->add_points(cluster,0);
            dpc->add_points(make_points_cluster(cluster, wpid_params));
            // SPDLOG_INFO("CTest dpcl->get_num_points() {} dpc->get_points().size() {}",
                        // dpcl->get_num_points(), dpc->get_points().size());
            const auto dir_hough = dpc->vhough_transform(extreme_points.first, extending_dis);
            SPDLOG_INFO("CTest Cluster {} dir_hough {} ", iclus, dir_hough);
            // const auto dir_hough_legacy = dpcl->vhough_transform(extreme_points.first, extending_dis);
        }
    }

}
#pragma GCC diagnostic pop
