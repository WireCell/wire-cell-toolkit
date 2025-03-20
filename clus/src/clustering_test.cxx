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


void WireCell::PointCloud::Facade::clustering_test(
    Grouping& live_grouping,
    const Grouping& dead_grouping,
    cluster_set_t& cluster_connected_dead,
    const IDetectorVolumes::pointer dv

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
        const auto& wpids = cluster->wpids();
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
            auto& kd3d = cluster->kd3d();
            SPDLOG_INFO("CTest Cluster {} kd3d.ndim() {} kd3d.npoints() {}", iclus, kd3d.ndim(), kd3d.npoints());
            auto& kd2dp0 = cluster->kd2d(0, WirePlaneId(kAllLayers, 0, 0));
            SPDLOG_INFO("CTest Cluster {} kd2dp0.ndim() {} kd2dp0.npoints() {}", iclus, kd2dp0.ndim(), kd2dp0.npoints());
            auto& kd2dp0_a0f1 = cluster->kd2d(0, WirePlaneId(kAllLayers, 0, 1));
            SPDLOG_INFO("CTest Cluster {} kd2dp0_a0f1.ndim() {} kd2dp0_a0f1.npoints() {}", iclus, kd2dp0_a0f1.ndim(), kd2dp0_a0f1.npoints());

            {
                auto& sv3d = cluster->sv3d();
                const auto fpc = sv3d.flat_pc("3d", {"uwire_index"});
                SPDLOG_INFO("CTest Cluster {} sv3d.keys().size() {} sv3d.size_major() {}",
                            iclus, fpc.keys().size(), fpc.size_major());
                const auto& uwire_index = fpc.get("uwire_index")->elements<int>();
                SPDLOG_INFO("CTest Cluster {} two calls uwire_index[0] {}", iclus, uwire_index[0]);
            }

            {
                const auto uwire_index = cluster->points_property<int>("uwire_index");
                SPDLOG_INFO("CTest Cluster {} one call uwire_index[0] {}", iclus, uwire_index[0]);
            }

            const auto x = cluster->points_property<double>("x");
            const auto y = cluster->points_property<double>("y");
            const auto z = cluster->points_property<double>("z");
            const auto wpid_ident = cluster->points_property<int>("wpid");
            for (size_t ipt=0; ipt!=x.size(); ipt++) {
                WirePlaneId wpid(wpid_ident[ipt]);
                SPDLOG_INFO("CTest Cluster {} wpid {} x {} y {} z {}", iclus, wpid.name(), x[ipt], y[ipt], z[ipt]);
                break; // only one point
            }
            break; // only one cluster
        }
    }

}
#pragma GCC diagnostic pop
