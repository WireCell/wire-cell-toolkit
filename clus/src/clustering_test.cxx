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
            SPDLOG_DEBUG("metadata(a0f0pA): {}", fastWriter.write(dv->metadata(a0f0pA)));
        }
    }

}
#pragma GCC diagnostic pop
