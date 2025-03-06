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
    cluster_set_t& cluster_connected_dead            // in/out
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

    for (size_t iclus = 0; iclus != live_clusters.size(); iclus++) {
        Cluster* cluster = live_clusters.at(iclus);
        for (size_t iblob = 0; iblob != cluster->children().size(); iblob++) {
            const auto* blob = cluster->children().at(iblob);
            SPDLOG_INFO("CTest Cluster {} Blob {} blob->wpid().name() {}", iclus, iblob, blob->wpid().name());
            break;
        }
        break;
    }
}
#pragma GCC diagnostic pop
