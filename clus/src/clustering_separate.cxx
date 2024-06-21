#include <WireCellClus/ClusteringFuncs.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;
void WireCell::PointCloud::Facade::clustering_separate(Grouping& live_grouping,
                                                       std::map<int, std::pair<double, double>>& dead_u_index,
                                                       std::map<int, std::pair<double, double>>& dead_v_index,
                                                       std::map<int, std::pair<double, double>>& dead_w_index)
{
    std::vector<Cluster*> live_clusters = live_grouping.children();  // copy
    for (size_t ilive = 0; ilive < live_clusters.size(); ++ilive) {
        const auto& live = live_clusters.at(ilive);
        live->Create_graph();
    }
}