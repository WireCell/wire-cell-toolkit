#include <WireCellClus/ClusteringFuncs.h>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;

void WireCell::PointCloud::Facade::clustering_examine_x_boundary(
    Grouping& live_grouping,
    const double low_limit,
    const double high_limit
    )
{
    std::vector<Cluster *> live_clusters = live_grouping.children();  // copy
    const auto &tp = live_grouping.get_params();
    // this is for 4 time slices
    double time_slice_width = tp.nticks_live_slice * tp.tick_drift;

    // std::vector<PR3DCluster *> new_clusters;
    // std::vector<PR3DCluster *> del_clusters;

    for (size_t i = 0; i != live_clusters.size(); i++) {
        Cluster *cluster = live_clusters.at(i);
        // only examine big clusters ...
        if (cluster->get_length() > 5 * units::cm && cluster->get_length() < 150 * units::cm) {
            // cluster->Create_point_cloud();
            std::unordered_map<int, Cluster*> id2clusters = cluster->examine_x_boundary(low_limit, high_limit);
            // if (clusters.size() != 0) {
            //     del_clusters.push_back(cluster);
            //     std::copy(clusters.begin(), clusters.end(), std::back_inserter(new_clusters));
            // }
        }
    }

    // for (auto it = new_clusters.begin(); it != new_clusters.end(); it++) {
    //     PR3DCluster *ncluster = (*it);
    //     // ncluster->Create_point_cloud();
    //     std::vector<int> range_v1 = ncluster->get_uvwt_range();
    //     double length_1 = sqrt(2. / 3. *
    //                                (pow(pitch_u * range_v1.at(0), 2) + pow(pitch_v * range_v1.at(1), 2) +
    //                                 pow(pitch_w * range_v1.at(2), 2)) +
    //                            pow(time_slice_width * range_v1.at(3), 2));
    //     cluster_length_map[ncluster] = length_1;
    //     live_clusters.push_back(ncluster);
    // }

    // for (auto it = del_clusters.begin(); it != del_clusters.end(); it++) {
    //     PR3DCluster *ocluster = (*it);
    //     cluster_length_map.erase(ocluster);
    //     live_clusters.erase(find(live_clusters.begin(), live_clusters.end(), ocluster));
    //     delete ocluster;
    // }
}