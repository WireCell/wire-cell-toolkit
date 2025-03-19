#include <WireCellClus/ClusteringFuncs.h>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;

void WireCell::PointCloud::Facade::clustering_examine_x_boundary(
    Grouping& live_grouping, const IDetectorVolumes::pointer dv
    )
{
    // Check that live_grouping has exactly one wpid
	if (live_grouping.wpids().size() != 1 ) {
		throw std::runtime_error("Live or Dead grouping must have exactly one wpid");
	}

    std::vector<Cluster *> live_clusters = live_grouping.children();  // copy
    // sort the clusters by length using a lambda function
    // std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster *cluster1, const Cluster *cluster2) {
    //     return cluster1->get_length() > cluster2->get_length();
    // });

    // const auto &tp = live_grouping.get_params();
    // this is for 4 time slices
    // double time_slice_width = tp.nticks_live_slice * tp.tick_drift;


    // std::cout << "Test: " << tp.FV_xmin << " " << tp.FV_xmax << " " << tp.FV_xmin_margin << " " << tp.FV_xmax_margin << std::endl;
    // std::cout << "Test: " << dv->metadata(*live_grouping.wpids().begin())["FV_xmin"].asDouble() << " " << dv->metadata(*live_grouping.wpids().begin())["FV_xmax"].asDouble() << " " << dv->metadata(*live_grouping.wpids().begin())["FV_xmin_margin"].asDouble() << " " << dv->metadata(*live_grouping.wpids().begin())["FV_xmax_margin"].asDouble() << std::endl;

    double FV_xmin = dv->metadata(*live_grouping.wpids().begin())["FV_xmin"].asDouble() ;
    double FV_xmax = dv->metadata(*live_grouping.wpids().begin())["FV_xmax"].asDouble() ;
    double FV_xmin_margin = dv->metadata(*live_grouping.wpids().begin())["FV_xmin_margin"].asDouble() ;
    double FV_xmax_margin = dv->metadata(*live_grouping.wpids().begin())["FV_xmax_margin"].asDouble() ;

    // std::vector<PR3DCluster *> new_clusters;
    // std::vector<PR3DCluster *> del_clusters;

    for (size_t i = 0; i != live_clusters.size(); i++) {
        Cluster *cluster = live_clusters.at(i);
        // only examine big clusters ...
        if (cluster->get_length() > 5 * units::cm && cluster->get_length() < 150 * units::cm) {
            // cluster->Create_point_cloud();
            // std::cout << "Cluster " << i << " old pointer " << cluster << " nchildren " << cluster->nchildren() << std::endl;
            auto b2groupid = cluster->examine_x_boundary(FV_xmin - FV_xmin_margin, FV_xmax + FV_xmax_margin);
            if (b2groupid.empty()) {
                continue;
            }
            live_grouping.separate(cluster, b2groupid, true);
            assert(cluster == nullptr);


            // std::cout << "Cluster " << i << " is seperated into " << id2clusters.size() << " clusters" << std::endl;
            // for (auto [id, ncluster] : id2clusters) {
            //     std::cout << "id " << id << " new pointer " << ncluster << " nchildren " << ncluster->nchildren() << std::endl;
            // }
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
