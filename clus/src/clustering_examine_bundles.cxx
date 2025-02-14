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

// #define __DEBUG__
#ifdef __DEBUG__
#define LogDebug(x) std::cout << "[yuhw]: " << __LINE__ << " : " << x << std::endl
#else
#define LogDebug(x)
#endif

void WireCell::PointCloud::Facade::clustering_examine_bundles(Grouping& live_grouping,
    const bool use_ctpc)
{
    std::cout << "Test Examine Bundles" << std::endl;

    std::vector<Cluster *> live_clusters = live_grouping.children();
    // for (size_t i = 0; i != live_clusters.size(); i++) {
    //    auto blobs = live_clusters.at(i)->kd_blobs();
    //    int nblobs = blobs.size();

    //    if(nblobs > 10){
    //     std::cout << "Test: " << nblobs << " " <<  std::endl;
    //     auto local_pcs = live_clusters.at(i)->local_pcs();
    //     for (auto it = local_pcs.begin(); it !=local_pcs.end(); it++){
    //         auto keys = it->second.keys();
    //         for (auto it1 = keys.begin(); it1 != keys.end(); it1++){
    //             std::cout << "Test: " << it->first << " " << *it1 << std::endl;
    //         }
    //     }
    //     auto flash = live_clusters.at(i)->get_scalar<int>("flash");
    //     std::cout << "Test: Flash: " << flash << std::endl;
    //    }
    // }


    for (size_t i=0;i!=live_clusters.size();i++){
        // if there is a cc component, record the main cluster as id of the blobs???

        // currently reset the cc component (todo: find the main component)

        // do the examine graph
        auto b2groupid = live_clusters.at(i)->examine_graph(true);
    
        // find the main cluster and set it to the cc tree ...

    }

}