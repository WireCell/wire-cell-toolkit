#include "WireCellClus/IClusteringMethod.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"

class ClusteringExamineXBoundary;
WIRECELL_FACTORY(ClusteringExamineXBoundary, ClusteringExamineXBoundary,
                 WireCell::IConfigurable, WireCell::Clus::IClusteringMethod)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;

static void clustering_examine_x_boundary(
    Grouping& live_grouping, 
    IDetectorVolumes::pointer dv,
    const Tree::Scope& scope
    );

class ClusteringExamineXBoundary : public IConfigurable, public Clus::IClusteringMethod, private NeedDV, private NeedScope {
public:
    ClusteringExamineXBoundary() {}
    virtual ~ClusteringExamineXBoundary() {}

    void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedScope::configure(config);
    }

    void clustering(Grouping& live_clusters, Grouping&) const {
        clustering_examine_x_boundary(live_clusters, m_dv, m_scope);
    }
    
};


// This function only handles Single APA/Face!
static void clustering_examine_x_boundary(
    Grouping& live_grouping, 
    const IDetectorVolumes::pointer dv,
    const Tree::Scope& scope
    )
{
    // Check that live_grouping has less than one wpid
    if (live_grouping.wpids().size() > 1) {
        for (const auto& wpid : live_grouping.wpids()) {
            std::cout << "Live grouping wpid: " << wpid.name() << std::endl;
        }
        raise<ValueError>("Live %d > 1", live_grouping.wpids().size());
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
        if (!cluster->get_scope_filter(scope)) continue;
        
        if (cluster->get_default_scope().hash() != scope.hash()) {
            cluster->set_default_scope(scope);
            // std::cout << "Test: Set default scope: " << pc_name << " " << coords[0] << " " << coords[1] << " " << coords[2] << " " << cluster->get_default_scope().hash() << " " << scope.hash() << std::endl;
        }
        // only examine big clusters ...
        if (cluster->get_length() > 5 * units::cm && cluster->get_length() < 150 * units::cm) {
            // cluster->Create_point_cloud();
            // std::cout << "Cluster " << i << " old pointer " << cluster << " nchildren " << cluster->nchildren() << std::endl;
            auto b2groupid = cluster->examine_x_boundary(FV_xmin - FV_xmin_margin, FV_xmax + FV_xmax_margin);
            if (b2groupid.empty()) {
                continue;
            }

            // Perform separation
            auto scope_transform = cluster->get_scope_transform(scope);
            auto id2clusters = live_grouping.separate(cluster, b2groupid, true);
            // Apply the scope filter settings to all new clusters
            for (auto& [id, new_cluster] : id2clusters) {
                new_cluster->set_default_scope(scope);
                new_cluster->set_scope_filter(scope, true);
                new_cluster->set_scope_transform(scope, scope_transform);
            }

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


    // {
    //     auto live_clusters = live_grouping.children(); // copy
    //      // Process each cluster
    //      for (size_t iclus = 0; iclus < live_clusters.size(); ++iclus) {
    //          Cluster* cluster = live_clusters.at(iclus);
    //          auto& scope = cluster->get_default_scope();
    //          std::cout << "Test: " << iclus << " " << cluster->nchildren() << " " << scope.pcname << " " << scope.coords[0] << " " << scope.coords[1] << " " << scope.coords[2] << " " << cluster->get_scope_filter(scope)<< " " << cluster->get_center() << std::endl;
    //      }
    //    }

    // set cluster id ... 
    int cluster_id = 1;
    for (auto* cluster : live_grouping.children()) {
        cluster->set_cluster_id(cluster_id++);
    }
    
}
