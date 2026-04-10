#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Spdlog.h"

class ClusteringSwitchScope;
WIRECELL_FACTORY(ClusteringSwitchScope, ClusteringSwitchScope,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;


static void clustering_switch_scope(
        Grouping& live_grouping,
        IPCTransformSet::pointer pcts,
        const std::string& correction_name
    );

class ClusteringSwitchScope : public IConfigurable, public Clus::IEnsembleVisitor, private NeedPCTS {
public:
    ClusteringSwitchScope() {}
    virtual ~ClusteringSwitchScope() {}

    void configure(const WireCell::Configuration& config) {
        NeedPCTS::configure(config);
        correction_name_ = convert<std::string>(config["correction_name"], "T0Correction");
    }

    void visit(Ensemble& ensemble) const {
        auto live_vec = ensemble.with_name("live");
        if (live_vec.empty()) {
            spdlog::warn("ClusteringSwitchScope: no 'live' grouping found, skipping");
            return;
        }
        clustering_switch_scope(*live_vec.at(0), m_pcts, correction_name_);
    }

private:
    std::string correction_name_{"T0Correction"};
};

static void clustering_switch_scope(
    Grouping& live_grouping,
    const Clus::IPCTransformSet::pointer pcts,
    const std::string& correction_name
)
{
    if (live_grouping.wpids().empty()) {
        throw std::runtime_error("Live grouping must have at least one wpid");
    }

    // Snapshot the cluster list before separation modifies the grouping.
    const std::vector<Cluster*> live_clusters = live_grouping.children();

    for (Cluster* cluster : live_clusters) {
        // Apply the named correction to each blob's 3d point cloud.
        // Returns a per-blob filter flag: 1 if any corrected point is inside the
        // active detector volume, 0 otherwise.
        // add_corrected_points() dispatches to the appropriate IPCTransform and
        // raises RuntimeError for unrecognised correction names.
        const std::vector<int> filter_results = cluster->add_corrected_points(pcts, correction_name);

        // Retrieve the correction scope registered by add_corrected_points().
        const auto correction_scope = cluster->get_scope(correction_name);

        // Set the correction scope as the cluster's default before separation so
        // that Grouping::separate() -> Cluster::from() propagates it to all children.
        cluster->set_default_scope(correction_scope);
        cluster->set_scope_transform(correction_scope, correction_name);

        // Split into two sub-clusters: id=0 (all blobs failed filter),
        // id=1 (at least one blob passed filter).
        // 'cluster' is destroyed by separate() (remove=true); do not use it after this line.
        auto separated_clusters = live_grouping.separate(cluster, filter_results, true);

        // from() (called inside separate()) already copied default_scope and
        // scope_transform to each separated cluster. Only the scope_filter is new.
        for (auto& [id, new_cluster] : separated_clusters) {
            new_cluster->set_scope_filter(correction_scope, id == 1);
        }
    }
}
