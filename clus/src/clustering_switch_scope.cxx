#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Spdlog.h"
#include "WireCellUtil/Logging.h"

class ClusteringSwitchScope;
WIRECELL_FACTORY(ClusteringSwitchScope, ClusteringSwitchScope,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;

// A named logger so these messages flow through WCT's configured sinks
// instead of relying on the process-global default logger.
static Log::logptr_t logger() {
    static Log::logptr_t l = Log::logger("clus.SwitchScope");
    return l;
}


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
            logger()->warn("ClusteringSwitchScope: no 'live' grouping found, skipping");
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
        logger()->warn("clustering_switch_scope: live grouping has no wpids, skipping");
        return;
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

        // Carry the flash-merge per-blob provenance across the rebuild:
        // separate() -> from() copies scalars/flags/scopes but NOT node-local
        // PCs, so the "perblob" real_cluster_id / real_cluster_main arrays
        // (persisted through the pctree tarball by the QL job's
        // save_real_cluster_id) would die here -- and they are what the Bee
        // writer's colors and TaggerCheckTGM main_component_mode="real" read.
        // Row-partition them by the same filter that partitions the blobs.
        // Absent arrays (legacy tarballs) => no-op.
        std::vector<int> src_rcid, src_rcmain;
        if (cluster->has_pcarray<int>("real_cluster_id", "perblob")) {
            auto sp = cluster->get_pcarray<int>("real_cluster_id", "perblob");
            src_rcid.assign(sp.begin(), sp.end());
        }
        if (cluster->has_pcarray<int>("real_cluster_main", "perblob")) {
            auto sp = cluster->get_pcarray<int>("real_cluster_main", "perblob");
            src_rcmain.assign(sp.begin(), sp.end());
        }

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
            // Re-attach the provenance rows of the blobs that went to this
            // part.  separate() keeps the blobs of each part in their source
            // order, so filtering the source array by the same group id
            // reproduces the parallel per-blob order.  Skip (fail open) on
            // any size mismatch.
            const int part_id = id;         // lambdas cannot capture a structured binding pre-C++20
            Cluster* part = new_cluster;
            auto carve = [&](const std::vector<int>& src, const char* aname) {
                if (src.empty() || src.size() != filter_results.size()) return;
                std::vector<int> sub;
                sub.reserve(src.size());
                for (size_t i = 0; i < src.size(); ++i) {
                    if (filter_results[i] == part_id) sub.push_back(src[i]);
                }
                if (sub.size() != (size_t)part->nchildren()) return;
                part->put_pcarray(sub, aname, "perblob");
            };
            carve(src_rcid, "real_cluster_id");
            carve(src_rcmain, "real_cluster_main");
        }
    }
}
