#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/Logger.h"

#include <unordered_map>


class ClusteringExamineBundles;
WIRECELL_FACTORY(ClusteringExamineBundles, ClusteringExamineBundles,
                 WireCell::INamed, WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)


using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;


static void clustering_examine_bundles(
        Grouping& live_grouping,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts,
        const Tree::Scope& scope,
        const bool use_ctpc,
        const std::string& graph_name,
        const bool use_flash_t0,
        const double flash_t0_window,
        WireCell::Log::logptr_t log);

class ClusteringExamineBundles : public IConfigurable, public Clus::IEnsembleVisitor, public Aux::Logger, private NeedDV, private NeedPCTS, private NeedScope {
public:
    ClusteringExamineBundles()
        : Aux::Logger("ClusteringExamineBundles", "clus")
    {}
    virtual ~ClusteringExamineBundles() {}
    
    void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config);
        NeedScope::configure(config);

        // If false, then DV and PCTS are not needed.
        use_ctpc_ = get<bool>(config, "use_ctpc", use_ctpc_);
        graph_name_ = get<std::string>(config, "graph_name", graph_name_);
        use_flash_t0_ = get<bool>(config, "use_flash_t0", use_flash_t0_);
        flash_t0_window_ = get<double>(config, "flash_t0_window", flash_t0_window_);
    }

    void visit(Ensemble& ensemble) const {
        auto& live = *ensemble.with_name("live").at(0);
        clustering_examine_bundles(live, m_dv, m_pcts, m_scope, use_ctpc_, graph_name_,
                                   use_flash_t0_, flash_t0_window_, log);
    }

private:
    bool use_ctpc_{true};
    std::string graph_name_{"relaxed"};
    bool use_flash_t0_{false};
    double flash_t0_window_{80*units::ns};
};


// The original developers do not care.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

// #define __DEBUG__
#ifdef __DEBUG__
#define LogDebug(x) std::cout << "[yuhw]: " << __LINE__ << " : " << x << std::endl
#else
#define LogDebug(x)
#endif

// All APA Faces 
static void clustering_examine_bundles(
    Grouping& live_grouping,
    IDetectorVolumes::pointer dv,
    IPCTransformSet::pointer pcts,
    const Tree::Scope& scope,
    const bool use_ctpc,
    const std::string& graph_name,
    const bool use_flash_t0,
    const double flash_t0_window,
    WireCell::Log::logptr_t log)
{
    // std::cout << "Test Examine Bundles" << std::endl;

    // Flash-aware pre-step: collapse each matched-flash-time group of clusters
    // into a single cluster, so the per-cluster labeling below partitions one
    // cluster per flash group (main sub-component = -1), like clustering_isolated
    // but grouped by flash time instead of geometry.  merge_clusters() (see
    // ClusteringFuncs.cxx) carries the longest contributing member's flash onto
    // the merged cluster.  With use_flash_t0=false this is skipped entirely.
    if (use_flash_t0) {
        auto pre_clusters = live_grouping.children();
        auto flash_t0_group = assign_flash_t0_groups(pre_clusters, flash_t0_window);

        cluster_connectivity_graph_t g;
        std::unordered_map<int, int> ilive2desc;
        std::unordered_map<const Cluster*, int> map_cluster_index;
        for (size_t ilive = 0; ilive < pre_clusters.size(); ++ilive) {
            auto live = pre_clusters[ilive];
            map_cluster_index[live] = ilive;
            ilive2desc[ilive] = boost::add_vertex(ilive, g);
        }
        // Edge between every pair of in-scope clusters sharing a flash-time group.
        // Unmatched clusters get unique singleton groups, so they are never linked.
        for (size_t i = 0; i < pre_clusters.size(); ++i) {
            auto c1 = pre_clusters[i];
            if (!c1->get_scope_filter(scope)) continue;
            for (size_t j = i + 1; j < pre_clusters.size(); ++j) {
                auto c2 = pre_clusters[j];
                if (!c2->get_scope_filter(scope)) continue;
                if (flash_t0_group.at(c1) != flash_t0_group.at(c2)) continue;
                boost::add_edge(ilive2desc[map_cluster_index[c1]],
                                ilive2desc[map_cluster_index[c2]], g);
            }
        }

        // Diagnostic: how many in-scope matched clusters fall into how many
        // distinct flash-time groups (each multi-member group is merged into one).
        std::set<int> matched_groups;
        int n_matched_inscope = 0;
        for (auto c : pre_clusters) {
            if (c->get_scope_filter(scope) && c->get_scalar<int>("flash", -1) >= 0) {
                ++n_matched_inscope;
                matched_groups.insert(flash_t0_group.at(c));
            }
        }
        log->debug("flash-t0 merge: {} in-scope matched clusters -> {} flash-time groups",
                   n_matched_inscope, matched_groups.size());

        merge_clusters(g, live_grouping);
    }

    std::vector<Cluster *> live_clusters = live_grouping.children();

    for (size_t i=0;i!=live_clusters.size();i++){
        if (!live_clusters.at(i)->get_scope_filter(scope)) continue; // move on if the cluster is not in the scope filter ...
        if (live_clusters.at(i)->get_default_scope().hash() != scope.hash()) {
            live_clusters.at(i)->set_default_scope(scope);
            // std::cout << "Test: Set default scope: " << pc_name << " " << coords[0] << " " << coords[1] << " " << coords[2] << " " << cluster->get_default_scope().hash() << " " << scope.hash() << std::endl;
        }

        // if there is a cc component, record the main cluster as id of the blobs???
        auto old_cc_array = live_clusters.at(i)->get_pcarray("isolated", "perblob");
        log->trace("old_cc_array: {} clusters", std::set<int>(old_cc_array.begin(), old_cc_array.end()).size());

        // currently reset the cc component (todo: find the main component)

        // do the examine graph
        auto b2groupid = live_clusters.at(i)->connected_blobs(dv, pcts, graph_name);
        log->trace("b2groupid: {} clusters", std::set<int>(b2groupid.begin(), b2groupid.end()).size());
        
        bool flag_largest = false;
        // Compare old and new cluster groupings
        if (old_cc_array.size() == b2groupid.size()) {
            // Check if main cluster exists in old array
            bool has_main = std::find(old_cc_array.begin(), old_cc_array.end(), -1) != old_cc_array.end();
            if (has_main) {
                // Count overlap between old main cluster and each new cluster
                std::map<int, int> overlap_counts;
                for (size_t j = 0; j < old_cc_array.size(); j++) {
                    if (old_cc_array[j] == -1) {  // this blob was in main cluster
                        overlap_counts[b2groupid[j]]++;
                    }
                }
                
                // Find cluster with maximum overlap
                int max_overlap = 0;
                int main_cluster_id = -1;
                for (const auto& pair : overlap_counts) {
                    if (pair.second > max_overlap) {
                        max_overlap = pair.second;
                        main_cluster_id = pair.first;
                    }
                }
                
                // Mark new main cluster
                for (auto& id : b2groupid) {
                    if (id == main_cluster_id) {
                        id = -1;
                    }
                }
            } else{
                flag_largest = true;
            }
        }else{
            flag_largest = true;
        }

        //run the longest ...
        if (flag_largest){
            // No main cluster existed, find longest cluster.
            // Single sweep: compute length for each unique non-(-1) id and track the maximum.
            std::map<int, double> cluster_lengths;
            for (const auto& id : b2groupid) {
                if (id >= 0 && cluster_lengths.find(id) == cluster_lengths.end()) {
                    cluster_lengths[id] = get_length(live_clusters.at(i), b2groupid, id);
                }
            }

            // Find cluster with maximum length
            double max_length = 0;
            int longest_cluster_id = -1;
            for (const auto& pair : cluster_lengths) {
                if (pair.second > max_length) {
                    max_length = pair.second;
                    longest_cluster_id = pair.first;
                }
            }
           // std::cout << "Test: " << longest_cluster_id << " " << max_length << std::endl;
            
            // Mark longest cluster as main
            for (auto& id : b2groupid) {
                if (id == longest_cluster_id) {
                    id = -1;
                }
            }
        }

        live_clusters.at(i)->put_pcarray(b2groupid, "isolated", "perblob");

    }






}
