#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"


class ClusteringExamineBundles;
WIRECELL_FACTORY(ClusteringExamineBundles, ClusteringExamineBundles,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)


using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;


static void clustering_examine_bundles(
        Grouping& live_grouping,
        IDetectorVolumes::pointer dv, 
        IPCTransformSet::pointer pcts,
        const Tree::Scope& scope,
        const bool use_ctpc);

class ClusteringExamineBundles : public IConfigurable, public Clus::IEnsembleVisitor, private NeedDV, private NeedPCTS, private NeedScope {
public:
    ClusteringExamineBundles() {}
    virtual ~ClusteringExamineBundles() {}
    
    void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config);
        NeedScope::configure(config);

        // If false, then DV and PCTS are not needed.
        use_ctpc_ = get<bool>(config, "use_ctpc", use_ctpc_);
    }

    void visit(Ensemble& ensemble) const {
        auto& live = *ensemble.with_name("live").at(0);
        clustering_examine_bundles(live, m_dv, m_pcts, m_scope, use_ctpc_);
    }
        
private:
    bool use_ctpc_{true};
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
    const bool use_ctpc)
{
    // std::cout << "Test Examine Bundles" << std::endl;

    std::vector<Cluster *> live_clusters = live_grouping.children();

    for (size_t i=0;i!=live_clusters.size();i++){
        if (!live_clusters.at(i)->get_scope_filter(scope)) continue; // move on if the cluster is not in the scope filter ...
        if (live_clusters.at(i)->get_default_scope().hash() != scope.hash()) {
            live_clusters.at(i)->set_default_scope(scope);
            // std::cout << "Test: Set default scope: " << pc_name << " " << coords[0] << " " << coords[1] << " " << coords[2] << " " << cluster->get_default_scope().hash() << " " << scope.hash() << std::endl;
        }

        // if there is a cc component, record the main cluster as id of the blobs???
        auto old_cc_array = live_clusters.at(i)->get_pcarray("isolated", "perblob");
        
        // currently reset the cc component (todo: find the main component)

        // do the examine graph
        auto b2groupid = live_clusters.at(i)->connected_blobs(dv, pcts);
        
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
            // No main cluster existed, find longest cluster
            std::map<int, double> cluster_lengths;
            std::set<int> unique_ids;
            for (const auto& id : b2groupid) {
                if (id >= 0) {  // skip any existing -1 values
                    unique_ids.insert(id);
                }
            }
            
            // Calculate length for each unique cluster ID
            for (const auto& id : unique_ids) {
                double length = get_length(live_clusters.at(i) , b2groupid, id);
                cluster_lengths[id] = length;
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
