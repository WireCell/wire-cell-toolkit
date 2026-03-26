#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/Logger.h"


class ClusteringRecoveringBundle;
WIRECELL_FACTORY(ClusteringRecoveringBundle, ClusteringRecoveringBundle,
                 WireCell::INamed, WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

/**
 * Clustering function that processes beam-flash flagged clusters and separates them
 * into individual bundles based on isolated blob components.
 * This function recovers separated clusters from over-clustered beam-flash events.
 */
class ClusteringRecoveringBundle : public IConfigurable, public Clus::IEnsembleVisitor, public Aux::Logger, private NeedDV, private NeedPCTS, private NeedScope {
public:
    ClusteringRecoveringBundle()
        : Aux::Logger("ClusteringRecoveringBundle", "clus")
    {}
    virtual ~ClusteringRecoveringBundle() {}

    virtual void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config);
        NeedScope::configure(config);
        m_grouping_name = get<std::string>(config, "grouping", "live");
        m_array_name = get<std::string>(config, "array_name", "isolated");
        m_pcarray_name = get<std::string>(config, "pcarray_name", "perblob");
        m_graph_name = get<std::string>(config, "graph_name", m_graph_name);
    }
    
    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["detector_volumes"] = "";
        cfg["pc_transforms"] = "";
        cfg["pc_name"] = "3d";
        cfg["coords"] = Json::arrayValue;
        cfg["grouping"] = m_grouping_name;
        cfg["array_name"] = m_array_name;
        cfg["pcarray_name"] = m_pcarray_name;
        cfg["graph_name"] = m_graph_name;
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const {
        // Get the specified grouping (default: "live")
        auto groupings = ensemble.with_name(m_grouping_name);
        if (groupings.empty()) {
            log->debug("No '{}' grouping found", m_grouping_name);
            return;
        }

        auto& grouping = *groupings.at(0);

        // Container to hold clusters after the initial filter
        std::vector<Cluster*> filtered_clusters;

        for (auto* cluster : grouping.children()) {
            if (cluster->get_flag(Flags::beam_flash)){
                filtered_clusters.push_back(cluster);
            }
        }

        log->debug("Found {} beam-flash flagged clusters", filtered_clusters.size());

        // Process each filtered cluster
        for (auto* cluster : filtered_clusters) {
            process_cluster(grouping, cluster, m_dv, m_pcts);
        }

        size_t nassoc = 0;
        for (auto* cluster : grouping.children()) {
            if (cluster->get_flag(Flags::associated_cluster)) {
                ++nassoc;
            }
        }
        log->debug("Associated clusters: {}", nassoc);
    }

private:
    std::string m_grouping_name{"live"};
    std::string m_array_name{"isolated"};
    std::string m_pcarray_name{"perblob"};
    std::string m_graph_name{"relaxed"};

    void process_cluster(Grouping& grouping, Cluster* cluster,
                         IDetectorVolumes::pointer dv,
                         IPCTransformSet::pointer pcts) const {

        std::vector<int> cc_vec;

        // Examine step: recompute connected components via graph algorithm
        if (cluster->get_scope_filter(m_scope)) {
            if (cluster->get_default_scope().hash() != m_scope.hash()) {
                cluster->set_default_scope(m_scope);
            }       

            // Determine which component is the "main" cluster
            auto old_cc_array = cluster->get_pcarray(m_array_name, m_pcarray_name);
            log->debug("old_cc_array: {} clusters", std::set<int>(old_cc_array.begin(), old_cc_array.end()).size());

            auto b2groupid = cluster->connected_blobs(dv, pcts, m_graph_name);
            log->debug("examine in recovering: b2groupid has {} components",
                       std::set<int>(b2groupid.begin(), b2groupid.end()).size());
            bool flag_largest = false;
            if (old_cc_array.size() == b2groupid.size()) {
                bool has_main = std::find(old_cc_array.begin(), old_cc_array.end(), -1) != old_cc_array.end();
                if (has_main) {
                    // Find the new component that overlaps most with the old main
                    std::map<int, int> overlap_counts;
                    for (size_t j = 0; j < old_cc_array.size(); j++) {
                        if (old_cc_array[j] == -1) overlap_counts[b2groupid[j]]++;
                    }
                    int max_overlap = 0, main_id = -1;
                    for (const auto& kv : overlap_counts) {
                        if (kv.second > max_overlap) { max_overlap = kv.second; main_id = kv.first; }
                    }
                    for (auto& id : b2groupid) { if (id == main_id) id = -1; }
                } else {
                    flag_largest = true;
                }
            } else {
                flag_largest = true;
            }

            if (flag_largest) {
                // Use the longest component as the main cluster
                std::set<int> unique_ids;
                for (const auto& id : b2groupid) { if (id >= 0) unique_ids.insert(id); }
                std::map<int, double> cluster_lengths;
                for (const auto& id : unique_ids) {
                    cluster_lengths[id] = get_length(cluster, b2groupid, id);
                }
                double max_length = 0; int longest_id = -1;
                for (const auto& kv : cluster_lengths) {
                    if (kv.second > max_length) { max_length = kv.second; longest_id = kv.first; }
                }
                for (auto& id : b2groupid) { if (id == longest_id) id = -1; }
            }

            cluster->put_pcarray(b2groupid, m_array_name, m_pcarray_name);
            cc_vec = b2groupid;
        } else {
            // Fall back to stored pcarray if scope is not available
            auto cc = cluster->get_pcarray(m_array_name, m_pcarray_name);
            cc_vec.assign(cc.begin(), cc.end());
        }

        // Skip if there are fewer than 2 blobs
        if (cc_vec.size() < 2) return;
        // Skip if all blobs belong to a single component
        if (std::set<int>(cc_vec.begin(), cc_vec.end()).size() < 2) return;

        // Perform the separation
        const int main_cluster_id = cluster->ident();
        auto splits = grouping.separate(cluster, cc_vec);
        cluster->set_flag(Flags::main_cluster);

        int sub_id = 1;
        for (auto& [id, new_cluster] : splits) {
            new_cluster->set_ident(main_cluster_id * 100 + sub_id);
            ++sub_id;
            new_cluster->set_flag(Flags::associated_cluster);
        }
    }
};