#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"


class ClusteringRecoveringBundle;
WIRECELL_FACTORY(ClusteringRecoveringBundle, ClusteringRecoveringBundle,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

/**
 * Clustering function that processes beam-flash flagged clusters and separates them
 * into individual bundles based on isolated blob components.
 * This function recovers separated clusters from over-clustered beam-flash events.
 */
class ClusteringRecoveringBundle : public IConfigurable, public Clus::IEnsembleVisitor {
public:
    ClusteringRecoveringBundle() {}
    virtual ~ClusteringRecoveringBundle() {}

    virtual void configure(const WireCell::Configuration& config) {
        m_grouping_name = get<std::string>(config, "grouping", "live");
        m_array_name = get<std::string>(config, "array_name", "isolated");
        m_pcarray_name = get<std::string>(config, "pcarray_name", "perblob");
    }
    
    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping_name;
        cfg["array_name"] = m_array_name;
        cfg["pcarray_name"] = m_pcarray_name;
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const {
        using spdlog::debug;
        
        // Get the specified grouping (default: "live")
        auto groupings = ensemble.with_name(m_grouping_name);
        if (groupings.empty()) {
            debug("ClusteringRecoveringBundle: No '{}' grouping found", m_grouping_name);
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

        debug("ClusteringRecoveringBundle: Found {} beam-flash flagged clusters", 
              filtered_clusters.size());

        // Process each filtered cluster
        for (auto* cluster : filtered_clusters) {
            process_cluster(grouping, cluster);
        }
        
        debug("ClusteringRecoveringBundle: Processing complete");
    }

private:
    std::string m_grouping_name{"live"};
    std::string m_array_name{"isolated"};
    std::string m_pcarray_name{"perblob"};

    void process_cluster(Grouping& grouping, Cluster* cluster) const {


        // std::cout << grouping.children().size() << " clusters before separation." << std::endl;

        // Separate the clusters into separated pieces
        auto cc = cluster->get_pcarray(m_array_name, m_pcarray_name);
        // Convert span to vector
        std::vector<int> cc_vec(cc.begin(), cc.end());
        
        // Skip if there are fewer than 2 components to separate
        if (cc_vec.size() < 2) return;
        
        // Perform the separation
        auto splits = grouping.separate(cluster, cc_vec);
        cluster->set_flag(Flags::main_cluster);

        // Apply the scope filter settings to all new clusters
        for (auto& [id, new_cluster] : splits) {
            // Store the split/group ID as a scalar
            // new_cluster->set_scalar<int>("split_id", id);
            // // Optionally also store the original parent's ident
            // new_cluster->set_scalar<int>("parent_ident", cluster->ident());
            new_cluster->set_flag(Flags::associated_cluster);
        }

        // std::cout << grouping.children().size() << " clusters after separation." << std::endl;
    }
};