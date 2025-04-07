#include <WireCellClus/ClusteringFuncs.h>
#include "WireCellUtil/ExecMon.h"
#include <set>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;

void WireCell::PointCloud::Facade::clustering_switch_scope(
    Grouping& live_grouping,
    const IDetectorVolumes::pointer dv,          // detector volumes
    const std::string& pc_name,                  // point cloud name
    const std::vector<std::string>& coords,      // coordinate names
    const std::string& correction_name          // name of correction to apply
)
{
    using spdlog::debug;
    using spdlog::info;

    // Check that live_grouping has at least one wpid
    if (live_grouping.wpids().empty()) {
        throw std::runtime_error("Live grouping must have at least one wpid");
    }

    // Get all clusters from the grouping
    std::vector<Cluster*> live_clusters = live_grouping.children(); // copy
    
    // Set the default scope for all clusters in the live grouping
    Tree::Scope default_scope{pc_name, coords};
    
    // Process each cluster
    for (size_t iclus = 0; iclus < live_clusters.size(); ++iclus) {
        Cluster* cluster = live_clusters.at(iclus);
                
        if (correction_name == "T0Correction") {
                
            // Get original bounds before correction
            info("Cluster {} original bounds:", iclus);
            const auto [earliest_orig, latest_orig] = cluster->get_earliest_latest_points();
            info("  earliest: {}", earliest_orig);
            info("  latest: {}", latest_orig);
            
            // Add corrected points - this returns filter values for each blob
            std::vector<int> filter_results = cluster->add_corrected_points(dv, correction_name);
            
            // Get the new scope with corrected points
            const auto correction_scope = cluster->get_scope(correction_name);
            
            // Set this as the default scope for viewing
            cluster->set_default_scope(correction_scope);
            
            // Get bounds after correction
            info("Cluster {} corrected bounds:", iclus);
            const auto [earliest_corr, latest_corr] = cluster->get_earliest_latest_points();
            info("  earliest: {}", earliest_corr);
            info("  latest: {}", latest_corr);
            
            // Get unique filter result values
            std::set<int> filter_result_set(filter_results.begin(), filter_results.end());
            info("Cluster {} has {} unique filter results:", iclus, filter_result_set.size());
            for (const auto& result : filter_result_set) {
                info("  filter result: {}", result);
            }
            
            // Separate the cluster based on filter results
            // This will create new clusters in the grouping
            auto separated_clusters = live_grouping.separate(cluster, filter_results, true);
            
            // Process each separated cluster
            for (auto& [id, new_cluster] : separated_clusters) {
                info("  Separated cluster filter={}, nchildren={}", id, new_cluster->nchildren());
                
                // Set the new scope as default for the separated cluster
                new_cluster->set_default_scope(correction_scope);
                
                // Get bounds of the separated cluster
                const auto [earliest_sep, latest_sep] = new_cluster->get_earliest_latest_points();
                info("    earliest: {}", earliest_sep);
                info("    latest: {}", latest_sep);
            }
        }
    }
    
    info("Completed scope switching with correction: {}", correction_name);
}

#pragma GCC diagnostic pop