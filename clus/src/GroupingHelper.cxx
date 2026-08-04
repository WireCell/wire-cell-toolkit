#include "WireCellClus/GroupingHelper.h"

namespace {
    /// Order clusters by their stable ident() instead of their address.
    struct ClusterIdentLess {
        bool operator()(const WireCell::Clus::Facade::Cluster* a,
                        const WireCell::Clus::Facade::Cluster* b) const {
            return a->ident() < b->ident();
        }
    };
}

std::map<WireCell::Clus::Facade::Cluster*, std::tuple<WireCell::Clus::Facade::Cluster*, int, WireCell::Clus::Facade::Cluster*>>
WireCell::Clus::Facade::process_groupings_helper(
    WireCell::Clus::Facade::Grouping& original,
    WireCell::Clus::Facade::Grouping& shadow,
    const std::string& aname,
    const std::string& pname)  // Removed const here
{
    // current cluster,  corresponding shadow_cluster, its id, the main cluster of this cluster ...
    std::map<Cluster*, std::tuple<Cluster*, int, Cluster*>> result;
    
    // Step 1: Map original clusters to shadow clusters.
    //
    // Ordered by ident(), not by address.  The step-2 loop below walks this map
    // and calls Grouping::separate() on each entry, and separate() MINTS NEW
    // CLUSTERS -- so the walk order sets the idents the split products receive,
    // and an address-ordered walk made those idents depend on heap layout.
    // ident() is unique among a grouping's children, so this is a total order.
    //
    // NOTE: as of this writing this whole function is UNREACHABLE -- its only
    // call site is commented out (clustering_retile.cxx:162) and no other
    // translation unit references it.  The fix is therefore inert by
    // construction and no A/B gate applies to it; it is made here so the defect
    // does not come back to life with the call site (doc pr/28 sec 15).
    std::map<Cluster*, Cluster*, ClusterIdentLess> orig_to_shadow;
    for (auto* orig_cluster : original.children()) {
        for (auto* shad_cluster : shadow.children()) {
            if (orig_cluster->ident() == shad_cluster->ident()) {
                orig_to_shadow[orig_cluster] = shad_cluster;
                break;
            }
        }
    }
    
    // std::cout << "haha: " << orig_to_shadow.size() << " " << original.children().size() << " " << shadow.children().size() << std::endl;

    // Step 2: Process each pair
    for (const auto& [orig_cluster, shad_cluster] : orig_to_shadow) {
        // std::cout << orig_cluster << " " << shad_cluster << std::endl;
        // Get cluster index array
        auto cc = orig_cluster->get_pcarray(aname, pname);
        std::vector<int> cc_vec(cc.begin(), cc.end());
        // Create a non-const pointer for separate()
        Cluster* mutable_cluster = orig_cluster;
        // Separate clusters
        auto scope_transform = mutable_cluster->get_scope_transform(mutable_cluster->get_default_scope());
        auto& scope = mutable_cluster->get_default_scope();
        mutable_cluster->get_scope_filter(scope);
        auto orig_splits = original.separate(mutable_cluster, cc_vec);
       


        // Get cluster index array
        auto shad_cc = shad_cluster->get_pcarray(aname, pname);
        std::vector<int> shad_cc_vec(shad_cc.begin(), shad_cc.end());
        // Create a non-const pointer for separate()
        Cluster* mutable_shad_cluster = shad_cluster;
        // Separate clusters
        mutable_shad_cluster->get_scope_filter(scope);
        auto shad_splits = shadow.separate(mutable_shad_cluster, shad_cc_vec);
       

        // fill in the main cluster information ...
        result[mutable_cluster] = std::make_tuple(mutable_shad_cluster, -1, mutable_cluster);

        for (const auto& [id1, cluster1] : orig_splits) {
            for (const auto& [id2, cluster2] : shad_splits){
                if (id2==id1){
                    result[cluster1] = std::make_tuple(cluster2, id1, mutable_cluster);
                    break;
                }
           }
        }
    }
    
    return result;
}
