#include "WireCellClus/GroupingHelper.h"

std::map<WireCell::Clus::Facade::Cluster*, std::tuple<WireCell::Clus::Facade::Cluster*, int, WireCell::Clus::Facade::Cluster*>> 
WireCell::Clus::Facade::process_groupings_helper(
    WireCell::Clus::Facade::Grouping& original,
    WireCell::Clus::Facade::Grouping& shadow,
    const std::string& aname,
    const std::string& pname)  // Removed const here
{
    // current cluster,  corresponding shadow_cluster, its id, the main cluster of this cluster ...
    std::map<Cluster*, std::tuple<Cluster*, int, Cluster*>> result;
    
    // Step 1: Map original clusters to shadow clusters
    std::map<Cluster*, Cluster*> orig_to_shadow;
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
