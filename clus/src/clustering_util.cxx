#include <WireCellClus/ClusteringFuncs.h>

#include <iostream>              // temp debug

using namespace WireCell::PointCloud::Facade;


void WireCell::PointCloud::Facade::merge_clusters(
    cluster_connectivity_graph_t& g,
    Grouping& grouping,
    cluster_set_t& known_clusters, // in/out
    const std::string& aname, const std::string& pcname)
{
    std::unordered_map<int, int> desc2id;
    std::unordered_map<int, std::set<int> > id2desc;
    /*int num_components =*/ boost::connected_components(g, boost::make_assoc_property_map(desc2id));
    for (const auto& [desc, id] : desc2id) {
        id2desc[id].insert(desc);
    }

    // Note, here we do an unusual thing and COPY the vector of children
    // facades.  In most simple access we would get the reference to the child
    // vector to save a little copy time.  We explicitly copy here as we must
    // preserve the original order of children facades even as we remove them
    // from the grouping.  As each child facade is removed, it's
    // unique_ptr<node> is returned which we ignore/drop and thus the child
    // facade dies along with its node.  This leaves the orig_clusters element
    // that was just holding the pointer to the doomed facade now holding
    // invalid memory.  But, it is okay as we never revisit the same cluster in
    // the grouping.  All that to explain a missing "&"! :)
    auto orig_clusters = grouping.children();

    const bool savecc = aname.size() > 0 && pcname.size() > 0;

    for (const auto& [id, descs] : id2desc) {
        if (descs.size() < 2) {
            continue;
        }

        // it starts with no cluster facade
        Cluster& fresh_cluster = grouping.make_child();

        std::vector<int> cc;
        int parent_id = 0;
        for (const auto& desc : descs) {
            const int idx = g[desc];
            if (idx < 0) {  // no need anymore ...
                continue;
            }

            auto live = orig_clusters[idx];
            fresh_cluster.take_children(*live, true);

            if (savecc) {
                cc.resize(fresh_cluster.nchildren(), parent_id);
                ++parent_id;
            }

            known_clusters.erase(live);
            grouping.destroy_child(live);
            assert(live == nullptr);
        }
        if (savecc) {
            fresh_cluster.put_pcarray(cc, aname, pcname);
        }
        known_clusters.insert(&fresh_cluster);
    }



    // fixme: sanity check / debugging.  remove this if you find it committed.
    for (const auto* cluster : grouping.children()) {
        if (!cluster) {
            std::cerr << "merge_clusters: null live cluster on output!\n";
            continue;
        }
        if (! cluster->nchildren()) {
            std::cerr << "merge_clusters: empty live cluster on output!\n";
            continue;
        }
        for (const auto* blob : cluster->children()) {
            if (!blob) {
                std::cerr << "merge_clusters: null live blob on output!\n";
                continue;
            }
        }
    }

}


