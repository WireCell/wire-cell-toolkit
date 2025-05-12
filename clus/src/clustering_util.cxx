#include <WireCellClus/ClusteringFuncs.h>

#include <iostream>              // temp debug

using namespace WireCell::Clus::Facade;


// Add this to your clustering_util.cxx file

std::tuple<geo_point_t, double, double, double> 
WireCell::Clus::Facade::extract_geometry_params(
    const Grouping& grouping,
    const IDetectorVolumes::pointer dv)
{
    geo_point_t drift_dir(1, 0, 0);  // initialize drift direction
    double angle_u = 0, angle_v = 0, angle_w = 0;  // initialize angles

    // Find the first valid WirePlaneId in the grouping
    for (const auto& gwpid : grouping.wpids()) {
        // Update drift direction based on face orientation
        int face_dirx = dv->face_dirx(gwpid);
        drift_dir.x(face_dirx);
        
        // Create wpids for all three planes with the same APA and face
        WirePlaneId wpid_u(kUlayer, gwpid.face(), gwpid.apa());
        WirePlaneId wpid_v(kVlayer, gwpid.face(), gwpid.apa());
        WirePlaneId wpid_w(kWlayer, gwpid.face(), gwpid.apa());
        
        // Get wire directions for all planes
        Vector wire_dir_u = dv->wire_direction(wpid_u);
        Vector wire_dir_v = dv->wire_direction(wpid_v);
        Vector wire_dir_w = dv->wire_direction(wpid_w);
        
        // Calculate angles
        angle_u = std::atan2(wire_dir_u.z(), wire_dir_u.y());
        angle_v = std::atan2(wire_dir_v.z(), wire_dir_v.y());
        angle_w = std::atan2(wire_dir_w.z(), wire_dir_w.y());
        
        // Only need to process the first valid WirePlaneId
        break;
    }
    
    return std::make_tuple(drift_dir, angle_u, angle_v, angle_w);
}

std::vector<Cluster*> WireCell::Clus::Facade::merge_clusters(
    cluster_connectivity_graph_t& g,
    Grouping& grouping,
    const std::string& aname, const std::string& pcname)
{
    std::unordered_map<int, int> desc2id;
    std::unordered_map<int, std::set<int> > id2desc;
    /*int num_components =*/ boost::connected_components(g, boost::make_assoc_property_map(desc2id));
    for (const auto& [desc, id] : desc2id) {
        id2desc[id].insert(desc);
    }

    std::vector<Cluster*> fresh;

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
            const Tree::Scope& default_scope = live->get_default_scope();
            bool flag = live->get_scope_filter(default_scope);
            auto scope_transform = live->get_scope_transform(live->get_default_scope());
            fresh_cluster.take_children(*live, true);
            // set scope filter assuming union of all clusters
            {   
                fresh_cluster.set_default_scope(default_scope);
                if (flag) fresh_cluster.set_scope_filter(default_scope, flag);
                fresh_cluster.set_scope_transform(default_scope, scope_transform);
            }

            if (savecc) {
                cc.resize(fresh_cluster.nchildren(), parent_id);
                ++parent_id;
            }

            grouping.destroy_child(live);
            assert(live == nullptr);
        }
        if (savecc) {
            fresh_cluster.put_pcarray(cc, aname, pcname);
        }

        // Normally, is is weird/wrong to store a pointer to a reference.  But,
        // we know the Cluster facade is held by the pc tree node that we just
        // added to the grouping node.  
        fresh.push_back(&fresh_cluster);
    }

    return fresh;
}


