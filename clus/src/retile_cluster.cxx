// This provides RetileCluster aka "IPCTreeMutate".  
//
// Warning: this lives up to its name.  It may change the input cluster.
// 
// It requires the input cluster node to have a (grouping) node parent.
//
// The retiling of a cluster follows this general sequence:
//
// 1) constructs layers of "activity" from input grouping.
// 2) applies "hacks" to the activity.
// 3) runs WCT tiling to create blobs.
// 4) runs blobs sampling to make point clouds
// 5) produces clusters such that the new blobs formed from an old cluster form a new "shadow" cluster.
// 6) forms a PC-tree
// 7) outputs the new grouping

#include "retile_cluster.h"  // Include the header instead of defining the class here

WIRECELL_FACTORY(RetileCluster, WireCell::Clus::RetileCluster,
                 WireCell::IConfigurable, WireCell::IPCTreeMutate)




// Segregate this weird choice for namespace.
namespace WCF = WireCell::Clus::Facade;

// Nick name for less typing.
namespace WRG = WireCell::RayGrid;

// static void debug_cluster(WCF::Cluster* cluster, const std::string& ctx)
// {
//     std::cout << "Cluster: " << cluster->node()->value.as_string(false) << "\n";
//     auto nblobs = cluster->kd_blobs().size();
//     auto nchild = cluster->nchildren();
//     if (nblobs == nchild) return;
//     std::cout << ctx << " n kd_blobs=" << nblobs << " nchildren=" << nchild << "\n";
// }




// Now can handle all APA/Faces 
void RetileCluster::configure(const WireCell::Configuration& cfg)
{
    NeedDV::configure(cfg);
    NeedPCTS::configure(cfg);

    if (cfg.isMember("samplers") && cfg["samplers"].isArray()) {
        // Process array of samplers
        for (const auto& sampler_cfg : cfg["samplers"]) {
            int apa = sampler_cfg["apa"].asInt();
            int face = sampler_cfg["face"].asInt();
            std::string sampler_name = sampler_cfg["name"].asString();
            
            if (sampler_name.empty()) {
                raise<ValueError>("RetileCluster requires an IBlobSampler name for APA %d face %d", apa, face);
            }
            // std::cout << "Test: " << apa << " " << face << " " << sampler_name << std::endl;
            auto sampler_ptr = Factory::find_tn<IBlobSampler>(sampler_name);
            m_samplers[apa][face] = sampler_ptr;
        }
    }

    std::vector<IAnodePlane::pointer> anodes_tn;
    for (const auto& aname : cfg["anodes"]) {
        auto anode = Factory::find_tn<IAnodePlane>(aname.asString());
        anodes_tn.push_back(anode);
        for (const auto& face1 : anode->faces()) {
            int apa = anode->ident();
            int face = face1->which();
            m_face[apa][face] = face1;
            const auto& coords = face1->raygrid();
            if (coords.nlayers() != 5) {
                raise<ValueError>("unexpected number of ray grid layers: %d", coords.nlayers());
            }
            // std::cout <<"Test: " << apa << " " << face << " " << coords.nlayers() << std::endl;
            // Get wire info for each plane
            m_plane_infos[apa][face].clear();
            m_plane_infos[apa][face].push_back(Aux::get_wire_plane_info(face1, kUlayer));
            m_plane_infos[apa][face].push_back(Aux::get_wire_plane_info(face1, kVlayer));
            m_plane_infos[apa][face].push_back(Aux::get_wire_plane_info(face1, kWlayer));
        }
    }
    
    // Add time cut configuration
    m_cut_time_low = get(cfg, "cut_time_low", -1e9);
    m_cut_time_high = get(cfg, "cut_time_high", 1e9);

}


// Step 0.  The RetileCluster only directly gets a cluster node but needs
// context from the parent grouping which likely can be cached between calls.
Facade::Cluster* RetileCluster::reinitialize(Points::node_type& node) const
{
    auto* cluster = node.value.facade<Cluster>();
    if (!cluster || !cluster->grouping()) {
        return nullptr;
    }
    if (m_grouping && m_grouping == cluster->grouping()) {
        return cluster;
    }
    m_grouping = cluster->grouping();

    m_wpid_angles.clear();
    for (const auto& gwpid : m_grouping->wpids()) {
        // gwpids are "all" type - no specific layer so we must remake per-layer wpids
        int apa = gwpid.apa();
        int face = gwpid.face();
        std::vector<double> angles(3);
        for (size_t ind=0; ind<3; ++ind) {
            // iplane2layer is in WirePlaneId.h
            WirePlaneId wpid(iplane2layer[ind], face, apa);
            Vector wire_dir = m_dv->wire_direction(wpid);
            angles[ind] = std::atan2(wire_dir.z(), wire_dir.y());
        }
        m_wpid_angles[gwpid] = angles;
    }
    return cluster;
}


// Step 1. Build activities from blobs in a cluster.
void RetileCluster::get_activity(const Cluster& cluster, std::map<std::pair<int, int>, std::vector<WRG::measure_t> >& map_slices_measures, int apa, int face) const
{
    const int nlayers = 2+3;

    // checkme: this assumes "iend" is the usual one-past-last aka [ibeg,iend)
    // forms a half-open range.  I'm not sure if PointTreeBuilding is following
    // this or not.

    

    // for (auto& info : plane_infos) {
    //     std::cout << "test1: " << info.start_index << " " << info.end_index << " " << info.total_wires << std::endl;
    // }

    int (WCF::Blob::*wmin[])(void) const = {
        &WCF::Blob::u_wire_index_min,
        &WCF::Blob::v_wire_index_min,
        &WCF::Blob::w_wire_index_min
    };

    int (WCF::Blob::*wmax[])(void) const = {
        &WCF::Blob::u_wire_index_max,
        &WCF::Blob::v_wire_index_max,
        &WCF::Blob::w_wire_index_max
    };
        
    const double hit=1.0;       // actual charge value does not matter to tiling.

    for (const auto* fblob : cluster.children()) {
        int tslice_beg = fblob->slice_index_min();
        int tslice_end = fblob->slice_index_max();

        // if blob is not consistent skip ...
        auto blob_wpid = fblob->wpid();
        if (blob_wpid.apa()!=apa || blob_wpid.face()!=face) continue;

        auto& measures = map_slices_measures[std::make_pair(tslice_beg, tslice_end)];
        
        // if (tslice_beg == tslice_end) {
        //     std::cout << "Test: Same: " << tslice_beg << " " << tslice_end << std::endl;
        // }

        if (measures.size()==0){
            measures.resize(nlayers);
            // what to do the first two views???
            measures[0].push_back(1);
            measures[1].push_back(1);
            measures[2].resize(m_plane_infos.at(apa).at(face)[0].total_wires, 0);
            measures[3].resize(m_plane_infos.at(apa).at(face)[1].total_wires, 0);
            measures[4].resize(m_plane_infos.at(apa).at(face)[2].total_wires, 0);
            // std::cout << measures[2].size() << " " << measures[3].size() << " " << measures[4].size() << std::endl;
        }

        // the three views ...
        for (int index=0; index<3; ++index) {
            const int layer = index + 2;
            WRG::measure_t& m = measures[layer];
            // Make each "wire" in each blob's bounds of this plane "hit".
            int ibeg = (fblob->*wmin[index])();
            int iend = (fblob->*wmax[index])();    
            while (ibeg < iend) {
                m[ibeg++] = hit;
            }
            //std::cout << ibeg << " " << iend << " " << index << " " << hit << std::endl;
        }
    }

    // std::cout << "Test: Org: " << map_slices_measures.size() << " " << cluster.children().size() << std::endl;

}


// Step 2. Modify activity to suit.
void RetileCluster::hack_activity(
    const Cluster& cluster,
    std::map<std::pair<int, int>, std::vector<WRG::measure_t> >& map_slices_measures,
    const std::vector<size_t>& path_wcps,
    int apa, int face) const
{

    // for (auto it = map_slices_measures.begin(); it!= map_slices_measures.end(); it++){
    //     std::cout << "Before: " << it->first.first << " " << it->first.second << " " << it->second.size() << std::endl;
    //     for (int i=0; i!=5; i++){
    //         std::cout << it->second[i].size() << " ";
    //     }
    //     std::cout << std::endl;
    // }


    const double low_dis_limit = 0.3 * units::cm;
    // Get path points
    // auto path_wcps = cluster.get_path_wcps();
    std::vector<std::pair<geo_point_t, WirePlaneId>> path_pts;

    // Convert list points to vector with interpolation
    for (const auto& wcp : path_wcps) {
        geo_point_t p= cluster.point3d_raw(wcp); // index ... // raw data points ...
        auto wpid_p = cluster.wire_plane_id(wcp); // wpid ...
        // std::cerr << "retile: path:" << wcp << " p:" << p << " wpid:" << wpid_p << "\n";
        if (path_pts.empty()) {
            path_pts.push_back(std::make_pair(p, wpid_p));
        } else {
            double dis = (p - path_pts.back().first).magnitude();
            if (dis < low_dis_limit) {
                path_pts.push_back(std::make_pair(p, wpid_p));
            } else {
                int ncount = int(dis/low_dis_limit) + 1;
                auto p2 = path_pts.back().first;
                auto wpid2 = path_pts.back().second;
                for (int i=0; i < ncount; i++) {
                    Point p1 = p2 + (p - p2) * (i+1)/ncount;
                    auto wpid_p1 = get_wireplaneid(p1, wpid_p, wpid2, m_dv);
                    path_pts.push_back(std::make_pair(p1, wpid_p1));
                }
            }
        }
    }


    std::vector<std::pair<int,int>> wire_limits;
    for (int i=0; i!=3; i++){
        wire_limits.push_back(std::make_pair(m_plane_infos.at(apa).at(face)[i].start_index, m_plane_infos.at(apa).at(face)[i].end_index));
        // std::cout << "Test: " << apa << " " << face << " " << wire_limits[i].first << " " << wire_limits[i].second << std::endl;
    }

    // this is to get the end of the time tick range = start_tick + tick_span
    const int tick_span = map_slices_measures.begin()->first.second -  map_slices_measures.begin()->first.first;

    // std::cout << "Test:  " << apa << " " << face << " " << tick_span << std::endl;

    // Flag points that have sufficient activity around them
    std::vector<bool> path_pts_flag(path_pts.size(), false);
    for (size_t i = 0; i < path_pts.size(); i++) {
        if (path_pts[i].second.apa() != apa || path_pts[i].second.face() != face) continue;
        auto [time_tick_u, u_wire] = cluster.grouping()->convert_3Dpoint_time_ch(path_pts[i].first, apa, m_face.at(apa).at(face)->which(), 0);
        auto [time_tick_v, v_wire] = cluster.grouping()->convert_3Dpoint_time_ch(path_pts[i].first, apa, m_face.at(apa).at(face)->which(), 1);
        auto [time_tick_w, w_wire] = cluster.grouping()->convert_3Dpoint_time_ch(path_pts[i].first, apa, m_face.at(apa).at(face)->which(), 2);
        //std::cout << time_tick_u <<  " " << u_wire << " " << v_wire << " " << w_wire << std::endl;

        int aligned_tick = std::round(time_tick_u *1.0/ tick_span) * tick_span;
        std::pair<int, int> tick_range = std::make_pair(aligned_tick, aligned_tick + tick_span);

        // Check for activity in neighboring wires/time
        // For each plane (U,V,W), count activity in current and adjacent wires
        std::vector<int> wire_hits = {0,0,0}; // counts for U,V,W planes
        std::vector<int> wires = {u_wire, v_wire, w_wire};
        
        for (size_t plane = 0; plane < 3; plane++) {
            // Check activity in current and adjacent wires
            for (int delta : {-1, 0, 1}) {
                int wire = wires[plane] + delta;
                if (wire < wire_limits[plane].first || wire > wire_limits[plane].second) 
                    continue;
                    
                int layer = plane + 2;
                if (map_slices_measures.find(tick_range) != map_slices_measures.end()) {
                    if (map_slices_measures[tick_range][layer][wire] > 0) {
                        wire_hits[plane] += (delta == 0) ? 1 : (delta == -1) ? 2 : 1;
                    }
                }
            }
        }
        
        // Set flag if sufficient activity found
        if (wire_hits[0] > 0 && wire_hits[1] > 0 && wire_hits[2] > 0 && 
            (wire_hits[0] + wire_hits[1] + wire_hits[2] >= 6)) {
            path_pts_flag[i] = true;
        }
        // std::cout << path_pts[i] << " " << wire_hits[0] << " " << wire_hits[1] << " " << wire_hits[2] << " " << path_pts_flag[i] << " " << aligned_tick/tick_span << " " << u_wire << " " << v_wire << " " << w_wire << " " << time_tick_u << " " << std::round(time_tick_u / tick_span) << std::endl;
        // std::cout << wire_hits[0] << " " << wire_hits[1] << " " << wire_hits[2] << " " << path_pts_flag[i] << std::endl;    
    }

    // Add missing activity based on path points
    for (size_t i = 0; i < path_pts.size(); i++) {
        if (path_pts[i].second.apa() != apa || path_pts[i].second.face() != face) continue;

        // Skip if point is well-covered by existing activity
        if (i == 0) {
            if (path_pts_flag[i] && path_pts_flag[i+1]) continue;
        } else if (i+1 == path_pts.size()) {
            if (path_pts_flag[i] && path_pts_flag[i-1]) continue;
        } else {
            if (path_pts_flag[i-1] && path_pts_flag[i] && path_pts_flag[i+1]) continue;
        }

        auto [time_tick_u, u_wire] = cluster.grouping()->convert_3Dpoint_time_ch(path_pts[i].first, apa, m_face.at(apa).at(face)->which(), 0);
        auto [time_tick_v, v_wire] = cluster.grouping()->convert_3Dpoint_time_ch(path_pts[i].first, apa, m_face.at(apa).at(face)->which(), 1);
        auto [time_tick_w, w_wire] = cluster.grouping()->convert_3Dpoint_time_ch(path_pts[i].first, apa, m_face.at(apa).at(face)->which(), 2);

        int aligned_tick = std::round(time_tick_u *1.0/ tick_span) * tick_span;

        // Add activity around this point
        for (int dt = -3; dt <= 3; dt++) {
            int time_slice = aligned_tick + dt * tick_span;
            if (time_slice < 0) continue;

            // Find or create time slice in measures map
            auto slice_key = std::make_pair(time_slice, time_slice+tick_span);  
            if (map_slices_measures.find(slice_key) == map_slices_measures.end()) {
                auto& measures = map_slices_measures[slice_key];
                measures = std::vector<WRG::measure_t>(5);  // 2+3 layers
                measures[0].push_back(1);  // First layer measurement 
                measures[1].push_back(1);  // Second layer measurement
                measures[2].resize(m_plane_infos.at(apa).at(face)[0].total_wires, 0);
                measures[3].resize(m_plane_infos.at(apa).at(face)[1].total_wires, 0); 
                measures[4].resize(m_plane_infos.at(apa).at(face)[2].total_wires, 0);
            }

            // Add activity for each plane
            std::vector<int> wires = {u_wire, v_wire, w_wire};
            for (size_t plane = 0; plane < 3; plane++) {
                auto& measures = map_slices_measures[slice_key][plane+2]; // +2 to skip first two layers
                
                for (int dw = -3; dw <= 3; dw++) {
                    int wire = wires[plane] + dw;
                    if (wire < wire_limits[plane].first || wire > wire_limits[plane].second ||
                        std::abs(dw) + std::abs(dt) > 3) 
                        continue;
    
                    measures.at(wire) = 1.0;  // Set activity
                }
            }

        }
    }


   // Loop through the map and remove slices with no activity in any plane view
    auto it = map_slices_measures.begin();
    while (it != map_slices_measures.end()) {
        bool missing_activity = false;
        
        // For each wire plane (U, V, W)
        for (int pind = 0; pind < 3; pind++) {
            const auto& measures = it->second[pind + 2]; // +2 to skip first two layers
            
            // Check if this plane has NO activity
            if (std::none_of(measures.begin(), measures.end(), [](double val) { return val > 0.0; })) {
                missing_activity = true;
                break;
            }
        }
        
        // If any plane has no activity, remove this slice
        if (missing_activity) {
            it = map_slices_measures.erase(it);
        } else {
            ++it;
        }
    }
   

}



// Step 3. Form IBlobs from activities.
std::vector<IBlob::pointer> RetileCluster::make_iblobs(std::map<std::pair<int, int>, std::vector<WRG::measure_t> >& map_slices_measures, int apa, int face) const
{
    std::vector<IBlob::pointer> ret;

    const auto& coords = m_face.at(apa).at(face)->raygrid();
    int blob_ident=0;
    int slice_ident = 0;
    for (auto it = map_slices_measures.begin(); it != map_slices_measures.end(); it++){
        // Do the actual tiling.
        WRG::activities_t activities = RayGrid::make_activities(m_face.at(apa).at(face)->raygrid(), it->second);
        auto bshapes = WRG::make_blobs(coords, activities);

    
        // {
        //     std::cerr << "abc: "
        //               << " s:"<<slice_ident<<" b:"<<blob_ident << " "
        //               << bshapes.size() << " " << activities.size() << " " << std::endl;
        //     for (const auto& activity : activities) {
        //         std::cerr <<"act: "
        //                   << " s:"<<slice_ident<<" b:"<<blob_ident << " "
        //                   << activity.as_string() << std::endl;
        //     }
        // }

        // Convert RayGrid blob shapes into IBlobs 
        const float blob_value = 0.0;  // tiling doesn't consider particular charge
        const float blob_error = 0.0;  // tiling doesn't consider particular charge
    
        for (const auto& bshape : bshapes) {
            IFrame::pointer sframe = nullptr;

            // {
            //     std::cerr << "blob: "
            //               << " s:"<<slice_ident<<" b:"<<blob_ident << " "
            //               << bshape << std::endl;
            //     for (const auto& strip : bshape.strips()) {
            //         std::cerr << "strip: "
            //                   << " s:"<<slice_ident<<" b:"<<blob_ident << " "
            //                   << strip << std::endl;
            //     }
            // }

            // 500 ns should be passed from outside?
            ISlice::pointer slice = std::make_shared<Aux::SimpleSlice>(sframe, slice_ident++, it->first.first*500*units::ns, (it->first.second - it->first.first)*500*units::ns);

            IBlob::pointer iblob = std::make_shared<Aux::SimpleBlob>(blob_ident++, blob_value,
                                                                 blob_error, bshape, slice, m_face.at(apa).at(face));
            ret.push_back(iblob);

        }
    }

    // std::cout << "Test: Blobs: " << ret.size() << std::endl;

    return ret;
}

std::set<const WireCell::Clus::Facade::Blob*> 
RetileCluster::remove_bad_blobs(const Cluster& cluster, Cluster& shad_cluster, int tick_span, int apa, int face) const
{
    const auto& wpids = cluster.grouping()->wpids();
    const auto& shad_wpids = shad_cluster.grouping()->wpids();
    // if (wpids.size() > 1 || shad_wpids.size() > 1) {
    //     throw std::runtime_error("Live or Dead grouping must have exactly one wpid: wpids.size()=" + 
    //                  std::to_string(wpids.size()) + ", shad_wpids.size()=" + 
    //                  std::to_string(shad_wpids.size()));
    // }
    
    // Since wpids is a set, we need to get the first element using an iterator
    // WirePlaneId wpid = *wpids.begin();
    // WirePlaneId shad_wpid = *shad_wpids.begin();
    // if (wpid != shad_wpid) {
    //     throw std::runtime_error("Live and Dead grouping must have the same wpid");
    // }
    // int apa = wpid.apa();
    // int face = wpid.face();

    // Implementation here
    // Get time-organized map of original blobs
    const auto& orig_time_blob_map = cluster.time_blob_map().at(apa).at(face);
    
    // Get time-organized map of newly created blobs
    const auto& new_time_blob_map = shad_cluster.time_blob_map().at(apa).at(face);
    
    // Track blobs that need to be removed
    std::set<const Blob*> blobs_to_remove;

    // Examine each new blob
    for (const auto& [time_slice, new_blobs] : new_time_blob_map) {
        // std::cout << time_slice << " " << new_blobs.size() << std::endl;

        for (const Blob* new_blob : new_blobs) {
            bool flag_good = false;
            
            // Check overlap with blobs in previous time slice
            if (orig_time_blob_map.find(time_slice - tick_span) != orig_time_blob_map.end()) {
                for (const Blob* orig_blob : orig_time_blob_map.at(time_slice - tick_span)) {
                    if (new_blob->overlap_fast(*orig_blob, 1)) {
                        flag_good = true;
                        break;
                    }
                }
            }
            
            // Check overlap with blobs in same time slice
            if (!flag_good && orig_time_blob_map.find(time_slice) != orig_time_blob_map.end()) {
                for (const Blob* orig_blob : orig_time_blob_map.at(time_slice)) {
                    if (new_blob->overlap_fast(*orig_blob, 1)) {
                        flag_good = true;
                        break;
                    }
                }
            }
            
            // Check overlap with blobs in next time slice
            if (!flag_good && orig_time_blob_map.find(time_slice + tick_span) != orig_time_blob_map.end()) {
                for (const Blob* orig_blob : orig_time_blob_map.at(time_slice + tick_span)) {
                    if (new_blob->overlap_fast(*orig_blob, 1)) {
                        flag_good = true;
                        break;
                    }
                }
            }
            
            // If no overlap found with original blobs in nearby time slices, mark for removal
            if (!flag_good) {
                blobs_to_remove.insert(new_blob);
            }
        }
    }
    
    // Remove the bad blobs
    return blobs_to_remove;
   
    
}

          
Points::node_ptr RetileCluster::mutate(Points::node_type& node) const
{
    auto* orig_cluster = reinitialize(node);
    if (!orig_cluster) {
        return nullptr;
    }

    // Only retile clusters with flashes that are in the window
    auto flash = orig_cluster->get_flash();
    if (!flash) {
        return nullptr;
    }
    double flash_time = flash.time();
    if (! (flash_time >= m_cut_time_low && flash_time <= m_cut_time_high)) {
        return nullptr;
    }
                
    // get the span of indices
    auto cc = orig_cluster->get_pcarray("isolated", "perblob");
    // convert span to vector
    std::vector<int> cc_vec(cc.begin(), cc.end());
    // for (const auto& val : cc_vec) {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;

    auto scope = orig_cluster->get_default_scope();
    auto scope_transform = orig_cluster->get_scope_transform(scope);
    // origi_cluster still have the original main cluster ... 
    // debug_cluster(orig_cluster, "Start:");

    // std::cout << "Xin1:  " << orig_cluster->get_scope_filter(scope) << " " << orig_cluster->get_default_scope() << " " << orig_cluster->get_scope_transform(scope) << std::endl;

    auto splits = m_grouping->separate(orig_cluster, cc_vec);
    // debug_cluster(orig_cluster, "Mid:");

    // std::cout << "Xin2:  " << orig_cluster->get_scope_filter(scope) << " " << orig_cluster->get_default_scope() << " " << orig_cluster->get_scope_transform(scope)<< std::endl;


    // orig_cluster->set_scope_filter(scope, true);
    orig_cluster->set_default_scope(scope); // need this to clear cache ... 
    // orig_cluster->set_scope_transform(scope,scope_transform);
                
    // std::cout << "Xin:  " << orig_cluster->get_scope_filter(scope) << " " << orig_cluster->get_default_scope() << " " << orig_cluster->get_scope_transform(scope)<< std::endl;

    std::map<int, Cluster*> map_id_cluster = splits;
    map_id_cluster[-1] = orig_cluster;

    Cluster *shadow_orig_cluster=nullptr;

    // A temporary node and grouping facade to hold separate clusters that we
    // the construct into a single "shadow" cluster for return.
    Points::node_t shad_node;
    auto* shad_grouping = shad_node.value.facade<Grouping>();
    shad_grouping->from(*m_grouping); // copies $#@%# state

    std::map<int, Cluster*> shadow_splits;
    for (auto& [id, cluster] : map_id_cluster) {

        auto& shad_cluster = shad_grouping->make_child();
        shad_cluster.set_ident(cluster->ident());                    

        if (id==-1) shadow_orig_cluster = &shad_cluster;
        else shadow_splits[id] = &shad_cluster;

        // Needed in hack_activity() but call it here to avoid call overhead.
        const auto& path_wcps = cluster_path_wcps(cluster);

        auto wpids = cluster->wpids_blob();
        std::set<WirePlaneId> wpid_set(wpids.begin(), wpids.end());
        for (auto it = wpid_set.begin(); it != wpid_set.end(); ++it) {
            int apa = it->apa();
            int face = it->face();
            const auto& angles = m_wpid_angles.at(*it);

            // Step 1.
            std::map<std::pair<int, int>, std::vector<WRG::measure_t> > map_slices_measures;
            get_activity(*cluster, map_slices_measures, apa, face);

            // Step 2.
            hack_activity(*cluster, map_slices_measures, path_wcps, apa, face); // may need more args

            // Check for time slices with same start and end
            // for (const auto& [time_range, measures] : map_slices_measures) {
            //     if (time_range.first == 480 or time_range.first == 1148) {
            //         std::cout << "Warning: Time slice with same start and end found: " 
            //                   << time_range.first << " " << time_range.second << std::endl;
            //     }
            // }
                        
            // Step 3.  Must make IBlobs for this is what the sampler takes.
            auto shad_iblobs = make_iblobs(map_slices_measures, apa, face); // may need more args

            // Steps 4-6.
            auto niblobs = shad_iblobs.size();

            // This is the 3rd generation of copy-paste for sampling.  Gen 2 is
            // in UbooneClusterSource.  OG is in PointTreeBuilding.  The reason
            // for the copy-pastes is insufficient attentino to proper code
            // factoring starting in PointTreeBuilding.  Over time, it is almost
            // guaranteed these copy-pastes become out-of-sync.  A 4th copy is
            // likely found in the steiner-related area.

            for (size_t bind=0; bind<niblobs; ++bind) {
                if (!m_samplers.at(apa).at(face)) {
                    shad_cluster.make_child();
                    continue;
                }
                const IBlob::pointer iblob = shad_iblobs[bind];
                auto sampler = m_samplers.at(apa).at(face);
                const double tick = 500*units::ns;

                auto pcs = Aux::sample_live(sampler, iblob, angles, tick, bind);
                /// DO NOT EXTEND FURTHER! see #426, #430

                if (pcs.empty()) {
                    SPDLOG_DEBUG("retile: skipping blob {} with no points", iblob->ident());
                    continue;
                }
                shad_cluster.node()->insert(Tree::Points(std::move(pcs)));
            }
            int tick_span = map_slices_measures.begin()->first.second -  map_slices_measures.begin()->first.first;
            // std::cout << "Test: " << shad_cluster.npoints() << " " << " " << shad_cluster.nchildren() << std::endl;

            // remove blobs after creating facade_blobs ... 
            auto blobs_to_remove = remove_bad_blobs(*cluster, shad_cluster, tick_span, apa, face);
            for (const Blob* blob : blobs_to_remove) {
                Blob& b = const_cast<Blob&>(*blob);
                shad_cluster.remove_child(b);
            }
            // shad_cluster.clear_cache();
            // std::cout << "Test: " << apa << " " << face << " " << blobs_to_remove.size() << std::endl;
        }
                    

        // add the new scope to the newly corrected shad_cluster ...
        auto& default_scope = cluster->get_default_scope();
        auto& raw_scope = cluster->get_raw_scope();

        if (default_scope.hash()!=raw_scope.hash()){
            auto correction_name = cluster->get_scope_transform(default_scope);
            // std::vector<int> filter_results = c
            shad_cluster.add_corrected_points(m_pcts, correction_name);
            // Get the new scope with corrected points
            const auto correction_scope = shad_cluster.get_scope(correction_name);
            // // Set this as the default scope for viewing
            shad_cluster.from(*cluster); // copy state from original cluster
            // std::cout << "Test: Same:" << default_scope.hash() << " " << raw_scope.hash() << std::endl; 
        }

    }


    // Restore input cluster
    auto cc2 = m_grouping->merge(splits,orig_cluster);

    // Record how we had split it.
    orig_cluster->put_pcarray(cc2, "isolated", "perblob");

    // Merge the separate shadow clusters into one.
    auto cc3 = shad_grouping->merge(shadow_splits, shadow_orig_cluster);

    // Record its splits
    shadow_orig_cluster->put_pcarray(cc3, "isolated", "perblob");


    // Send merged cluster node to caller.
    return shad_node.remove(shad_node.children().front());
}


