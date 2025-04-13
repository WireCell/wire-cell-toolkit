// Developers: see important comments in header file.
//
// A "FIXME" really must be fixed if you expect anything to work.  Lower case
// "fixme" are "normal" fixmes that everyone ignores.


#include "WireCellClus/ClusteringRetile.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SamplingHelpers.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PointTree.h"

#include "WireCellAux/SimpleSlice.h"
#include "WireCellClus/GroupingHelper.h"


using namespace WireCell;

// Segregate this weird choice for namespace.
namespace WCC = WireCell::PointCloud::Facade;

// Nick name for less typing.
namespace WRG = WireCell::RayGrid;


// Now can handle all APA/Faces 
WCC::ClusteringRetile::ClusteringRetile(const WireCell::Configuration& cfg)
{
    // auto sampler = get<std::string>(cfg, "sampler","");
    // if (sampler.empty()) {
    //     raise<ValueError>("ClusteringRetile requires an IBlobSampler type/name in 'sampler' parameter");
    // }
    // std::cout << "Test: " << sampler << std::endl;
    // m_sampler = Factory::find_tn<IBlobSampler>(sampler);

    if (cfg.isMember("samplers") && cfg["samplers"].isArray()) {
        // Process array of samplers
        for (const auto& sampler_cfg : cfg["samplers"]) {
            int apa = sampler_cfg["apa"].asInt();
            int face = sampler_cfg["face"].asInt();
            std::string sampler_name = sampler_cfg["name"].asString();
            
            if (sampler_name.empty()) {
                raise<ValueError>("ClusteringRetile requires an IBlobSampler name for APA %d face %d", apa, face);
            }
            // std::cout << "Test: " << apa << " " << face << " " << sampler_name << std::endl;
            auto sampler_ptr = Factory::find_tn<IBlobSampler>(sampler_name);
            m_samplers[apa][face] = sampler_ptr;
        }
    }

    // auto anode_tn = get<std::string>(cfg, "anode","");
    // if (anode_tn.empty()) {
    //     raise<ValueError>("ClusteringRetile requires an IAnodePlane type/name in 'anode' parameter");
    // }
    // int face_index = get(cfg, "face", 0);
    // auto anode = Factory::find_tn<IAnodePlane>(anode_tn);
    // m_face = anode->faces()[face_index];
    // if (!m_face) {
    //     raise<ValueError>("ClusteringRetile got null IAnodeFace at index=%d from %s", face_index, anode_tn);
    // }

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

    // Get the detector volumes pointer
    m_dv = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());
    if (m_dv == nullptr) {
        raise<ValueError>("failed to get IDetectorVolumes %s", cfg["detector_volumes"].asString());
    }
}


// Step 1. Build activities from blobs in a cluster.
void WCC::ClusteringRetile::get_activity(const Cluster& cluster, std::map<std::pair<int, int>, std::vector<WRG::measure_t> >& map_slices_measures, int apa, int face) const
{
    const int nlayers = 2+3;

    // checkme: this assumes "iend" is the usual one-past-last aka [ibeg,iend)
    // forms a half-open range.  I'm not sure if PointTreeBuilding is following
    // this or not.

    

    // for (auto& info : plane_infos) {
    //     std::cout << "test1: " << info.start_index << " " << info.end_index << " " << info.total_wires << std::endl;
    // }

    int (WCC::Blob::*wmin[])(void) const = {
        &WCC::Blob::u_wire_index_min,
        &WCC::Blob::v_wire_index_min,
        &WCC::Blob::w_wire_index_min
    };

    int (WCC::Blob::*wmax[])(void) const = {
        &WCC::Blob::u_wire_index_max,
        &WCC::Blob::v_wire_index_max,
        &WCC::Blob::w_wire_index_max
    };
        
    const double hit=1.0;       // actual charge value does not matter to tiling.

    for (const auto* fblob : cluster.children()) {
        int tslice_beg = fblob->slice_index_min();
        int tslice_end = fblob->slice_index_max();

        // if blob is not consistent skip ...
        auto blob_wpid = fblob->wpid();
        if (blob_wpid.apa()!=apa || blob_wpid.face()!=face) continue;

        auto& measures = map_slices_measures[std::make_pair(tslice_beg, tslice_end)];
        
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
void WCC::ClusteringRetile::hack_activity(const Cluster& cluster, std::map<std::pair<int, int>, std::vector<WRG::measure_t> >& map_slices_measures, int apa, int face) const
{
    const double low_dis_limit = 0.3 * units::cm;
    // Get path points
    auto path_wcps = cluster.get_path_wcps();
    std::vector<std::pair<geo_point_t, WirePlaneId>> path_pts;

    // Convert list points to vector with interpolation
    for (const auto& wcp : path_wcps) {
        geo_point_t p= cluster.point3d_raw(wcp); // index ... // raw data points ...
        auto wpid_p = cluster.wire_plane_id(wcp); // wpid ...
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
    }

    // this is to get the end of the time tick range = start_tick + tick_span
    int tick_span = map_slices_measures.begin()->first.second -  map_slices_measures.begin()->first.first;


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
                    measures[wire] = 1.0;  // Set activity
                }
            }
        }
    }


    // std::cout << "Test: Alt: " << map_slices_measures.size() << " " << cluster.children().size() << std::endl;

    // for (auto it = map_slices_measures.begin(); it!= map_slices_measures.end(); it++){
    //     std::cout << it->first.first << " " << it->first.second << " " << it->second.size() << std::endl;
    //     // for (int i=0; i!=5; i++){
    //     //     std::cout << it->second[i].size() << " ";
    //     // }
    //     // std::cout << std::endl;
    // }

}



// Step 3. Form IBlobs from activities.
std::vector<IBlob::pointer> WCC::ClusteringRetile::make_iblobs(std::map<std::pair<int, int>, std::vector<WRG::measure_t> >& map_slices_measures, int apa, int face) const
{
    std::vector<IBlob::pointer> ret;

    const auto& coords = m_face.at(apa).at(face)->raygrid();
    int blob_ident=0;
    int slice_ident = 0;
    for (auto it = map_slices_measures.begin(); it != map_slices_measures.end(); it++){
        // Do the actual tiling.
        WRG::activities_t activities = RayGrid::make_activities(m_face.at(apa).at(face)->raygrid(), it->second);
        auto bshapes = WRG::make_blobs(coords, activities);

        //std::cout << "abc: " << bshapes.size() << " " << activities.size() << " " << std::endl;
        // for (const auto& activity : activities) {
        //     std::cout << activity.as_string() << std::endl;
        // }

        // Convert RayGrid blob shapes into IBlobs 
        const float blob_value = 0.0;  // tiling doesn't consider particular charge
        const float blob_error = 0.0;  // tiling doesn't consider particular charge
    
        for (const auto& bshape : bshapes) {
            IFrame::pointer sframe = nullptr;

            // 500 ns should be passed from outside?
            ISlice::pointer slice = std::make_shared<Aux::SimpleSlice>(sframe, slice_ident++, it->first.first*500*units::ns, (it->first.second - it->first.first)*500*units::ns);

            //     ISlice::pointer slice = nullptr; // fixme: okay?
            IBlob::pointer iblob = std::make_shared<Aux::SimpleBlob>(blob_ident++, blob_value,
                                                                 blob_error, bshape, slice, m_face.at(apa).at(face));
            //     std::cout << "Test: " << iblob << std::endl;
            // FIXME: (maybe okay?) GridTiling produces an IBlobSet here which holds
            // ISlice info.  Are we losing anything important not including that
            // info?
            ret.push_back(iblob);
        }
    }

    // std::cout << "Test: Blobs: " << ret.size() << std::endl;

    return ret;
}

std::set<const WireCell::PointCloud::Facade::Blob*> 
WireCell::PointCloud::Facade::ClusteringRetile::remove_bad_blobs(const Cluster& cluster, Cluster& shad_cluster, int tick_span, int apa, int face) const
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



void WCC::ClusteringRetile::operator()(WCC::Grouping& original, WCC::Grouping& shadow, cluster_set_t&) const
{


    // FIXME: With #377 fixed, we would make the shadow grouping here from
    // scratch and add it to the input map by name, eg "shadow".  Instead, we
    // smash whatever is in the 2nd grouping to fill with "shadow" clusters.  In
    // other cluster functions this second grouping is interpreted as holding
    // "dead" clusters.

    // std::cout << "Test: " << original.wpids().size() << std::endl;
    // Check that live_grouping has exactly one wpid
    // if (original.wpids().size() > 1 ) {
    //     throw std::runtime_error("Live or Dead grouping must have exactly one wpid");
    // }
    auto wpids = original.wpids();

    // Example usage in clustering_parallel_prolong()
    // Find the first valid WirePlaneId in the grouping
    std::map<WirePlaneId , std::tuple<geo_point_t, double, double, double>> wpid_params;
    std::set<int> apas;
    for (const auto& gwpid : original.wpids()) {
        int apa = gwpid.apa();
        int face = gwpid.face();
        apas.insert(apa);

        // Create wpids for all three planes with this APA and face
        WirePlaneId wpid_u(kUlayer, face, apa);
        WirePlaneId wpid_v(kVlayer, face, apa);
        WirePlaneId wpid_w(kWlayer, face, apa);
     
        // Get drift direction based on face orientation
        int face_dirx = m_dv->face_dirx(wpid_u);
        geo_point_t drift_dir(face_dirx, 0, 0);
        
        // Get wire directions for all planes
        Vector wire_dir_u = m_dv->wire_direction(wpid_u);
        Vector wire_dir_v = m_dv->wire_direction(wpid_v);
        Vector wire_dir_w = m_dv->wire_direction(wpid_w);

        // Calculate angles
        double angle_u = std::atan2(wire_dir_u.z(), wire_dir_u.y());
        double angle_v = std::atan2(wire_dir_v.z(), wire_dir_v.y());
        double angle_w = std::atan2(wire_dir_w.z(), wire_dir_w.y());

        wpid_params[gwpid] = std::make_tuple(drift_dir, angle_u, angle_v, angle_w);
    }


    // reset the shadown clusters' content ... 
    // std::cout << shadow.children().size() << std::endl;
    shadow.local_pcs().clear();
    for (auto* fcluster : shadow.children()) {
        shadow.remove_child(*fcluster);
    }
    shadow.clear_cache();
    // std::cout << shadow.children().size() << std::endl;
    // const auto [angle_u,angle_v,angle_w] = original.wire_angles();


    for (auto* orig_cluster : original.children()) {
        // auto& scope = orig_cluster->get_default_scope();
        // auto& scope_raw = orig_cluster->get_raw_scope();

        // if (scope.hash()!=scope_raw.hash()){
        //     throw std::runtime_error("live grouping must have the raw points scope as the default");
        // }

        // find the flash time:
        auto flash = orig_cluster->get_flash();
        // int nblobs =
        // orig_cluster->kd_blobs().size();
        
        // Apply time cut
        if (flash) {
            double flash_time = flash.time();
            // std::cout << "Test: " << flash_time << " " << std::endl;

            if (flash_time >= m_cut_time_low && flash_time <= m_cut_time_high) {
                // std::cout << "Tests: " << nblobs << " at time " << flash_time/units::us << " " << m_cut_time_low/units::us << " " << m_cut_time_high/units::us << "\n";
                
                // get the span of indices
                auto cc = orig_cluster->get_pcarray("isolated", "perblob");
                // convert span to vector
                std::vector<int> cc_vec(cc.begin(), cc.end());
                // for (const auto& val : cc_vec) {
                //     std::cout << val << " ";
                // }
                // std::cout << std::endl;

                // use the vector for separate()
                // origi_cluster still have the original main cluster ... 
                auto splits = original.separate(orig_cluster, cc_vec);
             
                
                std::map<int, Cluster*> map_id_cluster = splits;
                map_id_cluster[-1] = orig_cluster;

                Cluster *shadow_orig_cluster;
                std::map<int, Cluster*> shadow_splits;

                for (auto& [id, cluster] : map_id_cluster) {

                    // make a shadow cluster, insert ID ...
                    auto& shad_cluster = shadow.make_child();
                    shad_cluster.set_ident(cluster->ident());                    
                    // std::cout <<"Test: bcd: " << cluster->ident() << " " << shad_cluster.ident() << std::endl;

                    if (id==-1) shadow_orig_cluster = &shad_cluster;
                    else shadow_splits[id] = &shad_cluster;

                    // find the highest and lowest points
                    std::pair<geo_point_t, geo_point_t> pair_points = cluster->get_highest_lowest_points();
                    //std::cout << pair_points.first << " " << pair_points.second << std::endl;
                    int high_idx = cluster->get_closest_point_index(pair_points.first);
                    int low_idx = cluster->get_closest_point_index(pair_points.second);
                    cluster->dijkstra_shortest_paths(high_idx, false);
                    cluster->cal_shortest_path(low_idx);

                    auto wpids = cluster->wpids_blob();
                    std::set<WirePlaneId> wpid_set(wpids.begin(), wpids.end());
                    for (auto it = wpid_set.begin(); it != wpid_set.end(); ++it) {
                        int apa = it->apa();
                        int face = it->face();
                        auto [drift_dir, angle_u, angle_v, angle_w] = wpid_params.at(*it);
                        // std::cout << "Test: " << apa << " " << face << " " << angle_u << " " << angle_v << " " << angle_w << std::endl;

                        // Step 1.
                        std::map<std::pair<int, int>, std::vector<WRG::measure_t> > map_slices_measures;
                        get_activity(*cluster, map_slices_measures, apa, face);

                        // Step 2.
                        hack_activity(*cluster, map_slices_measures, apa, face); // may need more args
                        
                        // Step 3.  Must make IBlobs for this is what the sampler takes.
                        auto shad_iblobs = make_iblobs(map_slices_measures, apa, face); // may need more args

                        // Steps 4-6.
                        auto niblobs = shad_iblobs.size();
                        // Forgive me (and small-f fixme), but this is now the 3rd generation of
                        // copy-paste.  Gen 2 is in UbooneClusterSource.  OG is in
                        // PointTreeBuilding.  The reason for the copy-pastes is insufficient
                        // factoring of the de-factor standard sampling code in PointTreeBuilding.
                        // Over time, it is almost guaranteed these copy-pastes become out-of-sync. 

                        for (size_t bind=0; bind<niblobs; ++bind) {
                            if (!m_samplers.at(apa).at(face)) {
                                shad_cluster.make_child();
                                continue;
                            }
                            const IBlob::pointer iblob = shad_iblobs[bind];

                            // Sample the iblob, make a new blob node.
                            PointCloud::Tree::named_pointclouds_t pcs;

                            auto [pc3d, aux] = m_samplers.at(apa).at(face)->sample_blob(iblob, bind);
                            
                            // how to sample points ... 
                            // std::cout << pc3d.size() << " " << aux.size() << " " <<  pc3d.get("x")->size_major() << " " << pc3d.get("y")->size_major() << " " << pc3d.get("z")->size_major() << std::endl;
                            // const auto& arr_x1 = pc3d.get("x")->elements<Point::coordinate_t>();

                            /// These seem unused and bring in yet more copy-paste code
                            // pcs.emplace("2dp0", WireCell::Aux::make2dds(pc3d, angle_u));
                            // pcs.emplace("2dp1", WireCell::Aux::make2dds(pc3d, angle_v));
                            // pcs.emplace("2dp2", WireCell::Aux::make2dds(pc3d, angle_w));
                            auto pc2dp0 = WireCell::Aux::make2dds(pc3d, angle_u);
                            auto pc2dp1 = WireCell::Aux::make2dds(pc3d, angle_v);
                            auto pc2dp2 = WireCell::Aux::make2dds(pc3d, angle_w);
                            pc3d.add("2dp0_x", *pc2dp0.get("x"));
                            pc3d.add("2dp0_y", *pc2dp0.get("y"));
                            pc3d.add("2dp1_x", *pc2dp1.get("x"));
                            pc3d.add("2dp1_y", *pc2dp1.get("y"));
                            pc3d.add("2dp2_x", *pc2dp2.get("x"));
                            pc3d.add("2dp2_y", *pc2dp2.get("y"));
                            pcs.emplace("3d", pc3d);
                            // std::cout << pcs["3d"].get("x")->size_major() << " " << pcs["3d"].get("y")->size_major() << " " << pcs["3d"].get("z")->size_major() << std::endl;
                            // const auto& arr_x = pcs["3d"].get("x")->elements<Point::coordinate_t>();
                            // std::cout << arr_x.size() << " " << arr_x1.size() << std::endl;
                            // std::cout << iblob->shape() << std::endl;
                            if (pc3d.get("x")->size_major() > 0){
                                const Point center = WireCell::Aux::calc_blob_center(pcs["3d"]);
                                auto scalar_ds = WireCell::Aux::make_scalar_dataset(iblob, center, pcs["3d"].get("x")->size_major(), 500*units::ns);
                                int max_wire_interval = aux.get("max_wire_interval")->elements<int>()[0];
                                int min_wire_interval = aux.get("min_wire_interval")->elements<int>()[0];
                                int max_wire_type = aux.get("max_wire_type")->elements<int>()[0];
                                int min_wire_type = aux.get("min_wire_type")->elements<int>()[0];
                                scalar_ds.add("max_wire_interval", Array({(int)max_wire_interval}));
                                scalar_ds.add("min_wire_interval", Array({(int)min_wire_interval}));
                                scalar_ds.add("max_wire_type", Array({(int)max_wire_type}));
                                scalar_ds.add("min_wire_type", Array({(int)min_wire_type}));
                                pcs.emplace("scalar", std::move(scalar_ds));

                                shad_cluster.node()->insert(Tree::Points(std::move(pcs)));
                            }else{
                                SPDLOG_WARN("blob {} has no points", iblob->ident());
                            }
                        }
                        int tick_span = map_slices_measures.begin()->first.second -  map_slices_measures.begin()->first.first;
                            // std::cout << "Test: " << shad_cluster.npoints() << " " << " " << shad_cluster.nchildren() << std::endl;

                        // remove blobs after creating facade_blobs ... 
                        auto blobs_to_remove = remove_bad_blobs(*cluster, shad_cluster, tick_span, apa, face);
                        for (const Blob* blob : blobs_to_remove) {
                            Blob& b = const_cast<Blob&>(*blob);
                            shad_cluster.remove_child(b);
                        }
                        shad_cluster.clear_cache();
                    }
                    
                    // // Reset cached data that depends on cluster contents
                    // shad_cluster.reset_pca();         // Reset PCA calculations
                    // // Force rebuild of time blob map by accessing it
                    // shad_cluster.time_blob_map();
                    // shad_cluster.point3d(0); // This will trigger PC tree rebuild
                    // std::cout << shad_cluster.npoints() << " " << shad_cluster.nbpoints() << " " << shad_cluster.nchildren() << std::endl;

                    // How to call overlap_fast ??? 
                    // for (auto* fblob : shad_cluster.children()) {
                    //     shad_cluster.remove_child(*fblob);
                    // }

                    // std::cout << "Test: remove: "  << " " << cluster->kd_blobs().size() << " " << shad_cluster.kd_blobs().size()  << std::endl;
               
                    // Example code to access shadown cluster information ...
                    // // shad cluster getting highest and lowest points and then do shortest path ... 
                    // // find the highest and lowest points
                    // std::pair<geo_point_t, geo_point_t> shad_pair_points = shad_cluster.get_highest_lowest_points();
                    // //std::cout << pair_points.first << " " << pair_points.second << std::endl;
                    // int shad_high_idx = shad_cluster.get_closest_point_index(shad_pair_points.first);
                    // int shad_low_idx = shad_cluster.get_closest_point_index(shad_pair_points.second);
                    // shad_cluster.dijkstra_shortest_paths(shad_high_idx, false);
                    // shad_cluster.cal_shortest_path(shad_low_idx);
                    // {
                    //     auto path_wcps = shad_cluster.get_path_wcps();                
                    //     // Convert list points to vector with interpolation
                    //     for (const auto& wcp : path_wcps) {
                    //         geo_point_t p= shad_cluster.point3d(wcp);
                    //         std::cout << p << std::endl;
                    //     }
                    // }

                    // add the new scope to the newly corrected shad_cluster ...
                    auto& default_scope = cluster->get_default_scope();
                    auto& raw_scope = cluster->get_raw_scope();

                    if (default_scope.hash()!=raw_scope.hash()){
                        auto correction_name = cluster->get_scope_transform(default_scope);
                        // std::vector<int> filter_results = c
                        shad_cluster.add_corrected_points(m_dv, correction_name);
                        // Get the new scope with corrected points
                        const auto correction_scope = shad_cluster.get_scope(correction_name);
                        // // Set this as the default scope for viewing
                        shad_cluster.set_default_scope(correction_scope);
                        shad_cluster.set_scope_transform(correction_scope, correction_name);
                    }

                }

                //     // FIXME: These two methods need to be added to the Cluster Facade.
                //     // They should set/get "cluster_id" from the "cluster_scalar" PC.
                //     // FIXME: above we add cluster_id to the "cluster_scalar" PC.  Do we
                //     // want to also try to copy over the "light" index entry?

                auto cc2 = original.merge(splits,orig_cluster);
                // for (const auto& val : cc2) {
                //     std::cout << val << " ";
                // }
                // std::cout << std::endl;
                orig_cluster->put_pcarray(cc2, "isolated", "perblob");

                auto cc3 = shadow.merge(shadow_splits,shadow_orig_cluster);
                // for (const auto& val : cc3) {
                //     std::cout << val << " ";
                // }
                // std::cout << std::endl;
                shadow_orig_cluster->put_pcarray(cc3, "isolated", "perblob");

            }
            // FIXME: do we need/want to copy over any PCs in the grouping?
            // Specifically optical light/flash/flashlight PCs?

        }
    }

     // set cluster id ... 
    int cluster_id = 1;
    for (auto* cluster : original.children()) {
        cluster->set_cluster_id(cluster_id++);
    }

    // // Process groupings after all shadow clusters are created
    // auto cluster_mapping = process_groupings(original, shadow);
    //     int cluster_id = 1;
    // for (const auto& [orig_cluster, tuple] : cluster_mapping) {
    //     std::cout << orig_cluster << " " << std::get<0>(tuple) << " " << std::get<1>(tuple) << " " << std::get<2>(tuple) << std::endl;
    // //     auto* shad_cluster = pair.first;
    // //     auto* main_cluster = pair.second;
    // shad_cluster->set_cluster_id(cluster_id);
    // main_cluster->set_cluster_id(cluster_id);
    // cluster_id ++;
    // //     // You can now use these mapped clusters for further processing
    // //     // For example, transfer any necessary properties or perform additional operations
    // }

   
}

std::map<WCC::Cluster*, std::tuple<WCC::Cluster*, int, WCC::Cluster*>> 
WCC::ClusteringRetile::process_groupings(
    WCC::Grouping& original,
    WCC::Grouping& shadow,
    const std::string& aname,
    const std::string& pname) const
{
  return process_groupings_helper(original, shadow, aname, pname);
}
