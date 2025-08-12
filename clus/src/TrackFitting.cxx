#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/TrackFitting_Util.h"
#include "WireCellUtil/Logging.h"

#include <Eigen/IterativeLinearSolvers>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

using geo_point_t = WireCell::Point;


TrackFitting::TrackFitting(FittingType fitting_type) 
    : m_fitting_type(fitting_type) 
{

}

void TrackFitting::add_segment(std::shared_ptr<PR::Segment> segment){
    m_segments.insert(segment);
    m_clusters.insert(segment->cluster());

    if (m_grouping == nullptr){
        m_grouping = segment->cluster()->grouping();

        BuildGeometry();
    }

    for (auto& cluster: m_clusters){
        for (auto& blob: cluster->children()){
            m_blobs.insert(blob);
        }
    }

    std::cout << "TrackFitting: Added segment with " << segment->wcpts().size() << " points." << " " << m_clusters.size() << " " << m_blobs.size() << std::endl;
}

void TrackFitting::BuildGeometry(){
    // Get all the wire plane IDs from the grouping
    const auto& wpids = m_grouping->wpids();    
    compute_wireplane_params(wpids, m_dv, wpid_params, wpid_U_dir, wpid_V_dir, wpid_W_dir, apas);

    // Clear existing maps
    wpid_offsets.clear();
    wpid_slopes.clear();
     // Get all unique APA/face combinations
    std::set<std::pair<int, int>> apa_face_combinations;

    // loop over wpids ...
    for (const auto& wpid : wpids) {
        double time_slice_width = //m_dv->metadata(wpid)["nticks_live_slice"].asDouble() *  
        m_dv->metadata(wpid)["tick_drift"].asDouble();

        WirePlaneId wpid_u(kUlayer, wpid.face(), wpid.apa());
        WirePlaneId wpid_v(kVlayer, wpid.face(), wpid.apa());
        WirePlaneId wpid_w(kWlayer, wpid.face(), wpid.apa());

        double pitch_u = m_dv->pitch_vector(wpid_u).magnitude();
        double pitch_v = m_dv->pitch_vector(wpid_v).magnitude();
        double pitch_w = m_dv->pitch_vector(wpid_w).magnitude();

        wpid_geoms[wpid] = std::make_tuple(time_slice_width, pitch_u, pitch_v, pitch_w);
        // std::cout << "Geometry: " << time_slice_width/units::cm << " " << pitch_u/units::cm << " " << pitch_v/units::cm << " " << pitch_w/units::cm << std::endl;

        apa_face_combinations.insert({wpid.apa(), wpid.face()});
    }
    
    // Process each APA/face combination
    for (const auto& [apa, face] : apa_face_combinations) {
        try {
            // Get anode interface for this APA/face
            auto anode = m_grouping->get_anode(apa);
            if (!anode) {
                std::cerr << "TrackFitting: Could not get anode for APA " << apa << std::endl;
                continue;
            }
            
            auto iface = anode->faces()[face];
            if (!iface) {
                std::cerr << "TrackFitting: Could not get face " << face << " for APA " << apa << std::endl;
                continue;
            }
            
            // Get geometry parameters from grouping
            const auto& pitch_mags = m_grouping->pitch_mags();
            const auto& proj_centers = m_grouping->proj_centers();
            
            // Get wire angles for this APA/face
            const auto [angle_u, angle_v, angle_w] = m_grouping->wire_angles(apa, face);
            std::vector<double> angles = {angle_u, angle_v, angle_w};
            
            // Get time/drift parameters from grouping cache
            double time_offset = m_grouping->get_time_offset().at(apa).at(face);
            double drift_speed = m_grouping->get_drift_speed().at(apa).at(face);
            double tick = m_grouping->get_tick().at(apa).at(face);
            
            // Get drift direction and origin from anode face
            double xsign = iface->dirx();
            double xorig = iface->planes()[2]->wires().front()->center().x();
            
            // Create WirePlaneId for this APA/face combination
            WirePlaneId wpid(kAllLayers, face, apa);
            
            // Calculate slopes and offsets for each plane
            std::pair<double, double> slope_yu_zu, slope_yv_zv, slope_yw_zw;
            double offset_u, offset_v, offset_w, offset_t;
            
            // U plane (plane index 0)
            double pitch_u = pitch_mags.at(apa).at(face).at(0);
            double center_u = proj_centers.at(apa).at(face).at(0);
            offset_u = -(center_u + 0.5 * pitch_u) / pitch_u;
            slope_yu_zu = {-sin(angles[0]) / pitch_u, cos(angles[0]) / pitch_u};
            
            // V plane (plane index 1)
            double pitch_v = pitch_mags.at(apa).at(face).at(1);
            double center_v = proj_centers.at(apa).at(face).at(1);
            offset_v = -(center_v + 0.5 * pitch_v) / pitch_v;
            slope_yv_zv = {-sin(angles[1]) / pitch_v, cos(angles[1]) / pitch_v};
            
            // W plane (plane index 2)
            double pitch_w = pitch_mags.at(apa).at(face).at(2);
            double center_w = proj_centers.at(apa).at(face).at(2);
            offset_w = -(center_w + 0.5 * pitch_w) / pitch_w;
            slope_yw_zw = {-sin(angles[2]) / pitch_w, cos(angles[2]) / pitch_w};
            
            // Time conversion parameters
            // From drift2time: time = (drift - xorig)/(xsign * drift_speed) - time_offset
            // tick_index = round(time / tick)
            double slope_t = 1.0 / (xsign * drift_speed * tick);
            offset_t = -(xorig / (xsign * drift_speed) + time_offset) / tick;
            
            // Store in maps
            wpid_offsets[wpid] = std::make_tuple(offset_t, offset_u, offset_v, offset_w);
            wpid_slopes[wpid] = std::make_tuple(
                slope_t,           // T slope (for x direction)
                slope_yu_zu,       // U plane slopes (y, z)
                slope_yv_zv,       // V plane slopes (y, z)
                slope_yw_zw        // W plane slopes (y, z)
            );
            
            // // Debug output (optional - can be removed)
            // std::cout << "TrackFitting: Initialized geometry for APA " << apa 
            //           << " Face " << face << std::endl;
            // std::cout << "  Offsets: T=" << offset_t << " U=" << offset_u 
            //           << " V=" << offset_v << " W=" << offset_w << std::endl;
            // std::cout << "  Slopes: T=" << slope_t 
            //           << " U=(" << slope_yu_zu.first << "," << slope_yu_zu.second << ")"
            //           << " V=(" << slope_yv_zv.first << "," << slope_yv_zv.second << ")"
            //           << " W=(" << slope_yw_zw.first << "," << slope_yw_zw.second << ")" << std::endl;
                      
        } catch (const std::exception& e) {
            std::cerr << "TrackFitting: Error initializing geometry for APA " << apa 
                      << " Face " << face << ": " << e.what() << std::endl;
        }
    }
    
    // std::cout << "TrackFitting: Geometry initialization complete. Processed " 
    //           << wpid_offsets.size() << " wire plane configurations." << std::endl;


}

geo_point_t TrackFitting::adjust_rough_path(PR::Segment& segment){
    return geo_point_t(0,0,0);
}

IAnodePlane::pointer TrackFitting::get_anode(int apa_ident) const {
    if (!m_grouping) {
        std::cerr << "TrackFitting: No grouping available to get anode" << std::endl;
        return nullptr;
    }
    
    try {
        return m_grouping->get_anode(apa_ident);
    } catch (const std::exception& e) {
        std::cerr << "TrackFitting: Error getting anode " << apa_ident << ": " << e.what() << std::endl;
        return nullptr;
    }
}

std::map<int, IAnodePlane::pointer> TrackFitting::get_all_anodes() const {
    std::map<int, IAnodePlane::pointer> result;
    
    if (!m_grouping) {
        return result;
    }
    
    // Get all unique APAs from the clusters
    std::set<int> apa_idents;
    // Extract APAs from cluster's wire plane IDs
    auto wpids = m_grouping->wpids();
    for (const auto& wpid : wpids) {
        apa_idents.insert(wpid.apa());
    }
    
    // Get anode for each APA
    for (int apa_ident : apa_idents) {
        auto anode = get_anode(apa_ident);
        if (anode) {
            result[apa_ident] = anode;
        }
    }
    
    return result;
}

int TrackFitting::get_channel_for_wire(int apa, int face, int plane, int wire) const {
    m_cache_stats.total_lookups++;
    
    PlaneKey plane_key = std::make_tuple(apa, face, plane);
    
    // Check hot cache first (O(1) for frequently accessed planes)
    auto hot_it = m_hot_cache.find(plane_key);
    if (hot_it != m_hot_cache.end()) {
        if (wire >= 0 && wire < static_cast<int>(hot_it->second.size())) {
            m_cache_stats.hot_hits++;
            return hot_it->second[wire];
        }
        return -1; // Wire index out of bounds
    }
    
    // Check cold cache (individual wire lookups)
    WireKey wire_key = std::make_tuple(apa, face, plane, wire);
    auto cold_it = m_cold_cache.find(wire_key);
    if (cold_it != m_cold_cache.end()) {
        m_cache_stats.cold_hits++;
        
        // Update access count for this plane
        m_access_count[plane_key]++;
        
        // Promote to hot cache if threshold reached
        if (m_access_count[plane_key] >= HOT_THRESHOLD) {
            cache_entire_plane(apa, face, plane);
        }
        
        return cold_it->second;
    }
    
    // Cache miss - fetch from anode and cache result
    int channel = fetch_channel_from_anode(apa, face, plane, wire);
    if (channel != -1) {
        m_cold_cache[wire_key] = channel;
        m_access_count[plane_key]++;
        m_cache_stats.cold_entries_count++;
    }
    
    return channel;
}

std::vector<std::tuple<int, int, int>> TrackFitting::get_wires_for_channel(int apa, int channel_number) const {
    std::vector<std::tuple<int, int, int>> result;
    
    auto anode = get_anode(apa);
    if (!anode) {
        return result;
    }
    
    // Get all wires for this channel (handles wrapped wires)
    auto wires = anode->wires(channel_number);
    
    for (const auto& wire : wires) {
        auto wpid = wire->planeid();
        result.emplace_back(wpid.face(), wpid.index(), wire->index());
    }
    
    return result;
}

void TrackFitting::clear_cache() const {
    m_hot_cache.clear();
    m_cold_cache.clear();
    m_access_count.clear();
    m_cache_stats = {0, 0, 0, 0, 0};
}

TrackFitting::CacheStats TrackFitting::get_cache_stats() const {
    auto stats = m_cache_stats;
    stats.hot_planes_count = m_hot_cache.size();
    stats.cold_entries_count = m_cold_cache.size();
    return stats;
}

void TrackFitting::cache_entire_plane(int apa, int face, int plane) const {
    auto anode = get_anode(apa);
    if (!anode) return;
    
    const auto& faces = anode->faces();
    if (face >= static_cast<int>(faces.size()) || !faces[face]) return;
    
    const auto& planes = faces[face]->planes();
    if (plane >= static_cast<int>(planes.size())) return;
    
    const auto& wires = planes[plane]->wires();
    PlaneKey plane_key = std::make_tuple(apa, face, plane);
    
    // Cache entire plane (this is the "hot" cache promotion)
    auto& hot_vec = m_hot_cache[plane_key];
    hot_vec.resize(wires.size());
    for (size_t i = 0; i < wires.size(); ++i) {
        hot_vec[i] = wires[i]->channel();
    }
    
    // Remove individual wire entries from cold cache to save memory
    for (size_t i = 0; i < wires.size(); ++i) {
        WireKey wire_key = std::make_tuple(apa, face, plane, static_cast<int>(i));
        if (m_cold_cache.erase(wire_key)) {
            m_cache_stats.cold_entries_count--;
        }
    }
    
    m_cache_stats.hot_planes_count++;
    
    std::cout << "TrackFitting: Promoted plane (" << apa << "," << face << "," << plane 
              << ") to hot cache with " << wires.size() << " wires" << std::endl;
}

int TrackFitting::fetch_channel_from_anode(int apa, int face, int plane, int wire) const {
    auto anode = get_anode(apa);
    if (!anode) return -1;
    
    const auto& faces = anode->faces();
    if (face >= static_cast<int>(faces.size()) || !faces[face]) return -1;
    
    const auto& planes = faces[face]->planes();
    if (plane >= static_cast<int>(planes.size())) return -1;
    
    const auto& wires = planes[plane]->wires();
    if (wire >= static_cast<int>(wires.size())) return -1;
    
    return wires[wire]->channel();
}



void TrackFitting::prepare_data() {
    if (m_charge_data.size()!=0) return;

    // Process every Facade::Cluster in m_clusters
    for (auto& cluster : m_clusters) {
        // Get boundary range using get_uvwt_range which returns map<WirePlaneId, tuple<int,int,int,int>>
        auto uvwt_ranges = cluster->get_uvwt_range();
        
        // Get the grouping from the cluster
        auto grouping = cluster->grouping();
        
        // Process each wpid (wire plane ID) separately
        for (const auto& [wpid, range_tuple] : uvwt_ranges) {
            int apa = wpid.apa();
            int face = wpid.face();
            
            // Get the ranges for this wpid
            // auto [u_size, v_size, w_size, t_size] = range_tuple;
            
            // Get min/max values for this specific apa/face
            auto [u_min, v_min, w_min, t_min] = cluster->get_uvwt_min(apa, face);
            auto [u_max, v_max, w_max, t_max] = cluster->get_uvwt_max(apa, face);
            
            // Process each plane (0=U, 1=V, 2=W)
            for (int plane = 0; plane < 3; ++plane) {
                int wire_min, wire_max, time_min, time_max;
                
                // Set the wire range based on plane
                switch (plane) {
                    case 0: wire_min = u_min; wire_max = u_max; break;
                    case 1: wire_min = v_min; wire_max = v_max; break;
                    case 2: wire_min = w_min; wire_max = w_max; break;
                }
                time_min = t_min;
                time_max = t_max;
                
                // Get charge information for this plane
                auto charge_map = grouping->get_overlap_good_ch_charge(
                    time_min, time_max, wire_min, wire_max, apa, face, plane);
                
                // Process each charge entry
                for (const auto& [time_wire, charge_data] : charge_map) {
                    int time_slice = time_wire.first;
                    int wire_index = time_wire.second;
                    double charge = charge_data.first;
                    double charge_err = charge_data.second;

                    int channel = fetch_channel_from_anode(apa, face, plane, wire_index);

                    // Create key for m_charge_data
                    CoordReadout data_key(apa, time_slice, channel);
                    
                    int flag = 1; // Default flag for all-live-channel case
                    
                    // Check for negative charge
                    if (charge < 0) {
                        charge = 0;
                        charge_err = 1000;
                        flag = 2;
                    }
                    
                    // Save to m_charge_data
                    m_charge_data[data_key] = {charge, charge_err, flag};
                }
            }
        }
    }
        
    for (auto& cluster : m_clusters) {
         // Get the grouping from the cluster
        auto grouping = cluster->grouping();
        // Handle dead channels - loop over all Facade::Blobs in cluster
        for (const auto* blob : cluster->children()) {
            auto wpid = blob->wpid();
            int apa = wpid.apa();
            int face = wpid.face();
            
            // Check each plane for dead channels
            for (int plane = 0; plane < 3; ++plane) {
                // Check if this plane is bad for this blob
                if (grouping->is_blob_plane_bad(blob, plane)) {
                    // Get blob properties
                    double blob_charge = blob->charge();
                    
                    // Get wire range for this plane
                    int wire_min, wire_max;
                    switch (plane) {
                        case 0: 
                            wire_min = blob->u_wire_index_min();
                            wire_max = blob->u_wire_index_max();
                            break;
                        case 1: 
                            wire_min = blob->v_wire_index_min();
                            wire_max = blob->v_wire_index_max();
                            break;
                        case 2: 
                            wire_min = blob->w_wire_index_min();
                            wire_max = blob->w_wire_index_max();
                            break;
                    }
                    
                    int num_wires = wire_max - wire_min;
                    if (num_wires <= 0) continue;
                    
                    // Get time range
                    int time_min = blob->slice_index_min();
                    // int time_max = blob->slice_index_max();
                    
                    // Process each dead pixel
                    int time_slice = time_min;
                    for (int wire_index = wire_min; wire_index < wire_max; ++wire_index) {
                        int channel = fetch_channel_from_anode(apa, face, plane, wire_index);
                        CoordReadout data_key(apa, time_slice, channel);
                            
                        // Check if content exists
                        auto it = m_charge_data.find(data_key);
                        
                        if (it == m_charge_data.end()) {
                            // No existing content
                            double charge = blob_charge / num_wires;
                            double charge_err = sqrt(pow(charge * 0.1, 2) + pow(600, 2));
                            m_charge_data[data_key] = {charge, charge_err, 0};
                        } else if (it->second.flag == 0) {
                            // Existing content with flag = 0
                            double new_charge = blob_charge / num_wires;
                            double new_charge_err = sqrt(pow(new_charge * 0.1, 2) + pow(600, 2));
                            
                            it->second.charge += new_charge;
                            it->second.charge_err = sqrt(pow(it->second.charge_err, 2) + pow(new_charge_err, 2));
                        }
                        // If flag != 0, do nothing
                    }
                }
            }
        }
    }


    std::cout << "Number of Measurements: " << m_charge_data.size() << std::endl;
}

void TrackFitting::fill_global_rb_map() {
    // Clear the global readout map first
    if (global_rb_map.size() != 0 ) return;

    auto clusters = m_grouping->children();
    // Loop over the m_grouping's clusters
    for (auto& cluster : clusters) {
        // For each cluster, loop over its blobs
        if (!cluster->get_scope_filter(cluster->get_default_scope())) continue;

        auto blobs = cluster->children();
        for (auto blob : blobs) {
            if (!blob) continue;
            
            // Get the wire plane ID for this blob to determine apa and face
            auto wpid = blob->wpid();
            int apa = wpid.apa();
            int face = wpid.face();
            
            // Get the time slice bounds for this blob
            int time_slice_min = blob->slice_index_min();
            // int time_slice_max = blob->slice_index_max();
            
            // For every blob, loop over its planes (U=0, V=1, W=2)
            for (int plane = 0; plane < 3; ++plane) {
                 if (m_grouping->is_blob_plane_bad(blob, plane)) continue;

                // Get wire bounds for this plane in the blob
                int wire_min, wire_max;
                switch (plane) {
                    case 0: // U plane
                        wire_min = blob->u_wire_index_min();
                        wire_max = blob->u_wire_index_max();
                        break;
                    case 1: // V plane
                        wire_min = blob->v_wire_index_min();
                        wire_max = blob->v_wire_index_max();
                        break;
                    case 2: // W plane
                        wire_min = blob->w_wire_index_min();
                        wire_max = blob->w_wire_index_max();
                        break;
                    default:
                        continue;
                }
                
                // Skip if no valid wire range
                if (wire_min >= wire_max) continue;
                
                // Loop over time slices in this blob
                int time_slice = time_slice_min; 
                // Loop over wire indices in this plane
                for (int wire_index = wire_min; wire_index < wire_max; ++wire_index) {
                    // Convert wire coordinates to channel using existing helper function
                    int channel = fetch_channel_from_anode(apa, face, plane, wire_index);
                    if (channel == -1) continue; // Skip invalid channels
                    
                    // Create CoordReadout key and find out its CoordReadout
                    CoordReadout coord_key(apa, time_slice, channel);
                    
                    // Fill in global_rb_map - add this blob to the set for this coordinate
                    global_rb_map[coord_key].insert(blob);
                    // std::cout << "Added blob to global_rb_map at " << coord_key.apa << " " << coord_key.time << " " << coord_key.channel  << std::endl;
                }
            }
        }
    }
    
    std::cout << "Global RB Map filled with " << global_rb_map.size() << " coordinate entries." << std::endl;
} 


std::vector<WireCell::Point> TrackFitting::organize_orig_path(std::shared_ptr<PR::Segment> segment, double low_dis_limit, double end_point_limit) {
    std::vector<WireCell::Point> pts;
    
    // Get the WCPoints from the segment
    const auto& segment_wcpts = segment->wcpts();
    if (segment_wcpts.empty()) {
        return pts;
    }
    
    // Convert WCPoints to vector for easier manipulation
    std::vector<WireCell::Point> temp_wcps_vec;
    for (const auto& wcp : segment_wcpts) {
        temp_wcps_vec.push_back(wcp.point);
    }
    
    // Fill in the beginning point ...
    {
        WireCell::Point p1 = temp_wcps_vec.front();
        WireCell::Point p2 = temp_wcps_vec.front();
        double dis1 = 0;
        for (auto it = temp_wcps_vec.begin(); it != temp_wcps_vec.end(); it++) {
            p2 = *it;
            dis1 = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
            if (dis1 > low_dis_limit) break;
        }
        if (dis1 != 0) {
            WireCell::Point extended_p1(
                p1.x() + (p1.x() - p2.x()) / dis1 * end_point_limit,
                p1.y() + (p1.y() - p2.y()) / dis1 * end_point_limit,
                p1.z() + (p1.z() - p2.z()) / dis1 * end_point_limit
            );
            pts.push_back(extended_p1);
        }
    }

    // std::cout << "Test: " <<  pts.size() << " " << pts.back() << " " << temp_wcps_vec.front() << std::endl;
    
    // Fill in the middle part
    for (size_t i = 0; i != temp_wcps_vec.size(); i++) {
        WireCell::Point p1 = temp_wcps_vec.at(i);
        
        double dis = low_dis_limit;
        if (pts.size() > 0) {
            dis = sqrt(pow(p1.x() - pts.back().x(), 2) + pow(p1.y() - pts.back().y(), 2) + pow(p1.z() - pts.back().z(), 2));
        }
        
        if (dis < low_dis_limit * 0.8) {
            continue;
        } else if (dis < low_dis_limit * 1.6) {
            pts.push_back(p1);
        } else {
            int npoints = std::round(dis / low_dis_limit);
            WireCell::Point p_save = pts.back();
            for (int j = 0; j != npoints; j++) {
                WireCell::Point p(
                    p_save.x() + (p1.x() - p_save.x()) / npoints * (j + 1),
                    p_save.y() + (p1.y() - p_save.y()) / npoints * (j + 1),
                    p_save.z() + (p1.z() - p_save.z()) / npoints * (j + 1)
                );
                pts.push_back(p);
            }
        }
    }
    
    // Fill in the end part
    {
        WireCell::Point p1 = temp_wcps_vec.back();
        WireCell::Point p2 = temp_wcps_vec.back();
        double dis1 = 0;
        for (auto it = temp_wcps_vec.rbegin(); it != temp_wcps_vec.rend(); it++) {
            p2 = *it;
            dis1 = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
            if (dis1 > low_dis_limit) break;
        }
        if (dis1 != 0) {
            WireCell::Point extended_p1(
                p1.x() + (p1.x() - p2.x()) / dis1 * end_point_limit,
                p1.y() + (p1.y() - p2.y()) / dis1 * end_point_limit,
                p1.z() + (p1.z() - p2.z()) / dis1 * end_point_limit
            );
            pts.push_back(extended_p1);
        }
    }

    // std::cout << "Test: " <<  pts.size() << " " << pts.back() << " " << temp_wcps_vec.back() << std::endl;

    
    return pts;
}

std::vector<WireCell::Point> TrackFitting::examine_end_ps_vec(std::shared_ptr<PR::Segment> segment,const std::vector<WireCell::Point>& pts, bool flag_start, bool flag_end) {
    std::list<WireCell::Point> ps_list(pts.begin(), pts.end());
    
    // get the cluster from the segment
    auto cluster = segment->cluster();
    const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
    double cluster_t0 = cluster->get_flash().time();

    if (flag_start) {
        // test start
        WireCell::Point temp_start = ps_list.front(); 
        while (ps_list.size() > 0) {
            // figure out the wpid for ps_list.front() ... 
            auto test_wpid = m_dv->contained_by(ps_list.front());

            if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                // this function takes the raw points ...
                auto temp_p_raw = transform->backward(ps_list.front(), cluster_t0, test_wpid.face(), test_wpid.apa());
                // std::cout << temp_p_raw << " " << ps_list.front() << " " << test_wpid.apa() << " " << test_wpid.face() << std::endl;
                if (m_grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*units::cm, 0, 0)) break;
            }
            temp_start = ps_list.front();
            ps_list.pop_front();
        }
        
        if (ps_list.size() > 0) {
            double dis_step = 0.2*units::cm;
            double temp_dis = sqrt(pow(temp_start.x() - ps_list.front().x(), 2) + pow(temp_start.y() - ps_list.front().y(), 2) + pow(temp_start.z() - ps_list.front().z(), 2));
            int ntest = std::round(temp_dis/dis_step);
            for (size_t i = 1; i < ntest; i++) {
                WireCell::Point test_p(temp_start.x() + (ps_list.front().x() - temp_start.x())/ntest * i,
                                       temp_start.y() + (ps_list.front().y() - temp_start.y())/ntest * i,
                                       temp_start.z() + (ps_list.front().z() - temp_start.z())/ntest * i);
                // figure out the wpid for the test_p ...
                auto test_wpid = m_dv->contained_by(test_p);
                if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                    // this function takes the raw points ...
                    auto temp_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                    if (m_grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*units::cm, 0, 0)) {
                        ps_list.push_front(test_p);
                        break;
                    }
                }
            }
        } else {
            ps_list.push_front(temp_start);
        }
    }
    
    if (flag_end) {
        WireCell::Point temp_end = ps_list.back();
        while (ps_list.size() > 0) {
            // figure out the wpid for the ps_list.back() ...
            auto test_wpid = m_dv->contained_by(ps_list.back());
            if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                //this function takes the raw points ...
                auto temp_p_raw = transform->backward(ps_list.back(), cluster_t0, test_wpid.face(), test_wpid.apa());
                if (m_grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*units::cm, 0, 0)) break;
            }
            temp_end = ps_list.back();
            ps_list.pop_back();
        }
        if (ps_list.size() > 0) {
            double dis_step = 0.2*units::cm;
            double temp_dis = sqrt(pow(temp_end.x() - ps_list.back().x(), 2) + pow(temp_end.y() - ps_list.back().y(), 2) + pow(temp_end.z() - ps_list.back().z(), 2));
            int ntest = std::round(temp_dis/dis_step);
            for (size_t i = 1; i < ntest; i++) {
                WireCell::Point test_p(temp_end.x() + (ps_list.back().x() - temp_end.x())/ntest * i,
                                       temp_end.y() + (ps_list.back().y() - temp_end.y())/ntest * i,
                                       temp_end.z() + (ps_list.back().z() - temp_end.z())/ntest * i);

                auto test_wpid = m_dv->contained_by(test_p);
                // figure out the wpid for the test_p ...
                if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                    // the following function takes raw points ...
                    auto temp_p_raw = transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                    if (m_grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2*units::cm, 0, 0)) {
                        ps_list.push_back(test_p);
                        break;
                    }
                }
            }
        } else {
            ps_list.push_back(temp_end);
        }
    }
    
    std::vector<WireCell::Point> tmp_pts(ps_list.begin(), ps_list.end());
    return tmp_pts;
}


void TrackFitting::organize_ps_path(std::shared_ptr<PR::Segment> segment, std::vector<WireCell::Point>& pts, double low_dis_limit, double end_point_limit) {
    
    std::vector<WireCell::Point> ps_vec = examine_end_ps_vec(segment, pts, true, true);
    if (ps_vec.size() <= 1) ps_vec = pts;
 
    pts.clear();
    // fill in the beginning part
    {
        WireCell::Point p1 = ps_vec.front();
        WireCell::Point p2 = ps_vec.front();
        double dis1 = 0;
        for (auto it = ps_vec.begin(); it != ps_vec.end(); it++) {
            p2 = *it;
            dis1 = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
            if (dis1 > low_dis_limit) break;
        }
        if (dis1 > low_dis_limit) {
            WireCell::Point extended_p1(
                p1.x() + (p1.x() - p2.x()) / dis1 * end_point_limit,
                p1.y() + (p1.y() - p2.y()) / dis1 * end_point_limit,
                p1.z() + (p1.z() - p2.z()) / dis1 * end_point_limit
            );
            pts.push_back(extended_p1);
        }
    }
    
    // fill in the middle part
    for (size_t i = 0; i != ps_vec.size(); i++) {
        WireCell::Point p1 = ps_vec.at(i);
        double dis;
        if (pts.size() != 0) {
            dis = sqrt(pow(p1.x() - pts.back().x(), 2) + pow(p1.y() - pts.back().y(), 2) + pow(p1.z() - pts.back().z(), 2));
        } else {
            dis = sqrt(pow(p1.x() - ps_vec.back().x(), 2) + pow(p1.y() - ps_vec.back().y(), 2) + pow(p1.z() - ps_vec.back().z(), 2));
        }
        
        if (dis < low_dis_limit * 0.8) {
            continue;
        } else if (dis < low_dis_limit * 1.6) {
            pts.push_back(p1);
        } else {
            int npoints = std::round(dis / low_dis_limit);
            WireCell::Point p_save = pts.back();
            for (int j = 0; j != npoints; j++) {
                WireCell::Point p(
                    p_save.x() + (p1.x() - p_save.x()) / npoints * (j + 1),
                    p_save.y() + (p1.y() - p_save.y()) / npoints * (j + 1),
                    p_save.z() + (p1.z() - p_save.z()) / npoints * (j + 1)
                );
                pts.push_back(p);
            }
        }
    }
    
    // fill in the end part
    if (end_point_limit != 0) {
        WireCell::Point p1 = ps_vec.back();
        WireCell::Point p2 = ps_vec.back();
        double dis1 = 0;
        for (auto it = ps_vec.rbegin(); it != ps_vec.rend(); it++) {
            p2 = *it;
            dis1 = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
            if (dis1 > low_dis_limit) break;
        }
        if (dis1 != 0) {
            WireCell::Point extended_p1(
                p1.x() + (p1.x() - p2.x()) / dis1 * end_point_limit,
                p1.y() + (p1.y() - p2.y()) / dis1 * end_point_limit,
                p1.z() + (p1.z() - p2.z()) / dis1 * end_point_limit
            );
            pts.push_back(extended_p1);
        }
    } else {
        WireCell::Point p1 = ps_vec.back();
        double dis1 = sqrt(pow(p1.x() - pts.back().x(), 2) + pow(p1.y() - pts.back().y(), 2) + pow(p1.z() - pts.back().z(), 2));
        if (dis1 >= 0.45*units::cm)
            pts.push_back(p1);
    }
    
    if (pts.size() <= 1)
        pts = ps_vec;
}


 void TrackFitting::form_point_association(std::shared_ptr<PR::Segment> segment,WireCell::Point &p, PlaneData& temp_2dut, PlaneData& temp_2dvt, PlaneData& temp_2dwt, double dis_cut, int nlevel, double time_tick_cut ){

     // Clear previous associations
    temp_2dut.associated_2d_points.clear();
    temp_2dvt.associated_2d_points.clear();
    temp_2dwt.associated_2d_points.clear();
    
    // Get cluster from segment
    auto cluster = segment->cluster();
    const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
    double cluster_t0 = cluster->get_flash().time();
    // find the raw point ...


    // Get closest point in cluster and find neighbors using graph
    auto closest_result = cluster->get_closest_wcpoint(geo_point_t(p.x(), p.y(), p.z()));
    size_t closest_point_index = closest_result.first;
    geo_point_t closest_point = closest_result.second;

    double temp_dis = sqrt(pow(closest_point.x() - p.x(), 2) + 
                        pow(closest_point.y() - p.y(), 2) + 
                        pow(closest_point.z() - p.z(), 2));
    // find the WPID for this point ...
    WirePlaneId wpid = cluster->wire_plane_id(closest_point_index);
    int apa = wpid.apa();
    int face = wpid.face();

    // Get geometry information from TrackFitting internal data
    auto paras = wpid_params.find(wpid);
    auto geoms = wpid_geoms.find(wpid);

    double angle_u = std::get<1>(paras->second);
    double angle_v = std::get<2>(paras->second);
    double angle_w = std::get<3>(paras->second);

    double time_tick_width = std::get<0>(geoms->second);
    double pitch_u = std::get<1>(geoms->second);
    double pitch_v = std::get<2>(geoms->second);
    double pitch_w = std::get<3>(geoms->second);

     // Get wire indices and time slice for current point
    WirePlaneId current_wpid = wpid;
    int cur_wire_u = cluster->wire_index(closest_point_index, 0);
    int cur_wire_v = cluster->wire_index(closest_point_index, 1);
    int cur_wire_w = cluster->wire_index(closest_point_index, 2);
    int cur_time_slice =  cluster->blob_with_point(closest_point_index)->slice_index_min();
    int cur_ntime_ticks = cluster->blob_with_point(closest_point_index)->slice_index_max() - cur_time_slice;

    if (temp_dis < dis_cut){
        
        // auto p_raw = transform->backward(p, cluster_t0, apa, face);
        // std::cout << "WirePlaneId: " << wpid << ", Angles: (" << angle_u << ", " << angle_v << ", " << angle_w << ")" << " " << time_tick_width/units::cm << " " << pitch_u/units::cm << " " << pitch_v/units::cm << " " << pitch_w/units::cm << std::endl;
        // Get graph algorithms interface
        // auto cached_gas = cluster->get_cached_graph_algorithms();
        // for (auto ga: cached_gas){
        //     std::cout << "GraphAlgorithm name: " << ga << std::endl;
        // }
        
        const auto& ga = cluster->graph_algorithms("basic_pid");
        //Find nearby points using graph traversal
        auto total_vertices_found = ga.find_neighbors_nlevel(closest_point_index, nlevel);
        // std::cout << "Neighbors: " << closest_point_index << " " << total_vertices_found.size() << std::endl;
        
        // Collect nearby blobs and their properties
        std::set<const Facade::Blob*> nearby_blobs_set;
        for (auto vertex_idx : total_vertices_found) {
            const Facade::Blob* blob = cluster->blob_with_point(vertex_idx);
            if (blob) {
                auto blob_wpid = blob->wpid();
                if (blob_wpid.apa()==apa && blob_wpid.face()==face)
                    nearby_blobs_set.insert(blob);
            }
            // // print out the distance between the vertex_idx and the original point 
            // geo_point_t vertex_point = cluster->point3d(vertex_idx);
            // double distance = sqrt(pow(vertex_point.x() - p.x(), 2) + 
            //                       pow(vertex_point.y() - p.y(), 2) + 
            //                       pow(vertex_point.z() - p.z(), 2));
            // std::cout << "Vertex " << vertex_idx << " distance to original point: " 
            //           << distance/units::cm << " cm" << std::endl;
        }
        
       
        
        // std::cout << "Cur: " << cur_time_slice << " " << cur_wire_u << " " << cur_wire_v << " " << cur_wire_w << std::endl;
        
        // Calculate adaptive distance cuts for each plane
        double dis_cut_u = dis_cut;
        double dis_cut_v = dis_cut;
        double dis_cut_w = dis_cut;
        double max_time_slice_u = 0, max_time_slice_v = 0, max_time_slice_w = 0;
        
        // Find maximum time slice differences for adaptive cuts
        for (const auto* blob : nearby_blobs_set) {
            int this_time_slice = blob->slice_index_min();
            
            // Check U plane
            if (cur_wire_u >= blob->u_wire_index_min()-1 && cur_wire_u < blob->u_wire_index_max() + 1) {
                max_time_slice_u = std::max(max_time_slice_u, static_cast<double>(abs(this_time_slice - cur_time_slice)));
            }
            
            // Check V plane
            if (cur_wire_v >= blob->v_wire_index_min()-1 && cur_wire_v < blob->v_wire_index_max() + 1) {
                max_time_slice_v = std::max(max_time_slice_v, static_cast<double>(abs(this_time_slice - cur_time_slice)));
            }
            
            // Check W plane
            if (cur_wire_w >= blob->w_wire_index_min()-1 && cur_wire_w < blob->w_wire_index_max() + 1) {
                max_time_slice_w = std::max(max_time_slice_w, static_cast<double>(abs(this_time_slice - cur_time_slice)));
            }
        }
        
        // Update distance cuts based on time slice spans
        if (max_time_slice_u * time_tick_width * 1.2 < dis_cut_u)
            dis_cut_u = max_time_slice_u * time_tick_width * 1.2;
        if (max_time_slice_v * time_tick_width * 1.2 < dis_cut_v)
            dis_cut_v = max_time_slice_v * time_tick_width * 1.2;
        if (max_time_slice_w * time_tick_width * 1.2 < dis_cut_w)
            dis_cut_w = max_time_slice_w * time_tick_width * 1.2;
        
        // Process each nearby blob for wire range calculations
        for (const auto* blob : nearby_blobs_set) {
            int this_time_slice = blob->slice_index_min();
            
            // Calculate remaining distance cuts accounting for time offset
            double rem_dis_cut_u = pow(dis_cut_u, 2) - pow((cur_time_slice - this_time_slice) * time_tick_width, 2);
            double rem_dis_cut_v = pow(dis_cut_v, 2) - pow((cur_time_slice - this_time_slice) * time_tick_width, 2);
            double rem_dis_cut_w = pow(dis_cut_w, 2) - pow((cur_time_slice - this_time_slice) * time_tick_width, 2);
            
            if ((rem_dis_cut_u > 0 || rem_dis_cut_v > 0 || rem_dis_cut_w > 0) && abs(cur_time_slice - this_time_slice) <= time_tick_cut) {
                
                // Calculate minimum wire distances
                float min_u_dis, min_v_dis, min_w_dis;
                
                // U wire distance
                if (cur_wire_u < blob->u_wire_index_min()) {
                    min_u_dis = blob->u_wire_index_min() - cur_wire_u;
                } else if (cur_wire_u <= blob->u_wire_index_max()) {
                    min_u_dis = 0;
                } else {
                    min_u_dis = cur_wire_u - blob->u_wire_index_max();
                }
                
                // V wire distance
                if (cur_wire_v < blob->v_wire_index_min()) {
                    min_v_dis = blob->v_wire_index_min() - cur_wire_v;
                } else if (cur_wire_v <= blob->v_wire_index_max()) {
                    min_v_dis = 0;
                } else {
                    min_v_dis = cur_wire_v - blob->v_wire_index_max();
                }
                
                // W wire distance
                if (cur_wire_w < blob->w_wire_index_min()) {
                    min_w_dis = blob->w_wire_index_min() - cur_wire_w;
                } else if (cur_wire_w <= blob->w_wire_index_max()) {
                    min_w_dis = 0;
                } else {
                    min_w_dis = cur_wire_w - blob->w_wire_index_max();
                }
                
                // Use the dedicated calculate_ranges_simplified function
                float range_u, range_v, range_w;
                WireCell::Clus::TrackFittingUtil::calculate_ranges_simplified(
                    angle_u, angle_v, angle_w,
                    rem_dis_cut_u, rem_dis_cut_v, rem_dis_cut_w,
                    min_u_dis, min_v_dis, min_w_dis,
                    pitch_u, pitch_v, pitch_w,
                    range_u, range_v, range_w);
                
                // If all ranges are positive, add wire indices to associations
                if (range_u > 0 && range_v > 0 && range_w > 0) {
                    // Calculate wire limits
                    float low_u_limit = cur_wire_u - sqrt(range_u) / pitch_u;
                    float high_u_limit = cur_wire_u + sqrt(range_u) / pitch_u;
                    float low_v_limit = cur_wire_v - sqrt(range_v) / pitch_v;
                    float high_v_limit = cur_wire_v + sqrt(range_v) / pitch_v;
                    float low_w_limit = cur_wire_w - sqrt(range_w) / pitch_w;
                    float high_w_limit = cur_wire_w + sqrt(range_w) / pitch_w;
                    
                    // Add U plane associations
                    for (int j = std::round(low_u_limit); j <= std::round(high_u_limit); j++) {
                        Coord2D coord(current_wpid.apa(), current_wpid.face(), 
                                    this_time_slice, j, 
                                    get_channel_for_wire(current_wpid.apa(), current_wpid.face(), 0, j), 
                                    WirePlaneLayer_t::kUlayer);
                        temp_2dut.associated_2d_points.insert(coord);
                    }
                    
                    // Add V plane associations
                    for (int j = std::round(low_v_limit); j <= std::round(high_v_limit); j++) {
                        Coord2D coord(current_wpid.apa(), current_wpid.face(), 
                                    this_time_slice, j, 
                                    get_channel_for_wire(current_wpid.apa(), current_wpid.face(), 1, j), 
                                    WirePlaneLayer_t::kVlayer);
                        temp_2dvt.associated_2d_points.insert(coord);
                    }
                    
                    // Add W plane associations
                    for (int j = std::round(low_w_limit); j <= std::round(high_w_limit); j++) {
                        Coord2D coord(current_wpid.apa(), current_wpid.face(), 
                                    this_time_slice, j, 
                                    get_channel_for_wire(current_wpid.apa(), current_wpid.face(), 2, j), 
                                    WirePlaneLayer_t::kWlayer);
                        temp_2dwt.associated_2d_points.insert(coord);
                    }
                }
            }
        }
    }

    // std::cout << "Pixels: " << temp_2dut.associated_2d_points.size() << " " << temp_2dvt.associated_2d_points.size() << " " << temp_2dwt.associated_2d_points.size() << std::endl;

    // Steiner Tree ... 
    if (cluster->has_graph("steiner_graph") && cluster->has_pc("steiner_pc")) {
        auto graph_name = "steiner_graph";
        auto pc_name = "steiner_pc";   
        const auto& steiner_pc = cluster->get_pc(pc_name);
        const auto& coords = cluster->get_default_scope().coords;
        const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
        const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
        const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
        const auto& wpid_array = steiner_pc.get("wpid")->elements<WirePlaneId>();


        auto steiner_search_result = cluster->kd_steiner_knn(1, p);
        auto steiner_search_point = cluster->kd_steiner_points(steiner_search_result);

        size_t closest_point_index = steiner_search_result.front().first;
        auto closest_point = steiner_search_point.front().first;
        auto closest_point_wpid = steiner_search_point.front().second.first;

        // 
        // std::cout << closest_point_index << " " << closest_point << " " << p << " " << closest_point_wpid << " " << test << std::endl;
        
        double temp_dis = sqrt(pow(closest_point.x() - p.x(), 2) + 
                               pow(closest_point.y() - p.y(), 2) + 
                               pow(closest_point.z() - p.z(), 2));

        if (temp_dis < dis_cut && apa == closest_point_wpid.apa() && face == closest_point_wpid.face()){
            // Get graph algorithms interface
            const auto& ga = cluster->graph_algorithms(graph_name);
            // Find nearby points using graph traversal (equivalent to original nested loop)
            auto total_vertices_found = ga.find_neighbors_nlevel(closest_point_index, nlevel);

            // find the raw point ...
            auto closest_point_raw = transform->backward(closest_point, cluster_t0, apa, face);


            auto cur_u = m_grouping->convert_3Dpoint_time_ch(closest_point_raw, apa, face, 0);
            auto cur_v = m_grouping->convert_3Dpoint_time_ch(closest_point_raw, apa, face, 1);
            auto cur_w = m_grouping->convert_3Dpoint_time_ch(closest_point_raw, apa, face, 2);

            int cur_time_slice = std::floor(std::get<0>(cur_u)/cur_ntime_ticks)*cur_ntime_ticks;
            int cur_wire_u = std::get<1>(cur_u);
            int cur_wire_v = std::get<1>(cur_v);
            int cur_wire_w = std::get<1>(cur_w);

            // std::cout << cur_time_slice << " " << cur_wire_u << " " << cur_wire_v << " " << cur_wire_w << std::endl;
        
            // Calculate adaptive distance cuts (equivalent to original max_time_slice_u/v/w calculation)
            double dis_cut_u = dis_cut;
            double dis_cut_v = dis_cut;
            double dis_cut_w = dis_cut;
            
            double max_time_slice_u = 0;
            double max_time_slice_v = 0;
            double max_time_slice_w = 0;

            std::map<int, std::tuple<int, int, int, int>> map_vertex_info;
            // Collect point indices 
            for (auto vertex_idx : total_vertices_found) {
                    auto vertex_wpid = wpid_array[vertex_idx];
                    // std::cout << vertex_wpid << " " << std::endl;

                    if (vertex_wpid.apa() != apa || vertex_wpid.face() != face) continue;

                    // Handle points not associated with blobs
                    geo_point_t vertex_point = {x_coords[vertex_idx], 
                                                y_coords[vertex_idx], 
                                                z_coords[vertex_idx]};

                    auto vertex_point_raw = transform->backward(vertex_point, cluster_t0, apa, face);
                    auto vertex_u = m_grouping->convert_3Dpoint_time_ch(vertex_point_raw, apa, face, 0);
                    auto vertex_v = m_grouping->convert_3Dpoint_time_ch(vertex_point_raw, apa, face, 1);
                    auto vertex_w = m_grouping->convert_3Dpoint_time_ch(vertex_point_raw, apa, face, 2);

                    int vertex_time_slice = std::floor(std::get<0>(vertex_u)/cur_ntime_ticks)*cur_ntime_ticks;
                    int vertex_wire_u = std::get<1>(vertex_u);
                    int vertex_wire_v = std::get<1>(vertex_v);
                    int vertex_wire_w = std::get<1>(vertex_w);

                    map_vertex_info[vertex_idx] = std::make_tuple(vertex_time_slice, vertex_wire_u, vertex_wire_v, vertex_wire_w);

                    // Check U, V, W wire proximity
                    if (abs(vertex_wire_u - cur_wire_u) <= 1) {
                        if (abs(vertex_time_slice - cur_time_slice) > max_time_slice_u)
                            max_time_slice_u = abs(vertex_time_slice - cur_time_slice);
                    }
                    if (abs(vertex_wire_v - cur_wire_v) <= 1) {
                        if (abs(vertex_time_slice - cur_time_slice) > max_time_slice_v)
                            max_time_slice_v = abs(vertex_time_slice - cur_time_slice);
                    }
                    if (abs(vertex_wire_w - cur_wire_w) <= 1) {
                        if (abs(vertex_time_slice - cur_time_slice) > max_time_slice_w)
                            max_time_slice_w = abs(vertex_time_slice - cur_time_slice);
                    }
            }

            // Apply adaptive cuts (equivalent to original adaptive cut logic)
            if (max_time_slice_u * time_tick_width * 1.2 < dis_cut_u)
                dis_cut_u = max_time_slice_u * time_tick_width * 1.2;
            if (max_time_slice_v * time_tick_width * 1.2 < dis_cut_v)
                dis_cut_v = max_time_slice_v * time_tick_width * 1.2;
            if (max_time_slice_w * time_tick_width * 1.2 < dis_cut_w)
                dis_cut_w = max_time_slice_w * time_tick_width * 1.2;
        
            // // Process each nearby blob for wire range calculations (equivalent to final loop)
            for (auto vertex_info : map_vertex_info){
                // auto vertex_idx = vertex_info.first;
                auto vertex_time_slice = std::get<0>(vertex_info.second);
                auto vertex_wire_u = std::get<1>(vertex_info.second);
                auto vertex_wire_v = std::get<2>(vertex_info.second);
                auto vertex_wire_w = std::get<3>(vertex_info.second);

                // Calculate remaining distance cuts accounting for time offset
                double rem_dis_cut_u = pow(dis_cut_u, 2) - pow((cur_time_slice - vertex_time_slice) * time_tick_width, 2);
                double rem_dis_cut_v = pow(dis_cut_v, 2) - pow((cur_time_slice - vertex_time_slice) * time_tick_width, 2);
                double rem_dis_cut_w = pow(dis_cut_w, 2) - pow((cur_time_slice - vertex_time_slice) * time_tick_width, 2);

                if ((rem_dis_cut_u > 0 || rem_dis_cut_v > 0 || rem_dis_cut_w > 0) && abs(cur_time_slice - vertex_time_slice) <= time_tick_cut) {
                
                // Calculate minimum wire distances (equivalent to original min_u/v/w_dis calculation)
                float min_u_dis, min_v_dis, min_w_dis;
                
                // U wire distance
                if (cur_wire_u < vertex_wire_u - 1) {
                    min_u_dis = vertex_wire_u - 1 - cur_wire_u;
                } else if (cur_wire_u > vertex_wire_u + 1) {
                    min_u_dis = cur_wire_u - vertex_wire_u - 1;
                } else {
                    min_u_dis = 0;
                }

                // V wire distance
                if (cur_wire_v < vertex_wire_v - 1) {
                    min_v_dis = vertex_wire_v - 1 - cur_wire_v;
                } else if (cur_wire_v > vertex_wire_v + 1) {
                    min_v_dis = cur_wire_v - vertex_wire_v - 1;
                } else {
                    min_v_dis = 0;
                }

                // W wire distance
                if (cur_wire_w < vertex_wire_w - 1) {
                    min_w_dis = vertex_wire_w - 1 - cur_wire_w;
                } else if (cur_wire_w > vertex_wire_w + 1) {
                    min_w_dis = cur_wire_w - vertex_wire_w - 1;
                } else {
                    min_w_dis = 0;
                }

                // Use the dedicated calculate_ranges_simplified function (replaces original range calculation)
                float range_u, range_v, range_w;
                WireCell::Clus::TrackFittingUtil::calculate_ranges_simplified(
                    angle_u, angle_v, angle_w,
                    rem_dis_cut_u, rem_dis_cut_v, rem_dis_cut_w,
                    min_u_dis, min_v_dis, min_w_dis,
                    pitch_u, pitch_v, pitch_w,
                    range_u, range_v, range_w);
                
                    // If all ranges are positive, add wire indices to associations
                    if (range_u > 0 && range_v > 0 && range_w > 0) {
                        // Calculate wire limits (equivalent to original low/high_limit calculations)
                        float low_u_limit = cur_wire_u - sqrt(range_u) / pitch_u;
                        float high_u_limit = cur_wire_u + sqrt(range_u) / pitch_u;
                        float low_v_limit = cur_wire_v - sqrt(range_v) / pitch_v;
                        float high_v_limit = cur_wire_v + sqrt(range_v) / pitch_v;
                        float low_w_limit = cur_wire_w - sqrt(range_w) / pitch_w;
                        float high_w_limit = cur_wire_w + sqrt(range_w) / pitch_w;
                        
                        // Add U plane associations (equivalent to temp_2dut.insert)
                        for (int j = std::round(low_u_limit); j <= std::round(high_u_limit); j++) {
                            Coord2D coord(apa, face, 
                                        vertex_time_slice, j, 
                                        get_channel_for_wire(apa, face, 0, j), 
                                        WirePlaneLayer_t::kUlayer);
                            temp_2dut.associated_2d_points.insert(coord);
                        }
                        // Add V
                        for (int j = std::round(low_v_limit); j <= std::round(high_v_limit); j++) {
                            Coord2D coord(apa, face, 
                                        vertex_time_slice, j, 
                                        get_channel_for_wire(apa, face, 1, j), 
                                        WirePlaneLayer_t::kVlayer);
                            temp_2dvt.associated_2d_points.insert(coord);
                        }

                        // Add W
                        for (int j = std::round(low_w_limit); j <= std::round(high_w_limit); j++) {
                            Coord2D coord(apa, face, 
                                        vertex_time_slice, j, 
                                        get_channel_for_wire(apa, face, 2, j), 
                                        WirePlaneLayer_t::kWlayer);
                            temp_2dwt.associated_2d_points.insert(coord);
                        }
                    }
                }
            }
        }
    }

    
    // std::cout << "Pixels 1: " << temp_2dut.associated_2d_points.size() << " " << temp_2dvt.associated_2d_points.size() << " " << temp_2dwt.associated_2d_points.size() << std::endl;

     // Fallback: simple projection if no associations were found from complex method
    if (temp_2dut.associated_2d_points.size() == 0 && 
        temp_2dvt.associated_2d_points.size() == 0 && 
        temp_2dwt.associated_2d_points.size() == 0) {
               
        // Convert time_tick_cut to integer for iteration
        int time_cut = static_cast<int>(time_tick_cut);
        int wire_cut = time_cut/cur_ntime_ticks;
        
        // Simple diamond pattern projection around current point
        for (int i = -wire_cut; i <= wire_cut; i++) {
            // loop over time dimension ...
            for (int j = -time_cut; j <= time_cut; j+=cur_ntime_ticks) {
                if (abs(i*cur_ntime_ticks) + abs(j) <= time_cut) {
                    // U plane projection
                    Coord2D coord_u(apa, face, cur_time_slice + j, cur_wire_u + i, 
                                   get_channel_for_wire(apa, face, 0, cur_wire_u + i), 
                                   kUlayer);
                    temp_2dut.associated_2d_points.insert(coord_u);
                    
                    // V plane projection  
                    Coord2D coord_v(apa, face, cur_time_slice + j, cur_wire_v + i, 
                                   get_channel_for_wire(apa, face, 1, cur_wire_v + i), 
                                   kVlayer);
                    temp_2dvt.associated_2d_points.insert(coord_v);
                    
                    // W plane projection
                    Coord2D coord_w(apa, face, cur_time_slice + j, cur_wire_w + i, 
                                   get_channel_for_wire(apa, face, 2, cur_wire_w + i), 
                                   kWlayer);
                    temp_2dwt.associated_2d_points.insert(coord_w);
                }
            }
        }
    }

 }


 void TrackFitting::examine_point_association(std::shared_ptr<PR::Segment> segment, WireCell::Point &p, PlaneData& temp_2dut, PlaneData& temp_2dvt, PlaneData& temp_2dwt, bool flag_end_point, double charge_cut){

    // Get cluster from segment
    auto cluster = segment->cluster();
    const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
    double cluster_t0 = cluster->get_flash().time();

    // find the apa and face ...
    auto wpid = m_dv->contained_by(p);
    int apa = wpid.apa();
    int face = wpid.face();
    
    if (apa == -1 || face == -1) return;

    // // Convert 3D point to wire/time coordinates for each plane
    geo_point_t p_raw = transform->backward(p, cluster_t0, apa, face);
    auto [time_u, wire_u] = m_grouping->convert_3Dpoint_time_ch(p_raw, apa, face, 0);
    auto [time_v, wire_v] = m_grouping->convert_3Dpoint_time_ch(p_raw, apa, face, 1);
    auto [time_w, wire_w] = m_grouping->convert_3Dpoint_time_ch(p_raw, apa, face, 2);

    std::set<int> temp_types_u;
    std::set<int> temp_types_v;
    std::set<int> temp_types_w;
    
    std::set<Coord2D> saved_2dut;
    std::set<Coord2D> saved_2dvt;
    std::set<Coord2D> saved_2dwt;

    std::vector<float> results;
    results.resize(3,0);
    
    // Process U plane
    for (auto it = temp_2dut.associated_2d_points.begin(); it != temp_2dut.associated_2d_points.end(); it++){
        CoordReadout coord_key(it->apa, it->time, it->channel);
        auto charge_it = m_charge_data.find(coord_key);
        if (charge_it != m_charge_data.end() && charge_it->second.charge > charge_cut) {
            temp_types_u.insert(charge_it->second.flag);
            if (charge_it->second.flag == 0) results.at(0)++;
            saved_2dut.insert(*it);
        }
    }

    // Process V plane
    for (auto it = temp_2dvt.associated_2d_points.begin(); it != temp_2dvt.associated_2d_points.end(); it++){
        CoordReadout coord_key(it->apa, it->time, it->channel);
        auto charge_it = m_charge_data.find(coord_key);
        if (charge_it != m_charge_data.end() && charge_it->second.charge > charge_cut) {
            temp_types_v.insert(charge_it->second.flag);
            if (charge_it->second.flag == 0) results.at(1)++;
            saved_2dvt.insert(*it);
        }
    }

    // Process W plane
    for (auto it = temp_2dwt.associated_2d_points.begin(); it != temp_2dwt.associated_2d_points.end(); it++){
        CoordReadout coord_key(it->apa, it->time, it->channel);
        auto charge_it = m_charge_data.find(coord_key);
        if (charge_it != m_charge_data.end() && charge_it->second.charge > charge_cut) {
            temp_types_w.insert(charge_it->second.flag);
            if (charge_it->second.flag == 0) results.at(2)++;
            saved_2dwt.insert(*it);
        }
    }

    // Calculate quality ratios
    if (temp_2dut.associated_2d_points.size() != 0)
        results.at(0) = (saved_2dut.size() - results.at(0)*1.0)/temp_2dut.associated_2d_points.size();
    else
        results.at(0) = 0;
    
    if (temp_2dvt.associated_2d_points.size() != 0)
        results.at(1) = (saved_2dvt.size() - results.at(1)*1.0)/temp_2dvt.associated_2d_points.size();
    else
        results.at(1) = 0;
    
    if (temp_2dwt.associated_2d_points.size() != 0)
        results.at(2) = (saved_2dwt.size() - results.at(2)*1.0)/temp_2dwt.associated_2d_points.size();
    else
        results.at(2) = 0;

    // Reset if only flag 0 found
    if (temp_types_u.find(0) != temp_types_u.end() && temp_types_u.size() == 1){
        saved_2dut.clear();
        results.at(0) = 0;
    }
    if (temp_types_v.find(0) != temp_types_v.end() && temp_types_v.size() == 1){
        saved_2dvt.clear();
        results.at(1) = 0;
    }
    if (temp_types_w.find(0) != temp_types_w.end() && temp_types_w.size() == 1){
        saved_2dwt.clear();
        results.at(2) = 0;
    }

    // Handle dead plane scenarios
    // U and V planes are dead ...
    if (saved_2dut.size() == 0 && saved_2dvt.size() == 0 && saved_2dwt.size() != 0){
        int channel_u = get_channel_for_wire(apa, face, 0, wire_u);
        int channel_v = get_channel_for_wire(apa, face, 1, wire_v);
        saved_2dut.insert(Coord2D(apa, face, time_u, wire_u, channel_u, kUlayer));
        saved_2dvt.insert(Coord2D(apa, face, time_v, wire_v, channel_v, kVlayer));
        
        // W plane check for outliers
        if (!flag_end_point && saved_2dwt.size() > 0)
        {
            std::pair<double, double> ave_pos = std::make_pair(0,0);
            double total_charge = 0;
            for (auto it1 = saved_2dwt.begin(); it1 != saved_2dwt.end(); it1++){
                CoordReadout coord_key(it1->apa, it1->time, it1->channel);
                auto charge_it = m_charge_data.find(coord_key);
                if (charge_it != m_charge_data.end()){
                    ave_pos.first += it1->wire * charge_it->second.charge;
                    ave_pos.second += it1->time * charge_it->second.charge;
                    total_charge += charge_it->second.charge;
                }
            }
            if (total_charge != 0){
                ave_pos.first /= total_charge;
                ave_pos.second /= total_charge;
            }
            double rms = 0;
            for (auto it1 = saved_2dwt.begin(); it1 != saved_2dwt.end(); it1++){
                rms += pow(it1->wire - ave_pos.first, 2) + pow(it1->time - ave_pos.second, 2);
            }
            rms = sqrt(rms/saved_2dwt.size());

            if (sqrt(pow(ave_pos.first - wire_w, 2) + pow(ave_pos.second - time_w, 2)) > 0.75*rms && 
                saved_2dwt.size() <= 5 && saved_2dwt.size() < 0.2 * temp_2dwt.associated_2d_points.size()){
                saved_2dwt.clear();
                int channel_w = get_channel_for_wire(apa, face, 2, wire_w);
                saved_2dwt.insert(Coord2D(apa, face, time_w, wire_w, channel_w, kWlayer));
                results.at(2) = 0;
            }
        }
    }
    else if (saved_2dut.size() == 0 && saved_2dwt.size() == 0 && saved_2dvt.size() != 0){
        // U and W planes are dead ...
        int channel_u = get_channel_for_wire(apa, face, 0, wire_u);
        int channel_w = get_channel_for_wire(apa, face, 2, wire_w);
        saved_2dut.insert(Coord2D(apa, face, time_u, wire_u, channel_u, kUlayer));
        saved_2dwt.insert(Coord2D(apa, face, time_w, wire_w, channel_w, kWlayer));
        
        // V plane check for outliers
        if (!flag_end_point && saved_2dvt.size() > 0)
        {
            std::pair<double, double> ave_pos = std::make_pair(0,0);
            double total_charge = 0;
            for (auto it1 = saved_2dvt.begin(); it1 != saved_2dvt.end(); it1++){
                CoordReadout coord_key(it1->apa, it1->time, it1->channel);
                auto charge_it = m_charge_data.find(coord_key);
                if (charge_it != m_charge_data.end()){
                    ave_pos.first += it1->wire * charge_it->second.charge;
                    ave_pos.second += it1->time * charge_it->second.charge;
                    total_charge += charge_it->second.charge;
                }
            }
            if (total_charge != 0){
                ave_pos.first /= total_charge;
                ave_pos.second /= total_charge;
            }
            double rms = 0;
            for (auto it1 = saved_2dvt.begin(); it1 != saved_2dvt.end(); it1++){
                rms += pow(it1->wire - ave_pos.first, 2) + pow(it1->time - ave_pos.second, 2);
            }
            rms = sqrt(rms/saved_2dvt.size());

            if (sqrt(pow(ave_pos.first - wire_v, 2) + pow(ave_pos.second - time_v, 2)) > 0.75*rms && 
                saved_2dvt.size() <= 5 && saved_2dvt.size() < 0.2 * temp_2dvt.associated_2d_points.size()){
                saved_2dvt.clear();
                int channel_v = get_channel_for_wire(apa, face, 1, wire_v);
                saved_2dvt.insert(Coord2D(apa, face, time_v, wire_v, channel_v, kVlayer));
                results.at(1) = 0;
            }
        }
    }
    else if (saved_2dvt.size() == 0 && saved_2dwt.size() == 0 && saved_2dut.size() != 0){
        // V and W planes are dead ...
        int channel_v = get_channel_for_wire(apa, face, 1, wire_v);
        int channel_w = get_channel_for_wire(apa, face, 2, wire_w);
        saved_2dvt.insert(Coord2D(apa, face, time_v, wire_v, channel_v, kVlayer));
        saved_2dwt.insert(Coord2D(apa, face, time_w, wire_w, channel_w, kWlayer));
        
        // U plane check for outliers
        if (!flag_end_point && saved_2dut.size() > 0)
        {
            std::pair<double, double> ave_pos = std::make_pair(0,0);
            double total_charge = 0;
            for (auto it1 = saved_2dut.begin(); it1 != saved_2dut.end(); it1++){
                CoordReadout coord_key(it1->apa, it1->time, it1->channel);
                auto charge_it = m_charge_data.find(coord_key);
                if (charge_it != m_charge_data.end()){
                    ave_pos.first += it1->wire * charge_it->second.charge;
                    ave_pos.second += it1->time * charge_it->second.charge;
                    total_charge += charge_it->second.charge;
                }
            }
            if (total_charge != 0){
                ave_pos.first /= total_charge;
                ave_pos.second /= total_charge;
            }
            double rms = 0;
            for (auto it1 = saved_2dut.begin(); it1 != saved_2dut.end(); it1++){
                rms += pow(it1->wire - ave_pos.first, 2) + pow(it1->time - ave_pos.second, 2);
            }
            rms = sqrt(rms/saved_2dut.size());

            if (sqrt(pow(ave_pos.first - wire_u, 2) + pow(ave_pos.second - time_u, 2)) > 0.75*rms && 
                saved_2dut.size() <= 5 && saved_2dut.size() < 0.2 * temp_2dut.associated_2d_points.size()){
                saved_2dut.clear();
                int channel_u = get_channel_for_wire(apa, face, 0, wire_u);
                saved_2dut.insert(Coord2D(apa, face, time_u, wire_u, channel_u, kUlayer));
                results.at(0) = 0;
            }
        }
    }
    // Handle partial dead plane scenarios (only one plane dead, check outliers in others)
    else if (saved_2dut.size() == 0 && saved_2dwt.size() != 0 && saved_2dvt.size() != 0){
        // Only U plane is dead, check W and V plane outliers
        auto check_outliers = [&](std::set<Coord2D>& saved_plane, std::vector<float>& results, int result_idx, 
                                 const std::set<Coord2D>& temp_plane, int expected_wire, int expected_time) {
            if (!flag_end_point && saved_plane.size() > 0)
            {
                std::pair<double, double> ave_pos = std::make_pair(0,0);
                double total_charge = 0;
                for (auto it1 = saved_plane.begin(); it1 != saved_plane.end(); it1++){
                    CoordReadout coord_key(it1->apa, it1->time, it1->channel);
                    auto charge_it = m_charge_data.find(coord_key);
                    if (charge_it != m_charge_data.end()){
                        ave_pos.first += it1->wire * charge_it->second.charge;
                        ave_pos.second += it1->time * charge_it->second.charge;
                        total_charge += charge_it->second.charge;
                    }
                }
                if (total_charge != 0){
                    ave_pos.first /= total_charge;
                    ave_pos.second /= total_charge;
                }
                double rms = 0;
                for (auto it1 = saved_plane.begin(); it1 != saved_plane.end(); it1++){
                    rms += pow(it1->wire - ave_pos.first, 2) + pow(it1->time - ave_pos.second, 2);
                }
                rms = sqrt(rms/saved_plane.size());

                if (sqrt(pow(ave_pos.first - expected_wire, 2) + pow(ave_pos.second - expected_time, 2)) > 0.75*rms && 
                    saved_plane.size() <= 5 && saved_plane.size() < 0.2 * temp_plane.size()){
                    saved_plane.clear();
                    int channel = get_channel_for_wire(apa, face, result_idx == 2 ? 2 : 1, expected_wire);
                    WirePlaneLayer_t plane_layer = (result_idx == 2) ? kWlayer : kVlayer;
                    saved_plane.insert(Coord2D(apa, face, expected_time, expected_wire, channel, plane_layer));
                    results.at(result_idx) = 0;
                }
            }
        };
        
        check_outliers(saved_2dwt, results, 2, temp_2dwt.associated_2d_points, wire_w, time_w);
        check_outliers(saved_2dvt, results, 1, temp_2dvt.associated_2d_points, wire_v, time_v);
    }
    else if (saved_2dvt.size() == 0 && saved_2dut.size() != 0 && saved_2dwt.size() != 0){
        // Only V plane is dead, check U and W plane outliers
        auto check_outliers = [&](std::set<Coord2D>& saved_plane, std::vector<float>& results, int result_idx, 
                                 const std::set<Coord2D>& temp_plane, int expected_wire, int expected_time) {
            if (!flag_end_point && saved_plane.size() > 0)
            {
                std::pair<double, double> ave_pos = std::make_pair(0,0);
                double total_charge = 0;
                for (auto it1 = saved_plane.begin(); it1 != saved_plane.end(); it1++){
                    CoordReadout coord_key(it1->apa, it1->time, it1->channel);
                    auto charge_it = m_charge_data.find(coord_key);
                    if (charge_it != m_charge_data.end()){
                        ave_pos.first += it1->wire * charge_it->second.charge;
                        ave_pos.second += it1->time * charge_it->second.charge;
                        total_charge += charge_it->second.charge;
                    }
                }
                if (total_charge != 0){
                    ave_pos.first /= total_charge;
                    ave_pos.second /= total_charge;
                }
                double rms = 0;
                for (auto it1 = saved_plane.begin(); it1 != saved_plane.end(); it1++){
                    rms += pow(it1->wire - ave_pos.first, 2) + pow(it1->time - ave_pos.second, 2);
                }
                rms = sqrt(rms/saved_plane.size());

                if (sqrt(pow(ave_pos.first - expected_wire, 2) + pow(ave_pos.second - expected_time, 2)) > 0.75*rms && 
                    saved_plane.size() <= 5 && saved_plane.size() < 0.2 * temp_plane.size()){
                    saved_plane.clear();
                    int channel = get_channel_for_wire(apa, face, result_idx == 0 ? 0 : 2, expected_wire);
                    WirePlaneLayer_t plane_layer = (result_idx == 0) ? kUlayer : kWlayer;
                    saved_plane.insert(Coord2D(apa, face, expected_time, expected_wire, channel, plane_layer));
                    results.at(result_idx) = 0;
                }
            }
        };
        
        check_outliers(saved_2dut, results, 0, temp_2dut.associated_2d_points, wire_u, time_u);
        check_outliers(saved_2dwt, results, 2, temp_2dwt.associated_2d_points, wire_w, time_w);
    }
    else if (saved_2dwt.size() == 0 && saved_2dut.size() != 0 && saved_2dvt.size() != 0){
        // Only W plane is dead, check U and V plane outliers  
        auto check_outliers = [&](std::set<Coord2D>& saved_plane, std::vector<float>& results, int result_idx, 
                                 const std::set<Coord2D>& temp_plane, int expected_wire, int expected_time) {
            if (!flag_end_point && saved_plane.size() > 0)
            {
                std::pair<double, double> ave_pos = std::make_pair(0,0);
                double total_charge = 0;
                for (auto it1 = saved_plane.begin(); it1 != saved_plane.end(); it1++){
                    CoordReadout coord_key(it1->apa, it1->time, it1->channel);
                    auto charge_it = m_charge_data.find(coord_key);
                    if (charge_it != m_charge_data.end()){
                        ave_pos.first += it1->wire * charge_it->second.charge;
                        ave_pos.second += it1->time * charge_it->second.charge;
                        total_charge += charge_it->second.charge;
                    }
                }
                if (total_charge != 0){
                    ave_pos.first /= total_charge;
                    ave_pos.second /= total_charge;
                }
                double rms = 0;
                for (auto it1 = saved_plane.begin(); it1 != saved_plane.end(); it1++){
                    rms += pow(it1->wire - ave_pos.first, 2) + pow(it1->time - ave_pos.second, 2);
                }
                rms = sqrt(rms/saved_plane.size());

                if (sqrt(pow(ave_pos.first - expected_wire, 2) + pow(ave_pos.second - expected_time, 2)) > 0.75*rms && 
                    saved_plane.size() <= 5 && saved_plane.size() < 0.2 * temp_plane.size()){
                    saved_plane.clear();
                    int channel = get_channel_for_wire(apa, face, result_idx, expected_wire);
                    WirePlaneLayer_t plane_layer = (result_idx == 0) ? kUlayer : kVlayer;
                    saved_plane.insert(Coord2D(apa, face, expected_time, expected_wire, channel, plane_layer));
                    results.at(result_idx) = 0;
                }
            }
        };
        
        check_outliers(saved_2dut, results, 0, temp_2dut.associated_2d_points, wire_u, time_u);
        check_outliers(saved_2dvt, results, 1, temp_2dvt.associated_2d_points, wire_v, time_v);
    }
    else if (saved_2dwt.size() == 0 && saved_2dut.size() == 0 && saved_2dvt.size() == 0){
        // All planes are dead, use fallback coordinates
        int channel_u = get_channel_for_wire(apa, face, 0, wire_u);
        int channel_v = get_channel_for_wire(apa, face, 1, wire_v);
        int channel_w = get_channel_for_wire(apa, face, 2, wire_w);
        saved_2dut.insert(Coord2D(apa, face, time_u, wire_u, channel_u, kUlayer));
        saved_2dvt.insert(Coord2D(apa, face, time_v, wire_v, channel_v, kVlayer));
        saved_2dwt.insert(Coord2D(apa, face, time_w, wire_w, channel_w, kWlayer));
    }
    
    // Update PlaneData with filtered results
    temp_2dut.associated_2d_points = saved_2dut;
    temp_2dvt.associated_2d_points = saved_2dvt;
    temp_2dwt.associated_2d_points = saved_2dwt;
    
    // Update quantity fields with calculated results
    temp_2dut.quantity = results.at(0);
    temp_2dvt.quantity = results.at(1);
    temp_2dwt.quantity = results.at(2);

 }

void TrackFitting::form_map(std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>>& ptss, double end_point_factor, double mid_point_factor, int nlevel, double time_tick_cut, double charge_cut) {
    // Implementation of form_map function

    m_3d_to_2d.clear();
    m_2d_to_3d.clear();

    std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> saved_pts;
    int count = 0;
    
    // Calculate distances between consecutive points
    std::vector<double> distances;
    for (size_t i = 0; i + 1 != ptss.size(); i++) {
        distances.push_back(sqrt(pow(ptss.at(i+1).first.x() - ptss.at(i).first.x(), 2) +
                               pow(ptss.at(i+1).first.y() - ptss.at(i).first.y(), 2) +
                               pow(ptss.at(i+1).first.z() - ptss.at(i).first.z(), 2)));
    }

    // Loop over the path
    for (size_t i = 0; i != ptss.size(); i++) {
        double dis_cut;
        if (i == 0) {
            dis_cut = std::min(distances.at(i) * end_point_factor, 4/3. * end_point_factor * units::cm);
        } else if (i + 1 == ptss.size()) {
            dis_cut = std::min(distances.back() * end_point_factor, 4/3. * end_point_factor * units::cm);
        } else {
            dis_cut = std::min(std::max(distances.at(i-1) * mid_point_factor, distances.at(i) * mid_point_factor), 
                              4/3. * mid_point_factor * units::cm);
        }

        // check point's apa and face ...
        // find the apa and face ...
        auto wpid = m_dv->contained_by(ptss.at(i).first);
        auto segment = ptss.at(i).second;
        int apa = wpid.apa();
        int face = wpid.face();
        
        if (apa != -1 && face != -1) {

            TrackFitting::PlaneData temp_2dut, temp_2dvt, temp_2dwt;
            form_point_association(segment, ptss.at(i).first, temp_2dut, temp_2dvt, temp_2dwt, dis_cut, nlevel, time_tick_cut);

            if (i == 0 || i == 1 || i + 1 == ptss.size() || i + 2 == ptss.size()) {
                examine_point_association(segment, ptss.at(i).first, temp_2dut, temp_2dvt, temp_2dwt, true, charge_cut);
            } else {
                examine_point_association(segment, ptss.at(i).first, temp_2dut, temp_2dvt, temp_2dwt, false, charge_cut);
            }

            // Fill the mapping data if we have valid associations
            if (temp_2dut.quantity + temp_2dvt.quantity + temp_2dwt.quantity > 0) {
                m_3d_to_2d[count].set_plane_data(WirePlaneLayer_t::kUlayer, temp_2dut);
                m_3d_to_2d[count].set_plane_data(WirePlaneLayer_t::kVlayer, temp_2dvt);
                m_3d_to_2d[count].set_plane_data(WirePlaneLayer_t::kWlayer, temp_2dwt);


                // Fill reverse mapping for U plane
                for (auto it = temp_2dut.associated_2d_points.begin(); it != temp_2dut.associated_2d_points.end(); it++) {
                    if (m_2d_to_3d.find(*it) == m_2d_to_3d.end()) {
                        std::set<int> temp_set;
                        temp_set.insert(count);
                        m_2d_to_3d[*it] = temp_set;
                    } else {
                        m_2d_to_3d[*it].insert(count);
                    }
                }

                for (auto it = temp_2dvt.associated_2d_points.begin(); it != temp_2dvt.associated_2d_points.end(); it++) {
                    if (m_2d_to_3d.find(*it) == m_2d_to_3d.end()) {
                        std::set<int> temp_set;
                        temp_set.insert(count);
                        m_2d_to_3d[*it] = temp_set;
                    } else {
                        m_2d_to_3d[*it].insert(count);
                    }
                }

                for (auto it = temp_2dwt.associated_2d_points.begin(); it != temp_2dwt.associated_2d_points.end(); it++) {
                    if (m_2d_to_3d.find(*it) == m_2d_to_3d.end()) {
                        std::set<int> temp_set;
                        temp_set.insert(count);
                        m_2d_to_3d[*it] = temp_set;
                    } else {
                        m_2d_to_3d[*it].insert(count);
                    }
                }

                saved_pts.push_back(std::make_pair(ptss.at(i).first, segment));
                count++;
            }
        }
    }

    // std::cout << ptss.size() << " " << saved_pts.size() << " " << m_2d_to_3d.size() << " " << m_3d_to_2d.size() << std::endl;
    
    ptss = saved_pts;


    // {
    //     int apa = 0, face = 0;
    //     auto cur_u = m_grouping->convert_3Dpoint_time_ch(pts.front(), apa, face, 0);
    //     auto cur_v = m_grouping->convert_3Dpoint_time_ch(pts.front(), apa, face, 1);
    //     auto cur_w = m_grouping->convert_3Dpoint_time_ch(pts.front(), apa, face, 2);

    //     std::cout << std::get<0>(cur_u) << " " << std::get<1>(cur_u) << " " << std::get<1>(cur_v) << " " << std::get<1>(cur_w)  << std::endl;

    //     WirePlaneId wpid(kAllLayers, face, apa);

    //     auto pt = std::get<0>(wpid_offsets[wpid]) + pts.front().x() * std::get<0>(wpid_slopes[wpid]);
    //     auto pu = std::get<1>(wpid_offsets[wpid]) + std::get<1>(wpid_slopes[wpid]).first * pts.front().y() + std::get<1>(wpid_slopes[wpid]).second * pts.front().z();
    //     auto pv = std::get<2>(wpid_offsets[wpid]) + std::get<2>(wpid_slopes[wpid]).first * pts.front().y() + std::get<2>(wpid_slopes[wpid]).second * pts.front().z();
    //     auto pw = std::get<3>(wpid_offsets[wpid]) + std::get<3>(wpid_slopes[wpid]).first * pts.front().y() + std::get<3>(wpid_slopes[wpid]).second * pts.front().z();

    //     std::cout << pt << " " << pu << " " << pv << " " << pw << std::endl;
    // }
}


// track trajectory fitting // should fit all APA ...
void TrackFitting::trajectory_fit(std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>>& pss_vec, int charge_div_method, double div_sigma){
    if (pss_vec.empty()) return;

    // Create charge division factor maps
    // apa/face --> time/wire, 3D indx --> fac
    std::map<std::pair<int, int>, std::map<std::tuple<int, int, int>, double>> map_Udiv_fac;
    std::map<std::pair<int, int>, std::map<std::tuple<int, int, int>, double>> map_Vdiv_fac;
    std::map<std::pair<int, int>, std::map<std::tuple<int, int, int>, double>> map_Wdiv_fac;

    // Charge division method
    // Equal division
    for (auto it = m_2d_to_3d.begin(); it != m_2d_to_3d.end(); it++) {
        for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++) {
            WirePlaneLayer_t plane = it->first.plane;
            int apa = it->first.apa;
            int face = it->first.face;
            int time = it->first.time;
            int wire = it->first.wire;
            if (plane == WirePlaneLayer_t::kUlayer) {
                map_Udiv_fac[std::make_pair(apa, face)][std::make_tuple(time, wire, *it1)] = 1.0 / it->second.size();
            } else if (plane == WirePlaneLayer_t::kVlayer) {
                map_Vdiv_fac[std::make_pair(apa, face)][std::make_tuple(time, wire, *it1)] = 1.0 / it->second.size();
            } else if (plane == WirePlaneLayer_t::kWlayer) {
                map_Wdiv_fac[std::make_pair(apa, face)][std::make_tuple(time, wire, *it1)] = 1.0 / it->second.size();
            }
        }
    }

    if (charge_div_method == 2) {
        // Use div_sigma for Gaussian weighting
        // Process each plane separately
        std::map<WirePlaneLayer_t, std::map<std::pair<int, int>, std::map<std::tuple<int, int, int>, double>>*> plane_maps = {
            {WirePlaneLayer_t::kUlayer, &map_Udiv_fac},
            {WirePlaneLayer_t::kVlayer, &map_Vdiv_fac},
            {WirePlaneLayer_t::kWlayer, &map_Wdiv_fac}
        };

        for (auto& [plane, div_fac_map] : plane_maps) {
            // Calculate Gaussian weights
            for (auto& [apa_face, coord_idx_fac] : *div_fac_map) {
                double sum = 0;
                int apa = apa_face.first;
                int face = apa_face.second; 

                WirePlaneId wpid(kAllLayers, face, apa);
                auto offset_t = std::get<0>(wpid_offsets[wpid]);
                auto offset_u = std::get<1>(wpid_offsets[wpid]);
                auto offset_v = std::get<2>(wpid_offsets[wpid]);
                auto offset_w = std::get<3>(wpid_offsets[wpid]);
                auto slope_x = std::get<0>(wpid_slopes[wpid]);
                auto slope_yu = std::get<1>(wpid_slopes[wpid]).first;
                auto slope_zu = std::get<1>(wpid_slopes[wpid]).second;
                auto slope_yv = std::get<2>(wpid_slopes[wpid]).first;
                auto slope_zv = std::get<2>(wpid_slopes[wpid]).second;
                auto slope_yw = std::get<3>(wpid_slopes[wpid]).first;
                auto slope_zw = std::get<3>(wpid_slopes[wpid]).second;

                auto time_tick_width = std::get<0>(wpid_geoms[wpid]);
                auto pitch_u = std::get<1>(wpid_geoms[wpid]);
                auto pitch_v = std::get<2>(wpid_geoms[wpid]);
                auto pitch_w = std::get<3>(wpid_geoms[wpid]);

                // Calculate weights
                for (auto& [coord_idx, fac] : coord_idx_fac) {
                    int time = std::get<0>(coord_idx);
                    int wire = std::get<1>(coord_idx);
                    int idx = std::get<2>(coord_idx);

                    double central_t = slope_x * pss_vec[idx].first.x() + offset_t;
                    double central_ch = 0;
                        
                    if (plane == WirePlaneLayer_t::kUlayer) {
                        central_ch = slope_yu * pss_vec[idx].first.y() + slope_zu * pss_vec[idx].first.z() + offset_u;
                    } else if (plane == WirePlaneLayer_t::kVlayer) {
                        central_ch = slope_yv * pss_vec[idx].first.y() + slope_zv * pss_vec[idx].first.z() + offset_v;
                    } else if (plane == WirePlaneLayer_t::kWlayer) {
                        central_ch = slope_yw * pss_vec[idx].first.y() + slope_zw * pss_vec[idx].first.z() + offset_w;
                    }
                        
                    double pitch = (plane == WirePlaneLayer_t::kUlayer) ? pitch_u :
                                    (plane == WirePlaneLayer_t::kVlayer) ? pitch_v : pitch_w;
                        
                    double factor = exp(-0.5 * (pow((central_t - time) * time_tick_width, 2) + pow((central_ch - wire) * pitch, 2)) / pow(div_sigma, 2));
                        
                        fac = factor;
                        sum += factor;
                }
                
                // Normalize weights
                if (sum > 0) {
                    for (auto& [coord_idx, fac] : coord_idx_fac) {
                        fac /= sum; 
                    }
                }
            }
        }
    }
    
    // Main fitting loop using Eigen
    Eigen::VectorXd pos_3D(3 * pss_vec.size());
    
    for (size_t i = 0; i < pss_vec.size(); i++) {
        // Get 2D associations for this 3D point
        const auto& point_info = m_3d_to_2d[i];

        auto segment = pss_vec.at(i).second;
        auto cluster = segment->cluster();
        const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
        double cluster_t0 = cluster->get_flash().time();
        
        auto plane_data_u = point_info.get_plane_data(WirePlaneLayer_t::kUlayer);
        auto plane_data_v = point_info.get_plane_data(WirePlaneLayer_t::kVlayer);
        auto plane_data_w = point_info.get_plane_data(WirePlaneLayer_t::kWlayer);
        
        int n_2D_u = 2 * plane_data_u.associated_2d_points.size();
        int n_2D_v = 2 * plane_data_v.associated_2d_points.size();
        int n_2D_w = 2 * plane_data_w.associated_2d_points.size();

        std::cout << i << " " << n_2D_u << " " << n_2D_v << " " << n_2D_w << std::endl;
        
        Eigen::VectorXd temp_pos_3D(3), data_u_2D(n_2D_u), data_v_2D(n_2D_v), data_w_2D(n_2D_w);
        Eigen::VectorXd temp_pos_3D_init(3);
        Eigen::SparseMatrix<double> RU(n_2D_u, 3);
        Eigen::SparseMatrix<double> RV(n_2D_v, 3);
        Eigen::SparseMatrix<double> RW(n_2D_w, 3);
        
        auto test_wpid = m_dv->contained_by(pss_vec[i].first);
        // Initialization with its raw position
        auto p_raw = transform->backward(pss_vec[i].first, cluster_t0, test_wpid.face(), test_wpid.apa());

        temp_pos_3D_init(0) = p_raw.x();
        temp_pos_3D_init(1) = p_raw.y();
        temp_pos_3D_init(2) = p_raw.z();

        // Initialize data vectors
        data_u_2D.setZero();
        data_v_2D.setZero();
        data_w_2D.setZero();
        
        // Fill U plane data
        int index = 0;
        for (auto it = plane_data_u.associated_2d_points.begin(); it != plane_data_u.associated_2d_points.end(); it++) {
            
            // Get charge measurement
            CoordReadout charge_key(it->apa, it->time, it->channel);
            double charge = 100, charge_err = 1000; // Default values
            
            auto charge_it = m_charge_data.find(charge_key);
            if (charge_it != m_charge_data.end()) {
                charge = charge_it->second.charge;
                charge_err = charge_it->second.charge_err;
            }
            
            if (charge < 100) {
                charge = 100;
                charge_err = 1000;
            }

            // Get division factor
            double div_factor = 1.0;
            auto apa_face_key = std::make_pair(it->apa, it->face);
            auto div_key = std::make_tuple(it->time, it->wire, (int)i);
            auto div_it1 = map_Udiv_fac.find(apa_face_key);
            if (div_it1 != map_Udiv_fac.end()) { 
                auto div_it2 = div_it1->second.find(div_key);
                if (div_it2 != div_it1->second.end()) {
                    div_factor = div_it2->second; 
                }
            }
  
            double scaling = (charge / charge_err) * div_factor;
            
            // Apply quality factor (simplified version)
            if (plane_data_u.quantity < 0.5) {
                if (plane_data_u.quantity != 0) {
                    scaling *= pow(plane_data_u.quantity / 0.5, 1);
                } else {
                    scaling *= 0.05;
                }
            } 

            WirePlaneId wpid(kAllLayers, it->face, it->apa);
            auto offset_t = std::get<0>(wpid_offsets[wpid]);
            auto offset_u = std::get<1>(wpid_offsets[wpid]);
            auto slope_x = std::get<0>(wpid_slopes[wpid]);
            auto slope_yu = std::get<1>(wpid_slopes[wpid]).first;
            auto slope_zu = std::get<1>(wpid_slopes[wpid]).second;
               
            if (scaling != 0) {
                data_u_2D(2 * index) = scaling * (it->wire - offset_u);
                data_u_2D(2 * index + 1) = scaling * (it->time - offset_t);
                
                RU.insert(2 * index, 1) = scaling * slope_yu;     // Y --> U
                RU.insert(2 * index, 2) = scaling * slope_zu;     // Z --> U
                RU.insert(2 * index + 1, 0) = scaling * slope_x;  // X --> T
            }
            // std::cout << index << " U " << n_2D_u << std::endl;
            index++;
        }
        
        // std::cout << "Fill V " << std::endl;
        // Fill V plane data (similar to U)
        index = 0;
        for (auto it = plane_data_v.associated_2d_points.begin(); it != plane_data_v.associated_2d_points.end(); it++) {
            
            // Get charge measurement
            CoordReadout charge_key(it->apa, it->time, it->channel);
            double charge = 100, charge_err = 1000; // Default values
            
            auto charge_it = m_charge_data.find(charge_key);
            if (charge_it != m_charge_data.end()) {
                charge = charge_it->second.charge;
                charge_err = charge_it->second.charge_err;
            }
            
            if (charge < 100) {
                charge = 100;
                charge_err = 1000;
            }

            // Get division factor
            double div_factor = 1.0;
            auto apa_face_key = std::make_pair(it->apa, it->face);
            auto div_key = std::make_tuple(it->time, it->wire, (int)i);
            auto div_it1 = map_Vdiv_fac.find(apa_face_key);
            if (div_it1 != map_Vdiv_fac.end()) { 
                auto div_it2 = div_it1->second.find(div_key);
                if (div_it2 != div_it1->second.end()) {
                    div_factor = div_it2->second; 
                }
            }
  
            double scaling = (charge / charge_err) * div_factor;
            
            // Apply quality factor (simplified version)
            if (plane_data_v.quantity < 0.5) {
                if (plane_data_v.quantity != 0) {
                    scaling *= pow(plane_data_v.quantity / 0.5, 1);
                } else {
                    scaling *= 0.05;
                }
            } 

            WirePlaneId wpid(kAllLayers, it->face, it->apa);
            auto offset_t = std::get<0>(wpid_offsets[wpid]);
            auto offset_v = std::get<2>(wpid_offsets[wpid]);
            auto slope_x = std::get<0>(wpid_slopes[wpid]);
            auto slope_yv = std::get<2>(wpid_slopes[wpid]).first;
            auto slope_zv = std::get<2>(wpid_slopes[wpid]).second;
            
            // std::cout << "Test: " << std::endl;
            if (scaling != 0) {
                data_v_2D(2 * index) = scaling * (it->wire - offset_v);
                data_v_2D(2 * index + 1) = scaling * (it->time - offset_t);
                
                RV.insert(2 * index, 1) = scaling * slope_yv;     // Y --> V
                RV.insert(2 * index, 2) = scaling * slope_zv;     // Z --> V
                RV.insert(2 * index + 1, 0) = scaling * slope_x;  // X --> T
            }
            // std::cout << index << " V " << n_2D_v << std::endl;

            index++;
        }
        
        // std::cout << "Fill W " << std::endl;

        // Fill W plane data (similar to U and V)
        index = 0;
        for (auto it = plane_data_w.associated_2d_points.begin(); it != plane_data_w.associated_2d_points.end(); it++) {
            
            // Get charge measurement
            CoordReadout charge_key(it->apa, it->time, it->channel);
            double charge = 100, charge_err = 1000; // Default values
            
            auto charge_it = m_charge_data.find(charge_key);
            if (charge_it != m_charge_data.end()) {
                charge = charge_it->second.charge;
                charge_err = charge_it->second.charge_err;
            }
            
            if (charge < 100) {
                charge = 100;
                charge_err = 1000;
            }

            // Get division factor
            double div_factor = 1.0;
            auto apa_face_key = std::make_pair(it->apa, it->face);
            auto div_key = std::make_tuple(it->time, it->wire, (int)i);
            auto div_it1 = map_Wdiv_fac.find(apa_face_key);
            if (div_it1 != map_Wdiv_fac.end()) { 
                auto div_it2 = div_it1->second.find(div_key);
                if (div_it2 != div_it1->second.end()) {
                    div_factor = div_it2->second; 
                }
            }
  
            double scaling = (charge / charge_err) * div_factor;
            
            // Apply quality factor (simplified version)
            if (plane_data_w.quantity < 0.5) {
                if (plane_data_w.quantity != 0) {
                    scaling *= pow(plane_data_w.quantity / 0.5, 1);
                } else {
                    scaling *= 0.05;
                }
            } 

            WirePlaneId wpid(kAllLayers, it->face, it->apa);
            auto offset_t = std::get<0>(wpid_offsets[wpid]);
            auto offset_w = std::get<3>(wpid_offsets[wpid]);
            auto slope_x = std::get<0>(wpid_slopes[wpid]);
            auto slope_yw = std::get<3>(wpid_slopes[wpid]).first;
            auto slope_zw = std::get<3>(wpid_slopes[wpid]).second;
            
            if (scaling != 0) {
                data_w_2D(2 * index) = scaling * (it->wire - offset_w);
                data_w_2D(2 * index + 1) = scaling * (it->time - offset_t);
                
                RW.insert(2 * index, 1) = scaling * slope_yw;     // Y --> W
                RW.insert(2 * index, 2) = scaling * slope_zw;     // Z --> W
                RW.insert(2 * index + 1, 0) = scaling * slope_x;  // X --> T
            }
            // std::cout << index << " W " << n_2D_w << std::endl;

            index++;
        }
        
        // Solve the least squares problem
        Eigen::SparseMatrix<double> RUT = RU.transpose();
        Eigen::SparseMatrix<double> RVT = RV.transpose();
        Eigen::SparseMatrix<double> RWT = RW.transpose();
        
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
        Eigen::VectorXd b = RUT * data_u_2D + RVT * data_v_2D + RWT * data_w_2D;
        Eigen::SparseMatrix<double> A = RUT * RU + RVT * RV + RWT * RW;
        
        solver.compute(A);
        temp_pos_3D = solver.solveWithGuess(b, temp_pos_3D_init);
        
        // Store result or use initial position if solver failed
        // these are raw positions ...
        if (std::isnan(solver.error())) {
            pos_3D(3 * i) = temp_pos_3D_init(0);
            pos_3D(3 * i + 1) = temp_pos_3D_init(1);
            pos_3D(3 * i + 2) = temp_pos_3D_init(2);
        } else {
            pos_3D(3 * i) = temp_pos_3D(0);
            pos_3D(3 * i + 1) = temp_pos_3D(1);
            pos_3D(3 * i + 2) = temp_pos_3D(2);
        }
    }
    
    // Clear and rebuild fine tracking path
    fine_tracking_path.clear();
    pu.clear();
    pv.clear();
    pw.clear();
    pt.clear();
    paf.clear();
    
    std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> temp_fine_tracking_path;
    std::vector<std::pair<int, int> > saved_paf;
    int skip_count = 0;
    
    for (size_t i = 0; i < pss_vec.size(); i++) {
        WireCell::Point p_raw(pos_3D(3 * i), pos_3D(3 * i + 1), pos_3D(3 * i + 2));
        auto segment = pss_vec.at(i).second;
        auto cluster = segment->cluster();
        const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
        double cluster_t0 = cluster->get_flash().time();
        auto test_wpid = m_dv->contained_by(pss_vec[i].first);

        auto p = transform->forward(p_raw, cluster_t0, test_wpid.face(), test_wpid.apa());
        auto apa_face = std::make_pair(test_wpid.face(), test_wpid.apa());
        // all corrected points ...
        bool flag_skip = skip_trajectory_point(p, apa_face, i, pss_vec, fine_tracking_path);
        // Protection against too many consecutive skips
        if (flag_skip) {
            skip_count++;
            if (skip_count <= 3) {
                continue;
            } else {
                skip_count = 0;
            }
        }

        // now all corrected points ... 
        temp_fine_tracking_path.push_back(pss_vec[i]);
        fine_tracking_path.push_back(std::make_pair(p, segment));
        saved_paf.push_back(std::make_pair(test_wpid.face(), test_wpid.apa()));
    }
    
    // Apply trajectory smoothing (simplified version of the area-based correction)
    for (size_t i = 0; i < fine_tracking_path.size(); i++) {
        bool flag_replace = false;
        
        // Check triangle area for smoothness (-1, +1 neighbors)
        if (i != 0 && i + 1 != fine_tracking_path.size()) {
            double a = sqrt(pow(fine_tracking_path[i-1].first.x() - fine_tracking_path[i].first.x(), 2) +
                          pow(fine_tracking_path[i-1].first.y() - fine_tracking_path[i].first.y(), 2) +
                          pow(fine_tracking_path[i-1].first.z() - fine_tracking_path[i].first.z(), 2));
            double b = sqrt(pow(fine_tracking_path[i+1].first.x() - fine_tracking_path[i].first.x(), 2) +
                          pow(fine_tracking_path[i+1].first.y() - fine_tracking_path[i].first.y(), 2) +
                          pow(fine_tracking_path[i+1].first.z() - fine_tracking_path[i].first.z(), 2));
            double c = sqrt(pow(fine_tracking_path[i-1].first.x() - fine_tracking_path[i+1].first.x(), 2) +
                          pow(fine_tracking_path[i-1].first.y() - fine_tracking_path[i+1].first.y(), 2) +
                          pow(fine_tracking_path[i-1].first.z() - fine_tracking_path[i+1].first.z(), 2));
            
            if (c > 0) {
                double s = (a + b + c) / 2.0;
                double area1 = sqrt(s * (s - a) * (s - b) * (s - c));
                
                // Compare with original point
                a = sqrt(pow(fine_tracking_path[i-1].first.x() - temp_fine_tracking_path[i].first.x(), 2) +
                       pow(fine_tracking_path[i-1].first.y() - temp_fine_tracking_path[i].first.y(), 2) +
                       pow(fine_tracking_path[i-1].first.z() - temp_fine_tracking_path[i].first.z(), 2));
                b = sqrt(pow(fine_tracking_path[i+1].first.x() - temp_fine_tracking_path[i].first.x(), 2) +
                       pow(fine_tracking_path[i+1].first.y() - temp_fine_tracking_path[i].first.y(), 2) +
                       pow(fine_tracking_path[i+1].first.z() - temp_fine_tracking_path[i].first.z(), 2));
                
                s = (a + b + c) / 2.0;
                double area2 = sqrt(s * (s - a) * (s - b) * (s - c));
                
                if (area1 > 1.8 * units::mm * c && area1 > 1.7 * area2) {
                    flag_replace = true;
                }
            }
        }
        
        if (flag_replace) {
            fine_tracking_path[i] = temp_fine_tracking_path[i];
        }
    }
    
    // Generate 2D projections
    for (size_t i = 0; i < fine_tracking_path.size(); i++) {
        WireCell::Point p = fine_tracking_path[i].first;
        auto segment = fine_tracking_path[i].second;
        auto cluster = segment->cluster();
        const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
        double cluster_t0 = cluster->get_flash().time();

        int apa = saved_paf.at(i).first;
        int face = saved_paf.at(i).second;

        auto p_raw = transform->backward(p, cluster_t0, apa, face);
        WirePlaneId wpid(kAllLayers, face, apa);
        auto offset_t = std::get<0>(wpid_offsets[wpid]);
        auto offset_u = std::get<1>(wpid_offsets[wpid]);
        auto offset_v = std::get<2>(wpid_offsets[wpid]);
        auto offset_w = std::get<3>(wpid_offsets[wpid]);
        auto slope_x = std::get<0>(wpid_slopes[wpid]);
        auto slope_yu = std::get<1>(wpid_slopes[wpid]).first;
        auto slope_zu = std::get<1>(wpid_slopes[wpid]).second;
        auto slope_yv = std::get<2>(wpid_slopes[wpid]).first;
        auto slope_zv = std::get<2>(wpid_slopes[wpid]).second;
        auto slope_yw = std::get<3>(wpid_slopes[wpid]).first;
        auto slope_zw = std::get<3>(wpid_slopes[wpid]).second;

        pu.push_back(offset_u + (slope_yu * p_raw.y() + slope_zu * p_raw.z()));
        pv.push_back(offset_v + (slope_yv * p_raw.y() + slope_zv * p_raw.z()));
        pw.push_back(offset_w + (slope_yw * p_raw.y() + slope_zw * p_raw.z()));
        pt.push_back(offset_t + slope_x * p_raw.x());
        paf.push_back(std::make_pair(apa, face));

    }
    
    // Update the input vector with the fitted results
    pss_vec = fine_tracking_path;
}

bool TrackFitting::skip_trajectory_point(WireCell::Point& p, std::pair<int, int>& apa_face, int i, std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>>& pss_vec,  std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>>& fine_tracking_path){
    return false;
}