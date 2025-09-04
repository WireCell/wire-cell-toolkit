#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/TrackFitting_Util.h"
#include "WireCellUtil/Logging.h"


using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

using geo_point_t = WireCell::Point;


TrackFitting::TrackFitting(FittingType fitting_type) 
    : m_fitting_type(fitting_type) 
{

}

// ============================================================================
// Parameter management methods
// ============================================================================

void TrackFitting::set_parameter(const std::string& name, double value) {
    // Map parameter names to struct members
    if (name == "DL") {
        m_params.DL = value;
    } else if (name == "DT") {
        m_params.DT = value;
    } else if (name == "col_sigma_w_T") {
        m_params.col_sigma_w_T = value;
    } else if (name == "ind_sigma_u_T") {
        m_params.ind_sigma_u_T = value;
    } else if (name == "ind_sigma_v_T") {
        m_params.ind_sigma_v_T = value;
    } else if (name == "rel_uncer_ind") {
        m_params.rel_uncer_ind = value;
    } else if (name == "rel_uncer_col") {
        m_params.rel_uncer_col = value;
    } else if (name == "add_uncer_ind") {
        m_params.add_uncer_ind = value;
    } else if (name == "add_uncer_col") {
        m_params.add_uncer_col = value;
    } else if (name == "add_sigma_L") {
        m_params.add_sigma_L = value;
    } else if (name == "low_dis_limit") {
        m_params.low_dis_limit = value;
    } else if (name == "end_point_limit") {
        m_params.end_point_limit = value;
    } else if (name == "time_tick_cut") {
        m_params.time_tick_cut = value;
    } else if (name == "rel_charge_uncer") {
        m_params.rel_charge_uncer = value;
    } else if (name == "add_charge_uncer") {
        m_params.add_charge_uncer = value;
    } else if (name == "default_charge_th") {
        m_params.default_charge_th = value;
    } else if (name == "default_charge_err") {
        m_params.default_charge_err = value;
    } else if (name == "scaling_quality_th") {
        m_params.scaling_quality_th = value;
    } else if (name == "scaling_ratio") {
        m_params.scaling_ratio = value;
    } else if (name == "area_ratio1") {
        m_params.area_ratio1 = value;
    } else if (name == "area_ratio2") {
        m_params.area_ratio2 = value;
    } else if (name == "skip_default_ratio_1") {
        m_params.skip_default_ratio_1 = value;
    } else if (name == "skip_ratio_cut") {
        m_params.skip_ratio_cut = value;
    } else if (name == "skip_ratio_1_cut") {
        m_params.skip_ratio_1_cut = value;
    } else if (name == "skip_angle_cut_1") {
        m_params.skip_angle_cut_1 = value;
    } else if (name == "skip_angle_cut_2") {
        m_params.skip_angle_cut_2 = value;
    } else if (name == "skip_angle_cut_3") {
        m_params.skip_angle_cut_3 = value;
    } else if (name == "skip_dis_cut") {
        m_params.skip_dis_cut = value;
    } else if (name == "default_dQ_dx") {
        m_params.default_dQ_dx = value;
    } else if (name == "end_point_factor") {
        m_params.end_point_factor = value;
    } else if (name == "mid_point_factor") {
        m_params.mid_point_factor = value;
    } else if (name == "nlevel") {
        m_params.nlevel = static_cast<int>(value);
    } else if (name == "charge_cut") {
        m_params.charge_cut = value;
    } else if (name == "share_charge_err") {
        m_params.share_charge_err = value;
    } else if (name == "min_drift_time") {
        m_params.min_drift_time = value;
    } else if (name == "search_range") {
        m_params.search_range = value;
    } else if (name == "dead_ind_weight") {
        m_params.dead_ind_weight = value;
    } else if (name == "dead_col_weight") {
        m_params.dead_col_weight = value;
    } else if (name == "close_ind_weight") {
        m_params.close_ind_weight = value;
    } else if (name == "close_col_weight") {
        m_params.close_col_weight = value;
    } else if (name == "overlap_th") {
        m_params.overlap_th = value;
    } else if (name == "dx_norm_length") {
        m_params.dx_norm_length = value;
    } else if (name == "lambda") {
        m_params.lambda = value;
    } else if (name == "div_sigma") {
        m_params.div_sigma = value;
    } else {
        raise<ValueError>("TrackFitting: Unknown parameter name '%s'", name.c_str());
    }
}

double TrackFitting::get_parameter(const std::string& name) const {
    // Map parameter names to struct members
    if (name == "DL") {
        return m_params.DL;
    } else if (name == "DT") {
        return m_params.DT;
    } else if (name == "col_sigma_w_T") {
        return m_params.col_sigma_w_T;
    } else if (name == "ind_sigma_u_T") {
        return m_params.ind_sigma_u_T;
    } else if (name == "ind_sigma_v_T") {
        return m_params.ind_sigma_v_T;
    } else if (name == "rel_uncer_ind") {
        return m_params.rel_uncer_ind;
    } else if (name == "rel_uncer_col") {
        return m_params.rel_uncer_col;
    } else if (name == "add_uncer_ind") {
        return m_params.add_uncer_ind;
    } else if (name == "add_uncer_col") {
        return m_params.add_uncer_col;
    } else if (name == "add_sigma_L") {
        return m_params.add_sigma_L;
    } else if (name == "low_dis_limit") {
        return m_params.low_dis_limit;
    } else if (name == "end_point_limit") {
        return m_params.end_point_limit;
    } else if (name == "time_tick_cut") {
        return m_params.time_tick_cut;
    } else if (name == "rel_charge_uncer") {
        return m_params.rel_charge_uncer;
    } else if (name == "add_charge_uncer") {
        return m_params.add_charge_uncer;
    } else if (name == "default_charge_th") {
        return m_params.default_charge_th;
    } else if (name == "default_charge_err") {
        return m_params.default_charge_err;
    } else if (name == "scaling_quality_th") {
        return m_params.scaling_quality_th;
    } else if (name == "scaling_ratio") {
        return m_params.scaling_ratio;
    } else if (name == "area_ratio1") {
        return m_params.area_ratio1;
    } else if (name == "area_ratio2") {
        return m_params.area_ratio2;
    } else if (name == "skip_default_ratio_1") {
        return m_params.skip_default_ratio_1;
    } else if (name == "skip_ratio_cut") {
        return m_params.skip_ratio_cut;
    } else if (name == "skip_ratio_1_cut") {
        return m_params.skip_ratio_1_cut;
    } else if (name == "skip_angle_cut_1") {
        return m_params.skip_angle_cut_1;
    } else if (name == "skip_angle_cut_2") {
        return m_params.skip_angle_cut_2;
    } else if (name == "skip_angle_cut_3") {
        return m_params.skip_angle_cut_3;
    } else if (name == "skip_dis_cut") {
        return m_params.skip_dis_cut;
    } else if (name == "default_dQ_dx") {
        return m_params.default_dQ_dx;
    } else if (name == "end_point_factor") {
        return m_params.end_point_factor;
    } else if (name == "mid_point_factor") {
        return m_params.mid_point_factor;
    } else if (name == "nlevel") {
        return static_cast<double>(m_params.nlevel);
    } else if (name == "charge_cut") {
        return m_params.charge_cut;
    } else if (name == "share_charge_err") {
        return m_params.share_charge_err;
    } else if (name == "min_drift_time") {
        return m_params.min_drift_time;
    } else if (name == "search_range") {
        return m_params.search_range;
    } else if (name == "dead_ind_weight") {
        return m_params.dead_ind_weight;
    } else if (name == "dead_col_weight") {
        return m_params.dead_col_weight;
    } else if (name == "close_ind_weight") {
        return m_params.close_ind_weight;
    } else if (name == "close_col_weight") {
        return m_params.close_col_weight;
    } else if (name == "overlap_th") {
        return m_params.overlap_th;
    } else if (name == "dx_norm_length") {
        return m_params.dx_norm_length;
    } else if (name == "lambda") {
        return m_params.lambda;
    } else if (name == "div_sigma") {
        return m_params.div_sigma;
    } else {
        raise<ValueError>("TrackFitting: Unknown parameter name '%s'", name.c_str());
        return 0;
    }
}

void TrackFitting::clear_segments(){
    m_segments.clear();
    m_clusters.clear();
    m_blobs.clear(); 
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
    
    // std::cout << "TrackFitting: Promoted plane (" << apa << "," << face << "," << plane 
            //   << ") to hot cache with " << wires.size() << " wires" << std::endl;
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

            u_min -= 5; v_min -=5; w_min-=5;
            u_max += 5; v_max +=5; w_max+=5;
            t_min -= 20;
            t_max += 20;
            // std::cout << "U Limits: " << u_min << " " << u_max << std::endl;
            // std::cout << "V Limits: " << v_min << " " << v_max << std::endl;
            // std::cout << "W Limits: " << w_min << " " << w_max << std::endl;

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
                            double charge_err = sqrt(pow(charge * m_params.rel_charge_uncer, 2) + pow(m_params.add_charge_uncer, 2));
                            m_charge_data[data_key] = {charge, charge_err, 0};
                        } else if (it->second.flag == 0) {
                            // Existing content with flag = 0
                            double new_charge = blob_charge / num_wires;
                            double new_charge_err = sqrt(pow(new_charge * m_params.rel_charge_uncer, 2) + pow(m_params.add_charge_uncer, 2));
                            
                            it->second.charge += new_charge;
                            it->second.charge_err = sqrt(pow(it->second.charge_err, 2) + pow(new_charge_err, 2));
                        }
                        // If flag != 0, do nothing
                    }
                }
            }
        }
    }


    // std::cout << "Number of Measurements: " << m_charge_data.size() << std::endl;
    // for (const auto& [coord_key, charge_measurement] : m_charge_data) {
    //     std::cout << "CoordReadout: (APA=" << coord_key.apa
    //               << ", Time=" << coord_key.time
    //               << ", Channel=" << coord_key.channel
    //               << ") -> Charge=" << charge_measurement.charge
    //               << ", ChargeErr=" << charge_measurement.charge_err
    //               << ", Flag=" << charge_measurement.flag << std::endl;
    // }
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

    // std::cout << "Test b: " <<  pts.size() << " " << pts.back() << " " << temp_wcps_vec.front() << std::endl;
    
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

    // std::cout << "Test m: " <<  pts.size() << " " << pts.back() << " " << temp_wcps_vec.back() << std::endl;

    
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

    // std::cout << "Test e: " <<  pts.size() << " " << pts.back() << " " << temp_wcps_vec.back() << std::endl;

    
    return pts;
}

std::vector<WireCell::Point> TrackFitting::examine_end_ps_vec(std::shared_ptr<PR::Segment> segment,const std::vector<WireCell::Point>& pts, bool flag_start, bool flag_end) {
    std::list<WireCell::Point> ps_list(pts.begin(), pts.end());
    
    // get the cluster from the segment
    auto cluster = segment->cluster();
    const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
    double cluster_t0 = cluster->get_cluster_t0();

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
        
        // std::cout << i << " " << dis << " " << low_dis_limit * 0.8 << " " << low_dis_limit * 1.6 << std::endl;

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
    double cluster_t0 = cluster->get_cluster_t0();
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
        
       
        // std::cout << nearby_blobs_set.size() << " nearby blobs found for point " <<  cur_time_slice << " " << cur_wire_u << " " << cur_wire_v << " " << cur_wire_w << " " << dis_cut <<  std::endl;
        
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

        // std::cout << dis_cut_u << " " << dis_cut_v << " " << dis_cut_w << std::endl;
        
        // Process each nearby blob for wire range calculations
        for (const auto* blob : nearby_blobs_set) {
            int this_time_slice = blob->slice_index_min();

            // std::cout << "Blob info: " << blob->u_wire_index_min() << " " << blob->u_wire_index_max() << " " << blob->v_wire_index_min() << " " << blob->v_wire_index_max() << " " << blob->w_wire_index_min() << " " << blob->w_wire_index_max() << " " << this_time_slice << std::endl; 
            
            // Calculate remaining distance cuts accounting for time offset
            double rem_dis_sq_cut_u = pow(dis_cut_u, 2) - pow((cur_time_slice - this_time_slice) * time_tick_width, 2);
            double rem_dis_sq_cut_v = pow(dis_cut_v, 2) - pow((cur_time_slice - this_time_slice) * time_tick_width, 2);
            double rem_dis_sq_cut_w = pow(dis_cut_w, 2) - pow((cur_time_slice - this_time_slice) * time_tick_width, 2);

            // std::cout << rem_dis_cut_u << " " << rem_dis_cut_v << " " << rem_dis_cut_w << " " << cur_time_slice << " " <<this_time_slice << " " << time_tick_cut << std::endl;

            if ((rem_dis_sq_cut_u > 0 || rem_dis_sq_cut_v > 0 || rem_dis_sq_cut_w > 0) && abs(cur_time_slice - this_time_slice) <= time_tick_cut) {

                // Calculate minimum wire distances
                float min_u_dis, min_v_dis, min_w_dis;
                
                // U wire distance
                if (cur_wire_u < blob->u_wire_index_min()) {
                    min_u_dis = blob->u_wire_index_min() - cur_wire_u;
                } else if (cur_wire_u < blob->u_wire_index_max()) {
                    min_u_dis = 0;
                } else {
                    min_u_dis = cur_wire_u - blob->u_wire_index_max()+1;
                }
                
                // V wire distance
                if (cur_wire_v < blob->v_wire_index_min()) {
                    min_v_dis = blob->v_wire_index_min() - cur_wire_v;
                } else if (cur_wire_v < blob->v_wire_index_max()) {
                    min_v_dis = 0;
                } else {
                    min_v_dis = cur_wire_v - blob->v_wire_index_max()+1;
                }
                
                // W wire distance
                if (cur_wire_w < blob->w_wire_index_min()) {
                    min_w_dis = blob->w_wire_index_min() - cur_wire_w;
                } else if (cur_wire_w < blob->w_wire_index_max()) {
                    min_w_dis = 0;
                } else {
                    min_w_dis = cur_wire_w - blob->w_wire_index_max()+1;
                }
                
                // Use the dedicated calculate_ranges_simplified function
                float range_sq_u, range_sq_v, range_sq_w;
                WireCell::Clus::TrackFittingUtil::calculate_ranges_simplified(
                    angle_u, angle_v, angle_w,
                    rem_dis_sq_cut_u, rem_dis_sq_cut_v, rem_dis_sq_cut_w,
                    min_u_dis, min_v_dis, min_w_dis,
                    pitch_u, pitch_v, pitch_w,
                    range_sq_u, range_sq_v, range_sq_w);

                // std::cout << "Cuts: " << range_sq_u << " " << range_sq_v << " " << range_sq_w << std::endl;


                // If all ranges are positive, add wire indices to associations
                if (range_sq_u > 0 && range_sq_v > 0 && range_sq_w > 0) {
                    // Calculate wire limits
                    float low_u_limit = cur_wire_u - sqrt(range_sq_u) / pitch_u;
                    float high_u_limit = cur_wire_u + sqrt(range_sq_u) / pitch_u;
                    float low_v_limit = cur_wire_v - sqrt(range_sq_v) / pitch_v;
                    float high_v_limit = cur_wire_v + sqrt(range_sq_v) / pitch_v;
                    float low_w_limit = cur_wire_w - sqrt(range_sq_w) / pitch_w;
                    float high_w_limit = cur_wire_w + sqrt(range_sq_w) / pitch_w;

                    // std::cout << low_u_limit << " " << high_u_limit << " "
                    //           << low_v_limit << " " << high_v_limit << " "
                    //           << low_w_limit << " " << high_w_limit << " " << this_time_slice << std::endl;

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

        // std::cout << "Steiner " << temp_dis << " " << dis_cut << " " <<apa << " " << closest_point_wpid.apa() << " " << face << " " << closest_point_wpid.face() << std::endl;

        if (temp_dis < dis_cut && apa == closest_point_wpid.apa() && face == closest_point_wpid.face()){
            // Get graph algorithms interface
            const auto& ga = cluster->graph_algorithms(graph_name);

            // // Print Steiner graph statistics
            // const auto& graph_steiner = cluster->get_graph("steiner_graph");
            // std::cout << "Steiner graph: vertices = " << boost::num_vertices(graph_steiner)
            //           << ", edges = " << boost::num_edges(graph_steiner) << std::endl;
            
            // Find nearby points using graph traversal (equivalent to original nested loop)
            auto total_vertices_found = ga.find_neighbors_nlevel(closest_point_index, nlevel);

            // find the raw point ...
            auto closest_point_raw = transform->backward(closest_point, cluster_t0, apa, face);

            // std::cout << p << " " << closest_point << " "  << closest_point_raw << std::endl;

            auto cur_u = m_grouping->convert_3Dpoint_time_ch(closest_point_raw, apa, face, 0);
            auto cur_v = m_grouping->convert_3Dpoint_time_ch(closest_point_raw, apa, face, 1);
            auto cur_w = m_grouping->convert_3Dpoint_time_ch(closest_point_raw, apa, face, 2);

            int cur_time_slice = std::floor(std::get<0>(cur_u)/cur_ntime_ticks)*cur_ntime_ticks;
            int cur_wire_u = std::get<1>(cur_u);
            int cur_wire_v = std::get<1>(cur_v);
            int cur_wire_w = std::get<1>(cur_w);

            // std::cout << "B: " << cluster_t0 << " " << std::get<0>(cur_u) << " " << cur_time_slice << " " << cur_wire_u << " " << cur_wire_v << " " << cur_wire_w << " " << total_vertices_found.size() << " " << nlevel << std::endl;
        
            // Calculate adaptive distance cuts (equivalent to original max_time_slice_u/v/w calculation)
            double dis_cut_u = dis_cut;
            double dis_cut_v = dis_cut;
            double dis_cut_w = dis_cut;
            
            double max_time_slice_u = 0;
            double max_time_slice_v = 0;
            double max_time_slice_w = 0;

            std::map<int, std::tuple<int, int, int, int, int, int> > map_time_wires;

            // std::map<int, std::tuple<int, int, int, int>> map_vertex_info;
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

                    int umin = vertex_wire_u, umax = vertex_wire_u+1;
                    int vmin = vertex_wire_v, vmax = vertex_wire_v+1;
                    int wmin = vertex_wire_w, wmax = vertex_wire_w+1;
                    auto it = map_time_wires.find(vertex_time_slice);
                    if (it == map_time_wires.end()) {
                        // No entry yet, insert new boundaries
                        map_time_wires[vertex_time_slice] = std::make_tuple(umin, umax, vmin, vmax, wmin, wmax);
                    } else {
                        // Update boundaries if needed
                        auto& tup = it->second;
                        std::get<0>(tup) = std::min(std::get<0>(tup), umin);
                        std::get<1>(tup) = std::max(std::get<1>(tup), umax);
                        std::get<2>(tup) = std::min(std::get<2>(tup), vmin);
                        std::get<3>(tup) = std::max(std::get<3>(tup), vmax);
                        std::get<4>(tup) = std::min(std::get<4>(tup), wmin);
                        std::get<5>(tup) = std::max(std::get<5>(tup), wmax);
                    }
            }

            for (const auto& [vertex_time_slice, wire_ranges] : map_time_wires) {
                int umin = std::get<0>(wire_ranges);
                int umax = std::get<1>(wire_ranges);
                int vmin = std::get<2>(wire_ranges);
                int vmax = std::get<3>(wire_ranges);
                int wmin = std::get<4>(wire_ranges);
                int wmax = std::get<5>(wire_ranges);

                for (auto vertex_wire_u = umin; vertex_wire_u < umax; ++vertex_wire_u) {
                    // Check U, V, W wire proximity
                    if (abs(vertex_wire_u - cur_wire_u)*pitch_u<=dis_cut) {
                        if (abs(vertex_time_slice - cur_time_slice) > max_time_slice_u)
                            max_time_slice_u = abs(vertex_time_slice - cur_time_slice);
                    }
                }
                for (auto vertex_wire_v = vmin; vertex_wire_v < vmax; ++vertex_wire_v) {
                    // Check U, V, W wire proximity
                    if (abs(vertex_wire_v - cur_wire_v)*pitch_v<=dis_cut) {
                        if (abs(vertex_time_slice - cur_time_slice) > max_time_slice_v)
                            max_time_slice_v = abs(vertex_time_slice - cur_time_slice);
                    }
                }
                for (auto vertex_wire_w = wmin; vertex_wire_w < wmax; ++vertex_wire_w) {
                    // Check U, V, W wire proximity
                    if (abs(vertex_wire_w - cur_wire_w)*pitch_w<=dis_cut) {
                        if (abs(vertex_time_slice - cur_time_slice) > max_time_slice_w)
                            max_time_slice_w = abs(vertex_time_slice - cur_time_slice);
                    }
                }
            }

            // Apply adaptive cuts (equivalent to original adaptive cut logic)
            if (max_time_slice_u * time_tick_width * 1.2 < dis_cut_u)
                dis_cut_u = max_time_slice_u * time_tick_width * 1.2;
            if (max_time_slice_v * time_tick_width * 1.2 < dis_cut_v)
                dis_cut_v = max_time_slice_v * time_tick_width * 1.2;
            if (max_time_slice_w * time_tick_width * 1.2 < dis_cut_w)
                dis_cut_w = max_time_slice_w * time_tick_width * 1.2;

            // std::cout << "Steiner dis: " << dis_cut_u << " " << dis_cut_v << " " << dis_cut_w << " " << std::endl;

            // Process each nearby blob for wire range calculations (equivalent to final loop)
            for (const auto& [vertex_time_slice, wire_ranges] : map_time_wires) {
                int umin = std::get<0>(wire_ranges);
                int umax = std::get<1>(wire_ranges);
                int vmin = std::get<2>(wire_ranges); 
                int vmax = std::get<3>(wire_ranges);
                int wmin = std::get<4>(wire_ranges);
                int wmax = std::get<5>(wire_ranges);
        

                // Calculate remaining distance cuts accounting for time offset
                double rem_dis_cut_u = pow(dis_cut_u, 2) - pow((cur_time_slice - vertex_time_slice) * time_tick_width, 2);
                double rem_dis_cut_v = pow(dis_cut_v, 2) - pow((cur_time_slice - vertex_time_slice) * time_tick_width, 2);
                double rem_dis_cut_w = pow(dis_cut_w, 2) - pow((cur_time_slice - vertex_time_slice) * time_tick_width, 2);

                // std::cout << rem_dis_cut_u << " " << rem_dis_cut_v << " " << rem_dis_cut_w << " " << std::endl;

                if ((rem_dis_cut_u > 0 || rem_dis_cut_v > 0 || rem_dis_cut_w > 0) && abs(cur_time_slice - vertex_time_slice) <= time_tick_cut) {
                
                // Calculate minimum wire distances (equivalent to original min_u/v/w_dis calculation)
                float min_u_dis, min_v_dis, min_w_dis;
                
                // U wire distance
                if (cur_wire_u < umin) {
                    min_u_dis = umin - cur_wire_u;
                } else if (cur_wire_u >= umax) {
                    min_u_dis = cur_wire_u - umax + 1;
                } else {
                    min_u_dis = 0;
                }

                // V wire distance
                if (cur_wire_v < vmin) {
                    min_v_dis = vmin - cur_wire_v;
                } else if (cur_wire_v >= vmax) {
                    min_v_dis = cur_wire_v - vmax + 1;
                } else {
                    min_v_dis = 0;
                }

                // W wire distance
                if (cur_wire_w < wmin) {
                    min_w_dis = wmin - cur_wire_w;
                } else if (cur_wire_w >= wmax) {
                    min_w_dis = cur_wire_w - wmax + 1;
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

                // std::cout << min_u_dis << " " << min_v_dis << " " << min_w_dis << " " << range_u << " " << range_v << " " << range_w << std::endl;
                
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

    // std::cout << "Pixels 2: " << temp_2dut.associated_2d_points.size() << " " << temp_2dvt.associated_2d_points.size() << " " << temp_2dwt.associated_2d_points.size() << std::endl;


 }


 void TrackFitting::examine_point_association(std::shared_ptr<PR::Segment> segment, WireCell::Point &p, PlaneData& temp_2dut, PlaneData& temp_2dvt, PlaneData& temp_2dwt, bool flag_end_point, double charge_cut){

    // Get cluster from segment
    auto cluster = segment->cluster();
    const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
    double cluster_t0 = cluster->get_cluster_t0();

    auto first_blob = segment->cluster()->children()[0];
    int cur_ntime_ticks = first_blob->slice_index_max() - first_blob->slice_index_min();

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
        // std::cout << "V: " << it->time/4 << " " << it->channel << " " << std::endl;
        if (charge_it != m_charge_data.end() && charge_it->second.charge > charge_cut) {
        // std::cout << "V: " << it->time/4 << " " << it->channel << " " << charge_it->second.charge << " " << charge_cut << std::endl;

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

    // std::cout << saved_2dut.size() << " " << saved_2dvt.size() << " " << saved_2dwt.size() << " " << charge_cut << std::endl;


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
                rms += pow(it1->wire - ave_pos.first, 2) + pow((it1->time - ave_pos.second)/cur_ntime_ticks, 2);
            }
            rms = sqrt(rms/saved_2dwt.size());

            if (sqrt(pow(ave_pos.first - wire_w, 2) + pow((ave_pos.second - time_w)/cur_ntime_ticks, 2)) > 0.75*rms && 
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
                rms += pow(it1->wire - ave_pos.first, 2) + pow((it1->time - ave_pos.second)/cur_ntime_ticks, 2);
            }
            rms = sqrt(rms/saved_2dvt.size());

            if (sqrt(pow(ave_pos.first - wire_v, 2) + pow((ave_pos.second - time_v)/cur_ntime_ticks, 2)) > 0.75*rms && 
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
                rms += pow(it1->wire - ave_pos.first, 2) + pow((it1->time - ave_pos.second)/cur_ntime_ticks, 2);
            }
            rms = sqrt(rms/saved_2dut.size());

            if (sqrt(pow(ave_pos.first - wire_u, 2) + pow((ave_pos.second - time_u)/cur_ntime_ticks, 2)) > 0.75*rms && 
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
                    rms += pow(it1->wire - ave_pos.first, 2) + pow((it1->time - ave_pos.second)/cur_ntime_ticks, 2);
                }
                rms = sqrt(rms/saved_plane.size());

                if (sqrt(pow(ave_pos.first - expected_wire, 2) + pow((ave_pos.second - expected_time)/cur_ntime_ticks, 2)) > 0.75*rms && 
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
                    rms += pow(it1->wire - ave_pos.first, 2) + pow((it1->time - ave_pos.second)/cur_ntime_ticks, 2);
                }
                rms = sqrt(rms/saved_plane.size());

                if (sqrt(pow(ave_pos.first - expected_wire, 2) + pow((ave_pos.second - expected_time)/cur_ntime_ticks, 2)) > 0.75*rms && 
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
                    rms += pow(it1->wire - ave_pos.first, 2) + pow((it1->time - ave_pos.second)/cur_ntime_ticks, 2);
                }
                rms = sqrt(rms/saved_plane.size());

                if (sqrt(pow(ave_pos.first - expected_wire, 2) + pow((ave_pos.second - expected_time)/cur_ntime_ticks, 2)) > 0.75*rms && 
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

        // std::cout << i << " " << distances.at(i) << " " << end_point_factor << " " << dis_cut << std::endl;

        // check point's apa and face ...
        // find the apa and face ...
        auto wpid = m_dv->contained_by(ptss.at(i).first);
        auto segment = ptss.at(i).second;
        int apa = wpid.apa();
        int face = wpid.face();
        
        if (apa != -1 && face != -1) {

            TrackFitting::PlaneData temp_2dut, temp_2dvt, temp_2dwt;
            form_point_association(segment, ptss.at(i).first, temp_2dut, temp_2dvt, temp_2dwt, dis_cut, nlevel, time_tick_cut);

            // std::cout << i << " " << ptss.at(i).first << " " << temp_2dut.associated_2d_points.size() << " " << temp_2dvt.associated_2d_points.size() << " " << temp_2dwt.associated_2d_points.size() << std::endl;

            if (i == 0 || i == 1 || i + 1 == ptss.size() || i + 2 == ptss.size()) {
                examine_point_association(segment, ptss.at(i).first, temp_2dut, temp_2dvt, temp_2dwt, true, charge_cut);
            } else {
                examine_point_association(segment, ptss.at(i).first, temp_2dut, temp_2dvt, temp_2dwt, false, charge_cut);
            }
            // std::cout << i << " E " << ptss.at(i).first << " " << temp_2dut.associated_2d_points.size() << " " << temp_2dvt.associated_2d_points.size() << " " << temp_2dwt.associated_2d_points.size() << std::endl;


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

    // std::cout << "Form Map: " << ptss.size() << " " << saved_pts.size() << " " << m_2d_to_3d.size() << " " << m_3d_to_2d.size() << std::endl;
    
    ptss = saved_pts;


    // {
    //     int apa = 0, face = 0;
    //     auto cur_u = m_grouping->convert_3Dpoint_time_ch(ptss.back().first, apa, face, 0);
    //     auto cur_v = m_grouping->convert_3Dpoint_time_ch(ptss.back().first, apa, face, 1);
    //     auto cur_w = m_grouping->convert_3Dpoint_time_ch(ptss.back().first, apa, face, 2);

    //     std::cout << std::get<0>(cur_u) << " " << std::get<1>(cur_u) << " " << std::get<1>(cur_v) << " " << std::get<1>(cur_w)  << std::endl;

    //     WirePlaneId wpid(kAllLayers, face, apa);

    //     auto pt = std::get<0>(wpid_offsets[wpid]) + ptss.back().first.x() * std::get<0>(wpid_slopes[wpid]);
    //     auto pu = std::get<1>(wpid_offsets[wpid]) + std::get<1>(wpid_slopes[wpid]).first * ptss.back().first.y() + std::get<1>(wpid_slopes[wpid]).second * ptss.back().first.z();
    //     auto pv = std::get<2>(wpid_offsets[wpid]) + std::get<2>(wpid_slopes[wpid]).first * ptss.back().first.y() + std::get<2>(wpid_slopes[wpid]).second * ptss.back().first.z();
    //     auto pw = std::get<3>(wpid_offsets[wpid]) + std::get<3>(wpid_slopes[wpid]).first * ptss.back().first.y() + std::get<3>(wpid_slopes[wpid]).second * ptss.back().first.z();

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
                int apa = apa_face.first;
                int face = apa_face.second; 

                WirePlaneId wpid(kAllLayers, face, apa);
                auto offset_it = wpid_offsets.find(wpid);
                auto slope_it = wpid_slopes.find(wpid);
                auto geom_it = wpid_geoms.find(wpid);

                auto offset_t = std::get<0>(offset_it->second);
                auto offset_u = std::get<1>(offset_it->second);
                auto offset_v = std::get<2>(offset_it->second);
                auto offset_w = std::get<3>(offset_it->second);
                auto slope_x = std::get<0>(slope_it->second);
                auto slope_yu = std::get<1>(slope_it->second).first;
                auto slope_zu = std::get<1>(slope_it->second).second;
                auto slope_yv = std::get<2>(slope_it->second).first;
                auto slope_zv = std::get<2>(slope_it->second).second;
                auto slope_yw = std::get<3>(slope_it->second).first;
                auto slope_zw = std::get<3>(slope_it->second).second;

                auto time_tick_width = std::get<0>(geom_it->second);
                auto pitch_u = std::get<1>(geom_it->second);
                auto pitch_v = std::get<2>(geom_it->second);
                auto pitch_w = std::get<3>(geom_it->second);

               
                // double sum = 0;
                std::map<std::pair<int, int>, double> map_tw_sum;

                // Calculate weights
                for (auto& [coord_idx, fac] : coord_idx_fac) {
                    int time = std::get<0>(coord_idx);
                    int wire = std::get<1>(coord_idx);
                    int idx = std::get<2>(coord_idx);

                    

                    auto segment = pss_vec.at(idx).second;
                    auto cluster = segment->cluster();
                    const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
                    double cluster_t0 = cluster->get_cluster_t0();

                    auto test_wpid = m_dv->contained_by(pss_vec[idx].first);
                    // Initialization with its raw position
                    auto p_raw = transform->backward(pss_vec[idx].first, cluster_t0, test_wpid.face(), test_wpid.apa());

                    double central_t = slope_x * p_raw.x() + offset_t;
                    double central_ch = 0;
                        
                    if (plane == WirePlaneLayer_t::kUlayer) {
                        central_ch = slope_yu * p_raw.y() + slope_zu * p_raw.z() + offset_u;
                    } else if (plane == WirePlaneLayer_t::kVlayer) {
                        central_ch = slope_yv * p_raw.y() + slope_zv * p_raw.z() + offset_v;
                    } else if (plane == WirePlaneLayer_t::kWlayer) {
                        central_ch = slope_yw * p_raw.y() + slope_zw * p_raw.z() + offset_w;
                    }
                        
                    double pitch = (plane == WirePlaneLayer_t::kUlayer) ? pitch_u :
                                    (plane == WirePlaneLayer_t::kVlayer) ? pitch_v : pitch_w;
                        
                    double factor = exp(-0.5 * (pow((central_t - time) * time_tick_width, 2) + pow((central_ch - wire) * pitch, 2)) / pow(div_sigma, 2));

                    // std::cout << plane << " " << time << " " << wire << " " << idx << " " << central_t << " " << central_ch << " " << factor << std::endl;

                    fac = factor;
                    // sum += factor;

                    auto [it, inserted] = map_tw_sum.try_emplace(std::make_pair(time, wire), factor);
                    if (!inserted) {
                        it->second += factor;
                    }

                }


                
                // Normalize weights
                for (auto& [coord_idx, fac] : coord_idx_fac) {
                    double sum = map_tw_sum[std::make_pair(std::get<0>(coord_idx), std::get<1>(coord_idx))];
                    fac /= sum;
                    // std::cout << plane << " " << std::get<0>(coord_idx) << " " << std::get<1>(coord_idx) << " " << fac << " " << sum << std::endl;
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
        double cluster_t0 = cluster->get_cluster_t0();

        auto plane_data_u = point_info.get_plane_data(WirePlaneLayer_t::kUlayer);
        auto plane_data_v = point_info.get_plane_data(WirePlaneLayer_t::kVlayer);
        auto plane_data_w = point_info.get_plane_data(WirePlaneLayer_t::kWlayer);
        
        int n_2D_u = 2 * plane_data_u.associated_2d_points.size();
        int n_2D_v = 2 * plane_data_v.associated_2d_points.size();
        int n_2D_w = 2 * plane_data_w.associated_2d_points.size();

        // std::cout << i << " " << n_2D_u << " " << n_2D_v << " " << n_2D_w << std::endl;
        
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
            double charge = m_params.default_charge_th, charge_err = m_params.default_charge_err; // Default values

            auto charge_it = m_charge_data.find(charge_key);
            if (charge_it != m_charge_data.end()) {
                charge = charge_it->second.charge;
                charge_err = charge_it->second.charge_err;
            }

            if (charge < m_params.default_charge_th) {
                charge = m_params.default_charge_th;
                charge_err = m_params.default_charge_err;
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
            if (plane_data_u.quantity < m_params.scaling_quality_th) {
                if (plane_data_u.quantity != 0) {
                    scaling *= pow(plane_data_u.quantity / m_params.scaling_quality_th, 1);
                } else {
                    scaling *= m_params.scaling_ratio;
                }
            } 

            WirePlaneId wpid(kAllLayers, it->face, it->apa);
            auto offset_it = wpid_offsets.find(wpid);
            auto slope_it = wpid_slopes.find(wpid);

            auto offset_t = std::get<0>(offset_it->second);
            auto offset_u = std::get<1>(offset_it->second);
            auto slope_x = std::get<0>(slope_it->second);
            auto slope_yu = std::get<1>(slope_it->second).first;
            auto slope_zu = std::get<1>(slope_it->second).second;
               
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
            double charge = m_params.default_charge_th, charge_err = m_params.default_charge_err; // Default values

            auto charge_it = m_charge_data.find(charge_key);
            if (charge_it != m_charge_data.end()) {
                charge = charge_it->second.charge;
                charge_err = charge_it->second.charge_err;
            }

            if (charge < m_params.default_charge_th) {
                charge = m_params.default_charge_th;
                charge_err = m_params.default_charge_err;
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
            if (plane_data_v.quantity < m_params.scaling_quality_th) {
                if (plane_data_v.quantity != 0) {
                    scaling *= pow(plane_data_v.quantity / m_params.scaling_quality_th, 1);
                } else {
                    scaling *= m_params.scaling_ratio;
                }
            } 

            WirePlaneId wpid(kAllLayers, it->face, it->apa);
            auto offset_it = wpid_offsets.find(wpid);
            auto slope_it = wpid_slopes.find(wpid);

            auto offset_t = std::get<0>(offset_it->second);
            auto offset_v = std::get<2>(offset_it->second);
            auto slope_x = std::get<0>(slope_it->second);
            auto slope_yv = std::get<2>(slope_it->second).first;
            auto slope_zv = std::get<2>(slope_it->second).second;
            
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
            double charge = m_params.default_charge_th, charge_err = m_params.default_charge_err; // Default values

            auto charge_it = m_charge_data.find(charge_key);
            if (charge_it != m_charge_data.end()) {
                charge = charge_it->second.charge;
                charge_err = charge_it->second.charge_err;
            }

            if (charge < m_params.default_charge_th) {
                charge = m_params.default_charge_th;
                charge_err = m_params.default_charge_err;
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
            if (plane_data_w.quantity < m_params.scaling_quality_th) {
                if (plane_data_w.quantity != 0) {
                    scaling *= pow(plane_data_w.quantity / m_params.scaling_quality_th, 1);
                } else {
                    scaling *= m_params.scaling_ratio;
                }
            } 

            WirePlaneId wpid(kAllLayers, it->face, it->apa);
            auto offset_it = wpid_offsets.find(wpid);
            auto slope_it = wpid_slopes.find(wpid);

            auto offset_t = std::get<0>(offset_it->second);
            auto offset_w = std::get<3>(offset_it->second);
            auto slope_x = std::get<0>(slope_it->second);
            auto slope_yw = std::get<3>(slope_it->second).first;
            auto slope_zw = std::get<3>(slope_it->second).second;
            
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
        // std::cout << "Track Fitting: " << i << " " << temp_pos_3D(0) << " " << temp_pos_3D(1) << " " << temp_pos_3D(2) << " " << temp_pos_3D_init(0) << " " << temp_pos_3D_init(1) << " " << temp_pos_3D_init(2) << std::endl;
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
        double cluster_t0 = cluster->get_cluster_t0();
        auto test_wpid = m_dv->contained_by(pss_vec[i].first);

        auto p = transform->forward(p_raw, cluster_t0, test_wpid.face(), test_wpid.apa());
        auto apa_face = std::make_pair(test_wpid.face(), test_wpid.apa());
        // all corrected points ...
        bool flag_skip = skip_trajectory_point(p, apa_face, i, pss_vec, fine_tracking_path);

        // std::cout << "Skip: " << i << " " << flag_skip << std::endl;
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

            if (area1 > m_params.area_ratio1 * c && area1 > m_params.area_ratio2 * area2) {
                flag_replace = true;
            }
        }

        //-2, +1
        if ((!flag_replace) && i>=2 && i+1 != fine_tracking_path.size()){
            double a = sqrt(pow(fine_tracking_path.at(i-2).first.x() - fine_tracking_path.at(i).first.x(),2)
                    +pow(fine_tracking_path.at(i-2).first.y() - fine_tracking_path.at(i).first.y(),2)
                    +pow(fine_tracking_path.at(i-2).first.z() - fine_tracking_path.at(i).first.z(),2));
            double b = sqrt(pow(fine_tracking_path.at(i+1).first.x() - fine_tracking_path.at(i).first.x(),2)
                    +pow(fine_tracking_path.at(i+1).first.y() - fine_tracking_path.at(i).first.y(),2)
                    +pow(fine_tracking_path.at(i+1).first.z() - fine_tracking_path.at(i).first.z(),2));
            double c = sqrt(pow(fine_tracking_path.at(i-2).first.x() - fine_tracking_path.at(i+1).first.x(),2)
                    +pow(fine_tracking_path.at(i-2).first.y() - fine_tracking_path.at(i+1).first.y(),2)
                    +pow(fine_tracking_path.at(i-2).first.z() - fine_tracking_path.at(i+1).first.z(),2));
            double s = (a+b+c)/2.;
            double area1 = sqrt(s*(s-a)*(s-b)*(s-c));

            a = sqrt(pow(fine_tracking_path.at(i-2).first.x() - temp_fine_tracking_path.at(i).first.x(),2)
                +pow(fine_tracking_path.at(i-2).first.y() - temp_fine_tracking_path.at(i).first.y(),2)
                +pow(fine_tracking_path.at(i-2).first.z() - temp_fine_tracking_path.at(i).first.z(),2));
            b = sqrt(pow(fine_tracking_path.at(i+1).first.x() - temp_fine_tracking_path.at(i).first.x(),2)
                +pow(fine_tracking_path.at(i+1).first.y() - temp_fine_tracking_path.at(i).first.y(),2)
                +pow(fine_tracking_path.at(i+1).first.z() - temp_fine_tracking_path.at(i).first.z(),2));
            s = (a+b+c)/2.;
            double area2 = sqrt(s*(s-a)*(s-b)*(s-c));
            //std::cout << i << " B " << area1/c << " " << area2/c  << std::endl;
            if (area1 > 1.8*units::mm * c && area1 > 1.7 * area2) flag_replace = true;	
        }
        //-1, +2
        if ((!flag_replace) && i>0 && i+2<fine_tracking_path.size()){
            double a = sqrt(pow(fine_tracking_path.at(i-1).first.x() - fine_tracking_path.at(i).first.x(),2)
                    +pow(fine_tracking_path.at(i-1).first.y() - fine_tracking_path.at(i).first.y(),2)
                    +pow(fine_tracking_path.at(i-1).first.z() - fine_tracking_path.at(i).first.z(),2));
            double b = sqrt(pow(fine_tracking_path.at(i+2).first.x() - fine_tracking_path.at(i).first.x(),2)
                    +pow(fine_tracking_path.at(i+2).first.y() - fine_tracking_path.at(i).first.y(),2)
                    +pow(fine_tracking_path.at(i+2).first.z() - fine_tracking_path.at(i).first.z(),2));
            double c = sqrt(pow(fine_tracking_path.at(i-1).first.x() - fine_tracking_path.at(i+2).first.x(),2)
                    +pow(fine_tracking_path.at(i-1).first.y() - fine_tracking_path.at(i+2).first.y(),2)
                    +pow(fine_tracking_path.at(i-1).first.z() - fine_tracking_path.at(i+2).first.z(),2));
            double s = (a+b+c)/2.;
            double area1 = sqrt(s*(s-a)*(s-b)*(s-c));

            a = sqrt(pow(fine_tracking_path.at(i-1).first.x() - temp_fine_tracking_path.at(i).first.x(),2)
                +pow(fine_tracking_path.at(i-1).first.y() - temp_fine_tracking_path.at(i).first.y(),2)
                +pow(fine_tracking_path.at(i-1).first.z() - temp_fine_tracking_path.at(i).first.z(),2));
            b = sqrt(pow(fine_tracking_path.at(i+2).first.x() - temp_fine_tracking_path.at(i).first.x(),2)
                +pow(fine_tracking_path.at(i+2).first.y() - temp_fine_tracking_path.at(i).first.y(),2)
                +pow(fine_tracking_path.at(i+2).first.z() - temp_fine_tracking_path.at(i).first.z(),2));
            s = (a+b+c)/2.;
            double area2 = sqrt(s*(s-a)*(s-b)*(s-c));

            //      std::cout << i << " C " << area1/c << " " << area2/c  << std::endl;
            
            if (area1 > 1.8*units::mm * c && area1 > 1.7 * area2) flag_replace = true;
        }
        
        
        if (flag_replace) {
            fine_tracking_path[i] = temp_fine_tracking_path[i];
        }

        // std::cout << i << " " << flag_replace << " " << std::endl;
    }
    
    // Generate 2D projections
    pu.clear();
    pv.clear();
    pw.clear();
    pt.clear();
    paf.clear();
    for (size_t i = 0; i < fine_tracking_path.size(); i++) {
        WireCell::Point p = fine_tracking_path[i].first;
        auto segment = fine_tracking_path[i].second;
        auto cluster = segment->cluster();
        const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
        double cluster_t0 = cluster->get_cluster_t0();

        int apa = saved_paf.at(i).first;
        int face = saved_paf.at(i).second;

        auto p_raw = transform->backward(p, cluster_t0, apa, face);
        WirePlaneId wpid(kAllLayers, face, apa);
        auto offset_it = wpid_offsets.find(wpid);
        auto slope_it = wpid_slopes.find(wpid);

        auto offset_t = std::get<0>(offset_it->second);
        auto offset_u = std::get<1>(offset_it->second);
        auto offset_v = std::get<2>(offset_it->second);
        auto offset_w = std::get<3>(offset_it->second);
        auto slope_x = std::get<0>(slope_it->second);
        auto slope_yu = std::get<1>(slope_it->second).first;
        auto slope_zu = std::get<1>(slope_it->second).second;
        auto slope_yv = std::get<2>(slope_it->second).first;
        auto slope_zv = std::get<2>(slope_it->second).second;
        auto slope_yw = std::get<3>(slope_it->second).first;
        auto slope_zw = std::get<3>(slope_it->second).second;

        pu.push_back(offset_u + (slope_yu * p_raw.y() + slope_zu * p_raw.z()));
        pv.push_back(offset_v + (slope_yv * p_raw.y() + slope_zv * p_raw.z()));
        pw.push_back(offset_w + (slope_yw * p_raw.y() + slope_zw * p_raw.z()));
        pt.push_back(offset_t + slope_x * p_raw.x());
        paf.push_back(std::make_pair(apa, face));

    }

        // std::cout << m_params.DL << std::endl;

    
    // Update the input vector with the fitted results
    pss_vec = fine_tracking_path;
}

bool TrackFitting::skip_trajectory_point(WireCell::Point& p, std::pair<int, int>& apa_face, int i, std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>>& pss_vec,  std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>>& fine_tracking_path){
      // Extract APA and face information
    int apa = apa_face.first;
    int face = apa_face.second;
    
    // Get geometry parameters for this APA/face
    WirePlaneId wpid(kAllLayers, face, apa);
    
    auto offset_it = wpid_offsets.find(wpid);
    auto slope_it = wpid_slopes.find(wpid);
    
    if (offset_it == wpid_offsets.end() || slope_it == wpid_slopes.end()) {
        return false; // Can't process without geometry info
    }
    
    // Extract offsets and slopes
    double offset_t = std::get<0>(offset_it->second);
    double offset_u = std::get<1>(offset_it->second);
    double offset_v = std::get<2>(offset_it->second);
    double offset_w = std::get<3>(offset_it->second);
    
    double slope_x = std::get<0>(slope_it->second);
    double slope_yu = std::get<1>(slope_it->second).first;
    double slope_zu = std::get<1>(slope_it->second).second;
    double slope_yv = std::get<2>(slope_it->second).first;
    double slope_zv = std::get<2>(slope_it->second).second;
    double slope_yw = std::get<3>(slope_it->second).first;
    double slope_zw = std::get<3>(slope_it->second).second;

    auto segment = pss_vec.at(i).second;
    auto cluster = segment->cluster();
    const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
    double cluster_t0 = cluster->get_cluster_t0();

     auto first_blob = cluster->children()[0];
    int cur_ntime_ticks = first_blob->slice_index_max() - first_blob->slice_index_min();

    // Initialization with its raw position
    auto p_raw = transform->backward(p, cluster_t0, face, apa);
    // Calculate 2D projections for current point p
    int t1 = std::round((offset_t + slope_x * p_raw.x())/cur_ntime_ticks) * cur_ntime_ticks;  // this needs some fix ... 
    int u1 = std::round(offset_u + (slope_yu * p_raw.y() + slope_zu * p_raw.z()));
    int v1 = std::round(offset_v + (slope_yv * p_raw.y() + slope_zv * p_raw.z()));
    int w1 = std::round(offset_w + (slope_yw * p_raw.y() + slope_zw * p_raw.z()));

    // // test ...
    // auto cur_u = m_grouping->convert_3Dpoint_time_ch(p_raw, 0, 0, 0);
    // auto cur_v = m_grouping->convert_3Dpoint_time_ch(p_raw, 0, 0, 1);
    // auto cur_w = m_grouping->convert_3Dpoint_time_ch(p_raw, 0, 0, 2);
    // // std::cout << t1 << " " << u1 << " " << v1 << " " << w1 << " " << std::get<0>(cur_u) << " " << std::get<1>(cur_u) << " " << std::get<1>(cur_v) << " " << std::get<1>(cur_w) << std::endl;



    // Calculate 2D projections for comparison point pss_vec[i]
    WireCell::Point ps_point = pss_vec.at(i).first;
    auto ps_point_raw = transform->backward(ps_point, cluster_t0, face, apa);
    int t2 = std::round((offset_t + slope_x * ps_point_raw.x())/cur_ntime_ticks)*cur_ntime_ticks; // this needs some fix ...
    int u2 = std::round(offset_u + (slope_yu * ps_point_raw.y() + slope_zu * ps_point_raw.z()));
    int v2 = std::round(offset_v + (slope_yv * ps_point_raw.y() + slope_zv * ps_point_raw.z()));
    int w2 = std::round(offset_w + (slope_yw * ps_point_raw.y() + slope_zw * ps_point_raw.z()));

    // Helper lambda to get charge from nearby coordinates
    auto get_charge_sum = [&](int wire, int time, WirePlaneLayer_t plane) -> double {
        double charge_sum = 0.0;
        
        // std::cout << m_charge_data.size() << std::endl;
        // for (const auto& [coord_key, charge_measurement] : m_charge_data) {
        //     int apa = coord_key.apa;
        //     int time = coord_key.time;
        //     int channel = coord_key.channel;
        //     // Get wires for this channel
        //     std::cout << "apa: " << apa << ", time: " << time << ", channel: " << channel << std::endl;
        // }

        // Convert WirePlaneLayer_t to plane number: kUlayer(1)->0, kVlayer(2)->1, kWlayer(4)->2
        int plane_num = -1;
        if (plane == kUlayer) plane_num = 0;
        else if (plane == kVlayer) plane_num = 1;
        else if (plane == kWlayer) plane_num = 2;
        
        // Search in a 3x3 neighborhood (±1 in wire and time)
        // Only consider the center, wire ±1, and time ±1 (total 5 combinations)
        for (int dw = -1; dw <= 1; dw++) {
            int channel = get_channel_for_wire(apa, face, plane_num, wire + dw);
            if (channel < 0) continue;
            // Center (dt = 0)
            {
            CoordReadout charge_key(apa, time, channel);
            auto charge_it = m_charge_data.find(charge_key);
            if (charge_it != m_charge_data.end() && charge_it->second.flag != 0) {
                charge_sum += charge_it->second.charge;
            }
            }
        }
        // Time -1
        {
            int channel = get_channel_for_wire(apa, face, plane_num, wire);
            if (channel >= 0) {
            CoordReadout charge_key(apa, time - cur_ntime_ticks, channel);
            auto charge_it = m_charge_data.find(charge_key);
            if (charge_it != m_charge_data.end() && charge_it->second.flag != 0) {
                charge_sum += charge_it->second.charge;
            }
            }
        }
        // Time +1
        {
            int channel = get_channel_for_wire(apa, face, plane_num, wire);
            if (channel >= 0) {
            CoordReadout charge_key(apa, time + cur_ntime_ticks, channel);
            auto charge_it = m_charge_data.find(charge_key);
            if (charge_it != m_charge_data.end() && charge_it->second.flag != 0) {
                charge_sum += charge_it->second.charge;
            }
            }
        }
              
        return charge_sum;
    };
    
    // Get charges for point p (c1)
    double c1_u = get_charge_sum(u1, t1, kUlayer);
    double c1_v = get_charge_sum(v1, t1, kVlayer);
    double c1_w = get_charge_sum(w1, t1, kWlayer);
    
    // Get charges for comparison point (c2)
    double c2_u = get_charge_sum(u2, t2, kUlayer);
    double c2_v = get_charge_sum(v2, t2, kVlayer);
    double c2_w = get_charge_sum(w2, t2, kWlayer);
    
    // std::cout << "Skip inside: " << t1 << " " << u1 << " " << v1 << " " << w1 << " | " << t2 << " " << u2 << " " << v2 << " " << w2 << " | " << c1_u << " " << c1_v << " " << c1_w << " " << c2_u << " " << c2_v << " " << c2_w << std::endl;

    // Calculate charge ratios
    double ratio = 0;
    double ratio_1 = 1;
    
    // U plane ratio
    if (c2_u != 0) {
        ratio += c1_u / c2_u;
        if (c1_u != 0)
            ratio_1 *= c1_u / c2_u;
        else
            ratio_1 *= m_params.skip_default_ratio_1;
    } else {
        ratio += 1;
    }
    
    // V plane ratio
    if (c2_v != 0) {
        ratio += c1_v / c2_v;
        if (c1_v != 0)
            ratio_1 *= c1_v / c2_v;
        else
            ratio_1 *= m_params.skip_default_ratio_1;
    } else {
        ratio += 1;
    }
    
    // W plane ratio
    if (c2_w != 0) {
        ratio += c1_w / c2_w;
        if (c1_w != 0)
            ratio_1 *= c1_w / c2_w;
        else
            ratio_1 *= m_params.skip_default_ratio_1;
    } else {
        ratio += 1;
    }

    // std::cout << "Inside: " << ratio << " " << ratio_1 << std::endl;
    
    // Apply charge-based correction
    if (ratio / 3.0 < m_params.skip_ratio_cut || ratio_1 < m_params.skip_ratio_1_cut) {
        p = ps_point;
    }
    
    // Angle constraint checking
    if (fine_tracking_path.size() >= 2) {
        // Get direction vectors for angle calculation
        auto& last_point = fine_tracking_path[fine_tracking_path.size()-1].first;
        auto& second_last_point = fine_tracking_path[fine_tracking_path.size()-2].first;
        
        // Vector from second-to-last to last point in fine tracking path
        WireCell::Vector v1(last_point.x() - second_last_point.x(),
                           last_point.y() - second_last_point.y(),
                           last_point.z() - second_last_point.z());
        
        // Vector from last point to current point p
        WireCell::Vector v2(p.x() - last_point.x(),
                           p.y() - last_point.y(),
                           p.z() - last_point.z());

        // std::cout << ratio << " " << ratio_1 << " (" << p.x() << " " << p.y() << " " <<p.z() <<") (" << last_point.x() << " " << last_point.y() << " " << last_point.z() << ") (" << second_last_point.x() << " " << second_last_point.y() << " " << second_last_point.z() << ")" << std::endl;

        // Calculate angle between vectors (in degrees)
        double dot_product = v1.dot(v2);
        double mag1 = sqrt(v1.x()*v1.x() + v1.y()*v1.y() + v1.z()*v1.z());
        double mag2 = sqrt(v2.x()*v2.x() + v2.y()*v2.y() + v2.z()*v2.z());
        
        double angle = 180.0;  // Default to large angle if vectors are too small
        if (mag1 > 0 && mag2 > 0) {
            double cos_angle = dot_product / (mag1 * mag2);
            // Clamp to [-1, 1] to handle numerical errors
            cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
            angle = acos(cos_angle) * 180.0 / M_PI;
        }
        
        // Calculate angle between consecutive segments in original path for comparison
        double angle1 = 180.0;
        if (i >= 2) {
            auto& p_i = pss_vec[i].first;
            auto& p_i1 = pss_vec[i-1].first;
            auto& p_i2 = pss_vec[i-2].first;
            
            WireCell::Vector v3(p_i1.x() - p_i2.x(), p_i1.y() - p_i2.y(), p_i1.z() - p_i2.z());
            WireCell::Vector v4(p_i.x() - p_i1.x(), p_i.y() - p_i1.y(), p_i.z() - p_i1.z());
            
            double dot_34 = v3.dot(v4);
            double mag3 = sqrt(v3.x()*v3.x() + v3.y()*v3.y() + v3.z()*v3.z());
            double mag4 = sqrt(v4.x()*v4.x() + v4.y()*v4.y() + v4.z()*v4.z());
            
            if (mag3 > 0 && mag4 > 0) {
                double cos_angle1 = dot_34 / (mag3 * mag4);
                cos_angle1 = std::max(-1.0, std::min(1.0, cos_angle1));
                angle1 = acos(cos_angle1) * 180.0 / M_PI;
            }
        }
        
        // Get hit information for dead channel detection
        // Check if we have valid 3D to 2D mapping for current point index
        bool has_u_hits = false, has_v_hits = false, has_w_hits = false;
        if (m_3d_to_2d.find(i) != m_3d_to_2d.end()) {
            const auto& point_info = m_3d_to_2d.at(i);
            has_u_hits = point_info.get_plane_data(kUlayer).quantity > 0;
            has_v_hits = point_info.get_plane_data(kVlayer).quantity > 0;
            has_w_hits = point_info.get_plane_data(kWlayer).quantity > 0;
        }
        
        // Check for dead channel conditions
        int dead_plane_count = 0;
        if (!has_u_hits) dead_plane_count++;
        if (!has_v_hits) dead_plane_count++;
        if (!has_w_hits) dead_plane_count++;

        // std::cout << "Inside: " << angle << " " << angle1 << " " << dead_plane_count << " " << mag2 << std::endl;
        
        if (angle > m_params.skip_angle_cut_3 && dead_plane_count >= 2) {
            return true;
        }
        
        // Check for fold-back or extreme angles
        if (angle > m_params.skip_angle_cut_1 || angle > angle1 + m_params.skip_angle_cut_2) {
            return true;
        }
        
        // Check last point protection
        if (i + 1 == pss_vec.size() && angle > m_params.skip_angle_cut_3 && mag2 < m_params.skip_dis_cut) {
            return true;
        }
    }


    
    return false;

 }


 double TrackFitting::cal_gaus_integral(int tbin, int wbin, double t_center, double t_sigma, 
                                       double w_center, double w_sigma, int flag, double nsigma, int cur_ntime_ticks) {
    // flag = 0: no boundary effect, pure Gaussian, time or collection plane
    // flag = 1: taking into account boundary effect for induction plane
    // flag = 2: more complex induction plane response
    
    double result = 0;
    
    // *** COORDINATE SYSTEM CLARIFICATION ***
    // In this toolkit convention:
    // - w_center = offset_u + (slope_yu * p.y + slope_zu * p.z)  [continuous coordinate]
    // - wbin = std::round(w_center)  [bin index - nearest integer]
    // - t_center = offset_t + slope_t * p.x  [continuous coordinate] 
    // - tbin = std::round(t_center)  [bin index - nearest integer]
    // Therefore: compare tbin vs t_center and wbin vs w_center DIRECTLY
    
    // Check if we're within nsigma range of both time and wire centers
    if (fabs(tbin - t_center) <= nsigma * t_sigma &&              // Direct comparison: bin index vs continuous center
        fabs(wbin - w_center) <= nsigma * w_sigma) {              // Direct comparison: bin index vs continuous center
        
        // Time dimension integration 
        // If tbin = std::round(t_center), then bin spans [tbin-0.5, tbin+0.5]
        result = 0.5 * (std::erf((tbin + 0.5*cur_ntime_ticks - t_center) / sqrt(2.) / t_sigma) - 
                       std::erf((tbin - 0.5*cur_ntime_ticks - t_center) / sqrt(2.) / t_sigma));
        
        if (flag == 0) {
            // Pure Gaussian case - simple wire dimension integration
            // If wbin = std::round(w_center), then bin spans [wbin-0.5, wbin+0.5]
            result *= 0.5 * (std::erf((wbin + 0.5 - w_center) / sqrt(2.) / w_sigma) - 
                            std::erf((wbin - 0.5 - w_center) / sqrt(2.) / w_sigma));

        // std::cout << tbin << " " << t_center << " " << t_sigma << " " << 0.5 * (std::erf((tbin + 0.5*cur_ntime_ticks - t_center) / sqrt(2.) / t_sigma) - 
        //                std::erf((tbin - 0.5*cur_ntime_ticks - t_center) / sqrt(2.) / t_sigma)) << " | " <<  0.5 * (std::erf((wbin + 0.5 - w_center) / sqrt(2.) / w_sigma) - 
        //                     std::erf((wbin - 0.5 - w_center) / sqrt(2.) / w_sigma)) << std::endl;

        } else if (flag == 1) {
            // Induction plane with bipolar response
            // All boundaries shift by -0.5 due to bin convention change
            
            // First part: positive lobe (was wbin+0.5 to wbin+1.5, now wbin+0.0 to wbin+1.0)
            double x2 = wbin + 1.0;     // was wbin + 1.5, shifted by -0.5
            double x1 = wbin + 0.0;     // was wbin + 0.5, shifted by -0.5 (bin center)
            double x0 = w_center;
            
            double content1 = 0.5 * (std::erf((x2 - x0) / sqrt(2.) / w_sigma) - 
                                    std::erf((x1 - x0) / sqrt(2.) / w_sigma));
            
            // Weight calculation for positive lobe
            double w1 = -pow(w_sigma, 2) / (-1) / sqrt(2. * 3.1415926) / w_sigma *
                       (exp(-pow(x0 - x2, 2) / 2. / pow(w_sigma, 2)) - 
                        exp(-pow(x0 - x1, 2) / 2. / pow(w_sigma, 2))) /
                       (0.5 * std::erf((x2 - x0) / sqrt(2.) / w_sigma) - 
                        0.5 * std::erf((x1 - x0) / sqrt(2.) / w_sigma)) +
                       (x0 - x2) / (-1);
            
            // Second part: negative lobe (was wbin-0.5 to wbin+0.5, now wbin-1.0 to wbin+0.0)
            x2 = wbin + 0.0;            // was wbin + 0.5, shifted by -0.5 (bin center)
            x1 = wbin - 1.0;            // was wbin - 0.5, shifted by -0.5
            
            double content2 = 0.5 * (std::erf((x2 - x0) / sqrt(2.) / w_sigma) - 
                                    std::erf((x1 - x0) / sqrt(2.) / w_sigma));
            
            // Weight calculation for negative lobe
            double w2 = -pow(w_sigma, 2) / (-1) / sqrt(2. * 3.1415926) / w_sigma *
                       (exp(-pow(x0 - x2, 2) / 2. / pow(w_sigma, 2)) - 
                        exp(-pow(x0 - x1, 2) / 2. / pow(w_sigma, 2))) /
                       (0.5 * std::erf((x2 - x0) / sqrt(2.) / w_sigma) - 
                        0.5 * std::erf((x1 - x0) / sqrt(2.) / w_sigma)) +
                       (x0 - x2) / (-1);
            
            // Combine positive and negative contributions
            result *= (content1 * w1 + content2 * (1 - w2));
            
        } else if (flag == 2) {
            // More complex induction response with multiple components
            // All boundaries shift by -0.5 due to bin convention change
            double sum = 0;
            
            // Component 1: (was wbin+0.5 to wbin+1.0, now wbin+0.0 to wbin+0.5)
            double x2 = wbin + 0.5;     // was wbin + 1.0, shifted by -0.5
            double x1 = wbin + 0.0;     // was wbin + 0.5, shifted by -0.5 (bin center)
            double x0 = w_center;
            
            double content1 = 0.5 * (std::erf((x2 - x0) / sqrt(2.) / w_sigma) - 
                                    std::erf((x1 - x0) / sqrt(2.) / w_sigma));
            double w1 = -pow(w_sigma, 2) / (-1) / sqrt(2. * 3.1415926) / w_sigma *
                       (exp(-pow(x0 - x2, 2) / 2. / pow(w_sigma, 2)) - 
                        exp(-pow(x0 - x1, 2) / 2. / pow(w_sigma, 2))) /
                       (0.5 * std::erf((x2 - x0) / sqrt(2.) / w_sigma) - 
                        0.5 * std::erf((x1 - x0) / sqrt(2.) / w_sigma)) +
                       (x0 - x2) / (-1);
            
            sum += content1 * (0.545 + 0.697 * w1);
            
            // Component 2: (was wbin+1.0 to wbin+1.5, now wbin+0.5 to wbin+1.0)
            x2 = wbin + 1.0;            // was wbin + 1.5, shifted by -0.5
            x1 = wbin + 0.5;            // was wbin + 1.0, shifted by -0.5
            
            content1 = 0.5 * (std::erf((x2 - x0) / sqrt(2.) / w_sigma) - 
                             std::erf((x1 - x0) / sqrt(2.) / w_sigma));
            w1 = -pow(w_sigma, 2) / (-1) / sqrt(2. * 3.1415926) / w_sigma *
                (exp(-pow(x0 - x2, 2) / 2. / pow(w_sigma, 2)) - 
                 exp(-pow(x0 - x1, 2) / 2. / pow(w_sigma, 2))) /
                (0.5 * std::erf((x2 - x0) / sqrt(2.) / w_sigma) - 
                 0.5 * std::erf((x1 - x0) / sqrt(2.) / w_sigma)) +
                (x0 - x2) / (-1);
            
            sum += content1 * (0.11364 + 0.1 * w1);
            
            // Component 3: (was wbin+0.0 to wbin+0.5, now wbin-0.5 to wbin+0.0)
            x2 = wbin + 0.0;            // was wbin + 0.5, shifted by -0.5 (bin center)
            x1 = wbin - 0.5;            // was wbin + 0.0, shifted by -0.5
            
            content1 = 0.5 * (std::erf((x2 - x0) / sqrt(2.) / w_sigma) - 
                             std::erf((x1 - x0) / sqrt(2.) / w_sigma));
            w1 = -pow(w_sigma, 2) / (-1) / sqrt(2. * 3.1415926) / w_sigma *
                (exp(-pow(x0 - x2, 2) / 2. / pow(w_sigma, 2)) - 
                 exp(-pow(x0 - x1, 2) / 2. / pow(w_sigma, 2))) /
                (0.5 * std::erf((x2 - x0) / sqrt(2.) / w_sigma) - 
                 0.5 * std::erf((x1 - x0) / sqrt(2.) / w_sigma)) +
                (x0 - x2) / (-1);
            
            sum += content1 * (0.545 + 0.697 * (1 - w1));
            
            // Component 4: (was wbin-0.5 to wbin+0.0, now wbin-1.0 to wbin-0.5)
            x2 = wbin - 0.5;            // was wbin + 0.0, shifted by -0.5
            x1 = wbin - 1.0;            // was wbin - 0.5, shifted by -0.5
            
            content1 = 0.5 * (std::erf((x2 - x0) / sqrt(2.) / w_sigma) - 
                             std::erf((x1 - x0) / sqrt(2.) / w_sigma));
            w1 = -pow(w_sigma, 2) / (-1) / sqrt(2. * 3.1415926) / w_sigma *
                (exp(-pow(x0 - x2, 2) / 2. / pow(w_sigma, 2)) - 
                 exp(-pow(x0 - x1, 2) / 2. / pow(w_sigma, 2))) /
                (0.5 * std::erf((x2 - x0) / sqrt(2.) / w_sigma) - 
                 0.5 * std::erf((x1 - x0) / sqrt(2.) / w_sigma)) +
                (x0 - x2) / (-1);
            
            sum += content1 * (0.11364 + 0.1 * (1 - w1));
            
            result *= sum;
        }
    }
    
    return result;
}


double TrackFitting::cal_gaus_integral_seg(int tbin, int wbin, std::vector<double>& t_centers, std::vector<double>& t_sigmas, std::vector<double>& w_centers, std::vector<double>& w_sigmas, std::vector<double>& weights, int flag, double nsigma, int cur_ntime_ticks){
  double result = 0;
  double result1 = 0;

  for (size_t i=0;i!=t_centers.size();i++){
    result += cal_gaus_integral(tbin,wbin,t_centers.at(i), t_sigmas.at(i), w_centers.at(i), w_sigmas.at(i),flag,nsigma,cur_ntime_ticks) * weights.at(i);
    result1 += weights.at(i);

    // std::cout << cal_gaus_integral(tbin,wbin,t_centers.at(i), t_sigmas.at(i), w_centers.at(i), w_sigmas.at(i),flag,nsigma,cur_ntime_ticks) << " " << weights.at(i) << std::endl;
  }

  result /= result1;
  
  return result;
}


void TrackFitting::update_dQ_dx_data() {
    // Step 1: Loop over m_clusters to collect all track blobs
    std::set<Facade::Blob*> track_blobs_set;
    for (auto cluster : m_clusters) {
        // Collect blobs from each cluster using toolkit convention
        for (auto blob : cluster->children()) {
            track_blobs_set.insert(blob);
        }
    }
    
    // Step 2: Check each measurement in global_rb_map
    for (const auto& [coord_key, blob_set] : global_rb_map) {
        // coord_key is of type CoordReadout
        // blob_set is of type std::set<const Facade::Blob*>

        bool is_shared = false;

        for (auto blob : blob_set) {
            if (track_blobs_set.find(blob) == track_blobs_set.end()) {
                // Found a blob not belonging to our track clusters
                is_shared = true;
                break;
            }
        }
            
        if (is_shared) {
            // Find and modify the measurement if it exists
            auto charge_it = m_charge_data.find(coord_key);
            if (charge_it != m_charge_data.end()) {
                // Get current measurement and increase charge error for shared measurements
                ChargeMeasurement& measurement = charge_it->second;
                m_orig_charge_data[coord_key] = measurement;
                measurement.charge_err = m_params.share_charge_err; // High penalty for shared wires
            }
        }
    }
}

void TrackFitting::recover_original_charge_data(){
    for (const auto& [coord_key, measurement] : m_orig_charge_data) {
        m_charge_data[coord_key] = measurement;
    }
    m_orig_charge_data.clear();
}

std::vector<std::pair<double, double>> TrackFitting::calculate_compact_matrix(
    Eigen::SparseMatrix<double>& weight_matrix,
    const Eigen::SparseMatrix<double>& response_matrix_transpose,
    int n_2d_measurements,
    int n_3d_positions,
    double cut_position){
    std::vector<std::pair<double,double> > results(n_3d_positions, std::make_pair(0,0));

    // Initialize count vector for 2D measurements
    std::vector<int> count_2d(n_2d_measurements, 1);
    
    // Maps for storing relationships between 2D and 3D indices
    std::map<int, std::set<int>> map_2d_to_3d;
    std::map<int, std::set<int>> map_3d_to_2d;
    std::map<std::pair<int, int>, double> map_pair_values;
    
    // Build mapping structures by iterating through sparse matrix
    for (int k = 0; k < response_matrix_transpose.outerSize(); ++k) {
        int count = 0;
        
        for (Eigen::SparseMatrix<double>::InnerIterator it(response_matrix_transpose, k); it; ++it) {
            int row = it.row();
            int col = it.col();
            double value = it.value();
            
            // Build 2D to 3D mapping
            if (map_2d_to_3d.find(col) != map_2d_to_3d.end()) {
                map_2d_to_3d[col].insert(row);
            } else {
                std::set<int> temp_set;
                temp_set.insert(row);
                map_2d_to_3d[col] = temp_set;
            }
            
            // Build 3D to 2D mapping  
            if (map_3d_to_2d.find(row) != map_3d_to_2d.end()) {
                map_3d_to_2d[row].insert(col);
            } else {
                std::set<int> temp_set;
                temp_set.insert(col);
                map_3d_to_2d[row] = temp_set;
            }
            
            // Store pair values for later lookup
            map_pair_values[std::make_pair(row, col)] = value;
            count++;
        }
        
        count_2d.at(k) = count;
    }
    
    // Calculate average count for 3D positions
    std::vector<std::pair<double, int>> average_count(n_3d_positions);
    for (auto it = map_3d_to_2d.begin(); it != map_3d_to_2d.end(); ++it) {
        int row = it->first;
        double sum1 = 0.0;
        double sum2 = 0.0;
        int flag = 0;
        
        for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1) {
            int col = *it1;
            double val = map_pair_values[std::make_pair(row, col)];
            sum1 += count_2d[col] * val;
            sum2 += val;
            if (count_2d[col] > 2) {
                flag = 1;
            }
        }
        
        average_count.at(row) = std::make_pair(sum1 / sum2, flag);
    }
    
    // Update 2D measurement weights based on 3D position sharing
    for (auto it = map_2d_to_3d.begin(); it != map_2d_to_3d.end(); ++it) {
        int col = it->first;
        double sum1 = 0.0;
        double sum2 = 0.0;
        int flag = 0;
        
        for (auto it1 = it->second.begin(); it1 != it->second.end(); ++it1) {
            int row = *it1;
            double val = map_pair_values[std::make_pair(row, col)];
            if (average_count.at(row).second == 1) {
                flag = 1;
            }
            sum1 += average_count.at(row).first * val;
            sum2 += val;
        }
        
        // Adjust weight matrix coefficients based on sharing criteria
        if (flag == 1 && weight_matrix.coeffRef(col, col) == 1 && sum1 > cut_position * sum2) {
            weight_matrix.coeffRef(col, col) = std::pow(1.0 / (sum1 / sum2 - cut_position + 1), 2);
        }
    }
    
    // Calculate sharing ratios between neighboring 3D positions
    
    for (auto it = map_3d_to_2d.begin(); it != map_3d_to_2d.end(); ++it) {
        int row = it->first;
        auto it_prev = map_3d_to_2d.find(row - 1);
        auto it_next = map_3d_to_2d.find(row + 1);
        
        double sum[3] = {0.0, 0.0, 0.0};
        
        // Count total connections for current 3D position
        for (auto it3 = it->second.begin(); it3 != it->second.end(); ++it3) {
            sum[0] += 1.0;  // Total count
        }
        
        // Count shared connections with previous neighbor
        if (it_prev != map_3d_to_2d.end()) {
            std::vector<int> common_results(it->second.size());
            auto it3 = std::set_intersection(
                it->second.begin(), it->second.end(),
                it_prev->second.begin(), it_prev->second.end(),
                common_results.begin()
            );
            common_results.resize(it3 - common_results.begin());
            
            for (auto it4 = common_results.begin(); it4 != common_results.end(); ++it4) {
                sum[1] += 1.0;  // Shared with previous
            }
        }
        
        // Count shared connections with next neighbor
        if (it_next != map_3d_to_2d.end()) {
            std::vector<int> common_results(it->second.size());
            auto it3 = std::set_intersection(
                it->second.begin(), it->second.end(),
                it_next->second.begin(), it_next->second.end(),
                common_results.begin()
            );
            common_results.resize(it3 - common_results.begin());
            
            for (auto it4 = common_results.begin(); it4 != common_results.end(); ++it4) {
                sum[2] += 1.0;  // Shared with next
            }
        }
        
        // Calculate overlap ratios
        if (sum[0] > 0) {
            results.at(row).first = sum[1] / sum[0];   // Previous neighbor ratio
            results.at(row).second = sum[2] / sum[0];  // Next neighbor ratio
        }
    }

    return results;
}

void TrackFitting::dQ_dx_fill(double dis_end_point_ext) {
    if (fine_tracking_path.size() <= 1) return;
    
    // Resize vectors to match fine_tracking_path size
    dQ.resize(fine_tracking_path.size(), 0);
    dx.resize(fine_tracking_path.size(), 0);
    reduced_chi2.resize(fine_tracking_path.size(), 0);
    
    // Loop through each point in the fine tracking path
    for (size_t i = 0; i != fine_tracking_path.size(); i++) {
        WireCell::Point curr_rec_pos = fine_tracking_path.at(i).first;
        WireCell::Point prev_rec_pos, next_rec_pos;
        
        if (i == 0) {
            // First point: extrapolate backward from the direction to next point
            next_rec_pos = WireCell::Point(
                (fine_tracking_path.at(i).first.x() + fine_tracking_path.at(i+1).first.x()) / 2.0,
                (fine_tracking_path.at(i).first.y() + fine_tracking_path.at(i+1).first.y()) / 2.0,
                (fine_tracking_path.at(i).first.z() + fine_tracking_path.at(i+1).first.z()) / 2.0
            );
            
            double length = sqrt(
                pow(fine_tracking_path.at(i+1).first.x() - fine_tracking_path.at(i).first.x(), 2) +
                pow(fine_tracking_path.at(i+1).first.y() - fine_tracking_path.at(i).first.y(), 2) +
                pow(fine_tracking_path.at(i+1).first.z() - fine_tracking_path.at(i).first.z(), 2)
            );
            
            if (length == 0) {
                prev_rec_pos = fine_tracking_path.at(i).first;
            } else {
                prev_rec_pos = WireCell::Point(
                    fine_tracking_path.at(i).first.x() - (fine_tracking_path.at(i+1).first.x() - fine_tracking_path.at(i).first.x()) / length * dis_end_point_ext,
                    fine_tracking_path.at(i).first.y() - (fine_tracking_path.at(i+1).first.y() - fine_tracking_path.at(i).first.y()) / length * dis_end_point_ext,
                    fine_tracking_path.at(i).first.z() - (fine_tracking_path.at(i+1).first.z() - fine_tracking_path.at(i).first.z()) / length * dis_end_point_ext
                );
            }
        } else if (i + 1 == fine_tracking_path.size()) {
            // Last point: extrapolate forward from the direction from previous point
            prev_rec_pos = WireCell::Point(
                (fine_tracking_path.at(i).first.x() + fine_tracking_path.at(i-1).first.x()) / 2.0,
                (fine_tracking_path.at(i).first.y() + fine_tracking_path.at(i-1).first.y()) / 2.0,
                (fine_tracking_path.at(i).first.z() + fine_tracking_path.at(i-1).first.z()) / 2.0
            );
            
            double length = sqrt(
                pow(fine_tracking_path.at(i-1).first.x() - fine_tracking_path.at(i).first.x(), 2) +
                pow(fine_tracking_path.at(i-1).first.y() - fine_tracking_path.at(i).first.y(), 2) +
                pow(fine_tracking_path.at(i-1).first.z() - fine_tracking_path.at(i).first.z(), 2)
            );
            
            if (length == 0) {
                next_rec_pos = fine_tracking_path.at(i).first;
            } else {
                next_rec_pos = WireCell::Point(
                    fine_tracking_path.at(i).first.x() - (fine_tracking_path.at(i-1).first.x() - fine_tracking_path.at(i).first.x()) / length * dis_end_point_ext,
                    fine_tracking_path.at(i).first.y() - (fine_tracking_path.at(i-1).first.y() - fine_tracking_path.at(i).first.y()) / length * dis_end_point_ext,
                    fine_tracking_path.at(i).first.z() - (fine_tracking_path.at(i-1).first.z() - fine_tracking_path.at(i).first.z()) / length * dis_end_point_ext
                );
            }
        } else {
            // Middle points: use midpoints to neighboring points
            prev_rec_pos = WireCell::Point(
                (fine_tracking_path.at(i).first.x() + fine_tracking_path.at(i-1).first.x()) / 2.0,
                (fine_tracking_path.at(i).first.y() + fine_tracking_path.at(i-1).first.y()) / 2.0,
                (fine_tracking_path.at(i).first.z() + fine_tracking_path.at(i-1).first.z()) / 2.0
            );
            
            next_rec_pos = WireCell::Point(
                (fine_tracking_path.at(i).first.x() + fine_tracking_path.at(i+1).first.x()) / 2.0,
                (fine_tracking_path.at(i).first.y() + fine_tracking_path.at(i+1).first.y()) / 2.0,
                (fine_tracking_path.at(i).first.z() + fine_tracking_path.at(i+1).first.z()) / 2.0
            );
        }
        
        // Calculate dx as sum of distances to previous and next positions
        dx.at(i) = sqrt(
            pow(curr_rec_pos.x() - prev_rec_pos.x(), 2) +
            pow(curr_rec_pos.y() - prev_rec_pos.y(), 2) +
            pow(curr_rec_pos.z() - prev_rec_pos.z(), 2)
        ) + sqrt(
            pow(curr_rec_pos.x() - next_rec_pos.x(), 2) +
            pow(curr_rec_pos.y() - next_rec_pos.y(), 2) +
            pow(curr_rec_pos.z() - next_rec_pos.z(), 2)
        );
        
        // Set placeholder dQ value (5000 * dx as in original)
        dQ.at(i) = m_params.default_dQ_dx * dx.at(i);
        
        // Initialize reduced_chi2 to 0
        reduced_chi2.at(i) = 0;
    }

}

void WireCell::Clus::TrackFitting::dQ_dx_fit(double dis_end_point_ext, bool flag_dQ_dx_fit_reg) {
    if (fine_tracking_path.size() <= 1) return;
    
    // Clear output vectors
    dQ.clear();
    dx.clear();
    reduced_chi2.clear();
    
    // Update charge data for shared wires (uses existing toolkit function)
    update_dQ_dx_data();
    
    const double DL = m_params.DL;                    // WAS: const double DL = 6.4e-7;
    const double DT = m_params.DT;                    // WAS: const double DT = 9.8e-7;
    const double col_sigma_w_T = m_params.col_sigma_w_T;  // WAS: const double col_sigma_w_T = 0.188060 * 0.2;
    const double ind_sigma_u_T = m_params.ind_sigma_u_T;  // WAS: const double ind_sigma_u_T = 0.402993 * 0.3;
    const double ind_sigma_v_T = m_params.ind_sigma_v_T;  // WAS: const double ind_sigma_v_T = 0.402993 * 0.5;
    const double rel_uncer_ind = m_params.rel_uncer_ind; // WAS: const double rel_uncer_ind = 0.075;
    const double rel_uncer_col = m_params.rel_uncer_col; // WAS: const double rel_uncer_col = 0.05;
    const double add_uncer_ind = m_params.add_uncer_ind; // WAS: const double add_uncer_ind = 0.0;
    const double add_uncer_col = m_params.add_uncer_col; // WAS: const double add_uncer_col = 300.0;
    const double add_sigma_L = m_params.add_sigma_L;     // WAS: const double add_sigma_L = 1.428249 * 0.5;
    
    std::map<CoordReadout, std::pair<ChargeMeasurement, std::set<Coord2D>>> map_U_charge_2D, map_V_charge_2D, map_W_charge_2D;
    // Fill the maps from m_charge_data
    // Fill the maps from m_charge_data
    for (const auto& [coord_readout, charge_measurement] : m_charge_data) {
        int apa = coord_readout.apa;
        int time = coord_readout.time;
        int channel = coord_readout.channel;
        
        // Get wires for this channel using the dedicated function
        auto wires_info = get_wires_for_channel(apa, channel);
        if (wires_info.empty()) continue; // Skip if no wire mapping found
        
        std::set<TrackFitting::Coord2D> associated_coords;
        int plane = -1; // asssuming all wires are from the same plane name ...
        // Process each wire associated with this channel
        for (const auto& wire_info : wires_info) {
            int face = std::get<0>(wire_info);
            plane = std::get<1>(wire_info);
            int wire = std::get<2>(wire_info);
            
            // Convert plane int to WirePlaneLayer_t
            WirePlaneLayer_t plane_layer = (plane == 0) ? kUlayer : 
                                        (plane == 1) ? kVlayer : kWlayer;
            
            // Create TrackFitting::Coord2D with all fields filled
            TrackFitting::Coord2D coord_2d(apa, face, time, wire, channel, plane_layer);
            associated_coords.insert(coord_2d);
        }

        // Create the pair for storage
        std::pair<ChargeMeasurement, std::set<TrackFitting::Coord2D>> charge_coord_pair = std::make_pair(charge_measurement, associated_coords);

        // Store in appropriate plane map
        switch (plane) {
            case 0: // U plane
                map_U_charge_2D[coord_readout] = charge_coord_pair;
                break;
            case 1: // V plane  
                map_V_charge_2D[coord_readout] = charge_coord_pair;
                break;
            case 2: // W plane
                map_W_charge_2D[coord_readout] = charge_coord_pair;
                break;
        }
    }


    std::cout << "dQ/dx: " << map_U_charge_2D.size() << " " << map_V_charge_2D.size() << " " << map_W_charge_2D.size() << std::endl;
    // for (const auto& [coord_key, result] : map_U_charge_2D) {
    //     std::cout << "CoordReadout: APA=" << coord_key.apa
    //               << ", Time=" << coord_key.time
    //               << ", Channel=" << coord_key.channel << std::endl;
    //     const auto& measurement = result.first;
    //     std::cout << "  Charge: " << measurement.charge
    //               << ", ChargeErr: " << measurement.charge_err
    //               << ", Flag: " << measurement.flag << std::endl;
    //     std::cout << "  Associated Coord2D set size: " << result.second.size() << std::endl;
    //     for (const auto& coord2d : result.second) {
    //         std::cout << "    Coord2D: APA=" << coord2d.apa
    //                   << ", Face=" << coord2d.face
    //                   << ", Time=" << coord2d.time
    //                   << ", Wire=" << coord2d.wire
    //                   << ", Channel=" << coord2d.channel
    //                   << ", Plane=" << coord2d.plane << std::endl;
    //     }
    // }


    int n_3D_pos = fine_tracking_path.size();
    // need to separate measurements into U, V, W and form separate matrices ... 
    // need to store measurement --> U, V, W --> measurements
    int n_2D_u = map_U_charge_2D.size();
    int n_2D_v = map_V_charge_2D.size();
    int n_2D_w = map_W_charge_2D.size();
    
    if (n_2D_u == 0 && n_2D_v == 0 && n_2D_w == 0) return;
    
    // // Initialize Eigen matrices and vectors
    Eigen::VectorXd pos_3D(n_3D_pos), data_u_2D(n_2D_u), data_v_2D(n_2D_v), data_w_2D(n_2D_w);
    Eigen::VectorXd pred_data_u_2D(n_2D_u), pred_data_v_2D(n_2D_v), pred_data_w_2D(n_2D_w);
    
    Eigen::SparseMatrix<double> RU(n_2D_u, n_3D_pos);
    Eigen::SparseMatrix<double> RV(n_2D_v, n_3D_pos);
    Eigen::SparseMatrix<double> RW(n_2D_w, n_3D_pos);
    
    Eigen::VectorXd pos_3D_init(n_3D_pos);
    std::vector<int> reg_flag_u(n_3D_pos, 0), reg_flag_v(n_3D_pos, 0), reg_flag_w(n_3D_pos, 0);
    
    
    // Initialize solution vector
    for (int i = 0; i < n_3D_pos; i++) {
        pos_3D_init(i) = 50000.0; // Initial guess
    }
    
    // Fill data vectors with charge/uncertainty ratios
    {
        int n_u = 0;
        for (const auto& [coord_key, result] : map_U_charge_2D) {
            const auto& measurement = result.first;
            if (measurement.charge >0) {
                double charge = measurement.charge;
                double charge_err = measurement.charge_err;
                double total_err = sqrt(pow(charge_err, 2) + pow(charge * rel_uncer_ind, 2) + pow(add_uncer_ind, 2));
                data_u_2D(n_u) = charge / total_err;
            } else {
                data_u_2D(n_u) = 0;
            }
            // std::cout << coord_key.time << " " << coord_key.channel << " " << measurement.charge << " " << measurement.charge_err << " " << data_u_2D(n_u) << std::endl;
            n_u++;
        }
        int n_v = 0;
        for (const auto& [coord_key, result] : map_V_charge_2D) {
            const auto& measurement = result.first;
            if (measurement.charge >0) {
                double charge = measurement.charge;
                double charge_err = measurement.charge_err;
                double total_err = sqrt(pow(charge_err, 2) + pow(charge * rel_uncer_ind, 2) + pow(add_uncer_ind, 2));
                data_v_2D(n_v) = charge / total_err;
            } else {
                data_v_2D(n_v) = 0;
            }
            n_v++;
        }
        int n_w = 0;
        for (const auto& [coord_key, result] : map_W_charge_2D) {
            const auto& measurement = result.first;
            if (measurement.charge >0) {
                double charge = measurement.charge;
                double charge_err = measurement.charge_err;
                double total_err = sqrt(pow(charge_err, 2) + pow(charge * rel_uncer_col, 2) + pow(add_uncer_col, 2));
                data_w_2D(n_w) = charge / total_err;
            } else {
                data_w_2D(n_w) = 0;
            }
            n_w++;
        }
    }
    
    // Calculate dx values (path segment lengths)
    dx.resize(n_3D_pos);
    for (int i = 0; i < n_3D_pos; i++) {
        WireCell::Point prev_rec_pos, next_rec_pos;
        WireCell::Point curr_rec_pos = fine_tracking_path.at(i).first;
        
        if (i == 0) {
            // First point: extrapolate backward
            if (n_3D_pos > 1) {
                WireCell::Point next_point = fine_tracking_path.at(i+1).first;
                WireCell::Vector dir = next_point - curr_rec_pos;
                double length = dir.magnitude();
                if (length > 0) {
                    prev_rec_pos = curr_rec_pos - (dir / length) * dis_end_point_ext;
                } else {
                    prev_rec_pos = curr_rec_pos;
                }
                next_rec_pos = (curr_rec_pos + next_point) * 0.5;
            } else {
                prev_rec_pos = curr_rec_pos;
                next_rec_pos = curr_rec_pos;
            }
        } else if (i == n_3D_pos - 1) {
            // Last point: extrapolate forward
            WireCell::Point prev_point = fine_tracking_path.at(i-1).first;
            WireCell::Vector dir = curr_rec_pos - prev_point;
            double length = dir.magnitude();
            if (length > 0) {
                next_rec_pos = curr_rec_pos + (dir / length) * dis_end_point_ext;
            } else {
                next_rec_pos = curr_rec_pos;
            }
            prev_rec_pos = (curr_rec_pos + prev_point) * 0.5;
        } else {
            // Middle point
            prev_rec_pos = (curr_rec_pos + fine_tracking_path.at(i-1).first) * 0.5;
            next_rec_pos = (curr_rec_pos + fine_tracking_path.at(i+1).first) * 0.5;
        }
        
        dx[i] = (curr_rec_pos - prev_rec_pos).magnitude() + (curr_rec_pos - next_rec_pos).magnitude();

        // std::cout << i << " " << dx[i] << std::endl;
    }
    
    // Build response matrices using geometry information
    for (int i = 0; i < n_3D_pos; i++) {
        WireCell::Point curr_rec_pos = fine_tracking_path.at(i).first;
        auto segment = fine_tracking_path.at(i).second;
        auto cluster = segment->cluster();
        const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
        double cluster_t0 = cluster->get_cluster_t0();
        auto first_blob = cluster->children()[0];
        int cur_ntime_ticks = first_blob->slice_index_max() - first_blob->slice_index_min();


        int apa = paf.at(i).first;
        int face = paf.at(i).second;
                
        WirePlaneId wpid_key(kAllLayers, face, apa);
        
        // Get geometry parameters from wpid_offsets and wpid_slopes
        auto offset_it = wpid_offsets.find(wpid_key);
        auto slope_it = wpid_slopes.find(wpid_key);
        auto geom_it = wpid_geoms.find(wpid_key);
        
        if (offset_it == wpid_offsets.end() || slope_it == wpid_slopes.end() || geom_it == wpid_geoms.end()) continue;
        
        double offset_t = std::get<0>(offset_it->second);
        double offset_u = std::get<1>(offset_it->second);
        double offset_v = std::get<2>(offset_it->second);
        double offset_w = std::get<3>(offset_it->second);
        
        double slope_x = std::get<0>(slope_it->second);
        auto slope_yu = std::get<1>(slope_it->second).first;
        auto slope_zu = std::get<1>(slope_it->second).second;
        auto slope_yv = std::get<2>(slope_it->second).first;
        auto slope_zv = std::get<2>(slope_it->second).second;
        auto slope_yw = std::get<3>(slope_it->second).first;
        auto slope_zw = std::get<3>(slope_it->second).second;
        
        double time_tick_width = std::get<0>(geom_it->second);
        double pitch_u = std::get<1>(geom_it->second);
        double pitch_v = std::get<2>(geom_it->second);
        double pitch_w = std::get<3>(geom_it->second);

        
        // Calculate previous and next positions for Gaussian integration
        WireCell::Point prev_rec_pos, next_rec_pos;
        if (i == 0) {
            if (n_3D_pos > 1) {
                WireCell::Point next_point = fine_tracking_path.at(i+1).first;
                next_rec_pos = (curr_rec_pos + next_point) * 0.5;
                WireCell::Vector dir = next_point - curr_rec_pos;
                double length = dir.magnitude();
                if (length > 0) {
                    prev_rec_pos = curr_rec_pos - (dir / length) * dis_end_point_ext;
                } else {
                    prev_rec_pos = curr_rec_pos;
                }
            } else {
                prev_rec_pos = next_rec_pos = curr_rec_pos;
            }
        } else if (i == n_3D_pos - 1) {
            WireCell::Point prev_point = fine_tracking_path.at(i-1).first;
            prev_rec_pos = (curr_rec_pos + prev_point) * 0.5;
            WireCell::Vector dir = curr_rec_pos - prev_point;
            double length = dir.magnitude();
            if (length > 0) {
                next_rec_pos = curr_rec_pos + (dir / length) * dis_end_point_ext;
            } else {
                next_rec_pos = curr_rec_pos;
            }
        } else {
            prev_rec_pos = (curr_rec_pos + fine_tracking_path.at(i-1).first) * 0.5;
            next_rec_pos = (curr_rec_pos + fine_tracking_path.at(i+1).first) * 0.5;
        }
        
        // Create Gaussian integration points and weights
        std::vector<double> centers_U, centers_V, centers_W, centers_T;
        std::vector<double> sigmas_U, sigmas_V, sigmas_W, sigmas_T;
        std::vector<double> weights;
        
        // Sample 5 points along each half-segment
        for (int j = 0; j < 5; j++) {
            // First half (prev to curr)
            WireCell::Point reco_pos = prev_rec_pos + (curr_rec_pos - prev_rec_pos) * (j + 0.5) / 5.0;
            // find out the raw position ...
            auto reco_pos_raw = transform->backward(reco_pos, cluster_t0, apa, face);

            double central_T = offset_t + slope_x * reco_pos_raw.x();
            double central_U = offset_u + (slope_yu * reco_pos_raw.y() + slope_zu * reco_pos_raw.z());
            double central_V = offset_v + (slope_yv * reco_pos_raw.y() + slope_zv * reco_pos_raw.z());
            double central_W = offset_w + (slope_yw * reco_pos_raw.y() + slope_zw * reco_pos_raw.z());
            double weight = (curr_rec_pos - prev_rec_pos).magnitude();
            
            // Calculate drift time and diffusion
            double drift_time = std::max(m_params.min_drift_time, reco_pos.x() / time_tick_width * 0.5*units::us );
            double diff_sigma_L = sqrt(2 * DL * drift_time);
            double diff_sigma_T = sqrt(2 * DT * drift_time);
            
            double sigma_L = sqrt(pow(diff_sigma_L, 2) + pow(add_sigma_L, 2)) / time_tick_width;
            double sigma_T_u = sqrt(pow(diff_sigma_T, 2) + pow(ind_sigma_u_T, 2)) / pitch_u;
            double sigma_T_v = sqrt(pow(diff_sigma_T, 2) + pow(ind_sigma_v_T, 2)) / pitch_v;
            double sigma_T_w = sqrt(pow(diff_sigma_T, 2) + pow(col_sigma_w_T, 2)) / pitch_w;
            
            centers_U.push_back(central_U);
            centers_V.push_back(central_V);
            centers_W.push_back(central_W);
            centers_T.push_back(central_T);
            weights.push_back(weight);
            sigmas_U.push_back(sigma_T_u);
            sigmas_V.push_back(sigma_T_v);
            sigmas_W.push_back(sigma_T_w);
            sigmas_T.push_back(sigma_L);
            
            // Second half (curr to next)
            reco_pos = next_rec_pos + (curr_rec_pos - next_rec_pos) * (j + 0.5) / 5.0;
            reco_pos_raw = transform->backward(reco_pos, cluster_t0, apa, face);

            central_T = offset_t + slope_x * reco_pos_raw.x();
            central_U = offset_u + (slope_yu * reco_pos_raw.y() + slope_zu * reco_pos_raw.z());
            central_V = offset_v + (slope_yv * reco_pos_raw.y() + slope_zv * reco_pos_raw.z());
            central_W = offset_w + (slope_yw * reco_pos_raw.y() + slope_zw * reco_pos_raw.z());
            weight = (curr_rec_pos - next_rec_pos).magnitude();

            drift_time = std::max(m_params.min_drift_time, reco_pos.x() / time_tick_width * 0.5*units::us );
            diff_sigma_L = sqrt(2 * DL * drift_time);
            diff_sigma_T = sqrt(2 * DT * drift_time);

            // std::cout << drift_time << " " << DL << " " << DT << " " << diff_sigma_L << " " << diff_sigma_T << std::endl;


            sigma_L = sqrt(pow(diff_sigma_L, 2) + pow(add_sigma_L, 2)) / time_tick_width;
            sigma_T_u = sqrt(pow(diff_sigma_T, 2) + pow(ind_sigma_u_T, 2)) / pitch_u;
            sigma_T_v = sqrt(pow(diff_sigma_T, 2) + pow(ind_sigma_v_T, 2)) / pitch_v;
            sigma_T_w = sqrt(pow(diff_sigma_T, 2) + pow(col_sigma_w_T, 2)) / pitch_w;
            
            centers_U.push_back(central_U);
            centers_V.push_back(central_V);
            centers_W.push_back(central_W);
            centers_T.push_back(central_T);
            weights.push_back(weight);
            sigmas_U.push_back(sigma_T_u);
            sigmas_V.push_back(sigma_T_v);
            sigmas_W.push_back(sigma_T_w);
            sigmas_T.push_back(sigma_L);
        }

        // std::cout << i << " U ";
        // for (size_t idx = 0; idx < centers_U.size(); ++idx) {
        //     std::cout << centers_U[idx] << " ";
        // }
        // std::cout << std::endl;

        // std::cout << i << " V ";
        // for (size_t idx = 0; idx < centers_V.size(); ++idx) {
        //     std::cout << centers_V[idx] << " ";
        // }
        // std::cout << std::endl;

        // std::cout << i << " W ";
        // for (size_t idx = 0; idx < centers_W.size(); ++idx) {
        //     std::cout << centers_W[idx] << " ";
        // }
        // std::cout << std::endl;

        // std::cout << i << " T ";
        // for (size_t idx = 0; idx < centers_T.size(); ++idx) {
        //     std::cout << centers_T[idx] << " ";
        // }
        // std::cout << std::endl;

        // std::cout << i << " Weights ";
        // for (size_t idx = 0; idx < weights.size(); ++idx) {
        //     std::cout << weights[idx] << " ";
        // }
        // std::cout << std::endl;

        // std::cout <<i << " SU ";
        // for (size_t idx = 0; idx < sigmas_U.size(); ++idx) {
        //     std::cout << sigmas_U[idx] << " ";
        // }
        // std::cout << std::endl;

        // std::cout << i << " SV ";
        // for (size_t idx = 0; idx < sigmas_V.size(); ++idx) {
        //     std::cout << sigmas_V[idx] << " ";
        // }
        // std::cout << std::endl;

        // std::cout << i << " SW ";
        // for (size_t idx = 0; idx < sigmas_W.size(); ++idx) {
        //     std::cout << sigmas_W[idx] << " ";
        // }
        // std::cout << std::endl;

        // std::cout << i << " ST ";
        // for (size_t idx = 0; idx < sigmas_T.size(); ++idx) {
        //     std::cout << sigmas_T[idx] << " ";
        // }
        // std::cout << std::endl;

        // Fill response matrices using Gaussian integration
        int n_u = 0;
        std::set<std::pair<int, int>> set_UT;
        for (const auto& [coord_key, result] : map_U_charge_2D) {
            const auto& measurement = result.first;
            const auto& Coord2D_set = result.second;    

            for (const auto& coord2d : Coord2D_set) {
                // coord2d: TrackFitting::Coord2D
                // Only process if plane matches U
                if (coord2d.plane != kUlayer || coord2d.apa != apa || coord2d.face != face) continue;
                int wire = coord2d.wire;
                int time = coord2d.time;

                set_UT.insert(std::make_pair(wire, time));

                // if (wire !=938 || time != 7176) continue;

                if (abs(wire - centers_U.front()) <= m_params.search_range && abs(time - centers_T.front()) <= m_params.search_range * cur_ntime_ticks) {
                    double value = cal_gaus_integral_seg(time, wire, centers_T, sigmas_T, centers_U, sigmas_U, weights, 0, 4, cur_ntime_ticks);


                    if (measurement.flag == 0 && value > 0) reg_flag_u[i] = 1; // Dead channel
                    
                    if (value > 0 && measurement.charge > 0 && measurement.flag != 0) {
                        double charge = measurement.charge;
                        double charge_err = measurement.charge_err;
                        double total_err = sqrt(pow(charge_err, 2) + pow(charge * rel_uncer_ind, 2) + pow(add_uncer_ind, 2));
                        RU.insert(n_u, i) = value / total_err;

                        // std::cout << time << " " << wire << " " << i << " " << value << " " << total_err << std::endl;
                    }
                }
            }
            n_u++;
        }
        int n_v = 0;
        std::set<std::pair<int, int>> set_VT;
        for (const auto& [coord_key, result] : map_V_charge_2D) {
            const auto& measurement = result.first;
            const auto& Coord2D_set = result.second;    

            for (const auto& coord2d : Coord2D_set) {
                // coord2d: TrackFitting::Coord2D
                // Only process if plane matches V
                if (coord2d.plane != kVlayer || coord2d.apa != apa || coord2d.face != face) continue;
                int wire = coord2d.wire;
                int time = coord2d.time;
                set_VT.insert(std::make_pair(wire, time));

                if (abs(wire - centers_V.front()) <= m_params.search_range && abs(time - centers_T.front()) <= m_params.search_range * cur_ntime_ticks) {
                    double value = cal_gaus_integral_seg(time, wire, centers_T, sigmas_T, centers_V, sigmas_V, weights, 0, 4, cur_ntime_ticks);

                    if (measurement.flag == 0 && value > 0) reg_flag_v[i] = 1; // Dead channel

                    if (value > 0 && measurement.charge > 0 && measurement.flag != 0) {
                        double charge = measurement.charge;
                        double charge_err = measurement.charge_err;
                        double total_err = sqrt(pow(charge_err, 2) + pow(charge * rel_uncer_ind, 2) + pow(add_uncer_ind, 2));
                        RV.insert(n_v, i) = value / total_err;
                    }

                }
            }
            n_v++;
        }
        int n_w = 0;
        std::set<std::pair<int, int>> set_WT;
        for (const auto& [coord_key, result] : map_W_charge_2D) {
            const auto& measurement = result.first;
            const auto& Coord2D_set = result.second;    

            for (const auto& coord2d : Coord2D_set) {
                // coord2d: TrackFitting::Coord2D
                // Only process if plane matches W
                if (coord2d.plane != kWlayer || coord2d.apa != apa || coord2d.face != face) continue;
                int wire = coord2d.wire;
                int time = coord2d.time;
                set_WT.insert(std::make_pair(wire, time));
                if (abs(wire - centers_W.front()) <= m_params.search_range && abs(time - centers_T.front()) <= m_params.search_range * cur_ntime_ticks) {
                    double value = cal_gaus_integral_seg(time, wire, centers_T, sigmas_T, centers_W, sigmas_W, weights, 0, 4, cur_ntime_ticks);

                    if (measurement.flag == 0 && value > 0) reg_flag_w[i] = 1; // Dead channel

                    if (value > 0 && measurement.charge > 0 && measurement.flag != 0) {
                        double charge = measurement.charge;
                        double charge_err = measurement.charge_err;
                        double total_err = sqrt(pow(charge_err, 2) + pow(charge * rel_uncer_col, 2) + pow(add_uncer_col, 2));
                        RW.insert(n_w, i) = value / total_err;
                    }
                }
            }
            n_u++;
        }


        // Additional dead channel checks
        if (reg_flag_u[i] == 0) { // apa, face
            for (size_t kk = 0; kk < centers_U.size(); kk++) {
                if (set_UT.find(std::make_pair(std::round(centers_U[kk]), std::round(centers_T[kk]/cur_ntime_ticks)*cur_ntime_ticks)) == set_UT.end()) {
                    reg_flag_u[i] = 1;
                    break;
                }
            }
        }
        if (reg_flag_v[i] == 0) { // apa, face
            for (size_t kk = 0; kk < centers_V.size(); kk++) {
                if (set_VT.find(std::make_pair(std::round(centers_V[kk]), std::round(centers_T[kk]/cur_ntime_ticks)*cur_ntime_ticks)) == set_VT.end()) {
                    reg_flag_v[i] = 1;
                    break;
                }
            }
        }
        if (reg_flag_w[i] == 0) { // apa, face
            for (size_t kk = 0; kk < centers_W.size(); kk++) {
                if (set_WT.find(std::make_pair(std::round(centers_W[kk]), std::round(centers_T[kk]/cur_ntime_ticks)*cur_ntime_ticks)) == set_WT.end()) {
                    reg_flag_w[i] = 1;
                    break;
                }
            }
        }
        // std::cout << i << " " << reg_flag_u[i] << " " << reg_flag_v[i] << " " << reg_flag_w[i] << std::endl;

    }
    
    // Calculate compact matrices for overlap analysis
    Eigen::SparseMatrix<double> RUT = RU.transpose();
    Eigen::SparseMatrix<double> RVT = RV.transpose();
    Eigen::SparseMatrix<double> RWT = RW.transpose();
    
    Eigen::SparseMatrix<double> MU(n_2D_u, n_2D_u), MV(n_2D_v, n_2D_v), MW(n_2D_w, n_2D_w);
    for (int k = 0; k < n_2D_u; k++) MU.insert(k, k) = 1;
    for (int k = 0; k < n_2D_v; k++) MV.insert(k, k) = 1;
    for (int k = 0; k < n_2D_w; k++) MW.insert(k, k) = 1;
    
    std::vector<std::pair<double, double>> overlap_u = calculate_compact_matrix(MU, RUT, n_2D_u, n_3D_pos, 3);
    std::vector<std::pair<double, double>> overlap_v = calculate_compact_matrix(MV, RVT, n_2D_v, n_3D_pos, 3);
    std::vector<std::pair<double, double>> overlap_w = calculate_compact_matrix(MW, RWT, n_2D_w, n_3D_pos, 2);
    
    // Add regularization based on dead channels and overlaps
    Eigen::SparseMatrix<double> FMatrix(n_3D_pos, n_3D_pos);

    const double dead_ind_weight = m_params.dead_ind_weight;
    const double dead_col_weight = m_params.dead_col_weight;
    const double close_ind_weight = m_params.close_ind_weight;
    const double close_col_weight = m_params.close_col_weight;

    for (int i = 0; i < n_3D_pos; i++) {
        bool flag_u = reg_flag_u[i];
        bool flag_v = reg_flag_v[i];
        bool flag_w = reg_flag_w[i];
        
        if (n_3D_pos != 1) {
            double weight = 0;
            if (flag_u) weight += dead_ind_weight;
            if (flag_v) weight += dead_ind_weight;
            if (flag_w) weight += dead_col_weight;
            
            if (i==0){
                if (overlap_u[i].second > m_params.overlap_th) weight += close_ind_weight * pow(2 * overlap_u[i].second - 1, 2);
                if (overlap_v[i].second > m_params.overlap_th) weight += close_ind_weight * pow(2 * overlap_v[i].second - 1, 2);
                if (overlap_w[i].second > m_params.overlap_th) weight += close_col_weight * pow(2 * overlap_w[i].second - 1, 2);
            }else if (i==n_3D_pos-1){
                if (overlap_u[i].second > m_params.overlap_th) weight += close_ind_weight * pow(2 * overlap_u[i].first - 1, 2);
                if (overlap_v[i].second > m_params.overlap_th) weight += close_ind_weight * pow(2 * overlap_v[i].first - 1, 2);
                if (overlap_w[i].second > m_params.overlap_th) weight += close_col_weight * pow(2 * overlap_w[i].first - 1, 2);
            }else{
                if (overlap_u.at(i).first + overlap_u.at(i).second > 2*m_params.overlap_th) weight += close_ind_weight * pow(overlap_u.at(i).first + overlap_u.at(i).second - 1,2);
                if (overlap_v.at(i).first + overlap_v.at(i).second > 2*m_params.overlap_th) weight += close_ind_weight * pow(overlap_v.at(i).first + overlap_v.at(i).second - 1,2);
                if (overlap_w.at(i).first + overlap_w.at(i).second > 2*m_params.overlap_th) weight += close_col_weight * pow(overlap_w.at(i).first + overlap_w.at(i).second - 1,2);

            }
            
            double dx_norm = (dx[i] + 0.001*units::cm) / m_params.dx_norm_length; // Normalize by 0.6 mm
            
            if (i == 0) {
                FMatrix.insert(0, 0) = -weight / dx_norm;
                if (n_3D_pos > 1) FMatrix.insert(0, 1) = weight / dx_norm;
            } else if (i == n_3D_pos - 1) {
                FMatrix.insert(i, i) = -weight / dx_norm; 
                FMatrix.insert(i, i-1) = weight / dx_norm;
            } else {
                FMatrix.insert(i, i) = -2.0 * weight / dx_norm;
                FMatrix.insert(i, i+1) = weight / dx_norm;
                FMatrix.insert(i, i-1) = weight / dx_norm;
            }
        }
    }
    
    // Apply regularization scaling
    double lambda = m_params.lambda;
    FMatrix *= lambda;
    if (!flag_dQ_dx_fit_reg) FMatrix *= 0.01;
    
    // Solve the linear system
    Eigen::SparseMatrix<double> FMatrixT = FMatrix.transpose();
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    
    Eigen::VectorXd b = RUT * MU * data_u_2D + RVT * MV * data_v_2D + RWT * MW * data_w_2D;
    Eigen::SparseMatrix<double> A = RUT * MU * RU + RVT * MV * RV + RWT * MW * RW + FMatrixT * FMatrix;
    
    solver.compute(A);
    pos_3D = solver.solveWithGuess(b, pos_3D_init);
    
    if (std::isnan(solver.error())) {
        pos_3D = solver.solve(b);
    }
    
    // Extract dQ values and apply corrections
    dQ.resize(n_3D_pos);
    for (int i=0;i!=n_3D_pos;i++){
        dQ[i] = pos_3D(i);
    }
 
    // Calculate predictions and reduced chi-squared
    pred_data_u_2D = RU * pos_3D;
    pred_data_v_2D = RV * pos_3D;
    pred_data_w_2D = RW * pos_3D;
    
    // Calculate reduced chi-squared for each 3D point
    reduced_chi2.resize(n_3D_pos);
    for (int k = 0; k < n_3D_pos; k++) {
        double sum[3] = {0, 0, 0};
        double sum1[3] = {0, 0, 0};
        
        for (Eigen::SparseMatrix<double>::InnerIterator it(RU, k); it; ++it) {
            if (pred_data_u_2D(it.row()) > 0) {
                sum[0] += pow(data_u_2D(it.row()) - pred_data_u_2D(it.row()), 2) * 
                          (it.value() * pos_3D(k)) / pred_data_u_2D(it.row());
                sum1[0] += (it.value() * pos_3D(k)) / pred_data_u_2D(it.row());
            }
        }
        
        for (Eigen::SparseMatrix<double>::InnerIterator it(RV, k); it; ++it) {
            if (pred_data_v_2D(it.row()) > 0) {
                sum[1] += pow(data_v_2D(it.row()) - pred_data_v_2D(it.row()), 2) * 
                          (it.value() * pos_3D(k)) / pred_data_v_2D(it.row());
                sum1[1] += (it.value() * pos_3D(k)) / pred_data_v_2D(it.row());
            }
        }
        
        for (Eigen::SparseMatrix<double>::InnerIterator it(RW, k); it; ++it) {
            if (pred_data_w_2D(it.row()) > 0) {
                sum[2] += pow(data_w_2D(it.row()) - pred_data_w_2D(it.row()), 2) * 
                          (it.value() * pos_3D(k)) / pred_data_w_2D(it.row());
                sum1[2] += (it.value() * pos_3D(k)) / pred_data_w_2D(it.row());
            }
        }
        
        double total_chi2 = sum[0] + sum[1] + sum[2] / 4.0; // Weight collection plane differently
        double total_weight = sum1[0] + sum1[1] + sum1[2];
        
        reduced_chi2[k] = (total_weight > 0) ? sqrt(total_chi2 / total_weight) : 0;
    }
    
    // Restore original charge data
    recover_original_charge_data();
}

void TrackFitting::do_single_tracking(std::shared_ptr<PR::Segment> segment, bool flag_dQ_dx_fit_reg, bool flag_dQ_dx_fit, bool flag_force_load_data) {
      // Clear all internal tracking vectors
    fine_tracking_path.clear();
    dQ.clear();
    dx.clear();
    pu.clear();
    pv.clear();
    pw.clear();
    pt.clear();
    paf.clear();
    reduced_chi2.clear();
    
    bool flag_1st_tracking = true;
    bool flag_2nd_tracking = true;
    bool flag_dQ_dx = flag_dQ_dx_fit;
    
    // Prepare the data for the fit - collect charge information from 2D projections
    if (flag_force_load_data || global_rb_map.size() == 0){
        prepare_data();
        fill_global_rb_map();
    }

    // std::cout << "Global Blob Map: " << global_rb_map.size() << std::endl;

    // First round of organizing the path from the path_wcps (shortest path)
    double low_dis_limit = m_params.low_dis_limit;
    double end_point_limit = m_params.end_point_limit;

    if (m_segments.find(segment) == m_segments.end()) {
        // Handle empty segments case - could log warning or return
        return;
    }
    // auto segment = *m_segments.begin();
    
    auto pts = organize_orig_path(segment, low_dis_limit, end_point_limit); 
    if (pts.size() == 0) return;
    else if (pts.size() == 1) {
        const auto& segment_wcpts = segment->wcpts();
        if (!segment_wcpts.empty()) {
            const auto& last_segment_point = segment_wcpts.back().point;
            if (sqrt(pow(last_segment_point.x() - pts.back().x(), 2) + 
                     pow(last_segment_point.y() - pts.back().y(), 2) + 
                     pow(last_segment_point.z() - pts.back().z(), 2)) < 0.01*units::cm) {
                return;
            } else {
                WireCell::Point p2(last_segment_point.x(), last_segment_point.y(), last_segment_point.z());
                pts.push_back(p2);
            }
        }
    }

    std::cout << "After organization " << pts.size() << std::endl;

    std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> ptss;
    for (const auto& pt : pts) {
        ptss.emplace_back(pt, segment);
    }

    if (flag_1st_tracking) {
        form_map(ptss, m_params.end_point_factor, m_params.mid_point_factor, m_params.nlevel, m_params.time_tick_cut, m_params.charge_cut);
        trajectory_fit(ptss, 1, m_params.div_sigma);
    }
    // Check for very close start/end points and reset if needed
    if (ptss.size() == 2) {
        if (sqrt(pow(ptss.front().first.x() - ptss.back().first.x(), 2) + 
                 pow(ptss.front().first.y() - ptss.back().first.y(), 2) + 
                 pow(ptss.front().first.z() - ptss.back().first.z(), 2)) < 0.1*units::cm) {
            ptss.clear();
            const auto& segment_wcpts = segment->wcpts();
            WireCell::Point p1(segment_wcpts.front().point.x(), segment_wcpts.front().point.y(), segment_wcpts.front().point.z());
            ptss.push_back(std::make_pair(p1,segment));
            WireCell::Point p2(segment_wcpts.back().point.x(), segment_wcpts.back().point.y(), segment_wcpts.back().point.z());
            ptss.push_back(std::make_pair(p2,segment));
        }
    } 
    
    if (ptss.size() <= 1) return;
    
    if (flag_2nd_tracking) {
        // Second round trajectory fit with tighter parameters
        low_dis_limit = m_params.low_dis_limit/2.;
        end_point_limit = m_params.end_point_limit/2.;

        pts.clear();
        for (const auto& pt_pair : ptss) {
            pts.push_back(pt_pair.first);
        }
        

        // // hack pts 
        // pts.resize(18);
        // pts[0] = WireCell::Point(2192.13, -873.682, 2094.73);
        // pts[1] = WireCell::Point(2190.05, -877.433, 2095.44);
        // pts[2] = WireCell::Point(2187.79, -882.693, 2096.37);
        // pts[3] = WireCell::Point(2181.57, -896.034, 2099.36);
        // pts[4] = WireCell::Point(2178.83, -904.709, 2100.54);
        // pts[5] = WireCell::Point(2173.46, -917.087, 2103.09);
        // pts[6] = WireCell::Point(2168.32, -927.231, 2105.47);
        // pts[7] = WireCell::Point(2162.27, -938.16, 2108.49);
        // pts[8] = WireCell::Point(2158.9, -946.832, 2110.14);
        // pts[9] = WireCell::Point(2153.28, -958.607, 2113.05);
        // pts[10] = WireCell::Point(2147.69, -966.479, 2115.99);
        // pts[11] = WireCell::Point(2143.06, -977.06, 2118.02);
        // pts[12] = WireCell::Point(2139.68, -984.335, 2119.6);
        // pts[13] = WireCell::Point(2134.63, -997.786, 2121.84);
        // pts[14] = WireCell::Point(2127.28, -1006.52, 2125.44);
        // pts[15] = WireCell::Point(2123.03, -1016.8, 2127.48);
        // pts[16] = WireCell::Point(2121.65, -1024.26, 2128.3);
        // pts[17] = WireCell::Point(2122.02, -1026.47, 2128.18);
        // //

        organize_ps_path(segment, pts, low_dis_limit, end_point_limit);
        
        // std::cout << pts.size() << std::endl;
        // for (size_t i = 0; i < pts.size(); ++i) {
        //     std::cout << "pts[" << i << "] = (" << pts[i].x() << ", " << pts[i].y() << ", " << pts[i].z() << ")" << std::endl;
        // }

        ptss.clear();
        for (const auto& pt : pts) {
            ptss.emplace_back(pt, segment);
        }
        form_map(ptss, m_params.end_point_factor, m_params.mid_point_factor, m_params.nlevel, m_params.time_tick_cut, m_params.charge_cut);
        trajectory_fit(ptss, 2, m_params.div_sigma);

        pts.clear();
        for (const auto& pt_pair : ptss) {
            pts.push_back(pt_pair.first);
        }
        
        // Final path organization
        organize_ps_path(segment, pts, low_dis_limit, 0);

        // Check for very close start/end points and reset if needed
        if (pts.size() == 2) {
            if (sqrt(pow(pts.front().x() - pts.back().x(), 2) + 
                 pow(pts.front().y() - pts.back().y(), 2) + 
                 pow(pts.front().z() - pts.back().z(), 2)) < 0.1*units::cm) {
            pts.clear();
            const auto& segment_wcpts = segment->wcpts();
            WireCell::Point p1(segment_wcpts.front().point.x(), segment_wcpts.front().point.y(), segment_wcpts.front().point.z());
            pts.push_back(p1);
            WireCell::Point p2(segment_wcpts.back().point.x(), segment_wcpts.back().point.y(), segment_wcpts.back().point.z());
            pts.push_back(p2);
            }
        }

        // std::cout << pts.size() << std::endl;
        // for (size_t i = 0; i < pts.size(); ++i) {
        //     std::cout << "pts[" << i << "] = (" << pts[i].x() << ", " << pts[i].y() << ", " << pts[i].z() << ")" << std::endl;
        // }

        // hack pts ...
        pts.clear();
        pts.push_back(WireCell::Point(2192.09, -872.848, 2094.53));
        pts.push_back(WireCell::Point(2189.97, -878.182, 2095.57));
        pts.push_back(WireCell::Point(2188.06, -882.53, 2096.41));
        pts.push_back(WireCell::Point(2186.15, -886.879, 2097.25));
        pts.push_back(WireCell::Point(2183.37, -891.079, 2098.36));
        pts.push_back(WireCell::Point(2180.51, -898.565, 2099.78));
        pts.push_back(WireCell::Point(2177.7, -905.317, 2101.2));
        pts.push_back(WireCell::Point(2174.72, -912.331, 2102.5));
        pts.push_back(WireCell::Point(2170.58, -918.449, 2104.41));
        pts.push_back(WireCell::Point(2168.22, -925.268, 2105.48));
        pts.push_back(WireCell::Point(2166.1, -929.509, 2106.5));
        pts.push_back(WireCell::Point(2163.08, -936.788, 2108.1));
        pts.push_back(WireCell::Point(2159.92, -942.301, 2109.39));
        pts.push_back(WireCell::Point(2157.91, -947.826, 2110.62));
        pts.push_back(WireCell::Point(2155.32, -951.674, 2111.93));
        pts.push_back(WireCell::Point(2152.74, -958.492, 2113.46));
        pts.push_back(WireCell::Point(2150.38, -963.231, 2114.67));
        pts.push_back(WireCell::Point(2147.41, -967.608, 2116.04));
        pts.push_back(WireCell::Point(2144.44, -971.985, 2117.41));
        pts.push_back(WireCell::Point(2142.82, -978.976, 2118.24));
        pts.push_back(WireCell::Point(2140.26, -983.316, 2119.43));
        pts.push_back(WireCell::Point(2136.82, -990.357, 2120.85));
        pts.push_back(WireCell::Point(2133.68, -997.199, 2122.16));
        pts.push_back(WireCell::Point(2131.01, -1001.64, 2123.41));
        pts.push_back(WireCell::Point(2127.42, -1006.39, 2125.27));
        pts.push_back(WireCell::Point(2125.14, -1011.79, 2126.47));
        pts.push_back(WireCell::Point(2121.17, -1018.17, 2128.31));
        pts.push_back(WireCell::Point(2119.77, -1024.55, 2128.91));
        //


        // Generate 2D projections
        pu.clear();
        pv.clear();
        pw.clear();
        pt.clear();
        ptss.clear();
        paf.clear();
        for (const auto& p : pts) {
            auto cluster = segment->cluster();
            const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
            double cluster_t0 = cluster->get_cluster_t0();

            auto test_wpid = m_dv->contained_by(p);

            if (test_wpid.apa()==-1) continue;
            int apa = test_wpid.apa();
            int face = test_wpid.face();

            auto p_raw = transform->backward(p, cluster_t0, apa, face);
            WirePlaneId wpid(kAllLayers, face, apa);
            auto offset_it = wpid_offsets.find(wpid);
            auto slope_it = wpid_slopes.find(wpid);

            auto offset_t = std::get<0>(offset_it->second);
            auto offset_u = std::get<1>(offset_it->second);
            auto offset_v = std::get<2>(offset_it->second);
            auto offset_w = std::get<3>(offset_it->second);
            auto slope_x = std::get<0>(slope_it->second);
            auto slope_yu = std::get<1>(slope_it->second).first;
            auto slope_zu = std::get<1>(slope_it->second).second;
            auto slope_yv = std::get<2>(slope_it->second).first;
            auto slope_zv = std::get<2>(slope_it->second).second;
            auto slope_yw = std::get<3>(slope_it->second).first;
            auto slope_zw = std::get<3>(slope_it->second).second;

            ptss.emplace_back(p, segment);
            pu.push_back(offset_u + (slope_yu * p_raw.y() + slope_zu * p_raw.z()));
            pv.push_back(offset_v + (slope_yv * p_raw.y() + slope_zv * p_raw.z()));
            pw.push_back(offset_w + (slope_yw * p_raw.y() + slope_zw * p_raw.z()));
            pt.push_back(offset_t + slope_x * p_raw.x());
            paf.push_back(std::make_pair(apa, face));

        }
    }
    
    
    fine_tracking_path = ptss;

    if (flag_dQ_dx) {
        // Store the fine tracking path as pairs of (Point, Segment)
        // Perform dQ/dx fit using the prepared charge data
        dQ_dx_fit(end_point_limit, flag_dQ_dx_fit_reg);
    } else {
        dQ_dx_fill(end_point_limit);
    }

    // Now put the results back into the
    // Create vector of Fit objects from the internal tracking results
    std::vector<PR::Fit> segment_fits;
    segment_fits.reserve(fine_tracking_path.size());
    
    // Check that all vectors have consistent sizes
    size_t npoints = fine_tracking_path.size();
    if (dQ.size() != npoints || dx.size() != npoints || 
        pu.size() != npoints || pv.size() != npoints || 
        pw.size() != npoints || pt.size() != npoints || 
        reduced_chi2.size() != npoints) {
        throw std::runtime_error("TrackFitting::do_single_tracking: inconsistent vector sizes for fit output!");
    }
    
    // Calculate cumulative range (distance along track)
    std::vector<double> cumulative_range(npoints, 0.0);
    if (npoints > 1) {
        for (size_t i = 1; i < npoints; ++i) {
            const auto& p1 = fine_tracking_path[i-1].first;
            const auto& p2 = fine_tracking_path[i].first;
            double step_distance = sqrt(pow(p2.x() - p1.x(), 2) + 
                                      pow(p2.y() - p1.y(), 2) + 
                                      pow(p2.z() - p1.z(), 2));
            cumulative_range[i] = cumulative_range[i-1] + step_distance;
        }
    }
    
    // Convert internal results to PR::Fit objects
    for (size_t i = 0; i < npoints; ++i) {
        PR::Fit fit;
        
        // Set the fitted 3D point
        fit.point = fine_tracking_path[i].first;
        
        // Set physics quantities
        fit.dQ =  dQ[i];
        fit.dx = dx[i];
        fit.pu = pu[i];
        fit.pv = pv[i];
        fit.pw = pw[i];
        fit.pt = pt[i];
        fit.paf = paf[i];
        fit.reduced_chi2 = reduced_chi2[i];

        // Set trajectory information
        fit.index = static_cast<int>(i);
        fit.range = cumulative_range[i];
        
        // Set fix flags (typically fix endpoints for track fitting)
        fit.flag_fix = false;        
        segment_fits.push_back(fit);
    }
    
    // Assign the fits to the segment
    segment->fits(segment_fits);
    
}
