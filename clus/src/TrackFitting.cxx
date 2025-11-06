#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/TrackFitting_Util.h"
#include "WireCellClus/PRSegmentFunctions.h"

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



void TrackFitting::clear_graph(){
    m_graph = nullptr;
    m_clusters.clear();
    m_blobs.clear(); 
}


void TrackFitting::clear_segments(){
    m_segments.clear();
    m_clusters.clear();
    m_blobs.clear(); 
}

void TrackFitting::add_graph(std::shared_ptr<PR::Graph> graph){
    m_graph = graph;

    if (!m_graph){
        return;
    }

    // Get edges from the graph and extract segments from them
    auto edge_range = boost::edges(*m_graph);
    std::set<std::shared_ptr<PR::Segment>> segments_set;
    
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        auto& edge_bundle = (*m_graph)[*e_it];
        if (edge_bundle.segment) {
            segments_set.insert(edge_bundle.segment);
            m_clusters.insert(edge_bundle.segment->cluster());
        }
    }

    if (m_grouping == nullptr && !segments_set.empty()){
        auto first_segment = *segments_set.begin();
        m_grouping = first_segment->cluster()->grouping();
        BuildGeometry();
    }

    for (auto& cluster: m_clusters){
        for (auto& blob: cluster->children()){
            m_blobs.insert(blob);
        }
    }

    std::cout << "TrackFitting: Added graph with " << segments_set.size() << " segments." << " " << m_clusters.size() << " " << m_blobs.size() << std::endl;
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

// ============================================================================
// Helper functions for organize_segments_path methods
// ============================================================================

void TrackFitting::check_and_reset_close_vertices() {
    if (!m_graph) return;
    
    auto edge_range = boost::edges(*m_graph);
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        auto& edge_bundle = (*m_graph)[*e_it];
        auto segment = edge_bundle.segment;
        if (!segment) continue;
        
        // Get vertices connected to this segment
        auto vd1 = boost::source(*e_it, *m_graph);
        auto vd2 = boost::target(*e_it, *m_graph);
        auto& v1_bundle = (*m_graph)[vd1];
        auto& v2_bundle = (*m_graph)[vd2];
        auto start_v = v1_bundle.vertex;
        auto end_v = v2_bundle.vertex;
        
        if (!start_v || !end_v) continue;
        
        // Determine which vertex corresponds to which end of the segment
        const auto& segment_wcpts = segment->wcpts();
        if (segment_wcpts.empty()) continue;
        
        // Check vertex ordering by comparing with segment endpoints
        // Use wcpt().point to determine which vertex matches the segment's front/back
        double dist_start_front = (start_v->wcpt().point - segment_wcpts.front().point).magnitude();
        double dist_start_back  = (start_v->wcpt().point - segment_wcpts.back().point).magnitude();

        // std::cout << "Vertex Distance Check: " << start_v->wcpt().point << " to front " << segment_wcpts.front().point << " = " << dist_start_front/units::cm << " cm, to back " << segment_wcpts.back().point << " = " << dist_start_back/units::cm << " cm" << std::endl;

        // If start_v is closer to the back, swap
        if (dist_start_back < dist_start_front) {
            std::swap(start_v, end_v);
            std::swap(vd1, vd2);
        }
        
        // Check if vertices are too close together
        double vertex_distance = sqrt(
            pow(start_v->fit().point.x() - end_v->fit().point.x(), 2) +
            pow(start_v->fit().point.y() - end_v->fit().point.y(), 2) +
            pow(start_v->fit().point.z() - end_v->fit().point.z(), 2)
        );
        
        if (vertex_distance < 0.01 * units::cm) {
            // Reset vertices to original points if they are endpoints (degree 1)
            if (boost::degree(vd1, *m_graph) == 1) {
                PR::Fit start_fit = start_v->fit();
                start_fit.point = start_v->wcpt().point;
                start_v->fit(start_fit);
            }
            if (boost::degree(vd2, *m_graph) == 1) {
                PR::Fit end_fit = end_v->fit();
                end_fit.point = end_v->wcpt().point;
                end_v->fit(end_fit);
            }
        }
    }
}

bool TrackFitting::get_ordered_segment_vertices(
    std::shared_ptr<PR::Segment> segment,
    const PR::edge_descriptor& ed,
    std::shared_ptr<PR::Vertex>& start_v,
    std::shared_ptr<PR::Vertex>& end_v,
    PR::node_descriptor& vd1,
    PR::node_descriptor& vd2
) {
    if (!m_graph || !segment) return false;
    
    // Get vertices connected to this segment
    vd1 = boost::source(ed, *m_graph);
    vd2 = boost::target(ed, *m_graph);
    auto& v1_bundle = (*m_graph)[vd1];
    auto& v2_bundle = (*m_graph)[vd2];
    start_v = v1_bundle.vertex;
    end_v = v2_bundle.vertex;
    
    if (!start_v || !end_v) return false;
    
    // Determine which vertex corresponds to which end of the segment
    const auto& segment_wcpts = segment->wcpts();
    if (segment_wcpts.empty()) return false;
    
    // Use wcpt().point to determine which vertex matches the segment's front/back
    double dist_start_front = (start_v->wcpt().point - segment_wcpts.front().point).magnitude();
    double dist_start_back  = (start_v->wcpt().point - segment_wcpts.back().point).magnitude();

    // std::cout << "Vertex Distance Check 1: " << start_v->wcpt().point << " to front " << segment_wcpts.front().point << " = " << dist_start_front/units::cm << " cm, to back " << segment_wcpts.back().point << " = " << dist_start_back/units::cm << " cm" << std::endl;

    // If start_v is closer to the back, swap
    if (dist_start_back < dist_start_front) {
        std::swap(start_v, end_v);
        std::swap(vd1, vd2);
    }
    
    return true;
}

std::vector<PR::Fit> TrackFitting::generate_fits_with_projections(
    std::shared_ptr<PR::Segment> segment,
    const std::vector<WireCell::Point>& pts
) {
    std::vector<PR::Fit> fits;
    if (!segment || pts.empty()) return fits;
    
    auto cluster = segment->cluster();
    const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
    double cluster_t0 = cluster->get_cluster_t0();
    
    for (size_t i = 0; i != pts.size(); i++) {
        PR::Fit fit;
        fit.point = pts.at(i);
        fit.dQ = 0;
        fit.dx = -1;
        fit.reduced_chi2 = 0;
        
        // Generate 2D projections
        auto test_wpid = m_dv->contained_by(pts.at(i));
        if (test_wpid.apa() != -1 && test_wpid.face() != -1) {
            int apa = test_wpid.apa();
            int face = test_wpid.face();
            
            auto p_raw = transform->backward(pts.at(i), cluster_t0, apa, face);
            WirePlaneId wpid(kAllLayers, face, apa);
            auto offset_it = wpid_offsets.find(wpid);
            auto slope_it = wpid_slopes.find(wpid);
            
            if (offset_it != wpid_offsets.end() && slope_it != wpid_slopes.end()) {
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
                
                fit.pu = offset_u + (slope_yu * p_raw.y() + slope_zu * p_raw.z());
                fit.pv = offset_v + (slope_yv * p_raw.y() + slope_zv * p_raw.z());
                fit.pw = offset_w + (slope_yw * p_raw.y() + slope_zw * p_raw.z());
                fit.pt = offset_t + slope_x * p_raw.x();
                fit.paf = std::make_pair(apa, face);
            }
        }
        
        fits.push_back(fit);
    }
    
    return fits;
}

void TrackFitting::organize_segments_path_3rd(double step_size){
    if (!m_graph) return;
    
    // First pass: check for vertices that are too close together
    check_and_reset_close_vertices();
    
    // Second pass: organize segments path with uniform step size
    auto edge_range = boost::edges(*m_graph);
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        auto& edge_bundle = (*m_graph)[*e_it];
        auto segment = edge_bundle.segment;
        if (!segment) continue;
        
        // Get ordered vertices
        std::shared_ptr<PR::Vertex> start_v, end_v;
        PR::node_descriptor vd1, vd2;
        if (!get_ordered_segment_vertices(segment, *e_it, start_v, end_v, vd1, vd2)) continue;
        
        // Check if vertices are endpoints (degree == 1)
        bool flag_startv_end = (boost::degree(vd1, *m_graph) == 1);
        bool flag_endv_end = (boost::degree(vd2, *m_graph) == 1);
        
        std::vector<WireCell::Point> pts, curr_pts;
        
        // Get current fitted path from the segment
        if (!segment->fits().empty()) {
            for (const auto& fit : segment->fits()) {
                curr_pts.push_back(fit.point);
            }
        } else {
            // If no fits, use original wcpts
            for (const auto& wcpt : segment->wcpts()) {
                curr_pts.push_back(wcpt.point);
            }
        }
        
        // Examine end points
        curr_pts = examine_end_ps_vec(segment, curr_pts, flag_startv_end, flag_endv_end);
        
        WireCell::Point start_p = curr_pts.front();
        WireCell::Point end_p = curr_pts.back();
        
        // Build points with uniform step size
        pts.push_back(start_p);
        double extra_dis = 0;
        
        for (size_t i = 0; i != curr_pts.size(); i++) {
            WireCell::Point p1 = curr_pts.at(i);
            
            double dis_end = sqrt(pow(p1.x() - end_p.x(), 2) + pow(p1.y() - end_p.y(), 2) + pow(p1.z() - end_p.z(), 2));
            if (dis_end < step_size) continue;
            
            double dis_prev = sqrt(pow(p1.x() - pts.back().x(), 2) + pow(p1.y() - pts.back().y(), 2) + pow(p1.z() - pts.back().z(), 2));
            
            if (dis_prev + extra_dis > step_size) {
                extra_dis += dis_prev;
                while (extra_dis > step_size) {
                    WireCell::Point tmp_p(
                        pts.back().x() + (p1.x() - pts.back().x()) / dis_prev * step_size,
                        pts.back().y() + (p1.y() - pts.back().y()) / dis_prev * step_size,
                        pts.back().z() + (p1.z() - pts.back().z()) / dis_prev * step_size
                    );
                    pts.push_back(tmp_p);
                    dis_prev = sqrt(pow(p1.x() - pts.back().x(), 2) + pow(p1.y() - pts.back().y(), 2) + pow(p1.z() - pts.back().z(), 2));
                    extra_dis -= step_size;
                }
            } else if (dis_prev + extra_dis < step_size) {
                extra_dis += dis_prev;
                continue;
            } else {
                pts.push_back(p1);
                extra_dis = 0;
            }
        }
        
        // Handle end point properly
        {
            double dis1 = sqrt(pow(pts.back().x() - end_p.x(), 2) + pow(pts.back().y() - end_p.y(), 2) + pow(pts.back().z() - end_p.z(), 2));
            
            if (dis1 < step_size * 0.6) {
                if (pts.size() <= 1) {
                    // Do nothing
                } else {
                    double dis2 = sqrt(pow(pts.back().x() - pts.at(pts.size()-2).x(), 2) + 
                                     pow(pts.back().y() - pts.at(pts.size()-2).y(), 2) + 
                                     pow(pts.back().z() - pts.at(pts.size()-2).z(), 2));
                    double dis3 = (dis1 + dis2) / 2.0;
                    
                    WireCell::Point tmp_p(
                        pts.at(pts.size()-2).x() + (pts.back().x() - pts.at(pts.size()-2).x()) / dis2 * dis3,
                        pts.at(pts.size()-2).y() + (pts.back().y() - pts.at(pts.size()-2).y()) / dis2 * dis3,
                        pts.at(pts.size()-2).z() + (pts.back().z() - pts.at(pts.size()-2).z()) / dis2 * dis3
                    );
                    pts.pop_back();
                    pts.push_back(tmp_p);
                }
            } else if (dis1 > step_size * 1.6) {
                int npoints = std::round(dis1 / step_size);
                WireCell::Point p_save = pts.back();
                for (int j = 0; j + 1 < npoints; j++) {
                    WireCell::Point p(
                        p_save.x() + (end_p.x() - p_save.x()) / npoints * (j + 1),
                        p_save.y() + (end_p.y() - p_save.y()) / npoints * (j + 1),
                        p_save.z() + (end_p.z() - p_save.z()) / npoints * (j + 1)
                    );
                    pts.push_back(p);
                }
            }
            
            pts.push_back(end_p);
        }
        
        // Ensure there is no single point
        if (pts.size() == 1) {
            pts.push_back(end_p);
        }
        
        // Generate 2D projections and store fit points in the segment
        segment->fits(generate_fits_with_projections(segment, pts));
    }
}

void TrackFitting::organize_segments_path_2nd(double low_dis_limit, double end_point_limit){
    if (!m_graph) return;
    
    // First pass: check for vertices that are too close together
    check_and_reset_close_vertices();
    
    // Second pass: organize segments path with 2D projection
    auto edge_range = boost::edges(*m_graph);
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        auto& edge_bundle = (*m_graph)[*e_it];
        auto segment = edge_bundle.segment;
        if (!segment) continue;
        
        // Get ordered vertices
        std::shared_ptr<PR::Vertex> start_v, end_v;
        PR::node_descriptor vd1, vd2;
        if (!get_ordered_segment_vertices(segment, *e_it, start_v, end_v, vd1, vd2)) continue;
        
        // Check if vertices are endpoints (degree == 1)
        bool flag_startv_end = (boost::degree(vd1, *m_graph) == 1);
        bool flag_endv_end = (boost::degree(vd2, *m_graph) == 1);
        
        std::vector<WireCell::Point> pts, curr_pts;
        
        // Get current fitted path from the segment
        if (!segment->fits().empty()) {
            for (const auto& fit : segment->fits()) {
                curr_pts.push_back(fit.point);
            }
        } else {
            // If no fits, use original wcpts
            for (const auto& wcpt : segment->wcpts()) {
                curr_pts.push_back(wcpt.point);
            }
        }
        
        // Examine end points
        curr_pts = examine_end_ps_vec(segment, curr_pts, flag_startv_end, flag_endv_end);
        
        WireCell::Point start_p, end_p;
        
        // Process start vertex
        if (!start_v->fit().flag_fix) {
            start_p = curr_pts.front();
            
            if (flag_startv_end) {
                WireCell::Point p2 = curr_pts.front();
                double dis1 = 0;
                for (auto it = curr_pts.begin(); it != curr_pts.end(); it++) {
                    p2 = *it;
                    dis1 = sqrt(pow(start_p.x() - p2.x(), 2) + pow(start_p.y() - p2.y(), 2) + pow(start_p.z() - p2.z(), 2));
                    if (dis1 > low_dis_limit) break;
                }
                if (dis1 != 0) {
                    start_p = WireCell::Point(
                        start_p.x() + (start_p.x() - p2.x()) / dis1 * end_point_limit,
                        start_p.y() + (start_p.y() - p2.y()) / dis1 * end_point_limit,
                        start_p.z() + (start_p.z() - p2.z()) / dis1 * end_point_limit
                    );
                }
            }
            
            // Set fit point for start vertex
            PR::Fit start_fit = start_v->fit();
            start_fit.point = start_p;
            start_v->fit(start_fit);
        } else {
            start_p = start_v->fit().point;
        }
        
        // Process end vertex
        if (!end_v->fit().flag_fix) {
            end_p = curr_pts.back();
            
            if (flag_endv_end) {
                WireCell::Point p2 = curr_pts.back();
                double dis1 = 0;
                for (auto it = curr_pts.rbegin(); it != curr_pts.rend(); it++) {
                    p2 = *it;
                    dis1 = sqrt(pow(end_p.x() - p2.x(), 2) + pow(end_p.y() - p2.y(), 2) + pow(end_p.z() - p2.z(), 2));
                    if (dis1 > low_dis_limit) break;
                }
                if (dis1 != 0) {
                    end_p = WireCell::Point(
                        end_p.x() + (end_p.x() - p2.x()) / dis1 * end_point_limit,
                        end_p.y() + (end_p.y() - p2.y()) / dis1 * end_point_limit,
                        end_p.z() + (end_p.z() - p2.z()) / dis1 * end_point_limit
                    );
                }
            }
            
            // Set fit point for end vertex
            PR::Fit end_fit = end_v->fit();
            end_fit.point = end_p;
            end_v->fit(end_fit);
        } else {
            end_p = end_v->fit().point;
        }
        
        // Build the middle points
        pts.push_back(start_p);
        for (size_t i = 0; i != curr_pts.size(); i++) {
            WireCell::Point p1 = curr_pts.at(i);
            double dis = low_dis_limit;
            double dis1 = sqrt(pow(p1.x() - end_p.x(), 2) + pow(p1.y() - end_p.y(), 2) + pow(p1.z() - end_p.z(), 2));
            if (pts.size() > 0) {
                dis = sqrt(pow(p1.x() - pts.back().x(), 2) + pow(p1.y() - pts.back().y(), 2) + pow(p1.z() - pts.back().z(), 2));
            }
            
            if (dis1 < low_dis_limit * 0.8) {
                continue;
            } else if (dis < low_dis_limit * 0.8) {
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
        
        // Handle final connection to end point
        {
            double dis1 = sqrt(pow(pts.back().x() - end_p.x(), 2) + pow(pts.back().y() - end_p.y(), 2) + pow(pts.back().z() - end_p.z(), 2));
            if (dis1 < low_dis_limit * 0.2) {
                if (pts.size() > 1) pts.pop_back();
            } else if (dis1 > low_dis_limit * 1.6) {
                int npoints = std::round(dis1 / low_dis_limit);
                WireCell::Point p_save = pts.back();
                for (int j = 0; j + 1 < npoints; j++) {
                    WireCell::Point p(
                        p_save.x() + (end_p.x() - p_save.x()) / npoints * (j + 1),
                        p_save.y() + (end_p.y() - p_save.y()) / npoints * (j + 1),
                        p_save.z() + (end_p.z() - p_save.z()) / npoints * (j + 1)
                    );
                    pts.push_back(p);
                }
            }
            pts.push_back(end_p);
        }
        
        // Handle case where only one point exists
        if (pts.size() == 1) {
            pts.push_back(end_p);
        }
        
        // Generate 2D projections and store fit points in the segment
        segment->fits(generate_fits_with_projections(segment, pts));
    }
}


void TrackFitting::organize_segments_path(double low_dis_limit, double end_point_limit){
    if (!m_graph) return;
    
    // Iterate over all edges (segments) in the graph
    auto edge_range = boost::edges(*m_graph);
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        auto& edge_bundle = (*m_graph)[*e_it];
        auto segment = edge_bundle.segment;
        if (!segment) continue;
        
        // Get ordered vertices
        std::shared_ptr<PR::Vertex> start_v, end_v;
        PR::node_descriptor vd1, vd2;
        if (!get_ordered_segment_vertices(segment, *e_it, start_v, end_v, vd1, vd2)) continue;
        
        // Check if vertices are endpoints (degree == 1)
        bool flag_startv_end = (boost::degree(vd1, *m_graph) == 1);
        bool flag_endv_end = (boost::degree(vd2, *m_graph) == 1);
        
        std::vector<WireCell::Point> pts;
        std::vector<WireCell::Point> temp_wcps_vec;
        
        // Convert WCPoints to Points
        for (const auto& wcp : segment->wcpts()) {
            temp_wcps_vec.push_back(wcp.point);
        }
        
        WireCell::Point start_p, end_p;
        
        start_v->fit().index = -1;
        end_v->fit().index = -1;

        // Process start vertex
        if (!start_v->fit().flag_fix) {
            start_p = temp_wcps_vec.front();
            
            if (flag_startv_end) {
                WireCell::Point p2 = temp_wcps_vec.front();
                double dis1 = 0;
                for (auto it = temp_wcps_vec.begin(); it != temp_wcps_vec.end(); it++) {
                    p2 = *it;
                    dis1 = sqrt(pow(start_p.x() - p2.x(), 2) + pow(start_p.y() - p2.y(), 2) + pow(start_p.z() - p2.z(), 2));
                    if (dis1 > low_dis_limit) break;
                }
                if (dis1 != 0) {
                    start_p = WireCell::Point(
                        start_p.x() + (start_p.x() - p2.x()) / dis1 * end_point_limit,
                        start_p.y() + (start_p.y() - p2.y()) / dis1 * end_point_limit,
                        start_p.z() + (start_p.z() - p2.z()) / dis1 * end_point_limit
                    );
                }
            }
            
            // Set fit point for start vertex
            PR::Fit start_fit = start_v->fit();
            start_fit.point = start_p;
            start_v->fit(start_fit);
        } else {
            start_p = start_v->fit().point;
        }
        
        // Process end vertex
        if (!end_v->fit().flag_fix) {
            end_p = temp_wcps_vec.back();
            
            if (flag_endv_end) {
                WireCell::Point p2 = temp_wcps_vec.back();
                double dis1 = 0;
                for (auto it = temp_wcps_vec.rbegin(); it != temp_wcps_vec.rend(); it++) {
                    p2 = *it;
                    dis1 = sqrt(pow(end_p.x() - p2.x(), 2) + pow(end_p.y() - p2.y(), 2) + pow(end_p.z() - p2.z(), 2));
                    if (dis1 > low_dis_limit) break;
                }
                if (dis1 != 0) {
                    end_p = WireCell::Point(
                        end_p.x() + (end_p.x() - p2.x()) / dis1 * end_point_limit,
                        end_p.y() + (end_p.y() - p2.y()) / dis1 * end_point_limit,
                        end_p.z() + (end_p.z() - p2.z()) / dis1 * end_point_limit
                    );
                }
            }
            
            // Set fit point for end vertex
            PR::Fit end_fit = end_v->fit();
            end_fit.point = end_p;
            end_v->fit(end_fit);
        } else {
            end_p = end_v->fit().point;
        }
        
        // Build the middle points
        pts.push_back(start_p);
        for (size_t i = 0; i != temp_wcps_vec.size(); i++) {
            WireCell::Point p1 = temp_wcps_vec.at(i);
            double dis = low_dis_limit;
            double dis1 = sqrt(pow(p1.x() - end_p.x(), 2) + pow(p1.y() - end_p.y(), 2) + pow(p1.z() - end_p.z(), 2));
            if (pts.size() > 0) {
                dis = sqrt(pow(p1.x() - pts.back().x(), 2) + pow(p1.y() - pts.back().y(), 2) + pow(p1.z() - pts.back().z(), 2));
            }
            
            if (dis1 < low_dis_limit * 0.8) {
                continue;
            } else if (dis < low_dis_limit * 0.8) {
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
        
        // Handle final connection to end point
        {
            double dis1 = sqrt(pow(pts.back().x() - end_p.x(), 2) + pow(pts.back().y() - end_p.y(), 2) + pow(pts.back().z() - end_p.z(), 2));
            if (dis1 < low_dis_limit * 0.2) {
                if (pts.size() > 1) pts.pop_back();
            } else if (dis1 > low_dis_limit * 1.6) {
                int npoints = std::round(dis1 / low_dis_limit);
                WireCell::Point p_save = pts.back();
                for (int j = 0; j + 1 < npoints; j++) {
                    WireCell::Point p(
                        p_save.x() + (end_p.x() - p_save.x()) / npoints * (j + 1),
                        p_save.y() + (end_p.y() - p_save.y()) / npoints * (j + 1),
                        p_save.z() + (end_p.z() - p_save.z()) / npoints * (j + 1)
                    );
                    pts.push_back(p);
                }
            }
            pts.push_back(end_p);
        }
        
        // Generate 2D projections and store fit points in the segment
        segment->fits(generate_fits_with_projections(segment, pts));
    }



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

void TrackFitting::update_association(std::shared_ptr<PR::Segment> segment, PlaneData& temp_2dut, PlaneData& temp_2dvt, PlaneData& temp_2dwt){
    if (!m_graph || !segment) return;
    
    // Get cluster and transformation info
    auto cluster = segment->cluster();
    const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
    // double cluster_t0 = cluster->get_cluster_t0();
    
    // Collect all segments from the graph for comparison
    std::vector<std::shared_ptr<PR::Segment>> all_segments;
    auto edge_range = boost::edges(*m_graph);
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        auto& edge_bundle = (*m_graph)[*e_it];
        if (edge_bundle.segment && edge_bundle.segment != segment) {
            all_segments.push_back(edge_bundle.segment);
        }
    }
    
    // Process U plane (plane 0)
    std::set<Coord2D> save_2dut;
    for (auto it = temp_2dut.associated_2d_points.begin(); it != temp_2dut.associated_2d_points.end(); it++) {
        const auto& coord = *it;
        
        // Convert 2D coordinates to 3D point using existing geometry infrastructure
        // Get APA/face info from the coordinate
        int apa = coord.apa;
        int face = coord.face;
        
        // Create a test point in raw coordinates
        WirePlaneId wpid(kUlayer, face, apa);
        auto offset_it = wpid_offsets.find(WirePlaneId(kAllLayers, face, apa));
        auto slope_it = wpid_slopes.find(WirePlaneId(kAllLayers, face, apa));
        
        if (offset_it == wpid_offsets.end() || slope_it == wpid_slopes.end()) continue;
        
        auto offset_t = std::get<0>(offset_it->second);
        auto offset_u = std::get<1>(offset_it->second);
        auto slope_x = std::get<0>(slope_it->second);
        auto slope_yu = std::get<1>(slope_it->second).first;
        // auto slope_zu = std::get<1>(slope_it->second).second;
        
        // Convert wire/time coordinates to 3D position
        double raw_x = (coord.time - offset_t) / slope_x;
        double raw_y = (coord.wire - offset_u)/slope_yu;
        
        // Create a test point (simplified - may need proper wire geometry)
        WireCell::Point test_point(raw_x, raw_y, 0); // Y,Z will be determined by wire geometry
        
        // Get distances to main segment and all other segments
        auto main_distances = segment_get_closest_2d_distances(segment, test_point, apa, face, "fit");
        double min_dis_track = std::get<0>(main_distances); // U plane distance
        
        double min_dis1_track = 1e9;
        for (const auto& other_seg : all_segments) {
            auto other_distances = segment_get_closest_2d_distances(other_seg, test_point, apa, face, "fit");
            double temp_dis = std::get<0>(other_distances); // U plane distance
            if (temp_dis < min_dis1_track) {
                min_dis1_track = temp_dis;
            }
        }
        
        // Apply selection criteria
        if (min_dis_track < min_dis1_track || min_dis_track < 0.3 * units::cm) {
            save_2dut.insert(*it);
        }
    }
    
    // Process V plane (plane 1)
    std::set<Coord2D> save_2dvt;
    for (auto it = temp_2dvt.associated_2d_points.begin(); it != temp_2dvt.associated_2d_points.end(); it++) {
        const auto& coord = *it;
        
        int apa = coord.apa;
        int face = coord.face;
        
        WirePlaneId wpid(kVlayer, face, apa);
        auto offset_it = wpid_offsets.find(WirePlaneId(kAllLayers, face, apa));
        auto slope_it = wpid_slopes.find(WirePlaneId(kAllLayers, face, apa));
        
        if (offset_it == wpid_offsets.end() || slope_it == wpid_slopes.end()) continue;
        
        auto offset_t = std::get<0>(offset_it->second);
        auto offset_v = std::get<2>(offset_it->second);
        auto slope_x = std::get<0>(slope_it->second);
        auto slope_yv = std::get<2>(slope_it->second).first;
        
        double raw_x = (coord.time - offset_t) / slope_x;
        double raw_y = (coord.wire - offset_v)/slope_yv;
        WireCell::Point test_point(raw_x, raw_y, 0);
        
        auto main_distances = segment_get_closest_2d_distances(segment, test_point, apa, face, "fit");
        double min_dis_track = std::get<1>(main_distances); // V plane distance
        
        double min_dis1_track = 1e9;
        for (const auto& other_seg : all_segments) {
            auto other_distances = segment_get_closest_2d_distances(other_seg, test_point, apa, face, "fit");
            double temp_dis = std::get<1>(other_distances); // V plane distance
            if (temp_dis < min_dis1_track) {
                min_dis1_track = temp_dis;
            }
        }
        
        if (min_dis_track < min_dis1_track || min_dis_track < 0.3 * units::cm) {
            save_2dvt.insert(*it);
        }
    }
    
    // Process W plane (plane 2)
    std::set<Coord2D> save_2dwt;
    for (auto it = temp_2dwt.associated_2d_points.begin(); it != temp_2dwt.associated_2d_points.end(); it++) {
        const auto& coord = *it;
        
        int apa = coord.apa;
        int face = coord.face;
        
        WirePlaneId wpid(kWlayer, face, apa);
        auto offset_it = wpid_offsets.find(WirePlaneId(kAllLayers, face, apa));
        auto slope_it = wpid_slopes.find(WirePlaneId(kAllLayers, face, apa));
        
        if (offset_it == wpid_offsets.end() || slope_it == wpid_slopes.end()) continue;
        
        auto offset_t = std::get<0>(offset_it->second);
        auto offset_w = std::get<3>(offset_it->second);
        auto slope_x = std::get<0>(slope_it->second);
        auto slope_yw = std::get<3>(slope_it->second).first;
        
        double raw_x = (coord.time - offset_t) / slope_x;
        double raw_y = (coord.wire - offset_w)/slope_yw;
        WireCell::Point test_point(raw_x, raw_y, 0);
        
        auto main_distances = segment_get_closest_2d_distances(segment, test_point, apa, face, "fit");
        double min_dis_track = std::get<2>(main_distances); // W plane distance
        
        double min_dis1_track = 1e9;
        for (const auto& other_seg : all_segments) {
            auto other_distances = segment_get_closest_2d_distances(other_seg, test_point, apa, face, "fit");
            double temp_dis = std::get<2>(other_distances); // W plane distance
            if (temp_dis < min_dis1_track) {
                min_dis1_track = temp_dis;
            }
        }
        
        if (min_dis_track < min_dis1_track || min_dis_track < 0.3 * units::cm) {
            save_2dwt.insert(*it);
        }
    }
    
    // Update the input plane data with filtered results
    temp_2dut.associated_2d_points = save_2dut;
    temp_2dvt.associated_2d_points = save_2dvt;
    temp_2dwt.associated_2d_points = save_2dwt;
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

void TrackFitting::form_map_graph(bool flag_exclusion, double end_point_factor, double mid_point_factor, int nlevel, double time_tick_cut, double charge_cut){
    // Clear existing mappings
    m_3d_to_2d.clear();
    m_2d_to_3d.clear();

    // Reset fit properties for all vertices first
    for (auto vp = boost::vertices(*m_graph); vp.first != vp.second; ++vp.first) {
        auto vd = *vp.first;
        auto& v_bundle = (*m_graph)[vd];
        if (v_bundle.vertex) {
            bool flag_fix = v_bundle.vertex->flag_fix();
            v_bundle.vertex->reset_fit_prop();
            v_bundle.vertex->flag_fix(flag_fix);
        }
    }

    // Collect segments and reset their fit properties
    std::vector<std::shared_ptr<PR::Segment>> segments;
    auto edge_range = boost::edges(*m_graph);
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        auto& edge_bundle = (*m_graph)[*e_it];
        if (edge_bundle.segment) {
            segments.push_back(edge_bundle.segment);
            edge_bundle.segment->reset_fit_prop();
        }
    }

    int count = 0;
    
    // Process each segment
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        auto& edge_bundle = (*m_graph)[*e_it];
        if (!edge_bundle.segment) continue;
        
        auto segment = edge_bundle.segment;
        auto& fits = segment->fits();
        if (fits.empty()) continue;

        // Get start and end vertices for this segment
        auto vd1 = boost::source(*e_it, *m_graph);
        auto vd2 = boost::target(*e_it, *m_graph);
        auto& v_bundle1 = (*m_graph)[vd1];
        auto& v_bundle2 = (*m_graph)[vd2];
        
        std::shared_ptr<PR::Vertex> start_v = nullptr, end_v = nullptr;
        
        if (v_bundle1.vertex && v_bundle2.vertex) {
            // Determine which vertex is start and which is end by comparing with first/last fit points   
            double dist1_first = (v_bundle1.vertex->wcpt().point - segment->wcpts().front().point).magnitude();
            double dist1_last = (v_bundle1.vertex->wcpt().point - segment->wcpts().back().point).magnitude();

            if (dist1_first < dist1_last) {
                start_v = v_bundle1.vertex;
                end_v = v_bundle2.vertex;
            } else {
                start_v = v_bundle2.vertex;
                end_v = v_bundle1.vertex;
            }
        }

        // Calculate distances between consecutive fit points
        std::vector<double> distances;
        for (size_t i = 0; i + 1 < fits.size(); i++) {
            distances.push_back((fits[i+1].point - fits[i].point).magnitude());
        }

        // std::vector<WireCell::Point> saved_pts;
        // std::vector<int> saved_index;
        // std::vector<bool> saved_skip;
        std::vector<PR::Fit> saved_fits;

        // Process each fit point
        for (size_t i = 0; i < fits.size(); i++) {
            double dis_cut;
            
            if (i == 0) {
                dis_cut = std::min(distances.empty() ? 0.0 : distances[0] * end_point_factor, 
                                  4.0/3.0 * end_point_factor * units::cm);
                if (start_v && start_v->fit_range() < 0) {
                    start_v->fit_range(dis_cut);
                } else if (start_v && dis_cut < start_v->fit_range()) {
                    start_v->fit_range(dis_cut);
                }
            } else if (i + 1 == fits.size()) {
                dis_cut = std::min(distances.empty() ? 0.0 : distances.back() * end_point_factor, 
                                  4.0/3.0 * end_point_factor * units::cm);
                if (end_v && end_v->fit_range() < 0) {
                    end_v->fit_range(dis_cut);
                } else if (end_v && dis_cut < end_v->fit_range()) {
                    end_v->fit_range(dis_cut);
                }
            } else {
                double dist_prev = i > 0 ? distances[i-1] : 0.0;
                double dist_next = i < distances.size() ? distances[i] : 0.0;
                dis_cut = std::min(std::max(dist_prev * mid_point_factor, dist_next * mid_point_factor), 
                                  4.0/3.0 * mid_point_factor * units::cm);
            }

            // Not the first and last point - process middle points
            if (i != 0 && i + 1 != fits.size()) {
                TrackFitting::PlaneData temp_2dut, temp_2dvt, temp_2dwt;
                form_point_association(segment, fits[i].point, temp_2dut, temp_2dvt, temp_2dwt, dis_cut, nlevel, time_tick_cut);

                if (flag_exclusion) {
                    update_association(segment, temp_2dut, temp_2dvt, temp_2dwt);
                }

                // Examine point association
                bool is_end_point = (i == 1 || i + 2 == fits.size());
                examine_point_association(segment, fits[i].point, temp_2dut, temp_2dvt, temp_2dwt, is_end_point, charge_cut);

                if (temp_2dut.quantity + temp_2dvt.quantity + temp_2dwt.quantity > 0) {
                    // Store in mapping structures
                    m_3d_to_2d[count].set_plane_data(WirePlaneLayer_t::kUlayer, temp_2dut);
                    m_3d_to_2d[count].set_plane_data(WirePlaneLayer_t::kVlayer, temp_2dvt);
                    m_3d_to_2d[count].set_plane_data(WirePlaneLayer_t::kWlayer, temp_2dwt);

                    // Fill reverse mappings
                    for (const auto& coord : temp_2dut.associated_2d_points) {
                        m_2d_to_3d[coord].insert(count);
                    }
                    for (const auto& coord : temp_2dvt.associated_2d_points) {
                        m_2d_to_3d[coord].insert(count);
                    }
                    for (const auto& coord : temp_2dwt.associated_2d_points) {
                        m_2d_to_3d[coord].insert(count);
                    }

                    saved_fits.push_back(fits[i]);
                    saved_fits.back().index = count;
                    saved_fits.back().flag_fix = false;

                    count++;
                }
            } else if (i == 0) {
                // First point
                saved_fits.push_back(fits[i]);
                saved_fits.back().flag_fix = true;
                if (start_v->fit_index() == -1) {
                    start_v->fit_index(count);
                    saved_fits.back().index = count;
                    count++;
                }else{
                    saved_fits.back().index = start_v->fit_index();
                    // saved_index.push_back(start_v->fit_index());

                }
            } else if (i + 1 == fits.size()) {
                // Last point
                saved_fits.push_back(fits[i]);
                saved_fits.back().flag_fix = true;
                if (end_v->fit_index() == -1) {
                    end_v->fit_index(count);
                    saved_fits.back().index = count;
                    count++;
                }else{
                    saved_fits.back().index = end_v->fit_index();
                }
            }
        }

        // std::cout << "Form Map: "  << " " << saved_pts.size() << " " << m_2d_to_3d.size() << " " << m_3d_to_2d.size() << std::endl;


        // Set fit associate vector for the segment
        segment->set_fit_associate_vec(saved_fits, m_dv, "fit");
    }

    // Deal with all vertices again
    for (auto vp = boost::vertices(*m_graph); vp.first != vp.second; ++vp.first) {
        auto vd = *vp.first;
        auto& v_bundle = (*m_graph)[vd];
        if (!v_bundle.vertex) continue;
        
        auto vertex = v_bundle.vertex;
        double dis_cut = vertex->fit_range();
        int vertex_count = vertex->fit_index();
        
        if (dis_cut > 0 && vertex_count >= 0) {
            WireCell::Point pt = vertex->fit().point;

            TrackFitting::PlaneData temp_2dut, temp_2dvt, temp_2dwt;
            
            // For vertex, we need to pass a dummy segment - use first available segment
            std::shared_ptr<PR::Segment> dummy_segment = nullptr;
            if (!segments.empty()) {
                dummy_segment = segments[0];
            }
            
            if (dummy_segment) {
                form_point_association(dummy_segment, pt, temp_2dut, temp_2dvt, temp_2dwt, dis_cut, nlevel, time_tick_cut);
                examine_point_association(dummy_segment, pt, temp_2dut, temp_2dvt, temp_2dwt, true, charge_cut);

                // Store vertex associations
                m_3d_to_2d[vertex_count].set_plane_data(WirePlaneLayer_t::kUlayer, temp_2dut);
                m_3d_to_2d[vertex_count].set_plane_data(WirePlaneLayer_t::kVlayer, temp_2dvt);
                m_3d_to_2d[vertex_count].set_plane_data(WirePlaneLayer_t::kWlayer, temp_2dwt);

                // Fill reverse mappings
                for (const auto& coord : temp_2dut.associated_2d_points) {
                    m_2d_to_3d[coord].insert(vertex_count);
                }
                for (const auto& coord : temp_2dvt.associated_2d_points) {
                    m_2d_to_3d[coord].insert(vertex_count);
                }
                for (const auto& coord : temp_2dwt.associated_2d_points) {
                    m_2d_to_3d[coord].insert(vertex_count);
                }
            }
        }
    }

    //  std::cout << "Form Map: "  << " I " <<  " " << m_2d_to_3d.size() << " " << m_3d_to_2d.size() << std::endl;
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


WireCell::Point TrackFitting::fit_point(WireCell::Point& init_p, int i, std::shared_ptr<PR::Segment> segment, std::map<std::pair<int, int>, std::map<std::tuple<int, int, int>, double>>& map_Udiv_fac, std::map<std::pair<int, int>, std::map<std::tuple<int, int, int>, double>>& map_Vdiv_fac, std::map<std::pair<int, int>, std::map<std::tuple<int, int, int>, double>>& map_Wdiv_fac, double offset_t, double slope_x, double offset_u, double slope_yu, double slope_zu, double offset_v, double slope_yv, double slope_zv, double offset_w, double slope_yw, double slope_zw){
    // Check if the point index exists in the 3D to 2D mapping
    auto point_it = m_3d_to_2d.find(i);
    if (point_it == m_3d_to_2d.end()) {
        return init_p; // Return original point if no 2D associations found
    }
    
    const auto& point_info = point_it->second;
    
    // Get plane data for each wire plane
    auto plane_data_u = point_info.get_plane_data(WirePlaneLayer_t::kUlayer);
    auto plane_data_v = point_info.get_plane_data(WirePlaneLayer_t::kVlayer);
    auto plane_data_w = point_info.get_plane_data(WirePlaneLayer_t::kWlayer);
    
    int n_2D_u = 2 * plane_data_u.associated_2d_points.size();
    int n_2D_v = 2 * plane_data_v.associated_2d_points.size();
    int n_2D_w = 2 * plane_data_w.associated_2d_points.size();
    
    // Set up Eigen matrices and vectors
    Eigen::VectorXd temp_pos_3D(3), data_u_2D(n_2D_u), data_v_2D(n_2D_v), data_w_2D(n_2D_w);
    Eigen::VectorXd temp_pos_3D_init(3);
    Eigen::SparseMatrix<double> RU(n_2D_u, 3);
    Eigen::SparseMatrix<double> RV(n_2D_v, 3);
    Eigen::SparseMatrix<double> RW(n_2D_w, 3);
    
    auto test_wpid = m_dv->contained_by(init_p);
    auto cluster = segment->cluster();
    auto cluster_t0 = cluster->get_cluster_t0();
    const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
    // Initialization with its raw position
    auto p_raw = transform->backward(init_p, cluster_t0, test_wpid.face(), test_wpid.apa());

    // Initialize with input point
    temp_pos_3D_init(0) = p_raw.x();
    temp_pos_3D_init(1) = p_raw.y();
    temp_pos_3D_init(2) = p_raw.z();
    
    // Initialize data vectors to zero
    data_u_2D.setZero();
    data_v_2D.setZero();
    data_w_2D.setZero();
    
    // Fill U plane data
    int index = 0;
    for (auto it = plane_data_u.associated_2d_points.begin(); it != plane_data_u.associated_2d_points.end(); it++) {
        // Get charge measurement
        CoordReadout charge_key(it->apa, it->time, it->channel);
        double charge = m_params.default_charge_th, charge_err = m_params.default_charge_err; // Default values (100, 1000 from WCP)

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
        auto div_key = std::make_tuple(it->time, it->wire, i);
        auto div_it1 = map_Udiv_fac.find(apa_face_key);
        if (div_it1 != map_Udiv_fac.end()) { 
            auto div_it2 = div_it1->second.find(div_key);
            if (div_it2 != div_it1->second.end()) {
                div_factor = div_it2->second; 
            }
        }

        double scaling = (charge / charge_err) * div_factor;
        
        if (plane_data_u.quantity < m_params.scaling_quality_th) {
            if (plane_data_u.quantity != 0) {
                scaling *= pow(plane_data_u.quantity / m_params.scaling_quality_th, 1);
            } else {
                scaling *= m_params.scaling_ratio;
            }
        } 
       
        
        if (scaling != 0) {
            data_u_2D(2 * index) = scaling * (it->wire - offset_u);
            data_u_2D(2 * index + 1) = scaling * (it->time - offset_t);
            
            RU.insert(2 * index, 1) = scaling * slope_yu;     // Y --> U
            RU.insert(2 * index, 2) = scaling * slope_zu;     // Z --> U
            RU.insert(2 * index + 1, 0) = scaling * slope_x;  // X --> T
        }
        index++;
    }
    
    // Fill V plane data
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
        auto div_key = std::make_tuple(it->time, it->wire, i);
        auto div_it1 = map_Vdiv_fac.find(apa_face_key);
        if (div_it1 != map_Vdiv_fac.end()) { 
            auto div_it2 = div_it1->second.find(div_key);
            if (div_it2 != div_it1->second.end()) {
                div_factor = div_it2->second; 
            }
        }

        double scaling = (charge / charge_err) * div_factor;
        
        if (plane_data_v.quantity < m_params.scaling_quality_th) {
            if (plane_data_v.quantity != 0) {
                scaling *= pow(plane_data_v.quantity / m_params.scaling_quality_th, 1);
            } else {
                scaling *= m_params.scaling_ratio;
            }
        } 
        
        if (scaling != 0) {
            data_v_2D(2 * index) = scaling * (it->wire - offset_v);
            data_v_2D(2 * index + 1) = scaling * (it->time - offset_t);
            
            RV.insert(2 * index, 1) = scaling * slope_yv;     // Y --> V
            RV.insert(2 * index, 2) = scaling * slope_zv;     // Z --> V
            RV.insert(2 * index + 1, 0) = scaling * slope_x;  // X --> T
        }
        index++;
    }
    
    // Fill W plane data
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
        auto div_key = std::make_tuple(it->time, it->wire, i);
        auto div_it1 = map_Wdiv_fac.find(apa_face_key);
        if (div_it1 != map_Wdiv_fac.end()) { 
            auto div_it2 = div_it1->second.find(div_key);
            if (div_it2 != div_it1->second.end()) {
                div_factor = div_it2->second; 
            }
        }

        double scaling = (charge / charge_err) * div_factor;
        
        if (plane_data_w.quantity < m_params.scaling_quality_th) {
            if (plane_data_w.quantity != 0) {
                scaling *= pow(plane_data_w.quantity / m_params.scaling_quality_th, 1);
            } else {
                scaling *= m_params.scaling_ratio;
            }
        } 
        
        if (scaling != 0) {
            data_w_2D(2 * index) = scaling * (it->wire - offset_w);
            data_w_2D(2 * index + 1) = scaling * (it->time - offset_t);
            
            RW.insert(2 * index, 2) = scaling * slope_zw;     // Z --> W
            RW.insert(2 * index + 1, 0) = scaling * slope_x;  // X --> T
        }
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
    
    WireCell::Point final_p_raw, final_p;
    
    // Check if solver succeeded
    if (std::isnan(solver.error())) {
        // Return initialization if solver failed
        final_p_raw = WireCell::Point(temp_pos_3D_init(0), temp_pos_3D_init(1), temp_pos_3D_init(2));
    } else {
        // Return fitted result
        final_p_raw = WireCell::Point(temp_pos_3D(0), temp_pos_3D(1), temp_pos_3D(2));
    }

    final_p = transform->forward(final_p_raw, cluster_t0, test_wpid.face(), test_wpid.apa());

    return final_p;
}

void TrackFitting::multi_trajectory_fit(int charge_div_method, double div_sigma){
    if (!m_graph) return;
    
    // create pss_vec ... 
    std::map<int, std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> map_index_pss;
    auto edge_range = boost::edges(*m_graph);
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        auto& edge_bundle = (*m_graph)[*e_it];
        if (!edge_bundle.segment) continue;
        auto segment = edge_bundle.segment;
        const auto& fits = segment->fits();
        for (const auto& fit : fits) {
            int idx = fit.index;
            map_index_pss[idx] = std::make_pair(fit.point, segment);
        }
    }


    // Create charge division factor maps for each APA/face combination
    // Structure: [apa/face] -> [time/wire/3D_index] -> factor
    std::map<std::pair<int, int>, std::map<std::tuple<int, int, int>, double>> map_Udiv_fac;
    std::map<std::pair<int, int>, std::map<std::tuple<int, int, int>, double>> map_Vdiv_fac;
    std::map<std::pair<int, int>, std::map<std::tuple<int, int, int>, double>> map_Wdiv_fac;

    // Equal division - same as in trajectory_fit
    for (auto it = m_2d_to_3d.begin(); it != m_2d_to_3d.end(); it++) {
        for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++) {
            WirePlaneLayer_t plane = it->first.plane;
            int apa = it->first.apa;
            int face = it->first.face;
            int time = it->first.time;
            int wire = it->first.wire;
            
            auto apa_face_key = std::make_pair(apa, face);
            auto div_key = std::make_tuple(time, wire, *it1);
            
            if (plane == WirePlaneLayer_t::kUlayer) {
                map_Udiv_fac[apa_face_key][div_key] = 1.0 / it->second.size();
            } else if (plane == WirePlaneLayer_t::kVlayer) {
                map_Vdiv_fac[apa_face_key][div_key] = 1.0 / it->second.size();
            } else if (plane == WirePlaneLayer_t::kWlayer) {
                map_Wdiv_fac[apa_face_key][div_key] = 1.0 / it->second.size();
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

                    auto pss_vec_it = map_index_pss.find(idx);
                    if (pss_vec_it == map_index_pss.end()) continue;

                    auto segment = pss_vec_it->second.second;
                    auto cluster = segment->cluster();
                    const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
                    double cluster_t0 = cluster->get_cluster_t0();

                    auto test_wpid = m_dv->contained_by(pss_vec_it->second.first);
                    // Initialization with its raw position
                    auto p_raw = transform->backward(pss_vec_it->second.first, cluster_t0, test_wpid.face(), test_wpid.apa());

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
    
    // Process vertices first
    for (auto vp = boost::vertices(*m_graph); vp.first != vp.second; ++vp.first) {
        auto vd = *vp.first;
        auto& v_bundle = (*m_graph)[vd];
        if (!v_bundle.vertex) continue;
        
        auto vertex = v_bundle.vertex;
        int i = vertex->fit_index();
        auto segment = map_index_pss[i].second;
        auto cluster = segment->cluster();
        auto cluster_t0 = cluster->get_cluster_t0();
        const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));

        bool flag_fit_fix = vertex->flag_fix();
        WireCell::Point init_p = vertex->fit().point;
        
        if (!flag_fit_fix && i >= 0) {
            // Get geometry parameters for this point
            auto test_wpid = m_dv->contained_by(init_p);
            if (test_wpid.apa() == -1 || test_wpid.face() == -1) continue;
            
            WirePlaneId wpid(kAllLayers, test_wpid.face(), test_wpid.apa());
            auto offset_it = wpid_offsets.find(wpid);
            auto slope_it = wpid_slopes.find(wpid);
            
            if (offset_it == wpid_offsets.end() || slope_it == wpid_slopes.end()) continue;
            
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
            
            // Fit the vertex point
            WireCell::Point fitted_p = fit_point(init_p, i, segment, map_Udiv_fac, map_Vdiv_fac, map_Wdiv_fac,
                                                offset_t, slope_x, offset_u, slope_yu, slope_zu,
                                                offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw);
            
            // std::cout << "Vertex Fit " << i << " " << init_p << " --> " << fitted_p << std::endl;

            // Update vertex with fitted position and projections
            auto& vertex_fit = vertex->fit();
            vertex_fit.point = fitted_p;
            vertex_fit.index = -1;

            // Store 2D projections (following WCP pattern with offsets)
            auto fitted_p_raw = transform->backward(fitted_p, cluster_t0, test_wpid.apa(), test_wpid.face());
            vertex_fit.pu = offset_u  + (slope_yu * fitted_p_raw.y() + slope_zu * fitted_p_raw.z());
            vertex_fit.pv = offset_v  + (slope_yv * fitted_p_raw.y() + slope_zv * fitted_p_raw.z());
            vertex_fit.pw = offset_w  + (slope_yw * fitted_p_raw.y() + slope_zw * fitted_p_raw.z());
            vertex_fit.pt = offset_t  + slope_x * fitted_p_raw.x();

            vertex_fit.paf = std::make_pair(test_wpid.apa(), test_wpid.face());
        }
    }
    
    // Process segments (tracks)
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        auto& edge_bundle = (*m_graph)[*e_it];
        if (!edge_bundle.segment) continue;
        
        auto segment = edge_bundle.segment;
        auto cluster = segment->cluster();
        auto cluster_t0 = cluster->get_cluster_t0();
        const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
        auto& fits = segment->fits();
        if (fits.empty()) continue;
        
        // Get start and end vertices
        auto vd1 = boost::source(*e_it, *m_graph);
        auto vd2 = boost::target(*e_it, *m_graph);
        auto& v_bundle1 = (*m_graph)[vd1];
        auto& v_bundle2 = (*m_graph)[vd2];
        
        std::shared_ptr<PR::Vertex> start_v = nullptr, end_v = nullptr;
        if (v_bundle1.vertex && v_bundle2.vertex) {
            // Determine which vertex is start and which is end by comparing with first/last fit points   
            double dist1_first = (v_bundle1.vertex->wcpt().point - segment->wcpts().front().point).magnitude();
            double dist1_last = (v_bundle1.vertex->wcpt().point - segment->wcpts().back().point).magnitude();

            if (dist1_first < dist1_last) {
                start_v = v_bundle1.vertex;
                end_v = v_bundle2.vertex;
            } else {
                start_v = v_bundle2.vertex;
                end_v = v_bundle1.vertex;
            }
        }
        
        // Get initial points and fit info
        std::vector<WireCell::Point> init_ps;
        std::vector<int> init_indices;
        std::vector<bool> init_fit_skip;
        
        for (const auto& fit : fits) {
            init_ps.push_back(fit.point);
            init_indices.push_back(fit.index);
            init_fit_skip.push_back(fit.flag_fix);
        }
        
        if (init_ps.empty()) continue;
        
        // Set endpoint positions from fitted vertices
        if (start_v) init_ps.front() = start_v->fit().point;
        if (end_v) init_ps.back() = end_v->fit().point;
        
        std::vector<WireCell::Point> final_ps;
        
        // Fit each point in the segment
        for (size_t i = 0; i < init_ps.size(); i++) {
            if (i == 0) {
                // Use start vertex position
                final_ps.push_back(start_v ? start_v->fit().point : init_ps[i]);
            } else if (i + 1 == init_ps.size()) {
                // Use end vertex position  
                final_ps.push_back(end_v ? end_v->fit().point : init_ps[i]);
            } else {
                WireCell::Point temp_p = init_ps[i];
                
                if (!init_fit_skip[i] && init_indices[i] >= 0) {
                    // Get geometry parameters
                    auto test_wpid = m_dv->contained_by(init_ps[i]);
                    if (test_wpid.apa() != -1 && test_wpid.face() != -1) {
                        WirePlaneId wpid(kAllLayers, test_wpid.face(), test_wpid.apa());
                        auto offset_it = wpid_offsets.find(wpid);
                        auto slope_it = wpid_slopes.find(wpid);
                        
                        if (offset_it != wpid_offsets.end() && slope_it != wpid_slopes.end()) {
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
                            
                            // Fit the point
                            temp_p = fit_point(init_ps[i], init_indices[i], segment, map_Udiv_fac, map_Vdiv_fac, map_Wdiv_fac,
                                             offset_t, slope_x, offset_u, slope_yu, slope_zu,
                                             offset_v, slope_yv, slope_zv, offset_w, slope_yw, slope_zw);
                        }
                    }
                }
                // std::cout << "Fit point Track " << i << ": (" << init_ps[i].x() << ", " << init_ps[i].y() << ", " << init_ps[i].z() << ")" <<  " : (" << temp_p.x() << ", " << temp_p.y() << ", " << temp_p.z() << ")" << std::endl;
                final_ps.push_back(temp_p);
            }
        }
        
        // Apply trajectory examination/smoothing
        std::vector<WireCell::Point> examined_ps = examine_segment_trajectory(segment, final_ps, init_ps);
        
        // std::cout  << " fitted with " << examined_ps.size() << " " << init_ps.size() << " " << final_ps.size() << " points." << std::endl;

        // Update segment with fitted results
        std::vector<PR::Fit> new_fits;
        for (size_t i = 0; i < examined_ps.size(); i++) {
            PR::Fit fit;
            fit.point = examined_ps[i];
            fit.index = -1;
            fit.flag_fix = false;
            
            // Calculate 2D projections
            auto test_wpid = m_dv->contained_by(examined_ps[i]);


            if (test_wpid.apa() != -1 && test_wpid.face() != -1) {
                WirePlaneId wpid(kAllLayers, test_wpid.face(), test_wpid.apa());
                auto offset_it = wpid_offsets.find(wpid);
                auto slope_it = wpid_slopes.find(wpid);
                
                auto examined_p_raw = transform->backward(examined_ps[i], cluster_t0, test_wpid.apa(), test_wpid.face());
                fit.paf = std::make_pair(test_wpid.apa(), test_wpid.face());

                if (offset_it != wpid_offsets.end() && slope_it != wpid_slopes.end()) {
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

                    fit.pu = offset_u + (slope_yu * examined_p_raw.y() + slope_zu * examined_p_raw.z());
                    fit.pv = offset_v + (slope_yv * examined_p_raw.y() + slope_zv * examined_p_raw.z());
                    fit.pw = offset_w + (slope_yw * examined_p_raw.y() + slope_zw * examined_p_raw.z());
                    fit.pt = offset_t + slope_x * examined_p_raw.x();
                }
            }
            
            new_fits.push_back(fit);
        }
        
        // Update segment fits
        segment->fits() = new_fits;
    }
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
            if (area1 > m_params.area_ratio1 * c && area1 > m_params.area_ratio2 * area2) flag_replace = true;	
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
            
            if (area1 > m_params.area_ratio1 * c && area1 > m_params.area_ratio2 * area2) flag_replace = true;
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

std::vector<WireCell::Point> TrackFitting::examine_segment_trajectory(std::shared_ptr<PR::Segment> segment, std::vector<WireCell::Point>& final_ps_vec, std::vector<WireCell::Point>& init_ps_vec){
    // Create local trajectory data structures
    std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> pss_vec;
    std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> fine_tracking_path;
    std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> temp_fine_tracking_path;
    std::vector<std::pair<int, int>> saved_paf;
    
    // Initialize pss_vec from input vectors
    if (final_ps_vec.size() != init_ps_vec.size()) {
        return std::vector<WireCell::Point>(); // Return empty if sizes don't match
    }
    
    for (size_t i = 0; i < final_ps_vec.size(); i++) {
        // std::cout << i << " " << final_ps_vec[i].x() << " " << final_ps_vec[i].y() << " " << final_ps_vec[i].z() << " : " << init_ps_vec[i].x() << " " << init_ps_vec[i].y() << " " << init_ps_vec[i].z() << std::endl;

        pss_vec.push_back(std::make_pair(final_ps_vec[i], segment));
    }
    
    // First pass: apply skip_trajectory_point logic
    int skip_count = 0;
    for (size_t i = 0; i < pss_vec.size(); i++) {
        WireCell::Point p = final_ps_vec[i];
        
        // Get APA and face information
        auto test_wpid = m_dv->contained_by(p);
        auto apa_face = std::make_pair(test_wpid.apa(), test_wpid.face());
        
        // Apply skip trajectory point check
        bool flag_skip = skip_trajectory_point(p, apa_face, i, pss_vec, fine_tracking_path);
        
        // Vertex points (first and last) should not be skipped
        if (i == 0 || i + 1 == final_ps_vec.size()) {
            flag_skip = false;
        }
        
        // Protection against too many consecutive skips
        if (flag_skip) {
            skip_count++;
            if (skip_count <= 3) {
                continue;
            } else {
                skip_count = 0;
            }
        }
        
        // Store points for trajectory smoothing
        temp_fine_tracking_path.push_back(std::make_pair(init_ps_vec[i], segment));
        fine_tracking_path.push_back(std::make_pair(p, segment));
        saved_paf.push_back(std::make_pair(test_wpid.apa(), test_wpid.face()));
    }
    
    // Second pass: Apply trajectory smoothing (area-based correction)
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

        // -2, +1 neighbor check
        if ((!flag_replace) && i >= 2 && i + 1 != fine_tracking_path.size()) {
            double a = sqrt(pow(fine_tracking_path[i-2].first.x() - fine_tracking_path[i].first.x(), 2) +
                          pow(fine_tracking_path[i-2].first.y() - fine_tracking_path[i].first.y(), 2) +
                          pow(fine_tracking_path[i-2].first.z() - fine_tracking_path[i].first.z(), 2));
            double b = sqrt(pow(fine_tracking_path[i+1].first.x() - fine_tracking_path[i].first.x(), 2) +
                          pow(fine_tracking_path[i+1].first.y() - fine_tracking_path[i].first.y(), 2) +
                          pow(fine_tracking_path[i+1].first.z() - fine_tracking_path[i].first.z(), 2));
            double c = sqrt(pow(fine_tracking_path[i-2].first.x() - fine_tracking_path[i+1].first.x(), 2) +
                          pow(fine_tracking_path[i-2].first.y() - fine_tracking_path[i+1].first.y(), 2) +
                          pow(fine_tracking_path[i-2].first.z() - fine_tracking_path[i+1].first.z(), 2));
            double s = (a + b + c) / 2.0;
            double area1 = sqrt(s * (s - a) * (s - b) * (s - c));

            a = sqrt(pow(fine_tracking_path[i-2].first.x() - temp_fine_tracking_path[i].first.x(), 2) +
                    pow(fine_tracking_path[i-2].first.y() - temp_fine_tracking_path[i].first.y(), 2) +
                    pow(fine_tracking_path[i-2].first.z() - temp_fine_tracking_path[i].first.z(), 2));
            b = sqrt(pow(fine_tracking_path[i+1].first.x() - temp_fine_tracking_path[i].first.x(), 2) +
                    pow(fine_tracking_path[i+1].first.y() - temp_fine_tracking_path[i].first.y(), 2) +
                    pow(fine_tracking_path[i+1].first.z() - temp_fine_tracking_path[i].first.z(), 2));
            s = (a + b + c) / 2.0;
            double area2 = sqrt(s * (s - a) * (s - b) * (s - c));
            
            if (area1 > m_params.area_ratio1 * c && area1 > m_params.area_ratio2 * area2) {
                flag_replace = true;	
            }
        }
        
        // -1, +2 neighbor check
        if ((!flag_replace) && i > 0 && i + 2 < fine_tracking_path.size()) {
            double a = sqrt(pow(fine_tracking_path[i-1].first.x() - fine_tracking_path[i].first.x(), 2) +
                          pow(fine_tracking_path[i-1].first.y() - fine_tracking_path[i].first.y(), 2) +
                          pow(fine_tracking_path[i-1].first.z() - fine_tracking_path[i].first.z(), 2));
            double b = sqrt(pow(fine_tracking_path[i+2].first.x() - fine_tracking_path[i].first.x(), 2) +
                          pow(fine_tracking_path[i+2].first.y() - fine_tracking_path[i].first.y(), 2) +
                          pow(fine_tracking_path[i+2].first.z() - fine_tracking_path[i].first.z(), 2));
            double c = sqrt(pow(fine_tracking_path[i-1].first.x() - fine_tracking_path[i+2].first.x(), 2) +
                          pow(fine_tracking_path[i-1].first.y() - fine_tracking_path[i+2].first.y(), 2) +
                          pow(fine_tracking_path[i-1].first.z() - fine_tracking_path[i+2].first.z(), 2));
            double s = (a + b + c) / 2.0;
            double area1 = sqrt(s * (s - a) * (s - b) * (s - c));

            a = sqrt(pow(fine_tracking_path[i-1].first.x() - temp_fine_tracking_path[i].first.x(), 2) +
                    pow(fine_tracking_path[i-1].first.y() - temp_fine_tracking_path[i].first.y(), 2) +
                    pow(fine_tracking_path[i-1].first.z() - temp_fine_tracking_path[i].first.z(), 2));
            b = sqrt(pow(fine_tracking_path[i+2].first.x() - temp_fine_tracking_path[i].first.x(), 2) +
                    pow(fine_tracking_path[i+2].first.y() - temp_fine_tracking_path[i].first.y(), 2) +
                    pow(fine_tracking_path[i+2].first.z() - temp_fine_tracking_path[i].first.z(), 2));
            s = (a + b + c) / 2.0;
            double area2 = sqrt(s * (s - a) * (s - b) * (s - c));
            
            if (area1 > m_params.area_ratio1 * c && area1 > m_params.area_ratio2 * area2) {
                flag_replace = true;
            }
        }
        
        // Replace with original point if flagged
        if (flag_replace) {
            fine_tracking_path[i] = temp_fine_tracking_path[i];
        }
    }
    
    // Extract the final trajectory points
    std::vector<WireCell::Point> result_ps;
    for (const auto& point_pair : fine_tracking_path) {
        result_ps.push_back(point_pair.first);
    }
    
    return result_ps;
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

std::vector<std::vector<double>> TrackFitting::calculate_compact_matrix_multi(std::vector<std::vector<int> >& connected_vec,Eigen::SparseMatrix<double>& weight_matrix, const Eigen::SparseMatrix<double>& response_matrix_transpose, int n_2d_measurements, int n_3d_positions, double cut_position){
    // Initialize results vector - returns sharing ratios for each 3D position
    std::vector<std::vector<double>> results(n_3d_positions);

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
    
    // Calculate sharing ratios between connected 3D positions (key difference from regular version)
    for (auto it = map_3d_to_2d.begin(); it != map_3d_to_2d.end(); ++it) {
        int row = it->first;
        
        // Skip if row is out of bounds for connected_vec
        if (row >= static_cast<int>(connected_vec.size())) continue;
        
        // For each connected neighbor defined in connected_vec
        for (size_t i = 0; i < connected_vec.at(row).size(); i++) {
            double sum[2] = {0.0, 0.0};
            
            // Find the connected neighbor
            auto it1 = map_3d_to_2d.find(connected_vec.at(row).at(i));
            
            // Count total connections for current 3D position
            for (auto it3 = it->second.begin(); it3 != it->second.end(); ++it3) {
                sum[0] += 1.0;  // Total count (using 1 instead of val as in WCP)
            }
            
            // Count shared connections with this connected neighbor
            if (it1 != map_3d_to_2d.end()) {
                std::vector<int> common_results(it->second.size());
                auto it3 = std::set_intersection(
                    it->second.begin(), it->second.end(),
                    it1->second.begin(), it1->second.end(),
                    common_results.begin()
                );
                common_results.resize(it3 - common_results.begin());
                
                for (auto it4 = common_results.begin(); it4 != common_results.end(); ++it4) {
                    sum[1] += 1.0;  // Shared count (using 1 instead of val as in WCP)
                }
            }
            
            // Calculate sharing ratio for this connected neighbor
            results.at(row).push_back(sum[1] / (sum[0] + 1e-9));
        }
    }
    
    // Ensure all result vectors have the correct size
    for (size_t i = 0; i < results.size(); i++) {
        if (i < connected_vec.size() && results.at(i).size() != connected_vec.at(i).size()) {
            results.at(i).resize(connected_vec.at(i).size(), 0.0);
        }
    }
    
    // Return the 2D vector directly (consistent with WCPPID)
    return results;
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

            // std::cout << "Row: " << row << ", Col: " << col << ", Value: " << value << std::endl;

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
        // std::cout << row << " " << sum1 << " " << sum2 << " " << flag << std::endl;
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
        
        // std::cout << row << " " << sum[0] << " " << sum[1] << " " << sum[2] << std::endl;

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

void TrackFitting::dQ_dx_multi_fit(double dis_end_point_ext, bool flag_dQ_dx_fit_reg){
    if (!m_graph) return;
    
    // Clear output vectors
    dQ.clear();
    dx.clear();
    reduced_chi2.clear();

    // Update charge data for shared wires
    update_dQ_dx_data();
    
    // Use parameters from member variable
    const double DL = m_params.DL;
    const double DT = m_params.DT;
    const double col_sigma_w_T = m_params.col_sigma_w_T;
    const double ind_sigma_u_T = m_params.ind_sigma_u_T;
    const double ind_sigma_v_T = m_params.ind_sigma_v_T;
    const double rel_uncer_ind = m_params.rel_uncer_ind;
    const double rel_uncer_col = m_params.rel_uncer_col;
    const double add_uncer_ind = m_params.add_uncer_ind;
    const double add_uncer_col = m_params.add_uncer_col;
    const double add_sigma_L = m_params.add_sigma_L;
    
    // Prepare charge data maps similar to dQ_dx_fit
    std::map<CoordReadout, std::pair<ChargeMeasurement, std::set<Coord2D>>> map_U_charge_2D, map_V_charge_2D, map_W_charge_2D;
    
    // Fill the maps from m_charge_data
    for (const auto& [coord_readout, charge_measurement] : m_charge_data) {
        int apa = coord_readout.apa;
        int time = coord_readout.time;
        int channel = coord_readout.channel;
        
        auto wires_info = get_wires_for_channel(apa, channel);
        if (wires_info.empty()) continue;
        
        std::set<TrackFitting::Coord2D> associated_coords;
        int plane = -1;
        
        for (const auto& wire_info : wires_info) {
            int face = std::get<0>(wire_info);
            plane = std::get<1>(wire_info);
            int wire = std::get<2>(wire_info);
            
            WirePlaneLayer_t plane_layer = (plane == 0) ? kUlayer : 
                                          (plane == 1) ? kVlayer : kWlayer;
            
            TrackFitting::Coord2D coord_2d(apa, face, time, wire, channel, plane_layer);
            associated_coords.insert(coord_2d);
        }
        
        std::pair<ChargeMeasurement, std::set<TrackFitting::Coord2D>> charge_coord_pair = 
            std::make_pair(charge_measurement, associated_coords);
        
        switch (plane) {
            case 0: map_U_charge_2D[coord_readout] = charge_coord_pair; break;
            case 1: map_V_charge_2D[coord_readout] = charge_coord_pair; break;
            case 2: map_W_charge_2D[coord_readout] = charge_coord_pair; break;
        }
    }
    
    // Count total 3D positions from all segments and vertices
    int n_3D_pos = 0;
    // reset fit_index for all vertices and segment points
    
    std::map<std::shared_ptr<PR::Vertex>, int> vertex_index_map;
    std::map<std::pair<std::shared_ptr<PR::Segment>, int>, int> segment_point_index_map;
    
    // First pass: assign indices to vertices and segments
    auto edge_range = boost::edges(*m_graph);
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        auto& edge_bundle = (*m_graph)[*e_it];
        if (!edge_bundle.segment) continue;
        
        auto segment = edge_bundle.segment;
        auto& fits = segment->fits();
        if (fits.empty()) continue;
        
        // Get start and end vertices
        auto vd1 = boost::source(*e_it, *m_graph);
        auto vd2 = boost::target(*e_it, *m_graph);
        auto& v_bundle1 = (*m_graph)[vd1];
        auto& v_bundle2 = (*m_graph)[vd2];
        
        std::shared_ptr<PR::Vertex> start_v = nullptr, end_v = nullptr;
        // std::cout << "Check: " << v_bundle1.vertex->fit().index << " " << v_bundle2.vertex->fit().index << " " << fits.front().index << " " << fits.back().index << " " << fits.size() << std::endl;
        if (v_bundle1.vertex && v_bundle2.vertex) {
            if (v_bundle1.vertex->fit().index == fits.front().index) {
                start_v = v_bundle1.vertex;
                end_v = v_bundle2.vertex;
            } else {
                start_v = v_bundle2.vertex;
                end_v = v_bundle1.vertex;
            }
        }
        
        // Assign indices to points
        for (size_t i = 0; i < fits.size(); i++) {
            int idx = fits[i].index;
            if (i == 0) {
                // Start vertex
                if (start_v && vertex_index_map.find(start_v) == vertex_index_map.end()) {
                    vertex_index_map[start_v] = start_v->fit().index;
                    idx = start_v->fit().index;
                }
                if (start_v) {
                    segment_point_index_map[std::make_pair(segment, i)] = vertex_index_map[start_v];
                } else {
                    segment_point_index_map[std::make_pair(segment, i)] = idx;
                }
            } else if (i + 1 == fits.size()) {
                // End vertex
                if (end_v && vertex_index_map.find(end_v) == vertex_index_map.end()) {
                    vertex_index_map[end_v] = end_v->fit().index;
                    idx = end_v->fit().index;
                }
                if (end_v) {
                    segment_point_index_map[std::make_pair(segment, i)] = vertex_index_map[end_v];
                } else {
                    segment_point_index_map[std::make_pair(segment, i)] = idx;
                }
            } else {
                // Middle points
                segment_point_index_map[std::make_pair(segment, i)] = idx;
            }
            // Track maximum index to determine n_3D_pos
            if (idx + 1 > n_3D_pos) n_3D_pos = idx + 1;
        }
    }
    
    // std::cout << n_3D_pos << " 3D positions identified for dQ/dx fitting." << std::endl;

    if (n_3D_pos == 0) return;
    
    int n_2D_u = map_U_charge_2D.size();
    int n_2D_v = map_V_charge_2D.size();
    int n_2D_w = map_W_charge_2D.size();
    
    if (n_2D_u == 0 && n_2D_v == 0 && n_2D_w == 0) return;
    
    // Initialize Eigen matrices and vectors
    Eigen::VectorXd pos_3D(n_3D_pos), data_u_2D(n_2D_u), data_v_2D(n_2D_v), data_w_2D(n_2D_w);
    Eigen::VectorXd pred_data_u_2D(n_2D_u), pred_data_v_2D(n_2D_v), pred_data_w_2D(n_2D_w);
    Eigen::SparseMatrix<double> RU(n_2D_u, n_3D_pos);
    Eigen::SparseMatrix<double> RV(n_2D_v, n_3D_pos);
    Eigen::SparseMatrix<double> RW(n_2D_w, n_3D_pos);
    
    std::vector<WireCell::Point> traj_pts(n_3D_pos);
    std::vector<double> local_dx(n_3D_pos, 0);
    std::vector<double> traj_reduced_chi2(n_3D_pos, 0);
    std::vector<int> reg_flag_u(n_3D_pos, 0), reg_flag_v(n_3D_pos, 0), reg_flag_w(n_3D_pos, 0);
    
    // Initialize solution vector
    Eigen::VectorXd pos_3D_init(n_3D_pos);
    for (int i = 0; i < n_3D_pos; i++) {
        pos_3D_init(i) = 5000.0*6; // Initial guess for single MIP
    }
    
    // Fill data vectors with charge/uncertainty ratios
    {
        int n_u = 0;
        for (const auto& [coord_key, result] : map_U_charge_2D) {
            const auto& measurement = result.first;
            if (measurement.charge > 0) {
                double charge = measurement.charge;
                double charge_err = measurement.charge_err;
                double total_err = sqrt(pow(charge_err, 2) + pow(charge * rel_uncer_ind, 2) + pow(add_uncer_ind, 2));
                data_u_2D(n_u) = charge / total_err;
            } else {
                data_u_2D(n_u) = 0;
            }
            n_u++;
        }
        
        int n_v = 0;
        for (const auto& [coord_key, result] : map_V_charge_2D) {
            const auto& measurement = result.first;
            if (measurement.charge > 0) {
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
            if (measurement.charge > 0) {
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

    // Build response matrices using cal_gaus_integral_seg
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        auto& edge_bundle = (*m_graph)[*e_it];
        if (!edge_bundle.segment) continue;
        
        auto segment = edge_bundle.segment;
        auto& fits = segment->fits();
        if (fits.empty()) continue;
        
        // Get time ticks from cluster
        int cur_ntime_ticks = 10; // Default value, should be calculated from cluster
        if (edge_bundle.segment && edge_bundle.segment->cluster()) {
            auto cluster = edge_bundle.segment->cluster();
            auto first_blob = cluster->children()[0];
            cur_ntime_ticks = first_blob->slice_index_max() - first_blob->slice_index_min();
        }

          // Fill trajectory points
        for (size_t i = 0; i < fits.size(); i++) {
            int idx = segment_point_index_map[std::make_pair(segment, i)];
            traj_pts[idx] = fits[i].point;
        }
        
        // Calculate dx values for middle points
        for (size_t i = 1; i + 1 < fits.size(); i++) {
            int idx = segment_point_index_map[std::make_pair(segment, i)];
            
            WireCell::Point prev_pos = fits[i-1].point;
            WireCell::Point curr_pos = fits[i].point;
            WireCell::Point next_pos = fits[i+1].point;
            
            WireCell::Point prev_mid = 0.5 * (prev_pos + curr_pos);
            WireCell::Point next_mid = 0.5 * (next_pos + curr_pos);
            
            double dx = (curr_pos - prev_mid).magnitude() + (curr_pos - next_mid).magnitude();
            local_dx[idx] = dx;
        }

        // Process middle points
        for (size_t i = 1; i + 1 < fits.size(); i++) {
            int idx = segment_point_index_map[std::make_pair(segment, i)];
            
            WireCell::Point prev_pos = fits[i-1].point;
            WireCell::Point curr_pos = fits[i].point;
            WireCell::Point next_pos = fits[i+1].point;
        
            // Create sampling points for Gaussian integration
            std::vector<double> centers_U, centers_V, centers_W, centers_T;
            std::vector<double> sigmas_T, sigmas_U, sigmas_V, sigmas_W;
            std::vector<double> weights;
            
            // Get geometry parameters
            auto test_wpid = m_dv->contained_by(curr_pos);
            if (test_wpid.apa() == -1 || test_wpid.face() == -1) continue;
            int apa = test_wpid.apa();
            int face = test_wpid.face();

            WirePlaneId wpid(kAllLayers, test_wpid.face(), test_wpid.apa());
            auto offset_it = wpid_offsets.find(wpid);
            auto slope_it = wpid_slopes.find(wpid);
            auto geom_it = wpid_geoms.find(wpid);
            
            if (offset_it == wpid_offsets.end() || slope_it == wpid_slopes.end() || geom_it == wpid_geoms.end()) continue;
            auto cluster = segment->cluster();
            const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
            double cluster_t0 = cluster->get_cluster_t0();

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

            // Get anode and tick information for drift time calculation
            auto anode = m_grouping->get_anode(test_wpid.apa());
            auto iface = anode->faces()[test_wpid.face()];
            // double xsign = iface->dirx();
            double xorig = iface->planes()[2]->wires().front()->center().x();  // Anode plane position
            // double tick_size = m_grouping->get_tick().at(test_wpid.apa()).at(test_wpid.face());
            double drift_speed = m_grouping->get_drift_speed().at(test_wpid.apa()).at(test_wpid.face());

            // Sample 5 points each from prev->curr and curr->next
            for (int j = 0; j < 5; j++) {
                // First half: prev -> curr
                WireCell::Point reco_pos = prev_pos + (curr_pos - prev_pos) * (j + 0.5) / 5.0;
                auto reco_pos_raw = transform->backward(reco_pos, cluster_t0, test_wpid.face(), test_wpid.apa());
                
                double central_T = offset_t + slope_x * reco_pos_raw.x();
                double central_U = offset_u + (slope_yu * reco_pos_raw.y() + slope_zu * reco_pos_raw.z());
                double central_V = offset_v + (slope_yv * reco_pos_raw.y() + slope_zv * reco_pos_raw.z());
                double central_W = offset_w + (slope_yw * reco_pos_raw.y() + slope_zw * reco_pos_raw.z());
                
                double weight = (prev_pos - curr_pos).magnitude();
                
                // Calculate drift time from drift distance
                double drift_distance = std::abs(reco_pos.x() - xorig);
                double drift_time = std::max(m_params.min_drift_time, drift_distance / drift_speed);
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
                
                // Second half: curr -> next
                reco_pos = next_pos + (curr_pos - next_pos) * (j + 0.5) / 5.0;
                reco_pos_raw = transform->backward(reco_pos, cluster_t0, test_wpid.face(), test_wpid.apa()); 

                central_T = offset_t + slope_x * reco_pos_raw.x();
                central_U = offset_u + (slope_yu * reco_pos_raw.y() + slope_zu * reco_pos_raw.z());
                central_V = offset_v + (slope_yv * reco_pos_raw.y() + slope_zv * reco_pos_raw.z());
                central_W = offset_w + (slope_yw * reco_pos_raw.y() + slope_zw * reco_pos_raw.z());
                
                weight = (next_pos - curr_pos).magnitude();
                
                // Calculate drift time from drift distance
                drift_distance = std::abs(reco_pos.x() - xorig);
                drift_time = std::max(m_params.min_drift_time, drift_distance / drift_speed);
                
                diff_sigma_L = sqrt(2 * DL * drift_time);
                diff_sigma_T = sqrt(2 * DT * drift_time);
                
            
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
            
            // Fill response matrices using Gaussian integrals - U plane
            int n_u = 0;
            std::set<std::pair<int, int>> set_UT;

            for (const auto& [coord_key, result] : map_U_charge_2D) {
                // const auto& measurement = result.first;
                const auto& coord_2d_set = result.second;
                for (const auto& coord_2d : coord_2d_set) {
                    if (coord_2d.plane != kUlayer || coord_2d.apa != apa || coord_2d.face != face) continue;
                    int wire = coord_2d.wire;
                    int time = coord_2d.time;
                    set_UT.insert(std::make_pair(wire, time));

                    if (abs(wire - centers_U.front()) <= m_params.search_range && abs(time - centers_T.front()) <= m_params.search_range * cur_ntime_ticks) {
                        
                        double value = cal_gaus_integral_seg(coord_2d.time, coord_2d.wire, centers_T, sigmas_T, 
                                                           centers_U, sigmas_U, weights, 0, 4, cur_ntime_ticks);
                        
                        if (result.first.flag == 0 && value > 0) reg_flag_u[idx] = 1;
                        
                        if (value > 0 && result.first.charge > 0 && result.first.flag != 0) {
                            double total_err = sqrt(pow(result.first.charge_err, 2) + 
                                                  pow(result.first.charge * rel_uncer_ind, 2) + 
                                                  pow(add_uncer_ind, 2));
                            RU.insert(n_u, idx) = value / total_err;
                        }
                    }
                }
                n_u++;
            }
            
            // Fill response matrices using Gaussian integrals - V plane
            int n_v = 0;
            std::set<std::pair<int, int>> set_VT;
            for (const auto& [coord_key, result] : map_V_charge_2D) {
                const auto& coord_2d_set = result.second;
                for (const auto& coord_2d : coord_2d_set) {
                    if (coord_2d.plane != kVlayer || coord_2d.apa != apa || coord_2d.face != face) continue;
                    int wire = coord_2d.wire;
                    int time = coord_2d.time;
                    set_VT.insert(std::make_pair(wire, time));
                      if (abs(wire - centers_V.front()) <= m_params.search_range && abs(time - centers_T.front()) <= m_params.search_range * cur_ntime_ticks) {
                        
                        double value = cal_gaus_integral_seg(coord_2d.time, coord_2d.wire, centers_T, sigmas_T, 
                                                           centers_V, sigmas_V, weights, 0, 4, cur_ntime_ticks);
                        
                        if (result.first.flag == 0 && value > 0) reg_flag_v[idx] = 1;
                        
                        if (value > 0 && result.first.charge > 0 && result.first.flag != 0) {
                            double total_err = sqrt(pow(result.first.charge_err, 2) + 
                                                  pow(result.first.charge * rel_uncer_ind, 2) + 
                                                  pow(add_uncer_ind, 2));
                            RV.insert(n_v, idx) = value / total_err;
                        }
                    }
                }
                n_v++;
            }
            
            // Fill response matrices using Gaussian integrals - W plane
            int n_w = 0;
            std::set<std::pair<int, int>> set_WT;
            for (const auto& [coord_key, result] : map_W_charge_2D) {
                const auto& coord_2d_set = result.second;
                for (const auto& coord_2d : coord_2d_set) {
                    if (coord_2d.plane != kWlayer || coord_2d.apa != apa || coord_2d.face != face) continue;
                    int wire = coord_2d.wire;
                    int time = coord_2d.time;
                    set_WT.insert(std::make_pair(wire, time));
                    if (abs(wire - centers_W.front()) <= m_params.search_range && abs(time - centers_T.front()) <= m_params.search_range * cur_ntime_ticks) {

                        double value = cal_gaus_integral_seg(time, wire, centers_T, sigmas_T,
                                                           centers_W, sigmas_W, weights, 0, 4, cur_ntime_ticks);
                        
                        if (result.first.flag == 0 && value > 0) reg_flag_w[idx] = 1;
                        
                        if (value > 0 && result.first.charge > 0 && result.first.flag != 0) {
                            double total_err = sqrt(pow(result.first.charge_err, 2) + 
                                                  pow(result.first.charge * rel_uncer_col, 2) + 
                                                  pow(add_uncer_col, 2));
                            RW.insert(n_w, idx) = value / total_err;
                        }
                    }
                }
                n_w++;
            }
            
            // Additional checks on dead channels for segments
            if (reg_flag_u[idx] == 0) {
                for (size_t kk = 0; kk < centers_U.size(); kk++) {
                    if (set_UT.find(std::make_pair(std::round(centers_U[kk]), std::round(centers_T[kk]/cur_ntime_ticks)*cur_ntime_ticks)) == set_UT.end()) {
                        reg_flag_u[idx] = 1;
                        break;
                    }
                }
            }
            
            if (reg_flag_v[idx] == 0) {
                for (size_t kk = 0; kk < centers_V.size(); kk++) {
                     if (set_VT.find(std::make_pair(std::round(centers_V[kk]), std::round(centers_T[kk]/cur_ntime_ticks)*cur_ntime_ticks)) == set_VT.end()) {
                        reg_flag_v[idx] = 1;
                        break;
                    }
                }
            }
            
            if (reg_flag_w[idx] == 0) {
                for (size_t kk = 0; kk < centers_W.size(); kk++) {
                    if (set_WT.find(std::make_pair(std::round(centers_W[kk]), std::round(centers_T[kk]/cur_ntime_ticks)*cur_ntime_ticks)) == set_WT.end()) {
                        reg_flag_w[idx] = 1;
                        break;
                    }
                }
            }
        }
    }


     // Calculate dx for vertices (endpoints)
    for (const auto& [vertex, vertex_idx] : vertex_index_map) {
        std::vector<WireCell::Point> connected_pts;
        
        WireCell::Clus::Facade::Cluster* cluster = nullptr;
        // Find connected segments
        auto vertex_desc = vertex->get_descriptor();
        if (vertex_desc != PR::Graph::null_vertex()) {
            auto adj_edges = boost::adjacent_vertices(vertex_desc, *m_graph);
            for (auto v_it = adj_edges.first; v_it != adj_edges.second; ++v_it) {
                auto edge_desc = boost::edge(vertex_desc, *v_it, *m_graph);
                if (edge_desc.second) {
                    auto& edge_bundle = (*m_graph)[edge_desc.first];
                    if (edge_bundle.segment && !edge_bundle.segment->fits().empty()) {
                        auto& fits = edge_bundle.segment->fits();
                        cluster = edge_bundle.segment->cluster();

                        if (fits.size() > 1){
                        
                            // std::cout << vertex->fit().point << " " << vertex->fit().index << " " << fits.front().index << " " << fits.back().index << " " << vertex->fit().paf << " " << fits.front().paf << " " << fits.back().paf << " " << (vertex->fit().paf == fits[1].paf) << " " << (vertex->fit().paf == fits[fits.size() - 2].paf) << std::endl;

                            if (vertex->fit().index == fits.front().index) {
                                // Start point
                                if (vertex->fit().paf == fits[1].paf) connected_pts.push_back((fits[1].point + vertex->fit().point) / 2.0);
                            } else {
                                // End point
                                if (vertex->fit().paf == fits[fits.size() - 2].paf) connected_pts.push_back((fits[fits.size() - 2].point + vertex->fit().point) / 2.0);
                            }
                        }
                    }
                }
            }
        }

        // If only one connection, extend endpoint
        if (connected_pts.size() == 1) {
            WireCell::Point curr_pos = vertex->fit().point;
            WireCell::Vector dir = (connected_pts[0] - curr_pos).norm();
            WireCell::Point extended = curr_pos - dir * dis_end_point_ext;
            connected_pts.push_back(extended);

            std::cout << (extended - curr_pos).magnitude()/units::cm << " " << (connected_pts[0] - curr_pos).magnitude()/units::cm << std::endl;
        }
        
        // Calculate total dx
        double total_dx = 0;
        for (const auto& pt : connected_pts) {
            // std::cout << "Vertex connected point: (" << pt.x() << ", " << pt.y() << ", " << pt.z() << ")" << " " << (pt - vertex->fit().point).magnitude()/units::cm << std::endl;
            total_dx += (pt - vertex->fit().point).magnitude();
        }
        local_dx[vertex_idx] = total_dx;

        WireCell::Point curr_pos = vertex->fit().point;

        // Create sampling points for Gaussian integration
        std::vector<double> centers_U, centers_V, centers_W, centers_T;
        std::vector<double> sigmas_T, sigmas_U, sigmas_V, sigmas_W;
        std::vector<double> weights;
        
        // Get geometry parameters
        auto test_wpid = m_dv->contained_by(curr_pos);
        if (test_wpid.apa() == -1 || test_wpid.face() == -1) continue;
        int apa = test_wpid.apa();
        int face = test_wpid.face();

        WirePlaneId wpid(kAllLayers, test_wpid.face(), test_wpid.apa());
        auto offset_it = wpid_offsets.find(wpid);
        auto slope_it = wpid_slopes.find(wpid);
        auto geom_it = wpid_geoms.find(wpid);
        
        if (offset_it == wpid_offsets.end() || slope_it == wpid_slopes.end() || geom_it == wpid_geoms.end()) continue;
        const auto transform = m_pcts->pc_transform(cluster->get_scope_transform(cluster->get_default_scope()));
        double cluster_t0 = cluster->get_cluster_t0();

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

        // Get anode and tick information for drift time calculation
        auto anode = m_grouping->get_anode(test_wpid.apa());
        auto iface = anode->faces()[test_wpid.face()];
        // double xsign = iface->dirx();
        double xorig = iface->planes()[2]->wires().front()->center().x();  // Anode plane position
        // double tick_size = m_grouping->get_tick().at(test_wpid.apa()).at(test_wpid.face());
        double drift_speed = m_grouping->get_drift_speed().at(test_wpid.apa()).at(test_wpid.face());

        auto first_blob = cluster->children()[0];
        int cur_ntime_ticks = first_blob->slice_index_max() - first_blob->slice_index_min();

        for (size_t k = 0; k < connected_pts.size(); k++) {
            for (int j = 0; j < 5; j++) {
                WireCell::Point reco_pos = connected_pts[k] + (curr_pos - connected_pts[k]) * (j + 0.5) / 5.0;
            
                auto reco_pos_raw = transform->backward(reco_pos, cluster_t0, test_wpid.face(), test_wpid.apa());
                
                double central_T = offset_t + slope_x * reco_pos_raw.x();
                double central_U = offset_u + (slope_yu * reco_pos_raw.y() + slope_zu * reco_pos_raw.z());
                double central_V = offset_v + (slope_yv * reco_pos_raw.y() + slope_zv * reco_pos_raw.z());
                double central_W = offset_w + (slope_yw * reco_pos_raw.y() + slope_zw * reco_pos_raw.z());
                
                double weight = (connected_pts[k] - curr_pos).magnitude();
                
                // Calculate drift time from drift distance
                double drift_distance = std::abs(reco_pos.x() - xorig);
                double drift_time = std::max(m_params.min_drift_time, drift_distance / drift_speed);
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
            }
        }

        int n_u = 0;
        std::set<std::pair<int, int>> set_UT;
        for (const auto& [coord_key, result] : map_U_charge_2D) {
            // const auto& measurement = result.first;
            const auto& coord_2d_set = result.second;
            for (const auto& coord_2d : coord_2d_set) {
                if (coord_2d.plane != kUlayer || coord_2d.apa != apa || coord_2d.face != face) continue;
                int wire = coord_2d.wire;
                int time = coord_2d.time;
                set_UT.insert(std::make_pair(wire, time));

                if (abs(wire - centers_U.front()) <= m_params.search_range && abs(time - centers_T.front()) <= m_params.search_range * cur_ntime_ticks) {
                    
                    double value = cal_gaus_integral_seg(coord_2d.time, coord_2d.wire, centers_T, sigmas_T, 
                                                        centers_U, sigmas_U, weights, 0, 4, cur_ntime_ticks);
                    
                    if (result.first.flag == 0 && value > 0) reg_flag_u[vertex_idx] = 1;
                    
                    if (value > 0 && result.first.charge > 0 && result.first.flag != 0) {
                        double total_err = sqrt(pow(result.first.charge_err, 2) + 
                                                pow(result.first.charge * rel_uncer_ind, 2) + 
                                                pow(add_uncer_ind, 2));
                        RU.insert(n_u, vertex_idx) = value / total_err;
                    }
                }
            }
            n_u++;
        }

        // Fill response matrices using Gaussian integrals - V plane
        int n_v = 0;
        std::set<std::pair<int, int>> set_VT;
        for (const auto& [coord_key, result] : map_V_charge_2D) {
            const auto& coord_2d_set = result.second;
            for (const auto& coord_2d : coord_2d_set) {
                if (coord_2d.plane != kVlayer || coord_2d.apa != apa || coord_2d.face != face) continue;
                int wire = coord_2d.wire;
                int time = coord_2d.time;
                set_VT.insert(std::make_pair(wire, time));
                    if (abs(wire - centers_V.front()) <= m_params.search_range && abs(time - centers_T.front()) <= m_params.search_range * cur_ntime_ticks) {
                    
                    double value = cal_gaus_integral_seg(coord_2d.time, coord_2d.wire, centers_T, sigmas_T, 
                                                        centers_V, sigmas_V, weights, 0, 4, cur_ntime_ticks);
                    
                    if (result.first.flag == 0 && value > 0) reg_flag_v[vertex_idx] = 1;
                    
                    if (value > 0 && result.first.charge > 0 && result.first.flag != 0) {
                        double total_err = sqrt(pow(result.first.charge_err, 2) + 
                                                pow(result.first.charge * rel_uncer_ind, 2) + 
                                                pow(add_uncer_ind, 2));
                        RV.insert(n_v, vertex_idx) = value / total_err;
                    }
                }
            }
            n_v++;
        }
        
        // Fill response matrices using Gaussian integrals - W plane
        int n_w = 0;
        std::set<std::pair<int, int>> set_WT;
        for (const auto& [coord_key, result] : map_W_charge_2D) {
            const auto& coord_2d_set = result.second;
            for (const auto& coord_2d : coord_2d_set) {
                if (coord_2d.plane != kWlayer || coord_2d.apa != apa || coord_2d.face != face) continue;
                int wire = coord_2d.wire;
                int time = coord_2d.time;
                set_WT.insert(std::make_pair(wire, time));
                if (abs(wire - centers_W.front()) <= m_params.search_range && abs(time - centers_T.front()) <= m_params.search_range * cur_ntime_ticks) {

                    double value = cal_gaus_integral_seg(time, wire, centers_T, sigmas_T,
                                                        centers_W, sigmas_W, weights, 0, 4, cur_ntime_ticks);
                    
                    if (result.first.flag == 0 && value > 0) reg_flag_w[vertex_idx] = 1;
                    
                    if (value > 0 && result.first.charge > 0 && result.first.flag != 0) {
                        double total_err = sqrt(pow(result.first.charge_err, 2) + 
                                                pow(result.first.charge * rel_uncer_col, 2) + 
                                                pow(add_uncer_col, 2));
                        RW.insert(n_w, vertex_idx) = value / total_err;
                    }
                }
            }
            n_w++;
        }

        // Additional checks on dead channels for segments
        if (reg_flag_u[vertex_idx] == 0) {
            for (size_t kk = 0; kk < centers_U.size(); kk++) {
                if (set_UT.find(std::make_pair(std::round(centers_U[kk]), std::round(centers_T[kk]/cur_ntime_ticks)*cur_ntime_ticks)) == set_UT.end()) {
                    reg_flag_u[vertex_idx] = 1;
                    break;
                }
            }
        }
        
        if (reg_flag_v[vertex_idx] == 0) {
            for (size_t kk = 0; kk < centers_V.size(); kk++) {
                    if (set_VT.find(std::make_pair(std::round(centers_V[kk]), std::round(centers_T[kk]/cur_ntime_ticks)*cur_ntime_ticks)) == set_VT.end()) {
                    reg_flag_v[vertex_idx] = 1;
                    break;
                }
            }
        }
        
        if (reg_flag_w[vertex_idx] == 0) {
            for (size_t kk = 0; kk < centers_W.size(); kk++) {
                if (set_WT.find(std::make_pair(std::round(centers_W[kk]), std::round(centers_T[kk]/cur_ntime_ticks)*cur_ntime_ticks)) == set_WT.end()) {
                    reg_flag_w[vertex_idx] = 1;
                    break;
                }
            }
        }
    }
        
        
    // Build connected_vec for regularization
    std::vector<std::vector<int>> connected_vec(n_3D_pos);
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        auto& edge_bundle = (*m_graph)[*e_it];
        if (!edge_bundle.segment) continue;
        
        auto segment = edge_bundle.segment;
        auto& fits = segment->fits();
        if (fits.empty()) continue;
        
        for (size_t i = 1; i + 1 < fits.size(); i++) {
            int idx = segment_point_index_map[std::make_pair(segment, i)];
            int prev_idx = segment_point_index_map[std::make_pair(segment, i-1)];
            int next_idx = segment_point_index_map[std::make_pair(segment, i+1)];
            
            connected_vec[idx].push_back(prev_idx);
            connected_vec[idx].push_back(next_idx);
        }
    }
    
    // Add vertex connections
    for (const auto& [vertex, vertex_idx] : vertex_index_map) {
        // Find connected segments
        auto vertex_desc = vertex->get_descriptor();
        if (vertex_desc != PR::Graph::null_vertex()) {
            auto adj_edges = boost::adjacent_vertices(vertex_desc, *m_graph);
            for (auto v_it = adj_edges.first; v_it != adj_edges.second; ++v_it) {
                auto edge_desc = boost::edge(vertex_desc, *v_it, *m_graph);
                if (edge_desc.second) {
                    auto& edge_bundle = (*m_graph)[edge_desc.first];
                    if (edge_bundle.segment && !edge_bundle.segment->fits().empty()) {
                        auto& fits = edge_bundle.segment->fits();
                        // Find connected segment points
                        if (vertex_idx == fits.front().index) {
                           connected_vec[vertex_idx].push_back(fits[1].index);
                        } else if (vertex_idx == fits.back().index) {
                           connected_vec[vertex_idx].push_back(fits[fits.size() - 2].index);
                        }
                    }
                }
            }
        }
    }
    
    // Build weight matrices and apply compact matrix analysis
    Eigen::SparseMatrix<double> MU(n_2D_u, n_2D_u), MV(n_2D_v, n_2D_v), MW(n_2D_w, n_2D_w);
    for (int k = 0; k < n_2D_u; k++) MU.insert(k, k) = 1;
    for (int k = 0; k < n_2D_v; k++) MV.insert(k, k) = 1;
    for (int k = 0; k < n_2D_w; k++) MW.insert(k, k) = 1;
    
    Eigen::SparseMatrix<double> RUT = RU.transpose();
    Eigen::SparseMatrix<double> RVT = RV.transpose();
    Eigen::SparseMatrix<double> RWT = RW.transpose();
    
    // Apply compact matrix regularization
    auto overlap_u = calculate_compact_matrix_multi(connected_vec, MU, RUT, n_2D_u, n_3D_pos, 3.0);
    auto overlap_v = calculate_compact_matrix_multi(connected_vec, MV, RVT, n_2D_v, n_3D_pos, 3.0);
    auto overlap_w = calculate_compact_matrix_multi(connected_vec, MW, RWT, n_2D_w, n_3D_pos, 2.0);
    
    // Build regularization matrix
    Eigen::SparseMatrix<double> FMatrix(n_3D_pos, n_3D_pos);
    
    const double dead_ind_weight = m_params.dead_ind_weight;
    const double dead_col_weight = m_params.dead_col_weight;
    const double close_ind_weight = m_params.close_ind_weight;
    const double close_col_weight = m_params.close_col_weight;
    
    for (size_t i = 0; i < n_3D_pos; i++) {
        if (i >= connected_vec.size()) continue;
        
        bool flag_u = reg_flag_u[i];
        bool flag_v = reg_flag_v[i];
        bool flag_w = reg_flag_w[i];
        
        double weight = 0;
        if (flag_u) weight += dead_ind_weight;
        if (flag_v) weight += dead_ind_weight;
        if (flag_w) weight += dead_col_weight;
        
        double scaling = (connected_vec[i].size() > 2) ? 2.0 / connected_vec[i].size() : 1.0;
        
        for (size_t j = 0; j < connected_vec[i].size(); j++) {
            double weight1 = weight;
            int row = i;
            int col = connected_vec[i][j];
            
            if (overlap_u[i][j] > m_params.overlap_th) weight1 += close_ind_weight * pow(overlap_u[i][j] - 0.5, 2);
            if (overlap_v[i][j] > m_params.overlap_th) weight1 += close_ind_weight * pow(overlap_v[i][j] - 0.5, 2);
            if (overlap_w[i][j] > m_params.overlap_th) weight1 += close_col_weight * pow(overlap_w[i][j] - 0.5, 2);
            
            double dx_norm = (local_dx[row] + 0.001 * units::cm) / m_params.dx_norm_length;
            FMatrix.coeffRef(row, row) += -weight1 * scaling / dx_norm;
            FMatrix.coeffRef(row, col) += weight1 * scaling / dx_norm;
        }
    }
    
    // Apply regularization strength
    double lambda = m_params.lambda*8.0/5.0; // adjusted for multi-track fitting ...
    if (!flag_dQ_dx_fit_reg) lambda *= 0.01;
    FMatrix *= lambda;
    
    Eigen::SparseMatrix<double> FMatrixT = FMatrix.transpose();
    
    // Solve the system
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    Eigen::VectorXd b = RUT * MU * data_u_2D + RVT * MV * data_v_2D + RWT * MW * data_w_2D;
    Eigen::SparseMatrix<double> A = RUT * MU * RU + RVT * MV * RV + RWT * MW * RW + FMatrixT * FMatrix;
    
    solver.compute(A);
    pos_3D = solver.solveWithGuess(b, pos_3D_init);
    
    if (std::isnan(solver.error())) {
        pos_3D = solver.solve(b);
    }
    
    // Calculate predictions
    pred_data_u_2D = RU * pos_3D;
    pred_data_v_2D = RV * pos_3D;
    pred_data_w_2D = RW * pos_3D;
    
    // Calculate reduced chi2
    traj_reduced_chi2.clear();
    // std::vector<int> traj_ndf(n_3D_pos, 0);
    
    // Calculate chi-squared contributions from U plane
    for (int k = 0; k < RU.outerSize(); ++k) {
        double sum[3]={0,0,0};
        double sum1[3] = {0,0,0};
        for (Eigen::SparseMatrix<double>::InnerIterator it(RU, k); it; ++it) {
            sum[0] += pow(data_u_2D(it.row()) - pred_data_u_2D(it.row()),2) * (it.value() * pos_3D(k) )/pred_data_u_2D(it.row());
            sum1[0] += (it.value() * pos_3D(k) )/pred_data_u_2D(it.row());      
        }
        for (Eigen::SparseMatrix<double>::InnerIterator it(RV, k); it; ++it) {
            sum[1] += pow(data_v_2D(it.row()) - pred_data_v_2D(it.row()),2) * (it.value() * pos_3D(k))/pred_data_v_2D(it.row());
            sum1[1] += (it.value() * pos_3D(k))/pred_data_v_2D(it.row());
        }
        for (Eigen::SparseMatrix<double>::InnerIterator it(RW, k); it; ++it) {
            sum[2] += pow(data_w_2D(it.row()) - pred_data_w_2D(it.row()),2) * (it.value() * pos_3D(k))/pred_data_w_2D(it.row());
            sum1[2] += (it.value()*pos_3D(k))/pred_data_w_2D(it.row());
        }
        traj_reduced_chi2.push_back(sqrt((sum[0] + sum[1] + sum[2]/4.)/(sum1[0]+sum1[1]+sum1[2])));
    }
    
   
    
    // Update vertex and segment fit results
    for (const auto& [vertex, vertex_idx] : vertex_index_map) {
        if (vertex_idx >= n_3D_pos) continue;
        
        double dQ = pos_3D(vertex_idx);
        double dx = local_dx[vertex_idx];
        double reduced_chi2 = traj_reduced_chi2[vertex_idx];
        
        // Update vertex fit information
        auto& vertex_fit = vertex->fit();
        vertex_fit.dQ = dQ;
        vertex_fit.dx = dx;
        vertex_fit.reduced_chi2 = reduced_chi2;
    }
    
    for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        auto& edge_bundle = (*m_graph)[*e_it];
        if (!edge_bundle.segment) continue;
        
        auto segment = edge_bundle.segment;
        auto& fits = segment->fits();
        if (fits.empty()) continue;
        
        // Update segment fit information
        for (size_t i = 0; i < fits.size(); i++) {
            int idx = segment_point_index_map[std::make_pair(segment, i)];
            if (idx >= n_3D_pos) continue;
            
            fits[i].dQ = pos_3D(idx);
            fits[i].dx = local_dx[idx];
            fits[i].reduced_chi2 = traj_reduced_chi2[idx];
        }
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

        // if (i!=0) std::cout << i << " TTT " << (paf.at(i) == paf.at(i-1)) << std::endl;

        if (i == 0 || (i!=0 && paf.at(i) != paf.at(i-1))) {
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
        } else if (i == n_3D_pos - 1 || (i!=n_3D_pos-1 && paf.at(i) != paf.at(i+1))) {
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
        } else if (paf.at(i) == paf.at(i-1) && paf.at(i) == paf.at(i+1)){
            // Middle point
            prev_rec_pos = (curr_rec_pos + fine_tracking_path.at(i-1).first) * 0.5;
            next_rec_pos = (curr_rec_pos + fine_tracking_path.at(i+1).first) * 0.5;
        }else {
            // Default case (should not happen)
            prev_rec_pos = curr_rec_pos;
            next_rec_pos = curr_rec_pos;
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

        // Get anode and tick information for drift time calculation
        auto anode = m_grouping->get_anode(apa);
        auto iface = anode->faces()[face];
        // double xsign = iface->dirx();
        double xorig = iface->planes()[2]->wires().front()->center().x();  // Anode plane position
        // double tick_size = m_grouping->get_tick().at(apa).at(face);
        double drift_speed = m_grouping->get_drift_speed().at(apa).at(face);

        
        // Calculate previous and next positions for Gaussian integration
        WireCell::Point prev_rec_pos, next_rec_pos;
        if (i == 0 || (i!=0 && paf.at(i) != paf.at(i-1))) {
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
        } else if (i == n_3D_pos - 1 || (i!=n_3D_pos-1 && paf.at(i) != paf.at(i+1))) {
            WireCell::Point prev_point = fine_tracking_path.at(i-1).first;
            prev_rec_pos = (curr_rec_pos + prev_point) * 0.5;
            WireCell::Vector dir = curr_rec_pos - prev_point;
            double length = dir.magnitude();
            if (length > 0) {
                next_rec_pos = curr_rec_pos + (dir / length) * dis_end_point_ext;
            } else {
                next_rec_pos = curr_rec_pos;
            }
        } else if (paf.at(i) == paf.at(i-1) && paf.at(i) == paf.at(i+1)) {
            prev_rec_pos = (curr_rec_pos + fine_tracking_path.at(i-1).first) * 0.5;
            next_rec_pos = (curr_rec_pos + fine_tracking_path.at(i+1).first) * 0.5;
        }else{
            prev_rec_pos = curr_rec_pos;
            next_rec_pos = curr_rec_pos;
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
            
            // Calculate drift time from drift distance
            double drift_distance = std::abs(reco_pos.x() - xorig);
            double drift_time = std::max(m_params.min_drift_time, drift_distance / drift_speed);
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

            // Calculate drift time from drift distance
            drift_distance = std::abs(reco_pos.x() - xorig);
            drift_time = std::max(m_params.min_drift_time, drift_distance / drift_speed);
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

                        // std::cout << n_u << " " << i << " " << time << " " << wire << " " << i << " " << value / total_err << std::endl;
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

                        // std::cout << n_w << " " << i << " " << time << " " << wire << " " << i << " " << value / total_err << std::endl;
                    
                    }
                }
            }
            n_w++;
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
    
    // std::cout << "U: " << std::endl;
    std::vector<std::pair<double, double>> overlap_u = calculate_compact_matrix(MU, RUT, n_2D_u, n_3D_pos, 3);
    // std::cout << "V: " <<std::endl;
    std::vector<std::pair<double, double>> overlap_v = calculate_compact_matrix(MV, RVT, n_2D_v, n_3D_pos, 3);
    // std::cout << "W: " << std::endl;
    std::vector<std::pair<double, double>> overlap_w = calculate_compact_matrix(MW, RWT, n_2D_w, n_3D_pos, 2);
    
//     for (size_t i=0;i!=n_3D_pos;i++){
//         // std::cout << i << " " << reg_flag_u.at(i) << " " << reg_flag_v.at(i) << " " << reg_flag_w.at(i) << std::endl;
//         std::cout << i << " " << (overlap_u.at(i).first + overlap_u.at(i).second)/2. << " " 
//             << (overlap_v.at(i).first + overlap_v.at(i).second)/2. << " " 
//             << (overlap_w.at(i).first + overlap_w.at(i).second)/2. << " " << std::endl;
// //            << MU.coeffRef(i,i) << " " << MV.coeffRef(i,i) << " " << MW.coeffRef(i,i) << std::endl;
//     } 

        // int n_w = 0;
        // for (const auto& [coord_key, result] : map_W_charge_2D) {
        //     const auto& measurement = result.first;
        //     const auto& Coord2D_set = result.second;    
        //     int wire, time;
        //     for (const auto& coord2d : Coord2D_set) {
        //         wire = coord2d.wire;
        //         time = coord2d.time;
        //     }
        //     std::cout << n_w << " " << wire << " " << time << " " << MW.coeffRef(n_w,n_w) << " " << measurement.charge << " " << measurement.charge_err << std::endl;
        //     n_w++;
        // }
    

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
                if (overlap_u[i].first > m_params.overlap_th) weight += close_ind_weight * pow(2 * overlap_u[i].first - 1, 2);
                if (overlap_v[i].first > m_params.overlap_th) weight += close_ind_weight * pow(2 * overlap_v[i].first - 1, 2);
                if (overlap_w[i].first > m_params.overlap_th) weight += close_col_weight * pow(2 * overlap_w[i].first - 1, 2);
            }else{
                if (overlap_u.at(i).first + overlap_u.at(i).second > 2*m_params.overlap_th) weight += close_ind_weight * pow(overlap_u.at(i).first + overlap_u.at(i).second - 1,2);
                if (overlap_v.at(i).first + overlap_v.at(i).second > 2*m_params.overlap_th) weight += close_ind_weight * pow(overlap_v.at(i).first + overlap_v.at(i).second - 1,2);
                if (overlap_w.at(i).first + overlap_w.at(i).second > 2*m_params.overlap_th) weight += close_col_weight * pow(overlap_w.at(i).first + overlap_w.at(i).second - 1,2);

            }
            
            double dx_norm = (dx[i] + 0.001*units::cm) / m_params.dx_norm_length; // Normalize by 0.6 mm
            
            if (i == 0) {
                FMatrix.insert(0, 0) = -weight / dx_norm;
                if (n_3D_pos > 1) FMatrix.insert(0, 1) = weight / ((dx[1] + 0.001*units::cm) / m_params.dx_norm_length) ;
            } else if (i == n_3D_pos - 1) {
                FMatrix.insert(i, i) = -weight / dx_norm; 
                FMatrix.insert(i, i-1) = weight / ((dx[i-1] + 0.001*units::cm) / m_params.dx_norm_length);
            } else {
                FMatrix.insert(i, i) = -2.0 * weight / dx_norm;
                FMatrix.insert(i, i+1) = weight / ((dx[i+1] + 0.001*units::cm) / m_params.dx_norm_length);
                FMatrix.insert(i, i-1) = weight / ((dx[i-1] + 0.001*units::cm) / m_params.dx_norm_length);
            }
            // std::cout << i << " " << flag_u << " " << flag_v << " " << flag_w << " " << overlap_u.at(i).first << " | " <<  overlap_u.at(i).second << " " << overlap_v.at(i).first << " | " <<  overlap_v.at(i).second << " " << overlap_w.at(i).first << " | " << overlap_w.at(i).second << " " << weight << " " << dx_norm << std::endl;
            // std::cout << i << " " << FMatrix.coeff(i, i) << std::endl;
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
    
    // for (int i = 0; i < b.size(); ++i) {
    //     // Example: print or process each element of b
    //     std::cout << "b[" << i << "] = " << b[i] << " " << A.coeff(i, i) << " " << lambda <<  " " << flag_dQ_dx_fit_reg << std::endl;
    //     // You can add your processing logic here
    // }

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

void TrackFitting::do_multi_tracking(bool flag_dQ_dx_fit_reg, bool flag_dQ_dx_fit, bool flag_force_load_data, bool flag_exclusion, bool flag_hack){

    // Reset fit properties for all vertices first
    for (auto vp = boost::vertices(*m_graph); vp.first != vp.second; ++vp.first) {
        auto vd = *vp.first;
        auto& v_bundle = (*m_graph)[vd];
        if (v_bundle.vertex) {
            bool flag_fix = v_bundle.vertex->flag_fix();
            if (!flag_fix){
                WireCell::Point p = v_bundle.vertex->wcpt().point;
                auto& vertex_fit = v_bundle.vertex->fit();
                vertex_fit.point = p;
            }

            // std::cout << "Vertex fit point before reset: " << flag_fix << " " << v_bundle.vertex->fit().point << " " << v_bundle.vertex->wcpt().point << std::endl;
            // v_bundle.vertex->reset_fit_prop();
            // v_bundle.vertex->flag_fix(flag_fix);
        }
    }

    bool flag_1st_tracking = true;
    bool flag_2nd_tracking = true;
    bool flag_dQ_dx = flag_dQ_dx_fit;
    
    // Prepare the data for the fit - collect charge information from 2D projections
    if (flag_force_load_data || global_rb_map.size() == 0){
        prepare_data();
        fill_global_rb_map();
    }

       

    // auto edge_range = boost::edges(*m_graph);
    //  for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
    //         auto& edge_bundle = (*m_graph)[*e_it];
    //         if (edge_bundle.segment) {
    //             std::cout << "Segment fits size: " << edge_bundle.segment->fits().size() << std::endl;
    //             for (const auto& fit : edge_bundle.segment->fits()) {
    //                 std::cout << "  Fit point: " << fit.point << " " << fit.index << std::endl;
    //             }
    //         }
    //     }
    //     for (auto vp = boost::vertices(*m_graph); vp.first != vp.second; ++vp.first) {
    //         auto vd = *vp.first;
    //         auto& v_bundle = (*m_graph)[vd];
    //         if (v_bundle.vertex) {
    //             std::cout << v_bundle.vertex->fit().point << " VTX " << v_bundle.vertex->fit().index << std::endl;
    //         }
    //     }

    // for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
    //     auto& edge_bundle = (*m_graph)[*e_it];
    //     if (edge_bundle.segment) {
    //         std::cout << "Segment fits size: " << edge_bundle.segment->fits().size() << std::endl;
    //     }
    // }

    // First round of organizing the path from the path_wcps (shortest path)
    double low_dis_limit = m_params.low_dis_limit;
    double end_point_limit = m_params.end_point_limit;
    organize_segments_path(low_dis_limit, end_point_limit);

    // std::cout << "After first organization " << std::endl;

    // for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
    //     auto& edge_bundle = (*m_graph)[*e_it];
    //     if (edge_bundle.segment) {
    //         std::cout << "Segment fits size: " << edge_bundle.segment->fits().size() << std::endl;
    //     }
    // }
    //   for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
    //         auto& edge_bundle = (*m_graph)[*e_it];
    //         if (edge_bundle.segment) {
    //             std::cout << "Segment fits size: " << edge_bundle.segment->fits().size() << std::endl;
    //             for (const auto& fit : edge_bundle.segment->fits()) {
    //                 std::cout << "  Fit point: " << fit.point << " " << fit.index << std::endl;
    //             }
    //         }
    //     }
    //     for (auto vp = boost::vertices(*m_graph); vp.first != vp.second; ++vp.first) {
    //         auto vd = *vp.first;
    //         auto& v_bundle = (*m_graph)[vd];
    //         if (v_bundle.vertex) {
    //             std::cout << v_bundle.vertex->fit().point << " VTX " << v_bundle.vertex->fit().index << std::endl;
    //         }
    //     }
  
    if (flag_1st_tracking){

        // auto edge_range = boost::edges(*m_graph);
        // for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        //     auto& edge_bundle = (*m_graph)[*e_it];
        //     if (edge_bundle.segment) {
        //         std::cout << "Segment fits size: " << edge_bundle.segment->fits().size() << std::endl;
        //         for (const auto& fit : edge_bundle.segment->fits()) {
        //             std::cout << "  Fit point: " << fit.point << " " << fit.index << std::endl;
        //         }
        //     }
        // }
        // for (auto vp = boost::vertices(*m_graph); vp.first != vp.second; ++vp.first) {
        //     auto vd = *vp.first;
        //     auto& v_bundle = (*m_graph)[vd];
        //     if (v_bundle.vertex) {
        //         std::cout << v_bundle.vertex->fit().point << " VTX " << v_bundle.vertex->fit().index << std::endl;
        //     }
        // }
        form_map_graph(flag_exclusion, m_params.end_point_factor, m_params.mid_point_factor, m_params.nlevel, m_params.time_tick_cut, m_params.charge_cut);
        //  for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        //     auto& edge_bundle = (*m_graph)[*e_it];
        //     if (edge_bundle.segment) {
        //         std::cout << "Segment fits size: " << edge_bundle.segment->fits().size() << std::endl;
        //         for (const auto& fit : edge_bundle.segment->fits()) {
        //             std::cout << "  Fit point: " << fit.point << " " << fit.index << std::endl;
        //         }
        //     }
        // }
        // for (auto vp = boost::vertices(*m_graph); vp.first != vp.second; ++vp.first) {
        //     auto vd = *vp.first;
        //     auto& v_bundle = (*m_graph)[vd];
        //     if (v_bundle.vertex) {
        //         std::cout << v_bundle.vertex->fit().point << " VTX " << v_bundle.vertex->fit().index << std::endl;
        //     }
        // }
        multi_trajectory_fit(1, m_params.div_sigma);
    }

  
    if (flag_2nd_tracking){
        // second round trajectory fit ...
        low_dis_limit = m_params.low_dis_limit/2.;
        end_point_limit = m_params.end_point_limit/2.;
        
        // auto edge_range = boost::edges(*m_graph);
        // for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        //     auto& edge_bundle = (*m_graph)[*e_it];
        //     if (edge_bundle.segment) {
        //         std::cout << "Segment fits size: " << edge_bundle.segment->fits().size() << std::endl;
        //         for (const auto& fit : edge_bundle.segment->fits()) {
        //             std::cout << "  Fit point: " << fit.point << " " << fit.index << std::endl;
        //         }
        //     }
        // }
        // for (auto vp = boost::vertices(*m_graph); vp.first != vp.second; ++vp.first) {
        //     auto vd = *vp.first;
        //     auto& v_bundle = (*m_graph)[vd];
        //     if (v_bundle.vertex) {
        //         std::cout << v_bundle.vertex->fit().point << " VTX " << v_bundle.vertex->fit().index << std::endl;
        //     }
        // }

        // organize path
        organize_segments_path_2nd(low_dis_limit, end_point_limit);    
        
        // std::cout << "After second organization " << std::endl;
        //  for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        //     auto& edge_bundle = (*m_graph)[*e_it];
        //     if (edge_bundle.segment) {
        //         std::cout << "Segment fits size: " << edge_bundle.segment->fits().size() << std::endl;
        //         for (const auto& fit : edge_bundle.segment->fits()) {
        //             std::cout << "  Fit point: " << fit.point << " " << fit.index << std::endl;
        //         }
        //     }
        // }
        // for (auto vp = boost::vertices(*m_graph); vp.first != vp.second; ++vp.first) {
        //     auto vd = *vp.first;
        //     auto& v_bundle = (*m_graph)[vd];
        //     if (v_bundle.vertex) {
        //         std::cout << v_bundle.vertex->fit().point << " VTX " << v_bundle.vertex->fit().index << std::endl;
        //     }
        // }

        form_map_graph(flag_exclusion, m_params.end_point_factor, m_params.mid_point_factor, m_params.nlevel, m_params.time_tick_cut, m_params.charge_cut);
        multi_trajectory_fit(1, m_params.div_sigma);

        //  for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        //     auto& edge_bundle = (*m_graph)[*e_it];
        //     if (edge_bundle.segment) {
        //         std::cout << "Segment fits size: " << edge_bundle.segment->fits().size() << std::endl;
        //         for (const auto& fit : edge_bundle.segment->fits()) {
        //             std::cout << "  Fit point: " << fit.point << " " << fit.index << " " << fit.paf.first << " " << fit.paf.second << std::endl;
        //         }
        //     }
        // }
        // for (auto vp = boost::vertices(*m_graph); vp.first != vp.second; ++vp.first) {
        //     auto vd = *vp.first;
        //     auto& v_bundle = (*m_graph)[vd];
        //     if (v_bundle.vertex) {
        //         std::cout << v_bundle.vertex->fit().point << " VTX " << v_bundle.vertex->fit().index << " " << v_bundle.vertex->fit().paf.first << " " << v_bundle.vertex->fit().paf.second << std::endl;
        //     }
        // }

        // organize path
        low_dis_limit = 0.6*units::cm;
        organize_segments_path_3rd(low_dis_limit);

        //  for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
        //     auto& edge_bundle = (*m_graph)[*e_it];
        //     if (edge_bundle.segment) {
        //         std::cout << "Segment fits size: " << edge_bundle.segment->fits().size() << std::endl;
        //         for (const auto& fit : edge_bundle.segment->fits()) {
        //             std::cout << "  Fit point: " << fit.point << " " << fit.index << " " << fit.paf.first << " " << fit.paf.second << std::endl;
        //         }
        //     }
        // }
        // for (auto vp = boost::vertices(*m_graph); vp.first != vp.second; ++vp.first) {
        //     auto vd = *vp.first;
        //     auto& v_bundle = (*m_graph)[vd];
        //     if (v_bundle.vertex) {
        //         std::cout << v_bundle.vertex->fit().point << " VTX " << v_bundle.vertex->fit().index << " " << v_bundle.vertex->fit().paf.first << " " << v_bundle.vertex->fit().paf.second << std::endl;
        //     }
        // }

    }
  
  
    if (flag_dQ_dx){
        for (auto vp = boost::vertices(*m_graph); vp.first != vp.second; ++vp.first) {
            auto vd = *vp.first;
            auto& v_bundle = (*m_graph)[vd];
            if (v_bundle.vertex) {
                bool flag_fix = v_bundle.vertex->flag_fix();
                v_bundle.vertex->reset_fit_prop();
                v_bundle.vertex->flag_fix(flag_fix);
            }
        }

        auto edge_range = boost::edges(*m_graph);    
        for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
            auto& edge_bundle = (*m_graph)[*e_it];
            if (edge_bundle.segment) {
                edge_bundle.segment->reset_fit_prop();
            }
        }
        form_map_graph(flag_exclusion, m_params.end_point_factor, m_params.mid_point_factor, m_params.nlevel, m_params.time_tick_cut, m_params.charge_cut);

    
        dQ_dx_multi_fit(end_point_limit, flag_dQ_dx_fit_reg);

        for (auto e_it = edge_range.first; e_it != edge_range.second; ++e_it) {
            auto& edge_bundle = (*m_graph)[*e_it];
            if (edge_bundle.segment) {
                std::cout << "Segment fits size: " << edge_bundle.segment->fits().size() << std::endl;
                for (const auto& fit : edge_bundle.segment->fits()) {
                    std::cout << "  Fit point: " << fit.point << " " << fit.index << " " << fit.dQ << " " << fit.dx/units::cm << " " << fit.reduced_chi2 << std::endl;
                }
            }
        }
        for (auto vp = boost::vertices(*m_graph); vp.first != vp.second; ++vp.first) {
            auto vd = *vp.first;
            auto& v_bundle = (*m_graph)[vd];
            if (v_bundle.vertex) {
                std::cout << v_bundle.vertex->fit().point << " VTX " << v_bundle.vertex->fit().index << " " << v_bundle.vertex->fit().dQ << " " << v_bundle.vertex->fit().dx/units::cm << " " << v_bundle.vertex->fit().reduced_chi2 << std::endl;
            }
        }

    }
}


void TrackFitting::do_single_tracking(std::shared_ptr<PR::Segment> segment, bool flag_dQ_dx_fit_reg, bool flag_dQ_dx_fit, bool flag_force_load_data, bool flag_hack) {
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

        // if (flag_hack){
        //     // hack pts ...
        //     pts.clear();
        //     pts.push_back(WireCell::Point(219.209*units::cm, -87.2848*units::cm, 209.453*units::cm));
        //     pts.push_back(WireCell::Point(219.011*units::cm, -87.8189*units::cm, 209.55*units::cm));
        //     pts.push_back(WireCell::Point(218.613*units::cm, -88.663*units::cm, 209.722*units::cm));
        //     pts.push_back(WireCell::Point(218.329*units::cm, -89.169*units::cm, 209.853*units::cm));
        //     pts.push_back(WireCell::Point(218.09*units::cm, -89.8885*units::cm, 209.969*units::cm));
        //     pts.push_back(WireCell::Point(217.858*units::cm, -90.3128*units::cm, 210.076*units::cm));
        //     pts.push_back(WireCell::Point(217.627*units::cm, -90.7371*units::cm, 210.184*units::cm));
        //     pts.push_back(WireCell::Point(217.423*units::cm, -91.4211*units::cm, 210.266*units::cm));
        //     pts.push_back(WireCell::Point(217.111*units::cm, -91.8551*units::cm, 210.423*units::cm));
        //     pts.push_back(WireCell::Point(216.84*units::cm, -92.3369*units::cm, 210.55*units::cm));
        //     pts.push_back(WireCell::Point(216.248*units::cm, -92.7898*units::cm, 210.798*units::cm));
            
        //     // pts.clear();
        //     // pts.push_back(WireCell::Point(216.165*units::cm, -93.8529*units::cm, 210.894*units::cm));
        //     // pts.push_back(WireCell::Point(215.996*units::cm, -94.3005*units::cm, 210.942*units::cm));
        //     // pts.push_back(WireCell::Point(215.888*units::cm, -94.8903*units::cm, 211.002*units::cm));
        //     // pts.push_back(WireCell::Point(215.38*units::cm, -95.4503*units::cm, 211.276*units::cm));
        //     // pts.push_back(WireCell::Point(215.128*units::cm, -96.2677*units::cm, 211.414*units::cm));
        //     // pts.push_back(WireCell::Point(214.801*units::cm, -96.7483*units::cm, 211.583*units::cm));
        //     // pts.push_back(WireCell::Point(214.429*units::cm, -97.4961*units::cm, 211.745*units::cm));
        //     // pts.push_back(WireCell::Point(214.225*units::cm, -98.0668*units::cm, 211.862*units::cm));
        //     // pts.push_back(WireCell::Point(213.766*units::cm, -98.833*units::cm, 212.054*units::cm));
        //     // pts.push_back(WireCell::Point(213.401*units::cm, -99.5325*units::cm, 212.199*units::cm));
        //     // pts.push_back(WireCell::Point(213.182*units::cm, -100.042*units::cm, 212.303*units::cm));
        //     // pts.push_back(WireCell::Point(212.763*units::cm, -100.63*units::cm, 212.518*units::cm));
        //     // pts.push_back(WireCell::Point(212.414*units::cm, -101.107*units::cm, 212.697*units::cm));
        //     // pts.push_back(WireCell::Point(212.25*units::cm, -101.563*units::cm, 212.785*units::cm));
        //     // pts.push_back(WireCell::Point(212.086*units::cm, -102.018*units::cm, 212.872*units::cm));
        //     // pts.push_back(WireCell::Point(212.224*units::cm, -102.854*units::cm, 212.839*units::cm));
        //     //

        //     // pts.clear();
        //     // pts.push_back(WireCell::Point(219.209 * units::cm, -87.2848 * units::cm, 209.453 * units::cm));
        //     // pts.push_back(WireCell::Point(218.997 * units::cm, -87.8182 * units::cm, 209.557 * units::cm));
        //     // pts.push_back(WireCell::Point(218.806 * units::cm, -88.253 * units::cm, 209.641 * units::cm));
        //     // pts.push_back(WireCell::Point(218.615 * units::cm, -88.6879 * units::cm, 209.725 * units::cm));
        //     // pts.push_back(WireCell::Point(218.337 * units::cm, -89.1079 * units::cm, 209.836 * units::cm));
        //     // pts.push_back(WireCell::Point(218.051 * units::cm, -89.8565 * units::cm, 209.978 * units::cm));
        //     // pts.push_back(WireCell::Point(217.77 * units::cm, -90.5317 * units::cm, 210.12 * units::cm));
        //     // pts.push_back(WireCell::Point(217.472 * units::cm, -91.2331 * units::cm, 210.25 * units::cm));
        //     // pts.push_back(WireCell::Point(217.058 * units::cm, -91.8449 * units::cm, 210.441 * units::cm));
        //     // pts.push_back(WireCell::Point(216.822 * units::cm, -92.5268 * units::cm, 210.548 * units::cm));
        //     // pts.push_back(WireCell::Point(216.61 * units::cm, -92.9509 * units::cm, 210.65 * units::cm));
        //     // pts.push_back(WireCell::Point(216.308 * units::cm, -93.6788 * units::cm, 210.81 * units::cm));
        //     // pts.push_back(WireCell::Point(215.992 * units::cm, -94.2301 * units::cm, 210.939 * units::cm));
        //     // pts.push_back(WireCell::Point(215.791 * units::cm, -94.7826 * units::cm, 211.062 * units::cm));
        //     // pts.push_back(WireCell::Point(215.532 * units::cm, -95.1674 * units::cm, 211.193 * units::cm));
        //     // pts.push_back(WireCell::Point(215.274 * units::cm, -95.8492 * units::cm, 211.346 * units::cm));
        //     // pts.push_back(WireCell::Point(215.038 * units::cm, -96.3231 * units::cm, 211.467 * units::cm));
        //     // pts.push_back(WireCell::Point(214.741 * units::cm, -96.7608 * units::cm, 211.604 * units::cm));
        //     // pts.push_back(WireCell::Point(214.444 * units::cm, -97.1985 * units::cm, 211.741 * units::cm));
        //     // pts.push_back(WireCell::Point(214.282 * units::cm, -97.8976 * units::cm, 211.824 * units::cm));
        //     // pts.push_back(WireCell::Point(214.026 * units::cm, -98.3316 * units::cm, 211.943 * units::cm));
        //     // pts.push_back(WireCell::Point(213.682 * units::cm, -99.0357 * units::cm, 212.085 * units::cm));
        //     // pts.push_back(WireCell::Point(213.368 * units::cm, -99.7199 * units::cm, 212.216 * units::cm));
        //     // pts.push_back(WireCell::Point(213.101 * units::cm, -100.164 * units::cm, 212.341 * units::cm));
        //     // pts.push_back(WireCell::Point(212.742 * units::cm, -100.639 * units::cm, 212.527 * units::cm));
        //     // pts.push_back(WireCell::Point(212.514 * units::cm, -101.179 * units::cm, 212.647 * units::cm));
        //     // pts.push_back(WireCell::Point(212.117 * units::cm, -101.817 * units::cm, 212.831 * units::cm));
        //     // pts.push_back(WireCell::Point(211.977 * units::cm, -102.455 * units::cm, 212.891 * units::cm));
        // }    

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
    std::vector<geo_point_t> path_points;
    path_points.reserve(fine_tracking_path.size());
    for (size_t i = 0; i < npoints; ++i) {
        PR::Fit fit;
        
        // Set the fitted 3D point
        fit.point = fine_tracking_path[i].first;
        path_points.push_back(fit.point);

        // Set physics quantities
        fit.dQ =  dQ[i];
        fit.dx = dx[i];
        fit.pu = pu[i];
        fit.pv = pv[i];
        fit.pw = pw[i];
        fit.pt = pt[i];
        fit.paf = paf[i];
        // std::cout <<"test " << fit.paf.first << " " << fit.paf.second << " " << paf[i].first << " " << paf[i].second << std::endl;
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

    // replace point cloud after track fitting ...
    PR::create_segment_point_cloud(segment, path_points, m_dv, "main");
    PR::create_segment_fit_point_cloud(segment, m_dv, "fit");

}
