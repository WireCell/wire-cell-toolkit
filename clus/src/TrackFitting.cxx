#include "WireCellClus/TrackFitting.h"
#include "WireCellUtil/Logging.h"

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
    m_grouping = segment->cluster()->grouping();

    for (auto& cluster: m_clusters){
        for (auto& blob: cluster->children()){
            m_blobs.insert(blob);
        }
    }

    std::cout << "TrackFitting: Added segment with " << segment->wcpts().size() << " points." << " " << m_clusters.size() << " " << m_blobs.size() << std::endl;
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