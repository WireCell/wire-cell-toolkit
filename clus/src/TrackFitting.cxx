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

void TrackFitting::add_segment(PR::Segment* segment){
    m_segments.insert(segment);
    m_clusters.insert(segment->cluster());
    m_grouping = segment->cluster()->grouping();

    // std::cout << "TrackFitting: Added segment with " << segment->wcpts().size() << " points." << std::endl;
}

geo_point_t TrackFitting::adjust_rough_path(PR::Segment& segment){
    return geo_point_t(0,0,0);
}

void TrackFitting::prepare_data(){

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