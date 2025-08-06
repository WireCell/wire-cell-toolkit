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

void TrackFitting::collect_charge(double dis_cut, double range_cut){
    // collect charge within dis_cut and range_cut
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