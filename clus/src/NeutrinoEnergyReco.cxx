#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRShowerFunctions.h"

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;


double PatternAlgorithms::cal_corr_factor(WireCell::Point& pt, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    double corr_factor = 1.0;
    // So far this is an empty class that needs to be filled with actual logic ...

    // Example 1: Find APA and face using detector volumes
    // The WirePlaneId contains apa and face information
    WirePlaneId wpid = dv->contained_by(pt);
    int apa = wpid.apa();
    int face = wpid.face();
    int plane = wpid.index();  // 0=U, 1=V, 2=W
    
    // std::cout << "Point at x=" << pt.x()/units::cm << " y=" << pt.y()/units::cm 
    //           << " z=" << pt.z()/units::cm << " cm" << std::endl;
    // std::cout << "  APA=" << apa << " Face=" << face << " Plane=" << plane << std::endl;
    
    // Example 2: Find the grouping from track_fitter
    // The track_fitter contains a reference to the grouping
    auto grouping = track_fitter.grouping();
    
    // if (grouping) {
    //     std::cout << "  Found grouping with " << grouping->children().size() << " clusters" << std::endl;
        
    //     // Example: Access detector parameters from grouping cache
    //     double tick = grouping->get_tick().at(apa).at(face);
    //     double time_offset = grouping->get_time_offset().at(apa).at(face);
    //     double drift_speed = grouping->get_drift_speed().at(apa).at(face);
        
    //     std::cout << "  Tick=" << tick/units::us << " us, "
    //               << "TimeOffset=" << time_offset/units::us << " us, "
    //               << "DriftSpeed=" << drift_speed/(units::mm/units::us) << " mm/us" << std::endl;
        
    //     // Example: Access channel information for a wire index
    //     int wire_index = 100;  // example wire index
    //     WirePlaneLayer_t layer = static_cast<WirePlaneLayer_t>(plane);
        
    //     try {
    //         auto channel = grouping->get_plane_channel_wind(apa, face, layer, wire_index);
    //         int channel_ident = channel->ident();
    //         std::cout << "  Wire " << wire_index << " -> Channel " << channel_ident << std::endl;
    //     } catch (...) {
    //         std::cout << "  Wire " << wire_index << " not found in this plane" << std::endl;
    //     }
        
    //     // Example: Access charge data in a region around the point
    //     // Convert 3D point to time slice and wire index
    //     auto [time_slice, wire_idx] = grouping->convert_3Dpoint_time_ch(pt, apa, face, plane);
    //     std::cout << "  3D point converts to: time_slice=" << time_slice 
    //               << ", wire_index=" << wire_idx << std::endl;
        
    //     // Get charge data in a small window around this point
    //     int time_window = 5;  // +/- 5 time slices
    //     int wire_window = 3;  // +/- 3 wires
        
    //     auto charge_map = grouping->get_overlap_good_ch_charge(
    //         time_slice - time_window, time_slice + time_window,
    //         wire_idx - wire_window, wire_idx + wire_window,
    //         apa, face, plane
    //     );
        
    //     std::cout << "  Found " << charge_map.size() << " charge measurements nearby" << std::endl;
        
    //     // Example: Iterate through charge measurements
    //     double total_charge = 0;
    //     for (const auto& [key, value] : charge_map) {
    //         int t_slice = key.first;
    //         int w_idx = key.second;
    //         double charge = value.first;
    //         double uncertainty = value.second;
    //         total_charge += charge;
            
    //         // Optionally print details
    //         // std::cout << "    t=" << t_slice << " w=" << w_idx 
    //         //           << " Q=" << charge << " +/- " << uncertainty << std::endl;
    //     }
    //     std::cout << "  Total charge in window: " << total_charge << std::endl;
    // } else {
    //     std::cout << "  Warning: No grouping found in track_fitter" << std::endl;
    // }

    (void)apa;
    (void)face;
    (void)plane;
    (void)grouping;
    
    return corr_factor;
}
