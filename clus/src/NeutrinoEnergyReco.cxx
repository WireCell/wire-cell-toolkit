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


double PatternAlgorithms::cal_kine_charge(ShowerPtr shower, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    (void)graph; // Unused for now
    
    if (!shower) return 0.0;
    
    auto grouping = track_fitter.grouping();
    if (!grouping) return 0.0;
    
    double kine_energy = 0;
    
    // Recombination and fudge factors based on particle type
    double fudge_factor = 0.95;
    double recom_factor = 0.7;
    
    if (shower->get_flag_shower()) {
        recom_factor = 0.5; // assume shower
        fudge_factor = 0.8; // shower ...
    } else if (std::abs(shower->get_particle_type()) == 2212) {
        recom_factor = 0.35; // proton
    }
    
    // Collect 2D charge data divided by plane
    std::map<TrackFitting::CoordReadout, TrackFitting::ChargeMeasurement> charge_2d_u;
    std::map<TrackFitting::CoordReadout, TrackFitting::ChargeMeasurement> charge_2d_v;
    std::map<TrackFitting::CoordReadout, TrackFitting::ChargeMeasurement> charge_2d_w;
    std::map<std::pair<int, int>, std::vector<std::tuple<int, int, int>>> map_apa_ch_plane_wires;
    
    track_fitter.collect_2D_charge(charge_2d_u, charge_2d_v, charge_2d_w, map_apa_ch_plane_wires);
    
    // Get point clouds from shower
    auto pcloud1 = shower->get_pcloud("associate_points"); // associated points
    auto pcloud2 = shower->get_pcloud("fit");              // fit points
    
    if (!pcloud1 && !pcloud2) return 0;
    if (!pcloud1 && pcloud2) pcloud1 = pcloud2;
    if (!pcloud2 && pcloud1) pcloud2 = pcloud1;
    
    double sum_u_charge = 0;
    double sum_v_charge = 0;
    double sum_w_charge = 0;
    
    double dis_cut = 0.6 * units::cm;
    
    // Process U plane charges
    for (const auto& [coord_key, charge_data] : charge_2d_u) {
        int time_slice = coord_key.time;
        int channel = coord_key.channel;
        int apa = coord_key.apa;
        
        // Get wire info for this channel
        auto apa_ch_key = std::make_pair(apa, channel);
        auto wire_it = map_apa_ch_plane_wires.find(apa_ch_key);
        if (wire_it == map_apa_ch_plane_wires.end()) continue;
        
        int face = -1;
        for (const auto& [f, plane, wire] : wire_it->second) {
            if (plane == 0) { // U plane
                face = f;
                break;
            }
        }
        if (face < 0) continue;
        
        // Convert time and channel to 2D point
        auto p2d = grouping->convert_time_ch_2Dpoint(time_slice, channel, apa, face, 0);
        WireCell::Point test_p2d(p2d.first, p2d.second, 0);
        
        // Find closest 3D point in point clouds
        double dis = 1e9;
        size_t point_index = 0;
        const Facade::Cluster* closest_cluster = nullptr;
        
        if (pcloud1) {
            auto result = pcloud1->get_closest_2d_point_info(test_p2d, 0, face, apa);
            dis = std::get<0>(result);
            closest_cluster = std::get<1>(result);
            point_index = std::get<2>(result);
        }
        
        if (dis < dis_cut && closest_cluster) {
            const auto& points = pcloud1->get_points();
            if (point_index < points.size()) {
                WireCell::Point test_p(points[point_index].x, points[point_index].y, points[point_index].z);
                double factor = cal_corr_factor(test_p, track_fitter, dv);
                sum_u_charge += charge_data.charge * factor;
            }
        } else if (pcloud2) {
            // Try second point cloud
            auto result = pcloud2->get_closest_2d_point_info(test_p2d, 0, face, apa);
            dis = std::get<0>(result);
            closest_cluster = std::get<1>(result);
            point_index = std::get<2>(result);
            
            if (dis < dis_cut && closest_cluster) {
                const auto& points = pcloud2->get_points();
                if (point_index < points.size()) {
                    WireCell::Point test_p(points[point_index].x, points[point_index].y, points[point_index].z);
                    double factor = cal_corr_factor(test_p, track_fitter, dv);
                    sum_u_charge += charge_data.charge * factor;
                }
            }
        }
    }
    
    // Process V plane charges
    for (const auto& [coord_key, charge_data] : charge_2d_v) {
        int time_slice = coord_key.time;
        int channel = coord_key.channel;
        int apa = coord_key.apa;
        
        // Get wire info for this channel
        auto apa_ch_key = std::make_pair(apa, channel);
        auto wire_it = map_apa_ch_plane_wires.find(apa_ch_key);
        if (wire_it == map_apa_ch_plane_wires.end()) continue;
        
        int face = -1;
        for (const auto& [f, plane, wire] : wire_it->second) {
            if (plane == 1) { // V plane
                face = f;
                break;
            }
        }
        if (face < 0) continue;
        
        // Convert time and channel to 2D point
        auto p2d = grouping->convert_time_ch_2Dpoint(time_slice, channel, apa, face, 1);
        WireCell::Point test_p2d(p2d.first, p2d.second, 0);
        
        // Find closest 3D point in point clouds
        double dis = 1e9;
        size_t point_index = 0;
        const Facade::Cluster* closest_cluster = nullptr;
        
        if (pcloud1) {
            auto result = pcloud1->get_closest_2d_point_info(test_p2d, 1, face, apa);
            dis = std::get<0>(result);
            closest_cluster = std::get<1>(result);
            point_index = std::get<2>(result);
        }
        
        if (dis < dis_cut && closest_cluster) {
            const auto& points = pcloud1->get_points();
            if (point_index < points.size()) {
                WireCell::Point test_p(points[point_index].x, points[point_index].y, points[point_index].z);
                double factor = cal_corr_factor(test_p, track_fitter, dv);
                sum_v_charge += charge_data.charge * factor;
            }
        } else if (pcloud2) {
            // Try second point cloud
            auto result = pcloud2->get_closest_2d_point_info(test_p2d, 1, face, apa);
            dis = std::get<0>(result);
            closest_cluster = std::get<1>(result);
            point_index = std::get<2>(result);
            
            if (dis < dis_cut && closest_cluster) {
                const auto& points = pcloud2->get_points();
                if (point_index < points.size()) {
                    WireCell::Point test_p(points[point_index].x, points[point_index].y, points[point_index].z);
                    double factor = cal_corr_factor(test_p, track_fitter, dv);
                    sum_v_charge += charge_data.charge * factor;
                }
            }
        }
    }
    
    // Process W plane charges
    for (const auto& [coord_key, charge_data] : charge_2d_w) {
        int time_slice = coord_key.time;
        int channel = coord_key.channel;
        int apa = coord_key.apa;
        
        // Get wire info for this channel
        auto apa_ch_key = std::make_pair(apa, channel);
        auto wire_it = map_apa_ch_plane_wires.find(apa_ch_key);
        if (wire_it == map_apa_ch_plane_wires.end()) continue;
        
        int face = -1;
        for (const auto& [f, plane, wire] : wire_it->second) {
            if (plane == 2) { // W plane
                face = f;
                break;
            }
        }
        if (face < 0) continue;
        
        // Convert time and channel to 2D point
        auto p2d = grouping->convert_time_ch_2Dpoint(time_slice, channel, apa, face, 2);
        WireCell::Point test_p2d(p2d.first, p2d.second, 0);
        
        // Find closest 3D point in point clouds
        double dis = 1e9;
        size_t point_index = 0;
        const Facade::Cluster* closest_cluster = nullptr;
        
        if (pcloud1) {
            auto result = pcloud1->get_closest_2d_point_info(test_p2d, 2, face, apa);
            dis = std::get<0>(result);
            closest_cluster = std::get<1>(result);
            point_index = std::get<2>(result);
        }
        
        if (dis < dis_cut && closest_cluster) {
            const auto& points = pcloud1->get_points();
            if (point_index < points.size()) {
                WireCell::Point test_p(points[point_index].x, points[point_index].y, points[point_index].z);
                double factor = cal_corr_factor(test_p, track_fitter, dv);
                sum_w_charge += charge_data.charge * factor;
            }
        } else if (pcloud2) {
            // Try second point cloud
            auto result = pcloud2->get_closest_2d_point_info(test_p2d, 2, face, apa);
            dis = std::get<0>(result);
            closest_cluster = std::get<1>(result);
            point_index = std::get<2>(result);
            
            if (dis < dis_cut && closest_cluster) {
                const auto& points = pcloud2->get_points();
                if (point_index < points.size()) {
                    WireCell::Point test_p(points[point_index].x, points[point_index].y, points[point_index].z);
                    double factor = cal_corr_factor(test_p, track_fitter, dv);
                    sum_w_charge += charge_data.charge * factor;
                }
            }
        }
    }
    
    // Calculate overall charge using weighted average
    double charge[3] = {sum_u_charge, sum_v_charge, sum_w_charge};
    double weight[3] = {0.25, 0.25, 1.0};
    
    // Find min, max, and median charges
    int min_index = 0, max_index = 0, med_index = 0;
    double min_charge = 1e9, max_charge = -1e9;
    for (int i = 0; i < 3; i++) {
        if (min_charge > charge[i]) {
            min_charge = charge[i];
            min_index = i;
        }
        if (max_charge < charge[i]) {
            max_charge = charge[i];
            max_index = i;
        }
    }
    
    if (min_index != max_index) {
        for (int i = 0; i < 3; i++) {
            if (i == min_index) continue;
            if (i == max_index) continue;
            med_index = i;
        }
    } else {
        min_index = 0;
        med_index = 1;
        max_index = 2;
    }
    
    // Calculate asymmetries
    double max_asy = 0;
    if (charge[med_index] + charge[max_index] > 0) {
        max_asy = std::abs(charge[med_index] - charge[max_index]) / 
                  (charge[med_index] + charge[max_index]);
    }
    
    // Calculate overall charge
    double overall_charge = (weight[0]*charge[0] + weight[1]*charge[1] + weight[2]*charge[2]) /
                           (weight[0] + weight[1] + weight[2]);
    
    // Exclude maximal charge if asymmetry is too large
    if (max_asy > 0.04) {
        overall_charge = (weight[med_index] * charge[med_index] + 
                         weight[min_index] * charge[min_index]) /
                        (weight[med_index] + weight[min_index]);
    }
    
    // Convert charge to kinetic energy
    // Using W-value of 23.6 eV per electron-ion pair
    kine_energy = overall_charge / recom_factor / fudge_factor * 23.6 / 1e6 * units::MeV;
    
    return kine_energy;
}