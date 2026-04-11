
#include "improvecluster_1.h"
#include <chrono>

WIRECELL_FACTORY(ImproveCluster_1, WireCell::Clus::ImproveCluster_1,
                 WireCell::INamed, WireCell::IConfigurable, WireCell::IPCTreeMutate)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;

// Segregate this weird choice for namespace.
namespace WCF = WireCell::Clus::Facade;

// Nick name for less typing.
namespace WRG = WireCell::RayGrid;

namespace WireCell::Clus {

    ImproveCluster_1::ImproveCluster_1()
        : Aux::Logger("ImproveCluster_1", "clus")
    {
    }

    ImproveCluster_1::~ImproveCluster_1() 
    {
    }

    void ImproveCluster_1::configure(const WireCell::Configuration& cfg)
    {
        // Base class configure() handles NeedDV, NeedPCTS, samplers, and anodes.
        RetileCluster::configure(cfg);
    }

    Configuration ImproveCluster_1::default_configuration() const
    {
        Configuration cfg = RetileCluster::default_configuration();
        
      
        
        return cfg;
    }


    std::unique_ptr<ImproveCluster_1::node_t> ImproveCluster_1::mutate(node_t& node) const
    {
        using Clock = std::chrono::steady_clock;
        using MS = std::chrono::duration<double, std::milli>;
        auto t_mutate_start = Clock::now();
        auto t0 = Clock::now();

        // get the original cluster
        auto* orig_cluster = reinitialize(node);
        SPDLOG_LOGGER_DEBUG(log, "timing: reinitialize took {} ms", MS(Clock::now()-t0).count());
        

        // std::cout << m_grouping->get_name() << " " << m_wpid_angles.size() << std::endl;

        const auto wpid_set = orig_cluster->wpids_blob_set();

        // // Needed in hack_activity() but call it here to avoid call overhead.
        // // find the highest and lowest points
        // auto pair_points = orig_cluster->get_two_boundary_wcps();
        // auto first_index  =   orig_cluster->get_closest_point_index(pair_points.first);
        // auto second_index =   orig_cluster->get_closest_point_index(pair_points.second);
        // std::vector<size_t> path_wcps = orig_cluster->graph_algorithms("basic_pid").shortest_path(first_index, second_index);


        // make a new node from the existing grouping
        auto& new_cluster = m_grouping->make_child(); // make a new cluster inside the existing grouping ...


        // std::cout << "Xin3: " << path_wcps.size() << " " << pair_points.first.x() << " " 
        //         << pair_points.first.y() << " " 
        //         << pair_points.first.z() << " | "
        //         << pair_points.second.x() << " " 
        //         << pair_points.second.y() << " " 
        //         << pair_points.second.z() << std::endl;
            

        for (auto it = wpid_set.begin(); it != wpid_set.end(); ++it) {
            int apa = it->apa();
            int face = it->face();
            const auto& angles = m_wpid_angles.at(*it);

            std::map<std::pair<int, int>, std::vector<WRG::measure_t> > map_slices_measures;
            
            // std::map<std::pair<int, int>, std::vector<WRG::measure_t> > map_slices_measures_orig;
            // get_activity(*orig_cluster, map_slices_measures_orig, apa, face);

            t0 = Clock::now();
            get_activity_improved(*orig_cluster, map_slices_measures, apa, face);
            SPDLOG_LOGGER_DEBUG(log, "timing: get_activity_improved (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());

            // Step 2 (hack_activity_improved) is intentionally NOT called here.
            // The prototype's Improve_PR3DCluster_1 adds dead/good channels only — it
            // has no path-tube logic.  Path-tube hacking is performed exclusively by
            // ImproveCluster_2::mutate, which calls hack_activity_improved twice.

            // test ...
            // std::cout << "Test: Improved: " << map_slices_measures.size() << " " << orig_cluster->children().size() << std::endl;
            // for (const auto& [slice_key, measures] : map_slices_measures) {
            //     std::cout << "Slice: [" << slice_key.first << ", " << slice_key.second << ") ";
            //     for (size_t i = 2; i < 5; ++i) {
            //         bool in_range = false;
            //         int start_idx = -1;
            //         std::cout << "Layer " << i << " ";
            //         for (size_t idx = 0; idx < measures[i].size(); ++idx) {
            //             if (measures[i][idx] > 0.0) {
            //                 if (!in_range) {
            //                     start_idx = idx;
            //                     in_range = true;
            //                 }
            //             } else {
            //                 if (in_range) {
            //                     std::cout << "[" << start_idx << ", " << idx-1 << ") ";
            //                     in_range = false;
            //                 }
            //             }
            //         }
            //         if (in_range) {
            //             std::cout << "[" << start_idx << ", " << measures[i].size()-1 << ") ";
            //         }
            //     }
            //     std::cout << std::endl;
            // }



            // Step 3.
            t0 = Clock::now();
            auto iblobs = make_iblobs_improved(map_slices_measures, apa, face);
            SPDLOG_LOGGER_DEBUG(log, "timing: make_iblobs_improved (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());
            SPDLOG_LOGGER_DEBUG(log, "{} blobs -> {} iblobs for apa {} face {}", orig_cluster->nchildren(), iblobs.size(), apa, face);

            auto niblobs = iblobs.size();
            
            // start to sampling points 
            int npoints = 0;
            t0 = Clock::now();
            for (size_t bind=0; bind<niblobs; ++bind) {
          
                const IBlob::pointer iblob = iblobs[bind];
                auto sampler = m_samplers.at(apa).at(face);
                const double tick = m_grouping->get_tick().at(apa).at(face);

                auto pcs = Aux::sample_live(sampler, iblob, angles, tick, bind);
                // DO NOT EXTEND FURTHER! see #426, #430

                if (pcs["3d"].size()==0) continue; // no points ...
                // Access 3D coordinates
                auto pc3d = pcs["3d"];  // Get the 3D point cloud dataset
                auto x_coords = pc3d.get("x")->elements<double>();  // Get X coordinates
                // auto y_coords = pc3d.get("y")->elements<double>();  // Get Y coordinates  
                // auto z_coords = pc3d.get("z")->elements<double>();  // Get Z coordinates
                // auto ucharge_val = pc3d.get("ucharge_val")->elements<double>();  // Get U charge
                // auto vcharge_val = pc3d.get("vcharge_val")->elements<double>();  // Get V charge
                // auto wcharge_val = pc3d.get("wcharge_val")->elements<double>();  // Get W charge
                // auto ucharge_err = pc3d.get("ucharge_unc")->elements<double>();  // Get U charge error
                // auto vcharge_err = pc3d.get("vcharge_unc")->elements<double>();  // Get V charge error
                // auto wcharge_err = pc3d.get("wcharge_unc")->elements<double>();  // Get W charge error

                // std::cout << "ImproveCluster_1 PCS: " << pcs.size() << " " 
                //           << pcs["3d"].size() << " " 
                //           << x_coords.size() << std::endl;     

                npoints +=x_coords.size();
                if (pcs.empty()) {
                    SPDLOG_DEBUG("ImproveCluster_1: skipping blob {} with no points", iblob->ident());
                    continue;
                }
                new_cluster.node()->insert(Tree::Points(std::move(pcs)));

            }
            SPDLOG_LOGGER_DEBUG(log, "timing: sample_live loop (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());
            SPDLOG_LOGGER_DEBUG(log, "{} points sampled for apa {} face {} Blobs {}", npoints, apa, face, niblobs);


            // remove bad blobs ...
            t0 = Clock::now();
            if (map_slices_measures.empty()) continue; // no tiled blobs for this face
            int tick_span = map_slices_measures.begin()->first.second -  map_slices_measures.begin()->first.first;
            auto blobs_to_remove = remove_bad_blobs(*orig_cluster, new_cluster, tick_span, apa, face);
            for (const Blob* blob : blobs_to_remove) {
                Blob& b = const_cast<Blob&>(*blob);
                new_cluster.remove_child(b);
            }
            SPDLOG_LOGGER_DEBUG(log, "timing: remove_bad_blobs (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());
            SPDLOG_LOGGER_DEBUG(log, "{} blobs removed for apa {} face {} remaining {}", blobs_to_remove.size(), apa, face, new_cluster.children().size());
        }


        auto& default_scope = orig_cluster->get_default_scope();
        auto& raw_scope = orig_cluster->get_raw_scope();

        SPDLOG_LOGGER_DEBUG(log, "Scope: {} {}", default_scope.hash(), raw_scope.hash());
        if (default_scope.hash()!=raw_scope.hash()){
            t0 = Clock::now();
            auto correction_name = orig_cluster->get_scope_transform(default_scope);
            // std::vector<int> filter_results = c
            new_cluster.add_corrected_points(m_pcts, correction_name);
            // Set this as the default scope for viewing
            new_cluster.from(*orig_cluster); // copy state from original cluster
            SPDLOG_LOGGER_DEBUG(log, "timing: add_corrected_points took {} ms", MS(Clock::now()-t0).count());
            // std::cout << "Test: Same:" << default_scope.hash() << " " << raw_scope.hash() << std::endl; 
        }


        // auto retiled_node = new_cluster.node();

        // std::cout << m_grouping->get_name() << " " << m_grouping->children().size() << std::endl;

        SPDLOG_LOGGER_DEBUG(log, "timing: mutate() TOTAL took {} ms", MS(Clock::now()-t_mutate_start).count());
        return m_grouping->remove_child(new_cluster);
    }






void ImproveCluster_1::get_activity_improved(const Cluster& cluster, std::map<std::pair<int, int>,std::vector<WireCell::RayGrid::measure_t>>& map_slices_measures, int apa, int face) const{

    auto uvwt_min = cluster.get_uvwt_min(apa, face);
    auto uvwt_max = cluster.get_uvwt_max(apa, face);
    // Track the bounds for optimization
    int min_time = std::get<3>(uvwt_min);
    int max_time = std::get<3>(uvwt_max) + 1;
    int min_uch = std::get<0>(uvwt_min), max_uch = std::get<0>(uvwt_max) + 1;
    int min_vch = std::get<1>(uvwt_min), max_vch = std::get<1>(uvwt_max) + 1;
    int min_wch = std::get<2>(uvwt_min), max_wch = std::get<2>(uvwt_max) + 1;

    // get grouping information
    auto grouping = cluster.grouping();

    // Note: In toolkit, grouping provides methods to get dead channels
    // although this is cham actually wire index ...
    auto dead_uchs_range  = grouping->get_overlap_dead_chs(min_time, max_time, min_uch, max_uch, apa, face, 0, 0);
    auto dead_vchs_range  = grouping->get_overlap_dead_chs(min_time, max_time, min_vch, max_vch, apa, face, 1, 0);
    auto dead_wchs_range  = grouping->get_overlap_dead_chs(min_time, max_time, min_wch, max_wch, apa, face, 2, 0);

    // auto dead_uchs_all = grouping->get_all_dead_chs(apa, face, 0);
    // std::cout << "dead_uchs_all: ";
    // for (const auto& [ch, ranges ]: dead_uchs_all) {
    //     std::cout << ch << " " << ranges.first << " " << ranges.second << std::endl;
    // }
    // std::cout << std::endl;

    // std::cout << "dead_uch_ranges: " << std::endl;
    // for (const auto& [start, end] : dead_uchs_range) {
    //     std::cout << "[" << start << ", " << end << ") " << std::endl;
    // }
    // std::cout << "dead_vch_ranges: " << std::endl;
    // for (const auto& [start, end] : dead_vchs_range) {
    //     std::cout << "[" << start << ", " << end << ") " << std::endl;
    // }
    // std::cout << "dead_wch_ranges: " << std::endl;
    // for (const auto& [start, end] : dead_wchs_range) {
    //     std::cout << "[" << start << ", " << end << ") " << std::endl;
    // }


    // althoguh ch, but wire index ...
    std::map<std::pair<int,int>, std::pair<double,double>> map_u_tcc = grouping->get_overlap_good_ch_charge(min_time, max_time, min_uch, max_uch, apa, face, 0);
    std::map<std::pair<int,int>, std::pair<double,double>> map_v_tcc = grouping->get_overlap_good_ch_charge(min_time, max_time, min_vch, max_vch, apa, face, 1);
    std::map<std::pair<int,int>, std::pair<double,double>> map_w_tcc = grouping->get_overlap_good_ch_charge(min_time, max_time, min_wch, max_wch, apa, face, 2);


    // for (const auto& [time_ch, charge_info] : map_u_tcc) {
    //     std::cout << "U plane: time_slice=" << time_ch.first 
    //               << ", ch=" << time_ch.second 
    //               << ", charge=" << charge_info.first 
    //               << ", error=" << charge_info.second << std::endl;
    // }

    // //print out for debug ...
    // std::cout << min_time << " " << max_time << " "
    //           << min_uch << " " << max_uch << " "
    //           << min_vch << " " << max_vch << " "
    //           << min_wch << " " << max_wch << " " << dead_uchs_range.size() << " " << dead_vchs_range.size() << " " << dead_wchs_range.size() << " " << map_u_tcc.size() << " " << map_v_tcc.size() << " " << map_w_tcc.size() << std::endl;


    // Maps for tracking time slices and channels for each wire plane
    std::map<int, std::set<int>> u_time_chs; // U plane time-channel map
    std::map<int, std::set<int>> v_time_chs; // V plane time-channel map  
    std::map<int, std::set<int>> w_time_chs; // W plane time-channel map

    // Derive tick_span from face metadata rather than blob geometry so that it
    // is stable even if no blobs exist yet for this face, and to avoid the
    // "last-blob wins" bug that occurred before Bug #3 was fixed.
    const int tick_span = m_grouping->get_nticks_per_slice().at(apa).at(face);

    // Step 1: Fill maps according to existing blobs in cluster (this face only).
    auto children = cluster.children();
    for (auto child : children) {
        auto blob = child->value().facade<Blob>();
        if (!blob) continue;

        // Skip blobs belonging to a different (apa, face) — wire indices are
        // face-local and must not be mixed across faces.
        auto blob_wpid = blob->wpid();
        if (blob_wpid.apa() != apa || blob_wpid.face() != face) continue;

        // Get the time slice bounds for this blob
        int time_slice_min = blob->slice_index_min();
        int time_slice_max = blob->slice_index_max();
        
        // Process each time slice in the blob
        for (int time_slice = time_slice_min; time_slice < time_slice_max; time_slice = time_slice + tick_span) {
            // Initialize channel sets if not present
            if (u_time_chs.find(time_slice) == u_time_chs.end()) {
                u_time_chs[time_slice] = std::set<int>();
                v_time_chs[time_slice] = std::set<int>();
                w_time_chs[time_slice] = std::set<int>();
            }
            // Process each wire plane (U=0, V=1, W=2)
            for (int plane = 0; plane < 3; ++plane) {
                // Get wire bounds for this plane in the blob
                int wire_min = (plane == 0) ? blob->u_wire_index_min() :
                              (plane == 1) ? blob->v_wire_index_min() : blob->w_wire_index_min();
                int wire_max = (plane == 0) ? blob->u_wire_index_max() :
                              (plane == 1) ? blob->v_wire_index_max() : blob->w_wire_index_max();
                
                // Process each wire in the range
                for (int wire_ch = wire_min; wire_ch < wire_max; ++wire_ch) {
                    
                    // Store in appropriate plane map
                    if (plane == 0) { // U plane
                        u_time_chs[time_slice].insert(wire_ch);
                    } else if (plane == 1) { // V plane  
                        v_time_chs[time_slice].insert(wire_ch);
                    } else { // W plane
                        w_time_chs[time_slice].insert(wire_ch);
                    }
                }
            }
        }
    }


    //  std::cout << u_time_chs.size() << " " << v_time_chs.size() << " " << w_time_chs.size() << " " << u_time_chs.begin()->second.size() << " " << v_time_chs.begin()->second.size() << " " << w_time_chs.begin()->second.size() << std::endl;

    // Distance cut for dead channel inclusion (20 cm as in original code)
    const double dis_cut = 20 * units::cm;

    // NOTE: ch values here are wire indices, not channel IDs.  On MicroBooNE
    // they coincide; on detectors with non-trivial wire→channel mappings (e.g.
    // wrapped wires) the dead/good-channel fill will be incorrect until
    // get_overlap_dead_chs / get_overlap_good_ch_charge are updated to accept
    // wire indices natively.  See §6.5 of the port review.

    // Convenience aliases indexed by plane (0=U, 1=V, 2=W).
    const std::vector<std::pair<int,int>>* dead_ch_ranges[3] = {
        &dead_uchs_range, &dead_vchs_range, &dead_wchs_range
    };
    std::map<int, std::set<int>>* time_chs[3] = {
        &u_time_chs, &v_time_chs, &w_time_chs
    };
    const std::map<std::pair<int,int>, std::pair<double,double>>* tcc_maps[3] = {
        &map_u_tcc, &map_v_tcc, &map_w_tcc
    };

    // Step 2: Handle dead channels — one loop over planes replaces three identical blocks.
    for (int pl = 0; pl < 3; ++pl) {
        for (const auto& [start, end] : *dead_ch_ranges[pl]) {
            for (int ch = start; ch < end; ++ch) {
                for (int time_slice = min_time; time_slice < max_time; time_slice += tick_span) {
                    auto [x_pos, y_pos] = grouping->convert_time_wire_2Dpoint(time_slice, ch, apa, face, pl);
                    std::vector<float_t> query_point = {static_cast<float_t>(x_pos), static_cast<float_t>(y_pos)};
                    auto ret_matches = cluster.kd2d(apa, face, pl).knn(1, query_point);
                    if (sqrt(ret_matches[0].second) < dis_cut) (*time_chs[pl])[time_slice].insert(ch);
                }
            }
        }
    }

    // Step 3: Handle good channels from CTPC — one loop over planes.
    for (int pl = 0; pl < 3; ++pl) {
        for (const auto& [time_ch, charge_info] : *tcc_maps[pl]) {
            int time_slice = time_ch.first;
            int ch = time_ch.second;
            auto [x_pos, y_pos] = grouping->convert_time_wire_2Dpoint(time_slice, ch, apa, face, pl);
            std::vector<float_t> query_point = {static_cast<float_t>(x_pos), static_cast<float_t>(y_pos)};
            auto ret_matches = cluster.kd2d(apa, face, pl).knn(1, query_point);
            if (sqrt(ret_matches[0].second) > dis_cut) continue;
            (*time_chs[pl])[time_slice].insert(ch);
        }
    }

    // Step 4: Convert to toolkit activity format (RayGrid measures).
    // Layers 0 and 1 are the geometric (non-wire) ray layers in the RayGrid
    // scheme; layers 2, 3, 4 carry U, V, W wire activity respectively.
    const int nlayers = 2 + 3;
    for (int pl = 0; pl < 3; ++pl) {
        for (const auto& [time_slice, ch_set] : *time_chs[pl]) {
            auto slice_key = std::make_pair(time_slice, time_slice + tick_span);
            auto& measures = map_slices_measures[slice_key];
            if (measures.empty()) {
                measures.resize(nlayers);
                measures[0].push_back(1);
                measures[1].push_back(1);
                for (int i = 0; i < 3; ++i)
                    measures[2+i].resize(m_plane_infos.at(apa).at(face)[i].total_wires, 0);
            }
            WRG::measure_t& m = measures[2 + pl];
            for (int ch : ch_set) {
                double charge = 1e-3; // sentinel: dead/forced channel, converted to (0, large_err) in make_iblobs_improved
                auto it = tcc_maps[pl]->find(std::make_pair(time_slice, ch));
                if (it != tcc_maps[pl]->end()) charge = it->second.first;
                m[ch] = charge;
            }
        }
    }

 
}


// Step 2. Modify activity to suit.
void ImproveCluster_1::hack_activity_improved(const Cluster& cluster, std::map<std::pair<int, int>, std::vector<WRG::measure_t> >& map_slices_measures, const std::vector<size_t>& path_wcps, int apa, int face) const
{

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

    // Guard: if activity is empty there is nothing to hack.
    if (map_slices_measures.empty()) return;

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
                        wire_hits[plane] += (delta == 0) ? 2 : 1;
                    }
                }
            }
        }
        
        // Set flag if sufficient activity found
        if (wire_hits[0] >=2 && wire_hits[1] >=2 && wire_hits[2] >=2) {
            path_pts_flag[i] = true;
        }
        
        // std::cout << i << " " << path_pts[i].first << " " << path_pts_flag[i] << std::endl;

        //std::cout << path_pts[i] << " " << wire_hits[0] << " " << wire_hits[1] << " " << wire_hits[2] << " " << path_pts_flag[i] << " " << aligned_tick/tick_span << " " << u_wire << " " << v_wire << " " << w_wire << " " << time_tick_u << " " << std::round(time_tick_u / tick_span) << std::endl;
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
                         pow(dw,2) + pow(dt,2)>3*3) 
                        continue;
                    if (measures.at(wire) > 0.0) continue; // Already has activity
                    measures.at(wire) = 1.0e-3;  // Set activity
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


std::vector<const WireCell::Clus::Facade::Blob*>
ImproveCluster_1::remove_bad_blobs(const Cluster& cluster, Cluster& shad_cluster, int tick_span, int apa, int face) const
{
    // Get time-organized maps of original and new blobs
    const auto& orig_time_blob_map = cluster.time_blob_map().at(apa).at(face);
    const auto& new_time_blob_map = shad_cluster.time_blob_map().at(apa).at(face);

    // Build index mappings for new blobs.  Sort by (time_slice, blob->ident())
    // so vertex IDs are deterministic across runs regardless of heap layout.
    std::map<int, const Blob*> map_index_blob;
    std::map<const Blob*, int> map_blob_index;
    std::vector<const Blob*> all_new_blobs;

    int index = 0;
    for (const auto& [time_slice, new_blobs] : new_time_blob_map) {
        // Sort within each time slice by blob ident for determinism.
        std::vector<const Blob*> sorted_blobs(new_blobs.begin(), new_blobs.end());
        // Sort by wire-range corners for a deterministic, geometry-based order that
        // is independent of heap-allocated pointer values.
        std::sort(sorted_blobs.begin(), sorted_blobs.end(),
                  [](const Blob* a, const Blob* b) {
                      if (a->u_wire_index_min() != b->u_wire_index_min()) return a->u_wire_index_min() < b->u_wire_index_min();
                      if (a->v_wire_index_min() != b->v_wire_index_min()) return a->v_wire_index_min() < b->v_wire_index_min();
                      return a->w_wire_index_min() < b->w_wire_index_min();
                  });
        for (const Blob* blob : sorted_blobs) {
            map_index_blob[index] = blob;
            map_blob_index[blob]  = index;
            all_new_blobs.push_back(blob);
            index++;
        }
    }
    
    // If no new blobs or only one blob, nothing to filter.
    if (all_new_blobs.size() <= 1) {
        return {};
    }
    
    // Create graph for new blobs - establish connectivity between adjacent time slices
    const int N = all_new_blobs.size();
    boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                         boost::no_property, boost::property<boost::edge_weight_t, double>>
        temp_graph(N);
    
    // Build graph edges between blobs in adjacent time slices that overlap spatially
    for (const auto& [time_slice, current_blobs] : new_time_blob_map) {
        // Connect to next time slice
        auto next_it = new_time_blob_map.find(time_slice + tick_span);
        if (next_it != new_time_blob_map.end()) {
            for (const Blob* blob1 : current_blobs) {
                for (const Blob* blob2 : next_it->second) {
                    int index1 = map_blob_index[blob1];
                    int index2 = map_blob_index[blob2];
                    
                    // Add edge if blobs overlap spatially (similar to prototype's Overlap_fast)
                    if (blob1->overlap_fast(*blob2, 1)) {
                        add_edge(index1, index2, 1.0, temp_graph);
                    }
                }
            }
        }
    }
    
    // Find connected components (groups of spatially/temporally connected blobs)
    std::vector<int> component(num_vertices(temp_graph));
    const int num_components = connected_components(temp_graph, &component[0]);
    
    std::vector<const Blob*> blobs_to_remove;

    // If we have multiple disconnected components, validate each component
    if (num_components > 1) {
        std::set<int> good_components;
        
        // Examine each connected component to determine if it's "good"
        for (int i = 0; i < static_cast<int>(component.size()); ++i) {
            int comp_id = component[i];
            
            // Skip if we've already validated this component
            if (good_components.find(comp_id) != good_components.end()) {
                continue;
            }
            
            const Blob* blob = map_index_blob[i];
            int time_slice = blob->slice_index_min(); // Get time slice for this blob
            bool flag_good = false;
            
            // Check overlap with original blobs in previous time slice
            if (!flag_good) {
                auto prev_it = orig_time_blob_map.find(time_slice - tick_span);
                if (prev_it != orig_time_blob_map.end()) {
                    for (const Blob* orig_blob : prev_it->second) {
                        if (blob->overlap_fast(*orig_blob, 1)) {
                            flag_good = true;
                            break;
                        }
                    }
                }
            }
            
            // Check overlap with original blobs in same time slice
            if (!flag_good) {
                auto same_it = orig_time_blob_map.find(time_slice);
                if (same_it != orig_time_blob_map.end()) {
                    for (const Blob* orig_blob : same_it->second) {
                        if (blob->overlap_fast(*orig_blob, 1)) {
                            flag_good = true;
                            break;
                        }
                    }
                }
            }
            
            // Check overlap with original blobs in next time slice
            if (!flag_good) {
                auto next_it = orig_time_blob_map.find(time_slice + tick_span);
                if (next_it != orig_time_blob_map.end()) {
                    for (const Blob* orig_blob : next_it->second) {
                        if (blob->overlap_fast(*orig_blob, 1)) {
                            flag_good = true;
                            break;
                        }
                    }
                }
            }
            
            // If this component representative blob has good overlap, mark entire component as good
            if (flag_good) {
                good_components.insert(comp_id);
            }
        }
        
        // Collect blobs from bad components for removal.
        // Iterate in deterministic vertex-ID order (assigned from sorted blobs above).
        for (int i = 0; i < static_cast<int>(component.size()); ++i) {
            int comp_id = component[i];
            if (good_components.find(comp_id) == good_components.end()) {
                const Blob* blob = map_index_blob[i];
                blobs_to_remove.push_back(blob);
            }
        }
    }

    return blobs_to_remove;
}

std::vector<IBlob::pointer> ImproveCluster_1::make_iblobs_improved(std::map<std::pair<int, int>, std::vector<WRG::measure_t> >& map_slices_measures, int apa, int face) const
{
    std::vector<IBlob::pointer> ret;

    const auto& coords = m_face.at(apa).at(face)->raygrid();
    int blob_ident=0;
    int slice_ident = 0;

    const double tick = m_grouping->get_tick().at(apa).at(face);


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
    
        // Convert measures to ISlice activity map
        // Layers 2, 3, 4 correspond to U, V, W wire planes

        IFrame::pointer sframe = nullptr;

         // Create the slice with activity
        auto sslice = std::make_shared<Aux::SimpleSlice>(sframe, slice_ident++, it->first.first*tick, (it->first.second - it->first.first)*tick);
        // Copy the prepared activity map into the slice
        auto& slice_activity = sslice->activity();

        for (int plane_idx = 0; plane_idx < 3; ++plane_idx) {
            const int layer = plane_idx + 2;
            const auto& plane_measures = it->second[layer];
            
            // Get the wire plane for this face and plane
            auto face_ptr = m_face.at(apa).at(face);
            auto planes = face_ptr->planes();
            if (static_cast<size_t>(plane_idx) >= planes.size()) continue;
            
            auto wire_plane = planes[plane_idx];
            const auto& channels = wire_plane->channels();
            
            // Map wire indices to channels and populate activity
            for (size_t wire_idx = 0; wire_idx < plane_measures.size(); ++wire_idx) {
                if (plane_measures[wire_idx] > 0.0) {
                    // Find the channel corresponding to this wire index
                    if (wire_idx < channels.size()) {
                        auto ichan = channels[wire_idx];
                        if (ichan) {
                            // Set activity with value and zero uncertainty
                            if (plane_measures[wire_idx]==1e-3){
                                slice_activity[ichan] = ISlice::value_t(0.0, 1e12);
                            } else {
                                slice_activity[ichan] = ISlice::value_t(plane_measures[wire_idx], 0.0);
                            }
                        }
                    }
                }
            }
        }


        for (const auto& bshape : bshapes) {

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
            // ISlice::pointer slice = sslice;

            IBlob::pointer iblob = std::make_shared<Aux::SimpleBlob>(blob_ident++, blob_value,
                                                                 blob_error, bshape, sslice, m_face.at(apa).at(face));
            ret.push_back(iblob);

        }
    }

    // std::cout << "Test: Blobs: " << ret.size() << std::endl;

    return ret;
}


} // namespace WireCell::Clus

