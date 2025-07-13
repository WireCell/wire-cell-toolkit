
    #include "improvecluster_1.h"

WIRECELL_FACTORY(ImproveCluster_1, WireCell::Clus::ImproveCluster_1,
                 WireCell::IConfigurable, WireCell::IPCTreeMutate)

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
    {
    }

    ImproveCluster_1::~ImproveCluster_1() 
    {
    }

    void ImproveCluster_1::configure(const WireCell::Configuration& cfg)
    {
        // Configure base class first
        RetileCluster::configure(cfg);
        
        NeedDV::configure(cfg);
        NeedPCTS::configure(cfg);

        if (cfg.isMember("samplers") && cfg["samplers"].isArray()) {
            // Process array of samplers
            for (const auto& sampler_cfg : cfg["samplers"]) {
                int apa = sampler_cfg["apa"].asInt();
                int face = sampler_cfg["face"].asInt();
                std::string sampler_name = sampler_cfg["name"].asString();
                
                if (sampler_name.empty()) {
                    raise<ValueError>("RetileCluster requires an IBlobSampler name for APA %d face %d", apa, face);
                }
                // std::cout << "Test: " << apa << " " << face << " " << sampler_name << std::endl;
                auto sampler_ptr = Factory::find_tn<IBlobSampler>(sampler_name);
                m_samplers[apa][face] = sampler_ptr;
            }
        }

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
    }

    Configuration ImproveCluster_1::default_configuration() const
    {
        Configuration cfg = RetileCluster::default_configuration();
        
      
        
        return cfg;
    }


    std::unique_ptr<ImproveCluster_1::node_t> ImproveCluster_1::mutate(node_t& node) const
    {
        // get the original cluster
        auto* orig_cluster = reinitialize(node);
        

        // std::cout << m_grouping->get_name() << " " << m_wpid_angles.size() << std::endl;

        auto wpids = orig_cluster->wpids_blob();
        std::set<WirePlaneId> wpid_set(wpids.begin(), wpids.end());

        // Needed in hack_activity() but call it here to avoid call overhead.
        // find the highest and lowest points
        auto pair_points = orig_cluster->get_two_boundary_wcps();
        auto first_index  =   orig_cluster->get_closest_point_index(pair_points.first);
        auto second_index =   orig_cluster->get_closest_point_index(pair_points.second);
        std::vector<size_t> path_wcps = orig_cluster->graph_algorithms("basic_pid").shortest_path(first_index, second_index);


        // make a new node from the existing grouping
        auto& new_cluster = m_grouping->make_child(); // make a new cluster inside the existing grouping ...


        std::cout << "Xin3: " << path_wcps.size() << " " << pair_points.first.x() << " " 
                << pair_points.first.y() << " " 
                << pair_points.first.z() << " | "
                << pair_points.second.x() << " " 
                << pair_points.second.y() << " " 
                << pair_points.second.z() << std::endl;
            

        for (auto it = wpid_set.begin(); it != wpid_set.end(); ++it) {
            int apa = it->apa();
            int face = it->face();
            const auto& angles = m_wpid_angles.at(*it);

            std::map<std::pair<int, int>, std::vector<WRG::measure_t> > map_slices_measures;
            
            // std::map<std::pair<int, int>, std::vector<WRG::measure_t> > map_slices_measures_orig;
            // get_activity(*orig_cluster, map_slices_measures_orig, apa, face);

            get_activity_improved(*orig_cluster, map_slices_measures, apa, face);

            // Step 2.
            hack_activity_improved(*orig_cluster, map_slices_measures, path_wcps, apa, face); // may need more args


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
            auto iblobs = make_iblobs(map_slices_measures, apa, face);

            std::cout << orig_cluster->nchildren() << " " << iblobs.size() << " iblobs for apa " << apa << " face " << face << std::endl;

            auto niblobs = iblobs.size();
            
            // start to sampling points 

            for (size_t bind=0; bind<niblobs; ++bind) {
          
                const IBlob::pointer iblob = iblobs[bind];
                auto sampler = m_samplers.at(apa).at(face);
                const double tick = m_grouping->get_tick().at(apa).at(face);

                auto pcs = Aux::sample_live(sampler, iblob, angles, tick, bind);
                // DO NOT EXTEND FURTHER! see #426, #430

                // if (pcs["3d"].size()==0) continue; // no points ...
                // // Access 3D coordinates
                // auto pc3d = pcs["3d"];  // Get the 3D point cloud dataset
                // auto x_coords = pc3d.get("x")->elements<double>();  // Get X coordinates
                // auto y_coords = pc3d.get("y")->elements<double>();  // Get Y coordinates  
                // auto z_coords = pc3d.get("z")->elements<double>();  // Get Z coordinates
                // auto ucharge_val = pc3d.get("ucharge_val")->elements<double>();  // Get U charge
                // auto vcharge_val = pc3d.get("vcharge_val")->elements<double>();  // Get V charge
                // auto wcharge_val = pc3d.get("wcharge_val")->elements<double>();  // Get W charge
                // auto ucharge_err = pc3d.get("ucharge_unc")->elements<double>();  // Get U charge error
                // auto vcharge_err = pc3d.get("vcharge_unc")->elements<double>();  // Get V charge error
                // auto wcharge_err = pc3d.get("wcharge_unc")->elements<double>();  // Get W charge error

                // std::cout << "Xin4: " << pcs.size() << " " 
                //         << pcs["3d"].size() << " " 
                //         << x_coords.size() << " " 
                //         << y_coords.size() << " "
                //         << z_coords.size() << " "
                //         << ucharge_val.size() << " "
                //         << vcharge_val.size() << " "
                //         << wcharge_val.size() << " " 
                //         << ucharge_err.size() << " " 
                //         << vcharge_err.size() << " " 
                //         << wcharge_err.size() << " " 
                //         << std::endl;

                if (pcs.empty()) {
                    SPDLOG_DEBUG("ImproveCluster_1: skipping blob {} with no points", iblob->ident());
                    continue;
                }
                new_cluster.node()->insert(Tree::Points(std::move(pcs)));

            }
 

            // remove bad blobs ...
            int tick_span = map_slices_measures.begin()->first.second -  map_slices_measures.begin()->first.first;
            auto blobs_to_remove = remove_bad_blobs(*orig_cluster, new_cluster, tick_span, apa, face);
            for (const Blob* blob : blobs_to_remove) {
                Blob& b = const_cast<Blob&>(*blob);
                new_cluster.remove_child(b);
            }
            std::cout << "Xin5: " << blobs_to_remove.size() << " blobs removed for apa " << apa << " face " << face << " " << new_cluster.children().size() << std::endl;
        }


        auto& default_scope = orig_cluster->get_default_scope();
        auto& raw_scope = orig_cluster->get_raw_scope();

        std::cout << "Xin6: " << default_scope.hash() << " " << raw_scope.hash() << std::endl;
        if (default_scope.hash()!=raw_scope.hash()){
            auto correction_name = orig_cluster->get_scope_transform(default_scope);
            // std::vector<int> filter_results = c
            new_cluster.add_corrected_points(m_pcts, correction_name);
            // Get the new scope with corrected points
            const auto& correction_scope = new_cluster.get_scope(correction_name);
            // Set this as the default scope for viewing
            new_cluster.from(*orig_cluster); // copy state from original cluster
            // std::cout << "Test: Same:" << default_scope.hash() << " " << raw_scope.hash() << std::endl; 
        }


        auto retiled_node = new_cluster.node();

        // std::cout << m_grouping->get_name() << " " << m_grouping->children().size() << std::endl;

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


    // //print out for debug ...
    // std::cout << min_time << " " << max_time << " "
    //           << min_uch << " " << max_uch << " "
    //           << min_vch << " " << max_vch << " "
    //           << min_wch << " " << max_wch << " " << dead_uchs_range.size() << " " << dead_vchs_range.size() << " " << dead_wchs_range.size() << " " << map_u_tcc.size() << " " << map_v_tcc.size() << " " << map_w_tcc.size() << std::endl;


    // Maps for tracking time slices and channels for each wire plane
    std::map<int, std::set<int>> u_time_chs; // U plane time-channel map
    std::map<int, std::set<int>> v_time_chs; // V plane time-channel map  
    std::map<int, std::set<int>> w_time_chs; // W plane time-channel map

    int tick_span;

    // Step 1: Fill maps according to existing blobs in cluster
    auto children = cluster.children();
    for (auto child : children) {
        auto blob = child->value().facade<Blob>();
        if (!blob) continue;
        
        // Get the time slice bounds for this blob
        int time_slice_min = blob->slice_index_min();
        int time_slice_max = blob->slice_index_max();
        tick_span = time_slice_max - time_slice_min;
        
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
         
    // Step 2: Handle dead channels from CTPC (using grouping interface)
    for (const auto& [start, end]: dead_uchs_range) {
        for (int ch = start; ch < end; ++ch) {
            for (int time_slice = min_time; time_slice < max_time; time_slice+=tick_span) {
                auto [x_pos, y_pos] = grouping->convert_time_ch_2Dpoint(time_slice, ch, apa, face, 0);
                std::vector<float_t> query_point = {static_cast<float_t>(x_pos), static_cast<float_t>(y_pos)};
                const auto& skd = cluster.kd2d(apa, face, 0);
                auto ret_matches = skd.knn(1, query_point);
                // std::cout << ret_matches[0].first << " " << sqrt(ret_matches[0].second) / units::cm << " " << dis_cut / units::cm << std::endl;
                if (sqrt(ret_matches[0].second) < dis_cut) u_time_chs[time_slice].insert(ch);
            }
        }
    }
    for (const auto& [start, end]: dead_vchs_range) {
        for (int ch = start; ch < end; ++ch) {
            for (int time_slice = min_time; time_slice < max_time; time_slice+=tick_span) {
                auto [x_pos, y_pos] = grouping->convert_time_ch_2Dpoint(time_slice, ch, apa, face, 1);
                std::vector<float_t> query_point = {static_cast<float_t>(x_pos), static_cast<float_t>(y_pos)};
                const auto& skd = cluster.kd2d(apa, face, 1);
                auto ret_matches = skd.knn(1, query_point);
                if (sqrt(ret_matches[0].second) < dis_cut) v_time_chs[time_slice].insert(ch);
            }
        }
    }
    for (const auto& [start, end]: dead_wchs_range) {
        for (int ch = start; ch < end; ++ch) {
            for (int time_slice = min_time; time_slice < max_time; time_slice+=tick_span) {
                auto [x_pos, y_pos] = grouping->convert_time_ch_2Dpoint(time_slice, ch, apa, face, 2);
                std::vector<float_t> query_point = {static_cast<float_t>(x_pos), static_cast<float_t>(y_pos)};
                const auto& skd = cluster.kd2d(apa, face, 2);
                auto ret_matches = skd.knn(1, query_point);
                if (sqrt(ret_matches[0].second) < dis_cut) w_time_chs[time_slice].insert(ch);
            }
        }
    }
   
   
    
    // Step 3: Deal with good channels from CTPC
    std::map<std::pair<int, int>, double> time_ch_charge_map;

    // Process U plane good channels
    for (const auto& [time_ch, charge_info] : map_u_tcc) {
        int time_slice = time_ch.first;
        int ch = time_ch.second;
        auto [x_pos, y_pos] = grouping->convert_time_ch_2Dpoint(time_slice, ch, apa, face, 0);
        std::vector<float_t> query_point = {static_cast<float_t>(x_pos), static_cast<float_t>(y_pos)};
        const auto& skd = cluster.kd2d(apa, face, 0);
        auto ret_matches = skd.knn(1, query_point);
        double temp_min_dis = sqrt(ret_matches[0].second);
        if (temp_min_dis > dis_cut) continue;
        u_time_chs[time_slice].insert(ch);
    }
    // Process V plane good channels
    for (const auto& [time_ch, charge_info] : map_v_tcc) {
        int time_slice = time_ch.first;
        int ch = time_ch.second;
        auto [x_pos, y_pos] = grouping->convert_time_ch_2Dpoint(time_slice, ch, apa, face, 1);
        std::vector<float_t> query_point = {static_cast<float_t>(x_pos), static_cast<float_t>(y_pos)};
        const auto& skd = cluster.kd2d(apa, face, 1);
        auto ret_matches = skd.knn(1, query_point);
        double temp_min_dis = sqrt(ret_matches[0].second); 
        if (temp_min_dis > dis_cut) continue;
        v_time_chs[time_slice].insert(ch);
    }
    // Process W plane good channels
    for (const auto& [time_ch, charge_info] : map_w_tcc) {
        int time_slice = time_ch.first;
        int ch = time_ch.second;
        auto [x_pos, y_pos] = grouping->convert_time_ch_2Dpoint(time_slice, ch, apa, face, 2);
        std::vector<float_t> query_point = {static_cast<float_t>(x_pos), static_cast<float_t>(y_pos)};
        const auto& skd = cluster.kd2d(apa, face, 2);
        auto ret_matches = skd.knn(1, query_point);
        double temp_min_dis = sqrt(ret_matches[0].second);
        if (temp_min_dis > dis_cut) continue;
        w_time_chs[time_slice].insert(ch);
    }

    // Step 4: Convert to toolkit activity format (RayGrid measures)
    const int nlayers = 2+3;
    for (const auto& [time_slice, ch_set] : u_time_chs) {
        auto slice_key = std::make_pair(time_slice, time_slice + tick_span);

        auto& measures = map_slices_measures[slice_key];
         if (measures.size()==0){
            measures.resize(nlayers);
            // what to do the first two views???
            measures[0].push_back(1);
            measures[1].push_back(1);
            measures[2].resize(m_plane_infos.at(apa).at(face)[0].total_wires, 0);
            measures[3].resize(m_plane_infos.at(apa).at(face)[1].total_wires, 0);
            measures[4].resize(m_plane_infos.at(apa).at(face)[2].total_wires, 0);
            
            // std::cout << "Test2: " << measures[2].size() << " " << measures[3].size() << " " << measures[4].size() << std::endl;
        }

        WRG::measure_t& m = measures[2+0]; // U plane is layer 2
        for (int ch : ch_set) {
            double charge = 1e-3; // Default charge value
            auto it = map_u_tcc.find(std::make_pair(time_slice, ch));
            if (it != map_u_tcc.end()) {
                charge = it->second.second; // Use the charge from the map
            }
            m[ch] = charge;
        }
    }
    for (const auto& [time_slice, ch_set] : v_time_chs) {
        auto slice_key = std::make_pair(time_slice, time_slice + tick_span);

        auto& measures = map_slices_measures[slice_key];
         if (measures.size()==0){
            measures.resize(nlayers);
            // what to do the first two views???
            measures[0].push_back(1);
            measures[1].push_back(1);
            measures[2].resize(m_plane_infos.at(apa).at(face)[0].total_wires, 0);
            measures[3].resize(m_plane_infos.at(apa).at(face)[1].total_wires, 0);
            measures[4].resize(m_plane_infos.at(apa).at(face)[2].total_wires, 0);
            
            // std::cout << "Test3: " << measures[2].size() << " " << measures[3].size() << " " << measures[4].size() << std::endl;
        }

        WRG::measure_t& m = measures[2+1]; // V plane is layer 3
        for (int ch : ch_set) {
            double charge = 1e-3; // Default charge value
            auto it = map_v_tcc.find(std::make_pair(time_slice, ch));
            if (it != map_v_tcc.end()) {
                charge = it->second.second; // Use the charge from the map
            }
            m[ch] = charge;
        }
    }
    for (const auto& [time_slice, ch_set] : w_time_chs) {
        auto slice_key = std::make_pair(time_slice, time_slice + tick_span);
        auto& measures = map_slices_measures[slice_key];
         if (measures.size()==0){
            measures.resize(nlayers);
            // what to do the first two views???
            measures[0].push_back(1);
            measures[1].push_back(1);
            measures[2].resize(m_plane_infos.at(apa).at(face)[0].total_wires, 0);
            measures[3].resize(m_plane_infos.at(apa).at(face)[1].total_wires, 0);
            measures[4].resize(m_plane_infos.at(apa).at(face)[2].total_wires, 0);
            
            // std::cout << "Test4: " << measures[2].size() << " " << measures[3].size() << " " << measures[4].size() << std::endl;
        }
        WRG::measure_t& m = measures[2+2]; // W plane is layer 4
        for (int ch : ch_set) {
            double charge = 1e-3; // Default charge value
            auto it = map_w_tcc.find(std::make_pair(time_slice, ch));
            if (it != map_w_tcc.end()) {
                charge = it->second.second; // Use the charge from the map
            }
            m[ch] = charge;
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
                        wire_hits[plane] += (delta == 0) ? 1 : (delta == -1) ? 2 : 1;
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


std::set<const WireCell::Clus::Facade::Blob*> 
ImproveCluster_1::remove_bad_blobs(const Cluster& cluster, Cluster& shad_cluster, int tick_span, int apa, int face) const
{
     // Get time-organized maps of original and new blobs
    const auto& orig_time_blob_map = cluster.time_blob_map().at(apa).at(face);
    const auto& new_time_blob_map = shad_cluster.time_blob_map().at(apa).at(face);
    
    // Build index mappings for new blobs (similar to prototype's mcell indexing)
    std::map<int, const Blob*> map_index_blob;
    std::map<const Blob*, int> map_blob_index;
    std::vector<const Blob*> all_new_blobs;
    
    int index = 0;
    for (const auto& [time_slice, new_blobs] : new_time_blob_map) {
        for (const Blob* blob : new_blobs) {
            map_index_blob[index] = blob;
            map_blob_index[blob] = index;
            all_new_blobs.push_back(blob);
            index++;
        }
    }
    
    // If no new blobs or only one blob, return empty set (no graph needed)
    if (all_new_blobs.size() <= 1) {
        return std::set<const Blob*>();
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
    
    std::set<const Blob*> blobs_to_remove;
    
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
        
        // Collect blobs from bad components for removal
        for (int i = 0; i < static_cast<int>(component.size()); ++i) {
            int comp_id = component[i];
            if (good_components.find(comp_id) == good_components.end()) {
                // This component is not good, mark its blobs for removal
                const Blob* blob = map_index_blob[i];
                blobs_to_remove.insert(blob);
            }
        }
    }
    
    return blobs_to_remove;
}


} // namespace WireCell::Clus

