#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/FiducialUtils.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/TrackFitting.h"  
#include "WireCellClus/TrackFittingPresets.h"
#include "WireCellClus/PRSegmentFunctions.h"

#include "WireCellIface/IScalarFunction.h"
#include "WireCellUtil/KSTest.h"




class TaggerCheckSTM;
WIRECELL_FACTORY(TaggerCheckSTM, TaggerCheckSTM,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

struct edge_base_t {
    typedef boost::edge_property_tag kind;
};

/**
 * Clustering function that checks the main cluster from clustering_recovering_bundle
 * for Short Track Muon (STM) characteristics and sets the STM flag when conditions are met.
 * This function works on clusters that have already been processed by clustering_recovering_bundle.
 */
class TaggerCheckSTM : public IConfigurable, public Clus::IEnsembleVisitor, private Clus::NeedDV, private Clus::NeedPCTS {
public:
    TaggerCheckSTM() {
        // Initialize with default preset
        m_track_fitter = TrackFittingPresets::create_with_current_values();
    }
    virtual ~TaggerCheckSTM() {}

    virtual void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config); 
        m_grouping_name = get<std::string>(config, "grouping", "live");

        m_trackfitting_config_file = get<std::string>(config, "trackfitting_config_file", "");
    
        if (!m_trackfitting_config_file.empty()) {
            std::cout << "TaggerCheckSTM: Loading TrackFitting config from: " << m_trackfitting_config_file << std::endl;
            load_trackfitting_config(m_trackfitting_config_file);
        } else {
            std::cout << "TaggerCheckSTM: No TrackFitting config file specified, using defaults" << std::endl;
        }

        // Configure the LinterpFunction - similar to how drifter is configured
        auto linterp_name = get<std::string>(config, "linterp_function", "Muon");
        if (!linterp_name.empty()) {
            m_linterp_function = Factory::find_tn<IScalarFunction>(linterp_name);
            if (!m_linterp_function) {
                std::cout << "TaggerCheckSTM: Failed to find LinterpFunction: " <<  linterp_name << std::endl;
                THROW(ValueError() << errmsg{"Failed to find LinterpFunction: " + linterp_name});
            }
            std::cout << "TaggerCheckSTM: Successfully configured LinterpFunction: " << linterp_name << std::endl;
        } else {
            std::cout << "TaggerCheckSTM: No LinterpFunction configured" << std::endl;
        }

    }
    
    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping_name;
        cfg["detector_volumes"] = "DetectorVolumes";
        cfg["pc_transforms"] = "PCTransformSet";  

        cfg["trackfitting_config_file"] = ""; 
        cfg["linterp_function"] = "";  // empty means user must provide

        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const {
        
        // Configure the track fitter with detector volume
        m_track_fitter.set_detector_volume(m_dv);
        m_track_fitter.set_pc_transforms(m_pcts); 

        // Get the specified grouping (default: "live")
        auto groupings = ensemble.with_name(m_grouping_name);
        if (groupings.empty()) {
            return;
        }
        
        auto& grouping = *groupings.at(0);
        
        // Find clusters that have the main_cluster flag (set by clustering_recovering_bundle)
        Cluster* main_cluster = nullptr;

        for (auto* cluster : grouping.children()) {
            if (cluster->get_flag(Flags::main_cluster)) {
                main_cluster = cluster;
            }
        }

        std::cout << "TaggerCheckSTM: Found " << (main_cluster ? 1 : 0)
                  << " main clusters to check for STM conditions." << std::endl;

        // For each main cluster, find its associated clusters
        std::map<Cluster*, std::vector<Cluster*>> main_to_associated;
        if (main_cluster) {
            std::vector<Cluster*> associated_clusters;
            
            // Find all clusters with the associated_cluster flag
            for (auto* cluster : grouping.children()) {
                if (cluster->get_flag(Flags::associated_cluster)) {
                    associated_clusters.push_back(cluster);
                }
            }
            
            main_to_associated[main_cluster] = associated_clusters;
            
            // std::cout << "TaggerCheckSTM: Main cluster " << main_cluster->ident() 
            //           << " has " << associated_clusters.size() << " associated clusters: ";
            // for (auto* assoc : associated_clusters) {
            //     std::cout << assoc->ident() << " ";
            // }
            // std::cout << std::endl;
        }

        // Process each main cluster
        size_t stm_count = 0;

        // validation check ... temporary ...
        {
            auto boundary_indices = main_cluster->get_two_boundary_steiner_graph_idx("steiner_graph", "steiner_pc", true);

            const auto& steiner_pc = main_cluster->get_pc("steiner_pc");
            const auto& coords = main_cluster->get_default_scope().coords;
            const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
            const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
            const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();

         // Add the two boundary points as additional extreme point groups
            geo_point_t boundary_point_first(x_coords[boundary_indices.first], 
                                        y_coords[boundary_indices.first], 
                                        z_coords[boundary_indices.first]);
            geo_point_t boundary_point_second(x_coords[boundary_indices.second], 
                                        y_coords[boundary_indices.second], 
                                        z_coords[boundary_indices.second]);
            geo_point_t first_wcp = boundary_point_first;
            geo_point_t last_wcp = boundary_point_second;

            std::cout << "End Points: " << first_wcp << " " << last_wcp << std::endl;
            auto path_points = do_rough_path(*main_cluster, first_wcp, last_wcp);

            // hack the path_points according to WCP ...
            path_points.clear();
            path_points.emplace_back(2195.39, -869.317, 2090.5);
            path_points.emplace_back(2193.19, -873.647, 2092);
            path_points.emplace_back(2190.99, -878.843, 2095);
            path_points.emplace_back(2188.79, -882.307, 2095);
            path_points.emplace_back(2186.59, -885.771, 2095);
            path_points.emplace_back(2184.38, -889.235, 2095);
            path_points.emplace_back(2182.18, -894.431, 2098);
            path_points.emplace_back(2182.18, -897.896, 2098);
            path_points.emplace_back(2179.98, -901.36, 2098);
            path_points.emplace_back(2177.78, -906.556, 2101);
            path_points.emplace_back(2177.78, -910.02, 2101);
            path_points.emplace_back(2173.37, -913.484, 2101);
            path_points.emplace_back(2171.17, -918.68, 2104);
            path_points.emplace_back(2168.97, -922.144, 2104);
            path_points.emplace_back(2168.97, -925.608, 2104);
            path_points.emplace_back(2166.77, -929.073, 2104);
            path_points.emplace_back(2164.57, -928.206, 2105.5);
            path_points.emplace_back(2164.57, -933.402, 2108.5);
            path_points.emplace_back(2162.36, -936.867, 2108.5);
            path_points.emplace_back(2160.16, -940.331, 2108.5);
            path_points.emplace_back(2160.16, -943.795, 2108.5);
            path_points.emplace_back(2157.96, -948.991, 2111.5);
            path_points.emplace_back(2155.76, -952.455, 2111.5);
            path_points.emplace_back(2153.56, -951.589, 2113);
            path_points.emplace_back(2153.56, -955.053, 2113);
            path_points.emplace_back(2153.56, -958.517, 2113);
            path_points.emplace_back(2151.35, -961.982, 2113);
            path_points.emplace_back(2149.15, -961.982, 2113);
            path_points.emplace_back(2146.95, -967.178, 2116);
            path_points.emplace_back(2146.95, -970.642, 2116);
            path_points.emplace_back(2144.75, -974.106, 2116);
            path_points.emplace_back(2142.55, -977.57, 2116);
            path_points.emplace_back(2140.34, -982.766, 2119);
            path_points.emplace_back(2138.14, -986.23, 2119);
            path_points.emplace_back(2135.94, -987.096, 2117.5);
            path_points.emplace_back(2135.94, -992.292, 2120.5);
            path_points.emplace_back(2133.74, -994.025, 2120.5);
            path_points.emplace_back(2133.74, -997.489, 2120.5);
            path_points.emplace_back(2131.54, -998.355, 2122);
            path_points.emplace_back(2131.54, -1001.82, 2122);
            path_points.emplace_back(2129.33, -1004.42, 2123.5);
            path_points.emplace_back(2127.13, -1005.28, 2122);
            path_points.emplace_back(2124.93, -1010.48, 2125);
            path_points.emplace_back(2124.93, -1012.21, 2128);
            path_points.emplace_back(2120.53, -1014.81, 2129.5);
            path_points.emplace_back(2120.53, -1018.27, 2129.5);
            path_points.emplace_back(2120.53, -1032.13, 2129.5);

            std::cout << path_points.size() << " Path Points: " << std::endl;
            // for (const auto& point : path_points) {
            //     std::cout << "Path point: x=" << point.x() << " y=" << point.y() << " z=" << point.z() << std::endl;
            // }
            // std::cout << std::endl;

            // Create segment for tracking
            auto segment = create_segment_for_cluster(*main_cluster, path_points);

            geo_point_t test_p(10,10,10);
            const auto& fit_seg_dpc = segment->dpcloud("main");
            auto closest_result = fit_seg_dpc->kd3d().knn(1, test_p);    
            double closest_3d_distance = sqrt(closest_result[0].second);
            auto closest_2d_u = fit_seg_dpc->get_closest_2d_point_info(test_p, 0, 0, 0);
            auto closest_2d_v = fit_seg_dpc->get_closest_2d_point_info(test_p, 1, 0, 0);
            auto closest_2d_w = fit_seg_dpc->get_closest_2d_point_info(test_p, 2, 0, 0);
            std::cout << closest_3d_distance << " " <<  std::get<0>(closest_2d_u) << " " << std::get<0>(closest_2d_v) << " " << std::get<0>(closest_2d_w) << std::endl;
            std::cout << std::get<2>(closest_2d_u) << " " << std::get<2>(closest_2d_v) << " " << std::get<2>(closest_2d_w) << std::endl;

            m_track_fitter.add_segment(segment);
            m_track_fitter.do_single_tracking(segment, false);
            // Extract fit results from the segment
            const auto& fits = segment->fits();
            
            // Print position, dQ, and dx for each fit point
            std::cout << "Fit results for " << fits.size() << " points:" << std::endl;
            for (size_t i = 0; i < fits.size(); ++i) {
                const auto& fit = fits[i];
                std::cout << "  Point " << i << ": position=(" 
                         << fit.point.x()/units::cm << ", " << fit.point.y()/units::cm << ", " << fit.point.z()/units::cm
                         << "), dQ=" << fit.dQ << ", dx=" << fit.dx/units::cm << std::endl;
            }
            std::cout << std::endl;

        }

        // if (check_stm_conditions(*main_cluster, main_to_associated[main_cluster] )) {
        //     main_cluster->set_flag(Flags::STM);
        //     stm_count++;
        // }
        
        (void)stm_count;
    }

private:
    std::string m_grouping_name{"live"};
    std::string m_trackfitting_config_file;  // Path to TrackFitting config file
    mutable TrackFitting m_track_fitter; 

    WireCell::IScalarFunction::pointer m_linterp_function;

    void load_trackfitting_config(const std::string& config_file) {
        try {
            // Load JSON file
            std::ifstream file(config_file);
            if (!file.is_open()) {
                std::cerr << "TaggerCheckSTM: Cannot open config file: " << config_file << std::endl;
                return;
            }
            
            Json::Value root;
            Json::CharReaderBuilder builder;
            std::string errs;
            
            if (!Json::parseFromStream(builder, file, &root, &errs)) {
                std::cerr << "TaggerCheckSTM: Failed to parse JSON: " << errs << std::endl;
                return;
            }
            
            // Apply each parameter from the JSON file
            for (const auto& param_name : root.getMemberNames()) {
                if (param_name.substr(0, 1) == "_") continue;  // Skip comments
                
                try {
                    double value = root[param_name].asDouble();
                    m_track_fitter.set_parameter(param_name, value);
                    std::cout << "TaggerCheckSTM: Set " << param_name << " = " << value << std::endl;
                } catch (const std::exception& e) {
                    std::cerr << "TaggerCheckSTM: Failed to set parameter " << param_name 
                            << ": " << e.what() << std::endl;
                }
            }
            
            std::cout << "TaggerCheckSTM: Successfully loaded TrackFitting configuration" << std::endl;
            
        } catch (const std::exception& e) {
            std::cerr << "TaggerCheckSTM: Exception loading config: " << e.what() << std::endl;
            std::cerr << "TaggerCheckSTM: Using default TrackFitting parameters" << std::endl;
        }
    }

    std::vector<geo_point_t> do_rough_path(const Cluster& cluster,geo_point_t& first_point, geo_point_t& last_point) const{
         // 1. Get Steiner point cloud and graph
        // const auto& steiner_pc = cluster.get_pc("steiner_pc");
        // const auto& steiner_graph = cluster.get_graph("steiner_graph");
        
        // 2. Find the closest point indices in the Steiner point cloud
  
        // Find closest indices in the steiner point cloud
        auto first_knn_results = cluster.kd_steiner_knn(1, first_point, "steiner_pc");
        auto last_knn_results = cluster.kd_steiner_knn(1, last_point, "steiner_pc");
        
        auto first_index = first_knn_results[0].first;  // Get the index from the first result
        auto last_index = last_knn_results[0].first;   // Get the index from the first result
 
        // 4. Use Steiner graph to find the shortest path
        const std::vector<size_t>& path_indices = 
            cluster.graph_algorithms("steiner_graph").shortest_path(first_index, last_index);
            
        std::vector<geo_point_t> path_points;
        const auto& steiner_pc = cluster.get_pc("steiner_pc");
        const auto& coords = cluster.get_default_scope().coords;
        const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
        const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
        const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();

        for (size_t idx : path_indices) {
            path_points.emplace_back(x_coords[idx], y_coords[idx], z_coords[idx]);
        }
        return path_points;
    }

    // return a vector of point, also the mid_p is also a return point ...
    std::vector<geo_point_t> adjust_rough_path(const Cluster& cluster, geo_point_t& mid_p) const{

        const geo_point_t drift_dir_abs(1,0,0); 
        // use the m_track_fitter ...
        auto fine_tracking_path = m_track_fitter.get_fine_tracking_path();
        auto dQ = m_track_fitter.get_dQ();
        auto dx = m_track_fitter.get_dx();

        mid_p.at(0) = fine_tracking_path.at(0).first.x();
        mid_p.at(1) = fine_tracking_path.at(0).first.y();
        mid_p.at(2) = fine_tracking_path.at(0).first.z();

        // Initialize variables
        int save_i = 0;
        bool flag_crawl = false;

        // Initialize angle vectors
        std::vector<double> refl_angles(fine_tracking_path.size(), 0);
        std::vector<double> para_angles(fine_tracking_path.size(), 0);

        // First part: Calculate reflection and parallel angles for each point
        for (size_t i = 0; i != fine_tracking_path.size(); i++) {
            double angle1 = 0;  // reflection angle
            double angle2 = 0;  // parallel angle
            
            // Calculate angles using vectors to neighboring points at different distances
            for (int j = 0; j != 6; j++) {
                WireCell::Vector v10(0, 0, 0);  // Vector from previous point
                WireCell::Vector v20(0, 0, 0);  // Vector to next point
                
                // Backward vector (from point i-j-1 to point i)
                if (i > j) {
                    v10 = WireCell::Vector(fine_tracking_path.at(i).first.x() - fine_tracking_path.at(i-j-1).first.x(),
                                        fine_tracking_path.at(i).first.y() - fine_tracking_path.at(i-j-1).first.y(),
                                        fine_tracking_path.at(i).first.z() - fine_tracking_path.at(i-j-1).first.z());
                }
                
                // Forward vector (from point i to point i+j+1)
                if (i + j + 1 < fine_tracking_path.size()) {
                    v20 = WireCell::Vector(fine_tracking_path.at(i+j+1).first.x() - fine_tracking_path.at(i).first.x(),
                                        fine_tracking_path.at(i+j+1).first.y() - fine_tracking_path.at(i).first.y(),
                                        fine_tracking_path.at(i+j+1).first.z() - fine_tracking_path.at(i).first.z());
                }
                
                if (j == 0) {
                    // For the first iteration, set initial values
                    if (v10.magnitude() > 0 && v20.magnitude() > 0) {
                        angle1 = std::acos(v10.dot(v20) / (v10.magnitude() * v20.magnitude())) / 3.1415926 * 180.0;
                    }                
                    // Calculate angles with drift direction
                    if (v10.magnitude() > 0) {
                        double angle_v10 = std::acos(v10.dot(drift_dir_abs) / v10.magnitude()) / 3.1415926 * 180.0;
                        angle2 = std::abs(angle_v10 - 90.0);
                    }
                    if (v20.magnitude() > 0) {
                        double angle_v20 = std::acos(v20.dot(drift_dir_abs) / v20.magnitude()) / 3.1415926 * 180.0;
                        angle2 = std::max(angle2, std::abs(angle_v20 - 90.0));
                    }
                } else {
                    // For subsequent iterations, take minimum values
                    if (v10.magnitude() != 0 && v20.magnitude() != 0) {
                        double temp_angle1 = std::acos(v10.dot(v20) / (v10.magnitude() * v20.magnitude())) / 3.1415926 * 180.0;
                        angle1 = std::min(temp_angle1, angle1);
                        
                        double angle_v10 = std::acos(v10.dot(drift_dir_abs) / v10.magnitude()) / 3.1415926 * 180.0;
                        double angle_v20 = std::acos(v20.dot(drift_dir_abs) / v20.magnitude()) / 3.1415926 * 180.0;
                        double temp_angle2 = std::max(std::abs(angle_v10 - 90.0), std::abs(angle_v20 - 90.0));
                        angle2 = std::min(temp_angle2, angle2);
                    }
                }
            }
            
            refl_angles.at(i) = angle1;
            para_angles.at(i) = angle2;
        }

        // Second part: Analyze charge and find breakpoints
        for (int i = 0; i != fine_tracking_path.size(); i++) {
            double min_dQ_dx = dQ.at(i) / dx.at(i);
            
            // Find minimum dQ/dx in the next 5 points
            for (size_t j = 1; j != 6; j++) {
                if (i + j < fine_tracking_path.size()) {
                    if (dQ.at(i + j) / dx.at(i + j) < min_dQ_dx) {
                        min_dQ_dx = dQ.at(i + j) / dx.at(i + j);
                    }
                }
            }
            
            // Calculate sum of reflection angles in a local window
            double sum_angles = 0;
            double nsum = 0;
            
            for (int j = -2; j != 3; j++) {
                if (i + j >= 0 && i + j < fine_tracking_path.size()) {
                    if (para_angles.at(i + j) > 10) {
                        sum_angles += pow(refl_angles.at(i + j), 2);
                        nsum++;
                    }
                }
            }
            
            if (nsum != 0) {
                sum_angles = sqrt(sum_angles / nsum);
            }
            
            // First breakpoint condition: Low charge with significant angles
            if (min_dQ_dx < 1000 && para_angles.at(i) > 10 && refl_angles.at(i) > 25) {
                std::cout << "Mid_Point_Break: " << i << " " << refl_angles.at(i) << " " 
                        << para_angles.at(i) << " " << min_dQ_dx << " " 
                        << fine_tracking_path.at(i).first.x() << " " 
                        << fine_tracking_path.at(i).first.y() << " " 
                        << fine_tracking_path.at(i).first.z() << std::endl;
                flag_crawl = true;
                save_i = i;
                break;
            }
            // Second breakpoint condition: Higher angle thresholds with geometric constraints
            else if (para_angles.at(i) > 15 && refl_angles.at(i) > 27 && sum_angles > 12.5) {
                // Calculate angle between vectors from start to point and point to end
                WireCell::Vector v10(fine_tracking_path.at(i).first.x() - fine_tracking_path.front().first.x(),
                                    fine_tracking_path.at(i).first.y() - fine_tracking_path.front().first.y(),
                                    fine_tracking_path.at(i).first.z() - fine_tracking_path.front().first.z());

                WireCell::Vector v20(fine_tracking_path.back().first.x() - fine_tracking_path.at(i).first.x(),
                                    fine_tracking_path.back().first.y() - fine_tracking_path.at(i).first.y(),
                                    fine_tracking_path.back().first.z() - fine_tracking_path.at(i).first.z());

                double angle3 = 0;
                if (v10.magnitude() > 0 && v20.magnitude() > 0) {
                    angle3 = std::acos(v10.dot(v20) / (v10.magnitude() * v20.magnitude())) / 3.1415926 * 180.0;
                }
                
                // Skip if the angle is too small (nearly straight line)
                if (angle3 < 20) continue;
                
                std::cout << "Mid_Point_Break: " << i << " " << refl_angles.at(i) << " " 
                        << para_angles.at(i) << " " << angle3 << " " << min_dQ_dx << " " 
                        << fine_tracking_path.at(i).first.x() << " " 
                        << fine_tracking_path.at(i).first.y() << " " 
                        << fine_tracking_path.at(i).first.z() << std::endl;
                flag_crawl = true;
                save_i = i;
                break;
            }
        }

        std::vector<geo_point_t> out_path_points;

        if (flag_crawl){
            // Start to Crawl
            const double step_dis = 1.0 * units::cm;
            
            // Get point clouds and coordinate arrays from cluster
            const auto& steiner_pc = cluster.get_pc("steiner_pc");
            const auto& coords = cluster.get_default_scope().coords;
            const auto& steiner_x = steiner_pc.get(coords.at(0))->elements<double>();
            const auto& steiner_y = steiner_pc.get(coords.at(1))->elements<double>();
            const auto& steiner_z = steiner_pc.get(coords.at(2))->elements<double>();
      
            
            // Initialization - get starting point from fine tracking path
            geo_point_t p(fine_tracking_path.at(save_i).first.x(), 
                        fine_tracking_path.at(save_i).first.y(), 
                        fine_tracking_path.at(save_i).first.z());
            
            // Find closest point in steiner point cloud
            auto curr_knn_results = cluster.kd_steiner_knn(1, p, "steiner_pc");
            size_t curr_index = curr_knn_results[0].first;
            geo_point_t curr_wcp(steiner_x[curr_index], steiner_y[curr_index], steiner_z[curr_index]);
            
            // Calculate previous point direction
            geo_point_t prev_p(0, 0, 0);
            int num_p = 0;
            for (size_t i = 1; i != 6; i++) {
                if (save_i >= i) {
                    prev_p.at(0) += fine_tracking_path.at(save_i - i).first.x();
                    prev_p.at(1) += fine_tracking_path.at(save_i - i).first.y();
                    prev_p.at(2) += fine_tracking_path.at(save_i - i).first.z();
                    num_p++;
                }
            }
            prev_p.at(0) /= num_p;
            prev_p.at(1) /= num_p;
            prev_p.at(2) /= num_p;
            
            // Calculate initial direction
            WireCell::Vector dir(p.at(0) - prev_p.at(0), p.at(1) - prev_p.at(1), p.at(2) - prev_p.at(2));
            dir = dir.norm();
            
            bool flag_continue = true;
            while (flag_continue) {
                flag_continue = false;
                
                for (int i = 0; i != 3; i++) {
                    // Calculate test point
                    geo_point_t test_p(curr_wcp.at(0) + dir.x() * step_dis * (i + 1),
                                       curr_wcp.at(1) + dir.y() * step_dis * (i + 1),
                                       curr_wcp.at(2) + dir.z() * step_dis * (i + 1));
                    // Try normal point cloud first
                    auto search_result = cluster.get_closest_wcpoint(test_p);
                    geo_point_t next_wcp = search_result.second;
                    WireCell::Vector dir1(next_wcp.at(0) - curr_wcp.at(0),
                                          next_wcp.at(1) - curr_wcp.at(1),
                                          next_wcp.at(2) - curr_wcp.at(2));

                    // Check angle constraint (30 degrees)
                    if (dir1.magnitude() != 0 && (std::acos(dir1.dot(dir) / dir1.magnitude()) / 3.1415926 * 180.0 < 30)) {
                        flag_continue = true;
                        curr_wcp = next_wcp;
                        dir = dir1 + dir * 5.0 * units::cm; // momentum trick
                        dir = dir.norm();
                        break;
                    }
                    
                    // Try steiner point cloud
                    auto next_knn_steiner = cluster.kd_steiner_knn(1, test_p, "steiner_pc");
                    size_t next_index = next_knn_steiner[0].first;
                    next_wcp.at(0) = steiner_x[next_index];
                    next_wcp.at(1) = steiner_y[next_index];
                    next_wcp.at(2) = steiner_z[next_index];
                    WireCell::Vector dir2(next_wcp.at(0) - curr_wcp.at(0),
                                          next_wcp.at(1) - curr_wcp.at(1),
                                          next_wcp.at(2) - curr_wcp.at(2));
                    
                    // Check angle constraint (30 degrees)
                    if (dir2.magnitude() != 0 && (std::acos(dir2.dot(dir) / dir2.magnitude()) / 3.1415926 * 180.0 < 30)) {
                        flag_continue = true;
                        curr_wcp = next_wcp;
                        dir = dir2 + dir * 5.0 * units::cm; // momentum trick
                        dir = dir.norm();
                        break;
                    }
                }
            }
            
            // Find first and last points in steiner point cloud
            geo_point_t first_p(fine_tracking_path.front().first.x(), 
                            fine_tracking_path.front().first.y(), 
                            fine_tracking_path.front().first.z());
            auto first_knn_results = cluster.kd_steiner_knn(1, first_p, "steiner_pc");
            size_t first_index = first_knn_results[0].first;
            
            geo_point_t last_p(fine_tracking_path.back().first.x(), 
                            fine_tracking_path.back().first.y(), 
                            fine_tracking_path.back().first.z());
            auto last_knn_results = cluster.kd_steiner_knn(1, last_p, "steiner_pc");
            size_t last_index = last_knn_results[0].first;
            
            // Update current point to closest steiner point
            auto curr_steiner_knn = cluster.kd_steiner_knn(1, curr_wcp, "steiner_pc");
            curr_index = curr_steiner_knn[0].first;
            
            mid_p = curr_wcp;

            std::cout << "First, Center: " << steiner_x[first_index] << " " << steiner_y[first_index] 
                    << " " << steiner_z[first_index] << " " << steiner_x[curr_index] << " " 
                    << steiner_y[curr_index] << " " << steiner_z[curr_index] << std::endl;
            
            // Calculate distance from current to last point
            double dis = std::sqrt(std::pow(steiner_x[curr_index] - steiner_x[last_index], 2) + 
                                std::pow(steiner_y[curr_index] - steiner_y[last_index], 2) + 
                                std::pow(steiner_z[curr_index] - steiner_z[last_index], 2));
            
            if (dis > 1.0 * units::cm) {
                // Find path from first to current point
                const std::vector<size_t>& path1_indices = 
                    cluster.graph_algorithms("steiner_graph").shortest_path(first_index, curr_index);
                
                // Find path from current to last point
                const std::vector<size_t>& path2_indices = 
                    cluster.graph_algorithms("steiner_graph").shortest_path(curr_index, last_index);
                
                // Combine paths, removing duplicate middle point
                // Copy first path to temporary storage
                std::vector<size_t> temp_path_indices = path1_indices;
                
                // Find overlapping portion between end of path1 and beginning of path2
                int count = 0;
                auto it1 = temp_path_indices.rbegin();  // reverse iterator for temp path
                for (auto it = path2_indices.begin(); it != path2_indices.end() && it1 != temp_path_indices.rend(); ++it, ++it1) {
                    if (*it == *it1) {
                        count++;
                    } else {
                        break;
                    }
                }
                
                // Remove from end of temp_path_indices
                for (int i = 0; i != count; i++) {
                    temp_path_indices.pop_back();
                }
                
                // Add first path (without overlapping end)
                for (size_t idx : temp_path_indices) {
                    out_path_points.emplace_back(steiner_x[idx], steiner_y[idx], steiner_z[idx]);
                }
                
                // Add second path (without overlapping beginning)
                for (size_t idx : path2_indices) {
                    out_path_points.emplace_back(steiner_x[idx], steiner_y[idx], steiner_z[idx]);
                }
                
            }
        }

        return out_path_points;
    }

    int find_first_kink(std::shared_ptr<PR::Segment> segment) const{
        // Implement your logic to find the first kink in the cluster

        auto& cluster = *segment->cluster();

        // Get FiducialUtils from the grouping
        auto fiducial_utils = cluster.grouping()->get_fiducialutils();
        if (!fiducial_utils) {
            std::cout << "TaggerCheckSTM: No FiducialUtils available in find_first_kink" << std::endl;
            return -1;
        }
        const auto transform = m_pcts->pc_transform(cluster.get_scope_transform(cluster.get_default_scope()));
        double cluster_t0 = cluster.get_flash().time();
        
        // Extract fit results from the segment
        const auto& fits = segment->fits();
        
        // Convert fit data to vectors matching the TrackFitting interface
        std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> fine_tracking_path;
        std::vector<double> dQ, dx, pu, pv, pw, pt;
        std::vector<std::pair<int,int>> paf;
        
        for (const auto& fit : fits) {
            fine_tracking_path.emplace_back(fit.point, segment);
            dQ.push_back(fit.dQ);
            dx.push_back(fit.dx);
            pu.push_back(fit.pu);
            pv.push_back(fit.pv);
            pw.push_back(fit.pw);
            pt.push_back(fit.pt);
            paf.push_back(fit.paf);
        }

        if (fine_tracking_path.empty()) {
            return -1;
        }

        // Define drift direction (X direction in detector coordinates)
        WireCell::Vector drift_dir_abs(1.0, 0.0, 0.0);
        
        // Initialize angle vectors
        std::vector<double> refl_angles(fine_tracking_path.size(), 0);
        std::vector<double> para_angles(fine_tracking_path.size(), 0);
        std::vector<double> ave_angles(fine_tracking_path.size(), 0);
        std::vector<int> max_numbers(fine_tracking_path.size(), -1);
        
        // Calculate reflection and parallel angles for each point
        for (size_t i = 0; i != fine_tracking_path.size(); i++) {
            double angle1 = 0;
            double angle2 = 0;
            
            for (int j = 0; j != 6; j++) {
                WireCell::Vector v10(0, 0, 0);
                WireCell::Vector v20(0, 0, 0);
                
                // Backward vector (from point i-j-1 to point i)
                if (i > j) {
                    v10 = WireCell::Vector(fine_tracking_path.at(i).first.x() - fine_tracking_path.at(i-j-1).first.x(),
                                        fine_tracking_path.at(i).first.y() - fine_tracking_path.at(i-j-1).first.y(),
                                        fine_tracking_path.at(i).first.z() - fine_tracking_path.at(i-j-1).first.z());
                }
                
                // Forward vector (from point i to point i+j+1)
                if (i + j + 1 < fine_tracking_path.size()) {
                    v20 = WireCell::Vector(fine_tracking_path.at(i+j+1).first.x() - fine_tracking_path.at(i).first.x(),
                                        fine_tracking_path.at(i+j+1).first.y() - fine_tracking_path.at(i).first.y(),
                                        fine_tracking_path.at(i+j+1).first.z() - fine_tracking_path.at(i).first.z());
                }
                
                if (j == 0) {
                    // For the first iteration, set initial values
                    if (v10.magnitude() > 0 && v20.magnitude() > 0) {
                        angle1 = std::acos(v10.dot(v20) / (v10.magnitude() * v20.magnitude())) / 3.1415926 * 180.0;
                    }                
                    // Calculate angles with drift direction
                    if (v10.magnitude() > 0) {
                        double angle_v10 = std::acos(v10.dot(drift_dir_abs) / v10.magnitude()) / 3.1415926 * 180.0;
                        angle2 = std::abs(angle_v10 - 90.0);
                    }
                    if (v20.magnitude() > 0) {
                        double angle_v20 = std::acos(v20.dot(drift_dir_abs) / v20.magnitude()) / 3.1415926 * 180.0;
                        angle2 = std::max(angle2, std::abs(angle_v20 - 90.0));
                    }
                } else {
                    // For subsequent iterations, take minimum values
                    if (v10.magnitude() != 0 && v20.magnitude() != 0) {
                        double temp_angle1 = std::acos(v10.dot(v20) / (v10.magnitude() * v20.magnitude())) / 3.1415926 * 180.0;
                        angle1 = std::min(temp_angle1, angle1);
                        
                        double angle_v10 = std::acos(v10.dot(drift_dir_abs) / v10.magnitude()) / 3.1415926 * 180.0;
                        double angle_v20 = std::acos(v20.dot(drift_dir_abs) / v20.magnitude()) / 3.1415926 * 180.0;
                        double temp_angle2 = std::max(std::abs(angle_v10 - 90.0), std::abs(angle_v20 - 90.0));
                        angle2 = std::min(temp_angle2, angle2);
                    }
                }
            }
            
            refl_angles.at(i) = angle1;
            para_angles.at(i) = angle2;
        }
        
        // Calculate average angles in a 5-point window
        for (int i = 0; i != fine_tracking_path.size(); i++) {
            double sum_angles = 0;
            double nsum = 0;
            double max_angle = 0;
            int max_num = -1;
            
            for (int j = -2; j != 3; j++) {
                if (i + j >= 0 && i + j < fine_tracking_path.size()) {
                    if (para_angles.at(i + j) > 12) {
                        sum_angles += pow(refl_angles.at(i + j), 2);
                        nsum++;
                        if (refl_angles.at(i + j) > max_angle) {
                            max_angle = refl_angles.at(i + j);
                            max_num = i + j;
                        }
                    }
                }
            }
            
            if (nsum != 0) sum_angles = sqrt(sum_angles / nsum);
            ave_angles.at(i) = sum_angles;
            max_numbers.at(i) = max_num;
        }
        
        // Look for kink candidates
        for (int i = 0; i != fine_tracking_path.size(); i++) {
            geo_point_t current_point(fine_tracking_path.at(i).first.x(),
                                    fine_tracking_path.at(i).first.y(),
                                    fine_tracking_path.at(i).first.z());
            
            // Check basic angle conditions and fiducial volume
            if ((refl_angles.at(i) > 20 && ave_angles.at(i) > 10) && 
                fiducial_utils->inside_fiducial_volume(current_point)) {
                
                // Calculate angle between start-to-kink and kink-to-end vectors
                WireCell::Vector v10(fine_tracking_path.at(i).first.x() - fine_tracking_path.front().first.x(),
                                    fine_tracking_path.at(i).first.y() - fine_tracking_path.front().first.y(),
                                    fine_tracking_path.at(i).first.z() - fine_tracking_path.front().first.z());
                WireCell::Vector v20(fine_tracking_path.back().first.x() - fine_tracking_path.at(i).first.x(),
                                    fine_tracking_path.back().first.y() - fine_tracking_path.at(i).first.y(),
                                    fine_tracking_path.back().first.z() - fine_tracking_path.at(i).first.z());
                
                double angle3 = 0;
                if (v10.magnitude() > 0 && v20.magnitude() > 0) {
                    angle3 = std::acos(v10.dot(v20) / (v10.magnitude() * v20.magnitude())) / 3.1415926 * 180.0;
                }
                
                double angle3p = angle3;
                if (i + 1 != fine_tracking_path.size()) {
                    WireCell::Vector v11(fine_tracking_path.at(i+1).first.x() - fine_tracking_path.front().first.x(),
                                        fine_tracking_path.at(i+1).first.y() - fine_tracking_path.front().first.y(),
                                        fine_tracking_path.at(i+1).first.z() - fine_tracking_path.front().first.z());
                    WireCell::Vector v21(fine_tracking_path.back().first.x() - fine_tracking_path.at(i+1).first.x(),
                                        fine_tracking_path.back().first.y() - fine_tracking_path.at(i+1).first.y(),
                                        fine_tracking_path.back().first.z() - fine_tracking_path.at(i+1).first.z());
                    if (v11.magnitude() > 0 && v21.magnitude() > 0) {
                        angle3p = std::acos(v11.dot(v21) / (v11.magnitude() * v21.magnitude())) / 3.1415926 * 180.0;
                    }
                }
                
                // need to calculate current_point_raw ...
                WireCell::Point current_point_raw= transform->backward(current_point, cluster_t0, paf.at(i).second, paf.at(i).first);

                // Apply selection criteria
                if ((angle3 < 20 && ave_angles.at(i) < 20) || 
                    (angle3 < 12.5 && fiducial_utils->inside_dead_region(current_point_raw, paf.at(i).first, paf.at(i).second, 2)) || 
                    angle3 < 7.5 || i <= 4) continue;
                
                if ((angle3 > 30 && (refl_angles.at(i) > 25.5 && ave_angles.at(i) > 12.5)) ||
                    (angle3 > 40 && angle3 > angle3p && v10.magnitude() > 5*units::cm && v20.magnitude() > 5*units::cm)) {
                    
        //             // Special handling for shortened Y region
        //             if (pw.at(i) > 7135-5 && pw.at(i) < 7264+5) {
        //                 bool flag_bad = false;
        //                 // Check for dead channels around this position
        //                 // This would need access to ch_mcell_set_map equivalent in toolkit
        //                 // For now, apply stricter angle cuts in this region
        //                 if (pw.at(i) > 7135 && pw.at(i) < 7264) {
        //                     if (refl_angles.at(i) < 27 || ave_angles.at(i) < 15) continue;
        //                 }
        //             }
                    
                    // Calculate charge density before and after kink
                    double sum_fQ = 0;
                    double sum_fx = 0;
                    double sum_bQ = 0;
                    double sum_bx = 0;
                    
                    for (int k = 0; k != 10; k++) {
                        if (i >= k + 1) {
                            sum_fQ += dQ.at(i - k - 1);
                            sum_fx += dx.at(i - k - 1);
                        }
                        if (i + k + 1 < dQ.size()) {
                            sum_bQ += dQ.at(i + k + 1);
                            sum_bx += dx.at(i + k + 1);
                        }
                    }
                    
                    sum_fQ /= (sum_fx / units::cm + 1e-9) * 50e3;
                    sum_bQ /= (sum_bx / units::cm + 1e-9) * 50e3;
                    
                    // Final selection criteria
                    if ((sum_fQ > 0.6 && sum_bQ > 0.6) || 
                        (sum_fQ + sum_bQ > 1.4 && (sum_fQ > 0.8 || sum_bQ > 0.8) && 
                        v10.magnitude() > 10*units::cm && v20.magnitude() > 10*units::cm)) {
                        
                        if (i + 2 < dQ.size()) {
                            std::cout << "Kink: " << i << " " << refl_angles.at(i) << " " << para_angles.at(i) 
                                    << " " << ave_angles.at(i) << " " << max_numbers.at(i) << " " << angle3 
                                    << " " << dQ.at(i)/dx.at(i)*units::cm/50e3 << " " << pu.at(i) 
                                    << " " << pv.at(i) << " " << pw.at(i) << std::endl;
                            return max_numbers.at(i);
                        }
                    }
                }
            }
        }

        for (int i=0;i!=fine_tracking_path.size();i++){
            // std::cout << i << " " << refl_angles.at(i) << " " << ave_angles.at(i) << " " << inside_fiducial_volume(fine_tracking_path.at(i)) << std::endl;
            
            geo_point_t current_point(fine_tracking_path.at(i).first.x(),
                                    fine_tracking_path.at(i).first.y(),
                                    fine_tracking_path.at(i).first.z());
            
            if ((refl_angles.at(i) > 20 && ave_angles.at(i) > 15) && 
                fiducial_utils->inside_fiducial_volume(current_point)) {

                WireCell::Vector v10(fine_tracking_path.at(i).first.x() - fine_tracking_path.front().first.x(),
                                fine_tracking_path.at(i).first.y() - fine_tracking_path.front().first.y(),
                                fine_tracking_path.at(i).first.z() - fine_tracking_path.front().first.z());
                WireCell::Vector v20(fine_tracking_path.back().first.x() - fine_tracking_path.at(i).first.x(),
                                fine_tracking_path.back().first.y() - fine_tracking_path.at(i).first.y(),
                                fine_tracking_path.back().first.z() - fine_tracking_path.at(i).first.z());
                
                double angle3 = 0;
                if (v10.magnitude() > 0 && v20.magnitude() > 0) {
                    angle3 = std::acos(v10.dot(v20) / (v10.magnitude() * v20.magnitude())) / 3.1415926 * 180.0;
                }
                
                // Convert to raw coordinates for dead region check
                WireCell::Point current_point_raw = transform->backward(current_point, cluster_t0, paf.at(i).second, paf.at(i).first);
                
                if ((angle3 < 20 && ave_angles.at(i) < 20) || 
                    (angle3 < 12.5 && fiducial_utils->inside_dead_region(current_point_raw, paf.at(i).first, paf.at(i).second, 2)) || 
                    angle3 < 7.5 || i <= 4) continue;
                
                if (angle3 > 30){
            //         // shorted Y ...
            //         if (pw.at(i) > 7135 && pw.at(i) < 7264){
            //             bool flag_bad = false;
            //             for (int k=-1;k!=2;k++){
            //                 if (grouping->is_wire_dead(paf.at(i).first, paf.at(i).second, 1, std::round(pv.at(i)+k), std::round(pt.at(i)))){
            //                     flag_bad = true;
            //                     break;
            //                 }
            //             }
            //             if (flag_bad) continue;
            //         }
                    bool flag_bad_u = false;
                    {
                        for (int k=-1;k!=2;k++){
                            if (cluster.grouping()->is_wire_dead(paf.at(i).first, paf.at(i).second, 0, std::round(pu.at(i)+k), std::round(pt.at(i)))){
                                flag_bad_u = true;
                                break;
                            }
                        }
                    }
                    bool flag_bad_v = false;
                    {
                        for (int k=-1;k!=2;k++){
                            if (cluster.grouping()->is_wire_dead(paf.at(i).first, paf.at(i).second, 1, std::round(pv.at(i)+k), std::round(pt.at(i)))){
                                flag_bad_v = true;
                                break;
                            }
                        }
                    }
                    bool flag_bad_w = false;
                    {
                        for (int k=-1;k!=2;k++){
                            if (cluster.grouping()->is_wire_dead(paf.at(i).first, paf.at(i).second, 2, std::round(pw.at(i)+k), std::round(pt.at(i)))){
                                flag_bad_w = true;
                                break;
                            }
                        }
                    }
                    
                    double sum_fQ = 0;
                    double sum_fx = 0;
                    double sum_bQ = 0;
                    double sum_bx = 0;
                    for (int k=0;k!=10;k++){
                        if (i>=k+1){
                            sum_fQ += dQ.at(i-k-1);
                            sum_fx += dx.at(i-k-1);
                        }
                        if (i+k+1 < dQ.size()){
                            sum_bQ += dQ.at(i+k+1);
                            sum_bx += dx.at(i+k+1);
                        }
                    }
                    sum_fQ /= (sum_fx/units::cm+1e-9)*50e3;
                    sum_bQ /= (sum_bx/units::cm+1e-9)*50e3;
                    //std::cout << sum_fQ << " " << sum_bQ << std::endl;
                    if (std::abs(sum_fQ-sum_bQ) < 0.07*(sum_fQ+sum_bQ) && (flag_bad_u||flag_bad_v||flag_bad_w)) continue;
                    
                    if (sum_fQ > 0.6 && sum_bQ > 0.6 ){
                        if (i+2<dQ.size()){
                            std::cout << "Kink: " << i << " " << refl_angles.at(i) << " " << para_angles.at(i) << " " << ave_angles.at(i) << " " << max_numbers.at(i) << " " << angle3 << " " << dQ.at(i)/dx.at(i)*units::cm/50e3 << std::endl;
                            return max_numbers.at(i);
                        }
                    }
                }
            }
        }


        return fine_tracking_path.size();  // Placeholder return value
    }

    bool detect_proton(std::shared_ptr<PR::Segment> segment, int kink_num, std::vector<std::shared_ptr<PR::Segment>>& fitted_segments) const{
        auto& cluster = *segment->cluster();
        // Get FiducialUtils from the grouping
        auto fiducial_utils = cluster.grouping()->get_fiducialutils();
        if (!fiducial_utils) {
            std::cout << "TaggerCheckSTM: No FiducialUtils available in find_first_kink" << std::endl;
            return -1;
        }
        const auto transform = m_pcts->pc_transform(cluster.get_scope_transform(cluster.get_default_scope()));
        // double cluster_t0 = cluster.get_flash().time();
        
        // Extract fit results from the segment
        const auto& fits = segment->fits();
        
        // Convert fit data to vectors matching the TrackFitting interface
        std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> fine_tracking_path;
        std::vector<double> dQ, dx;
        std::vector<std::pair<int,int>> paf;
        
        for (const auto& fit : fits) {
            fine_tracking_path.emplace_back(fit.point, segment);
            dQ.push_back(fit.dQ);
            dx.push_back(fit.dx);
            paf.push_back(fit.paf);
        }

        if (fine_tracking_path.empty()) {
            return false;
        }

        // Extract points vector for compatibility with prototype algorithm
        std::vector<WireCell::Point> pts;
        for (const auto& path_point : fine_tracking_path) {
            pts.push_back(path_point.first);
        }

        // Determine end point
        WireCell::Point end_p;
        if (kink_num == pts.size()) {
            end_p = pts.back();
        } else {
            end_p = pts.at(kink_num);
        }

        // Calculate main track direction vector
        geo_point_t p1(pts.front().x() - pts.back().x(),
                    pts.front().y() - pts.back().y(),
                    pts.front().z() - pts.back().z());

        // Check fitted segments for Michel electrons and delta rays
        for (size_t i = 1; i < fitted_segments.size(); i++) {
            const auto& seg_fits = fitted_segments[i]->fits();
            if (seg_fits.empty()) continue;

            double dis = sqrt(pow(end_p.x() - seg_fits.front().point.x(), 2) +
                            pow(end_p.y() - seg_fits.front().point.y(), 2) +
                            pow(end_p.z() - seg_fits.front().point.z(), 2));

            // Protection against Michel electron
            double seg_length = segment_track_length(fitted_segments[i]);
            double seg_dQ_dx = segment_median_dQ_dx(fitted_segments[i]) * units::cm / 50000;
            
            if (dis < 1*units::cm && seg_length > 4*units::cm && seg_dQ_dx > 0.5) {
                return false;
            }

            // Check for delta rays
            if (dis > 10*units::cm && seg_length > 4*units::cm) {
                geo_point_t p2(seg_fits.back().point.x() - seg_fits.front().point.x(),
                        seg_fits.back().point.y() - seg_fits.front().point.y(),
                        seg_fits.back().point.z() - seg_fits.front().point.z()); 
                
                // Get closest point on the segment to point p
                const auto& fit_seg_dpc = segment->dpcloud("main");
                auto closest_result_back = fit_seg_dpc->kd3d().knn(1, seg_fits.back().point);
                size_t closest_index_back = closest_result_back[0].first;

                auto closest_result_front = fit_seg_dpc->kd3d().knn(1, seg_fits.front().point);
                size_t closest_index_front = closest_result_front[0].first;
                
                // Access the actual 3D points at the found indices
                const auto& dpc_points = fit_seg_dpc->get_points();
                const auto& closest_point_back = dpc_points[closest_index_back];
                const auto& closest_point_front = dpc_points[closest_index_front];
                
                // Access 3D coordinates
                geo_point_t back_point(closest_point_back.x, closest_point_back.y, closest_point_back.z);
                geo_point_t front_point(closest_point_front.x, closest_point_front.y, closest_point_front.z);
                
                geo_point_t p3(0,0,0);
                if (pow(end_p.x() - back_point.x(), 2) + pow(end_p.y() - back_point.y(), 2) + pow(end_p.z() - back_point.z(), 2) < 
                    pow(end_p.x() - front_point.x(), 2) + pow(end_p.y() - front_point.y(), 2) + pow(end_p.z() - front_point.z(), 2)) {
                    p3.set(front_point.x() - back_point.x(), front_point.y() - back_point.y(), front_point.z() - back_point.z());
                } else {
                    p3.set(back_point.x() - front_point.x(), back_point.y() - front_point.y(), back_point.z() - front_point.z());
                }

                // Judge direction for delta ray detection
                if ((p2.angle(p1)/3.1415926*180. < 20 || 
                    (p3.angle(p1)/3.1415926*180. < 15 && p3.magnitude() > 4*units::cm && p2.angle(p1)/3.1415926*180. < 35)) &&
                    seg_dQ_dx > 0.8) {
                    std::cout << "Delta Ray Dir: " << p2.angle(p1)/3.1415926*180. << " " 
                            << p3.angle(p1)/3.1415926*180. << " " << p3.magnitude()/units::cm << " " 
                            << dis/units::cm << " " << seg_length/units::cm << std::endl;
                    return true;
                }
            }
        }

        // Calculate cumulative distances and dQ/dx along track
        std::vector<double> L(pts.size(), 0);
        std::vector<double> dQ_dx(pts.size(), 0);
        double dis = 0;
        L[0] = dis;
        dQ_dx[0] = dQ[0] / (dx[0] / units::cm + 1e-9);
        
        for (size_t i = 1; i != pts.size(); i++) {
            dis += sqrt(pow(pts[i].x() - pts[i-1].x(), 2) + pow(pts[i].y() - pts[i-1].y(), 2) + pow(pts[i].z() - pts[i-1].z(), 2));
            L[i] = dis;
            dQ_dx[i] = dQ[i] / (dx[i] / units::cm + 1e-9);
        }

        double end_L;
        double max_num;
        if (kink_num == pts.size()) {
            end_L = L.back();
            max_num = L.size();
        } else {
            end_L = L[kink_num] - 0.5*units::cm;
            max_num = kink_num;
        }

        // Find the maximum bin
        double max_bin = -1;
        double max_sum = 0;
        for (size_t i = 0; i != L.size(); i++) {
            double sum = 0;
            double nsum = 0;
            double temp_max_bin = i;
            double temp_max_val = dQ_dx[i];
            
            if (L[i] < end_L + 0.5*units::cm && L[i] > end_L - 40*units::cm && i < max_num) {
                sum += dQ_dx[i]; nsum++;
                if (i >= 2) {
                    sum += dQ_dx[i-2]; nsum++;
                    if (dQ_dx[i-2] > temp_max_val && i-2 < max_num) {
                        temp_max_val = dQ_dx[i-2];
                        temp_max_bin = i-2;
                    }
                }
                if (i >= 1) {
                    sum += dQ_dx[i-1]; nsum++;
                    if (dQ_dx[i-1] > temp_max_val && i-1 < max_num) {
                        temp_max_val = dQ_dx[i-1];
                        temp_max_bin = i-1;
                    }
                }
                if (i+1 < L.size()) {
                    sum += dQ_dx[i+1]; nsum++;
                    if (dQ_dx[i+1] > temp_max_val && i+1 < max_num) {
                        temp_max_val = dQ_dx[i+1];
                        temp_max_bin = i+1;
                    }
                }
                if (i+2 < L.size()) {
                    sum += dQ_dx[i+2]; nsum++;
                    if (dQ_dx[i+2] > temp_max_val && i+2 < max_num) {
                        temp_max_val = dQ_dx[i+2];
                        temp_max_bin = i+2;
                    }
                }
                sum /= nsum;
                if (sum > max_sum) {
                    max_sum = sum;
                    max_bin = temp_max_bin;
                }
            }
        }

        end_L = L[max_bin] + 0.2*units::cm;
        int ncount = 0, ncount_p = 0;
        std::vector<double> vec_x, vec_xp;
        std::vector<double> vec_y, vec_yp;

        for (size_t i = 0; i != L.size(); i++) {
            if (end_L - L[i] < 35*units::cm && end_L - L[i] > 3*units::cm) {
                vec_x.push_back(end_L - L[i]);
                vec_y.push_back(dQ_dx[i]);
                ncount++;
            }

            if (end_L - L[i] < 20*units::cm) {
                vec_xp.push_back(end_L - L[i]);
                vec_yp.push_back(dQ_dx[i]);
                ncount_p++;
            }
        }

        if (ncount >= 5) {
            // Create reference vectors for comparison
            std::vector<double> muon_ref(ncount);
            std::vector<double> const_ref(ncount, 50e3);
            std::vector<double> muon_ref_p(ncount_p);

            for (size_t i = 0; i != ncount; i++) {
                muon_ref[i] = m_linterp_function->scalar_function((vec_x[i])/units::cm);
            }
            for (size_t i = 0; i != ncount_p; i++) {
                muon_ref_p[i] = m_linterp_function->scalar_function((vec_xp[i])/units::cm);
            }

            // Perform KS-like tests using kslike_compare
            double ks1 = WireCell::kslike_compare(vec_y, muon_ref);
            double ratio1 = std::accumulate(muon_ref.begin(), muon_ref.end(), 0.0) / 
                        (std::accumulate(vec_y.begin(), vec_y.end(), 0.0) + 1e-9);
            double ks2 = WireCell::kslike_compare(vec_y, const_ref);
            double ratio2 = std::accumulate(const_ref.begin(), const_ref.end(), 0.0) / 
                        (std::accumulate(vec_y.begin(), vec_y.end(), 0.0) + 1e-9);
            double ks3 = WireCell::kslike_compare(vec_yp, muon_ref_p);
            double ratio3 = std::accumulate(vec_yp.begin(), vec_yp.end(), 0.0) / 
                        (std::accumulate(muon_ref_p.begin(), muon_ref_p.end(), 0.0) + 1e-9);

            std::cout << "End proton detection: " << ks1 << " " << ks2 << " " << ratio1 << " " << ratio2 
                    << " " << ks3 << " " << ratio3 << " " << ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 
                    << " " << dQ_dx[max_bin]/50e3 << " " << dQ_dx.size() - max_bin << " " << std::endl;

            if (ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 > 0.02 && dQ_dx[max_bin]/50e3 > 2.3 &&
                (dQ_dx.size() - max_bin <= 3 || (ks2 < 0.05 && dQ_dx.size() - max_bin <= 12))) {
                if (dQ_dx.size()-max_bin <= 1 && dQ_dx[max_bin]/50e3 > 2.5 && ks2 < 0.035 && fabs(ratio2-1) < 0.1)
                return true;
                if (dQ_dx.size()-max_bin <= 1 && ((dQ_dx[max_bin]/50e3 < 3.0 && ((ks1 < 0.06 && ks2 > 0.03) || (ks1 < 0.065 && ks2 > 0.04))) || (ks1 < 0.035 && dQ_dx[max_bin]/50e3 < 4.0)))
                return false;
                if (ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 > 0.027)
                return true;
            }

            // Check for proton with very high dQ_dx
            double track_medium_dQ_dx = segment_median_dQ_dx(fitted_segments[0]) * units::cm / 50000.;
            std::cout << "End proton detection1: " << track_medium_dQ_dx << " " << dQ_dx[max_bin]/50e3 
                    << " " << ks3 << " " << ratio3 << std::endl;
                    
            if (track_medium_dQ_dx < 1.0 && dQ_dx.at(max_bin)/50e3 > 3.5){
                if ((ks3 > 0.06 && ratio3 > 1.1 && ks1 > 0.045) || (ks3 > 0.1 && ks2 < 0.19) || (ratio3 > 1.3)) return true;
                if ((ks2 < 0.045 && ks3 > 0.03) || (dQ_dx.at(max_bin)/50e3 > 4.3 && ks3 > 0.03)) return true;
            }else if (track_medium_dQ_dx < 1 && dQ_dx.at(max_bin)/50e3 > 3.0){
                if (ks3 > 0.12 && ks1 > 0.03) return true;
            }
        }

        return false;
    }

    bool eval_stm(std::shared_ptr<PR::Segment> segment, int kink_num, double peak_range = 40*units::cm, double offset_length = 0*units::cm, double com_range = 35*units::cm, bool flag_strong_check = false) const{
        auto& cluster = *segment->cluster();
        // Get FiducialUtils from the grouping
        auto fiducial_utils = cluster.grouping()->get_fiducialutils();
        if (!fiducial_utils) {
            std::cout << "TaggerCheckSTM: No FiducialUtils available in find_first_kink" << std::endl;
            return -1;
        }
        const auto transform = m_pcts->pc_transform(cluster.get_scope_transform(cluster.get_default_scope()));
        // double cluster_t0 = cluster.get_flash().time();
        
        // Extract fit results from the segment
        const auto& fits = segment->fits();
        
        // Convert fit data to vectors matching the TrackFitting interface
        std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> fine_tracking_path;
        std::vector<double> dQ, dx;
        std::vector<std::pair<int,int>> paf;
        
        for (const auto& fit : fits) {
            fine_tracking_path.emplace_back(fit.point, segment);
            dQ.push_back(fit.dQ);
            dx.push_back(fit.dx);
            paf.push_back(fit.paf);
        }

        if (fine_tracking_path.empty()) {
            return false;
        }

        // Extract points vector for compatibility with prototype algorithm
        std::vector<WireCell::Point> pts;
        for (const auto& path_point : fine_tracking_path) {
            pts.push_back(path_point.first);
        }
        
        std::vector<double> L(pts.size(), 0);
        std::vector<double> dQ_dx(pts.size(), 0);
        double dis = 0;
        L[0] = dis;
        dQ_dx[0] = dQ[0] / (dx[0] / units::cm + 1e-9);

        for (size_t i = 1; i != pts.size(); i++) {
            dis += sqrt(pow(pts[i].x() - pts[i-1].x(), 2) + 
                        pow(pts[i].y() - pts[i-1].y(), 2) + 
                        pow(pts[i].z() - pts[i-1].z(), 2));
            L[i] = dis;
            dQ_dx[i] = dQ[i] / (dx[i] / units::cm + 1e-9);
        }

        double end_L;
        double max_num;
        if (kink_num == pts.size()) {
            end_L = L.back();
            max_num = L.size();
        } else {
            end_L = L[kink_num] - 0.5 * units::cm;
            max_num = kink_num;
        }

        double max_bin = -1;
        double max_sum = 0;
        for (size_t i = 0; i != L.size(); i++) {
            double sum = 0;
            double nsum = 0;
            double temp_max_bin = i;
            double temp_max_val = dQ_dx[i];
            
            if (L[i] < end_L + 0.5 * units::cm && L[i] > end_L - peak_range && i < max_num) {
                sum += dQ_dx[i]; nsum++;
                if (i >= 2) {
                    sum += dQ_dx[i-2]; nsum++;
                    if (dQ_dx[i-2] > temp_max_val && i-2 < max_num) {
                        temp_max_val = dQ_dx[i-2];
                        temp_max_bin = i-2;
                    }
                }
                if (i >= 1) {
                    sum += dQ_dx[i-1]; nsum++;
                    if (dQ_dx[i-1] > temp_max_val && i-1 < max_num) {
                        temp_max_val = dQ_dx[i-1];
                        temp_max_bin = i-1;
                    }
                }
                if (i+1 < L.size()) {
                    sum += dQ_dx[i+1]; nsum++;
                    if (dQ_dx[i+1] > temp_max_val && i+1 < max_num) {
                        temp_max_val = dQ_dx[i+1];
                        temp_max_bin = i+1;
                    }
                }
                if (i+2 < L.size()) {
                    sum += dQ_dx[i+2]; nsum++;
                    if (dQ_dx[i+2] > temp_max_val && i+2 < max_num) {
                        temp_max_val = dQ_dx[i+2];
                        temp_max_bin = i+2;
                    }
                }
                sum /= nsum;
                if (sum > max_sum) {
                    max_sum = sum;
                    max_bin = temp_max_bin;
                }
            }
        }

        if (max_bin == -1)
            max_bin = max_num;

        end_L = L[max_bin] + 0.2 * units::cm;
        int ncount = 0;
        std::vector<double> vec_x;
        std::vector<double> vec_y;
        std::vector<double> vec_res_x;
        std::vector<double> vec_res_y;

        for (size_t i = 0; i != L.size(); i++) {
            if (end_L - L[i] < com_range && end_L - L[i] > 0) {
                vec_x.push_back(end_L - L[i]);
                vec_y.push_back(dQ_dx[i]);
                ncount++;
            } else if (L[i] > end_L) {
                vec_res_x.push_back(L[i] - end_L);
                vec_res_y.push_back(dQ_dx[i]);
            }
        }

        double ave_res_dQ_dx = 0;
        double res_length = 0;
        for (size_t i = 0; i != vec_res_y.size(); i++) {
            ave_res_dQ_dx += vec_res_y[i];
        }

        if (vec_res_y.size() > 0) {
            res_length = vec_res_x.back();
            ave_res_dQ_dx /= 1. * vec_res_y.size();
        }

        double res_length1 = 0, res_dis1 = 0;
        if (max_bin + 3 < L.size()) {
            res_length1 = L.back() - L[max_bin + 3];
            res_dis1 = sqrt(pow(pts.back().x() - pts[max_bin + 3].x(), 2) +
                        pow(pts.back().y() - pts[max_bin + 3].y(), 2) +
                        pow(pts.back().z() - pts[max_bin + 3].z(), 2));
        }

        // Create vectors for KS test instead of histograms
        std::vector<double> test_data(ncount);
        std::vector<double> ref_muon(ncount);
        std::vector<double> ref_flat(ncount);

        for (size_t i = 0; i != ncount; i++) {
            test_data[i] = vec_y[i];
            ref_muon[i] = m_linterp_function->scalar_function((vec_x[i] + offset_length) / units::cm);
            ref_flat[i] = 50e3;
        }

        double ks1 = WireCell::kslike_compare(test_data, ref_muon);
        double ratio1 = std::accumulate(ref_muon.begin(), ref_muon.end(), 0.0) / 
                        (std::accumulate(test_data.begin(), test_data.end(), 0.0) + 1e-9);
        double ks2 = WireCell::kslike_compare(test_data, ref_flat);
        double ratio2 = std::accumulate(ref_flat.begin(), ref_flat.end(), 0.0) / 
                        (std::accumulate(test_data.begin(), test_data.end(), 0.0) + 1e-9);

        std::cout << "KS value: " << flag_strong_check << " " << ks1 << " " << ks2 << " " << ratio1 << " " << ratio2 << " " << ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 << " "  << res_dis1/(res_length1+1e-9) << " " << res_length /units::cm << " " << ave_res_dQ_dx/50000. << std::endl;

        if (ks1 - ks2 >= 0.0) return false;
        if (sqrt(pow(ks2/0.06, 2) + pow((ratio2-1)/0.06, 2)) < 1.4 && 
            ks1 - ks2 + (fabs(ratio1-1) - fabs(ratio2-1))/1.5*0.3 > -0.02) return false;

        if (((res_length > 8*units::cm && ave_res_dQ_dx/50000. > 0.9 && res_length1 > 5*units::cm) ||
            (res_length1 > 1.5*units::cm && ave_res_dQ_dx/50000. > 2.3)) && res_dis1/(res_length1+1e-9) > 0.99)
            return false;

        // If residual does not look like a michel electron
        if ((res_length > 20 * units::cm && ave_res_dQ_dx/50000. > 1.2 && 
            ks1 - ks2 + (fabs(ratio1-1) - fabs(ratio2-1))/1.5*0.3 > -0.02) ||
            (res_length > 16 * units::cm && ave_res_dQ_dx > 72500) || 
            (res_length > 10 * units::cm && ave_res_dQ_dx > 72500 && 
            ks1 - ks2 + (fabs(ratio1-1) - fabs(ratio2-1))/1.5*0.3 > -0.05) ||
            (res_length > 10 * units::cm && ave_res_dQ_dx > 85000) ||
            (res_length > 6 * units::cm && ave_res_dQ_dx > 92500) ||
            (res_length > 6 * units::cm && ave_res_dQ_dx > 72500 && 
            ks1 - ks2 + (fabs(ratio1-1) - fabs(ratio2-1))/1.5*0.3 > -0.05) ||
            (res_length > 4 * units::cm && ave_res_dQ_dx/50000. > 1.4 && 
            ks1 - ks2 + (fabs(ratio1-1) - fabs(ratio2-1))/1.5*0.3 > 0.02) ||
            (res_length > 2*units::cm && ave_res_dQ_dx/50000. > 4.5))
            return false;

        if (!flag_strong_check) {
            if (ks1 - ks2 < -0.02 && ((ks2 > 0.09 && fabs(ratio2-1) > 0.1) || ratio2 > 1.5 || ks2 > 0.2)) 
                return true;
            if (ks1 - ks2 + (fabs(ratio1-1) - fabs(ratio2-1))/1.5*0.3 < 0) 
                return true;
        } else {
            if (ks1 - ks2 < -0.02 && (ks2 > 0.09 || ratio2 > 1.5) && ks1 < 0.05 && fabs(ratio1-1) < 0.1) 
                return true;
            if (ks1 - ks2 + (fabs(ratio1-1) - fabs(ratio2-1))/1.5*0.3 < 0 && ks1 < 0.05 && fabs(ratio1-1) < 0.1) 
                return true;
        }


        return false;
    }

    std::shared_ptr<PR::Segment> create_segment_for_cluster(WireCell::Clus::Facade::Cluster& cluster, 
                               const std::vector<geo_point_t>& path_points) const{
    
        // Step 3: Prepare segment data
        std::vector<PR::WCPoint> wcpoints;
        // const auto transform = m_pcts->pc_transform(cluster.get_scope_transform(cluster.get_default_scope()));
        // double cluster_t0 = cluster.get_flash().time();
        // Step 4: Create segment connecting the vertices
        auto segment = PR::make_segment();
        
        // create and associate Dynamic Point Cloud
        for (const auto& point : path_points) {
            PR::WCPoint wcp;
            wcp.point = point; 
            wcpoints.push_back(wcp);
        }

        // Step 5: Configure the segment
        segment->wcpts(wcpoints).cluster(&cluster).dirsign(1); // direction: +1, 0, or -1
               
        // auto& wcpts = segment->wcpts();
        // for (size_t i=0;i!=path_points.size(); i++){
        //     std::cout << "A: " << i << " " << path_points.at(i) << " " << wcpts.at(i).point << std::endl;
        // }
        create_segment_point_cloud(segment, path_points, m_dv, "main");

        return segment;
    }

    void search_other_tracks(Cluster& cluster, std::vector<std::shared_ptr<PR::Segment>>& fitted_segments, double search_range = 1.5*units::cm, double scaling_2d = 0.8) const{
        std::vector<std::shared_ptr<PR::Segment> > other_segments;

        // Early return if no existing segment
        if (fitted_segments.empty()) return;

        const auto& steiner_pc = cluster.get_pc("steiner_pc");
        const auto& coords = cluster.get_default_scope().coords;
        const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
        const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
        const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
        const auto& wpid_array = steiner_pc.get("wpid")->elements<WirePlaneId>();

        const size_t N = x_coords.size();
        if (N == 0) return;
        
        std::vector<bool> flag_tagged(N, false);
        int num_tagged = 0;

        const auto transform = m_pcts->pc_transform(cluster.get_scope_transform(cluster.get_default_scope()));
        double cluster_t0 = cluster.get_flash().time();
        
        // Step 1: Tag points within search_range of existing tracks        
        for (size_t i = 0; i < N; i++) {
            geo_point_t p(x_coords[i], y_coords[i], z_coords[i]);
            double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;
            
            // Get closest 2D distances for each plane using DynamicPointCloud
            // Get wire plane parameters for 2D projections
            WirePlaneId wpid = wpid_array[i];
            int apa = wpid.apa();
            int face = wpid.face();

            for (const auto& fit_seg : fitted_segments) {
                // Get closest point on the segment to point p
                const auto& fit_seg_dpc = fit_seg->dpcloud("main");

                auto closest_result = fit_seg_dpc->kd3d().knn(1, p);
                if (closest_result.empty()) continue;
                
                // size_t closest_index = closest_result[0].first;
                double closest_3d_distance = sqrt(closest_result[0].second);
                
                if (closest_3d_distance < search_range) {
                    flag_tagged[i] = true;
                    num_tagged++;
                    break;
                }    
    
                // Check distances in each plane (U, V, W)
                for (int plane = 0; plane < 3; plane++) {
                    auto closest_2d = fit_seg_dpc->get_closest_2d_point_info(p, plane, face, apa);
                    double dist_2d = std::get<0>(closest_2d);
                    
                    if (plane == 0 && dist_2d < min_dis_u) min_dis_u = dist_2d;
                    else if (plane == 1 && dist_2d < min_dis_v) min_dis_v = dist_2d;
                    else if (plane == 2 && dist_2d < min_dis_w) min_dis_w = dist_2d;
                }
            }
            
            // Additional tagging based on 2D projections and dead channels
            if (!flag_tagged[i]) {
                // Check if point should be tagged based on 2D distances or dead channels
                // figure out the raw_point ...
                auto p_raw= transform->backward(p, cluster_t0, face, apa);

                // Note: Dead channel checking would require access to detector status
                bool u_ok = (min_dis_u < scaling_2d * search_range || cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0));  // U plane
                bool v_ok = (min_dis_v < scaling_2d * search_range || cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1));  // V plane  
                bool w_ok = (min_dis_w < scaling_2d * search_range || cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2));  // W plane
                
                if (u_ok && v_ok && w_ok) {
                    flag_tagged[i] = true;
                }
            }
        }
        (void) num_tagged;

        // Step 2: Get Steiner_Graph and its terminals ...
        const auto& steiner_graph = cluster.get_graph("steiner_graph");
        const auto& flag_steiner_terminal = steiner_pc.get("flag_steiner_terminal")->elements<int>();
        
        // Step 3: Identify terminal vertices from cluster boundary points
        std::vector<size_t> terminals;
        std::map<size_t, size_t> map_oindex_tindex;
        for (size_t i = 0;i!=flag_steiner_terminal.size();i++){
            if (flag_steiner_terminal[i]){
                map_oindex_tindex[i] = terminals.size();
                terminals.push_back(i);
            }
        }

        // Step 4: Use cluster's graph for Voronoi computation
        using namespace WireCell::Clus::Graphs::Weighted;
        auto vor = voronoi(steiner_graph, terminals);
        // Access nearest terminal for any vertex like this:
        // auto nearest_terminal_for_vertex_i = vor.terminal[i];

        // Step 5: Build terminal graph and find MST
        using Base = boost::property<edge_base_t, edge_type>;  // Not edge_name_t!
        using WeightProperty = boost::property<boost::edge_weight_t, double, Base>;
        using TerminalGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                                boost::no_property, WeightProperty>;
        
        TerminalGraph terminal_graph(N);
        std::map<std::pair<size_t, size_t>, std::pair<double, edge_type>> map_saved_edge;
        
        // Build terminal graph from Voronoi regions
        auto edge_weight = get(boost::edge_weight, steiner_graph);

        for (auto w : boost::make_iterator_range(edges(steiner_graph))) {
            size_t nearest_to_source = vor.terminal[source(w, steiner_graph)];
            size_t nearest_to_target = vor.terminal[target(w, steiner_graph)];
            
            if (nearest_to_source != nearest_to_target) {
                double weight = vor.distance[source(w, steiner_graph)] + 
                            vor.distance[target(w, steiner_graph)] + 
                            edge_weight[w];  // Don't forget the edge weight!
                
                // Convert terminal indices back to actual terminal vertices
                auto edge_pair1 = std::make_pair(nearest_to_source, nearest_to_target);
                auto edge_pair2 = std::make_pair(nearest_to_target, nearest_to_source);

                auto it1 = map_saved_edge.find(edge_pair1);
                auto it2 = map_saved_edge.find(edge_pair2);

                if (it1 != map_saved_edge.end()) {
                    // Update (A,B) if better
                    if (weight < it1->second.first) {
                        it1->second = std::make_pair(weight, w);
                    }
                } else if (it2 != map_saved_edge.end()) {
                    // Update (B,A) if better  
                    if (weight < it2->second.first) {
                        it2->second = std::make_pair(weight, w);
                    }
                } else {
                    // Create new entry
                    map_saved_edge[edge_pair1] = std::make_pair(weight, w);
                }
            }
        }
        
        // Add edges with compound properties
        for (const auto& [edge_pair, weight_info] : map_saved_edge) {
            boost::add_edge(edge_pair.first, edge_pair.second, 
                        WeightProperty(weight_info.first, Base(weight_info.second)),
                        terminal_graph);
        }
        
        // Step 6: Find minimum spanning tree
        std::vector<boost::graph_traits<TerminalGraph>::edge_descriptor> mst_edges;
        boost::kruskal_minimum_spanning_tree(terminal_graph, std::back_inserter(mst_edges));
        
        // Step 7: Create cluster graph and find connected components
        TerminalGraph terminal_graph_cluster(terminals.size());
        std::map<size_t, std::set<size_t>> map_connection;
        
        for (const auto& edge : mst_edges) {
            size_t source_idx = boost::source(edge, terminal_graph);
            size_t target_idx = boost::target(edge, terminal_graph);
            
            if (flag_tagged[source_idx] == flag_tagged[target_idx]) {
                boost::add_edge(map_oindex_tindex[source_idx], map_oindex_tindex[target_idx], terminal_graph_cluster);
            } else { 
                if (map_connection.find(source_idx)==map_connection.end()){
                    std::set<size_t> temp_results;
                    temp_results.insert(target_idx);
                    map_connection[source_idx] = temp_results;
                }else{
                    map_connection[source_idx].insert(target_idx);
                }
                if (map_connection.find(target_idx)==map_connection.end()){
                    std::set<size_t> temp_results;
                    temp_results.insert(source_idx);
                    map_connection[target_idx] = temp_results;
                }else{
                    map_connection[target_idx].insert(source_idx);
                }
            }
        }
        
        // Step 8: Find connected components
        std::vector<int> component(boost::num_vertices(terminal_graph_cluster));
        const int num_components = boost::connected_components(terminal_graph_cluster, &component[0]);
        std::vector<int> ncounts(num_components, 0);
        std::vector<std::vector<size_t>> sep_clusters(num_components);
        
        for (size_t i = 0; i < component.size(); ++i) {
            ncounts[component[i]]++;
            sep_clusters[component[i]].push_back(terminals[i]);
        }
        
        // Step 9: Filter and create new segments for valid clusters
        for (int comp_idx = 0; comp_idx < num_components; comp_idx++) {
            // Skip if inside original track or just one point
            if (flag_tagged[sep_clusters[comp_idx].front()] || ncounts[comp_idx] == 1) continue;    
            // Find connection point to existing track
            size_t special_A = SIZE_MAX;
            for (size_t j = 0; j < ncounts[comp_idx]; j++) {
                if (map_connection.find(sep_clusters[comp_idx][j]) != map_connection.end()) {
                    special_A = sep_clusters[comp_idx][j];
                    break;
                }
            }
            if (special_A == SIZE_MAX) continue;
            
            // Find furthest point from special_A
            size_t special_B = special_A;
            double max_dis = 0;
            int number_not_faked = 0;
            double max_dis_u = 0, max_dis_v = 0, max_dis_w = 0;
            
            for (size_t j = 0; j < ncounts[comp_idx]; j++) {
                size_t point_idx = sep_clusters[comp_idx][j];
                geo_point_t p1(x_coords[special_A], y_coords[special_A], z_coords[special_A]);
                geo_point_t p2(x_coords[point_idx], y_coords[point_idx], z_coords[point_idx]);

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                if (dis > max_dis) {
                    max_dis = dis;
                    special_B = point_idx;
                }
                
                // Check if this track segment is "fake" (too close to existing tracks)
                double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;
                WirePlaneId wpid = wpid_array[point_idx];
                int apa = wpid.apa();
                int face = wpid.face();

                for (const auto& fit_seg : fitted_segments) {
                    const auto& fit_seg_dpc = fit_seg->dpcloud("main");
                    for (int plane = 0; plane < 3; plane++) {
                        auto closest_2d = fit_seg_dpc->get_closest_2d_point_info(p2, plane, face, apa);
                        double dist_2d = std::get<0>(closest_2d);
                        
                        if (plane == 0 && dist_2d < min_dis_u) min_dis_u = dist_2d;
                        else if (plane == 1 && dist_2d < min_dis_v) min_dis_v = dist_2d;
                        else if (plane == 2 && dist_2d < min_dis_w) min_dis_w = dist_2d;
                    }
                }
                

                auto p_raw= transform->backward(p2, cluster_t0, face, apa);

                int flag_num = 0;
                if (min_dis_u > scaling_2d * search_range && (!cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0))) flag_num++;
                if (min_dis_v > scaling_2d * search_range && (!cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1))) flag_num++;
                if (min_dis_w > scaling_2d * search_range && (!cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2))) flag_num++;

                if (min_dis_u > max_dis_u && (!cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0))) max_dis_u = min_dis_u;
                if (min_dis_v > max_dis_v && (!cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1))) max_dis_v = min_dis_v;
                if (min_dis_w > max_dis_w && (!cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2))) max_dis_w = min_dis_w;

                if (flag_num >= 2) number_not_faked++;
            }
            
            // Apply quality cuts (from prototype)
            if (number_not_faked < 4 && (number_not_faked < 0.15 * ncounts[comp_idx] || number_not_faked == 1)) continue;
            
            bool quality_check = ((max_dis_u/units::cm > 4 || max_dis_v/units::cm > 4 || max_dis_w/units::cm > 4) && 
                                max_dis_u + max_dis_v + max_dis_w > 7*units::cm) ||
                                (number_not_faked > 4 && number_not_faked >= 0.75*ncounts[comp_idx]);
            
            if (!quality_check) continue;
            
            // Step 10: Create new segment for this cluster
            std::vector<size_t> path_indices;
            
            // Use cluster's shortest path algorithm to connect special_A to special_B
            auto path_wcps = cluster.graph_algorithms("steiner_graph", m_dv, m_pcts).shortest_path(special_A, special_B);

            // Convert path to points
            std::vector<geo_point_t> path_points;
            for (size_t idx : path_wcps) {
                path_points.push_back({x_coords[idx],y_coords[idx],z_coords[idx]});
            }

            auto new_segment = create_segment_for_cluster(cluster, path_points);
        
            m_track_fitter.add_segment(new_segment);
            m_track_fitter.do_single_tracking(new_segment);

            fitted_segments.push_back(new_segment);
        }

        // std::cout << "Fitted segments: " << fitted_segments.size() << std::endl;

    }

    bool check_other_clusters(Cluster& main_cluster, std::vector<Cluster*> associated_clusters) const {
        int number_clusters = 0;
        double total_length = 0;
        
        // Iterate through all associated clusters
        for (auto it = associated_clusters.begin(); it != associated_clusters.end(); it++) {
            Cluster* cluster = *it;
            
            // Get the two boundary points (equivalent to get_two_boundary_wcps in prototype)
            std::pair<geo_point_t, geo_point_t> boundary_points = cluster->get_two_boundary_wcps();
            
            // Calculate coverage_x (difference in x coordinates)
            double coverage_x = boundary_points.first.x() - boundary_points.second.x();
            
            // Get cluster length using the toolkit's built-in function
            double length = cluster->get_length();
            
            // Get closest points between main cluster and current cluster
            std::tuple<int, int, double> results = main_cluster.get_closest_points(*cluster);
            
            // Apply the same conditions as in the prototype:
            // - Distance between clusters < 25 cm
            // - Absolute coverage in X > 0.75 cm  
            // - Length > 5 cm
            if (std::get<2>(results) < 25*units::cm && 
                fabs(coverage_x) > 0.75*units::cm && 
                length > 5*units::cm) {
                number_clusters++;
                total_length += length;
            }
        }
        
        // Apply the final condition from prototype
        if (number_clusters > 0 && 
            (number_clusters/3. + total_length/(35*units::cm)/number_clusters) >= 1) {
            std::cout << "Other clusters: " << number_clusters << " " 
                    << (number_clusters/3. + total_length/(35*units::cm)/number_clusters) 
                    << std::endl;
            return true;
        }

      return false;
    }

    bool check_other_tracks(Cluster& cluster, std::vector<std::shared_ptr<PR::Segment>>& fitted_segments) const {
        if (fitted_segments.size() <= 1) return false;
    
        int ntracks = 0;
        
        geo_point_t drift_dir_abs(1, 0, 0);
        
        // Loop through segments starting from index 1 (skip main segment)
        for (size_t i = 1; i < fitted_segments.size(); i++) {
            auto segment = fitted_segments[i];
            // Use helper functions from PRSegmentFunctions.h
            double track_length1 = segment_track_length(segment) / units::cm;
            double track_medium_dQ_dx = segment_median_dQ_dx(segment) * units::cm / 50000.;
            double track_length_threshold = segment_track_length_threshold(segment, 75000./units::cm) / units::cm;
            
            // Calculate direction vector (geometric path from front to back)
            const auto& fits = segment->fits();
            if (fits.empty()) continue;  // Skip if no fits available
            
            const auto& front_pt = fits.front().point;
            const auto& back_pt = fits.back().point;
            geo_point_t dir1(
                front_pt.x() - back_pt.x(),
                front_pt.y() - back_pt.y(), 
                front_pt.z() - back_pt.z()
            );
            
            // Calculate geometric length using helper function (maps to get_track_length(2))
            double straightness_ratio = segment_geometric_length(segment) / segment_track_length(segment);
            
            // Main logic from prototype
            if (track_length1 > 5 && track_medium_dQ_dx > 0.4) {
                ntracks++;
            }
            if (track_length1 > 40 && track_medium_dQ_dx > 0.8) return true;
            
            double angle_deg = fabs(dir1.angle(drift_dir_abs) - 3.1415926 / 2.) * 180. / 3.1415926;
            if (fabs(angle_deg - 90.0) < 7.5) continue;  // Skip tracks nearly parallel to drift
            
            // Complex condition from prototype
            if (((track_length1 > 5 && track_medium_dQ_dx > 0.7) &&
                ((track_medium_dQ_dx - 0.7)/0.1 > (19 - track_length1)/7.) &&
                straightness_ratio > 0.99) ||
                (track_length1 > 4 && track_medium_dQ_dx > 1.5 && straightness_ratio > 0.975)) {
                return true;
            }
            
            if (track_medium_dQ_dx > 1.5 && track_length1 > 8 && straightness_ratio < 0.9) continue;
            
            if ((track_medium_dQ_dx > 1.5 && track_length1 > 3) ||
                (track_medium_dQ_dx > 2.5 && track_length1 > 2.5) ||
                (track_length_threshold > 5 && ((track_length_threshold > 0.6 * track_length1) || track_length1 > 20))) {

                if (track_length1 < 5 && track_medium_dQ_dx < 2) continue;
                else if (track_length1 < 25 && track_medium_dQ_dx < 1) continue;
                else if (track_length1 < 10 && track_medium_dQ_dx < 85/50.) continue;
                else if (track_length1 < 3.5 && track_medium_dQ_dx < 110/50.) continue;
                else return true;
            }
        }
        
        if (ntracks >= 3) return true;
        return false;
    }

    /**
     * Check if a cluster meets the conditions for STM (Short Track Muon) tagging.
     * This is where you'll implement your specific STM detection algorithm.
     * 
     * @param cluster The main cluster to analyze
     * @return true if cluster should be flagged as STM
     */
    bool check_stm_conditions(Cluster& cluster, std::vector<Cluster*> associated_clusters) const {
        // get all the angles ...

        // Get all the wire plane IDs from the grouping
        const auto& wpids = cluster.grouping()->wpids();
        // Key: pair<APA, face>, Value: drift_dir, angle_u, angle_v, angle_w
        std::map<WirePlaneId , std::tuple<geo_point_t, double, double, double>> wpid_params;
        std::map<WirePlaneId, std::pair<geo_point_t, double> > wpid_U_dir;
        std::map<WirePlaneId, std::pair<geo_point_t, double> > wpid_V_dir;
        std::map<WirePlaneId, std::pair<geo_point_t, double> > wpid_W_dir;
        std::set<int> apas;
        compute_wireplane_params(wpids, m_dv, wpid_params, wpid_U_dir, wpid_V_dir, wpid_W_dir, apas);

        std::cout << "TaggerCheckSTM: Checking cluster with " << wpids.size() 
                  << " wire plane IDs and " << apas.size() << " APAs." << std::endl;

        // Early exit if no steiner graph points
        if (!cluster.has_pc("steiner_pc") || cluster.get_pc("steiner_pc").size() == 0) {
            return false;
        }

        // Get the main PCA axis direction
        const auto& pca = cluster.get_pca();
        geo_vector_t main_dir = pca.axis.at(0);

        // need to set later accoding to APA and face
        geo_point_t drift_dir, U_dir, V_dir, W_dir;

        std::vector<geo_point_t> candidate_exit_wcps;
        std::set<int> temp_set;
        std::pair<int, int> boundary_indices;

        // First round check - get boundary points from steiner graph
        boundary_indices = cluster.get_two_boundary_steiner_graph_idx("steiner_graph", "steiner_pc", true);

        // Get extreme points
        std::vector<std::vector<geo_point_t>> out_vec_wcps = cluster.get_extreme_wcps();

        // Get the steiner_pc to access actual points using boundary_indices
        const auto& steiner_pc = cluster.get_pc("steiner_pc");
        const auto& coords = cluster.get_default_scope().coords;
        const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
        const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
        const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
        
        // Add the two boundary points as additional extreme point groups
        geo_point_t boundary_point_first(x_coords[boundary_indices.first], 
                                     y_coords[boundary_indices.first], 
                                     z_coords[boundary_indices.first]);
        geo_point_t boundary_point_second(x_coords[boundary_indices.second], 
                                     y_coords[boundary_indices.second], 
                                     z_coords[boundary_indices.second]);
        {
            std::vector<geo_point_t> temp_wcps;
            temp_wcps.push_back(boundary_point_first);
            out_vec_wcps.push_back(temp_wcps);
        }
        {
            std::vector<geo_point_t> temp_wcps;          
            temp_wcps.push_back(boundary_point_second);
            out_vec_wcps.push_back(temp_wcps);
        }

        // Get FiducialUtils from the grouping
        auto fiducial_utils = cluster.grouping()->get_fiducialutils();
        if (!fiducial_utils) {
            std::cout << "TaggerCheckSTM: No FiducialUtils available" << std::endl;
            return false;
        }

        // Boundary check
        for (size_t i = 0; i != out_vec_wcps.size(); i++) {
            bool flag_save = false;
            
            // Check all the points in this extreme point group
            for (size_t j = 0; j != out_vec_wcps.at(i).size(); j++) {
                geo_point_t p1 = out_vec_wcps.at(i).at(j);
                
                if (!fiducial_utils->inside_fiducial_volume(p1)) {
                    candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
                    flag_save = true;
                    break;
                }
            }
            
            if (!flag_save) {
                // Check direction using vhough_transform
                geo_point_t p1 = out_vec_wcps.at(i).at(0);
                geo_vector_t dir_vec = cluster.vhough_transform(p1, 30*units::cm);
                dir_vec = dir_vec * (-1.0);  // Reverse direction
                
                // Convert to geo_point_t for angle calculations
                geo_point_t dir(dir_vec.x(), dir_vec.y(), dir_vec.z());
                
                // Check U, V, and W angles
                geo_point_t dir_1(0, dir.y(), dir.z());

                // get apa and face from the point p1 ... 
                auto wpid_p1 = m_dv->contained_by(p1);
                // Fill drift_dir, U_dir, V_dir, W_dir from the maps using wpid_p1
                auto it_params = wpid_params.find(wpid_p1);
                if (it_params != wpid_params.end()) {
                    drift_dir = std::get<0>(it_params->second);
                } else {
                    std::cerr << "TaggerCheckSTM: wpid_params not found for wpid_p1" << std::endl;
                }

                auto it_U = wpid_U_dir.find(wpid_p1);
                if (it_U != wpid_U_dir.end()) {
                    U_dir = it_U->second.first;
                } else {
                    std::cerr << "TaggerCheckSTM: wpid_U_dir not found for wpid_p1" << std::endl;
                }

                auto it_V = wpid_V_dir.find(wpid_p1);
                if (it_V != wpid_V_dir.end()) {
                    V_dir = it_V->second.first;
                } else {
                    std::cerr << "TaggerCheckSTM: wpid_V_dir not found for wpid_p1" << std::endl;
                }

                auto it_W = wpid_W_dir.find(wpid_p1);
                if (it_W != wpid_W_dir.end()) {
                    W_dir = it_W->second.first;
                } else {
                    std::cerr << "TaggerCheckSTM: wpid_W_dir not found for wpid_p1" << std::endl;
                }

                // std::cout << "TaggerCheckSTM: Checking angles for point " 
                //           << p1 << " with wpid " << wpid_p1 << " " << drift_dir << " " << U_dir << " " << V_dir << " " << W_dir << std::endl;
                
                // Calculate angles with wire directions
                double angle1 = acos(dir_1.dot(U_dir) / (dir_1.magnitude() * U_dir.magnitude()));
                geo_point_t tempV1(fabs(dir.x()), 
                                sqrt(dir.y()*dir.y() + dir.z()*dir.z()) * sin(angle1), 0);
                double angle1_1 = acos(tempV1.dot(drift_dir) / (tempV1.magnitude() * drift_dir.magnitude())) / 3.1415926 * 180.;
                
                double angle2 = acos(dir_1.dot(V_dir) / (dir_1.magnitude() * V_dir.magnitude()));
                geo_point_t tempV2(fabs(dir.x()), 
                                sqrt(dir.y()*dir.y() + dir.z()*dir.z()) * sin(angle2), 0);
                double angle2_1 = acos(tempV2.dot(drift_dir) / (tempV2.magnitude() * drift_dir.magnitude())) / 3.1415926 * 180.;
                
                double angle3 = acos(dir_1.dot(W_dir) / (dir_1.magnitude() * W_dir.magnitude()));
                geo_point_t tempV3(fabs(dir.x()), 
                                sqrt(dir.y()*dir.y() + dir.z()*dir.z()) * sin(angle3), 0);
                double angle3_1 = acos(tempV3.dot(drift_dir) / (tempV3.magnitude() * drift_dir.magnitude())) / 3.1415926 * 180.;
                
                if ((angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)) {
                    if (!fiducial_utils->check_signal_processing(cluster, p1, dir_vec, 1*units::cm)) {
                        flag_save = true;
                        candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
                    }
                }
                
                if (!flag_save) {
                    // Calculate angle between direction and main axis
                    double main_angle = acos(dir_vec.dot(main_dir) / (dir_vec.magnitude() * main_dir.magnitude()));
                    double angle_deg = fabs((3.1415926/2. - main_angle) / 3.1415926 * 180.);
                    
                    if (angle_deg > 60) {
                        if (!fiducial_utils->check_dead_volume(cluster, p1, dir_vec, 1*units::cm)) {
                            flag_save = true;
                            candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
                        }
                    }
                }
            }
        }

        // Determine which boundary points are exit candidates
        for (size_t i = 0; i != candidate_exit_wcps.size(); i++) {
            double dis1 = (candidate_exit_wcps.at(i) - boundary_point_first).magnitude();
            double dis2 = (candidate_exit_wcps.at(i) - boundary_point_second).magnitude();

            // Check if essentially one of the extreme points
            if (dis1 < dis2) {
                if (dis1 < 1.0*units::cm) temp_set.insert(0);
            } else {
                if (dis2 < 1.0*units::cm) temp_set.insert(1);
            }
        }

        // Protection against two end point situation
        if (temp_set.size() == 2) {
            geo_point_t tp1 = boundary_point_first;
            geo_point_t tp2 = boundary_point_second;
            
            temp_set.clear();
            
            if (!fiducial_utils->inside_fiducial_volume(tp1)) temp_set.insert(0);
            if (!fiducial_utils->inside_fiducial_volume(tp2)) temp_set.insert(1);
            if (temp_set.size() == 0) {
                temp_set.insert(0);
                temp_set.insert(1);
            }
        }

        // Second round check if no candidates found
        if (temp_set.size() == 0) {
            candidate_exit_wcps.clear();
            
            // Repeat the process with flag_cosmic = false
            boundary_indices = cluster.get_two_boundary_steiner_graph_idx("steiner_graph", "steiner_pc", false);
            out_vec_wcps = cluster.get_extreme_wcps();
            
            // Get the steiner_pc to access actual points using boundary_indices
            const auto& steiner_pc = cluster.get_pc("steiner_pc");
            const auto& coords = cluster.get_default_scope().coords;
            const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
            const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
            const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
            
            // Add the two boundary points as additional extreme point groups
            geo_point_t boundary_point_first(x_coords[boundary_indices.first], 
                                        y_coords[boundary_indices.first], 
                                        z_coords[boundary_indices.first]);
            geo_point_t boundary_point_second(x_coords[boundary_indices.second], 
                                        y_coords[boundary_indices.second], 
                                        z_coords[boundary_indices.second]);

            // Add boundary points again
            {
                std::vector<geo_point_t> temp_wcps;
                temp_wcps.push_back(boundary_point_first);
                out_vec_wcps.push_back(temp_wcps);
            }
            {
                std::vector<geo_point_t> temp_wcps;
                temp_wcps.push_back(boundary_point_second);
                out_vec_wcps.push_back(temp_wcps);
            }
            
            // Repeat boundary check (same logic as above)
            for (size_t i = 0; i != out_vec_wcps.size(); i++) {
                bool flag_save = false;
                
                for (size_t j = 0; j != out_vec_wcps.at(i).size(); j++) {
                    geo_point_t p1 = out_vec_wcps.at(i).at(j);
                    
                    if (!fiducial_utils->inside_fiducial_volume(p1)) {
                        candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
                        flag_save = true;
                        break;
                    }
                }
                
                if (!flag_save) {
                    geo_point_t p1 = out_vec_wcps.at(i).at(0);
                    geo_vector_t dir_vec = cluster.vhough_transform(p1, 30*units::cm);
                    dir_vec = dir_vec * (-1.0);
                    
                    geo_point_t dir(dir_vec.x(), dir_vec.y(), dir_vec.z());
                    geo_point_t dir_1(0, dir.y(), dir.z());
                    
                    // get apa and face from the point p1 ... 
                    auto wpid_p1 = m_dv->contained_by(p1);
                    // Fill drift_dir, U_dir, V_dir, W_dir from the maps using wpid_p1
                    auto it_params = wpid_params.find(wpid_p1);
                    if (it_params != wpid_params.end()) {
                        drift_dir = std::get<0>(it_params->second);
                    } else {
                        std::cerr << "TaggerCheckSTM: wpid_params not found for wpid_p1" << std::endl;
                    }

                    auto it_U = wpid_U_dir.find(wpid_p1);
                    if (it_U != wpid_U_dir.end()) {
                        U_dir = it_U->second.first;
                    } else {
                        std::cerr << "TaggerCheckSTM: wpid_U_dir not found for wpid_p1" << std::endl;
                    }

                    auto it_V = wpid_V_dir.find(wpid_p1);
                    if (it_V != wpid_V_dir.end()) {
                        V_dir = it_V->second.first;
                    } else {
                        std::cerr << "TaggerCheckSTM: wpid_V_dir not found for wpid_p1" << std::endl;
                    }

                    auto it_W = wpid_W_dir.find(wpid_p1);
                    if (it_W != wpid_W_dir.end()) {
                        W_dir = it_W->second.first;
                    } else {
                        std::cerr << "TaggerCheckSTM: wpid_W_dir not found for wpid_p1" << std::endl;
                    }


                    double angle1 = acos(dir_1.dot(U_dir) / (dir_1.magnitude() * U_dir.magnitude()));
                    geo_point_t tempV1(fabs(dir.x()), 
                                    sqrt(dir.y()*dir.y() + dir.z()*dir.z()) * sin(angle1), 0);
                    double angle1_1 = acos(tempV1.dot(drift_dir) / (tempV1.magnitude() * drift_dir.magnitude())) / 3.1415926 * 180.;
                    
                    double angle2 = acos(dir_1.dot(V_dir) / (dir_1.magnitude() * V_dir.magnitude()));
                    geo_point_t tempV2(fabs(dir.x()), 
                                    sqrt(dir.y()*dir.y() + dir.z()*dir.z()) * sin(angle2), 0);
                    double angle2_1 = acos(tempV2.dot(drift_dir) / (tempV2.magnitude() * drift_dir.magnitude())) / 3.1415926 * 180.;
                    
                    double angle3 = acos(dir_1.dot(W_dir) / (dir_1.magnitude() * W_dir.magnitude()));
                    geo_point_t tempV3(fabs(dir.x()), 
                                    sqrt(dir.y()*dir.y() + dir.z()*dir.z()) * sin(angle3), 0);
                    double angle3_1 = acos(tempV3.dot(drift_dir) / (tempV3.magnitude() * drift_dir.magnitude())) / 3.1415926 * 180.;
                    
                    if ((angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)) {
                        if (!fiducial_utils->check_signal_processing(cluster, p1, dir_vec, 1*units::cm)) {
                            flag_save = true;
                            candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
                        }
                    }
                    
                    if (!flag_save) {
                        double main_angle = acos(dir_vec.dot(main_dir) / (dir_vec.magnitude() * main_dir.magnitude()));
                        double angle_deg = fabs((3.1415926/2. - main_angle) / 3.1415926 * 180.);
                        
                        if (angle_deg > 60) {
                            if (!fiducial_utils->check_dead_volume(cluster, p1, dir_vec, 1*units::cm)) {
                                flag_save = true;
                                candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
                            }
                        }
                    }
                }
            }
            
            // Determine boundary points again
            for (size_t i = 0; i != candidate_exit_wcps.size(); i++) {
                double dis1 = (candidate_exit_wcps.at(i) - boundary_point_first).magnitude();
                double dis2 = (candidate_exit_wcps.at(i) - boundary_point_second).magnitude();

                if (dis1 < dis2) {
                    if (dis1 < 1.0*units::cm) temp_set.insert(0);
                } else {
                    if (dis2 < 1.0*units::cm) temp_set.insert(1);
                }
            }
            
            // Protection against two end point situation
            if (temp_set.size() == 2) {
                geo_point_t tp1 = boundary_point_first;
                geo_point_t tp2 = boundary_point_second;

                temp_set.clear();
                
                if (!fiducial_utils->inside_fiducial_volume(tp1)) temp_set.insert(0);
                if (!fiducial_utils->inside_fiducial_volume(tp2)) temp_set.insert(1);
                if (temp_set.size() == 0) {
                    temp_set.insert(0);
                    temp_set.insert(1);
                }
            }
        }

        // Fully contained, so not a STM
        if (candidate_exit_wcps.size() == 0) {
            std::cout << "STMTagger: Mid Point: A" << std::endl;
            return false;
        }

        // std::cout << "end_point: " << temp_set.size() << " " << candidate_exit_wcps.size() << std::endl;

        // Determine first and last points for further analysis
        geo_point_t first_wcp, last_wcp;
        bool flag_double_end = false;

        if (temp_set.size() != 0) {
            if (*temp_set.begin() == 0) {
                first_wcp = boundary_point_first;
                last_wcp = boundary_point_second;
            } else {
                first_wcp = boundary_point_second;
                last_wcp = boundary_point_first;
            }
            if (temp_set.size() == 2) flag_double_end = true;
        } else {
            if (candidate_exit_wcps.size() == 1) {
                first_wcp = candidate_exit_wcps.at(0);

                geo_vector_t dir1 = boundary_point_first - candidate_exit_wcps.at(0);
                geo_vector_t dir2 = boundary_point_second - candidate_exit_wcps.at(0);
                double dis1 = dir1.magnitude();
                double dis2 = dir2.magnitude();
                
                double angle_between = acos(dir1.dot(dir2) / (dis1 * dis2));
                
                if (angle_between > 120./180.*3.1415926 && dis1 > 20*units::cm && dis2 > 20*units::cm) {
                    std::cout << "Mid Point: B" << std::endl;
                    return false;
                } else {
                    if (dis1 < dis2) {
                        last_wcp = boundary_point_second;
                    } else {
                        last_wcp = boundary_point_first;
                    }
                }
            } else {
                std::cout << "Mid Point: C" << std::endl;
                return false;
            }
        }

        bool flag_other_clusters = check_other_clusters(cluster, associated_clusters);

        std::cout << "STM analysis: flag_double_end=" << flag_double_end 
                  << ", flag_other_clusters=" << flag_other_clusters << std::endl;

        // Forward check
        {
            if (flag_double_end) std::cout << "Forward check!" << std::endl;
            
            // Do rough path tracking
            auto path_points = do_rough_path(cluster, first_wcp, last_wcp);
            
            {
                // Create segment for tracking
                auto segment = create_segment_for_cluster(cluster, path_points);
                m_track_fitter.add_segment(segment);
                m_track_fitter.do_single_tracking(segment, false);
                // Extract fit results from the segment
                const auto& fits = segment->fits();
                if (fits.size() <=3) return false;
            }

            std::cout << "Finish first round of fitting" << std::endl;

            geo_point_t mid_point(0,0,0);
            auto adjusted_path_points = adjust_rough_path(cluster, mid_point);
            auto adjusted_segment = create_segment_for_cluster(cluster, adjusted_path_points);
            m_track_fitter.clear_segments();
            m_track_fitter.add_segment(adjusted_segment);
            m_track_fitter.do_single_tracking(adjusted_segment);

            std::cout << "Finish second round of fitting" << std::endl;

            std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> fine_tracking_path;
            std::vector<double> dQ, dx;
            std::vector<std::pair<int,int>> paf;

            const auto& fits = adjusted_segment->fits();
            for (const auto& fit : fits) {
                fine_tracking_path.emplace_back(fit.point, adjusted_segment);
                dQ.push_back(fit.dQ);
                dx.push_back(fit.dx);
                paf.push_back(fit.paf);
            }

            // Extract points for compatibility
            std::vector<WireCell::Point> pts;
            for (const auto& path_point : fine_tracking_path) {
                pts.push_back(path_point.first);
            }

            std::cout << "Collect points " << pts.size() << std::endl;

            int kink_num = find_first_kink(adjusted_segment);

            std::cout << "Kink Number: " << kink_num << std::endl;

            double left_L = 0; 
            double left_Q = 0;
            double exit_L = 0; 
            double exit_Q = 0;
            
            for (size_t i=0; i != kink_num && i < dx.size(); i++){
                exit_L += dx.at(i);
                exit_Q += dQ.at(i);
            }
            for (size_t i = kink_num; i < dx.size(); i++){
                left_L += dx.at(i);
                left_Q += dQ.at(i);
            }
            
            std::cout << "Left: " << exit_L/units::cm << " " << left_L/units::cm << " " 
                      << (left_Q/(left_L/units::cm+1e-9))/50e3 << " " 
                      << (exit_Q/(exit_L/units::cm+1e-9)/50e3) << std::endl;

            // TGM (Through-Going Muon) check
            if ((!fiducial_utils->inside_fiducial_volume(pts.front())) && (!fiducial_utils->inside_fiducial_volume(pts.back()))){

                bool flag_TGM_anode = false;
                if ((pts.back().x() < 2*units::cm || pts.front().x() < 2*units::cm) && 
                    kink_num >= 0 && kink_num < pts.size()) {
                    if (pts.at(kink_num).x() < 6*units::cm){
                        geo_point_t v10(pts.back().x()-pts.at(kink_num).x(),
                                       pts.back().y()-pts.at(kink_num).y(),
                                       pts.back().z()-pts.at(kink_num).z());
                        geo_point_t v20(pts.front().x()-pts.at(kink_num).x(),
                                       pts.front().y()-pts.at(kink_num).y(),
                                       pts.front().z()-pts.at(kink_num).z());
                        
                        if ((fabs(v10.angle(drift_dir)/3.1415926*180.-90)<12.5 && v10.magnitude()>15*units::cm) || 
                            (fabs(v20.angle(drift_dir)/3.1415926*180.-90)<12.5 && v20.magnitude()>15*units::cm)) {
                            flag_TGM_anode = true;
                        }
                    }
                }
                
                if ((exit_L < 3*units::cm || left_L < 3*units::cm) || flag_TGM_anode){
                    std::cout << "TGM: " << pts.front() << " " << pts.back() << std::endl;
                    cluster.set_flag(Flags::TGM);
                    return true;
                }
                
            } 
            else if ((!fiducial_utils->inside_fiducial_volume(pts.front())) && left_L < 3*units::cm){
                // Check dead volume
                WireCell::Point p1 = pts.back();
                geo_point_t dir_vec = cluster.vhough_transform(p1, 30*units::cm);
                dir_vec *= -1;
                
                if (!fiducial_utils->check_dead_volume(cluster, p1, dir_vec, 1*units::cm)){
                    if (exit_L < 3*units::cm || left_L < 3*units::cm){
                        std::cout << "TGM: " << pts.front() << " " << pts.back() << std::endl;
                        cluster.set_flag(Flags::TGM);
                        return true;
                    }
                }
            } 

            // STM evaluation logic
            if (left_L > 40*units::cm || (left_L > 7.5*units::cm && (left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0)){
                if (!flag_double_end){
                    std::cout << "Mid Point A " << " Fid "
                             << " " << mid_point << " " << left_L << " " 
                             << (left_Q/(left_L/units::cm+1e-9)/50e3) << std::endl;
                    return false;
                }
            } else {
                bool flag_fix_end = false;
                if (exit_L < 35*units::cm || ((left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0 && left_L > 2*units::cm)) {
                    flag_fix_end = true;
                }
                
                // Readjust parameters for short tracks
                if ((left_L < 8*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3)< 1.5) ||
                    (left_L < 6*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.7) ||
                    (left_L < 5*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.8) ||
                    (left_L < 3*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.9)){
                    left_L = 0;
                    kink_num = dQ.size();
                    exit_L = 40*units::cm;
                    flag_fix_end = false;
                }

                bool flag_pass = false;

                if (!flag_other_clusters){
                    if (left_L < 40*units::cm) {
                        if (flag_fix_end){
                            flag_pass = eval_stm(adjusted_segment, kink_num, 5*units::cm, 0., 35*units::cm) ||
                                       eval_stm(adjusted_segment, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
                        } else {
                            flag_pass = eval_stm(adjusted_segment, kink_num, 40*units::cm - left_L, 0., 35*units::cm) ||
                                       eval_stm(adjusted_segment, kink_num, 40*units::cm - left_L, 3.*units::cm, 35*units::cm);
                        }
                        
                        if (!flag_pass){
                            if (flag_fix_end){
                                flag_pass = eval_stm(adjusted_segment, kink_num, 5*units::cm, 0., 15*units::cm) ||
                                           eval_stm(adjusted_segment, kink_num, 5*units::cm, 3.*units::cm, 15*units::cm);
                            } else {
                                flag_pass = eval_stm(adjusted_segment, kink_num, 40*units::cm - left_L, 0., 15*units::cm) ||
                                           eval_stm(adjusted_segment, kink_num, 40*units::cm - left_L, 3.*units::cm, 15*units::cm);
                            }
                        }
                    }
                    
                    if (left_L < 20*units::cm){
                        if (!flag_pass){
                            if (flag_fix_end){
                                flag_pass = eval_stm(adjusted_segment, kink_num, 5*units::cm, 0., 35*units::cm) ||
                                           eval_stm(adjusted_segment, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
                            } else {
                                flag_pass = eval_stm(adjusted_segment, kink_num, 20*units::cm - left_L, 0., 35*units::cm) ||
                                           eval_stm(adjusted_segment, kink_num, 20*units::cm - left_L, 3.*units::cm, 35*units::cm);
                            }
                        }
                        
                        if (!flag_pass){
                            if (flag_fix_end){
                                flag_pass = eval_stm(adjusted_segment, kink_num, 5*units::cm, 0., 15*units::cm) ||
                                           eval_stm(adjusted_segment, kink_num, 5*units::cm, 3.*units::cm, 15*units::cm);
                            } else {
                                flag_pass = eval_stm(adjusted_segment, kink_num, 20*units::cm - left_L, 0., 15*units::cm) ||
                                           eval_stm(adjusted_segment, kink_num, 20*units::cm - left_L, 3.*units::cm, 15*units::cm);
                            }
                        }
                    }
                } else {
                    if (flag_fix_end) {
                        flag_pass = eval_stm(adjusted_segment, kink_num, 5*units::cm, 0., 35*units::cm, true);
                    } else {
                        flag_pass = eval_stm(adjusted_segment, kink_num, 40*units::cm, 0., 35*units::cm, true);
                    }
                }
                
                if (flag_pass) {
                    std::vector<std::shared_ptr<PR::Segment>> fitted_segments;
                    fitted_segments.push_back(adjusted_segment);
                    search_other_tracks(cluster, fitted_segments);

                    if (check_other_tracks(cluster, fitted_segments)){
                        std::cout << "Mid Point Tracks" << std::endl;
                        return false;
                    }
                    
                    if (!detect_proton(adjusted_segment, kink_num, fitted_segments)) return true;
                }
            }
        }
        
        // Backward check (if double-ended)
        if (flag_double_end){
            std::cout << "Backward check!" << std::endl;
            
            {
                m_track_fitter.clear_segments();
                // Do backward path tracking
                auto path_points = do_rough_path(cluster, last_wcp, first_wcp);
                auto segment = create_segment_for_cluster(cluster, path_points);
                m_track_fitter.add_segment(segment);
                m_track_fitter.do_single_tracking(segment, false);

            }
            geo_point_t mid_point(0,0,0);
            auto adjusted_path_points = adjust_rough_path(cluster, mid_point);
            auto adjusted_segment = create_segment_for_cluster(cluster, adjusted_path_points);
            m_track_fitter.clear_segments();
            m_track_fitter.add_segment(adjusted_segment);
            m_track_fitter.do_single_tracking(adjusted_segment);

            std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> fine_tracking_path;
            std::vector<double> dQ, dx;
            std::vector<std::pair<int,int>> paf;

            const auto& fits = adjusted_segment->fits();
            for (const auto& fit : fits) {
                fine_tracking_path.emplace_back(fit.point, adjusted_segment);
                dQ.push_back(fit.dQ);
                dx.push_back(fit.dx);
                paf.push_back(fit.paf);
            }

            // Extract points for compatibility
            std::vector<WireCell::Point> pts;
            for (const auto& path_point : fine_tracking_path) {
                pts.push_back(path_point.first); 
            }

            int kink_num = find_first_kink(adjusted_segment);
            
            double left_L = 0;
            double left_Q = 0;
            double exit_L = 0;
            double exit_Q = 0;
            
            for (size_t i=0; i != kink_num && i < dx.size(); i++){
                exit_L += dx.at(i);
                exit_Q += dQ.at(i);
            }
            for (size_t i = kink_num; i != dx.size(); i++){
                left_L += dx.at(i);
                left_Q += dQ.at(i);
            }
            
            std::cout << "Left: " << exit_L/units::cm << " " << left_L/units::cm << " " 
                      << (left_Q/(left_L/units::cm+1e-9))/50e3 << " " 
                      << (exit_Q/(exit_L/units::cm+1e-9)/50e3) << std::endl;
            
            // TGM check for backward direction
            if ((!fiducial_utils->inside_fiducial_volume(pts.front())) && 
                (!fiducial_utils->inside_fiducial_volume(pts.back()))){
                
                bool flag_TGM_anode = false;
                
                if ((pts.back().x() < 2*units::cm || pts.front().x() < 2*units::cm) && 
                    kink_num >= 0 && kink_num < pts.size()) {
                    if (pts.at(kink_num).x() < 6*units::cm){
                        geo_point_t v10(pts.back().x()-pts.at(kink_num).x(),
                                       pts.back().y()-pts.at(kink_num).y(),
                                       pts.back().z()-pts.at(kink_num).z());
                        geo_point_t v20(pts.front().x()-pts.at(kink_num).x(),
                                       pts.front().y()-pts.at(kink_num).y(),
                                       pts.front().z()-pts.at(kink_num).z());
                        
                        if ((fabs(v10.angle(drift_dir)/3.1415926*180.-90)<12.5 && v10.magnitude()>15*units::cm) || 
                            (fabs(v20.angle(drift_dir)/3.1415926*180.-90)<12.5 && v20.magnitude()>15*units::cm)) {
                            flag_TGM_anode = true;
                        }
                    }
                }
                
                if ((exit_L < 3*units::cm || left_L < 3*units::cm) || flag_TGM_anode){
                    std::cout << "TGM: " << pts.front() << " " << pts.back() << std::endl;
                    cluster.set_flag(Flags::TGM);
                    return true;
                }
            }

            if (left_L > 40*units::cm || (left_L > 7.5*units::cm && (left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0)){
                std::cout << "Mid Point A " << " Fid" 
                         << " " << mid_point << " " << left_L << " " 
                         << (left_Q/(left_L/units::cm+1e-9)/50e3) << std::endl;
                return false;
            } else {
                bool flag_fix_end = false;
                if (exit_L < 35*units::cm || ((left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0 && left_L > 2*units::cm)) {
                    flag_fix_end = true;
                }
                
                if ((left_L < 8*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3)< 1.5) ||
                    (left_L < 6*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.7) ||
                    (left_L < 3*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.9)){
                    left_L = 0;
                    kink_num = dQ.size();
                    exit_L = 40*units::cm;
                    flag_fix_end = false;
                }

                bool flag_pass = false;
                if (!flag_other_clusters){
                    if (left_L < 40*units::cm) {
                        if (flag_fix_end){
                            flag_pass = eval_stm(adjusted_segment, kink_num, 5*units::cm, 0., 35*units::cm) ||
                                       eval_stm(adjusted_segment, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
                        } else {
                            flag_pass = eval_stm(adjusted_segment, kink_num, 40*units::cm - left_L, 0., 35*units::cm) ||
                                       eval_stm(adjusted_segment, kink_num, 40*units::cm - left_L, 3.*units::cm, 35*units::cm);
                        }
                        
                        if (!flag_pass){
                            if (flag_fix_end){
                                flag_pass = eval_stm(adjusted_segment, kink_num, 5*units::cm, 0., 15*units::cm) ||
                                           eval_stm(adjusted_segment, kink_num, 5*units::cm, 3.*units::cm, 15*units::cm);
                            } else {
                                flag_pass = eval_stm(adjusted_segment, kink_num, 40*units::cm - left_L, 0., 15*units::cm) ||
                                           eval_stm(adjusted_segment, kink_num, 40*units::cm - left_L, 3.*units::cm, 15*units::cm);
                            }
                        }
                    }
                    
                    if (left_L < 20*units::cm){
                        if (!flag_pass){
                            if (flag_fix_end){
                                flag_pass = eval_stm(adjusted_segment, kink_num, 5*units::cm, 0., 35*units::cm) ||
                                           eval_stm(adjusted_segment, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
                            } else {
                                flag_pass = eval_stm(adjusted_segment, kink_num, 20*units::cm - left_L, 0., 35*units::cm) ||
                                           eval_stm(adjusted_segment, kink_num, 20*units::cm - left_L, 3.*units::cm, 35*units::cm);
                            }
                        }
                        
                        if (!flag_pass){
                            if (flag_fix_end){
                                flag_pass = eval_stm(adjusted_segment, kink_num, 5*units::cm, 0., 15*units::cm) ||
                                           eval_stm(adjusted_segment, kink_num, 5*units::cm, 3.*units::cm, 15*units::cm);
                            } else {
                                flag_pass = eval_stm(adjusted_segment, kink_num, 20*units::cm - left_L, 0., 15*units::cm) ||
                                           eval_stm(adjusted_segment, kink_num, 20*units::cm - left_L, 3.*units::cm, 15*units::cm);
                            }
                        }
                    }
                } else {
                    if (flag_fix_end) {
                        flag_pass = eval_stm(adjusted_segment, kink_num, 5*units::cm, 0., 35*units::cm, true);
                    } else {
                        flag_pass = eval_stm(adjusted_segment, kink_num, 40*units::cm, 0., 35*units::cm, true);
                    }
                }

                if (flag_pass) {
                    std::vector<std::shared_ptr<PR::Segment>> fitted_segments;
                    fitted_segments.push_back(adjusted_segment);
                    search_other_tracks(cluster, fitted_segments);
                    
                    if (check_other_tracks(cluster, fitted_segments)){
                        std::cout << "Mid Point Tracks" << std::endl;
                        return false;
                    }
                    
                    if (!detect_proton(adjusted_segment, kink_num, fitted_segments)) return true;
                }
            }
        }
        
        std::cout << "Mid Point " << std::endl;







        // // std::cout << "STMTagger tracking " << first_wcp << " " << last_wcp << std::endl;

        // // temporary tracking implementation ...
        // auto path_points = do_rough_path(cluster, first_wcp, last_wcp);
        // // Optional: Print path info for debugging
        // std::cout << "TaggerCheckSTM: Steiner path: " << path_points.size() << " points from index " << first_wcp << " " <<path_points.front() << " " << last_wcp << " " << path_points.back() << std::endl;

        // auto segment = create_segment_for_cluster(cluster, path_points);
        // // auto& wcpts = segment->wcpts();
        // // for (size_t i=0;i!=path_points.size(); i++){
        // //     std::cout << i << " " << path_points.at(i) << " " << wcpts.at(i).point << std::endl;
        // // }
        // m_track_fitter.add_segment(segment);

        // auto ch = m_track_fitter.get_channel_for_wire(0,0,1,50);
        // auto test_results = m_track_fitter.get_wires_for_channel(0,ch);
        // std::cout << ch << " " << test_results.size() << " wires. " << " " << std::get<0>(test_results.front()) << " " << std::get<1>(test_results.front()) << " " << std::get<2>(test_results.front()) << std::endl;

        // m_track_fitter.do_single_tracking(segment);


        // geo_point_t mid_point(0,0,0);
        // adjust_rough_path(cluster, mid_point);

        // auto kink_num = find_first_kink(segment);

        // check_other_clusters(cluster, associated_clusters);

        // std::vector<std::shared_ptr<PR::Segment>> fitted_segments;
        // fitted_segments.push_back(segment);
        // search_other_tracks(cluster, fitted_segments);

        // detect_proton(segment, kink_num, fitted_segments);

        // eval_stm(segment, kink_num);

        // // missing check other tracks ...
        // m_track_fitter.prepare_data();
        // m_track_fitter.fill_global_rb_map();
        // auto organized_path = m_track_fitter.organize_orig_path(segment);
        // // auto test = m_track_fitter.examine_end_ps_vec(segment, organized_path, true, true);
        // auto test_path = organized_path;
        // m_track_fitter.organize_ps_path(segment, test_path, 1.2*units::cm, 0.6*units::cm);
        // std::cout << "TaggerCheckSTM: Organized path: " << organized_path.size() << " points." << " original " << segment->wcpts().size() << " " << test_path.size() << std::endl;
        // // std::cout << m_track_fitter.get_pc_transforms() << " " << m_track_fitter.get_detector_volume() << std::endl;
        // // WireCell::Point p = organized_path.front();
        // // TrackFitting::PlaneData temp_2dut, temp_2dvt, temp_2dwt;
        // // m_track_fitter.form_point_association(segment, p, temp_2dut, temp_2dvt, temp_2dwt, 1.0*units::cm, 3, 20);
        // // m_track_fitter.examine_point_association(segment, p, temp_2dut, temp_2dvt, temp_2dwt, true);
        // // std::cout << "2D Association: " << temp_2dut.associated_2d_points.size() << " " << temp_2dut.quantity << " " << temp_2dvt.associated_2d_points.size() << " " << temp_2dvt.quantity << " " << temp_2dwt.associated_2d_points.size() << " " << temp_2dwt.quantity << std::endl;
        // std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> ptss;
        // for (const auto& p : organized_path) {
        //     ptss.emplace_back(p, segment);
        // }
        // m_track_fitter.form_map(ptss);
        // m_track_fitter.trajectory_fit(ptss);
        // m_track_fitter.dQ_dx_fit();
        // std::cout << m_track_fitter.get_parameter("DL") << std::endl;
        // std::cout << "TaggerCheckSTM: Formed map with " << organized_path.size() << " points." << std::endl;

        return false;
    }
    
 
};
