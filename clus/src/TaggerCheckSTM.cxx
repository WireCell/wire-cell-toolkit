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


class TaggerCheckSTM;
WIRECELL_FACTORY(TaggerCheckSTM, TaggerCheckSTM,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

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

    }
    
    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping_name;
        cfg["detector_volumes"] = "DetectorVolumes";
        cfg["pc_transforms"] = "PCTransformSet";  

        cfg["trackfitting_config_file"] = ""; 
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
        std::vector<Cluster*> main_clusters;

        for (auto* cluster : grouping.children()) {
            if (cluster->get_flag(Flags::main_cluster)) {
                main_clusters.push_back(cluster);
            }
        }

        std::cout << "TaggerCheckSTM: Found " << main_clusters.size() 
                  << " main clusters to check for STM conditions." << std::endl;

        // Process each main cluster
        size_t stm_count = 0;
        for (auto* cluster : main_clusters) {
            if (check_stm_conditions(*cluster)) {
                cluster->set_flag(Flags::STM);
                stm_count++;
            }
        }
        
        (void)stm_count;
    }

private:
    std::string m_grouping_name{"live"};
    std::string m_trackfitting_config_file;  // Path to TrackFitting config file
    mutable TrackFitting m_track_fitter; 

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

    int find_first_kink(const WireCell::Clus::Facade::Cluster& cluster) const{
        // Implement your logic to find the first kink in the cluster
        auto fine_tracking_path = m_track_fitter.get_fine_tracking_path();
        // auto dQ = m_track_fitter.get_dQ();
        // auto dx = m_track_fitter.get_dx();

        return -1;  // Placeholder return value
    }

    std::shared_ptr<PR::Segment> create_segment_for_cluster(WireCell::Clus::Facade::Cluster& cluster, 
                               const std::vector<geo_point_t>& path_points) const{
    
        // Step 3: Prepare segment data
        std::vector<PR::WCPoint> wcpoints;
        
        for (const auto& point : path_points) {
            PR::WCPoint wcp;
            wcp.point = point; 
            wcpoints.push_back(wcp);
        }
        
        // Step 4: Create segment connecting the vertices
        auto segment = PR::make_segment();
        
        // Step 5: Configure the segment
        segment->wcpts(wcpoints).cluster(&cluster).dirsign(1); // direction: +1, 0, or -1
               
        // auto& wcpts = segment->wcpts();
        // for (size_t i=0;i!=path_points.size(); i++){
        //     std::cout << "A: " << i << " " << path_points.at(i) << " " << wcpts.at(i).point << std::endl;
        // }

        return segment;
    }
   
    /**
     * Check if a cluster meets the conditions for STM (Short Track Muon) tagging.
     * This is where you'll implement your specific STM detection algorithm.
     * 
     * @param cluster The main cluster to analyze
     * @return true if cluster should be flagged as STM
     */
    bool check_stm_conditions(Cluster& cluster) const {
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


        // std::cout << "STMTagger tracking " << first_wcp << " " << last_wcp << std::endl;

        // temporary tracking implementation ...
        auto path_points = do_rough_path(cluster, first_wcp, last_wcp);
        // Optional: Print path info for debugging
        std::cout << "TaggerCheckSTM: Steiner path: " << path_points.size() << " points from index " << first_wcp << " " <<path_points.front() << " " << last_wcp << " " << path_points.back() << std::endl;

        auto segment = create_segment_for_cluster(cluster, path_points);
        // auto& wcpts = segment->wcpts();
        // for (size_t i=0;i!=path_points.size(); i++){
        //     std::cout << i << " " << path_points.at(i) << " " << wcpts.at(i).point << std::endl;
        // }
        m_track_fitter.add_segment(segment);

        auto ch = m_track_fitter.get_channel_for_wire(0,0,1,50);
        auto test_results = m_track_fitter.get_wires_for_channel(0,ch);
        std::cout << ch << " " << test_results.size() << " wires. " << " " << std::get<0>(test_results.front()) << " " << std::get<1>(test_results.front()) << " " << std::get<2>(test_results.front()) << std::endl;

        m_track_fitter.do_single_tracking();


        geo_point_t mid_point(0,0,0);
        adjust_rough_path(cluster, mid_point);

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