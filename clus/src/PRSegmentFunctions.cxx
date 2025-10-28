#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/DynamicPointCloud.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/KSTest.h"
#include <cmath>
#include <numeric>
#include <algorithm>

namespace WireCell::Clus::PR {
    void create_segment_point_cloud(SegmentPtr segment,
                                const std::vector<geo_point_t>& path_points,
                                const IDetectorVolumes::pointer& dv,
                                const std::string& cloud_name)
    {
        if (!segment || !segment->cluster()) {
            raise<RuntimeError>("create_segment_point_cloud: invalid segment or missing cluster");
        }
        
        auto& cluster = *segment->cluster();
        
        // Create point-plane pairs
        std::vector<std::pair<geo_point_t, WirePlaneId>> point_plane_pairs;
        for (const auto& point : path_points) {
            WirePlaneId wpid = dv->contained_by(point);
            point_plane_pairs.emplace_back(point, wpid);
        }
        
        // Get wpid_params (from detector configuration)
        const auto& wpids = cluster.grouping()->wpids();
        std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> wpid_params;
        std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_U_dir;
        std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_V_dir;
        std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_W_dir;
        std::set<int> apas;
        Facade::compute_wireplane_params(wpids, dv, wpid_params, wpid_U_dir, wpid_V_dir, wpid_W_dir, apas);
        
        // Create DynamicPointCloud
        auto dpc = std::make_shared<Facade::DynamicPointCloud>(wpid_params);
        
        // Create DPCPoints using factory function
        auto dpc_points = Facade::make_points_direct(&cluster, dv, wpid_params, point_plane_pairs, true);
        
        // Add points to DynamicPointCloud
        dpc->add_points(dpc_points);
        
        // Remove existing point cloud if it exists
        if (segment->dpcloud(cloud_name)) {
            segment->dpcloud(cloud_name, nullptr);
        }
        
        // Associate with segment
        segment->dpcloud(cloud_name, dpc);
    }

    void create_segment_fit_point_cloud(SegmentPtr segment,
                                    const IDetectorVolumes::pointer& dv,
                                    const std::string& cloud_name){
        std::vector<geo_point_t> fit_points;
        
        if (!segment || !segment->cluster()) {
            raise<RuntimeError>("create_segment_fit_point_cloud: invalid segment or missing cluster");
        }
        
        // Extract points from segment fits
        const auto& fits = segment->fits();
        fit_points.reserve(fits.size());
        for (const auto& fit : fits) {
            if (fit.valid()) {
                fit_points.push_back(fit.point);
            }
        }
        create_segment_point_cloud(segment, fit_points, dv, cloud_name);
  
    }


    std::pair<double, WireCell::Point> segment_get_closest_point(SegmentPtr seg, const WireCell::Point& point, const std::string& cloud_name){
        double min_dist = 1e9;
        WireCell::Point closest_point(0,0,0);
        
        if (!seg) {
            raise<RuntimeError>("get_closest_point: invalid segment");
        }
        
        auto dpc = seg->dpcloud(cloud_name);
        if (!dpc) {
            raise<RuntimeError>("get_closest_point: segment missing DynamicPointCloud with name " + cloud_name);
        }
        
        const auto& points = dpc->get_points();
        if (points.empty()) {
            raise<RuntimeError>("get_closest_point: DynamicPointCloud has no points");
        }
        
        // Use KD-tree to find the closest point
        auto& kd_tree = dpc->kd3d();
        auto knn_results = kd_tree.knn(1, point);
        
        if (!knn_results.empty()) {
            size_t closest_index = knn_results[0].first;
            min_dist = std::sqrt(knn_results[0].second); // knn returns squared distance
            
            // Get the actual point from the DynamicPointCloud
            const auto& dpc_point = points[closest_index];
            closest_point = WireCell::Point(dpc_point.x, dpc_point.y, dpc_point.z);
        }
        
        return {min_dist, closest_point};
    }

    std::tuple<double, double, double> segment_get_closest_2d_distances(SegmentPtr seg, const WireCell::Point& point, int apa, int face, const std::string& cloud_name) {
        if (!seg) {
            raise<RuntimeError>("segment_get_closest_2d_distances: invalid segment");
        }

        auto dpc = seg->dpcloud(cloud_name);
        if (!dpc) {
            raise<RuntimeError>("segment_get_closest_2d_distances: segment missing DynamicPointCloud with name 'fit'");
        }
        
        const auto& points = dpc->get_points();
        if (points.empty()) {
            raise<RuntimeError>("segment_get_closest_2d_distances: DynamicPointCloud has no points");
        }
        
        // Use DynamicPointCloud's optimized method to get 2D distances for each plane
        auto closest_2d_u = dpc->get_closest_2d_point_info(point, 0, face, apa);  // U plane
        auto closest_2d_v = dpc->get_closest_2d_point_info(point, 1, face, apa);  // V plane  
        auto closest_2d_w = dpc->get_closest_2d_point_info(point, 2, face, apa);  // W plane
        
        // Extract distances for each plane (U=0, V=1, W=2)
        double min_dist_u = std::get<0>(closest_2d_u);
        double min_dist_v = std::get<0>(closest_2d_v);
        double min_dist_w = std::get<0>(closest_2d_w);
        
        return std::make_tuple(min_dist_u, min_dist_v, min_dist_w);
    }

    std::tuple<WireCell::Point, WireCell::Vector, WireCell::Vector, bool> segment_search_kink(SegmentPtr seg, WireCell::Point& start_p, const std::string& cloud_name, double dQ_dx_threshold){
        auto tmp_results = segment_get_closest_point(seg, start_p, cloud_name);
        WireCell::Point test_p = tmp_results.second;

        WireCell::Vector drift_dir_abs(1,0,0);
        
        const auto& fits = seg->fits();
        if (fits.empty()) {
            WireCell::Point p1 = WireCell::Point(0,0,0);
            WireCell::Vector dir(0,0,0);
            return std::make_tuple(p1, dir, dir, false);
        }

        std::vector<double> refl_angles(fits.size(), 0);
        std::vector<double> para_angles(fits.size(), 0);
        
        // Start the angle search
        for (size_t i = 0; i < fits.size(); i++) {
            double angle1 = 0;
            double angle2 = 0;
            
            for (int j = 0; j < 6; j++) {
                WireCell::Vector v10(0,0,0);
                WireCell::Vector v20(0,0,0);
                
                if (i >= (j+1)*2) {
                    v10 = fits[i].point - fits[i-(j+1)*2].point;
                } else {
                    v10 = fits[i].point - fits.front().point;
                }
                
                if (i+(j+1)*2 < fits.size()) {
                    v20 = fits[i+(j+1)*2].point - fits[i].point;
                } else {
                    v20 = fits.back().point - fits[i].point;
                }
                
                if (j == 0) {
                    double dot_product = v10.dot(v20);
                    double mag_product = v10.magnitude() * v20.magnitude();
                    if (mag_product > 0) {
                        angle1 = std::acos(std::max(-1.0, std::min(1.0, dot_product / mag_product))) / M_PI * 180.0;
                    }
                    
                    double drift_dot1 = v10.dot(drift_dir_abs);
                    double drift_dot2 = v20.dot(drift_dir_abs);
                    double drift_mag1 = v10.magnitude();
                    double drift_mag2 = v20.magnitude();
                    
                    double drift_angle1 = 90.0, drift_angle2 = 90.0;
                    if (drift_mag1 > 0) drift_angle1 = std::acos(std::max(-1.0, std::min(1.0, drift_dot1 / drift_mag1))) / M_PI * 180.0;
                    if (drift_mag2 > 0) drift_angle2 = std::acos(std::max(-1.0, std::min(1.0, drift_dot2 / drift_mag2))) / M_PI * 180.0;
                    
                    angle2 = std::max(std::abs(drift_angle1 - 90.0), std::abs(drift_angle2 - 90.0));
                } else {
                    if (v10.magnitude() != 0 && v20.magnitude() != 0) {
                        double dot_product = v10.dot(v20);
                        double mag_product = v10.magnitude() * v20.magnitude();
                        double current_angle1 = std::acos(std::max(-1.0, std::min(1.0, dot_product / mag_product))) / M_PI * 180.0;
                        angle1 = std::min(current_angle1, angle1);
                        
                        double drift_dot1 = v10.dot(drift_dir_abs);
                        double drift_dot2 = v20.dot(drift_dir_abs);
                        double drift_mag1 = v10.magnitude();
                        double drift_mag2 = v20.magnitude();
                        
                        double drift_angle1 = 90.0, drift_angle2 = 90.0;
                        if (drift_mag1 > 0) drift_angle1 = std::acos(std::max(-1.0, std::min(1.0, drift_dot1 / drift_mag1))) / M_PI * 180.0;
                        if (drift_mag2 > 0) drift_angle2 = std::acos(std::max(-1.0, std::min(1.0, drift_dot2 / drift_mag2))) / M_PI * 180.0;
                        
                        double current_angle2 = std::max(std::abs(drift_angle1 - 90.0), std::abs(drift_angle2 - 90.0));
                        angle2 = std::min(current_angle2, angle2);
                    }
                }
            }
            
            refl_angles[i] = angle1;
            para_angles[i] = angle2;
        }

        bool flag_check = false;
        int save_i = -1;
        bool flag_switch = false;
        bool flag_search = false;
        
        for (size_t i = 0; i < fits.size(); i++) {
            // Check if close to test point
            double dist_to_test = (test_p - fits[i].point).magnitude();
            if (dist_to_test < 0.1 * units::cm) flag_check = true;
            
            // Check distance constraints
            double dist_to_front = (fits[i].point - fits.front().point).magnitude();
            double dist_to_back = (fits[i].point - fits.back().point).magnitude();
            double dist_to_start = (fits[i].point - start_p).magnitude();
            
            if (dist_to_front < 1*units::cm || 
                dist_to_back < 1*units::cm || 
                dist_to_start < 1*units::cm) continue;
            
            if (flag_check) {
                // Calculate average and max dQ/dx in local region
                double ave_dQ_dx = 0; 
                int ave_count = 0;
                double max_dQ_dx = fits[i].dQ / (fits[i].dx + 1e-9);
                
                for (int j = -2; j <= 2; j++) {
                    int idx = i + j;
                    if (idx >= 0 && idx < static_cast<int>(fits.size())) {
                        double local_dQ_dx = fits[idx].dQ / (fits[idx].dx + 1e-9);
                        ave_dQ_dx += local_dQ_dx;
                        ave_count++;
                        if (local_dQ_dx > max_dQ_dx) max_dQ_dx = local_dQ_dx;
                    }
                }
                if (ave_count != 0) ave_dQ_dx /= ave_count;
                
                // Calculate angle sums
                double sum_angles = 0;
                double nsum = 0;
                double sum_angles1 = 0;
                double nsum1 = 0;
                
                for (int j = -2; j <= 2; j++) {
                    int idx = i + j;
                    if (idx >= 0 && idx < static_cast<int>(fits.size())) {
                        if (para_angles[idx] > 10) {
                            sum_angles += pow(refl_angles[idx], 2);
                            nsum++;
                        }
                        if (para_angles[idx] > 7.5) {
                            sum_angles1 += pow(refl_angles[idx], 2);
                            nsum1++;
                        }
                    }
                }
                if (nsum != 0) sum_angles = sqrt(sum_angles / nsum);
                if (nsum1 != 0) sum_angles1 = sqrt(sum_angles1 / nsum1);
                
                // Apply kink detection criteria
                if (para_angles[i] > 10 && refl_angles[i] > 30 && sum_angles > 15) {
                    save_i = i;
                    break;
                } else if (para_angles[i] > 7.5 && refl_angles[i] > 45 && sum_angles1 > 25) {
                    save_i = i;
                    break;
                } else if (para_angles[i] > 15 && refl_angles[i] > 27 && sum_angles > 12.5) {
                    save_i = i;
                    break;
                } else if (para_angles[i] > 15 && refl_angles[i] > 22 && sum_angles > 19 && 
                          max_dQ_dx > dQ_dx_threshold*1.5 && ave_dQ_dx > dQ_dx_threshold) {
                    save_i = i;
                    flag_search = true;
                    break;
                }
            }
        }
        
        // Return results
        if (save_i > 0 && save_i+1 < static_cast<int>(fits.size())) {
            WireCell::Point p = fits[save_i].point;
            
            WireCell::Point prev_p(0,0,0);
            WireCell::Point next_p(0,0,0);
            int num_p = 0;
            int num_p1 = 0;
            
            double length1 = 0;
            double length2 = 0;
            WireCell::Point last_p1, last_p2;
            
            // Calculate direction vectors by averaging nearby points
            for (int i = 1; i < 10; i++) {
                if (save_i >= i) {
                    length1 += (fits[save_i-i].point - fits[save_i-i+1].point).magnitude();
                    prev_p = prev_p + fits[save_i-i].point;
                    last_p1 = fits[save_i-i].point;
                    num_p++;
                }
                if (save_i+i < static_cast<int>(fits.size())) {
                    length2 += (fits[save_i+i].point - fits[save_i+i-1].point).magnitude();
                    next_p = next_p + fits[save_i+i].point;
                    last_p2 = fits[save_i+i].point;
                    num_p1++;
                }
            }
            
            double length1_1 = (last_p1 - fits[save_i].point).magnitude();
            double length2_1 = (last_p2 - fits[save_i].point).magnitude();
            
            // Check for direction switch
            if (std::abs(length2 - length2_1) < 0.03 * length2_1 && length1 * length2_1 > 1.06 * length2 * length1_1) {
                flag_switch = true;
                flag_search = true;
            } else if (std::abs(length1 - length1_1) < 0.03 * length1_1 && length2 * length1_1 > 1.06 * length1 * length2_1) {
                flag_search = true;
            }
            
            prev_p = prev_p * (1.0/num_p);
            next_p = next_p * (1.0/num_p1);
            
            WireCell::Vector dir = (p - prev_p).norm();
            WireCell::Vector dir1 = (p - next_p).norm();
            
            // Calculate local charge density
            double sum_dQ = 0, sum_dx = 0;
            for (int i = -2; i <= 2; i++) {
                int idx = save_i + i;
                if (idx >= 0 && idx < static_cast<int>(fits.size())) {
                    sum_dQ += fits[idx].dQ;
                    sum_dx += fits[idx].dx;
                }
            }
            
            if (flag_search) {
                if (flag_switch) {
                    return std::make_tuple(p, dir1, dir, true);
                } else {
                    return std::make_tuple(p, dir, dir1, true);
                }
            } else if (sum_dQ / (sum_dx + 1e-9) > 25000/units::cm) { //not too low ...
                if (flag_switch) {
                    return std::make_tuple(p, dir1, dir, false);
                } else {
                    return std::make_tuple(p, dir, dir1, false);
                }
            } else {
                if (flag_switch) {
                    return std::make_tuple(p, dir1, dir, true);
                } else {
                    return std::make_tuple(p, dir, dir1, true);
                }
            }
        } else {
            WireCell::Point p1 = fits.back().point;
            WireCell::Vector dir(0,0,0);
            return std::make_tuple(p1, dir, dir, false);
        }
    }




    bool break_segment(Graph& graph, SegmentPtr seg, Point point, double max_dist/*=1e9*/)
    {
        /// sanity checks
        if (! seg->descriptor_valid()) {
            raise<RuntimeError>("break_segment: segment has invalid descriptor\n");
        }
        auto ed = seg->get_descriptor();
        auto vd1 = boost::source(ed, graph);
        auto vd2 = boost::target(ed, graph);
        auto [_, ingraph] = boost::edge(vd1, vd2, graph);
        if (! ingraph) {
            raise<RuntimeError>("break_segment: segment not in graph\n");
        }

        const auto& fits = seg->fits();
        auto itfits = closest_point(fits, point, owp_to_point<Fit>);

        // reject if test point is at begin or end of fits.
        if (itfits == fits.begin() || itfits+1 == fits.end()) {
            return false;
        }

        const auto& wcpts = seg->wcpts();        
        auto itwcpts = closest_point(wcpts, point, owp_to_point<WCPoint>);

        // clamp the wcpts to not be first/last
        if (itwcpts == wcpts.begin()) {
            ++itwcpts;
        }
        else if (itwcpts+1 == wcpts.end()) {
            --itwcpts;
        }

        
        // update graph
        remove_segment(graph, seg);

        auto vtx1 = graph[vd1].vertex;
        auto vtx2 = graph[vd2].vertex;
        auto vtx = make_vertex(graph);

        // WARNING there is no "direction" in the graph.  You can not assume the
        // "source()" of a segment is closest to the segments first point.  As
        // of now, at least...
        auto seg1 = make_segment(graph, vtx, vtx1);
        auto seg2 = make_segment(graph, vtx, vtx2);


        // fill in the new objects.  All three get the middle thing

        seg1->wcpts(std::vector<WCPoint>(wcpts.begin(), itwcpts+1));
        seg2->wcpts(std::vector<WCPoint>(itwcpts, wcpts.end()));
        vtx->wcpt(*itwcpts);

        seg1->fits(std::vector<Fit>(fits.begin(), itfits+1));
        seg2->fits(std::vector<Fit>(itfits, fits.end()));
        vtx->fit(*itfits);

        //.... more for segment
        // dir_weak
        // flags (dir, shower traj, shower topo)
        // particle type and mass and score
        // points clouds
            
        return true;
    }


    double segment_track_length(SegmentPtr seg, int flag, int n1, int n2, WireCell::Vector dir_perp)
    {
        double length = 0;
        
        if (flag == 1) {
            // Sum dx values from fits (equivalent to original flag==1 case)
            auto& fits = seg->fits();
            if (n1>=0 && n2 >=0){
                n1 = std::max(0, n1);
                n2 = std::min(static_cast<int>(fits.size())-1, n2);
                for (int i = n1; i+1 <= n2; ++i) {
                    auto& fit = fits[i];
                    if (fit.valid() && fit.dx > 0) {
                        length += fit.dx;
                    }
                }
            }else{
                for (auto& fit : fits) {
                    if (fit.valid() && fit.dx > 0) {
                        length += fit.dx;
                    }
                }
            }
        } else {
            // Calculate geometric length from fits (equivalent to original flag==0 case)
            const auto& fits = seg->fits();
            if (fits.size() < 2) {
                return 0.0;
            }
            if (n1 >=0 && n2 >=0){
                n1 = std::max(0, n1);
                n2 = std::min(static_cast<int>(fits.size())-1, n2);
                for (int i = n1; i + 1 <= n2; i++) {
                    const Point& p1 = fits[i].point;
                    const Point& p2 = fits[i + 1].point;
                    WireCell::Vector segment_vec = p2 - p1;
                    if (dir_perp.magnitude() > 0) {
                        double mag_sq = segment_vec.magnitude() * segment_vec.magnitude();
                        double dot_sq = std::pow(segment_vec.dot(dir_perp.norm()), 2);
                        length += std::sqrt(mag_sq - dot_sq);
                    }else{
                        length += segment_vec.magnitude();
                    }
                }
            }else{
                for (size_t i = 0; i + 1 < fits.size(); i++) {
                    const Point& p1 = fits[i].point;
                    const Point& p2 = fits[i + 1].point;
                    WireCell::Vector segment_vec = p2 - p1;
                    if (dir_perp.magnitude() > 0) {
                        double mag_sq = segment_vec.magnitude() * segment_vec.magnitude();
                        double dot_sq = std::pow(segment_vec.dot(dir_perp.norm()), 2);
                        length += std::sqrt(mag_sq - dot_sq);
                    }else{
                        length += segment_vec.magnitude();
                    }
                }
            }
        }
        
        return length;
    }

    double segment_track_direct_length(SegmentPtr seg, int n1, int n2, WireCell::Vector dir_perp){
        double length = 0;
        
        const auto& fits = seg->fits();
        if (fits.empty()) {
            return 0.0;
        }

        if (n1<0 && n2 <0){
            n1 = 0;
            n2 = static_cast<int>(fits.size()) - 1;
        }
        
        // Clamp indices to valid range (following WCPPID logic)
        if (n1 < 0) n1 = 0;
        if (n1 >= static_cast<int>(fits.size())) n1 = static_cast<int>(fits.size()) - 1;
        if (n2 < 0) n2 = 0;
        if (n2 >= static_cast<int>(fits.size())) n2 = static_cast<int>(fits.size()) - 1;
        
        const Point& p1 = fits[n1].point;
        const Point& p2 = fits[n2].point;
        WireCell::Vector temp_dir = p1 - p2;
        
        if (dir_perp.magnitude() > 0) {
            // Calculate length with perpendicular direction subtracted
            double mag_sq = temp_dir.magnitude() * temp_dir.magnitude();
            double dot_sq = std::pow(temp_dir.dot(dir_perp.norm()), 2);
            length = std::sqrt(mag_sq - dot_sq);
        } else {
            // Simple direct distance
            length = temp_dir.magnitude();
        }

        return length;
    }

    double segment_track_max_deviation(SegmentPtr seg, int n1, int n2){
        double max_deviation = 0.0;

        const auto& fits = seg->fits();
        if (fits.empty()) {
            return 0.0;
        }
        
        if (n1<0 && n2 <0){
            n1 = 0;
            n2 = static_cast<int>(fits.size()) - 1;
        }

        // Handle default values and clamp indices (following WCPPID logic)
        if (n1 < 0) n1 = 0;
        if (n1 >= static_cast<int>(fits.size())) n1 = static_cast<int>(fits.size()) - 1;
        if (n2 < 0) n2 = static_cast<int>(fits.size()) - 1;
        if (n2 >= static_cast<int>(fits.size())) n2 = static_cast<int>(fits.size()) - 1;
        
        // Ensure n1 <= n2
        if (n1 > n2) std::swap(n1, n2);
        
        if (n1 != n2) {
            const Point& p1 = fits[n1].point;
            const Point& p2 = fits[n2].point;
            WireCell::Vector line_vec = p2 - p1;
            double line_length_sq = line_vec.magnitude2();
            
            for (int i = n1; i <= n2; i++) {
                const Point& test_point = fits[i].point;
                
                if (line_length_sq > 0) {
                    // Calculate distance from point to line using projection
                    WireCell::Vector point_vec = test_point - p1;
                    double projection = point_vec.dot(line_vec) / line_length_sq;
                    
                    // Clamp projection to line segment bounds
                    projection = std::max(0.0, std::min(1.0, projection));
                    
                    // Find closest point on line segment
                    Point closest_on_line = p1 + line_vec * projection;
                    
                    // Calculate distance
                    double distance = (test_point - closest_on_line).magnitude();
                    
                    if (distance > max_deviation) {
                        max_deviation = distance;
                    }
                } else {
                    // Line has zero length, distance is just point-to-point distance
                    double distance = (test_point - p1).magnitude();
                    if (distance > max_deviation) {
                        max_deviation = distance;
                    }
                }
            }
        }
        
        return max_deviation;
    }
        
        


    double segment_median_dQ_dx(SegmentPtr seg, int n1, int n2)
    {
        auto& fits = seg->fits();
        if (fits.empty()) {
            return 0.0;
        }
        
        // Handle default parameters (equivalent to get_medium_dQ_dx())
        if (n1 < 0 && n2 < 0) {
            n1 = 0;
            n2 = static_cast<int>(fits.size());
        }
        
        // Clamp indices to valid range (equivalent to WCPPID bounds checking)
        if (n1 < 0) n1 = 0;
        if (n1 + 1 > static_cast<int>(fits.size())) n1 = static_cast<int>(fits.size()) - 1;
        if (n2 < 0) n2 = 0;
        if (n2 + 1 > static_cast<int>(fits.size())) n2 = static_cast<int>(fits.size()) - 1;
        
        std::vector<double> vec_dQ_dx;
        vec_dQ_dx.reserve(n2 - n1 + 1);
        
        // Loop over specified range [n1, n2] (inclusive, matching WCPPID)
        for (int i = n1; i <= n2 && i < static_cast<int>(fits.size()); i++) {
            auto& fit = fits[i];
            if (fit.valid() && fit.dx > 0 && fit.dQ >= 0) {
                // Add small epsilon to avoid division by zero (same as WCPPID: 1e-9)
                vec_dQ_dx.push_back(fit.dQ / (fit.dx + 1e-9));
            }
        }
        
        if (vec_dQ_dx.empty()) {
            return 0.0;
        }
        
        // Use nth_element to find median (exact WCPPID algorithm)
        size_t median_index = vec_dQ_dx.size() / 2;
        std::nth_element(vec_dQ_dx.begin(), 
                        vec_dQ_dx.begin() + median_index, 
                        vec_dQ_dx.end());
        
        return vec_dQ_dx[median_index];
    }
    
    double segment_rms_dQ_dx(SegmentPtr seg)
    {
        auto& fits = seg->fits();
        if (fits.empty()) {
            return 0.0;
        }
        
        std::vector<double> vec_dQ_dx;
        vec_dQ_dx.reserve(fits.size());
        
        for (auto& fit : fits) {
            if (fit.valid() && fit.dx > 0 && fit.dQ >= 0) {
                // Add small epsilon to avoid division by zero (same as original)
                vec_dQ_dx.push_back(fit.dQ / (fit.dx + 1e-9));
            }
        }
        
        if (vec_dQ_dx.empty()) {
            return 0.0;
        }
        
        // Calculate mean
        double sum = std::accumulate(vec_dQ_dx.begin(), vec_dQ_dx.end(), 0.0);
        double mean = sum / vec_dQ_dx.size();
        
        // Calculate variance
        double sq_sum = std::inner_product(vec_dQ_dx.begin(), vec_dQ_dx.end(), vec_dQ_dx.begin(), 0.0);
        double variance = sq_sum / vec_dQ_dx.size() - mean * mean;
        
        return std::sqrt(variance);
    }


    double segment_track_length_threshold(SegmentPtr seg, double threshold)
    {
        auto& fits = seg->fits();
        if (fits.empty()) {
            return 0.0;
        }
        
        double length = 0;
        for (auto& fit : fits) {
            if (fit.valid() && fit.dx > 0 ) {
                double dQ_dx = fit.dQ / (fit.dx + 1e-9); // Add epsilon to avoid division by zero
                if (dQ_dx > threshold || threshold == 0) {
                    length += fit.dx;
                }
            }
        }
        
        return length;
    }




   

    double segment_geometric_length(SegmentPtr seg, int n1, int n2, WireCell::Vector dir_perp)
    {
        return segment_track_length(seg, 0, n1, n2, dir_perp); // Always use geometric calculation
    }


    bool eval_ks_ratio(double ks1, double ks2, double ratio1, double ratio2){
        //  std::cout << ks1 << " " << ks2 << " " << ratio1 << " " << ratio2 << " " << sqrt(pow(ks2/0.06,2)+pow((ratio2-1)/0.06,2)) << " " << ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 << " " <<  ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 << " " << std::endl;
        if (ks1-ks2 >= 0.0) return false;
        if (sqrt(pow(ks2/0.06,2)+pow((ratio2-1)/0.06,2))< 1.4 && ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 > -0.02) return false;

        if (ks1 - ks2 < -0.02 && ((ks2 > 0.09 && fabs(ratio2-1) >0.1) || ratio2 > 1.5 || ks2 > 0.2)) return true;
        if ( ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 < 0) return true;

        return false;
    }

    bool segment_is_shower_trajectory(SegmentPtr seg, double step_size, double mip_dQ_dx){
        bool flag_shower_trajectory = false;
        double length = segment_track_length(seg, 0);

        // Too long
        if (length > 50 * units::cm) return flag_shower_trajectory;
        
        const auto& fits = seg->fits();
        if (fits.empty()) return flag_shower_trajectory;
        
        int ncount = std::round(length / step_size);
        if (ncount == 0) ncount = 1;
        
        std::vector<std::pair<int,int>> sections(ncount);
        for (int i = 0; i < ncount; i++) {
            sections[i] = std::make_pair(
                std::round(fits.size() / ncount * i),
                std::round(fits.size() / ncount * (i + 1))
            );
        }
        sections.back().second = fits.size() - 1;
        
        int n_shower_like = 0;
        WireCell::Vector drift_dir(1, 0, 0);
        
        for (size_t j = 0; j < ncount; j++) {
            int first_idx = sections[j].first;
            int second_idx = sections[j].second;
            
            if (first_idx >= static_cast<int>(fits.size())) first_idx = fits.size() - 1;
            if (second_idx >= static_cast<int>(fits.size())) second_idx = fits.size() - 1;
            
            WireCell::Vector dir_1 = fits[first_idx].point - fits[second_idx].point;
            if (dir_1.magnitude() > 0) {
                dir_1 = dir_1.norm();
            }
            
            double tmp_dQ_dx = segment_median_dQ_dx(seg) / (mip_dQ_dx);
            
            // Calculate angle difference
            double dot_product = drift_dir.dot(dir_1);
            double angle_rad = std::acos(std::max(-1.0, std::min(1.0, dot_product)));
            double angle_diff = std::abs(angle_rad / M_PI * 180.0 - 90.0);
            
            if (angle_diff > 10) { // Not parallel case
                double direct_length = segment_track_direct_length(seg, first_idx, second_idx, WireCell::Vector(0,0,0));
                double integrated_length = segment_track_length(seg, 0, first_idx, second_idx, WireCell::Vector(0,0,0));
                double max_dev = segment_track_max_deviation(seg, first_idx, second_idx);
                
                double length_ratio;
                if (direct_length == 0) length_ratio = 1;
                else length_ratio = direct_length / integrated_length;
                
                if (tmp_dQ_dx * 0.11 + 2 * length_ratio < 2.03 && 
                    tmp_dQ_dx < 2 && 
                    length_ratio < 0.95 && 
                    (angle_diff < 60 || integrated_length < 10 * units::cm || 
                     (integrated_length >= 10 * units::cm && max_dev > 0.75 * units::cm))) {
                    n_shower_like++;
                }
            } else { // Parallel case
                WireCell::Vector dir_2 = drift_dir.cross(dir_1);
                if (dir_2.magnitude() > 0) {
                    dir_2 = dir_2.norm();
                }
                
                double direct_length = segment_track_direct_length(seg, first_idx, second_idx, dir_2);
                double integrated_length = segment_track_length(seg, 0, first_idx, second_idx, dir_2);
                double max_dev = segment_track_max_deviation(seg, first_idx, second_idx);
                
                double length_ratio;
                if (direct_length == 0) length_ratio = 1;
                else length_ratio = direct_length / integrated_length;
                
                if (tmp_dQ_dx * 0.11 + 2 * length_ratio < 2.06 && 
                    tmp_dQ_dx < 2 && 
                    length_ratio < 0.97 && 
                    (integrated_length < 10 * units::cm || 
                     (integrated_length >= 10 * units::cm && max_dev > 0.75 * units::cm))) {
                    n_shower_like++;
                }
            }
        }
        
        if (n_shower_like >= 0.5 * sections.size()) {
            flag_shower_trajectory = true;
        }
        
        // Set the flag on the segment if it's identified as shower trajectory
        if (flag_shower_trajectory) {
            seg->set_flags(SegmentFlags::kShowerTrajectory);
        }
        
        return flag_shower_trajectory;
    }

    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg){
        const auto& fits = seg->fits();
        if (fits.size() < 2) {
            return WireCell::Vector(0, 0, 0);
        }
        
        WireCell::Point p(0, 0, 0);
        int flag_dir = seg->dirsign();
        
        if (flag_dir == 1) {
            // Forward direction: from first point using next few points
            for (size_t i = 1; i < 5 && i < fits.size(); i++) {
                p = p + (fits[i].point - fits[0].point);
            }
        } else if (flag_dir == -1) {
            // Backward direction: from last point using previous few points
            for (size_t i = 1; i < 5 && (fits.size() - i - 1) < fits.size(); i++) {
                if (fits.size() - i - 1 < fits.size()) {
                    p = p + (fits[fits.size() - i - 1].point - fits.back().point);
                }
            }
        } else {
            // Default case (flag_dir == 0): use forward direction
            for (size_t i = 1; i < 5 && i < fits.size(); i++) {
                p = p + (fits[i].point - fits[0].point);
            }
        }
        
        WireCell::Vector v1(p.x(), p.y(), p.z());
        if (v1.magnitude() > 0) {
            v1 = v1.norm();
        }
        return v1;
    }
    
    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg, WireCell::Point& p, double dis_cut){
        const auto& fits = seg->fits();
        if (fits.empty()) {
            return WireCell::Vector(0, 0, 0);
        }
        
        WireCell::Point p1(0, 0, 0);
        int ncount = 0;
        
        for (size_t i = 0; i < fits.size(); i++) {
            double dis = (fits[i].point - p).magnitude();
            if (dis < dis_cut) {
                p1 = p1 + fits[i].point;
                ncount++;
            }
        }
        
        if (ncount == 0) {
            return WireCell::Vector(0, 0, 0);
        }
        
        WireCell::Point avg_point = p1 * (1.0 / ncount);
        WireCell::Vector v1 = avg_point - p;
        if (v1.magnitude() > 0) {
            v1 = v1.norm();
        }
        return v1;
    }
    
    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg, int direction, int num_points, int start){
        const auto& fits = seg->fits();
        if (fits.empty() || start >= static_cast<int>(fits.size()) || start <= 0) {
            std::cout << "bad start point in segment_cal_dir_3vector" << std::endl;
            return WireCell::Vector(0, 0, 0);
        }
        
        WireCell::Point p(0, 0, 0);
        
        if (direction == 1) {
            // Forward direction
            for (int i = start; i < start + num_points - 1 && i < static_cast<int>(fits.size()); i++) {
                p = p + (fits[i].point - fits[start - 1].point);
            }
        } else if (direction == -1) {
            // Backward direction
             for (int i = start; i < start + num_points - 1; i++) {
                // WCPPID's bounds check
                if (i + 1 > static_cast<int>(fits.size())) break;
                
                // Ensure backward indices are valid
                int back_idx = fits.size() - i - 1;
                int ref_idx = fits.size() - start;
                
                if (back_idx >= 0 && back_idx < static_cast<int>(fits.size()) && 
                    ref_idx >= 0 && ref_idx < static_cast<int>(fits.size())) {
                    p = p + (fits[back_idx].point - fits[ref_idx].point);
                }
            }
        }
        
        WireCell::Vector v1(p.x(), p.y(), p.z());
        if (v1.magnitude() > 0) {
            v1 = v1.norm();
        }
        return v1;
    }

    double segment_cal_kine_dQdx(SegmentPtr seg, const IRecombinationModel::pointer& recomb_model){
        if (!seg || !recomb_model) {
            return 0.0;
        }
        
        auto& fits = seg->fits();
        if (fits.empty()) {
            return 0.0;
        }
        
        double kine_energy = 0.0;
        
        
        for (size_t i = 0; i < fits.size(); i++) {
            if (!fits[i].valid() || fits[i].dx <= 0) continue;
            
            double dX = fits[i].dx;
            double dQ = fits[i].dQ;
            if (i == 0 && fits.size() > 1) {
                // First point: check against distance to next point
                double dis = (fits[1].point - fits[0].point).magnitude();
                if (dX> dis * 1.5) {
                    dX = dis;
                }
            } else if (i + 1 == fits.size() && fits.size() > 1) {
                // Last point: check against distance to previous point
                double dis = (fits[i].point - fits[i-1].point).magnitude();
                if (dX > dis * 1.5) {
                    dX = dis;
                }
            }
            // std::cout << i << " " << fits[i].dQ << " " << fits[i].dx/units::cm << " " << dX/units::cm << std::endl;
            // Filter out unreasonable values (same threshold as original)
            if (dQ/dX / (43e3/units::cm) > 1000) dQ = 0;
            
            // Calculate dE/dx using Box model inverse formula from original code
            double dE = recomb_model->dE(dQ, dX);

            // std::cout << dQ << " " << dX << " " << dE << std::endl;

            // Apply bounds (same as original)
            if (dE < 0) dE = 0;
            if (dE > 50 * units::MeV / units::cm * dX) dE = 50 * units::MeV / units::cm * dX;

            // Calculate path length with special handling for first and last points
            kine_energy += dE;
        }
        
        return kine_energy;
    }
    
    double cal_kine_dQdx(std::vector<double>& vec_dQ, std::vector<double>& vec_dx, const IRecombinationModel::pointer& recomb_model){
        if (vec_dQ.size() != vec_dx.size() || vec_dQ.empty() || !recomb_model) {
            return 0.0;
        }
        
        double kine_energy = 0.0;
        
        for (size_t i = 0; i < vec_dQ.size(); i++) {
              // Calculate dQ/dx with units conversion (same as original)
            double dQ = vec_dQ[i];
            double dx = vec_dx[i];
            
            // Filter out unreasonable values (same threshold as original)
            if (dQ/dx / (43e3/units::cm) > 1000) dQ = 0;
            
            // Calculate dE/dx using Box model inverse formula from original code
            double dE = recomb_model->dE(dQ, dx);

            // double dQp = (*recomb_model)(dE, dx);

            // std::cout << dQ << " " << dx << " " << dE << " " << units::MeV << " " << dQp << std::endl;
            
            // Apply bounds (same as original)
            if (dE < 0) dE = 0;
            if (dE > 50 * units::MeV / units::cm * dx) dE = 50 * units::MeV / units::cm * dx;

            // Calculate path length with special handling for first and last points
            kine_energy += dE;
        }
        
        return kine_energy;
    }

    std::vector<double> do_track_comp(std::vector<double>& L , std::vector<double>& dQ_dx, double compare_range, double offset_length, const Clus::ParticleDataSet::pointer& particle_data,  double MIP_dQdx){
        
        double end_L = L.back() + 0.15*units::cm - offset_length;
        
        int ncount = 0;
        std::vector<double> vec_x;
        std::vector<double> vec_y;
        
        for (size_t i = 0; i != L.size(); i++) {
            if (end_L - L.at(i) < compare_range && end_L - L.at(i) > 0) { // check up to compared range
                vec_x.push_back(end_L - L.at(i));
                vec_y.push_back(dQ_dx.at(i));
                ncount++;
            }
        }
        
        // Create reference vectors for different particles
        std::vector<double> muon_ref(ncount);
        std::vector<double> const_ref(ncount, MIP_dQdx);  // MIP-like constant
        std::vector<double> proton_ref(ncount);
        std::vector<double> electron_ref(ncount);
        
        for (size_t i = 0; i != ncount; i++) {
            muon_ref[i] = particle_data->get_dEdx_function("muon")->scalar_function((vec_x[i])/units::cm) /units::cm;
            proton_ref[i] = particle_data->get_dEdx_function("proton")->scalar_function((vec_x[i])/units::cm) / units::cm;
            electron_ref[i] = particle_data->get_dEdx_function("electron")->scalar_function((vec_x[i])/units::cm) / units::cm;
        }
        
        // Perform KS-like tests using kslike_compare
        double ks1 = WireCell::kslike_compare(vec_y, muon_ref);
        double ratio1 = std::accumulate(muon_ref.begin(), muon_ref.end(), 0.0) / 
                        (std::accumulate(vec_y.begin(), vec_y.end(), 0.0) + 1e-9);
        
        double ks2 = WireCell::kslike_compare(vec_y, const_ref);
        double ratio2 = std::accumulate(const_ref.begin(), const_ref.end(), 0.0) / 
                        (std::accumulate(vec_y.begin(), vec_y.end(), 0.0) + 1e-9);
        
        double ks3 = WireCell::kslike_compare(vec_y, proton_ref);
        double ratio3 = std::accumulate(proton_ref.begin(), proton_ref.end(), 0.0) / 
                        (std::accumulate(vec_y.begin(), vec_y.end(), 0.0) + 1e-9);
        
        double ks4 = WireCell::kslike_compare(vec_y, electron_ref);
        double ratio4 = std::accumulate(electron_ref.begin(), electron_ref.end(), 0.0) / 
                        (std::accumulate(vec_y.begin(), vec_y.end(), 0.0) + 1e-9);
        
        // std::cout << ks1 << " " << ratio1 << " " << ks2 << " " << ratio2 << " " << ks3 << " " << ratio3 << " " << ks4 << " " << ratio4 << std::endl;

        std::vector<double> results;
        // Convert bool result to double (1.0 for true, 0.0 for false)
        results.push_back(eval_ks_ratio(ks1, ks2, ratio1, ratio2) ? 1.0 : 0.0); // direction metric
        results.push_back(sqrt(pow(ks1, 2) + pow(ratio1-1, 2))); // muon information
        results.push_back(sqrt(pow(ks3, 2) + pow(ratio3-1, 2))); // proton information  
        results.push_back(sqrt(pow(ks4, 2) + pow(ratio4-1, 2))); // electron information
        
        return results;
    }
   
    double cal_kine_range(double L, int pdg_code, const Clus::ParticleDataSet::pointer& particle_data){

        IScalarFunction::pointer range_function = nullptr;

        if (abs(pdg_code) == 11) {        // electron
            range_function = particle_data->get_range_function("electron");
        }
        else if (abs(pdg_code) == 13) {   // muon
            range_function = particle_data->get_range_function("muon");
        }
        else if (abs(pdg_code) == 211) {  // pion
            range_function = particle_data->get_range_function("pion");
        }
        else if (abs(pdg_code) == 321) {  // kaon
            range_function = particle_data->get_range_function("kaon");
        }
        else if (abs(pdg_code) == 2212) { // proton
            range_function = particle_data->get_range_function("proton");
        }
        
        if (!range_function) {
            // Default to muon if particle type not recognized
            range_function = particle_data->get_range_function("muon"); 
        }
        double kine_energy = range_function->scalar_function(L/units::cm) * units::MeV;
        return kine_energy;
    }

    // success, flag_dir, particle_type, particle_score
    std::tuple<bool, int, int, double> segment_do_track_pid(SegmentPtr segment, std::vector<double>& L , std::vector<double>& dQ_dx, double compare_range , double offset_length, bool flag_force, const Clus::ParticleDataSet::pointer& particle_data, double MIP_dQdx){
        
        if (L.size() != dQ_dx.size() || L.empty() || !segment) {
            return std::make_tuple(false, 0, 0, 0.0);
        }
        
        std::vector<double> rL(L.size(), 0);
        std::vector<double> rdQ_dx(L.size(), 0);
        
        // Get reverse vectors
        for (size_t i = 0; i != L.size(); i++) {
            rL.at(i) = L.back() - L.at(L.size() - 1 - i);
            rdQ_dx.at(i) = dQ_dx.at(L.size() - 1 - i);
        }
        
        std::vector<double> result_forward = do_track_comp(L, dQ_dx, compare_range, offset_length, particle_data, MIP_dQdx);
        std::vector<double> result_backward = do_track_comp(rL, rdQ_dx, compare_range, offset_length, particle_data, MIP_dQdx);
        
        // Direction determination
        bool flag_forward = static_cast<bool>(std::round(result_forward.at(0)));
        bool flag_backward = static_cast<bool>(std::round(result_backward.at(0)));
        
        // Calculate length from path (total walk length over fits or wcpts)
        double length = segment_track_length(segment, 0);


        // // Calculate straight-line distance between endpoints (length1 equivalent)
        // double length1 = 0.0;
        // auto& fits = segment->fits();
        // length1 = (fits.front().point - fits.back().point).magnitude();
        
        // Forward particle type determination
        int forward_particle_type = 13; // default muon
        double min_forward_val = result_forward.at(1);
        if (result_forward.at(2) < min_forward_val) {
            min_forward_val = result_forward.at(2);
            forward_particle_type = 2212; // proton
        }
        if (result_forward.at(3) < min_forward_val && length < 20*units::cm) {
            min_forward_val = result_forward.at(3);
            forward_particle_type = 11; // electron
        }
        
        // Backward particle type determination  
        int backward_particle_type = 13; // default muon
        double min_backward_val = result_backward.at(1);
        if (result_backward.at(2) < min_backward_val) {
            min_backward_val = result_backward.at(2);
            backward_particle_type = 2212; // proton
        }
        if (result_backward.at(3) < min_backward_val && length < 20*units::cm) {
            min_backward_val = result_backward.at(3);
            backward_particle_type = 11; // electron
        }
        
        // Decision logic
        int flag_dir = 0;
        int particle_type = 0;
        double particle_score = 0.0;
        
        if (flag_forward == 1 && flag_backward == 0) {
            flag_dir = 1;
            particle_type = forward_particle_type;
            particle_score = min_forward_val;
            return std::make_tuple(true, flag_dir, particle_type, particle_score);
        }
        else if (flag_forward == 0 && flag_backward == 1) {
            flag_dir = -1;
            particle_type = backward_particle_type;
            particle_score = min_backward_val;
            return std::make_tuple(true, flag_dir, particle_type, particle_score);
        }
        else if (flag_forward == 1 && flag_backward == 1) {
            if (min_forward_val < min_backward_val) {
                flag_dir = 1;
                particle_type = forward_particle_type;
                particle_score = min_forward_val;
            }
            else {
                flag_dir = -1;
                particle_type = backward_particle_type;
                particle_score = min_backward_val;
            }
            return std::make_tuple(true, flag_dir, particle_type, particle_score);
        }
        else if (flag_forward == 0 && flag_backward == 0 && flag_force) {
            if (min_forward_val < min_backward_val) {
                particle_score = min_forward_val;
                particle_type = forward_particle_type;
                flag_dir = 1;
            }
            else {
                particle_score = min_backward_val;
                particle_type = backward_particle_type;
                flag_dir = -1;
            }
            return std::make_tuple(true, flag_dir, particle_type, particle_score);
        }
        
        // Reset before return - failure case
        return std::make_tuple(false, 0, 0, 0.0);
    }

    // 4-momentum: E, px, py, pz, 
    WireCell::D4Vector<double> segment_cal_4mom(SegmentPtr segment, int pdg_code, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx){
        double length = segment_track_length(segment, 0);
        double kine_energy = 0;

        WireCell::D4Vector<double> results(0.0, 0.0, 0.0, 0.0); // 4-momentum: E, px, py, pz

        if (length < 4*units::cm){
            kine_energy = segment_cal_kine_dQdx(segment, recomb_model); // short track 
        }else if (segment->flags_any(PR::SegmentFlags::kShowerTrajectory)){
            kine_energy = segment_cal_kine_dQdx(segment, recomb_model);
        }else{
            kine_energy = cal_kine_range(length, pdg_code, particle_data);
        }
        // results[4] = kine_energy;

        double particle_mass = particle_data->get_particle_mass(pdg_code);

        results[0]= kine_energy + particle_mass;
        double mom = sqrt(pow(results[3],2) - pow(particle_mass,2));
        auto v1 = segment_cal_dir_3vector(segment);
        results[1] = mom * v1.x();
        results[2] = mom * v1.y();
        results[3] = mom * v1.z();

        return results;
    }

    void segment_determine_dir_track(SegmentPtr segment, int start_n, int end_n, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx, bool flag_print) {
        if (!segment || !particle_data) {
            return;
        }
        
        // Reset direction flag
        segment->dirsign(0);
        
        const auto& fits = segment->fits();
        int npoints = fits.size();
        int start_n1 = 0, end_n1 = npoints - 1;
        
        // If more than one point, exclude the vertex
        if (end_n != 1) {
            end_n1 = npoints - 2;
            npoints -= 1;
        }
        if (start_n != 1) {
            npoints -= 1;
            start_n1 = 1;
        }
        
        if (npoints == 0 || end_n1 < start_n1) return;
        
        std::vector<double> L(npoints, 0);
        std::vector<double> dQ_dx(npoints, 0);
        
        double dis = 0;
        for (int i = start_n1; i <= end_n1; i++) {
            L.at(i - start_n1) = dis;
            if (fits[i].dx > 0) {
                dQ_dx.at(i - start_n1) = fits[i].dQ / (fits[i].dx + 1e-9);
            }
            if (i + 1 < static_cast<int>(fits.size())) {
                dis += (fits[i+1].point - fits[i].point).magnitude();
            }
        }
        
        int pdg_code = 0;
        double particle_score = 0.0;
        
        if (npoints >= 2) { // reasonably long
            bool tmp_flag_pid = false;
            
            if (start_n == 1 && end_n == 1 && npoints >= 15) {
                // Can use the dQ/dx to do PID and direction
                auto result = segment_do_track_pid(segment, L, dQ_dx, 35*units::cm, 1*units::cm, true, particle_data);
                tmp_flag_pid = std::get<0>(result);
                if (tmp_flag_pid) {
                    segment->dirsign(std::get<1>(result));
                    pdg_code = std::get<2>(result);
                    particle_score = std::get<3>(result);
                }
                
                if (!tmp_flag_pid) {
                    result = segment_do_track_pid(segment, L, dQ_dx, 15*units::cm, 1*units::cm, true, particle_data);
                    tmp_flag_pid = std::get<0>(result);
                    if (tmp_flag_pid) {
                        segment->dirsign(std::get<1>(result));
                        pdg_code = std::get<2>(result);
                        particle_score = std::get<3>(result);
                    }
                }
            } else {
                // Can use the dQ/dx to do PID and direction
                auto result = segment_do_track_pid(segment, L, dQ_dx, 35*units::cm, 0*units::cm, false, particle_data);
                tmp_flag_pid = std::get<0>(result);
                if (tmp_flag_pid) {
                    segment->dirsign(std::get<1>(result));
                    pdg_code = std::get<2>(result);
                    particle_score = std::get<3>(result);
                }
                
                if (!tmp_flag_pid) {
                    result = segment_do_track_pid(segment, L, dQ_dx, 15*units::cm, 0*units::cm, false, particle_data);
                    tmp_flag_pid = std::get<0>(result);
                    if (tmp_flag_pid) {
                        segment->dirsign(std::get<1>(result));
                        pdg_code = std::get<2>(result);
                        particle_score = std::get<3>(result);
                    }
                }
            }
        }
        
        double length = segment_track_length(segment, 0);
        
        // Short track what to do???
        if (pdg_code == 0) {
            // Calculate median dQ/dx
            double medium_dQ_dx = segment_median_dQ_dx(segment);
            if (medium_dQ_dx > MIP_dQdx * 1.75) {
                pdg_code = 2212; // proton
            } else if (medium_dQ_dx < MIP_dQdx * 1.2) {
                pdg_code = 13; // muon
            } else if (medium_dQ_dx < MIP_dQdx * 1.5 && length < 4*units::cm) {
                pdg_code = 13;
            }
        }
        
        // Electron and both end contain stuff
        if (abs(pdg_code) == 11 && (start_n > 1 && end_n > 1)) {
            segment->dir_weak(true);
            segment->dirsign(0);
            if (particle_score < 0.15) pdg_code = 13;
        } else if (abs(pdg_code) == 11 && ((start_n > 1 && segment->dirsign() == -1) || (end_n > 1 && segment->dirsign() == 1))) {
            segment->dir_weak(true);
            segment->dirsign(0);
            if (particle_score < 0.15) pdg_code = 13;
        } else if (length < 1.5*units::cm) {
            segment->dir_weak(true);
        }
        
        // Vertex activities
        if (length < 1.5*units::cm && (start_n == 1 || end_n == 1)) {
            if (start_n == 1 && end_n > 2) {
                segment->dirsign(-1);
                double medium_dQ_dx = segment_median_dQ_dx(segment);
                if (medium_dQ_dx > MIP_dQdx * 1.75) {
                    pdg_code = 2212;
                } else if (medium_dQ_dx < MIP_dQdx * 1.2) {
                    pdg_code = 211;
                }
            } else if (end_n == 1 && start_n > 2) {
                segment->dirsign(1);
                double medium_dQ_dx = segment_median_dQ_dx(segment);
                if (medium_dQ_dx > MIP_dQdx * 1.75) {
                    pdg_code = 2212;
                } else if (medium_dQ_dx < MIP_dQdx * 1.2) {
                    pdg_code = 211;
                }
            }
        }
        
        // If the particle score is really bad, make it a shower
        if (length > 10*units::cm && particle_score > 1.0 && particle_score < 100) {
            pdg_code = 11;
            particle_score = 200;
            segment->dirsign(0);
        }
        
        // Set particle mass and calculate 4-momentum
        if (pdg_code != 0) {
            // Calculate 4-momentum using the identified particle type
            auto four_momentum = segment_cal_4mom(segment, pdg_code, particle_data, recomb_model);

            // Create ParticleInfo with the identified particle
            auto pinfo = std::make_shared<Aux::ParticleInfo>(
                pdg_code,                    // PDG code
                particle_data->get_particle_mass(pdg_code), // mass
                particle_data->pdg_to_name(pdg_code),       // name
                four_momentum                     // 4-momentum
            );
            
            // Set additional properties if available
            pinfo->set_particle_score(particle_score); // This method would need to be added
            
            // Store particle info in segment (this would require adding particle_info to Segment class)
            segment->particle_info(pinfo);
        }
                
        if (flag_print && pdg_code != 0) {
            std::cout << "Segment PID: PDG=" << pdg_code 
                      << ", Score=" << particle_score 
                      << ", Length=" << length / units::cm << " cm"
                      << ", Direction=" << segment->dirsign() 
                      << (segment->dir_weak() ? " (weak)" : "") 
                      << ", Medium dQ/dx=" << segment_median_dQ_dx(segment) / (MIP_dQdx) 
                      << " MIP"
                      << std::endl;
        }
    }

     void segment_determine_shower_direction_trajectory(SegmentPtr segment, int start_n, int end_n, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx, bool flag_print){
        segment->dirsign(0);
        double length = segment_track_length(segment, 0);
        
        // hack for now ...
        int pdg_code = 11;
        
        if (start_n==1 && end_n >1){
            segment->dirsign(-1);
        }else if (start_n > 1 && end_n == 1){
            segment->dirsign(1);
        }else{
            segment_determine_dir_track(segment, start_n, end_n, particle_data, recomb_model, MIP_dQdx, false);
            if (segment->particle_info()->pdg() != 11){
                segment->dirsign(0);
            }
        }
                
        auto four_momentum = segment_cal_4mom(segment, pdg_code, particle_data, recomb_model);

        // Create ParticleInfo with the identified particle
        auto pinfo = std::make_shared<Aux::ParticleInfo>(
            pdg_code,                    // PDG code
            particle_data->get_particle_mass(pdg_code), // mass
            particle_data->pdg_to_name(pdg_code),       // name
            four_momentum                     // 4-momentum
        );
                
        // Store particle info in segment (this would require adding particle_info to Segment class)
        segment->particle_info(pinfo);

         if (flag_print ) {
            std::cout << "Segment PID: PDG=" << pdg_code 
                      << ", Length=" << length / units::cm << " cm"
                      << ", Direction=" << segment->dirsign() 
                      << (segment->dir_weak() ? " (weak)" : "") 
                      << ", Medium dQ/dx=" << segment_median_dQ_dx(segment) / (MIP_dQdx) 
                      << " MIP"
                      << std::endl;
        }


     }

    void clustering_points_segments(std::set<SegmentPtr> segments, const IDetectorVolumes::pointer& dv, const std::string& cloud_name, double search_range, double scaling_2d){
        std::map<Facade::Cluster*, std::set<SegmentPtr> > map_cluster_segs;
        for (auto seg : segments){
            if (seg->cluster()){
                map_cluster_segs[seg->cluster()].insert(seg); 
            }
        }
        
        for (auto it : map_cluster_segs){
            auto clus = it.first;
            auto& segs = it.second;
            
            // get the default point cloud from cluster
            const auto& points = clus->points();

             // Get the graph directly
            const auto& graph = clus->find_graph("basic_pid");

            std::map<SegmentPtr, std::vector<geo_point_t>> map_segment_points;
            std::map<int, std::pair<SegmentPtr, double>> map_pindex_segment;
            //std::cout << "Cluster has " << npoints << " points and " << segs.size() << " segments." << std::endl;
            //std::cout << "Number of vertices in the graph: " << boost::num_vertices(graph) << std::endl;
            
            // core algorithms
             
            // define steiner terminal for segments ...
            for (auto seg: segs){
                auto& fits = seg->fits();
                if (fits.size() > 2){
                    for (size_t i = 1; i+1 < fits.size(); i++){
                        geo_point_t gp = {fits[i].point.x(), fits[i].point.y(), fits[i].point.z()};
                        // use cluster to get the indices of the closest 5 points
                        auto closest_results = clus->kd_knn(5, gp);
                        for (const auto& [point_index, distance] : closest_results) {
                            if (map_pindex_segment.find(point_index) == map_pindex_segment.end()) {
                                map_pindex_segment[point_index] = std::make_pair(seg, distance);
                                break;
                            }
                        }
                    }
                }else{
                    geo_point_t gp = {(fits[0].point.x()+fits[1].point.x())/2., (fits[0].point.y()+fits[1].point.y())/2., (fits[0].point.z()+fits[1].point.z())/2.};
                    auto closest_results = clus->kd_knn(5, gp);
                    for (const auto& [point_index, distance] : closest_results) {
                        if (map_pindex_segment.find(point_index) == map_pindex_segment.end()) {
                            map_pindex_segment[point_index] = std::make_pair(seg, distance);
                            break;
                        }
                    }
                }
            }


            // these are terminals ...
            if (map_pindex_segment.size()>0){
               
                // Convert terminals from int to vertex_type
                std::vector<WireCell::Clus::Graphs::Weighted::vertex_type> terminals;
                for (auto it = map_pindex_segment.begin(); it!=map_pindex_segment.end(); it++){
                    terminals.push_back(static_cast<WireCell::Clus::Graphs::Weighted::vertex_type>(it->first));
                }

                auto vor = WireCell::Clus::Graphs::Weighted::voronoi(graph, terminals);
                
                // Now we can find the nearest terminal for every vertex in the graph
                // The Voronoi diagram provides:
                // - vor.terminal[v]: the nearest terminal vertex for vertex v
                // - vor.distance[v]: the distance to the nearest terminal for vertex v
                std::map<int, std::pair<int, double>> vertex_to_nearest_terminal;
                // Iterate through all vertices in the graph
                const int num_graph_vertices = boost::num_vertices(graph);
                for (int vertex_idx = 0; vertex_idx < num_graph_vertices; ++vertex_idx) {
                    // Get the nearest terminal for this vertex
                    int nearest_terminal_idx = vor.terminal[vertex_idx];
                    double distance_to_terminal = vor.distance[vertex_idx];
                    // Store the mapping
                    vertex_to_nearest_terminal[vertex_idx] = std::make_pair(nearest_terminal_idx, distance_to_terminal);
                }
                // std::cout << "Debug: Number of graph vertices: " << num_graph_vertices << std::endl;
                // now examine to remove ghost points ....
                for (size_t i=0;i!=num_graph_vertices;i++){
                    if (map_pindex_segment.find(vertex_to_nearest_terminal.at(i).first) == map_pindex_segment.end()) continue;
                    geo_point_t gp(points[0][i], points[1][i], points[2][i]);
                    auto main_sg = map_pindex_segment[vertex_to_nearest_terminal.at(i).first].first;

                    auto point_wpid = clus->wire_plane_id(i);
                    auto apa = point_wpid.apa();
                    auto face = point_wpid.face();

                    // use the dynamic point cloud of fit, and then derive distances ... 
                    // Get 3D closest point using the "fit" point cloud
                    std::pair<double, WireCell::Point> closest_dis_point = segment_get_closest_point(main_sg, gp, "fit");

                    // Calculate 2D distances for each wire plane (U, V, W) using APA/face information
                    std::tuple<double, double, double> closest_2d_dis = segment_get_closest_2d_distances(main_sg, gp, apa, face, "fit");
                    
                    std::tuple<double, double, double> min_2d_dis = closest_2d_dis;
                
                    // check against main_sg;
                    bool flag_change = true;
                
                    // Compare against all segments in the cluster to find minimum 2D distances
                    for (auto seg : segs) {
                       if (main_sg == seg) continue;
                    
                        // Get 2D distances for this segment
                        std::tuple<double, double, double> temp_2d_dis = segment_get_closest_2d_distances(seg, gp, apa, face, "fit");
                        // Update minimum distances for each plane
                        if (std::get<0>(temp_2d_dis) < std::get<0>(min_2d_dis)) std::get<0>(min_2d_dis) = std::get<0>(temp_2d_dis);
                        if (std::get<1>(temp_2d_dis) < std::get<1>(min_2d_dis)) std::get<1>(min_2d_dis) = std::get<1>(temp_2d_dis);
                        if (std::get<2>(temp_2d_dis) < std::get<2>(min_2d_dis)) std::get<2>(min_2d_dis) = std::get<2>(temp_2d_dis);
                    }
                
                    if (std::get<0>(min_2d_dis) == std::get<0>(closest_2d_dis) && std::get<1>(min_2d_dis) == std::get<1>(closest_2d_dis) && std::get<2>(min_2d_dis) == std::get<2>(closest_2d_dis)) // all closest
                    flag_change = false;
                    else if (std::get<0>(min_2d_dis) == std::get<0>(closest_2d_dis) && std::get<1>(min_2d_dis) == std::get<1>(closest_2d_dis) ) //&& (std::get<2>(closest_2d_dis) < scaling_2d * search_range || closest_dis_point.first < search_range)) // 2 closest
                    flag_change = false;
                    else if (std::get<0>(min_2d_dis) == std::get<0>(closest_2d_dis) && std::get<2>(min_2d_dis) == std::get<2>(closest_2d_dis) ) //&& (std::get<1>(closest_2d_dis) < scaling_2d * search_range || closest_dis_point.first < search_range))
                    flag_change = false;
                    else if (std::get<1>(min_2d_dis) == std::get<1>(closest_2d_dis) && std::get<2>(min_2d_dis) == std::get<2>(closest_2d_dis) ) //&& (std::get<0>(closest_2d_dis) < scaling_2d * search_range || closest_dis_point.first < search_range))
                    flag_change = false;
                    else if (std::get<0>(min_2d_dis) == std::get<0>(closest_2d_dis) && (closest_dis_point.first < search_range || (std::get<1>(closest_2d_dis) < scaling_2d * search_range &&  std::get<2>(closest_2d_dis) < scaling_2d * search_range)) )
                    flag_change = false;
                    else if (std::get<1>(min_2d_dis) == std::get<1>(closest_2d_dis) && (closest_dis_point.first < search_range || (std::get<0>(closest_2d_dis) < scaling_2d * search_range &&  std::get<2>(closest_2d_dis) < scaling_2d * search_range) ))
                    flag_change = false;
                    else if (std::get<2>(min_2d_dis) == std::get<2>(closest_2d_dis) && (closest_dis_point.first < search_range || (std::get<1>(closest_2d_dis) < scaling_2d * search_range &&  std::get<0>(closest_2d_dis) < scaling_2d * search_range) ))
                    flag_change = false;

                    // deal with dead channels ...
                    if (!flag_change){
                        auto grouping = clus->grouping();
                        int ch_range = 0; // Default channel range for dead channel checking
                        
                        // Check U plane (pind=0) for dead channels
                        if (grouping->get_closest_dead_chs(gp, ch_range, apa, face, 0) && std::get<0>(closest_2d_dis) > scaling_2d * search_range){
                            if (std::get<1>(closest_2d_dis) < scaling_2d * search_range ||  std::get<2>(closest_2d_dis) < scaling_2d * search_range)
                                flag_change = true;
                        // Check V plane (pind=1) for dead channels
                        }else if (grouping->get_closest_dead_chs(gp, ch_range, apa, face, 1) && std::get<1>(closest_2d_dis) > scaling_2d * search_range){
                            if (std::get<0>(closest_2d_dis) < scaling_2d * search_range ||  std::get<2>(closest_2d_dis) < scaling_2d * search_range)
                                flag_change = true;
                        // Check W plane (pind=2) for dead channels
                        }else if (grouping->get_closest_dead_chs(gp, ch_range, apa, face, 2) && std::get<2>(closest_2d_dis) > scaling_2d * search_range){
                            if (std::get<1>(closest_2d_dis) < scaling_2d * search_range ||  std::get<0>(closest_2d_dis) < scaling_2d * search_range)
                                flag_change = true;
                        }
                   } 
                
                    // change the point's clustering ...
                    if (!flag_change){
                        map_segment_points[main_sg].push_back(gp);
                    }
                }
            }


            // convert points to geo_point_t format
            // add points to segments ... 
            for (const auto& [seg, geo_points] : map_segment_points) {
                create_segment_point_cloud(seg, geo_points, dv, cloud_name);
            }
        }
    }

    bool segment_determine_shower_direction(SegmentPtr segment, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, const std::string& cloud_name, double MIP_dQdx, double rms_cut){
        segment->dirsign(0);
        const auto& fits = segment->fits();
        
        if (fits.empty()) return false;
        
        // Get the fit point cloud for KD-tree queries
        auto dpcloud_fit = segment->dpcloud("fit");
        if (!dpcloud_fit) return false;
        
        // Get the associated point cloud 
        auto dpcloud_assoc = segment->dpcloud(cloud_name);
        if (!dpcloud_assoc) return false;
        
        const auto& assoc_points = dpcloud_assoc->get_points();
        if (assoc_points.empty()) return false;
        
        // Initialize vectors to store analysis results for each fit point
        std::vector<std::vector<WireCell::Point>> local_points_vec(fits.size());
        std::vector<std::tuple<double, double, double>> vec_rms_vals(fits.size(), std::make_tuple(0,0,0));
        std::vector<double> vec_dQ_dx(fits.size(), 0);
        std::vector<WireCell::Vector> vec_dir(fits.size());
        
        // Build KD-tree index for fit points and associate points with nearest fit point
        auto& kd_tree_fit = dpcloud_fit->kd3d();
        
        for (const auto& pt : assoc_points) {
            WireCell::Point test_p(pt.x, pt.y, pt.z);
            auto results = kd_tree_fit.knn(1, test_p);
            if (!results.empty()) {
                size_t closest_fit_idx = results.front().first;
                local_points_vec.at(closest_fit_idx).push_back(test_p);
            }
        }
        
        WireCell::Vector drift_dir_abs(1, 0, 0);  // Drift direction
        
        // Calculate local directions and RMS spreads for each fit point
        for (size_t i = 0; i < local_points_vec.size(); i++) {
            // Calculate local direction from neighboring fit points
            WireCell::Vector v1(0, 0, 0);
            for (size_t j = 1; j < 3; j++) {
                if (i + j < fits.size()) {
                    v1 += WireCell::Vector(
                        fits[i+j].point.x() - fits[i].point.x(),
                        fits[i+j].point.y() - fits[i].point.y(),
                        fits[i+j].point.z() - fits[i].point.z()
                    );
                }
                if (i >= j) {
                    v1 += WireCell::Vector(
                        fits[i].point.x() - fits[i-j].point.x(),
                        fits[i].point.y() - fits[i-j].point.y(),
                        fits[i].point.z() - fits[i-j].point.z()
                    );
                }
            }
            
            WireCell::Vector dir_1 = v1.magnitude() > 0 ? v1.norm() : WireCell::Vector(1, 0, 0);
            vec_dir.at(i) = dir_1;
            
            // Set up orthogonal coordinate system
            WireCell::Vector dir_2, dir_3;
            double angle_deg = std::acos(dir_1.dot(drift_dir_abs)) * 180.0 / M_PI;
            
            if (angle_deg < 7.5) {
                dir_1 = WireCell::Vector(1, 0, 0);
                dir_2 = WireCell::Vector(0, 1, 0);
                dir_3 = WireCell::Vector(0, 0, 1);
            } else {
                dir_2 = drift_dir_abs.cross(dir_1).norm();
                dir_3 = dir_1.cross(dir_2);
            }
            
            // Project associated points onto the local coordinate system
            std::vector<std::tuple<double, double, double>> vec_projs;
            for (const auto& pt : local_points_vec.at(i)) {
                double proj_1 = dir_1.dot(pt);
                double proj_2 = dir_2.dot(pt);
                double proj_3 = dir_3.dot(pt);
                vec_projs.push_back(std::make_tuple(proj_1, proj_2, proj_3));
            }
            
            // Calculate RMS spread in each direction
            int ncount = local_points_vec.at(i).size();
            if (ncount > 1) {
                WireCell::Point fit_pt(fits[i].point.x(), fits[i].point.y(), fits[i].point.z());
                std::tuple<double, double, double> means = std::make_tuple(
                    dir_1.dot(fit_pt),
                    dir_2.dot(fit_pt),
                    dir_3.dot(fit_pt)
                );
                
                for (const auto& proj : vec_projs) {
                    std::get<0>(vec_rms_vals.at(i)) += std::pow(std::get<0>(proj) - std::get<0>(means), 2);
                    std::get<1>(vec_rms_vals.at(i)) += std::pow(std::get<1>(proj) - std::get<1>(means), 2);
                    std::get<2>(vec_rms_vals.at(i)) += std::pow(std::get<2>(proj) - std::get<2>(means), 2);
                }
                
                std::get<0>(vec_rms_vals.at(i)) = std::sqrt(std::get<0>(vec_rms_vals.at(i)) / ncount);
                std::get<1>(vec_rms_vals.at(i)) = std::sqrt(std::get<1>(vec_rms_vals.at(i)) / ncount);
                std::get<2>(vec_rms_vals.at(i)) = std::sqrt(std::get<2>(vec_rms_vals.at(i)) / ncount);
            }
            
            // Calculate dQ/dx
            vec_dQ_dx.at(i) = fits[i].dQ / (fits[i].dx + 1e-9) / MIP_dQdx;
        }
        
        // Analyze spread characteristics
        double max_spread = 0;
        double large_spread_length = 0;
        double total_effective_length = 0;
        double total_length = 0;
        
        // bool flag_prev = false;
        for (size_t i = 0; i + 1 < local_points_vec.size(); i++) {
            double length = std::sqrt(
                std::pow(fits[i+1].point.x() - fits[i].point.x(), 2) +
                std::pow(fits[i+1].point.y() - fits[i].point.y(), 2) +
                std::pow(fits[i+1].point.z() - fits[i].point.z(), 2)
            );
            total_length += length;
            
            if (std::get<2>(vec_rms_vals.at(i)) != 0) {
                total_effective_length += length;
                if (std::get<2>(vec_rms_vals.at(i)) > rms_cut) {
                    large_spread_length += length;
                    // flag_prev = true;
                }
                if (std::get<2>(vec_rms_vals.at(i)) > max_spread) {
                    max_spread = std::get<2>(vec_rms_vals.at(i));
                }
            }
        }
        
        // Determine direction based on spread analysis
        int flag_dir = 0;
        
        // Check if this looks like a shower based on spread
        bool is_shower_like = (
            (max_spread > 0.7*units::cm && large_spread_length > 0.2 * total_effective_length && 
             total_effective_length > 3*units::cm && total_effective_length < 15*units::cm && 
             (large_spread_length > 2.7*units::cm || large_spread_length > 0.35 * total_effective_length)) ||
            (max_spread > 0.8*units::cm && large_spread_length > 0.3 * total_effective_length && 
             total_effective_length >= 15*units::cm) ||
            (max_spread > 1.0*units::cm && large_spread_length > 0.4 * total_effective_length)
        );
        
        if (is_shower_like) {
            WireCell::Vector main_dir1, main_dir2;
            bool flag_skip_angle1 = false;
            bool flag_skip_angle2 = false;
            
            // Create copies of points since segment_cal_dir_3vector expects non-const reference
            WireCell::Point front_pt = fits.front().point;
            WireCell::Point back_pt = fits.back().point;
            
            if (fits.front().point.z() < fits.back().point.z()) {
                main_dir1 = segment_cal_dir_3vector(segment, front_pt, 15*units::cm);
                main_dir2 = segment_cal_dir_3vector(segment, back_pt, 6*units::cm);
                double angle1 = std::acos(main_dir1.dot(drift_dir_abs)) * 180.0 / M_PI;
                if (std::fabs(angle1 - 90) < 10) flag_skip_angle1 = true;
            } else {
                main_dir1 = segment_cal_dir_3vector(segment, front_pt, 6*units::cm);
                main_dir2 = segment_cal_dir_3vector(segment, back_pt, 15*units::cm);
                double angle2 = std::acos(main_dir2.dot(drift_dir_abs)) * 180.0 / M_PI;
                if (std::fabs(angle2 - 90) < 10) flag_skip_angle2 = true;
            }
            
            // Group consecutive segments with large spread in forward direction
            // Each tuple: (start_index, end_index, max_rms_in_range)
            std::vector<std::tuple<int, int, double>> threshold_segs;
            
            for (size_t i = 0; i < vec_dQ_dx.size(); i++) {
                double angle = std::acos(main_dir1.dot(vec_dir.at(i))) * 180.0 / M_PI;
                if ((angle < 30 || (flag_skip_angle1 && angle < 60)) && 
                    (std::get<2>(vec_rms_vals.at(i))/units::cm > 0.4 || vec_dQ_dx.at(i) > 1.6)) {
                    
                    if (threshold_segs.empty()) {
                        // Start new segment group
                        threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                    } else {
                        // Check if continuous with previous group
                        if (i == std::get<1>(threshold_segs.back()) + 1) {
                            // Extend existing group
                            std::get<1>(threshold_segs.back()) = i;
                            if (std::get<2>(threshold_segs.back()) < std::get<2>(vec_rms_vals.at(i))) {
                                std::get<2>(threshold_segs.back()) = std::get<2>(vec_rms_vals.at(i));
                            }
                        } else {
                            // Start new group (gap detected)
                            threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                        }
                    }
                }
            }
            
            // Calculate total and max continuous length for forward direction
            double total_length1 = 0, max_length1 = 0;
            for (const auto& seg : threshold_segs) {
                int start_n = std::get<0>(seg);
                if (start_n > 0) start_n--;  // Include one segment before
                int end_n = std::get<1>(seg);
                
                double tmp_length = 0;
                for (int i = start_n; i < end_n && i + 1 < (int)fits.size(); i++) {
                    tmp_length += std::sqrt(
                        std::pow(fits[i+1].point.x() - fits[i].point.x(), 2) +
                        std::pow(fits[i+1].point.y() - fits[i].point.y(), 2) +
                        std::pow(fits[i+1].point.z() - fits[i].point.z(), 2)
                    );
                }
                total_length1 += tmp_length;
                if (tmp_length > max_length1) max_length1 = tmp_length;
            }
            
            // Group consecutive segments with large spread in backward direction
            threshold_segs.clear();
            
            for (int i = vec_dQ_dx.size() - 1; i >= 0; i--) {
                double angle = 180 - std::acos(main_dir2.dot(vec_dir.at(i))) * 180.0 / M_PI;
                if ((angle < 30 || (flag_skip_angle2 && angle < 60)) && 
                    (std::get<2>(vec_rms_vals.at(i))/units::cm > 0.4 || vec_dQ_dx.at(i) > 1.6)) {
                    
                    if (threshold_segs.empty()) {
                        // Start new segment group
                        threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                    } else {
                        // Check if continuous with previous group (decrementing)
                        if (i == std::get<1>(threshold_segs.back()) - 1) {
                            // Extend existing group
                            std::get<1>(threshold_segs.back()) = i;
                            if (std::get<2>(threshold_segs.back()) < std::get<2>(vec_rms_vals.at(i))) {
                                std::get<2>(threshold_segs.back()) = std::get<2>(vec_rms_vals.at(i));
                            }
                        } else {
                            // Start new group (gap detected)
                            threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                        }
                    }
                }
            }
            
            // Calculate total and max continuous length for backward direction
            double total_length2 = 0, max_length2 = 0;
            for (const auto& seg : threshold_segs) {
                int start_n = std::get<0>(seg);
                if (start_n < (int)fits.size() - 1) start_n++;  // Include one segment after
                int end_n = std::get<1>(seg);
                
                double tmp_length = 0;
                for (int i = start_n; i > end_n && i > 0; i--) {
                    tmp_length += std::sqrt(
                        std::pow(fits[i-1].point.x() - fits[i].point.x(), 2) +
                        std::pow(fits[i-1].point.y() - fits[i].point.y(), 2) +
                        std::pow(fits[i-1].point.z() - fits[i].point.z(), 2)
                    );
                }
                total_length2 += tmp_length;
                if (tmp_length > max_length2) max_length2 = tmp_length;
            }
            
            // Compare using both total and max continuous lengths
            if (total_length1 + max_length1 > 1.1 * (total_length2 + max_length2)) {
                flag_dir = 1;
            } else if (1.1 * (total_length1 + max_length1) < total_length2 + max_length2) {
                flag_dir = -1;
            }
        } else {
            // Not shower-like, use simpler direction determination
            if (total_length < 5*units::cm) {
                if (!segment_is_shower_trajectory(segment)) segment_determine_dir_track(segment, 0, fits.size(), particle_data, recomb_model);
                // For short segments, could call determine_dir_track here if needed
            } else {
                // Count consistent directions at each end
                WireCell::Point front_pt = fits.front().point;
                WireCell::Point back_pt = fits.back().point;
                WireCell::Vector main_dir_front = segment_cal_dir_3vector(segment, front_pt, 6*units::cm); 
                int ncount_front = 0;
                for (size_t i = 0; i < vec_dQ_dx.size(); i++) {
                    double angle = std::acos(main_dir_front.dot(vec_dir.at(i))) * 180.0 / M_PI;
                    if (angle < 30) ncount_front++;
                }

                WireCell::Vector main_dir_back = segment_cal_dir_3vector(segment, back_pt, 6*units::cm);
                int ncount_back = 0;
                for (int i = vec_dQ_dx.size() - 1; i >= 0; i--) {
                    double angle = 180 - std::acos(main_dir_back.dot(vec_dir.at(i))) * 180.0 / M_PI;
                    if (angle < 30) ncount_back++;
                }
                
                if (1.2 * ncount_front < ncount_back) {
                    flag_dir = -1;
                } else if (ncount_front > 1.2 * ncount_back) {
                    flag_dir = 1;
                }
            }
        }
        
        segment->dirsign(flag_dir);
        return (flag_dir != 0);
    }

    bool segment_is_shower_topology(SegmentPtr segment, bool tmp_val, double MIP_dQ_dx){
        int flag_dir = 0;
        bool flag_shower_topology = tmp_val; 
        const auto& fits = segment->fits();

        if (fits.empty()) return false;
        
        // Get the fit point cloud for KD-tree queries
        auto dpcloud_fit = segment->dpcloud("fit");
        if (!dpcloud_fit) return false;
        
        // Get the associated point cloud 
        auto dpcloud_assoc = segment->dpcloud("associated");
        if (!dpcloud_assoc) return false;
        
        const auto& assoc_points = dpcloud_assoc->get_points();
        if (assoc_points.empty()) return false;

        // Initialize vectors to store analysis results for each fit point
        std::vector<std::vector<WireCell::Point>> local_points_vec(fits.size());
        std::vector<std::tuple<double, double, double>> vec_rms_vals(fits.size(), std::make_tuple(0,0,0));
        std::vector<double> vec_dQ_dx(fits.size(), 0);
        
        // Build KD-tree index for fit points and associate points with nearest fit point
        auto& kd_tree_fit = dpcloud_fit->kd3d();
        
        for (const auto& pt : assoc_points) {
            WireCell::Point test_p(pt.x, pt.y, pt.z);
            auto results = kd_tree_fit.knn(1, test_p);
            if (!results.empty()) {
                size_t closest_fit_idx = results.front().first;
                local_points_vec.at(closest_fit_idx).push_back(test_p);
            }
        }
        
        WireCell::Vector drift_dir_abs(1, 0, 0);  // Drift direction
        
        // Calculate local directions and RMS spreads for each fit point
        for (size_t i = 0; i < local_points_vec.size(); i++) {
            // Calculate local direction from neighboring fit points
            WireCell::Vector v1(0, 0, 0);
            for (size_t j = 1; j < 3; j++) {
                if (i + j < fits.size()) {
                    v1 += WireCell::Vector(
                        fits[i+j].point.x() - fits[i].point.x(),
                        fits[i+j].point.y() - fits[i].point.y(),
                        fits[i+j].point.z() - fits[i].point.z()
                    );
                }
                if (i >= j) {
                    v1 += WireCell::Vector(
                        fits[i].point.x() - fits[i-j].point.x(),
                        fits[i].point.y() - fits[i-j].point.y(),
                        fits[i].point.z() - fits[i-j].point.z()
                    );
                }
            }
            
            WireCell::Vector dir_1 = v1.magnitude() > 0 ? v1.norm() : WireCell::Vector(1, 0, 0);
            
            // Set up orthogonal coordinate system
            WireCell::Vector dir_2, dir_3;
            double angle_deg = std::acos(dir_1.dot(drift_dir_abs)) * 180.0 / M_PI;
            
            if (angle_deg < 7.5) {
                dir_1 = WireCell::Vector(1, 0, 0);
                dir_2 = WireCell::Vector(0, 1, 0);
                dir_3 = WireCell::Vector(0, 0, 1);
            } else {
                dir_2 = drift_dir_abs.cross(dir_1).norm();
                dir_3 = dir_1.cross(dir_2);
            }
            
            // Project associated points onto the local coordinate system
            std::vector<std::tuple<double, double, double>> vec_projs;
            for (const auto& pt : local_points_vec.at(i)) {
                double proj_1 = dir_1.dot(pt);
                double proj_2 = dir_2.dot(pt);
                double proj_3 = dir_3.dot(pt);
                vec_projs.push_back(std::make_tuple(proj_1, proj_2, proj_3));
            }
            
            // Calculate RMS spread in each direction
            int ncount = local_points_vec.at(i).size();
            if (ncount > 1) {
                WireCell::Point fit_pt(fits[i].point.x(), fits[i].point.y(), fits[i].point.z());
                std::tuple<double, double, double> means = std::make_tuple(
                    dir_1.dot(fit_pt),
                    dir_2.dot(fit_pt),
                    dir_3.dot(fit_pt)
                );
                
                for (const auto& proj : vec_projs) {
                    std::get<0>(vec_rms_vals.at(i)) += std::pow(std::get<0>(proj) - std::get<0>(means), 2);
                    std::get<1>(vec_rms_vals.at(i)) += std::pow(std::get<1>(proj) - std::get<1>(means), 2);
                    std::get<2>(vec_rms_vals.at(i)) += std::pow(std::get<2>(proj) - std::get<2>(means), 2);
                }
                
                std::get<0>(vec_rms_vals.at(i)) = std::sqrt(std::get<0>(vec_rms_vals.at(i)) / ncount);
                std::get<1>(vec_rms_vals.at(i)) = std::sqrt(std::get<1>(vec_rms_vals.at(i)) / ncount);
                std::get<2>(vec_rms_vals.at(i)) = std::sqrt(std::get<2>(vec_rms_vals.at(i)) / ncount);
            }
            
            // Calculate dQ/dx
            vec_dQ_dx.at(i) = fits[i].dQ / (fits[i].dx + 1e-9) / MIP_dQ_dx;
        }
        
        // Analyze spread characteristics
        double max_spread = 0;
        double large_spread_length = 0;
        double total_effective_length = 0;
        
        double max_cont_length = 0;
        double max_cont_weighted_length = 0;
        double cont_length = 0;
        double cont_weighted_length = 0;
        bool flag_prev = false;
        
        for (size_t i = 0; i + 1 < local_points_vec.size(); i++) {
            double length = std::sqrt(
                std::pow(fits[i+1].point.x() - fits[i].point.x(), 2) +
                std::pow(fits[i+1].point.y() - fits[i].point.y(), 2) +
                std::pow(fits[i+1].point.z() - fits[i].point.z(), 2)
            );
            
            if (std::get<2>(vec_rms_vals.at(i)) != 0) {
                total_effective_length += length;
                if (std::get<2>(vec_rms_vals.at(i)) > 0.4 * units::cm) {
                    large_spread_length += length;
                    cont_length += length;
                    cont_weighted_length += length * std::get<2>(vec_rms_vals.at(i));
                    flag_prev = true;
                } else {
                    if (flag_prev && cont_length > max_cont_length) {
                        max_cont_length = cont_length;
                        max_cont_weighted_length = cont_weighted_length;
                    }
                    cont_length = 0;
                    cont_weighted_length = 0;
                    flag_prev = false;
                }
                if (std::get<2>(vec_rms_vals.at(i)) > max_spread) {
                    max_spread = std::get<2>(vec_rms_vals.at(i));
                }
            }
        }
        (void)max_cont_weighted_length; // Currently unused
        
        // Determine if this is shower topology based on spread patterns
        if ((max_spread > 0.7*units::cm && large_spread_length > 0.2 * total_effective_length && 
             total_effective_length > 3*units::cm && total_effective_length < 15*units::cm && 
             (large_spread_length > 2.7*units::cm || large_spread_length > 0.35 * total_effective_length)) ||
            (max_spread > 0.8*units::cm && large_spread_length > 0.3 * total_effective_length && 
             total_effective_length >= 15*units::cm) ||
            (max_spread > 0.8*units::cm && large_spread_length > 8*units::cm && 
             large_spread_length > 0.18 * total_effective_length) ||
            (max_spread > 1.0*units::cm && large_spread_length > 0.4 * total_effective_length) ||
            (max_spread > 1.0*units::cm && large_spread_length > 5*units::cm && 
             large_spread_length > 0.23 * total_effective_length)) {
            
            flag_shower_topology = true;
        }
        
        // If identified as shower topology, determine direction
        if (flag_shower_topology) {
            // Group consecutive segments with large spread in forward direction
            std::vector<std::tuple<int, int, double>> threshold_segs;
            
            for (size_t i = 0; i < vec_dQ_dx.size(); i++) {
                if (std::get<2>(vec_rms_vals.at(i))/units::cm > 0.4) {
                    if (threshold_segs.empty()) {
                        threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                    } else {
                        if (i == std::get<1>(threshold_segs.back()) + 1) {
                            // Extend existing group
                            std::get<1>(threshold_segs.back()) = i;
                            if (std::get<2>(threshold_segs.back()) < std::get<2>(vec_rms_vals.at(i))) {
                                std::get<2>(threshold_segs.back()) = std::get<2>(vec_rms_vals.at(i));
                            }
                        } else {
                            // Start new group
                            threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                        }
                    }
                }
            }
            
            // Calculate total and max continuous length for forward direction
            double total_length1 = 0, max_length1 = 0;
            for (const auto& seg : threshold_segs) {
                int start_n = std::get<0>(seg);
                if (start_n > 0) start_n--;
                int end_n = std::get<1>(seg);
                
                double tmp_length = 0;
                for (int i = start_n; i < end_n && i + 1 < (int)fits.size(); i++) {
                    tmp_length += std::sqrt(
                        std::pow(fits[i+1].point.x() - fits[i].point.x(), 2) +
                        std::pow(fits[i+1].point.y() - fits[i].point.y(), 2) +
                        std::pow(fits[i+1].point.z() - fits[i].point.z(), 2)
                    );
                }
                total_length1 += tmp_length;
                if (tmp_length > max_length1) max_length1 = tmp_length;
            }
            
            // Group consecutive segments with large spread in backward direction
            threshold_segs.clear();
            
            for (int i = vec_dQ_dx.size() - 1; i >= 0; i--) {
                if (std::get<2>(vec_rms_vals.at(i))/units::cm > 0.4) {
                    if (threshold_segs.empty()) {
                        threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                    } else {
                        if (i == std::get<1>(threshold_segs.back()) - 1) {
                            // Extend existing group
                            std::get<1>(threshold_segs.back()) = i;
                            if (std::get<2>(threshold_segs.back()) < std::get<2>(vec_rms_vals.at(i))) {
                                std::get<2>(threshold_segs.back()) = std::get<2>(vec_rms_vals.at(i));
                            }
                        } else {
                            // Start new group
                            threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                        }
                    }
                }
            }
            
            // Calculate total and max continuous length for backward direction
            double total_length2 = 0, max_length2 = 0;
            for (const auto& seg : threshold_segs) {
                int start_n = std::get<0>(seg);
                if (start_n < (int)fits.size() - 1) start_n++;
                int end_n = std::get<1>(seg);
                
                double tmp_length = 0;
                for (int i = start_n; i > end_n && i > 0; i--) {
                    tmp_length += std::sqrt(
                        std::pow(fits[i-1].point.x() - fits[i].point.x(), 2) +
                        std::pow(fits[i-1].point.y() - fits[i].point.y(), 2) +
                        std::pow(fits[i-1].point.z() - fits[i].point.z(), 2)
                    );
                }
                total_length2 += tmp_length;
                if (tmp_length > max_length2) max_length2 = tmp_length;
            }
            
            // Determine direction based on spread comparison
            if (total_length1 + max_length1 > 1.1 * (total_length2 + max_length2)) {
                flag_dir = 1;
            } else if (1.1 * (total_length1 + max_length1) < total_length2 + max_length2) {
                flag_dir = -1;
            }
            
            // Override shower topology for very long segments with little spread
            double tmp_total_length = segment_track_length(segment, 0);
            if (tmp_total_length > 50*units::cm && 
                total_length1 < 0.25 * tmp_total_length && 
                total_length2 < 0.25 * tmp_total_length) {
                flag_dir = 0;
                flag_shower_topology = false;
            }
        }

        if (flag_shower_topology) segment->set_flags(SegmentFlags::kShowerTopology);
        segment->dirsign(flag_dir);
        return flag_shower_topology;
    }

}
