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
        
        


    double segment_median_dQ_dx(SegmentPtr seg)
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
        
        // Use nth_element to find median (same algorithm as original)
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
            for (int i = start; i < start + num_points - 1 && (fits.size() - i - 1) < fits.size(); i++) {
                if (fits.size() - start < fits.size()) {
                    p = p + (fits[fits.size() - i - 1].point - fits[fits.size() - start].point);
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
            proton_ref[i] = particle_data->get_dEdx_function("proton")->scalar_function((vec_x[i])/units::cm)/ units::cm;
            electron_ref[i] = particle_data->get_dEdx_function("electron")->scalar_function((vec_x[i])/units::cm)/ units::cm;
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
        
        std::vector<double> results;
        // Convert bool result to double (1.0 for true, 0.0 for false)
        results.push_back(eval_ks_ratio(ks1, ks2, ratio1, ratio2) ? 1.0 : 0.0); // direction metric
        results.push_back(sqrt(pow(ks1, 2) + pow(ratio1-1, 2))); // muon information
        results.push_back(sqrt(pow(ks3, 2) + pow(ratio3-1, 2))); // proton information  
        results.push_back(sqrt(pow(ks4, 2) + pow(ratio4-1, 2))); // electron information
        
        return results;
    }
   
    double cal_kine_range(double L, int particle_type, const Clus::ParticleDataSet::pointer& particle_data){

        IScalarFunction::pointer range_function = nullptr;
        
        if (abs(particle_type) == 11) {        // electron
            range_function = particle_data->get_range_function("electron");
        }
        else if (abs(particle_type) == 13) {   // muon
            range_function = particle_data->get_range_function("muon");
        }
        else if (abs(particle_type) == 211) {  // pion
            range_function = particle_data->get_range_function("pion");
        }
        else if (abs(particle_type) == 321) {  // kaon
            range_function = particle_data->get_range_function("kaon");
        }
        else if (abs(particle_type) == 2212) { // proton
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

    // 4-momentum: px, py, pz, E and the kine_energy ...
    std::vector<double> segment_cal_4mom(SegmentPtr segment, int particle_type, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx){
        double length = segment_track_length(segment, 0);
        double kine_energy = 0;

        std::vector<double> results(5,0.0); // 4-momentum: px, py, pz, E and the kine_energy ...

        if (length < 4*units::cm){
            kine_energy = segment_cal_kine_dQdx(segment, recomb_model); // short track 
        }else if (segment->flags_any(PR::SegmentFlags::kShowerTrajectory)){
            kine_energy = segment_cal_kine_dQdx(segment, recomb_model);
        }else{
            kine_energy = cal_kine_range(length, particle_type, particle_data);
        }
        results[4] = kine_energy;

        double particle_mass = particle_data->get_particle_mass(particle_type);

        results[3]= kine_energy + particle_mass;
        double mom = sqrt(pow(results[3],2) - pow(particle_mass,2));
        auto v1 = segment_cal_dir_3vector(segment);
        results[0] = mom * v1.x();
        results[1] = mom * v1.y();
        results[2] = mom * v1.z();

        return results;
    }

}
