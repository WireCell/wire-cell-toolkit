#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/DynamicPointCloud.h"
#include "WireCellClus/ClusteringFuncs.h"

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



}
