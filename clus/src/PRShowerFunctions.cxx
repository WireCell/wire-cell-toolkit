#include "WireCellClus/PRShowerFunctions.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/DynamicPointCloud.h"

using namespace WireCell::Clus::PR;

namespace WireCell::Clus::PR {

    std::pair<double, WireCell::Point> shower_get_closest_point(Shower& shower, const WireCell::Point& point, const std::string& cloud_name /* = "fit" */){
        // Get the dynamic point cloud from the shower
        auto pcloud = shower.get_pcloud(cloud_name);
        
        // If no point cloud exists, return invalid result
        if (!pcloud) {
            return std::make_pair(-1.0, WireCell::Point(0, 0, 0));
        }
        
        // Get the 3D KD-tree
        auto& kd3d = pcloud->kd3d();
        
        // Prepare query point
        std::vector<double> query = {point.x(), point.y(), point.z()};
        
        // Find the nearest neighbor
        auto results = kd3d.knn(1, query);
        
        // Check if we found a point
        if (results.empty()) {
            return std::make_pair(-1.0, WireCell::Point(0, 0, 0));
        }
        
        // Get the result
        const size_t idx = results[0].first;
        const double distance = sqrt(results[0].second);  // KD-tree returns squared distance
        
        // Get the actual point from the point cloud
        const auto& points = pcloud->get_points();
        const auto& closest_pt = points[idx];
        
        return std::make_pair(distance, WireCell::Point(closest_pt.x, closest_pt.y, closest_pt.z));
    }

    double shower_get_closest_dis(Shower& shower, SegmentPtr seg, const std::string& cloud_name /* = "fit" */){
        // Get the first point from segment's DynamicPointCloud
        auto seg_dpc = seg->dpcloud(cloud_name);
        if (!seg_dpc) {
            return -1.0;
        }
        
        const auto& seg_points = seg_dpc->get_points();
        if (seg_points.empty()) {
            return -1.0;
        }
        
        // Use the first point from segment's point cloud
        WireCell::Point first_point(seg_points.front().x, seg_points.front().y, seg_points.front().z);
        
        // Step 1: Get closest point in shower to the first point of the segment
        WireCell::Point test_p = shower_get_closest_point(shower, first_point, cloud_name).second;
        
        // Step 2: Get closest point in segment to that point (uses segment's DPC with KD-tree)
        test_p = segment_get_closest_point(seg, test_p, cloud_name).second;
        
        // Step 3: Get closest point in shower to that point
        test_p = shower_get_closest_point(shower, test_p, cloud_name).second;
        
        // Step 4: Get closest distance in segment to that point
        return segment_get_closest_point(seg, test_p, cloud_name).first;
    }

    double shower_get_dis(Shower& shower, SegmentPtr seg, const std::string& cloud_name /* = "fit" */){
        double min_dis = 1e9;
        WireCell::Point min_point;
        
        // Get the first point from the input segment's DynamicPointCloud
        auto seg_dpc = seg->dpcloud(cloud_name);
        if (!seg_dpc) {
            return -1.0;
        }
        
        const auto& seg_points = seg_dpc->get_points();
        if (seg_points.empty()) {
            return -1.0;
        }
        
        WireCell::Point test_p(seg_points.front().x, seg_points.front().y, seg_points.front().z);
        
        // Get the view graph to access segments
        const auto& view = shower.view_graph();
        
        // First iteration: find the closest point in any shower segment to test_p
        for (auto edesc : shower.edges()) {
            SegmentPtr sg = view[edesc].segment;
            if (!sg) continue;
            
            auto results = segment_get_closest_point(sg, test_p, cloud_name);
            if (results.first < min_dis) {
                min_dis = results.first;
                min_point = results.second;
            }
        }
        
        // Get closest point in input segment to that minimum point
        auto results1 = segment_get_closest_point(seg, min_point, cloud_name);
        test_p = results1.second;
        
        // Second iteration: find the closest distance from any shower segment to the new test point
        for (auto edesc : shower.edges()) {
            SegmentPtr sg = view[edesc].segment;
            if (!sg) continue;
            
            auto results = segment_get_closest_point(sg, test_p, cloud_name);
            if (results.first < min_dis) {
                min_dis = results.first;
                min_point = results.second;
            }
        }
        
        return min_dis;
    }

    WireCell::Vector shower_cal_dir_3vector(Shower& shower, const WireCell::Point& p, double dis_cut /* = 15*units::cm */){
        WireCell::Point p_sum(0, 0, 0);
        int ncount = 0;
        
        // Get the view graph to access segments
        const auto& view = shower.view_graph();
        
        // Loop through all segments in the shower
        for (auto edesc : shower.edges()) {
            SegmentPtr seg = view[edesc].segment;
            if (!seg) continue;
            
            // Get the segment's fits
            const auto& fits = seg->fits();
            
            // Check each fit point in the segment
            for (const auto& fit : fits) {
                // Calculate distance from this fit point to p
                double dis = sqrt(pow(fit.point.x() - p.x(), 2) + 
                                 pow(fit.point.y() - p.y(), 2) + 
                                 pow(fit.point.z() - p.z(), 2));
                
                // If within distance cutoff, accumulate the point
                if (dis < dis_cut) {
                    p_sum = WireCell::Point(p_sum.x() + fit.point.x(), 
                                           p_sum.y() + fit.point.y(), 
                                           p_sum.z() + fit.point.z());
                    ncount++;
                }
            }
        }
        
        // If no points were found, return zero vector
        if (ncount == 0) {
            return WireCell::Vector(0, 0, 0);
        }
        
        // Calculate the average point and direction vector
        WireCell::Point p_avg(p_sum.x() / ncount, 
                             p_sum.y() / ncount, 
                             p_sum.z() / ncount);
        
        // Direction vector from p to average point
        WireCell::Vector dir(p_avg.x() - p.x(), 
                            p_avg.y() - p.y(), 
                            p_avg.z() - p.z());
        
        // Normalize the vector
        double magnitude = dir.magnitude();
        if (magnitude > 0) {
            dir = dir.norm();
        }
        
        return dir;
    }



} // namespace WireCell::Clus::PR
