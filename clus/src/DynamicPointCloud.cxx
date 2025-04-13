#include "WireCellClus/DynamicPointCloud.h"
#include "WireCellClus/Facade.h"
#include "WireCellClus/Facade_Util.h"
#include "WireCellUtil/Logging.h"

#include <boost/histogram.hpp>
#include <boost/histogram/algorithm/sum.hpp>

using namespace WireCell;
using namespace WireCell::PointCloud::Facade;
using spdlog::debug;

#ifdef __DEBUG__
#define LogDebug(x) std::cout << "[DPC]: " << __LINE__ << " : " << x << std::endl
#else
#define LogDebug(x)
#endif

DynamicPointCloud::nfkd_t &DynamicPointCloud::kd3d() const
{
    if (!m_kd3d) {
        m_kd3d = std::make_unique<nfkd_t>(3);
    }
    return *m_kd3d;
}

DynamicPointCloud::nfkd_t &DynamicPointCloud::kd2d(const int plane, const int face, const int apa) const
{
    WirePlaneId wpid(iplane2layer[plane], face, apa);
    // SPDLOG_DEBUG("DynamicPointCloud: kd2d {} {} {} wpid {}", plane, face, apa, wpid.name());
    auto iter = m_kd2d.find(wpid.ident());
    if (iter == m_kd2d.end()) {
        m_kd2d[wpid.ident()] = std::make_unique<nfkd_t>(2);
    }
    return *m_kd2d[wpid.ident()];
}

const std::unordered_map<size_t, size_t> &DynamicPointCloud::kd2d_l2g(const int plane, const int face,
                                                                      const int apa) const
{
    WirePlaneId wpid(iplane2layer[plane], face, apa);
    auto iter = m_kd2d_index_l2g.find(wpid.ident());

    // Debug logging for the requested wpid information
    // std::cout << "DynamicPointCloud: kd2d_l2g requested for wpid " << wpid.name() 
    //           << " (ident: " << wpid.ident() << ", face: " << wpid.face() 
    //           << ", apa: " << wpid.apa() << ")" << std::endl;
    // for (auto it = m_kd2d_index_l2g.begin(); it != m_kd2d_index_l2g.end(); ++it) {
    //     std::cout << "DynamicPointCloud: kd2d_l2g available wpid " << WirePlaneId(it->first).name()
    //           << " (ident: " << it->first << ", face: " << WirePlaneId(it->first).face()
    //           << ", apa: " << WirePlaneId(it->first).apa() << ")" << std::endl;
    // }

    if (iter == m_kd2d_index_l2g.end()) {
        // Create empty mapping for this wpid instead of raising an error
        m_kd2d_index_l2g[wpid.ident()] = std::unordered_map<size_t, size_t>();
        iter = m_kd2d_index_l2g.find(wpid.ident());
        SPDLOG_DEBUG("DynamicPointCloud: created empty 2D index l2g for wpid {}", wpid.name());
    }
    return iter->second;
}

const std::unordered_map<size_t, std::vector<size_t> > &DynamicPointCloud::kd2d_g2l(const int plane, const int face,
                                                                      const int apa) const
{
    WirePlaneId wpid(iplane2layer[plane], face, apa);
    auto iter = m_kd2d_index_g2l.find(wpid.ident());
    if (iter == m_kd2d_index_g2l.end()) {
        raise<RuntimeError>("DynamicPointCloud: missing 2D index g2l for wpid %s", wpid.name());
    }
    return iter->second;
}


void DynamicPointCloud::add_points(const std::vector<DPCPoint> &points) {
    if (points.empty()) {
        return;
    }
    
    // Preallocate memory and get original size
    size_t original_size = m_points.size();
    m_points.reserve(original_size + points.size());
    
    // Move the points instead of copying
    m_points.insert(m_points.end(),
                   std::make_move_iterator(points.begin()),
                   std::make_move_iterator(points.end()));
    
    // Process 3D KD tree
    auto &kd3d = this->kd3d();
    
    // Prepare batch data for 3D KD tree
    NFKDVec::Tree<double>::points_type pts3d(3);
    for (size_t i = 0; i < 3; ++i) {
        pts3d[i].reserve(points.size());
    }
    
    // Prepare maps to store 2D points for each plane and track local-to-global mappings
    std::map<int, NFKDVec::Tree<double>::points_type> planes_pts;
    std::map<int, std::vector<size_t>> planes_global_indices;
    
    // Extract data and prepare for batch processing
    for (size_t i = 0; i < points.size(); ++i) {
        size_t global_idx = original_size + i;
        const auto &pt = points[i];
        
        // Add to 3D points
        pts3d[0].push_back(pt.x);
        pts3d[1].push_back(pt.y);
        pts3d[2].push_back(pt.z);
        
        // Check 2D projection validity
        if (pt.x_2d.size() != 3 || pt.y_2d.size() != 3) {
            raise<RuntimeError>("DynamicPointCloud: unexpected 2D projection size x_2d %d y_2d %d", 
                               pt.x_2d.size(), pt.y_2d.size());
        }
        
        // Skip 2D KD if wpid is not valid
        WirePlaneId wpid_volume(pt.wpid);
        if (wpid_volume.face() == -1 || wpid_volume.apa() == -1) {
            continue;
        }
        
        // Process 2D points for each plane
        for (size_t pindex = 0; pindex < 3; ++pindex) {

            // std::cout << "Test: " << pindex << " " << pt.x_2d[pindex].size() << " " << pt.y_2d[pindex].size() << " " << pt.wpid_2d[pindex].size() << std::endl;

            // Add 2D point to plane data
            for (size_t j = 0; j < pt.x_2d[pindex].size(); ++j) {
                WirePlaneId wpid_2d(pt.wpid_2d[pindex].at(j));
                WirePlaneId wpid_plane(iplane2layer[pindex], wpid_2d.face(), wpid_2d.apa());
                int key = wpid_plane.ident();
                // Initialize plane data structures if not exists
                if (planes_pts.find(key) == planes_pts.end()) {
                    planes_pts[key] = NFKDVec::Tree<double>::points_type(2);
                    planes_pts[key][0].reserve(points.size());
                    planes_pts[key][1].reserve(points.size());
                    planes_global_indices[key].reserve(points.size());
                }
                planes_pts[key][0].push_back(pt.x_2d[pindex][j]);
                planes_pts[key][1].push_back(pt.y_2d[pindex][j]);
                planes_global_indices[key].push_back(global_idx);
            }
        }

    }
    
    // Batch append 3D points
    kd3d.append(pts3d);
    
    
    // Create a reverse mapping from layer to iplane based on the existing iplane2layer array
    std::unordered_map<WirePlaneLayer_t, int> layer2iplane;
    for (int i = 0; i < 3; ++i) {
        layer2iplane[iplane2layer[i]] = i;
    }


    // Batch append 2D points for each plane
    for (const auto& [key, pts] : planes_pts) {
        WirePlaneId wpid_plane(key);
        int pindex = layer2iplane[wpid_plane.layer()];
        auto& kd2d = this->kd2d(pindex, wpid_plane.face(), wpid_plane.apa());
        
        // Get the starting index for the new points
        size_t start_idx = kd2d.npoints();
        
        // Batch append 2D points
        kd2d.append(pts);
        
        // Update index mappings
        const auto& indices = planes_global_indices[key];
        for (size_t i = 0; i < indices.size(); ++i) {
            size_t local_idx = start_idx + i;
            size_t global_idx = indices[i];
            m_kd2d_index_l2g[key][local_idx] = global_idx;
            m_kd2d_index_g2l[key][global_idx].push_back(local_idx); // save things to a vector
        }
    }
}




std::vector<std::tuple<double, const Cluster *, size_t>>
DynamicPointCloud::get_2d_points_info(const geo_point_t &p, const double radius, const int plane, const int face,
                                     const int apa) const
{
    // Create WirePlaneId once
    WirePlaneId wpid_volume(kAllLayers, face, apa);
    
    // Get KD tree and mapping
    auto &kd2d = this->kd2d(plane, face, apa);
    auto &l2g = this->kd2d_l2g(plane, face, apa);
    
    // Get angle parameters - lookup once
    if (m_wpid_params.find(wpid_volume) == m_wpid_params.end()) {
        raise<RuntimeError>("DynamicPointCloud: missing wpid params for wpid %s", wpid_volume.name());
    }
    const auto [_, angle_u, angle_v, angle_w] = m_wpid_params.at(wpid_volume);
    
    // Compute projected point
    const double angle = (plane == 0) ? angle_u : ((plane == 1) ? angle_v : angle_w);
    const double projected_y = cos(angle) * p.z() - sin(angle) * p.y();
    
    // Prepare query point
    std::vector<double> query = {p.x(), projected_y};

    // Perform radius search
    auto results = kd2d.radius(radius * radius, query);
    
    // Optimize for empty results case
    if (results.empty()) {
        return {};
    }
    
    // Pre-allocate return vector
    std::vector<std::tuple<double, const Cluster *, size_t>> return_results;
    return_results.reserve(results.size());
    
    // Process results
    for (const auto &[local_idx, dist_squared] : results) {
        const size_t global_idx = l2g.at(local_idx);
        const auto &pt = m_points[global_idx];
        
        // Option 1: Return squared distances if caller doesn't need actual distances
        // return_results.emplace_back(dist_squared, pt.cluster, global_idx);
        
        // Option 2: Compute sqrt only when needed
        return_results.emplace_back(sqrt(dist_squared), pt.cluster, global_idx);
    }

    return return_results;
}



std::tuple<double, const Cluster *, size_t>
DynamicPointCloud::get_closest_2d_point_info(const geo_point_t &p, const int plane, const int face, const int apa) const
{
    // Create WirePlaneId only once
    const WirePlaneId wpid_volume(kAllLayers, face, apa);
    
    // Get KD tree and mapping - only need l2g here, g2l isn't used
    auto &kd2d = this->kd2d(plane, face, apa);
    auto &l2g = this->kd2d_l2g(plane, face, apa);
    
    // Check and get angle parameters
    auto wpid_iter = m_wpid_params.find(wpid_volume);
    if (wpid_iter == m_wpid_params.end()) {
        raise<RuntimeError>("DynamicPointCloud: missing wpid params for wpid %s", wpid_volume.name());
    }
    
    // Calculate angle more directly based on plane parameter
    const auto &[_, angle_u, angle_v, angle_w] = wpid_iter->second;
    const double angle = (plane == 0) ? angle_u : ((plane == 1) ? angle_v : angle_w);
    
    // Prepare query point more efficiently
    const double projected_y = cos(angle) * p.z() - sin(angle) * p.y();
    const std::vector<double> query = {p.x(), projected_y};

    // Perform nearest neighbor search
    auto results = kd2d.knn(1, query);

    // Early return for empty results
    if (results.empty()) {
        return std::make_tuple(-1.0, nullptr, static_cast<size_t>(-1));
    }
    
    // Process the single result
    const size_t local_idx = results[0].first;
    const double distance = sqrt(results[0].second);  // Only compute sqrt once
    const size_t global_idx = l2g.at(local_idx);
    const auto &pt = m_points[global_idx];
    
    return std::make_tuple(distance, pt.cluster, global_idx);
}

std::pair<double, double> DynamicPointCloud::hough_transform(const geo_point_t &origin, const double dis) const
{
    auto &kd3d = this->kd3d();

    // Preallocate with reasonable capacity to reduce reallocations
    std::vector<geo_point_t> pts;
    std::vector<const Blob *> blobs;
    pts.reserve(100);  // Reasonable starting size, adjust based on typical usage
    blobs.reserve(100);

    // Create query point
    std::vector<double> query = {origin.x(), origin.y(), origin.z()};

    // Perform radius search
    auto results = kd3d.radius(dis * dis, query);

    // Preallocate exact size now that we know it
    pts.reserve(results.size());
    blobs.reserve(results.size());

    for (const auto &[idx, _] : results) {
        const auto &pt = m_points[idx];
        pts.push_back({pt.x, pt.y, pt.z});
        blobs.push_back(pt.blob);
    }

    namespace bh = boost::histogram;
    namespace bha = boost::histogram::algorithm;

    constexpr double pi = 3.141592653589793;
    const double ten_cm = 10.0 * units::cm;

    // Parameter axis 1 is theta angle
    const int nbins1 = 180;
    auto theta_param = [](const Vector &dir) {
        const Vector Z(0, 0, 1);
        return acos(Z.dot(dir));
    };
    double min1 = 0, max1 = pi;

    // Parameter axis 2 is phi angle
    const int nbins2 = 360;
    const double min2 = -pi;
    const double max2 = +pi;
    auto phi_param = [](const Vector &dir) {
        const Vector X(1, 0, 0);
        const Vector Y(0, 1, 0);
        return atan2(Y.dot(dir), X.dot(dir));
    };

    auto hist = bh::make_histogram(bh::axis::regular<>(nbins1, min1, max1), bh::axis::regular<>(nbins2, min2, max2));

    // Early return if no points found
    if (pts.empty()) {
        return {0.0, 0.0};
    }

    for (size_t ind = 0; ind < blobs.size(); ++ind) {
        const auto *blob = blobs[ind];
        auto charge = blob ? blob->charge() : 1.0;

        if (charge <= 0) continue;

        const auto npoints = blob ? blob->npoints() : 1;
        const auto &pt = pts[ind];

        // Use the original normalization method
        const Vector dir = (pt - origin).norm();
        const double r = (pt - origin).magnitude();

        const double p1 = theta_param(dir);
        const double p2 = phi_param(dir);

        if (r < ten_cm) {
            hist(p1, p2, bh::weight(charge / npoints));
        }
        else {
            // Use original formula
            hist(p1, p2, bh::weight(charge / npoints * pow(ten_cm / r, 2)));
        }
    }

    // Use the original max finding approach
    auto indexed = bh::indexed(hist);
    auto it = std::max_element(indexed.begin(), indexed.end());
    const auto &cell = *it;
    return {cell.bin(0).center(), cell.bin(1).center()};
}




geo_point_t DynamicPointCloud::vhough_transform(const geo_point_t &origin, const double dis) const
{
    const auto [th, phi] = hough_transform(origin, dis);
    return {sin(th) * cos(phi), sin(th) * sin(phi), cos(th)};
}



std::vector<DynamicPointCloud::DPCPoint> PointCloud::Facade::make_points_cluster(
    const Cluster *cluster, const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params, bool flag_wrap)
{
    if (!cluster) {
        SPDLOG_WARN("make_points_cluster: null cluster return empty points");
        return {};
    }

    const size_t num_points = cluster->npoints();
    std::vector<DynamicPointCloud::DPCPoint> dpc_points;
    dpc_points.reserve(num_points);

    const auto &winds = cluster->wire_indices();
    
    // Cache commonly referenced WPIDs and their params to avoid map lookups
    std::unordered_map<int, std::tuple<geo_point_t, double, double, double>> cached_params;
    
    for (size_t ipt = 0; ipt < num_points; ++ipt) {
        geo_point_t pt = cluster->point3d(ipt);
        const auto wpid = cluster->wire_plane_id(ipt);
        // std::cout << " DEBUG PointCloud::Facade::make_points_cluster wpid " << wpid.name() << std::endl;
        int wpid_ident = wpid.ident();
        
        // Check cache first, then populate if needed
        auto param_it = cached_params.find(wpid_ident);
        if (param_it == cached_params.end()) {
            auto wpid_it = wpid_params.find(wpid);
            if (wpid_it == wpid_params.end()) {
                raise<RuntimeError>("make_points_cluster: missing wpid params for wpid %s", wpid.name());
            }
            param_it = cached_params.emplace(wpid_ident, wpid_it->second).first;
        }
        
        const auto &[drift_dir, angle_u, angle_v, angle_w] = param_it->second;
        const double angle_uvw[3] = {angle_u, angle_v, angle_w};

        DynamicPointCloud::DPCPoint point;
        point.x = pt.x();
        point.y = pt.y();
        point.z = pt.z();
        point.wpid = wpid.ident();
        point.cluster = cluster;
        point.blob = cluster->blob_with_point(ipt);
        
        // Pre-allocate vectors with correct size
        point.x_2d.resize(3);
        point.y_2d.resize(3);
        point.wpid_2d.resize(3);
        
        if (flag_wrap){
            fill_wrap_points(cluster, pt, WirePlaneId(wpid), point.x_2d, point.y_2d, point.wpid_2d);
        }else{
            for (size_t pindex = 0; pindex < 3; ++pindex) {
                point.x_2d[pindex].push_back(point.x);
                point.y_2d[pindex].push_back(cos(angle_uvw[pindex]) * point.z - sin(angle_uvw[pindex]) * point.y);
                point.wpid_2d[pindex].push_back(wpid.ident());
            }
        }

        
        point.wind = {winds[0][ipt], winds[1][ipt], winds[2][ipt]};
        point.dist_cut = {-1e12, -1e12, -1e12};

        dpc_points.push_back(std::move(point));
    }

    return dpc_points;
}



std::vector<DynamicPointCloud::DPCPoint>
PointCloud::Facade::make_points_cluster_skeleton(const Cluster *cluster, const IDetectorVolumes::pointer dv,
                             const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params, bool flag_wrap,
                             const double step)
{
    std::vector<DynamicPointCloud::DPCPoint> dpc_points;

    if (!cluster) {
        SPDLOG_WARN("make_points_cluster_skeleton: null cluster return empty points");
        return dpc_points;
    }

    // Estimate capacity to avoid reallocations
    const std::list<size_t> &path_wcps = cluster->get_path_wcps();
    size_t estimated_capacity = path_wcps.size() * 2; // Rough estimate
    dpc_points.reserve(estimated_capacity);

    // Cache for angle values per wpid to avoid repeated tuple unpacking
    std::unordered_map<int, std::array<double, 3>> wpid_angles_cache;
    
    geo_point_t prev_wcp = cluster->point3d(path_wcps.front());
    auto prev_wpid = cluster->wire_plane_id(path_wcps.front());

    // Pre-computed constants
    const double dist_cut_value = 2.4 * units::cm;
    
    for (auto it = path_wcps.begin(); it != path_wcps.end(); it++) {
        geo_point_t test_point = cluster->point3d(*it);
        double dis = (test_point - prev_wcp).magnitude();
        auto wpid_test_point = cluster->wire_plane_id(*it);

        if (wpid_params.find(wpid_test_point) == wpid_params.end()) {
            raise<RuntimeError>("make_points_cluster: missing wpid params for wpid %s", wpid_test_point.name());
        }
        
        // Get or compute angle values for this wpid
        std::array<double, 3> angle_uvw;
        auto cache_it = wpid_angles_cache.find(wpid_test_point.ident());
        if (cache_it == wpid_angles_cache.end()) {
            const auto& [drift_dir, angle_u, angle_v, angle_w] = wpid_params.at(wpid_test_point);
            angle_uvw = {angle_u, angle_v, angle_w};
            wpid_angles_cache[wpid_test_point.ident()] = angle_uvw;
        } else {
            angle_uvw = cache_it->second;
        }

        if (dis <= step) {
            DynamicPointCloud::DPCPoint point;
            point.wpid = wpid_test_point.ident(); // Direct assignment without recreation
            point.cluster = cluster;
            point.blob = nullptr;
            
            // Pre-allocate and initialize vectors
            point.x_2d.resize(3);
            point.y_2d.resize(3);
            point.wpid_2d.resize(3);
            point.wind = {-1e12, -1e12, -1e12};
            point.dist_cut = {dist_cut_value, dist_cut_value, dist_cut_value};
            
            point.x = test_point.x();
            point.y = test_point.y();
            point.z = test_point.z();

            if (flag_wrap){
                fill_wrap_points(cluster, test_point, WirePlaneId(wpid_test_point), point.x_2d, point.y_2d, point.wpid_2d);
            }else{
                for (size_t pindex = 0; pindex < 3; ++pindex) {
                    point.x_2d[pindex].push_back(test_point.x());
                    point.y_2d[pindex].push_back(cos(angle_uvw[pindex]) * test_point.z() - sin(angle_uvw[pindex]) * test_point.y());
                    point.wpid_2d[pindex].push_back(wpid_test_point.ident());
                }
            }
            
            dpc_points.push_back(std::move(point));
        }
        else {
            int num_points = int(dis / step) + 1;
            
            // Pre-compute direction vectors to avoid recalculation in loop
            double dx = (test_point.x() - prev_wcp.x()) / num_points;
            double dy = (test_point.y() - prev_wcp.y()) / num_points;
            double dz = (test_point.z() - prev_wcp.z()) / num_points;
            
            for (int k = 0; k != num_points; k++) {
                DynamicPointCloud::DPCPoint point;
                
                // Faster interpolation with pre-computed increments
                double t = (k + 1.0);
                point.x = prev_wcp.x() + t * dx;
                point.y = prev_wcp.y() + t * dy;
                point.z = prev_wcp.z() + t * dz;
                
                geo_point_t temp_point(point.x, point.y, point.z);
                auto temp_wpid = WirePlaneId(get_wireplaneid(temp_point, prev_wpid, wpid_test_point, 
                                               cluster->grouping()->get_detector_volumes()));
                
                point.wpid = temp_wpid.ident(); // Direct assignment
                point.cluster = cluster;
                point.blob = nullptr;
                point.wind = {-1e12, -1e12, -1e12};
                point.dist_cut = {dist_cut_value, dist_cut_value, dist_cut_value};
                point.x_2d.resize(3);
                point.y_2d.resize(3);
                point.wpid_2d.resize(3);

                
                if (temp_wpid.apa() != -1) {
                    // Get cached angles if available
                    std::array<double, 3> temp_angle_uvw;
                    auto cache_it = wpid_angles_cache.find(temp_wpid.ident());
                    if (cache_it == wpid_angles_cache.end()) {
                        const auto& [drift_dir, angle_u, angle_v, angle_w] = wpid_params.at(temp_wpid);
                        temp_angle_uvw = {angle_u, angle_v, angle_w};
                        wpid_angles_cache[temp_wpid.ident()] = temp_angle_uvw;
                    } else {
                        temp_angle_uvw = cache_it->second;
                    }
                    
                    if (flag_wrap){
                        fill_wrap_points(cluster, temp_point, temp_wpid,  point.x_2d, point.y_2d, point.wpid_2d);
                    }else{
                        for (size_t pindex = 0; pindex < 3; ++pindex) {
                            point.x_2d[pindex].push_back(point.x);
                            point.y_2d[pindex].push_back(cos(temp_angle_uvw[pindex]) * point.z - 
                                                sin(temp_angle_uvw[pindex]) * point.y);
                            point.wpid_2d[pindex].push_back(temp_wpid.ident());

                        }
                    }
                }
                // } else {
                //     // point.x_2d = {-1e12, -1e12, -1e12};
                //     // point.y_2d = {-1e12, -1e12, -1e12};
                // }
                
                dpc_points.push_back(std::move(point));
            }
        }

        prev_wcp = test_point;
        prev_wpid = wpid_test_point;
    }
    
    return dpc_points;
}



std::vector<DynamicPointCloud::DPCPoint> PointCloud::Facade::make_points_linear_extrapolation(
    const Cluster *cluster, const geo_point_t &p_test, const geo_point_t &dir_unmorm, const double range,
    const double step, const double angle, const IDetectorVolumes::pointer dv,
    const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params)
{
    std::vector<DynamicPointCloud::DPCPoint> dpc_points;

    if (!cluster) {
        SPDLOG_WARN("make_points_linear_extrapolation: null cluster return empty points");
        return dpc_points;
    }

    // Get wpid once from the grouping
    const auto wpid = *(cluster->grouping()->wpids().begin());
    
    // Check wpid early
    if (wpid_params.find(wpid) == wpid_params.end()) {
        raise<RuntimeError>("make_points_cluster: missing wpid params for wpid %s", wpid.name());
    }
    
    // Pre-compute constants
    const double DEG_TO_RAD = 3.1415926/180.0;
    const double sin_angle_rad = sin(angle * DEG_TO_RAD);
    const double MIN_DIS_CUT = 2.4 * units::cm;
    const double MAX_DIS_CUT = 13.0 * units::cm;
    
    // Normalize direction once
    geo_point_t dir = dir_unmorm.norm();

    // Calculate segment count and distances
    int num_points = int(range / step) + 1;
    double dis_seg = range / num_points;
    
    // Pre-compute direction scaling
    double dx_step = dir.x() * dis_seg;
    double dy_step = dir.y() * dis_seg;
    double dz_step = dir.z() * dis_seg;
    
    // Get angle values once
    const auto [drift_dir, angle_u, angle_v, angle_w] = wpid_params.at(wpid);
    const double cos_angle_uvw[3] = {cos(angle_u), cos(angle_v), cos(angle_w)};
    const double sin_angle_uvw[3] = {sin(angle_u), sin(angle_v), sin(angle_w)};
    
    // Resize once
    dpc_points.resize(num_points);
    
    for (int k = 0; k < num_points; k++) {
        // Calculate position once
        double k_dis = k * dis_seg;
        double x = p_test.x() + k * dx_step;
        double y = p_test.y() + k * dy_step;
        double z = p_test.z() + k * dz_step;
        
        // Calculate distance cut
        double dis_cut = std::min(std::max(MIN_DIS_CUT, k_dis * sin_angle_rad), MAX_DIS_CUT);
        int dis_cut_int = static_cast<int>(dis_cut);
        
        // Set up point
        DynamicPointCloud::DPCPoint& point = dpc_points[k];
        point.cluster = cluster;
        point.blob = nullptr;
        point.x = x;
        point.y = y;
        point.z = z;
        point.wpid = wpid.ident();
        
        // Initialize arrays
        point.x_2d.resize(3); 
        point.y_2d.resize(3);
        point.wpid_2d.resize(3);
        point.wind = {-1e12, -1e12, -1e12};
        point.dist_cut = {dis_cut_int, dis_cut_int, dis_cut_int};
        
        // Calculate 2D projections
        for (size_t pindex = 0; pindex < 3; ++pindex) {
            point.x_2d[pindex].push_back(x);
            point.y_2d[pindex].push_back(cos_angle_uvw[pindex] * z - sin_angle_uvw[pindex] * y);
            point.wpid_2d[pindex].push_back(wpid.ident());
        }
    }

    return dpc_points;
}


void PointCloud::Facade::fill_wrap_points(const Cluster *cluster, const geo_point_t &point, const WirePlaneId& wpid, std::vector<std::vector<double>>& p_x, std::vector<std::vector<double>>& p_y, std::vector<std::vector<int>>& p_wpid){
    int apa = wpid.apa();
    int face = wpid.face();
    auto grouping = cluster->grouping();
    std::map<int, std::vector<double>> map_angles; // face -->angles 
    // std::cout << "fill_wrap_points: apa " << apa << " face " << face << std::endl;
    const auto wire_angles = grouping->wire_angles(apa, face);
    auto& angles = map_angles[face];
    angles.push_back(std::get<0>(wire_angles));
    angles.push_back(std::get<1>(wire_angles));
    angles.push_back(std::get<2>(wire_angles));

    // find the drift time ...
    const auto map_time_offset = grouping->get_time_offset().at(apa);
    const auto map_drift_speed = grouping->get_drift_speed().at(apa);
    double time_offset = map_time_offset.at(face);
    double drift_speed = map_drift_speed.at(face);

    // std::cout << "Test: " << map_time_offset.size() <<  " " << map_time_offset.begin()->first << " " <<  std::endl;

    auto anode = grouping->get_anode(apa);
    const auto iface = anode->faces()[face];
    const double time = drift2time(iface, time_offset, drift_speed, point.x());

    const auto map_pitch_mags = grouping->pitch_mags().at(apa);
    const auto map_proj_centers = grouping->proj_centers().at(apa);

    for (size_t pind = 0; pind < 3; ++pind) {
        // find the wire index ...
        const double angle = map_angles.at(face)[pind];
        const double pitch = map_pitch_mags.at(face).at(pind);
        const double center = map_proj_centers.at(face).at(pind);
        int wind = point2wind(point, angle, pitch, center);
        if (wind < 0) wind = 0;
        size_t max_wind = grouping->get_plane_channels(apa, face, iplane2layer[pind]).size() - 1;
        if (wind > max_wind) wind = max_wind;
        // get channel ...
        auto channel = grouping->get_plane_channel_wind(apa, face, iplane2layer[pind], wind);

        // get all wires
        auto wires = anode->wires(channel->ident());
        for (const auto &wire : wires) {
            auto wire_wpid = wire->planeid();
            // std::cout << "Test: " << map_time_offset.size() <<  " " << map_time_offset.begin()->first << " " << wire_wpid.face() << std::endl;
            p_x[pind].push_back(time2drift(anode->faces()[wire_wpid.face()], map_time_offset.at(wire_wpid.face()), map_drift_speed.at(wire_wpid.face()), time));
            if (map_angles.find(wire_wpid.face()) == map_angles.end()) {
                const auto wire_angles1 = grouping->wire_angles(apa, wire_wpid.face());
                auto& angles = map_angles[wire_wpid.face()];
                angles.push_back(std::get<0>(wire_angles1));
                angles.push_back(std::get<1>(wire_angles1));
                angles.push_back(std::get<2>(wire_angles1));
            }
            p_y[pind].push_back(wind2point2dproj(wind, map_angles.at(wire_wpid.face()).at(pind), map_pitch_mags.at(wire_wpid.face()).at(pind), map_proj_centers.at(wire_wpid.face()).at(pind)));
            p_wpid[pind].push_back(WirePlaneId(kAllLayers, wire_wpid.face(), wire_wpid.apa()).ident());
        }
    }
    
}
