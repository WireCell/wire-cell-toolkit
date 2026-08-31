#include "WireCellClus/DynamicPointCloud.h"
#include "WireCellClus/Facade.h"
#include "WireCellClus/Facade_Util.h"
#include "WireCellUtil/Logging.h"

#include <boost/histogram.hpp>
#include <boost/histogram/algorithm/sum.hpp>

#include <climits>

using namespace WireCell;
using namespace WireCell::Clus::Facade;
using spdlog::debug;

#ifdef __DEBUG__
#define LogDebug(x) std::cout << "[DPC]: " << __LINE__ << " : " << x << std::endl
#else
#define LogDebug(x)
#endif

const static std::array<int, 3> wind_bogus = {INT_MIN, INT_MIN, INT_MIN};
const static std::array<double, 3> dist_cut_bogus = {-1e12, -1e12, -1e12};


void DPCBatch::append(const DPCBatch& other)
{
    const size_t pbase = p2d_x.size();
    x.insert(x.end(), other.x.begin(), other.x.end());
    y.insert(y.end(), other.y.begin(), other.y.end());
    z.insert(z.end(), other.z.begin(), other.z.end());
    wpid.insert(wpid.end(), other.wpid.begin(), other.wpid.end());
    cluster.insert(cluster.end(), other.cluster.begin(), other.cluster.end());
    blob.insert(blob.end(), other.blob.begin(), other.blob.end());
    wind.insert(wind.end(), other.wind.begin(), other.wind.end());
    dist_cut.insert(dist_cut.end(), other.dist_cut.begin(), other.dist_cut.end());
    p2d_x.insert(p2d_x.end(), other.p2d_x.begin(), other.p2d_x.end());
    p2d_y.insert(p2d_y.end(), other.p2d_y.begin(), other.p2d_y.end());
    p2d_wpid.insert(p2d_wpid.end(), other.p2d_wpid.begin(), other.p2d_wpid.end());
    p2d_off.reserve(p2d_off.size() + other.p2d_off.size() - 1);
    for (size_t k = 1; k < other.p2d_off.size(); ++k) {
        p2d_off.push_back(pbase + other.p2d_off[k]);
    }
}

void DPCBatch::append(const DPCBatch& other, const std::vector<size_t>& rows)
{
    reserve(size() + rows.size());
    for (const size_t r : rows) {
        add_point(other.x[r], other.y[r], other.z[r], other.wpid[r],
                  other.cluster[r], other.blob[r], other.wind[r], other.dist_cut[r]);
        for (size_t p = 0; p < 3; ++p) {
            const auto [b, e] = other.proj_range(r, p);
            for (size_t k = b; k < e; ++k) {
                add_proj(other.p2d_x[k], other.p2d_y[k], other.p2d_wpid[k]);
            }
            end_plane();
        }
    }
}


DynamicPointCloud::nfkd_t &DynamicPointCloud::kd3d() const
{
    if (!m_kd3d) {
        m_kd3d = std::make_unique<nfkd_t>(3);
        // Zero-copy: the 3D tree reads the batch columns in place (their
        // member addresses are stable; vector growth is announced via
        // append_external in index_new_points).
        m_kd3d->bind_external({&m_pts.x, &m_pts.y, &m_pts.z});
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


void DynamicPointCloud::add_points(const DPCBatch &points) {
    if (points.empty()) {
        return;
    }
    size_t original_size = m_pts.size();
    m_pts.append(points);
    index_new_points(original_size);
}

void DynamicPointCloud::add_points(DPCBatch &&points) {
    if (points.empty()) {
        return;
    }
    size_t original_size = m_pts.size();
    if (original_size == 0) {
        // fresh cloud: take the whole storage in O(1)
        m_pts = std::move(points);
    }
    else {
        m_pts.append(points);
    }
    index_new_points(original_size);
}

void DynamicPointCloud::add_points(const DynamicPointCloud &other) {
    add_points(other.m_pts);
}

void DynamicPointCloud::add_points(const DynamicPointCloud &other, const std::vector<size_t> &rows) {
    if (rows.empty()) {
        return;
    }
    size_t original_size = m_pts.size();
    m_pts.append(other.m_pts, rows);
    index_new_points(original_size);
}

// Index m_pts[original_size..end) into the 3D and per-plane 2D k-d
// trees.  Reads from m_pts (NOT the caller's batch, which may have
// been moved from).
void DynamicPointCloud::index_new_points(size_t original_size) {
    const size_t nnew = m_pts.size() - original_size;

    // CSR structural invariant (replaces the old per-point 2D
    // projection size check).
    if (m_pts.p2d_off.size() != 3*m_pts.size() + 1) {
        raise<RuntimeError>("DynamicPointCloud: corrupt 2D projection CSR: %d offsets for %d points",
                            m_pts.p2d_off.size(), m_pts.size());
    }

    // Process 3D KD tree.  The tree is bound to the batch columns
    // (zero-copy); it only needs to be told how many points arrived.
    auto &kd3d = this->kd3d();

    // Prepare maps to store 2D points for each plane and track local-to-global mappings
    std::map<int, NFKDVec::Tree<double>::points_type> planes_pts;
    std::map<int, std::vector<size_t>> planes_global_indices;

    // Extract data and prepare for batch processing
    for (size_t i = 0; i < nnew; ++i) {
        size_t global_idx = original_size + i;

        // Skip 2D KD if wpid is not valid
        WirePlaneId wpid_volume(m_pts.wpid[global_idx]);
        if (wpid_volume.face() == -1 || wpid_volume.apa() == -1) {
            continue;
        }

        // Process 2D points for each plane
        for (size_t pindex = 0; pindex < 3; ++pindex) {
            const auto [jb, je] = m_pts.proj_range(global_idx, pindex);
            for (size_t j = jb; j < je; ++j) {
                WirePlaneId wpid_2d(m_pts.p2d_wpid[j]);
                WirePlaneId wpid_plane(iplane2layer[pindex], wpid_2d.face(), wpid_2d.apa());
                int key = wpid_plane.ident();
                // Initialize plane data structures if not exists
                if (planes_pts.find(key) == planes_pts.end()) {
                    planes_pts[key] = NFKDVec::Tree<double>::points_type(2);
                    planes_pts[key][0].reserve(nnew);
                    planes_pts[key][1].reserve(nnew);
                    planes_global_indices[key].reserve(nnew);
                }
                planes_pts[key][0].push_back(m_pts.p2d_x[j]);
                planes_pts[key][1].push_back(m_pts.p2d_y[j]);
                planes_global_indices[key].push_back(global_idx);
            }
        }

    }

    // Batch append 3D points
    kd3d.append_external(nnew);


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

geo_point_t DynamicPointCloud::get_center_point_radius(const geo_point_t &p_test, const double radius) const{
    auto &kd3d = this->kd3d();

    // Create query point
    std::vector<double> query = {p_test.x(), p_test.y(), p_test.z()};

    // Perform radius search (NFKDVec uses squared distance)
    auto results = kd3d.radius(radius * radius, query);

    // Calculate center point
    geo_point_t center(0, 0, 0);
    int ncount = 0;

    for (const auto &[idx, _] : results) {
        center.set(center.x() + m_pts.x[idx], center.y() + m_pts.y[idx], center.z() + m_pts.z[idx]);
        ncount++;
    }

    if (ncount > 0) {
        center.set(center.x() / ncount, center.y() / ncount, center.z() / ncount);
    }

    return center;
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
        return_results.emplace_back(sqrt(dist_squared), m_pts.cluster[global_idx], global_idx);
    }

    return return_results;
}



std::array<double, 3> DynamicPointCloud::get_angles(int face, int apa) const
{
    WirePlaneId wpid_volume(kAllLayers, face, apa);
    auto it = m_wpid_params.find(wpid_volume);
    if (it == m_wpid_params.end()) {
        return {0.0, 0.0, 0.0};
    }
    const auto& [_, angle_u, angle_v, angle_w] = it->second;
    return {angle_u, angle_v, angle_w};
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

    return std::make_tuple(distance, m_pts.cluster[global_idx], global_idx);
}

std::tuple<double, const Cluster *, size_t>
DynamicPointCloud::get_closest_2d_point_info_direct(
    double drift, double wire_perp, const int plane, const int face, const int apa) const
{
    auto &kd2d = this->kd2d(plane, face, apa);
    auto &l2g  = this->kd2d_l2g(plane, face, apa);

    // Query directly with (drift, wire_perp) — no angle projection needed because
    // convert_time_wire_2Dpoint already returns coordinates in the wire-perpendicular space
    // that matches the KD2D tree's storage format.
    const std::vector<double> query = {drift, wire_perp};
    auto results = kd2d.knn(1, query);

    if (results.empty()) {
        return std::make_tuple(-1.0, nullptr, static_cast<size_t>(-1));
    }

    const size_t local_idx  = results[0].first;
    const double distance   = sqrt(results[0].second);
    const size_t global_idx = l2g.at(local_idx);

    return std::make_tuple(distance, m_pts.cluster[global_idx], global_idx);
}


std::pair<double, double> DynamicPointCloud::hough_transform(const geo_point_t &origin, const double dis) const
{
    auto &kd3d = this->kd3d();

    // Create query point
    std::vector<double> query = {origin.x(), origin.y(), origin.z()};

    // Perform radius search
    auto results = kd3d.radius(dis * dis, query);

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
    if (results.empty()) {
        return {0.0, 0.0};
    }

    for (const auto &[idx, _] : results) {
        const auto *blob = m_pts.blob[idx];
        auto charge = blob ? blob->charge() : 1.0;

        if (charge <= 0) continue;

        const auto npoints = blob ? blob->npoints() : 1;
        const geo_point_t pt(m_pts.x[idx], m_pts.y[idx], m_pts.z[idx]);

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
    // Compare cell values explicitly: the accessor's own comparison operators
    // are ambiguous under Apple clang (boost::histogram detail::operators).
    auto it = std::max_element(indexed.begin(), indexed.end(), [](const auto& a, const auto& b) {
        return static_cast<double>(*a) < static_cast<double>(*b);
    });
    const auto &cell = *it;
    return {cell.bin(0).center(), cell.bin(1).center()};
}




geo_point_t DynamicPointCloud::vhough_transform(const geo_point_t &origin, const double dis) const
{
    const auto [th, phi] = hough_transform(origin, dis);
    return {sin(th) * cos(phi), sin(th) * sin(phi), cos(th)};
}



DPCBatch Clus::Facade::make_points_cluster(
    const Cluster *cluster, const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params, bool flag_wrap)
{
    DPCBatch batch;
    if (!cluster) {
        SPDLOG_WARN("make_points_cluster: null cluster return empty points");
        return batch;
    }

    const size_t num_points = cluster->npoints();
    batch.reserve(num_points);

    const auto &winds = cluster->wire_indices();

    // Cache commonly referenced WPIDs and their params to avoid map lookups
    std::unordered_map<int, std::tuple<geo_point_t, double, double, double>> cached_params;

    for (size_t ipt = 0; ipt < num_points; ++ipt) {
        geo_point_t pt = cluster->point3d(ipt);
        const auto wpid = cluster->wire_plane_id(ipt);
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

        batch.add_point(pt.x(), pt.y(), pt.z(), wpid_ident,
                        cluster, cluster->blob_with_point(ipt),
                        {winds[0][ipt], winds[1][ipt], winds[2][ipt]}, dist_cut_bogus);

        if (flag_wrap){
            fill_wrap_points(cluster, pt, WirePlaneId(wpid), batch);
        }else{
            for (size_t pindex = 0; pindex < 3; ++pindex) {
                batch.add_proj(pt.x(), cos(angle_uvw[pindex]) * pt.z() - sin(angle_uvw[pindex]) * pt.y(),
                               wpid_ident);
                batch.end_plane();
            }
        }
    }

    return batch;
}

DPCBatch Clus::Facade::make_points_cluster_steiner(const Cluster *cluster, const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params, bool flag_wrap){
    DPCBatch batch;
    if (!cluster) {
        SPDLOG_WARN("make_points_cluster_steiner: null cluster return empty points");
        return batch;
    }

    // Check if steiner point cloud exists
    if (!cluster->has_pc("steiner_pc")) {
        SPDLOG_WARN("make_points_cluster_steiner: cluster has no steiner_pc");
        return batch;
    }

    const auto& steiner_pc = cluster->get_pc("steiner_pc");
    const auto& coords = cluster->get_default_scope().coords;
    auto x_ptr = steiner_pc.get(coords.at(0));
    auto y_ptr = steiner_pc.get(coords.at(1));
    auto z_ptr = steiner_pc.get(coords.at(2));
    auto wpid_ptr = steiner_pc.get("wpid");
    if (!x_ptr || !y_ptr || !z_ptr || !wpid_ptr) {
        SPDLOG_WARN("make_points_cluster_steiner: steiner_pc missing coordinate arrays, returning empty");
        return batch;
    }
    const auto& x_coords = x_ptr->elements<double>();
    const auto& y_coords = y_ptr->elements<double>();
    const auto& z_coords = z_ptr->elements<double>();
    const auto& wpid_array = wpid_ptr->elements<WirePlaneId>();

    const size_t num_points = x_coords.size();
    batch.reserve(num_points);

    // Cache commonly referenced WPIDs and their params to avoid map lookups
    std::unordered_map<int, std::tuple<geo_point_t, double, double, double>> cached_params;

    for (size_t ipt = 0; ipt < num_points; ++ipt) {
        geo_point_t pt(x_coords[ipt], y_coords[ipt], z_coords[ipt]);
        const auto wpid = wpid_array[ipt];
        int wpid_ident = wpid.ident();

        // Check cache first, then populate if needed
        auto param_it = cached_params.find(wpid_ident);
        if (param_it == cached_params.end()) {
            auto wpid_it = wpid_params.find(wpid);
            if (wpid_it == wpid_params.end()) {
                raise<RuntimeError>("make_points_cluster_steiner: missing wpid params for wpid %s", wpid.name());
            }
            param_it = cached_params.emplace(wpid_ident, wpid_it->second).first;
        }

        const auto &[drift_dir, angle_u, angle_v, angle_w] = param_it->second;
        const double angle_uvw[3] = {angle_u, angle_v, angle_w};

        // Steiner points don't have blob associations nor wire indices
        batch.add_point(pt.x(), pt.y(), pt.z(), wpid_ident,
                        cluster, nullptr, wind_bogus, dist_cut_bogus);

        if (flag_wrap){
            fill_wrap_points(cluster, pt, wpid, batch);
        }else{
            for (size_t pindex = 0; pindex < 3; ++pindex) {
                batch.add_proj(pt.x(), cos(angle_uvw[pindex]) * pt.z() - sin(angle_uvw[pindex]) * pt.y(),
                               wpid_ident);
                batch.end_plane();
            }
        }
    }

    return batch;
}


DPCBatch Clus::Facade::make_points_direct(const Cluster *cluster, const IDetectorVolumes::pointer dv, const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params, std::vector<std::pair<geo_point_t,WirePlaneId>>& points_info, bool flag_wrap){
    DPCBatch batch;

    if (!cluster) {
        SPDLOG_WARN("make_points_cluster_skeleton: null cluster return empty points");
        return batch;
    }
    batch.reserve(points_info.size());

    // Cache for angle values per wpid to avoid repeated tuple unpacking
    std::unordered_map<int, std::array<double, 3>> wpid_angles_cache;

    for (auto& [test_point, wpid_test_point] : points_info) {
        // Skip points outside the detector volume (apa=-1) or with unknown wpid
        if (wpid_test_point.apa() == -1) continue;
        if (wpid_params.find(wpid_test_point) == wpid_params.end()) {
            raise<RuntimeError>("make_points_cluster: missing wpid params for wpid %s", wpid_test_point.name());
        }

        batch.add_point(test_point.x(), test_point.y(), test_point.z(),
                        wpid_test_point.ident(), cluster, nullptr,
                        wind_bogus, dist_cut_bogus);

        // (the apa() != -1 test always holds here: -1 was skipped above)
        {
            // Get cached angles if available
            std::array<double, 3> temp_angle_uvw;
            auto cache_it = wpid_angles_cache.find(wpid_test_point.ident());
            if (cache_it == wpid_angles_cache.end()) {
                const auto& [drift_dir, angle_u, angle_v, angle_w] = wpid_params.at(wpid_test_point);
                temp_angle_uvw = {angle_u, angle_v, angle_w};
                wpid_angles_cache[wpid_test_point.ident()] = temp_angle_uvw;
            } else {
                temp_angle_uvw = cache_it->second;
            }


            if (flag_wrap){
                fill_wrap_points(cluster, test_point, wpid_test_point, batch);
            }else{
                for (size_t pindex = 0; pindex < 3; ++pindex) {
                    batch.add_proj(test_point.x(),
                                   cos(temp_angle_uvw[pindex]) * test_point.z() -
                                       sin(temp_angle_uvw[pindex]) * test_point.y(),
                                   wpid_test_point.ident());
                    batch.end_plane();
                }
            }
        }
    }


    return batch;

}


DPCBatch
Clus::Facade::make_points_cluster_skeleton(
    const Cluster *cluster, const IDetectorVolumes::pointer dv,
    const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params,
    const std::vector<size_t>& path_wcps,
    bool flag_wrap,
    const double step)
{
    DPCBatch batch;

    if (!cluster) {
        SPDLOG_WARN("make_points_cluster_skeleton: null cluster return empty points");
        return batch;
    }

    // Estimate capacity to avoid reallocations
    size_t estimated_capacity = path_wcps.size() * 2; // Rough estimate
    batch.reserve(estimated_capacity);

    // Cache for angle values per wpid to avoid repeated tuple unpacking
    std::unordered_map<int, std::array<double, 3>> wpid_angles_cache;

    geo_point_t prev_wcp = cluster->point3d(path_wcps.front());
    auto prev_wpid = cluster->wire_plane_id(path_wcps.front());

    // Pre-computed constants
    const double dist_cut_value = 2.4 * units::cm;
    const std::array<double, 3> dist_cut_skel = {dist_cut_value, dist_cut_value, dist_cut_value};

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
            batch.add_point(test_point.x(), test_point.y(), test_point.z(),
                            wpid_test_point.ident(), cluster, nullptr,
                            wind_bogus, dist_cut_skel);

            if (flag_wrap){
                fill_wrap_points(cluster, test_point, WirePlaneId(wpid_test_point), batch);
            }else{
                for (size_t pindex = 0; pindex < 3; ++pindex) {
                    batch.add_proj(test_point.x(),
                                   cos(angle_uvw[pindex]) * test_point.z() - sin(angle_uvw[pindex]) * test_point.y(),
                                   wpid_test_point.ident());
                    batch.end_plane();
                }
            }
        }
        else {
            int num_points = int(dis / step) + 1;

            // Pre-compute direction vectors to avoid recalculation in loop
            double dx = (test_point.x() - prev_wcp.x()) / num_points;
            double dy = (test_point.y() - prev_wcp.y()) / num_points;
            double dz = (test_point.z() - prev_wcp.z()) / num_points;

            for (int k = 0; k != num_points; k++) {
                // Faster interpolation with pre-computed increments
                double t = (k + 1.0);
                const double px = prev_wcp.x() + t * dx;
                const double py = prev_wcp.y() + t * dy;
                const double pz = prev_wcp.z() + t * dz;

                geo_point_t temp_point(px, py, pz);
                auto temp_wpid = WirePlaneId(get_wireplaneid(temp_point, prev_wpid, wpid_test_point, dv));

                batch.add_point(px, py, pz, temp_wpid.ident(), cluster, nullptr,
                                wind_bogus, dist_cut_skel);

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
                        fill_wrap_points(cluster, temp_point, temp_wpid, batch);
                    }else{
                        for (size_t pindex = 0; pindex < 3; ++pindex) {
                            batch.add_proj(px,
                                           cos(temp_angle_uvw[pindex]) * pz -
                                               sin(temp_angle_uvw[pindex]) * py,
                                           temp_wpid.ident());
                            batch.end_plane();
                        }
                    }
                }
                else {
                    // out-of-volume interpolated point: no 2D projections
                    batch.end_plane();
                    batch.end_plane();
                    batch.end_plane();
                }
            }
        }

        prev_wcp = test_point;
        prev_wpid = wpid_test_point;
    }

    return batch;
}



DPCBatch Clus::Facade::make_points_linear_extrapolation(
    const Cluster *cluster, const geo_point_t &p_test, const geo_point_t &dir_unmorm, const double range,
    const double step, const double angle, const IDetectorVolumes::pointer dv,
    const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params,
    const WirePlaneId &seed_wpid)
{
    DPCBatch batch;

    if (!cluster) {
        SPDLOG_WARN("make_points_linear_extrapolation: null cluster return empty points");
        return batch;
    }

    // Fallback wpid: the extrapolation seed's volume if given, else the
    // grouping's first wpid (legacy single-volume behavior).
    const auto wpid = (seed_wpid.apa() >= 0) ? seed_wpid : *(cluster->grouping()->wpids().begin());

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

    // Multi-volume groupings (drift-side groups): bucket each synthetic point
    // into the volume containing it so it lands in the correct per-(face,apa)
    // 2D KD trees; out-of-volume points keep the seed wpid.  Per-wpid
    // projection constants are cached on first use.
    const bool multi_volume = (wpid_params.size() > 1) && dv;
    std::map<int, std::array<double, 6>> trig_cache;  // wpid ident -> cos_uvw, sin_uvw

    batch.reserve(num_points);

    for (int k = 0; k < num_points; k++) {
        // Calculate position once
        double k_dis = k * dis_seg;
        double x = p_test.x() + k * dx_step;
        double y = p_test.y() + k * dy_step;
        double z = p_test.z() + k * dz_step;

        // Calculate distance cut
        const double dis_cut = std::floor(std::min(std::max(MIN_DIS_CUT, k_dis * sin_angle_rad), MAX_DIS_CUT));

        WirePlaneId pt_wpid = wpid;
        const double* cos_uvw = cos_angle_uvw;
        const double* sin_uvw = sin_angle_uvw;
        if (multi_volume) {
            const auto vol = dv->contained_by(Point(x, y, z));
            if (vol.face() >= 0) {
                const WirePlaneId cand(kAllLayers, vol.face(), vol.apa());
                if (wpid_params.find(cand) != wpid_params.end()) pt_wpid = cand;
            }
            auto [tit, fresh] = trig_cache.try_emplace(pt_wpid.ident());
            if (fresh) {
                const auto& [dd, au, av, aw] = wpid_params.at(pt_wpid);
                tit->second = {cos(au), cos(av), cos(aw), sin(au), sin(av), sin(aw)};
            }
            cos_uvw = tit->second.data();
            sin_uvw = tit->second.data() + 3;
        }

        batch.add_point(x, y, z, pt_wpid.ident(), cluster, nullptr,
                        wind_bogus, {dis_cut, dis_cut, dis_cut});

        // Calculate 2D projections
        for (size_t pindex = 0; pindex < 3; ++pindex) {
            batch.add_proj(x, cos_uvw[pindex] * z - sin_uvw[pindex] * y, pt_wpid.ident());
            batch.end_plane();
        }
    }

    return batch;
}


void Clus::Facade::fill_wrap_points(const Cluster *cluster, const geo_point_t &point, const WirePlaneId& wpid, DPCBatch& batch){
    int apa = wpid.apa();
    int face = wpid.face();
    auto grouping = cluster->grouping();
    std::map<int, std::vector<double>> map_angles; // face -->angles
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
        auto plane_ptr =iface->plane(pind);
        const auto& wires_all = plane_ptr->wires();
        size_t max_wind = wires_all.size() - 1;
        if ((size_t)wind > max_wind) wind = max_wind;
        // get channel ...
        auto wire = wires_all[wind];
        int channel_number = wire->channel();

        // get all wires
        auto wires = anode->wires(channel_number);
        for (const auto &wire : wires) {
            auto wire_wpid = wire->planeid();

            const double wx = time2drift(anode->faces()[wire_wpid.face()], map_time_offset.at(wire_wpid.face()), map_drift_speed.at(wire_wpid.face()), time);
            if (map_angles.find(wire_wpid.face()) == map_angles.end()) {
                const auto wire_angles1 = grouping->wire_angles(apa, wire_wpid.face());
                auto& angles = map_angles[wire_wpid.face()];
                angles.push_back(std::get<0>(wire_angles1));
                angles.push_back(std::get<1>(wire_angles1));
                angles.push_back(std::get<2>(wire_angles1));
            }

            double wy;
            // Check if this wire is the same as the original wire (wire index, apa, face are all the same)
            if (wire_wpid.apa() == wpid.apa() && wire_wpid.face() == wpid.face() && wire->index() == wind) {
                // Use the original wire's angles to calculate the projection
                wy = cos(angles[pind]) * point.z() - sin(angles[pind]) * point.y();
            } else {
                // Use the current algorithm
                wy = wind2point2dproj(wind, map_angles.at(wire_wpid.face()).at(pind), map_pitch_mags.at(wire_wpid.face()).at(pind), map_proj_centers.at(wire_wpid.face()).at(pind));
            }
            batch.add_proj(wx, wy, WirePlaneId(kAllLayers, wire_wpid.face(), wire_wpid.apa()).ident());
        }
        batch.end_plane();
    }

}
