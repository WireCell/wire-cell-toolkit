#include "WireCellClus/DynamicPointCloud.h"
#include "WireCellClus/Facade_Util.h"
#include "WireCellClus/Facade.h"
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

DynamicPointCloud::nfkd_t& DynamicPointCloud::kd3d() const {
    if (!m_kd3d) {
        m_kd3d = std::make_unique<nfkd_t>(3);
    }
    return *m_kd3d;
}

DynamicPointCloud::nfkd_t& DynamicPointCloud::kd2d(const int plane, const int face, const int apa) const {
    WirePlaneId wpid(iplane2layer[plane], face, apa);
    auto iter = m_kd2d.find(wpid);
    if (iter == m_kd2d.end()) {
        m_kd2d[wpid] = std::make_unique<nfkd_t>(2);
    }
    return *m_kd2d[wpid];
}

const std::unordered_map<size_t, size_t>& DynamicPointCloud::kd2d_l2g(const int plane, const int face, const int apa) const {
    WirePlaneId wpid(iplane2layer[plane], face, apa);
    auto iter = m_kd2d_index_l2g.find(wpid);
    if (iter == m_kd2d_index_l2g.end()) {
        raise<RuntimeError>("DynamicPointCloud: missing 2D index l2g for wpid %s", wpid.name());
    }
    return iter->second;
}

const std::unordered_map<size_t, size_t>& DynamicPointCloud::kd2d_g2l(const int plane, const int face, const int apa) const {
    WirePlaneId wpid(iplane2layer[plane], face, apa);
    auto iter = m_kd2d_index_g2l.find(wpid);
    if (iter == m_kd2d_index_g2l.end()) {
        raise<RuntimeError>("DynamicPointCloud: missing 2D index g2l for wpid %s", wpid.name());
    }
    return iter->second;
}


void DynamicPointCloud::add_points(const std::vector<DPCPoint> &points) {
    // move data to self
    m_points = std::move(points);

    // process KD trees
    auto &kd3d = this->kd3d();
    for (size_t ipt = 0; ipt < m_points.size(); ++ipt) {
        const auto &pt = m_points.at(ipt);
        kd3d.append({{pt.x}, {pt.y}, {pt.z}});
        WirePlaneId wpid_volume(pt.wpid);
        if (!wpid_volume or wpid_volume.layer() != kAllLayers) {
            SPDLOG_TRACE("DynamicPointCloud: !wpid {} skipping 2D KD", wpid_volume.name());
            continue;
        }
        if (pt.x_2d.size() != 3 or pt.y_2d.size() != 3) {
            raise<RuntimeError>("DynamicPointCloud: unexpected 2D projection size x_2d %d y_2d %d", pt.x_2d.size(), pt.y_2d.size());
        }
        for (size_t pindex = 0; pindex < 3; ++pindex) {
            auto &kd2d = this->kd2d(pindex, wpid_volume.face(), wpid_volume.apa());
            kd2d.append({{pt.x_2d[pindex]}, {pt.y_2d[pindex]}});
            WirePlaneId wpid_plane(iplane2layer[pindex], wpid_volume.face(), wpid_volume.apa());
            m_kd2d_index_l2g[wpid_plane][kd2d.npoints() - 1] = ipt;
            m_kd2d_index_g2l[wpid_plane][ipt] = kd2d.npoints() - 1;
        }
    }
}

std::vector<std::tuple<double, const Cluster *, size_t>> 
DynamicPointCloud::get_2d_points_info(const geo_point_t &p, const double radius, const int plane, const int face, const int apa) const {
    WirePlaneId wpid(iplane2layer[plane], face, apa);
    auto iter = m_kd2d.find(wpid);
    if (iter == m_kd2d.end()) {
        raise<RuntimeError>("DynamicPointCloud: missing 2D KD for wpid %s", wpid.name());
    }
    auto &kd2d = *iter->second;
    auto &l2g = this->kd2d_l2g(plane, face, apa);
    auto &g2l = this->kd2d_g2l(plane, face, apa);

    if (m_wpid_params.find(wpid) == m_wpid_params.end()) {
        raise<RuntimeError>("DynamicPointCloud: missing wpid params for wpid %s", wpid.name());
    }
    const auto [_, angle_u, angle_v, angle_w] = m_wpid_params.at(wpid);
    double angle_uvw[3] = {angle_u, angle_v, angle_w};

    // Prepare query point
    std::vector<double> query = {p.x(), cos(angle_uvw[plane]) * p.z() - sin(angle_uvw[plane]) * p.y()};
    
    // Perform radius search
    auto results = kd2d.radius(radius * radius, query);
    
    std::vector<std::tuple<double, const Cluster *, size_t>> return_results;
    for (const auto& [local_idx, dist_squared] : results) {
        size_t global_idx = l2g.at(local_idx);
        const auto& pt = m_points[global_idx];
        return_results.push_back(std::make_tuple(sqrt(dist_squared), pt.cluster, global_idx));
    }
    
    return return_results;
}

std::tuple<double, const Cluster *, size_t> 
DynamicPointCloud::get_closest_2d_point_info(const geo_point_t &p, const int plane, const int face, const int apa) const {
    WirePlaneId wpid(iplane2layer[plane], face, apa);
    auto iter = m_kd2d.find(wpid);
    if (iter == m_kd2d.end()) {
        raise<RuntimeError>("DynamicPointCloud: missing 2D KD for wpid %s", wpid.name());
    }
    auto &kd2d = *iter->second;
    auto &l2g = this->kd2d_l2g(plane, face, apa);
    auto &g2l = this->kd2d_g2l(plane, face, apa);

    if (m_wpid_params.find(wpid) == m_wpid_params.end()) {
        raise<RuntimeError>("DynamicPointCloud: missing wpid params for wpid %s", wpid.name());
    }
    const auto [_, angle_u, angle_v, angle_w] = m_wpid_params.at(wpid);
    double angle_uvw[3] = {angle_u, angle_v, angle_w};

    // Prepare query point
    std::vector<double> query = {p.x(), cos(angle_uvw[plane]) * p.z() - sin(angle_uvw[plane]) * p.y()};
    
    // Perform radius search
    auto results = kd2d.knn(1, query);
    
    if (results.size() == 1) {
        size_t global_idx = l2g.at(results[0].first);
        const auto& pt = m_points[global_idx];
        return std::make_tuple(sqrt(results[0].second), pt.cluster, global_idx);
    }
    else {
        return std::make_tuple(-1, nullptr, -1);
    }
}

std::pair<double, double> DynamicPointCloud::hough_transform(const geo_point_t &origin, const double dis) const
{
    
    auto &kd3d = this->kd3d();
    
    // Find points within radius
    std::vector<size_t> point_indices;
    std::vector<const Blob*> blobs;
    std::vector<geo_point_t> pts;
    
    // Create query point
    std::vector<double> query = {origin.x(), origin.y(), origin.z()};
    
    // Perform radius search
    auto results = kd3d.radius(dis * dis, query);
    
    for (const auto& [idx, _] : results) {
        const auto& pt = m_points[idx];
        pts.push_back({pt.x, pt.y, pt.z});
        blobs.push_back(pt.blob);
    }
    
    namespace bh = boost::histogram;
    namespace bha = boost::histogram::algorithm;
    
    constexpr double pi = 3.141592653589793;
    
    // Parameter axis 1 is theta angle
    const int nbins1 = 180;
    auto theta_param = [](const Vector& dir) {
        const Vector Z(0, 0, 1);
        return acos(Z.dot(dir));
    };
    double min1 = 0, max1 = pi;
    
    // Parameter axis 2 is phi angle
    const int nbins2 = 360;
    const double min2 = -pi;
    const double max2 = +pi;
    auto phi_param = [](const Vector& dir) {
        const Vector X(1, 0, 0);
        const Vector Y(0, 1, 0);
        return atan2(Y.dot(dir), X.dot(dir));
    };
    
    auto hist = bh::make_histogram(
        bh::axis::regular<>(nbins1, min1, max1), 
        bh::axis::regular<>(nbins2, min2, max2)
    );
    
    for (size_t ind = 0; ind < blobs.size(); ++ind) {
        const auto* blob = blobs[ind];
        auto charge = blob ? blob->charge() : 1.0;
        
        if (charge <= 0) continue;
        
        const auto npoints = blob ? blob->npoints() : 1;
        const auto& pt = pts[ind];
        
        const Vector dir = (pt - origin).norm();
        const double r = (pt - origin).magnitude();
        
        const double p1 = theta_param(dir);
        const double p2 = phi_param(dir);
        
        if (r < 10 * units::cm) {
            hist(p1, p2, bh::weight(charge / npoints));
        }
        else {
            hist(p1, p2, bh::weight(charge / npoints * pow(10 * units::cm / r, 2)));
        }
    }
    
    auto indexed = bh::indexed(hist);
    auto it = std::max_element(indexed.begin(), indexed.end());
    const auto& cell = *it;
    return {cell.bin(0).center(), cell.bin(1).center()};
}

geo_point_t DynamicPointCloud::vhough_transform(const geo_point_t &origin, const double dis) const
{
    const auto [th, phi] = hough_transform(origin, dis);
    return {sin(th) * cos(phi), sin(th) * sin(phi), cos(th)};
}

std::vector<DynamicPointCloud::DPCPoint>
    make_points_cluster(const Cluster *cluster,
                        const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params) {
    std::vector<DynamicPointCloud::DPCPoint> dpc_points;
    
    if (!cluster) {
        SPDLOG_WARN("make_points_cluster: null cluster return empty points");
        return dpc_points;
    }
    // use this so we can add a scope in the future
    const auto& points_3d = cluster->points();
    const auto& winds = cluster->wire_indices();
    const auto& wpids = cluster->points_property<int>("wpid");

    dpc_points.resize(cluster->npoints());
    for (size_t ipt=0; ipt<cluster->npoints(); ++ipt) {
        const auto x = points_3d[0][ipt];
        const auto y = points_3d[1][ipt];
        const auto z = points_3d[2][ipt];
        const auto& wpid = WirePlaneId(wpids[ipt]);
        if (wpid_params.find(wpid) == wpid_params.end()) {
            raise<RuntimeError>("make_points_cluster: missing wpid params for wpid %s", wpid.name());
        }
        const auto [drift_dir, angle_u, angle_v, angle_w] = wpid_params.at(wpid);
        const double angle_uvw[3] = {angle_u, angle_v, angle_w};
        
        DynamicPointCloud::DPCPoint point;
        point.x = x;
        point.y = y;
        point.z = z;
        point.wpid = wpid;
        point.cluster = cluster;
        point.blob = cluster->blob_with_point(ipt);
        point.x_2d.resize(3);
        point.y_2d.resize(3);
        point.wind = {winds[0][ipt], winds[1][ipt], winds[2][ipt]};
        point.dist_cut = {-1e12, -1e12, -1e12};
        
        for (size_t pindex = 0; pindex < 3; ++pindex) {
            point.x_2d[pindex] = x;
            point.y_2d[pindex] = cos(angle_uvw[pindex]) * z - sin(angle_uvw[pindex]) * y;
        }
        
        dpc_points[ipt] = std::move(point);
    }

    return dpc_points;
}