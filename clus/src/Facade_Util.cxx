#include "WireCellClus/Facade.h"
#include "WireCellClus/Facade_Cluster.h"
#include <boost/container_hash/hash.hpp>

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Facade;
// using WireCell::PointCloud::Dataset;
using namespace WireCell::PointCloud::Tree;  // for "Points" node value type
// using WireCell::PointCloud::Tree::named_pointclouds_t;

#include "WireCellUtil/Logging.h"
using spdlog::debug;

// int global_counter_get_closest_wcpoint = 0;

// #define __DEBUG__
#ifdef __DEBUG__
#define LogDebug(x) std::cout << "[yuhw]: " << __LINE__ << " : " << x << std::endl
#else
#define LogDebug(x)
#endif

/// unused
#if 0
namespace {
    // helper to dump a dataset
    std::string dump_ds(const WireCell::PointCloud::Dataset& ds) {
        std::stringstream ss;
        for (const auto& key : ds.keys()) {;
            const auto& arr = ds.get(key);
            ss << " {" << key << ":" << arr->dtype() << ":" << arr->shape()[0] << "} ";
            // const auto& arr = ds.get(key)->elements<float>();
            // for(auto elem : arr) {
            //     ss << elem << " ";
            // }
        }
        return ss.str();
    }
    std::string dump_pcs(const ScopedBase::pointclouds_t& pcs) {
        std::stringstream ss;
        for (const auto& pc : pcs) {
            ss << dump_ds(pc) << std::endl;
        }
        return ss.str();
    }
}
#endif

void Facade::Simple3DPointCloud::add(const point_type& new_pt) {
    if (new_pt.size() != 3) {
        raise<ValueError>("points must be 3D");
    }
    for (size_t ind=0; ind<3; ++ind) {
        points()[ind].push_back(new_pt[ind]);
    }
}
const Facade::Simple3DPointCloud::nfkd_t& Facade::Simple3DPointCloud::kd(bool rebuild) const
{
    if (rebuild) m_kd = nullptr;
    if (m_kd) return *m_kd;
    m_kd = std::make_unique<nfkd_t>(points());
    return *m_kd;
}
Facade::Simple3DPointCloud::results_type Facade::Simple3DPointCloud::get_closest_index(const geo_point_t& p, const size_t N) const {
    return kd().knn(N, p);
}
std::pair<size_t, geo_point_t> Facade::Simple3DPointCloud::get_closest_wcpoint(const geo_point_t& p) const {
    // global_counter_get_closest_wcpoint++;
    const auto knn_res = kd().knn(1, p);
    if (knn_res.size() != 1) {
        raise<ValueError>("no points found");
    }
    const auto ind = knn_res[0].first;
    geo_point_t pt = {points()[0][ind], points()[1][ind], points()[2][ind]};
    // std::cout << "get_closest_wcpoint: " << p << " " << ind << " " << pt << " " << knn_res[0].second << std::endl;
    return std::make_pair(ind, pt);
}
std::pair<int, double> Facade::Simple3DPointCloud::get_closest_point_along_vec(const geo_point_t& p_test1,
                                                                          const geo_point_t& dir, double test_dis,
                                                                          double dis_step, double angle_cut,
                                                                          double dis_cut) const
{
    geo_point_t p_test;

    double min_dis = 1e9;
    double min_dis1 = 1e9;
    int min_index = -1;

    for (int i = 0; i != int(test_dis / dis_step) + 1; i++) {
        p_test.set(p_test1.x() + dir.x() * i * dis_step, p_test1.y() + dir.y() * i * dis_step,
                   p_test1.z() + dir.z() * i * dis_step);

        auto knn_res = kd().knn(1, p_test);
        if (knn_res.size() == 0) {
            raise<ValueError>("no points found");
        }
        auto ind = knn_res[0].first;
        geo_point_t pts = {points()[0][ind], points()[1][ind], points()[2][ind]};

        double dis = sqrt(pow(p_test.x() - pts.x(), 2) + pow(p_test.y() - pts.y(), 2) +
                          pow(p_test.z() - pts.z(), 2));
        double dis1 = sqrt(pow(p_test1.x() - pts.x(), 2) + pow(p_test1.y() - pts.y(), 2) +
                           pow(p_test1.z() - pts.z(), 2));

        if (dis < std::min(dis1 * tan(angle_cut / 180. * 3.1415926), dis_cut)) {
            if (dis < min_dis) {
                min_dis = dis;
                min_index = ind;
                min_dis1 = dis1;
            }
            /// HARDCODED:
            if (dis < 3 * units::cm) return std::make_pair(ind, dis1);
        }
    }

    return std::make_pair(min_index, min_dis1);
}
std::tuple<int, int, double> Facade::Simple3DPointCloud::get_closest_points(const Simple3DPointCloud& other) const {

    int p1_index = 0;
    int p2_index = 0;
    geo_point_t p1 = point(p1_index);
    geo_point_t p2 = other.point(p2_index);
    int p1_save = 0;
    int p2_save = 0;
    double min_dis = 1e9;

    int prev_index1 = -1;
    int prev_index2 = -1;
    while (p1_index != prev_index1 || p2_index != prev_index2) {
        prev_index1 = p1_index;
        prev_index2 = p2_index;
        std::tie(p2_index, p2) = other.get_closest_wcpoint(p1);
        std::tie(p1_index, p1) = get_closest_wcpoint(p2);
    }
    // std::cout << "get_closest_points: " << p1_index << " " << p2_index << std::endl;
    double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
    if (dis < min_dis) {
        min_dis = dis;
        p1_save = p1_index;
        p2_save = p2_index;
    }

    prev_index1 = -1;
    prev_index2 = -1;
    p1_index = points()[0].size() - 1;
    p2_index = 0;
    p1 = point(p1_index);
    p2 = other.point(p2_index);
    while (p1_index != prev_index1 || p2_index != prev_index2) {
        prev_index1 = p1_index;
        prev_index2 = p2_index;
        std::tie(p2_index, p2) = other.get_closest_wcpoint(p1);
        std::tie(p1_index, p1) = get_closest_wcpoint(p2);
    }
    // std::cout << "get_closest_points: " << p1_index << " " << p2_index << std::endl;
    dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
    if (dis < min_dis) {
        min_dis = dis;
        p1_save = p1_index;
        p2_save = p2_index;
    }

    prev_index1 = -1;
    prev_index2 = -1;
    p1_index = 0;
    p2_index = other.points()[0].size() - 1;
    p1 = point(p1_index);
    p2 = other.point(p2_index);
    while (p1_index != prev_index1 || p2_index != prev_index2) {
        prev_index1 = p1_index;
        prev_index2 = p2_index;
        std::tie(p2_index, p2) = other.get_closest_wcpoint(p1);
        std::tie(p1_index, p1) = get_closest_wcpoint(p2);
    }
    // std::cout << "get_closest_points: " << p1_index << " " << p2_index << std::endl;
    dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
    if (dis < min_dis) {
        min_dis = dis;
        p1_save = p1_index;
        p2_save = p2_index;
    }

    prev_index1 = -1;
    prev_index2 = -1;
    p1_index = points()[0].size() - 1;
    p2_index = other.points()[0].size() - 1;
    p1 = point(p1_index);
    p2 = other.point(p2_index);
    while (p1_index != prev_index1 || p2_index != prev_index2) {
        prev_index1 = p1_index;
        prev_index2 = p2_index;
        std::tie(p2_index, p2) = other.get_closest_wcpoint(p1);
        std::tie(p1_index, p1) = get_closest_wcpoint(p2);
    }
    // std::cout << "get_closest_points: " << p1_index << " " << p2_index << std::endl;
    dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
    if (dis < min_dis) {
        min_dis = dis;
        p1_save = p1_index;
        p2_save = p2_index;
    }

    return std::make_tuple(p1_save, p2_save, min_dis);
}

std::ostream& Facade::operator<<(std::ostream& os, const Simple3DPointCloud& s3dpc)
{
    const size_t npts = s3dpc.points()[0].size();
    os << "Simple3DPointCloud " << " npts=" << npts;
    for (size_t ind=0; ind<npts; ++ind) {
        os << " " << s3dpc.point(ind);
    }
    return os;
}

Facade::Multi2DPointCloud::Multi2DPointCloud(const double angle_u, const double angle_v, const double angle_w) : angle_uvw{angle_u, angle_v, angle_w} {
    for (size_t plane = 0; plane < 3; ++plane) {
        points(plane).resize(2);
    }
}

void Facade::Multi2DPointCloud::add(const geo_point_t& new_pt) {
    for (size_t plane = 0; plane < 3; ++plane) {
        double x = new_pt[0];
        double y = cos(angle_uvw[plane]) * new_pt[2] - sin(angle_uvw[plane]) * new_pt[1];
        points(plane)[0].push_back(x);
        points(plane)[1].push_back(y);
    }
}

const Facade::Multi2DPointCloud::nfkd_t& Facade::Multi2DPointCloud::kd(const size_t plane, const bool rebuild) const
{
    if (rebuild) m_kd[plane] = nullptr;
    if (m_kd[plane]) return *m_kd[plane];
    m_kd[plane] = std::make_unique<nfkd_t>(points(plane));
    return *m_kd[plane];
}
std::pair<int, double> Facade::Multi2DPointCloud::get_closest_2d_dis(const geo_point_t& p, size_t plane) const
{
    double x = p[0];
    double y = cos(angle_uvw[plane]) * p.z() - sin(angle_uvw[plane]) * p.y();
    std::vector<double> query = {x, y};
    const auto& res = kd(plane).knn(1, query);

    if (res.size() == 1)
        return std::make_pair(res[0].first, sqrt(res[0].second)); /// Note, res.second is the square of the distance!
    else
        return std::make_pair(-1, 1e9);
}
std::vector<std::pair<size_t, double>> Facade::Multi2DPointCloud::get_closest_2d_index(const geo_point_t &p, const double radius, size_t plane) const
{
    double x = p[0];
    double y = cos(angle_uvw[plane]) * p.z() - sin(angle_uvw[plane]) * p.y();
    std::vector<double> query = {x, y};
    // const auto& res = kd(plane).knn(1, query);
    const auto& res = kd(plane).radius(radius * radius, query);
    std::vector<std::pair<size_t, double>> ret;
    for (const auto& r : res) {
        ret.push_back(std::make_pair(r.first, sqrt(r.second)));
    }
    return ret;
}

std::ostream& Facade::operator<<(std::ostream& os, const Multi2DPointCloud& m2dpc)
{
    os << "Multi2DPointCloud " << " get_num_points " << m2dpc.get_num_points();
    return os;
}

Facade::DynamicPointCloud::DynamicPointCloud(const double angle_u, const double angle_v, const double angle_w)
  : m_pc2d(angle_u, angle_v, angle_w)
{
}

void Facade::DynamicPointCloud::add_points(const Cluster* cluster, const int flag, const double step)
{
    size_t current_size = get_num_points();
    const auto& winds = cluster->wire_indices();

    if (flag == 0) {
        // add actual points in
        // WCP::WCPointCloud<double>& pcloud = cluster->get_point_cloud()->get_cloud();
        // WCP::WC2DPointCloud<double>& pcloud_u = cluster->get_point_cloud()->get_cloud_u();
        // WCP::WC2DPointCloud<double>& pcloud_v = cluster->get_point_cloud()->get_cloud_v();
        // WCP::WC2DPointCloud<double>& pcloud_w = cluster->get_point_cloud()->get_cloud_w();

        // cloud.pts.resize(current_size + pcloud.pts.size());
        // cloud_u.pts.resize(current_size + pcloud.pts.size());
        // cloud_v.pts.resize(current_size + pcloud.pts.size());
        // cloud_w.pts.resize(current_size + pcloud.pts.size());
        // vec_index_cluster.resize(current_size + pcloud.pts.size());

        for (size_t i = 0; i != cluster->npoints(); i++) {
            // vec_index_cluster.at(current_size + i) = cluster;
            m_clusters.push_back(cluster);
            m_pc3d.add({cluster->point3d(i).x(), cluster->point3d(i).y(), cluster->point3d(i).z()});
            m_pc2d.add(cluster->point3d(i));
            for (size_t plane = 0; plane < 3; ++plane) {
                m_winds[plane].push_back(winds[plane][i]);
            }
            m_blobs.push_back(cluster->blob_with_point(i));

            // cloud.pts[current_size + i].x = pcloud.pts.at(i).x;
            // cloud.pts[current_size + i].y = pcloud.pts.at(i).y;
            // cloud.pts[current_size + i].z = pcloud.pts.at(i).z;
            // cloud.pts[current_size + i].index_u = pcloud.pts.at(i).index_u;
            // cloud.pts[current_size + i].index_v = pcloud.pts.at(i).index_v;
            // cloud.pts[current_size + i].index_w = pcloud.pts.at(i).index_w;
            // cloud.pts[current_size + i].mcell = pcloud.pts.at(i).mcell;
            // cloud.pts[current_size + i].index = current_size + i;
            // cloud_u.pts[current_size + i].x = pcloud_u.pts.at(i).x;
            // cloud_u.pts[current_size + i].y = pcloud_u.pts.at(i).y;
            // cloud_u.pts[current_size + i].index = current_size + i;
            // cloud_v.pts[current_size + i].x = pcloud_v.pts.at(i).x;
            // cloud_v.pts[current_size + i].y = pcloud_v.pts.at(i).y;
            // cloud_v.pts[current_size + i].index = current_size + i;
            // cloud_w.pts[current_size + i].x = pcloud_w.pts.at(i).x;
            // cloud_w.pts[current_size + i].y = pcloud_w.pts.at(i).y;
            // cloud_w.pts[current_size + i].index = current_size + i;
        }
        // if (pcloud.pts.size() > 0) {
        //     index->addPoints(current_size, current_size + pcloud.pts.size() - 1);
        //     index_u->addPoints(current_size, current_size + pcloud.pts.size() - 1);
        //     index_v->addPoints(current_size, current_size + pcloud.pts.size() - 1);
        //     index_w->addPoints(current_size, current_size + pcloud.pts.size() - 1);
        // }
    }
    else {
        // add skeleton points in
        const std::list<size_t>& path_wcps = cluster->get_path_wcps();

        std::vector<geo_point_t> pts;
        geo_point_t prev_wcp = cluster->point3d(path_wcps.front());
        for (auto it = path_wcps.begin(); it != path_wcps.end(); it++) {
            geo_point_t test_point = cluster->point3d(*it);
            double dis =
                sqrt(pow(test_point.x() - prev_wcp.x(), 2) + pow(test_point.y() - prev_wcp.y(), 2) + pow(test_point.z() - prev_wcp.z(), 2));
            if (dis <= step) {
                // geo_point_t current_pt((*it).x(), (*it).y(), (*it).z());
                pts.push_back(test_point);
            }
            else {
                int num_points = int(dis / (step)) + 1;
                double dis_seg = dis / num_points;
                for (int k = 0; k != num_points; k++) {
                    geo_point_t current_pt(prev_wcp.x() + (k + 1.) / num_points * (test_point.x() - prev_wcp.x()),
                                           prev_wcp.y() + (k + 1.) / num_points * (test_point.y() - prev_wcp.y()),
                                           prev_wcp.z() + (k + 1.) / num_points * (test_point.z() - prev_wcp.z()));
                    pts.push_back(current_pt);
                }
            }
            prev_wcp = test_point;
        }

        // cloud.pts.resize(current_size + pts.size());
        // cloud_u.pts.resize(current_size + pts.size());
        // cloud_v.pts.resize(current_size + pts.size());
        // cloud_w.pts.resize(current_size + pts.size());
        // vec_index_cluster.resize(current_size + pts.size());
        // int i = 0;
        for (auto it = pts.begin(); it != pts.end(); it++) {
            m_clusters.push_back(cluster);
            m_blobs.push_back(nullptr);
            m_pc3d.add({(*it).x(), (*it).y(), (*it).z()});
            m_pc2d.add((*it));
            for (size_t plane = 0; plane < 3; ++plane) {
                m_winds[plane].push_back(2.4 * units::cm);
            }

            // vec_index_cluster.at(current_size + i) = cluster;
            // cloud.pts[current_size + i].x = (*it).x;
            // cloud.pts[current_size + i].y = (*it).y;
            // cloud.pts[current_size + i].z = (*it).z;
            // cloud.pts[current_size + i].index_u = 2.4 * units::cm;
            // cloud.pts[current_size + i].index_v = 2.4 * units::cm;
            // cloud.pts[current_size + i].index_w = 2.4 * units::cm;
            // cloud.pts[current_size + i].mcell = 0;
            // cloud.pts[current_size + i].index = current_size + i;
            // cloud_u.pts[current_size + i].x = (*it).x;
            // cloud_u.pts[current_size + i].y = cos(angle_u) * (*it).z - sin(angle_u) * (*it).y;
            // cloud_u.pts[current_size + i].index = current_size + i;
            // cloud_v.pts[current_size + i].x = (*it).x;
            // cloud_v.pts[current_size + i].y = cos(angle_v) * (*it).z - sin(angle_v) * (*it).y;
            // cloud_v.pts[current_size + i].index = current_size + i;
            // cloud_w.pts[current_size + i].x = (*it).x;
            // cloud_w.pts[current_size + i].y = cos(angle_w) * (*it).z - sin(angle_w) * (*it).y;
            // cloud_w.pts[current_size + i].index = current_size + i;

            // i++;
        }
        // if (pts.size() > 0) {
        //     index->addPoints(current_size, current_size + pts.size() - 1);
        //     index_u->addPoints(current_size, current_size + pts.size() - 1);
        //     index_v->addPoints(current_size, current_size + pts.size() - 1);
        //     index_w->addPoints(current_size, current_size + pts.size() - 1);
        // }
    }
}

void Facade::DynamicPointCloud::add_points(const Cluster* cluster, const geo_point_t& p_test,
                                                              const geo_point_t& dir_unmorm, const double range,
                                                              const double step, const double angle)
{
    size_t current_size = get_num_points();
    geo_point_t dir = dir_unmorm.norm();

    int num_points = int(range / (step)) + 1;
    double dis_seg = range / num_points;

    /// TODO: resize is faster, but needs more interface implementation
    for (int k = 0; k != num_points; k++) {
        // 13 cm  = 75 * sin(10/180.*3.1415926)
        double dis_cut =
            std::min(std::max(2.4 * units::cm, k * dis_seg * sin(angle / 180. * 3.1415926)), 13 * units::cm);
        m_clusters.push_back(cluster);
        m_blobs.push_back(nullptr);
        m_pc3d.add({p_test.x() + k * dir.x() * dis_seg, p_test.y() + k * dir.y() * dis_seg,
                    p_test.z() + k * dir.z() * dis_seg});
        m_winds[0].push_back(int(dis_cut));
        m_winds[1].push_back(int(dis_cut));
        m_winds[2].push_back(int(dis_cut));
        m_pc2d.add({p_test.x() + k * dir.x() * dis_seg, p_test.y() + k * dir.y() * dis_seg,
                    p_test.z() + k * dir.z() * dis_seg});
    }
}

std::vector<std::tuple<double, const Cluster*, size_t>> Facade::DynamicPointCloud::get_2d_points_info(
    const geo_point_t& p, const double radius, const int plane)
{
    std::vector<std::pair<size_t, double>> results = m_pc2d.get_closest_2d_index(p, radius, plane);
    std::vector<std::tuple<double, const Cluster*, size_t>> return_results;

    for (size_t i = 0; i != results.size(); i++) {
        return_results.push_back(std::make_tuple(results.at(i).second, m_clusters.at(results.at(i).first),
                                                 (size_t)results.at(i).first));
    }

    return return_results;
}



#include <boost/histogram.hpp>
#include <boost/histogram/algorithm/sum.hpp>
namespace bh = boost::histogram;
namespace bha = boost::histogram::algorithm;

// Example parameter calculating functions used by directional hough
// transforms.
static double theta_angle(const Vector& dir)
{
    const Vector Z(0, 0, 1);
    return acos(Z.dot(dir));
}
static double theta_cosine(const Vector& dir)
{
    const Vector Z(0, 0, 1);
    return Z.dot(dir);
}
static double phi_angle(const Vector& dir)
{
    const Vector X(1, 0, 0);
    const Vector Y(0, 1, 0);
    return atan2(Y.dot(dir), X.dot(dir));
}

std::pair<double, double> Facade::DynamicPointCloud::hough_transform(const geo_point_t& origin, const double dis) const
{
    std::vector<geo_point_t> pts;
    std::vector<const Blob*> blobs;
    auto results = m_pc3d.kd().radius(dis * dis, origin);
    for (const auto& [point_index, _] : results) {
        pts.push_back(m_pc3d.point(point_index));
        blobs.push_back(m_blobs.at(point_index));
    }

    constexpr double pi = 3.141592653589793;

    using direction_parameter_function_f = std::function<double(const Vector& dir)>;

    // Parameter axis 1 is some measure of theta angle (angle or cosine)
    const int nbins1 = 180;
    // param_space == costh_phi
    direction_parameter_function_f theta_param = theta_angle;
    double min1 = 0, max1 = pi;

    // Parameter axis 2 is only supported by phi angle
    const int nbins2 = 360;
    const double min2 = -pi;
    const double max2 = +pi;
    direction_parameter_function_f phi_param = phi_angle;

    auto hist = bh::make_histogram(bh::axis::regular<>(nbins1, min1, max1), bh::axis::regular<>(nbins2, min2, max2));

    for (size_t ind = 0; ind < blobs.size(); ++ind) {
        const auto* blob = blobs[ind];
        auto charge = blob->charge();
        // protection against the charge=0 case ...
        // if (charge == 0) charge = 1;
        if (charge <= 0) continue;

        const auto npoints = blob->npoints();
        const auto& pt = pts[ind];

        const Vector dir = (pt - origin).norm();
        const double r = (pt - origin).magnitude();

        const double p1 = theta_param(dir);
        const double p2 = phi_param(dir);
        if (r < 10 * units::cm) {
            hist(p1, p2, bh::weight(charge / npoints));
        }
        else {
            // hough->Fill(vec.Theta(), vec.Phi(), q * pow(10 * units::cm / r, 2));
            hist(p1, p2, bh::weight(charge / npoints * pow(10 * units::cm / r, 2)));
        }
    }

    auto indexed = bh::indexed(hist);
    auto it = std::max_element(indexed.begin(), indexed.end());
    const auto& cell = *it;
    return {cell.bin(0).center(), cell.bin(1).center()};
}


geo_point_t Facade::DynamicPointCloud::vhough_transform(const geo_point_t& origin, const double dis) const
{
    // TODO: only support theta_phi
    const auto [th, phi] = hough_transform(origin, dis);
    return {sin(th) * cos(phi), sin(th) * sin(phi), cos(th)};
}

// dirft = xorig + xsign * (time + m_time_offset) * m_drift_speed
double Facade::time2drift(const IAnodeFace::pointer anodeface, const double time_offset, const double drift_speed, double time) {
    const Pimpos* colpimpos = anodeface->planes()[2]->pimpos();
    double xsign = colpimpos->axis(0)[0];
    double xorig = anodeface->planes()[2]->wires().front()->center().x();
    const double drift = (time + time_offset)*drift_speed;
    /// TODO: how to determine xsign?
    return xorig + xsign*drift;
}

// time = (drift - xorig) / (xsign * m_drift_speed) - m_time_offset
double Facade::drift2time(const IAnodeFace::pointer anodeface, const double time_offset, const double drift_speed, double drift) {
    const Pimpos* colpimpos = anodeface->planes()[2]->pimpos();
    double xsign = colpimpos->axis(0)[0];
    double xorig = anodeface->planes()[2]->wires().front()->center().x();
    return (drift - xorig) / (xsign * drift_speed) - time_offset;
}

int Facade::point2wind(const geo_point_t& point, const double angle, const double pitch, const double center)
{
    // double y = cos(angles[pind]) * point[2] - sin(angles[pind]) * point[1];
    // y = mag * wind + center
    double y = cos(angle) * point[2] - sin(angle) * point[1];
    double wind = (y - center) / pitch;
    return std::round(wind);
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
