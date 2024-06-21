#include "WireCellClus/Facade.h"
#include <boost/container_hash/hash.hpp>

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Facade;
// using WireCell::PointCloud::Dataset;
using namespace WireCell::PointCloud::Tree;  // for "Points" node value type
// using WireCell::PointCloud::Tree::named_pointclouds_t;

#include "WireCellUtil/Logging.h"
using spdlog::debug;

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
    const auto knn_res = kd().knn(1, p);
    if (knn_res.size() != 1) {
        raise<ValueError>("no points found");
    }
    const auto ind = knn_res[0].first;
    geo_point_t pt = {points()[0][ind], points()[1][ind], points()[2][ind]};
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
    int p1_save;
    int p2_save;
    double min_dis = 1e9;

    int prev_index1 = -1;
    int prev_index2 = -1;
    while (p1_index != prev_index1 || p2_index != prev_index2) {
        prev_index1 = p1_index;
        prev_index2 = p2_index;
        std::tie(p2_index, p2) = other.get_closest_wcpoint(p1);
        std::tie(p1_index, p1) = get_closest_wcpoint(p2);
    }
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
