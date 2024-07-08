#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Grouping.h"
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

std::ostream& Facade::operator<<(std::ostream& os, const Facade::Grouping& grouping)
{
    os << "<Grouping [" << (void*) &grouping << "]:" << " nclusters=" << grouping.nchildren() << ">";
    return os;
}

std::string Facade::dump(const Facade::Grouping& grouping, int level)
{
    std::stringstream ss;

    ss << grouping;
    if (level == 0) {
        return ss.str();
    }
    ss << "\n";
    size_t nc = 0;
    for (const auto* cluster : grouping.children()) {
        ss << nc++ << "\t" << *cluster << "\n";
        if (level == 1) {
            continue;
        }
        size_t nb = 0;
        for (const auto* blob : cluster->children()) {
            ss << nb++ << "\t\t" << *blob << "\n";
        }
    }
    return ss.str();
}

size_t Grouping::hash() const
{
    std::size_t h = 0;
    boost::hash_combine(h, m_tp.pitch_u);
    boost::hash_combine(h, m_tp.pitch_v);
    boost::hash_combine(h, m_tp.pitch_w);
    boost::hash_combine(h, m_tp.angle_u);
    boost::hash_combine(h, m_tp.angle_v);
    boost::hash_combine(h, m_tp.angle_w);
    boost::hash_combine(h, m_tp.tick_drift);
    auto clusters = children();  // copy vector
    sort_clusters(clusters);
    for (const Cluster* cluster : clusters) {
        boost::hash_combine(h, cluster->hash());
    }
    return h;
}

const Grouping::kd2d_t& Grouping::kd2d(const int face, const int pind) const
{
    const auto sname = String::format("ctpc_f%dp%d", face, pind);
    Tree::Scope scope = {sname, {"x", "y"}, 1};
    const auto& sv = m_node->value.scoped_view(scope);
    // std::cout << "sname: " << sname << " npoints: " << sv.kd().npoints() << std::endl;
    return sv.kd();
}

void Grouping::fill_proj_centers_pitch_mags() const
{
    const int ndummy_layers = 2;
    if (!m_anode) {
        raise<ValueError>("anode is null");
    }
    for (const auto& face : m_anode->faces()) {
        const auto& coords = face->raygrid();
        // skip dummy layers so the vector matches 0, 1, 2 plane order
        for (int layer=ndummy_layers; layer<coords.nlayers(); ++layer) {
            const auto& pitch_dir = coords.pitch_dirs()[layer];
            const auto& center = coords.centers()[layer];
            double proj_center = center.dot(pitch_dir);
            m_proj_centers[face->which()][layer-ndummy_layers] = proj_center;
            m_pitch_mags[face->which()][layer-ndummy_layers] = coords.pitch_mags()[layer];
        }
    }
}

const Facade::mapfp_t<double>& Grouping::proj_centers() const
{
    if (!m_proj_centers.empty()) return m_proj_centers;
    fill_proj_centers_pitch_mags();
    return m_proj_centers;
}

const Facade::mapfp_t<double>& Grouping::pitch_mags() const
{
    if (!m_pitch_mags.empty()) return m_pitch_mags;
    fill_proj_centers_pitch_mags();
    return m_pitch_mags;
}


bool Grouping::is_good_point(const geo_point_t& point, const int face, double radius, int ch_range, int allowed_bad) const {
    const int nplanes = 3;
    int matched_planes = 0;
    for (int pind = 0; pind < nplanes; ++pind) {
        if (get_closest_points(point, radius, face, pind).size() > 0) {
            matched_planes++;
        } else if (get_closest_dead_chs(point, ch_range, face, pind)) {
            matched_planes++;
        }
    }
    // std::cout << "matched_planes: " << matched_planes << std::endl;
    if (matched_planes >= nplanes - allowed_bad) {
        return true;
    }
    return false;
}

Grouping::kd_results_t Grouping::get_closest_points(const geo_point_t& point, const double radius, const int face,
                                                    int pind) const
{
    double x = point[0];
    const auto [angle_u,angle_v,angle_w] = wire_angles();
    std::vector<double> angles = {angle_u, angle_v, angle_w};
    double y = cos(angles[pind]) * point[2] - sin(angles[pind]) * point[1];
    const auto& skd = kd2d(face, pind);
    return skd.radius<std::vector<double>>(radius * radius, {x, y});
}

bool Grouping::get_closest_dead_chs(const geo_point_t& point, const int ch_range, const int face, int pind) const {
    const auto [tind, wind] = convert_3Dpoint_time_ch(point, face, pind);
    const auto& ch2xrange = get_dead_winds(face, pind);
    for (int ch = wind - ch_range; ch <= wind + ch_range; ++ch) {
        if (ch2xrange.find(ch) ==  ch2xrange.end()) continue;
        const auto [xmin, xmax] = ch2xrange.at(ch);
        if (point[0] >= xmin && point[0] <= xmax) {
            std::cout << "ch " << ch << " x " << point[0] << " xmin " << xmin << " xmax " << xmax << std::endl;
            return true;
        }
    }
    return false;
}

std::tuple<int, int> Grouping::convert_3Dpoint_time_ch(const geo_point_t& point, const int face, const int pind) const {
    if (m_anode == nullptr) {
        raise<ValueError>("Anode is null");
    }
    const auto& iface = m_anode->face(face);
    if (iface == nullptr) {
        raise<ValueError>("Face is null");
    }

    const auto [angle_u,angle_v,angle_w] = wire_angles();
    std::vector<double> angles = {angle_u, angle_v, angle_w};
    const double angle = angles[pind];
    const double pitch = pitch_mags().at(face).at(pind);
    const double center = proj_centers().at(face).at(pind);
    const int wind = point2wind(point, angle, pitch, center);

    const auto params = get_params();
    const double time = drift2time(iface, params.time_offset, params.drift_speed, point[0]);
    const int tind = std::round(time / params.tick);

    return {tind, wind};
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
