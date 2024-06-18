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
