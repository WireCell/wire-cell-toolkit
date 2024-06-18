/** A facade over a PC tree giving semantics to otherwise nodes.
 *
 */

#ifndef WIRECELL_CLUS_FACADEUTIL
#define WIRECELL_CLUS_FACADEUTIL

#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Spdlog.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"

#include "WireCellUtil/Graph.h"


// using namespace WireCell;  NO!  do not open up namespaces in header files!

namespace WireCell::PointCloud::Facade {
    using points_t = Tree::Points;
    using node_t = Tree::Points::node_t;
    using node_ptr = std::unique_ptr<node_t>;
    using geo_point_t = WireCell::Point;
    using geo_vector_t = WireCell::Vector;

    // map for face, plane to something
    /// TODO: face (ident? which?) -> plane (index) -> Dataset
    template<typename T> using mapfp_t = std::unordered_map<int, std::unordered_map<int, T>>;

    // FIXME: why define these out?
    // Haiwang: these are to make changing the underlying types easier and more consistent.
    // Probably need to move this to a better place.
    using float_t = double;
    using int_t = int;

    // FIXME: refactor to vector<pitch>, etc?  or vector<TPCPlane> with ::pitch/::angle?
    struct TPCParams {
        float_t pitch_u{3 * units::mm};
        float_t pitch_v{3 * units::mm};
        float_t pitch_w{3 * units::mm};
        float_t angle_u{1.0472};   // 60 degrees    uboone geometry ...
        float_t angle_v{-1.0472};  //-60 degrees   uboone geometry ...
        float_t angle_w{0};        // 0 degrees    uboone geometry ...
        float_t drift_speed{1.101 * units::mm / units::us};
        float_t tick{0.5 * units::us};           // 0.5 mm per tick
        float_t tick_drift{drift_speed * tick};  // tick * speed
        float_t time_offset{-1600 * units::us};
    };

    double time2drift(const IAnodeFace::pointer anodeface, const double time_offset, const double drift_speed,
                      const double time);
    double drift2time(const IAnodeFace::pointer anodeface, const double time_offset, const double drift_speed,
                      const double drift);
    int point2wind(const geo_point_t& point, const double angle, const double pitch, const double center);

    // fixme: why do we inline these?
    inline double cal_proj_angle_diff(const geo_vector_t& dir1, const geo_vector_t& dir2, double plane_angle)
    {
        geo_vector_t temp_dir1;
        geo_vector_t temp_dir2;

        temp_dir1.set(dir1.x(), 0, -sin(plane_angle) * dir1.y() + cos(plane_angle) * dir1.z());
        temp_dir2.set(dir2.x(), 0, -sin(plane_angle) * dir2.y() + cos(plane_angle) * dir2.z());

        return temp_dir1.angle(temp_dir2);
    }

    inline bool is_angle_consistent(const geo_vector_t& dir1, const geo_vector_t& dir2, bool same_direction,
                                    double angle_cut, double uplane_angle, double vplane_angle, double wplane_angle,
                                    int num_cut = 2)
    {
        double angle_u = cal_proj_angle_diff(dir1, dir2, uplane_angle);
        double angle_v = cal_proj_angle_diff(dir1, dir2, vplane_angle);
        double angle_w = cal_proj_angle_diff(dir1, dir2, wplane_angle);
        int num = 0;
        // input is degrees ...
        angle_cut *= 3.1415926 / 180.;

        if (same_direction) {
            if (angle_u <= angle_cut) num++;
            if (angle_v <= angle_cut) num++;
            if (angle_w <= angle_cut) num++;
        }
        else {
            if ((3.1415926 - angle_u) <= angle_cut) num++;
            if ((3.1415926 - angle_v) <= angle_cut) num++;
            if ((3.1415926 - angle_w) <= angle_cut) num++;
        }

        if (num >= num_cut) return true;
        return false;
    }

}  // namespace WireCell::PointCloud::Facade

#endif
