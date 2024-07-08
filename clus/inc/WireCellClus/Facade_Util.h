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


#include "WCPQuickhull/QuickHull.h"

// extern int global_counter_get_closest_wcpoint;


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
        int nticks_live_slice{4};
    };

    class Simple3DPointCloud {
       public:
        using nfkd_t = NFKDVec::Tree<double, NFKDVec::IndexStatic>;
        // these derived types do not depend on static/dynamic
        using points_type = nfkd_t::points_type;
        using results_type = nfkd_t::results_type;
        using point_type = std::vector<double>;
        inline const points_type& points() const { return m_points; }
        inline points_type& points() { return m_points; }
        inline geo_point_t point(const size_t ind) const {
            if (ind >= m_points[0].size()) {
                raise<IndexError>("point index %d out of range %d", ind, m_points[0].size());
            }
            return {m_points[0][ind], m_points[1][ind], m_points[2][ind]};
        }
        void add(const point_type& new_pt);
        size_t get_num_points() const { return m_points[0].size(); }
        const nfkd_t& kd(bool rebuild=false) const;
        results_type get_closest_index(const geo_point_t& p, const size_t N) const;
        /// @return index, geo_point_t
        std::pair<size_t, geo_point_t> get_closest_wcpoint(const geo_point_t& p) const;

        /// @param p_test1 is the point to start from
        /// @param dir is the direction to search along
        /// @param test_dis is the distance to search
        /// @param dis_step is the step size to sample along dir
        /// @param angle_cut in degrees
        /// @param dis_cut is the maximum distance to search
        std::pair<int, double> get_closest_point_along_vec(
            const geo_point_t& p_test1,
            const geo_point_t& dir,
            double test_dis,
            double dis_step,
            double angle_cut,
            double dis_cut) const;

        /// @brief  return local indices instead of global
        std::tuple<int, int, double> get_closest_points(const Simple3DPointCloud& other) const;
       private:
        points_type m_points{3};
        mutable std::unique_ptr<nfkd_t> m_kd{nullptr}; // lazy
    };
    std::ostream& operator<<(std::ostream& os, const Simple3DPointCloud& s3dpc);

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
