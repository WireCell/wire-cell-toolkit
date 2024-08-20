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
#include "WireCellUtil/Graph.h"
// #include "WireCellUtil/D2Vector.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"


#include "WCPQuickhull/QuickHull.h"

// extern int global_counter_get_closest_wcpoint;


// using namespace WireCell;  NO!  do not open up namespaces in header files!

namespace WireCell::PointCloud::Facade {
    using points_t = Tree::Points;
    using node_t = Tree::Points::node_t;
    using node_ptr = std::unique_ptr<node_t>;
    using geo_point_t = WireCell::Point;
    using geo_vector_t = WireCell::Vector;
    using geo_point2d_t = D3Vector<double>;
    using geo_vector2d_t = D3Vector<double>;

    class Cluster;

    // map for face, plane to something
    /// TODO: face (ident? which?) -> plane (index) -> Dataset
    template<typename T> using mapfp_t = std::unordered_map<int, std::unordered_map<int, T>>;

    // FIXME: why define these out?
    // Haiwang: these are to make changing the underlying types easier and more consistent.
    // Probably need to move this to a better place.
    using float_t = double;
    using int_t = int;


    using namespace boost;
    struct VertexProp {
        int index;
        // WCPointCloud<double>::WCPoint wcpoint;
        //  add pointer to merged cell
    };
    struct EdgeProp {
        float dist;  // edge distance
    };
    typedef adjacency_list<vecS, vecS, undirectedS, VertexProp, EdgeProp> MCUGraph;
    typedef graph_traits<MCUGraph>::vertex_descriptor vertex_descriptor;
    typedef graph_traits<MCUGraph>::edge_descriptor edge_descriptor;

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


    class Multi2DPointCloud {
       public:
        Multi2DPointCloud(double angle_u, double angle_v, double angle_w);
        using nfkd_t = NFKDVec::Tree<double, NFKDVec::IndexDynamic>;
        // these derived types do not depend on static/dynamic
        using points_type = nfkd_t::points_type;
        using results_type = nfkd_t::results_type;
        using point_type = std::vector<double>;
        inline const points_type& points(const size_t plane) const { return m_points[plane]; }
        inline points_type& points(const size_t plane) { return m_points[plane]; }
        inline point_type point(const size_t plane, const size_t ind) const {
            if (ind >= m_points[plane][0].size()) {
                raise<IndexError>("point index %d out of range %d", ind, m_points[0].size());
            }
            return {m_points[plane][0][ind], m_points[plane][1][ind]};
        }
        void add(const geo_point_t& new_pt);
        size_t get_num_points() const { return m_points[0][0].size(); }
        const nfkd_t& kd(const size_t plane, const bool rebuild=false) const;
        std::pair<int, double> get_closest_2d_dis(const geo_point_t &p, size_t plane) const;
       private:
        points_type m_points[3];
        mutable std::unique_ptr<nfkd_t> m_kd[3]; // lazy
        double angle_uvw[3];
    };
    std::ostream& operator<<(std::ostream& os, const Multi2DPointCloud& s3dpc);


    class DynamicPointCloud {
       public:
        DynamicPointCloud(double angle_u, double angle_v, double angle_w);
        using points3d_type = Simple3DPointCloud::points_type;
        using points2d_type = Multi2DPointCloud::points_type;
        using point_type = std::vector<double>;
        inline size_t get_num_points() const { return m_pc3d.get_num_points(); }
        inline point_type point2d(const size_t plane, const size_t ind) const {
            return m_pc2d.point(plane, ind);
        }
        inline geo_point_t point3d(const size_t ind) const { return m_pc3d.point(ind); }

        /// @brief add points from p_test along dir with range and step
        void add_points(const Cluster* cluster, const geo_point_t& p_test, const geo_point_t& dir_unmorm, const double range,
                        const double step, const double angle);

        /// @return: dist, Cluster, point_index
        std::vector<std::tuple<double, Cluster*, size_t>> get_2d_points_info(const geo_point_t& p, const double radius,
                                                                             const int plane);

       private:
        Multi2DPointCloud m_pc2d;
        Simple3DPointCloud m_pc3d;
        std::vector<int> m_winds[3]; // u, v, w
        std::vector<const Cluster*> m_index_cluster;
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
