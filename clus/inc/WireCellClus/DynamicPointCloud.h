#ifndef WIRECELLCLUS_DYNAMICPOINTCLOUD_H
#define WIRECELLCLUS_DYNAMICPOINTCLOUD_H

#include "WireCellUtil/Graph.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Spdlog.h"
#include "WireCellUtil/Units.h"
// #include "WireCellUtil/D2Vector.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDetectorVolumes.h"

#include <array>

namespace WireCell::Clus::Facade {
    using points_t = PointCloud::Tree::Points;
    using node_t = PointCloud::Tree::Points::node_t;
    using node_ptr = std::unique_ptr<node_t>;
    using geo_point_t = WireCell::Point;
    using geo_vector_t = WireCell::Vector;
    using geo_point2d_t = D3Vector<double>;
    using geo_vector2d_t = D3Vector<double>;
    // using wpid_params_t = std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>>;

    class Cluster;
    class Blob;

    /// Columnar (SoA) point storage shared by the make_points_* builders
    /// (as the transfer format) and DynamicPointCloud (as the internal
    /// storage).  This replaces the former per-point DPCPoint struct
    /// whose five nested vectors cost ~a dozen heap allocations per
    /// point.  One entry per point in the scalar columns; the per-plane
    /// 2D projections are CSR-encoded: for point i and plane p the
    /// projection entries occupy [p2d_off[3*i+p], p2d_off[3*i+p+1]) of
    /// the flat p2d_* arrays.
    struct DPCBatch {
        // scalar columns, one entry per point
        std::vector<double> x, y, z;
        std::vector<int> wpid;
        std::vector<const Cluster*> cluster;
        std::vector<const Blob*> blob;
        std::vector<std::array<int, 3>> wind;          // bogus row when absent
        std::vector<std::array<double, 3>> dist_cut;   // bogus row when absent
        // 2D projections (CSR, see above)
        std::vector<size_t> p2d_off{0};
        std::vector<double> p2d_x, p2d_y;
        std::vector<int> p2d_wpid;

        size_t size() const { return x.size(); }
        bool empty() const { return x.empty(); }

        void reserve(size_t npts)
        {
            x.reserve(npts); y.reserve(npts); z.reserve(npts);
            wpid.reserve(npts);
            cluster.reserve(npts); blob.reserve(npts);
            wind.reserve(npts); dist_cut.reserve(npts);
            p2d_off.reserve(3*npts + 1);
            p2d_x.reserve(3*npts); p2d_y.reserve(3*npts);
            p2d_wpid.reserve(3*npts);
        }

        /// Builder protocol: push one point's scalars with add_point(),
        /// then for each of the 3 planes in order push that plane's
        /// projection entries with add_proj() and seal the plane with
        /// end_plane().  (3 end_plane() calls per add_point.)
        void add_point(double px, double py, double pz, int w,
                       const Cluster* c, const Blob* b,
                       const std::array<int, 3>& wd,
                       const std::array<double, 3>& dc)
        {
            x.push_back(px); y.push_back(py); z.push_back(pz);
            wpid.push_back(w);
            cluster.push_back(c); blob.push_back(b);
            wind.push_back(wd); dist_cut.push_back(dc);
        }
        void add_proj(double xx, double yy, int wp)
        {
            p2d_x.push_back(xx); p2d_y.push_back(yy); p2d_wpid.push_back(wp);
        }
        void end_plane() { p2d_off.push_back(p2d_x.size()); }

        /// CSR range of point i, plane p.
        std::pair<size_t, size_t> proj_range(size_t i, size_t p) const
        {
            return {p2d_off[3*i + p], p2d_off[3*i + p + 1]};
        }

        /// Append all rows of another batch (CSR offsets rebased).
        void append(const DPCBatch& other);
        /// Append the selected rows of another batch, in the given order.
        void append(const DPCBatch& other, const std::vector<size_t>& rows);
    };

    class DynamicPointCloud {
       public:
        using nfkd_t = NFKDVec::Tree<double, NFKDVec::IndexDynamic>;

        DynamicPointCloud(const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params)
          : m_wpid_params(wpid_params)
        {
        }
        ~DynamicPointCloud() = default;

        void add_points(const DPCBatch &points);
        // Rvalue overload: moves the batch in (column moves; for a fresh
        // cloud this is the whole storage in O(1)).
        void add_points(DPCBatch &&points);
        // Append all (or selected) points of another cloud.
        void add_points(const DynamicPointCloud &other);
        void add_points(const DynamicPointCloud &other, const std::vector<size_t> &rows);

        /// Columnar access to the stored points.
        const DPCBatch& points() const { return m_pts; }
        size_t npoints() const { return m_pts.size(); }
        geo_point_t point3d(size_t i) const { return {m_pts.x[i], m_pts.y[i], m_pts.z[i]}; }
        const Cluster* cluster(size_t i) const { return m_pts.cluster[i]; }
        const Blob* blob(size_t i) const { return m_pts.blob[i]; }
        int wpid(size_t i) const { return m_pts.wpid[i]; }
        const std::array<int, 3>& wind(size_t i) const { return m_pts.wind[i]; }
        double dist_cut(size_t i, int plane) const { return m_pts.dist_cut[i][plane]; }

        DynamicPointCloud::nfkd_t &kd3d() const;
        DynamicPointCloud::nfkd_t &kd2d(const int plane, const int face, const int apa) const;
        const std::unordered_map<size_t, size_t> &kd2d_l2g(const int plane, const int face, const int apa) const;
        const std::unordered_map<size_t, std::vector<size_t>> &kd2d_g2l(const int plane, const int face, const int apa) const;

        geo_point_t get_center_point_radius(const geo_point_t &p_test, const double radius) const;

        /// @brief: kd2d().radius(radius)
        /// @return: [dist, Cluster, global point_index]
        std::vector<std::tuple<double, const Cluster *, size_t>> get_2d_points_info(const geo_point_t &p,
                                                                                    const double radius,
                                                                                    const int plane, const int face,
                                                                                    const int apa) const;
        /// @brief: kd2d().knn(1)
        /// @brief: dist, Cluster, global point_index
        std::tuple<double, const Cluster *, size_t> get_closest_2d_point_info(const geo_point_t &p, const int plane,
                                                                              const int face, const int apa) const;

        /// @brief Like get_closest_2d_point_info but takes pre-projected (drift, wire_perp) coordinates
        /// directly from Grouping::convert_time_wire_2Dpoint(), bypassing the internal angle projection.
        /// Use this when the 2D coordinates are already in the wire-perpendicular space.
        std::tuple<double, const Cluster *, size_t> get_closest_2d_point_info_direct(
            double drift, double wire_perp, const int plane, const int face, const int apa) const;

        /// @brief Return (angle_u, angle_v, angle_w) for the given face/apa from the stored wpid_params.
        /// Returns {0,0,0} if the face/apa combination is not found.
        std::array<double, 3> get_angles(int face, int apa) const;

        /// @brief Read-only access to the full wpid_params map.
        const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>>& get_wpid_params() const {
            return m_wpid_params;
        }

        /// @brief Merge in any wpid_params entries from @p other not already present.
        /// Used when building a shower DPC that spans segments from multiple APAs/faces.
        void merge_wpid_params(const DynamicPointCloud& other) {
            for (const auto& [wpid, params] : other.m_wpid_params) {
                m_wpid_params.emplace(wpid, params);  // no-op if key already present
            }
        }


        std::pair<double, double> hough_transform(const geo_point_t &origin, const double dis) const;
        geo_point_t vhough_transform(const geo_point_t &origin, const double dis) const;

       private:
        // Index m_pts[original_size..end) into the 3D/2D k-d trees.
        void index_new_points(size_t original_size);

        // main data (columnar)
        DPCBatch m_pts;

        // for 3D only consider all apa for now
        mutable std::unique_ptr<nfkd_t> m_kd3d{nullptr};

        std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> m_wpid_params;

        // for 2D, wpid to kd
        mutable std::map<int, std::unique_ptr<nfkd_t>> m_kd2d;
        mutable std::unordered_map<int, std::unordered_map<size_t, size_t>> m_kd2d_index_l2g;
        mutable std::unordered_map<int, std::unordered_map<size_t, std::vector<size_t> >> m_kd2d_index_g2l;
    };

    DPCBatch
    make_points_cluster(const Cluster *cluster,
                        const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params, bool flag_wrap = false);

    DPCBatch
    make_points_cluster_steiner(const Cluster *cluster,
                        const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params, bool flag_wrap = false);

    DPCBatch make_points_cluster_skeleton(
        const Cluster *cluster, const IDetectorVolumes::pointer dv,
        const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params,
        const std::vector<size_t>& path_wcps,
        bool flag_wrap = false,
        const double step = 0.6 * units::cm);

    DPCBatch make_points_direct(const Cluster *cluster, const IDetectorVolumes::pointer dv, const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params, std::vector<std::pair<geo_point_t,WirePlaneId>>& points_info, bool flag_wrap = false);

    /// @brief add points from p_test along dir with range and step
    /// @attention: the index_uvw is hacked to store the distance cut
    /// seed_wpid: drift volume of the extrapolation seed point.  With a
    /// multi-volume wpid_params each synthetic point is bucketed into the
    /// volume containing it (dv->contained_by), falling back to seed_wpid
    /// (or the grouping's first wpid when seed_wpid is invalid) when the ray
    /// exits all sensitive volumes.
    DPCBatch make_points_linear_extrapolation(
        const Cluster *cluster, const geo_point_t &p_test, const geo_point_t &dir_unmorm, const double range,
        const double step, const double angle, const IDetectorVolumes::pointer dv,
        const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params,
        const WirePlaneId &seed_wpid = WirePlaneId(kUnknownLayer, -1, -1));

    /// Append the wrapped-wire 2D projections of one point into the
    /// batch: for each of the 3 planes, push that plane's entries and
    /// seal it with end_plane().  (The batch's add_point for this point
    /// must already have been called.)
    void fill_wrap_points(const Cluster *cluster, const geo_point_t &point, const WirePlaneId &wpid_point, DPCBatch& batch);

}  // namespace WireCell::Clus::Facade

#endif  // WIRECELLCLUS_DYNAMICPOINTCLOUD_H
