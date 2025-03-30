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

namespace WireCell::PointCloud::Facade {
    using points_t = Tree::Points;
    using node_t = Tree::Points::node_t;
    using node_ptr = std::unique_ptr<node_t>;
    using geo_point_t = WireCell::Point;
    using geo_vector_t = WireCell::Vector;
    using geo_point2d_t = D3Vector<double>;
    using geo_vector2d_t = D3Vector<double>;
    // using wpid_params_t = std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>>;

    class Cluster;
    class Blob;

    class DynamicPointCloud {
       public:
        struct DPCPoint {
            double x, y, z;
            // this could be a duplicated from Blob
            // if wpid then x/y_2d length must be 3
            int wpid;
            const Cluster *cluster;
            const Blob *blob;
            std::vector<std::vector<double>> x_2d;
            std::vector<std::vector<double>> y_2d;
            std::vector<int> wpid_2d;
            std::vector<int> wind;      // length 3 or 0
            std::vector<int> dist_cut;  // length 3 or 0
        };
        using nfkd_t = NFKDVec::Tree<double, NFKDVec::IndexDynamic>;

        DynamicPointCloud(const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params)
          : m_wpid_params(wpid_params)
        {
        }
        ~DynamicPointCloud() = default;

        void add_points(const std::vector<DPCPoint> &points);
        const std::vector<DPCPoint>& get_points() const { return m_points; }

        DynamicPointCloud::nfkd_t &kd3d() const;
        DynamicPointCloud::nfkd_t &kd2d(const int plane, const int face, const int apa) const;
        const std::unordered_map<size_t, size_t> &kd2d_l2g(const int plane, const int face, const int apa) const;
        const std::unordered_map<size_t, std::vector<size_t>> &kd2d_g2l(const int plane, const int face, const int apa) const;

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

        std::pair<double, double> hough_transform(const geo_point_t &origin, const double dis) const;
        geo_point_t vhough_transform(const geo_point_t &origin, const double dis) const;

       private:
        // main data
        std::vector<DPCPoint> m_points;

        // for 3D only consider all apa for now
        mutable std::unique_ptr<nfkd_t> m_kd3d{nullptr};

        std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> m_wpid_params;

        // for 2D, wpid to kd
        mutable std::map<int, std::unique_ptr<nfkd_t>> m_kd2d;
        std::unordered_map<int, std::unordered_map<size_t, size_t>> m_kd2d_index_l2g;
        std::unordered_map<int, std::unordered_map<size_t, std::vector<size_t> >> m_kd2d_index_g2l;
    };

    std::vector<DynamicPointCloud::DPCPoint>
    make_points_cluster(const Cluster *cluster,
                        const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params, bool flag_wrap = false);

    std::vector<DynamicPointCloud::DPCPoint> make_points_cluster_skeleton(
        const Cluster *cluster, const IDetectorVolumes::pointer dv,
        const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params, bool flag_wrap = false,
        const double step = 0.6 * units::cm);

    /// @brief add points from p_test along dir with range and step
    /// @attention: the index_uvw is hacked to store the distance cut
    std::vector<DynamicPointCloud::DPCPoint> make_points_linear_extrapolation(
        const Cluster *cluster, const geo_point_t &p_test, const geo_point_t &dir_unmorm, const double range,
        const double step, const double angle, const IDetectorVolumes::pointer dv,
        const std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> &wpid_params);

    void fill_wrap_points(const Cluster *cluster, const geo_point_t &point, int wpid_point, std::vector<std::vector<double>>& p_x, std::vector<std::vector<double>>& p_y, std::vector<int>& p_wpid);

}  // namespace WireCell::PointCloud::Facade

#endif  // WIRECELLCLUS_DYNAMICPOINTCLOUD_H
