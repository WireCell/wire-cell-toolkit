/** A facade over a PC tree giving semantics to otherwise nodes.
 *
 */

#ifndef WIRECELL_CLUS_FACADECLUSTER
#define WIRECELL_CLUS_FACADECLUSTER

#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Spdlog.h"
#include "WireCellUtil/Graph.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"

#include "WireCellClus/Facade_Util.h"


// using namespace WireCell;  NO!  do not open up namespaces in header files!

namespace WireCell::PointCloud::Facade {
    class Blob;
    class Grouping;

    // Give a node "Cluster" semantics.  A cluster node's children are blob nodes.
    class Cluster : public NaryTree::FacadeParent<Blob, points_t> {
        // The expected scope.
        const Tree::Scope scope = {"3d", {"x", "y", "z"}};
        const Tree::Scope scope_wire_index = {"3d", {"uwire_index", "vwire_index", "wwire_index"}};
        Tree::Scope scope2ds[3] = {
            {"2dp0", {"x", "y"}},
            {"2dp1", {"x", "y"}},
            {"2dp2", {"x", "y"}}
        };

       public:
        Cluster() = default;
        virtual ~Cluster() {}

        // Return the grouping to which this cluster is a child.  May be nullptr.
        Grouping* grouping();
        const Grouping* grouping() const;

        // Get the scoped view for the "3d" point cloud (x,y,z)
        using sv3d_t = Tree::ScopedView<double>;
        const sv3d_t& sv3d() const;

        // Access the k-d tree for "3d" point cloud (x,y,z).
        // Note, this may trigger the underlying k-d build.
        using kd3d_t = sv3d_t::nfkd_t;
        const kd3d_t& kd3d() const;

        using kd_results_t = kd3d_t::results_type;
        // Perform a k-d tree radius query.  This radius is linear distance
        kd_results_t kd_radius(double radius_not_squared, const geo_point_t& query_point) const;
        // Perform a k-d tree NN query.
        kd_results_t kd_knn(int nnearest, const geo_point_t& query_point) const;

        std::vector<geo_point_t> kd_points(const kd_results_t& res);
        std::vector<geo_point_t> kd_points(const kd_results_t& res) const;

        // print all blob information
        void print_blobs_info() const;

        std::string dump() const;
        std::string dump_graph() const;

        // Get all blobs in k-d tree order.  This is different than children()
        // order and different that sort_blobs() order.
        std::vector<Blob*> kd_blobs();
        std::vector<const Blob*> kd_blobs() const;

        // Return the blob with the point at the given k-d tree point index.
        Blob* blob_with_point(size_t point_index);
        const Blob* blob_with_point(size_t point_index) const;

        // Return blobs from k-d tree result set.
        std::vector<Blob*> blobs_with_points(const kd_results_t& res);
        std::vector<const Blob*> blobs_with_points(const kd_results_t& res) const;

        // Return the 3D point at the k-d tree point index.  Calling this in a
        // tight loop should probably be avoided.  Instead get the full points() array.
        geo_point_t point3d(size_t point_index) const;
        // alias for point3d to match the Simple3DPointCloud interface
        geo_point_t point(size_t point_index) const;

        // Return vector is size 3 holding vectors of size npoints providing k-d tree coordinate points.
        using points_type = kd3d_t::points_type;
        const points_type& points() const;

        // Return charge-weighted average position of points of blobs within distance of point.
        geo_point_t calc_ave_pos(const geo_point_t& origin, const double dis) const;

        // Return blob containing the returned point that is closest to the given point.
        using point_blob_map_t = std::map<geo_point_t, const Blob*>;

        // WCP: get_closest_wcpoint
        std::pair<geo_point_t, const Blob*> get_closest_point_blob(const geo_point_t& point) const;
        /// WCP: get_closest_wcpoint: index, geo_point_t
        std::pair<size_t, geo_point_t> get_closest_wcpoint(const geo_point_t& p) const;

        // 
        size_t get_closest_point_index(const geo_point_t& point) const;

        // WCP: get_closest_dis
        double get_closest_dis(const geo_point_t& point) const;
        
        /// @return idx from self, idx from other and distance
        std::tuple<int, int, double> get_closest_points(const Cluster& other) const;

        // Return set of blobs each with a corresponding point.  The set
        // includes blobs with at least one point within the given radius of the
        // given point.  The point is one in the blob and that is closest to the
        // given point.
        //
        // Note: radius must provide a LINEAR distance measure.
        using const_blob_point_map_t = std::map<const Blob*, geo_point_t>;
        const_blob_point_map_t get_closest_blob(const geo_point_t& point, double radius) const;

        std::pair<geo_point_t, double> get_closest_point_along_vec(geo_point_t& p_test, geo_point_t dir,
                                                                   double test_dis, double dis_step, double angle_cut,
                                                                   double dis_cut) const;

        // Return the number of points in the k-d tree
        int npoints() const;

        // Number of points according to sum of Blob::nbpoints()
        // WCP: int get_num_points()
        size_t nbpoints() const;

        // Return the number of points within radius of the point.  Note, radius
        // is a LINEAR distance through the L2 metric is used internally.
        // WCP: int get_num_points(Point& p_test, double dis);
        int nnearby(const geo_point_t& point, double radius) const;

        // Return the number of points in the k-d tree partitioned into pair
        // (#forward,#backward) based on given direction of view from the given
        // point.
        // WCP: std::pair<int,int> get_num_points(Point& p, TVector3& dir, double dis);
        // if dis < 0, then no cut on distance
        std::pair<int, int> ndipole(const geo_point_t& point, const geo_point_t& dir, const double dis=-1) const;

        // Return the number of points with in the radius of the given point in
        // the k-d tree partitioned into pair (#forward,#backward) based on
        // given direction of view from the given point.
        //
        // Note: the radius is a LINEAR distance measure.
        // std::pair<int, int> nprojection(const geo_point_t& point, const geo_point_t& dir, double radius) const;

        // WCP: get_two_extreme_points
        // TODO: configurable dist cut?
        // 1, determines the most extreme points along the y, x, z axes
        // 2, calculates which pair of these points has the greatest distance between them
        // 3, adjusted using local averaging, calc_ave_pos
        std::pair<geo_point_t,geo_point_t> get_two_extreme_points() const;

        // Return the points at the extremes of the given Cartesian axis.  Default is Y-axis.
        //
        // Note: the two points are in DESCENDING order!
        // WCP: axis=1: get_highest_lowest_wcps, 2: get_front_back_wcps, 0: get_earliest_latest_wcps (reverse order)
        std::pair<geo_point_t, geo_point_t> get_highest_lowest_points(size_t axis = 1) const; 

        // Return the points at the extremes of the X-axis.
        // WCP: get_earliest_latest_wcps
        // Note: the two points are in ASCENDING order!
        std::pair<geo_point_t, geo_point_t> get_earliest_latest_points() const;

        // WCP: get_front_back_wcps
        std::pair<geo_point_t, geo_point_t> get_front_back_points() const;

        /// TODO: old_wcp and dir are used as local vars inside the function, make the IO more clear?
        geo_point_t get_furthest_wcpoint(geo_point_t old_wcp, geo_point_t dir, const double step = 5*units::cm, const int allowed_nstep = 12) const;

        /// @brief adjusts the positions of two points (start and end points)
        /// to be more in line with the overall point cloud.
        void adjust_wcpoints_parallel(size_t& start_idx, size_t& end_idx) const;
        
        /// WCP: Construct_skeleton
        bool construct_skeleton(const bool use_ctpc);

        /// section for 2D PC

        // Get the scoped view for the "3d" point cloud (x,y,z)
        using sv2d_t = Tree::ScopedView<double>;
        const sv2d_t& sv2d(const size_t plane) const;
        using kd2d_t = sv2d_t::nfkd_t;
        const kd2d_t& kd2d(const size_t plane) const;

        /// 
        std::vector<size_t> get_closest_2d_index(const geo_point_t& p, const double search_radius, const int plane) const;

        std::vector<const Blob*> is_connected(const Cluster& c, const int offset) const;

        // The Hough-based direction finder works in a 2D parameter space.  The
        // first dimension is associated with theta (angle w.r.t. Z-axis) and
        // can use an angle or a cosine measure.  Theta angle measure has a
        // non-uniform metric space, especially near the poles.  Perhaps this
        // bias is useful.  The cosine(theta) metric space is uniform.
        enum HoughParamSpace {
            costh_phi,  // (cos(theta), phi)
            theta_phi   // (theta, phi)
        };

        // Return parameter values characterizing the points within radius of
        // given point.
        //
        // Note: radius must provide a LINEAR distance measure.
        std::pair<double, double> hough_transform(const geo_point_t& point, const double radius,
                                                  HoughParamSpace param_space = HoughParamSpace::theta_phi,
                                                  std::shared_ptr<const Simple3DPointCloud> s3dpc = nullptr,
                                                  const std::vector<size_t>& global_indices = {}) const;

        // Call hough_transform() and transform result as to a directional vector representation.
        //
        // Note: radius must provide a LINEAR distance measure.
        // VHoughTrans in WCP
        geo_vector_t vhough_transform(const geo_point_t& point, const double radius,
                                      HoughParamSpace param_space = HoughParamSpace::theta_phi,
                                      std::shared_ptr<const Simple3DPointCloud> = nullptr,
                                      const std::vector<size_t>& global_indices = {}) const;

        // Return a quasi geometric size of the cluster based on its transverse
        // extents in each view and in time.
        double get_length() const;

        // Return blob at the front of the time blob map.  Raises ValueError if cluster is empty.
        const Blob* get_first_blob() const;

        // Return blob at the back of the time blob map.  Raises ValueError if cluster is empty.
        const Blob* get_last_blob() const;

        // number of unique slice times, i.e. time_blob_map().size()
        size_t get_num_time_slices() const;

        // Return a value representing the content of this cluster.
        size_t hash() const;

        // Check facade consistency between blob view and k-d tree view.
        bool sanity(Log::logptr_t log = nullptr) const;

        inline MCUGraph* get_graph() { return m_graph.get(); }
        inline const MCUGraph* get_graph() const { return m_graph.get(); }
        void Create_graph(const bool use_ctpc = true) const;

        /// @brief edges inside blobs and between overlapping blobs
        /// @attention has distance-based cuts
        void Establish_close_connected_graph() const;
        /// @attention some distance-based cuts
        void Connect_graph(const bool use_ctpc = false) const;

        ///
        void dijkstra_shortest_paths(const size_t pt_idx, const bool use_ctpc = true) const;

        ///
        void cal_shortest_path(const size_t dest_wcp_index) const;


        ///
        inline const std::list<size_t>& get_path_wcps() const { return m_path_wcps; }
        
        // TODO: relying on scoped_view to do the caching?
        using wire_indices_t = std::vector<std::vector<int_t>>;
        const wire_indices_t& wire_indices() const;

        std::vector<geo_point_t> get_hull() const;

        geo_point_t get_center() const;
        geo_vector_t get_pca_axis(int axis) const;
        double get_pca_value(int axis) const;

        // start slice index (tick number) to blob facade pointer can be
        // duplicated, example usage:
        // https://github.com/HaiwangYu/learn-cpp/blob/main/test-multimap.cxx
        // WCP: get_time_cells_set_map
        using BlobSet = std::set<const Blob*>;
        using time_blob_map_t = std::map<int, BlobSet>;
        const time_blob_map_t& time_blob_map() const;

        /// @brief break a cluster if it crosses the boundary
        /// @return if breaking, return broken clusters, if not, return empty map
        std::unordered_map<int, Cluster*> examine_x_boundary(const double low_limit = -1*units::cm, const double high_limit = 257*units::cm);

        /// @brief get_mcell_indices in WCP
        /// TODO: currently return copy, return a const reference?
        std::vector<int> get_blob_indices(const Blob*) const;

        /// @brief to assess whether a given point (p_test) in a cluster is a vertex, or endpoint, based on asymmetry and occupancy criteria.
        /// @note p_test will be updated
        bool judge_vertex(geo_point_t& p_test, const double asy_cut = 1. / 3., const double occupied_cut = 0.85);

       private:
        mutable time_blob_map_t m_time_blob_map;  // lazy, do not access directly.

        mutable std::map<const Blob*, std::vector<int>> m_map_mcell_indices; // lazy, do not access directly.


        // Cached and lazily calculated in get_length().
        // Getting a new node invalidates by setting to 0.
        mutable double m_length{0};
        // Cached and lazily calculated in npoints()
        mutable int m_npoints{0};

        void Calc_PCA() const;
        mutable bool m_pca_calculated{false};
        // lazy, do not access directly.
        mutable geo_point_t m_center;
        mutable geo_vector_t m_pca_axis[3];
        mutable double m_pca_values[3];

        // m_graph
        mutable std::unique_ptr<MCUGraph> m_graph;
        // create things for Dijkstra
        mutable std::vector<vertex_descriptor> m_parents;
        mutable std::vector<int> m_distances;
        mutable int m_source_pt_index{-1};
        mutable std::list<size_t> m_path_wcps;
        mutable std::list<const Blob*> m_path_mcells;

       public:  // made public only for debugging
        // Return the number of unique wires or ticks.
        std::tuple<int, int, int, int> get_uvwt_range() const;
        std::tuple<int, int, int, int> get_uvwt_min() const;
        std::tuple<int, int, int, int> get_uvwt_max() const;
    };
    std::ostream& operator<<(std::ostream& os, const Cluster& cluster);

    // Return true if a is less than b.  May be used as 3rd arg in std::sort to
    // get ascending order.  For descending, pass to sort() rbegin()/rend()
    // instead of begin()/end()..
    bool cluster_less(const Cluster* a, const Cluster* b);
    // Apply standard sort to put clusters in descending order.
    void sort_clusters(std::vector<const Cluster*>& clusters);
    void sort_clusters(std::vector<Cluster*>& clusters);


}  // namespace WireCell::PointCloud::Facade

template <> struct fmt::formatter<WireCell::PointCloud::Facade::Cluster> : fmt::ostream_formatter {};

#endif
