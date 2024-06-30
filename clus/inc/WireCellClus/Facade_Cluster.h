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

    struct VertexProp {
        int index;
        // WCPointCloud<double>::WCPoint wcpoint;
        //  add pointer to merged cell
    };
    struct EdgeProp {
        float dist;  // edge distance
    };

    using namespace boost;
    typedef adjacency_list<vecS, vecS, undirectedS, VertexProp, EdgeProp> MCUGraph;
    typedef graph_traits<MCUGraph>::vertex_descriptor vertex_descriptor;
    typedef graph_traits<MCUGraph>::edge_descriptor edge_descriptor;

    // Give a node "Cluster" semantics.  A cluster node's children are blob nodes.
    class Cluster : public NaryTree::FacadeParent<Blob, points_t> {
        // The expected scope.
        const Tree::Scope scope = {"3d", {"x", "y", "z"}};
        const Tree::Scope scope_wire_index = {"3d", {"uwire_index", "vwire_index", "wwire_index"}};

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

        // Return vector is size 3 holding vectors of size npoints providing k-d tree coordinate points.
        using points_type = kd3d_t::points_type;
        const points_type& points() const;

        // Return charge-weighted average position of points of blobs within distance of point.
        geo_point_t calc_ave_pos(const geo_point_t& origin, const double dis) const;

        // Return blob containing the returned point that is closest to the given point.
        using point_blob_map_t = std::map<geo_point_t, const Blob*>;
        std::pair<geo_point_t, const Blob*> get_closest_point_blob(const geo_point_t& point) const;

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
        size_t nbpoints() const;

        // Return the number of points within radius of the point.  Note, radius
        // is a LINEAR distance through the L2 metric is used internally.
        int nnearby(const geo_point_t& point, double radius) const;

        // Return the number of points in the k-d tree partitioned into pair
        // (#forward,#backward) based on given direction of view from the given
        // point.
        std::pair<int, int> ndipole(const geo_point_t& point, const geo_point_t& dir) const;

        // Return the number of points with in the radius of the given point in
        // the k-d tree partitioned into pair (#forward,#backward) based on
        // given direction of view from the given point.
        //
        // Note: the radius is a LINEAR distance measure.
        // std::pair<int, int> nprojection(const geo_point_t& point, const geo_point_t& dir, double radius) const;

        // Return the points at the extremes of the X-axis.
        //
        // Note: the two points are in ASCENDING order!
        std::pair<geo_point_t, geo_point_t> get_earliest_latest_points() const;

        // Return the points at the extremes of the given Cartesian axis.  Default is Y-axis.
        //
        // Note: the two points are in DESCENDING order!
        std::pair<geo_point_t, geo_point_t> get_highest_lowest_points(size_t axis = 1) const;

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

        // Return a value representing the content of this cluster.
        size_t hash() const;

        // Check facade consistency between blob view and k-d tree view.
        bool sanity(Log::logptr_t log = nullptr) const;
        
        // FIXME: move to private after debugging
        // graph
        MCUGraph* graph;
        void Create_graph(const bool use_ctpc = false);

        /// @brief edges inside blobs and between overlapping blobs
        void Establish_close_connected_graph();
        void Connect_graph(const bool use_ctpc = false);
        
        // TODO: relying on scoped_view to do the caching?
        using wire_indices_t = std::vector<std::vector<int_t>>;
        const wire_indices_t& wire_indices();

        std::vector<geo_point_t> get_hull() const;

        geo_point_t get_center() const;
        geo_vector_t get_pca_axis(int axis) const;
        double get_pca_value(int axis) const;

       private:
        // start slice index (tick number) to blob facade pointer can be
        // duplicated, example usage:
        // https://github.com/HaiwangYu/learn-cpp/blob/main/test-multimap.cxx

        using BlobSet = std::set<const Blob*>;
        using time_blob_map_t = std::map<int, BlobSet>;
        const time_blob_map_t& time_blob_map() const;
        mutable time_blob_map_t m_time_blob_map;  // lazy, do not access directly.

        std::map<const Blob*, std::vector<int>> m_map_mcell_indices; // lazy, do not access directly.
        std::vector<int> get_blob_indices(const Blob*);


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
