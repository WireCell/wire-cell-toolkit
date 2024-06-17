/** A facade over a PC tree giving semantics to otherwise nodes.
 *
 */

#ifndef WIRECELLIMG_POINTCLOUDFACADE
#define WIRECELLIMG_POINTCLOUDFACADE

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

    class Cluster;

    /// Give a node "Blob" semantics
    class Blob : public NaryTree::Facade<points_t> {
       public:
        Blob() = default;
        virtual ~Blob() {}

        // Return the cluster to which this blob is a child.  May be nullptr.
        Cluster* cluster();
        const Cluster* cluster() const;

        geo_point_t center_pos() const;

        bool overlap_fast(const Blob& b, const int offset) const;

        float_t charge() const { return charge_; }
        float_t center_x() const { return center_x_; }
        float_t center_y() const { return center_y_; }
        float_t center_z() const { return center_z_; }
        int_t npoints() const { return npoints_; }

        // units are number of ticks
        int_t slice_index_min() const { return slice_index_min_; }
        int_t slice_index_max() const { return slice_index_max_; }

        int_t u_wire_index_min() const { return u_wire_index_min_; }
        int_t u_wire_index_max() const { return u_wire_index_max_; }
        int_t v_wire_index_min() const { return v_wire_index_min_; }
        int_t v_wire_index_max() const { return v_wire_index_max_; }
        int_t w_wire_index_min() const { return w_wire_index_min_; }
        int_t w_wire_index_max() const { return w_wire_index_max_; }

        int get_max_wire_interval() const { return m_max_wire_interval;}
        int get_min_wire_interval() const { return m_min_wire_interval;}
        int get_max_wire_type() const { return m_max_wire_type;}
        int get_min_wire_type() const { return m_min_wire_type;}

        // Return a value representing the content of this blob.
        size_t hash() const;

        // Return the scope points.
        std::vector<geo_point_t> points() const;
        size_t nbpoints() const;

        // Check facade consistency
        bool sanity(Log::logptr_t log = nullptr) const;

       private:
        float_t charge_{0};
        float_t center_x_{0};
        float_t center_y_{0};
        float_t center_z_{0};
        int_t npoints_{0};

        int_t slice_index_min_{0};  // unit: tick
        int_t slice_index_max_{0};

        int_t u_wire_index_min_{0};
        int_t u_wire_index_max_{0};
        int_t v_wire_index_min_{0};
        int_t v_wire_index_max_{0};
        int_t w_wire_index_min_{0};
        int_t w_wire_index_max_{0};

        // FIXME: dummy values for now
        int m_max_wire_interval{2}; 
        int m_min_wire_interval{1};
        int m_max_wire_type{2}; // 0: u, 1: v, 2: w
        int m_min_wire_type{0}; // 0: u, 1: v, 2: w

       protected:
        // Receive notification when this facade is created on a node.
        virtual void on_construct(node_type* node);
    };
    std::ostream& operator<<(std::ostream& os, const Blob& blob);

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
                                                  HoughParamSpace param_space = HoughParamSpace::theta_phi) const;

        // Call hough_transform() and transform result as to a directional vector representation.
        //
        // Note: radius must provide a LINEAR distance measure.
        geo_vector_t vhough_transform(const geo_point_t& point, const double radius,
                                      HoughParamSpace param_space = HoughParamSpace::theta_phi) const;

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

       public:  // made public only for debugging
        // Return the number of unique wires or ticks.
        std::tuple<int, int, int, int> get_uvwt_range() const;
        std::tuple<int, int, int, int> get_uvwt_min() const;
        std::tuple<int, int, int, int> get_uvwt_max() const;
    };
    std::ostream& operator<<(std::ostream& os, const Cluster& cluster);

    // Give a node "Grouping" semantics.  A grouping node's children are cluster
    // nodes that are related in some way.
    class Grouping : public NaryTree::FacadeParent<Cluster, points_t> {

        TPCParams m_tp{};  // use default value by default.
        /// TODO: replace TPCParams with this in the future?
        IAnodePlane::pointer m_anode{nullptr};

       public:
        // MUST call this sometimes after construction if non-default value needed.
        void set_params(const TPCParams& tp) { m_tp = tp; }
        const TPCParams& get_params() const { return m_tp; }
        std::tuple<double, double, double> wire_angles() const { return {m_tp.angle_u, m_tp.angle_v, m_tp.angle_w}; }

        void set_anode(const IAnodePlane::pointer anode) { m_anode = anode; }

        // Return a value representing the content of this grouping.
        size_t hash() const;

        std::map<int, std::pair<double, double>>& get_dead_winds(const int face, const int pind) const
        {
            // make one if not exist
            return m_dead_winds[face][pind];
        }
        using sv2d_t = Tree::ScopedView<float_t>;
        using kd2d_t = sv2d_t::nfkd_t;
        using kd_results_t = kd2d_t::results_type;

        const kd2d_t& kd2d(const int face, const int pind) const;

        const mapfp_t<double>& proj_centers() const; // lazy, do not access directly.
        const mapfp_t<double>& pitch_mags() const;   // lazy, do not access directly.

        bool is_good_point(const geo_point_t& point, const int face, const double radius = 0.6 * units::cm, const int ch_range = 1,
                           const int allowed_bad = 1) const;

        /// @brief
        /// @param point
        /// @param radius
        /// @param face
        /// @param pind plane index
        /// @return
        kd_results_t get_closest_points(const geo_point_t& point, const double radius, const int face, int pind) const;

        /// @brief true if the point is within the dead region, [wind+ch_range, wind-ch_range] and [xmin, xmax]
        bool get_closest_dead_chs(const geo_point_t& point, const int ch_range, const int face, int pind) const;

        /// @brief convert_3Dpoint_time_ch
        std::tuple<int, int> convert_3Dpoint_time_ch(const geo_point_t& point, const int face, const int pind) const;

       private:
        void fill_proj_centers_pitch_mags() const;
        mutable mapfp_t<double> m_proj_centers;
        mutable mapfp_t<double> m_pitch_mags;
        mutable mapfp_t< std::map<int, std::pair<double, double>> > m_dead_winds;
    };
    std::ostream& operator<<(std::ostream& os, const Grouping& grouping);

    double time2drift(const IAnodeFace::pointer anodeface, const double time_offset, const double drift_speed,
                      const double time);
    double drift2time(const IAnodeFace::pointer anodeface, const double time_offset, const double drift_speed,
                      const double drift);
    int point2wind(const geo_point_t& point, const double angle, const double pitch, const double center);

    // Return true if a is less than b.  May be used as 3rd arg in std::sort to
    // get ascending order.  For descending, pass to sort() rbegin()/rend()
    // instead of begin()/end()..
    bool blob_less(const Blob* a, const Blob* b);
    // Apply standard sort to put blobs in descending order.
    void sort_blobs(std::vector<const Blob*>& blobs);
    void sort_blobs(std::vector<Blob*>& blobs);

    // Return true if a is less than b.  May be used as 3rd arg in std::sort to
    // get ascending order.  For descending, pass to sort() rbegin()/rend()
    // instead of begin()/end()..
    bool cluster_less(const Cluster* a, const Cluster* b);
    // Apply standard sort to put clusters in descending order.
    void sort_clusters(std::vector<const Cluster*>& clusters);
    void sort_clusters(std::vector<Cluster*>& clusters);

    // Dump the grouping to a string eg for debugging.  Level 0 dumps info about
    // just the grouping, 1 includes info about the clusters and 2 includes info
    // about the blobs.
    std::string dump(const Grouping& grouping, int level = 0);

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

template <> struct fmt::formatter<WireCell::PointCloud::Facade::Blob> : fmt::ostream_formatter {};
template <> struct fmt::formatter<WireCell::PointCloud::Facade::Cluster> : fmt::ostream_formatter {};
template <> struct fmt::formatter<WireCell::PointCloud::Facade::Grouping> : fmt::ostream_formatter {};

#endif
