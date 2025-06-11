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

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IDetectorVolumes.h"

#include "WireCellClus/Facade_Mixins.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_ClusterCache.h"
#include "WireCellClus/IPCTransform.h"
#include "WireCellClus/Graphs.h"

#include <functional>

namespace WireCell::Clus::Facade {

    using IPCTransformSet = Clus::IPCTransformSet;
    using namespace WireCell::PointCloud;

    class Blob;
    class Grouping;

    // Give a node "Cluster" semantics.  A cluster node's children are blob nodes.
    class Cluster : public NaryTree::FacadeParent<Blob, points_t>
                  , public Mixins::Cached<Cluster, ClusterCache>
                  , public Mixins::Graphs
    {
      public:
        Cluster() : Mixins::Cached<Cluster, ClusterCache>(*this, "cluster_scalar") {}
        virtual ~Cluster() {}

        // return raw pc information ...
        void set_default_scope(const Tree::Scope& scope);
        const Tree::Scope& get_default_scope() const {return m_default_scope;}
        const Tree::Scope& get_raw_scope() const {return m_scope_3d_raw;}

        // set, get scope filter ...
        void set_scope_filter(const Tree::Scope& scope, bool flag=true);
        bool get_scope_filter(const Tree::Scope& scope) const;

        void set_scope_transform(const Tree::Scope& scope, const std::string& transform_name);

        // If no scope given, will use default scope. 
        std::string get_scope_transform(Tree::Scope scope = {}) const;

        const Tree::Scope& get_scope(const std::string& scope_name) const;

        /// Set other's default scope, filter and transform to this.  See also from().
        void default_scope_from(const Cluster& other);

        /// Set on this various meta information from the other including
        /// default scope and any flags in the default flags_ prefix.
        ///
        /// This may be appropriate to call when this cluster is made from
        /// others such as in a cluster merge.
        ///
        /// See set_flag(), get_flag(), flag_names() and flag_from() provided by
        /// base class.
        void from(const Cluster& other);

        double get_cluster_t0() const;
        void set_cluster_t0(double cluster_t0);

        // scopes_from() and from()

        void set_cluster_id(int id);
        int get_cluster_id() const;

        /// @param correction_name: T0Correction
        std::vector<int> add_corrected_points(const Clus::IPCTransformSet::pointer pcts,
                                              const std::string &correction_name);


        /// Return the grouping to which this cluster is a child.  May be nullptr.
        Grouping* grouping();
        const Grouping* grouping() const;

        /// Order is synchronized with children().
        std::vector<WireCell::WirePlaneId> wpids_blob() const;

        /// return the wpid given a point ...
        WirePlaneId wpid(const geo_point_t& point) const;

        /// Get an arbitrary scoped view.  Make sure type T matches the type of
        /// the scope's coords arrays.
        template <typename T=double>
        const Tree::ScopedView<T>& sv(const Tree::Scope& sc) const
        {
            return m_node->value.scoped_view<T>(sc);
        }

        /// Default scoped view is a view of the default scope.
        template <typename T=double>
        const Tree::ScopedView<T>& sv() const
        {
            return m_node->value.scoped_view<T>(m_default_scope);
        }


        /// Return a vector of values from given array name "key" that spans the
        /// points in nodes of the 3D RAW scoped view.  Note, the "key" name
        /// need not be in the RAW scope.
        template <typename T=double>
        const std::vector<T> points_property(const std::string& key) const
        {
            return sv().template flat_vector<T>(key);
        }

        /// Return a vector of values from given array name "key" that spans the
        /// points in nodes of the given scoped.  Note, the "key" name need not
        /// be in the RAW scope.  If the pcname is not given, the scope pcname
        /// is used.
        template <typename T=double, typename ST=double>
        const std::vector<T> points_property(const std::string& key,
                                             const Tree::Scope& scope,
                                             std::string pcname="") const
        {
            return sv<ST>(scope).template flat_vector<T>(key, pcname);
        }

        // Get a 3D scopeed view scoped view for the "3d" point cloud (x,y,z)
        using sv3d_t = Tree::ScopedView<double>;
        const sv3d_t& sv3d() const;
        const sv3d_t& sv3d_raw() const;

        // Access the k-d tree for "3d" point cloud (x,y,z).
        // Note, this may trigger the underlying k-d build.
        using kd3d_t = sv3d_t::nfkd_t;
        const kd3d_t& kd3d() const;
        const kd3d_t& kd() const;

        const kd3d_t& kd3d_raw() const;

    public:


        using kd_results_t = kd3d_t::results_type;
        // Perform a k-d tree radius query.  This radius is linear distance
        kd_results_t kd_radius(double radius_not_squared, const geo_point_t& query_point) const;
        // Perform a k-d tree NN query.
        kd_results_t kd_knn(int nnearest, const geo_point_t& query_point) const;

        /// Return 
        std::vector<geo_point_t> kd_points(const kd_results_t& res);
        std::vector<geo_point_t> kd_points(const kd_results_t& res) const;

        // std::vector<geo_point_t> kd_points_raw(const kd_results_t& res);
        // std::vector<geo_point_t> kd_points_raw(const kd_results_t& res) const;

        // print all blob information
        void print_blobs_info() const;

        std::string dump() const;
        // std::string dump_graph() const;

        // Get all blobs in k-d tree order.  This is different than children()
        // order and different that sort_blobs() order.
        std::vector<Blob*> kd_blobs();
        std::vector<const Blob*> kd_blobs() const;
        // // Return the number of blobs from the k-d tree
        // size_t nkd_blobs() const;

        // Return the blob with the point at the given k-d tree point index.
        Blob* blob_with_point(size_t point_index);
        const Blob* blob_with_point(size_t point_index) const;

        // Return blobs from k-d tree result set.
        std::vector<Blob*> blobs_with_points(const kd_results_t& res);
        std::vector<const Blob*> blobs_with_points(const kd_results_t& res) const;

        // Return the 3D point at the k-d tree point index.  Calling this in a
        // tight loop should probably be avoided.  Instead get the full points() array.
        geo_point_t point3d(size_t point_index) const;
        geo_point_t point3d_raw(size_t point_index) const;

        // return WirePlaneId for an index ...
        WirePlaneId wire_plane_id(size_t point_index) const;

        // Return vector is size 3 holding vectors of size npoints providing k-d tree coordinate points.
        using points_type = kd3d_t::points_type;
        // Return points in a scope in point order.  If no scope is given, use default_scope.
        const points_type& points() const;
        const points_type& points_raw() const;

        // Return charge-weighted average position of points of blobs within distance of point.
        geo_point_t calc_ave_pos(const geo_point_t& origin, const double dis) const;
        geo_point_t calc_ave_pos(const geo_point_t& origin, int N) const;


        // In the public section of the Cluster class:
        geo_vector_t calc_dir(const geo_point_t& p_test, const geo_point_t& p, double dis) const;

        // Return blob containing the returned point that is closest to the given point.
        using point_blob_map_t = std::map<geo_point_t, const Blob*>;

        // WCP: get_closest_wcpoint
        std::pair<geo_point_t, const Blob*> get_closest_point_blob(const geo_point_t& point) const;
        /// WCP: get_closest_wcpoint: index, geo_point_t
        std::pair<size_t, geo_point_t> get_closest_wcpoint(const geo_point_t& p) const;

        /// PCType needs to have point() and get_closest_wcpoint() methods
        template <typename PCType>
        std::tuple<int, int, double> get_closest_points(const PCType& two) const
        {
            return get_closest_points(*this, two);
        }

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
        const_blob_point_map_t get_closest_blob(const geo_point_t& point, int N) const;


        std::pair<geo_point_t, double> get_closest_point_along_vec(geo_point_t& p_test, geo_point_t dir,
                                                                   double test_dis, double dis_step, double angle_cut,
                                                                   double dis_cut) const;

        // Return the number of points in the k-d tree
        int npoints() const;

        // Number of points according to sum of Blob::nbpoints()
        // WCP: int get_num_points()
        // size_t nbpoints() const;

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

        std::pair<geo_point_t, geo_point_t> get_main_axis_points() const;

        /// TODO: old_wcp and dir are used as local vars inside the function, make the IO more clear?
        geo_point_t get_furthest_wcpoint(geo_point_t old_wcp, geo_point_t dir, const double step = 5*units::cm, const int allowed_nstep = 12) const;

        /// @brief adjusts the positions of two points (start and end points)
        /// to be more in line with the overall point cloud.
        void adjust_wcpoints_parallel(size_t& start_idx, size_t& end_idx) const;
        

        /// section for 2D PC

        // Get the scoped view for the "3d" point cloud (x,y,z)
        using sv2d_t = Tree::ScopedView<double>;
        /// @param plane 0, 1, 2
        /// @param wpid currently provides the apa and face
        const sv2d_t& sv2d(const int apa, const int face, const size_t plane) const;
        using kd2d_t = sv2d_t::nfkd_t;
        const kd2d_t& kd2d(const int apa, const int face, const size_t plane) const;

        /// 
        std::vector<size_t> get_closest_2d_index(const geo_point_t& p, const double search_radius, const int apa, const int face, const int plane) const;

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

        // // Return blob at the front of the time blob map.  Raises ValueError if cluster is empty.
        // const Blob* get_first_blob() const;

        // // Return blob at the back of the time blob map.  Raises ValueError if cluster is empty.
        // const Blob* get_last_blob() const;

        // number of unique slice times, i.e. time_blob_map().size()
        // size_t get_num_time_slices() const;

        // Return a value representing the content of this cluster.
        size_t hash() const;

        // Check facade consistency between blob view and k-d tree view.
        bool sanity(Log::logptr_t log = nullptr) const;


        ///
        /// Blob-level connected components
        /// 
        /// Return a connected components array that is aligned with the blob
        /// node children list.  Each element gives a "group number" identifying
        /// a connected subgraph in which the corresponding blob resides.  The
        /// special group number -1 indicates the corresponding blob does not
        /// contribute points to the graph.
        ///
        /// The "relaxed" graph (ne' "overclustering protection") is used.  See
        /// graph_algorithms() family of methods.
        ///
        /// Note, this method used to be called "examine_graph()".
        ///
        std::vector<int> connected_blobs(IDetectorVolumes::pointer dv, IPCTransformSet::pointer pcts) const;
        
        ///
        /// Graph algorithms hold a graph and are cached
        ///
        /// Get a graph algorithms by its "flavor" and throw KeyError if not
        /// found.
        /// 
        /// Note, as a special case, the default graph flavor "basic" can and
        /// will be produced on the fly.
        const WireCell::Clus::Graphs::Weighted::GraphAlgorithms& 
        graph_algorithms(const std::string& flavor = "basic") const;

        /// Get and construct if needed a GA for a graph of a known flavor and
        /// that uses detector information in its construction.  Known flavors
        /// include:
        /// 
        /// - "ctpc" :: likely used for Djikstra's shortest paths
        /// - "relaxed" :: likely used for connected blobs
        /// 
        /// If the flavor is not in this known set, KeyError is
        /// thrown.
        ///
        /// Note, "relaxed" used to be known as "overclustering protection").
        const WireCell::Clus::Graphs::Weighted::GraphAlgorithms& 
        graph_algorithms(const std::string& flavor,
                         IDetectorVolumes::pointer dv, 
                         IPCTransformSet::pointer pcts) const;



        // Return 3D points for given indices in the 3d PC.
        std::vector<geo_point_t> indices_to_points(const std::vector<size_t>& path_indices) const;

        // void organize_points_path_vec(std::vector<geo_point_t>& path_points, double low_dis_limit) const;
        // void organize_path_points(std::vector<geo_point_t>& path_points, double low_dis_limit) const;

        // TODO: relying on scoped_view to do the caching?
        using wire_indices_t = std::vector<std::vector<int_t>>;
        const wire_indices_t& wire_indices() const;

        std::vector<geo_point_t> get_hull() const;

        // Return PCA calculated on blob children sample points
        // PCA has attributes: {center,axis,values}
        using PCA = ClusterCache::PCA;
        PCA& get_pca() const;

        // start slice index (tick number) to blob facade pointer can be
        // duplicated, example usage:
        // https://github.com/HaiwangYu/learn-cpp/blob/main/test-multimap.cxx
        // WCP: get_time_cells_set_map
        using time_blob_map_t = ClusterCache::time_blob_map_t;
        const time_blob_map_t& time_blob_map() const;
   
        /// @brief Determine if a cluster may be separated due to crossing the boundary.
        /// @return connected components array or empty if separation is not warranted.
        std::vector<int>
        examine_x_boundary(const double low_limit = -1*units::cm, const double high_limit = 257*units::cm);

        /// @brief get_mcell_indices 
        /// WCP: get_cell_times_set_map
        /// TODO: currently return copy, return a const reference?
        std::vector<int> get_blob_indices(const Blob*) const;

        // Return the number of unique wires or ticks.
        std::map<WirePlaneId, std::tuple<int, int, int, int> > get_uvwt_range() const;
        std::tuple<int, int, int, int> get_uvwt_min() const;
        std::tuple<int, int, int, int> get_uvwt_max() const;


        /// @brief to assess whether a given point (p_test) in a cluster is a vertex, or endpoint, based on asymmetry and occupancy criteria.
        /// @note p_test will be updated
        bool judge_vertex(geo_point_t& p_test, const IDetectorVolumes::pointer dv, const double asy_cut = 1. / 3., const double occupied_cut = 0.85);

        class Flash {
            friend class Cluster;
            bool m_valid{false};
            double m_time{0}, m_value{0};
            int m_ident{-1}, m_type{-1};
            std::vector<double> m_times, m_values, m_errors;
        public:

            /// A "false" means there was not "flash" PC array and all values
            /// are invalid.  A "true" does not guarantee all values are valid.
            explicit operator bool() const { return m_valid;}

            /// Any "singular" methods are about the flash itself.

            /// Get the time of the flash.  
            double time() const { return m_time; }

            /// Get the measurement of the flash
            double value() const {return m_value; }

            /// The ID of the flash
            int ident() const { return m_ident; }

            /// The type of the flash.
            int type() const { return m_type; }

            /// Any "plural" methods return per-optical-detector quantities.
            /// They will be empty() if the "light" and "flashlight" arrays are
            /// missing.  These vectors have the same size.

            // keep these return-by-value.

            /// Times of individual optical detector readouts.
            std::vector<double> times() const { return m_times; }
            /// Measurement values from optical detectors.
            std::vector<double> values() const { return m_values; }
            /// Measurement uncertainty from optical detectors.
            std::vector<double> errors() const { return m_errors; }
        };
        // Return a flash.  If there is none, it will hold default values.
        Flash get_flash() const;

       private:

        // default scope for all points with raw x,y,z as coordinates
        std::map<std::string, Tree::Scope> m_scopes = {
            {"scope_3d_raw", {"3d", {"x", "y", "z"}}}
        };
        // FIXME: shoud we remove this in the future?
        const Tree::Scope& m_scope_3d_raw = m_scopes.at("scope_3d_raw");
        const Tree::Scope m_scope_wire_index = {"3d", {"uwire_index", "vwire_index", "wwire_index"}};
        std::string m_scope2ds_prefix[3] = {"2dp0", "2dp1", "2dp2"};
        Tree::Scope m_default_scope = m_scope_3d_raw;
        std::map<size_t, bool> m_map_scope_filter={{m_scope_3d_raw.hash(), true}};
        std::map<size_t, std::string> m_map_scope_transform={{m_scope_3d_raw.hash(), "Unity"}};


       protected:

        //
        // Caching.
        //
        // See the ClusterCache struct in Facade_ClusterCache.h.
        //
        // DO NOT PUT BARE CACHE ITEMS DIRECTLY IN THE Cluster class.
        //
        virtual void fill_cache(ClusterCache& cache) const;
    };                          // Cluster
    std::ostream& operator<<(std::ostream& os, const Cluster& cluster);



    // Return true if a is less than b.  May be used as 3rd arg in std::sort to
    // get ascending order.  For descending, pass to sort() rbegin()/rend()
    // instead of begin()/end()..
    bool cluster_less(const Cluster* a, const Cluster* b);
    // Apply standard sort to put clusters in descending order.
    void sort_clusters(std::vector<const Cluster*>& clusters);
    void sort_clusters(std::vector<Cluster*>& clusters);

    std::map<WirePlaneId, std::tuple<int, int, int, int> > get_uvwt_range(const Cluster* cluster, const std::vector<int>& b2id, const int id);
    double get_length(const Cluster* cluster, const std::vector<int>& b2id, const int id);

    struct cluster_less_functor {
        bool operator()(const Cluster* a, const Cluster* b) const { return cluster_less(a, b); }
    };

   struct ComponentInfo {
        int component_id;
        std::vector<size_t> vertex_indices;
        size_t min_vertex;
        double total_weight;  // Add this member

        ComponentInfo(int id) 
            : component_id(id)
            , min_vertex(std::numeric_limits<size_t>::max())
            , total_weight(0.0)  // Initialize in constructor
        {}

        void add_vertex(size_t vertex_idx) {
            vertex_indices.push_back(vertex_idx);
            min_vertex = std::min(min_vertex, vertex_idx);
        }
    };

}  // namespace WireCell::Clus::Facade

template <> struct fmt::formatter<WireCell::Clus::Facade::Cluster> : fmt::ostream_formatter {};

#endif
