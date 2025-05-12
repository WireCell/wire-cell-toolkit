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
#include "WireCellIface/IDetectorVolumes.h"



#include "WCPQuickhull/QuickHull.h"

#include <string>


// extern int global_counter_get_closest_wcpoint;


// using namespace WireCell;  NO!  do not open up namespaces in header files!

namespace WireCell::Clus::Facade {

    struct DummyCache{};

    /// The Grouping/Cluster/Facade inherit from this to gain additional methods
    /// that are common to all three facade types.  The mixin itself needs to
    /// know its facade type and value but specifically does not include anything
    /// that requires parent or child types or values.
    ///
    /// It provides helper functions to deal with local PCs and an optional
    /// caching mechanism.  See comments on cache() and fill_cache() and
    /// clear_cache().  Note, using the cache mechanism does not preclude facade
    /// doing DIY caching.
    template<typename SelfType, typename CacheType=DummyCache>
    class Mixin {
        SelfType& self;
        std::string scalar_pc_name, ident_array_name;
        mutable std::unique_ptr<CacheType> m_cache;
    public:
        Mixin(SelfType& self, const std::string& scalar_pc_name, const std::string& ident_array_name = "ident")
            : self(self)
            , scalar_pc_name(scalar_pc_name)
            , ident_array_name(ident_array_name) {
            
        }

    protected:
        /// Facade cache management has three simple rules:
        ///
        /// Cache rule 1: The SelfType may call this to access a full and const cache.
        const CacheType& cache() const
        {
            if (! m_cache) {
                m_cache = std::make_unique<CacheType>();
                fill_cache(* const_cast<CacheType*>(m_cache.get()));
            }
            return *m_cache.get();
        }

        /// Cache rule 2:
        ///
        /// The SelfType overrides this method to fill an empty cache.  This is
        /// the only place where the cache object can be accessed by Self in
        /// mutable form.
        virtual void fill_cache(CacheType& cache) const {}
        
    public:
        /// Cache rule 3:
        ///
        /// The SelfType may override clear_cache(), for example to clear cached
        /// data not in the CacheType.  An override must then forward-call the
        /// Mixin::clear_cache().  The Mixin, the SelfType implementation and/or
        /// SelfType users may all call this method thought the goal is to make
        /// clear_cache() called in response to a tree notification.
        virtual void clear_cache() const
        {
            m_cache = nullptr;
        }

    public:

        /// Clear my node of all children nodes and purge my local PCs.
        /// Invalidates any cache.
        void clear()
        {
            // node level:
            self.node()->remove_children();
            // value level:
            self.local_pcs().clear();
            // facade cache level:
            clear_cache();
        }

        // Get the map from name to PC for all local PCs.
        WireCell::PointCloud::Tree::named_pointclouds_t& local_pcs()
        {
            return self.value().local_pcs();
        }
        const WireCell::PointCloud::Tree::named_pointclouds_t& local_pcs() const
        {
            return self.value().local_pcs();
        }

        // Return an "identifying number" from the "scalar" PC of the node.  As
        // with all "ident" values in WCT, there is no meaning ascribed to the
        // actual value (by WCT).  It is meant to refer to some external
        // identity.  If the scalar PC or the ident array are not found, the
        // default is returned.
        //
        // This is a special case method that merely delegates to get_scalar().
        int ident(int def = -1) const
        {
            return get_scalar<int>(ident_array_name, def);
        }

        // Set an ident number, delegating to set_scalar().
        void set_ident(int id)
        {
            set_scalar<int>(ident_array_name, id);
        }

        template <typename T = int>
        T get_element(const std::string& pcname, const std::string& aname, size_t index, T def = 0) const {
            const auto& lpcs = local_pcs();
            auto it = lpcs.find(pcname);
            if (it == lpcs.end()) {
                return def;
            }
            const auto arr = it->second.get(aname);
            if (!arr) {
                return def;
            }
            // std::cout << "test1 " << pcname << " " << aname << " " << index << " " << arr->template element<T>(index) << std::endl;
            return arr->template element<T>(index);
        }

        // Return a value from the scalar PC
        template <typename T = int>
        T get_scalar(const std::string& aname, T def = 0) const {
            return get_element(scalar_pc_name, aname, 0, def);
        }
        
        // Set a value on the scalar PC
        template <typename T = int>
        void set_scalar(const std::string& aname, T val = 0) {
            auto& lpcs = local_pcs();
            auto& cs = lpcs[scalar_pc_name]; // create if not existing
            auto arr = cs.get(aname);
            if (!arr) {
                cs.add(aname, PointCloud::Array({(T)val}));
                return;
            }
            arr->template element<T>(0) = (T)val;
        }

        /// A flag is a name that can be "set" on a facade.  It is simply an
        /// entry in the scalar PC.  Most imply, a flag is Boolean false if
        /// unset (not defined) or has value 0 and set if defined with non-zero
        /// value.  Non-boolean values are allowed.  The flag name has a prefix
        /// (default "flag_") to provide a namespace.
        void set_flag(const std::string& name, int value=1, const std::string& prefix="flag_") {
            set_scalar<int>(prefix + name, value);
        }

        /// Get the value of a flag.  If the flag is unset, return the
        /// default_value.  See set_flag().
        int get_flag(const std::string& name, int default_value=0, const std::string& prefix="flag_") const {
            return get_scalar<int>(prefix + name, default_value);
        }

        /// Get all set flag names with a given prefix.
        std::vector<std::string> flag_names(const std::string& prefix="flag_") const {
            std::vector<std::string> ret;
            const auto& spc = get_pc(scalar_pc_name);
            for (const auto& key : spc.keys()) {
                if (String::startswith(key, prefix)) {
                    ret.push_back(key.substr(0, prefix.size()));
                }
            }
            return ret;
        }

        // Any flag set on the other will be set on this.
        void flags_from(const SelfType& other, const std::string& prefix="flag_") {
            for (const auto& fname : other.flag_names(prefix)) {
                set_flag(fname, other.get_flag(fname, 0, prefix), prefix);
            }
        }


        bool has_pc(const std::string& pcname) const
        {
            static PointCloud::Dataset dummy;
            const auto& lpcs = local_pcs();
            auto it = lpcs.find(pcname);
            if (it == lpcs.end()) {
                return false;
            }
            return true;
        }

        // Const access to a local PC/Dataset.  If pcname is missing return
        // reference to an empty dataset.
        const PointCloud::Dataset& get_pc(const std::string& pcname) const
        {
            static PointCloud::Dataset dummy;
            const auto& lpcs = local_pcs();
            auto it = lpcs.find(pcname);
            if (it == lpcs.end()) {
                return dummy;
            }
            return it->second;
        }
        // Mutable access to a local PC/Dataset.  If pcname is missing, a new
        // dataset of that name will be created.
        PointCloud::Dataset& get_pc(const std::string& pcname)
        {
            static PointCloud::Dataset dummy;
            const auto& lpcs = local_pcs();
            return lpcs[pcname];
        }

        // Return true if this cluster has a PC array and PC of given names and type.
        template<typename ElementType=int>
        bool has_pcarray(const std::string& aname, const std::string& pcname) const {
            auto& lpc = local_pcs();
            auto lit = lpc.find(pcname);
            if (lit == lpc.end()) {
                return false;
            }

            auto arr = lit->second.get(aname);
            if (!arr) {
                return false;
            }
            return arr->template is_type<ElementType>();
        }

        // Return as a span an array named "aname" stored in clusters PC named
        // "pcname".  Returns default span if PC or array not found or there is
        // a type mismatch.  Note, span uses array data in place.
        template<typename ElementType=int>
        PointCloud::Array::span_t<ElementType>
        get_pcarray(const std::string& aname, const std::string& pcname) {

            auto& lpc = local_pcs();
            auto lit = lpc.find(pcname);
            if (lit == lpc.end()) {
                return {};
            }

            auto arr = lit->second.get(aname);
            if (!arr) {
                return {};
            }
            return arr->template elements<ElementType>();
        }
        template<typename ElementType=int>
        const PointCloud::Array::span_t<ElementType>
        get_pcarray(const std::string& aname, const std::string& pcname) const {

            auto& lpc = local_pcs();
            auto lit = lpc.find(pcname);
            if (lit == lpc.end()) {
                return {};
            }

            auto arr = lit->second.get(aname);
            if (!arr) {
                return {};
            }
            return arr->template elements<ElementType>();
        }

        // Store vector as an array named "aname" into this cluster's PC named "pcname".
        // Reminder, all arrays in a PC must have same major size.
        template<typename ElementType=int>
        void
        put_pcarray(const std::vector<ElementType>& vec,
                    const std::string& aname, const std::string& pcname) {

            auto &lpc = local_pcs();
            auto& pc = lpc[pcname];

            PointCloud::Array::shape_t shape = {vec.size()};

            auto arr = pc.get(aname);
            if (arr) {
                arr->template assign(vec.data(), shape, false);
            }
            else {
                pc.add(aname, PointCloud::Array(vec, shape, false));
            }
        }
    };


    using points_t = WireCell::PointCloud::Tree::Points;
    using node_t = WireCell::PointCloud::Tree::Points::node_t;
    using node_ptr = std::unique_ptr<node_t>;
    using geo_point_t = WireCell::Point;
    using geo_vector_t = WireCell::Vector;
    using geo_point2d_t = D3Vector<double>;
    using geo_vector2d_t = D3Vector<double>;

    class Cluster;
    class Blob;

    // map for face, plane to something
    /// TODO: face (ident? which?) -> plane (index) -> Dataset
    template<typename T> using mapfp_t = std::unordered_map<int, std::unordered_map<int, T>>;

    // FIXME: why define these out?
    // Haiwang: these are to make changing the underlying types easier and more consistent.
    // Probably need to move this to a better place.
    using float_t = double;
    using int_t = int;


    // AVOID DOING THIS in headers!!!  In this case it causes conflict between
    // boost::units and WireCell::Units in imp files that #include this one.
    //
    // If typing the namespace:: is too much, then one can do select "using
    // namespace::symbol".
    // 
    // using namespace boost;

    struct VertexProp {
        int index;
        // WCPointCloud<double>::WCPoint wcpoint;
        //  add pointer to merged cell
    };
    using EdgeProp = boost::property<boost::edge_weight_t, float>; 
    // struct EdgeProp {
    //     float dist;  // edge distance
    // };
    typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, VertexProp, EdgeProp> MCUGraph;
    typedef boost::graph_traits<MCUGraph>::vertex_descriptor vertex_descriptor;
    typedef boost::graph_traits<MCUGraph>::edge_descriptor edge_descriptor;

    // FIXME: refactor to vector<pitch>, etc?  or vector<TPCPlane> with ::pitch/::angle?
    struct TPCParams {
        int face{0};
        float_t pitch_u{3 * units::mm};
        float_t pitch_v{3 * units::mm};
        float_t pitch_w{3 * units::mm};
        float_t angle_u{1.0472};   // 60 degrees    uboone geometry ...
        float_t angle_v{-1.0472};  //-60 degrees   uboone geometry ...
        float_t angle_w{0};        // 0 degrees    uboone geometry ...
        float_t drift_speed{1.101 * units::mm / units::us};
        float_t tick{0.5 * units::us};           // 0.5 mm per tick
        float_t tick_drift{drift_speed * tick};  // tick * speed
        float_t time_offset{-1600 * units::us };
        int nticks_live_slice{4};

        float_t FV_xmin{1 * units::cm};
        float_t FV_xmax{255 * units::cm};
        float_t FV_ymin{-99.5 * units::cm};
        float_t FV_ymax{101.5 * units::cm};
        float_t FV_zmin{15 * units::cm};
        float_t FV_zmax{1022 * units::cm};
        float_t FV_xmin_margin{2 * units::cm};
        float_t FV_xmax_margin{2 * units::cm};
        float_t FV_ymin_margin{2.5 * units::cm};
        float_t FV_ymax_margin{2.5 * units::cm};
        float_t FV_zmin_margin{3 * units::cm};
        float_t FV_zmax_margin{3 * units::cm};
    };


    /// @brief  return local indices instead of global
    template <typename PCType1, typename PCType2>
    std::tuple<int, int, double> get_closest_points(const PCType1& one, const PCType2& two)
    {
        // improved algorithm ...
        double min_dis = 1e9;
        int p1_save = 0, p2_save = 0;
        
        // Sample points from first cloud at regular intervals
        int stride = std::max(1, (int)(one.points()[0].size() / 20)); // Sample ~20 points
        
        for(size_t i = 0; i < one.points()[0].size(); i += stride) {
            // Get K nearest neighbors from second cloud
            auto p1 = one.point(i);
            auto knn = two.kd().knn(5, p1); // Get 5 nearest neighbors
            
            // Refine search around these neighbors
            for(auto [idx2, dist] : knn) {
                auto p2 = two.point(idx2);
                
                // Local refinement by checking neighboring points
                int curr_idx1 = i;  // Keep track of current point index from first cloud
                std::tie(idx2, p2) = two.get_closest_wcpoint(p1);
                std::tie(curr_idx1, p1) = one.get_closest_wcpoint(p2);
                
                double dis = sqrt(pow(p1.x()-p2.x(),2) + pow(p1.y()-p2.y(),2) + pow(p1.z()-p2.z(),2));
                if(dis < min_dis) {
                    min_dis = dis;
                    p1_save = curr_idx1;  // Save the refined index from first cloud
                    p2_save = idx2;
                }
            }
        }
        
        return std::make_tuple(p1_save, p2_save, min_dis);
   
    }

    class Simple3DPointCloud {
       public:
        using nfkd_t = NFKDVec::Tree<double, NFKDVec::IndexDynamic>;
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
        nfkd_t& kd(bool rebuild=false);
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
        template <typename PCType>
        std::tuple<int, int, double> get_closest_points(const PCType& two) const
        {
            return Clus::Facade::get_closest_points(*this, two);
        
        }

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
        using coordinates_type = nfkd_t::coordinates_type;
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
        size_t get_num_points(const size_t plane) const { return m_points[plane][0].size(); }
        const nfkd_t& kd(const size_t plane, const bool rebuild=false) const;
        nfkd_t& kd(const size_t plane, const bool rebuild=false);
        std::pair<int, double> get_closest_2d_dis(const geo_point_t &p, size_t plane) const;
        std::vector<std::pair<size_t, double>> get_closest_2d_index_radius(const geo_point_t &p, const double radius, size_t plane) const;
        std::vector<std::pair<size_t, double>> get_closest_2d_index_knn(const geo_point_t &p, const int N, size_t plane) const;
       private:
        points_type m_points[3];
        mutable std::unique_ptr<nfkd_t> m_kd[3]; // lazy
        double angle_uvw[3];
    };
    std::ostream& operator<<(std::ostream& os, const Multi2DPointCloud& s3dpc);


    // class DynamicPointCloudLegacy {
    //    public:
    //     DynamicPointCloudLegacy(double angle_u, double angle_v, double angle_w);
    //     using points3d_type = Simple3DPointCloud::points_type;
    //     using points2d_type = Multi2DPointCloud::points_type;
    //     using point_type = std::vector<double>;
    //     inline size_t get_num_points() const { return m_pc3d.get_num_points(); }
    //     inline point_type point2d(const size_t plane, const size_t ind) const {
    //         return m_pc2d.point(plane, ind);
    //     }
    //     inline geo_point_t point3d(const size_t ind) const { return m_pc3d.point(ind); }

    //     // useful when hacking the winds with dist_cut
    //     inline int dist_cut(const size_t plane, const size_t ind) const { return m_winds[plane].at(ind); }

    //     /// @brief flag 0 points, flag 1 skeleton
    //     void add_points(const Cluster* cluster, const int flag=0, const double step = 0.6*units::cm); // flag 0 points, flag 1 scheleton

    //     /// @brief add points from p_test along dir with range and step
    //     /// @attention: the index_uvw is hacked to store the distance cut
    //     void add_points(const Cluster* cluster, const geo_point_t& p_test, const geo_point_t& dir_unmorm, const double range,
    //                     const double step, const double angle);

    //     /// @return: dist, Cluster, point_index
    //     std::vector<std::tuple<double, const Cluster*, size_t>> get_2d_points_info(const geo_point_t& p, const double radius,
    //                                                                          const int plane);
    //     /// @brief 
    //     std::tuple<double, const Cluster*, size_t> get_closest_2d_point_info(const geo_point_t& p, const int plane);

    //     std::pair<double, double> hough_transform(const geo_point_t& origin, const double dis) const;
    //     geo_point_t vhough_transform(const geo_point_t& origin, const double dis) const;
    //    private:
    //     Multi2DPointCloud m_pc2d;
    //     Simple3DPointCloud m_pc3d;
    //     std::vector<int> m_winds[3]; // u, v, w
    //     std::vector<const Cluster*> m_clusters;
    //     std::vector<const Blob*> m_blobs;
    // };

    void process_mst_deterministically(
            const boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
            boost::no_property, boost::property<boost::edge_weight_t, double>>& temp_graph,
            std::vector<std::vector<std::tuple<int,int,double>>>& index_index_dis,
            std::vector<std::vector<std::tuple<int,int,double>>>& index_index_dis_mst) ;

    double time2drift(const IAnodeFace::pointer anodeface, const double time_offset, const double drift_speed,
                      const double time);
    double drift2time(const IAnodeFace::pointer anodeface, const double time_offset, const double drift_speed,
                      const double drift);
    int point2wind(const geo_point_t& point, const double angle, const double pitch, const double center);
    double wind2point2dproj(const int wind, const double angle, const double pitch, const double center);

    WirePlaneId get_wireplaneid(const geo_point_t& point, const WirePlaneId& wpid1, const WirePlaneId& wpid2, IDetectorVolumes::pointer dv);
    WirePlaneId get_wireplaneid(const geo_point_t& p1, const WirePlaneId& wpid1, const geo_point_t& p2, const WirePlaneId& wpid2, IDetectorVolumes::pointer dv);

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



}  // namespace WireCell::Clus::Facade

#endif
