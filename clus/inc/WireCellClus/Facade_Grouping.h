/** A facade over a PC tree giving semantics to otherwise nodes.
 *
 */

#ifndef WIRECELL_CLUS_FACADEGROUPING
#define WIRECELL_CLUS_FACADEGROUPING

#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Spdlog.h"
#include "WireCellUtil/Graph.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IDetectorVolumes.h"

#include "WireCellClus/Facade_Util.h"


// using namespace WireCell;  NO!  do not open up namespaces in header files!

namespace WireCell::PointCloud::Facade {
    class Cluster;

    struct GroupingCache {

        mapfp_t<double> proj_centers;
        mapfp_t<double> pitch_mags;

        // what wpids the grouping has.
        std::set<WireCell::WirePlaneId> wpids;

        // #381 if you give a crap about dead_winds.  

    };

    // Give a node "Grouping" semantics.  A grouping node's children are cluster
    // nodes that are related in some way.
    class Grouping : public NaryTree::FacadeParent<Cluster, points_t>, public Mixin<Grouping, GroupingCache> {

        std::map<int, IAnodePlane::pointer> m_anodes;
        IDetectorVolumes::pointer m_dv{nullptr};

        /// TODO: remove these in the future
        IAnodePlane::pointer m_anode{nullptr};
        TPCParams m_tp{};  // use default value by default.

       public:

        Grouping() : Mixin<Grouping, GroupingCache>(*this, "grouping_scalar") {}

        // TODO: remove this in the future
        void set_params(const TPCParams& tp) { m_tp = tp; }
        void set_params(const WireCell::Configuration& cfg);
        const TPCParams& get_params() const { return m_tp; }
        std::tuple<double, double, double> wire_angles() const { return {m_tp.angle_u, m_tp.angle_v, m_tp.angle_w}; }

        /// TODO: remove this in the future
        // void set_anode(const IAnodePlane::pointer anode) { m_anode = anode; }
        const IAnodePlane::pointer get_anode() const { return m_anode; }

        void set_anodes(const std::vector<IAnodePlane::pointer>& anodes);
        const IAnodePlane::pointer get_anode(const int ident) const;

        void set_detector_volumes(const IDetectorVolumes::pointer dv) { m_dv = dv; }
        const IDetectorVolumes::pointer get_detector_volumes() const { return m_dv; }

        // Return a value representing the content of this grouping.
        size_t hash() const;

        std::set<WireCell::WirePlaneId> wpids() const { return cache().wpids; }

        const std::map< int, mapfp_t< std::map<int, std::pair<double, double>> > >& all_dead_winds() const {
            // this is added in order that we may dump it in json_summary() for debugging.
            return m_dead_winds;
        }

        // FIXME: need to remove apa=0
        std::map<int, std::pair<double, double>>& get_dead_winds(const int face, const int pind, const int apa=0) const
        {
            // make one if not exist
            return m_dead_winds[apa][face][pind];

            // This is utter garbage.  #381.
        }
        using sv2d_t = Tree::ScopedView<float_t>;
        using kd2d_t = sv2d_t::nfkd_t;
        using kd_results_t = kd2d_t::results_type;

        const kd2d_t& kd2d(const int face, const int pind, const int apa = 0) const;

        const mapfp_t<double>& proj_centers() const {
            return cache().proj_centers;
        }
        const mapfp_t<double>& pitch_mags() const {
            return cache().pitch_mags;
        }

        bool is_good_point(const geo_point_t& point, const int face, const double radius = 0.6 * units::cm, const int ch_range = 1,
                           const int allowed_bad = 1) const;
        // In Facade_Grouping.h, inside the Grouping class declaration
        bool is_good_point_wc(const geo_point_t& point, const int face, const double radius = 0.6 * units::cm, 
                            const int ch_range = 1, const int allowed_bad = 1) const;
        // In the Grouping class declaration in Facade_Grouping.h
        std::vector<int> test_good_point(const geo_point_t& point, const int face, 
                                double radius = 0.6 * units::cm, int ch_range = 1) const;


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
        // In class Grouping definition
        std::pair<double,double> convert_time_ch_2Dpoint(const int timeslice, const int channel, const int face, const int plane) const;

        /// @brief Get number of points for a given plane
        /// @param plane The plane index (0=U, 1=V, 2=W)
        /// @return Number of points in the specified plane
        size_t get_num_points(const int face, const int pind, const int apa = 0) const;

        // In Facade_Grouping.h, add to public section:
        double get_ave_3d_charge(const geo_point_t& point, const double radius = 0.3 * units::cm, const int face = 0) const;
        double get_ave_charge(const geo_point_t& point, const double radius = 0.3 * units::cm, const int face = 0, const int pind = 0, const int apa = 0) const;

        /// @brief Get ranges of dead channels that overlap with given time and channel window
        /// @param min_time Minimum time
        /// @param max_time Maximum time  
        /// @param min_ch Minimum channel
        /// @param max_ch Maximum channel
        /// @param face Face number
        /// @param pind Plane index
        /// @param flag_ignore_time If true, ignore time window check
        /// @return Vector of pairs representing ranges of dead channels
        std::vector<std::pair<int, int>> get_overlap_dead_chs(const int min_time, const int max_time, 
            const int min_ch, const int max_ch, const int face, const int pind, 
            const bool flag_ignore_time=false) const;

        // In Facade_Grouping.h, inside the Grouping class public section:
        std::map<int, std::pair<int, int>> get_all_dead_chs(const int face, const int pind, int expand = 12) const;
        // Get overlapping good channel charges in a time-channel window
        std::map<std::pair<int,int>, std::pair<double,double>> get_overlap_good_ch_charge(
            int min_time, int max_time, int min_ch, int max_ch, 
            const int face, const int pind, const int apa = 0) const;

        // We override this from Mixin in order to inject propagation of the
        // utter garbage handling of dead_winds.  If someone fixes that, this
        // method may be removed.  #381.
        virtual void clear_cache() const;

      private:

        // This "cache" is utterly abused.  Someone else fix it.  #381.
        mutable std::map< int, mapfp_t< std::map<int, std::pair<double, double>> > > m_dead_winds;

       protected:
        // Receive notification when this facade is created on a node. #381.
        virtual void on_construct(node_type* node);

        virtual void fill_cache(GroupingCache& cache) const;
        
    };
    std::ostream& operator<<(std::ostream& os, const Grouping& grouping);

    // Dump the grouping to a string eg for debugging.  Level 0 dumps info about
    // just the grouping, 1 includes info about the clusters and 2 includes info
    // about the blobs.
    std::string dump(const Grouping& grouping, int level = 0);

}  // namespace WireCell::PointCloud::Facade

template <> struct fmt::formatter<WireCell::PointCloud::Facade::Grouping> : fmt::ostream_formatter {};

#endif
