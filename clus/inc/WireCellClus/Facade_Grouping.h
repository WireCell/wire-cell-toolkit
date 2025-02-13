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

#include "WireCellClus/Facade_Util.h"


// using namespace WireCell;  NO!  do not open up namespaces in header files!

namespace WireCell::PointCloud::Facade {
    class Cluster;

    // Give a node "Grouping" semantics.  A grouping node's children are cluster
    // nodes that are related in some way.
    class Grouping : public NaryTree::FacadeParent<Cluster, points_t>, public Mixin<Grouping> {

        TPCParams m_tp{};  // use default value by default.
        /// TODO: replace TPCParams with this in the future?
        IAnodePlane::pointer m_anode{nullptr};

       public:

        Grouping() : Mixin<Grouping>(*this, "grouping_scalar") {}

        // MUST call this sometimes after construction if non-default value needed.
        // FIXME: TPCParams should be moved out of the facade!
        void set_params(const TPCParams& tp) { m_tp = tp; }
        void set_params(const WireCell::Configuration& cfg);
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
        size_t get_num_points(const int face, const int pind) const;

        // In Facade_Grouping.h, add to public section:
        double get_ave_3d_charge(const geo_point_t& point, const double radius = 0.3 * units::cm, const int face = 0) const;
        double get_ave_charge(const geo_point_t& point, const double radius = 0.3 * units::cm, const int face = 0, const int pind = 0) const;

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
            const int face, const int pind) const;

       private:
        void fill_proj_centers_pitch_mags() const;
        mutable mapfp_t<double> m_proj_centers;
        mutable mapfp_t<double> m_pitch_mags;
        mutable mapfp_t< std::map<int, std::pair<double, double>> > m_dead_winds;

       protected:
        // Receive notification when this facade is created on a node.
        virtual void on_construct(node_type* node);
    };
    std::ostream& operator<<(std::ostream& os, const Grouping& grouping);

    // Dump the grouping to a string eg for debugging.  Level 0 dumps info about
    // just the grouping, 1 includes info about the clusters and 2 includes info
    // about the blobs.
    std::string dump(const Grouping& grouping, int level = 0);

}  // namespace WireCell::PointCloud::Facade

template <> struct fmt::formatter<WireCell::PointCloud::Facade::Grouping> : fmt::ostream_formatter {};

#endif
