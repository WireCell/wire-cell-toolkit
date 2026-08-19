#ifndef WIRECELL_CLUS_PR_SHOWER_FUNCTIONS
#define WIRECELL_CLUS_PR_SHOWER_FUNCTIONS

#include "WireCellClus/PRShower.h"
#include "WireCellUtil/Units.h"


namespace WireCell::Clus::PR {

    /** Modify shower assuming shower kinematics.
     *
     * This free function is is equivalent to the method of WCP's
     * WCShower::calculate_kinematics().
     */

     std::pair<double, WireCell::Point> shower_get_closest_point(Shower& shower, const WireCell::Point& point, const std::string& cloud_name = "fit");

     double shower_get_closest_dis(Shower& shower, SegmentPtr seg, const std::string& cloud_name = "fit");

     double shower_get_dis(Shower& shower, SegmentPtr seg, const std::string& cloud_name = "fit");

     WireCell::Vector shower_cal_dir_3vector(Shower& shower, const WireCell::Point& p, double dis_cut = 15*units::cm);

     /** doc sbnd_xin/docs/pr/92 -- true iff the shower's start segment is the
      * collinear continuation of a straight long track that is NOT part of
      * the shower itself (e.g. the walled-off tail of an overclustered
      * cosmic: SBND 350935 shower 11001, 321371 shower 18004).
      *
      * Deliberately NOT segment_is_straight_long_track_or_continuation:
      * (a) that predicate's first arm (the segment itself is straight-long)
      * false-fires on legitimate EM trunks -- a 1.8 GeV NCpi0 gamma's 11 cm
      * trunk passes the direct > 0.93*length test (ncpi0 506114 shower
      * 19016); (b) its continuation arm accepts siblings that are members
      * of this same shower, which every dense EM tree has.  Here the
      * collinear sibling must be OUTSIDE the shower (membership tested
      * against the shower's own segment set via fill_sets, not shower-id
      * fields) and itself straight-long (or the combined seg+sib chain
      * qualifies).  Same-cluster siblings only, matching the round-9
      * continuation rule.
      */
     bool shower_start_is_track_continuation(Graph& graph, Shower& shower,
                                             double max_kink_deg = 25.0);

}

#endif
