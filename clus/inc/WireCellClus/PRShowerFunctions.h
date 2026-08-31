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

}

#endif
