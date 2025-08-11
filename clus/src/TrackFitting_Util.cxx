#include "WireCellClus/TrackFitting_Util.h"
#include <cmath>

using namespace WireCell::Clus::TrackFittingUtil;

void WireCell::Clus::TrackFittingUtil::calculate_ranges_simplified(
    double angle_u, double angle_v, double angle_w,
    double rem_dis_cut_u, double rem_dis_cut_v, double rem_dis_cut_w,
    double min_u_dis, double min_v_dis, double min_w_dis,
    double pitch_u, double pitch_v, double pitch_w,
    float& range_u, float& range_v, float& range_w) {
    
    // Geometric coupling coefficients
    double coupling_uv = fabs(cos(angle_u - angle_v));
    double coupling_uw = fabs(cos(angle_u - angle_w));
    double coupling_vw = fabs(cos(angle_v - angle_w));
    
    // Cost from other planes (weighted by coupling)
    double cost_u_from_v = coupling_uv * pow(min_v_dis * pitch_v, 2);
    double cost_u_from_w = coupling_uw * pow(min_w_dis * pitch_w, 2);
    
    double cost_v_from_u = coupling_uv * pow(min_u_dis * pitch_u, 2);
    double cost_v_from_w = coupling_vw * pow(min_w_dis * pitch_w, 2);
    
    double cost_w_from_u = coupling_uw * pow(min_u_dis * pitch_u, 2);
    double cost_w_from_v = coupling_vw * pow(min_v_dis * pitch_v, 2);
    
    // Calculate available ranges
    double available_u = rem_dis_cut_u * rem_dis_cut_u - cost_u_from_v - cost_u_from_w;
    double available_v = rem_dis_cut_v * rem_dis_cut_v - cost_v_from_u - cost_v_from_w;
    double available_w = rem_dis_cut_w * rem_dis_cut_w - cost_w_from_u - cost_w_from_v;
    
    range_u = (available_u > 0) ? sqrt(available_u) : 0;
    range_v = (available_v > 0) ? sqrt(available_v) : 0;
    range_w = (available_w > 0) ? sqrt(available_w) : 0;
}