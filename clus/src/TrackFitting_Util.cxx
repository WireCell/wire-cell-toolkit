#include "WireCellClus/TrackFitting_Util.h"
#include <cmath>
#include <iostream>

using namespace WireCell::Clus::TrackFittingUtil;

void WireCell::Clus::TrackFittingUtil::calculate_ranges_simplified(
    double angle_u, double angle_v, double angle_w,
    double rem_dis_sq_cut_u, double rem_dis_sq_cut_v, double rem_dis_sq_cut_w,
    double min_u_dis, double min_v_dis, double min_w_dis,
    double pitch_u, double pitch_v, double pitch_w,
    float& range_sq_u, float& range_sq_v, float& range_sq_w) {


    // Geometric coupling coefficients
    double coupling_uv = fabs(cos(angle_u - angle_v));
    double coupling_uw = fabs(cos(angle_u - angle_w));
    double coupling_vw = fabs(cos(angle_v - angle_w));

    // std::cout << "Angles: " << rem_dis_cut_u << " " << rem_dis_cut_v << " " << rem_dis_cut_w << " " << coupling_uv << " " << coupling_uw << " " << coupling_vw << " | " << angle_u << " " << angle_v << " " << angle_w << " " << min_u_dis << " " << min_v_dis << " " << min_w_dis << std::endl;


    // Cost from other planes (weighted by coupling)
    double cost_u_from_v = coupling_uv * pow(min_v_dis * pitch_v, 2)/coupling_vw;
    double cost_u_from_w = coupling_uw * pow(min_w_dis * pitch_w, 2)/coupling_vw;
    
    double cost_v_from_u = coupling_uv * pow(min_u_dis * pitch_u, 2)/coupling_uw;
    double cost_v_from_w = coupling_vw * pow(min_w_dis * pitch_w, 2)/coupling_uw;
    
    double cost_w_from_u = coupling_uw * pow(min_u_dis * pitch_u, 2)/coupling_uv;
    double cost_w_from_v = coupling_vw * pow(min_v_dis * pitch_v, 2)/coupling_uv;
    
    // Calculate available ranges
    double available_u = rem_dis_sq_cut_u*(coupling_uv + coupling_uw + coupling_vw)  - cost_u_from_v - cost_u_from_w;
    double available_v = rem_dis_sq_cut_v*(coupling_uv + coupling_uw + coupling_vw)  - cost_v_from_u - cost_v_from_w;
    double available_w = rem_dis_sq_cut_w*(coupling_uv + coupling_uw + coupling_vw)  - cost_w_from_u - cost_w_from_v;

    range_sq_u = (available_u > 0) ? available_u : 0;
    range_sq_v = (available_v > 0) ? available_v : 0;
    range_sq_w = (available_w > 0) ? available_w : 0;
}

WireCell::Point WireCell::Clus::TrackFittingUtil::recenter_point_transverse(
    const WireCell::Point& p0, const WireCell::Vector& dir,
    const std::vector<WireCell::Point>& pts, const std::vector<double>& weights,
    double sigma, double half_slab, int max_iter, double max_move)
{
    const double dmag = dir.magnitude();
    if (!(dmag > 0) || !(sigma > 0) || pts.size() != weights.size()) return p0;
    const WireCell::Vector d = dir / dmag;
    const double rmax2 = 9.0 * sigma * sigma;
    const double inv2s2 = 1.0 / (2.0 * sigma * sigma);
    WireCell::Point p = p0;
    for (int it = 0; it < max_iter; ++it) {
        double sw = 0, sx = 0, sy = 0, sz = 0;
        for (size_t i = 0; i != pts.size(); ++i) {
            if (!(weights[i] > 0)) continue;
            const WireCell::Vector v = pts[i] - p;
            const double along = v.dot(d);
            if (std::abs(along) > half_slab) continue;
            const WireCell::Vector perp = v - d * along;
            const double r2 = perp.dot(perp);
            if (r2 > rmax2) continue;
            const double w = weights[i] * std::exp(-r2 * inv2s2);
            sw += w; sx += w * perp.x(); sy += w * perp.y(); sz += w * perp.z();
        }
        if (!(sw > 0)) break;
        const WireCell::Vector step(sx / sw, sy / sw, sz / sw);
        p = p + step;
        if (step.magnitude() < 1e-3 * sigma) break;
    }
    if ((p - p0).magnitude() > max_move) return p0;
    return p;
}
