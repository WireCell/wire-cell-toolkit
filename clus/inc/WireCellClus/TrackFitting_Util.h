#ifndef WIRECELLCLUS_TRACKFITTING_UTIL_H
#define WIRECELLCLUS_TRACKFITTING_UTIL_H

#include "WireCellUtil/Point.h"

#include <vector>

namespace WireCell::Clus::TrackFittingUtil {

    /** doc pdvd/111 -- move a seed point TRANSVERSELY onto the local charge ridge.
     *
     * Gaussian-kernel mean shift restricted to the plane perpendicular to `dir`:
     * each step adds sum_i w_i k_i perp_i / sum_i w_i k_i, where perp_i is the
     * component of (pts[i] - p) perpendicular to dir, k_i = exp(-|perp_i|^2 / 2 sigma^2),
     * and only points with |along_i| <= half_slab and |perp_i| <= 3 sigma count.
     * The along-dir coordinate of the point never changes.  Stops after max_iter
     * steps or when a step is below 1e-3 sigma.
     *
     * Returns p0 unchanged when dir is zero, sigma <= 0, no point carries weight,
     * or the total move would exceed max_move (a far mode is not this track's ridge).
     * The sums run in the given order; callers wanting run-to-run bit identity
     * pass the points in a stable order.
     */
    WireCell::Point recenter_point_transverse(const WireCell::Point& p0, const WireCell::Vector& dir,
                                              const std::vector<WireCell::Point>& pts,
                                              const std::vector<double>& weights,
                                              double sigma, double half_slab, int max_iter, double max_move);

    /** Calculate ranges for track fitting using simplified coupling coefficients.
     *
     * This function computes available ranges for each wire plane (U, V, W) based on
     * geometric coupling between planes and minimum distance constraints.
     *
     * @param angle_u Wire angle for U plane
     * @param angle_v Wire angle for V plane  
     * @param angle_w Wire angle for W plane
     * @param rem_dis_cut_u Remaining distance cut for U plane
     * @param rem_dis_cut_v Remaining distance cut for V plane
     * @param rem_dis_cut_w Remaining distance cut for W plane
     * @param min_u_dis Minimum distance for U plane
     * @param min_v_dis Minimum distance for V plane
     * @param min_w_dis Minimum distance for W plane
     * @param pitch_u Wire pitch for U plane
     * @param pitch_v Wire pitch for V plane
     * @param pitch_w Wire pitch for W plane
     * @param range_u [out] Calculated range for U plane
     * @param range_v [out] Calculated range for V plane
     * @param range_w [out] Calculated range for W plane
     */
    void calculate_ranges_simplified(double angle_u, double angle_v, double angle_w,
                                   double rem_dis_cut_u, double rem_dis_cut_v, double rem_dis_cut_w,
                                   double min_u_dis, double min_v_dis, double min_w_dis,
                                   double pitch_u, double pitch_v, double pitch_w,
                                   float& range_u, float& range_v, float& range_w);

} // namespace WireCell::Clus::TrackFittingUtil

#endif