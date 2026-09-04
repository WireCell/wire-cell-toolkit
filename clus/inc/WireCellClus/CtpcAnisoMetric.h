/** Anisotropic (lattice-normalised) radius queries on the ctpc k-d trees.
 *
 * doc pdvd/34 and pdvd/36.  The ctpc ("charge-time point cloud") that
 * Grouping::kd2d() indexes is an EXACT 2-D lattice: x is the drift
 * coordinate quantised by the time slice (nticks_live_slice * tick *
 * drift_speed, 2.96 mm on PDVD, 3.13 mm on SBND) and y is the projected
 * wire coordinate quantised by the pitch (7.65 mm PDVD U/V, 5.10 mm PDVD W,
 * 3.00 mm SBND).  Every ctpc distance in the toolkit was an isotropic
 * Euclidean sqrt of the quadrature sum, which is harmless where the two cell
 * edges are nearly equal (SBND 0.96) and not where they are not (PDVD U/V
 * 2.58): a 0.2 cm radius spans 0.68 of a drift step but only 0.26 of a U/V
 * pitch, and an interior trajectory point passes the strict three-plane
 * good-point test 17.5 % of the time on PDVD against 62.4 % on SBND.
 *
 * The metric used here scales ONLY the pitch axis:
 *
 *     s   = min(1, drift_step / pitch)            (per apa, face, plane)
 *     d^2 = dx^2 + (s * dy)^2                     (< r^2 accepts)
 *
 * so the radius keeps its exact drift-axis meaning and the pitch reach, in
 * units of a pitch, becomes r / drift_step on every plane -- the same number
 * of lattice cells as on SBND.  The clamp at 1 makes a detector whose pitch is
 * finer than its drift step (SBND, 1.04) the identity: s == 1 is the legacy
 * isotropic metric exactly, not approximately.
 *
 * The tree stays isotropic (NFKDVec exposes one scalar squared radius and
 * the stored (x, y) are not touched -- see doc 34 sec 7 for why scaling the
 * stored y at build time was rejected).  The ellipse is answered exactly in
 * two levels:
 *
 *   level 1  an isotropic k-d query with the CIRCUMSCRIBING radius r / s
 *            (the ellipse's long semi-axis, along the pitch);
 *   level 2  the per-candidate ellipse test above.
 *
 * Level 1 must query with the LARGER radius, not r: the ellipse sticks out
 * past the circle of radius r along the pitch, and querying r would silently
 * drop points that should match.  Since every point of the ellipse lies
 * inside the circumscribing circle, no point can escape, and the filter is
 * exact.  For the existence-only form two brackets keep the common cases
 * allocation-free: a hit inside the INSCRIBED circle (radius r) is inside the
 * ellipse, and an empty circumscribing circle means an empty ellipse.  Both
 * are exact, so no correctness rests on them.
 *
 * TRAP.  Only an all-candidates query (NFKDVec::Tree::radius) may be
 * filtered.  A nearest-only query (knn(1), knn1) cannot: the Euclidean-nearest
 * point can lie outside the ellipse while a Euclidean-farther point lies
 * inside it.  exists_within() returns a boolean and is used here only as the
 * two exact brackets.
 *
 * The returned pairs keep the ISOTROPIC squared distance the tree computed
 * (physical length squared, what callers sqrt into a length); only the
 * membership decision is anisotropic.
 *
 * These are function templates over the tree type so that they can be
 * exercised on a bare NFKDVec::Tree without any detector geometry
 * (clus/test/doctest_ctpc_aniso_metric.cxx).  Their only production consumer
 * is Facade::Grouping::get_closest_points / has_closest_point.
 */

#ifndef WIRECELL_CLUS_CTPCANISOMETRIC
#define WIRECELL_CLUS_CTPCANISOMETRIC

#include <algorithm>
#include <array>
#include <cmath>

namespace WireCell::Clus::Facade {

    /// The pitch-axis scale s = min(1, drift_step / pitch).  Returns 1 (the
    /// identity, i.e. the legacy isotropic metric) whenever either input is
    /// non-finite or non-positive, so a grouping whose DetectorVolumes
    /// metadata carries no nticks_live_slice can never produce a nonsense
    /// metric -- the caller is expected to log that case loudly.
    inline double ctpc_yscale(double drift_step, double pitch)
    {
        if (!std::isfinite(drift_step) || !std::isfinite(pitch)) return 1.0;
        if (drift_step <= 0.0 || pitch <= 0.0) return 1.0;
        return std::min(1.0, drift_step / pitch);
    }

    /// The level-2 test: is the tree point (px, py) inside the ellipse of
    /// drift semi-axis r and pitch semi-axis r/s centred on the query?  Strict
    /// '<', matching nanoflann's radius() / exists_within() comparison.
    inline bool ctpc_in_ellipse(double px, double py, double qx, double qy, double radius, double yscale)
    {
        const double dx = px - qx;
        const double dy = (py - qy) * yscale;
        return dx * dx + dy * dy < radius * radius;
    }

    /// All tree points inside the ellipse.  Equivalent to kd.radius(r*r, q)
    /// when yscale >= 1 (and dispatches to it verbatim); otherwise the two
    /// level query above.  Result pairs are (flat index, isotropic squared
    /// distance), in the order the circumscribing query returned them.
    template <typename KD, typename Query>
    typename KD::results_type ctpc_radius_aniso(const KD& kd, const Query& query, double radius, double yscale)
    {
        if (yscale >= 1.0) {
            return kd.template radius<Query>(radius * radius, query);
        }
        const double r_out = radius / yscale;
        auto cands = kd.template radius<Query>(r_out * r_out, query);
        typename KD::results_type ret;
        ret.reserve(cands.size());
        for (const auto& c : cands) {
            const double px = kd.kdtree_get_pt(c.first, 0);
            const double py = kd.kdtree_get_pt(c.first, 1);
            if (ctpc_in_ellipse(px, py, query[0], query[1], radius, yscale)) {
                ret.push_back(c);
            }
        }
        return ret;
    }

    /// True iff any tree point lies inside the ellipse.  Boolean-identical
    /// to !ctpc_radius_aniso(...).empty(); the two exact brackets make the
    /// common cases a single early-terminating k-d descent with no
    /// allocation, and only the in-between case enumerates candidates.
    template <typename KD, typename Query>
    bool ctpc_exists_within_aniso(const KD& kd, const Query& query, double radius, double yscale)
    {
        if (yscale >= 1.0) {
            return kd.template exists_within<Query>(radius * radius, query);
        }
        // Inscribed circle (radius r) is contained in the ellipse => accept.
        if (kd.template exists_within<Query>(radius * radius, query)) return true;
        // Ellipse is contained in the circumscribing circle (radius r/s) => an
        // empty circle means an empty ellipse.
        const double r_out = radius / yscale;
        if (!kd.template exists_within<Query>(r_out * r_out, query)) return false;
        // In between: enumerate the circumscribing circle and test each.
        const auto cands = kd.template radius<Query>(r_out * r_out, query);
        for (const auto& c : cands) {
            const double px = kd.kdtree_get_pt(c.first, 0);
            const double py = kd.kdtree_get_pt(c.first, 1);
            if (ctpc_in_ellipse(px, py, query[0], query[1], radius, yscale)) return true;
        }
        return false;
    }

}  // namespace WireCell::Clus::Facade

#endif
