/** Gap measurement for the trajectory END TRIM (doc pdvd/38).
 *
 * `TrackFitting::examine_end_ps_vec` is a TIP test: it pops points off each end
 * of a fitted path while the strict three-plane good-point test fails, and
 * BREAKS at the first point that passes.  It never examines the interior.  That
 * is safe only while the good-point test is loose enough to be wrong in the
 * helpful direction -- doc pdvd/36 sec 10.1 measured the isotropic test failing
 * 82 % of points that sit on real charge, which is what used to walk a
 * 55 cm chord back off PDVD 039252/2 cluster 109 by accident.  With the
 * anisotropic ctpc metric (CtpcAnisoMetric.h) the same trajectory passes at
 * 0.91 on charge, the tip of a detached fragment correctly passes, and nothing
 * is popped: the chord survives as a 46 cm tail.
 *
 * The retune this enables: after the pop loop stops at a passing tip, ask how
 * far inward the next SUPPORTED point is.  On real charge under the metric the
 * failing runs are 1-6 points (0.6-3.6 cm at the ~0.6 cm fit-point spacing);
 * across the cluster-109 chord the run is 90 points (54 cm).  A threshold in
 * between separates a detached island from a track.  The same rule cannot be
 * written against the isotropic test, where failing runs are everywhere.
 *
 * Kept here, header-only and templated on the iterator and the predicate, so
 * that the rule can be doctested on a synthetic path with a synthetic
 * "supported" predicate -- `examine_end_ps_vec` itself needs a Grouping, a
 * DetectorVolumes and a PCTransformSet and is not reachable from a unit test.
 */
#ifndef WIRECELLCLUS_ENDTRIMGAP
#define WIRECELLCLUS_ENDTRIMGAP

#include <cstddef>
#include <iterator>
#include <cmath>

namespace WireCell::Clus {

    /// Default bound on how far inward to look for a supported point.  Not a
    /// physics knob: it caps the cost of the scan on a pathological path where
    /// nothing passes.  At the ~0.6 cm spacing of fitted trajectory points this
    /// is ~6 m of path, i.e. longer than any PDVD cluster.
    static constexpr std::size_t end_trim_gap_max_scan = 1000;

    /// What the scan found at one end of a fitted path.
    struct TipIsland {
        /// Points in the SUPPORTED stretch at the tip (>= 1 when the caller has
        /// already stopped its pop loop on a passing point).
        std::size_t npoints{0};
        /// Path length of that supported stretch.
        double supported_len{0.0};
        /// Path length of the unsupported run immediately inward of it.
        double gap_len{0.0};
        /// Path length of the supported stretch found immediately BEYOND that
        /// run -- the body the tip would be detached from.  Measured only until
        /// it exceeds `supported_len`, which is all the caller needs to know.
        double body_len{0.0};
        /// True when a supported point was found beyond that run, i.e. there is
        /// a body at all.  False (the conservative answer) when the scan ran
        /// off the end of the path or hit max_scan.
        bool detached{false};
    };

    /**
     * Measure the supported stretch at the tip `*first` and the unsupported run
     * that follows it, walking `[first, last)` in order (use `rbegin()/rend()`
     * to walk a path inward from its far end).
     *
     * The rule this exists for -- doc pdvd/36 sec 10.1 -- is that a trajectory
     * end can be a DETACHED ISLAND: a short stretch that really does sit on
     * charge (so the tip test passes and the legacy trim stops there), reached
     * from the body of the track only by a long chord through empty space.  On
     * PDVD 039252/2 cluster 109 the island is 1 point and the chord 53.7 cm; on
     * 039349/48 cluster 53 the island is 2 points and the chord 65.8 cm.  A
     * caller decides with
     *
     *     isl.detached && isl.gap_len > L
     *                  && isl.gap_len  > isl.supported_len
     *                  && isl.body_len > isl.supported_len
     *
     * i.e. the island is the SMALLEST of the three: shorter than the gap that
     * separates it and shorter than the body it is separated from.  Both
     * comparison clauses carry physics.
     *   - gap > island protects nothing on its own but states the geometry: a
     *     stretch holding on across something longer than itself is suspect.
     *   - body > island is what protects a real track end lying beyond a
     *     dead-channel region.  The strict three-plane test has no bad-plane
     *     allowance, so a dead region reads exactly like empty space; but there
     *     the SUPPORTED stretch at the tip is the track's own end (tens of cm)
     *     and the remainder beyond the gap is smaller, so nothing is popped.
     *     Without this clause the rule chews: a path whose supported runs are
     *     all short is eaten one island at a time (measured, doc 38 sec 2).
     *
     * `is_good` is called at most `max_scan` times.
     */
    template <typename Iter, typename GoodFn>
    TipIsland scan_tip_island(Iter first, Iter last, GoodFn is_good,
                              std::size_t max_scan = end_trim_gap_max_scan)
    {
        TipIsland out;
        if (first == last) return out;
        auto dist = [](Iter a, Iter b) {
            const double dx = b->x() - a->x();
            const double dy = b->y() - a->y();
            const double dz = b->z() - a->z();
            return std::sqrt(dx*dx + dy*dy + dz*dz);
        };
        std::size_t used = 0;
        Iter prev = first;
        Iter it = first;
        out.npoints = 1;                       // the tip itself
        // (1) the supported stretch at the tip
        for (++it; it != last && used < max_scan; ++it) {
            ++used;
            const double step = dist(prev, it);
            if (!is_good(*it)) {
                out.gap_len = step;            // first step into the gap
                prev = it;
                ++it;
                break;
            }
            out.supported_len += step;
            ++out.npoints;
            prev = it;
        }
        if (out.gap_len == 0.0) return out;    // ran out while still supported
        // (2) the unsupported run, and whether a body follows it
        for (; it != last && used < max_scan; ++it) {
            ++used;
            out.gap_len += dist(prev, it);
            if (is_good(*it)) { out.detached = true; prev = it; ++it; break; }
            prev = it;
        }
        if (!out.detached) return out;         // no body found
        // (3) enough of that body to compare against the island.  Stop as soon
        // as it wins: the caller only asks whether body_len > supported_len.
        for (; it != last && used < max_scan; ++it) {
            ++used;
            if (!is_good(*it)) break;
            out.body_len += dist(prev, it);
            prev = it;
            if (out.body_len > out.supported_len) break;
        }
        return out;
    }

}  // namespace WireCell::Clus

#endif
