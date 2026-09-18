/** Blank-plane admission policy for Steiner terminal candidates.

    doc pdvd/114.  This is the pure core of the `terminal_blank_plane_mode`
    knob of CreateSteinerGraph / Steiner::Grapher, split out of the
    (src-private) SteinerGrapher.h so a doctest can exercise it without a
    cluster -- the same reason SteinerThinning.h and CtpcAnisoMetric.h exist.

    Background.  Phase 1 of create_steiner_tree admits a point as a terminal
    candidate through Cluster::calc_charge_wcp with disable_dead_mix_cell =
    false: a plane at charge EXACTLY 0 passes the per-plane test and is left
    out of the charge RMS, so a point that sees charge on two planes and sits
    on an empty cell in the third is graded on the two bright planes alone.
    The retiler writes its painted and dead cells as charge 0, so such points
    exist in a halo beside every track, and doc 112 showed the Steiner seed
    leaves the image through them.  Doc 114 measured that in about three of
    four such off-image terminals the SAME blob also holds a three-plane,
    on-image candidate that lost the charge ordering -- and that the per-blob
    peak finder only suppresses ADJACENT candidates, so re-scoring alone does
    not reach them.

    The policy is a filter on a blob's candidate set, applied before the peak
    search, and it never touches a blob whose candidates ALL have a zero
    plane: inside a dead or inefficient region every point has one, so the
    graph still bridges the gap exactly as before.

      wcp             today: no filter (the C++ default; byte-identical)
      prefer3         drop the candidates with a zero plane when the blob holds
                      at least one candidate with charge on all three planes
      nearby          drop a candidate with a zero plane when a three-plane
                      candidate lies within `radius` of it (the caller supplies
                      the predicate; it is a kd radius query on the cluster)
      prefer3+nearby  both

    A dead plane is a zero plane too.  Neither rule needs the dead-channel
    registry: a blob wholly inside a dead region has no three-plane candidate
    and is left alone; at the region's edge a two-plane candidate may yield to
    a three-plane neighbour, which moves a terminal by at most a blob width.
 */

#ifndef WIRECELLCLUS_STEINERBLANKPLANE
#define WIRECELLCLUS_STEINERBLANKPLANE

#include <cstddef>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace WireCell::Clus::Steiner {

    enum class BlankPlaneMode { wcp = 0, prefer3, nearby, prefer3_nearby };

    /// Parse the config string.  Unknown strings return false and leave
    /// `mode` untouched, so a caller can refuse a typo instead of silently
    /// running the legacy path (the doc pdvd/113 retile_mode precedent).
    inline bool parse_blank_plane_mode(const std::string& s, BlankPlaneMode& mode)
    {
        if (s == "wcp") { mode = BlankPlaneMode::wcp; return true; }
        if (s == "prefer3") { mode = BlankPlaneMode::prefer3; return true; }
        if (s == "nearby") { mode = BlankPlaneMode::nearby; return true; }
        if (s == "prefer3+nearby") { mode = BlankPlaneMode::prefer3_nearby; return true; }
        return false;
    }

    /// Apply the policy to one blob's candidate set.
    ///
    /// @param candidates  (charge, point index) pairs of the blob, any order
    ///                    (the caller's std::set with std::greater at the
    ///                    production site).
    /// @param mode        the policy.
    /// @param nzero       nzero(idx) -> number of planes at charge exactly 0
    ///                    for point idx.
    /// @param near3       near3(idx) -> true when a three-plane candidate lies
    ///                    within the configured radius of point idx (only
    ///                    called in the nearby modes, and only for candidates
    ///                    that have a zero plane).
    /// @return the surviving candidates, same type, same order.  Under `wcp`
    ///         the input is returned unchanged.
    template <class Candidates, class NZero, class Near3>
    Candidates apply_blank_plane_policy(const Candidates& candidates, BlankPlaneMode mode, NZero nzero, Near3 near3)
    {
        if (mode == BlankPlaneMode::wcp || candidates.empty()) {
            return candidates;
        }
        const bool use_prefer = (mode == BlankPlaneMode::prefer3 || mode == BlankPlaneMode::prefer3_nearby);
        const bool use_near = (mode == BlankPlaneMode::nearby || mode == BlankPlaneMode::prefer3_nearby);
        bool blob_has_3 = false;
        if (use_prefer) {
            for (const auto& c : candidates) {
                if (nzero(c.second) == 0) { blob_has_3 = true; break; }
            }
        }
        Candidates out;
        for (const auto& c : candidates) {
            const int nz = nzero(c.second);
            bool keep = true;
            if (nz > 0) {
                if (use_prefer && blob_has_3) keep = false;
                if (keep && use_near && near3(c.second)) keep = false;
            }
            if (keep) out.insert(out.end(), c);
        }
        return out;
    }

}  // namespace WireCell::Clus::Steiner

#endif
