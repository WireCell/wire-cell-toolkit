/** Cross-blob minimum-separation thinning for Steiner terminals.

    doc pdvd/37 R1.  This is the pure-geometry core of
    Steiner::Grapher::thin_terminals_by_separation, split out of the
    (src-private) SteinerGrapher.h so it can be exercised directly by a
    doctest -- the same reason CtpcAnisoMetric.h exists.

    Background.  find_peak_point_indices runs once per blob, and its
    map_index_charge holds only that blob's points, so the local-maximum test
    skips every out-of-blob neighbour and can never remove a blob's LAST
    candidate.  Steiner terminal density is therefore floored at the
    candidate-bearing blob density -- 1.02 terminals per such blob measured in
    doc pdvd/31 round 6 -- which makes the spacing a function of the time-slice
    pitch and the track's drift alignment rather than of anything physical.
    Doc pdvd/37 sec.3.2 confirms that law holds on PDVD, SBND and uBooNE alike.

    The remedy here does not touch selection.  It filters the selected set: walk
    the terminals in decreasing charge order and admit one only when no
    already-admitted terminal lies strictly within min_separation of it.
 */

#ifndef WIRECELLCLUS_STEINERTHINNING
#define WIRECELLCLUS_STEINERTHINNING

#include "WireCellUtil/Point.h"

#include <cstddef>
#include <utility>
#include <vector>

namespace WireCell::Clus::Steiner {

    /// Greedy minimum-separation thinning.
    ///
    /// @param ordered  (position, id) pairs ALREADY in admission-priority
    ///                 order -- decreasing charge at the production call site.
    ///                 The caller owns the ordering, so this function is a pure
    ///                 function of its arguments and carries no tie-break of
    ///                 its own.
    /// @param min_separation  the spacing to enforce, in the WCT system of
    ///                 units.  <= 0 disables thinning and every id is returned.
    ///
    /// @return the ids kept, in the order they were admitted.
    ///
    /// The exclusion test is STRICT: a terminal exactly @p min_separation from
    /// an already-kept one is KEPT, so the guarantee the survivors satisfy is
    /// "no two closer than min_separation", not "no two at min_separation".
    /// The doctest pins this; it is not free to be "simplified".
    ///
    /// Cost is O(N) in the number of terminals via a uniform grid of cell edge
    /// min_separation, so only 27 cells are ever consulted per candidate.  The
    /// result is independent of that grid: cell membership is used to shrink
    /// the candidate list, never to decide admission.
    std::vector<size_t> thin_by_min_separation(
        const std::vector<std::pair<Point, size_t>>& ordered,
        double min_separation);

}  // namespace WireCell::Clus::Steiner

#endif
