/** 2-D charge support along a straight Steiner-graph edge.

    doc pdvd/111 round 2.  The STM rough path is a Dijkstra walk on
    "steiner_graph", whose edge weight prices charge at the two endpoints only
    (SteinerGrapher.cxx create_enhanced_steiner_graph).  Whether the straight
    segment BETWEEN two Steiner points sits on measured 2-D charge is therefore
    never looked at.  This header samples that segment and counts, per sample,
    how many wire planes see charge (or a dead channel) there, so the log-only
    WCT_STEINER_GRAPH_DUMP can record it per edge.

    The sampling is the one gap_edge_bad_fraction uses (PRSegmentFunctions.h,
    doc sbnd_xin/pr/51 round 5): round(length / step) intervals, both endpoints
    included.  The classifier is injected, so this core needs no cluster, no
    grouping and no geometry, and the doctest can pin it directly.
 */

#ifndef WIRECELLCLUS_STEINEREDGESUPPORT
#define WIRECELLCLUS_STEINEREDGESUPPORT

#include "WireCellUtil/Point.h"

#include <algorithm>
#include <cmath>

namespace WireCell::Clus::Steiner {

    /// Per-edge sample counts.  A plane is "ok" at a sample when it is live
    /// (charge within the radius) or, failing that, dead.
    struct EdgeSupport {
        int n{0};      ///< samples taken, both endpoints included
        int ok3{0};    ///< samples ok in all three planes
        int ok2{0};    ///< samples ok in at least two planes
        int ok1{0};    ///< samples ok in at least one plane
        int live1{0};  ///< samples live in at least one plane
    };

    /// Sample the segment a -> b every ~step and classify each sample.
    ///
    /// @param classify  bool(const Point&, int planes[6]).  It fills planes
    ///                  in the Grouping::test_good_point layout {live_u,
    ///                  live_v, live_w, dead_u, dead_v, dead_w} (counts, > 0
    ///                  = yes) and returns false when the sample lies outside
    ///                  every TPC.  Such a sample is counted in n and in no
    ///                  other field.
    /// @param step      sampling step in the WCT system of units; <= 0
    ///                  returns an all-zero record.
    template <class Classify>
    EdgeSupport edge_support(const Point& a, const Point& b, double step, Classify&& classify)
    {
        EdgeSupport s;
        if (!(step > 0)) return s;
        const double len = (b - a).magnitude();
        const int nsteps = std::max(1, static_cast<int>(std::round(len / step)));
        for (int k = 0; k <= nsteps; ++k) {
            const double f = static_cast<double>(k) / nsteps;
            const Point p(a.x() + (b.x() - a.x()) * f, a.y() + (b.y() - a.y()) * f,
                          a.z() + (b.z() - a.z()) * f);
            int planes[6] = {0, 0, 0, 0, 0, 0};
            const bool inside = classify(p, planes);
            int nok = 0;
            bool live = false;
            if (inside) {
                for (int pl = 0; pl < 3; ++pl) {
                    if (planes[pl] > 0) {
                        ++nok;
                        live = true;
                    }
                    else if (planes[pl + 3] > 0) {
                        ++nok;
                    }
                }
            }
            ++s.n;
            if (nok >= 3) ++s.ok3;
            if (nok >= 2) ++s.ok2;
            if (nok >= 1) ++s.ok1;
            if (live) ++s.live1;
        }
        return s;
    }

}  // namespace WireCell::Clus::Steiner

#endif
