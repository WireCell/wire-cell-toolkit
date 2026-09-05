/** BadBlobRuns -- the pure graph core of ImproveCluster_1::remove_bad_blobs.

    doc pdvd/40 round 3.  The retiler (ImproveCluster_2) fabricates blobs along
    a whole-cluster shortest path (hack_activity_improved) and then filters the
    result with remove_bad_blobs.  The historical filter -- a faithful port of
    the prototype's Improve_PR3DCluster_2 (ImprovePR3DCluster.cxx:515-620) --
    builds a graph over the NEW blobs with adjacent-slice overlap edges and,
    ONLY when that graph has more than one connected component, judges each
    component by whether its FIRST blob overlaps an ORIGINAL blob within one
    slice.  Two holes let a 122 cm fabricated column through on PDVD 039252/2
    cluster 119: a single component is never examined at all, and one
    representative blob speaks for a whole component.

    This header holds the decision logic as a pure function over indices so it
    can be unit-tested without a Grouping.  The caller (remove_bad_blobs)
    supplies the adjacency, the per-blob support flags and the blob centers.

    Semantics, in order:
      1. connected components over `edges`;
      2. the historical vote, UNCHANGED: with ncomp > 1, a component is good
         iff ANY of its blobs is supported (the historical loop tests a
         component's blobs in index order and stops at the first supported
         one -- `if (good.count(comp)) continue;` skips only components
         already found good); failing components are removed whole.
         NOTE a first reading of that loop as a "first blob only" vote was
         wrong and was caught by the 120-event hash gate (doc pdvd/40 round 3);
      3. run bound: inside good components, the connected pieces of the
         UNSUPPORTED blobs are "runs"; a run whose bounding-box diagonal of
         blob centers exceeds `max_run` is removed whole.  This is what closes
         both holes: a lone component is now examined (a), and a fabricated
         column attached to a real track no longer inherits its verdict (b).
         max_run <= 0 disables step 3.
    `legacy_component_vote` is the historical filter as a pure function; the
    unit test asserts analyze(...).removed_by_vote equals it.
 */
#ifndef WIRECELLCLUS_BADBLOBRUNS_H
#define WIRECELLCLUS_BADBLOBRUNS_H

#include "WireCellUtil/Point.h"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace WireCell::Clus::BadBlobRuns {

    struct Run {
        std::vector<int> blobs;      // indices, ascending
        int component{-1};
        double span{0};              // bbox diagonal of the blob centers
        Point center{0, 0, 0};       // bbox center
        int nslices{0};              // distinct slice keys among the blobs
    };

    struct Result {
        int ncomp{0};
        std::vector<int> component;           // per blob
        std::vector<bool> good_component;     // per component id
        std::vector<Run> runs;                // every unsupported run in a good component, longest first
        std::vector<int> removed_by_vote;     // indices, ascending
        std::vector<int> removed_by_run;      // indices, ascending
    };

    // Union-find labels, 0..ncomp-1, in order of first appearance by index.
    inline int label_components(int n, const std::vector<std::pair<int, int>>& edges,
                                std::vector<int>& comp)
    {
        std::vector<int> parent(n);
        for (int i = 0; i < n; ++i) parent[i] = i;
        auto find = [&](int a) {
            while (parent[a] != a) { parent[a] = parent[parent[a]]; a = parent[a]; }
            return a;
        };
        for (const auto& [a, b] : edges) {
            int ra = find(a), rb = find(b);
            if (ra != rb) parent[std::max(ra, rb)] = std::min(ra, rb);
        }
        comp.assign(n, -1);
        std::vector<int> relabel(n, -1);
        int ncomp = 0;
        for (int i = 0; i < n; ++i) {
            int r = find(i);
            if (relabel[r] < 0) relabel[r] = ncomp++;
            comp[i] = relabel[r];
        }
        return ncomp;
    }

    /// The historical filter: with >1 component, keep a component iff ANY of
    /// its blobs is supported.  Returns the indices to remove, ascending.
    inline std::vector<int> legacy_component_vote(int n, const std::vector<std::pair<int, int>>& edges,
                                                  const std::vector<bool>& supported)
    {
        std::vector<int> comp;
        const int ncomp = label_components(n, edges, comp);
        std::vector<int> out;
        if (ncomp <= 1) return out;
        std::vector<char> good(ncomp, 0);
        for (int i = 0; i < n; ++i)
            if (supported[i]) good[comp[i]] = 1;
        for (int i = 0; i < n; ++i)
            if (!good[comp[i]]) out.push_back(i);
        return out;
    }

    /// The round-3 filter.  `slice` may be empty (nslices then reports 0).
    inline Result analyze(int n, const std::vector<std::pair<int, int>>& edges,
                          const std::vector<bool>& supported, const std::vector<Point>& centers,
                          double max_run, const std::vector<int>& slice = {})
    {
        Result r;
        r.ncomp = label_components(n, edges, r.component);
        // The historical vote: a component is good iff any blob is supported;
        // a lone component is never voted on.
        r.good_component.assign(r.ncomp, true);
        if (r.ncomp > 1) {
            std::fill(r.good_component.begin(), r.good_component.end(), false);
            for (int i = 0; i < n; ++i)
                if (supported[i]) r.good_component[r.component[i]] = true;
        }
        for (int i = 0; i < n; ++i)
            if (!r.good_component[r.component[i]]) r.removed_by_vote.push_back(i);

        // Runs: components of the unsupported blobs inside good components.
        std::vector<int> keep(n, 0);
        for (int i = 0; i < n; ++i) keep[i] = (!supported[i] && r.good_component[r.component[i]]) ? 1 : 0;
        std::vector<std::pair<int, int>> uedges;
        for (const auto& e : edges)
            if (keep[e.first] && keep[e.second]) uedges.push_back(e);
        std::vector<int> ucomp;
        const int nu = label_components(n, uedges, ucomp);
        std::vector<Run> runs(nu);
        for (int i = 0; i < n; ++i)
            if (keep[i]) runs[ucomp[i]].blobs.push_back(i);
        for (auto& run : runs) {
            if (run.blobs.empty()) continue;
            run.component = r.component[run.blobs.front()];
            Point lo = centers[run.blobs.front()], hi = lo;
            std::vector<int> sl;
            for (int i : run.blobs) {
                const auto& c = centers[i];
                lo = Point(std::min(lo.x(), c.x()), std::min(lo.y(), c.y()), std::min(lo.z(), c.z()));
                hi = Point(std::max(hi.x(), c.x()), std::max(hi.y(), c.y()), std::max(hi.z(), c.z()));
                if (!slice.empty()) sl.push_back(slice[i]);
            }
            run.span = (hi - lo).magnitude();
            run.center = (hi + lo) * 0.5;
            std::sort(sl.begin(), sl.end());
            run.nslices = int(std::unique(sl.begin(), sl.end()) - sl.begin());
            r.runs.push_back(std::move(run));
        }
        std::stable_sort(r.runs.begin(), r.runs.end(),
                         [](const Run& a, const Run& b) { return a.span > b.span; });
        if (max_run > 0) {
            for (const auto& run : r.runs)
                if (run.span > max_run)
                    r.removed_by_run.insert(r.removed_by_run.end(), run.blobs.begin(), run.blobs.end());
            std::sort(r.removed_by_run.begin(), r.removed_by_run.end());
        }
        return r;
    }

}  // namespace WireCell::Clus::BadBlobRuns

#endif
