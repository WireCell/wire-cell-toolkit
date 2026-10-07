#include "WireCellClus/BeamDeghostFunctions.h"

#include <algorithm>
#include <map>

namespace WireCell::Clus::PR {

    static long long span(int a0, int a1, int b0, int b1)
    {
        const int lo = std::max(a0, b0), hi = std::min(a1, b1);
        return hi > lo ? hi - lo : 0;
    }

    long long beam_deghost_overlap(const DeghostBox& a, const DeghostBox& b)
    {
        if (a.wpid != b.wpid) return 0;
        const long long t = span(a.t0, a.t1, b.t0, b.t1);
        if (!t) return 0;
        const long long u = span(a.u0, a.u1, b.u0, b.u1);
        if (!u) return 0;
        const long long v = span(a.v0, a.v1, b.v0, b.v1);
        if (!v) return 0;
        const long long w = span(a.w0, a.w1, b.w0, b.w1);
        return t * u * v * w;
    }

    std::vector<int> beam_deghost_assign(const std::vector<DeghostBox>& cells,
                                         const std::vector<DeghostBox>& owners)
    {
        // Owners per face, sorted by (t0, input index); the longest owner span
        // per face bounds how far back a cell has to look.
        struct Face {
            std::vector<size_t> idx;
            int max_span{0};
        };
        std::map<int, Face> faces;   // int-keyed => sorted
        for (size_t i = 0; i < owners.size(); ++i) {
            auto& f = faces[owners[i].wpid];
            f.idx.push_back(i);
            f.max_span = std::max(f.max_span, owners[i].t1 - owners[i].t0);
        }
        for (auto& [wpid, f] : faces) {
            std::stable_sort(f.idx.begin(), f.idx.end(),
                             [&](size_t a, size_t b) { return owners[a].t0 < owners[b].t0; });
        }

        std::vector<int> ret(cells.size(), -1);
        std::map<int, long long> sums;   // cluster_id -> summed overlap
        for (size_t ic = 0; ic < cells.size(); ++ic) {
            const auto& c = cells[ic];
            auto fit = faces.find(c.wpid);
            if (fit == faces.end()) continue;
            const auto& f = fit->second;
            sums.clear();
            auto it = std::lower_bound(f.idx.begin(), f.idx.end(), c.t0 - f.max_span,
                                       [&](size_t a, int t) { return owners[a].t0 < t; });
            for (; it != f.idx.end() && owners[*it].t0 < c.t1; ++it) {
                const long long o = beam_deghost_overlap(c, owners[*it]);
                if (o > 0) sums[owners[*it].cluster_id] += o;
            }
            long long best = 0;
            for (const auto& [cid, s] : sums) {   // ascending cluster_id: a tie keeps the smaller
                if (s > best) { best = s; ret[ic] = cid; }
            }
        }
        return ret;
    }

}  // namespace WireCell::Clus::PR
