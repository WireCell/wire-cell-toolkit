#include "WireCellClus/ImageDeghostFunctions.h"

#include <cmath>
#include <map>

using namespace WireCell::Clus::PR;

std::vector<ImageDeghostT0> WireCell::Clus::PR::image_deghost_cluster_t0(const std::vector<int>& cell_dgc,
                                                                         const std::vector<int>& owner_of,
                                                                         const std::vector<char>& selected,
                                                                         const std::vector<double>& owner_t0,
                                                                         int n_dgc, double tol,
                                                                         std::vector<char>& take)
{
    const size_t ncell = cell_dgc.size();
    std::vector<ImageDeghostT0> out(n_dgc > 0 ? n_dgc : 0);
    take.assign(ncell, 0);
    auto is_taker = [&](size_t i) {
        const int k = owner_of[i];
        return k >= 0 && k < (int) selected.size() && selected[k];
    };
    // cells per (cascade cluster, owner); int keys => sorted, deterministic
    std::vector<std::map<int, size_t>> counts(out.size());
    for (size_t i = 0; i < ncell; ++i) {
        if (!is_taker(i)) continue;
        const int d = cell_dgc[i];
        if (d < 0 || d >= n_dgc) continue;
        ++counts[d][owner_of[i]];
    }
    for (int d = 0; d < n_dgc; ++d) {
        if (counts[d].empty()) continue;
        int best = -1; size_t nbest = 0;
        for (const auto& [k, n] : counts[d]) {   // ascending key: strict > keeps the smaller key on a tie
            if (n > nbest) { best = k; nbest = n; }
        }
        out[d].touched = true;
        out[d].t0 = owner_t0[best];
        for (const auto& [k, n] : counts[d]) {
            if (std::abs(owner_t0[k] - out[d].t0) > tol) out[d].conflict = true;
        }
    }
    for (size_t i = 0; i < ncell; ++i) {
        if (!is_taker(i)) continue;
        const int d = cell_dgc[i];
        if (d < 0 || d >= n_dgc) continue;
        if (std::abs(owner_t0[owner_of[i]] - out[d].t0) <= tol) take[i] = 1;
    }
    return out;
}
