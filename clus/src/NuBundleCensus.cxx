#include "WireCellClus/NuBundleCensus.h"

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace WireCell::Clus::PR;

std::vector<int> WireCell::Clus::PR::group_flashes(const std::vector<int>& gid,
                                                   const std::vector<int>& tpc,
                                                   const std::vector<double>& time_us,
                                                   double dt_us)
{
    const size_t n = std::min({gid.size(), tpc.size(), time_us.size()});
    std::vector<size_t> parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&parent](size_t i) {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    };
    if (dt_us > 0) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                if (tpc[i] == tpc[j]) continue;
                // Written as !(x < dt) so a NaN time never groups.
                if (!(std::abs(time_us[i] - time_us[j]) < dt_us)) continue;
                const size_t ri = find(i), rj = find(j);
                if (ri != rj) parent[std::max(ri, rj)] = std::min(ri, rj);
            }
        }
    }
    // Group id = the smallest gid among the members (independent of input order).
    std::vector<int> min_gid(n, 0);
    std::vector<bool> seen(n, false);
    for (size_t i = 0; i < n; ++i) {
        const size_t r = find(i);
        if (!seen[r] || gid[i] < min_gid[r]) {
            min_gid[r] = gid[i];
            seen[r] = true;
        }
    }
    std::vector<int> out(n);
    for (size_t i = 0; i < n; ++i) out[i] = min_gid[find(i)];
    return out;
}
