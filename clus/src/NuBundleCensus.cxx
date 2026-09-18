#include "WireCellClus/NuBundleCensus.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>

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

std::vector<std::size_t> WireCell::Clus::PR::dedup_flash_groups(const std::vector<int>& gid,
                                                                const std::map<int, int>& gid_group)
{
    std::set<int> kept_groups;
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < gid.size(); ++i) {
        auto it = gid_group.find(gid[i]);
        const int grp = (it == gid_group.end()) ? gid[i] : it->second;
        if (!kept_groups.insert(grp).second) continue;   // a longer candidate already holds this flash
        keep.push_back(i);
    }
    return keep;
}
