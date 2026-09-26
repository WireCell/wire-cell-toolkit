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

std::vector<std::pair<int, int>> WireCell::Clus::PR::cathode_pair_candidates(
    const std::vector<int>& gids, const std::map<int, int>& gid_group, const std::map<int, int>& gid_tpc)
{
    std::vector<int> sorted(gids.begin(), gids.end());
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    auto group_of = [&gid_group](int g) {
        auto it = gid_group.find(g);
        return it == gid_group.end() ? g : it->second;
    };
    auto tpc_of = [&gid_tpc](int g) {
        auto it = gid_tpc.find(g);
        return it == gid_tpc.end() ? -1 : it->second;
    };
    std::vector<std::pair<int, int>> out;
    for (size_t i = 0; i < sorted.size(); ++i) {
        for (size_t j = i + 1; j < sorted.size(); ++j) {
            const int a = sorted[i], b = sorted[j];
            if (group_of(a) != group_of(b)) continue;
            const int ta = tpc_of(a), tb = tpc_of(b);
            if (ta < 0 || tb < 0 || ta == tb) continue;
            out.emplace_back(a, b);
        }
    }
    return out;
}

bool WireCell::Clus::PR::bundle_contact(double d, double xa, double xb, double cathode_x, double xcut, double gap)
{
    // Written as !(x < cut) so a NaN never contacts.
    if (!(d < gap)) return false;
    if (xcut > 0) {
        if (!(std::abs(xa - cathode_x) < xcut)) return false;
        if (!(std::abs(xb - cathode_x) < xcut)) return false;
    }
    return true;
}

std::map<int, int> WireCell::Clus::PR::merge_bundles(const std::vector<int>& gids,
                                                     const std::vector<std::pair<int, int>>& contacts)
{
    std::map<int, int> parent;
    for (int g : gids) parent[g] = g;
    auto find = [&parent](int g) {
        while (parent[g] != g) {
            parent[g] = parent[parent[g]];
            g = parent[g];
        }
        return g;
    };
    for (const auto& [a, b] : contacts) {
        if (!parent.count(a) || !parent.count(b)) continue;
        const int ra = find(a), rb = find(b);
        if (ra == rb) continue;
        // Root = the smaller gid, so the merged bundle's key is order-independent.
        parent[std::max(ra, rb)] = std::min(ra, rb);
    }
    std::map<int, int> out;
    for (int g : gids) out[g] = find(g);
    return out;
}
