#include "WireCellClus/ImageDeghostFunctions.h"

#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"

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

size_t WireCell::Clus::PR::harmonize_blob_arrays(WireCell::Clus::Facade::Grouping& grouping, const std::string& pcname,
                                                 std::vector<std::string>& partial)
{
    // keys every blob's PC carries; std::map of strings => sorted, deterministic
    std::map<std::string, size_t> nkey;
    size_t nblob = 0;
    for (auto* c : grouping.children()) {
        for (auto* b : c->children()) {
            auto it = b->local_pcs().find(pcname);
            if (it == b->local_pcs().end()) continue;
            ++nblob;
            for (const auto& name : it->second.keys()) ++nkey[name];
        }
    }
    partial.clear();
    for (const auto& [name, n] : nkey) if (n != nblob) partial.push_back(name);
    size_t nerased = 0;
    if (partial.empty()) return 0;
    for (auto* c : grouping.children()) nerased += WireCell::Clus::Facade::erase_blob_arrays(*c, partial, pcname);
    return nerased;
}
