// Diagnostic-only integrity checker for the per-blob provenance arrays
// ("perblob" node-local PC: isolated / assoc_cluster_* / real_cluster_*).
//
// Doc 52 §11: the assoc pair came out of the per-APA -> all-APA handoff with
// the right group sizes attached to the WRONG blobs, and no existing check
// could see it -- every other per-blob array is constant within a pre-merge
// member, so a within-member permutation is invisible to them.  This checker
// makes the corruption observable at any pipeline boundary:
//
//  * length check: every array in a cluster's "perblob" PC must have exactly
//    nchildren rows (put_pcarray() can silently violate this: it assign()s
//    over an existing array with a new shape, bypassing Dataset::add()'s
//    major-axis check).
//  * scatter check: an assoc_cluster_main == 0 group is, by construction, one
//    small pre-merge cluster (clustering_isolated small/big length_cut <= 20
//    cm), so its blob centers must be spatially compact.  A group whose
//    max distance from its own centroid exceeds ~25 cm can only be rows
//    attached to the wrong blobs.
//
// Enabled ONLY by env WCT_PROV_CHECK=1 -- production runs never pay for it and
// never log from it.  Works on the raw pctree node so it can run both where
// facades exist (MultiAlgBlobClustering) and where they do not
// (PointTreeMerging, serialize/deserialize round-trips).
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Units.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <map>
#include <set>
#include <vector>

using namespace WireCell;
using WireCell::PointCloud::Tree::Points;

size_t WireCell::Clus::Facade::check_perblob_provenance(
    const Points::node_t& groot, const std::string& tag)
{
    static const bool enabled = []() {
        const char* env = std::getenv("WCT_PROV_CHECK");
        return env && env[0] && env[0] != '0';
    }();
    if (!enabled) return 0;

    static Log::logptr_t log = Log::logger("clus.provchk");
    const double scatter_cut = 25 * units::cm;

    size_t nclus = 0, nwith = 0, nbad = 0;
    double worst = -1;
    int worst_ic = -1, worst_gid = -2;
    std::array<double, 3> worst_ctr{0, 0, 0};

    int ic = -1;
    for (const auto* cnode : groot.children()) {
        ++ic;
        ++nclus;
        const auto& lpcs = cnode->value.local_pcs();
        auto it = lpcs.find("perblob");
        if (it == lpcs.end()) continue;
        ++nwith;
        const auto& ds = it->second;
        const size_t nb = cnode->nchildren();

        bool len_ok = true;
        for (const auto& key : ds.keys()) {
            const auto arr = ds.get(key);
            if (!arr) continue;
            if (arr->size_major() != nb) {
                log->error("[provchk:{}] cluster#{}: key '{}' has {} rows for {} blobs",
                           tag, ic, key, arr->size_major(), nb);
                ++nbad;
                len_ok = false;
            }
        }
        if (!len_ok) continue;

        const auto aid = ds.get("assoc_cluster_id");
        const auto amain = ds.get("assoc_cluster_main");
        if (!aid || !amain) continue;
        const auto vid = aid->elements<int>();
        const auto vmain = amain->elements<int>();
        if (vid.size() != nb || vmain.size() != nb) continue;

        // Blob centers in child (== row) order, from each blob's scalar PC.
        std::vector<std::array<double, 3>> ctr(nb);
        bool have_ctr = (nb > 0);
        const auto kids = cnode->children();
        for (size_t i = 0; i < nb && have_ctr; ++i) {
            const auto& bl = kids[i]->value.local_pcs();
            auto sit = bl.find("scalar");
            if (sit == bl.end()) { have_ctr = false; break; }
            const auto cx = sit->second.get("center_x");
            const auto cy = sit->second.get("center_y");
            const auto cz = sit->second.get("center_z");
            if (!cx || !cy || !cz) { have_ctr = false; break; }
            ctr[i] = {cx->elements<double>()[0],
                      cy->elements<double>()[0],
                      cz->elements<double>()[0]};
        }
        if (!have_ctr) continue;

        std::map<int, std::vector<size_t>> grp;
        for (size_t i = 0; i < nb; ++i) {
            if (vmain[i] == 0) grp[vid[i]].push_back(i);
        }
        for (const auto& [gid, rows] : grp) {
            if (rows.size() < 2) continue;
            std::array<double, 3> c{0, 0, 0};
            for (size_t r : rows) {
                for (int k = 0; k < 3; ++k) c[k] += ctr[r][k];
            }
            for (int k = 0; k < 3; ++k) c[k] /= rows.size();
            double r2max = 0;
            for (size_t r : rows) {
                double r2 = 0;
                for (int k = 0; k < 3; ++k) {
                    const double d = ctr[r][k] - c[k];
                    r2 += d * d;
                }
                if (r2 > r2max) r2max = r2;
            }
            const double rmax = std::sqrt(r2max);
            if (rmax > worst) {
                worst = rmax;
                worst_ic = ic;
                worst_gid = gid;
                worst_ctr = c;
            }
            if (rmax > scatter_cut) {
                log->error("[provchk:{}] cluster#{}: assoc group {} n={} rmax={:.1f} cm "
                           "centroid=({:.1f},{:.1f},{:.1f}) cm",
                           tag, ic, gid, rows.size(), rmax / units::cm,
                           c[0] / units::cm, c[1] / units::cm, c[2] / units::cm);
                ++nbad;
            }
        }
    }

    log->info("[provchk:{}] {} clusters ({} with perblob): {} problem(s); "
              "worst assoc-group rmax {:.1f} cm (cluster#{} grp {} at ({:.1f},{:.1f},{:.1f}) cm)",
              tag, nclus, nwith, nbad, worst / units::cm, worst_ic, worst_gid,
              worst_ctr[0] / units::cm, worst_ctr[1] / units::cm, worst_ctr[2] / units::cm);
    return nbad;
}

bool WireCell::Clus::Facade::realign_perblob_after_regroup(
    Cluster& cluster, const std::vector<int>& cc, const std::string& pcname)
{
    const size_t nb = cc.size();
    if (!nb) return false;

    auto& lpcs = cluster.value().local_pcs();
    auto it = lpcs.find(pcname);
    if (it == lpcs.end()) return false;
    if (it->second.size_major() != nb) return false;
    if ((size_t) cluster.nchildren() != nb) return false;

    // The post-round-trip children order: kept (cc<0) blobs first in original
    // relative order, then each non-negative group id ascending, rows in
    // original relative order within a group.
    std::vector<size_t> rows;
    rows.reserve(nb);
    for (size_t i = 0; i < nb; ++i) {
        if (cc[i] < 0) rows.push_back(i);
    }
    std::set<int> gids;
    for (int v : cc) {
        if (v >= 0) gids.insert(v);
    }
    for (int gid : gids) {
        for (size_t i = 0; i < nb; ++i) {
            if (cc[i] == gid) rows.push_back(i);
        }
    }

    it->second = it->second.subset(rows);
    return true;
}
