/** Reco <-> truth matching by shared reco points (ai-helper issue 41).
 *
 * The SPINE MatchProcessor rule, as used by M. Vicenzi's truthmatch.py, on
 * objects that live on the same reco points:
 *
 *   - every reco point carries a reco object and (through its blob's truth
 *     label) at most one truth object; an unlabeled point belongs to no
 *     truth object but still counts in its reco object's size;
 *   - IoU = |r & t| / (|r| + |t| - |r & t|), in points (n) and in charge (q);
 *   - reco -> truth takes each reco object's partners by decreasing point IoU,
 *     truth -> reco each truth object's, independently; a partner counts only
 *     above min_overlap (0 = any shared point); "mutual" where the two best
 *     partners name each other.
 *
 * Pure computation: plain structs, no ROOT, no PR graph.  The caller builds the
 * objects and the shared-point cells (SbndPrMagnifyTrackingVisitor, from the
 * T_overlap rows it writes, so a reader can redo the same match from the file).
 */
#ifndef WIRECELL_CLUS_TRUTHMATCH_H
#define WIRECELL_CLUS_TRUTHMATCH_H

#include <algorithm>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace WireCell::Clus::TruthMatch {

    /// Size of one object: points and charge.
    struct Size {
        double n{0};
        double q{0};
    };

    /// Points and charge shared by reco object `reco` and truth object `truth`.
    struct Cell {
        int reco{-1};
        int truth{-1};
        double n{0};
        double q{0};
    };

    /// One partner of an object.
    struct Partner {
        int id{-1};          // index of the partner on the other side
        double iou_n{0};     // point IoU
        double iou_q{0};     // charge IoU
        double shared_n{0};
        double shared_q{0};
    };

    struct Result {
        std::vector<std::vector<Partner>> reco;   // per reco object, by decreasing iou_n
        std::vector<std::vector<Partner>> truth;  // per truth object, by decreasing iou_n
        std::vector<char> reco_mutual;            // best partner's best partner is this object
        std::vector<char> truth_mutual;
    };

    inline double iou(double shared, double a, double b)
    {
        const double u = a + b - shared;
        return u > 0 ? shared / u : 0.0;
    }

    /// Match `nreco` x `ntruth` objects.  Cells with the same (reco, truth)
    /// are summed.  `topk` partners are kept per object (0 = all).  Ties in
    /// iou_n are broken by iou_q, then by the lower partner index (numpy's
    /// argmax keeps the first).
    inline Result match(const std::vector<Size>& reco, const std::vector<Size>& truth,
                        const std::vector<Cell>& cells, std::size_t topk = 3, double min_overlap = 0.0)
    {
        std::map<std::pair<int, int>, Size> shared;
        for (const auto& c : cells) {
            if (c.reco < 0 || c.truth < 0) continue;
            if (c.reco >= (int) reco.size() || c.truth >= (int) truth.size()) continue;
            auto& s = shared[{c.reco, c.truth}];
            s.n += c.n;
            s.q += c.q;
        }
        Result res;
        res.reco.resize(reco.size());
        res.truth.resize(truth.size());
        for (const auto& [rt, s] : shared) {
            const auto [r, t] = rt;
            const double in = iou(s.n, reco[r].n, truth[t].n);
            const double iq = iou(s.q, reco[r].q, truth[t].q);
            if (!(in > min_overlap)) continue;
            res.reco[r].push_back({t, in, iq, s.n, s.q});
            res.truth[t].push_back({r, in, iq, s.n, s.q});
        }
        auto order = [](const Partner& a, const Partner& b) {
            if (a.iou_n != b.iou_n) return a.iou_n > b.iou_n;
            if (a.iou_q != b.iou_q) return a.iou_q > b.iou_q;
            return a.id < b.id;
        };
        for (auto* side : {&res.reco, &res.truth}) {
            for (auto& v : *side) {
                std::sort(v.begin(), v.end(), order);
                if (topk > 0 && v.size() > topk) v.resize(topk);
            }
        }
        res.reco_mutual.assign(reco.size(), 0);
        res.truth_mutual.assign(truth.size(), 0);
        for (std::size_t r = 0; r < reco.size(); ++r) {
            if (res.reco[r].empty()) continue;
            const int t = res.reco[r].front().id;
            if (!res.truth[t].empty() && res.truth[t].front().id == (int) r) {
                res.reco_mutual[r] = 1;
                res.truth_mutual[t] = 1;
            }
        }
        return res;
    }

}  // namespace WireCell::Clus::TruthMatch

#endif
