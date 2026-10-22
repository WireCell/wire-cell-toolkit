/** Matching: reco object <-> truth node overlap metrics from one
    intersection table.

    Header-only (owner decision O1).  The table shared[reco][truth] =
    sum of the overlap (charge in electrons, or a count of blobs/points
    when every entry is added with weight 1) gives in one pass

    - purity = shared / reco total, efficiency (completeness) = shared /
      truth total, IoU = shared / (reco + truth - shared)
      (SPINE math/match.py);
    - the Pandora validation cut: shared >= min_shared, purity >= 0.5,
      completeness >= 0.1 (EventValidationBaseAlgorithm);
    - a one-to-one assignment, best IoU first, IoU >= iou_min, and the
      panoptic quality PQ = SQ x RQ on it (pimm instance metrics);
    - the adjusted Rand index of two labelings over a mask.

    wcp-porting-img pr_ml/docs/19 sec 3.6 and sec 7.4.
 */

#ifndef WIRECELLTRUTH_MATCHING
#define WIRECELLTRUTH_MATCHING

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <vector>

namespace WireCell::Truth {

    struct MatchPair {
        int reco{-1};
        int truth{-1};
        double shared{0};
        double purity{0};
        double efficiency{0};
        double iou{0};
    };

    struct PanopticQuality {
        double pq{0}, sq{0}, rq{0};
        size_t tp{0}, fp{0}, fn{0};
    };

    class Intersection {
      public:
        void add(int reco, int truth, double w)
        {
            if (w == 0.0) return;
            m_shared[reco][truth] += w;
            m_reco_tot[reco] += w;
            m_truth_tot[truth] += w;
        }
        /// A reco object with no truth overlap (ghost) or a truth node with no reco (lost).
        void add_reco_only(int reco, double w) { m_reco_tot[reco] += w; m_shared[reco]; }
        void add_truth_only(int truth, double w) { m_truth_tot[truth] += w; }

        size_t nreco() const { return m_reco_tot.size(); }
        size_t ntruth() const { return m_truth_tot.size(); }
        double reco_total(int reco) const { auto it = m_reco_tot.find(reco); return it == m_reco_tot.end() ? 0 : it->second; }
        double truth_total(int truth) const { auto it = m_truth_tot.find(truth); return it == m_truth_tot.end() ? 0 : it->second; }
        const std::map<int, double>& reco_totals() const { return m_reco_tot; }
        const std::map<int, double>& truth_totals() const { return m_truth_tot; }

        MatchPair pair(int reco, int truth) const
        {
            MatchPair p;
            p.reco = reco;
            p.truth = truth;
            auto rit = m_shared.find(reco);
            if (rit != m_shared.end()) {
                auto tit = rit->second.find(truth);
                if (tit != rit->second.end()) p.shared = tit->second;
            }
            const double rt = reco_total(reco), tt = truth_total(truth);
            p.purity = rt > 0 ? p.shared / rt : 0.0;
            p.efficiency = tt > 0 ? p.shared / tt : 0.0;
            const double uni = rt + tt - p.shared;
            p.iou = uni > 0 ? p.shared / uni : 0.0;
            return p;
        }

        /// Every non-zero pair, in (reco, truth) order.
        std::vector<MatchPair> pairs() const
        {
            std::vector<MatchPair> out;
            for (const auto& [r, ts] : m_shared) {
                for (const auto& [t, w] : ts) {
                    if (w > 0) out.push_back(pair(r, t));
                }
            }
            return out;
        }

        /// The truth node with the largest shared weight for a reco
        /// object (ties: the smaller truth id); reco -1 when none.
        MatchPair best_for_reco(int reco) const
        {
            MatchPair best;
            auto rit = m_shared.find(reco);
            if (rit == m_shared.end()) return best;
            for (const auto& [t, w] : rit->second) {
                if (w > best.shared) best = pair(reco, t);
            }
            return best;
        }
        MatchPair best_for_truth(int truth) const
        {
            MatchPair best;
            for (const auto& [r, ts] : m_shared) {
                auto tit = ts.find(truth);
                if (tit != ts.end() && tit->second > best.shared) best = pair(r, truth);
            }
            return best;
        }

        /// Pandora validation: a good match.
        static bool pandora_good(const MatchPair& p, double min_shared = 5, double purity_min = 0.5, double completeness_min = 0.1)
        {
            return p.shared >= min_shared && p.purity >= purity_min && p.efficiency >= completeness_min;
        }

        /// One-to-one, best IoU first, greedy; ties by (reco, truth) id.
        std::vector<MatchPair> one_to_one(double iou_min = 0.5) const
        {
            auto all = pairs();
            std::stable_sort(all.begin(), all.end(), [](const MatchPair& a, const MatchPair& b) {
                if (a.iou != b.iou) return a.iou > b.iou;
                if (a.reco != b.reco) return a.reco < b.reco;
                return a.truth < b.truth;
            });
            std::set<int> used_r, used_t;
            std::vector<MatchPair> out;
            for (const auto& p : all) {
                if (p.iou < iou_min) break;
                if (used_r.count(p.reco) || used_t.count(p.truth)) continue;
                used_r.insert(p.reco);
                used_t.insert(p.truth);
                out.push_back(p);
            }
            return out;
        }

        /// Panoptic quality on the one-to-one assignment.
        PanopticQuality panoptic(double iou_min = 0.5) const
        {
            PanopticQuality q;
            const auto m = one_to_one(iou_min);
            q.tp = m.size();
            q.fp = nreco() - q.tp;
            q.fn = ntruth() - q.tp;
            double siou = 0;
            for (const auto& p : m) siou += p.iou;
            q.sq = q.tp ? siou / q.tp : 0.0;
            const double den = q.tp + 0.5 * q.fp + 0.5 * q.fn;
            q.rq = den > 0 ? q.tp / den : 0.0;
            q.pq = q.sq * q.rq;
            return q;
        }

        /// Adjusted Rand index of two labelings (same length), over the
        /// elements where mask is non-zero (or all when mask is empty).
        static double ari(const std::vector<int>& a, const std::vector<int>& b, const std::vector<int>& mask = {})
        {
            std::map<std::pair<int,int>, double> nij;
            std::map<int, double> ai, bj;
            double n = 0;
            for (size_t i = 0; i < a.size() && i < b.size(); ++i) {
                if (!mask.empty() && i < mask.size() && mask[i] == 0) continue;
                nij[{a[i], b[i]}] += 1;
                ai[a[i]] += 1;
                bj[b[i]] += 1;
                n += 1;
            }
            auto c2 = [](double x) { return 0.5 * x * (x - 1); };
            double sum_ij = 0, sum_a = 0, sum_b = 0;
            for (const auto& [k, v] : nij) sum_ij += c2(v);
            for (const auto& [k, v] : ai) sum_a += c2(v);
            for (const auto& [k, v] : bj) sum_b += c2(v);
            const double cn = c2(n);
            if (cn <= 0) return 1.0;
            const double expected = sum_a * sum_b / cn;
            const double maxi = 0.5 * (sum_a + sum_b);
            if (maxi == expected) return 1.0;
            return (sum_ij - expected) / (maxi - expected);
        }

      private:
        std::map<int, std::map<int, double>> m_shared;
        std::map<int, double> m_reco_tot, m_truth_tot;
    };
}

#endif
