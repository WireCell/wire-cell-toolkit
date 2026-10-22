/** DepoPolyline: the ordered trajectory of one G4 track from its depos.

    Header-only (owner decision O1: clus visitors may include it without
    a link edge to WireCellTruth).

    The polyline vertices are the positions of the undrifted priors of a
    track's depos, transformed into the reconstructed frame by the
    caller (x from the drifted depo's arrival time through time2drift,
    y and z of the depo), ordered by the priors' true time.  It answers
    the per-point questions of doc 15 sec 2.1: distance to the
    trajectory, local tangent, arc-length position, nearness to an end
    point, and the kinks (direction changes over a +-window arc length),
    which become the v_heat_* keypoint fields.

    Distances use a k-d tree over the segment midpoints
    (WireCellUtil/NFKDVec.h) followed by the exact point-to-segment
    distance on the k nearest segments; k = 8 by default, exact for any
    k when the polyline has at most k segments.

    wcp-porting-img pr_ml/docs/19 sec 3.6.
 */

#ifndef WIRECELLTRUTH_DEPOPOLYLINE
#define WIRECELLTRUTH_DEPOPOLYLINE

#include "WireCellUtil/NFKDVec.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

namespace WireCell::Truth {

    class DepoPolyline {
      public:
        struct Vertex {
            Point pos;
            double t{0};        // ordering key (the prior's time)
            double energy{0};   // energy of the depo at this vertex (any unit; dedx uses it per length)
        };
        enum class KeyKind : int { start = 0, end = 1, kink = 2 };
        struct Keypoint {
            KeyKind kind;
            Point pos;
            double s{0};        // arc length from the start
            double angle{0};    // kink angle (rad), 0 for start/end
        };
        struct Hit {
            bool valid{false};
            double dist{std::numeric_limits<double>::max()};
            size_t seg{0};      // segment index: vertices seg, seg+1
            double s{0};        // arc length of the closest point
            Vector tangent;     // unit direction of the segment (start -> end)
            Point closest;
            double dedx{0};     // energy per length of the segment
        };

        DepoPolyline() = default;

        /// Add a vertex; build() orders them.
        void add(const Point& pos, double t, double energy = 0.0)
        {
            m_verts.push_back({pos, t, energy});
            m_built = false;
        }
        size_t nvertices() const { return m_verts.size(); }
        size_t nsegments() const { return m_verts.size() > 1 ? m_verts.size() - 1 : 0; }
        double length() const { return m_s.empty() ? 0.0 : m_s.back(); }
        const std::vector<Vertex>& vertices() const { return m_verts; }
        const Point& start() const { return m_verts.front().pos; }
        const Point& end() const { return m_verts.back().pos; }

        /// Order by time (stable), drop zero-length steps, compute arc
        /// lengths, tangents and the midpoint tree.  Idempotent.
        void build()
        {
            if (m_built) return;
            std::stable_sort(m_verts.begin(), m_verts.end(),
                             [](const Vertex& a, const Vertex& b) { return a.t < b.t; });
            // merge coincident consecutive vertices, keeping their energy
            std::vector<Vertex> kept;
            for (const auto& v : m_verts) {
                if (!kept.empty() && (kept.back().pos - v.pos).magnitude() <= 0.0) {
                    kept.back().energy += v.energy;
                    continue;
                }
                kept.push_back(v);
            }
            m_verts.swap(kept);
            const size_t n = m_verts.size();
            m_s.assign(n, 0.0);
            m_tan.assign(nsegments(), Vector(0, 0, 0));
            m_seglen.assign(nsegments(), 0.0);
            for (size_t i = 0; i + 1 < n; ++i) {
                const Vector d = m_verts[i + 1].pos - m_verts[i].pos;
                const double l = d.magnitude();
                m_seglen[i] = l;
                m_tan[i] = l > 0 ? d / l : Vector(0, 0, 0);
                m_s[i + 1] = m_s[i] + l;
            }
            m_tree.reset();
            if (nsegments() > 0) {
                NFKDVec::Tree<double, NFKDVec::IndexStatic>::points_type mids(3);
                for (size_t i = 0; i < nsegments(); ++i) {
                    const Point m = 0.5 * (m_verts[i].pos + m_verts[i + 1].pos);
                    mids[0].push_back(m.x());
                    mids[1].push_back(m.y());
                    mids[2].push_back(m.z());
                }
                m_tree = std::make_unique<NFKDVec::Tree<double, NFKDVec::IndexStatic>>(mids);
            }
            m_built = true;
        }

        /// Nearest point on the polyline; a one-vertex polyline answers
        /// with the vertex (seg 0, zero tangent).
        Hit nearest(const Point& p, size_t kcand = 8) const
        {
            Hit hit;
            if (m_verts.empty() || !m_built) return hit;
            if (nsegments() == 0) {
                hit.valid = true;
                hit.closest = m_verts[0].pos;
                hit.dist = (p - hit.closest).magnitude();
                hit.dedx = 0;
                return hit;
            }
            std::vector<size_t> cands;
            if (nsegments() <= kcand) {
                cands.resize(nsegments());
                std::iota(cands.begin(), cands.end(), 0);
            }
            else {
                const std::vector<double> q{p.x(), p.y(), p.z()};
                for (const auto& [ind, d2] : m_tree->knn(kcand, q)) cands.push_back(ind);
            }
            for (size_t seg : cands) {
                const auto& a = m_verts[seg].pos;
                const double l = m_seglen[seg];
                double u = 0.0;
                if (l > 0) {
                    u = (p - a).dot(m_tan[seg]) / l;
                    u = std::min(1.0, std::max(0.0, u));
                }
                const Point c = a + m_tan[seg] * (u * l);
                const double d = (p - c).magnitude();
                // ties: the earlier segment wins, whatever the candidate order
                if (d < hit.dist || (d == hit.dist && hit.valid && seg < hit.seg)) {
                    hit.valid = true;
                    hit.dist = d;
                    hit.seg = seg;
                    hit.s = m_s[seg] + u * l;
                    hit.tangent = m_tan[seg];
                    hit.closest = c;
                    hit.dedx = l > 0 ? (m_verts[seg].energy + m_verts[seg + 1].energy) * 0.5 / l : 0.0;
                }
            }
            return hit;
        }

        /// Start, end and the kinks: a vertex where the direction of
        /// the chord `window` behind and the chord `window` ahead differ
        /// by more than angle_min and the angle is the local maximum
        /// within +-window.  Positions are floating point.
        std::vector<Keypoint> keypoints(double window = 3 * units::cm, double angle_min = 15 * units::deg) const
        {
            std::vector<Keypoint> out;
            if (m_verts.empty() || !m_built) return out;
            out.push_back({KeyKind::start, m_verts.front().pos, 0.0, 0.0});
            const size_t n = m_verts.size();
            if (n < 3 || length() <= 0) {
                if (n > 1) out.push_back({KeyKind::end, m_verts.back().pos, length(), 0.0});
                return out;
            }
            std::vector<double> ang(n, 0.0);
            for (size_t i = 1; i + 1 < n; ++i) {
                size_t j = i, k = i;
                while (j > 0 && m_s[i] - m_s[j] < window) --j;
                while (k + 1 < n && m_s[k] - m_s[i] < window) ++k;
                if (m_s[i] - m_s[j] < 0.5 * window || m_s[k] - m_s[i] < 0.5 * window) continue; // too close to an end
                const Vector b = (m_verts[i].pos - m_verts[j].pos).norm();
                const Vector a = (m_verts[k].pos - m_verts[i].pos).norm();
                const double c = std::max(-1.0, std::min(1.0, a.dot(b)));
                ang[i] = std::acos(c);
            }
            for (size_t i = 1; i + 1 < n; ++i) {
                if (ang[i] < angle_min) continue;
                bool best = true;
                for (size_t j = 1; j + 1 < n && best; ++j) {
                    if (j == i || std::abs(m_s[j] - m_s[i]) > window) continue;
                    if (ang[j] > ang[i] || (ang[j] == ang[i] && j < i)) best = false;
                }
                if (best) out.push_back({KeyKind::kink, m_verts[i].pos, m_s[i], ang[i]});
            }
            out.push_back({KeyKind::end, m_verts.back().pos, length(), 0.0});
            return out;
        }

        /// Gaussian "heat" of a point for a keypoint kind: max over the
        /// keypoints of that kind of exp(-d^2 / 2 sigma^2).
        static double heat(const Point& p, const std::vector<Keypoint>& kps, KeyKind kind, double sigma = 1 * units::cm)
        {
            double best = 0.0;
            for (const auto& kp : kps) {
                if (kp.kind != kind) continue;
                const double d2 = (p - kp.pos).magnitude2();
                best = std::max(best, std::exp(-0.5 * d2 / (sigma * sigma)));
            }
            return best;
        }

      private:
        std::vector<Vertex> m_verts;
        std::vector<double> m_s, m_seglen;
        std::vector<Vector> m_tan;
        std::unique_ptr<NFKDVec::Tree<double, NFKDVec::IndexStatic>> m_tree;
        bool m_built{false};
    };
}

#endif
