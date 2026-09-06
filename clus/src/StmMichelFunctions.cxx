// doc pdvd/48 -- see the header.
#include "WireCellClus/StmMichelFunctions.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRTrajectoryView.h"
#include "WireCellClus/PRCommon.h"
#include "WireCellAux/ParticleInfo.h"

#include <boost/graph/graph_traits.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>

using namespace WireCell;
using namespace WireCell::Clus::PR;

WireCell::Point WireCell::Clus::PR::stm_michel_vertex_point(const VertexPtr& vtx)
{
    if (!vtx) return WireCell::Point(0, 0, 0);
    return vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
}

namespace {
    // Dijkstra from `from` over the whole reachable component, keyed on the
    // graph index (VertexIndexCmp), tie-broken on the edge index, so the
    // result does not depend on allocation order.  `stop_at` (may be null)
    // ends the search early once that vertex is settled.
    struct StmDijkstra {
        std::map<VertexPtr, double, VertexIndexCmp> dist;
        std::map<VertexPtr, std::pair<SegmentPtr, VertexPtr>, VertexIndexCmp> prev;
    };
    StmDijkstra stm_dijkstra(Graph& g, VertexPtr from, VertexPtr stop_at)
    {
        StmDijkstra out;
        if (!from || !from->descriptor_valid()) return out;
        // frontier: (distance, vertex index) -> vertex
        std::set<std::pair<double, size_t>> frontier;
        std::map<size_t, VertexPtr> by_index;

        out.dist[from] = 0.0;
        frontier.insert({0.0, from->get_graph_index()});
        by_index[from->get_graph_index()] = from;

        while (!frontier.empty()) {
            auto it = frontier.begin();
            const double d = it->first;
            VertexPtr v = by_index.at(it->second);
            frontier.erase(it);
            if (stop_at && v == stop_at) break;
            auto dv_it = out.dist.find(v);
            if (dv_it != out.dist.end() && d > dv_it->second) continue;  // stale entry

            for (auto e : sorted_out_edges(v->get_descriptor(), g)) {
                SegmentPtr seg = g[e].segment;
                if (!seg) continue;
                VertexPtr other = find_other_vertex(g, seg, v);
                if (!other) continue;
                const double w = std::max(0.0, segment_track_length(seg));
                const double nd = d + w;
                auto od = out.dist.find(other);
                if (od == out.dist.end() || nd < od->second) {
                    if (od != out.dist.end()) frontier.erase({od->second, other->get_graph_index()});
                    out.dist[other] = nd;
                    out.prev[other] = {seg, v};
                    frontier.insert({nd, other->get_graph_index()});
                    by_index[other->get_graph_index()] = other;
                }
            }
        }
        return out;
    }
}

VertexPtr WireCell::Clus::PR::stm_michel_farthest_vertex(Graph& g, VertexPtr from,
                                                         const std::function<bool(const VertexPtr&)>& accept)
{
    if (!from || !from->descriptor_valid()) return nullptr;
    auto dj = stm_dijkstra(g, from, nullptr);
    VertexPtr best; double best_d = -1;
    for (const auto& [v, d] : dj.dist) {      // index-ordered map: ties -> lowest index
        if (v == from) continue;
        if (accept && !accept(v)) continue;
        if (d > best_d) { best_d = d; best = v; }
    }
    return best;
}

std::vector<SegmentPtr> WireCell::Clus::PR::stm_michel_shortest_chain(Graph& g, VertexPtr from, VertexPtr to)
{
    std::vector<SegmentPtr> out;
    if (!from || !to || !from->descriptor_valid() || !to->descriptor_valid()) return out;
    if (from == to) return out;

    auto dj = stm_dijkstra(g, from, to);
    auto& dist = dj.dist;
    auto& prev = dj.prev;

    if (!dist.count(to)) return out;
    // Unwind.
    VertexPtr cur = to;
    while (cur != from) {
        auto p = prev.find(cur);
        if (p == prev.end()) { out.clear(); return out; }
        out.push_back(p->second.first);
        cur = p->second.second;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

std::vector<VertexPtr> WireCell::Clus::PR::stm_michel_chain_vertices(Graph& g, const std::vector<SegmentPtr>& chain, VertexPtr entry)
{
    std::vector<VertexPtr> vtxs;
    if (!entry) return vtxs;
    vtxs.push_back(entry);
    VertexPtr cur = entry;
    for (const auto& seg : chain) {
        VertexPtr next = find_other_vertex(g, seg, cur);
        if (!next) { vtxs.clear(); return vtxs; }
        vtxs.push_back(next);
        cur = next;
    }
    return vtxs;
}

StmMichelProfile WireCell::Clus::PR::stm_michel_profile(Graph& g, const std::vector<SegmentPtr>& chain, VertexPtr entry)
{
    StmMichelProfile prof;
    if (!entry || chain.empty()) return prof;

    VertexPtr cur = entry;
    double L = 0;
    bool have_prev = false;
    WireCell::Point prev_pt;
    for (size_t si = 0; si < chain.size(); ++si) {
        const auto& seg = chain[si];
        VertexPtr next = find_other_vertex(g, seg, cur);
        const auto& fits = seg->fits();
        if (!fits.empty()) {
            const auto vp = stm_michel_vertex_point(cur);
            const double d_front = (fits.front().point - vp).magnitude();
            const double d_back  = (fits.back().point - vp).magnitude();
            const bool forward = d_front <= d_back;
            const int n = static_cast<int>(fits.size());
            for (int k = 0; k < n; ++k) {
                const Fit& f = fits[forward ? k : (n - 1 - k)];
                if (f.dx <= 0 || f.dQ < 0) continue;
                if (have_prev) L += (f.point - prev_pt).magnitude();
                prev_pt = f.point;
                have_prev = true;
                prof.L.push_back(L);
                prof.dQdx.push_back(f.dQ / (f.dx / units::cm));
                prof.pts.push_back(f.point);
                prof.seg_idx.push_back(static_cast<int>(si));
            }
        }
        if (!next) break;
        cur = next;
    }
    prof.total_length = L;
    prof.rr.resize(prof.L.size());
    for (size_t i = 0; i < prof.L.size(); ++i) prof.rr[i] = L - prof.L[i];
    return prof;
}

StmMichelProfile WireCell::Clus::PR::stm_michel_profile_live(const StmMichelProfile& prof, double min_dqdx, int& n_dead)
{
    StmMichelProfile out;
    out.total_length = prof.total_length;
    n_dead = 0;
    for (size_t i = 0; i < prof.L.size(); ++i) {
        if (prof.dQdx[i] < min_dqdx) { ++n_dead; continue; }
        out.L.push_back(prof.L[i]); out.dQdx.push_back(prof.dQdx[i]); out.rr.push_back(prof.rr[i]);
        out.pts.push_back(prof.pts[i]); out.seg_idx.push_back(prof.seg_idx[i]);
    }
    return out;
}

double WireCell::Clus::PR::stm_michel_median(std::vector<double> v)
{
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    if (n % 2) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

StmMichelBragg WireCell::Clus::PR::stm_michel_bragg_contrast(const StmMichelProfile& prof,
                                                              const std::function<double(double)>& mu_dqdx_at_rr_cm,
                                                              double tail_lo, double tail_hi,
                                                              double pl_lo, double pl_hi)
{
    StmMichelBragg b;
    if (prof.empty()) return b;
    if (prof.total_length < pl_hi) {
        b.short_track = true;
        pl_lo *= 0.5;
        pl_hi *= 0.5;
    }
    std::vector<double> tail, plateau, tail_ref, plateau_ref;
    for (size_t i = 0; i < prof.rr.size(); ++i) {
        const double rr = prof.rr[i];
        const double ref = mu_dqdx_at_rr_cm ? mu_dqdx_at_rr_cm(rr / units::cm) : 0.0;
        if (rr >= tail_lo && rr <= tail_hi) { tail.push_back(prof.dQdx[i]); tail_ref.push_back(ref); }
        if (rr >= pl_lo && rr <= pl_hi)     { plateau.push_back(prof.dQdx[i]); plateau_ref.push_back(ref); }
    }
    b.n_tail = static_cast<int>(tail.size());
    b.n_plateau = static_cast<int>(plateau.size());
    b.tail_med = stm_michel_median(tail);
    b.plateau_med = stm_michel_median(plateau);
    if (b.n_tail >= 3 && b.n_plateau >= 3 && b.plateau_med > 0) {
        b.valid = true;
        b.contrast = b.tail_med / b.plateau_med;
        const double pref = stm_michel_median(plateau_ref);
        b.expected = pref > 0 ? stm_michel_median(tail_ref) / pref : 0.0;
    }
    return b;
}

bool WireCell::Clus::PR::stm_michel_seg_is_shower(const SegmentPtr& seg)
{
    if (!seg) return false;
    return seg->flags_any(SegmentFlags::kShowerTrajectory) ||
           seg->flags_any(SegmentFlags::kShowerTopology) ||
           (seg->has_particle_info() && std::abs(seg->particle_info()->pdg()) == 11);
}

namespace {
    StmMichelArm measure_arm(Graph& g, SegmentPtr ref_seg, SegmentPtr arm, VertexPtr vtx,
                             const StmMichelArmThresholds& th, double far_cap)
    {
        StmMichelArm a;
        a.seg = arm;
        if (!arm) return a;
        a.len = std::max(0.0, segment_track_length(arm));
        a.mip = th.mip_dqdx_median > 0 ? segment_median_dQ_dx(arm) / th.mip_dqdx_median : 0.0;
        a.shower_like = stm_michel_seg_is_shower(arm);
        const auto vp = stm_michel_vertex_point(vtx);
        a.kink_deg = ref_seg ? segment_pair_kink_deg(ref_seg, arm, vp, th.dir_window) : -1.0;
        VertexPtr far = find_other_vertex(g, arm, vtx);
        if (far && far->descriptor_valid()) {
            a.terminal = boost::out_degree(far->get_descriptor(), g) <= 1;
            a.far_len = a.terminal ? 0.0 : segment_far_subtree_track_length(g, far, arm, far_cap);
        }
        else {
            a.terminal = true;
        }
        return a;
    }
}

StmMichelArm WireCell::Clus::PR::stm_michel_classify_stop_arm(Graph& g, SegmentPtr last_muon, SegmentPtr arm,
                                                              VertexPtr stop, const StmMichelArmThresholds& th)
{
    StmMichelArm a = measure_arm(g, last_muon, arm, stop, th, th.michel_max_len);
    if (!arm) return a;
    const bool kink_ok = a.kink_deg >= 0;
    // A collinear MIP-like arm past the stop is the muon going on, whatever
    // the track/shower separation stamped on it: doc pdhd/03 sec 6 found the
    // shower flag on 20-24 cm arms at 2-3 deg and 1.2-1.3 MIP (029107/21
    // cluster 116, 029107/28 cluster 35), which doc pdvd/48's
    // "!shower_like" guard then handed to the Michel clause as 43-57 MeV
    // electrons.  The shower flag is not consulted here.
    if (kink_ok && a.kink_deg < th.continuation_max_angle_deg &&
        a.len > th.continuation_min_len &&
        a.mip >= th.continuation_mip_lo && a.mip <= th.continuation_mip_hi) {
        a.kind = StmMichelArm::kContinuation;
        return a;
    }
    // A Michel is a TURN (or a shower-flagged arm) of MIP-or-below charge and
    // Michel-sized reach.  It must NOT qualify through low dQ/dx alone: doc
    // pdvd/42 sec 4.4 measured the PDVD leftover past the tagger's stop as
    // collinear muon continuation at ~0.9 MIP, which a dQ/dx-only clause would
    // call a Michel (first census: a fifth of the "Michels" had kink < 30 deg).
    // doc pdhd/03 sec 6.8: the shower flag alone admitted the muon's own last
    // 5 cm (a collinear 1.6 MIP stub at 4 deg, 029107/1 cluster 113) as a
    // "Michel"; with michel_shower_min_kink_deg >= 0 a shower-flagged arm
    // whose kink is measurable must also turn by at least that much.
    const bool shower_admits = a.shower_like &&
        (th.michel_shower_min_kink_deg < 0 || !kink_ok || a.kink_deg >= th.michel_shower_min_kink_deg);
    if (a.len + a.far_len <= th.michel_max_len && a.mip > th.michel_mip_lo && a.mip < th.michel_mip_hi &&
        (shower_admits || (kink_ok && a.kink_deg >= th.michel_min_kink_deg))) {
        a.kind = StmMichelArm::kMichel;
        return a;
    }
    a.kind = StmMichelArm::kOther;
    return a;
}

StmMichelArm WireCell::Clus::PR::stm_michel_classify_chain_arm(Graph& g, SegmentPtr in_seg, SegmentPtr arm,
                                                               VertexPtr vtx, const StmMichelArmThresholds& th)
{
    StmMichelArm a = measure_arm(g, in_seg, arm, vtx, th, th.delta_max_len);
    if (!arm) return a;
    if (a.len <= th.delta_max_len && a.far_len <= th.delta_max_len && a.terminal) {
        a.kind = StmMichelArm::kDelta;
        return a;
    }
    if (a.len > th.delta_max_len && a.mip > th.hadron_mip) {
        a.kind = StmMichelArm::kHadron;
        return a;
    }
    a.kind = StmMichelArm::kOther;
    return a;
}
