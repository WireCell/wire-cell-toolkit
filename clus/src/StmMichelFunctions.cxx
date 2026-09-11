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

std::vector<VertexPtr> WireCell::Clus::PR::stm_michel_reachable_vertices(Graph& g, VertexPtr from)
{
    std::vector<VertexPtr> out;
    if (!from || !from->descriptor_valid()) return out;
    auto dj = stm_dijkstra(g, from, nullptr);
    for (const auto& [v, d] : dj.dist) { (void)d; out.push_back(v); }
    // explicit, whatever the map's own order
    std::sort(out.begin(), out.end(),
              [](const VertexPtr& a, const VertexPtr& b) { return a->get_graph_index() < b->get_graph_index(); });
    return out;
}

std::pair<VertexPtr, double> WireCell::Clus::PR::stm_michel_closest_vertex_of(const std::vector<VertexPtr>& cands,
                                                                              const WireCell::Point& pt,
                                                                              const std::function<bool(const VertexPtr&)>& accept)
{
    VertexPtr best;
    double best_d = 1e9;
    for (const auto& v : cands) {
        if (!v) continue;
        if (accept && !accept(v)) continue;
        const double d = (v->wcpt().point - pt).magnitude();
        if (d < best_d) { best_d = d; best = v; }
    }
    return {best, best_d};
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

// doc pdhd/15.  See the header for the citation and the calibration.
double WireCell::Clus::PR::stm_michel_charge_to_energy(double dQ_electrons, double recom_factor,
                                                       double fudge_factor, double w_value_ev)
{
    if (!std::isfinite(dQ_electrons) || dQ_electrons <= 0) return 0.0;
    if (!(recom_factor > 0) || !(fudge_factor > 0) || !(w_value_ev > 0)) return 0.0;
    return dQ_electrons / recom_factor / fudge_factor * w_value_ev / 1e6 * units::MeV;
}

// doc pdhd/17 sec 9.  See the header for the MIP-equivalence caveat.
double WireCell::Clus::PR::stm_michel_charge_to_energy_model(
    double dQ_electrons, const IRecombinationModel::pointer& model, double dedx_mev_per_cm)
{
    if (!std::isfinite(dQ_electrons) || dQ_electrons <= 0) return 0.0;
    if (!model) return 0.0;
    if (!std::isfinite(dedx_mev_per_cm) || dedx_mev_per_cm <= 0) return 0.0;
    // dx cancels between the two factors below; 1 cm keeps every intermediate
    // in the range the practical-unit models were written for (doc 88).
    const double dx = 1.0 * units::cm;
    const double dE = dedx_mev_per_cm * units::MeV / units::cm * dx;
    const double dQ_model = (*model)(dE, dx);          // electrons for that dE over dx
    if (!std::isfinite(dQ_model) || dQ_model <= 0) return 0.0;
    return dQ_electrons * dE / dQ_model;
}

// doc pdvd/81.  See the header; the arithmetic is kine_charge_from_maps's,
// statement for statement, so a symmetric input reproduces it bit for bit.
double WireCell::Clus::PR::stm_michel_combine_planes(const std::array<double, 3>& sums,
                                                     const std::array<double, 3>& weights,
                                                     double asym_switch, int* dropped_plane)
{
    if (dropped_plane) *dropped_plane = -1;
    int min_idx = 0, max_idx = 0, med_idx = 0;
    double min_q = 1e9, max_q = -1e9;
    for (int i = 0; i < 3; ++i) {
        if (sums[i] < min_q) { min_q = sums[i]; min_idx = i; }
        if (sums[i] > max_q) { max_q = sums[i]; max_idx = i; }
    }
    if (min_idx != max_idx) {
        for (int i = 0; i < 3; ++i) {
            if (i != min_idx && i != max_idx) { med_idx = i; break; }
        }
    }
    else {
        min_idx = 0; med_idx = 1; max_idx = 2;
    }
    const double weight_sum = weights[0] + weights[1] + weights[2];
    double max_asy = 0;
    if (sums[med_idx] + sums[max_idx] > 0)
        max_asy = std::abs(sums[med_idx] - sums[max_idx]) / (sums[med_idx] + sums[max_idx]);
    double overall = 0;
    if (weight_sum > 0)
        overall = (weights[0]*sums[0] + weights[1]*sums[1] + weights[2]*sums[2]) / weight_sum;
    if (max_asy > asym_switch) {
        const double pair_sum = weights[med_idx] + weights[min_idx];
        if (pair_sum > 0) {
            overall = (weights[med_idx]*sums[med_idx] + weights[min_idx]*sums[min_idx]) / pair_sum;
            if (dropped_plane) *dropped_plane = max_idx;
        }
    }
    return overall;
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

namespace {
    // 3-point running median over consecutive entries of v, same convention
    // as pdhd/stm_michel_scan/census_lib.py:running_median (which this
    // mirrors so the offline probe and the C++ agree on one definition).
    std::vector<double> running_median3(const std::vector<double>& v)
    {
        std::vector<double> out(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            const size_t lo = (i >= 1) ? i - 1 : 0;
            const size_t hi = std::min(v.size() - 1, i + 1);
            std::vector<double> w(v.begin() + lo, v.begin() + hi + 1);
            std::sort(w.begin(), w.end());
            out[i] = w[w.size() / 2];
        }
        return out;
    }

    // The reference plateau: median dQ/dx over [lo, hi] (halved for a chain
    // shorter than hi, same short-track convention stm_michel_bragg_contrast
    // uses), live points only.  Shared by stm_michel_stop_retreat and
    // stm_michel_stop_split so the two mechanisms cannot drift on what
    // "collapsed" and "Bragg rise" are judged against.  Returns 0 (never
    // negative) when there are too few live points to judge.
    double profile_plateau(const StmMichelProfile& prof, double lo, double hi, double min_dqdx_live)
    {
        if (prof.total_length < hi) { lo *= 0.5; hi *= 0.5; }
        std::vector<double> pl;
        for (size_t i = 0; i < prof.rr.size(); ++i) {
            if (prof.dQdx[i] < min_dqdx_live) continue;
            if (prof.rr[i] >= lo && prof.rr[i] <= hi) pl.push_back(prof.dQdx[i]);
        }
        if (pl.size() < 3) return 0;
        const double plateau = stm_michel_median(pl);
        return plateau > 0 ? plateau : 0;
    }
}

StmMichelRetreat WireCell::Clus::PR::stm_michel_stop_retreat(const StmMichelProfile& prof, int n_chain_segs,
                                                              const StmMichelRetreatThresholds& th)
{
    StmMichelRetreat out;
    if (th.max_drop <= 0 || prof.empty() || n_chain_segs <= 1) return out;

    // The reference plateau: SAME window and short-track halving as
    // stm_michel_bragg_contrast, computed ONCE from the full chain so the
    // reference does not chase the shrinking frame.  Live points only (a
    // dead stretch must not fake a low plateau or a fake collapse).
    const double plateau = profile_plateau(prof, th.plateau_lo, th.plateau_hi, th.min_dqdx_live);
    if (!(plateau > 0)) return out;   // cannot judge
    out.plateau = plateau;

    for (int n_drop = 1; n_drop <= th.max_drop && n_chain_segs - n_drop > 0; ++n_drop) {
        const int cutoff = n_chain_segs - n_drop;   // segments < cutoff are kept
        // the boundary L: earliest point of the segment about to be dropped
        // (ALL points, live or dead -- the geometry does not depend on charge)
        double boundary_L = -1;
        size_t boundary_i = 0;       // doc pdvd/82: the row itself, for the bend at it
        for (size_t i = 0; i < prof.L.size(); ++i) {
            if (prof.seg_idx[i] >= cutoff) { boundary_L = prof.L[i]; boundary_i = i; break; }
        }
        if (boundary_L < 0) break;   // no point carries that seg_idx -> nothing to drop
        const double drop_len = prof.total_length - boundary_L;
        if (drop_len > th.max_drop_len) break;   // only grows with n_drop; further tries are worse

        std::vector<double> tail;
        for (size_t i = 0; i < prof.L.size(); ++i) {
            if (prof.L[i] < boundary_L) continue;
            if (th.tail_strict && !(prof.L[i] > boundary_L)) continue;              // doc pdvd/74 (P3)
            if (!th.tail_sublive && prof.dQdx[i] < th.min_dqdx_live) continue;
            tail.push_back(prof.dQdx[i]);
        }
        if (static_cast<int>(tail.size()) < th.min_tail_pts) break;   // cannot judge this candidate
        const double tail_med = stm_michel_median(tail);
        // doc pdvd/82: the doc 57 reading is the FIRST clause and exits at the
        // same point it always did.  Only with tail_peak_frac > 0 does control
        // reach the peak-relative alternative below, so the knob-off path runs
        // the doc 57 test alone -- byte-identical by construction, not by
        // re-derivation.
        const bool collapsed_plateau = tail_med < th.collapse_frac * plateau;
        if (!collapsed_plateau && !(th.tail_peak_frac > 0)) break;   // not a collapsed tail

        // The profile that would SURVIVE this drop: live points with L <=
        // boundary_L, in ascending-L order (new_rr = boundary_L - L is then
        // descending).  The running median is taken over the whole surviving
        // sequence, exactly as census_lib.shape() computes it over the whole
        // array before selecting a window -- so a peak just outside the
        // window still borrows its neighbour correctly.
        std::vector<double> kept_q; std::vector<double> kept_new_rr;
        for (size_t i = 0; i < prof.L.size(); ++i) {
            if (prof.L[i] > boundary_L) continue;
            if (prof.dQdx[i] < th.min_dqdx_live) continue;
            kept_q.push_back(prof.dQdx[i]);
            kept_new_rr.push_back(boundary_L - prof.L[i]);
        }
        if (kept_q.size() < 3) break;
        const auto rm = running_median3(kept_q);
        double peak = -1;
        int n_win = 0;
        for (size_t i = 0; i < rm.size(); ++i) {
            if (kept_new_rr[i] > th.peak_window) continue;
            ++n_win;
            peak = std::max(peak, rm[i]);
        }
        if (n_win < 3) break;   // window too sparse to judge
        if (!(peak >= th.peak_frac * plateau)) break;   // nothing to retreat TO

        // doc pdvd/82: reached with collapsed_plateau false only when the knob
        // is on.  The tail is then judged against the PEAK, and the row's own
        // bend must carry the move -- charge shape alone cannot tell a Michel
        // continuing past the Bragg peak from a muon that simply keeps going
        // (the tail sits at 0.5-1.9 x plateau in both cases).
        double kink = -1;
        if (!collapsed_plateau) {
            kink = stm_michel_row_kink_deg(prof, boundary_i, th.dir_window);
            if (!(tail_med <= th.tail_peak_frac * peak && kink >= th.tail_peak_kink_min)) break;
        }

        out.n_drop = n_drop;
        out.drop_len = drop_len;
        out.last_tail_med = tail_med;
        out.last_peak = peak;
        out.by_tail_peak = !collapsed_plateau;
        out.last_kink_deg = kink;
    }
    return out;
}

double WireCell::Clus::PR::stm_michel_row_kink_deg(const StmMichelProfile& prof, size_t i, double window)
{
    if (i >= prof.L.size()) return -1;
    const double Li = prof.L[i];
    // prof.L is non-decreasing along the array index (it is the walked
    // polyline arclength), so a fixed-arclength arm on either side of `i` is
    // found by walking the index outward, not by searching rr (rr runs the
    // opposite way but is not otherwise guaranteed monotonic here).
    size_t j_lo = i;
    bool have_lo = false;
    while (j_lo > 0) {
        --j_lo;
        if (Li - prof.L[j_lo] >= window) { have_lo = true; break; }
    }
    if (!have_lo) return -1;
    size_t j_hi = i;
    bool have_hi = false;
    while (j_hi + 1 < prof.L.size()) {
        ++j_hi;
        if (prof.L[j_hi] - Li >= window) { have_hi = true; break; }
    }
    if (!have_hi) return -1;
    const Vector d_in = prof.pts[i] - prof.pts[j_lo];    // direction of travel arriving at i
    const Vector d_out = prof.pts[j_hi] - prof.pts[i];   // direction of travel leaving i
    if (!(d_in.magnitude() > 0) || !(d_out.magnitude() > 0)) return -1;
    return d_in.angle(d_out) * 180.0 / M_PI;
}

StmMichelSplit WireCell::Clus::PR::stm_michel_stop_split(const StmMichelProfile& prof, int n_chain_segs,
                                                          const StmMichelSplitThresholds& th)
{
    StmMichelSplit out;
    if (th.max_split <= 0 || prof.empty() || n_chain_segs <= 0) return out;

    // SAME reference the retreat uses -- the two mechanisms must not drift
    // on what "collapsed" and "Bragg rise" mean.
    const double plateau = profile_plateau(prof, th.plateau_lo, th.plateau_hi, th.min_dqdx_live);
    if (!(plateau > 0)) return out;   // cannot judge

    const int last_seg = n_chain_segs - 1;
    double best_kink = -1;
    for (size_t i = 0; i < prof.L.size(); ++i) {
        // Only rows INSIDE the last chain segment: anywhere else a graph
        // vertex already exists and belongs to stm_michel_stop_retreat, not
        // this function -- that is the whole division of labor between them.
        if (prof.seg_idx[i] != last_seg) continue;
        if (prof.rr[i] < th.min_drop || prof.rr[i] > th.max_drop_len) continue;

        const double kink = stm_michel_row_kink_deg(prof, i, th.dir_window);
        if (kink < th.kink_min_deg) continue;   // also excludes the -1 "unmeasurable" case

        const double boundary_L = prof.L[i];

        std::vector<double> tail;
        for (size_t j = 0; j < prof.L.size(); ++j) {
            if (prof.L[j] < boundary_L) continue;
            if (prof.dQdx[j] < th.min_dqdx_live) continue;
            tail.push_back(prof.dQdx[j]);
        }
        if (static_cast<int>(tail.size()) < th.min_tail_pts) continue;
        const double tail_med = stm_michel_median(tail);
        // doc pdvd/82: see stm_michel_stop_retreat -- the doc 58 reading is the
        // first clause and exits where it always did.
        const bool collapsed_plateau = tail_med < th.collapse_frac * plateau;
        if (!collapsed_plateau && !(th.tail_peak_frac > 0)) continue;   // not a collapsed tail

        std::vector<double> kept_q; std::vector<double> kept_new_rr;
        for (size_t j = 0; j < prof.L.size(); ++j) {
            if (prof.L[j] > boundary_L) continue;
            if (prof.dQdx[j] < th.min_dqdx_live) continue;
            kept_q.push_back(prof.dQdx[j]);
            kept_new_rr.push_back(boundary_L - prof.L[j]);
        }
        if (kept_q.size() < 3) continue;
        const auto rm = running_median3(kept_q);
        double peak = -1;
        int n_win = 0;
        for (size_t j = 0; j < rm.size(); ++j) {
            if (kept_new_rr[j] > th.peak_window) continue;
            ++n_win;
            peak = std::max(peak, rm[j]);
        }
        if (n_win < 3) continue;   // window too sparse to judge
        if (!(peak >= th.peak_frac * plateau)) continue;   // nothing to retreat TO
        // doc pdvd/82: the peak-relative tail, when it is what admits this row.
        // `kink` is the row's own bend, already measured above; this asks it to
        // clear the higher tail_peak_kink_min bar as well.
        if (!collapsed_plateau && !(tail_med <= th.tail_peak_frac * peak && kink >= th.tail_peak_kink_min)) continue;

        if (kink > best_kink) {
            best_kink = kink;
            out.ok = true;
            out.index = i;
            out.cut_rr = prof.rr[i];
            out.drop_len = prof.total_length - boundary_L;
            out.kink_deg = kink;
            out.plateau = plateau;
            out.tail_med = tail_med;
            out.peak = peak;
            out.by_tail_peak = !collapsed_plateau;
        }
    }
    return out;
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
    const bool p2_on = th.michel_mip_lo_turned >= 0 || th.michel_far_len_shower_max >= 0 || th.michel_kink_window > 0;
    if (!p2_on) {
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
    // doc pdvd/73 (P2).  (b): measure_arm capped the far walk at michel_max_len,
    // so a shower-flagged arm's subtree is re-measured up to its own cap, with
    // the stop vertex fenced off (a loop back to the stop must not add the
    // muon chain).  (c): the turn over the shorter window.  Only the Michel
    // clause reads either; the continuation test above and the persisted
    // kink keep the michel/continuation window.
    if (th.michel_far_len_shower_max >= 0 && a.shower_like && !a.terminal) {
        VertexPtr far = find_other_vertex(g, arm, stop);
        a.far_len = stm_michel_far_subtree_len(g, far, arm, stop,
                                               std::max(th.michel_max_len, th.michel_far_len_shower_max));
    }
    if (th.michel_kink_window > 0 && last_muon)
        a.kink_w_deg = segment_pair_kink_deg(last_muon, arm, stm_michel_vertex_point(stop), th.michel_kink_window);
    if (stm_michel_michel_gate(a.len, a.far_len, a.mip, a.kink_deg, a.kink_w_deg, a.shower_like, th)) {
        a.kind = StmMichelArm::kMichel;
        return a;
    }
    a.kind = StmMichelArm::kOther;
    return a;
}

bool WireCell::Clus::PR::stm_michel_michel_gate(double len, double far_len, double mip, double kink_deg,
                                                double kink_w_deg, bool shower_like,
                                                const StmMichelArmThresholds& th)
{
    const bool kink_ok = kink_deg >= 0;
    const bool reach = (th.michel_far_len_shower_max >= 0 && shower_like)
        ? (len <= th.michel_max_len && far_len <= th.michel_far_len_shower_max)    // (b)
        : (len + far_len <= th.michel_max_len);
    // (a): a diluted PDVD Michel is admitted below michel_mip_lo by its
    // TOPOLOGY -- a hard turn -- never by low dQ/dx alone (the comment above
    // the continuation test keeps its meaning).
    const bool charge_lo = mip > th.michel_mip_lo ||
        (th.michel_mip_lo_turned >= 0 && kink_ok && kink_deg >= th.michel_mip_lo_turned_kink_deg &&
         mip > th.michel_mip_lo_turned);
    const bool charge = charge_lo && mip < th.michel_mip_hi;
    const bool shower_admits = shower_like &&
        (th.michel_shower_min_kink_deg < 0 || !kink_ok || kink_deg >= th.michel_shower_min_kink_deg);
    const bool turn = (kink_ok && kink_deg >= th.michel_min_kink_deg) ||
        (th.michel_kink_window > 0 && kink_w_deg >= 0 && kink_w_deg >= th.michel_min_kink_deg);   // (c)
    return reach && charge && (shower_admits || turn);
}

double WireCell::Clus::PR::stm_michel_far_subtree_len(Graph& g, VertexPtr far_vtx, SegmentPtr stem,
                                                      VertexPtr stop_vtx, double cap)
{
    double total = 0;
    if (!far_vtx || !far_vtx->descriptor_valid()) return total;
    std::set<SegmentPtr> used_segs{stem};      // membership only, never iterated
    std::set<VertexPtr> used_vtx{far_vtx};
    std::vector<VertexPtr> stack{far_vtx};
    while (!stack.empty()) {
        VertexPtr v = stack.back();
        stack.pop_back();
        if (!v || !v->descriptor_valid()) continue;
        for (auto edesc : sorted_out_edges(v->get_descriptor(), g)) {
            SegmentPtr sg = g[edesc].segment;
            if (!sg || used_segs.count(sg)) continue;
            used_segs.insert(sg);
            VertexPtr ov = find_other_vertex(g, sg, v);
            if (stop_vtx && ov == stop_vtx) continue;   // back into the stop: the muon's side
            total += segment_track_length(sg);
            if (total > cap) return total;
            if (ov && !used_vtx.count(ov)) {
                used_vtx.insert(ov);
                stack.push_back(ov);
            }
        }
    }
    return total;
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

WireCell::Clus::PR::StmMichelNearArms
WireCell::Clus::PR::stm_michel_near_stop_arms(Graph& g, const std::vector<SegmentPtr>& chain,
                                              const std::vector<VertexPtr>& chain_vtxs, double max_dist,
                                              const StmMichelArmThresholds& th,
                                              const std::function<bool(const SegmentPtr&)>& skip)
{
    StmMichelNearArms out;
    if (!(max_dist > 0) || chain.size() < 2 || chain_vtxs.size() != chain.size() + 1) return out;
    std::set<SegmentPtr> chain_set(chain.begin(), chain.end());   // membership only, never iterated
    double dist = 0;
    for (size_t vi = chain.size() - 1; vi >= 1; --vi) {
        dist += std::max(0.0, segment_track_length(chain[vi]));
        if (dist > max_dist) break;
        VertexPtr v = chain_vtxs[vi];
        SegmentPtr in_seg = chain[vi - 1];
        if (!v || !v->descriptor_valid()) continue;
        std::vector<StmMichelArm> found;
        for (auto e : sorted_out_edges(v->get_descriptor(), g)) {
            auto arm = g[e].segment;
            if (!arm || chain_set.count(arm)) continue;
            if (skip && skip(arm)) continue;
            const auto body = stm_michel_classify_chain_arm(g, in_seg, arm, v, th);
            if (body.kind == StmMichelArm::kDelta || body.kind == StmMichelArm::kHadron) continue;
            ++out.n_examined;
            auto a = stm_michel_classify_stop_arm(g, in_seg, arm, v, th);
            if (a.kind == StmMichelArm::kMichel) found.push_back(a);
        }
        if (!found.empty()) {
            std::sort(found.begin(), found.end(), [](const StmMichelArm& a, const StmMichelArm& b) {
                if (a.len != b.len) return a.len > b.len;
                return a.seg->get_graph_index() < b.seg->get_graph_index();
            });
            out.vtx_index = static_cast<int>(vi);
            out.dist = dist;
            out.michel = std::move(found);
            return out;
        }
    }
    return out;
}

bool WireCell::Clus::PR::stm_michel_stop_gamma_ring(double d_stop, double len,
                                                    double inner_cm, double outer_cm,
                                                    double max_len_cm)
{
    if (!std::isfinite(d_stop) || !std::isfinite(len)) return false;
    if (d_stop <= inner_cm) return false;   // the Michel's, not the gamma's
    if (d_stop > outer_cm) return false;
    if (len > max_len_cm) return false;     // a blob, not another cosmic
    return true;
}

bool WireCell::Clus::PR::stm_michel_stop_gamma_energy(double ke_mev, double lo_mev, double hi_mev)
{
    if (!std::isfinite(ke_mev)) return false;
    return ke_mev >= lo_mev && ke_mev <= hi_mev;
}

double WireCell::Clus::PR::stm_michel_admit_radius(double michel_cm,
                                                   double gamma_cm, bool gamma_on,
                                                   double survey_cm, bool survey_on)
{
    double r = std::isfinite(michel_cm) ? michel_cm : 0.0;
    if (gamma_on && std::isfinite(gamma_cm)) r = std::max(r, gamma_cm);
    if (survey_on && std::isfinite(survey_cm)) r = std::max(r, survey_cm);
    return r;
}

unsigned WireCell::Clus::PR::stm_michel_topology_clear(unsigned reject_bits, int michel_found, int conn_type,
                                                       double ke_mev, double len_cm,
                                                       double ke_min_mev, double len_min_cm, bool clears_sparse)
{
    if (!michel_found) return 0;
    if (conn_type != 1 && conn_type != 2) return 0;   // attached or bridged; a charge-only object (3) is not topology
    if (!std::isfinite(ke_mev) || !std::isfinite(len_cm)) return 0;
    if (ke_mev < ke_min_mev || len_cm < len_min_cm) return 0;
    unsigned clearable = R_NO_BRAGG | R_SHAPE_FLAT;
    if (clears_sparse) clearable |= R_PROFILE_SPARSE;
    return reject_bits & clearable;
}

int WireCell::Clus::PR::stm_michel_gamma_gate(double d_stop, double len, double cos_dir, double d_mich,
                                              double d_body, double ke_mev, double radius, double max_len,
                                              double cos_min, double max_ke_mev)
{
    for (double v : {d_stop, len, cos_dir, d_mich, d_body, ke_mev})
        if (!std::isfinite(v)) return 6;
    if (d_stop > radius) return 1;       // near the stop
    if (len > max_len) return 2;         // a dot, not a track
    if (cos_dir < cos_min) return 3;     // along the Michel electron
    if (d_body <= d_mich) return 4;      // the Michel's, not the muon body's
    if (ke_mev > max_ke_mev) return 5;   // a gamma's energy, not an over-clustered lump's
    return 0;
}

std::vector<int> WireCell::Clus::PR::stm_michel_gamma_take(double core_ke_mev, const std::vector<double>& ke_mev,
                                                           double total_max_mev)
{
    std::vector<int> take(ke_mev.size(), 0);
    if (!std::isfinite(core_ke_mev)) return take;
    double total = core_ke_mev;
    for (size_t i = 0; i < ke_mev.size(); ++i) {
        const double e = ke_mev[i];
        if (!std::isfinite(e) || e < 0) continue;
        if (total + e > total_max_mev) continue;   // skipped; a smaller, farther blob may still fit
        take[i] = 1;
        total += e;
    }
    return take;
}

WireCell::Clus::PR::StmMichelMovedStopSpare
WireCell::Clus::PR::stm_michel_moved_stop_spare(double kink_deg, double reach_cm,
                                                double kink_min_deg, double reach_min_cm)
{
    // The kink test is doc pdvd/72's expression verbatim, first, so with the
    // reach test off this returns exactly what the T2c site decided before
    // doc pdvd/84.  A NaN fails both comparisons and is vetoed.
    if (kink_min_deg >= 0 && kink_deg >= kink_min_deg) return StmMichelMovedStopSpare::kKink;
    if (reach_min_cm >= 0 && reach_cm >= reach_min_cm) return StmMichelMovedStopSpare::kReach;
    return StmMichelMovedStopSpare::kVeto;
}

bool WireCell::Clus::PR::stm_michel_stop_gamma_withhold(bool require_stm, unsigned reject_bits)
{
    return require_stm && reject_bits != 0;
}

std::vector<char> WireCell::Clus::PR::stm_michel_rows_keep(const std::vector<int>& roles, int drop_role)
{
    std::vector<char> keep(roles.size(), 1);
    for (size_t i = 0; i < roles.size(); ++i)
        if (roles[i] == drop_role) keep[i] = 0;
    return keep;
}
