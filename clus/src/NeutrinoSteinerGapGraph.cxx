// doc sbnd_xin/docs/pr/51 round 5 -- steiner_gap_penalty.
//
// The round-4 probe (NeutrinoRoughPathProbe.cxx) measured the owner's
// near-vertex "short-cut" class as a PATH-COST problem: do_rough_path runs
// plain Dijkstra on "steiner_graph", whose edge weight prices charge only
// at the two edge endpoints with a hard [0.8, 1.2] dynamic range
// (SteinerGrapher.cxx create_enhanced_steiner_graph), so a gap-spanning
// chord beats following charge around any turn sharper than ~150 deg.  The
// probe's P3 counterfactual ladder validated the fix shipped here: derive a
// per-cluster "steiner_graph_gap" flavor, identical topology and vertex
// indexing, every edge longer than m_sgp_min_edge re-weighted
//     w' = w * (1 + m_steiner_gap_penalty * bad_fraction)
// with bad_fraction the sampled unsupported(+alpha*dead) fraction of the
// edge interior, and let do_rough_path route on it.  Built lazily, ONCE per
// cluster, only when a rough-path query actually happens (owner runtime
// requirement: the round-4 scan cost is 5-70 ms per cluster, and only
// PR-routed clusters ever pay it).  Nothing but do_rough_path selects the
// new flavor: TaggerCheckSTM's verbatim fork, TrackFitting's reader and
// every other "steiner_graph" consumer keep the base graph.
//
// m_steiner_gap_penalty <= 0 (C++ default 0) => first-line return => the
// flavor is never built => byte-identical production path.
//
// doc sbnd_xin/docs/pr/51 round 6 -- sgp_weak_scale / sgp_weak_qref.
// The unsupported-fraction penalty is blind to chords whose interior is
// image-SUPPORTED but charge-poor: 18259-131357's trunk chord (point q
// ~674-1276 vs ~1900-5400 on the true corner route) and 18255-506746's
// hairpin connector (first points q = 18/238/633) have P3 ladders that are
// literally scale-invariant 0..10.  When m_sgp_weak_scale > 0 each scanned
// edge additionally pays a thresholded charge-deficit term:
//     w' = w * (1 + gap_scale*bad + weak_scale*deficit)
//     deficit = 0.5*(max(0,1-q_s/qref) + max(0,1-q_t/qref))
// with endpoint charges recovered once per steiner vertex via
// calc_charge_wcp(idx, 4000, false) -- the same call/flags the production
// steiner edge weighting used (CreateSteinerGraph.cxx:262 passes
// disable_dead_mix_cell=false, forwarded by the pr/29 D2 knob, SBND ON).
// m_sgp_weak_scale <= 0 (C++ default 0) => the round-5 reweight statements
// run verbatim => the gap flavor is byte-identical to round-5 production.
//
// Toolkit-only; no prototype counterpart (the prototype has the same
// endpoint-only charge weighting; this is a deliberate, knob-gated
// divergence, not a port fix).

#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/Graphs.h"

#include "WireCellAux/Logger.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

namespace {

    // Classify one 3D point: 0=live, 1=dead, 2=unsupported (outside every
    // TPC counts as unsupported).  Verbatim from NeutrinoRoughPathProbe.cxx
    // (deliberately duplicated, not refactored: the probe is a shipped
    // diagnostic whose behavior must not move underneath its round-4
    // measurements).
    int classify_point(const WireCell::Point& p, IDetectorVolumes::pointer dv,
                       const WireCell::Clus::Facade::Grouping* grouping,
                       const WireCell::Clus::IPCTransform* transform, double cluster_t0,
                       double point_radius)
    {
        auto wpid = dv->contained_by(p);
        if (wpid.face() == -1 || wpid.apa() == -1) return 2;
        auto p_raw = transform->backward(p, cluster_t0, wpid.face(), wpid.apa());
        auto scores = grouping->test_good_point(p_raw, wpid.apa(), wpid.face(), point_radius, 0);
        if (scores[0] || scores[1] || scores[2]) return 0;
        if (scores[3] || scores[4] || scores[5]) return 1;
        return 2;
    }

}  // namespace

namespace WireCell::Clus::PR {

    double gap_edge_bad_fraction(const WireCell::Point& a, const WireCell::Point& b,
                                 double length, double step, double dead_alpha,
                                 const std::function<int(const WireCell::Point&)>& classify)
    {
        const int nsteps = std::max(1, static_cast<int>(std::round(length / step)));
        int n_live = 0, n_dead = 0, n_unsup = 0;
        for (int k = 0; k <= nsteps; ++k) {
            const double f = static_cast<double>(k) / nsteps;
            WireCell::Point p(a.x() + (b.x() - a.x()) * f, a.y() + (b.y() - a.y()) * f,
                              a.z() + (b.z() - a.z()) * f);
            const int cls = classify(p);
            if (cls == 0) n_live++;
            else if (cls == 1) n_dead++;
            else n_unsup++;
        }
        const int n = n_live + n_dead + n_unsup;
        return n ? (n_unsup + dead_alpha * n_dead) / n : 0.0;
    }

    double weak_charge_deficit(double qa, double qb, double qref)
    {
        if (qref <= 0) return 0.0;
        const double da = std::max(0.0, 1.0 - qa / qref);
        const double db = std::max(0.0, 1.0 - qb / qref);
        return 0.5 * (da + db);
    }

}  // namespace WireCell::Clus::PR

bool PatternAlgorithms::ensure_steiner_gap_graph(const Facade::Cluster& cluster)
{
    if (m_steiner_gap_penalty <= 0) return false;
    if (cluster.has_graph("steiner_graph_gap")) return true;

    // Prerequisites.  Any missing => one DEBUG line, fall back to the base
    // flavor (do_rough_path stays on "steiner_graph" for this cluster; the
    // next query retries, which is cheap -- the checks below are O(1)).
    if (!cluster.has_graph("steiner_graph") || !cluster.has_pc("steiner_pc") ||
        cluster.get_pc("steiner_pc").size_major() == 0) {
        SPDLOG_LOGGER_DEBUG(s_log, "sgp build: cluster {} missing steiner_graph/steiner_pc, using base flavor",
                            cluster.ident());
        return false;
    }
    const auto* grouping = cluster.grouping();
    auto transform = m_sgp_pcts ? m_sgp_pcts->pc_transform(cluster.get_scope_transform(cluster.get_default_scope()))
                                : nullptr;
    if (!m_sgp_dv || !grouping || !transform) {
        SPDLOG_LOGGER_DEBUG(s_log, "sgp build: cluster {} missing dv/grouping/transform, using base flavor",
                            cluster.ident());
        return false;
    }
    const double cluster_t0 = cluster.get_cluster_t0();

    const auto& base = cluster.get_graph("steiner_graph");
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& xs = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& ys = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& zs = steiner_pc.get(coords.at(2))->elements<double>();

    auto classify = [&](const WireCell::Point& p) {
        return classify_point(p, m_sgp_dv, grouping, transform.get(), cluster_t0, m_sgp_point_radius);
    };

    // doc pr/51 round 6: lazy per-vertex charge recovery for the weak-charge
    // deficit term.  steiner_pc rows are exact copies of default-pc points
    // (Dataset::subset), so an exact-match 1-NN into the default kd-tree
    // recovers the original point index; calc_charge_wcp(idx, 4000, false)
    // then mirrors the charges the production steiner edge weighting used
    // (CreateSteinerGraph.cxx:262 + pr/29 D2).  Computed at most once per
    // steiner vertex, only for endpoints of scanned edges, only when the
    // weak term is on.
    const bool weak_on = (m_sgp_weak_scale > 0);
    std::vector<double> qcache;
    if (weak_on) qcache.assign(steiner_pc.size_major(), -1.0);
    auto q_of = [&](size_t i) -> double {
        double& q = qcache[i];
        if (q < 0) {
            const WireCell::Point p(xs[i], ys[i], zs[i]);
            const size_t didx = cluster.get_closest_point_index(p);
            q = cluster.calc_charge_wcp(didx, 4000.0, false).second;
        }
        return q;
    };

    const auto t_start = std::chrono::steady_clock::now();

    // Copy-construct: preserves vertex indices, so kd_steiner_knn indices
    // into "steiner_pc" remain valid on the gap flavor.  boost::edges() on
    // vecS storage iterates in insertion order -- deterministic, and no
    // pointer-keyed container anywhere in this scan.
    Graphs::Weighted::Graph gap = base;
    const auto weight_map = boost::get(boost::edge_weight, gap);
    size_t edges_scanned = 0, edges_penalized = 0, edges_weak = 0;
    for (auto [ei, ei_end] = boost::edges(gap); ei != ei_end; ++ei) {
        const double w = boost::get(weight_map, *ei);
        if (w < m_sgp_min_edge) continue;
        edges_scanned++;
        const auto sidx = boost::source(*ei, gap);
        const auto tidx = boost::target(*ei, gap);
        const WireCell::Point a(xs[sidx], ys[sidx], zs[sidx]);
        const WireCell::Point b(xs[tidx], ys[tidx], zs[tidx]);
        // `w` is the charge-weighted length (geometric x [0.8, 1.2]) -- the
        // probe's P3 ladder used it for both the min-edge gate and the
        // sample count, and the shipped fix must match what was validated.
        const double bad = gap_edge_bad_fraction(a, b, w, m_sgp_sample_step, m_sgp_dead_alpha, classify);
        if (!weak_on) {
            // Round-5 statements, verbatim: weak off => byte-identical flavor.
            if (bad <= 0) continue;
            edges_penalized++;
            boost::put(weight_map, *ei, w * (1.0 + m_steiner_gap_penalty * bad));
            continue;
        }
        // doc pr/51 round 6: combined form.
        const double deficit = weak_charge_deficit(q_of(sidx), q_of(tidx), m_sgp_weak_qref);
        if (deficit > 0) edges_weak++;
        if (bad <= 0 && deficit <= 0) continue;
        edges_penalized++;
        boost::put(weight_map, *ei,
                   w * (1.0 + m_steiner_gap_penalty * bad + m_sgp_weak_scale * deficit));
    }

    const double scan_ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_start).count();

    // The gap flavor is a deterministic derived cache of "steiner_graph",
    // which has no writer during the tagger stage (producers run in
    // CreateSteinerGraph, before PR) -- the same justification as the
    // const_cast inside Facade::Cluster::graph_algorithms() const
    // ("we are caching, so const cast is okay", Facade_Cluster.cxx).  If a
    // future stage ever rewrites "steiner_graph" mid-PR it must invalidate
    // "steiner_graph_gap" alongside.
    const_cast<Facade::Cluster&>(cluster).give_graph("steiner_graph_gap", std::move(gap));

    if (!weak_on) {
        SPDLOG_LOGGER_DEBUG(s_log,
            "sgp build: cluster {} edges={} scanned={} penalized={} scale={:.1f} scan_ms={:.1f}",
            cluster.ident(), boost::num_edges(base), edges_scanned, edges_penalized,
            m_steiner_gap_penalty, scan_ms);
        return true;
    }

    // doc pr/51 round 6: extended sentinel.  Charge quartiles of the
    // recovered vertex charges tell the qref operating point where the
    // event's charge scale actually sits.
    std::vector<double> qs_seen;
    qs_seen.reserve(qcache.size());
    for (double q : qcache)
        if (q >= 0) qs_seen.push_back(q);
    std::sort(qs_seen.begin(), qs_seen.end());
    const auto quart = [&](double f) -> double {
        return qs_seen.empty() ? 0.0 : qs_seen[static_cast<size_t>(f * (qs_seen.size() - 1))];
    };
    SPDLOG_LOGGER_DEBUG(s_log,
        "sgp build: cluster {} edges={} scanned={} penalized={} scale={:.1f} scan_ms={:.1f}"
        " weak_scale={:.1f} qref={:.0f} weak_edges={} nq={} q25={:.0f} q50={:.0f} q75={:.0f}",
        cluster.ident(), boost::num_edges(base), edges_scanned, edges_penalized,
        m_steiner_gap_penalty, scan_ms, m_sgp_weak_scale, m_sgp_weak_qref, edges_weak,
        qs_seen.size(), quart(0.25), quart(0.50), quart(0.75));
    return true;
}
