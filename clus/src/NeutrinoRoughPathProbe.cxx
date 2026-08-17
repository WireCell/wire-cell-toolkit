// doc sbnd_xin/docs/pr/51 round 4 -- rough-path diagnostic probe.
//
// Owner Bee hand-scan of 131357/268067/285567/506746 rejected rounds 2-3
// (main_vertex_graph_audit): those fixed near-vertex GRAPH SHAPE; the
// reported symptom -- a short, low-dQ/dx "short-cut" that skips the true
// direct connection near the main vertex, sometimes with a visible jump --
// is a PATH-COST problem (see the m_rough_path_probe member block in
// NeutrinoPatternBase.h for the two candidate mechanisms, H1/H2).  This
// pass diagnoses which mechanism wrote the near-vertex segments and how
// large a charge-support penalty would change the Dijkstra answer.  It
// NEVER edits the graph, a segment, or a fit -- every line is
// SPDLOG_LOGGER_TRACE, gated entirely on m_rough_path_probe (C++ default
// false => immediate return => byte-identical).
//
// Toolkit-only; no prototype counterpart (this is instrumentation, not a
// ported algorithm).

#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/Graphs.h"

#include "WireCellAux/Logger.h"

#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <algorithm>
#include <chrono>
#include <map>
#include <sstream>

using namespace WireCell;
using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

namespace {

    double point_dis(const WireCell::Point& a, const WireCell::Point& b)
    {
        return (a - b).magnitude();
    }

    // Segments of `cluster` with at least one wcpt within `radius` of
    // `mv_pt` -- same scope idiom as main_vertex_graph_audit's
    // in_scope_segments, but read-only and with no `created` exemption
    // (the probe wants every near-vertex segment as it exists today).
    std::vector<SegmentPtr> near_vertex_segments(Graph& graph, const Facade::Cluster& cluster,
                                                  const WireCell::Point& mv_pt, double radius)
    {
        std::vector<SegmentPtr> out;
        for (const auto& ed : ordered_edges(graph)) {
            SegmentPtr sg = graph[ed].segment;
            if (!sg || sg->cluster() != &cluster) continue;
            bool near = false;
            for (const auto& wcp : sg->wcpts()) {
                if (point_dis(wcp.point, mv_pt) < radius) { near = true; break; }
            }
            if (near) out.push_back(sg);
        }
        return out;
    }

    // -- P1: path provenance -----------------------------------------
    //
    // There is no persistent per-segment origin tag in the data model (a
    // real one would need to be threaded through every wcpts()-writing
    // call site across NeutrinoStructureExaminer.cxx / NeutrinoOtherSegments.cxx /
    // NeutrinoGraphAudit.cxx / PRSegmentFunctions.cxx -- a much larger
    // change than a diagnostic probe warrants this round).  Instead this
    // reasons structurally from what is already on the segment: a
    // fresh do_rough_path result threads every intermediate point through
    // the Steiner point cloud (each wcpt has a real steiner_pc neighbor at
    // distance ~0), while examine_structure_1's straight-line replacement
    // literally IS the straight line between the two ends (each wcpt lies
    // near the chord).  "splice" is recognized structurally too (the mvga
    // op3 signature: an interior wcpt exactly matches the main vertex's
    // point).  This is a heuristic classifier, not ground truth -- stated
    // as such in every TRACE line.
    std::string classify_provenance(const SegmentPtr& sg, const WireCell::Point& mv_pt)
    {
        const auto& wcpts = sg->wcpts();
        if (wcpts.size() < 3) return "short";  // too few points to classify

        // "splice" signature: an interior point sits on the main vertex.
        for (size_t i = 1; i + 1 < wcpts.size(); ++i) {
            if (point_dis(wcpts[i].point, mv_pt) < 0.05 * units::cm) return "splice";
        }

        // Chord straightness: RMS perpendicular deviation of interior wcpts
        // from the straight line front->back, normalized by the chord
        // length.  A straight-line replacement (H2) reads ~0; a rough
        // Steiner path (H1) that actually follows charge around a turn
        // reads large.
        const WireCell::Point& p0 = wcpts.front().point;
        const WireCell::Point& p1 = wcpts.back().point;
        const double chord = point_dis(p0, p1);
        if (chord < 0.1 * units::cm) return "short";
        WireCell::Point dir = (p1 - p0) * (1.0 / chord);
        double sum_perp2 = 0;
        for (size_t i = 1; i + 1 < wcpts.size(); ++i) {
            WireCell::Point v = wcpts[i].point - p0;
            double along = v.dot(dir);
            WireCell::Point perp = v - dir * along;
            sum_perp2 += perp.dot(perp);
        }
        const double rms_perp = std::sqrt(sum_perp2 / (wcpts.size() - 2));
        const double straightness = rms_perp / chord;  // ~0 = straight line
        return (straightness < 0.02) ? "straighten(rough_or_straight)" : "rough(curved)";
    }

    // -- P2: per-segment support + charge profile ---------------------
    //
    // Classifies each consecutive wcpt-pair midpoint as live-supported,
    // dead-supported, or unsupported (Grouping::test_good_point at the
    // tight PR-level radius/ch_range NeutrinoStructureExaminer.cxx uses for
    // its own straight-line acceptance test, 0.2cm/0 -- NOT the looser
    // 0.6cm/1 graph-builder default), and reports the fitted dQ/dx ratio
    // (to m_mip_dqdx_median) for fit points falling in that span, flagging
    // whether those fit points were actually fitted or still carry
    // TrackFitting's default_dQ_dx seed (fit.valid() is the only signal
    // available at this level; the seed vs. measured distinction is a
    // known blind spot, called out explicitly in the doc).
    struct SupportTally {
        int n_live{0}, n_dead{0}, n_unsup{0};
    };

    // Classify one 3D point: 0=live, 1=dead, 2=unsupported (outside every
    // TPC counts as unsupported).
    int classify_point(const WireCell::Point& p, IDetectorVolumes::pointer dv,
                        WireCell::Clus::Facade::Grouping* grouping,
                        const WireCell::Clus::IPCTransform* transform, double cluster_t0)
    {
        auto wpid = dv->contained_by(p);
        if (wpid.face() == -1 || wpid.apa() == -1) return 2;
        auto p_raw = transform->backward(p, cluster_t0, wpid.face(), wpid.apa());
        auto scores = grouping->test_good_point(p_raw, wpid.apa(), wpid.face(), 0.2 * units::cm, 0);
        if (scores[0] || scores[1] || scores[2]) return 0;
        if (scores[3] || scores[4] || scores[5]) return 1;
        return 2;
    }

    void p2_support_profile(const SegmentPtr& sg, IDetectorVolumes::pointer dv,
                             WireCell::Clus::Facade::Grouping* grouping,
                             const WireCell::Clus::IPCTransform* transform, double cluster_t0,
                             double mip_dqdx_median)
    {
        const auto& wcpts = sg->wcpts();
        if (wcpts.size() < 2) return;
        SupportTally total;
        const double step = 0.3 * units::cm;
        for (size_t i = 0; i + 1 < wcpts.size(); ++i) {
            const WireCell::Point& a = wcpts[i].point;
            const WireCell::Point& b = wcpts[i + 1].point;
            const double d = point_dis(a, b);
            const int nsteps = std::max(1, static_cast<int>(std::round(d / step)));
            for (int k = 0; k <= nsteps; ++k) {
                const double f = static_cast<double>(k) / nsteps;
                WireCell::Point p(a.x() + (b.x() - a.x()) * f, a.y() + (b.y() - a.y()) * f,
                                   a.z() + (b.z() - a.z()) * f);
                const int cls = classify_point(p, dv, grouping, transform, cluster_t0);
                if (cls == 0) total.n_live++;
                else if (cls == 1) total.n_dead++;
                else total.n_unsup++;
            }
        }
        const int n_sampled = total.n_live + total.n_dead + total.n_unsup;
        int n_fitted = 0, n_default_seed = 0;
        double dqdx_sum = 0;
        for (const auto& fit : sg->fits()) {
            if (!fit.valid() || fit.dx <= 0) continue;
            const double ratio = (fit.dQ / fit.dx) / mip_dqdx_median;
            dqdx_sum += ratio;
            // TrackFitting::dQ_dx_fill seeds dQ = default_dQ_dx * dx before any
            // fit runs (TrackFitting.cxx:6068-6079); default_dQ_dx = 5000,
            // i.e. ratio ~= 5000/43000 ~= 0.116 (m_mip_dqdx_median = 43000/cm
            // by default).  A point sitting within 1% of that ratio is
            // flagged as "maybe still the seed, not a real fit" -- reported,
            // not filtered, because this is a heuristic, not a certainty.
            if (std::abs(ratio - 0.116) < 0.001) n_default_seed++;
            n_fitted++;
        }
        const double mean_ratio = n_fitted ? dqdx_sum / n_fitted : -1;
        SPDLOG_LOGGER_TRACE(s_log,
            "rough_path_probe P2: seg={} len={:.2f}cm nwcpt={} sampled={} live={} dead={} "
            "unsup={} unsup_frac={:.2f} nfits={} maybe_default_seed={} mean_dqdx_ratio={:.3f}",
            static_cast<const void*>(sg.get()), segment_track_length(sg) / units::cm, wcpts.size(),
            n_sampled, total.n_live, total.n_dead, total.n_unsup,
            n_sampled ? static_cast<double>(total.n_unsup) / n_sampled : -1.0, n_fitted,
            n_default_seed, mean_ratio);
    }

    // -- P3: counterfactual Dijkstra ladder ----------------------------
    //
    // Builds ONE re-weighted scratch copy of the cluster's "steiner_graph"
    // per scale (shared across every near-vertex segment, both for realism
    // -- Design A would build one penalized graph per cluster, not per
    // segment -- and for the cost measurement below), re-runs Dijkstra for
    // each flagged segment's endpoint pair, and reports whether the
    // resulting path changes and how much unsupported length it still
    // carries.  Read-only: the real "steiner_graph" is never touched.
    struct EdgeSupport { double bad_fraction{0}; bool bridge{false}; };

    void p3_ladder(Graph& /*graph*/, Facade::Cluster& cluster, VertexPtr main_vertex,
                   const std::vector<SegmentPtr>& segs, IDetectorVolumes::pointer dv,
                   WireCell::Clus::Facade::Grouping* grouping,
                   const WireCell::Clus::IPCTransform* transform, double cluster_t0)
    {
        if (!cluster.has_graph("steiner_graph") || !cluster.has_pc("steiner_pc")) {
            SPDLOG_LOGGER_TRACE(s_log, "rough_path_probe P3: no steiner_graph/steiner_pc on cluster {}, skipping",
                                cluster.ident());
            return;
        }
        const auto& base = cluster.get_graph("steiner_graph");
        const auto& steiner_pc = cluster.get_pc("steiner_pc");
        const auto& coords = cluster.get_default_scope().coords;
        const auto& xs = steiner_pc.get(coords.at(0))->elements<double>();
        const auto& ys = steiner_pc.get(coords.at(1))->elements<double>();
        const auto& zs = steiner_pc.get(coords.at(2))->elements<double>();
        auto pt_at = [&](size_t i) { return WireCell::Point(xs[i], ys[i], zs[i]); };

        const size_t n_edges = boost::num_edges(base);
        const double min_len = 0.5 * units::cm;

        // Per-edge support (computed once at scale-independent cost, reused
        // across the whole ladder): sample every min_len-scanned edge at
        // 0.3cm, classify each sample, and store the unsupported+dead
        // fraction (dead counts at alpha=0.25, per the H2 finding that
        // is_good_point/test_good_point treats dead as traversable but the
        // penalty must not be zero there either).
        const auto edge_weight_map = boost::get(boost::edge_weight, base);
        auto t_scan_start = std::chrono::steady_clock::now();
        std::vector<EdgeSupport> support(n_edges);
        std::vector<double> geom_len(n_edges);
        // (min-vertex, max-vertex) -> edge array index, for per-hop lookup
        // when walking a reconstructed path below (the graph is built with
        // vecS/vecS storage and no edge_index property, so boost::edges()
        // order is the only stable index this graph offers).
        std::map<std::pair<size_t, size_t>, size_t> edge_index_of;
        size_t edges_scanned = 0;
        {
            size_t idx = 0;
            for (auto [ei, ei_end] = boost::edges(base); ei != ei_end; ++ei, ++idx) {
                const double w = boost::get(edge_weight_map, *ei);
                geom_len[idx] = w;
                const auto sidx = boost::source(*ei, base);
                const auto tidx = boost::target(*ei, base);
                edge_index_of[{std::min(sidx, tidx), std::max(sidx, tidx)}] = idx;
                if (w < min_len) continue;
                edges_scanned++;
                WireCell::Point a = pt_at(sidx), b = pt_at(tidx);
                const double step = 0.3 * units::cm;
                const int nsteps = std::max(1, static_cast<int>(std::round(w / step)));
                int n_live = 0, n_dead = 0, n_unsup = 0;
                for (int k = 0; k <= nsteps; ++k) {
                    const double f = static_cast<double>(k) / nsteps;
                    WireCell::Point p(a.x() + (b.x() - a.x()) * f, a.y() + (b.y() - a.y()) * f,
                                       a.z() + (b.z() - a.z()) * f);
                    const int cls = classify_point(p, dv, grouping, transform, cluster_t0);
                    if (cls == 0) n_live++;
                    else if (cls == 1) n_dead++;
                    else n_unsup++;
                }
                const int n = n_live + n_dead + n_unsup;
                support[idx].bad_fraction = n ? (n_unsup + 0.25 * n_dead) / n : 0.0;
                // "bridge" here just means "long enough to matter"; the
                // finer bridge-vs-close-connect-vs-intra-blob distinction
                // needs the pre-Steiner-reduction graph and is reported as
                // a length/support pair instead (sufficient to answer
                // "does a gap-spanning edge exist on the winning path").
                support[idx].bridge = true;
            }
        }
        const double scan_ms =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t_scan_start).count();
        SPDLOG_LOGGER_TRACE(s_log,
            "rough_path_probe P3: cluster={} steiner_graph edges={} edges_scanned(>{:.1f}cm)={} scan_ms={:.1f}",
            cluster.ident(), n_edges, min_len / units::cm, edges_scanned, scan_ms);

        static const std::vector<double> kScales = {0, 1, 2, 3, 5, 10};

        for (const auto& sg : segs) {
            const auto& wcpts = sg->wcpts();
            if (wcpts.size() < 2) continue;
            WireCell::Point p_front = wcpts.front().point;
            WireCell::Point p_back = wcpts.back().point;
            auto first_knn = cluster.kd_steiner_knn(1, p_front, "steiner_pc");
            auto last_knn = cluster.kd_steiner_knn(1, p_back, "steiner_pc");
            if (first_knn.empty() || last_knn.empty()) continue;
            const size_t src = first_knn[0].first;
            const size_t dst = last_knn[0].first;

            // Today's path length + directly-sampled unsupported length, for
            // reference (sampled independently of the steiner_graph edge
            // set, since a straightened/spliced segment's wcpts need not
            // land on steiner_graph edges at all).
            double today_len = 0, today_unsup_len = 0;
            for (size_t i = 0; i + 1 < wcpts.size(); ++i) {
                const WireCell::Point& a = wcpts[i].point;
                const WireCell::Point& b = wcpts[i + 1].point;
                const double d = point_dis(a, b);
                today_len += d;
                const double step = 0.3 * units::cm;
                const int nsteps = std::max(1, static_cast<int>(std::round(d / step)));
                int n_unsup = 0;
                for (int k = 0; k <= nsteps; ++k) {
                    const double f = static_cast<double>(k) / nsteps;
                    WireCell::Point p(a.x() + (b.x() - a.x()) * f, a.y() + (b.y() - a.y()) * f,
                                       a.z() + (b.z() - a.z()) * f);
                    if (classify_point(p, dv, grouping, transform, cluster_t0) == 2) n_unsup++;
                }
                today_unsup_len += d * static_cast<double>(n_unsup) / (nsteps + 1);
            }

            std::ostringstream ladder;
            std::vector<size_t> baseline_path;
            for (double scale : kScales) {
                Graphs::Weighted::Graph scratch = base;
                if (scale > 0) {
                    auto sw = boost::get(boost::edge_weight, scratch);
                    size_t idx = 0;
                    for (auto [ei, ei_end] = boost::edges(scratch); ei != ei_end; ++ei, ++idx) {
                        if (!support[idx].bridge) continue;
                        const double w0 = geom_len[idx];
                        boost::put(sw, *ei, w0 * (1.0 + scale * support[idx].bad_fraction));
                    }
                }
                std::vector<Graphs::Weighted::vertex_type> preds(boost::num_vertices(scratch));
                std::vector<Graphs::Weighted::dijkstra_distance_type> dists(boost::num_vertices(scratch));
                boost::dijkstra_shortest_paths(
                    scratch, src,
                    boost::predecessor_map(&preds[0]).distance_map(&dists[0]));
                // Walk predecessors src<-...<-dst to reconstruct the path.
                std::vector<size_t> path;
                size_t cur = dst;
                bool ok = true;
                for (int guard = 0; guard < static_cast<int>(boost::num_vertices(scratch)) + 1; ++guard) {
                    path.push_back(cur);
                    if (cur == src) break;
                    if (preds[cur] == cur) { ok = false; break; }  // unreachable
                    cur = preds[cur];
                }
                if (!ok || path.empty() || path.back() != src) {
                    ladder << "scale=" << scale << ":unreachable ";
                    continue;
                }
                std::reverse(path.begin(), path.end());
                if (scale == 0) baseline_path = path;
                double plen = 0, unsup_len = 0;
                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    const double d = point_dis(pt_at(path[i]), pt_at(path[i + 1]));
                    plen += d;
                    // Reuse the per-edge bad_fraction if this hop is one of
                    // the scanned base-graph edges; else assume supported
                    // (intra-blob / sub-min_len hop -- on charge by
                    // construction, per SteinerGrapher.cxx).
                    auto key = std::make_pair(std::min(path[i], path[i + 1]), std::max(path[i], path[i + 1]));
                    auto it = edge_index_of.find(key);
                    if (it != edge_index_of.end()) {
                        unsup_len += d * support[it->second].bad_fraction;
                    }
                }
                const bool changed = (scale > 0) && (path != baseline_path);
                ladder << "scale=" << scale << ":len=" << plen / units::cm << "cm"
                       << ":unsup=" << unsup_len / units::cm << "cm"
                       << (changed ? ":CHANGED" : ":same") << " ";
            }
            SPDLOG_LOGGER_TRACE(s_log,
                "rough_path_probe P3: seg={} today_len={:.2f}cm today_unsup_len={:.2f}cm src={} dst={} ladder: {}",
                static_cast<const void*>(sg.get()), today_len / units::cm, today_unsup_len / units::cm, src, dst,
                ladder.str());
        }
    }

}  // namespace

void PatternAlgorithms::rough_path_probe(Graph& graph, Facade::Cluster& cluster, VertexPtr main_vertex,
                                          TrackFitting& track_fitter, IDetectorVolumes::pointer dv)
{
    if (!m_rough_path_probe) return;
    if (!main_vertex || !main_vertex->descriptor_valid()) return;

    const WireCell::Point mv_pt =
        main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;

    auto segs = near_vertex_segments(graph, cluster, mv_pt, m_mvga_radius);
    SPDLOG_LOGGER_TRACE(s_log, "rough_path_probe: cluster={} main_vertex=({:.2f},{:.2f},{:.2f}) "
                               "n_near_vertex_segments={} radius={:.1f}cm",
                        cluster.ident(), mv_pt.x() / units::cm, mv_pt.y() / units::cm,
                        mv_pt.z() / units::cm, segs.size(), m_mvga_radius / units::cm);

    // P1: provenance (heuristic classifier -- see classify_provenance doc).
    for (const auto& sg : segs) {
        SPDLOG_LOGGER_TRACE(s_log, "rough_path_probe P1: seg={} len={:.2f}cm nwcpt={} prov={}",
                            static_cast<const void*>(sg.get()), segment_track_length(sg) / units::cm,
                            sg->wcpts().size(), classify_provenance(sg, mv_pt));
    }

    auto transform = track_fitter.get_pc_transforms()->pc_transform(cluster.get_scope_transform(cluster.get_default_scope()));
    const double cluster_t0 = cluster.get_cluster_t0();
    auto grouping = cluster.grouping();
    if (!transform || !grouping) {
        SPDLOG_LOGGER_TRACE(s_log, "rough_path_probe: no transform/grouping on cluster {}, skipping P2/P3",
                            cluster.ident());
        return;
    }

    // P2: per-segment support + charge profile.
    for (const auto& sg : segs) {
        p2_support_profile(sg, dv, grouping, transform.get(), cluster_t0, m_mip_dqdx_median);
    }

    // P3: counterfactual Dijkstra ladder (shared scratch graph per scale).
    p3_ladder(graph, cluster, main_vertex, segs, dv, grouping, transform.get(), cluster_t0);
}
