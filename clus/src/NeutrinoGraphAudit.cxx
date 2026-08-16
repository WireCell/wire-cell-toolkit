// doc sbnd_xin/docs/pr/51 -- main-vertex graph audit.
//
// Near-vertex PR graph-shape repair on the main cluster, gated on
// m_main_vertex_graph_audit (C++ default false => immediate return =>
// byte-identical).  See the m_main_vertex_graph_audit member block in
// NeutrinoPatternBase.h for the full design, the measured per-event
// pathologies (131357 / 268067 / 360535 / 142421 / 285567) and the guard
// list.  Toolkit-only pass; prototype ancestry (not a port target):
// NeutrinoID_final_structure.h examine_structure_final_1 (lines 545-696)
// and _1p (lines 402-544), NeutrinoID_improve_vertex.h
// eliminate_short_vertex_activities (lines 365-1018).

#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"

#include "WireCellAux/Logger.h"

#include <algorithm>
#include <deque>

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

namespace {

    // Points of a segment for geometry tests: valid fit points when the
    // segment has been fitted, wcpts otherwise (a rough-path segment created
    // mid-pass has no fits until the op4 refit).
    std::vector<WireCell::Point> seg_points(const SegmentPtr& sg)
    {
        std::vector<WireCell::Point> pts;
        for (const auto& fit : sg->fits()) {
            if (fit.valid()) pts.push_back(fit.point);
        }
        if (pts.empty()) {
            for (const auto& wcp : sg->wcpts()) pts.push_back(wcp.point);
        }
        return pts;
    }

    size_t seg_valid_fits(const SegmentPtr& sg)
    {
        size_t n = 0;
        for (const auto& fit : sg->fits()) {
            if (fit.valid()) ++n;
        }
        return n;
    }

    double point_dis(const WireCell::Point& a, const WireCell::Point& b)
    {
        return (a - b).magnitude();
    }

    // All graph vertices of `cluster` reachable from `from` without
    // traversing `exclude`.  Membership only -- iteration order does not
    // leak into any output.
    std::set<Graph::vertex_descriptor> reachable_without(
        const Graph& graph, const Facade::Cluster& cluster,
        const VertexPtr& from, const SegmentPtr& exclude)
    {
        std::set<Graph::vertex_descriptor> seen;
        if (!from || !from->descriptor_valid()) return seen;
        std::deque<Graph::vertex_descriptor> todo{from->get_descriptor()};
        seen.insert(from->get_descriptor());
        while (!todo.empty()) {
            auto vd = todo.front();
            todo.pop_front();
            for (auto [ei, ei_end] = boost::out_edges(vd, graph); ei != ei_end; ++ei) {
                if (graph[*ei].segment == exclude) continue;
                auto nd = boost::target(*ei, graph);
                VertexPtr nv = graph[nd].vertex;
                if (!nv || nv->cluster() != &cluster) continue;
                if (seen.insert(nd).second) todo.push_back(nd);
            }
        }
        return seen;
    }

}  // namespace

bool PatternAlgorithms::main_vertex_graph_audit(Graph& graph, Facade::Cluster& cluster,
                                                VertexPtr main_vertex,
                                                TrackFitting& track_fitter,
                                                IDetectorVolumes::pointer dv)
{
    if (!m_main_vertex_graph_audit) return false;
    if (!main_vertex || !main_vertex->descriptor_valid()) return false;

    WireCell::Point mv_pt =
        main_vertex->fit().valid() ? main_vertex->fit().point : main_vertex->wcpt().point;

    // Segments created by this pass's own reconnects: exempt from every op
    // (no delete/recreate cycling; they carry no fits until the op4 refit).
    IndexedSegmentSet created;

    int n_op1 = 0, n_op2 = 0, n_op3 = 0;
    constexpr int kEditCap = 8;  // per-op paranoia cap; every real case needs <= 3

    // In-scope segments of the main cluster, in stable graph-index order.
    auto in_scope_segments = [&]() {
        std::vector<SegmentPtr> segs;
        for (const auto& vd : ordered_nodes(graph)) {
            VertexPtr vtx = graph[vd].vertex;
            if (!vtx || vtx->cluster() != &cluster) continue;
            for (auto edesc : sorted_out_edges(vd, graph)) {
                SegmentPtr sg = graph[edesc].segment;
                if (!sg) continue;
                if (std::find(segs.begin(), segs.end(), sg) != segs.end()) continue;
                if (created.count(sg)) continue;
                bool near = false;
                for (const auto& p : seg_points(sg)) {
                    if (point_dis(p, mv_pt) < m_mvga_radius) { near = true; break; }
                }
                if (near) segs.push_back(sg);
            }
        }
        std::sort(segs.begin(), segs.end(), SegmentIndexCmp{});
        return segs;
    };

    // Nearest in-cluster vertex to `p` among `allowed` descriptors,
    // excluding `not_this`.  ordered_nodes + strict < keeps the argmin
    // deterministic (first minimum in node-index order wins).
    auto nearest_vertex_in = [&](const WireCell::Point& p,
                                 const std::set<Graph::vertex_descriptor>& allowed,
                                 const VertexPtr& not_this) {
        double min_dis = 1e9;
        VertexPtr min_vtx = nullptr;
        for (const auto& vd : ordered_nodes(graph)) {
            VertexPtr vtx = graph[vd].vertex;
            if (!vtx || vtx->cluster() != &cluster) continue;
            if (vtx == not_this) continue;
            if (!allowed.empty() && !allowed.count(vd)) continue;
            WireCell::Point vp = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
            double dis = point_dis(p, vp);
            if (dis < min_dis) {
                min_dis = dis;
                min_vtx = vtx;
            }
        }
        return std::make_pair(min_vtx, min_dis);
    };

    // Direct rough-path edge v1 -> v2 unless one already exists.  Returns
    // the (possibly pre-existing) connecting segment, or nullptr when no
    // usable path exists.  New segments enter `created`.
    auto connect_direct = [&](VertexPtr v1, VertexPtr v2) -> SegmentPtr {
        if (!v1 || !v2 || v1 == v2) return nullptr;
        if (SegmentPtr ex = find_segment(graph, v1, v2)) return ex;
        Facade::geo_point_t p1 = v1->wcpt().point;
        Facade::geo_point_t p2 = v2->wcpt().point;
        auto path_points = do_rough_path(cluster, p1, p2);
        if (path_points.size() < 2) return nullptr;
        auto sg = create_segment_for_cluster(cluster, dv, path_points, 0);
        if (!sg) return nullptr;
        add_segment(graph, sg, v1, v2);
        created.insert(sg);
        return sg;
    };

    // Drop a now-degree-0 vertex (never the main vertex, never a protected
    // break -- pr/48 / pr/50 precedent).
    auto cleanup_vertex = [&](VertexPtr vtx) {
        if (!vtx || vtx == main_vertex) return;
        if (vtx->flags_any(VertexFlags::kProtectedBreak)) return;
        if (!vtx->descriptor_valid()) return;
        if (boost::degree(vtx->get_descriptor(), graph) == 0) {
            remove_vertex(graph, vtx);
        }
    };

    // ---- op1: duplicate-corridor merge -----------------------------------
    if (m_mvga_dup_frac > 0 && m_mvga_dup_tol > 0) {
        bool flag_continue = true;
        while (flag_continue && n_op1 < kEditCap) {
            flag_continue = false;
            auto segs = in_scope_segments();
            for (size_t i = 0; i + 1 < segs.size() && !flag_continue; ++i) {
                for (size_t j = i + 1; j < segs.size() && !flag_continue; ++j) {
                    SegmentPtr sa = segs[i];
                    SegmentPtr sb = segs[j];
                    double la = segment_track_length(sa);
                    double lb = segment_track_length(sb);
                    SegmentPtr shorter = (la <= lb) ? sa : sb;
                    SegmentPtr longer  = (la <= lb) ? sb : sa;
                    auto pts_s = seg_points(shorter);
                    auto pts_l = seg_points(longer);
                    // 3-4 point paths make the fraction degenerate -- op3's
                    // point-degeneracy gate owns those.
                    if (pts_s.size() <= static_cast<size_t>(m_mvga_stub_pts)) continue;
                    double frac = path_overlap_fraction(pts_s, pts_l, m_mvga_dup_tol);
                    // Probe line for threshold calibration (TRACE, cf. op2).
                    if (frac > 0.25) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op1 eval cluster={} pair len {:.2f}/{:.2f}cm npts {}/{} overlap={:.2f}",
                            cluster.ident(), segment_track_length(shorter)/units::cm,
                            segment_track_length(longer)/units::cm,
                            pts_s.size(), pts_l.size(), frac);
                    }
                    if (frac < m_mvga_dup_frac) continue;
                    // Near-parallel guard: a corridor duplicate runs
                    // (anti)parallel to its ribbon (268067 rider 13 deg,
                    // 360535 pair ~parallel, 285567 shorts 7-11 deg); a
                    // genuine small-opening V whose short prong hugs the
                    // long one within tol must NOT merge.  Chord-vs-chord,
                    // folded to [0,90] deg.  <= 0 disables.
                    if (m_mvga_dup_angle > 0 && pts_s.size() >= 2 && pts_l.size() >= 2) {
                        auto chord = [](const std::vector<WireCell::Point>& pts) {
                            return pts.back() - pts.front();
                        };
                        auto ca = chord(pts_s);
                        auto cb = chord(pts_l);
                        double den = ca.magnitude() * cb.magnitude();
                        if (den > 0) {
                            double cosang = std::abs(ca.dot(cb)) / den;
                            double ang = std::acos(std::clamp(cosang, 0.0, 1.0)) / 3.1415926 * 180.0;
                            if (ang > m_mvga_dup_angle) continue;
                        }
                    }

                    // Keep the higher-integrated-charge member (tie: keep longer).
                    double qa = segment_integrated_dQ(sa);
                    double qb = segment_integrated_dQ(sb);
                    SegmentPtr loser;
                    if (qa == qb) loser = shorter;
                    else loser = (qa < qb) ? sa : sb;
                    SegmentPtr survivor = (loser == sa) ? sb : sa;

                    auto [lv1, lv2] = find_vertices(graph, loser);
                    auto [sv1, sv2] = find_vertices(graph, survivor);
                    if (!lv1 || !lv2 || !sv1 || !sv2) continue;

                    // Reconnect plan: each loser endpoint that is not a
                    // survivor endpoint gets a direct edge to the nearest
                    // survivor endpoint (doc pr/51 op1; 268067: V_A -> MAIN,
                    // 360535: MAIN -> V3).  Pre-verify the rough paths so a
                    // failed reconnect never strands the graph mid-edit.
                    std::vector<std::pair<VertexPtr, VertexPtr>> plans;
                    bool feasible = true;
                    for (VertexPtr le : {lv1, lv2}) {
                        if (le == sv1 || le == sv2) continue;
                        WireCell::Point lp = le->fit().valid() ? le->fit().point : le->wcpt().point;
                        WireCell::Point p1 = sv1->fit().valid() ? sv1->fit().point : sv1->wcpt().point;
                        WireCell::Point p2 = sv2->fit().valid() ? sv2->fit().point : sv2->wcpt().point;
                        VertexPtr target = (point_dis(lp, p1) <= point_dis(lp, p2)) ? sv1 : sv2;
                        if (find_segment(graph, le, target)) continue;  // already linked
                        Facade::geo_point_t a = le->wcpt().point;
                        Facade::geo_point_t b = target->wcpt().point;
                        if (do_rough_path(cluster, a, b).size() < 2) { feasible = false; break; }
                        plans.emplace_back(le, target);
                    }
                    if (!feasible) continue;

                    remove_segment(graph, loser);
                    for (auto& [le, target] : plans) connect_direct(le, target);
                    cleanup_vertex(lv1);
                    cleanup_vertex(lv2);

                    SPDLOG_LOGGER_DEBUG(s_log,
                        "mvga: op1 dup-merge cluster={} removed seg len={:.2f}cm sumdQ={:.3g} "
                        "overlap={:.2f}@{:.1f}mm vs survivor len={:.2f}cm sumdQ={:.3g} reconnects={}",
                        cluster.ident(), segment_track_length(loser)/units::cm,
                        (loser == sa) ? qa : qb, frac, m_mvga_dup_tol/units::mm,
                        segment_track_length(survivor)/units::cm,
                        (loser == sa) ? qb : qa, plans.size());
                    ++n_op1;
                    flag_continue = true;
                }
            }
        }
    }

    // ---- op2: charge-less-bridge removal ---------------------------------
    if (m_mvga_bridge_mip > 0) {
        bool flag_continue = true;
        while (flag_continue && n_op2 < kEditCap) {
            flag_continue = false;
            for (SegmentPtr sg : in_scope_segments()) {
                if (seg_valid_fits(sg) == 0) continue;  // unfitted: no charge verdict possible
                auto [v1, v2] = find_vertices(graph, sg);
                if (!v1 || !v2 || !v1->descriptor_valid() || !v2->descriptor_valid()) continue;
                // Terminal stubs are op3's turf.
                if (boost::degree(v1->get_descriptor(), graph) < 2 ||
                    boost::degree(v2->get_descriptor(), graph) < 2) continue;

                double ratio = segment_median_dQ_dx(sg) / m_mip_dqdx_median;
                // Probe line for threshold calibration (TRACE: invisible at
                // the production debug level, so censuses stay sentinel-only).
                SPDLOG_LOGGER_TRACE(s_log,
                    "mvga: op2 eval cluster={} len={:.2f}cm dqdx_ratio={:.3f}",
                    cluster.ident(), segment_track_length(sg)/units::cm, ratio);
                if (ratio >= m_mvga_bridge_mip) continue;

                // Deletion feasibility: every side must stay connected to
                // the main vertex, or reconnect directly to a reachable
                // vertex within m_mvga_reconnect (268067: V_B at 4.2 cm from
                // MAIN reconnects; the roundabout cycle collapses).
                auto reach = reachable_without(graph, cluster, main_vertex, sg);
                std::vector<std::pair<VertexPtr, VertexPtr>> plans;
                bool feasible = true;
                for (VertexPtr ve : {v1, v2}) {
                    if (reach.count(ve->get_descriptor())) continue;  // still connected
                    WireCell::Point vp = ve->fit().valid() ? ve->fit().point : ve->wcpt().point;
                    auto [target, dis] = nearest_vertex_in(vp, reach, ve);
                    if (!target || dis > m_mvga_reconnect) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op2 skip cluster={} disconnected side has no reachable vertex "
                            "within {:.1f}cm (nearest {:.2f}cm)",
                            cluster.ident(), m_mvga_reconnect/units::cm,
                            target ? dis/units::cm : -1.0);
                        feasible = false; break;
                    }
                    Facade::geo_point_t a = ve->wcpt().point;
                    Facade::geo_point_t b = target->wcpt().point;
                    if (do_rough_path(cluster, a, b).size() < 2) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op2 skip cluster={} rough path failed for reconnect ({:.2f}cm)",
                            cluster.ident(), dis/units::cm);
                        feasible = false; break;
                    }
                    plans.emplace_back(ve, target);
                }
                if (!feasible) continue;

                double len = segment_track_length(sg);
                remove_segment(graph, sg);
                for (auto& [ve, target] : plans) connect_direct(ve, target);
                cleanup_vertex(v1);
                cleanup_vertex(v2);

                SPDLOG_LOGGER_DEBUG(s_log,
                    "mvga: op2 bridge-removal cluster={} len={:.2f}cm dqdx_ratio={:.3f} reconnects={}",
                    cluster.ident(), len/units::cm, ratio, plans.size());
                ++n_op2;
                flag_continue = true;
                break;
            }
        }
    }

    // ---- op3: micro-stub absorb + vertex re-seat -------------------------
    if (m_mvga_stub > 0) {
        bool flag_continue = true;
        while (flag_continue && n_op3 < kEditCap) {
            flag_continue = false;
            if (!main_vertex->descriptor_valid()) break;

            // Anchor list: the main vertex first (re-seat eligible), then --
            // when m_mvga_satellite > 0 (round 3, doc pr/51) -- every other
            // main-cluster vertex within m_mvga_satellite of the CURRENT
            // main-vertex position, absorb-only.  Re-derived every
            // iteration (mv_pt moves on a re-seat, so a hoisted list would
            // go stale); ordered_nodes keeps satellite order deterministic.
            // m_mvga_satellite == 0 ⇒ the list is exactly {main_vertex},
            // byte-identical to round 2.
            std::vector<std::pair<VertexPtr, bool>> anchors;
            anchors.emplace_back(main_vertex, true);
            if (m_mvga_satellite > 0) {
                for (const auto& vd : ordered_nodes(graph)) {
                    VertexPtr vtx = graph[vd].vertex;
                    if (!vtx || vtx->cluster() != &cluster) continue;
                    if (vtx == main_vertex) continue;
                    if (vtx->flags_any(VertexFlags::kProtectedBreak)) continue;
                    if (boost::degree(vd, graph) < 2) continue;  // nothing to shed without disconnecting
                    WireCell::Point vp = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
                    if (point_dis(vp, mv_pt) >= m_mvga_satellite) continue;
                    anchors.emplace_back(vtx, false);
                }
            }

            bool fired_here = false;
            for (auto& anchor_pair : anchors) {
                VertexPtr anchor = anchor_pair.first;
                bool is_main = anchor_pair.second;
                if (!anchor->descriptor_valid()) continue;
                auto a_vd = anchor->get_descriptor();

                // Snapshot the incident segments before editing (the
                // examine_vertices_4 pattern).
                std::vector<SegmentPtr> incident;
                for (auto edesc : sorted_out_edges(a_vd, graph)) {
                    SegmentPtr sg = graph[edesc].segment;
                    if (sg && !created.count(sg)) incident.push_back(sg);
                }
                // doc pr/86 P1b: "nothing to shed" is a terminal-absorb
                // argument -- the interposed SPLICE re-attaches every far
                // prong to the anchor, so a degree-1 main anchor whose
                // single edge is an interposed stub (evt67394's click shape;
                // 26 of the 86 sec 10.2 Class-B cases) is exactly its
                // target, not a disconnection risk.  The terminal branch
                // below re-imposes >= 2 (reason=deg1-terminal).  Knob off
                // (default) => the old gate verbatim.
                const bool deg1_splice = m_mvga_interposed_deg1 && is_main &&
                                         m_mvga_interposed && incident.size() == 1;
                if (incident.size() < 2 && !deg1_splice) continue;  // nothing to shed at this anchor

                double anchor_dis = is_main ? 0.0
                    : point_dis(anchor->fit().valid() ? anchor->fit().point : anchor->wcpt().point, mv_pt);

                // doc pr/86 P1: the interposed-splice branch may use a wider
                // candidate ceiling than the terminal absorb -- two thirds of
                // the merged-prong defect sits above m_mvga_stub (pr/86
                // sec 10.3) while the absorb, where the pr/85 sec 10.6
                // adverse movers live, keeps m_mvga_stub untouched.  Default
                // m_mvga_interposed_len = 0 keeps cand_ceiling == m_mvga_stub
                // so the knob-off control flow is byte-identical.
                const double cand_ceiling =
                    (is_main && m_mvga_interposed && m_mvga_interposed_len > m_mvga_stub)
                    ? m_mvga_interposed_len : m_mvga_stub;

                for (SegmentPtr stub : incident) {
                    double len = segment_track_length(stub);
                    if (len >= cand_ceiling) {
                        // doc pr/86 P3: op3's early gates used to decline
                        // silently; 3 of the pr/86 sec 3 top ten were
                        // unexplained purely for that reason.
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op3 decline cluster={} anchor={} len={:.2f}cm reason=ceiling",
                            cluster.ident(), is_main ? "main" : "sat", len/units::cm);
                        continue;
                    }
                    VertexPtr vf = find_other_vertex(graph, stub, anchor);
                    if (!vf || !vf->descriptor_valid()) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op3 decline cluster={} anchor={} len={:.2f}cm reason=far-invalid",
                            cluster.ident(), is_main ? "main" : "sat", len/units::cm);
                        continue;
                    }
                    if (vf == main_vertex) {  // never treat the main-vertex edge as a stub
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op3 decline cluster={} anchor={} len={:.2f}cm reason=far-is-main",
                            cluster.ident(), is_main ? "main" : "sat", len/units::cm);
                        continue;
                    }
                    // doc pr/85: an INTERPOSED stub's far vertex carries the
                    // real prong(s), so its degree is >= 2 -- the terminal
                    // absorb below can never reach it.  m_mvga_interposed
                    // opens that class at the main-vertex anchor only; off
                    // (default) this line is exactly the old terminal-only
                    // rejection.
                    const bool interposed = (boost::degree(vf->get_descriptor(), graph) != 1);
                    if (interposed && !(m_mvga_interposed && is_main)) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op3 decline cluster={} anchor={} len={:.2f}cm reason=interposed-not-main",
                            cluster.ident(), is_main ? "main" : "sat", len/units::cm);
                        continue;
                    }
                    if (vf->flags_any(VertexFlags::kProtectedBreak)) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op3 decline cluster={} anchor={} len={:.2f}cm reason=protected",
                            cluster.ident(), is_main ? "main" : "sat", len/units::cm);
                        continue;
                    }
                    // doc pr/86 P1: the wider ceiling admits candidates only
                    // into the interposed splice; the terminal absorb keeps
                    // the production m_mvga_stub ceiling.  Knob off this can
                    // never be the rejecting gate (every survivor of the
                    // cand_ceiling test already has len < m_mvga_stub).
                    if (!interposed && len >= m_mvga_stub) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op3 decline cluster={} anchor={} len={:.2f}cm reason=ceiling-terminal",
                            cluster.ident(), is_main ? "main" : "sat", len/units::cm);
                        continue;
                    }
                    // doc pr/86 P1b: a deg1_splice-admitted anchor may only
                    // take the interposed branch -- a terminal absorb of the
                    // anchor's single edge would disconnect it (and with no
                    // siblings best_overlap=0, the degeneracy gate alone
                    // could accept).  Knob off this is unreachable (the
                    // anchor gate already required incident >= 2).
                    if (!interposed && incident.size() < 2) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op3 decline cluster={} anchor={} len={:.2f}cm reason=deg1-terminal",
                            cluster.ident(), is_main ? "main" : "sat", len/units::cm);
                        continue;
                    }

                    if (interposed) {
                        // ---- op3 interposed-stub absorb (doc pr/85) ------
                        // Never falls through to the terminal overlap /
                        // degeneracy gates (an interposed stub is collinear
                        // with its continuation, not with a sibling at the
                        // anchor).  All far prongs are carried through the
                        // stub onto the main vertex; the far vertex dies.
                        std::vector<SegmentPtr> prongs;
                        for (auto edesc : sorted_out_edges(vf->get_descriptor(), graph)) {
                            SegmentPtr p = graph[edesc].segment;
                            if (p && p != stub) prongs.push_back(p);
                        }
                        if (prongs.empty()) continue;

                        // Far-end collinearity gate: the stub must continue
                        // a prong straight through vf (both directions taken
                        // FROM vf, so a continuation reads ~180 deg).
                        WireCell::Point vf_pt = vf->fit().valid() ? vf->fit().point : vf->wcpt().point;
                        WireCell::Vector dir_stub = segment_cal_dir_3vector(stub, vf_pt, 10*units::cm);
                        double best_angle = 0;
                        for (SegmentPtr p : prongs) {
                            WireCell::Vector dir_p = segment_cal_dir_3vector(p, vf_pt, 10*units::cm);
                            if (dir_stub.magnitude() == 0 || dir_p.magnitude() == 0) continue;
                            double angle = std::acos(std::clamp(
                                dir_stub.dot(dir_p) / (dir_stub.magnitude() * dir_p.magnitude()),
                                -1.0, 1.0)) / 3.1415926 * 180.0;
                            if (angle > best_angle) best_angle = angle;
                        }
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op3 eval-interposed cluster={} len={:.2f}cm vf_deg={} far_angle={:.1f}deg",
                            cluster.ident(), len/units::cm,
                            boost::degree(vf->get_descriptor(), graph), best_angle);
                        if (best_angle < m_mvga_interposed_angle) continue;

                        // All-or-nothing: every prong must pre-verify
                        // (endpoint wcpt matches, far-slot free -- B.7)
                        // before any is moved.
                        bool ok_all = true;
                        for (SegmentPtr p : prongs) {
                            if (!carry_prong_verify(graph, p, stub, vf, anchor)) { ok_all = false; break; }
                        }
                        if (!ok_all) continue;

                        const size_t vf_deg = boost::degree(vf->get_descriptor(), graph);
                        for (SegmentPtr p : prongs) {
                            carry_prong_execute(graph, p, stub, vf, anchor, dv);
                        }
                        remove_segment(graph, stub);
                        cleanup_vertex(vf);  // degree 0 now; not protected (gated above)

                        SPDLOG_LOGGER_DEBUG(s_log,
                            "mvga: op3 stub-interposed cluster={} len={:.2f}cm vf_deg={} carried={} far_angle={:.1f}deg",
                            cluster.ident(), len/units::cm, vf_deg, prongs.size(), best_angle);
                        ++n_op3;
                        fired_here = true;
                        flag_continue = true;
                        break;
                    }

                    // Gates: corridor overlap with a sibling prong at the
                    // same anchor, or point degeneracy (142421's 3-point
                    // stubs make overlap meaningless).
                    auto pts_stub = seg_points(stub);
                    size_t nfit = seg_valid_fits(stub);
                    double best_overlap = 0;
                    for (SegmentPtr sib : incident) {
                        if (sib == stub) continue;
                        double frac = path_overlap_fraction(pts_stub, seg_points(sib), m_mvga_dup_tol);
                        if (frac > best_overlap) best_overlap = frac;
                    }
                    // doc pr/86 P4: all four pr/85 sec 10.6 adverse absorbs
                    // sat at the MAIN anchor (d=0.00) while the wanted absorb
                    // the same guard blocked (evt30504) is at a SATELLITE
                    // anchor 1.26 cm away.  m_mvga_sat_dup_frac > 0 gives
                    // satellite anchors their own overlap threshold; 0
                    // (default) => m_mvga_dup_frac everywhere => byte-identical.
                    const double dup_frac_eff = (!is_main && m_mvga_sat_dup_frac > 0)
                        ? m_mvga_sat_dup_frac : m_mvga_dup_frac;
                    const bool overlap_gate = (dup_frac_eff > 0) && (best_overlap >= dup_frac_eff);
                    const bool degen_gate = (m_mvga_stub_pts > 0) &&
                                            (nfit <= static_cast<size_t>(m_mvga_stub_pts));

                    // Probe line for satellite-radius calibration (TRACE,
                    // cf. op1/op2's eval lines; op3 had none until round 3).
                    SPDLOG_LOGGER_TRACE(s_log,
                        "mvga: op3 eval cluster={} anchor={} d={:.2f}cm len={:.2f}cm nfit={} overlap={:.2f}",
                        cluster.ident(), is_main ? "main" : "sat", anchor_dis/units::cm,
                        len/units::cm, nfit, best_overlap);

                    if (!overlap_gate && !degen_gate) continue;

                    // Re-seat sub-case (131357): main-vertex anchor only,
                    // overlap-gated stub whose collinear-continuation
                    // sibling exists, and the main vertex itself is not a
                    // protected snap/break corner.  Satellite anchors are
                    // absorb-only (round 3) -- the endpoint-match / seat
                    // logic below names main_vertex explicitly.
                    SegmentPtr cont = nullptr;
                    double cont_angle = 0;
                    if (is_main && overlap_gate && m_mvga_reseat_angle > 0 &&
                        !main_vertex->flags_any(VertexFlags::kProtectedBreak)) {
                        WireCell::Point mvp = mv_pt;
                        WireCell::Vector dir_stub = segment_cal_dir_3vector(stub, mvp, 10*units::cm);
                        for (SegmentPtr sib : incident) {
                            if (sib == stub) continue;
                            WireCell::Vector dir_sib = segment_cal_dir_3vector(sib, mvp, 10*units::cm);
                            if (dir_stub.magnitude() == 0 || dir_sib.magnitude() == 0) continue;
                            double angle = std::acos(std::clamp(
                                dir_stub.dot(dir_sib) / (dir_stub.magnitude() * dir_sib.magnitude()),
                                -1.0, 1.0)) / 3.1415926 * 180.0;
                            if (angle >= m_mvga_reseat_angle && angle > cont_angle) {
                                cont_angle = angle;
                                cont = sib;
                            }
                        }
                    }

                    if (cont) {
                        // Extend the continuation sibling through the main
                        // vertex to the stub far end, then move the main
                        // vertex there (examine_structure_final_1p mechanics).
                        const auto& vec_sib = cont->wcpts();
                        const auto& vec_stub = stub->wcpts();
                        if (vec_sib.empty() || vec_stub.empty()) continue;
                        WireCell::Point main_wcpt_point = main_vertex->wcpt().point;
                        bool sib_front  = (point_dis(vec_sib.front().point,  main_wcpt_point) < 0.01*units::cm);
                        bool sib_back   = (point_dis(vec_sib.back().point,   main_wcpt_point) < 0.01*units::cm);
                        bool stub_front = (point_dis(vec_stub.front().point, main_wcpt_point) < 0.01*units::cm);
                        bool stub_back  = (point_dis(vec_stub.back().point,  main_wcpt_point) < 0.01*units::cm);
                        if (!(sib_front || sib_back) || !(stub_front || stub_back)) continue;

                        std::list<WCPoint> merged(vec_sib.begin(), vec_sib.end());
                        if (sib_front && stub_front) {
                            for (auto it = vec_stub.begin(); it != vec_stub.end(); ++it)
                                if (point_dis(it->point, merged.front().point) > 0.01*units::cm) merged.push_front(*it);
                        } else if (sib_front && !stub_front) {
                            for (auto it = vec_stub.rbegin(); it != vec_stub.rend(); ++it)
                                if (point_dis(it->point, merged.front().point) > 0.01*units::cm) merged.push_front(*it);
                        } else if (!sib_front && stub_front) {
                            for (auto it = vec_stub.begin(); it != vec_stub.end(); ++it)
                                if (point_dis(it->point, merged.back().point) > 0.01*units::cm) merged.push_back(*it);
                        } else {
                            for (auto it = vec_stub.rbegin(); it != vec_stub.rend(); ++it)
                                if (point_dis(it->point, merged.back().point) > 0.01*units::cm) merged.push_back(*it);
                        }

                        std::vector<WCPoint> new_wcpts(merged.begin(), merged.end());
                        cont->wcpts(new_wcpts);
                        std::vector<Facade::geo_point_t> cont_pts;
                        for (const auto& wcp : new_wcpts) cont_pts.push_back(wcp.point);
                        create_segment_point_cloud(cont, cont_pts, dv, "main");

                        WireCell::Point vf_pt = vf->fit().valid() ? vf->fit().point : vf->wcpt().point;
                        double reseat_dis = point_dis(vf_pt, mv_pt);
                        WCPoint vf_wcp = vf->wcpt();
                        main_vertex->wcpt(vf_wcp);
                        if (vf->fit().valid()) main_vertex->fit(vf->fit());

                        remove_segment(graph, stub);
                        cleanup_vertex(vf);
                        mv_pt = main_vertex->fit().valid() ? main_vertex->fit().point
                                                           : main_vertex->wcpt().point;

                        SPDLOG_LOGGER_DEBUG(s_log,
                            "mvga: op3 stub-reseat cluster={} anchor=main len={:.2f}cm nfit={} "
                            "overlap={:.2f} cont_angle={:.1f}deg reseat_dis={:.2f}cm",
                            cluster.ident(), len/units::cm, nfit, best_overlap, cont_angle,
                            reseat_dis/units::cm);
                    } else {
                        remove_segment(graph, stub);
                        cleanup_vertex(vf);
                        SPDLOG_LOGGER_DEBUG(s_log,
                            "mvga: op3 stub-absorb cluster={} anchor={} d={:.2f}cm len={:.2f}cm "
                            "nfit={} overlap={:.2f} gate={}",
                            cluster.ident(), is_main ? "main" : "sat", anchor_dis/units::cm,
                            len/units::cm, nfit, best_overlap,
                            overlap_gate ? "overlap" : "degenerate");
                    }
                    ++n_op3;
                    fired_here = true;
                    flag_continue = true;
                    break;
                }
                if (fired_here) break;
            }
        }
    }

    // ---- op4: one local refit --------------------------------------------
    const bool fired = (n_op1 + n_op2 + n_op3) > 0;
    if (fired) {
        track_fitter.do_multi_tracking(true, true, false, m_fit_exclusion, false, &cluster);
        SPDLOG_LOGGER_DEBUG(s_log,
            "mvga: fired cluster={} op1={} op2={} op3={} (refit done)",
            cluster.ident(), n_op1, n_op2, n_op3);
    }
    return fired;
}
