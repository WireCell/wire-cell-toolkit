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

    // doc sbnd_xin/docs/pr/86 sec 15 (round 2): straight, charge-vetoed,
    // Steiner-snapped chain from start_wcp to end_wcp.  Forked BY
    // DUPLICATION from examine_structure_2's inner recipe
    // (NeutrinoStructureExaminer.cxx lines 344-439; prototype
    // NeutrinoID_examine_structure.h lines 96-160) -- the production
    // examiner file stays byte-untouched.  Recipe: sample the straight
    // line every 0.6 cm; a sample outside every TPC or failing
    // is_good_point (0.2 cm) is bad; n_bad > 1 => return false (charge
    // genuinely bends -- caller keeps its old chain).  Survivors are
    // snapped to the steiner_pc with 0.01 cm consecutive dedup, endpoints
    // kept verbatim (graph-consistency: chains must terminate on vertex
    // wcpts).
    bool straight_steiner_chain(Facade::Cluster& cluster, TrackFitting& track_fitter,
                                const WireCell::IDetectorVolumes::pointer& dv,
                                const WCPoint& start_wcp, const WCPoint& end_wcp,
                                std::vector<WCPoint>& out,
                                double good_radius)
    {
        if (!cluster.has_pc("steiner_pc")) return false;
        const auto transform = track_fitter.get_pc_transforms()->pc_transform(
            cluster.get_scope_transform(cluster.get_default_scope()));
        auto grouping = cluster.grouping();
        if (!transform || !grouping) return false;
        const double cluster_t0 = cluster.get_cluster_t0();

        const auto start_p = start_wcp.point;
        const auto end_p = end_wcp.point;
        const double step_size = 0.6 * WireCell::units::cm;
        const double distance = point_dis(start_p, end_p);
        const int ncount = std::round(distance / step_size);

        std::vector<WireCell::Point> line_pts;
        int n_bad = 0;
        for (int i = 1; i < ncount; i++) {
            WireCell::Point test_p(
                start_p.x() + (end_p.x() - start_p.x()) / ncount * i,
                start_p.y() + (end_p.y() - start_p.y()) / ncount * i,
                start_p.z() + (end_p.z() - start_p.z()) / ncount * i);
            line_pts.push_back(test_p);
            auto test_wpid = dv->contained_by(test_p);
            if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                auto temp_p_raw = transform->backward(test_p, cluster_t0,
                                                      test_wpid.face(), test_wpid.apa());
                if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(),
                                             good_radius, 0, 0)) {
                    n_bad++;
                }
            }
            else {
                n_bad++;  // outside all TPCs => bad sample (es2 B.3)
            }
            if (n_bad > 1) return false;
        }

        const auto& steiner_pc = cluster.get_pc("steiner_pc");
        const auto& coords = cluster.get_default_scope().coords;
        const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
        const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
        const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();

        const double tol = 0.01 * WireCell::units::cm;
        std::vector<WCPoint> chain;
        chain.push_back(start_wcp);
        for (const auto& p : line_pts) {
            auto knn = cluster.kd_steiner_knn(1, p, "steiner_pc");
            if (knn.empty()) continue;
            size_t idx = knn[0].first;
            WCPoint wcp;
            wcp.point = Facade::geo_point_t(x_coords[idx], y_coords[idx], z_coords[idx]);
            if (point_dis(wcp.point, chain.back().point) > tol) chain.push_back(wcp);
        }
        if (point_dis(chain.back().point, end_wcp.point) > tol) chain.push_back(end_wcp);
        if (chain.size() < 2) return false;
        out = std::move(chain);
        return true;
    }

    // doc pr/86 sec 15 R1: re-derive the near-anchor stretch of a
    // just-carried prong.  The merged chain terminates on the anchor wcpt
    // at one end (carry_prong_execute contract); walk from that end to the
    // first point at arc length >= reach (the old stub length + the knob,
    // so the straightened span always crosses the dissolved junction), and
    // replace chain[anchor..P*] with the straight recipe.  Charge veto or
    // any precondition failure => chain untouched (concatenation stands).
    bool straighten_spliced_prong(SegmentPtr prong, const WCPoint& anchor_wcp,
                                  double reach, Facade::Cluster& cluster,
                                  TrackFitting& track_fitter,
                                  const WireCell::IDetectorVolumes::pointer& dv,
                                  double good_radius)
    {
        if (!prong) return false;
        const auto& chain = prong->wcpts();
        if (chain.size() < 3) return false;
        const double tol = 0.01 * WireCell::units::cm;
        const bool at_front = point_dis(chain.front().point, anchor_wcp.point) < tol;
        const bool at_back = point_dis(chain.back().point, anchor_wcp.point) < tol;
        if (!at_front && !at_back) return false;

        std::vector<WCPoint> fwd(chain.begin(), chain.end());
        if (!at_front) std::reverse(fwd.begin(), fwd.end());

        size_t istar = fwd.size() - 1;
        double arc = 0;
        for (size_t i = 1; i < fwd.size(); ++i) {
            arc += point_dis(fwd[i].point, fwd[i - 1].point);
            if (arc >= reach) {
                istar = i;
                break;
            }
        }
        if (istar < 1) return false;

        std::vector<WCPoint> straight;
        if (!straight_steiner_chain(cluster, track_fitter, dv, fwd.front(), fwd[istar],
                                    straight, good_radius)) {
            return false;
        }
        for (size_t i = istar + 1; i < fwd.size(); ++i) {
            if (point_dis(fwd[i].point, straight.back().point) > tol) straight.push_back(fwd[i]);
        }
        if (straight.size() < 2) return false;
        if (!at_front) std::reverse(straight.begin(), straight.end());

        prong->wcpts(straight);
        std::vector<Facade::geo_point_t> pts;
        for (const auto& wcp : straight) pts.push_back(wcp.point);
        create_segment_point_cloud(prong, pts, dv, "main");
        return true;
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
    // doc pr/83 r3: radius is now a parameter so op1 can consult its own
    // m_mvga_op1_radius (op2 keeps m_mvga_radius); radius < 0 => unscoped
    // (whole main cluster).  include_created lifts the `created` exemption
    // for the m_mvga_op1_post pass only -- op1/op2 pass false and see the
    // exact legacy set.
    auto in_scope_segments = [&](double radius, bool include_created = false) {
        std::vector<SegmentPtr> segs;
        for (const auto& vd : ordered_nodes(graph)) {
            VertexPtr vtx = graph[vd].vertex;
            if (!vtx || vtx->cluster() != &cluster) continue;
            for (auto edesc : sorted_out_edges(vd, graph)) {
                SegmentPtr sg = graph[edesc].segment;
                if (!sg) continue;
                if (std::find(segs.begin(), segs.end(), sg) != segs.end()) continue;
                if (!include_created && created.count(sg)) continue;
                bool near = (radius < 0);
                if (!near) {
                    for (const auto& p : seg_points(sg)) {
                        if (point_dis(p, mv_pt) < radius) { near = true; break; }
                    }
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

    // ---- op0: pass-through split (doc sbnd_xin/docs/pr/103, SBND 18255-405707)
    // A junction J within m_mvga_passthru of the main vertex M carries a
    // prong S whose fitted interior passes within m_mvga_passthru_tol of M
    // (405707: the muon's 2 cm back-extension to J, where the proton
    // attaches).  Production then reads the M-J connector as a DUPLICATE of
    // S (op1 overlap=1.00) and deletes it -- reconnects=0, because the
    // reconnect plan's "already linked" test sees the loser itself -- so J
    // is cut off from M until stitch_disconnected_main_cluster re-adds a
    // 4-point stub, and J's other prong (the proton) never reaches the
    // vertex ("takes a shortcut").  op0 makes the graph say what the fit
    // says: S is split at its wcpt closest to M, the remainder re-anchored
    // on M (SegmentPtr identity kept, fits stale until op4 -- the
    // carry_prong_execute contract), and the J->M piece becomes the
    // connecting stub (an existing M-J segment is reused), registered in
    // `passthru_stubs` so op3 treats it like a created stub: per-prong
    // charge-verified straighten splice of J's remaining prongs onto M,
    // admitted in pass 0.  Knob 0 (default) => block never runs and
    // `passthru_stubs` stays empty => byte-identical.
    IndexedSegmentSet passthru_stubs;
    int n_op0 = 0;
    if (m_mvga_passthru > 0 && main_vertex->descriptor_valid()) {
        const double minrem = 3*units::cm;   // S must continue at least this far beyond M
        const double tol = 0.01*units::cm;
        const WCPoint mv_wcp = main_vertex->wcpt();
        const auto mv_w = mv_wcp.point;
        bool flag_continue = true;
        while (flag_continue && n_op0 < kEditCap) {
            flag_continue = false;
            for (const auto& vd : ordered_nodes(graph)) {
                if (flag_continue) break;
                VertexPtr J = graph[vd].vertex;
                if (!J || J->cluster() != &cluster || J == main_vertex) continue;
                if (!J->descriptor_valid()) continue;
                if (J->flags_any(VertexFlags::kProtectedBreak)) continue;
                if (boost::degree(vd, graph) < 2) continue;
                WireCell::Point jp = J->fit().valid() ? J->fit().point : J->wcpt().point;
                const double dj = point_dis(jp, mv_pt);
                if (dj >= m_mvga_passthru || dj < 0.05*units::cm) continue;

                std::vector<SegmentPtr> inc_J;
                for (auto edesc : sorted_out_edges(vd, graph)) {
                    SegmentPtr sg = graph[edesc].segment;
                    if (sg) inc_J.push_back(sg);
                }
                SegmentPtr link = find_segment(graph, J, main_vertex);

                // best candidate S: smallest miss distance to M's wcpt
                SegmentPtr best_S;
                size_t best_k = 0;
                double best_miss = 1e9, best_arc = 0, best_rem = 0;
                bool best_front = true;
                for (SegmentPtr S : inc_J) {
                    if (S == link) continue;
                    if (created.count(S) || passthru_stubs.count(S)) continue;
                    VertexPtr vfar = find_other_vertex(graph, S, J);
                    if (!vfar || vfar == main_vertex || !vfar->descriptor_valid()) continue;
                    if (find_segment(graph, main_vertex, vfar)) continue;  // B.7: (M, far) slot must be free
                    const auto& w = S->wcpts();
                    if (w.size() < 3) continue;
                    const bool front = point_dis(w.front().point, J->wcpt().point) < tol;
                    if (!front && point_dis(w.back().point, J->wcpt().point) >= tol) continue;
                    // Orient J-first; test the closest point on the wcpt POLYLINE
                    // (wcpt chains are sparse -- 405707's muon has a 5.8 cm first
                    // hop -- so a nearest-wcpt test cannot see the pass-through).
                    std::vector<WCPoint> fw(w.begin(), w.end());
                    if (!front) std::reverse(fw.begin(), fw.end());
                    std::vector<double> arc(fw.size(), 0.0);
                    for (size_t i = 1; i < fw.size(); ++i) arc[i] = arc[i-1] + point_dis(fw[i].point, fw[i-1].point);
                    const double total = arc.back();
                    size_t kbest = fw.size();   // split hop: between fw[kbest] and fw[kbest+1]
                    double miss = 1e9, aj = -1;
                    for (size_t i = 0; i + 1 < fw.size(); ++i) {
                        const auto a = fw[i].point, b = fw[i+1].point;
                        const auto ab = b - a;
                        const double L2 = ab.dot(ab);
                        if (L2 <= 0) continue;
                        double t = (mv_w - a).dot(ab) / L2;
                        t = std::clamp(t, 0.0, 1.0);
                        const auto q = a + ab * t;
                        const double d = point_dis(q, mv_w);
                        const double aq = arc[i] + t * std::sqrt(L2);    // arc from J
                        if (aq < 0.3*units::cm || aq > m_mvga_passthru + m_mvga_passthru_tol) continue;
                        if (d < miss) { miss = d; kbest = i; aj = aq; }
                    }
                    bool other = false;   // J must keep a prong besides S and the link
                    for (SegmentPtr o : inc_J) if (o != S && o != link) { other = true; break; }
                    SPDLOG_LOGGER_TRACE(s_log,
                        "mvga: op0 eval cluster={} dJ={:.2f}cm S_len={:.2f}cm npts={} miss={:.2f}cm arc={:.2f}cm rem={:.2f}cm other={}",
                        cluster.ident(), dj/units::cm, total/units::cm, fw.size(),
                        (kbest == fw.size()) ? -1.0 : miss/units::cm, aj/units::cm,
                        (kbest == fw.size()) ? -1.0 : (total - aj)/units::cm, other ? 1 : 0);
                    if (kbest == fw.size() || miss >= m_mvga_passthru_tol) continue;
                    if (total - aj < minrem) continue;
                    if (!other) continue;
                    if (miss < best_miss) {
                        best_miss = miss; best_S = S; best_k = kbest; best_front = front;
                        best_arc = aj; best_rem = total - aj;
                    }
                }
                if (!best_S) continue;

                const auto& w = best_S->wcpts();
                std::vector<WCPoint> fwd(w.begin(), w.end());
                if (!best_front) std::reverse(fwd.begin(), fwd.end());
                const size_t k = best_k;   // M lies on the hop fwd[k] -> fwd[k+1]
                std::vector<WCPoint> stub_w(fwd.begin(), fwd.begin() + k + 1);   // J .. fwd[k]
                if (point_dis(stub_w.back().point, mv_w) > tol) stub_w.push_back(mv_wcp);
                std::vector<WCPoint> rem_w;
                rem_w.push_back(mv_wcp);                                             // M .. far
                for (size_t i = k + 1; i < fwd.size(); ++i)
                    if (point_dis(fwd[i].point, rem_w.back().point) > tol) rem_w.push_back(fwd[i]);
                if (rem_w.size() < 2 || stub_w.size() < 2) continue;
                VertexPtr vfar = find_other_vertex(graph, best_S, J);
                if (!vfar) continue;

                SegmentPtr stub = link;
                const bool stub_new = !stub;
                if (stub_new) {
                    std::vector<Facade::geo_point_t> spts;
                    for (const auto& q : stub_w) spts.push_back(q.point);
                    stub = create_segment_for_cluster(cluster, dv, spts, 0);
                    if (!stub) continue;
                    add_segment(graph, stub, J, main_vertex);
                }
                passthru_stubs.insert(stub);

                best_S->wcpts(rem_w);
                std::vector<Facade::geo_point_t> rpts;
                for (const auto& q : rem_w) rpts.push_back(q.point);
                create_segment_point_cloud(best_S, rpts, dv, "main");
                // Trim the (stale) fit points to the remainder: op1/op3 read
                // seg_points()/segment_track_length() from the fits until the
                // op4 refit, and the J->M stretch must not be claimed twice.
                {
                    auto& fits = best_S->fits();
                    if (fits.size() >= 3) {
                        size_t imin = 0; double dmin = 1e9;
                        for (size_t i = 0; i < fits.size(); ++i) {
                            const double d = point_dis(fits[i].point, mv_w);
                            if (d < dmin) { dmin = d; imin = i; }
                        }
                        const bool j_first = point_dis(fits.front().point, J->wcpt().point)
                                           < point_dis(fits.back().point, J->wcpt().point);
                        std::vector<Fit> kept;
                        if (j_first) kept.assign(fits.begin() + imin, fits.end());
                        else         kept.assign(fits.begin(), fits.begin() + imin + 1);
                        if (kept.size() >= 2) best_S->fits(kept);
                    }
                }
                remove_segment(graph, best_S);
                add_segment(graph, best_S, main_vertex, vfar);

                SPDLOG_LOGGER_DEBUG(s_log,
                    "mvga: op0 passthru-split cluster={} dJ={:.2f}cm miss={:.2f}cm arc={:.2f}cm rem={:.2f}cm "
                    "J_deg={} stub={} stub_npts={}",
                    cluster.ident(), dj/units::cm, best_miss/units::cm, best_arc/units::cm, best_rem/units::cm,
                    inc_J.size(), stub_new ? "created" : "existing", stub_w.size());
                ++n_op0;
                flag_continue = true;
            }
        }
    }

    // ---- op1: duplicate-corridor merge -----------------------------------
    // doc pr/83 r3 sec 9.2/9.3: op1's scope and threshold decouple from
    // op2/op3's via m_mvga_op1_radius / m_mvga_op1_dup_frac.  Both default
    // 0 => the shared members apply => byte-identical legacy.
    const double op1_radius = (m_mvga_op1_radius != 0) ? m_mvga_op1_radius : m_mvga_radius;
    const double op1_dup_frac = (m_mvga_op1_dup_frac > 0) ? m_mvga_op1_dup_frac : m_mvga_dup_frac;
    if (m_mvga_dup_frac > 0 && m_mvga_dup_tol > 0) {
        bool flag_continue = true;
        while (flag_continue && n_op1 < kEditCap) {
            flag_continue = false;
            auto segs = in_scope_segments(op1_radius);
            // pr/83 round 2 diagnostic: which segments op1 even considered,
            // and the main-vertex point the mvga_radius scope is centered on
            // -- a duplicate pair outside scope never reaches the eval TRACE
            // above and looks identical, in the log, to one that was in
            // scope but below the overlap/angle threshold.
            if (n_op1 == 0) {
                SPDLOG_LOGGER_TRACE(s_log, "mvga: op1 scope cluster={} mv=({:.2f},{:.2f},{:.2f})cm n_in_scope={}",
                    cluster.ident(), mv_pt.x()/units::cm, mv_pt.y()/units::cm, mv_pt.z()/units::cm, segs.size());
                for (SegmentPtr sg : segs) {
                    SPDLOG_LOGGER_TRACE(s_log, "mvga: op1 scope-member cluster={} len={:.2f}cm",
                        cluster.ident(), segment_track_length(sg)/units::cm);
                }
            }
            for (size_t i = 0; i + 1 < segs.size() && !flag_continue; ++i) {
                for (size_t j = i + 1; j < segs.size() && !flag_continue; ++j) {
                    SegmentPtr sa = segs[i];
                    SegmentPtr sb = segs[j];
                    // doc pr/103: the connector op0 just established is never a
                    // duplicate to be merged away (405707: with S's stale fits
                    // it reads overlap=1.00).  Empty set when the knob is off.
                    if (passthru_stubs.count(sa) || passthru_stubs.count(sb)) continue;
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
                    // doc pr/83 r3: the relaxed m_mvga_op1_dup_frac applies
                    // only when BOTH members are >= 10 cm -- the census's
                    // own min-length for a duplicate finding (pr83_dup_metric
                    // / pr83r2_census.py), i.e. exactly the class the knob
                    // exists to recover (404684: 15.8/61.5 cm at 0.74).
                    // Measured adversity is all short riders: 390842's
                    // 1.91 cm rider merges at 0.75 under a flat 0.7 and the
                    // kine walk then loses the whole 1.03 GeV muon chain
                    // (Enu 1148 -> 10 MeV); 285567 (1.98 cm @ 0.75) and
                    // 268067 (1.66 cm @ 0.75) flip nue_score the same way.
                    // A sub-10 cm rider's overlap fraction rides on a
                    // handful of points and stays gated at m_mvga_dup_frac
                    // (production 0.8).  Knob off => both gates equal =>
                    // byte-identical at any length.
                    const double frac_gate =
                        (std::min(la, lb) >= 10*units::cm) ? op1_dup_frac
                                                           : m_mvga_dup_frac;
                    if (frac < frac_gate) continue;
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
                            // pr/83 round 2 diagnostic: an overlap-gate pass
                            // that still declines needs its angle on record --
                            // op1's own eval TRACE (above) never reaches this
                            // far when frac < dup_frac, so a post-overlap
                            // decline was otherwise silent.
                            SPDLOG_LOGGER_TRACE(s_log,
                                "mvga: op1 angle cluster={} overlap={:.2f} angle={:.1f}deg gate={}",
                                cluster.ident(), frac, ang,
                                (ang > m_mvga_dup_angle) ? "decline" : "pass");
                            if (ang > m_mvga_dup_angle) continue;
                        } else {
                            SPDLOG_LOGGER_TRACE(s_log,
                                "mvga: op1 angle cluster={} overlap={:.2f} angle=zero-chord gate=pass",
                                cluster.ident(), frac);
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
                    if (!lv1 || !lv2 || !sv1 || !sv2) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op1 find-vertices-failed cluster={} overlap={:.2f} lv1={} lv2={} sv1={} sv2={}",
                            cluster.ident(), frac, (bool)lv1, (bool)lv2, (bool)sv1, (bool)sv2);
                        continue;
                    }

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
                    if (!feasible) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op1 reconnect-infeasible cluster={} loser_len={:.2f}cm",
                            cluster.ident(), segment_track_length(loser)/units::cm);
                        continue;
                    }

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
            for (SegmentPtr sg : in_scope_segments(m_mvga_radius)) {
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

    // doc pr/86 sec 15 R2: op3 and op3.5 interleave.  An op3.5 collapse can
    // turn the next junction up the approach chain into an exact
    // interposed-splice case (349945: merging the outer zigzag pair leaves
    // an elbow vertex whose main-incident stub + far prongs are op3's
    // shape), but op3 has already finished by then -- so re-run op3 after
    // any pass in which op3.5 fired.  Knob off => op3.5 never fires => a
    // single pass, op3 exactly once => byte-identical.  Caps: the shared
    // n_op3/n_op3b kEditCap counters persist across passes; the pass count
    // itself is bounded as belt-and-braces.
    int n_op3b = 0;
    for (int audit_pass = 0; audit_pass < 4; ++audit_pass) {
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
                    if (!sg) continue;
                    // doc pr/86 sec 15: op1/op2 reconnect products stay
                    // exempt in pass 0 (the delete/recreate cycling guard),
                    // but on an op3.5 re-entry pass they are legitimate
                    // SPLICE stubs (349945: the op1-created 3.06 cm elbow
                    // stub is the last hop of the owner's direct track and
                    // was starved forever).  op1/op2 never re-run in later
                    // passes, so no cycle; the terminal branch below still
                    // declines created stubs (reason=created-terminal).
                    // audit_pass > 0 requires an op3.5 fire => knob on =>
                    // pass 0 alone is byte-identical.
                    if (created.count(sg) && audit_pass == 0) continue;
                    incident.push_back(sg);
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
                    // doc pr/86 sec 15: a created reconnect admitted on a
                    // re-entry pass may only be SPLICED -- terminal-absorbing
                    // it would delete what op1/op2 just built.  Unreachable
                    // in pass 0 (created segs are filtered there).
                    if (!interposed && created.count(stub)) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op3 decline cluster={} anchor={} len={:.2f}cm reason=created-terminal",
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
                        // doc pr/86 sec 15: created reconnects carry no fits,
                        // so the fit-based direction is (0,0,0) and the angle
                        // gate silently declines the very stub the re-entry
                        // pass admitted (349945's elbow).  Mirror the
                        // centroid-direction arithmetic of
                        // segment_cal_dir_3vector(seg, p, cut) over
                        // seg_points() (wcpt fallback), for CREATED segments
                        // only -- production segments keep the round-1
                        // zero-dir decline semantics (284145).
                        auto dir_from_points = [&](SegmentPtr sg) {
                            WireCell::Point acc(0, 0, 0);
                            int nc = 0;
                            for (const auto& q : seg_points(sg)) {
                                if (point_dis(q, vf_pt) < 10*units::cm) { acc = acc + q; ++nc; }
                            }
                            if (!nc) return WireCell::Vector(0, 0, 0);
                            WireCell::Vector v = acc*(1.0/nc) - vf_pt;
                            if (v.magnitude() > 0) v = v.norm();
                            return v;
                        };
                        // Fallback is additionally gated on a round-2 knob:
                        // created PRONGS reach this loop even in pass 0
                        // (prongs are collected without the created filter),
                        // and an ungated fallback changed knob-off behavior
                        // on 284145/319611 (the round-1 zero-dir declines).
                        const bool r2_dir_fallback =
                            (m_mvga_splice_straighten > 0 || m_mvga_approach_collapse > 0);
                        WireCell::Vector dir_stub = segment_cal_dir_3vector(stub, vf_pt, 10*units::cm);
                        if (dir_stub.magnitude() == 0 && r2_dir_fallback && created.count(stub)) dir_stub = dir_from_points(stub);
                        double best_angle = 0;
                        for (SegmentPtr p : prongs) {
                            WireCell::Vector dir_p = segment_cal_dir_3vector(p, vf_pt, 10*units::cm);
                            if (dir_p.magnitude() == 0 && r2_dir_fallback && created.count(p)) dir_p = dir_from_points(p);
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

                        // doc pr/86 sec 15: created-stub splice (re-entry
                        // passes only -- created stubs are filtered from
                        // `incident` in pass 0, so knob-off never reaches
                        // this).  The far-angle gate asks "is the stub a
                        // collinear continuation", the wrong question for a
                        // connection op1 itself built: 349945's elbow reads
                        // 52.6 deg while the straight chord from the anchor
                        // past the junction is charge-covered.  Replace the
                        // angle gate with a PER-PRONG charge-verified
                        // straighten gate: carry exactly those prongs whose
                        // straightened approach passes the es2 charge veto;
                        // the junction and stub survive for the rest.
                        // doc pr/103: op0's pass-through stubs take the same
                        // per-prong charge-verified path (they are never in
                        // `created`, so pass 0 admits them; knob off => empty set).
                        // doc pr/103: `mvga_interposed_fallback` -- when the
                        // far-angle gate would decline (65289 118.5 deg,
                        // 345633 126.2, 400856 122.4, 287517 64.1 at a 0.64 cm
                        // stub: the same 3-track vertex split in two), take
                        // the same per-prong charge-verified path instead of
                        // declining; prongs whose straight chord from the
                        // anchor is off charge stay on the junction.  Off
                        // (default) => the old angle decline verbatim.
                        // Restricted to a DEGREE-2 far vertex (stub + ONE prong:
                        // a single prong's elbow near the vertex -- the op3.5
                        // shape without its whole-chain chord).  A far vertex
                        // carrying >= 2 long prongs is a genuine multi-track
                        // junction; when the click-anchored A/B was run on such
                        // cases (65289 / 66712 / 345633) the owner's vertex was
                        // J itself (2.4-3.7 cm from the reco main), and
                        // collapsing J into M pulled every prong away from it
                        // -- a vertex-placement question, not a topology one.
                        // Angle floor (doc pr/103 sec 6 ledger): best_angle == 0
                        // means NO measurable prong direction (fit-less prong)
                        // -- the legacy gate declines those and so must the
                        // fallback (235435: two 0.0-deg fires cut two prongs
                        // off the vertex, Enu 818 -> 555 MeV); and a prong
                        // leaving J nearly parallel to the stub (< ~45 deg) is a
                        // hairpin / back-fold whose carry only shortens the
                        // track (389588, 314705, 394532).  Production floor 45.
                        const bool angle_fallback = m_mvga_interposed_fallback &&
                                                    best_angle > 0 &&
                                                    best_angle >= m_mvga_interposed_fallback_min_angle &&
                                                    best_angle < m_mvga_interposed_angle &&
                                                    boost::degree(vf->get_descriptor(), graph) == 2 &&
                                                    !created.count(stub) && !passthru_stubs.count(stub);
                        if (created.count(stub) || passthru_stubs.count(stub) || angle_fallback) {
                            if (m_mvga_splice_straighten <= 0) continue;  // no charge gate available
                            // doc pr/83 r3 (sec 8.5 fallback): decline a
                            // multi-prong carry outright above the cap --
                            // the stub stays as the shared trunk.  0
                            // (default) => unlimited => byte-identical.
                            if (m_mvga_carry_max > 0 &&
                                static_cast<int>(prongs.size()) > m_mvga_carry_max) {
                                SPDLOG_LOGGER_TRACE(s_log,
                                    "mvga: op3 created-splice decline cluster={} reason=carry-max prongs={} cap={}",
                                    cluster.ident(), prongs.size(), m_mvga_carry_max);
                                continue;
                            }
                            double stub_arc = 0;
                            const auto& swc0 = stub->wcpts();
                            for (size_t i = 1; i < swc0.size(); ++i)
                                stub_arc += point_dis(swc0[i].point, swc0[i-1].point);
                            const double good_r = (m_mvga_straighten_radius > 0)
                                ? m_mvga_straighten_radius : 0.2*units::cm;
                            const WCPoint anchor_wcp = anchor->wcpt();
                            int carried = 0;
                            for (SegmentPtr p : prongs) {
                                if (!carry_prong_verify(graph, p, stub, vf, anchor)) continue;
                                const auto& pw = p->wcpts();
                                if (pw.size() < 2) continue;
                                std::vector<WCPoint> fwd(pw.begin(), pw.end());
                                if (point_dis(pw.front().point, vf->wcpt().point) >= 0.01*units::cm)
                                    std::reverse(fwd.begin(), fwd.end());
                                size_t istar = fwd.size() - 1;
                                double arc = 0;
                                for (size_t i = 1; i < fwd.size(); ++i) {
                                    arc += point_dis(fwd[i].point, fwd[i-1].point);
                                    if (arc >= m_mvga_splice_straighten) { istar = i; break; }
                                }
                                std::vector<WCPoint> chk;
                                if (!straight_steiner_chain(cluster, track_fitter, dv,
                                                            anchor_wcp, fwd[istar], chk, good_r)) {
                                    SPDLOG_LOGGER_TRACE(s_log,
                                        "mvga: op3 created-splice decline cluster={} reason=charge-veto",
                                        cluster.ident());
                                    continue;
                                }
                                carry_prong_execute(graph, p, stub, vf, anchor, dv);
                                straighten_spliced_prong(p, anchor_wcp,
                                                         stub_arc + m_mvga_splice_straighten,
                                                         cluster, track_fitter, dv, good_r);
                                ++carried;
                            }
                            if (!carried) continue;
                            if (vf->descriptor_valid() &&
                                boost::degree(vf->get_descriptor(), graph) == 1) {
                                remove_segment(graph, stub);
                                cleanup_vertex(vf);
                            }
                            SPDLOG_LOGGER_DEBUG(s_log,
                                "mvga: op3 created-splice cluster={} stub_arc={:.2f}cm carried={} vf_kept={} kind={} far_angle={:.1f}deg",
                                cluster.ident(), stub_arc/units::cm, carried,
                                vf->descriptor_valid() ? 1 : 0,
                                angle_fallback ? "angle-fallback" : (passthru_stubs.count(stub) ? "passthru" : "created"),
                                best_angle);
                            ++n_op3;
                            fired_here = true;
                            flag_continue = true;
                            break;
                        }

                        if (best_angle < m_mvga_interposed_angle) continue;

                        // doc pr/83 r3 (sec 8.5 fallback): decline a
                        // multi-prong carry outright above the cap.  0
                        // (default) => unlimited => byte-identical.
                        if (m_mvga_carry_max > 0 &&
                            static_cast<int>(prongs.size()) > m_mvga_carry_max) {
                            SPDLOG_LOGGER_TRACE(s_log,
                                "mvga: op3 stub-interposed decline cluster={} reason=carry-max prongs={} cap={}",
                                cluster.ident(), prongs.size(), m_mvga_carry_max);
                            continue;
                        }

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

                        // doc pr/86 sec 15 R1: straighten each carried
                        // prong's near-anchor stretch (stub length + reach,
                        // so the span crosses the dissolved junction) so the
                        // op4 refit seeds straight instead of around the
                        // concatenation kink.  Charge veto keeps genuinely
                        // bending prongs untouched.  Knob 0 (default) =>
                        // this block never runs => byte-identical.
                        if (m_mvga_splice_straighten > 0) {
                            const WCPoint anchor_wcp = anchor->wcpt();
                            // Fit-based len is 0 for a created (fitless)
                            // stub; the wcpt arc is the real junction
                            // distance the straighten span must cross.
                            double stub_arc = 0;
                            const auto& swc = stub->wcpts();
                            for (size_t i = 1; i < swc.size(); ++i)
                                stub_arc += point_dis(swc[i].point, swc[i-1].point);
                            const double reach = std::max(len, stub_arc) + m_mvga_splice_straighten;
                            const double good_r = (m_mvga_straighten_radius > 0)
                                ? m_mvga_straighten_radius : 0.2*units::cm;
                            int n_straight = 0;
                            for (SegmentPtr p : prongs) {
                                if (straighten_spliced_prong(p, anchor_wcp, reach,
                                                             cluster, track_fitter, dv,
                                                             good_r)) {
                                    ++n_straight;
                                }
                            }
                            SPDLOG_LOGGER_DEBUG(s_log,
                                "mvga: op3 splice-straighten cluster={} carried={} straightened={} reach={:.2f}cm",
                                cluster.ident(), prongs.size(), n_straight, reach/units::cm);
                        }

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

    // ---- op3.5: near-vertex approach collapse (doc pr/86 sec 15 R2) ------
    // Late re-run of the examine_structure_2 merge restricted to degree-2
    // junction vertices within m_mvga_approach_collapse of the main vertex:
    // improve_vertex builds the multi-segment zigzag approaches (349945's
    // 14 cm polyline over 5.8 cm) AFTER examine_structure_2 last ran, and
    // mvga is the last pass that edits this graph.  Iterated to a fixed
    // point; each fire replaces two chain segments with one straight,
    // charge-vetoed, Steiner-snapped segment (co-located-endpoint case
    // skipped -- divergence from es2's B.7 vertex merge, stated in the
    // doc).  Knob 0 (default) => the whole block is skipped =>
    // byte-identical.
    const int op3b_before = n_op3b;
    if (m_mvga_approach_collapse > 0) {
        bool go = true;
        while (go && n_op3b < kEditCap) {
            go = false;
            for (const auto& vd : ordered_nodes(graph)) {
                VertexPtr vtx = graph[vd].vertex;
                if (!vtx || vtx->cluster() != &cluster) continue;
                if (vtx == main_vertex) continue;
                if (boost::degree(vd, graph) != 2) continue;
                if (vtx->flags_any(VertexFlags::kProtectedBreak)) continue;
                WireCell::Point vp = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
                if (point_dis(vp, mv_pt) >= m_mvga_approach_collapse) continue;

                auto edges2 = sorted_out_edges(vd, graph);
                SegmentPtr sg1 = graph[edges2[0]].segment;
                SegmentPtr sg2 = graph[edges2[1]].segment;
                if (!sg1 || !sg2) continue;
                VertexPtr vtx1 = find_other_vertex(graph, sg1, vtx);
                VertexPtr vtx2 = find_other_vertex(graph, sg2, vtx);
                if (!vtx1 || !vtx2 || vtx1 == vtx2) continue;
                // co-located endpoints: es2 merges the vertices (B.7); a
                // near-vertex zigzag junction is never co-located -- skip.
                if (point_dis(vtx1->wcpt().point, vtx2->wcpt().point) < 0.01*units::cm) continue;
                // setS edge aliasing (examine_structure_review.md B.7): the
                // (vtx1, vtx2) slot must be free.
                if (find_segment(graph, vtx1, vtx2)) continue;

                // doc pr/99 round 2 guards (design block in
                // NeutrinoPatternBase.h at the m_mvga_ac_* members).  Each
                // knob 0/false => its test is skipped => byte-identical.
                if (m_mvga_ac_no_cascade && (created.count(sg1) || created.count(sg2))) {
                    SPDLOG_LOGGER_DEBUG(s_log,
                        "mvga: op3.5 decline cluster={} d={:.2f}cm reason=no-cascade",
                        cluster.ident(), point_dis(vp, mv_pt)/units::cm);
                    continue;
                }
                const double ac_chord =
                    point_dis(vtx1->wcpt().point, vtx2->wcpt().point);
                if (m_mvga_ac_chord_max > 0 && ac_chord > m_mvga_ac_chord_max) {
                    SPDLOG_LOGGER_DEBUG(s_log,
                        "mvga: op3.5 decline cluster={} d={:.2f}cm chord={:.2f}cm reason=chord-cap",
                        cluster.ident(), point_dis(vp, mv_pt)/units::cm,
                        ac_chord/units::cm);
                    continue;
                }

                std::vector<WCPoint> straight;
                // pr/99 round 2: the collapse chord gets its own veto radius
                // (prototype es2: 0.2 cm) when m_mvga_ac_veto_radius > 0;
                // legacy falls back to the R1 straighten radius rule.
                const double good_r = (m_mvga_ac_veto_radius > 0)
                    ? m_mvga_ac_veto_radius
                    : ((m_mvga_straighten_radius > 0)
                       ? m_mvga_straighten_radius : 0.2*units::cm);
                if (!straight_steiner_chain(cluster, track_fitter, dv,
                                            vtx1->wcpt(), vtx2->wcpt(), straight, good_r)) {
                    SPDLOG_LOGGER_TRACE(s_log,
                        "mvga: op3.5 decline cluster={} d={:.2f}cm reason=charge-veto",
                        cluster.ident(), point_dis(vp, mv_pt)/units::cm);
                    continue;
                }

                auto new_seg = make_segment();
                new_seg->wcpts(straight).cluster(&cluster).dirsign(0);
                std::vector<Facade::geo_point_t> mpts;
                for (const auto& wcp : straight) mpts.push_back(wcp.point);
                create_segment_point_cloud(new_seg, mpts, dv, "main");

                remove_segment(graph, sg1);
                remove_segment(graph, sg2);
                remove_vertex(graph, vtx);
                add_segment(graph, new_seg, vtx1, vtx2);
                // doc pr/86 sec 15: collapse products join `created` so the
                // op3 re-entry pass may SPLICE but never terminal-absorb
                // them (67394: the absorb deleted a fresh 0.56 cm collapse
                // product at the main anchor -- nfit=0, overlap=1.00 -- and
                // moved the vertex 1.02 cm off the click).
                created.insert(new_seg);

                SPDLOG_LOGGER_DEBUG(s_log,
                    "mvga: op3.5 approach-collapse cluster={} d={:.2f}cm chord={:.2f}cm npts={}",
                    cluster.ident(), point_dis(vp, mv_pt)/units::cm,
                    point_dis(vtx1->wcpt().point, vtx2->wcpt().point)/units::cm,
                    straight.size());
                ++n_op3b;
                go = true;
                break;
            }
        }
    }

    // No op3.5 fire this pass => nothing new for op3 => stop.
    if (n_op3b == op3b_before) break;
    }  // audit_pass (op3 <-> op3.5 interleave)

    // ---- op4: one local refit --------------------------------------------
    // doc pr/103: op0 re-anchors wcpt chains (fits stale) -> it needs the refit too; n_op0 == 0 when the knob is off.
    const bool fired_ops = (n_op0 + n_op1 + n_op2 + n_op3 + n_op3b) > 0;
    if (fired_ops) {
        track_fitter.do_multi_tracking(true, true, false, m_fit_exclusion, false, &cluster);
        SPDLOG_LOGGER_DEBUG(s_log,
            "mvga: fired cluster={} op1={} op2={} op3={} (refit done)",
            cluster.ident(), n_op1, n_op2, n_op3);
        if (n_op0 > 0) {
            SPDLOG_LOGGER_DEBUG(s_log, "mvga: op0 fired cluster={} splits={}", cluster.ident(), n_op0);
        }
        if (n_op3b > 0) {
            SPDLOG_LOGGER_DEBUG(s_log,
                "mvga: op3.5 fired cluster={} collapses={}", cluster.ident(), n_op3b);
        }
    }

    // ---- op1-post: duplicate-corridor pass over the REFITTED op3 products --
    // doc pr/83 r3 (sec 8, class A): the pr/85 carry + pr/86 straighten give
    // N carried prongs the SAME near-anchor trunk geometry, and op1 -- which
    // runs BEFORE op3 and skips the `created` set -- can never see it
    // (138009: six stacked prongs, 12 dup pairs, op1=0).  This pass re-runs
    // op1's exact metric and merge recipe AFTER the op4 refit, with the
    // `created` exemption lifted; benign carries overlap nothing and are
    // untouched.  Position is load-bearing: measured pre-refit (74544 /
    // 138009 TRACE), the carried prongs read overlap 0.30-0.50 -- it is the
    // refit that routes them onto the shared charge ridge (>= 0.7, the
    // census geometry) -- so this pass must see REFITTED points, and a
    // second refit (below) re-derives the survivor set.  Duplicated from
    // op1 above rather than shared (op1 stays byte-untouched); false
    // (default) => skipped => byte-identical.
    int n_op1_post = 0;
    if (m_mvga_op1_post && m_mvga_dup_frac > 0 && m_mvga_dup_tol > 0) {
        // Iterate to a FIXED POINT, refitting between rounds: the refit that
        // follows a merge round can itself pull the remaining prongs onto
        // the newly-consolidated charge ridge (138009 TRACE: the last pair
        // reads 0.50 before the first round's 5 merges + refit, 0.71 after
        // -- above threshold, invisible to a single pass).  kEditCap bounds
        // total merges, so the outer loop terminates.
        int post_rounds = 0;
        while (post_rounds < kEditCap && n_op1_post < kEditCap) {
        const int post_before = n_op1_post;
        bool flag_continue = true;
        while (flag_continue && n_op1_post < kEditCap) {
            flag_continue = false;
            auto segs = in_scope_segments(op1_radius, /*include_created=*/true);
            for (size_t i = 0; i + 1 < segs.size() && !flag_continue; ++i) {
                for (size_t j = i + 1; j < segs.size() && !flag_continue; ++j) {
                    SegmentPtr sa = segs[i];
                    SegmentPtr sb = segs[j];
                    // doc pr/103: the connector op0 just established is never a
                    // duplicate to be merged away (405707: with S's stale fits
                    // it reads overlap=1.00).  Empty set when the knob is off.
                    if (passthru_stubs.count(sa) || passthru_stubs.count(sb)) continue;
                    double la = segment_track_length(sa);
                    double lb = segment_track_length(sb);
                    SegmentPtr shorter = (la <= lb) ? sa : sb;
                    SegmentPtr longer  = (la <= lb) ? sb : sa;
                    auto pts_s = seg_points(shorter);
                    auto pts_l = seg_points(longer);
                    if (pts_s.size() <= static_cast<size_t>(m_mvga_stub_pts)) continue;
                    double frac = path_overlap_fraction(pts_s, pts_l, m_mvga_dup_tol);
                    if (frac > 0.25) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op1-post eval cluster={} pair len {:.2f}/{:.2f}cm npts {}/{} overlap={:.2f}",
                            cluster.ident(), segment_track_length(shorter)/units::cm,
                            segment_track_length(longer)/units::cm,
                            pts_s.size(), pts_l.size(), frac);
                    }
                    if (frac < op1_dup_frac) continue;
                    // doc pr/99 round 2 (m_mvga_dup_starved_mip design block
                    // in NeutrinoPatternBase.h): a non-null forced loser
                    // overrides the integrated-charge rule below -- set only
                    // on the angle-decline branch when exactly one member is
                    // charge-starved and the other healthy.
                    SegmentPtr forced_loser = nullptr;
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
                            if (ang > m_mvga_dup_angle) {
                                // Knobs 0 (default) => the decline stands
                                // exactly as before => byte-identical.
                                // Both tests required: the fitter SPLITS
                                // the shared corridor's charge across the
                                // pair (70084 measured 1.16/0.62 -- an
                                // absolute MIP floor can never separate
                                // them post-refit), so the discriminator
                                // is the op1-proj-style pair ASYMMETRY,
                                // plus an absolute cap on the loser so a
                                // genuine proton+MIP V (muon ~1.0) is
                                // never mistaken for a starved chord.
                                if (!(m_mvga_dup_starved_asym > 0) ||
                                    !(m_mvga_dup_starved_mip > 0)) continue;
                                if (seg_valid_fits(sa) == 0 || seg_valid_fits(sb) == 0) continue;
                                // Span comparability: a projective duplicate
                                // shares the corridor over its WHOLE span
                                // (70084: 15.7 vs 13.0 cm, ratio 0.83); a
                                // track paired with its own Bragg-peak stub
                                // or a short spur is NOT a duplicate (138009:
                                // 21.0 vs 3.2 cm, ratio 0.15 -- the first
                                // knob-on campaign deleted the 21 cm electron
                                // stem there and lost the nue selection).
                                if (m_mvga_dup_starved_span > 0 &&
                                    std::min(la, lb) / std::max(la, lb) < m_mvga_dup_starved_span) {
                                    SPDLOG_LOGGER_TRACE(s_log,
                                        "mvga: op1-post starved-decline cluster={} overlap={:.2f} "
                                        "span {:.2f}/{:.2f}cm (not comparable)",
                                        cluster.ident(), frac, std::min(la, lb)/units::cm,
                                        std::max(la, lb)/units::cm);
                                    continue;
                                }
                                const double ra = segment_median_dQ_dx(sa) / m_mip_dqdx_median;
                                const double rb = segment_median_dQ_dx(sb) / m_mip_dqdx_median;
                                if (ra <= 0 || rb <= 0) continue;
                                const double asym = std::min(ra, rb) / std::max(ra, rb);
                                const double lo = std::min(ra, rb);
                                const double hi = std::max(ra, rb);
                                // The single m_mvga_dup_starved_mip threshold
                                // separates the pair BOTH ways: loser at or
                                // below it, survivor at or above it (doc
                                // pr/99 sec 8: the survivor must carry the
                                // charge it claims -- 46363's 0.61-ratio
                                // "survivor" is itself starved, no verdict).
                                if (asym <= m_mvga_dup_starved_asym &&
                                    lo <= m_mvga_dup_starved_mip &&
                                    hi >= m_mvga_dup_starved_mip) {
                                    forced_loser = (ra <= rb) ? sa : sb;
                                }
                                if (!forced_loser) {
                                    SPDLOG_LOGGER_TRACE(s_log,
                                        "mvga: op1-post starved-decline cluster={} overlap={:.2f} "
                                        "angle={:.1f}deg ratios {:.2f}/{:.2f} asym={:.2f}",
                                        cluster.ident(), frac, ang, ra, rb, asym);
                                    continue;
                                }
                                SPDLOG_LOGGER_DEBUG(s_log,
                                    "mvga: op1-post starved-override cluster={} overlap={:.2f} "
                                    "angle={:.1f}deg starved_ratio={:.2f} healthy_ratio={:.2f} "
                                    "starved_len={:.2f}cm",
                                    cluster.ident(), frac, ang,
                                    (forced_loser == sa) ? ra : rb,
                                    (forced_loser == sa) ? rb : ra,
                                    segment_track_length(forced_loser)/units::cm);
                            }
                        }
                    }

                    // Keep the higher-integrated-charge member (tie: keep
                    // longer).  A created (fitless) member integrates 0 and
                    // loses to any fitted one -- the wanted outcome: the
                    // carried duplicate dies, the original survives.
                    // pr/99 round 2: a starved-override pair skips this rule
                    // -- the starved member dies regardless of totals.
                    double qa = segment_integrated_dQ(sa);
                    double qb = segment_integrated_dQ(sb);
                    SegmentPtr loser;
                    if (forced_loser) loser = forced_loser;
                    else if (qa == qb) loser = shorter;
                    else loser = (qa < qb) ? sa : sb;
                    SegmentPtr survivor = (loser == sa) ? sb : sa;

                    auto [lv1, lv2] = find_vertices(graph, loser);
                    auto [sv1, sv2] = find_vertices(graph, survivor);
                    if (!lv1 || !lv2 || !sv1 || !sv2) continue;

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
                    if (!feasible) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op1-post reconnect-infeasible cluster={} loser_len={:.2f}cm",
                            cluster.ident(), segment_track_length(loser)/units::cm);
                        continue;
                    }

                    remove_segment(graph, loser);
                    for (auto& [le, target] : plans) connect_direct(le, target);
                    cleanup_vertex(lv1);
                    cleanup_vertex(lv2);

                    SPDLOG_LOGGER_DEBUG(s_log,
                        "mvga: op1-post dup-merge cluster={} removed seg len={:.2f}cm sumdQ={:.3g} "
                        "overlap={:.2f}@{:.1f}mm vs survivor len={:.2f}cm sumdQ={:.3g} reconnects={}",
                        cluster.ident(), segment_track_length(loser)/units::cm,
                        (loser == sa) ? qa : qb, frac, m_mvga_dup_tol/units::mm,
                        segment_track_length(survivor)/units::cm,
                        (loser == sa) ? qb : qa, plans.size());
                    ++n_op1_post;
                    flag_continue = true;
                }
            }
        }
        if (n_op1_post == post_before) break;  // fixed point: nothing merged this round
        track_fitter.do_multi_tracking(true, true, false, m_fit_exclusion, false, &cluster);
        SPDLOG_LOGGER_DEBUG(s_log,
            "mvga: op1-post round={} cluster={} merges={} (refit done)",
            post_rounds, cluster.ident(), n_op1_post - post_before);
        ++post_rounds;
        }  // post_rounds fixed-point loop
    }
    if (n_op1_post > 0) {
        SPDLOG_LOGGER_DEBUG(s_log,
            "mvga: op1-post fired cluster={} merges={}",
            cluster.ident(), n_op1_post);
    }

    // ---- op1-proj: projective duplicate collapse at the main vertex ------
    // doc pr/83 r4: a 1-track-1-shower stem reported as TWO main-vertex
    // tracks that overlap in the 2D wire views while separating in 3D --
    // the fitter split one projective charge corridor across two 3D
    // interpretations, starving one of charge (stem dQ/dx ratios measured
    // 0.08-0.28 vs >= 0.7 for genuine two-prong vertices).  3D corridor
    // overlap reads 0.14-0.58, below every op1/op1-post gate, so rounds 1-3
    // never fire (138009 12094/12095 o3D=0.58, 168596 14168/14172 o3D=0.14,
    // 74544 12105/12107 o3D=0.46; each overlaps >= 0.75 in 2 of 3 views).
    // Candidates are ONLY segment pairs incident on the main vertex (the
    // shower-stem beginning; 390842-class near-vertex riders are NOT
    // incident and cannot enter).  Runs post-refit for the same reason as
    // op1-post: the dQ/dx starvation is a fit product.  Merge recipe is
    // op1's (keep higher integrated charge, pre-verified reconnects), in
    // its own fixed-point merge->refit loop.  m_mvga_proj_dup_frac == 0
    // (default) => skipped => byte-identical.
    int n_op1_proj = 0;
    if (m_mvga_proj_dup_frac > 0 && m_mvga_dup_tol > 0
        && main_vertex->descriptor_valid()) {
        auto grouping = cluster.grouping();
        // Per-view 2D overlap: fraction of A's fit points within
        // m_mvga_dup_tol of B's in view coords (x, cos(a)z - sin(a)y),
        // wire angles from the pair's (apa,face).  Returns the per-plane
        // fractions sorted descending, or empty on any inconsistency.
        auto view_overlaps = [&](const std::vector<Fit>& fa,
                                 const std::vector<Fit>& fb,
                                 const std::pair<int,int>& paf) {
            std::vector<double> out;
            if (!grouping) return out;
            const auto [au, av, aw] = grouping->wire_angles(paf.first, paf.second);
            for (double ang : {au, av, aw}) {
                const double c = std::cos(ang), s = std::sin(ang);
                int nin = 0, ntot = 0;
                for (const auto& pa : fa) {
                    if (!pa.valid()) continue;
                    ++ntot;
                    const double xa = pa.point.x();
                    const double wa = c*pa.point.z() - s*pa.point.y();
                    double best = 1e30;
                    for (const auto& pb : fb) {
                        if (!pb.valid()) continue;
                        const double dx = pb.point.x() - xa;
                        const double dw = c*pb.point.z() - s*pb.point.y() - wa;
                        best = std::min(best, dx*dx + dw*dw);
                    }
                    if (best < m_mvga_dup_tol*m_mvga_dup_tol) ++nin;
                }
                if (ntot == 0) return std::vector<double>{};
                out.push_back(static_cast<double>(nin)/ntot);
            }
            std::sort(out.rbegin(), out.rend());
            return out;
        };
        // dQ/dx over the first 8 cm of fitted path from the main vertex.
        auto stem_dqdx = [&](const SegmentPtr& sg, const WireCell::Point& mvp) {
            std::vector<const Fit*> fits;
            for (const auto& f : sg->fits())
                if (f.valid() && f.dx > 0 && f.dQ >= 0) fits.push_back(&f);
            if (fits.size() < 2) return -1.0;
            if (point_dis(fits.back()->point, mvp) < point_dis(fits.front()->point, mvp))
                std::reverse(fits.begin(), fits.end());
            double sq = 0, sx = 0;
            for (const Fit* f : fits) {
                sq += f->dQ;
                sx += f->dx;
                if (sx >= 8*units::cm) break;
            }
            return (sx > 0.5*units::cm) ? sq/sx : -1.0;
        };
        int proj_rounds = 0;
        while (proj_rounds < kEditCap && n_op1_proj < kEditCap) {
        const int proj_before = n_op1_proj;
        bool flag_continue = true;
        while (flag_continue && n_op1_proj < kEditCap) {
            flag_continue = false;
            if (!main_vertex->descriptor_valid()) break;
            const WireCell::Point mvp = main_vertex->fit().valid()
                ? main_vertex->fit().point : main_vertex->wcpt().point;
            std::vector<SegmentPtr> segs;
            for (auto edesc : sorted_out_edges(main_vertex->get_descriptor(), graph)) {
                SegmentPtr sg = graph[edesc].segment;
                if (!sg) continue;
                if (std::find(segs.begin(), segs.end(), sg) != segs.end()) continue;
                segs.push_back(sg);
            }
            std::sort(segs.begin(), segs.end(), SegmentIndexCmp{});
            for (size_t i = 0; i + 1 < segs.size() && !flag_continue; ++i) {
                for (size_t j = i + 1; j < segs.size() && !flag_continue; ++j) {
                    SegmentPtr sa = segs[i];
                    SegmentPtr sb = segs[j];
                    // doc pr/103: the connector op0 just established is never a
                    // duplicate to be merged away (405707: with S's stale fits
                    // it reads overlap=1.00).  Empty set when the knob is off.
                    if (passthru_stubs.count(sa) || passthru_stubs.count(sb)) continue;
                    double la = segment_track_length(sa);
                    double lb = segment_track_length(sb);
                    SegmentPtr shorter = (la <= lb) ? sa : sb;
                    SegmentPtr longer  = (la <= lb) ? sb : sa;
                    auto pts_s = seg_points(shorter);
                    auto pts_l = seg_points(longer);
                    if (pts_s.size() <= static_cast<size_t>(m_mvga_stub_pts)) continue;
                    // dQ/dx needs fits on both members.
                    if (seg_valid_fits(shorter) < 2 || seg_valid_fits(longer) < 2) continue;
                    // Near-parallel gate first (cheap): op1's chord test,
                    // with op1-proj's own ceiling when set (doc pr/83 r4b,
                    // 284206: residual pair at 22 deg; 0 = op1's shared
                    // m_mvga_dup_angle => byte-identical).
                    const double proj_angle =
                        (m_mvga_proj_angle > 0) ? m_mvga_proj_angle : m_mvga_dup_angle;
                    if (proj_angle > 0 && pts_s.size() >= 2 && pts_l.size() >= 2) {
                        auto ca = pts_s.back() - pts_s.front();
                        auto cb = pts_l.back() - pts_l.front();
                        double den = ca.magnitude() * cb.magnitude();
                        if (den > 0) {
                            double cosang = std::abs(ca.dot(cb)) / den;
                            double ang = std::acos(std::clamp(cosang, 0.0, 1.0)) / 3.1415926 * 180.0;
                            if (ang > proj_angle) continue;
                        }
                    }
                    // Same (apa,face) required -- a cross-face pair has no
                    // single view frame to compare in.
                    std::pair<int,int> paf{-1,-1};
                    for (const auto& f : longer->fits())
                        if (f.valid() && f.paf.first >= 0) { paf = f.paf; break; }
                    bool same_paf = (paf.first >= 0);
                    for (const auto& f : shorter->fits())
                        if (f.valid() && f.paf.first >= 0 && f.paf != paf) { same_paf = false; break; }
                    if (!same_paf) continue;
                    auto ov = view_overlaps(shorter->fits(), longer->fits(), paf);
                    if (ov.size() < 3) continue;
                    if (ov[1] > 0.25) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op1-proj eval cluster={} pair len {:.2f}/{:.2f}cm views={:.2f}/{:.2f}/{:.2f}",
                            cluster.ident(), la/units::cm, lb/units::cm, ov[0], ov[1], ov[2]);
                    }
                    // 2-of-3 views must read duplicate.
                    if (ov[1] < m_mvga_proj_dup_frac) continue;
                    // Stem dQ/dx asymmetry: the projective ghost is charge-
                    // starved; a genuine collinear two-prong is not.
                    const double da = stem_dqdx(sa, mvp);
                    const double db = stem_dqdx(sb, mvp);
                    if (da <= 0 || db <= 0) continue;
                    const double ratio = std::min(da, db) / std::max(da, db);
                    SPDLOG_LOGGER_TRACE(s_log,
                        "mvga: op1-proj dqdx cluster={} stem {:.3g}/{:.3g} ratio={:.2f} gate={}",
                        cluster.ident(), da, db, ratio,
                        (ratio < m_mvga_proj_dqdx_ratio) ? "pass" : "decline");
                    if (ratio >= m_mvga_proj_dqdx_ratio) continue;

                    // Keep the higher-integrated-charge member (op1 rule).
                    double qa = segment_integrated_dQ(sa);
                    double qb = segment_integrated_dQ(sb);
                    SegmentPtr loser;
                    if (qa == qb) loser = shorter;
                    else loser = (qa < qb) ? sa : sb;
                    SegmentPtr survivor = (loser == sa) ? sb : sa;

                    auto [lv1, lv2] = find_vertices(graph, loser);
                    auto [sv1, sv2] = find_vertices(graph, survivor);
                    if (!lv1 || !lv2 || !sv1 || !sv2) continue;

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
                    if (!feasible) {
                        SPDLOG_LOGGER_TRACE(s_log,
                            "mvga: op1-proj reconnect-infeasible cluster={} loser_len={:.2f}cm",
                            cluster.ident(), segment_track_length(loser)/units::cm);
                        continue;
                    }

                    remove_segment(graph, loser);
                    for (auto& [le, target] : plans) connect_direct(le, target);
                    cleanup_vertex(lv1);
                    cleanup_vertex(lv2);

                    SPDLOG_LOGGER_DEBUG(s_log,
                        "mvga: op1-proj dup-merge cluster={} removed seg len={:.2f}cm sumdQ={:.3g} "
                        "views={:.2f}/{:.2f}/{:.2f}@{:.1f}mm dqdx_ratio={:.2f} "
                        "vs survivor len={:.2f}cm sumdQ={:.3g} reconnects={}",
                        cluster.ident(), segment_track_length(loser)/units::cm,
                        (loser == sa) ? qa : qb, ov[0], ov[1], ov[2],
                        m_mvga_dup_tol/units::mm, ratio,
                        segment_track_length(survivor)/units::cm,
                        (loser == sa) ? qb : qa, plans.size());
                    ++n_op1_proj;
                    flag_continue = true;
                }
            }
        }
        if (n_op1_proj == proj_before) break;  // fixed point
        track_fitter.do_multi_tracking(true, true, false, m_fit_exclusion, false, &cluster);
        SPDLOG_LOGGER_DEBUG(s_log,
            "mvga: op1-proj round={} cluster={} merges={} (refit done)",
            proj_rounds, cluster.ident(), n_op1_proj - proj_before);
        ++proj_rounds;
        }  // proj_rounds fixed-point loop
    }
    if (n_op1_proj > 0) {
        SPDLOG_LOGGER_DEBUG(s_log,
            "mvga: op1-proj fired cluster={} merges={}",
            cluster.ident(), n_op1_proj);
    }
    return fired_ops || n_op1_post > 0 || n_op1_proj > 0;
}

// doc pr/83 r3 (sec 9.5, Mechanism C): one unscoped duplicate-corridor pass
// over a cluster being ABANDONED by swap_main_cluster.  Fork-by-duplication
// of op1's metric and merge recipe (main_vertex_graph_audit above stays
// byte-untouched): no radius (the cluster has no main vertex to center one
// on -- that absence is the entire defect), no `created` exemption (nothing
// here is mid-pass state), same overlap/angle gates via the op1-effective
// members, same keep-the-higher-charge merge, same pre-verified reconnects,
// one per-cluster refit if anything merged so the survivor geometry is what
// the bundle/nusel layer ultimately reports.  Caller gates on
// m_swap_orphan_dup_audit; this function itself is knob-free.
bool PatternAlgorithms::orphan_dup_audit(Graph& graph, Facade::Cluster& cluster,
                                         TrackFitting& track_fitter,
                                         IDetectorVolumes::pointer dv)
{
    if (!(m_mvga_dup_frac > 0 && m_mvga_dup_tol > 0)) return false;
    const double dup_frac = (m_mvga_op1_dup_frac > 0) ? m_mvga_op1_dup_frac : m_mvga_dup_frac;
    constexpr int kEditCap = 8;  // same paranoia cap as mvga

    auto cluster_segments = [&]() {
        std::vector<SegmentPtr> segs;
        for (const auto& vd : ordered_nodes(graph)) {
            VertexPtr vtx = graph[vd].vertex;
            if (!vtx || vtx->cluster() != &cluster) continue;
            for (auto edesc : sorted_out_edges(vd, graph)) {
                SegmentPtr sg = graph[edesc].segment;
                if (!sg) continue;
                if (std::find(segs.begin(), segs.end(), sg) != segs.end()) continue;
                segs.push_back(sg);
            }
        }
        std::sort(segs.begin(), segs.end(), SegmentIndexCmp{});
        return segs;
    };

    // Rough-path edge v1 -> v2 unless one already exists (op1's
    // connect_direct without the `created` bookkeeping -- no later op reads
    // it here).
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
        return sg;
    };

    // Drop a now-degree-0 vertex (no main vertex exists on an abandoned
    // cluster; protected breaks still survive -- pr/48 / pr/50 precedent).
    auto cleanup_vertex = [&](VertexPtr vtx) {
        if (!vtx) return;
        if (vtx->flags_any(VertexFlags::kProtectedBreak)) return;
        if (!vtx->descriptor_valid()) return;
        if (boost::degree(vtx->get_descriptor(), graph) == 0) {
            remove_vertex(graph, vtx);
        }
    };

    int n_merged = 0;
    bool flag_continue = true;
    while (flag_continue && n_merged < kEditCap) {
        flag_continue = false;
        auto segs = cluster_segments();
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
                if (pts_s.size() <= static_cast<size_t>(m_mvga_stub_pts)) continue;
                double frac = path_overlap_fraction(pts_s, pts_l, m_mvga_dup_tol);
                if (frac > 0.25) {
                    SPDLOG_LOGGER_TRACE(s_log,
                        "mvga: orphan-dup eval cluster={} pair len {:.2f}/{:.2f}cm npts {}/{} overlap={:.2f}",
                        cluster.ident(), segment_track_length(shorter)/units::cm,
                        segment_track_length(longer)/units::cm,
                        pts_s.size(), pts_l.size(), frac);
                }
                if (frac < dup_frac) continue;
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

                double qa = segment_integrated_dQ(sa);
                double qb = segment_integrated_dQ(sb);
                SegmentPtr loser;
                if (qa == qb) loser = shorter;
                else loser = (qa < qb) ? sa : sb;
                SegmentPtr survivor = (loser == sa) ? sb : sa;

                auto [lv1, lv2] = find_vertices(graph, loser);
                auto [sv1, sv2] = find_vertices(graph, survivor);
                if (!lv1 || !lv2 || !sv1 || !sv2) continue;

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
                if (!feasible) {
                    SPDLOG_LOGGER_TRACE(s_log,
                        "mvga: orphan-dup reconnect-infeasible cluster={} loser_len={:.2f}cm",
                        cluster.ident(), segment_track_length(loser)/units::cm);
                    continue;
                }

                remove_segment(graph, loser);
                for (auto& [le, target] : plans) connect_direct(le, target);
                cleanup_vertex(lv1);
                cleanup_vertex(lv2);

                SPDLOG_LOGGER_DEBUG(s_log,
                    "mvga: orphan-dup merge cluster={} removed seg len={:.2f}cm sumdQ={:.3g} "
                    "overlap={:.2f}@{:.1f}mm vs survivor len={:.2f}cm sumdQ={:.3g} reconnects={}",
                    cluster.ident(), segment_track_length(loser)/units::cm,
                    (loser == sa) ? qa : qb, frac, m_mvga_dup_tol/units::mm,
                    segment_track_length(survivor)/units::cm,
                    (loser == sa) ? qb : qa, plans.size());
                ++n_merged;
                flag_continue = true;
            }
        }
    }

    if (n_merged > 0) {
        track_fitter.do_multi_tracking(true, true, false, m_fit_exclusion, false, &cluster);
        SPDLOG_LOGGER_DEBUG(s_log,
            "mvga: orphan-dup fired cluster={} merges={} (refit done)",
            cluster.ident(), n_merged);
    }
    return n_merged > 0;
}

// doc sbnd_xin/docs/pr/84 round 2 (F3 = pr/84 P2 "conn3_stitch").  A main
// cluster is one contiguous lump of charge, yet its segment graph can be
// disconnected (pr/54 keep-isolated residuals arrive as -- the source
// comment's own words -- "a disconnected piece of this cluster's graph";
// snap_main_vertex_to_kink can strand the main vertex on a tiny component,
// SBND evt 283713).  Downstream, shower_conn3_unreachable promotes the
// unreachable pieces to conn-3 "association" showers while logging anchor
// distances of millimetres, and the Bee PF writer then hangs them under
// synthetic gamma/neutron carriers.  This pass runs AFTER mvga and BEFORE
// clustering_points: bridge each component whose closest approach to the
// reachable side is within m_conn3_stitch_max with a real rough-path
// segment, then refit once, so every later pass sees a connected graph and
// the piece is classified conn-1 naturally.  Wider gaps still fall through
// to the conn3_unreachable backstop.
bool PatternAlgorithms::stitch_disconnected_main_cluster(Graph& graph, Facade::Cluster& cluster,
                                                         VertexPtr main_vertex,
                                                         TrackFitting& track_fitter,
                                                         IDetectorVolumes::pointer dv)
{
    if (m_conn3_stitch_max <= 0) return false;
    if (!main_vertex || !main_vertex->descriptor_valid()) return false;

    constexpr int kEditCap = 8;  // same paranoia cap as mvga / orphan_dup_audit

    // All in-cluster segments in stable index order (orphan_dup_audit recipe).
    auto cluster_segments = [&]() {
        std::vector<SegmentPtr> segs;
        for (const auto& vd : ordered_nodes(graph)) {
            VertexPtr vtx = graph[vd].vertex;
            if (!vtx || vtx->cluster() != &cluster) continue;
            for (auto edesc : sorted_out_edges(vd, graph)) {
                SegmentPtr sg = graph[edesc].segment;
                if (!sg) continue;
                if (std::find(segs.begin(), segs.end(), sg) != segs.end()) continue;
                segs.push_back(sg);
            }
        }
        std::sort(segs.begin(), segs.end(), SegmentIndexCmp{});
        return segs;
    };

    // Rough-path edge v1 -> v2 unless one already exists (orphan_dup_audit's
    // connect_direct; no `created` bookkeeping -- nothing later reads it).
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
        return sg;
    };

    int n_stitched = 0;
    bool flag_continue = true;
    while (flag_continue && n_stitched < kEditCap) {
        flag_continue = false;

        const auto reachable = reachable_without(graph, cluster, main_vertex, nullptr);

        // Global argmin over (unreachable segment, reachable vertex) of the
        // segment's closest fitted approach to the vertex -- the same metric
        // shower_conn3_unreachable logs as anchor_dis.  Stable iteration
        // (segments by index, vertices by node order) + strict < keeps the
        // argmin deterministic.
        double best_dis = 1e9;
        SegmentPtr best_seg = nullptr;
        VertexPtr best_vtx = nullptr;
        for (const auto& sg : cluster_segments()) {
            auto [sv1, sv2] = find_vertices(graph, sg);
            if (!sv1 || !sv2) continue;
            // An edge's two endpoints share a component: testing one suffices.
            if (!sv1->descriptor_valid() || reachable.count(sv1->get_descriptor())) continue;
            for (const auto& vd : ordered_nodes(graph)) {
                if (!reachable.count(vd)) continue;
                VertexPtr vtx = graph[vd].vertex;
                if (!vtx || vtx->cluster() != &cluster) continue;
                WireCell::Point vp = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
                const double dis = segment_get_closest_point(sg, vp).first;
                if (dis >= 0 && dis < best_dis) {
                    best_dis = dis;
                    best_seg = sg;
                    best_vtx = vtx;
                }
            }
        }

        if (!best_seg || best_dis > m_conn3_stitch_max) break;

        // Component-side anchor: the bridged segment's endpoint vertex nearer
        // the chosen reachable vertex.  (The closest approach may be
        // mid-segment; vertex-to-vertex bridging is the accepted
        // approximation at these <~3 cm gaps -- breaking the segment is out
        // of scope, recorded in the round doc.)
        auto [cv1, cv2] = find_vertices(graph, best_seg);
        if (!cv1 || !cv2) break;
        WireCell::Point rp = best_vtx->fit().valid() ? best_vtx->fit().point : best_vtx->wcpt().point;
        auto vtx_pt = [](const VertexPtr& v) {
            return v->fit().valid() ? v->fit().point : v->wcpt().point;
        };
        VertexPtr comp_vtx = (point_dis(vtx_pt(cv1), rp) <= point_dis(vtx_pt(cv2), rp)) ? cv1 : cv2;

        SegmentPtr bridge = connect_direct(best_vtx, comp_vtx);
        if (!bridge) {
            // do_rough_path could not cross the gap (genuinely disconnected
            // charge): leave the component to the conn3_unreachable backstop.
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr84 conn3_stitch: decline cluster={} gap={:.2f}cm (no rough path)",
                cluster.ident(), best_dis/units::cm);
            break;
        }

        SPDLOG_LOGGER_DEBUG(s_log,
            "pr84 conn3_stitch: bridge cluster={} gap={:.2f}cm seg_len={:.2f}cm comp_seg_len={:.2f}cm",
            cluster.ident(), best_dis/units::cm,
            segment_track_length(bridge)/units::cm,
            segment_track_length(best_seg)/units::cm);
        ++n_stitched;
        flag_continue = true;
    }

    if (n_stitched > 0) {
        track_fitter.do_multi_tracking(true, true, false, m_fit_exclusion, false, &cluster);
        SPDLOG_LOGGER_DEBUG(s_log,
            "pr84 conn3_stitch: fired cluster={} bridges={} (refit done)",
            cluster.ident(), n_stitched);
    }
    return n_stitched > 0;
}
