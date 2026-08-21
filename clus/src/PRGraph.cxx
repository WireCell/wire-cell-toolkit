#include "WireCellClus/PRGraph.h"

#include <cstdio>
#include <cstdlib>
#include <deque>
#include <execinfo.h>

namespace WireCell::Clus::PR {

    // doc pr/30 §11 instrumentation (see PRGraph.h).
    PortAuditCounters   g_port_audit;
    GraphEndpointPolicy g_graph_endpoint_policy;
    // doc pr/36 §10 tagger-stage instrumentation (see PRGraph.h).
    Pr36AuditCounters   g_pr36_audit;
    // doc pr/33 §10 EM-shower-clustering instrumentation (see PRGraph.h).
    Pr33AuditCounters   g_pr33_audit;

    // Determinism-debug (WCT_DET_DEBUG=2): log every segment creation with
    // its stable index and a mini backtrace so creation-order divergence
    // between two runs can be diffed and attributed to the calling function.
    static bool det_dbg_graph() {
        static const char* v = std::getenv("WCT_DET_DEBUG");
        static const bool on = (v && v[0] == '2');
        return on;
    }

    bool add_vertex(Graph& g, VertexPtr vtx)
    {
        if (vtx->descriptor_valid()) {
            return false;
        }
        auto& gb = g[boost::graph_bundle];
        const size_t index = gb.num_node_indices;
        auto desc = boost::add_vertex(NodeBundle{vtx, index}, g);
        ++ gb.num_node_indices;
        vtx->set_descriptor(desc);
        vtx->set_graph_index(index);
        return true;
    }

    bool remove_vertex(Graph& graph, VertexPtr vtx)
    {
        if (! vtx->descriptor_valid()) { return false; }
        auto desc = vtx->get_descriptor();
        boost::remove_vertex(desc, graph);
        vtx->invalidate_descriptor();
        return true;
    }

    // doc pr/30 §11, P8 -- the endpoint-consistency check the port dropped.
    //
    // Prototype NeutrinoID::add_proto_connection refuses the connection and
    // prints "Error! Vertex and Segment does not match" unless the vertex's
    // wcpt INDEX equals the segment's front or back index
    // (NeutrinoID_proto_vertex.h:1952-1956).  The toolkit carries no Steiner
    // point-cloud index on a PR::Vertex, so the check is restated
    // positionally -- which is the owner-accepted substitution (doc pr/30
    // §10, clarification 1) -- rather than left out, which is what it was.
    //
    // The tolerance is a stand-in for index equality, not a physics cut: a
    // consistent vertex holds a COPY of the segment's terminal wcpt, so the
    // distance is exactly 0 in the normal case and the tolerance only absorbs
    // a later fit nudge.  0.3 cm is far under Steiner point spacing (it
    // cannot merge two distinct ends) and far over FP noise.
    //
    // Returns true when the pair is consistent.  Counting is unconditional;
    // only the REFUSAL is knobbed.
    //
    // WHAT THIS ACTUALLY CATCHES ON SBND -- read before turning strict on.
    // 108 firings over 22 of the 48 nueCC events (doc pr/30 §11.4), and every
    // one traced with WCT_DET_DEBUG=2 comes from PatternAlgorithms::
    // crawl_segment (NeutrinoStructureExaminer.cxx:892 and :944).  That
    // function rebuilds each of the vertex's segments so the new path
    // TERMINATES at vtx_new_point, attaches them, and only then assigns
    // `vertex->wcpt().point = vtx_new_point` (:948).  So the inconsistency is
    // TRANSIENT and self-healing: by the time crawl_segment returns, every
    // one of those connections is consistent.  It is not a lost invariant.
    // Consequently m_graph_endpoint_strict is a FALSE POSITIVE as placed --
    // it refuses legitimate connections mid-repair, and measurably damages
    // the result (22/48 events change, 5 nue candidates lost vs 1 gained).
    // A check for a PERSISTENT violation would have to run at end-of-stage,
    // not here.  This one is retained as a tripwire, not as a fix.
    static bool endpoint_consistent(const SegmentPtr& seg,
                                    const VertexPtr& vtx1, const VertexPtr& vtx2)
    {
        if (seg->wcpts().empty()) { return true; }   // nothing to check against
        const auto& front = seg->wcpts().front().point;
        const auto& back  = seg->wcpts().back().point;
        const double tol  = g_graph_endpoint_policy.tol;

        auto at_an_end = [&](const VertexPtr& v) {
            if (!v) { return true; }                 // null handled by callers
            const auto& p = v->wcpt().point;
            return ray_length(Ray{p, front}) <= tol || ray_length(Ray{p, back}) <= tol;
        };
        // The prototype tests each connection separately (one call per
        // vertex), so a pair is consistent exactly when BOTH ends are.
        return at_an_end(vtx1) && at_an_end(vtx2);
    }

    bool add_segment(Graph& g, SegmentPtr seg, VertexPtr vtx1, VertexPtr vtx2)
    {
        g_port_audit.add_segment_calls.fetch_add(1, std::memory_order_relaxed);

        // The check sits AFTER the descriptor_valid() early return, because the
        // prototype's counterpart guards add_proto_connection -- a connection
        // being MADE -- not a call that turns out to be a no-op.  On the 48
        // nueCC events this placement is measured to make no difference
        // (add_segment_reentry == 0 in every event: no call ever arrives with
        // the segment already in the graph), but the ordering is the correct
        // one on the argument and the doctest pins it.
        if (seg && seg->descriptor_valid()) {
            g_port_audit.add_segment_reentry.fetch_add(1, std::memory_order_relaxed);
        }

        bool changed = false;
        changed = add_vertex(g, vtx1) || changed;
        changed = add_vertex(g, vtx2) || changed;

        if (seg->descriptor_valid()) {
            return changed;
        }

        if (seg && vtx1 && vtx2 && !endpoint_consistent(seg, vtx1, vtx2)) {
            const auto n = g_port_audit.endpoint_mismatch.fetch_add(1, std::memory_order_relaxed);
            if (n < 20) {   // bound the log; the counter carries the full rate
                auto log = Log::logger("clus");
                // DEBUG, not WARN (doc pr/30 §12.10).  This fires ~108 times
                // per 48 SBND events and every one of those is the benign
                // crawl_segment transient documented above.  At WARN it is
                // noise that trains the reader to ignore the one channel that
                // would matter if it ever fired from another call site.
                const auto& f = seg->wcpts().front().point;
                const auto& b = seg->wcpts().back().point;
                SPDLOG_LOGGER_DEBUG(log,
                    "PR::add_segment: vertex/segment endpoint mismatch (doc pr/30 P8): "
                    "seg nwcpts={} front=({:.3f},{:.3f},{:.3f}) back=({:.3f},{:.3f},{:.3f}) "
                    "v1=({:.3f},{:.3f},{:.3f}) v2=({:.3f},{:.3f},{:.3f}) tol={:.3f} cm strict={}",
                    seg->wcpts().size(),
                    f.x()/units::cm, f.y()/units::cm, f.z()/units::cm,
                    b.x()/units::cm, b.y()/units::cm, b.z()/units::cm,
                    vtx1->wcpt().point.x()/units::cm, vtx1->wcpt().point.y()/units::cm,
                    vtx1->wcpt().point.z()/units::cm,
                    vtx2->wcpt().point.x()/units::cm, vtx2->wcpt().point.y()/units::cm,
                    vtx2->wcpt().point.z()/units::cm,
                    g_graph_endpoint_policy.tol/units::cm,
                    g_graph_endpoint_policy.strict);
            }
            if (g_graph_endpoint_policy.strict) {
                // Prototype behaviour: the connection is NOT made.  Most
                // prototype callers ignore add_proto_connection's return
                // value, so this reproduces the effective outcome as well as
                // the nominal one.  Default OFF -- this changes the graph.
                // NOTE the two vertices have already been added above, exactly
                // as the prototype leaves any vertex it had already recorded;
                // only the connection is withheld.
                g_port_audit.endpoint_refused.fetch_add(1, std::memory_order_relaxed);
                return changed;
            }
        }

        // doc sbnd_xin/docs/pr/31 §10.10 (F9): a self-loop connection about to
        // be made.  boost::add_edge accepts it and boost::degree then counts
        // it TWICE, where the prototype's map_vertex_segments set counts one
        // segment -- so every start_n/end_n "== 1" free-end test downstream
        // scores such a segment differently in the two trees.  Counted, not
        // fixed: if this never fires the divergence is vacuous.
        if (vtx1 == vtx2) {
            const auto n = g_port_audit.selfloop_segment.fetch_add(1, std::memory_order_relaxed);
            if (n < 20) {
                auto log = Log::logger("clus");
                SPDLOG_LOGGER_DEBUG(log,
                    "PR::add_segment: self-loop connection (doc pr/31 F9): "
                    "seg nwcpts={} vtx=({:.3f},{:.3f},{:.3f})",
                    seg ? seg->wcpts().size() : 0,
                    vtx1->wcpt().point.x()/units::cm, vtx1->wcpt().point.y()/units::cm,
                    vtx1->wcpt().point.z()/units::cm);
            }
        }

        auto& gb = g[boost::graph_bundle];
        const size_t index = gb.num_edge_indices;
        auto [desc,added] = boost::add_edge(vtx1->get_descriptor(),
                                            vtx2->get_descriptor(),
                                            EdgeBundle{seg, index}, g);

        seg->set_descriptor(desc);

        // Edge was added
        if (added) {
            ++ gb.num_edge_indices;
            seg->set_descriptor(desc);
            seg->set_graph_index(index);
            if (det_dbg_graph()) {
                const auto& p1 = vtx1->wcpt().point;
                const auto& p2 = vtx2->wcpt().point;
                double wf[3] = {0,0,0}, wb[3] = {0,0,0};
                if (!seg->wcpts().empty()) {
                    const auto& f = seg->wcpts().front().point;
                    const auto& b = seg->wcpts().back().point;
                    wf[0]=f.x(); wf[1]=f.y(); wf[2]=f.z();
                    wb[0]=b.x(); wb[1]=b.y(); wb[2]=b.z();
                }
                fprintf(stderr, "WCT_DETA seg idx=%zu nw=%zu v1=(%.4f,%.4f,%.4f) v2=(%.4f,%.4f,%.4f) wf=(%.4f,%.4f,%.4f) wb=(%.4f,%.4f,%.4f)\n",
                        index, seg->wcpts().size(), p1.x(), p1.y(), p1.z(), p2.x(), p2.y(), p2.z(),
                        wf[0], wf[1], wf[2], wb[0], wb[1], wb[2]);
                void* bt[8];
                int nbt = backtrace(bt, 8);
                char** syms = backtrace_symbols(bt, nbt);
                for (int i = 1; i < nbt && i < 6; ++i)
                    fprintf(stderr, "WCT_DETA   bt[%d] %s\n", i, syms ? syms[i] : "?");
                free(syms);
            }
            return true;
        }

        // Edge already existed, assure its object is this one.
        // Inherit the existing edge's graph index so m_graph_index is not left at SIZE_MAX.
        // doc sbnd_xin/docs/pr/31 §10.10 (F9): when a DIFFERENT segment already
        // owns this edge, the reassignment below orphans it from graph
        // iteration while it still holds a descriptor it believes valid; the
        // prototype's map_vertex_segments holds BOTH segments (a 2-cycle).
        // Counted before the overwrite, not fixed -- relaxing the graph type
        // (multisetS) or rejecting the call is its own escalation.
        if (seg && g[desc].segment && g[desc].segment != seg) {
            const auto n = g_port_audit.edge_aliased.fetch_add(1, std::memory_order_relaxed);
            if (n < 20) {
                auto log = Log::logger("clus");
                SPDLOG_LOGGER_DEBUG(log,
                    "PR::add_segment: parallel segment aliased onto existing edge (doc pr/31 F9): "
                    "edge idx={} old nwcpts={} new nwcpts={}",
                    g[desc].index, g[desc].segment->wcpts().size(), seg->wcpts().size());
            }
        }
        g[desc].segment = seg;
        seg->set_graph_index(g[desc].index);

        return changed;
    }

    // doc pr/67 round 2, P6.  See the declaration in PRGraph.h for why this is a
    // file-static mirror rather than a PatternAlgorithms member.
    static bool s_traj_cover_probe = false;
    void set_traj_cover_probe(bool enable) { s_traj_cover_probe = enable; }

    // doc sbnd_xin/docs/pr/96: WCT_PR96_REMSEG_DEBUG -- which CALLER deletes a
    // fitted segment.  pr/67's P6 sentinel above locates a removal
    // geometrically but not by call site, and on SBND 18255-279955 the branch
    // that covers the owner's uncovered charge is created and then removed
    // twice, so the call site is the whole question.  Same mini-backtrace
    // idiom as the WCT_DET_DEBUG=2 creation probe above.  Log-only and unset by
    // default => byte-identical.
    static bool remseg_dbg() {
        static const bool on = std::getenv("WCT_PR96_REMSEG_DEBUG") != nullptr;
        return on;
    }

    bool remove_segment(Graph& graph, SegmentPtr seg)
    {
        if (! seg->descriptor_valid()) { return false; }
        auto desc = seg->get_descriptor();
        if (s_traj_cover_probe) {
            // Owner hypothesis (b) at the SEGMENT level.  Round 1 tested only
            // point-level trimming (examine_end_ps_vec, P3); a whole fitted
            // segment deleted by one of remove_segment's many callers would not
            // have shown up there.  Endpoints locate the removal geometrically,
            // which is enough to match it against a hand-scanned coordinate.
            static auto s_log = WireCell::Log::logger("clus.PRGraph");
            const auto& fits = seg->fits();
            SPDLOG_LOGGER_DEBUG(s_log,
                "pr67 remove_segment: seg={} nfits={} front=({:.2f},{:.2f},{:.2f}) "
                "back=({:.2f},{:.2f},{:.2f}) cm",
                seg->id(), fits.size(),
                fits.empty() ? 0.0 : fits.front().point.x() / units::cm,
                fits.empty() ? 0.0 : fits.front().point.y() / units::cm,
                fits.empty() ? 0.0 : fits.front().point.z() / units::cm,
                fits.empty() ? 0.0 : fits.back().point.x() / units::cm,
                fits.empty() ? 0.0 : fits.back().point.y() / units::cm,
                fits.empty() ? 0.0 : fits.back().point.z() / units::cm);
        }
        if (remseg_dbg()) {
            const auto& fits = seg->fits();
            fprintf(stderr, "PR96REMSEG nfits=%zu front=(%.2f,%.2f,%.2f) back=(%.2f,%.2f,%.2f)\n",
                    fits.size(),
                    fits.empty() ? 0.0 : fits.front().point.x() / units::cm,
                    fits.empty() ? 0.0 : fits.front().point.y() / units::cm,
                    fits.empty() ? 0.0 : fits.front().point.z() / units::cm,
                    fits.empty() ? 0.0 : fits.back().point.x() / units::cm,
                    fits.empty() ? 0.0 : fits.back().point.y() / units::cm,
                    fits.empty() ? 0.0 : fits.back().point.z() / units::cm);
            void* bt[12];
            int nbt = backtrace(bt, 12);
            char** syms = backtrace_symbols(bt, nbt);
            for (int i = 1; i < nbt && i < 9; ++i)
                fprintf(stderr, "PR96REMSEG   bt[%d] %s\n", i, syms ? syms[i] : "?");
            free(syms);
        }
        boost::remove_edge(desc, graph);
        seg->invalidate_descriptor();
        return true;
    }


    std::pair<VertexPtr, VertexPtr> find_vertices(Graph& graph, SegmentPtr seg)
    {
        if (! seg->descriptor_valid()) { return std::pair<VertexPtr, VertexPtr>{}; }

        auto ed = seg->get_descriptor();

        auto vd1 = boost::source(ed, graph);
        auto vd2 = boost::target(ed, graph);

        auto vtx1 = graph[vd1].vertex;
        auto vtx2 = graph[vd2].vertex;

        // Every caller defends against null vertices in the returned pair
        // (`if (!v1 || !v2) ...`), so a null graph vertex is an anticipated
        // state -- but the wcpt() ordering distance below dereferenced them.
        if (!vtx1 || !vtx2) {
            return std::make_pair(vtx1, vtx2);
        }

        // Use wcpts if available, fall back to fits, then return unordered pair.
        WireCell::Point ept;
        if (!seg->wcpts().empty()) {
            ept = seg->wcpts().front().point;
        } else if (!seg->fits().empty()) {
            ept = seg->fits().front().point;
        } else {
            return std::make_pair(vtx1, vtx2);
        }

        double d1 = ray_length(Ray{vtx1->wcpt().point, ept});
        double d2 = ray_length(Ray{vtx2->wcpt().point, ept});

        if (d1 < d2) {
            return std::make_pair(vtx1, vtx2);
        }
        return std::make_pair(vtx2, vtx1);        
    }

    VertexPtr find_other_vertex(Graph& graph, SegmentPtr seg, VertexPtr vertex)
    {
        if (! seg->descriptor_valid()) { return nullptr; }
        if (! vertex->descriptor_valid()) { return nullptr; }

        auto ed = seg->get_descriptor();
        auto vd = vertex->get_descriptor();

        auto vd1 = boost::source(ed, graph);
        auto vd2 = boost::target(ed, graph);

        auto [ed2,ingraph] = boost::edge(vd1, vd2, graph);
        if (!ingraph)  { return nullptr; }

        // Check if the given vertex is connected to this segment
        if (vd == vd1) {
            return graph[vd2].vertex;
        }
        else if (vd == vd2) {
            return graph[vd1].vertex;
        }

        // Given vertex is not connected to this segment
        return nullptr;
    }

    SegmentPtr find_segment(Graph& graph, VertexPtr vtx1, VertexPtr vtx2)
    {
        if (! vtx1->descriptor_valid()) { return nullptr; }
        if (! vtx2->descriptor_valid()) { return nullptr; }

        auto vd1 = vtx1->get_descriptor();
        auto vd2 = vtx2->get_descriptor();

        // Check if edge exists between the two vertices
        auto [ed, exists] = boost::edge(vd1, vd2, graph);
        if (!exists) { return nullptr; }

        // Return the segment associated with this edge
        return graph[ed].segment;
    }

    IndexedSegmentSet unreachable_segments(const Graph& graph, VertexPtr root)
    {
        // Membership-only BFS over node descriptors: no iteration order can
        // leak into the result, so raw boost::out_edges is safe here.
        std::set<Graph::vertex_descriptor> seen;
        if (root && root->descriptor_valid()) {
            std::deque<Graph::vertex_descriptor> todo{root->get_descriptor()};
            seen.insert(root->get_descriptor());
            while (!todo.empty()) {
                auto vd = todo.front();
                todo.pop_front();
                for (auto [ei, eend] = boost::out_edges(vd, graph); ei != eend; ++ei) {
                    auto nd = boost::target(*ei, graph);
                    if (seen.insert(nd).second) todo.push_back(nd);
                }
            }
        }

        IndexedSegmentSet unreachable;
        for (const auto& ed : ordered_edges(graph)) {
            // Undirected graph: an edge is reachable iff either endpoint is.
            if (seen.count(boost::source(ed, graph)) || seen.count(boost::target(ed, graph))) {
                continue;
            }
            SegmentPtr seg = graph[ed].segment;
            if (seg) unreachable.insert(seg);
        }
        return unreachable;
    }

    IndexedVertexSet reachable_vertices(const Graph& graph, VertexPtr root)
    {
        // Membership-only BFS over node descriptors -- same construction as
        // unreachable_segments() above, which is deliberately duplicated
        // rather than factored out (production path stays byte-for-byte).
        std::set<Graph::vertex_descriptor> seen;
        if (root && root->descriptor_valid()) {
            std::deque<Graph::vertex_descriptor> todo{root->get_descriptor()};
            seen.insert(root->get_descriptor());
            while (!todo.empty()) {
                auto vd = todo.front();
                todo.pop_front();
                for (auto [ei, eend] = boost::out_edges(vd, graph); ei != eend; ++ei) {
                    auto nd = boost::target(*ei, graph);
                    if (seen.insert(nd).second) todo.push_back(nd);
                }
            }
        }

        IndexedVertexSet reachable;
        for (auto nd : ordered_nodes(graph)) {
            if (!seen.count(nd)) continue;
            if (VertexPtr vtx = graph[nd].vertex) reachable.insert(vtx);
        }
        return reachable;
    }

}

    
