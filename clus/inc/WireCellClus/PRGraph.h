/** Define functions that can operate on a "trajectory" graph.

    Functions provided should be preferred over equivalent graph related
    functions in `boost::` as the trajectory graph requires properly bookkeeping
    when adding/removing nodes and/or edges.  It is undefined behavior if you
    call `boost::` graph mutators and fail to follow the bookkeeping
    conventions.

 */
#ifndef WIRECELL_CLUS_PR_GRAPH
#define WIRECELL_CLUS_PR_GRAPH

#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRSegment.h"

#include "WireCellUtil/Graph.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Units.h"
#include <atomic>
#include <cstdint>
#include <limits>
#include <set>

namespace WireCell::Clus::PR {

    /// doc sbnd_xin/docs/pr/30 §11 -- port-fidelity audit instrumentation.
    ///
    /// Process-wide counters.  Every SBND runner drives ONE event per
    /// wire-cell process, so a cumulative total read at the end of
    /// TaggerCheckNeutrino::visit() IS the per-event total; nothing here
    /// needs resetting.  Atomic because clusters may be visited from more
    /// than one thread.
    ///
    /// These are pure diagnostics: incrementing them cannot change any
    /// reconstruction output, so they are unconditional and stay on in the
    /// knob-off arm.  That is the point -- the knob-off arm is what measures
    /// how often each divergence is reachable at all.
    struct PortAuditCounters {
        // F2 / P9 -- times a point outside every TPC hit each guard.
        std::atomic<uint64_t> oov_isochronous{0};   // modify_segment_isochronous
        std::atomic<uint64_t> oov_dead_scan{0};     // examine_vertices_1p
        std::atomic<uint64_t> oov_unique_scan{0};   // examine_vertices_3
        // P8 -- add_segment endpoint consistency.
        std::atomic<uint64_t> add_segment_calls{0};
        std::atomic<uint64_t> add_segment_reentry{0}; // seg already in the graph => no-op
        std::atomic<uint64_t> endpoint_mismatch{0}; // a connection being MADE is inconsistent
        std::atomic<uint64_t> endpoint_refused{0};  // ... and strict mode dropped it
        // P2 -- local-PCA first-segment endpoint refinement.
        std::atomic<uint64_t> pca_refine_calls{0};
        std::atomic<uint64_t> pca_refine_moved{0};
        std::atomic<uint64_t> pca_move_um_sum{0};   // summed |move| in microns
        std::atomic<uint64_t> pca_move_um_max{0};
        // P4 -- find_other_segments acceptance attribution.
        std::atomic<uint64_t> oseg_accept_proto{0};   // a prototype clause fired
        std::atomic<uint64_t> oseg_accept_relaxed{0}; // ONLY the toolkit clause fired
        std::atomic<uint64_t> oseg_reject{0};
        // doc sbnd_xin/docs/pr/54 -- the isolated-residual discard path had no
        // counter at all (§5a there): a candidate whose endpoints touch no
        // existing segment gets a real do_single_tracking fit and is then
        // silently thrown away.  keep fires only under the default-OFF
        // other_seg_keep_isolated knob; drop is unconditional visibility.
        std::atomic<uint64_t> oseg_isolated_drop{0};
        std::atomic<uint64_t> oseg_isolated_keep{0};
        // doc sbnd_xin/docs/pr/31 §10.10 (F9, was P12) -- the two degenerate
        // topologies the prototype's map_vertex_segments can hold and the PR
        // adjacency_list cannot: a self-loop counts 2 in boost::degree where
        // the prototype's set counts 1, and a second segment on an existing
        // vertex pair silently takes over the first segment's edge
        // (add_segment's !added branch).  Both PERMITTED by add_segment (no
        // guard), reachability unmeasured until now.  If both counters stay 0
        // across the nueCC48 manifest, F9 closes as vacuous.
        std::atomic<uint64_t> selfloop_segment{0};  // add_segment with vtx1 == vtx2
        std::atomic<uint64_t> edge_aliased{0};      // !added and a DIFFERENT live segment owned the edge
    };
    /// The one instance.  Defined in PRGraph.cxx.
    extern PortAuditCounters g_port_audit;

    /// doc sbnd_xin/docs/pr/32 §11 -- stage-4 (neutrino vertex ID) audit
    /// instrumentation, same contract as PortAuditCounters above: process-wide,
    /// one event per wire-cell process, never reset, pure diagnostics.
    struct Pr32AuditCounters {
        // F1 -- how often the fit point and the snapped wcpt actually differ,
        // and by how much (microns, so an integer atomic suffices).
        std::atomic<uint64_t> f1_reads{0};        // vertex-position reads at the 11 sites
        std::atomic<uint64_t> f1_fit_valid{0};    // ... of which had a valid fit
        std::atomic<uint64_t> f1_moved_um_sum{0}; // summed |fit - wcpt| over valid reads
        std::atomic<uint64_t> f1_moved_um_max{0};
        // F2 -- how often the stored kShowerTrajectory flag and a fresh 10 cm
        // recomputation disagree at the two outer gates, and how often the
        // recheck body actually runs.
        std::atomic<uint64_t> f2_gate_calls{0};
        std::atomic<uint64_t> f2_gate_disagree{0};
        std::atomic<uint64_t> f2_body_runs{0};
        std::atomic<uint64_t> f2_flag_cleared{0}; // refresh actually demoted a segment
        // F3 -- candidates the descriptor filter would drop.  If this is 0
        // everywhere, P7 is confirmed dead code and the knob can be retired.
        std::atomic<uint64_t> f3_candidates{0};
        std::atomic<uint64_t> f3_dropped{0};
        // F4 -- vertices flagged as main-vertex candidates.
        std::atomic<uint64_t> f4_flagged{0};
    };
    extern Pr32AuditCounters g_pr32_audit;

    /// doc sbnd_xin/docs/pr/36 §10 -- tagger-stage audit instrumentation,
    /// same contract as PortAuditCounters above: process-wide, one event per
    /// wire-cell process, never reset, pure diagnostics (unconditional, on in
    /// the knob-off arm).
    struct Pr36AuditCounters {
        // F1 (match_isFC): both-ways containment check, counted only when a
        // fiducial is configured (the knob-off arm cannot measure this).
        std::atomic<uint64_t> f1_fc_checks{0};
        std::atomic<uint64_t> f1_fc_disagree{0};  // legacy FiducialUtils vs configured fiducial
        // F2 (PDG-gate population): segments that the prototype's
        // get_particle_type() coercion (ProtoSegment.cxx:10-15) would process
        // as electrons and the toolkit's has_particle_info() gates skip.
        // Sweep = one pass over the graph before the tagger block; per-gate =
        // skips at each gate (0 = NuE:396 bare info gate, 1-9 = the NuE
        // pdg-11 gates in file order, 10 = SinglePhoton's).
        std::atomic<uint64_t> f2_sweep_segments{0};  // segments swept
        std::atomic<uint64_t> f2_sweep_hits{0};      // shower-flagged && !has_particle_info()
        static constexpr int f2_ngates = 11;
        std::atomic<uint64_t> f2_gate_skip[f2_ngates]{};
        // F5 (stem endpoint rule): per-call-site attribution across the 18
        // seg_endpoint_near sites (0-12 = NeutrinoTaggerNuE.cxx in file
        // order, 13-17 = NeutrinoTaggerSinglePhoton.cxx).  disagree = the
        // prototype wcpt rule and the legacy proximity rule pick different
        // fit endpoints; neither = neither wcpt endpoint equals the vertex
        // wcpt EXACTLY (fires => the skeleton-node coincidence premise of
        // doc pr/36 §10.15b is wrong and F5 needs redesign, not tolerance).
        static constexpr int f5_nsites = 18;
        std::atomic<uint64_t> f5_calls[f5_nsites]{};
        std::atomic<uint64_t> f5_disagree[f5_nsites]{};
        std::atomic<uint64_t> f5_neither[f5_nsites]{};
        // F6 (broken_muon cluster count): distinct-cluster-id count differs
        // from distinct-cluster-pointer count (id non-injectivity, doc 53).
        std::atomic<uint64_t> f6_id_vs_ptr_disagree{0};
    };
    extern Pr36AuditCounters g_pr36_audit;

    /// doc sbnd_xin/docs/pr/33 §10 -- EM-shower-clustering audit
    /// instrumentation, same contract as PortAuditCounters above:
    /// process-wide, one event per wire-cell process, never reset, pure
    /// diagnostics (unconditional, on in the knob-off arm).
    struct Pr33AuditCounters {
        // F1 (daughter-count callee): both callees computed at both sites.
        // differ = the consumed value (.first at the main-vertex site,
        // .second at the examine_showers site) differs between
        // calculate_num_daughter_showers (current) and _tracks (prototype).
        // gate_flip = the :303 proton-skip verdict itself would change.
        std::atomic<uint64_t> f1_mv_calls{0};
        std::atomic<uint64_t> f1_mv_differ{0};
        std::atomic<uint64_t> f1_mv_gate_flip{0};
        std::atomic<uint64_t> f1_ex_calls{0};
        std::atomic<uint64_t> f1_ex_differ{0};
        // F2 (whose PDG): per-site skip-verdict disagreement between the
        // legacy reading and the full prototype reading.  Site ids:
        // 0=:170 (in_main_cluster), 1=:525 (from_main_cluster, inverted),
        // 2=:1247 (examine_merge_showers), 3=:2193 (examine_showers),
        // 4=:2911 5=:2927 (id_pi0_without_vertex).
        static constexpr int f2_nsites = 6;
        std::atomic<uint64_t> f2_calls[f2_nsites]{};
        std::atomic<uint64_t> f2_disagree[f2_nsites]{};
        // F3 (pi0 id allocator): pi0s minted per finder.  Both nonzero in
        // one event = the collision the by-value seed makes possible.
        std::atomic<uint64_t> f3_pi0_with_vertex{0};
        std::atomic<uint64_t> f3_pi0_without_vertex{0};
        // F4 (get_flag_shower's abs(pdg)==11 term): segments at the :797
        // site where appending the term would flip is_shower.
        std::atomic<uint64_t> f4_flip{0};
        // F5 (shower_less tie-break): entries into the same-index
        // pointer-address fallback at :2848.  0 = unreachable on this
        // manifest (pr/32 P7 precedent: 0/2219).
        std::atomic<uint64_t> f5_fallback_hits{0};
    };
    extern Pr33AuditCounters g_pr33_audit;

    /// P8 knob transport.  add_segment() is a free function reached from ~30
    /// call sites and has no access to component configuration, so the two
    /// values it needs are written once per event by TaggerCheckNeutrino
    /// before any graph is built.  Read-mostly; never written concurrently
    /// with graph construction.
    struct GraphEndpointPolicy {
        bool   strict{false};              // refuse an inconsistent connection
        double tol{0.3 * units::cm};       // positional stand-in for wcpt-index equality
    };
    extern GraphEndpointPolicy g_graph_endpoint_policy;


    /// Make a vertex and add it to the graph.
    ///
    /// This method will properly adhere to the indexing policy and is
    /// recommended instead of creating a PR::Vertex by hand.
    ///
    /// This is templated to allow for non-default construction of the
    /// PR::Vertex.  Normal usage:
    ///
    /// @code{.cpp}
    /// auto my_vtx = make_vertex(my_graph);
    /// @endcode
    template <typename... Args>
    VertexPtr make_vertex(Graph& g, Args&&... args) {
        auto ptr = std::make_shared<Vertex>(std::forward<Args>(args)...);
        bool ok = add_vertex(g, ptr);
        if (!ok) {
            raise<RuntimeError>("make_vertex failed to add vertex which should never happen!");
        }
        return ptr;
    }

    /// Add an existing vertex to the graph.
    ///
    /// Return true if graph modified.
    ///
    /// If the PR::Vertex already has a valid descriptor, false is returned.
    ///
    /// Otherwise, a new graph node is added with its bundle holding this vertex
    /// and a unique index.  The resulting descriptor is set on the vertex and
    /// true is returned.
    ///
    /// User must take responsibility to provide a vertex that is consistent
    /// with the indexing policy.  For a safer function, use `make_vertex()`.
    bool add_vertex(Graph& g, VertexPtr vtx);


    /// Remove a vertex from the graph.
    ///
    /// Return true if graph modified.
    ///
    /// The descriptor held by the object is left invalid.
    bool remove_vertex(Graph& graph, VertexPtr vtx);


    /// Make a PR::Segment
    ///
    /// This makes an "orphaned" segment not associated to any graph.  Call
    /// `add_segment()` to add it to a graph, or better, use `make_segment()`
    /// to do both at once.
    template <typename... Args>
    SegmentPtr make_segment(Args&&... args) {
        return std::make_shared<Segment>(std::forward<Args>(args)...);
    }

    /// Add a PR::Segment to the graph as an edge between the nodes of two
    /// PR::Vertex instances.
    ///
    /// Return true if the graph was modified.  Modification means that a new
    /// vertex or edge was added.  In the case that vtx1 and vtx2 already have
    /// an edge, it is updated with the segment and that does not influence
    /// true/false return value.
    ///
    /// If the PR::Segment already has a valid descriptor then it is not added.
    ///
    /// If either PR::Vertex do not have a valid descriptor then they are added
    /// to the graph and will cause a true return value irrespective of the
    /// status of the segment.  But, users are urged to only use `make_index()`
    /// to construct a PR::Vertex which will assure graph addition.
    bool add_segment(Graph& g, SegmentPtr seg, VertexPtr vtx1, VertexPtr vtx2);


    /// Make a PR::Segment and add it to the graph as an edge between two vertices.
    ///
    template <typename... Args>
    SegmentPtr make_segment(Graph& g, VertexPtr vtx1, VertexPtr vtx2, Args&&... args) {
        SegmentPtr seg = std::make_shared<Segment>(std::forward<Args>(args)...);
        add_segment(g, seg, vtx1, vtx2);
        return seg;
    }


    /// Remove a segment from the graph.
    ///
    /// Return true if graph modified.
    ///
    /// The descriptor held by the object is left invalid.
    bool remove_segment(Graph& graph, SegmentPtr seg);

    /// doc pr/67 round 2, P6 -- enable the log-only removal sentinel inside
    /// remove_segment().  remove_segment is a free function with ~55 call sites
    /// across six files, so it cannot read PatternAlgorithms::m_traj_cover_probe;
    /// TaggerCheckNeutrino mirrors the knob here at configure time instead of
    /// tagging every call site.  Default false => byte-identical when off.
    void set_traj_cover_probe(bool enable);

    
    /// Return the end-point vertices of a segment.
    ///
    /// The pair will be nullptr if segment edge not in graph.
    ///
    /// The two Vertex objects are those associated with the source/target nodes
    /// of the segment's edge.  The pair is ordered.  The first Vertex is the
    /// one with a "wcpoint" closest to the segment's initial "wcpoint".
    std::pair<VertexPtr, VertexPtr> find_vertices(Graph& graph, SegmentPtr seg);

    /// Return the other vertex connected to a segment, given one known vertex.
    ///
    /// Returns nullptr if:
    /// - The segment edge is not in the graph
    /// - The given vertex is not connected to the segment
    ///
    /// This is useful when you have one vertex of a segment and need to find
    /// the vertex at the other end.
    VertexPtr find_other_vertex(Graph& graph, SegmentPtr seg, VertexPtr vertex);

    /// Find the segment connecting two vertices.
    ///
    /// Returns nullptr if:
    /// - Either vertex is not in the graph
    /// - No edge exists between the two vertices
    ///
    /// This function searches for an edge between the two given vertices and
    /// returns the associated segment if found.
    SegmentPtr find_segment(Graph& graph, VertexPtr vtx1, VertexPtr vtx2);

    /** Comparators for VertexPtr and SegmentPtr that order by the stable graph
     *  index stored directly in the object rather than by raw pointer address.
     *  No graph reference is required.
     *
     *  Use these (via IndexedVertexSet / IndexedSegmentSet) wherever a
     *  std::set<VertexPtr> or std::set<SegmentPtr> must produce a deterministic
     *  iteration order across runs.
     */
    /// The value Vertex::m_graph_index / Segment::m_graph_index hold before the
    /// object has ever been handed to PR::add_vertex / PR::add_segment.
    ///
    /// This is the index-keyed containers' one real hazard, and it is much
    /// broader than the rare "edge already existed, inherit its index" path in
    /// add_segment(): EVERY not-yet-added object carries this same value, so
    /// under VertexIndexCmp/SegmentIndexCmp they all compare EQUAL to one
    /// another.  An IndexedSegmentSet handed two such segments silently keeps
    /// one and drops the other, and find() on either returns whichever is held.
    /// Probing six SBND events found zero live occurrences, so this is latent
    /// rather than firing -- the warning below is what keeps it that way.
    constexpr size_t kUnindexed = std::numeric_limits<size_t>::max();

    /// Emit a one-time warning that an index-ordered container was asked to
    /// compare an object with no graph index.  noinline+cold on purpose: these
    /// comparators run O(n log n) inside clustering, so the hot path must stay
    /// at the two integer loads and compares below.  Defined here rather than in
    /// PRGraph.cxx because wcdoctest-clus compiles its own subset of clus/src
    /// and does not link libWireCellClus.
    [[gnu::noinline, gnu::cold]]
    inline void warn_unindexed(const char* what)
    {
        auto log = Log::logger("clus.PRGraph");
        log->warn("{} with no graph index reached an index-ordered container.  "
                  "All such objects compare EQUAL, so a set/map may silently drop "
                  "one and find() may return the wrong one.  Add the object to a "
                  "graph first, or key the container on the pointer.  Reported "
                  "once per process (doc pr/28 sec 14.7).",
                  what);
    }

    struct VertexIndexCmp {
        bool operator()(const VertexPtr& a, const VertexPtr& b) const {
            const size_t ia = a->get_graph_index();
            const size_t ib = b->get_graph_index();
            if (ia == kUnindexed || ib == kUnindexed) [[unlikely]] {
                // Function-local so the atomic is shared across every
                // instantiation and the check stays thread-safe; wire-cell runs
                // these comparators from several worker threads.
                static std::atomic<bool> warned{false};
                if (!warned.exchange(true)) warn_unindexed("PR::Vertex");
            }
            return ia < ib;
        }
    };
    struct SegmentIndexCmp {
        bool operator()(const SegmentPtr& a, const SegmentPtr& b) const {
            const size_t ia = a->get_graph_index();
            const size_t ib = b->get_graph_index();
            if (ia == kUnindexed || ib == kUnindexed) [[unlikely]] {
                static std::atomic<bool> warned{false};
                if (!warned.exchange(true)) warn_unindexed("PR::Segment");
            }
            return ia < ib;
        }
    };
    using IndexedVertexSet  = std::set<VertexPtr,  VertexIndexCmp>;
    using IndexedSegmentSet = std::set<SegmentPtr, SegmentIndexCmp>;

    /// Return every segment whose edge is NOT reachable from `root` by a
    /// breadth-first walk of the graph.
    ///
    /// Pure topology: no cluster or flag filtering (callers filter).  If
    /// `root` is null or not in the graph, every segment in the graph is
    /// returned.  The walk keys on membership only, so iteration order cannot
    /// leak into the result; the returned set is ordered by stable graph
    /// index.
    ///
    /// doc sbnd_xin/docs/pr/65 round 3: the shower absorbers' cluster()==
    /// main_cluster guards encode "already claimed by the main_vertex graph
    /// walk", which other_seg_keep_isolated (doc pr/54) broke by adding
    /// disconnected components to the main cluster's graph.  This helper is
    /// what lets those guards test reachability instead of cluster identity.
    IndexedSegmentSet unreachable_segments(const Graph& graph, VertexPtr root);

};
#endif
