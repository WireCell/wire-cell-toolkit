#ifndef WIRECELL_CLUS_MULTIALGBLOBCLUSTERING
#define WIRECELL_CLUS_MULTIALGBLOBCLUSTERING

#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/IClusGeomHelper.h"
#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/IBeeSink.h"
#include "WireCellClus/Facade.h"

#include "WireCellAux/Logger.h"

#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/ITerminal.h"

#include "WireCellUtil/Bee.h"

#include <vector>
#include <set>
#include <map>
#include <string>

namespace WireCell::Clus {

    class MultiAlgBlobClustering
        : public Aux::Logger, public ITensorSetFilter, public IConfigurable, public ITerminal
    {
       public:
        MultiAlgBlobClustering();
        virtual ~MultiAlgBlobClustering() = default;

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_pointer& in, output_pointer& out);

        virtual void finalize();

       private:

        /** Config: bee_zip
         *
         * The path to a zip file to fill with Bee JSON
         */
        Bee::Sink m_sink;
        int m_last_ident{-1};
        int m_initial_index{0};  // Default to 0 for backward compatibility

        // Optional shared Bee sink (config "bee_sink").  When set, all Bee
        // writes go to this shared single-zip sink (at an explicit per-event
        // index) instead of the per-node m_sink above, so a multi-node chain
        // produces one .zip.  Default: unset -> own m_sink (back-compat).
        IBeeSink::pointer m_shared_sink{nullptr};
        bool m_use_shared_sink{false};
        size_t m_bee_event_index{0};
        
        // Drift-side / APA grouping shared by the clustering bee output and the
        // dead-area output.  When configured, each cluster (or dead blob) is
        // routed by its APA into the matching group and dumped as a single bee
        // instance named "<algorithm>-<group name>".  Unset -> behavior unchanged.
        struct ApaGroup {
            std::string name;       // bee instance suffix, e.g. "group02"
            std::set<int> apas;     // APA idents that belong to this group
        };

        // Replace the existing bee points structures with a more flexible approach
        struct BeePointsConfig {
            // special name "img" dumps "live" before clustering
            // prepipeline=true dumps this set at that SAME pre-clustering point
            // (the raw per-anode matched clusters) without needing name=="img",
            // so a second raw set (e.g. per-drift-side img-side-* split by
            // apa_groups) can coexist with img-global.  Default false =>
            // set follows the legacy name=="img" pre/post routing (byte-identical).
            bool prepipeline{false};
            std::string name;   // bee type name
            std::string detector;  // bee geom name
            std::string algorithm; // bee alg name, defaults to type
            std::string pcname;    // PC to take
            std::string grouping;  // grouping to take (default to "live")
            std::string visitor;   // if given, dump just after this visitor runs and any cluster ID enumeration
            std::vector<std::string> coords;
            bool individual;
            int filter{1};// 1 for on, 0 for off, -1 for inverse filter
            double dQdx_scale{1.0};
            double dQdx_offset{0.0};
            bool use_associate_points{false};  // use dpcloud("associate_points") + shower-based charge
            bool use_graph_vertices{false};    // dump graph vertices; charge=15000 for main (kNeutrinoVertex), 0 otherwise
            // This set is PR output only: when the bound visitor produced no PR
            // graph (e.g. TaggerCheckNeutrino took its "no main cluster
            // selected" exit) leave the set EMPTY instead of falling back to the
            // generic cluster dump.  Default false = legacy fallback, which is
            // what the uBooNE "regular"/"steiner" sets bound to CreateSteinerGraph
            // rely on.  See sbnd_xin/docs/pr/3 sec. 9.
            bool require_pr_graph{false};
            // Prototype-parity options (default false => legacy output, byte-identical):
            // particle_ids: with use_associate_points, real_cluster_id follows the
            // prototype per-particle convention (NeutrinoID::fill_point_info):
            // every non-shower segment gets cluster_id*1000 + seg id instead of
            // the plain cluster_id collapse.
            bool particle_ids{false};
            // pseudo_shower_track_paint: with use_associate_points, a segment
            // belonging to a shower whose cached particle type is +-13 (a
            // long-muon pseudo-shower, seeded from segments_in_long_muon) is
            // painted as TRACK (q=0) instead of shower (q=15000).  The PF tree
            // (make_shower_leaf) displays the same shower as "mu-" from the
            // same cached field, so the legacy membership-first paint is
            // provably inconsistent for this class (doc sbnd_xin/docs/pr/45,
            // SBND 18255-56463: 411 cm cathode-crossing muon painted red).
            // Default false => legacy paint, byte-identical.
            bool pseudo_shower_track_paint{false};
            // include_vertex_points: PR-graph edge dumps (track_fit) also append
            // each vertex fit point with real_cluster_id=-1 (prototype
            // fill_skeleton_info_magnify vertex rows).
            bool include_vertex_points{false};

            // require_flag: restrict this set to clusters carrying the named
            // tagger flag, e.g. "STM" (Facade flags are flag_<NAME> scalars,
            // Facade_Mixins.h set_flag/get_flag; names in ClusteringFuncs.h's
            // Flags namespace).  Empty => every cluster, i.e. the legacy dump,
            // so configs without the key are byte-identical.  Used by PDVD's
            // stm / steiner_graph / steiner_terminals layers (doc pdvd/39).
            std::string require_flag;
            // steiner_terminals_only: for a pcname=="steiner_pc" set, dump only
            // the points flagged flag_steiner_terminal.  Default false => the
            // whole Steiner node cloud, byte-identical.
            bool steiner_terminals_only{false};

            // Optional per-set drift-side / APA grouping (see ApaGroup above).
            // Non-empty -> route this set's clusters into group buckets.
            std::vector<ApaGroup> apa_groups;
        };

        // Vector to store configurations for multiple bee points sets
        std::vector<BeePointsConfig> m_bee_points_configs;
        
         // Nested structure to store bee points objects for each configuration, by APA and face
        // First key: bee points set name, second key: "anode_id-face_id" string
        struct ApaBeePoints {
             // Default constructor (add this)
            ApaBeePoints() = default;
            
            // Global points (used when individual == false)
            Bee::Points global;
            
            // Individual points (used when individual == true)
            // Key is "anode_id-face_id" string
            std::map<int, std::map<int , Bee::Points> > by_apa_face; // apa, face

            // Grouped points (used when apa_groups is non-empty).
            // Key is the group name.
            std::map<std::string, Bee::Points> by_group;
        };
    
        Facade::Grouping& load_grouping(
            Facade::Ensemble& ensemble,
            const std::string& name,
            const std::string& path,
            const ITensorSet::pointer ints);

        std::map<std::string, ApaBeePoints> m_bee_points;

        // New helper function to fill bee points
        // doc 53: re-stamp "real_cluster_id" into ONE globally unique ident
        // epoch.  Called once per event after the clustering pipeline and
        // BEFORE the Bee fills, so the Bee per-blob labels and the saved pctree
        // carry the same ids.
        void restamp_real_cluster_id(Facade::Grouping& grouping) const;
        void fill_bee_points(const std::string& name, const Facade::Grouping& grouping);
        void fill_bee_points_from_cluster(
            Bee::Points& bpts, const Facade::Cluster& cluster,
            const std::string& pcname, const std::vector<std::string>& coords,
            int filter, double dQdx_scale = 1.0, double dQdx_offset = 0.0,
            bool steiner_terminals_only = false);
        // doc pr/94 Phase 4b: `tf_in` selects WHICH per-bundle TrackFitting to
        // render (null = the unnamed slot, i.e. the legacy single-candidate
        // behavior).  `do_reset` must be false on every call after the first of
        // a multi-bundle sequence: both functions reset the Bee::Points object
        // at entry, so resetting per bundle would leave only the LAST bundle's
        // points -- a bug that hides itself, since the symptom being fixed
        // ("the second candidate has no points") would still look cured.
        void fill_bee_points_from_pr_graph(const std::string& name, const Facade::Grouping& grouping,
                                           std::shared_ptr<WireCell::Clus::TrackFitting> tf_in = nullptr,
                                           bool do_reset = true);
        void fill_bee_vertices_from_pr_graph(const std::string& name, const Facade::Grouping& grouping,
                                             std::shared_ptr<WireCell::Clus::TrackFitting> tf_in = nullptr,
                                             bool do_reset = true);

        void fill_bee_patches_from_grouping(const Facade::Grouping& grouping);
        void fill_bee_patches_from_cluster(const Facade::Cluster& cluster);

        // ---- Particle-flow Bee output ----
        // Triggered after TaggerCheckNeutrino (or any configured visitor) runs.
        // Produces one file per event named "mc" (bare JSON array), matching the
        // prototype "mc" format read by the Bee viewer.
        // (public so doctest_clus_knob_defaults can pin the in-class defaults)
       public:
        struct BeePFConfig {
            std::string name{"mc"};          // Bee file name (default "mc")
            std::string visitor;             // dump after this visitor runs
            std::string grouping{"live"};    // grouping to read PR graph from
            // Prototype-parity options (defaults => legacy output, byte-identical):
            // prototype_names: node text uses the prototype's TDatabasePDG-style
            // names ("proton" not "p") and integer MeV ("proton  76 MeV").
            bool prototype_names{false};
            // KeepMC-style kinetic-energy floors (system units; 0 = keep all):
            // em_ke_min applies to gamma/e+-, np_ke_min to n/p/nuclei
            // (prototype WCReader::KeepMC: 5 MeV / 10 MeV).
            double em_ke_min{0.0};
            double np_ke_min{0.0};
            // ---- doc sbnd_xin/docs/pr/34 §10 port-fidelity knobs ----
            // pf_track_main_cluster_only: the track BFS skips segments whose
            // cluster ident differs from the main vertex's, restoring the
            // prototype's dropped guard (NeutrinoID.cxx:1488).  Ident, not
            // pointer: split products inherit the parent's ident.
            bool pf_track_main_cluster_only{false};
            // doc sbnd_xin/docs/pr/40 round 9 B2: let the track BFS traverse
            // clusters recorded by TrackFitting::get_bridged_cluster_ids()
            // (connected to the main cluster by an nv_bridge_track bridge
            // segment) despite pf_track_main_cluster_only.  Off, or an
            // empty bridged set => byte-identical to legacy.
            bool pf_track_bridged_clusters{false};
            // pf_shower_vertex_barrier: pre-seed visited_vtxs from every
            // shower's vertex set so the BFS does not expand THROUGH a shower
            // (prototype used_vertices seed, NeutrinoID.cxx:1597-1602).
            // Corrected semantics (doc pr/38): the barrier set excludes each
            // shower's START vertex -- the prototype's map_vtx_segs never
            // holds it (WCShower.cxx:547, :708-716, :733-745), so main-track
            // attachment junctions stay traversable and only shower-INTERIOR
            // vertices block.  The same knob also enables the orphan safety
            // net: BFS-unreached non-shower main-cluster segments are emitted
            // as root-level leaves (prototype flat-loop mc_mother=0 parity,
            // NeutrinoID.cxx:1485-1489) instead of silently vanishing.
            bool pf_shower_vertex_barrier{false};
            // pf_shower_parent_precedence: record the parent shower for
            // track-attached shower vertices and consult it before the
            // incoming track segment when a shower picks its parent
            // (prototype map_vertex_in_shower-first order, :1655/:1680/:1720).
            bool pf_shower_parent_precedence{false};
            // pf_pi0_node_per_id: one pi0 node per pi0 id instead of one per
            // parent (prototype map_pio_id_saved_pair, :1326/:1361).  The
            // merged node's home is the HIGHEST-ENERGY daughter's parent --
            // owner decision 2026-08-04, deliberately NOT the prototype's
            // first-writer-wins (a jsTree node has exactly one parent).
            bool pf_pi0_node_per_id{false};
            // pf_pdg_name_prototype_fallback: pi0/nuclei entries in the PDG
            // name table plus the prototype's numeric fallback
            // (WCReader.cc:529-547) instead of "particle".
            bool pf_pdg_name_prototype_fallback{false};
            // ---- doc sbnd_xin/docs/pr/38 Round 4 ----
            // pf_orphan_track_parentage: upgrade the pf_shower_vertex_barrier
            // orphan safety net from flat root-leaves to graph-faithful
            // parentage.  An orphan whose endpoint vertex carries a claimed
            // incoming track segment attaches as that segment's child; else an
            // endpoint inside a shower's view attaches as a child of that
            // shower's displayed leaf; orphan-of-orphan chains attach to each
            // other (pi+ -> proton).  Only orphans with no anchor at all fall
            // back to the legacy flat root emission.  This state is UNREACHABLE
            // in the prototype (it has no shower_absorb_track_guard, so such
            // tracks are absorbed into the shower) -- designed divergence, see
            // porting_dictionary.  Inert unless pf_shower_vertex_barrier is
            // also on.  C++ default false => byte-identical legacy output.
            bool pf_orphan_track_parentage{false};
            // ---- doc sbnd_xin/docs/pr/65 round 3 ----
            // pf_orphan_audit_only: rung 3 of the pr/65 ladder.  The flat
            // orphan safety net fabricates a root-level PF particle for every
            // BFS-unreached main-cluster segment -- a faithful port of the
            // prototype's mc_mother=0 default (NeutrinoID.cxx:1485-1489), but
            // a state the prototype itself can never reach (its residual
            // segments are attach-or-discard).  When on, keep the net's
            // VISIBILITY but drop the FABRICATION: emit one log line per
            // still-unclaimed segment (no display filters -- dirsign==0 and
            // empty-fit segments are named too, today they vanish silently)
            // and append no node.  Load-bearing only together with rung 1
            // (shower_absorb_unreachable_main): alone it would just hide the
            // charge.  Inert unless pf_shower_vertex_barrier is also on.
            // C++ default false => byte-identical legacy output.
            bool pf_orphan_audit_only{false};
            // ---- doc sbnd_xin/docs/pr/84 round 2 (display-only family) ----
            // pf_direct_when_touching: a conn-2/3 shower whose fitted charge
            // comes within pf_touch_max of the main vertex is rendered as a
            // direct leaf instead of under a synthetic gamma/neutron carrier.
            // The carrier means "the PR graph could not walk there"; when the
            // charge demonstrably touches the vertex the neutral parent is a
            // graph artifact, not physics (evts 283713/316025/407280).
            // start_connection_type is NOT modified -- kinematics,
            // mc_included and every tagger are untouched; only the rendered
            // tree (mabc zip 0-mc.json) can change.  pi0 daughter gammas are
            // exempt (their carrier is correct by construction).
            // C++ default false => byte-identical legacy output.
            bool pf_direct_when_touching{false};
            double pf_touch_max{3.0 * units::cm};   // read only when F1 on
            // doc 77 round 1 (2026-08-24): pf_touch_cross_main/_max removed --
            // zero movers on all 7 census candidates, F1.0 probe failure
            // (pr/84 sec 607/622).  See sbnd_xin/docs/77_knob-ledger.tsv.
            // pf_pseudo_gap_from_main (= pr/84 P3): in the "start_vtx not in
            // BFS tree" fallback, anchor the synthetic carrier at the MAIN
            // vertex instead of the shower's own start vertex, so a genuinely
            // remote association draws its real gap instead of collapsing to
            // a zero-length node (pr/84 sec 4: 15 such, median 38.8 cm).
            // C++ default false => byte-identical legacy output.
            bool pf_pseudo_gap_from_main{false};
            // pf_unique_node_ids (doc pr/84 round 3, G1): mc.json is jsTree
            // input and jsTree keys its model by node id, so a repeated id is
            // invalid output -- observed on SBND 394532, where a node and its
            // own descendant carried id 8033 and selecting it blanked the
            // whole PF panel.  Natural ids follow the prototype convention
            // (`cluster_id*1000 + seg_id`, NeutrinoID.cxx:1268) and are NOT
            // unique: two showers can share a start segment (see
            // m_shower_dedup_start_seg, which fixes that at the source), a
            // shower leaf can collide with the track node of the same
            // segment, and the pseudo/pi0 counter starts at 1 inside the same
            // number space.  When on, a colliding node is re-issued a fresh
            // unused id and the collision is logged -- the id is used for
            // nothing but jsTree identity (bee3 events/static/js/bee/physics/
            // mc.js draws from data.start/data.end only).  C++ default false
            // => byte-identical legacy output.
            bool pf_unique_node_ids{false};
            // pf_drop_stray_satellites (doc pr/92): skip showers whose id is
            // in TrackFitting::get_dropped_satellite_shower_ids() -- the
            // conn-2/3 satellites fill_kine_tree's kine_drop_stray_satellites
            // gate removed from kine_reco_Enu (overclustered cosmics, second
            // neutrinos).  Effective only when that kine knob is also on;
            // otherwise the id set is empty and this knob is inert.  C++
            // default false => byte-identical legacy output.
            bool pf_drop_stray_satellites{false};
            // pf_orphan_confident_track (doc pr/93 round 4): in the
            // pf_orphan_audit_only branch, EMIT a root PF node for an
            // unclaimed main-cluster segment that passes
            // segment_orphan_confident_track (confident non-electron
            // template PID + length > pf_orphan_track_min + straight-long).
            // pr/65 rung 4: the audit chose visibility-without-fabrication
            // for the general population; this narrow class (e.g. SBND
            // 18255-315167's 150.7cm score-0.101 proton, freed from shower
            // membership by shower_cone_absorb_guard but graph-disconnected
            // from the main vertex) is a real particle the owner wants in
            // the PF.  Audit log lines are unchanged for every segment.
            // C++ default false => byte-identical legacy output.
            bool pf_orphan_confident_track{false};
            double pf_orphan_track_min{50.0 * units::cm};  // read only when the bool is on
            // pf_orphan_guard_freed (doc pr/123 round 2): emit a root PF
            // node for an unclaimed segment carrying SegmentFlags::
            // kPass4GuardFreed -- a track the pass4 long-track guard
            // declined (SBND 18255-171572's 125cm muon), which is outside
            // the pf_orphan_confident_track scope (cross-cluster,
            // score-100 sentinel PID).  The flag IS the predicate: the
            // guard's own decline set, nothing wider (the 120-segment
            // any-cluster unclaimed population -- largely cosmics -- stays
            // invisible).  false = legacy = byte-identical.
            bool pf_orphan_guard_freed{false};
            // ---- doc sbnd_xin/docs/pr/128 (PF/kine completeness) ----
            // pf_orphan_near_cross_cluster: emit a PF node for an unclaimed
            // CROSS-CLUSTER track segment whose fit points come within
            // pf_orphan_near_gap of the emitted candidate.  Every PF orphan
            // pool -- and the pr/65 audit line itself -- is same_cluster
            // gated (prototype NeutrinoID.cxx:1488), so this class reaches no
            // output and is not even counted: 29 objects in 18 of 239 SBND
            // events, 21 within 10 cm of displayed content, many at gap 0.00,
            // against a main-cluster control of 1 (doc pr/128 §1).  The
            // 137238 class of doc pr/127, generalised.  Admission is
            // segment_near_candidate_track + the proximity test; rendering
            // follows the pr/123 displaced-object convention (pseudo-neutron
            // carrier, owner correction 2026-08-28).  PF ONLY -- the energy
            // moves with the kine twin kine_count_near_cross_cluster.
            // false = legacy = byte-identical.
            bool pf_orphan_near_cross_cluster{false};
            double pf_orphan_near_gap{5.0 * units::cm};      // read only when the bool is on
            double pf_orphan_near_min_len{30.0 * units::cm}; // read only when the bool is on
            // Continuation terms (doc pr/128 §3.1).  Proximity alone admits
            // cosmics that brush the far end of a displayed track: SBND
            // 18255-72786 gained +1151 MeV on a 701 MeV candidate that way.
            // A continuation joins END to END and runs STRAIGHT ON, the same
            // discriminator doc pr/127's sccc fix uses.
            double pf_orphan_near_end_tol{10.0 * units::cm};
            double pf_orphan_near_kink_deg{30.0};
            // pf_conn4_near_candidate: stop skipping a conn-4 shower whose
            // material is the candidate's own -- i.e. whose closest approach
            // to the main cluster is within pf_conn4_near_gap.  conn-4 means
            // "cluster >80 cm from the candidate" (NeutrinoShowerClustering.
            // cxx:3733) and its skip is correct for 490 of the 514 conn-4
            // showers measured (2815 of 3514 MeV, all >=50 cm away -- the
            // far-away over-clustered activity the owner ruled must NOT be
            // counted, 2026-08-29).  But three of our OWN passes stamp conn-4
            // on material they just rescued or shed -- pr/74
            // conn3_unreachable (:3858), pr/123 pass4_prune (:6435), pr/124
            // pass4_prune2 (:6645) -- and conn-4 is the one label that means
            // show nowhere and count nowhere.  That is 481.5 MeV in 4 events
            // (doc pr/128 §2), e.g. SBND 18255-105074's two pdg-13 showers,
            // 215.1 + 162.0 MeV, 82.9 and 58.2 cm, in the main cluster.
            // false = legacy = byte-identical.
            bool pf_conn4_near_candidate{false};
            double pf_conn4_near_gap{20.0 * units::cm};      // read only when the bool is on
            // pf_track_owns_loose_vertex (doc pr/93 round 4): in the F3a
            // root branch, a root-anchored shower's fill_sets() vertex VIEW
            // overrides the track BFS wherever the two disagree, and
            // pf_shower_parent_precedence then hangs everything anchored
            // there under the shower.  The view can hold a vertex none of
            // whose incident segments the shower owns -- the F12 absorb
            // guard add_vertex()es the frontier BEFORE refusing the segment
            // beyond it (PRShower.cxx guard_excludes walk termination), and
            // add_shower()/add_segment(seg,true) can do the same.  When on,
            // skip the claim when BOTH (a) the REAL track BFS walked a
            // segment to the vertex, and (b) the vertex is not an endpoint
            // of any MEMBER segment of the claiming shower (pure loose
            // association).  This is the general "track BFS beats shower
            // set" fix the pr/74-r4 comment deferred; superset of the
            // kMuonStemGuard protection for the loose-association case only
            // (both are kept).  SBND 18264-69314: the 151.9cm muon's far
            // endpoint (deg 2) is claimed by the 595 MeV root shower whose
            // nearest member is 35cm away, stealing the muon's own 67 MeV
            // conn-1 daughter shower.  Render-only (mc.json): parentage
            // changes only, never membership; kine reads none of these
            // maps.  C++ default false => byte-identical legacy output.
            bool pf_track_owns_loose_vertex{false};
        };
       private:
        std::vector<BeePFConfig> m_bee_pf_configs;

        // Storage: flushed at end of each event (same lifecycle as m_bee_points)
        std::map<std::string, WireCell::Bee::ParticleTree> m_bee_pf_trees;

        /// Render one neutrino candidate's particle flow into the named Bee
        /// tree.  The three trailing arguments are doc pr/94 Phase 4 and all
        /// default to the pre-pr/94 single-candidate behaviour:
        ///   tf_in            - render THIS TrackFitting instead of resolving
        ///                      the grouping's single unnamed slot implicitly
        ///                      (with N per-bundle fitters the unnamed slot is
        ///                      only ever bundle 0).
        ///   shared_used_ids  - the pf_unique_node_ids reissue set, hoisted to
        ///                      the caller's scope so reissued ids cannot
        ///                      collide BETWEEN bundles (each call otherwise
        ///                      restarts its own set at 1000000).
        ///   out_particles    - when non-null, this bundle's roots are wrapped
        ///                      in one synthetic node and APPENDED here instead
        ///                      of being handed to set_particles(), which is a
        ///                      plain overwrite (Bee.cxx:549-551) and would
        ///                      otherwise make bundle i erase bundle i-1.
        void fill_bee_pf_tree(const BeePFConfig& cfg, const Facade::Grouping& grouping, bool flag_print = false,
                              std::shared_ptr<WireCell::Clus::TrackFitting> tf_in = nullptr,
                              std::set<int>* shared_used_ids = nullptr,
                              Configuration* out_particles = nullptr);

        std::map<int, std::map<int, Bee::Patches>> m_bee_dead_patches;
        // Bee::Patches m_bee_dead; // dead region ...

        // Dead-area drift-side grouping (mirrors the clustering grouping).  When
        // m_dead_apa_groups is non-empty, dead blobs are routed by APA into one
        // Bee::Patches per group ("channel-deadarea-<group name>") instead of one
        // per (apa,face).  Empty -> per-(apa,face) output (unchanged).
        std::vector<ApaGroup> m_dead_apa_groups;
        std::map<std::string, Bee::Patches> m_bee_dead_groups;

        // ---- Optical flash / charge-light "op" Bee dump ----
        // When m_save_opflash is set, the merged-grouping root carries a
        // self-contained per-flash "opflash" point cloud (gid, time, ch, pe)
        // written by the upstream QLMatching, and matched clusters carry a
        // "matched_flash_gid" scalar + a "flashpred" PC (predicted per-channel
        // PE).  At the pre-pipeline "img" dump point we read these and emit the
        // "op" display so the flash/QL-matching result lands in the same zip
        // as the charge clusters.  Default OFF (no flash dump for detectors
        // that don't attach an "opflash" PC, e.g. uboone).
        bool m_save_opflash{false};
        // When set, fill_bee_flashes emits one "op" row per flash carrying ALL
        // matched cluster ids (op_cluster_ids array) with element-wise summed
        // predicted PE, instead of one row per (flash, cluster).  Default OFF so
        // existing output is bit-identical; enabled for the SBND all-APA match.
        bool m_bee_flash_per_flash{false};
        // doc pr/94 round 3.  Minimum total predicted light (PE) a matched
        // cluster must carry to appear in the "op" display's op_cluster_ids.
        // The legacy dump_light used a hard-coded 100 PE, which is why SBND
        // 18255/73038's 26.5 cm cathode activity -- genuinely matched to the
        // beam flash gid 14, but predicting only 3.6 PE of that flash's 602.6
        // PE -- was drawn as matched to NO flash while the PR chain happily
        // reconstructed it (owner Bee scan, doc pr/94 sec 9.9).  Default 100 =
        // the legacy filter, byte-identical output; 0 shows every genuine
        // match.
        double m_bee_flash_pred_min{100.0};
        // When > 0, group the root opflash flashes across both TPC sides by this
        // ±time window (stored as a per-flash "group" array on the root opflash
        // PC, pre-pipeline) so the Bee viewer can show a TPC0/TPC1 coincidence
        // together.  Default 0 = off, output bit-identical; set for SBND all-APA.
        double m_flash_group_window{0.0};
        // Flash-pair grouping policy when several matched flashes coincide within
        // the window.  false (default) = ONE closest cross-side pair per
        // coincidence neighborhood, the rest left solo.  true = GREEDY disjoint
        // pairs: repeatedly take the next-closest available cross-side pair, so
        // two distinct close crossers in one neighborhood both group.  Default
        // false keeps the grouped output bit-identical to the one-pair policy.
        bool m_flash_group_greedy{false};
        Bee::Flashes m_bee_flash;
        void fill_bee_flashes(const Facade::Grouping& grouping);

        // Add new member variables for run/subrun/event
        int m_runNo{0};
        int m_subRunNo{0};
        int m_eventNo{0};
        bool m_use_config_rse{false};  // Flag to determine if we use configured RSE
        // When set, take the event number from the per-event tensor-set ident
        // (m_eventNo = ident, run/subrun = 0).  Used by the bundled standalone
        // chain whose ident already carries the real event id.  Default off keeps
        // the existing use_config_rse / auto-increment behavior unchanged.
        bool m_rse_from_ident{false};

        // Like m_rse_from_ident for the EVENT number, but the configured run
        // and subrun are kept.  A group of events streamed through one process
        // can span several runs (SBND nueCC48 spans 12), so a per-group job
        // supplies the run/subrun of each event in m_rse_map; an ident with no
        // entry keeps the configured pair.  Default off => nothing changes.
        bool m_event_from_ident{false};
        // ident -> (run, subrun), consulted only when m_event_from_ident.
        std::map<int, std::pair<int,int>> m_rse_map;

        // Own (non-shared) Bee zip.  A printf conversion in the configured
        // "bee_zip" means one zip per event: the open zip is closed and the
        // next opened when the event number changes.  Without one the zip is
        // opened in configure() exactly as before.
        std::string m_bee_zip{"mabc.zip"};
        bool m_bee_zip_templated{false};
        /** doc 87: an EMPTY "bee_zip" (and no "bee_sink") means write no Bee
            zip at all.  Bee::Sink::reset("") would raise IOError -- an empty
            name was never a disable -- so the sink is simply never built and
            write_obj() becomes a no-op.  Default "mabc.zip" => unchanged. */
        bool m_bee_disabled{false};
        int m_bee_zip_open_evt{-1};
        /// Open the own-sink zip for the current event, if that is still needed.
        void ensure_own_sink();
        /// Set m_eventNo from the tensor ident, and run/subrun from m_rse_map.
        void apply_event_from_ident(int ident);

        void flush(int ident = -1);
        void flush(WireCell::Bee::Points& bpts, int ident);

        // Write one Bee object: routes to the shared sink (at the current
        // per-event index) when m_use_shared_sink, else to the per-node m_sink.
        size_t write_obj(const WireCell::Bee::Object& obj);

        bool m_save_deadarea{false};
        // save_real_cluster_id (default false = byte-identical legacy tarball):
        // before serializing the point-cloud trees, give EVERY cluster that
        // carries a "perblob" PC the "real_cluster_id" / "real_cluster_main"
        // arrays (fill-in: own ident / own "main_cluster" flag).  The tensor
        // serializer (TensorDM as_tensors) concatenates same-named local PCs
        // across nodes and silently drops arrays whose key is absent from the
        // first-seen node, so the flash-merge provenance written by
        // examine_bundles (real_cluster_id on merged clusters only) never
        // survived the tarball.  Homogenizing the key set at save time lets it
        // through; the fill-in values reproduce exactly what a reader would
        // assume for an unmerged cluster, so nothing downstream changes except
        // that the arrays now exist after a load.
        //
        // SCOPE INVARIANT: "real_cluster_id" is meaningful only WITHIN one
        // cluster -- see the merge_clusters() comment in ClusteringFuncs.h.
        // The recorded values come from the ident numbering in force when
        // examine_bundles ran; the fill-in below uses the CURRENT numbering
        // (enumerate_idents re-runs after every visitor).  The two epochs
        // overlap, so joining on this value across clusters is wrong.
        bool m_save_real_cluster_id{false};
        // real_cluster_id_global (default TRUE since doc 53, owner decision):
        // re-stamp "real_cluster_id" at save time into ONE globally unique
        // epoch -- the representative rows (real_cluster_main != 0) take the
        // cluster's own current ident, every other pre-merge group takes a
        // fresh id above the largest ident in the grouping.  Without it the
        // array mixes the ident numbering merge_clusters recorded with the one
        // enumerate_idents has installed since, and 31% of values name two
        // clusters (SBND d52ron 30-event set).
        //
        // Group membership is untouched, so consumers that only compare rows
        // within a cluster (ClusteringUnmergeBundle, TaggerCheckTGM
        // main_component_mode="real") are behaviourally unchanged -- measured:
        // partition and real_cluster_main byte-equal over 179 clusters, nusel
        // verdict tables row-for-row identical.  What changes is that the value
        // becomes a valid event-wide key.  A cluster that was never merged is
        // rewritten to exactly the values it already had.
        //
        // SCOPE: restamp_real_cluster_id() runs once per event right after the
        // clustering pipeline and BEFORE the Bee fills, so the Bee per-blob
        // label and the saved pctree carry the SAME ids.  Per-visitor Bee dumps
        // (trace_bee) are mid-pipeline snapshots and necessarily keep whatever
        // ids existed at that step.
        //
        // NOT gated on save_real_cluster_id: examine_bundles writes the array in
        // memory whether or not the tarball is saved, so the Bee labels are
        // affected either way.  Detectors that never write it (PDHD, PDVD --
        // their examine_bundles is disabled) see a structural no-op.
        //
        // Gated by save_real_cluster_id, so it is a STRUCTURAL no-op wherever
        // that is off -- i.e. every detector but SBND.  Set false only to
        // reproduce the two-epoch values for A/B archaeology.
        bool m_real_cluster_id_global{true};
        // save_assoc_cluster_id (default false = byte-identical legacy tarball):
        // the same homogenization for the isolated grouping's provenance pair
        // "assoc_cluster_id" / "assoc_cluster_main" (doc 52 Stage 1/2), written
        // by clustering_isolated's save_assoc_id and carried across merges by
        // merge_clusters.  Fill-in for a cluster the grouping never touched: own
        // ident, and main=1 -- an ungrouped cluster IS a main, which is the
        // sentinel the un-merge relies on to keep cathode crossers whole.
        bool m_save_assoc_cluster_id{false};
        // 1 = legacy bare-array channel-deadarea-*.json (default; back-compat for
        //     single-TPC viewers like the original Bee).  2 = wire-cell-bee3 v2
        //     wrapper {"version":2,"tpc":<apa>,"polygons":[...]} that places the
        //     dead-area slab on the per-TPC anode face.  See wire-cell-bee3
        //     docs/dead-area.md.  Currently we use the WCT anode ident as the
        //     bee TPC index, which is correct for single-face anode detectors
        //     (e.g. SBND) and may need a mapping table for multi-face anodes.
        int m_dead_area_version{1};


        // Count how many times we are called
        size_t m_count{0};

        /** Config: "groupings"
         *
         * List of groupings to select for processing.
         *
         * Default: ["live","dead"].
         */
        std::vector<std::string> m_groupings = {"live","dead"};

        /** Config: "inpath"
         *
         * The BASE datapath for the input pc tree data.  This may be a regular
         * expression which will be applied in a first-match basis against the
         * input tensor datapaths.  If the matched tensor is a pcdataset it is
         * interpreted as providing the nodes dataset.  Otherwise the matched
         * tensor must be a pcgraph.
         *
         * A "%d" will be interpolated with the ident number of the input tensor
         * set.
         * 
         * See also "insubpaths".
         */
        std::string m_inpath{".*"};

        /** Config: "outpath"
         *
         * The BASE datapath for the resulting pc tree data.
         *
         * A "%d" will be interpolated with the ident number of the input tensor
         * set.
         *
         * See outsubpaths.
         */
        std::string m_outpath{""};

        /** Config: insubpaths, outsubpaths.
         *
         * By default, a grouping of a given NAME is located at an input or
         * output path spelled as: "{inpath,outpath}/NAME".
         *
         * If a grouping NAME is found in either insubpath or outsubpath then
         * this default is overridden.  Both parameters are array of objects,
         * each object has keys "name" and "subpath".  The subpath is a simple
         * string suffix and thus should include a leading "/" if the user
         * wishes to locate the grouping in a "subdirectory".
         */
        // See issue #375 and #416.  
        std::map<std::string, std::string> m_insubpaths, m_outsubpaths;
        std::string inpath(const std::string& name, int ident);
        std::string outpath(const std::string& name, int ident);

        /** Config: "perf"
         *
         * If true, emit time/memory performance measures.  Default is false.
         */
        bool m_perf{false};

        /** Config: "grouping2file_prefix"
         *
         * If not empty, dump the final grouping to a file with this prefix + potential RSE info + .npz.
         */

        std::string m_grouping2file_prefix{};

        /** Config: "dump_json"
         *
         * If true, dump to files like {live,dead}-summary-<ident>.json a
         * summary of the groupings.  Default is false.  The dumps are rather
         * large despite omitting point cloud point data.  Use of jq or other
         * tool is expected.  These are intended only for debugging / validating.
         */
        bool m_dump_json{false};

        /** Config: "cluster_id_order"
         *
         * The various operations can lead to redundantly or non-sequentially
         * numbered cluster idents.  The application of merge_clusters() will
         * cause a resulting cluster to have the ID of its first contributing
         * constituent.  The application of separate() will cause all clusters
         * to have the ID of the original.
         *
         * When this parameter is given, the cluster IDs will be reset after
         * each component operation, on a per-grouping basis, so that they
         * represent a specific sort order.
         *
         * - "tree" :: use the as-is child-order which represents
         *   child-insertion order unless some operation has explicitly
         *   reordered the underlying PC tree.
         *
         * - "size" :: use the "size" of the cluster as determined by
         *   cluster_less() to order the cluster IDs.  This considers the
         *   cluster length, number of children, number of points followed by
         *   per-view min then max bounds and finally cluster center.
         *
         * - "" :: cluster IDs are not modified.
         *
         * By default, no ID rewriting is performed.
         *
         * Notes:
         *
         * - When an ID cluster order is applied, the ID counting starts from 1.
         * - The default (unset) ID is -1.
         */
        std::string m_clusters_id_order;

        // configurable parameters for dead-live clustering
        int m_dead_live_overlap_offset{2};

        /** Config: "ctpc_aniso_metric"
         *
         * doc pdvd/36.  If true, every grouping this node loads answers its
         * ctpc radius queries (Grouping::get_closest_points /
         * has_closest_point, hence every good-point test, connector and
         * charge average on it) with the lattice-normalised anisotropic
         * metric of CtpcAnisoMetric.h instead of the isotropic Euclidean one.
         * Job-wide on purpose: the callers must agree on the metric.
         * Default false = legacy path, byte-identical.  Even when on it is the
         * identity on SBND (pitch finer than the drift step, s clamps to 1);
         * it loosens the pitch axis on PDVD (U/V s = 0.387, W 0.581), PDHD
         * (~0.67) and uBooNE (0.73), so it must stay OFF outside PDVD.
         */
        bool m_ctpc_aniso_metric{false};

        // Keep track of configured clustering methods with their metadata to
        // assist in debugging/logging.
        struct EnsembleVisitor {
            std::string name;
            IEnsembleVisitor::pointer meth;
        };
        /** Config: pipeline
         *
         * Array of type/name of instances of IEnsembleVisitor to execute in the pipeline.
         */
        std::vector<EnsembleVisitor> m_pipeline;

        // the anode to be processed
        std::vector<IAnodePlane::pointer> m_anodes;

        IDetectorVolumes::pointer m_dv;

        // the face to be processed
        // int m_face{0};

        // the geometry helper
        // IClusGeomHelper::pointer m_geomhelper;
    };
}  // namespace WireCell::Clus

#endif
