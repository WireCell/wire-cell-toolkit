#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/ParticleDataSet.h"
#include "WireCellClus/FiducialUtils.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/Logger.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/TrackFitting.h"  
#include "WireCellClus/TrackFittingPresets.h"
#include "WireCellClus/PRSegmentFunctions.h"

#include "WireCellIface/IScalarFunction.h"
#include "WireCellUtil/KSTest.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

class TaggerCheckNeutrino : public Aux::Logger, public IConfigurable, public Clus::IEnsembleVisitor, private Clus::NeedDV, private Clus::NeedPCTS, private Clus::NeedRecombModel, private Clus::NeedParticleData, private Clus::NeedClusGeomHelper, private Clus::NeedFiducial {
public:
    TaggerCheckNeutrino() : Aux::Logger("TaggerCheckNeutrino", "clus") {
        // Initialize with default preset
        m_track_fitter = std::make_shared<TrackFitting>(TrackFittingPresets::create_with_current_values());
    }
    virtual ~TaggerCheckNeutrino() {}
    virtual void configure(const WireCell::Configuration& config) ;

    virtual Configuration default_configuration() const ;

    virtual void visit(Ensemble& ensemble) const;


    private:
        std::string m_grouping_name{"live"};
        std::string m_trackfitting_config_file;  // Path to TrackFitting config file
        bool m_perf{false};  // if true, print per-step timing to stdout
        // Route direction-weakness reads through segment_is_dir_weak() (the
        // faithful port of prototype ProtoSegment::is_dir_weak(), score
        // thresholds included).  C++ default false = legacy raw-flag reads,
        // byte-identical.  See PatternAlgorithms::m_dir_weak_use_score.
        bool m_dir_weak_use_score{false};
        // MIP dQ/dx scales in e/cm (uBooNE legacy defaults => absent keys are
        // byte-identical).  Two roles: mip_dqdx = flat-template amplitude
        // (legacy 50e3); mip_dqdx_median = median-ratio threshold scale
        // (legacy 43e3).  See PatternAlgorithms::m_mip_dqdx{,_median}.
        double m_mip_dqdx{50000.0};
        double m_mip_dqdx_median{43000.0};
        // Proton-template direction vote (doc sbnd_xin/docs/pr/8; default
        // false = legacy abstention).  Thresholds pending pr/8 calibration.
        bool   m_proton_dir_vote{false};
        double m_proton_dir_score_max{0.25};
        double m_proton_dir_asym_min{1.3};
        // Endpoint-trim retry (doc sbnd_xin/docs/pr/9 sec. 6 F1; default false
        // = legacy).  On abstention only, retry the dQ/dx PID once with 1
        // sample excluded at the hypothesized stopping end (ill-defined
        // endpoints dilute or inflate the tip dQ/dx).  Runs before the vote.
        bool   m_endpoint_trim_retry{false};
        // Minimum wcpt-path length (cm) for a segment to count as a leg in the
        // fit_vertex position fit (doc sbnd_xin/docs/pr/9 sec. 11 F3c; default
        // 0 = legacy include-all, byte-identical).  Stops vertex-activity
        // stubs from dragging the vertex: >=3 surviving legs => fit on the
        // survivors; <=2 => skip the fit (two-leg position already fit).
        double m_fit_vertex_min_seg_length{0};
        // Robust vertex fit (doc sbnd_xin/docs/pr/51 round 7; design at
        // NeutrinoPatternBase.h m_mvfit_robust and MyFCN.h RobustParams).
        // Lengths in cm here (cm -> internal at the push block); angle deg,
        // frac/aniso/pts unitless.  Default false = legacy byte-identical.
        bool   m_mvfit_robust{false};
        bool   m_mvfit_main_only{true};
        double m_mvfit_min_len{10.0};      // cm
        double m_mvfit_rin_margin{2.0};    // cm
        double m_mvfit_rout_frac{0.5};
        double m_mvfit_rout_min{9.0};      // cm
        double m_mvfit_rout_max{18.0};     // cm
        double m_mvfit_angle{20.0};        // deg
        int    m_mvfit_min_pts{5};
        double m_mvfit_min_aniso{3.0};
        double m_mvfit_prior_range{1.0};   // cm
        // Cathode kink veto (doc sbnd_xin/docs/pr/20 Part II B0), both in cm.
        // cathode_kink_xcut = 0 => segment_search_kink sees every fit point =>
        // byte-identical to the pre-pr/20 behavior.  cathode_x is the cathode
        // plane in the T0-corrected frame (SBND 0; the same convention
        // ClusteringCathodeConnect's cathode_x uses).
        double m_cathode_x{0};
        double m_cathode_kink_xcut{0};
        // Wide-baseline cathode kink accept (doc sbnd_xin/docs/pr/47 sec 8,
        // O1): fifth segment_search_kink accept path at cathode-crossing fit
        // indices, keyed on the skirt-excluded PCA turn angle across the
        // crossing (angle in deg; skirt/baseline in cm).  angle = 0 => path
        // never evaluated => byte-identical legacy kink search.
        double m_cathode_wide_kink_angle{0};     // deg
        double m_cathode_wide_kink_skirt{3};     // cm
        double m_cathode_wide_kink_baseline{15}; // cm
        // ---- doc sbnd_xin/docs/pr/48: back-to-back track fixes -------------
        // Two-end residual-range break (18255-51513/56211/57903/57485): a
        // single-non-stub-segment main cluster, both endpoints in the FV,
        // with dQ/dx rising at BOTH ends gets broken at the argmin of the
        // joint two-arm stopping-template score.  false => the pass never
        // runs => byte-identical.  Full rationale and per-knob meaning:
        // PatternAlgorithms::m_two_end_break / PR::TwoEndBreakOptions.
        // Lengths cm, angle deg, ratios/caps dimensionless.
        bool   m_two_end_break{false};
        double m_teb_min_len{10};        // cm
        double m_teb_min_arm{1.8};       // cm
        int    m_teb_min_arm_pts{4};
        double m_teb_stub_max{4};        // cm
        double m_teb_accept_range{15};   // cm
        double m_teb_rise_r1{1.3};
        double m_teb_rise_r2{1.15};
        double m_teb_abs_end_min{1.7};   // x mip_dqdx_median
        double m_teb_dip_floor{0.6};     // x mip_dqdx_median; dips below = instrumental
        double m_teb_score_cap_r1{0.6};
        double m_teb_score_cap_r2{0.9};
        double m_teb_turn_angle{25};     // deg; <= 0 disables route R2
        double m_teb_turn_baseline{35};  // cm
        double m_teb_turn_skirt{3};      // cm
        // doc sbnd_xin/docs/pr/90 round 2.  turn_min_arm_frac: route R2's
        // turn argmax PREFERS an index whose PCA arms can each span this
        // fraction of teb_turn_baseline, when such an index clears
        // teb_turn_angle on its own; otherwise the legacy argmax (starved
        // near-end candidates included) stands.  0 = legacy.
        // second_max (cm): entry gate tolerates extra long prongs when
        // exactly one segment exceeds this cap; 0 = legacy strict gate.
        double m_teb_turn_min_arm_frac{0.0};
        double m_teb_second_max{0};      // cm
        // doc pr/90 round 4 (sec 9.5 D1/D3/D4); all default OFF =>
        // byte-identical.  See the PatternAlgorithms member block.
        bool   m_teb_chain_topology{false};
        double m_teb_r3_turn{0.0};       // deg; <= 0 disables route R3
        double m_teb_r3_hot{0.0};        // x mip_dqdx_median; <= 0 disables route R3
        double m_teb_bragg_veto_turn{0.0}; // deg; <= 0 disables the R2 bragg veto
        // 59335 fix (a): the local-dQ/dx walk gate also stops the C4 /
        // straightness (flag_search) accepts.  false => byte-identical.
        bool   m_kink_walk_dqdx_stop{false};
        // Shared Bragg-hot ratio (x mip_dqdx_median) gating BOTH 59335
        // fixes to genuinely hot kinks/breaks.  Inert while both are false.
        double m_kink_dqdx_hot_ratio{1.7};
        // 59335 fix (b): a break born from a C4/A0 accept gets
        // VertexFlags::kProtectedBreak so examine_vertices_4's < 2 cm
        // absorption floor cannot erase it.  false => byte-identical.
        bool   m_kink_break_protect{false};
        // doc sbnd_xin/docs/pr/50 -- main-vertex kink-consistency snap
        // (172230-class near-vertex robustness).  Mirrors of the
        // PatternAlgorithms::m_vertex_kink_snap / m_vks_* members (see the
        // design comment there); numerics here in cm/deg, converted at the
        // visit() copy.  All defaults = the pass never fires =>
        // byte-identical.
        bool   m_vertex_kink_snap{false};
        double m_vks_radius{5.0};     // cm
        double m_vks_min_dis{0.5};    // cm
        double m_vks_angle{25.0};     // deg
        double m_vks_margin{10.0};    // deg
        double m_vks_collinear{30.0}; // deg
        double m_vks_skirt{0.3};      // cm
        double m_vks_baseline{2.0};   // cm
        double m_vks_min_arm{1.5};    // cm
        double m_vks_fit_miss{0.35};  // cm; snap only when the fit misses the image corner by >= this
        double m_vks_hot_ratio{0};    // x mip_dqdx_median; 0 = veto off (default: misfires on misassigned charge)
        double m_vks_carry_prong{0};  // cm; carry the old vertex's arms through the snap residual below this arc (doc pr/85); 0 = off, byte-identical
        // doc sbnd_xin/docs/pr/51 -- main-vertex graph audit (near-vertex
        // graph-shape repair: duplicate-corridor merge / charge-less-bridge
        // removal / micro-stub absorb + re-seat / one refit).  Mirrors of
        // the PatternAlgorithms::m_main_vertex_graph_audit / m_mvga_*
        // members (see the design comment there); numerics here in cm/deg,
        // converted at the visit() copy.  All defaults = the pass never
        // fires => byte-identical.
        // sbnd_xin/docs/73 sec 12 (round 3): mirror of
        // PatternAlgorithms::m_esva_ignore_empty_2d (see the design comment
        // there -- the empty-2D-index sentinel in
        // eliminate_short_vertex_activities case 5 must not read as
        // "covered" on cathode-crossing clusters).  false = legacy.
        bool   m_esva_ignore_empty_2d{false};
        bool   m_main_vertex_graph_audit{false};
        double m_mvga_radius{15.0};       // cm
        double m_mvga_dup_tol{1.4};       // cm
        double m_mvga_dup_frac{0.7};      // fraction
        double m_mvga_dup_angle{20.0};    // deg; <= 0 disables the op1 near-parallel guard
        double m_mvga_bridge_mip{0.5};    // x mip_dqdx_median (measured: 268067 bridge 0.436, real track 1.29)
        double m_mvga_reconnect{5.0};     // cm
        double m_mvga_stub{2.0};          // cm
        int    m_mvga_stub_pts{4};        // valid fit points
        double m_mvga_reseat_angle{150.0}; // deg
        double m_mvga_satellite{0};       // cm; 0 = main-vertex-only op3 scope (round 2), byte-identical (round 3, doc pr/51)
        bool   m_mvga_interposed{false};  // op3 interposed-stub absorb at the main-vertex anchor (doc pr/85); false = terminal-only, byte-identical
        double m_mvga_interposed_angle{150.0}; // deg; far-end collinearity gate for the interposed absorb
        double m_mvga_interposed_len{0.0};  // cm; interposed-splice candidate ceiling (doc pr/86); 0 = use mvga_stub, byte-identical
        double m_mvga_sat_dup_frac{0.0};    // fraction; satellite-anchor op3 overlap threshold (doc pr/86); 0 = use mvga_dup_frac, byte-identical
        bool   m_mvga_interposed_deg1{false}; // op3 interposed splice at degree-1 main anchors (doc pr/86); false = byte-identical
        double m_mvga_splice_straighten{0.0}; // cm; op3 post-carry straighten reach past the junction (doc pr/86 round 2); 0 = concatenation verbatim, byte-identical
        double m_mvga_approach_collapse{0.0}; // cm; op3.5 junction-collapse radius around the main vertex (doc pr/86 round 2); 0 = pass skipped, byte-identical
        double m_mvga_straighten_radius{0.0}; // cm; R1/R2 straight-chain charge-veto radius (doc pr/86 round 2); 0 = prototype 0.2 cm; inert unless straighten/collapse on
        double m_mvga_op1_radius{0.0};    // cm; op1-only scope radius (doc pr/83 r3); 0 = use mvga_radius, -1 = unscoped, byte-identical at 0
        double m_mvga_op1_dup_frac{0.0};  // fraction; op1-only overlap threshold (doc pr/83 r3); 0 = use mvga_dup_frac, byte-identical
        bool   m_mvga_op1_post{false};    // post-op3 duplicate-corridor pass incl. created segments (doc pr/83 r3 class A); false = byte-identical
        int    m_mvga_carry_max{0};       // op3 interposed-carry prong-count ceiling (doc pr/83 r3); 0 = unlimited, byte-identical
        bool   m_swap_orphan_dup_audit{false}; // dup-audit the abandoned main cluster inside swap_main_cluster (doc pr/83 r3 Mechanism C); false = byte-identical
        double m_mvga_proj_dup_frac{0.0};  // 2nd-best per-view overlap threshold for the projective dup collapse (doc pr/83 r4); 0 = disabled, byte-identical
        double m_mvga_proj_dqdx_ratio{0.4}; // stem dQ/dx asymmetry gate for the same pass (doc pr/83 r4); inert while frac == 0
        double m_mvga_proj_angle{0.0};    // deg; op1-proj chord-angle ceiling (doc pr/83 r4b); 0 = use mvga_dup_angle, byte-identical
        // Long shower-topology demote length, cm (doc sbnd_xin/docs/pr/25
        // sec 3).  0 => the guard never fires => byte-identical.  50 is the
        // scan-supported operating point (9/10 owner-scanned events; ~45
        // covers all 10).  Threaded to PatternAlgorithms and thence to every
        // segment_is_shower_topology call site -- both NeutrinoTrackShowerSep
        // and NeutrinoVertexFinder, so one segment cannot be classified two
        // ways within one event.
        double m_shower_topo_demote_len{0};

        // doc sbnd_xin/docs/pr/49 -- cross-cluster projection-ghost
        // deweighting of the trajectory fit's 2D charge association
        // (18255-57441 V-plane ghost).  A live cell outside the fitted
        // cluster's own blob coverage that sits inside a 3D-distant foreign
        // cluster's keeps its measurement at reduced weight (deweight, not
        // drop -- dead-channel single-view charge must stay usable).
        // Pushed to TrackFitting::Parameters::fit_blob_coverage per visit;
        // -1 = off (legacy, byte-identical), >= 0 = on with value =
        // wire/slice tolerance in cells (0 = strict).  NOTE: this jsonnet
        // key is the single source of truth -- it overrides any
        // trackfitting_config_file value at visit time.
        double m_fit_blob_coverage{-1};

        // doc sbnd_xin/docs/pr/50 (172230-class near-vertex robustness):
        // when true AND fit_blob_coverage >= 0, the deweighting is
        // SUSPENDED while find_proto_vertex forms the MAIN cluster's break
        // partition (the recursive kink walk refits after every break and
        // re-reads the fits, so any fit perturbation -- however far from
        // the vertex -- can reshuffle the whole partition and lose the
        // proto-vertex at the true kink; 18255-172230 lost its main vertex
        // to a 2.7 cm neighbor this way).  The partition then forms on
        // legacy fits, byte-identical to knob-off by construction, while
        // every later fitting stage (clustering_points onward, main-vertex
        // determination, improve_vertex, the final trajectory + dQ/dx)
        // keeps the pr/49 ghost protection.  Main cluster ONLY: a
        // non-main cluster's final trajectory is essentially its
        // find_proto_vertex fit (no later full refit), so deferring there
        // un-fixes the pr/49 ghosts (57441 cid 20 measured 1.12 -> 1.23 cm
        // under a global defer).  Default false = pr/49 behavior,
        // byte-identical.  Inert while fit_blob_coverage < 0.
        bool m_fit_blob_coverage_defer{false};

        // ---- doc sbnd_xin/docs/pr/30 §11 port-fidelity knobs ----------------
        // Threaded verbatim to PatternAlgorithms / PR::g_graph_endpoint_policy.
        // See PatternAlgorithms::m_fit_exclusion etc. for the full rationale.
        // The first three default to today's behaviour by being FALSE (new
        // behaviour is opt-in); the last two default to TRUE because the
        // behaviour they gate is already production.  Either way, defaults =>
        // byte-identical to the pre-pr/30 tree.
        bool   m_fit_exclusion{false};             // P1
        bool   m_graph_endpoint_strict{false};     // P8 (the WARN is unconditional)
        double m_graph_endpoint_tol{0.3};          // cm
        bool   m_oov_prototype_parity{false};      // F2 (all three sites)
        bool   m_first_seg_local_pca{true};        // P2
        bool   m_other_seg_relaxed_accept{true};   // P4
        bool   m_other_seg_empty_2d_guard{false};  // doc pr/45 -- -1 empty-2D-tree sentinel guard
        // doc sbnd_xin/docs/pr/54 -- keep well-supported isolated residual
        // segments instead of discarding them after the fit (18255-142421
        // separated EM shower).  All default to the legacy discard.
        bool   m_other_seg_keep_isolated{false};
        int    m_other_seg_keep_isolated_min_points{25};
        double m_other_seg_keep_isolated_min_length{3.0}; // cm; scaled at copy
        // doc sbnd_xin/docs/pr/67 round 3 (S2) -- size gate on the isochronous
        // snap in find_other_segments.  Legacy 10 cm; lowering it lets a short
        // isochronously-displaced branch reach modify_vertex/segment_isochronous.
        double m_iso_snap_min_dir_mag{10.0}; // cm; scaled at copy
        // doc sbnd_xin/docs/pr/59 round 2 -- when a segment created after its
        // cluster's association pass would otherwise be left with a null
        // associate_points cloud (18255-142421 seg 20; 116944-71372 segs
        // 19052/19053/136199), clear and re-run the WHOLE cluster's
        // association + re-classify just the rescued segments.  Default false
        // = legacy (orphan stays null).  See
        // PatternAlgorithms::reassociate_cluster_orphans.
        bool   m_assoc_full_recluster{false};
        // doc sbnd_xin/docs/pr/64 round 7 -- reassign points Stage C of
        // clustering_points_segments would otherwise drop to the same-cluster
        // segment that actually wins the global 2D projection contest,
        // instead of discarding them.  Default false = legacy (drop).  See
        // PatternAlgorithms::m_assoc_reassign_orphans.
        bool   m_assoc_reassign_orphans{false};
        // doc sbnd_xin/docs/pr/64 round 8 -- clears a merge survivor's
        // associate_points when examine_structure_final_1/_1p/_3 deletes a
        // segment that had non-empty associate_points, so pr/59's
        // reassociate_cluster_orphans any_orphan trigger correctly re-fires.
        // Default false = legacy (survivor keeps its stale cloud). See
        // PatternAlgorithms::m_assoc_clear_on_merge.
        bool   m_assoc_clear_on_merge{false};
        // doc sbnd_xin/docs/pr/31 §11 -- F2 (was P2).  true => skip the
        // stage-3 segment_determine_shower_direction call, so a topology
        // shower keeps the direction segment_is_shower_topology set, which is
        // what the prototype does.  false => today's path, byte-identical.
        // Rationale, the mutation audit and the F3 residual: see
        // NeutrinoPatternBase.h's m_shower_topo_proto_dir.
        bool   m_shower_topo_proto_dir{false};     // pr/31 F2
        // doc sbnd_xin/docs/pr/32 §11 -- the four kept findings of the stage-4
        // (neutrino vertex identification) audit.  All C++ default false =
        // today's path = byte-identical; the SBND operating point turns them
        // on.  Rationale and prototype anchors: NeutrinoPatternBase.h.
        bool   m_vertex_dir_use_fit_point{false};       // pr/32 F1 (was P1)
        bool   m_shower_traj_recheck_parity{false};     // pr/32 F2 (was P3)
        bool   m_main_vertex_require_descriptor{false}; // pr/32 F3 (was P7)
        bool   m_main_vertex_candidate_flag{false};     // pr/32 F4 (was P12)
        // doc sbnd_xin/docs/pr/31 §12 -- the §10.12 port-fidelity round
        // (topology/PID/direction audit, owner-filtered survivors).  All C++
        // default false = today's path = byte-identical; the SBND operating
        // point decides which are on.  Rationale and prototype anchors:
        // NeutrinoPatternBase.h.
        bool   m_cont_muon_dir3_30cm{false};            // pr/31 F5 (was P6)
        bool   m_track_comp_empty_abstain{false};       // pr/31 F6 (was P7)
        bool   m_shower_topo_reset{false};              // pr/31 F3 (was P13)
        bool   m_reclass_preserve_4mom{false};          // pr/31 F1 (was P1+P3a+P4)
        bool   m_reclass_never_computed_ke_floor{false}; // pr/40 round 2 F6
        bool   m_dir_track_median_local{false};         // pr/31 F4 (was P8)
        bool   m_examine_showers_vertex_by_index{false}; // pr/31 F7 (was P5) -- stays OFF pending pr/30 F4
        // Isochronous first-segment endpoint finding (doc sbnd_xin/docs/pr/24
        // round 2, SBND evt 271851); lengths in cm.  false => legacy
        // wire-footprint boundary endpoints, byte-identical.  When on, a
        // long cluster whose quantile-trimmed drift-x extent is small (a
        // filled 2-D sheet) gets its first-segment endpoints from the sheet's
        // principal axis instead, and the local-PCA endpoint refinement is
        // skipped on that branch.
        bool   m_iso_endpoint{false};
        double m_iso_endpoint_min_length{40};   // cm
        double m_iso_endpoint_max_xext{25};     // cm
        double m_iso_endpoint_xext_frac{0.35};
        double m_iso_endpoint_xext_quantile{0.02};
        double m_iso_endpoint_tube_radius{4};   // cm
        double m_iso_endpoint_min_aspect{0.12};
        // Round 5 (doc sbnd_xin/docs/pr/24 sec. 18).  examine_vertices_3's
        // get_local_extension recovery step can retract a track endpoint
        // instead of extending it (worst at an isochronous sheet's axial
        // extreme, where m_iso_endpoint picks its seed).  false => today's
        // unconditional accept, byte-identical.  m_v3_extension_min_gain (cm)
        // is the minimum required increase in distance to the segment's far
        // endpoint; a small negative default tolerates legacy's few-mm
        // legitimate retreat.
        // doc sbnd_xin/docs/pr/67: log-only trajectory-coverage probe, and the
        // counterfactual override for find_proto_vertex's hardcoded
        // branch-search round budget.  false / 0 => byte-identical.
        bool   m_traj_cover_probe{false};
        int    m_pr_find_other_rounds{0};
        bool   m_v3_extension_guard{false};
        double m_v3_extension_min_gain{-1.0};   // cm
        // doc sbnd_xin/docs/pr/72 round 2 -- guard examine_structure_3
        // against merging a genuine near-vertex track stub into an
        // unrelated shower/track trunk (18255-196649).  Default false =
        // legacy (unconditional merge on angle alone), byte-identical.
        // Numeric defaults fitted from a 117-event merge census; see
        // PatternAlgorithms::m_es3_stub_guard / es3_stub_suppress
        // (PRSegmentFunctions.h) for the predicate and the doc for the fit.
        bool   m_es3_stub_guard{false};
        double m_es3sg_stub_max{7.0};      // cm; scaled at copy
        double m_es3sg_len_ratio{2.0};
        double m_es3sg_ang3_min{15.0};     // degrees
        double m_es3sg_ang_ratio{1.0};
        bool   m_es3sg_require_terminal{true};
        // Detector-extent literals (docs/pr/2 sec. 2e(iv)), all in cm.
        // Defaults = the uBooNE prototype values (active volume y in
        // [-116,+117], z in [0,1037]) => absent keys are byte-identical.
        // cosmic_tagger()'s four "reaches the top" cuts, quoted as a distance
        // below the top face (117 - value): 17 / 15 / 37 / 67 cm.
        double m_cosmic_y_top_main{100};
        double m_cosmic_y_top_strict{102};
        double m_cosmic_y_top_loose{80};
        double m_cosmic_y_small_piece{50};
        // Denominator of the upstream-z vertex-score penalty (z-min_z)/scale in
        // compare_main_vertices{,_global}.  See PatternAlgorithms.
        double m_vertex_z_prior_scale{200};
        // SSM beam-line reference directions in the detector frame, {x,y,z}.
        // Defaults = the prototype's uBooNE BNB-target / NuMI-absorber vectors
        // (absent keys => byte-identical).  SBND has no value for either yet
        // (docs/pr/2 sec. 2e(i)).  See PatternAlgorithms::m_ssm_target_dir.
        std::vector<double> m_ssm_target_dir{0.46, 0.05, 0.885};
        std::vector<double> m_ssm_absorber_dir{0.33, 0.75, -0.59};
        // Charge -> kinetic-energy calibration constants of NeutrinoEnergyReco
        // (docs/pr/2 sec. 2e(iii)).  Defaults = the uBooNE literals they
        // replaced => absent keys are byte-identical.  Dimensionless except
        // m_kine_w_value (eV).  See PR::KineChargeOptions.
        double m_kine_fudge_factor{0.95};
        double m_kine_recom_factor{0.7};
        double m_kine_shower_fudge_factor{0.8};
        double m_kine_shower_recom_factor{0.5};
        double m_kine_proton_recom_factor{0.35};
        std::vector<double> m_kine_plane_weights{0.25, 0.25, 1.0};  // {U,V,W}
        double m_kine_plane_asym_switch{0.04};
        double m_kine_w_value{23.6};  // eV per electron-ion pair
        // doc pr/35 §10.2 (F1): read the shower PDG live from the start
        // segment at the fill_kine_tree sites instead of the cached field.
        bool m_kine_shower_pdg_live{false};
        // Muon median-dQ/dx-vs-length envelope {c0, c1, pivot_cm, power}:
        // cut = c0 + c1*pow(pivot/L, power), a multiple of mip_dqdx_median,
        // used at nine tagger sites (docs/pr/2 sec. 2e(iv)).  Defaults = the
        // prototype's empirical uBooNE stopping-muon refit => absent key is
        // byte-identical.  See PatternAlgorithms::m_muon_dqdx_curve.
        std::vector<double> m_muon_dqdx_curve{0.8866, 0.9533, 18.0, 0.4234};
        // Single-photon stem dE/dx conversion (docs/pr/2 sec. 2e(i) third
        // correctness item).  Default false = the inline uBooNE-field inverse
        // Modified-Box, byte-identical; true routes shw_sp_vec_*_dedx through
        // the configured recombination model.  The mean-dedx threshold
        // (MeV/cm) is coupled to that choice; default = the legacy 2.3.
        bool   m_sp_dedx_use_recomb_model{false};
        double m_sp_mean_dedx_cut{2.3};
        std::string m_dl_weights;              // path to SCN vertex .pth weights file (empty = DL disabled)
        // doc sbnd_xin/docs/pr/75: true iff dl_weights was configured non-empty
        // but Persist::resolve failed.  A failed resolve empties m_dl_weights,
        // so without this the scoreboard cannot tell doc pr/52 route 3
        // ("weights path not found") from route 1 ("DL never configured").
        bool m_dl_weights_missing{false};
        double m_dl_vtx_cut{25.0};             // max distance (mm) from DL prediction to accept candidate vertex (default 2.5 cm)
        double m_dQdx_scale{0.1};              // scale factor applied to dQ before passing to SCN network
        double m_dQdx_offset{-1000.0};         // offset applied after scaling: q_in = dQ * scale + offset
        bool   m_dl_vtx_rerank{true};              // if true, use top-K + soft re-rank; if false, use legacy single-best argmax
        int    m_dl_vtx_top_k{5};                // number of top voxels from DL inference to re-rank (only when rerank enabled)
        double m_dl_vtx_min_accept_score{4.0};    // minimum composite re-rank score to accept DL vertex (only when rerank enabled; matches the 4.0 advertised in default_configuration -- docs/pr/2 sec 8.4)
        double m_dl_vtx_score_scale{1000.0};      // scale factor applied to the raw DL score term in the composite re-rank score
        // doc sbnd_xin/docs/pr/51 (18255-506746): rerank-branch guard --
        // skip DL candidates hosted on a different cluster than the current
        // main cluster, so an accepted DL vertex can never swap the main
        // cluster (a confident wrong voxel's s_dl swamps every structural
        // term).  false = legacy = byte-identical.
        bool   m_dl_vtx_swap_guard{false};
        // doc sbnd_xin/docs/pr/89 Arm C (C2): rule-1 outgoing-prong topology
        // term in the DL rerank composite, s_topo = w * (frac - center) for
        // candidates with >= 1 decisive attached Bragg vote.  0 = legacy =
        // byte-identical (the term is never computed).  The offline C1
        // replay selected w = 3.0, center 0 (pr/89 sec 11.5).
        double m_dl_vtx_topo_weight{0.0};
        double m_dl_vtx_topo_center{0.0};
        // doc sbnd_xin/docs/pr/51 round 3: the traditional (non-DL)
        // determine_overall_main_vertex takes main_cluster and its vertex
        // map BY VALUE, yet internally decides cluster swaps via
        // examine_main_vertices / check_switch_main_cluster[_2], which call
        // swap_main_cluster -- a function that already mutates persistent
        // state (the Flags::main_cluster bit on both clusters, and
        // other_clusters by reference) even though the by-value pointer
        // return is discarded by the caller.  A firing swap today therefore
        // leaves the job in a half-applied state: flags/other_clusters say
        // the new cluster, the caller's main_cluster variable still says the
        // old one.  false (default) = legacy = the caller passes its own
        // copies and discards them exactly as before, byte-identical; true =
        // apply the decision, propagating both the new main_cluster and the
        // pruned vertex map back to the caller (matching the DL sibling's
        // by-reference contract).  A "mvsa:" DEBUG sentinel fires whenever a
        // swap is decided, in both states, so the off-arms self-census how
        // often the traditional path swaps in production today.
        bool   m_main_vertex_swap_apply{false};
        // doc sbnd_xin/docs/pr/51 round 4: diagnostic-only TRACE probe for the
        // near-vertex short-cut investigation (owner Bee scan of
        // 131357/268067/285567/506746 -- rounds 2-3's main_vertex_graph_audit
        // fixed graph SHAPE; the owner's report is that the rough-path
        // Dijkstra (do_rough_path on "steiner_graph") or the
        // examine_structure_1 straight-line replacement produce a
        // charge-unsupported near-vertex chord -- a path-COST problem, not a
        // graph-shape one).  Prints per near-vertex segment: (P1) the origin
        // tag of its current wcpts() -- "rough" (do_rough_path via
        // create_segment_for_cluster), "straighten" (examine_structure_1),
        // "splice" (mvga op3 re-seat), or "other"; (P2) a live/dead/
        // unsupported interior support profile against the fitted dQ/dx;
        // (P3) a counterfactual Dijkstra re-route on a scratch gap-penalized
        // copy of the Steiner graph across a scale ladder.  No graph, no
        // fit, and no segment content is ever changed -- every line is
        // SPDLOG_LOGGER_TRACE, gated entirely by this knob.  false (default)
        // => none of the probe code runs => byte-identical.
        bool   m_rough_path_probe{false};
        // doc sbnd_xin/docs/pr/51 round 5: steiner gap penalty -- the H1 fix
        // the round-4 probe validated (see the m_steiner_gap_penalty block
        // in NeutrinoPatternBase.h).  Scale of the per-edge unsupported-
        // fraction penalty on the lazily-derived "steiner_graph_gap" flavor
        // that do_rough_path routes on when this is > 0.  0 (default) =>
        // the flavor is never built => byte-identical.  Sub-knobs ride the
        // C++ defaults (cm where dimensional).
        double m_steiner_gap_penalty{0};
        double m_sgp_dead_alpha{0.25};   // dead-sample weight in bad_fraction
        double m_sgp_min_edge{0.5};      // cm; shorter edges never scanned
        double m_sgp_sample_step{0.3};   // cm; edge-interior sampling step
        double m_sgp_point_radius{0.2};  // cm; test_good_point radius
        // doc sbnd_xin/docs/pr/51 round 6: weak-charge deficit term on the
        // same gap flavor (see NeutrinoPatternBase.h).  Active only when
        // both m_steiner_gap_penalty > 0 and this scale > 0; 0 (default)
        // => round-5 reweight path verbatim => byte-identical.
        double m_sgp_weak_scale{0};      // weak-charge penalty scale; 0 = off
        double m_sgp_weak_qref{2000};    // charge ref (calc_charge_wcp units, no cm conversion)
        // doc sbnd_xin/docs/pr/73: per-edge DEBUG sentinel for the sgp
        // scan.  Log-only; false = legacy (never emits).
        bool   m_sgp_edge_probe{false};
        // doc sbnd_xin/docs/pr/73 round 2, F3a: cap (cm) on how far the
        // penalized route may stray from the unpenalized one in do_rough_path;
        // over the cap, the base route is kept.  OFF-TEST IS `< 0`, not the
        // `<= 0` the scale knobs above use -- 0 is a meaningful cap here.
        double m_sgp_max_sep{-1};        // cm; < 0 = off (unbounded, legacy)
        // doc sbnd_xin/docs/pr/83: orient break_segment() splits to the wcpt
        // path (find_vertices) instead of boost source/target, which do not
        // track orientation -- a reversed edge otherwise yields crossed
        // children with vertex fits on the wrong track ends and stacked
        // duplicate trajectories (mcp1k 283040/59899/72586).  Prototype
        // parity (NeutrinoID_proto_vertex.h:595-601).  false = legacy.
        bool   m_break_seg_orient{false};
        // doc sbnd_xin/docs/pr/75: record, per event, the numbers the two
        // vertex selectors compared (compare_main_vertices scores; DL top-K
        // voxels + the seven rerank terms + the accept decision) so a hand
        // scan can become tuning input.  Pure recording; false = legacy =>
        // the board stays empty and the compiled config omits the key.
        bool   m_vertex_scoreboard{false};
        // doc sbnd_xin/docs/pr/79 §10: live-feature harvest for DL-vertex
        // training -- the exact SCN input cloud plus the traditional-path
        // per-candidate features the scorers compute and discard.  Pure
        // recording (no decision reads it); REQUIRES vertex_scoreboard
        // (warned and inert otherwise).  false = legacy => key omitted from
        // compiled config, calib JSON schema unchanged.
        bool   m_dl_vtx_harvest{false};
        double m_beam_window_low{0};   // beam window [low, high) on cluster_t0 (matched flash time, WCT units).
        double m_beam_window_high{0};  // low >= high (default) disables the gate: uBooNE single-main behavior.
        bool m_nu_skip_cosmic{false};  // if true (beam-gate only), skip in-window mains already tagged
                                       // cosmic upstream: flag_TGM, flag_STM, or lm_flag > 0.
        bool m_nu_skip_cosmic_bundle{false};  // if true, that verdict vetoes the whole flash bundle
                                              // (every in-window main sharing matched_flash_gid with a
                                              // cosmic-tagged main), not just the tagged main itself.
                                              // Inert unless m_nu_skip_cosmic.
        double m_nu_skip_cosmic_bundle_min_length{0};  // cm.  > 0 spares an UNTAGGED in-window main at
                                                       // least this long from the bundle veto (its
                                                       // cosmic-tagged bundle-mate stays excluded --
                                                       // companions are associated-only).  0 = legacy:
                                                       // the veto removes every bundle-mate.
                                                       // Inert unless m_nu_skip_cosmic_bundle.
                                                       // docs/pr/16 design A (SBND 15).
        bool m_skip_cosmic_companions{false};  // doc pr/20 Part I P4.  If true, a COMPANION
                                               // (other_clusters member) that is TGM- or
                                               // STM-tagged and at least
                                               // m_cosmic_companion_min_length long is dropped
                                               // from the neutrino's companion list.  Nothing
                                               // tags a companion today unless the taggers run
                                               // with evaluate_demoted_mains (P3), so this is
                                               // inert without it.
        bool m_nu_fallback_demoted_mains{false};  // sbnd_xin/docs/73 sec 12 (round 3).  If true and
                                                  // the primary loop selected NO candidate, a second
                                                  // pass considers DEMOTED mains (Flags::demoted_main,
                                                  // set by ClusteringUnmergeBundle's
                                                  // restore_demoted_mains and scored by the taggers
                                                  // under evaluate_demoted_mains) with the same
                                                  // window / cosmic / bundle-veto gates.  Never runs
                                                  // when a main-cluster candidate exists, so such
                                                  // events are byte-identical.  false = legacy: a
                                                  // demoted main is never a candidate.
        bool m_sp_photon_flag{false};  // doc pr/26 sec. 8.2 port gap.  If true, the single-photon
                                       // tagger's verdict is stored in TaggerInfo::photon_flag,
                                       // as prototype NeutrinoID.cxx:271 does
                                       // (`if (flag_sp) tagger_info.photon_flag = true;`).  The
                                       // port already runs singlephoton_tagger() and fills its
                                       // ~90 shw_sp_* features; only the return value was
                                       // discarded.  C++ default false = legacy: photon_flag
                                       // stays at its init_tagger_info() 0, so the uBooNE tagger
                                       // ntuple branch is byte-identical.
        double m_cosmic_companion_min_length{0};  // cm.  A tagged companion SHORTER than this
                                                  // stays in regardless of verdict, so a
                                                  // mis-tagged short neutrino daughter can never
                                                  // be silently dropped.  0 = no floor: drop on
                                                  // the verdict alone.  Own tuning, deliberately
                                                  // NOT inheriting nu_skip_cosmic_bundle_min_length
                                                  // -- a different question.
        // ---- doc sbnd_xin/docs/pr/36 §10 tagger-stage knobs -----------------
        // F1 (= P1): fiducial volume for the match_isFC (cluster_fc_check)
        // recompute.  Unset by default => the historical FiducialUtils /
        // sensitive-volume-union fallback, i.e. an absent "fiducial" key is
        // byte-identical to the pre-knob code.  Same knob, semantics and
        // rationale as TaggerCheckSTM's (TaggerCheckSTM.cxx:69-116): the
        // fallback volume is NOT the one TaggerCheckTGM/FC/STM use, and
        // match_isFC is numu XGBoost input 70 -- a binary feature, so a
        // containment disagreement moves it full-range.
        bool m_use_fiducial{false};
        std::vector<double> m_fv_tolerance;
        // F3 (= P2): gate for the single-photon SCE correction.  A SEPARATE
        // bool from clus_geom_helper on purpose: the key is shared with
        // fill_kine_tree, and reusing it alone would silently couple the two
        // consumers the day someone configures the helper for kine (doc
        // pr/36 §10.4).  Owner decision 2026-08-04: ships OFF on SBND (no
        // SBND SCE helper exists), so this is plumbing, not behavior.
        bool m_sp_sce_correction{false};
        // F4/F5/F6/F7: threaded to PatternAlgorithms; see
        // NeutrinoPatternBase.h for the full rationale of each.
        bool m_tagger_ordered_segment_sets{false};
        bool m_stem_endpoint_wcpt_parity{false};
        bool m_broken_muon_cluster_id_count{false};
        bool m_neutrino_type_bitmask{false};
        // doc pr/33 §10 EM-shower-clustering knobs (F1a/F1b, F2a/F2b/F2c,
        // F3, F4, F5): threaded to PatternAlgorithms; see
        // NeutrinoPatternBase.h for the full rationale of each.
        bool m_daughter_count_proto_main_vertex{false};
        bool m_daughter_count_proto_examine_showers{false};
        bool m_shower_pdg_from_start_segment{false};
        bool m_shower_pdg_from_shower_type{false};
        bool m_shower_pdg_exact_muon_test{false};
        bool m_pi0_id_shared_allocator{false};
        bool m_shower_flag_pdg_electron{false};
        bool m_shower_less_id_tiebreak{false};
        // doc pr/39: prototype-parity exclusion of a shower's own start
        // vertex from the end_point farthest-vertex search.  See
        // PatternAlgorithms::m_shower_endpoint_exclude_start_vertex.
        bool m_shower_endpoint_exclude_start_vertex{false};
        // doc pr/91 round 1 F1: also skip a node no member segment touches.
        // See PatternAlgorithms::m_shower_endpoint_skip_orphan_vtx.
        bool m_shower_endpoint_skip_orphan_vtx{false};
        // doc pr/91 round 3: flood-fill frontier test = visited, not merely
        // present.  See PatternAlgorithms::m_shower_walk_visited_parity.
        bool m_shower_walk_visited_parity{false};
        // doc sbnd_xin/docs/pr/40 -- track (proton/pion/muon) mis-identified
        // as electron.  F1 restores prototype-faithful PID persistence
        // (threaded via track_pid_options()); F2/F3 guard the wholesale
        // track-to-electron conversion sites with the segment's own dQ/dx.
        // See NeutrinoPatternBase.h for the full rationale of each.
        bool m_track_pid_persist_dqdx{false};
        bool m_shower_reclass_dqdx_guard{false};
        bool m_shower_topo_dqdx_guard{false};
        // doc sbnd_xin/docs/pr/40 round 2 -- two follow-on defects measured on the
        // Bee display of pr/40's fix: F4 replaces the rest-mass-only 4-mom
        // stub with an unconditional segment_cal_4mom call (zero KE ->
        // correct KE); F5 relabels an electron-labelled segment to pion when
        // it emanates from the neutrino vertex and its far end is a
        // charge-confirmed PID'd proton (a daughter no electron can
        // produce); F6 fixes a related negative-KE construction in
        // reclass_pinfo.  See NeutrinoPatternBase.h for the full rationale.
        bool m_track_pid_persist_4mom{false};
        bool m_shower_proton_daughter_pion{false};
        // doc sbnd_xin/docs/pr/40 round 4 -- two follow-on defects measured on the
        // Bee display of pr/40 round 2/3's F5: F7 clears the shower flags when
        // F5 relabels a segment to pion, so it stops being wrapped as a Shower
        // (its proton daughter had been swallowed into the shower's segment
        // set); F8 relabels a muon segment to pion when its far (non-main-
        // vertex) end is a multi-proton (>=2, charge-confirmed) hadronic
        // vertex. See NeutrinoPatternBase.h for the full rationale.
        bool m_shower_proton_daughter_pion_dissolve{false};
        bool m_muon_multi_proton_pion{false};
        bool m_track_pid_persist_dqdx_electron_guard{false};        // doc pr/40 round 5 F9
        bool m_shower_connect_main_vertex_straight_guard{false};    // doc pr/40 round 5 F10
        bool m_shower_traj_straight_guard{false};                   // doc pr/40 round 5 F11
        bool m_shower_absorb_track_guard{false};                    // doc pr/40 round 6 F12
        bool m_shower_absorb_unreachable_main{false};               // doc pr/65 round 3
        bool m_shower_connect_protected_pion_guard{false};          // doc pr/40 round 6 F13
        bool m_michel_stem_muon_rescue{false};                      // doc pr/40 round 6 F14
        bool m_shower_in_cascade_guard{false};                      // doc pr/74 round 2 P1
        double m_shower_in_max_len{40};                             // cm; pr/74 P1 tunable
        double m_shower_in_mip_hi{1.3};                             // ratio; pr/74 P1 tunable
        // doc pr/40 round 9 -- the rounds-7+8 straight-track PID guard
        // family + the B2 cross-cluster bridge.  Rationale comments in
        // NeutrinoPatternBase.h.
        bool m_shower_connect_from_vertices_straight_guard{false};  // doc pr/40 round 9 (round 8 Part A)
        bool m_shower_connect_start_seg_straight_guard{false};      // doc pr/40 round 9 (round 7 c2c, D1 re-target)
        bool m_examine_direction_dirsign_shower_in_guard{false};    // doc pr/40 round 9 (round 7 c2a, D2 re-scope)
        bool m_daughter_shower_angle_reclass_straight_guard{false}; // doc pr/40 round 9 (round 7 c2b)
        bool m_shower_topo_reexam_straight_guard{false};            // doc pr/40 round 9 (round 7 c1 safety net)
        double m_sfv_kink_max{25.0};                                // degrees; continuation-arm tunable
        bool m_shower_nv_bridge_track{false};                       // doc pr/40 round 9 B2
        double m_shower_nv_bridge_max_gap{1.8};                     // cm; B2 gap cut (steiner-cloud closest approach)
        // doc pr/92 -- stray-satellite drop from kine/PF.  Rationale
        // comments in NeutrinoPatternBase.h (pr/92 block).
        bool m_kine_drop_stray_satellites{false};                   // doc pr/92 master
        double m_kine_sat_min_energy{20.0};                         // MeV; drop-candidate floor
        double m_kine_sat_prox_max{8.0};                            // cm; main-cluster proximity exemption
        double m_kine_sat_angle_bad{60.0};                          // degrees; Arm A attachment-angle cut
        double m_kine_sat_angle_main{45.0};                         // degrees; Arm B main-vertex-angle cut
        double m_kine_sat_far_dis{90.0};                            // cm; Arm B far-attachment trigger
        double m_kine_sat_axis_dis_cut{30.0};                       // cm; shower axis integration radius
        double m_kine_sat_cont_kink{25.0};                          // degrees; Arm C continuation kink
        bool m_michel_stem_michel_check{false};                     // doc pr/74 round 2 P2
        double m_michel_stem_max_far_len{40};                       // cm; pr/74 P2 tunable
        bool m_shower_stem_backfill{false};                         // doc pr/74 round 2 K4
        double m_stem_backfill_max_len{30};                         // cm; pr/74 K4 tunable
        double m_stem_backfill_mip_lo{0.75};                        // ratio; pr/74 K4 tunable
        double m_stem_backfill_mip_hi{3.5};                         // ratio; pr/74 K4 tunable
        double m_stem_backfill_min_shower_len{40};                  // cm; pr/74 K4 tunable
        bool m_shower_conn3_unreachable{false};                     // doc pr/74 round 2 K5 (pr/65 rung 2)
        double m_conn3_unreachable_min_len{10};                     // cm; pr/74 K5 tunable
        double m_conn3_stitch_max{0};                               // cm; doc pr/84 r2 F3; 0 = off = byte-identical
        bool m_shower_dedup_start_seg{false};                       // doc pr/84 r3 S1; false = off = byte-identical
        bool m_shower_traj_michel_stem{false};                      // doc pr/74 round 4 K6 (18255-506746 muon+Michel)
        double m_michel_stem_traj_min_len{15};                      // cm; pr/74 K6 tunable
        double m_michel_stem_traj_max_len{45};                      // cm; pr/74 K6 tunable
        double m_michel_stem_traj_mip_lo{1.3};                      // x MIP median; pr/74 K6 tunable
        double m_michel_stem_traj_max_far_len{40};                  // cm; pr/74 K6 tunable (own ceiling, not P2's)
        double m_michel_stem_traj_min_kink_deg{40.0};               // deg; pr/74 K6 tunable
        bool m_shower_long_muon_keep_type{false};                   // doc pr/44
        bool m_single_muon_proton_chain_veto{false};                // doc pr/43 round 2 K1
        bool m_single_muon_long_muon_claim{false};                  // doc pr/43 round 2 K2
        bool m_pid_flag_reconcile{false};                           // doc pr/43 round 2 K3
        bool m_long_muon_stub_bridge{false};                        // doc pr/46
        mutable std::shared_ptr<TrackFitting> m_track_fitter;

        void load_trackfitting_config(const std::string& config_file);

};