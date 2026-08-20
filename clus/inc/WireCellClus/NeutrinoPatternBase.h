#include "WireCellClus/PRGraph.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/PRShower.h"
#include "WireCellClus/NeutrinoTaggerInfo.h"
#include "WireCellClus/PRVertexScoreboard.h"
#include "WireCellClus/IClusGeomHelper.h"
#include "WireCellClus/PRSegmentFunctions.h"

#include <array>
#include <cmath>

namespace WireCell::Clus::PR {

    /** BDT input features for the best pi0 candidate found by id_pi0_with_vertex (flag==1)
     *  or id_pi0_without_vertex (flag==2).  All angles in radians, energies and distances
     *  in WireCell internal units (multiply by 1/units::MeV, 1/units::cm for output).
     *  Default-initialised to zero (flag==0 means no pi0 found).
     */
    /// Type aliases for the pre-collected 2D charge maps used by cal_kine_charge.
    /// Collect once per shower_clustering_with_nv call via TrackFitting::collect_2D_charge()
    /// and reuse across all showers/merges to avoid O(N_hits) re-collection per shower.
    using ChargeMap = std::map<TrackFitting::CoordReadout, TrackFitting::ChargeMeasurement>;
    using WireMap   = std::map<std::pair<int,int>, std::vector<std::tuple<int,int,int>>>;

    /// Calibration constants of the charge -> kinetic-energy conversion in
    /// NeutrinoEnergyReco.cxx (cal_kine_charge and calculate_shower_kinematics).
    /// Every default below is the uBooNE-tuned number that used to be a literal
    /// at the three duplicated fudge/recom blocks and inside
    /// kine_charge_from_maps -- absent config keys are therefore byte-identical.
    /// None of them has been re-derived for any other detector; they are exposed
    /// so that a calibration can move them without a code change (see
    /// wcp-porting-img sbnd_xin/docs/pr/2 sec. 2e(iii)).
    ///
    /// The conversion is
    ///     E = sum_plane(w_p * Q_p) / sum(w) / recom / fudge * w_value[eV] * 1e-6 MeV
    /// with the (med,max) plane-asymmetry switch described at plane_asym_switch.
    struct KineChargeOptions {
        /// Average recombination survival fraction and residual scale factor for
        /// a track-like object.  uBooNE values; field- and calibration-dependent.
        double fudge_factor{0.95};
        double recom_factor{0.7};
        /// Same pair for an object flagged shower-like (higher local dE/dx =>
        /// more recombination, hence the smaller survival fraction).
        double shower_fudge_factor{0.8};
        double shower_recom_factor{0.5};
        /// Recombination survival for |pdg| == 2212 (proton).  The fudge factor
        /// deliberately stays at fudge_factor for protons -- only recombination
        /// is specialised, as in the prototype.
        double proton_recom_factor{0.35};
        /// Per-plane weights {U, V, W} of the charge average.  The uBooNE
        /// 0.25/0.25/1.0 encodes "induction planes are a quarter as trustworthy
        /// as collection".  A zero weight is legitimate ("ignore this plane");
        /// only an all-zero triple is rejected.
        std::array<double, 3> plane_weights{{0.25, 0.25, 1.0}};
        /// If the two largest per-plane charges disagree by more than this
        /// relative asymmetry, the three-plane weighted average is replaced by
        /// the (median, minimum) pair -- i.e. the largest plane is dropped as
        /// likely contaminated.  Dimensionless, uBooNE-tuned.
        double plane_asym_switch{0.04};
        /// Argon W-value in eV per electron-ion pair.  Duplicated here rather
        /// than read from the recombination model's Wi (docs/pr/2 sec. 2e(iii),
        /// "low risk, note only") -- if the model's Wi is ever wired in, this
        /// knob is what it replaces.
        double w_value{23.6};
        /// doc pr/35 §10.2 (F1 = P1 + P8).  false = today's cached
        /// Shower::get_particle_type() read at the four fill_kine_tree sites.
        /// true = the prototype's live start-segment read (kine.h:53 :67 :175
        /// :187), correct independently of the cache's refresh schedule.
        /// Config key kine_shower_pdg_live; absent => legacy, byte-identical.
        bool shower_pdg_live{false};
        /// doc pr/99 round 3 (C1).  false = legacy: every shower's
        /// cal_kine_charge independently rescans the whole event's 2D charge
        /// maps, so two spatially interleaved showers each collect the SAME
        /// cells' full charge and the event Enu double-counts (SBND
        /// 18255-168596: one 2016 MeV cascade split into 1843+1265 MeV,
        /// Enu 2619->3533).  true = after all shower-structure passes,
        /// recompute every shower's kine_charge in ONE scan of the charge
        /// maps where each 2D cell is credited to exactly one shower (the
        /// nearest accepting cloud; tie -> lowest shower creation id).  The
        /// prototype is also ownership-free (NeutrinoID_energy_reco.h) --
        /// this is a deliberate divergence in both trees, not a port fix.
        /// Config key kine_charge_dedup; absent => legacy, byte-identical.
        bool dedup{false};
        /// doc pr/99 round 3 (C1b).  Prototype parity: WCP rebuilds a
        /// shower's point clouds from CURRENT members before every
        /// cal_kine_charge read (NeutrinoID_energy_reco.h:99); the toolkit
        /// clouds are add-only merges, so stale points from departed members
        /// keep pulling charge in.  true = the final recompute pass reads
        /// EPHEMERAL member-true rebuilt clouds (stored clouds untouched --
        /// taggers and the pi0 start_point derivation query them later, and
        /// a rebuild changes row order hence kNN tie-breaks).  Config key
        /// kine_charge_rebuild; absent => legacy, byte-identical.
        bool rebuild{false};
    };

    struct Pi0KineFeatures {
        int    flag{0};      ///< 0=none, 1=with_vertex, 2=without_vertex
        double mass{0};      ///< reconstructed pi0 invariant mass
        double vtx_dis{0};   ///< distance from pi0 decay vertex to main vertex
        double energy_1{0};  ///< kine_charge of first shower
        double theta_1{0};   ///< polar angle of first shower direction (from z-axis)
        double phi_1{0};     ///< azimuthal angle of first shower direction
        double dis_1{0};     ///< distance from pi0 vertex to first shower start point
        double energy_2{0};  ///< kine_charge of second shower
        double theta_2{0};   ///< polar angle of second shower direction
        double phi_2{0};     ///< azimuthal angle of second shower direction
        double dis_2{0};     ///< distance from pi0 vertex to second shower start point
        double angle{0};     ///< opening angle between the two shower directions
    };

    class PatternAlgorithms{
        public:
        bool m_perf{false};  // if true, print per-step timing to stdout

        // Direction-weakness read fidelity switch (default false = legacy port
        // behavior, byte-identical).  The prototype's ONLY public accessor is the
        // score-thresholded ProtoSegment::is_dir_weak() (the raw member is
        // private); the port read the raw flag at every one of the ~83 sites
        // instead.  When true, all reads route through segment_is_dir_weak()
        // (PRSegmentFunctions.cxx), restoring the prototype semantics: a mu/p
        // segment whose particle_score is poor (or still the sentinel 100) is
        // treated as weak.  See wcp-porting-img sbnd_xin/docs/pr/6.
        bool m_dir_weak_use_score{false};
        // All PatternAlgorithms code must read direction weakness through this,
        // never seg->dir_weak() directly (flag propagation in break_segment is
        // the one intentional raw read, and it lives in PRSegmentFunctions.cxx).
        bool seg_dir_weak(SegmentPtr seg) const;

        // MIP dQ/dx scales, INTERNAL units (charge per WCT length unit).
        // Two distinct uBooNE-legacy roles, previously hard-coded everywhere
        // (see wcp-porting-img sbnd_xin/docs/pr/7 sec. 5 and pr/8):
        //  - m_mip_dqdx        (legacy 50e3 e/cm): the flat-template amplitude in
        //    the do_track_comp KS comparison and the segment_cal_4mom /
        //    segment_is_shower_trajectory scale.
        //  - m_mip_dqdx_median (legacy 43e3 e/cm): the scale every median-dQ/dx
        //    ratio threshold (x1.2/1.3/1.4/1.75...) and BDT normalization is
        //    quoted against.
        // Defaults = the uBooNE numbers => absent config keys are byte-identical.
        // SBND production sets 56000 / 48000 e/cm (owner 2026-07-30; 48000 is a
        // placeholder scaled by the uBooNE 43/50 ratio, pending an SBND MIP
        // median measurement).
        double m_mip_dqdx{50000/units::cm};
        double m_mip_dqdx_median{43000/units::cm};

        // Long shower-topology demote length (doc sbnd_xin/docs/pr/25 sec 3).
        // Passed to every segment_is_shower_topology call site: a segment the
        // topology test would flag kShowerTopology whose geometric length
        // exceeds this is demoted to a track instead, so it receives real
        // track PID rather than the hard-coded pdg=11/score=100.
        // Motivation: the test's only measurement axis satisfies
        // dir_3.xhat = sin(angle-to-drift), so for a near-isochronous segment
        // it measures spread along the DRIFT axis, where points sit on a
        // 0.313 cm time-slice lattice -- against a 0.4 cm "large spread" cut.
        // 86 of 91 long firings across a 572-event manifest carry no evidence
        // above that noise floor, and an owner hand-scan of all 10 such
        // segments on a selected nu-candidate main cluster found 10/10 tracks.
        // C++ default 0 => the guard never fires => byte-identical.
        // 50*units::cm reproduces the scan-supported rule for 9 of the 10
        // scanned events; ~45*units::cm covers all 10 (evt 400504 measures
        // 49.1 cm by segment_track_length, the length this guard uses).
        double m_shower_topo_demote_len{0};

        // Isochronous first-segment endpoint finding (doc sbnd_xin/docs/pr/24
        // round 2, SBND evt 271851).  When ON and the cluster is sheet-like
        // (long, with a small quantile-trimmed drift-x extent), the first PR
        // segment's endpoints are picked as the quantile-trimmed extremes
        // along the sheet's principal axis (each snapped to a Steiner
        // terminal) instead of the wire-footprint boundary metric, and the
        // local-PCA endpoint refinement (degenerate on a filled 2-D sheet) is
        // skipped on that branch.  No prototype counterpart; the prototype's
        // boundary metric (PR3DCluster_path.h:530-536, wire-delta coefficients
        // *0.0) is preserved byte-identically when off.  Precedents:
        // adjust_wcpoints_parallel (data/src/PR3DCluster.cxx:428, separation
        // only) and search_for_connection_isochronous
        // (pid/src/PR3DCluster_graph.h:1445, call site disabled upstream).
        // C++ default false => legacy path byte-identical.
        bool   m_iso_endpoint{false};
        double m_iso_endpoint_min_length{40 * units::cm};
        double m_iso_endpoint_max_xext{25 * units::cm};
        double m_iso_endpoint_xext_frac{0.35};
        double m_iso_endpoint_xext_quantile{0.02};
        // Round 3 (doc pr/24 sec 15).  Round 2 picked each endpoint as the band
        // point nearest the centroid of a 3 cm band at the QUANTILE-TRIMMED
        // axial extreme -- inward-biased by construction (8-15 cm on a
        // 7000-point track), which left a tip stub for find_other_segments to
        // claim as its own segment and put a 0.9 deg "vertex" in a straight
        // track (SBND mcp1k evt 284794, 59899).  Instead: the endpoint is the
        // TRUE, untrimmed axial extreme over all qualified points (so it can
        // never move inward -- "leaves no stub" holds by construction), and
        // lateral centring is applied only WITHIN the 3 cm end band, which is
        // what still keeps the pick off a sheet's edge corner.
        // m_iso_endpoint_min_aspect additionally requires real 2-D sheet-ness
        // (trimmed transverse/axial extent ratio), so the branch is inert on
        // 1-D track-like clusters rather than merely harmless.
        // m_iso_endpoint_tube_radius is DIAGNOSTIC ONLY: the DEBUG line reports
        // whether the chosen endpoint lies within it of the axis line.  It was
        // a hard filter in a first draft; measured on the 38 gated clusters
        // that filter pulled endpoints up to 28.6 cm inward on long/curved
        // objects (doc pr/24 sec 15), the same failure it was meant to cure.
        double m_iso_endpoint_tube_radius{4 * units::cm};
        double m_iso_endpoint_min_aspect{0.12};

        // Round 5 (doc sbnd_xin/docs/pr/24 sec. 18, SBND 18259-42280 /
        // 18255-271851 / 18255-350186).  examine_vertices_3
        // (NeutrinoStructureExaminer.cxx) is the stage that revisits the main
        // cluster's two ORIGINAL init_first_segment endpoints and tries to
        // push each one further out via get_local_extension -- a 10 cm-radius
        // Hough-transform direction estimate (faithful port,
        // PR3DCluster_path.h:288-316).  Neither the toolkit nor the prototype
        // (NeutrinoID_proto_vertex.h:2412-2463) ever checks that the
        // "extension" actually moves the vertex FARTHER from the segment's
        // far endpoint -- only that it differs from both existing points and
        // that the rebuilt path isn't more than 2x longer. At the axial
        // extreme of an isochronous (drift-perpendicular) sheet -- exactly
        // where m_iso_endpoint picks its seed -- the local 10 cm neighbourhood
        // is dominated by the sheet's transverse spread, so the Hough
        // direction estimate is poorly conditioned and can point BACK into
        // the cluster.  Measured on all three events: the "extension" landed
        // 7.5-8.9 cm closer to the far endpoint than the original vertex,
        // amputating the delivered trajectory by exactly that much (round 4's
        // 8.4-10.9 cm undershoot). Legacy (non-iso) endpoints sit at a
        // different, better-conditioned point and lose only ~1 cm to the same
        // bug, which is why the owner's report is iso-only: the bug is
        // general, but iso_endpoint's picks are where it is worst.
        // This is a prototype limitation (M15), not a port error: fixed here
        // behind a knob rather than corrected unconditionally.
        // m_v3_extension_guard, when true, rejects get_local_extension's
        // result unless it increases the vertex's distance to the segment's
        // OTHER (far) endpoint by more than -m_v3_extension_min_gain (a small
        // negative tolerance permits the legacy arm's ~1 cm legitimate
        // retreat while rejecting the multi-cm amputation).
        // C++ default false => examine_vertices_3 keeps today's unconditional
        // accept => byte-identical.
        bool   m_v3_extension_guard{false};
        double m_v3_extension_min_gain{-1.0 * units::cm};

        // doc sbnd_xin/docs/pr/67 -- LOG-ONLY diagnostic probe for "the fitted
        // trajectory does not cover the image", worst in isochronous
        // topologies (owner's four cases: 18264-137238, 18259-42280,
        // 18345-21073, 18255-58717).  Emits, per main-cluster PR pass:
        //   P1  which gate of find_iso_first_segment_endpoints rejected a
        //       cluster.  Today only the ASPECT rejection logs, so a cluster
        //       thrown out by min_length / max_xext / xext_frac / degenerate
        //       axis is indistinguishable from one that was never tried --
        //       and 18255-58717 is exactly such a silent rejection.
        //   P2  get_local_extension's drift-angle early return.  That function
        //       returns the vertex UNCHANGED when the local Hough direction is
        //       within 7.5 deg of perpendicular-to-drift, i.e. precisely in
        //       the isochronous case, so the one stage meant to push an
        //       endpoint outward is a structural no-op there.
        //   P4  per-round find_other_segments census (segment count in/out),
        //       so "the 2-round branch-search budget was exhausted with work
        //       still to do" becomes visible rather than inferred.
        // C++ default false => no lines, no behavior change, byte-identical.
        bool   m_traj_cover_probe{false};

        // doc sbnd_xin/docs/pr/72 round 2 -- examine_structure_3
        // (NeutrinoStructureExaminer.cxx) merges any degree-2 junction
        // whose 10cm/3cm direction-agreement angles both clear lenient
        // thresholds, with no check for whether it is a genuine near-vertex
        // track stub meeting a shower trunk rather than one particle's
        // trajectory the tracker happened to split (18255-196649: a real
        // 6.28cm stub, terminal at the true neutrino vertex, silently
        // absorbed into a 33cm shower trunk).  m_es3_stub_guard true
        // suppresses the merge when es3_stub_suppress (PRSegmentFunctions.h)
        // says the junction looks like that case; see that function's doc
        // comment for the predicate and the 117-event census that fit the
        // five sub-parameters below. C++ default false => the guard in
        // examine_structure_3 is never evaluated => byte-identical to
        // before this round.
        bool   m_es3_stub_guard{false};
        double m_es3sg_stub_max{7 * units::cm};   // short-arm length ceiling
        double m_es3sg_len_ratio{2.0};            // long/short length ratio floor
        double m_es3sg_ang3_min{15.0};            // degrees; local-kink floor
        double m_es3sg_ang_ratio{1.0};            // require ang3 > ratio * ang10
        bool   m_es3sg_require_terminal{true};    // short arm's far end must be degree 1

        // doc pr/67 -- counterfactual for the owner's own 137238 hypothesis
        // ("is that limited by not sufficient rounds of doing the branch
        // searching?").  find_proto_vertex's nrounds_find_other_tracks is
        // HARDCODED at its three TaggerCheckNeutrino call sites (2 for the
        // main cluster and associated clusters, 1 for the third pass) with no
        // config surface at all.  When > 0 this overrides the main-cluster
        // count only, so the hypothesis can be measured instead of argued.
        // C++ default 0 => the hardcoded value stands => byte-identical.
        // DIAGNOSTIC: raising it changes reconstruction output by design.
        int    m_pr_find_other_rounds{0};

        // Cathode kink veto (doc sbnd_xin/docs/pr/20 Part II, B0).  Passed to
        // segment_search_kink from break_segments: a candidate fit point within
        // m_cathode_kink_xcut of the cathode plane at m_cathode_x is skipped, so
        // the ~2 cm transverse cathode mismatch cannot invent a vertex that breaks
        // a crossing cosmic into two particles.  C++ default 0 => no point is ever
        // skipped => byte-identical to the pre-pr/20 behavior.
        double m_cathode_x{0};
        double m_cathode_kink_xcut{0};

        // Wide-baseline cathode kink accept (doc sbnd_xin/docs/pr/47 sec 8,
        // O1).  Passed to segment_search_kink from break_segments: at a
        // cathode-crossing fit index the shipped index-windowed refl_angle
        // statistic is suppressed by the gap/distortion (52085: wide-baseline
        // 33-38 deg kink measures only ~23 deg locally), so a fifth accept
        // path fires when the skirt-excluded PCA turn angle across the
        // crossing is >= m_cathode_wide_kink_angle degrees.  Skirt/baseline
        // are internal-unit lengths (config takes cm).  C++ default angle 0
        // => path never evaluated => byte-identical legacy kink search.
        double m_cathode_wide_kink_angle{0};                 // deg; 0 = OFF
        double m_cathode_wide_kink_skirt{3*units::cm};
        double m_cathode_wide_kink_baseline{15*units::cm};

        // ---- doc sbnd_xin/docs/pr/83: oriented break_segment splits ---------
        //
        // m_break_seg_orient -- passed as break_segment()'s orient_split at
        // every PatternAlgorithms call site (break_two_end_dqdx,
        // shower start-segment break, snap_main_vertex_to_kink).
        // break_segment slices the parent's wcpts/fits into a front half and
        // a back half but hands them to (boost source, boost target), which
        // do NOT track path orientation; a reversed edge -- legal, and
        // routinely produced by the examine_vertices re-routes -- yields
        // crossed children: each carries the half that terminates at the
        // OTHER vertex, the vertex fit points land on the wrong track ends
        // (67-118 cm off on mcp1k 283040/59899/72586), the fitted
        // trajectories of both children stack on one arm, and the
        // fit-vs-wcpt divergence seeds runaway vertex-activity bridges.
        // The prototype cannot cross: it resolves (start_v, end_v) by
        // wcpt-index equality tested in both orientations
        // (NeutrinoID_proto_vertex.h:595-601), so true = prototype parity
        // via find_vertices().  C++ default false => byte-identical legacy.
        bool   m_break_seg_orient{false};

        // ---- doc sbnd_xin/docs/pr/48: back-to-back track fixes --------------
        //
        // m_two_end_break -- the two-end residual-range break pass
        // (break_two_end_dqdx, run inside find_proto_vertex after
        // examine_vertices): a main cluster whose topology is a single
        // non-stub segment, both endpoints in the fiducial volume, with dQ/dx
        // rising at BOTH ends (two Bragg ends = a junction at the interior
        // dip: 18255-51513/56211/57903/57485) gets broken at the argmin of
        // the joint two-arm stopping-template score
        // (segment_two_end_break_scan) and the new vertex protected
        // (VertexFlags::kProtectedBreak) from the examiner merge passes.
        // C++ default false => pass never runs => byte-identical.
        // The teb_* operating point mirrors TwoEndBreakOptions; lengths are
        // internal units here (config takes cm), angles degrees.
        bool   m_two_end_break{false};
        double m_teb_min_len{10*units::cm};
        double m_teb_min_arm{1.8*units::cm};
        int    m_teb_min_arm_pts{4};
        double m_teb_stub_max{4*units::cm};
        double m_teb_accept_range{15*units::cm};
        double m_teb_rise_r1{1.3};
        double m_teb_rise_r2{1.15};
        double m_teb_abs_end_min{1.7};
        double m_teb_dip_floor{0.6};
        double m_teb_score_cap_r1{0.6};
        double m_teb_score_cap_r2{0.9};
        double m_teb_turn_angle{25.0};                       // deg; <= 0 disables route R2
        double m_teb_turn_baseline{35*units::cm};
        double m_teb_turn_skirt{3*units::cm};
        // doc sbnd_xin/docs/pr/90 round 2 -- two default-OFF refinements:
        // m_teb_turn_min_arm_frac: route R2's turn argmax PREFERS an index
        // whose PCA arms can each span this fraction of teb_turn_baseline,
        // when one clears teb_turn_angle on its own; legacy argmax is the
        // fallback (mirrors TwoEndBreakOptions::turn_min_arm_frac); 0 =
        // legacy argmax, byte-identical.
        // m_teb_second_max: tolerate additional long (> teb_stub_max)
        // segments in the entry gate as long as exactly one segment exceeds
        // this cap (that one is the candidate); 0 = legacy strict
        // single-long-segment gate, byte-identical.  Internal length units
        // (config takes cm).
        double m_teb_turn_min_arm_frac{0.0};
        double m_teb_second_max{0};
        // doc sbnd_xin/docs/pr/90 round 4 (sec 9.5) -- three default-OFF
        // knobs for the round-3 residual classes:
        // m_teb_chain_topology (D1): when n_long > 1, admit iff the
        // cluster's segment graph is a simple path (the owner's "still a
        // line, no 3-track vertex") and the candidate is the unique longest
        // segment; chain-admitted candidates are scanned by route R3 only,
        // so admission also requires both R3 knobs.  false = legacy gate,
        // byte-identical.
        // m_teb_r3_turn (deg) / m_teb_r3_hot (x m_mip_dqdx_median) (D3):
        // route R3's local-turn threshold on the 10 cm-baseline turn and the
        // +-2 cm vertex-activity corroboration floor
        // (TwoEndBreakOptions::r3_turn / r3_hot).  <= 0 = R3 off.
        // m_teb_bragg_veto_turn (deg) (D4): veto an accepted R2 break whose
        // turn is below this when its short-arm end is not Bragg-consistent
        // (TwoEndBreakOptions::bragg_veto_turn; owner-calibrated keep/kill
        // rule, sec 9.4b).  <= 0 = off, byte-identical.
        bool   m_teb_chain_topology{false};
        double m_teb_r3_turn{0.0};
        double m_teb_r3_hot{0.0};
        double m_teb_bragg_veto_turn{0.0};

        // m_kink_walk_dqdx_stop -- 59335 fix (a): forwarded to
        // segment_search_kink so a dQ/dx-confident C4/straightness accept
        // stops the proto_extend_point walk AT the kink instead of
        // overshooting to the terminus (see the PRSegmentFunctions.h knob
        // comment).  C++ default false => byte-identical.
        bool   m_kink_walk_dqdx_stop{false};
        // Shared Bragg-hot ratio (x m_mip_dqdx_median) for BOTH 59335 fixes:
        // fix (a) stops the walk only at a kink whose 5-point local dQ/dx
        // exceeds it, fix (b) protects only a break sitting on charge above
        // it.  At the legacy 25/43 "not too low" scale nearly every C4
        // accept qualifies and the footprint is sample-wide; 1.7 (the same
        // scale as teb_abs_end_min) confines both fixes to genuine Bragg
        // stubs (59335 reads 2.5-6x).  Inert while both bools are false.
        double m_kink_dqdx_hot_ratio{1.7};

        // m_kink_break_protect -- 59335 fix (b): when break_segments breaks
        // at a kink accepted by C4 (dQ/dx-assisted) or A0 (wide-baseline
        // cathode accept), the new vertex gets VertexFlags::kProtectedBreak
        // so examine_vertices_4's unconditional < 2 cm stub-absorption floor
        // (NeutrinoStructureExaminer.cxx) cannot silently erase a
        // high-confidence break that produced a short arm (59335: a correct
        // C4 accept 0.28 cm from truth leaves a ~1 cm proton stub, absorbed
        // today).  C++ default false => flag never set => every examiner
        // check is a no-op => byte-identical.
        bool   m_kink_break_protect{false};

        // ---- doc sbnd_xin/docs/pr/50: main-vertex kink-consistency snap ----
        //
        // 172230-class near-vertex robustness: find_proto_vertex's recursive
        // break partition is globally sensitive to fit perturbations, so a
        // small distortion far from the vertex (there: 200 fit_blob_coverage
        // deweight firings ~90 cm away at a shower tip) can reshuffle the
        // partition and erase the proto-vertex at the TRUE image kink.
        // determine_main_vertex can only choose surviving graph vertices
        // (172230: nearest survivor 2.7 cm off), and improve_vertex/MyFCN is
        // a local optimizer (0.43 cm soft prior) that can never jump back --
        // the final trajectory rounds the corner ("jumps over the vertex")
        // with starved dQ/dx.  This pass runs ONCE, after the overall main
        // vertex is final and before the final improve_vertex: it scans the
        // WCPT paths (steiner/image-anchored -- TrackFitting never rewrites
        // wcpts, so they are immune to exactly the fit distortion that
        // caused the loss) of the segments incident to the main vertex for
        // a strong interior turn within m_vks_radius; when the image says
        // the track passes STRAIGHT THROUGH the vertex (an incident
        // segment's stub toward the turn is anti-parallel within
        // m_vks_collinear deg of another incident arm) and turns at the
        // interior point instead (turn >= m_vks_angle deg and >=
        // m_vks_margin deg stronger than the residual bend at the vertex),
        // it re-breaks the segment there (break_segment), protects the new
        // vertex (kProtectedBreak), re-seats the main vertex on it, and
        // refits once so improve_vertex polishes a corner-anchored
        // trajectory.  Guards: never fires on a kProtectedBreak vertex
        // (pr/48 TEB junctions), needs vertex degree >= 2, declines when
        // the collinear opposite arm is Bragg-hot at the vertex (median
        // dQ/dx within 3 cm > m_vks_hot_ratio x m_mip_dqdx_median -- a hot
        // arm ending at V is evidence V is a real junction), and when a
        // graph vertex already sits within m_vks_min_dis of the turn.
        // Single strongest winner, one snap per event, no recursion.
        // Toolkit-only (no prototype counterpart; nearest prototype
        // machinery is examine_structure_final_2/_3's 2.0/2.5 cm merges,
        // which cannot reach 2.7 cm).  C++ default false => the pass
        // returns immediately => byte-identical.
        bool   m_vertex_kink_snap{false};
        double m_vks_radius{5*units::cm};     ///< scan window outer edge (arclength from vertex)
        double m_vks_min_dis{0.5*units::cm};  ///< scan window inner edge + dedupe radius vs existing vertices
        double m_vks_angle{25};               ///< accept threshold on the interior turn, deg
        double m_vks_margin{10};              ///< interior turn must beat the vertex bend by this, deg
        double m_vks_collinear{30};           ///< pass-through test: stub vs other arm within this of 180 deg (172230 measures 23.5 at the true corner -- near-corner arms curve)
        double m_vks_skirt{0.3*units::cm};    ///< PCA window inner skirt (path_scan_vertex_kink)
        double m_vks_baseline{2*units::cm};   ///< PCA window length (path_scan_vertex_kink)
        double m_vks_min_arm{1.5*units::cm};  ///< outward arclength support required beyond the turn
        double m_vks_fit_miss{0.35*units::cm}; ///< snap only when the fit misses the image corner by at least this (a modeled kink has fit points on it)
        double m_vks_hot_ratio{0};            ///< OPTIONAL Bragg-hot veto scale, x m_mip_dqdx_median; default 0 = off (misfires on the failure class: misassigned charge reads hot)
        // doc sbnd_xin/docs/pr/85: carry the prongs through the snap stub.
        // The snap re-breaks the arm at the image corner K and re-seats the
        // main vertex there, but the OLD vertex keeps every other arm plus
        // the new residual edge to K -- that residual IS the "interposed
        // stub" of the pr/85 census (mode 1a-VIA: a >= 10 cm prong reaching
        // the neutrino vertex only through a sub-3 cm segment; SNAP arcs
        // 0.93 / 1.07 cm measured on evt59685 / evt280972 against interposed
        // stubs of 1.18 / 0.29 cm).  When best.arc (wcpt-space arclength old
        // vertex -> K, the same space break_segment cuts in) is below this
        // threshold, every arm at the old vertex is spliced through the
        // residual's wcpts onto the new main vertex (carry_prong_verify /
        // carry_prong_execute, PRSegmentFunctions.h: SegmentPtr identity
        // preserved, all-or-nothing -- any arm failing pre-verification
        // leaves the graph exactly as today), the residual is removed and
        // the old vertex, now degree 0, is dropped.  The carried arms ride
        // the snap's existing trailing do_multi_tracking refit.  0 (default)
        // => the block is unreachable => byte-identical inside the
        // production-ON snap pass.
        double m_vks_carry_prong{0};          ///< carry threshold on the snap arc, internal units; 0 = off

        // ---- doc sbnd_xin/docs/pr/51: main-vertex graph audit --------------
        //
        // Near-vertex PR *graph-shape* pathologies that survive every
        // point-level fit metric (the ribbons DO follow charge -- the graph
        // is wrong): (a) duplicated corridors -- two segments riding one
        // charge ribbon (360535's ~1 cm parallel pair splitting one MIP's
        // charge 934+1219; 268067's 15003 riding the proton at 86% overlap;
        // 285567's 8010/8032/8033 full-charge double-counting variant);
        // (b) charge-less bridges -- connectivity-only segments at a small
        // fraction of MIP median dQ/dx (268067's 15005 at ~0.1 MIP;
        // 285567's 8031 with 11/13 fit points at zero charge riding a
        // full-MIP track); (c) micro-stubs at the vertex inflating the
        // prong count, mostly carved by improve_vertex (131357's 1.56 cm
        // 9-point "proton"; 142421's 7081/7082; 285567's 81/82/83).
        // This pass runs ONCE on the main cluster, AFTER the final
        // improve_vertex (the stubs it must see are created there) and
        // BEFORE clustering_points/examine_direction, scoped to
        // m_mvga_radius around the main vertex.  Ordered operations, each
        // independently disabled by zeroing its knob, each edit printed as
        // one "mvga:" DEBUG sentinel with the measured quantities:
        //   op1 duplicate-corridor merge: pair overlap (shorter onto
        //       longer, path_overlap_fraction at m_mvga_dup_tol) >=
        //       m_mvga_dup_frac deletes the lower-integrated-charge member
        //       and reconnects its orphaned endpoint to the survivor's
        //       nearest endpoint by a direct rough-path edge;
        //   op2 charge-less-bridge removal: non-terminal segment with
        //       median dQ/dx < m_mvga_bridge_mip x m_mip_dqdx_median is
        //       deleted iff the far side stays connected to the main
        //       vertex or reconnects to a reachable vertex within
        //       m_mvga_reconnect;
        //   op3 micro-stub absorb + re-seat: terminal stub at an anchor
        //       vertex (the main vertex, or -- when m_mvga_satellite > 0 --
        //       any other main-cluster vertex within m_mvga_satellite of the
        //       main vertex, doc pr/51 round 3: 142421's 7082/7023 and
        //       285567's residual sit on such satellites, 1.2-1.5 cm out,
        //       outside the main-vertex-only round-2 scope) shorter than
        //       m_mvga_stub is absorbed when its corridor overlap with a
        //       sibling prong at the same anchor >= m_mvga_dup_frac OR it is
        //       point-degenerate (<= m_mvga_stub_pts valid fits); at the
        //       main vertex only, when the overlap gate passed and a sibling
        //       is the stub's collinear continuation (>= m_mvga_reseat_angle
        //       deg), the sibling is extended through the vertex and the
        //       main vertex re-seats at the stub far end (131357: 0.9 mm
        //       from the true image corner) -- satellite anchors are
        //       absorb-only, never re-seated;
        //   op4 one local do_multi_tracking refit (the examiner contract)
        //       so dQ/dx and the display layers reflect the audited graph.
        // Guards: kProtectedBreak vertices are never removed, re-seated, or
        // used as a satellite anchor (pr/48 + pr/50 snap precedent); the
        // main vertex is never removed; segments created by this pass's own
        // reconnects are exempt from every op, not just op1 (no
        // delete/recreate cycling); every op is edit-capped; no recursion.
        // Toolkit-only (no prototype counterpart; ancestry:
        // prototype NeutrinoID_final_structure.h examine_structure_final_1
        // lines 545-696 / _1p lines 402-544 and NeutrinoID_improve_vertex.h
        // eliminate_short_vertex_activities lines 365-1018).  C++ default
        // false => the pass returns immediately => byte-identical.
        // sbnd_xin/docs/73 sec 12 (round 3, SBND data evt 78242).
        // eliminate_short_vertex_activities case 5 deletes a non-pre-existing
        // segment when EVERY wcpt is within 0.45 cm of some pre-existing
        // segment in ALL THREE 2D views.  The per-view distance comes from
        // DynamicPointCloud::get_closest_2d_point_info, which returns the
        // SENTINEL -1 when the queried (plane,face,apa) 2D index is EMPTY --
        // i.e. when the pre-existing segment has no points in the query
        // point's APA.  -1 passes every "< 0.45 cm" test, so on a
        // CATHODE-CROSSING cluster (the only place a cluster's segments span
        // two APAs) one other-APA segment vacuously "covers" every point in
        // all views, n_good stays 0, and the junction segment the cathode
        // rescue exists to create is unconditionally deleted -- 78242's 71 cm
        // track_fit hole, and (downstream) its muon far half absorbed into an
        // EM shower via absorb_unreachable_main once the junction was gone.
        // The exemption for pre-existing segments does not save it because
        // examine_vertices_1's merge_two_segments_into_one creates a NEW
        // segment object.  When true, a -1 view distance is treated as "no
        // information" (1e9) instead of "covered".  Single-APA clusters have
        // no empty per-APA index, so uBooNE/prototype behavior (single TPC,
        // NeutrinoID_improve_vertex.h:365-1018) is unreachable by this knob.
        // C++ default false => byte-identical legacy.
        bool   m_esva_ignore_empty_2d{false};
        bool   m_main_vertex_graph_audit{false};
        double m_mvga_radius{15*units::cm};   ///< audit scope around the main vertex
        double m_mvga_dup_tol{1.4*units::cm}; ///< op1/op3 corridor-overlap point tolerance (must clear the fitter's ~1 cm ribbon separation: 360535's pair reads 0% at 0.6 cm, 77-80% at 1.4 cm)
        double m_mvga_dup_frac{0.7};          ///< op1/op3 overlap-fraction accept threshold
        double m_mvga_dup_angle{20};          ///< op1 near-parallel guard, deg (chord-vs-chord, folded to [0,90]); duplicates run (anti)parallel (13/7-11 deg measured), a genuine small-opening V must not merge; <= 0 disables
        double m_mvga_bridge_mip{0.5};        ///< op2 median-dQ/dx ceiling, x m_mip_dqdx_median.  Measured (268067 TRACE probe): the charge-less bridge 15005 reads 0.436 INTERNAL (its Bee display "0.1 MIP" was skewed by the affine q = dQ*0.1-1000 display transform), the genuine middle track 1.290 -- 0.5 separates them; the round-1 guess 0.33 missed the defining case.
        double m_mvga_reconnect{5*units::cm}; ///< op2 max direct-reconnect distance for a disconnected side (268067: V_B at 4.2 cm)
        double m_mvga_stub{2*units::cm};      ///< op3 terminal-stub length ceiling
        int    m_mvga_stub_pts{4};            ///< op3 point-degeneracy sub-gate: <= this many valid fits (overlap fractions are meaningless at 3-4 points)
        double m_mvga_reseat_angle{150};      ///< op3 re-seat collinearity threshold vs a sibling prong, deg (131357 measures ~180; 175 would be the final_1p analogue but near-corner arms curve)
        double m_mvga_satellite{0};           ///< op3 satellite-anchor radius around the main vertex (round 3, doc pr/51); 0 = main-vertex-only, byte-identical to round 2
        // doc sbnd_xin/docs/pr/85: op3 interposed-stub absorb.  op3 above is
        // TERMINAL-only -- its degree(far)==1 line can never reach an
        // INTERPOSED stub, whose far vertex carries the real prong(s) and so
        // has degree >= 2.  The pr/85 census over 462 hand-scanned events
        // measures that class at 11 of 78 near-vertex stubs (= the whole of
        // mode 1a-VIA, 21 events: a >= 10 cm prong connected to the clicked
        // neutrino vertex only through a sub-3 cm segment, created after
        // every examine_* pass by snap_main_vertex_to_kink and the final
        // improve_vertex).  With m_mvga_interposed, a main-vertex-anchored
        // stub under m_mvga_stub whose far vertex is degree >= 2, not
        // kProtectedBreak, and collinear within m_mvga_interposed_angle of a
        // prong at the far end (the stub is the missing last piece of a
        // through-going prong, not one arm of a genuine V) has ALL far
        // prongs spliced through the stub's wcpts onto the main vertex
        // (carry_prong_verify / carry_prong_execute, PRSegmentFunctions.h:
        // SegmentPtr identity preserved, all-or-nothing pre-verification),
        // the stub removed and the far vertex, now degree 0, dropped.  The
        // main vertex never moves; no new segment is created (the `created`
        // exemption set is untouched); one edit per restart under the same
        // kEditCap.  Deliberately out of reach: degree-1 main-vertex anchors
        // (evt360535's shape -- incident.size() < 2 gates the anchor loop).
        // false (default) => the degree(far)==1 line short-circuits exactly
        // as before => byte-identical.
        bool   m_mvga_interposed{false};      ///< op3 interposed-stub absorb at the main-vertex anchor (doc pr/85)
        double m_mvga_interposed_angle{150};  ///< far-end collinearity gate, deg (mirrors m_mvga_reseat_angle's 150)
        // doc sbnd_xin/docs/pr/86 P1: a separate candidate ceiling for the
        // interposed SPLICE only.  The pr/86 sec 10.3 merged-prong census
        // measures two thirds of the straight-track-misses-vertex defect as
        // an interposed segment of 2.5-10 cm -- above m_mvga_stub -- while
        // the interposed-segment length distribution of reachable orphans
        // stops at 5 cm (sec 4).  m_mvga_stub itself must NOT be raised: it
        // also gates the terminal absorb, which is where the pr/85 sec 10.6
        // adverse movers live.  > m_mvga_stub => interposed candidates are
        // admitted up to this length at the main-vertex anchor (the terminal
        // branch re-applies m_mvga_stub); 0 (default) => ceiling ==
        // m_mvga_stub => byte-identical.
        double m_mvga_interposed_len{0};      ///< interposed-splice candidate ceiling, internal length units (doc pr/86); 0 = use m_mvga_stub
        // doc pr/86 P4: pr/85 sec 10.6 raised m_mvga_dup_frac to 0.8 to stop
        // four adverse absorbs -- every one at the MAIN anchor with d=0.00.
        // The same guard blocks the wanted evt30504 absorb, which is at a
        // SATELLITE anchor (d=1.26 cm): a satellite absorb deletes a stub
        // away from the main vertex and cannot move it directly.  This knob
        // restores a separate (typically pre-pr/85 0.7) overlap threshold at
        // satellite anchors only.  0 (default) => m_mvga_dup_frac applies
        // everywhere => byte-identical.
        double m_mvga_sat_dup_frac{0};        ///< satellite-anchor op3 overlap threshold (doc pr/86); 0 = use m_mvga_dup_frac
        // doc pr/86 P1b: admit degree-1 MAIN anchors into op3 for the
        // interposed splice only.  The incident>=2 anchor gate is a
        // terminal-absorb argument (shedding the only edge disconnects the
        // anchor); the splice re-attaches all far prongs, and 26 of the 86
        // sec 10.2 Class-B cases -- including b1<=1 cm events where the
        // vertex is exactly right -- sit at a degree-1 anchor.  The
        // terminal branch re-imposes >= 2.  false (default) => the old
        // anchor gate verbatim => byte-identical.
        bool   m_mvga_interposed_deg1{false}; ///< op3 interposed splice at degree-1 main anchors (doc pr/86)
        // doc sbnd_xin/docs/pr/86 sec 15 (round 2): the interposed splice is
        // a graph relabel -- carry_prong_execute CONCATENATES the stub's and
        // prong's wcpt chains, and the op4 refit seeds from wcpts
        // (TrackFitting organize_segments_path), so the kink at the old far
        // vertex survives the fit and the rendered trajectory never changes
        // (owner: 38856's large-angle turn still there).  The prototype
        // straightens when it merges through a vertex
        // (NeutrinoID_examine_structure.h examine_structure_2 lines 76-165:
        // straight line, 0.6 cm steps, is_good_point charge veto, Steiner
        // snap) and only concatenates already-collinear merges.  > 0 =>
        // after each op3 carry, the merged chain from the anchor to the
        // first point >= (stub length + this reach) along the chain is
        // re-derived with that recipe; charge veto fails => keep the
        // concatenated chain.  0 (default) => concatenation verbatim =>
        // byte-identical.
        double m_mvga_splice_straighten{0};   ///< op3 post-carry straighten reach past the junction, internal length units (doc pr/86 round 2); 0 = off
        // doc sbnd_xin/docs/pr/86 sec 15 (round 2): near-vertex approach
        // collapse (op3.5).  The 349945 shape: a long prong reaches the main
        // vertex only via a multi-segment zigzag (14 cm of polyline over a
        // 5.8 cm straight line) built by improve_vertex AFTER
        // examine_structure_2 last ran, and nothing after mvga refits.  > 0
        // => degree-2, non-protected junction vertices within this radius of
        // the main vertex are re-tested with the examine_structure_2 merge
        // (straight, charge-vetoed, Steiner-snapped replacement segment;
        // co-located-endpoint case skipped -- divergence from es2's B.7
        // vertex merge, stated in doc), iterated to a fixed point; fires
        // count into the op4 refit trigger.  0 (default) => pass skipped
        // entirely => byte-identical.
        double m_mvga_approach_collapse{0};   ///< op3.5 junction-collapse radius around the main vertex, internal length units (doc pr/86 round 2); 0 = off
        // doc pr/86 sec 15: is_good_point radius for the R1/R2 straight-chain
        // charge veto.  The prototype es2 recipe uses 0.2 cm, which can
        // never approve cutting a corner whose charge ridge deviates ~1.6 cm
        // from the straight chord (the 349945 zigzag, the owner's "direct
        // track") -- Stage A measured every long-stub straighten and both
        // 349945 collapse junctions charge-veto'd at 0.2.  0 (default) =>
        // the prototype 0.2 cm => recipe-faithful; > 0 => that radius.
        // Inert unless splice_straighten/approach_collapse is on.
        double m_mvga_straighten_radius{0};   ///< straight-chain veto radius, internal length units (doc pr/86 round 2); 0 = prototype 0.2 cm
        // doc sbnd_xin/docs/pr/99 round 2 -- op3.5 approach-collapse guards.
        // The pr/99 firing census (mcp2k production: 46 fires, median chord
        // 15.7 cm, max 338.5 cm vs the 5.8 cm design case) showed op3.5
        // off-envelope: it replaced a true 51 cm V with straight chords on
        // 285567 (the flip regression) and built 315167's 146 cm fake EM
        // trunk.  Three independent guards, each 0/false = legacy:
        // - ac_veto_radius: dedicated is_good_point radius for the COLLAPSE
        //   chord's charge veto.  Production m_mvga_straighten_radius=1.0cm
        //   (relaxed for the pr/86 R1 splice-straighten purpose) leaked into
        //   the collapse veto; the prototype ancestor examine_structure_2
        //   (NeutrinoID_examine_structure.h line ~131) checks
        //   is_good_point(test_p, 0.2*units::cm, 0, 0).  > 0 => that radius
        //   for op3.5 only (R1 splice-straighten keeps m_mvga_straighten_radius).
        // - ac_chord_max: decline a collapse whose replacement chord
        //   |vtx1-vtx2| exceeds this; the prototype capped the removed
        //   segment lengths at 5 cm (commented out in the ancestor, but the
        //   design envelope is few-cm zigzags, not half-meter trunks).
        // - ac_no_cascade: never collapse a candidate whose sg1/sg2 is
        //   itself a `created` product (285567's second fire consumed its
        //   own nfits=0 chord, 53.6 -> 58.3 cm).
        double m_mvga_ac_veto_radius{0};      ///< op3.5-only charge-veto radius, internal length units (doc pr/99 round 2); 0 = legacy (straighten_radius rule)
        double m_mvga_ac_chord_max{0};        ///< op3.5 replacement-chord length cap, internal length units (doc pr/99 round 2); 0 = no cap
        bool   m_mvga_ac_no_cascade{false};   ///< op3.5: skip candidates touching `created` products (doc pr/99 round 2); false = legacy
        // doc sbnd_xin/docs/pr/99 round 2 -- op1-post charge second-opinion.
        // 70084: a 15.7 cm charge-starved chord rode a real prong at
        // overlap 0.87 but the pair's ~30 deg chord opening angle declined
        // the merge -- the angle guard exists so a genuine small-opening V
        // never merges, and here it protected a ghost.  The refit SPLITS
        // the shared corridor's charge across the pair (70084 measured
        // median-dQ/dx ratios 1.16/0.62 vs m_mip_dqdx_median), so the
        // discriminator is the op1-proj-style pair ASYMMETRY
        // (min/max <= starved_asym; 70084 reads 0.53, op1-proj production
        // gate 0.55) AND an absolute cap on the loser
        // (min <= starved_mip; a genuine proton+MIP V's muon reads ~1.0
        // and is never mistaken for a starved chord).  Both knobs > 0
        // required; the LOWER-charge member is deleted (never the healthy
        // one -- the pr/96 F3 safe direction).  Members without valid fits
        // are skipped (no charge verdict possible, the op2 rule).  Either
        // knob 0 (default) => the decline stands => byte-identical.
        double m_mvga_dup_starved_asym{0};    ///< op1-post starved-member override: pair min/max dQ/dx asymmetry gate (doc pr/99 round 2); 0 = off
        double m_mvga_dup_starved_mip{0};     ///< op1-post starved-member override: absolute cap on the loser, ratio vs m_mip_dqdx_median (doc pr/99 round 2); 0 = off
        double m_mvga_dup_starved_span{0};    ///< op1-post starved-member override: pair min/max LENGTH comparability floor (doc pr/99 round 2); 0 = no span test
        // doc sbnd_xin/docs/pr/83 round 3 (sec 9.6/sec 8.5): the duplicate-
        // corridor family.  op1 is a global graph-correctness check ("no two
        // segments double-book the same charge"), not inherently a
        // near-vertex one like op2/op3 -- but it inherits m_mvga_radius and
        // m_mvga_dup_frac, both tuned for other purposes (the radius bounds
        // the expensive near-vertex audit; 0.8 closed a pr/86 marginal-
        // overlap gap at the MAIN anchor).  Sec 9.2 measured 6 of 9 class-B
        // duplicates 5-7 cm OUTSIDE the 15 cm scope; sec 9.3 measured 1
        // declined at 0.74 by the raised 0.8.  These knobs decouple op1's
        // scope and threshold from op2/op3's without re-litigating either
        // tuning.  0 (default) => the shared member applies => byte-identical.
        double m_mvga_op1_radius{0};          ///< op1-only scope radius (doc pr/83 r3); 0 = use m_mvga_radius, -1 = unscoped (whole main cluster)
        double m_mvga_op1_dup_frac{0};        ///< op1-only overlap threshold (doc pr/83 r3); 0 = use m_mvga_dup_frac
        // doc pr/83 r3 (sec 8, class A): the pr/85 interposed carry re-
        // attaches N far prongs directly to the anchor and the pr/86
        // splice-straighten re-derives each carried prong's near-anchor
        // stretch over the same 26-35 cm reach -- N prongs acquire the SAME
        // trunk geometry (138009: six stacked electrons, 204 cm of fitted
        // length on a 43 cm trunk, kine_reco_Enu double-counted 1973 vs 1441
        // MeV).  op1 runs BEFORE op3 and its `created` exemption makes the
        // carried products structurally invisible (sec 8.3).  true => one
        // additional op1-style duplicate-corridor pass AFTER the op3/op3.5
        // interleave concludes (before the op4 refit), with the `created`
        // exemption lifted, at op1's effective radius/threshold.  Benign
        // carries overlap nothing and are untouched -- no pr/86 giveback at
        // the 33 non-stacking multi-carry sites (sec 8.5's cost analysis of
        // the blanket cap).  false (default) => pass skipped => byte-identical.
        bool   m_mvga_op1_post{false};        ///< post-op3 duplicate-corridor pass incl. created segments (doc pr/83 r3, class A)
        // doc pr/83 r3 (sec 8.5): the blanket fallback -- decline the
        // interposed absorb/splice outright when it would carry more than
        // this many far prongs (all 8 class-A events had carried >= 2; a cap
        // of 1 keeps the stub as the shared trunk, the physically correct
        // topology).  Costs pr/86 benefit at every benign multi-carry site,
        // so m_mvga_op1_post is the primary fix; this ships only if that
        // proves insufficient.  0 (default) => unlimited => byte-identical.
        int    m_mvga_carry_max{0};           ///< op3 interposed-carry prong-count ceiling (doc pr/83 r3); 0 = unlimited
        // doc pr/83 r3 (sec 9.5, Mechanism C + the 359980 follow-up): a
        // cluster that went through find_proto_vertex but is not the final
        // main cluster keeps its segments in the output yet never receives
        // a duplicate-corridor pass.  Two ways in: swap_main_cluster
        // re-points pattern recognition away from it (350935: an
        // overlap-1.00 dup born in find_proto_vertex survives to Bee
        // untouched), or it was a candidate that lost the main-cluster
        // contest with no swap at all (359980: dup on cluster 75, main is
        // 21).  true => one unscoped op1-style duplicate-corridor pass (no
        // vertex to center a radius on -- these clusters have no main
        // vertex), with a refit if anything merged: inside
        // swap_main_cluster on the abandoned cluster, and once more over
        // every non-main cluster before shower_clustering_with_nv consumes
        // them (TaggerCheckNeutrino.cxx).  false (default) => byte-identical.
        bool   m_swap_orphan_dup_audit{false}; ///< dup-audit non-main clusters: at swap + pre-shower sweep (doc pr/83 r3)

        // ---- doc sbnd_xin/docs/pr/83 round 4: projective duplicate collapse --
        // A 1-track-1-shower stem can be reported as TWO main-vertex tracks
        // that overlap in >=2 of the 3 wire views while separating in 3D
        // (the fitter places the same 2D charge on two 3D interpretations;
        // the starved one reads far-below-MIP stem dQ/dx).  3D corridor
        // overlap reads 0.14-0.58 (below every op1 gate), so round 3 never
        // fires (138009 12094/12095, 168596 14168/14172, 74544 12105/12107).
        // 0 = pass disabled = byte-identical legacy.  When > 0: minimum
        // 2nd-best per-view overlap fraction (views at wire_angles, coord
        // (x, cos(a)z - sin(a)y), tol = m_mvga_dup_tol).
        double m_mvga_proj_dup_frac{0};
        // Stem dQ/dx asymmetry gate: merge only when min/max of the two
        // members' dQ/dx over the first 8 cm from the main vertex is below
        // this ratio (measured: ghosts 0.08-0.28; real two-prong vertices
        // carry MIP-level charge on both prongs).  Only consulted when
        // m_mvga_proj_dup_frac > 0, so the default is inert.
        double m_mvga_proj_dqdx_ratio{0.4};
        // op1-proj's own chord-angle gate, degrees (doc pr/83 r4b, 284206:
        // the residual stem pair reads 22 deg -- just over op1's shared
        // 20 deg).  0 = use m_mvga_dup_angle (byte-identical legacy).
        double m_mvga_proj_angle{0};

        // ---- doc sbnd_xin/docs/pr/51 round 4: rough-path diagnostic probe --
        // Diagnostic-only TRACE instrumentation for the near-vertex
        // short-cut investigation (see TaggerCheckNeutrino.h for the full
        // rationale).  Never mutates a graph, segment, or fit; every line is
        // SPDLOG_LOGGER_TRACE.  false (default) => byte-identical.
        bool   m_rough_path_probe{false};

        // ---- doc sbnd_xin/docs/pr/51 round 5: steiner gap penalty -----------
        // The H1 fix the round-4 probe validated counterfactually: the
        // "steiner_graph" edge weight prices charge only at the two edge
        // endpoints with a hard [0.8, 1.2] dynamic range
        // (SteinerGrapher.cxx create_enhanced_steiner_graph), so a
        // gap-spanning chord beats following charge around any turn sharper
        // than ~150 deg.  When m_steiner_gap_penalty > 0, do_rough_path
        // lazily derives a per-cluster graph flavor "steiner_graph_gap" --
        // same topology and vertex indexing as "steiner_graph", every edge
        // longer than m_sgp_min_edge re-weighted
        //     w' = w * (1 + scale * bad_fraction)
        // with bad_fraction = (n_unsup + alpha*n_dead)/n from sampling the
        // edge interior at m_sgp_sample_step with
        // Grouping::test_good_point(p_raw, apa, face, m_sgp_point_radius, 0)
        // -- exactly the round-4 probe's P3 scan -- and routes on that
        // flavor.  Nothing else moves: TaggerCheckSTM's fork, TrackFitting
        // and every other "steiner_graph" consumer keep the base graph.
        // 0 (default) => first-line return => the flavor is never built =>
        // byte-identical.
        double m_steiner_gap_penalty{0};        ///< P3-ladder scale; 0 = off
        double m_sgp_dead_alpha{0.25};          ///< dead-sample weight in bad_fraction
        double m_sgp_min_edge{0.5*units::cm};   ///< edges shorter are never scanned/penalized
        double m_sgp_sample_step{0.3*units::cm};///< edge-interior sampling step
        double m_sgp_point_radius{0.2*units::cm};///< test_good_point radius (ch_range stays 0)
        // doc sbnd_xin/docs/pr/73: per-edge diagnostic sentinel for the scan
        // above.  When true, ensure_steiner_gap_graph emits one DEBUG line
        // per SCANNED edge -- endpoints, midpoint, w, bad, the two recovered
        // vertex charges, deficit -- so the penalized edges can be located in
        // space rather than counted.  Answers "are the weak-charge edges
        // inside the isochronous ribbon?", which doc pr/73 sec 4.8 leaves
        // open.  Log-only: reads nothing it does not already read and writes
        // no graph.  false (default) => not a single extra statement runs.
        bool   m_sgp_edge_probe{false};         ///< per-edge DEBUG sentinel; diagnostic only
        // doc sbnd_xin/docs/pr/51 round 6: weak-charge deficit term.  The
        // round-5 unsupported-fraction penalty is blind to chords whose
        // interior is image-SUPPORTED but charge-poor (18259-131357 trunk
        // chord, 18255-506746 hairpin connector: probe P3 ladders are
        // scale-invariant 0..10 on both).  When m_sgp_weak_scale > 0 (and
        // the gap flavor is enabled, i.e. m_steiner_gap_penalty > 0), each
        // scanned edge additionally pays
        //     w' = w * (1 + gap_scale*bad + weak_scale*deficit),
        //     deficit = 0.5*(max(0,1-q_s/qref) + max(0,1-q_t/qref)),
        // with q_* the endpoint steiner-vertex charges recovered via
        // calc_charge_wcp(idx, 4000, false) -- the same call the production
        // steiner edge weighting used (CreateSteinerGraph.cxx
        // create_steiner_tree(..., false, ...) + pr/29 D2 forwarding).
        // 0 (default) => the round-5 reweight path runs verbatim =>
        // byte-identical gap flavor.
        double m_sgp_weak_scale{0};             ///< weak-charge penalty scale; 0 = off
        double m_sgp_weak_qref{2000};           ///< charge ref (calc_charge_wcp RMS units); deficit=0 at/above
        // doc sbnd_xin/docs/pr/73 round 2, fix F3a: bound what the penalty may
        // do to the ROUTE.  On 18255-57903 the round-6 weighting talks
        // do_rough_path into a +4.87 cm excursion through an isochronous ghost
        // ribbon, relocating the main vertex and forcing a 5.9 cm bow that
        // nothing downstream can remove (multi_trajectory_fit PINS both segment
        // endpoints to the vertex fit points, TrackFitting.cxx:4246-4259).  The
        // two events round 6 was built to fix need 98 route-moving calls whose
        // largest excursion is 2.57 cm, so a cap separates them from the damage
        // (4.85 cm on the causal call) where capping the DETOUR does not:
        // in percent the detour criterion is inverted (the fixes routinely need
        // 30-40 % on short paths) and in cm its margin is only 1.17x.
        //
        // When the cap is >= 0 and the gap flavor is in use, do_rough_path also
        // routes on the untouched base flavor and, if the penalized route ever
        // gets further than the cap from it, KEEPS THE BASE ROUTE.
        //
        // NOTE the off-test is `< 0`, NOT the `<= 0` used by the sgp scale
        // knobs above: 0 is a meaningful cap here (reject any excursion at
        // all), so it cannot double as the off value.  -1*units::cm stays
        // negative through the cm->internal conversion in TaggerCheckNeutrino.
        // Default -1 => the base flavor is never routed and not one extra
        // statement runs => byte-identical.
        double m_sgp_max_sep{-1*units::cm};     ///< route excursion cap; < 0 = off (unbounded)
        IDetectorVolumes::pointer m_sgp_dv{nullptr};   ///< pushed from TaggerCheckNeutrino (NeedDV)
        IPCTransformSet::pointer  m_sgp_pcts{nullptr}; ///< pushed from TaggerCheckNeutrino (NeedPCTS)

        // ---- doc sbnd_xin/docs/pr/30 §11: four port-fidelity knobs ----------
        //
        // P1 / F1 -- `flag_exclusion` on do_multi_tracking.
        // 28 of the 30 live prototype call sites pass flag_exclusion = true
        // (NeutrinoID_proto_vertex.h:108,140,170,319,355,415,1363,1479,2358,
        //  2471,2537,2604,2724,2813,2995 + 13 more in improve_vertex /
        //  final_structure / examine_structure / track_shower); the toolkit
        // passes false at all 31.  Both signatures DEFAULT the flag to false
        // (prototype PR3DCluster.h:182, toolkit TrackFitting.h:446), so the
        // prototype's 28 are explicit overrides the port stopped writing --
        // the signature of a dropped argument, not a design choice.  With
        // exclusion on, form_map_graph calls update_association
        // (TrackFitting.cxx:3187, a faithful port of
        // PR3DCluster_multi_track_fitting.h:820/970-1075) which strips from
        // this segment's 2-D associations the (wire,tick) cells belonging to
        // OTHER segments.  Off, two segments crossing the same wire region
        // both claim the same charge.
        // C++ default false => every call site keeps passing false =>
        // byte-identical to the pre-pr/30 behaviour.
        // NOT applied at NeutrinoPatternBase.cxx:1491/1507 (break_segments --
        // the toolkit's correct match to the prototype's own two false sites,
        // NeutrinoID_proto_vertex.h:722/751, "fit dQ/dx here, do not exclude
        // others") nor at NeutrinoVertexFinder.cxx:500 (a single-segment local
        // graph, where exclusion has nothing to exclude).
        bool   m_fit_exclusion{false};

        // P8 -- endpoint-consistency check on PR::add_segment.
        // The prototype refuses a vertex/segment connection whose vertex wcpt
        // index is neither the segment's front nor its back
        // (NeutrinoID_proto_vertex.h:1952-1956, "Error! Vertex and Segment
        // does not match").  The toolkit has no index to compare, and no
        // positional analogue was written, so the check is absent rather than
        // translated.  It is the invariant find_vertices() silently depends on
        // (PRGraph.cxx:105-141 orders its pair by distance to the segment's
        // FIRST wcpt, which is only meaningful if both vertices really do sit
        // at the segment's two ends).
        // The WARN is unconditional and log-only.  m_graph_endpoint_strict
        // additionally REFUSES the connection, which changes the graph =>
        // C++ default false => byte-identical.
        // m_graph_endpoint_tol is the positional stand-in for wcpt-index
        // equality; 0.3 cm is well under the Steiner point spacing so it
        // cannot merge distinct ends, and well over FP noise in a copied point.
        bool   m_graph_endpoint_strict{false};
        double m_graph_endpoint_tol{0.3 * units::cm};

        // F2 / P9 -- polarity of the out-of-detector-volume point guard.
        // dv->contained_by(p) returns apa()/face() == -1 outside every TPC and
        // the downstream m_trigger_offsets.at(-1) / m_anodes.at(-1) then
        // throws, so the guard itself is REQUIRED (class-C crash, doc pr/11
        // §6.3; SBND MCP2025C evt 49951).  What was never chosen is what the
        // skipped point MEANS -- and the three in-scope sites silently picked
        // three different answers, all three the opposite of what the
        // prototype's own helper returns for a point with no readout:
        //   NeutrinoOtherSegments.cxx  modify_segment_isochronous
        //     toolkit: cannot increment n_bad          => "good/connected"
        //     prototype: is_good_point -> get_closest_points empty AND
        //       get_closest_dead_chs false => num_planes 0 => FALSE
        //       (ToyCTPointCloud.cxx:399-431)          => n_bad++  ("bad")
        //   NeutrinoStructureExaminer.cxx  examine_vertices_1p
        //     toolkit: flag_dead stays true            => "dead"
        //     prototype: get_closest_dead_chs -> channel absent from
        //       dead_uchs/vchs/wchs => FALSE           => flag_dead=false
        //   NeutrinoStructureExaminer.cxx  examine_vertices_3
        //     toolkit: cannot contribute num_unique    => "not unique" (the
        //       segment is then REMOVED from the graph)
        //     prototype: get_closest_2d_dis is a pure kd-tree 2-D distance
        //       with no volume check (ProtoSegment.cxx:1094-1103), so the
        //       point participates normally.
        // The first two restorations are EXACT (the helper demonstrably
        // returns false); the third is a DIRECTIONAL INFERENCE -- the
        // prototype computes a real distance rather than voting, and
        // "outside every TPC => far from every other segment in some view"
        // is very likely but not proven.  It is also the non-destructive
        // default, which is the right tiebreak for a branch whose other
        // outcome deletes a segment.  Flagged as such at the site.
        // C++ default false => all three keep today's polarity =>
        // byte-identical.
        bool   m_oov_prototype_parity{false};

        // P2 -- local-PCA refinement of the first segment's two boundary
        // endpoints (NeutrinoPatternBase.cxx, commit 1eb097a9 "fix a bug and
        // randomness").  Toolkit-only: the prototype takes
        // get_two_boundary_wcps(2) and uses it as-is
        // (NeutrinoID_proto_vertex.h:426).  Unlike the knobs above this is
        // ALREADY the production behaviour, so its default is TRUE and OFF is
        // the new option -- the point of the knob is that an unconditional,
        // un-knobbed departure in the single most load-bearing function of the
        // stage could not be measured at all.  Everything downstream of
        // init_first_segment derives from these two points.
        // C++ default true => byte-identical to the pre-pr/30 behaviour.
        // (Inert when m_iso_endpoint fires -- that branch bypasses the whole
        // legacy block including this refinement.)
        bool   m_first_seg_local_pca{true};

        // P4 -- the extra acceptance clause in find_other_segments.
        // Prototype NeutrinoID_proto_vertex.h:1368 accepts
        //   length > 30 cm || (direct/length < 0.78 && length > 10 cm &&
        //                      median_dQdx/MIP > 1.6)
        // the toolkit adds
        //   || (direct/length < 0.72 && length > 15 cm && median_dQdx/MIP > 1.05)
        // -- a strict widening (less curved, longer, less ionising), so it can
        // only ADD segments.  Same default logic as P2: this is production
        // today, so default TRUE, and OFF restores the prototype's clause.
        // C++ default true => byte-identical.
        bool   m_other_seg_relaxed_accept{true};

        // doc sbnd_xin/docs/pr/54 -- keep well-supported isolated residual
        // segments in find_other_segments (18255-142421 "missing gammas": a
        // 2930-point separated EM component of the main cluster is fragmented
        // by the terminal-graph partition, its pieces fit and then discarded
        // because neither endpoint touches the existing graph and the
        // isochronous snap does not apply).  The prototype has the same
        // discard: NeutrinoID_proto_vertex.h:1470-1475 pushes such candidates
        // into residual_segment_candidates, which is write-only -- never
        // consumed anywhere -- so this keep path is a toolkit-only extension
        // of an unfinished prototype feature, not a parity fix.  When ON, an
        // isolated candidate whose terminal-graph component has at least
        // min_points points AND whose fitted track is at least min_length
        // long is added to the graph as its own disconnected piece (own two
        // endpoint vertices) and refit jointly with the cluster, exactly like
        // the isochronous-accepted branch.  C++ default false => discard,
        // byte-identical.
        bool   m_other_seg_keep_isolated{false};
        int    m_other_seg_keep_isolated_min_points{25};
        double m_other_seg_keep_isolated_min_length{3.0 * units::cm}; // internal units

        // doc sbnd_xin/docs/pr/67 round 3 (S2) -- the size gate on the first
        // clause of the isochronous-snap guard in find_other_segments.  The
        // machinery behind that guard (modify_vertex_isochronous /
        // modify_segment_isochronous) is the only thing that ATTACHES an
        // isochronously-displaced branch to its parent; at the legacy 10 cm it
        // never runs on the short branches of doc pr/67 (dir_mag 4.3-4.7 cm).
        // Lowering the SIZE gate does not relax the ISOCHRONOUS requirement:
        // the vertex path keeps the caller's 15 deg perpendicular test and the
        // segment path keeps its own angle_cut test.  The second clause
        // (8 cm + 13 cm track length) and the >18 cm / >36 cm widening tiers
        // are deliberately untouched.  C++ default 10 cm => byte-identical.
        double m_iso_snap_min_dir_mag{10 * units::cm}; // internal units

        // doc sbnd_xin/docs/pr/59 round 2 -- gates reassociate_cluster_orphans
        // (checked inside that function, called unconditionally, matching the
        // m_main_vertex_graph_audit idiom).  C++ default false => legacy
        // (orphaned associate_points cloud stays null), byte-identical.
        bool   m_assoc_full_recluster{false};

        // doc sbnd_xin/docs/pr/64 round 7 -- passed as the trailing
        // reassign_orphans arg of clustering_points_segments (both live call
        // sites: clustering_points and reassociate_cluster_orphans's recluster
        // above, so a rescued segment gets the same association rule as the
        // main pass).  18259-18625: a 12-18 pt blob at PF segment 126042's own
        // fit endpoint (the photon-conversion vertex) is present in img charge
        // but absent from shower_track/associate_points -- Stage C's ghost
        // removal only ever DROPS a losing point, it never hands it to the
        // segment that actually wins the global 2D projection contest.  When
        // on, such a point is reassigned to the winning same-cluster segment
        // (never a different cluster's -- cross-cluster ghost rejection is
        // unaffected by construction).  C++ default false = legacy = drop,
        // byte-identical.  See PRSegmentFunctions.cxx's doc comment on
        // clustering_points_segments for the full mechanism and
        // WCT_PR64_ORPHAN_CENSUS for the log-only diagnostic.
        bool   m_assoc_reassign_orphans{false};

        // doc sbnd_xin/docs/pr/64 round 8 -- checked inside
        // examine_structure_final_1/_1p/_3 (NeutrinoStructureExaminer.cxx),
        // called unconditionally from determine_main_vertex.  Those passes
        // merge a short/duplicate/degenerate segment into a surviving
        // neighbor: the survivor's wcpts/fits/"main" cloud are rebuilt, but
        // its "associate_points" cloud (the actual charge/blob association
        // from clustering_points) is left untouched -- so if the DELETED
        // segment held associate_points, they are simply discarded with no
        // replacement (18259-18625: a 33-point, 6-wcpt segment absorbed at
        // the main vertex by examine_structure_final_1p, its points -- incl.
        // the reported 12-pt blob at (142.1,78.3,176.5) -- gone from the
        // final PF/Bee dump).  reassociate_cluster_orphans (pr/59, live in
        // SBND production) exists exactly to catch a stale/incomplete
        // association after determine_main_vertex, but its any_orphan
        // trigger only fires when some CURRENT segment has a completely
        // EMPTY associate_points cloud -- a survivor that already had SOME
        // points of its own never trips it, even though those points no
        // longer cover its newly-extended geometry.  When on: at each of
        // the four merge/delete sites, if the segment being removed had a
        // non-empty associate_points cloud, the designated survivor's
        // associate_points is cleared (set to null) too, so any_orphan
        // correctly sees a gap and pr/59's existing full-cluster
        // re-clustering re-derives a geometry-consistent association for
        // the whole cluster -- no new competition logic, reuses the
        // already-validated pr/59 machinery.  C++ default false = legacy =
        // survivor keeps its stale cloud, byte-identical.
        bool   m_assoc_clear_on_merge{false};

        // doc sbnd_xin/docs/pr/45 -- empty-2D-tree sentinel guard in
        // find_other_segments (SBND 18255-56463 cluster 14, the 30 cm
        // isochronous tail beyond segment 14006's end).
        //
        // DynamicPointCloud::get_closest_2d_point_info returns distance
        // -1.0 when the per-(plane,face,apa) 2D kd-tree is empty
        // (DynamicPointCloud.cxx knn empty-result branch), which happens for
        // any segment with zero fit points on the queried face -- routine in
        // SBND cathode-crossing clusters.  find_other_segments compares that
        // sentinel against thresholds as if it were a real distance:
        // -1 < scaling_2d*search_range "covers" all three planes at once, so
        // ONE far-TPC segment tags the ENTIRE near face as explained and no
        // residual component can ever form there (measured: all 194 tail-box
        // steiner points tagged with u=v=w=-1, 3D distance 3.9-28.9 cm).
        // The prototype cannot hit this: uBooNE is single-face, so
        // ProtoSegment::get_closest_2d_dis always has points.  Toolkit-only
        // port divergence, not a prototype limitation.
        // When ON, a negative 2D distance in find_other_segments' three
        // comparison sites (tagging loop, component not_faked census,
        // re-evaluation loop) is treated as "no projection information"
        // (1e9) instead of "distance zero".
        // C++ default false => byte-identical legacy behaviour.
        bool   m_other_seg_empty_2d_guard{false};

        // doc sbnd_xin/docs/pr/31 §11 -- F2 (was P2): the stage-3 shower
        // direction call site that has no prototype counterpart.
        //
        // determine_direction's kShowerTopology branch
        // (NeutrinoTrackShowerSep.cxx:121-135) calls
        // segment_determine_shower_direction -- 305 lines of associated-point
        // PCA, spread profiling and endpoint comparison
        // (PRSegmentFunctions.cxx:2208-2512) -- whose prototype original
        // ProtoSegment::determine_shower_direction() is called from exactly
        // ONE place in the whole prototype tree,
        // NeutrinoID_track_shower.h:1532, inside
        // compare_main_vertices_all_showers: stage 4, all-showers path only.
        // What the prototype runs HERE is determine_dir_shower_topology
        // (ProtoSegment.cxx:1677-1710), four live lines that set
        // particle_type=11 and particle_mass and DO NOT touch flag_dir -- so
        // in the prototype a topology shower leaves stage 3 with whatever
        // direction is_shower_topology's forward/backward large-spread
        // comparison left (ProtoSegment.cxx:523-527).
        //
        // The call mutates exactly one thing, twice: segment->dirsign(0) at
        // entry (:2209) and segment->dirsign(flag_dir) at the end (:2509).
        // Nothing else on the segment is written and the return value is
        // discarded at the call site, so suppressing it suppresses precisely
        // the direction overwrite.
        //
        // NOT filed as a clear porting bug: the prototype's own function
        // carries a "// hack for now" comment above particle_type = 11 and
        // has both of its direction blocks commented out, i.e. it is
        // self-declared provisional, and the toolkit's PCA may well be the
        // better physics.  What is certain is that the substitution is
        // unconditional, undeclared and reaches this stage's in/out maps and
        // stage 4's vertex scorer through dirsign.  This knob exists so the
        // question can be measured instead of argued.
        //
        // TRUE = prototype behaviour (skip the call).  C++ default FALSE =>
        // today's path => byte-identical.
        //
        // Residual when ON, stated because it is not full parity: the
        // prototype clears flag_dir at is_shower_topology's entry
        // (ProtoSegment.cxx:321, before its early returns) while the toolkit's
        // segment_is_shower_topology skips dirsign entirely on its four early
        // returns (PRSegmentFunctions.cxx:2518-2529).  Today that hole is
        // masked by this call's entry-side dirsign(0); with the knob ON it is
        // exposed, but only for a segment carrying a STALE kShowerTopology
        // flag -- i.e. only in combination with F3 (P13), which is what closes
        // it.  Do not fix that here.
        bool   m_shower_topo_proto_dir{false};

        // ---- doc sbnd_xin/docs/pr/32 §11 -- the four kept findings of the
        // stage-4 (neutrino vertex identification) port audit.  All four are
        // C++ default FALSE = today's path = byte-identical; the SBND
        // operating point turns them on in wct-pr-perevt.jsonnet.

        // F1 (was P1).  calc_conflict_maps and compare_main_vertices_all_showers
        // measure their direction vectors, PCA projection and z tie-breaks from
        // `vtx->wcpt().point` -- the vertex fit SNAPPED to the nearest Steiner
        // node -- where the prototype uses `get_fit_pt()`, the continuous fit
        // (NeutrinoID_track_shower.h:1804/1808 for the two direction vectors,
        // :1468/:1480/:1504-1505/:1542/:1551/:1567/:1574 in the all-showers
        // scorer).  Eleven expressions in two functions; the same file already
        // reads `fit().valid() ? fit().point : wcpt().point` at 23 other
        // expressions, so the file is inconsistent with itself and the
        // prototype has exactly one convention.
        //
        // The gap is Vertex::fit_distance() -- a quantisation onto the Steiner
        // lattice, not a fitted-vs-unfitted swap: MyFCN::UpdateInfo writes both
        // members (MyFCN.cxx:487 the fit, :496 the snap).  It matters because
        // these angles feed the <35/<70/<85/<110 conflict ladder whose rungs
        // are worth +5/+3/+1/+0.25 before the /4, against topology terms spaced
        // 0.125 apart.
        //
        // NOT byte-identical when ON, and known so in advance: fit_distance()
        // is nonzero on 127/127 vertices of evt 388 (doc pr/28 A4).
        //
        // The fallback to wcpt() is load-bearing and stays: PR::Vertex has no
        // constructor initialising m_fit.point from m_wcpt.point the way
        // ProtoVertex's does (ProtoVertex.cxx:17-19), so a vertex created after
        // the last fit and never fitted carries fit().point == (0,0,0).
        bool   m_vertex_dir_use_fit_point{false};

        // F2 (was P3).  The improve_vertex shower-trajectory recheck.  Three
        // changes that must move together:
        //   (a) prototype's two outer gates READ the stored flag
        //       (NeutrinoID_improve_vertex.h:248 and :287); the toolkit
        //       RECOMPUTES at both (NeutrinoVertexFinder.cxx:2337, :2409);
        //   (b) the inner test runs at 1.0 cm where the prototype uses the
        //       10 cm default (ProtoSegment.h:85);
        //   (c) the inner test passes m_mip_dqdx_median where the prototype's
        //       is_shower_trajectory divides by 50000/units::cm internally
        //       (ProtoSegment.cxx:566), i.e. the m_mip_dqdx analog.
        //
        // Fixing (b)+(c) ALONE makes the block dead code: with the outer gate
        // also recomputing at (10 cm, m_mip_dqdx), the inner
        // `!segment_is_shower_trajectory(...)` negates a condition established
        // one line earlier and never fires.  The 1.0 cm step is what keeps the
        // block alive after (a); it is not an independent typo.
        //
        // (a) is only possible once the flag can be cleared -- see
        // PR::g_shower_traj_refresh_flag, which this knob also drives.
        bool   m_shower_traj_recheck_parity{false};

        // F3 (was P7).  compare_main_vertices guards its proton-topology
        // pre-pass and its z-prior/per-segment ladder on descriptor_valid()
        // but NOT the min_z scan, the fiducial +0.5, the conflict penalty or
        // the argmax.  An invalid-descriptor candidate therefore reaches the
        // argmax carrying `0 + (0.5 if in FV) - conflicts/4`, and real
        // candidates routinely score negative, so it can win.  The prototype
        // has no descriptor concept and no such path.
        //
        // Believed unreachable today -- candidates come from ordered_nodes()
        // and nothing between collection and use removes a vertex -- but that
        // is a control-flow argument, not a measurement, which is why the drop
        // is counted (see the PR32AUDIT log line) rather than assumed zero.
        bool   m_main_vertex_require_descriptor{false};

        // F4 (was P12).  The prototype records the per-cluster candidate list
        // before filtering (NeutrinoID_track_shower.h:1332) and exposes it
        // (NeutrinoID.h:1720); the toolkit has no equivalent, which is why doc
        // pr/27 §6 records "candidate list, per-cluster and global -> nothing".
        // Every prototype consumer is an app-level output-tree filler, so this
        // is DIAGNOSTIC ONLY -- no algorithm reads it on either side.
        //
        // Ported as PR::VertexFlags::kMainCandidate rather than a new container
        // so it needs no plumbing and travels with the vertex; PrDisplayDump
        // emits it as "main_candidate".
        bool   m_main_vertex_candidate_flag{false};

        // ---- doc sbnd_xin/docs/pr/75 -- the vertex scoreboard doc pr/52 §5.1
        // asked for.  Records, per event, the numbers the two vertex selectors
        // actually compared: compare_main_vertices' additive score per
        // candidate, and the DL rerank's top-K voxels + seven composite terms
        // + accept/reject decision.  Pure recording -- no decision reads it,
        // and every read inside the guarded blocks must be const (see
        // compare_main_vertices: map_vertex_num is a std::map, so an
        // operator[] read on an unscored candidate would INSERT and mutate a
        // container the legacy path then walks).  C++ default false => the
        // board stays empty => byte-identical.
        bool   m_vertex_scoreboard{false};
        VertexScoreboard m_vtx_board;

        // ---- doc sbnd_xin/docs/pr/79 §10 -- dl_vtx_harvest.  ALWAYS the AND
        // vertex_scoreboard && dl_vtx_harvest (TaggerCheckNeutrino pushes the
        // conjunction), so every m_vtx_harvest fill site may assume the board
        // is active.  Recording only -- no decision reads it, and the same
        // const-read discipline as m_vertex_scoreboard applies (never
        // operator[] on map_vertex_num, never iterate pointer-keyed maps).
        bool   m_vtx_harvest{false};

        // ---- doc sbnd_xin/docs/pr/31 §12 -- the §10.12 port-fidelity round:
        // the five surviving bug-class findings of the topology/PID/direction
        // audit (F5, F6, F3, F1, F4) plus the deliberately-dormant F7.  All
        // C++ default FALSE = today's path = byte-identical; the SBND
        // operating point in wct-pr-perevt.jsonnet decides which are on.

        // F5 (was P6).  find_cont_muon_segment_nue's hoisted dir3 falls back
        // to the 15 cm dir1 when sg_length <= 30 cm; the prototype always
        // compares two 30 cm directions when either segment is long
        // (NeutrinoID_track_shower.h:2402-2408).  TRUE = unconditional 30 cm.
        bool   m_cont_muon_dir3_30cm{false};

        // F6 (was P7, value half).  do_track_comp's empty-comparison-window
        // (and missing-dEdx) return declares "direction confirmed" (1.0); the
        // prototype's degenerate answer is abstain (0.0) -- executed, not
        // inferred (zero-bin TH1F KolmogorovTest returns 0 for both templates,
        // eval_ks_ratio's first gate fails).  Travels via
        // TrackPidOptions::track_comp_empty_abstain.
        bool   m_track_comp_empty_abstain{false};

        // F3 (was P13).  segment_is_shower_topology never clears
        // kShowerTopology (set-only tail) and its four early returns skip
        // dirsign(); the prototype clears both at entry, before its early
        // returns (ProtoSegment.cxx:319-321).  Re-entry on the same segment is
        // the normal path: stage 3 plus three stage-4 sites.  TRUE = clear the
        // flag bit and zero dirsign at entry (unset_flags -- other flags
        // survive).  Also the closer of m_shower_topo_proto_dir's stated
        // residual (see its comment above): with BOTH on, a segment whose
        // clouds hit the early returns leaves undirected instead of keeping a
        // stale direction.  Interaction watched when ON:
        // m_shower_topo_demote_len's demotion is no longer undone by a stale
        // flag on the next pass.
        bool   m_shower_topo_reset{false};

        // F1 (was P1 + P3's 4-momentum half + P4) -- the largest-reach item.
        // Fifteen reclassification sites in NeutrinoTrackShowerSep.cxx rewrite
        // the 4-momentum where the prototype writes type and mass and guards
        // the recompute on a previously computed energy
        // (if (get_particle_4mom(3)>0) cal_4mom()).  See reclass_4mom's
        // comment in NeutrinoTrackShowerSep.cxx for the three shapes and the
        // accident evidence (the guard survives at exactly one site, WITH a
        // comment paraphrasing it).  TRUE = preserve the existing 4-momentum;
        // recompute only where the prototype's guard passes.  Moves
        // kine_reco_Enu directly -- validate alone (doc pr/31 §10.2).
        bool   m_reclass_preserve_4mom{false};

        // doc sbnd_xin/docs/pr/40 round 2 F6.  reclass_pinfo (NeutrinoTrackShowerSep.cxx)
        // is reachable, with m_reclass_preserve_4mom true, on a segment that
        // never had a particle_info (had==false): the legacy code finishes by
        // overwriting the just-built (mass,0,0,0) placeholder to an all-zero
        // 4-vector, so ParticleInfo::kinetic_energy() (= e() - mass) reads
        // -mass -- a NEGATIVE kinetic energy wherever a caller (the Bee
        // PF-tree writer among them) reads it.  TRUE = leave the (mass,0,0,0)
        // placeholder alone in that case (KE == 0) instead of zeroing further.
        // Independent of m_reclass_preserve_4mom's own value (only reached
        // through one of its sub-branches) and of F4/F5.  C++ default false =
        // legacy = byte-identical.
        bool   m_reclass_never_computed_ke_floor{false};

        // F4 (was P8).  segment_determine_dir_track's median dQ/dx comes from
        // segment_median_dQ_dx's FILTERED rebuild while the PID receives the
        // local unfiltered vector (zeros kept) -- the toolkit disagreeing with
        // itself about the same pathological point.  The prototype takes the
        // median over the very vector it hands to do_track_pid at all three of
        // its sites.  TRUE = median over the local vector.  Travels via
        // TrackPidOptions::dir_track_median_local (and a dedicated
        // median-only forward into segment_determine_shower_direction's
        // short-segment interior call, which otherwise passes default
        // options).  The filtered helper itself is unchanged for every other
        // caller.
        bool   m_dir_track_median_local{false};

        // F7 (was P5; pr/30 F4's sibling) -- IMPLEMENTED BUT DELIBERATELY NOT
        // FLIPPED.  examine_all_showers' asymmetric 165/150-degree acceptance
        // pair is keyed on find_vertices().first, which the toolkit orders by
        // proximity to the segment's first fit point where the prototype
        // orders by vertex id.  TRUE = order the pair by get_graph_index() at
        // this one call site -- A deterministic convention, not provably the
        // prototype's (creation orders differ between the trees).  Stays OFF
        // pending pr/30 F4's adjudication, which owns the global
        // find_vertices-ordering question; flipping this alone would fix one
        // of three known order-sensitive callers.
        bool   m_examine_showers_vertex_by_index{false};

        // Proton-template direction vote (doc pr/8; default false = legacy).
        // Thresholds are initial values pending the pr/8 sec. 6 calibration.
        bool   m_proton_dir_vote{false};
        double m_proton_dir_score_max{0.25};
        double m_proton_dir_asym_min{1.3};
        // Endpoint-trim retry (doc sbnd_xin/docs/pr/9 sec. 6 F1).  C++ default
        // false => legacy abstention path, byte-identical.
        bool   m_endpoint_trim_retry{false};
        // Minimum WCPT-path length for a segment to count as a leg in the
        // fit_vertex position fit (doc sbnd_xin/docs/pr/9 sec. 11 F3c, owner
        // 2026-07-30, deliberate prototype divergence): very short vertex-
        // activity stubs otherwise drag the fitted vertex (evt 172230: a
        // 0.62 cm drift-blur stub moved the true vertex 2.4 cm).  Measured in
        // wcpt space (a fresh stub's FIT cloud spreads past 1 cm).  If >=3
        // legs survive, the fit runs on the survivors; if the exclusion
        // leaves <=2, the vertex fit is skipped (an effectively two-leg
        // vertex was already fit by the plain two-leg pass, and refitting on
        // stub-contaminated re-tracked clouds reproduces the drag).  Excluded
        // segments stay in the graph and the particle flow.
        // C++ default 0 => no filtering, legacy byte-identical.
        double m_fit_vertex_min_seg_length{0};

        // doc sbnd_xin/docs/pr/51 round 7 (follow-ups 1-5): robust vertex
        // fit.  MyFCN::UpdateInfo rewrites each long leg's first 4 cm as a
        // straight line to the current vertex, and the next AddSegment reads
        // that leg's direction from the (1.5, 6] cm window -- mostly that
        // same rewrite -- so the fit self-confirms wherever the vertex is
        // (117-evt census: 5-8 % of main vertices carry a leg whose
        // re-seat-free outer axis disagrees > 20 deg; worst 74 deg with a
        // 5 cm outer impact parameter).  When on, fit_vertex hands MyFCN a
        // RobustParams block (see MyFCN.h for the algorithm): per-leg
        // DYNAMIC outer annulus (rin past the re-seat radius, rout scaling
        // with leg length -- a long track earns a longer lever), folded
        // disagreement gate, anisotropy + shower vetoes, substitution of
        // the leg's PCA/center only, and a relaxed prior iff a substituted
        // vertex has exactly 2 fittable legs (distorted 2-track vertices
        // need more corrective authority; >= 3 diverse tracks are already
        // precise and keep the 0.43 cm polish prior).  m_mvfit_main_only
        // restricts to the main (neutrino) vertex.  Offline prototype
        // (mvfit_proto.py, docs/pr/51_mvfit_census.tsv): 4/78 fittable main
        // vertices fire at these defaults, unfired vertices numerically
        // untouched; calibration case 18255-234638 (round-5 arm) lands
        // 0.32 cm from the round-6-confirmed charge tip.
        // C++ default false => AddSegment epilogue never runs, FitVertex
        // prior untouched: legacy byte-identical.
        bool m_mvfit_robust{false};
        bool m_mvfit_main_only{true};
        double m_mvfit_min_len{10 * units::cm};
        double m_mvfit_rin_margin{2 * units::cm};
        double m_mvfit_rout_frac{0.5};
        double m_mvfit_rout_min{9 * units::cm};
        double m_mvfit_rout_max{18 * units::cm};
        double m_mvfit_angle{20};       // deg, folded inner-vs-outer
        int m_mvfit_min_pts{5};
        double m_mvfit_min_aniso{3.0};
        double m_mvfit_prior_range{1.0 * units::cm};

        // doc sbnd_xin/docs/pr/40 F1 (= doc pr/7 sec 5 / pr/31 P14-F8).
        // Travels via TrackPidOptions::track_pid_persist_dqdx -- see its
        // comment in PRSegmentFunctions.h.  C++ default false = legacy =
        // byte-identical.
        bool   m_track_pid_persist_dqdx{false};

        // doc sbnd_xin/docs/pr/40 round 2 F4.  Travels via
        // TrackPidOptions::track_pid_persist_4mom -- see its comment in
        // PRSegmentFunctions.h.  C++ default false = legacy = byte-identical.
        bool   m_track_pid_persist_4mom{false};

        // doc sbnd_xin/docs/pr/40 F2: spare a segment from an unconditional
        // track-to-electron reclassification when its own median dQ/dx is
        // decisively proton- or muon-like (segment_dqdx_spares_electron_
        // reclass, same >1.75x / <1.2x MIP thresholds the short-track
        // fallback already trusts).  Applied at the three wholesale-
        // conversion sites Part 0 measured as this investigation's actual
        // writers when persistence (F1) either did not apply or was not
        // enough on its own: examine_all_showers' flag_change_showers loop,
        // both reclassify loops in improve_maps_shower_in_track_out (out_
        // tracks and no-direction segments), and improve_maps_no_dir_tracks'
        // Case E (muon-topology demotion).  This is a designed divergence
        // from the prototype, not a port-fidelity fix -- these sites'
        // wholesale conversion is deliberate prototype behaviour (doc pr/9
        // sec 4) -- see porting_dictionary.md.  C++ default false = legacy =
        // byte-identical.
        bool   m_shower_reclass_dqdx_guard{false};

        // doc sbnd_xin/docs/pr/40 F3.  Travels via segment_is_shower_
        // topology's dqdx_guard parameter -- see its comment in
        // PRSegmentFunctions.h.  Same designed-divergence status as F2.
        // C++ default false = legacy = byte-identical.
        bool   m_shower_topo_dqdx_guard{false};

        // doc sbnd_xin/docs/pr/40 round 2 F5: an electron cannot father a proton.
        // Unlike F2/F3 (which guard the stage-3 wholesale conversion sites in
        // NeutrinoTrackShowerSep.cxx), this fires inside
        // set_default_shower_particle_info -- the SINGLE stage-4 choke point
        // where any flag_shower segment still missing particle_info gets
        // defaulted to electron (NeutrinoPatternBase.cxx, called from
        // examine_direction).  The stage-3 sites have no main_vertex yet
        // (determine_main_vertex has not run), so the neutrino-vertex
        // topology test below cannot be evaluated there; this choke point is
        // also where SBND evt 256587 seg 11079 -- the actual owner-reported
        // case this closes -- was traced to (doc pr/40's Part-0 "topology
        // path; PID never ran").  Relabel to PION (211) instead of electron
        // when the candidate segment (a) emanates from the neutrino vertex
        // (graph identity, not a distance cut -- every measured case sits at
        // d=0.00 cm) and (b) its FAR end is a vertex whose out-edges include
        // a segment already PID'd proton (2212) AND independently
        // charge-confirmed (median dQ/dx > 1.75x MIP by its own points).
        // Measured (doc pr/40 round 2): 5/2209 electron-labelled segments in the
        // 48-event nueCC48 population satisfy both; a naive "any electron
        // segment with a >1.75x MIP neighbour" rule (no vertex requirement,
        // no charge confirmation) fires 348/2209 -- the neutrino-vertex
        // requirement is what excludes the ordinary, CORRECT nueCC topology
        // where an electron and a proton merely share the neutrino vertex
        // as siblings rather than parent/daughter.  This is a designed
        // divergence, not a port-fidelity fix -- the prototype has no
        // proton-daughter veto anywhere -- see porting_dictionary.md.  C++
        // default false = legacy = byte-identical.
        //
        // Round 2 found this override could be silently reverted by
        // Shower::update_particle_type (PRShower.cxx, called from 8 sites in
        // NeutrinoShowerClustering.cxx), which reasserts electron on a
        // shower's start segment whenever shower_length > track_length with
        // no PID/topology awareness (evt 256587 seg 11079: traced,
        // pdg 211 -> 11 at PRShower.cxx:801; only 1/2209 survived
        // end-to-end).  Round 3 (doc pr/40 round 3) closed this by threading
        // a matching guard into update_particle_type itself (see
        // PRShower.h's docstring on update_particle_type) -- gate-clean
        // (48/48 byte-identical off, evt 256587 now 211 end-to-end,
        // population census 2/2209, zero nusel verdict regression).  SBND
        // PRODUCTION DEFAULT ON (cfg-only flip, doc pr/40 round 3).  C++
        // default here stays false = legacy = byte-identical; the SBND
        // operating point sets it true in wct-pr-perevt.jsonnet.
        bool   m_shower_proton_daughter_pion{false};

        // doc sbnd_xin/docs/pr/40 round 4 F7 -- a pion is a track, not a
        // shower.  m_shower_proton_daughter_pion (F5, above) relabels a
        // shower-flagged segment's PDG to 211 but leaves kShowerTrajectory/
        // kShowerTopology set, so shower_clustering_with_nv still roots a
        // Shower there (NeutrinoShowerClustering.cxx: is_shower_seg tests
        // the flags, not the final pdg).  Two visible consequences on SBND
        // evt 256587 seg 11079: its charge-confirmed proton daughter
        // (seg 11080) never gets its own particle-flow node -- it is
        // pre-claimed into the shower's segment set
        // (MultiAlgBlobClustering.cxx fill_bee_pf_tree, used_segs =
        // shower_segs) -- and the pi+ Bee node's displayed endpoint is the
        // SHOWER's end point (a 0.35 cm fragment absorbed from a different,
        // non-main cluster), not segment 11079's own end.  When this knob
        // and m_shower_proton_daughter_pion are both on, the override also
        // clears the shower flags so the segment is treated as an ordinary
        // track from here on.  Designed divergence (the prototype never
        // reclassifies a shower as a track after F5-style relabeling); see
        // porting_dictionary.md.  Population (48-event nueCC48, on
        // work-pr40r3-on48): 2/2209 electron-labelled segments carry pdg 211
        // + flag_shower before this fix.  C++ default false = legacy =
        // byte-identical; requires m_shower_proton_daughter_pion also true
        // to have any effect.
        bool   m_shower_proton_daughter_pion_dissolve{false};

        // doc sbnd_xin/docs/pr/40 round 4 F8 -- a muon cannot terminate in a
        // multi-proton hadronic vertex.  SBND evt 489330: a mu- segment's
        // FAR end (not the neutrino vertex) has TWO charge-confirmed proton
        // daughters (same PID + charge-confirmation test as F5's
        // segment_has_proton_daughter, generalized to segment_at_multi_
        // proton_vertex with min_protons=2).  On a fire, relabel that one
        // segment to pion (211); no propagation across a degree-2 kink to
        // the parent muon segment further from the vertex (owner decision --
        // the sibling segment stays mu-).  A genuine numuCC muon-plus-two-
        // protons AT THE NEUTRINO VERTEX is the ordinary correct topology
        // and must not fire -- same main-vertex exclusion idiom as F5.
        // Designed divergence (the prototype has no proton-multiplicity
        // veto); see porting_dictionary.md.  Population (48-event nueCC48):
        // 1/N muon segments fires this rule (evt 489330 seg 4019); 6 more
        // muon segments have exactly one qualifying proton at a non-main
        // vertex and are deliberately NOT touched (owner's "two protons"
        // wording is read literally).  C++ default false = legacy =
        // byte-identical.
        bool   m_muon_multi_proton_pion{false};

        // doc sbnd_xin/docs/pr/40 round 5 F9 -- narrows F1
        // (m_track_pid_persist_dqdx).  Travels via
        // TrackPidOptions::track_pid_persist_dqdx_electron_guard -- see its
        // comment in PRSegmentFunctions.h.  SBND evt 84229 seg 19038: F1
        // persists a WEAK, UNDIRECTED (dirsign==0) electron guess on a
        // neighboring 4.9 cm segment that would otherwise have been
        // discarded (free_end_dir false, no confident direction at all).
        // That persisted pdg==11 is then read by a DIFFERENT segment's
        // flag_shower_in test (NeutrinoVertexFinder.cxx:1320: has_particle_
        // info() && |pdg|==11), which treats it as an established shower
        // neighbor and unconditionally demotes an independently-confident,
        // free-ended muon call (pdg 13, score 0.04, dirsign=-1) to electron
        // at :1370-1376.  F1's own design intent (persisting a computed
        // proton/muon identity that failed only the free-end topology
        // check) never needed an UNDIRECTED electron guess to survive -- of
        // F1's originally-measured rescue population (doc pr/40 Part 0:
        // 74544, 267597, 269774, 174637, 423981, 433451), none is pdg 11.
        // When this knob is on, F1 no longer rescues a pdg==11 conclusion
        // that lacks a free end (dirsign==0 or non-free topology); it still
        // rescues every non-electron case exactly as before.  Designed
        // narrowing of an already-shipped default-ON knob, not a port-
        // fidelity fix (prototype_base/pid/ProtoSegment.cxx:1637-1639 has
        // no electron-specific carve-out either, since WCPPID's electron
        // template competition works differently) -- see
        // porting_dictionary.md.  C++ default false = today's shipped F1
        // behaviour = byte-identical.
        // MEASURED (doc pr/40 round 5): fixes seg 19038's own pdg (13,
        // correct) but the Bee/mc.json node is unchanged -- 19038 is still
        // flood-filled into neighbor 19039's shower by
        // Shower::complete_structure_with_start_segment (PRShower.cxx:
        // 337-408), which has no per-segment test.  NOT flipped; G2 open.
        bool   m_track_pid_persist_dqdx_electron_guard{false};

        // doc sbnd_xin/docs/pr/40 round 5 F10.  Travels via
        // PatternAlgorithms::shower_clustering_connecting_to_main_vertex's
        // straight_guard parameter -- see its comment in
        // NeutrinoShowerClustering.cxx.  SBND evt 54341 seg 18005: this
        // function picks the single largest "EM-shower-like" candidate
        // connected to the main vertex using only geometric/topological
        // criteria (vertex multiplicity, total length, flag_good_track) --
        // none of which inspect the candidate's own dQ/dx or straightness --
        // and force-sets its start segment to pdg 11.  A 21.3 cm, 0.99-
        // straight stopping-muon stem satisfies every existing criterion.
        // When this knob is on, a candidate whose start segment is long and
        // straight (segment_is_straight_long_track, same shape as the
        // toolkit's own straightness demotion at NeutrinoVertexFinder.cxx:
        // 1432-1447, itself a byte-faithful port of prototype
        // NeutrinoID_track_shower.h:2042-2054) is excluded before the
        // max-length selection, the same way any other disqualified
        // candidate is skipped.  Designed divergence -- the prototype's
        // analogous main-vertex EM-shower selection
        // (NeutrinoID_shower_clustering.h) has no straightness veto either
        // -- see porting_dictionary.md.  C++ default false = legacy =
        // byte-identical.
        // MEASURED (doc pr/40 round 5): the split shape is achieved (seg
        // 18005 excludes from the shower), but once unshielded, ORDINARY
        // track PID -- not this guard -- calls it proton (2212) from its own
        // elevated dQ/dx, not the intended mu-.  Open physics question for
        // the owner, not a defect in this guard.  NOT flipped; G2 open.
        bool   m_shower_connect_main_vertex_straight_guard{false};

        // doc sbnd_xin/docs/pr/40 round 5 F11.  Travels via
        // segment_is_shower_trajectory's straight_guard parameter -- see
        // its comment in PRSegmentFunctions.h.  SBND evt 55715 seg 15007:
        // segment_is_shower_trajectory (the fallback shower test run only
        // when segment_is_shower_topology did not already set the flag,
        // NeutrinoTrackShowerSep.cxx:136-143) has no dQ/dx or straightness
        // guard at all -- unlike its sibling segment_is_shower_topology,
        // which F3 (m_shower_topo_dqdx_guard) already guards.  Once
        // kShowerTrajectory is set, determine_direction's trajectory branch
        // unconditionally assigns pdg 11 (segment_determine_shower_
        // direction_trajectory, PRSegmentFunctions.cxx).  This segment's own
        // median dQ/dx (1.70x MIP) sits just under the existing 1.75x
        // proton-like threshold segment_dqdx_spares_electron_reclass uses,
        // so that helper (reused unchanged by F2/F3) does not discriminate
        // here -- straightness does: 0.97 direct/arc ratio over 14.7 cm.
        // Same segment_is_straight_long_track helper and shape as F10.
        // Designed divergence -- the prototype's is_shower_trajectory
        // (ProtoSegment.cxx:543-613) has its own >50 cm hard length veto
        // but no direct/arc straightness veto below that -- see
        // porting_dictionary.md.  C++ default false = legacy =
        // byte-identical.
        // MEASURED (doc pr/40 round 5): fixes seg 15007's own pdg (13,
        // correct) but is a CONFIRMED REGRESSION -- clearing 15007's flag
        // makes shower_clustering_with_nv_in_main_cluster's is_shower_seg
        // (NeutrinoShowerClustering.cxx:116-119, untouched by F10) re-seed
        // the shower one segment further up, flipping seg 15005 pi+(211) ->
        // e-(11) against the owner's explicit "leave 15005 alone" decision.
        // Isolated to this knob alone via a clean single-knob A/B. NOT
        // flipped; G2 open.
        bool   m_shower_traj_straight_guard{false};

        // doc sbnd_xin/docs/pr/40 round 6 F12 -- the boundary-level fix
        // round 5's F9/F10/F11 measurement demanded.  Travels via
        // Shower::complete_structure_with_start_segment's absorb_track_guard
        // parameter (all 7 call sites in NeutrinoShowerClustering.cxx).  The
        // flood-fill has NO per-segment test (PRShower.cxx), so one
        // shower-flagged seed swallows every connected segment regardless of
        // its own PID; the absorbed muon's length then also biases
        // Shower::update_particle_type's majority vote toward electron, and
        // the merged Bee/mc.json node shows one "e-" (SBND evt 84229 seg
        // 19038, a confident mu- flood-filled into neighbor 19039's 4.9 cm
        // shower).  When on, a confidently PID'd non-electron (pdg != 0,
        // |pdg| != 11) that is long and straight (segment_is_straight_long_
        // track) is not absorbed; the walk terminates there; the segment
        // then gets its own PF node (view-keyed suppression in
        // fill_bee_pf_tree; the pf_shower_vertex_barrier orphan safety net
        // covers the BFS-unreachable case).  Long-muon pseudo-showers
        // (Shower::get_particle_type()==13) are exempt so broken muon
        // reassembly keeps working.  Designed divergence -- the prototype's
        // complete_structure_with_start_segment has no per-segment test
        // either -- see porting_dictionary.md.  C++ default false = legacy =
        // byte-identical.
        // MEASURED (doc pr/40 round 6): fixes BOTH remaining display cases.
        // 84229: seg 19038 gets its own mu- PF node (owner's split shape)
        // with 19040 as the separate e- child.  55715: 15006/15007 are no
        // longer absorbed into the 15005-seeded candidate, which then fails
        // EM acceptance -- 15007 gets its own mu- node.  Isolated per-case
        // via single-knob A/B arms.
        bool   m_shower_absorb_track_guard{false};

        // doc sbnd_xin/docs/pr/65 round 3 -- offer graph-unreachable
        // main-cluster segments to the proximity/angle shower absorbers.
        // The prototype's seg1->cluster()==main_cluster guards
        // (NeutrinoID_shower_clustering.h:1870/:1901 etc) encode "already
        // claimed by the main_vertex graph walk", an invariant that
        // other_seg_keep_isolated (doc pr/54) broke by adding kept residual
        // segments as DISCONNECTED components of the main cluster's graph.
        // When on, the six absorber guards in NeutrinoShowerClustering.cxx
        // test reachability instead of raw cluster identity, so those
        // fragments can be absorbed by the prototype's own already-tuned
        // 3.5 cm / cone thresholds instead of surfacing as fabricated
        // PF-root nodes.  On a connected main cluster the predicate is
        // identical to legacy.  C++ default false = legacy = byte-identical.
        bool   m_shower_absorb_unreachable_main{false};
        // Transient companion state for the knob above: the main-cluster
        // segments NOT reachable from main_vertex, recomputed unconditionally
        // at every shower_clustering_with_nv entry (empty when the knob is
        // off, so the guards' predicate is byte-identical to legacy).  Valid
        // for the span of that call only -- PR-graph topology is frozen there.
        IndexedSegmentSet m_absorb_unreachable_main_segs;

        // doc sbnd_xin/docs/pr/40 round 6 F13 -- closes round 5's F11
        // regression (SBND evt 55715 seg 15005 pi+ -> e-).  Travels into
        // shower_clustering_connecting_to_main_vertex's candidate skip chain
        // (NeutrinoShowerClustering.cxx).  Round-6 trace (probe tag
        // "connecting_to_main_vertex") showed the regression's writer is
        // THIS function, not the BFS re-seed round 5 hypothesized: once F11
        // declassifies daughter 15007 (correctly demoted to mu- by the
        // straightness demotion), the 6.1 cm pi+ parent 15005 -- UNDER
        // segment_is_straight_long_track's 10 cm floor, so F10's straight
        // branch cannot save it -- is selected as the EM candidate and
        // force-set to pdg 11.  When on, a pdg==211 candidate with a
        // charge-confirmed proton daughter (segment_has_proton_daughter,
        // the round-3 protected_pion predicate: "an electron cannot father
        // a proton") is skipped.  C++ default false = legacy =
        // byte-identical.
        // MEASURED (doc pr/40 round 6): DEAD as shaped, NEVER FLIP.  The
        // full-transition trace (WCT_PID_WRITE_DEBUG=2) showed 15005 is
        // already pdg 2212 by candidate-selection time -- its baseline 211
        // only ever existed downstream of 15007's WRONG e- label (Michel
        // rescue 2212->13 at NeutrinoVertexFinder.cxx:1804, then the
        // single-muon pion demotion 13->211 at :1679), a chain that F11
        // correctly eliminates.  A pdg==211 predicate therefore cannot fire
        // on the motivating case (confirmed: F11+F13 arm still shows the
        // merged e- node), and F12 alone yields the intended display.  Kept
        // as a documented negative result (doc pr/36 F2 precedent), default
        // false, excluded from the round-6 flip.
        bool   m_shower_connect_protected_pion_guard{false};

        // doc sbnd_xin/docs/pr/40 round 6 F14 -- closes round 5's F10
        // residual (SBND evt 54341 seg 18005: split achieved, stem labelled
        // proton).  Travels via override_michel_stem_muon (called next to F8
        // in examine_direction, last word before shower_clustering_with_nv).
        // The toolkit's own Michel rescue ("a stopped proton cannot produce
        // a Michel electron", examine_direction's weak-direction branch)
        // encodes exactly the right physics but only reaches seg_dir_weak
        // segments whose stopping vertex has EXACTLY 2 segments; 18005's
        // Bragg rise wins the proton template with a confident direction
        // (template competition has no absolute quality gate) and its
        // stopping vertex has degree 4.  When on, a pdg==2212 segment that
        // is long and straight, emanates from the main vertex, and whose far
        // (stopping) vertex carries >=1 shower-like sibling (kShowerTrajectory
        // or |pdg|==11 -- the existing Michel sibling test) is relabelled
        // mu- with a recomputed 4-momentum.  C++ default false = legacy =
        // byte-identical.
        // MEASURED (doc pr/40 round 6): 54341 node 18005 reads mu- 74 MeV
        // with children e- 19 MeV (18006) + mu- 11 MeV (18007) -- exactly
        // the owner-requested shape.  Confirmed no effect on the other two
        // cases (isolation arm F11+F13+F14 unchanged vs F11+F13).
        bool   m_michel_stem_muon_rescue{false};

        // doc sbnd_xin/docs/pr/74 round 2 P1 (SBND 18345 evt 53361 seg
        // 27004).  examine_direction's flag_shower_in cascade relabels any
        // downstream |pdg|==13/pdg==0 segment electron with no length
        // ceiling and no charge test (prototype-faithful:
        // NeutrinoID_track_shower.h:2004 has the same unconditional branch).
        // A 113.9 cm segment at 1.02x MIP -- as MIP-like as a track can be
        // -- was relabelled e- this way.  When on, the relabel is refused
        // for a segment that is BOTH long (> m_shower_in_max_len) AND
        // MIP-like (median dQ/dx < m_shower_in_mip_hi x mip_dqdx_median);
        // 90055's genuine shower bodies (0.75-0.90x MIP but 7-20 cm) are
        // spared by the length conjunct.  A zero/absent median never vetoes
        // (no evidence != MIP-like evidence, same convention as
        // segment_dqdx_spares_electron_reclass).  C++ default false =
        // legacy = byte-identical.
        bool   m_shower_in_cascade_guard{false};
        double m_shower_in_max_len{40*units::cm};
        double m_shower_in_mip_hi{1.3};

        // ------------------------------------------------------------------
        // doc sbnd_xin/docs/pr/40 round 9 -- the round-7/round-8 named guard
        // family: five sites that force pdg 11 onto a straight long track
        // with no geometry test, each declined (write skipped, structure
        // untouched) when segment_is_straight_long_track(_or_continuation)
        // fires.  All C++ default false = legacy = byte-identical.
        //
        // Site corrections found implementing (round-9 doc "D" findings):
        // D1 -- round 7's "candidate 2c at NeutrinoShowerClustering.cxx:1401"
        //   is the SAME write as round 8 Part A (line drift since the trace
        //   commits); the start_seg knob below is re-targeted at the
        //   accept-time set_pdg(11) in shower_clustering_connecting_to_main_
        //   vertex, the real same-cluster analogue.
        // D2 -- round 7's "dirsign()==0 branch at NeutrinoVertexFinder.cxx
        //   :1659" is dead code (dirsign is assigned 6 lines above); the
        //   real coverage gap is the pr/74 cascade veto's 40cm/1.3xMIP
        //   thresholds, so the examine_direction knob adds a GEOMETRY arm
        //   beside that veto in all three flag_shower_in branches.
        // ------------------------------------------------------------------

        // Round 8 Part A: shower_clustering_with_nv_from_vertices's pdg
        // block (the ONLY creator of a conn-2 shower from a fresh track,
        // SBND 286906 seg 9002 / 409546 seg 9000).  Uses the continuation-
        // aware predicate because PATH C hands the guard a broken sub-10cm
        // HALF of the track (286906: 8.68cm anchor at 4.9deg kink to the
        // 126.89cm body).  Also co-guards the update_particle_type call for
        // the same shower, which would otherwise redo the declined 13->11
        // via the majority vote (PRShower.cxx:981-1009).
        bool   m_shower_connect_from_vertices_straight_guard{false};

        // D1 re-target: the accept-time set_pdg(11) on the winning
        // main-vertex EM candidate's start segment (F10 vetoes only at seed
        // time on the seed's own geometry; a short anchor collinear with a
        // long straight sibling passes F10 and reaches this write).
        bool   m_shower_connect_start_seg_straight_guard{false};

        // doc sbnd_xin/docs/pr/93 round 3 -- four knobs for the five owner
        // events where an "electron" is really tracks or a hadronic-pi0
        // shower (SBND 18255-55595/348471/69314/292643/315167).  All C++
        // default false = legacy = byte-identical.
        //
        // Cause A (55595): improve_maps_no_dir_tracks Case B writes pdg 11
        // unconditionally when dirsign()==0, regardless of length (193.8cm
        // MIP muon here).  Same F2 dQ/dx spare-test its sibling Case E
        // already carries.
        bool   m_shower_reclass_case_b_dqdx_guard{false};
        // Cause B (348471, 69314): the new_shower_accepted and
        // merged_shower_start_segment sites force set_pdg(11) on a shower's
        // start segment with NO PID check -- decline the write when the
        // segment already carries a confident non-electron template PID
        // (segment_confident_nonelectron_pid).
        bool   m_shower_accept_pid_guard{false};
        // Shared minimum-length floor (internal units; 50cm default = the
        // scale of SBND's shower_topo_demote_len "a >50cm segment is not
        // EM-flaggable" rule) for the Cause A, Cause B, AND Cause D
        // declines of this family.
        // Below it, a confident non-electron template score is NOT reliable
        // evidence against electron: real nueCC48 electron stems of 22-47cm
        // carry 0.11-0.64 proton/muon scores (the template competition
        // never considers electron), and the un-floored guards regressed
        // 9+17 of 48 nueCC48 events on the attribution arms.  Inert while
        // both guards are off.
        double m_shower_pid_guard_min_len{50*units::cm};
        // Cause C (292643): update_particle_type's vote counts only
        // confirmed protons as track, so a muon/pion chain always votes
        // electron -- when on, any unflagged member with |pdg| in
        // {13,211,2212} counts as track (threaded as a trailing param to
        // Shower::update_particle_type).
        bool   m_shower_vote_track_pid_counts{false};
        // Cause D (315167): pass 3's direction-cone absorber
        // (shower_clustering_with_nv_from_main_cluster's angle/distance
        // sweep) has no PID or straightness check, and with
        // shower_absorb_unreachable_main ON (pr/65, SBND production) a
        // graph-unreachable MAIN-cluster segment is eligible there --
        // 315167's 150.7cm score-0.10 proton was cone-absorbed into a
        // 15.7cm EM stub's shower, whose energy is then computed from
        // total_length under the electron hypothesis (1046.7 MeV).  When
        // on, that absorber declines a confidently-PID'd non-electron
        // straight-long track (the flood-fill guard_excludes predicate,
        // pr/40 F12); the declined segment stays unclaimed.  Confirmed by
        // the WCT_SHOWER_ABSORB_DEBUG tape (site=pass3_cone seg=8001).
        bool   m_shower_cone_absorb_guard{false};

        // ------------------------------------------------------------------
        // doc sbnd_xin/docs/pr/93 round 4 -- PF-hierarchy fine-tunes.
        // ------------------------------------------------------------------
        // shower_detach_track_stem (348471 + 292643): a TRACK-HEADED shower
        // -- final start segment is a non-shower-flagged, non-zero,
        // non-+-11-pdg track -- is really "track + daughter EM shower" glued
        // into one object (348471: 53.5cm score-0.232 proton stem + 15 EM
        // satellites; 292643: 22.9cm pi+ stem [score-100 sentinel] + 13.9cm
        // mu continuation + 9 EM satellites).  When on, a post-pass after
        // shower dedup peels the main-cluster track prefix back into the
        // track pool and re-roots the shower at the prefix's far vertex
        // (conn 2), so PF renders "track -> pseudo-gamma -> EM shower", kine
        // costs the track by range (not shower charge + rest mass -- 348471
        // gained 938 MeV of proton mass on a charge-based 719 MeV aggregate
        // in round 3), and the pi0 pairing sees a real EM shower.
        // STRUCTURE-keyed on purpose (not the round-3 confident-score/50cm
        // family: 21 of the 50 non-trivial track-headed stems in the census
        // carry the score-100 sentinel, and 292643's stem is under the
        // floor).  Long-muon pseudo-showers (cached type +-13) are exempt;
        // pure track chains and single-blob protons self-exclude (the
        // detach refuses to empty a shower).
        bool   m_shower_detach_track_stem{false};
        // ------------------------------------------------------------------
        // doc sbnd_xin/docs/pr/99 round 2 -- shower_ghost_member_drop.
        // 395148: the 295 MeV "electron" held a 23.4 cm projective-ghost
        // member (median per-point charge -678, 93% of points dQ<=0, 2D
        // view overlap vs a sibling member 1.00/1.00/0.25 -- the pr/83
        // op1-proj class) 95 cm from the main vertex, outside every mvga
        // scope.  When on, a post-pass after shower dedup (and BEFORE the
        // pi0 finders, so pairing never sees the ghost) scans member PAIRS
        // of each shower with the op1-proj discriminator restated for
        // membership: same (apa,face), 2nd-best per-view 2D overlap >=
        // ghost_overlap_frac, candidate starved (median dQ/dx ratio <=
        // ghost_dqdx_ratio OR frac(dQ<=0) > 0.5) while the partner is
        // healthy (ratio >= 2x ghost_dqdx_ratio), candidate length >=
        // ghost_min_len.  The conjunction is mandatory -- the pr/99 census
        // measured 51/185 real multi-member showers with SOME charge-
        // starved fragment; overlap-with-a-healthy-sibling is what makes it
        // a ghost.  The ghost is dropped from the shower view (leaf-only,
        // Shower::drop_ghost_member contract) and DELETED from the graph so
        // it leaves the PF/Bee display; vote/kinematics/charge recomputed.
        // No refit at this seat -- acceptable for a charge-starved ghost.
        // false (default) => no pass => byte-identical; the numeric knobs
        // are inert while the bool is off.
        bool   m_shower_ghost_member_drop{false};   ///< doc pr/99 round 2; false = off
        double m_shower_ghost_overlap_frac{0.7};    ///< 2nd-best per-view overlap gate; inert while drop off
        double m_shower_ghost_dqdx_ratio{0.25};     ///< starved gate, ratio vs m_mip_dqdx_median; inert while drop off
        double m_shower_ghost_min_len{10*units::cm}; ///< candidate min length, internal units; inert while drop off
        // doc pr/99 round 3 (A5 hadronic-shower tag; owner: "many hadronic
        // showers / particle flow still labeled as electrons").  A claimed
        // EM shower (|particle_type|==11, conn type 1) whose surrounding
        // charge population SHRINKS downstream is a hadron, not an electron:
        // real EM cascades grow 3-6x in in-cylinder point population over
        // the first 30 cm (measured round 1: 168596/14153 5.7x, 360535/7060
        // 3.1x, adverse gamma control 506114 2.3x) while the pr/99 fakes
        // read 0.2x (315167, terminal Bragg rise 2930->10801) and 0.5x
        // (395148).  Verdict = smax >= min_len && (growth < growth_max ||
        // (bragg >= bragg_ratio && growth < growth_bragg)).  On verdict the
        // START SEGMENT's pdg is stamped 211 (pi+-) with mass + 4-momentum
        // refresh (the pi0 incoming-track stamp recipe) plus
        // shower->set_particle_type(211) -- the durable route, since
        // Shower::calculate_kinematics re-copies the start segment's pdg on
        // every recompute.  NOT 13: long-muon consumers route on |13| and
        // would misroute a shower absent from segments_in_long_muon.  The
        // re-typed shower is excluded from pi0 pairing via
        // m_hadronic_retyped_shower_ids.  A DEBUG census line is emitted for
        // EVERY evaluated shower (the offline calibration channel).  false
        // (default) => no pass => byte-identical; numerics inert while off.
        bool   m_shower_hadronic_tag{false};              ///< doc pr/99 round 3; false = off
        double m_shower_hadronic_min_len{10*units::cm};   ///< min trajectory extent to judge; inert while tag off
        double m_shower_hadronic_scan_len{30*units::cm};  ///< growth window along the trajectory; inert while tag off
        double m_shower_hadronic_bin{3*units::cm};        ///< arc-length bin; inert while tag off
        double m_shower_hadronic_r_cyl{8*units::cm};      ///< in-cylinder population radius; inert while tag off
        double m_shower_hadronic_r_core{1.2*units::cm};   ///< core radius, off-axis census only; inert while tag off
        double m_shower_hadronic_growth_max{0.8};         ///< growth below => hadronic; inert while tag off
        double m_shower_hadronic_growth_bragg{1.2};       ///< growth ceiling for Bragg-confirmed branch; inert while tag off
        double m_shower_hadronic_bragg_ratio{3.0};        ///< terminal/trunk median dQ/dx rise; inert while tag off
        double m_shower_hadronic_stem_ratio{0};           ///< proton-stem branch: stem median (MIP units) at/above this + growth < growth_bragg => hadronic; 0 = branch off
        /// Shower ids (stable per-run Shower::get_shower_id) re-typed by the
        /// A5 pass this event; guards the pi0 finders' candidate collection.
        /// Cleared at shower_clustering_with_nv pass entry.  Empty (tag off)
        /// => the guards are no-ops => byte-identical.
        std::set<int> m_hadronic_retyped_shower_ids;
        // kine_count_orphan_tracks (315167): fill_kine_tree counterpart of
        // the PF-side pf_orphan_confident_track knob (BeePFConfig).  A
        // confident straight-long main-cluster track that is graph-
        // disconnected from the main vertex (freed from shower membership by
        // shower_cone_absorb_guard) is reached by NEITHER the kine BFS nor
        // any shower: 315167's 150.7cm score-0.101 proton (~595 MeV KE)
        // silently vanished from kine_reco_Enu.  When on, a post-pass
        // pushes such segments via push_segment_kine (range KE + binding).
        // Shares segment_orphan_confident_track with the PF side so the two
        // outputs describe the same particle set.
        bool   m_kine_count_orphan_tracks{false};
        double m_kine_orphan_track_min{50*units::cm};

        // D2: geometry arm beside the pr/74 P1 cascade veto in examine_
        // direction's flag_shower_in branches (54629 seg 15007: 31cm,
        // 1.42xMIP -- fails BOTH of the P1 veto's conjuncts, but 0.97
        // direct/arc straightness catches it).
        bool   m_examine_direction_dirsign_shower_in_guard{false};

        // Round 7 candidate 2b: the num_daughter_showers>=4 / angle-
        // collinearity wholesale reclass write (NeutrinoVertexFinder.cxx,
        // 54629 seg 15011: 94.6cm, 0.980 straight).  Guards the write only,
        // NOT the outer condition (falsifying that would re-route control
        // into the pdg==11 else-if chain).
        bool   m_daughter_shower_angle_reclass_straight_guard{false};

        // Round 7 candidate 1 (320865 seg 13001: 48.1cm, 0.94 straight):
        // third arm in improve_vertex's topology re-exam decline condition
        // -- re-demote (unset kShowerTopology, keep the track PID) instead
        // of the set_flags+pdg-11 escape.  Framed as a defensible safety
        // net, NOT "the fix for 320865" (that segment only exists via a
        // pr/90 side effect; a pr/90-side fix remains the open alternative).
        // F11-displacement risk is highest here -- see the round-9 census.
        bool   m_shower_topo_reexam_straight_guard{false};

        // Max kink (DEGREES) for the collinear-continuation arm of
        // segment_is_straight_long_track_or_continuation (286906 measures
        // 4.9 deg; 25 matches the daughter-shower reclass's own >155-deg
        // opening-angle idiom, i.e. kink < 25).
        double m_sfv_kink_max{25.0};
        // Transient record of the anchors whose electron write the
        // from_vertices guard declined (same per-call lifecycle as
        // m_absorb_unreachable_main_segs below; empty when the knob is
        // off).  Kept as the marker for any later pass that would re-flip
        // the anchor -- and as the hook the B2 bridge consumes.
        IndexedSegmentSet m_sfv_declined_anchors;

        // ------------------------------------------------------------------
        // doc sbnd_xin/docs/pr/40 round 9 B2 -- cross-cluster track bridge.
        // When shower_clustering_with_nv_from_vertices's directional match
        // accepts a candidate in ANOTHER cluster, the candidate (or its
        // collinear continuation) is a straight long track, and the exact
        // steiner-cloud closest approach between the two clusters is below
        // m_shower_nv_bridge_max_gap: do NOT fabricate the conn-2 electron;
        // add a straight 2-point zero-charge bridge segment vertex->track
        // to the PR graph (stamped main_cluster) so PF/kinematics see one
        // continuous track (owner directive, 286906: gap 1.39cm bridges;
        // 521075: 2.92cm must not).  The gap is a signal-processing hole
        // with no charge -- no do_rough_path routing, no refit.  All
        // bookkeeping below is empty when the knob is off => byte-identical.
        // ------------------------------------------------------------------
        bool   m_shower_nv_bridge_track{false};
        double m_shower_nv_bridge_max_gap{1.8*units::cm};
        // ------------------------------------------------------------------
        // doc sbnd_xin/docs/pr/97 D1 -- shower_clustering_with_nv_from_vertices
        // builds `cluster_point_info main_pi` and sets ONLY .cluster and
        // .min_vertex; .min_angle/.min_dis/.min_point stay INDETERMINATE and
        // are filled in the vertex loop only if main_vertex is itself one of
        // main_cluster_vertices.  When the overall main vertex lives in a
        // different cluster (reachable, and the state of the [B] doctest
        // fixture) nothing fills them, and the "prefer main_pi over min_pi"
        // comparison a few lines later reads stale stack bytes -- a leftover
        // pointer, so the branch is decided by the address-space layout
        // (ASLR, or just the size of the environment).  The prototype has the
        // same hole (NeutrinoID_shower_clustering.h:1071-1073), so this is an
        // inherited defect, not a porting slip.
        // ON: sentinel-initialise so the "main vertex was never evaluated
        // against this cluster" case deterministically prefers min_pi.
        // OFF (default) = legacy indeterminate read, bit-for-bit.
        // ------------------------------------------------------------------
        bool   m_shower_nv_main_pi_init{false};
        // Transient per-call state (same lifecycle as
        // m_absorb_unreachable_main_segs above): the bridge segment plus the
        // bridged cluster's own segments, shielded from every shower
        // flood-fill/absorber so the rescued track is not re-swallowed
        // (the pr/54 lesson: cross-boundary graph edits break the
        // cluster()==main_cluster absorber invariant unless shielded).
        IndexedSegmentSet m_nv_bridge_shield_segs;
        std::set<int>     m_nv_bridge_cluster_ids;
        // Pre-seed a shower flood-fill's used_segments with the shield set
        // so complete_structure_with_start_segment never traverses the
        // bridge or absorbs the rescued cluster (a no-op when the set is
        // empty, i.e. whenever the bridge knob is off or never fired).
        IndexedSegmentSet& nv_bridge_seed(IndexedSegmentSet& used) {
            used.insert(m_nv_bridge_shield_segs.begin(), m_nv_bridge_shield_segs.end());
            return used;
        }

        // B2 helpers (NeutrinoPatternBase.cxx / NeutrinoShowerClustering.cxx)
        double cluster_steiner_gap(const Facade::Cluster& a, const Facade::Cluster& b) const;
        bool   nv_bridge_track(Graph& graph, Facade::Cluster* main_cluster,
                               Facade::Cluster* cluster, SegmentPtr sg1, VertexPtr vertex,
                               const WireCell::Point& point,
                               const std::vector<SegmentPtr>& cluster_segs,
                               TrackFitting& track_fitter, IDetectorVolumes::pointer dv,
                               const Clus::ParticleDataSet::pointer& particle_data,
                               const IRecombinationModel::pointer& recomb_model);
        // doc pr/93 round 4: the connect/register tail of nv_bridge_track,
        // factored out (byte-identical) so the sccc bridge replay below can
        // reuse it.  Returns the bridge segment or nullptr.
        SegmentPtr nv_bridge_connect(Graph& graph, Facade::Cluster* main_cluster,
                                     Facade::Cluster* cluster, VertexPtr vertex,
                                     VertexPtr far_vtx,
                                     const std::vector<SegmentPtr>& cluster_segs,
                                     TrackFitting& track_fitter,
                                     IDetectorVolumes::pointer dv);

        // ------------------------------------------------------------------
        // doc sbnd_xin/docs/pr/93 round 4 (straight_cont_cross_cluster,
        // SBND 18264-137238) -- a main-vertex kShowerTrajectory stem that is
        // the CROSS-CLUSTER continuation of a straight long track (a pr/57
        // W-plane-gap over-clustering split: 14.6cm flagged stem, degree-1
        // tip, 81cm straight muon body 3-4cm away in another cluster) seeds
        // a fake electron over the whole chain.  demote_cross_cluster_
        // straight_stems (called between examine_direction and
        // shower_clustering_with_nv) clears the flag + re-PIDs the stem as a
        // track when the cross-cluster arm of segment_is_straight_long_
        // track_or_continuation matches under the two-tier geometric gate:
        // base (max_gap, kink_max) OR aligned (gap_aligned, kink_tight) --
        // the owner-requested "very aligned buys a larger gap" tier;
        // pr/57:1336 records a 53-deg LOCAL kink at this break class, so
        // the tight-angle long-lever arm is the robust one.  The pr/57
        // owner-good W-gap pairs (21073 6-7, 122660 14-17, 61579 3-4) are
        // mandatory negative controls for any retune.
        // sccc_bridge_body (second rung): additionally RECORD a bridge
        // request {stem tip, continuation cluster, continuation vertex};
        // shower_clustering_with_nv replays it through nv_bridge_connect
        // AFTER its entry clears (building it earlier would be wiped), so
        // the muon body joins the PF/kine track chain through a real graph
        // edge.  All state empty when the master bool is off =>
        // byte-identical.  Gap params internal units; kinks DEGREES.
        // ------------------------------------------------------------------
        bool   m_straight_cont_cross_cluster{false};
        bool   m_sccc_bridge_body{false};
        double m_sccc_max_gap{5*units::cm};
        double m_sccc_kink_max{15.0};
        double m_sccc_gap_aligned{12*units::cm};
        double m_sccc_kink_tight{7.5};
        struct SCCCBridgeRequest {
            VertexPtr main_vtx;          // the demoted stem's degree-1 tip
            Facade::Cluster* cluster;    // the continuation's cluster
            VertexPtr far_vtx;           // the continuation's near endpoint vertex
        };
        // Insertion order = sorted_out_edges(main_vertex) order; a plain
        // vector, never iterated as a pointer-keyed container.
        std::vector<SCCCBridgeRequest> m_sccc_bridge_requests;
        // Demoted stems (per-event lifecycle: cleared at the demotion pass's
        // entry).  Consulted by examine_showers' retarget selection: a stem
        // this knob demoted to a muon head must not be re-adopted as a
        // shower start segment -- the retarget's merged_shower_start_segment
        // stamp writes pdg 11 with the score-100 sentinel, which
        // shower_accept_pid_guard cannot decline (measured on 18264-137238:
        // the demoted stem was re-captured as "e- 71 MeV").  Membership-
        // tested only; empty when the knob is off => byte-identical.
        IndexedSegmentSet m_sccc_shield_segs;
        // Clusters bridged by the sccc REPLAY (not by nv_bridge_track's own
        // Step-5 firing -- that population keeps its legacy behavior).
        // Consulted by from_vertices' Step-3 candidate analysis: a cluster
        // the replay already claimed as a muon-body chain must not be
        // re-analyzed as an EM-shower candidate (measured on 18264-137238:
        // Path C broke the bridged cluster's 3.2cm entry stub off as a
        // conn-2 electron that examine_shower_1 then spliced back onto the
        // demoted stem).  Cleared at shower_clustering_with_nv entry with
        // the rest of the bridge state; empty when the knob is off =>
        // byte-identical.
        std::set<int> m_sccc_bridged_cluster_ids;
        void demote_cross_cluster_straight_stems(Graph& graph, VertexPtr main_vertex,
                                                 const Clus::ParticleDataSet::pointer& particle_data,
                                                 const IRecombinationModel::pointer& recomb_model);

        // ------------------------------------------------------------------
        // doc sbnd_xin/docs/pr/92 -- drop stray satellite showers from the
        // kinematics tree (and, via TrackFitting transport, the Bee PF
        // tree).  fill_kine_tree's leftover pass admits every conn-2/3
        // shower with no direction or distance check (the prototype has the
        // identical hole, NeutrinoID_kine.h:209-255): overclustered cosmics
        // (SBND 350935 shower 11001, 449 MeV at 93 deg off its attachment
        // vertex; 321371 shower 18004, 98 MeV collinear tail of a 256 cm
        // dropped cosmic) and second neutrinos (389538 shower 19040,
        // 997 MeV attached 144.5 cm away, 69 deg off the main vertex) are
        // summed into kine_reco_Enu.  Candidates: BFS-unreached, conn 2/3,
        // start segment in a NON-main cluster, kine_best above the floor,
        // not pi0-paired, not within the proximity exemption of a
        // main-cluster attachment.  Drop arms (axis = fresh
        // shower_cal_dir_3vector from the start point -- the STORED
        // init_dir for conn 2/3 is exactly the vertex->start chord and
        // would always read 0 deg):
        //   A: angle(axis, start - attach_vtx) > m_kine_sat_angle_bad
        //   B: (attach farther than m_kine_sat_far_dis OR attach vertex not
        //      in the main cluster) AND angle(axis, start - main_vtx) >=
        //      m_kine_sat_angle_main
        //   C: shower_start_is_track_continuation (collinear straight-long
        //      sibling OUTSIDE the shower; see PRShowerFunctions.h)
        // All state empty / arms unreachable when the master bool is off =>
        // byte-identical.  Angles in DEGREES (sfv_kink_max precedent).
        // ------------------------------------------------------------------
        bool   m_kine_drop_stray_satellites{false};
        double m_kine_sat_min_energy{20*units::MeV};
        double m_kine_sat_prox_max{8*units::cm};
        double m_kine_sat_angle_bad{60.0};
        double m_kine_sat_angle_main{45.0};
        double m_kine_sat_far_dis{90*units::cm};
        double m_kine_sat_axis_dis_cut{30*units::cm};
        double m_kine_sat_cont_kink{25.0};
        // pr/92 round 2 (owner retune, 2026-08-18): the direction arms
        // A/B/C apply only to TRACK-like satellites (straight-long start
        // segment with <= max_nseg segments, or an out-of-shower track
        // continuation) -- for those, direction inconsistency means
        // overclustering.  EM-shower-like satellites are usually genuinely
        // detached (NCpi0-like) and are dropped only when FAR from the main
        // vertex (> em_far_dis; second neutrinos sit 169-250 cm, legit
        // detached fragments 18-119 cm on the survey samples) AND the
        // folded (sign-insensitive) main-vertex angle fails.
        int    m_kine_sat_track_max_nseg{3};
        double m_kine_sat_em_far_dis{150*units::cm};

        // doc sbnd_xin/docs/pr/74 round 2 P2 (SBND 18255 evt 90055 seg
        // 11045).  override_michel_stem_muon (F14 above) accepts ANY
        // shower-like sibling at the stem's far vertex as "the Michel
        // electron"; on a nueCC event the sibling is the 2020 MeV EM
        // shower's start segment and the rescue paints a muon at the
        // neutrino vertex.  A genuine Michel electron is TERMINAL: the
        // graph beyond the far vertex is a few cm of electron and nothing
        // else.  When on, the rescue additionally requires the total track
        // length reachable beyond the far vertex (excluding the stem) to be
        // below m_michel_stem_max_far_len; a shower trunk heading a
        // 155 cm+ downstream tree fails.  C++ default false = legacy =
        // byte-identical.
        bool   m_michel_stem_michel_check{false};
        double m_michel_stem_max_far_len{40*units::cm};

        // doc sbnd_xin/docs/pr/74 round 2 K4 (SBND 18255 evts 90055 +
        // 469665).  Shower formation walks outward from the main vertex and
        // starts a shower at the first shower-like segment; the track-typed
        // stem it walked PAST is structurally excluded, and no later step
        // pulls a vertex-attached stem into the shower beyond it (the
        // prototype has the same gap -- NeutrinoID_shower_clustering.h:
        // 1654-1706 has no backfill either, so this is a WCT improvement,
        // not a parity fix).  When on, a post-pass walks from each
        // substantial EM shower's attach vertex back toward the main vertex
        // and absorbs the chain while each segment is short
        // (< m_stem_backfill_max_len), not charge-hot (median dQ/dx <
        // m_stem_backfill_mip_hi x MIP median -- a Bragg proton at ~9x MIP
        // stops the walk, 469665's vertex proton survives), not in a long
        // muon, and not already claimed by a shower.  Membership then drives
        // the paint and the PF tree.  C++ default false = legacy =
        // byte-identical.
        bool   m_shower_stem_backfill{false};
        double m_stem_backfill_max_len{30*units::cm};
        double m_stem_backfill_mip_lo{0.75};
        double m_stem_backfill_mip_hi{3.5};
        double m_stem_backfill_min_shower_len{40*units::cm};

        // doc sbnd_xin/docs/pr/74 round 2 K5 = doc pr/65's deferred rung 2
        // (SBND 18255 evt 142421 seg 7013: 41.9 cm, 266 MeV, in a
        // disconnected main-cluster graph component, PF-invisible; the
        // pr/65 rung-3 audit prints it and by design fabricates nothing;
        // rung 1's absorbers are geometry-gated and never reach it).  When
        // on, shower_clustering_in_other_clusters's own leftover-cluster
        // branch -- the prototype's connection_type=3 pseudo-gamma path --
        // is extended to graph-unreachable, unclaimed main-cluster segments
        // above m_conn3_unreachable_min_len: same nearest-candidate-vertex
        // anchor, same <80 cm conn-3/conn-4 split, same
        // complete_structure_with_start_segment component claim.  Honors
        // pr/65 round 2's two hard requirements: no fabricated PF-root
        // particle, and arrival via the existing clustering algorithm.
        // C++ default false = legacy = byte-identical.
        bool   m_shower_conn3_unreachable{false};
        double m_conn3_unreachable_min_len{10*units::cm};

        // doc sbnd_xin/docs/pr/84 round 3.  Two PR::Showers can be built on
        // the SAME start segment, because the start-segment choice in
        // shower_clustering_in_other_clusters consults no claim at all (its
        // guards are cluster-level, and the segment it picks off the
        // cluster's vertex can belong to another cluster entirely -- the
        // prototype has the identical hole, NeutrinoID_shower_clustering.h
        // :1481-1495), and because the K5 conn3_unreachable branch checks
        // map_segment_in_shower, which update_shower_maps only refreshes at
        // the END of the function.  Measured on SBND 169626/174752/347129/
        // 394532: the twin renders a SECOND PF node carrying the same jsTree
        // id (`cluster_id*1000 + seg_id`), which breaks the Bee tree, and it
        // is counted a second time in kine_energy_particle -- 394532's
        // kine_reco_Enu 352.2 MeV against 255.5 de-duplicated.  When on, a
        // merge pass runs after examine_showers and before the pi0 finders:
        // each group of showers sharing a start segment collapses onto its
        // most-directly-connected member, which ABSORBS the others' segments
        // (Shower::add_shower, whose membership gate makes the overlap
        // idempotent), and the survivor's kinematics are recomputed by the
        // production calculate_shower_kinematics (flag_kinematics cleared;
        // every other shower keeps its flag and is untouched).  conn_type==4
        // members never participate: they are skipped by both the PF tree and
        // the kine tree, so folding their charge in would be a change with no
        // reported defect behind it.  false = no pass = byte-identical.
        bool   m_shower_dedup_start_seg{false};

        // doc sbnd_xin/docs/pr/84 round 2 (F3 = pr/84 P2 "conn3_stitch").
        // The main cluster's segment graph can end up disconnected (pr/54
        // keep-isolated residuals, snap-stranded vertices) even though the
        // cluster is one contiguous lump of charge; the unreachable pieces
        // are then promoted to conn-3 "association" showers at anchor
        // distances the log itself shows to be millimetres (pr/74
        // conn3_unreachable anchor_dis 0.3 cm).  When a component's closest
        // approach to a main-vertex-reachable vertex is within this radius,
        // bridge it with a real rough-path segment BEFORE clustering_points,
        // so the BFS reaches it and it is classified conn-1 naturally.
        // shower_conn3_unreachable stays on as the backstop for wider gaps.
        // Internal units (config surface takes cm).  C++ default 0 = off =
        // legacy = byte-identical.
        double m_conn3_stitch_max{0};

        // doc sbnd_xin/docs/pr/74 round 4 K6 (SBND 18255 evt 506746 seg
        // 21048).  A stopping muon that emits a Michel electron at the
        // neutrino vertex is reconstructed as ONE EM shower: track/shower
        // SEPARATION flags the muon kShowerTrajectory on its wiggliness
        // (segment_is_shower_trajectory), then
        // segment_determine_shower_direction_trajectory forces pdg 11 with
        // the sentinel score 100, and shower clustering absorbs the muon,
        // the Michel and its blob into a single "107 MeV electron" starting
        // at the neutrino vertex -- an electron where the charge says muon.
        //
        // The ordinary track PID cannot arbitrate: separate_track_shower
        // never runs it on a shower-flagged segment
        // (NeutrinoTrackShowerSep.cxx:298-317), and when it IS run (measured,
        // doc pr/74 round 4 Phase A) it ABSTAINS on 21048 -- pdg_code 0, no
        // store.  The discriminating evidence is topological, not template:
        // a Michel is emitted at a RANDOM angle off the stopping point and is
        // TERMINAL, whereas an electron trunk continues near-forward into its
        // own cascade.
        //
        // When on, a pass at the tail of the FINAL examine_direction (the
        // TaggerCheckNeutrino.cxx:1418 call, flag_final -- once, on the
        // neutrino main cluster, immediately before shower clustering)
        // demotes a main-vertex kShowerTrajectory segment to a stopping muon
        // when ALL of: length in [min_len, max_len]; median dQ/dx >=
        // mip_lo x MIP median (a stopping muon's last ~20 cm averages ~1.6x,
        // a single-electron trunk sits at ~1x); its far vertex has degree
        // exactly 2 (one stop, one Michel -- a real cascade branches); the
        // track length beyond that vertex is below max_far_len (terminal);
        // the one sibling there is shower-like; and the kink between them is
        // at least min_kink_deg.  The muon keeps kMuonStemGuard so
        // stem_backfill cannot absorb it straight back.  C++ default false =
        // legacy = byte-identical.
        bool   m_shower_traj_michel_stem{false};
        double m_michel_stem_traj_min_len{15*units::cm};
        double m_michel_stem_traj_max_len{45*units::cm};
        double m_michel_stem_traj_mip_lo{1.3};
        // Deliberately NOT m_michel_stem_max_far_len (P2 above): that member
        // is a VETO ceiling for the F14 rescue and this one is an ACCEPT
        // ceiling here.  Same number today, opposite roles -- sharing it
        // would move one pass silently when the owner tunes the other.
        double m_michel_stem_traj_max_far_len{40*units::cm};
        double m_michel_stem_traj_min_kink_deg{40.0};

        // doc sbnd_xin/docs/pr/44 -- a MULTI-segment long-muon pseudo-shower
        // seeded by shower_clustering_with_nv_in_main_cluster (cached
        // particle_type recorded 13 at the seed) must not have its start
        // segment majority-voted to electron by the update_particle_type call
        // that follows completion: the vote counts every non-proton member
        // (muons included) as shower_length, so a pure muon chain always trips
        // `shower_length > track_length` and the start segment is relabelled
        // 13 -> 11.  That update_particle_type call is a toolkit-only addition
        // (18f09178); the prototype goes straight to the deliberate
        // long-muon -> EM reclass loop (NeutrinoID_shower_clustering.h:
        // 1709-1717) and never re-types a long-muon start segment there.
        // When on, showers cached type +-13 skip the vote at that ONE site
        // (in_main_cluster seeding); every other update_particle_type site is
        // untouched.  SBND 18255 evt 142421: restores the ~143 cm collinear
        // MIP chain 7023->7024->7018 to muon; the fake "e- 163 MeV" (which
        // paired into the pi0) is no longer seeded because the 1.2 cm stub
        // candidate's shower shrinks to {7023, 7082} and fails the
        // n_multi_vtx acceptance.  C++ default false = legacy =
        // byte-identical.
        bool   m_shower_long_muon_keep_type{false};

        // doc sbnd_xin/docs/pr/40 round 10 -- segment_dqdx_spares_electron_
        // reclass (PRSegmentFunctions.cxx, doc pr/40 F2) already guards
        // examine_all_showers' cluster-wide "every non-shower segment here
        // becomes electron" reclassification (NeutrinoTrackShowerSep.cxx
        // ~2070-2105) with a flat median-dQ/dx test: ratio>1.75 (proton) or
        // ratio<1.2 (clean MIP) spares the segment.  A segment whose ratio
        // falls in [1.2, 1.75] gets no protection there -- and that gap is
        // exactly where a real, disconnected muon fragment can sit near
        // end-of-range.  SBND 18255 evt 314507 seg 17002 (32.3 cm, xMIP
        // 1.57x -- inside the gap): the flat guard doesn't fire,
        // examine_all_showers force-relabels it 13->11 (traced with
        // WCT_PID_WRITE_DEBUG=2 to NeutrinoTrackShowerSep.cxx:2091), even
        // though its own Bragg/dE-dx-template PID (segment_do_track_pid,
        // PRSegmentFunctions.cxx) already scored it 0.082 -- a confident
        // fit -- moments earlier.  Above 20 cm the electron template is
        // never in that PID's competition (PRSegmentFunctions.cxx ~2540),
        // so a real (<1.0, i.e. not the 100 "unscored" sentinel) score
        // there is unambiguous muon-or-proton evidence regardless of what
        // particle_info ended up attached.  When on, threaded into
        // examine_all_showers' reclassification test as an OR alongside
        // segment_reclass_dqdx_guard's own check (new helper
        // segment_bragg_spares_electron_reclass): a segment longer than
        // 20 cm with particle_score < 1.0 is spared the same way a
        // ratio-confirmed proton or clean MIP already is.  Purely
        // protective -- can only add a spare, never remove one the flat
        // guard already grants -- so it is additive to
        // shower_reclass_dqdx_guard, not a replacement.
        //
        // Restricted to is_main_cluster (examine_all_showers' own local
        // variable, no new plumbing).  Owner Bee review found a good Bragg
        // score is not reliable evidence inside a SATELLITE cluster (SBND
        // 18255-259542, cluster 124, disjoint from the main interaction):
        // a photon's early conversion stem can score well against the
        // muon template over the same 20-35 cm comparison window before
        // the cascade visibly multiplies, and satellite clusters already
        // have their own dedicated EM-vs-track classifier
        // (kine_drop_stray_satellites, NeutrinoKinematics.cxx, doc pr/92)
        // that correctly kept 259542 as EM.  314507 (the motivating case)
        // sits in the main cluster and is unaffected by this restriction.
        // C++ default false = legacy = byte-identical.
        bool   m_shower_bragg_protect_start_segment{false};

        // doc sbnd_xin/docs/pr/43 round 2 K1 -- the single-muon selection in
        // examine_direction vetoes a muon candidate only when a proton sits
        // at its IMMEDIATE far vertex (1-hop n_proton check); a proton
        // reached through a short degree-2 continuation stub is invisible.
        // SBND 18255 evt 54351: candidate 17007 (54.2 cm) wins over 17010
        // (42.6 cm) by length though its chain 17007 -> 17005 (2.7 cm stub)
        // terminates in charge-confirmed proton 17011 -- a muon cannot end
        // in a proton.  When on, the veto walks the bounded non-shower
        // degree-2 chain (segment_chain_has_proton, <=3 hops); a chain-vetoed
        // candidate is demoted to pion together with its continuation stubs
        // (segment_chain_continuation) and the selection re-picks among the
        // remaining candidates (17010 -> mu-).  Guard: the chain veto only
        // disqualifies when at least one chain-proton-free candidate exists
        // at the vertex; otherwise legacy selection stands (no new
        // demote-all-arms cases).  C++ default false = legacy.
        bool   m_single_muon_proton_chain_veto{false};

        // doc sbnd_xin/docs/pr/43 round 2 K2 -- the same selection loop
        // SKIPS out-edges that belong to a long-muon accumulation chain
        // (`segments_in_long_muon`), so the long muon neither competes for
        // nor claims the vertex muon slot, and a second, shorter pdg-13 arm
        // wins as "the" muon.  SBND 18255 evt 56463: the 411 cm chain
        // 14005+14007 is the muon (owner), yet 14006 (60.2 cm) also stays
        // mu- because 14005 was invisible to the competition.  When on, a
        // long-muon out-edge claims the muon slot with the chain's summed
        // length (deterministic: IndexedSegmentSet order); it is itself
        // never demoted (not in pion_sgs), and other pdg-13 arms demote to
        // pion through the existing conversion loop.  C++ default false.
        bool   m_single_muon_long_muon_claim{false};

        // doc sbnd_xin/docs/pr/46 -- long-muon stub bridge.  The long-muon
        // chain walk (find_cont_muon_segment) requires the junction angle
        // between 15 cm FITTED directions to be < 10 deg (15 deg when the
        // incoming segment is < 6 cm).  A broken long muon whose first piece
        // is a short vertex stub fails that test on the stub's own noisy
        // direction: SBND 18255 evt 55595, stub seg 7 (2.2-2.5 cm, seed
        // dqdx_ratio 0.486) vs the 192.9 cm MIP muon seg 5 (ratio 1.086)
        // measures 30.5-35.4 deg, so the chain never forms, the muon-slot
        // competition crowns the 2.7 cm sibling stub "mu- 48 MeV", and the
        // gateway stub demotes to "pi+ 5 MeV" with the 181.5 cm muon as its
        // child.  When on, the walk accepts a bridge candidate when ALL
        // hold: incoming segment < 6 cm, candidate > 35 cm (short downstream
        // muons keep genuine pi->mu decay kinematics), candidate median
        // dQ/dx ratio < 1.3 (unchanged MIP test), fitted angle < 45 deg
        // (Phase-0 separation: 55595 needs >= 35.4, genuine-pion evt 66118
        // measures 70-82 and must not merge), and the junction has NO other
        // track-like out-edge > 10 cm (not shower-flagged, pdg != +-11 --
        // owner criteria: multiple substantial outgoing tracks veto; a
        // delta-ray electron or tiny fragment does not).  Applied ONLY via
        // the flag_stub_bridge parameter from the formation walk in
        // examine_direction; examine_main_vertices_local and the NuMu
        // tagger keep legacy behavior.  Acceptance (45/35/size>1) and the
        // seed gate are unchanged.  Sample footprint: 2/206 diagnostic
        // events (55595, 61461).  C++ default false = legacy.
        bool   m_long_muon_stub_bridge{false};

        // doc sbnd_xin/docs/pr/43 round 2 K3 -- late particle-info/flag
        // reconciliation pass (reconcile_particle_flags), called from
        // TaggerCheckNeutrino AFTER shower_clustering_with_nv and BEFORE the
        // taggers, so taggers, kine, Bee PF tree and PR display all see one
        // consistent labeling.  Two rules: (1) forced-electron terminal
        // rescue -- a main-vertex proton's degree-2 continuation chain ending
        // in a segment that was FORCED pdg 11 with sentinel score 100 by
        // segment_determine_shower_direction_trajectory's unconditional
        // branches (never actually PID'd) gets ordinary track PID re-run and
        // a confident non-electron conclusion adopted (fresh ParticleInfo,
        // never in-place set_pdg); intermediate pdg-13 stubs (<=15 cm)
        // relabel pion (owner chain proton -> pi+ -> mu-, evt 57661).
        // (2) consistency guard -- a segment whose final pdg is a track type
        // (13/211/2212) but still carries kShowerTrajectory/kShowerTopology
        // has the stale flags cleared (pr/40 F7 precedent), and a
        // single-segment wrapper Shower whose segment is now a confirmed
        // track is dissolved so it renders as a live track node; long-muon
        // pseudo-showers (cached type +-13) are exempt from dissolution.
        // C++ default false = pass never runs = byte-identical.
        bool   m_pid_flag_reconcile{false};

        // ---- Detector-extent literals (doc sbnd_xin/docs/pr/2 sec. 2e(iv)) ----
        // The uBooNE active volume the prototype was written against is
        // y in [-116, +117] cm, z in [0, 1037] cm, x in [0, 256] cm
        // (prototype_base/pid/apps/wire-cell-prod-nue.cxx:417).  Every default
        // below is the prototype literal, so absent config keys are
        // byte-identical.  Values are INTERNAL units.
        //
        // cosmic_tagger() "reaches the top of the detector" tests, part D
        // (prototype NeutrinoID_cosmic_tagger.h:594-796).  Each is a fixed
        // DISTANCE BELOW THE TOP FACE in uBooNE terms (117 - cut), which is
        // what a downward cosmic entering the top satisfies regardless of
        // detector height -- so a taller detector moves them up, it does not
        // scale them.  SBND (top y = +200 cm) uses 183 / 185 / 163 / 133.
        double m_cosmic_y_top_main{100 * units::cm};    // :1175, 17 cm below top:
                                                        // the MAIN cluster's own highest
                                                        // point, relaxing the vertical-angle
                                                        // cut from 20 to 30 degrees.
        double m_cosmic_y_top_strict{102 * units::cm};  // :1190, 15 cm below top: event
                                                        // highest point, single-cosmic branch.
        double m_cosmic_y_top_loose{80 * units::cm};    // :1191, 37 cm below top: event
                                                        // highest point, global gate on the
                                                        // whole flagp_cosmic decision.
        double m_cosmic_y_small_piece{50 * units::cm};  // :1073, 67 cm below top: a <3 cm
                                                        // cluster counts as cosmic debris
                                                        // (acc_small_length) only if its PCA
                                                        // centre sits in the upper region.
        // Upstream-z preference in compare_main_vertices (:875) and
        // compare_main_vertices_global (:3001): each candidate is penalised by
        // (z - min_z) / m_vertex_z_prior_scale, competing with the +0.25-per-
        // track bonuses.  The uBooNE 200 cm gives ~5.2 penalty units across the
        // 1037 cm detector.  On a shorter detector the same 200 cm shrinks the
        // prior's dynamic range, which matters most in the _global comparison
        // (candidates from DIFFERENT clusters of the beam bundle, separations up
        // to the full detector length).  SBND uses 100 cm ~ 200 x 501/1037.
        double m_vertex_z_prior_scale{200 * units::cm};

        // Bundle the do_track_pid-related knobs for the free segment functions.
        TrackPidOptions track_pid_options() const {
            TrackPidOptions o;
            o.mip_dqdx = m_mip_dqdx;
            o.proton_dir_vote = m_proton_dir_vote;
            o.proton_dir_score_max = m_proton_dir_score_max;
            o.proton_dir_asym_min = m_proton_dir_asym_min;
            o.endpoint_trim_retry = m_endpoint_trim_retry;
            o.track_comp_empty_abstain = m_track_comp_empty_abstain;   // doc pr/31 §12 F6
            o.dir_track_median_local = m_dir_track_median_local;       // doc pr/31 §12 F4
            o.track_pid_persist_dqdx = m_track_pid_persist_dqdx;       // doc pr/40 F1
            o.track_pid_persist_4mom = m_track_pid_persist_4mom;       // doc pr/40 round 2 F4
            o.track_pid_persist_dqdx_electron_guard = m_track_pid_persist_dqdx_electron_guard; // doc pr/40 round 5 F9
            return o;
        }

        // Beam-line reference directions used by ssm_tagger, in the DETECTOR
        // frame.  They feed the 8 ssm_*_angle_{target,absorber} BDT features
        // via safe_acos(dir.dot(ref)) at NeutrinoTaggerSSM.cxx:908-909 (the
        // 10 cm initial-direction pair) and :1124-1125 (the nu/con_nu/prim_nu/
        // track momentum pairs).
        // Defaults = the prototype's hard-coded uBooNE numbers
        // (NeutrinoID_ssm_tagger.h): BNB target and NuMI absorber as seen from
        // MicroBooNE.  Absent config keys => byte-identical.  SBND needs its
        // own BNB-target vector and the NuMI-absorber features have no obvious
        // SBND meaning -- see wcp-porting-img sbnd_xin/docs/pr/2 sec. 2e(i);
        // the value itself is still an open gap, only its location moved.
        //
        // NOTE these are NOT unit vectors: |target| = 0.99866,
        // |absorber| = 1.00970.  safe_acos() clamps the dot product, so the
        // absorber angle saturates at exactly 0 for true angles out to ~7.9
        // deg.  Kept verbatim for prototype parity -- whoever supplies a
        // properly normalized SBND vector should expect these feature
        // distributions to shift even at the same physical direction.
        WireCell::Vector m_ssm_target_dir{0.46, 0.05, 0.885};
        WireCell::Vector m_ssm_absorber_dir{0.33, 0.75, -0.59};

        // Charge -> kinetic-energy calibration constants (NeutrinoEnergyReco.cxx).
        // Defaults = the uBooNE literals they replaced => byte-identical when the
        // config keys are absent.  See KineChargeOptions above.
        KineChargeOptions m_kine_charge{};

        // Muon median-dQ/dx-vs-length envelope {c0, c1, pivot, power}:
        //     dQ_dx_cut = c0 + c1 * pow(pivot / L, power)
        // a dimensionless multiple of m_mip_dqdx_median.  The default is the
        // prototype's empirical uBooNE stopping-muon refit (of the commented-out
        // 0.85 + 0.95*sqrt(25cm/L) family): short tracks are Bragg-dominated so
        // their MEDIAN dQ/dx sits well above MIP (2.5x at 5 cm), relaxing toward
        // the 0.8866 plateau asymptote.  Appears at NINE sites (numu x2,
        // vertex-finder, nue x4, ssm, cosmic), all as literals before this knob
        // (doc pr/2 sec 2e(iv) undercounted 3).  pivot is INTERNAL units.
        // Defaults = the uBooNE literals => absent config keys are
        // byte-identical.
        std::array<double, 4> m_muon_dqdx_curve{{0.8866, 0.9533, 18 * units::cm, 0.4234}};

        // The envelope, length in INTERNAL units.  Inline and shaped exactly
        // like the literal expression it replaced (c0 + c1*pow(180.0/L, c3)),
        // so the defaults are bit-identical.
        double muon_dqdx_cut(double length) const {
            return m_muon_dqdx_curve[0] +
                   m_muon_dqdx_curve[1] * std::pow(m_muon_dqdx_curve[2] / length, m_muon_dqdx_curve[3]);
        }
        // Same envelope for a length already in cm (the ssm_tagger convention,
        // whose sg_length values are divided by units::cm at the source).
        // pivot/units::cm = 18.0 exactly, so this too is bit-identical.
        double muon_dqdx_cut_cm(double length_cm) const {
            return m_muon_dqdx_curve[0] +
                   m_muon_dqdx_curve[1] * std::pow((m_muon_dqdx_curve[2] / units::cm) / length_cm, m_muon_dqdx_curve[3]);
        }

        // Single-photon stem dE/dx conversion (NeutrinoTaggerSinglePhoton.cxx).
        // Default false = the inline float inverse Modified-Box frozen at
        // uBooNE's field (A=1.0, B=0.255, rho=1.38, E=0.273 kV/cm), kept for
        // byte-identity.  When true, shw_sp_vec_{median,mean}_dedx go through
        // m_recomb_model->dE() instead -- the same configured model the rest of
        // the chain uses (docs/pr/2 sec 2e(i) third correctness item).
        bool m_sp_dedx_use_recomb_model{false};
        // The hard mip_quality cut on shw_sp_vec_mean_dedx (MeV/cm).  Stored as
        // float so the default compares bit-identically to the legacy 2.3f
        // literal.  uBooNE-tuned against the inline (0.273 kV/cm) dE/dx scale;
        // coupled to the knob above -- retune when switching models.
        float m_sp_mean_dedx_cut{2.3f};
        // The configured recombination model, set by TaggerCheckNeutrino::visit()
        // alongside the kine transfer.  Null when the component has none.
        IRecombinationModel::pointer m_recomb_model{};

        // ---- doc sbnd_xin/docs/pr/36 §10 tagger-stage knobs -----------------
        // All C++ default false = today's path = byte-identical; the SBND
        // operating point decides which are on.  Set by
        // TaggerCheckNeutrino::visit() alongside the other transfers.
        //
        // F4 (= P3 + P5 + track_overclustering's muon_segs): iterate the three
        // std::set<SegmentPtr>/std::set<ShowerPtr> accumulation loops in
        // NeutrinoTaggerNuE.cxx in graph-index order instead of
        // pointer-address order.  This is an M4 / CLAUDE.md §2 determinism
        // house-rule fix, NOT a fidelity fix: the prototype does the same
        // address-ordered thing (std::set<ProtoSegment*>,
        // NeutrinoID_nue_tagger.h:1036 iterated :1139), so turning this on
        // moves the toolkit FURTHER from the prototype's (unreproducible)
        // order while making our own runs order-stable.
        bool m_tagger_ordered_segment_sets{false};
        // F5 (= P6): pick the stem fit endpoint by the prototype's wcpt
        // identity rule (NeutrinoID_nue_tagger.h:71-75) instead of the
        // nearest-fit-endpoint proximity substitute, at the 18
        // seg_endpoint_near call sites.  WCPoint::index is not ported
        // (PRCommon.h:99) so identity is EXACT wcpt-position equality --
        // deliberately no tolerance (doc pr/36 §10.15b).
        bool m_stem_endpoint_wcpt_parity{false};
        // F6 (= P8): broken_muon_id's multi-cluster test counts distinct
        // cluster IDS (prototype NeutrinoID_nue_tagger.h ~:1183) instead of
        // distinct Facade::Cluster POINTERS.  Equal iff cluster<->id is
        // injective within the event -- doc 53's real_cluster_id epochs are
        // why this is a knob and not an assumption.
        bool m_broken_muon_cluster_id_count{false};
        // F7 (= P4): compute the prototype's per-tagger neutrino_type verdict
        // bitmask (see TaggerInfo::neutrino_type).  The matching T_tagger
        // branch is booked by UbooneTaggerOutputVisitor under the same
        // config key.
        bool m_neutrino_type_bitmask{false};

        // ---- doc sbnd_xin/docs/pr/33 §10 EM-shower-clustering knobs -------
        // All C++ default false = today's path = byte-identical; the SBND
        // operating point decides which are on.  Set by
        // TaggerCheckNeutrino::visit() alongside the pr/36 transfers.
        //
        // F1 (= P1): restore the prototype's calculate_num_daughter_tracks
        // callee (prototype NeutrinoID_shower_clustering.h:140 /
        // NeutrinoID_em_shower.h:17) where the port calls
        // calculate_num_daughter_showers.  Two knobs because the two sites
        // err in opposite directions and one knob could not attribute which
        // site moved a decision.
        bool m_daughter_count_proto_main_vertex{false};
        bool m_daughter_count_proto_examine_showers{false};
        // F2 (= P2): read the PDG off the object the prototype reads.
        // from_start_segment: four sites read shower->get_particle_type()
        // where the prototype reads get_start_segment()->get_particle_type()
        // (NeutrinoID_shower_clustering.h:1716,:387,:394,:497,:511).
        // from_shower_type: the one inverted site (:525 here, prototype
        // :1800) reads the start segment where the prototype reads the
        // shower.  exact_muon_test: drop std::abs at the two sites where
        // the prototype's muon test is exact (proto :1716, em_shower.h:10).
        // Prototype parity at the :170 site needs from_start_segment AND
        // exact_muon_test together; either alone is neither tree's behavior.
        bool m_shower_pdg_from_start_segment{false};
        bool m_shower_pdg_from_shower_type{false};
        bool m_shower_pdg_exact_muon_test{false};
        // F3 (= P3): the two pi0 finders share one id allocation stream
        // (prototype: member NeutrinoID.h:1982 mutated at :933 and :688),
        // so a with-vertex and a without-vertex pi0 in one event cannot
        // collide on pio_id.  Scoped to the two finders only --
        // shower_clustering_with_nv stays by value because the caller's
        // local is also seeded by-reference into ssm_tagger (doc pr/33
        // §10.10 amendment 1); the seeding-at-0 divergence is a recorded
        // gap, not part of this knob.
        bool m_pi0_id_shared_allocator{false};
        // F4 (= P6): the :797 is_shower mirror of the prototype's
        // get_flag_shower() gains the missing third disjunct
        // (ProtoSegment.cxx:1305-1312: fabs(particle_type)==11).
        bool m_shower_flag_pdg_electron{false};
        // F5 (= P12): shower_less's same-index tie-break compares
        // get_shower_id() instead of pointer addresses.  NOT a
        // prototype-fidelity fix (the prototype orders these by pointer):
        // a CLAUDE.md §2 determinism house-rule fix, shipped with an
        // unconditional fallback-hit counter (g_pr33_audit).
        bool m_shower_less_id_tiebreak{false};
        // doc pr/39: prototype-parity exclusion of a shower's own start
        // vertex from the farthest-vertex search that sets data.end_point
        // (Shower::calculate_kinematics{,_long_muon}).  Same rule as
        // Shower::fill_sets's exclude_start_vertex (doc pr/38 round 2,
        // WCShower.cxx:547 map_vtx_segs never holds start_vertex), extended
        // to the two calculate_kinematics searches and the long-muon search,
        // which pr/38 round 2 did not touch.  Without it, a detached
        // (conn_type 2/3) shower's end_point can land exactly on its own
        // start vertex (e.g. the neutrino vertex) instead of growing away
        // from it.  Default false = legacy search over every node,
        // byte-identical.
        bool m_shower_endpoint_exclude_start_vertex{false};

        // doc sbnd_xin/docs/pr/91 round 1 F1 -- the same farthest-vertex search
        // must also skip a node that NO member segment of the shower touches.
        //
        // Such orphan nodes are a toolkit-only artefact: set_start_vertex()
        // calls add_vertex() (PRShower.h fill_sets docstring, doc pr/38), so a
        // conn-2/3 shower's view carries a vertex from somebody else's cluster
        // -- routinely the nearest main-cluster or in-shower vertex chosen by
        // shower_clustering_in_other_clusters.  While that shower owns it,
        // m_shower_endpoint_exclude_start_vertex above hides it.  Then
        // Shower::add_shower's node loop imports it wholesale into an absorber,
        // where it is no longer the start vertex, the exclusion stops applying,
        // and it wins the search.  shower_dedup_start_seg (doc pr/84 round 3,
        // SBND ON) is what makes that routine: on SBND 169626/174752/347129/
        // 394532 it produced 6 orphan imports and 5 wrong end points -- 394532's
        // 30 MeV and 66 MeV showers report each other's charge as their end,
        // and a single-13.6cm-segment shower in 347129 reports an end 67.8 cm
        // away in another cluster.  Forcing the knob off reproduces 0 orphans /
        // 0 bad end points / 0 add_shower calls on the same four events.
        //
        // end_point feeds the Bee PF node's `data.end` and every
        // cal_dir_3vector / angle consumer; kine_charge, kine_energy_particle
        // and kine_reco_Enu sum over member segments and do NOT move.
        // Default false = legacy search, byte-identical.
        bool m_shower_endpoint_skip_orphan_vtx{false};

        // doc pr/91 round 3: Shower::complete_structure_with_start_segment's
        // frontier test is view MEMBERSHIP (!has_node), not visitation -- a
        // vertex added by set_start_vertex()/set_start_segment()/
        // add_segment()/add_shower() but never actually scanned by a
        // flood-fill worklist is permanently unreachable.  When true, the
        // frontier test switches to Shower::m_walked_nodes, the prototype's
        // map_vtx_segs equivalent (WCPPID::WCShower::set_start_vertex,
        // WCShower.cxx:529-532, never touches map_vtx_segs -- only the
        // pointer and connection type).  Instrumented on SBND nueCC evt
        // 168596: the 2039 MeV electron's original start vertex 14027 (seeded
        // by set_start_vertex, later superseded when examine_showers
        // re-seats the shower onto the main vertex) walls off a 4.74 cm
        // proton stub and, past it, a 7.7 cm electron stub 96% inside this
        // shower's own point cloud -- left as a separate shower that later
        // pairs into a spurious pi0 via the ownership-free kine_charge sum.
        // See PRShower.h complete_structure_with_start_segment for the full
        // mechanism.  Default false = legacy has_node()-gated frontier,
        // byte-identical.
        bool m_shower_walk_visited_parity{false};

        // 2D charge maps cached for the duration of shower_clustering_with_nv.
        // Populated once by collect_charge_maps(); reused by calculate_shower_kinematics
        // and all cal_kine_charge call sites to avoid O(N_hits) re-collection per shower.
        ChargeMap m_charge_2d_u, m_charge_2d_v, m_charge_2d_w;
        WireMap   m_map_apa_ch_plane_wires;
        // doc pr/35 §10.11a (F3) discriminating diagnostic: TrackFitting's
        // m_charge_data size at the moment the cache above was filled.  The
        // segment cal_kine_charge overload compares it against the size at its
        // own (fresh, per-call) collection; a mismatch means the cache and a
        // fresh collection would see different charge, i.e. caching the segment
        // path would be a behaviour change, not a perf change.  Log-only.
        size_t m_charge_data_size_at_collect{0};

        // Populate the cached charge maps from track_fitter.
        // Call once at the start of shower_clustering_with_nv; the maps are then
        // valid for the entire call tree beneath it.
        void collect_charge_maps(TrackFitting& track_fitter);
        std::vector<VertexPtr> find_cluster_vertices(Graph& graph, const Facade::Cluster& cluster);
        std::vector<SegmentPtr> find_cluster_segments(Graph& graph, const Facade::Cluster& cluster);
        bool clean_up_graph(Graph& graph, const Facade::Cluster& cluster);

        SegmentPtr init_first_segment(Graph& graph, Facade::Cluster& cluster, Facade::Cluster* main_cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, bool flag_back_search = true);

        // Isochronous first-segment endpoints (m_iso_endpoint; doc
        // sbnd_xin/docs/pr/24 round 2).  Returns true and fills p1/p2 (both
        // snapped to Steiner terminals of "steiner_pc") when the cluster
        // passes the sheet gate; false => caller uses the legacy boundary
        // search unchanged.
        bool find_iso_first_segment_endpoints(const Facade::Cluster& cluster, Facade::geo_point_t& p1, Facade::geo_point_t& p2) const;

        // find the shortest path using steiner graph
        std::vector<Facade::geo_point_t> do_rough_path(const Facade::Cluster& cluster,Facade::geo_point_t& first_point, Facade::geo_point_t& last_point);
        // doc sbnd_xin/docs/pr/51 round 5 (m_steiner_gap_penalty): make sure
        // the cluster carries the support-penalized "steiner_graph_gap"
        // flavor, building it lazily (once per cluster) from the installed
        // "steiner_graph".  Returns true iff the flavor exists (do_rough_path
        // then routes on it); false => caller stays on "steiner_graph"
        // unchanged (knob off, missing prerequisites, or empty graph).
        bool ensure_steiner_gap_graph(const Facade::Cluster& cluster);
        std::vector<Facade::geo_point_t> do_rough_path_reg_pc(const Facade::Cluster& cluster, Facade::geo_point_t& first_point, Facade::geo_point_t& last_point,  std::string graph_name = "relaxed_pid");
        // create a segment given a path
        SegmentPtr create_segment_for_cluster(WireCell::Clus::Facade::Cluster& cluster, IDetectorVolumes::pointer dv, const std::vector<Facade::geo_point_t>& path_points, int dir = 0);
        // create a segment given two vertices, null, if failed
        SegmentPtr create_segment_from_vertices(Graph& graph, Facade::Cluster& cluster, VertexPtr v1, VertexPtr v2, IDetectorVolumes::pointer dv);
        // replace a segment and vertex with another segment and vertex, assuming the original vertex only connect to this segment
        bool replace_segment_and_vertex(Graph& graph, SegmentPtr& seg, VertexPtr& vtx, std::list<Facade::geo_point_t>& path_point_list, Facade::geo_point_t& break_point, IDetectorVolumes::pointer dv);
        bool replace_segment_and_vertex(Graph& graph, SegmentPtr& seg, VertexPtr old_vertex, VertexPtr new_vertex, IDetectorVolumes::pointer dv);
        bool break_segment_into_two(Graph& graph, VertexPtr vtx1, SegmentPtr seg, VertexPtr vtx2, std::list<Facade::geo_point_t>& path_point_list1, Facade::geo_point_t& break_point, std::list<Facade::geo_point_t>& path_point_list2, IDetectorVolumes::pointer dv, SegmentPtr& out_seg2);


        // return the point and its index in the steiner tree as a pair
        std::pair<Facade::geo_point_t, size_t> proto_extend_point(const Facade::Cluster& cluster, Facade::geo_point_t& p, Facade::geo_vector_t& dir, Facade::geo_vector_t& dir_other, bool flag_continue, std::vector<Facade::geo_point_t>* walk_history = nullptr);
        // return Steiner Graph path in wcps_list1 and wcps_list2
        bool proto_break_tracks(const Facade::Cluster& cluster, const Facade::geo_point_t& first_wcp, Facade::geo_point_t& curr_wcp, const Facade::geo_point_t& last_wcp, std::list<Facade::geo_point_t>& wcps_list1, std::list<Facade::geo_point_t>& wcps_list2, bool flag_pass_check);
        // breaking segments ...
        bool break_segments(Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, std::vector<SegmentPtr>& remaining_segments, float dis_cut = 0);
        // merge vertices within 0.1 cm after break_segments, then refit if changed
        bool merge_nearby_vertices(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        // merge two segments to one
        bool merge_two_segments_into_one(Graph& graph, SegmentPtr& seg1, VertexPtr& vtx, SegmentPtr& seg2, IDetectorVolumes::pointer dv);
        // merge vertex into another
        bool merge_vertex_into_another(Graph& graph, VertexPtr& vtx_from, VertexPtr& vtx_to, IDetectorVolumes::pointer dv);

        // get direction with  distance cut ... 
        Facade::geo_vector_t vertex_get_dir(VertexPtr& vertex, Graph& graph, double dis_cut = 5*units::cm);
        Facade::geo_vector_t vertex_segment_get_dir(VertexPtr& vertex, SegmentPtr& segment, Graph& graph, double dis_cut = 2*units::cm);


        // Structure examination
        void examine_structure(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv); // call examine_structure_1 and examine_structure_2      
        bool examine_structure_1(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_2(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        bool examine_structure_3(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_4(VertexPtr vertex, bool flag_final_vertex, Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // identify other segments giving the graph ...
        void find_other_segments(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv , bool flag_break_track =true, double search_range = 1.5*units::cm, double scaling_2d = 0.8);
        VertexPtr find_vertex_other_segment(Graph& graph, Facade::Cluster& cluster, SegmentPtr seg, bool flag_forwrard, Facade::geo_point_t& wcp, TrackFitting& track_fitter, IDetectorVolumes::pointer dv );
        std::tuple<VertexPtr, SegmentPtr, Facade::geo_point_t> check_end_point(Graph& graph, Facade::Cluster& cluster, std::vector<Facade::geo_point_t>& tracking_path, bool flag_front = true, double vtx_cut1 = 0.9 * units::cm, double vtx_cut2 = 2.0 * units::cm, double sg_cut1  = 2.0 * units::cm, double sg_cut2  = 1.2 * units::cm);
        bool modify_vertex_isochronous(Graph& graph, Facade::Cluster& cluster, VertexPtr vtx, VertexPtr v1, SegmentPtr sg, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool modify_segment_isochronous(Graph& graph, Facade::Cluster& cluster, SegmentPtr sg1, VertexPtr v1, SegmentPtr sg, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, double dis_cut = 6*units::cm, double angle_cut = 15, double extend_cut = 15*units::cm);

        // examine segment
        void examine_segment(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool crawl_segment(Graph& graph, Facade::Cluster& cluster, SegmentPtr seg, VertexPtr vertex, TrackFitting& track_fitter, IDetectorVolumes::pointer dv );
        void examine_partial_identical_segments(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        //examine vertices
        void examine_vertices(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex = nullptr);
        bool examine_vertices_1(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex = nullptr);
        bool examine_vertices_1p(Graph&graph, VertexPtr v1, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_vertices_2(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex = nullptr);
        bool examine_vertices_4(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex = nullptr);
        bool examine_vertices_4p(Graph&graph, VertexPtr v1, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        Facade::geo_point_t get_local_extension(Facade::Cluster& cluster, const Facade::geo_point_t& wcp);
        void examine_vertices_3(Graph& graph, Facade::Cluster& main_cluster, std::pair<VertexPtr, VertexPtr> main_cluster_initial_pair_vertices, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // master pattern recognition function
        // particle_data (doc pr/48): the dE/dx stopping templates for the
        // two-end break pass (break_two_end_dqdx).  nullptr (the default)
        // makes that pass a no-op regardless of m_two_end_break.
        bool find_proto_vertex(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, bool flag_break_track = true, int nrounds_find_other_tracks = 2, bool flag_back_search = true, const Clus::ParticleDataSet::pointer& particle_data = nullptr);

        // doc sbnd_xin/docs/pr/48 sec 6 -- the two-end residual-range
        // back-to-back break pass.  Gated on m_two_end_break; runs inside
        // find_proto_vertex after examine_vertices on the main cluster only.
        // Topology gate: exactly one segment of this cluster longer than
        // m_teb_stub_max, both its endpoints inside the fiducial volume
        // (grouping's FiducialUtils; a missing FiducialUtils never fires --
        // conservative).  On a segment_two_end_break_scan accept, breaks the
        // segment at the located fit point (break_segment), marks the new
        // vertex VertexFlags::kProtectedBreak, and returns true (the caller
        // re-fits via its existing final do_multi_tracking).  Single break
        // per cluster by construction.
        bool break_two_end_dqdx(Graph& graph, Facade::Cluster& cluster, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data);
        
        void init_point_segment(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // examine the structure of the patterns ... 
        bool examine_structure_final_1(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_final_1p(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_final_2(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_final_3(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_final(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // EM shower related
        void clustering_points(Graph& graph, Facade::Cluster& cluster, const IDetectorVolumes::pointer& dv, const std::string& cloud_name = "associate_points", double search_range = 1.2*units::cm, double scaling_2d = 0.7);
        void separate_track_shower(Graph&graph, Facade::Cluster& cluster);
        // doc sbnd_xin/docs/pr/59 round 2: inert unless m_assoc_full_recluster.
        // A segment created after a cluster's association pass (e.g. examine_
        // structure_final*/examine_vertices_1's replacement polylines inside
        // determine_main_vertex) can end up with a null/empty associate_points
        // cloud -- see clustering_points_segments's doc pr/59 sentinels for the
        // full diagnosis.  If ANY of the cluster's current segments has a null
        // cloud, clears associate_points on every segment in the cluster and
        // re-runs clustering_points_segments over the whole cluster (a fresh
        // Voronoi + ghost-removal competition -- must include every segment,
        // never just the orphans, or the orphan would win points by default
        // with no sibling able to contest).  Then re-runs the
        // separate_track_shower classification (is_shower_topology, then
        // is_shower_trajectory if not topology-shower) ONLY on the segments
        // that were null before the clear, since those are the only ones a
        // classification pass could not previously reach and re-classifying an
        // already-correct segment is pure blast radius.  No-op (returns 0) if
        // no segment in the cluster is orphaned -- untouched clusters stay
        // byte-identical even with the knob on.  Returns the number of
        // segments rescued.
        size_t reassociate_cluster_orphans(Graph& graph, Facade::Cluster& cluster, const IDetectorVolumes::pointer& dv);
        // Direction
        void determine_direction(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        std::pair<int, double> calculate_num_daughter_showers(Graph& graph, VertexPtr vertex, SegmentPtr segment, bool flag_count_shower = true);
        std::pair<int, double> calculate_num_daughter_tracks(Graph& graph, VertexPtr vtx, SegmentPtr sg, bool flag_count_shower = false, double length_cut = 0);
        // find_cont_muon_segment_nue: find adjacent segment continuing in same direction (MIP-like)
        std::pair<SegmentPtr, VertexPtr> find_cont_muon_segment_nue(Graph& graph, SegmentPtr sg, VertexPtr vtx, bool flag_ignore_dQ_dx = false);
        void examine_good_tracks(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data);
        // about fix maps
        void fix_maps_multiple_tracks_in(Graph& graph, Facade::Cluster& cluster);
        void fix_maps_shower_in_track_out(Graph& graph, Facade::Cluster& cluster);
        void improve_maps_one_in(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_strong_check = true);
        void improve_maps_shower_in_track_out(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_strong_check = true);
        void improve_maps_no_dir_tracks(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void improve_maps_multiple_tracks_in(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void judge_no_dir_tracks_close_to_showers(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, IDetectorVolumes::pointer dv);
        bool examine_maps(Graph&graph, Facade::Cluster& cluster);
        void examine_all_showers(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data);
        void shower_determining_in_main_cluster(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, IDetectorVolumes::pointer dv);
        // main_vertex: doc sbnd_xin/docs/pr/40 round 2 F5 -- defaults to nullptr so
        // any hypothetical other call site (there is currently exactly one,
        // examine_direction, which always has a main_vertex) stays source-
        // compatible; segment_has_proton_daughter always returns false for a
        // null vertex, so an omitted argument is byte-identical to F5 off.
        void set_default_shower_particle_info(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, VertexPtr main_vertex = nullptr);

        // doc sbnd_xin/docs/pr/40 round 4 F8 -- see m_muon_multi_proton_pion's
        // docstring above.  Called right after set_default_shower_particle_info
        // (same call site, same per-cluster main_vertex, still the last word
        // before shower_clustering_with_nv runs).  No-op (returns immediately)
        // when m_muon_multi_proton_pion is false.
        void override_muon_multi_proton_pion(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, VertexPtr main_vertex = nullptr);

        // doc sbnd_xin/docs/pr/40 round 6 F14 -- see m_michel_stem_muon_
        // rescue's docstring above.  Same call site as F8 (right after
        // override_muon_multi_proton_pion, last word before
        // shower_clustering_with_nv), same per-cluster main_vertex.  No-op
        // (returns immediately) when m_michel_stem_muon_rescue is false.
        void override_michel_stem_muon(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, VertexPtr main_vertex = nullptr);

        // doc sbnd_xin/docs/pr/74 round 4 K6 -- see m_shower_traj_michel_stem's
        // docstring above.  A SEPARATE pass, not a branch inside
        // override_michel_stem_muon: F14 is SBND production ON and its body
        // stays byte-for-byte untouched (fork by duplication).  Called only
        // from the flag_final examine_direction, i.e. once, on the neutrino
        // main cluster, immediately before shower_clustering_with_nv.  No-op
        // (returns immediately) when m_shower_traj_michel_stem is false.
        void override_shower_traj_michel_stem(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, VertexPtr main_vertex);

        // doc sbnd_xin/docs/pr/43 round 2 K3 -- late particle-info/flag
        // reconciliation pass; see m_pid_flag_reconcile above.  Gated
        // internally on that knob (no-op when false).
        void reconcile_particle_flags(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ShowerIntMap& map_shower_pio_id, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);

        // PCA calculation
        std::pair<Facade::geo_point_t, Facade::geo_vector_t> calc_PCA_main_axis(std::vector<Facade::geo_point_t>& points);

        // vertex related functions 
        bool search_for_vertex_activities(Graph& graph, VertexPtr vertex, std::vector<SegmentPtr>& segments_set, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, double search_range = 1.5*units::cm);
        // existing_segments stays a plain std::set<SegmentPtr> (pointer-keyed)
        // ON PURPOSE.  It is iterated once, in the Case-5 block, where the body
        // only takes a running min over three distances -- order-insensitive.
        // Its find()/count() must stay POINTER identity: an index-ordered set
        // would compare by Segment::get_graph_index(), and that value is NOT
        // guaranteed unique across live Segment objects -- two sources exist.
        //   (1) Segment::m_graph_index defaults to SIZE_MAX (PRSegment.h:153),
        //       so EVERY segment not yet passed to PR::add_segment compares
        //       equal to every other such segment.
        //   (2) PR::add_segment on a vertex pair that already carries an edge
        //       takes the "edge already existed" path (PRGraph.cxx:86-89): it
        //       overwrites g[desc].segment and copies the existing edge index
        //       into the new segment, leaving the displaced one aliased.
        // MEASURED: making that swap moved kine_reco_Enu on SBND evt 239794
        // from 2930 to 1687 MeV -- that observation stands and is why this
        // container is not converted.  The CAUSE is NOT confirmed: probing this
        // exact set at its call site (NeutrinoVertexFinder.cxx:2289) on six
        // events INCLUDING evt 239794 found zero index collisions and zero
        // unindexed members, so neither (1) nor (2) reproduces here.  See
        // sbnd_xin/docs/pr/28 §14.8.  Do not cite the mechanism as established;
        // do keep the container pointer-keyed.
        bool eliminate_short_vertex_activities(Graph& graph, Facade::Cluster& cluster, VertexPtr main_vertex, std::set<SegmentPtr>& existing_segments, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        std::tuple<bool, int, int> examine_main_vertex_candidate(Graph& graph, VertexPtr vertex);
        VertexPtr compare_main_vertices_all_showers(Graph& graph, Facade::Cluster& cluster, std::vector<VertexPtr>& vertex_candidates, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        float calc_conflict_maps(Graph& graph, VertexPtr vertex);
        VertexPtr compare_main_vertices(Graph& graph, Facade::Cluster& cluster, std::vector<VertexPtr>& vertex_candidates);
        std::pair<SegmentPtr, VertexPtr> find_cont_muon_segment(Graph &graph, SegmentPtr sg, VertexPtr vtx, bool flag_ignore_dQ_dx = false, bool flag_stub_bridge = false);
        bool examine_direction(Graph& graph, VertexPtr vertex, VertexPtr main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_final = false);

        bool fit_vertex(Facade::Cluster& cluster, VertexPtr vertex, VertexPtr main_vertex, std::vector<SegmentPtr>& sg_set, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        void improve_vertex(Graph& graph, Facade::Cluster& cluster, VertexPtr& main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_search_vertex_activity = true , bool flag_final_vertex = false);
        // doc sbnd_xin/docs/pr/50 (see the m_vertex_kink_snap member block):
        // main-vertex kink-consistency snap.  Returns true iff it re-broke a
        // segment and re-seated main_vertex (by reference) on the new
        // kProtectedBreak vertex, followed by one full refit.  Gated on
        // m_vertex_kink_snap (default false => immediate return, no side
        // effects).
        bool snap_main_vertex_to_kink(Graph& graph, Facade::Cluster& cluster, VertexPtr& main_vertex, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        // doc sbnd_xin/docs/pr/51 (see the m_main_vertex_graph_audit member
        // block): near-vertex graph audit on the main cluster.  Returns true
        // iff at least one operation edited the graph (followed by one full
        // refit); may re-seat main_vertex's position in place (op3), never
        // replaces the pointer.  Gated on m_main_vertex_graph_audit (default
        // false => immediate return, no side effects).
        bool main_vertex_graph_audit(Graph& graph, Facade::Cluster& cluster, VertexPtr main_vertex, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        // doc pr/83 r3 (sec 9.5, Mechanism C): one unscoped duplicate-
        // corridor pass over `cluster` (op1's metric and merge recipe,
        // fork-by-duplication -- op1 itself is untouched), refitting the
        // cluster if anything merged.  Called from swap_main_cluster on the
        // ABANDONED main cluster when m_swap_orphan_dup_audit is set; that
        // gate lives in the caller so this function itself is knob-free.
        bool orphan_dup_audit(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        // doc sbnd_xin/docs/pr/84 round 2 (F3, see the m_conn3_stitch_max
        // member block): bridge disconnected main-cluster components whose
        // closest approach to the reachable side is within
        // m_conn3_stitch_max, then one full refit.  Returns true iff at
        // least one bridge was created.  Gated on m_conn3_stitch_max > 0
        // (default 0 => immediate return, no side effects).
        bool stitch_disconnected_main_cluster(Graph& graph, Facade::Cluster& cluster, VertexPtr main_vertex, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        // doc sbnd_xin/docs/pr/51 round 4 (see the m_rough_path_probe member
        // block): diagnostic-only TRACE probe for the near-vertex short-cut
        // investigation.  Never edits graph, segment, or fit content; always
        // returns without side effects.  Gated on m_rough_path_probe
        // (default false => immediate return).
        void rough_path_probe(Graph& graph, Facade::Cluster& cluster, VertexPtr main_vertex, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        void determine_main_vertex(Graph& graph, Facade::Cluster& cluster, VertexPtr& main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void change_daughter_type(Graph& graph, VertexPtr vertex, SegmentPtr segment, int particle_type, double mass, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void examine_main_vertices_local(Graph& graph, std::vector<VertexPtr>& vertices, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);

        // cluster functions ...
        Facade::geo_vector_t calc_dir_cluster(Graph& graph, Facade::Cluster& cluster, const Facade::geo_point_t& orig_p, double dis_cut);
        // doc pr/83 r3: the three trailing params exist only for
        // m_swap_orphan_dup_audit (the abandoned-cluster duplicate audit,
        // sec 9.5).  Defaulted null => every legacy call compiles and, knob
        // off, behaves byte-identically; knob on with any of them missing
        // logs a "skipped" TRACE rather than silently auditing nothing.
        Facade::Cluster* swap_main_cluster(Facade::Cluster& new_main_cluster, Facade::Cluster& old_main_cluster, std::vector<Facade::Cluster*>& other_clusters, Graph* graph = nullptr, TrackFitting* track_fitter = nullptr, IDetectorVolumes::pointer dv = nullptr);
        void examine_main_vertices(Graph& graph, ClusterVertexMap& map_cluster_main_vertices, Facade::Cluster*& main_cluster, std::vector<Facade::Cluster*>& other_clusters, TrackFitting* track_fitter = nullptr, IDetectorVolumes::pointer dv = nullptr);

        VertexPtr compare_main_vertices_global(Graph& graph, std::vector<VertexPtr>& vertex_candidates, Facade::Cluster& main_cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        Facade::Cluster* check_switch_main_cluster(Graph& graph, ClusterVertexMap map_cluster_main_vertices, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        Facade::Cluster* check_switch_main_cluster_2(Graph& graph, VertexPtr temp_main_vertex, Facade::Cluster* max_length_cluster, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, TrackFitting* track_fitter = nullptr, IDetectorVolumes::pointer dv = nullptr);
        // doc sbnd_xin/docs/pr/51 round 3: map_cluster_main_vertices and
        // main_cluster are BY REFERENCE (matching determine_overall_main_vertex_DL
        // below) so a cluster swap decided internally (examine_main_vertices /
        // check_switch_main_cluster[_2] -> swap_main_cluster, which already
        // mutates persistent Flags::main_cluster + the by-ref other_clusters)
        // is not silently discarded.  Callers that want the pre-round-3
        // discard behaviour pass throwaway local copies of both arguments --
        // see TaggerCheckNeutrino.cxx's m_main_vertex_swap_apply gate.
        VertexPtr determine_overall_main_vertex(Graph& graph, ClusterVertexMap& map_cluster_main_vertices, Facade::Cluster*& main_cluster, std::vector<Facade::Cluster*>& other_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_dev_chain = true);

        // Deep-learning vertex refinement.  Returns true if the DL network changed
        // the selected vertex (in which case the traditional determine_overall_main_vertex
        // should NOT be called).  Returns false if the DL network was unavailable, failed,
        // or did not improve on the candidate vertices (fall back to traditional).
        // dl_vtx_swap_guard (doc sbnd_xin/docs/pr/51, 18255-506746): when
        // true, the RERANK branch skips candidates hosted on a different
        // cluster than the current main cluster (each skip is one
        // "dl_swap_guard:" DEBUG sentinel), so an accepted DL vertex can
        // never swap the main cluster; if no candidate survives, the normal
        // traditional fallback runs.  506746: a single confident uBooNE-net
        // voxel (raw score 0.576 -> s_dl = +576 at score_scale 1000, which
        // swamps every +-2 structural term) moved the vertex 28 cm onto a
        // non-flash-matched cluster.  Default false => byte-identical.
        bool determine_overall_main_vertex_DL(Graph& graph, ClusterVertexMap& map_cluster_main_vertices, Facade::Cluster*& main_cluster, std::vector<Facade::Cluster*>& other_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, const std::string& dl_weights, double dl_vtx_cut, double dQdx_scale = 0.1, double dQdx_offset = -1000.0, bool flag_rerank = false, int dl_vtx_top_k = 5, double dl_vtx_min_accept_score = 4.0, double dl_vtx_score_scale = 1000.0, bool dl_vtx_swap_guard = false, double dl_vtx_topo_weight = 0.0, double dl_vtx_topo_center = 0.0);

        // global information transfer
        void transfer_info_from_segment_to_cluster(Graph& graph, Facade::Cluster& cluster, const std::string& cloud_name = "associated_points");

        // print information
        void print_segs_info(Graph& graph, Facade::Cluster& cluster, VertexPtr vertex= nullptr);

        // shower related functions
        void update_shower_maps(IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters);
        // doc pr/84 round 3 -- see m_shower_dedup_start_seg.  Returns the
        // number of showers absorbed (0 => nothing changed).
        int merge_showers_sharing_start_segment(IndexedShowerSet& showers);
        // doc sbnd_xin/docs/pr/74 round 2 K4 -- see m_shower_stem_backfill.
        void stem_backfill(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedSegmentSet& segments_in_long_muon, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void shower_clustering_with_nv_in_main_cluster(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void shower_clustering_connecting_to_main_vertex(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters);
        void shower_clustering_with_nv_from_main_cluster(Graph& graph, VertexPtr main_vertex, Facade::Cluster* main_cluster, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters);
        void shower_clustering_with_nv_from_vertices(Graph& graph, VertexPtr main_vertex, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void calculate_shower_kinematics(IndexedShowerSet& showers, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void examine_merge_showers(IndexedShowerSet& showers, VertexPtr main_vertex, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void shower_clustering_in_other_clusters(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_save = true);
        void examine_shower_1(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void examine_showers(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        // acc_segment_id is int& (doc pr/33 F3): the two pi0 finders share
        // one allocation stream when m_pi0_id_shared_allocator is on.  The
        // caller (shower_clustering_with_nv) binds a local copy, so nothing
        // propagates past it either way.
        void id_pi0_with_vertex(int& acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Pi0KineFeatures& pio_kine, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void id_pi0_without_vertex(int& acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Pi0KineFeatures& pio_kine, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void shower_clustering_with_nv(int acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Pi0KineFeatures& pio_kine, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);


        // deghost related functions ...
        void order_clusters(Graph& graph, std::vector<Facade::Cluster*>& ordered_clusters, std::map<Facade::Cluster*, std::vector<SegmentPtr> >& map_cluster_to_segments, std::map<Facade::Cluster*, double>& map_cluster_total_length);
        void deghost_clusters(Graph& graph, std::vector<Facade::Cluster*>& all_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        void order_segments(std::vector<SegmentPtr>& ordered_segments, std::vector<SegmentPtr>& segments);
        void deghost_segments(Graph& graph, ClusterVertexMap map_cluster_main_vertices, std::vector<Facade::Cluster*>& all_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        void deghosting(Graph& graph, ClusterVertexMap& map_cluster_main_vertices, std::vector<Facade::Cluster*>& all_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv );

        // energy calculation ...
        double cal_corr_factor(WireCell::Point& pt, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // kinematics output ...
        void init_tagger_info(TaggerInfo& ti);

        // cosmic tagger
        bool bad_reconstruction(Graph& graph, VertexPtr main_vertex, ShowerPtr shower,
                                bool flag_fill = false, TaggerInfo* ti = nullptr);
        bool cosmic_tagger(Graph& graph, VertexPtr main_vertex,
                           IndexedShowerSet& showers,
                           ShowerSegmentMap& map_segment_in_shower,
                           VertexShowerSetMap& map_vertex_to_shower,
                           IndexedSegmentSet& segments_in_long_muon,
                           Facade::Cluster* main_cluster,
                           std::vector<Facade::Cluster*>& all_clusters,
                           IDetectorVolumes::pointer dv,
                           TaggerInfo& ti);

        // numu tagger
        // count_daughters: counts track branches and total branches at the far end of a muon.
        // For a segment: BFS from the vertex closer to main_vertex through sg, count what is beyond.
        // For a shower:  finds the last segment of the long-muon chain, then counts from its near end.
        std::pair<int, int> count_daughters(Graph& graph, SegmentPtr sg, VertexPtr main_vertex);
        std::pair<int, int> count_daughters(Graph& graph, ShowerPtr shower, VertexPtr main_vertex,
                                            IndexedSegmentSet& segments_in_long_muon);
        // numu_tagger: fills TaggerInfo numu_cc_* BDT features and returns
        //   {flag_long_muon, max_muon_length}.
        // Prototype: WCPPID::NeutrinoID::numu_tagger() in NeutrinoID_numu_tagger.h.
        std::pair<bool, double> numu_tagger(Graph& graph,
                                            VertexPtr main_vertex,
                                            IndexedShowerSet& showers,
                                            IndexedSegmentSet& segments_in_long_muon,
                                            Facade::Cluster* main_cluster,
                                            TaggerInfo& ti);

        // nue_tagger: fills TaggerInfo NuE BDT features and returns flag_nue.
        // Prototype: WCPPID::NeutrinoID::nue_tagger() in NeutrinoID_nue_tagger.h.
        bool nue_tagger(Graph& graph,
                        Facade::Cluster* main_cluster,
                        VertexPtr main_vertex,
                        int apa, int face,
                        IndexedShowerSet& showers,
                        VertexShowerSetMap& map_vertex_to_shower,
                        IndexedShowerSet& pi0_showers,
                        ShowerIntMap& map_shower_pio_id,
                        std::map<int, std::vector<ShowerPtr>>& map_pio_id_showers,
                        std::map<int, std::pair<double,int>>& map_pio_id_mass,
                        IDetectorVolumes::pointer dv,
                        ParticleDataSet::pointer particle_data,
                        double muon_length,
                        TaggerInfo& ti);

        // singlephoton_tagger: fills TaggerInfo shw_sp_* BDT features and returns flag_sp.
        // Prototype: WCPPID::NeutrinoID::singlephoton_tagger() in NeutrinoID_singlephoton_tagger.h.
        // apa/face are derived internally from dv->contained_by(main_vertex_pt) so callers
        // do not need to know detector geometry details.
        // geom_helper (doc pr/36 §10.4, F3 = P2): the prototype SCE-corrects
        // every position feeding the shw_sp_* features
        // (func_pos_SCE_correction at NeutrinoID_singlephoton_tagger.h:13,
        // :103, :132, :317; :222 is commented out there).  Null (the legacy
        // and the knob-off value) => raw positions, byte-identical.
        bool singlephoton_tagger(Graph& graph,
                                 Facade::Cluster* main_cluster,
                                 VertexPtr main_vertex,
                                 IndexedShowerSet& showers,
                                 VertexShowerSetMap& map_vertex_to_shower,
                                 ShowerIntMap& map_shower_pio_id,
                                 std::map<int, std::vector<ShowerPtr>>& map_pio_id_showers,
                                 std::map<int, std::pair<double,int>>& map_pio_id_mass,
                                 IDetectorVolumes::pointer dv,
                                 WireCell::IClusGeomHelper::pointer geom_helper,
                                 TaggerInfo& ti);

        // ssm_tagger: fills TaggerInfo ssm_* and ssmsp_* features and returns flag_ssm.
        // Prototype: WCPPID::NeutrinoID::ssm_tagger() in NeutrinoID_ssm_tagger.h.
        bool ssm_tagger(Graph& graph,
                        VertexPtr main_vertex,
                        IndexedShowerSet& showers,
                        ShowerVertexMap& map_vertex_in_shower,
                        ShowerSegmentMap& map_segment_in_shower,
                        const Pi0KineFeatures& pio_kine,
                        int flag_ssmsp,
                        int& acc_segment_id,
                        const ParticleDataSet::pointer& particle_data,
                        const IRecombinationModel::pointer& recomb_model,
                        TaggerInfo& ti);

        KineInfo fill_kine_tree(VertexPtr main_vertex, IndexedShowerSet& showers,
                                const Pi0KineFeatures& pio_kine,
                                Graph& graph, TrackFitting& track_fitter,
                                IDetectorVolumes::pointer dv,
                                WireCell::IClusGeomHelper::pointer geom_helper,
                                const Clus::ParticleDataSet::pointer& particle_data,
                                const IRecombinationModel::pointer& recomb_model,
                                // doc pr/92: pi0-paired showers (protected from the
                                // stray-satellite drop; the TrackFitting copy is stale
                                // at call time) and the out-param collecting dropped
                                // shower ids for the Bee PF tree's matching gate.
                                const IndexedShowerSet& pi0_showers = IndexedShowerSet{},
                                std::set<int>* dropped_satellites = nullptr);
        // Convenience overloads: collect 2D charge maps internally (safe for isolated calls).
        double cal_kine_charge(ShowerPtr Shower, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        /// doc pr/99 round 3 (C1/C1b): final knob-gated recompute of every
        /// shower's kine_charge -- cross-shower 2D-cell ownership
        /// (m_kine_charge.dedup) and/or ephemeral member-true clouds
        /// (m_kine_charge.rebuild).  No-op unless one of the two flags is
        /// set.  Called from shower_clustering_with_nv after all shower-
        /// structure passes and BEFORE the pi0 finders, so pairing, taggers,
        /// fill_kine_tree and the PF/Bee display all read one consistent
        /// energy set while every mid-pipeline gate saw legacy values.
        void recompute_shower_kine_charge_final(IndexedShowerSet& showers, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        double cal_kine_charge(SegmentPtr segment, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        // Fast overload: reuse pre-collected 2D charge maps (avoids O(N_hits) collection per call).
        // Collect maps once with track_fitter.collect_2D_charge() and pass here when calling in a loop.
        double cal_kine_charge(ShowerPtr shower,
            const ChargeMap& charge_2d_u, const ChargeMap& charge_2d_v, const ChargeMap& charge_2d_w,
            const WireMap& map_apa_ch_plane_wires,
            TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

    };
}
