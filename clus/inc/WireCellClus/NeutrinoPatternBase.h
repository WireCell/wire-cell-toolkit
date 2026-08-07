#include "WireCellClus/PRGraph.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/PRShower.h"
#include "WireCellClus/NeutrinoTaggerInfo.h"
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

        // Cathode kink veto (doc sbnd_xin/docs/pr/20 Part II, B0).  Passed to
        // segment_search_kink from break_segments: a candidate fit point within
        // m_cathode_kink_xcut of the cathode plane at m_cathode_x is skipped, so
        // the ~2 cm transverse cathode mismatch cannot invent a vertex that breaks
        // a crossing cosmic into two particles.  C++ default 0 => no point is ever
        // skipped => byte-identical to the pre-pr/20 behavior.
        double m_cathode_x{0};
        double m_cathode_kink_xcut{0};

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

        // doc sbnd_xin/docs/pr/43 F1 -- extend the muon-candidate loop's
        // 1-hop proton veto (NeutrinoVertexFinder.cxx examine_direction,
        // `n_proton` at the candidate's immediate far vertex) to a bounded
        // multi-hop chain walk (segment_chain_has_proton,
        // PRSegmentFunctions.cxx, max_hops=3).  Run 18255 evt 54351:
        // muon-candidate seg 17007 (54.2cm) wins the single-muon-per-
        // cluster selection over seg 17010 (42.6cm) purely by length,
        // because 17007's own far vertex has no proton -- but two hops
        // further, through a 2.7cm muon-pdg continuation stub (seg 17005),
        // sits a charge-confirmed proton (seg 17011, median/MIP > 1.75).
        // 17007 cannot be the muon in a chain that terminates in a proton;
        // with this knob on, 17007 is excluded from candidacy (same
        // n_proton==0-style gate, now chain-aware), so 17010 (with no
        // proton anywhere in its own chain) wins uncontested; 17007 and the
        // stub 17005 both go through the existing "convert non-muon
        // candidates to pion" path unchanged.  Applies identically to the
        // pdg==0 branch just below (an undetermined-direction segment at
        // this vertex is also chain-walked before falling to the dqdx_ratio
        // proton/pion split).  Designed divergence (the prototype's muon
        // selection has no proton-chain veto at all, 1-hop or multi-hop);
        // see porting_dictionary.md.  C++ default false = legacy =
        // byte-identical.
        bool   m_muon_chain_proton_veto{false};

        // doc sbnd_xin/docs/pr/43 -- Shower::update_particle_type relabels
        // its start segment to electron (shower_length > track_length) but
        // leaves the Shower's OWN cached particle_type at a stale prior
        // value (e.g. 13, set when the segment was first detected as part
        // of a long-muon chain).  fill_bee_pf_tree's shower-leaf renderer
        // (MultiAlgBlobClustering.cxx make_shower_leaf) reads the CACHED
        // shower->get_particle_type(), so run 18255 evt 56463 displays a
        // Bee PF root node "mu- 903 MeV" for a segment (14005) whose own
        // particle_info is already pdg 11 -- the toolkit's own logic had
        // already reclassified it away from muon, but the display never
        // caught up, reading as a second, phantom muon at the vertex next
        // to the genuine one (14006).  Same divergence class as doc pr/35
        // F1 kine_shower_pdg_live, on the Bee side instead of the kine
        // side.  Threaded to Shower::update_particle_type's
        // refresh_type_cache parameter (PRShower.h) at every call site in
        // NeutrinoShowerClustering.cxx.  Designed as a cache-consistency
        // fix, not a new heuristic -- see porting_dictionary.md.  C++
        // default false = legacy = byte-identical.
        bool   m_shower_type_cache_refresh{false};

        // doc sbnd_xin/docs/pr/43 F3 -- see TrackPidOptions::
        // shower_traj_dqdx_guard's comment in PRSegmentFunctions.h.  Threaded
        // via track_pid_options() to segment_determine_shower_direction_
        // trajectory.  C++ default false = legacy = byte-identical.
        bool   m_shower_traj_dqdx_guard{false};

        // doc sbnd_xin/docs/pr/43 F3b -- run 18255 evt 57661: once
        // shower_traj_dqdx_guard rescues seg 18007 to a confident muon,
        // the short (5cm) stub 18005 between it and the main-vertex-
        // emanating proton (18003) is STILL pdg 13 (it was never itself a
        // shower-flagged segment, just an ordinary per-segment muon call
        // that nothing downstream revisited), giving proton -> muon ->
        // muon.  The owner's reading is proton -> pion -> muon: a
        // transitional stub between a proton and a confirmed muon reads as
        // the pion, not a second muon.  Threaded via override_shower_traj_
        // chain_pion (NeutrinoPatternBase.cxx), same call-site convention
        // as F8/F14: relabels every pdg-13 segment in a main-vertex
        // proton's short (<=15cm each) non-shower continuation chain to
        // pion EXCEPT the last (deepest) one, and only when every segment
        // in that chain is pdg 13 (so an already-mixed or single-segment
        // chain -- the ordinary proton -> muon topology -- is untouched).
        // The <=15cm-per-segment gate (same scale as the isolated_length_
        // cut convention) exists to avoid misreading a single long muon
        // track that pattern recognition merely fragmented into two
        // collinear pieces as a proton-pion-muon chain. Designed divergence
        // (the prototype has no such chain rule); see porting_dictionary.md.
        // C++ default false = legacy = byte-identical.
        bool   m_shower_traj_chain_pion{false};

        // doc sbnd_xin/docs/pr/43 F4 -- run 18255 evt 142421: fill_kine_tree
        // (NeutrinoKinematics.cxx) runs its OWN BFS with the same over-wide
        // shower barrier doc pr/38 already corrected on the Bee PF-tree
        // side (pf_shower_vertex_barrier).  `shower->fill_sets(...,
        // /*flag_exclude_start_segment=*/false)` at fill_kine_tree's
        // opening takes fill_sets's DEFAULT `exclude_start_vertex=false`,
        // so a detached (conn-type 2/3) shower's start vertex still blocks
        // this BFS even after pr/38 stopped it blocking the PF-tree BFS --
        // the vehicle (fill_sets's exclude_start_vertex parameter) exists,
        // was simply never threaded here.  Same fix, same knob shape as
        // pr/38's pf_shower_vertex_barrier: (a) pass exclude_start_vertex
        // = true at that fill_sets call, and (b) after the BFS, push every
        // still-unreached, non-shower, main-cluster segment with a
        // dirsign into kine_energy_particle via the existing
        // push_segment_kine lambda (prototype fill_kine_tree gives every
        // main-cluster segment a node in a flat pass before the BFS;
        // NeutrinoID_kine.h has no disconnected-fragment exclusion).  This
        // MOVES kine_reco_Enu and therefore numu_score/nue_score -- unlike
        // the PF-tree-only pr/38 fix, this is not display-only, so it
        // carries its own score-shift review (doc pr/43 G4).  Designed
        // parity with pr/38, not a prototype port; see
        // porting_dictionary.md.  C++ default false = legacy =
        // byte-identical.
        bool   m_kine_shower_vertex_barrier{false};

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
            o.shower_traj_dqdx_guard = m_shower_traj_dqdx_guard; // doc pr/43 F3
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
        bool find_proto_vertex(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, bool flag_break_track = true, int nrounds_find_other_tracks = 2, bool flag_back_search = true);
        
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

        // doc sbnd_xin/docs/pr/43 F3b -- see m_shower_traj_chain_pion's
        // docstring above.  No-op when m_shower_traj_chain_pion is false.
        void override_shower_traj_chain_pion(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, VertexPtr main_vertex = nullptr);

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
        std::pair<SegmentPtr, VertexPtr> find_cont_muon_segment(Graph &graph, SegmentPtr sg, VertexPtr vtx, bool flag_ignore_dQ_dx = false);
        bool examine_direction(Graph& graph, VertexPtr vertex, VertexPtr main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_final = false);

        bool fit_vertex(Facade::Cluster& cluster, VertexPtr vertex, VertexPtr main_vertex, std::vector<SegmentPtr>& sg_set, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        void improve_vertex(Graph& graph, Facade::Cluster& cluster, VertexPtr& main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_search_vertex_activity = true , bool flag_final_vertex = false);
        void determine_main_vertex(Graph& graph, Facade::Cluster& cluster, VertexPtr& main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void change_daughter_type(Graph& graph, VertexPtr vertex, SegmentPtr segment, int particle_type, double mass, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void examine_main_vertices_local(Graph& graph, std::vector<VertexPtr>& vertices, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);

        // cluster functions ...
        Facade::geo_vector_t calc_dir_cluster(Graph& graph, Facade::Cluster& cluster, const Facade::geo_point_t& orig_p, double dis_cut);
        Facade::Cluster* swap_main_cluster(Facade::Cluster& new_main_cluster, Facade::Cluster& old_main_cluster, std::vector<Facade::Cluster*>& other_clusters);
        void examine_main_vertices(Graph& graph, ClusterVertexMap& map_cluster_main_vertices, Facade::Cluster*& main_cluster, std::vector<Facade::Cluster*>& other_clusters);

        VertexPtr compare_main_vertices_global(Graph& graph, std::vector<VertexPtr>& vertex_candidates, Facade::Cluster& main_cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        Facade::Cluster* check_switch_main_cluster(Graph& graph, ClusterVertexMap map_cluster_main_vertices, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        Facade::Cluster* check_switch_main_cluster_2(Graph& graph, VertexPtr temp_main_vertex, Facade::Cluster* max_length_cluster, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters);
        VertexPtr determine_overall_main_vertex(Graph& graph, ClusterVertexMap map_cluster_main_vertices, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_dev_chain = true);

        // Deep-learning vertex refinement.  Returns true if the DL network changed
        // the selected vertex (in which case the traditional determine_overall_main_vertex
        // should NOT be called).  Returns false if the DL network was unavailable, failed,
        // or did not improve on the candidate vertices (fall back to traditional).
        bool determine_overall_main_vertex_DL(Graph& graph, ClusterVertexMap& map_cluster_main_vertices, Facade::Cluster*& main_cluster, std::vector<Facade::Cluster*>& other_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, const std::string& dl_weights, double dl_vtx_cut, double dQdx_scale = 0.1, double dQdx_offset = -1000.0, bool flag_rerank = false, int dl_vtx_top_k = 5, double dl_vtx_min_accept_score = 4.0, double dl_vtx_score_scale = 1000.0);

        // global information transfer
        void transfer_info_from_segment_to_cluster(Graph& graph, Facade::Cluster& cluster, const std::string& cloud_name = "associated_points");

        // print information
        void print_segs_info(Graph& graph, Facade::Cluster& cluster, VertexPtr vertex= nullptr);

        // shower related functions
        void update_shower_maps(IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters);
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
                                const IRecombinationModel::pointer& recomb_model);
        // Convenience overloads: collect 2D charge maps internally (safe for isolated calls).
        double cal_kine_charge(ShowerPtr Shower, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        double cal_kine_charge(SegmentPtr segment, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        // Fast overload: reuse pre-collected 2D charge maps (avoids O(N_hits) collection per call).
        // Collect maps once with track_fitter.collect_2D_charge() and pass here when calling in a loop.
        double cal_kine_charge(ShowerPtr shower,
            const ChargeMap& charge_2d_u, const ChargeMap& charge_2d_v, const ChargeMap& charge_2d_w,
            const WireMap& map_apa_ch_plane_wires,
            TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

    };
}
