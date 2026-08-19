#ifndef WIRECELL_CLUS_PRSEGMENTFUNCTIONS
#define WIRECELL_CLUS_PRSEGMENTFUNCTIONS

#include "WireCellClus/PRGraph.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/D4Vector.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IRecombinationModel.h"
#include "WireCellClus/ParticleDataSet.h"

#include <array>
#include <functional>

namespace WireCell::Clus::PR {

    using geo_point_t = WireCell::Point;

    /// Replace the segment in the graph with two new segments that meet at a
    /// new vertex nearest to the point.
    ///
    /// The input segment is removed from the graph.
    ///
    /// The point must be withing max_dist of the segment.
    ///
    /// orient_split=true resolves the parent's (front, back) vertices with
    /// find_vertices() before slicing, so each child carries the wcpt/fit
    /// half that actually terminates at its vertices.  false = legacy:
    /// slice by boost source/target, which do NOT track path orientation --
    /// on a reversed edge each child gets the WRONG half (doc
    /// sbnd_xin/docs/pr/83; C++ default false => byte-identical).
    ///
    /// Returns true if the graph was modified.
    std::tuple<bool, std::pair<SegmentPtr, SegmentPtr>, VertexPtr> break_segment(Graph& graph, SegmentPtr seg, Point point, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, const IDetectorVolumes::pointer& dv,
                       double max_dist=1e9*units::cm, bool orient_split=false);
    // patter recognition
    // cathode_kink_xcut > 0 makes the kink search skip candidate fit points within
    // that distance of the cathode plane at cathode_x, where the transverse
    // cathode mismatch fakes a turn and the para_angle guard (which exists to
    // protect isochronous imaging artifacts) is wide open because the crossing is
    // drift-x dominated -- see sbnd_xin/docs/pr/20 Part II.  Default 0 = no point
    // is ever skipped = byte-identical to the pre-pr/20 behavior.
    //
    // cathode_wide_kink_angle > 0 adds a fifth accept path at cathode-crossing
    // fit indices only (doc sbnd_xin/docs/pr/47 sec 8, option O1): the shipped
    // index-windowed refl_angle statistic (min over 2-12-point lever arms) is
    // suppressed by the cathode gap/distortion, so a genuine kink AT the
    // crossing can never reach the legacy thresholds there.  The new test
    // measures the skirt-excluded wide-baseline PCA turn angle across the
    // crossing (segment_cathode_wide_kink_accepts below) and accepts when it
    // is >= cathode_wide_kink_angle degrees.  Default 0 = path fully disabled
    // = byte-identical legacy search.  New algorithm, no prototype counterpart.
    //
    // kink_walk_dqdx_stop (doc sbnd_xin/docs/pr/48 sec 6, 59335 fix (a)): the
    // C4 / straightness accepts set flag_search and then return
    // flag_continue=true UNCONDITIONALLY, bypassing the local-dQ/dx walk gate
    // that every other accept honors -- so proto_extend_point walks past a
    // dQ/dx-confident kink toward the terminus (18255-59335: a correct C4
    // accept 0.28 cm from truth overshoots to 0.47 cm from the far end).
    // true = a BRAGG-HOT kink (local_dQdx > threshold * kink_hot_ratio)
    // stops the walk on the flag_search path too.  kink_hot_ratio defaults
    // to the legacy 25/43 "not too low" scale but callers should pass a
    // genuinely hot ratio (~1.7): at 25/43 nearly every C4 accept qualifies
    // and the footprint is sample-wide (doc pr/48 sec 9.6).  false = legacy
    // = byte-identical.
    //
    // accept_criterion (doc pr/48 sec 6, 59335 fix (b) transport): when
    // non-null, receives which accept fired for the returned kink -- -1 none,
    // 0 the wide-baseline cathode accept (A0), 1-4 the legacy criteria C1-C4.
    // Pure out-param, never read: nullptr (the default) is byte-identical.
    std::tuple<WireCell::Point, WireCell::Vector, WireCell::Vector, bool> segment_search_kink(SegmentPtr seg, WireCell::Point& start_p, const std::string& cloud_name = "fit", double dQ_dx_threshold = 43000/units::cm, double cathode_x = 0, double cathode_kink_xcut = 0, double cathode_wide_kink_angle = 0, double cathode_wide_kink_skirt = 3*units::cm, double cathode_wide_kink_baseline = 15*units::cm, bool kink_walk_dqdx_stop = false, int* accept_criterion = nullptr, double kink_hot_ratio = 25.0/43.0 );

    /// Wide-baseline cathode kink test (doc sbnd_xin/docs/pr/47 sec 8, O1).
    ///
    /// For every cathode crossing of the fit trajectory (consecutive fit
    /// points straddling x = cathode_x), fit a PCA direction to each arm over
    /// the arclength window [skirt, skirt+baseline] measured from the
    /// crossing (points within the skirt are EXCLUDED -- that is where the
    /// transverse cathode distortion lives), orient both directions along the
    /// direction of travel, and take the angle between them: ~0 deg for a
    /// straight through-going track, independent of any pure transverse
    /// offset between the two arms (a translation does not rotate a PCA
    /// axis).  A crossing with < 3 points in either arm's window never fires.
    /// Returns the accepted fit indices: per firing crossing, the
    /// crossing-adjacent index with the larger |x - cathode_x|, stepped
    /// outward (up to 2 indices) to clear the |x| < 0.45 cm active-volume
    /// slab hole, clamped to segment_search_kink's accept contract
    /// 0 < i < fits.size()-1.  Same statistic as the doc pr/47 census script
    /// (scripts/analysis/pr47/cathode_junction_census.py, skirt_turn_angle).
    /// All lengths in internal units; angle_cut_deg in degrees; angle_cut_deg
    /// <= 0 returns empty.
    std::vector<size_t> segment_cathode_wide_kink_accepts(
        const std::vector<Fit>& fits, double cathode_x, double angle_cut_deg,
        double skirt, double baseline);

    /// Wide-baseline PCA turn angle at one fit index (doc
    /// sbnd_xin/docs/pr/48 sec 6).  Same statistic as
    /// segment_cathode_wide_kink_accepts above, stripped of the cathode
    /// scoping: per-arm PCA direction over the arclength window
    /// [skirt, skirt+baseline] either side of fits[idx], each axis oriented
    /// along its own arm's chord, angle between them in degrees.  Returns 0
    /// when either arm has < 3 points in its window (a terminus-adjacent or
    /// short-segment index can never fire).  All lengths internal units.
    double segment_wide_turn_angle(const std::vector<Fit>& fits, size_t idx,
                                   double skirt, double baseline);

    /// doc sbnd_xin/docs/pr/50 (172230-class near-vertex robustness):
    /// result of path_scan_vertex_kink below.
    struct VertexKinkScanResult {
        bool found{false};
        int idx{-1};         ///< index into the oriented point list (0 = vertex end)
        double turn_deg{0};  ///< turn angle at that point, degrees (0 = straight)
        double arc{0};       ///< arclength from pts[0], internal units
    };

    /// doc sbnd_xin/docs/pr/50: ALL interior turns >= angle_cut of an
    /// oriented point path (pts[0] at a vertex) inside the arclength window
    /// [d_min, d_max], in index order.  Per candidate index the turn is
    /// segment_wide_turn_angle over a synthesized Fit vector (it reads only
    /// .point); when the vertex-side arm is too short for the symmetric PCA
    /// windows (< 3 points -- the 2.7 cm stub case) the vertex-side
    /// direction falls back to the chord pts[0] -> pts[i], against the
    /// outward arm's windowed PCA (itself falling back to the chord
    /// pts[i] -> pts.back()).  Candidates keep >= min_arm of outward
    /// arclength support beyond the turn.  The CALLER picks among the
    /// candidates with its own topology evidence -- the strongest turn is
    /// not necessarily the right corner (172230: a secondary wiggle at
    /// 4.9 cm out-turned the true corner at 2.4 cm by 3 deg).
    /// Deterministic integer index scan, no pointer-keyed containers.
    /// All lengths internal units, angles deg.
    std::vector<VertexKinkScanResult> path_scan_vertex_kink(const std::vector<WireCell::Point>& pts,
                                                            double d_min, double d_max,
                                                            double skirt, double baseline,
                                                            double angle_cut, double min_arm);

    /// doc sbnd_xin/docs/pr/54 (18255-142421 separated EM shower with no
    /// fitted trajectory): keep decision for an isolated residual candidate
    /// in find_other_segments -- one whose freshly-fit endpoints touch no
    /// existing segment and whose isochronous snap failed.  Legacy (and
    /// prototype: NeutrinoID_proto_vertex.h:1470-1475 pushes such candidates
    /// into the never-consumed residual_segment_candidates) behaviour is to
    /// discard; `keep_isolated=false` reproduces that unconditionally.  When
    /// the knob is on, keep candidates whose terminal-graph component has at
    /// least min_points points AND whose fitted track length is at least
    /// min_length -- separating "isolated but well-supported" (a separated
    /// shower) from "isolated and sparse" (noise fragments).  Pure
    /// predicate, all lengths internal units.
    bool other_seg_keep_isolated_ok(bool keep_isolated, int component_points,
                                    double track_length, int min_points,
                                    double min_length);

    /// doc sbnd_xin/docs/pr/51 (main-vertex graph audit): fraction of pts_a
    /// whose distance to the NEAREST point of pts_b is <= tol.  The
    /// duplicate-corridor test: two segments riding one charge ribbon give a
    /// high fraction for the shorter path against the longer one (360535:
    /// 0.77-0.80 at 1.4 cm for a ~1 cm-separated parallel pair; 268067:
    /// 0.86 at 0.6 cm for a corridor rider).  Returns 0 when either input
    /// is empty.  Brute force O(|a|*|b|) -- callers pass near-vertex fit
    /// paths (tens of points).  Deterministic, pure geometry.
    double path_overlap_fraction(const std::vector<WireCell::Point>& pts_a,
                                 const std::vector<WireCell::Point>& pts_b,
                                 double tol);

    /// doc sbnd_xin/docs/pr/51 round 5 (steiner_gap_penalty): unsupported
    /// fraction of the straight span a->b of length `length`, sampled every
    /// `step` (inclusive endpoints, nsteps = max(1, round(length/step)) --
    /// the round-4 probe's P3 loop verbatim).  `classify` returns 0=live,
    /// 1=dead, 2=unsupported for a 3D point; the result is
    /// (n_unsup + dead_alpha*n_dead) / n in [0, 1].  Pure function -- the
    /// caller injects the classifier, so doctests need no cluster.
    /// Implemented in NeutrinoSteinerGapGraph.cxx.
    double gap_edge_bad_fraction(const WireCell::Point& a, const WireCell::Point& b,
                                 double length, double step, double dead_alpha,
                                 const std::function<int(const WireCell::Point&)>& classify);

    /// doc sbnd_xin/docs/pr/51 round 6: thresholded weak-charge deficit of
    /// an edge from its endpoint charges,
    ///     0.5 * (max(0, 1 - qa/qref) + max(0, 1 - qb/qref)),
    /// in [0, 1]: 0 when both endpoints are at/above qref, 1 when both are
    /// chargeless.  The saturating base-weight form (Q0/(q+Q0)) is provably
    /// too weak to reroute a short weak-charge chord -- multiplicative
    /// penalty then hurts the longer strong-charge route more -- hence the
    /// hard threshold.  Pure function; doctested without a cluster.
    /// Implemented in NeutrinoSteinerGapGraph.cxx.
    double weak_charge_deficit(double qa, double qb, double qref);

    /// doc sbnd_xin/docs/pr/51 round 7 (robust vertex fit): annulus PCA of a
    /// point path around a vertex, replicating MyFCN::AddSegment's math on an
    /// arbitrary (rin, rout] shell: kept points are those with
    /// rin <= |p - vtx| <= rout, center = kept point closest to vtx,
    /// covariance about center divided by the kept count, eigenvalues
    /// descending with the (0.15 cm)^2 floor added, axes normalized.
    /// npts < 2 marks an unusable window (axes/vals zero).
    struct MvfitAnnulusPca {
        int npts{0};
        WireCell::Point center{0, 0, 0};
        std::array<WireCell::Point, 3> axes{WireCell::Point(0, 0, 0),
                                            WireCell::Point(0, 0, 0),
                                            WireCell::Point(0, 0, 0)};
        std::array<double, 3> vals{0, 0, 0};
    };
    MvfitAnnulusPca mvfit_annulus_pca(const std::vector<WireCell::Point>& pts,
                                      const WireCell::Point& vtx,
                                      double rin, double rout);

    /// doc sbnd_xin/docs/pr/51 round 7: dynamic outer direction-window radius
    /// for a leg of the given endpoint chord -- clamp(frac*chord, rmin, rmax).
    /// A longer track earns a longer direction lever arm.
    double mvfit_rout_dyn(double chord, double frac, double rmin, double rmax);

    /// doc sbnd_xin/docs/pr/51 round 7: per-leg disagreement gate.  True iff
    /// the folded angle (acos|dot|, eigenvector signs are arbitrary) between
    /// the production inner axis and the outer window's leading axis exceeds
    /// angle_deg, the outer window holds >= min_pts points, and its
    /// anisotropy sqrt(vals[0]/vals[1]) >= min_aniso (a fanning shower's
    /// outer window fails this intrinsically).  Pure function; doctested
    /// without a cluster.  Implemented in MyFCN.cxx.
    bool mvfit_leg_disagrees(const WireCell::Point& inner_axis,
                             const MvfitAnnulusPca& outer,
                             double angle_deg, int min_pts, double min_aniso);

    /// doc sbnd_xin/docs/pr/51: total associated charge of a segment --
    /// the sum of fit.dQ over valid fit points (fit.valid() && fit.dx > 0
    /// && fit.dQ >= 0, the same filter as segment_median_dQ_dx).  The
    /// duplicate-corridor merge keeps the higher-integrated-charge member.
    /// Returns 0 for a segment with no valid fits.
    double segment_integrated_dQ(SegmentPtr seg);

    /// Options for segment_two_end_break_scan (doc sbnd_xin/docs/pr/48 sec 6
    /// -- the two-end residual-range back-to-back break).  Every value here
    /// only matters once the owning component's two_end_break knob is ON;
    /// defaults are the doc pr/48 starting operating point.
    /// - mip_dqdx: flat-template amplitude handed to do_track_comp (the
    ///   50e3-e/cm role, same as TrackPidOptions::mip_dqdx).
    /// - mip_dqdx_median: the 43e3-e/cm median-threshold role; scales
    ///   abs_end_min.
    /// - min_len: shortest segment ever scanned.
    /// - min_arm / min_arm_pts: both arms of a candidate break index must
    ///   have at least this arclength AND this many fit points (4 points - 1
    ///   skip_stop_sample = 3 KS samples, the minimum for a meaningful
    ///   shape comparison).
    /// - accept_range: do_track_comp window for the ACCEPTANCE tier at the
    ///   located index (15 cm = the tier eval_ks_ratio's constants were
    ///   tuned at).  LOCALIZATION does not use the templates at all: a
    ///   stopping-template score is flat in the break index whenever one
    ///   arm is a template-perfect prefix/suffix (any prefix of a Bragg is
    ///   a Bragg), so the junction is located by its measured markers
    ///   instead -- the 3-point-mean dQ/dx DIP (route R1) and the
    ///   wide-baseline turn maximum (route R2).
    /// - rise_r1 / rise_r2: two-end dQ/dx rise floors (end-window median /
    ///   interior median, max over 2/4/8 cm end windows clamped to L/4) for
    ///   the dQ/dx-only route R1 and the angle-assisted route R2.
    /// - abs_end_min: R1's absolute end-median floor in units of
    ///   mip_dqdx_median (dead-region fake-dip compensator).
    /// - dip_floor: R1's candidate-dip floor in units of mip_dqdx_median --
    ///   a 3-point-mean dip BELOW this is instrumental (dead region / gap:
    ///   missing charge), not a junction (both particles sit at maximal
    ///   residual range there, so a genuine junction dip stays near MIP),
    ///   and is never a break candidate.
    /// - score_cap_r1 / score_cap_r2: acceptance-tier cap on
    ///   max over arms of min(muon score, proton score).
    /// - turn_angle (deg) / turn_baseline / turn_skirt: R2's wide-baseline
    ///   PCA turn requirement at the located index; turn_angle <= 0 disables
    ///   route R2 entirely.
    /// - turn_min_arm_frac: when > 0, R2's turn argmax runs a PREFERENCE
    ///   pass first over indices where each PCA arm's achievable arclength
    ///   (bounded by the segment end beyond the skirt) reaches this fraction
    ///   of turn_baseline; that winner is kept only if it clears turn_angle
    ///   on its own, else the legacy unrestricted argmax stands.  A starved
    ///   near-end window reads PCA jitter, not heading (doc
    ///   sbnd_xin/docs/pr/90 sec 3b: 4 pts / 1.94 cm of a 35 cm baseline
    ///   outscored the true 33 deg corner), but a hard filter is wrong too --
    ///   genuine corners do sit 4-5 cm from an end (pr/90 sec 8.6), so the
    ///   starved candidate stays as the fallback.  0 = legacy argmax,
    ///   byte-identical.
    /// - bragg_veto_turn (deg): doc sbnd_xin/docs/pr/90 sec 9.4b -- the
    ///   owner-calibrated keep/kill rule for accepted R2 breaks.  A genuine
    ///   back-to-back junction turns hard (all 5 owner-KEEP events >= 32.5
    ///   deg) while spurious accepts hug the 25 deg threshold (all 4
    ///   owner-KILL events 26.5-27.4 deg), and their end brightness is
    ///   vertex activity / track overlap, not a Bragg (dim-extended; a
    ///   genuine Bragg concentrates its charge).  An accepted R2 break with
    ///   turn < bragg_veto_turn whose SHORT-arm end is not Bragg-consistent
    ///   (peak >= 2.0 x mip_dqdx_median AND contiguous > 1.5x hot extent
    ///   from that end <= (peak - 1) x 1 cm/MIP, over an 8 cm end window)
    ///   is vetoed.  Route R1 (dip) accepts are untouched.  <= 0 = off,
    ///   byte-identical.  Calibrated operating point 30.0 (margins -8.7% /
    ///   +8.3% against the 11 owner verdicts; the -1 cm extent offset is the
    ///   sec 10.3 in-code recalibration that keeps 64503 vetoed).
    /// - r3_turn (deg) / r3_hot (x mip_dqdx_median): options for
    ///   segment_chain_turn_break_scan (route R3, doc pr/90 sec 9.5 D3);
    ///   unused by segment_two_end_break_scan itself.  Both > 0 enables R3.
    struct TwoEndBreakOptions {
        double mip_dqdx{50000/units::cm};
        double mip_dqdx_median{43000/units::cm};
        double min_len{10*units::cm};
        double min_arm{1.8*units::cm};
        int    min_arm_pts{4};
        double accept_range{15*units::cm};
        double rise_r1{1.3};
        double rise_r2{1.15};
        double abs_end_min{1.7};
        double dip_floor{0.6};
        double score_cap_r1{0.6};
        double score_cap_r2{0.9};
        double turn_angle{25.0};
        double turn_baseline{35*units::cm};
        double turn_skirt{3*units::cm};
        double turn_min_arm_frac{0.0};
        double bragg_veto_turn{0.0};
        double r3_turn{0.0};
        double r3_hot{0.0};
    };

    /// Result of segment_two_end_break_scan.  `found` is the overall accept
    /// (route1 || route2); `break_idx` is the accepted fit index (route R1's
    /// dip index, else route R2's turn-maximum index; -1 when nothing
    /// accepted).  idx_dip / idx_turn are the two candidate indices (-1 when
    /// ineligible / not computed); turn_deg is the turn at idx_turn.  The
    /// remaining members expose every intermediate the routes consulted, for
    /// logging and doctests.
    struct TwoEndBreakResult {
        bool   found{false};
        bool   route1{false};
        bool   route2{false};
        /// Route R3 (segment_chain_turn_break_scan) accept; exclusive with
        /// route1/route2 (a scan call runs either R1/R2 or R3, never both).
        bool   route3{false};
        int    break_idx{-1};
        int    idx_dip{-1};
        int    idx_turn{-1};
        double joint_score{1e9};
        double arm_a_len{0}, arm_b_len{0};
        double sA{1e9}, sB{1e9};
        bool   flagA{false}, flagB{false};
        double ratio_lo{0}, ratio_hi{0};
        double absmed_lo{0}, absmed_hi{0};
        double turn_deg{0};
        /// R2 bragg-veto diagnostics (doc pr/90 sec 9.4b): set when the veto
        /// evaluated an accepted R2 candidate.  veto_peak / veto_extent are
        /// the short-arm end peak (x mip_dqdx_median) and contiguous hot
        /// extent; bragg_vetoed means the accept was suppressed.
        bool   bragg_vetoed{false};
        double veto_peak{0}, veto_extent{0};
        /// Every candidate accept attempt, in evaluation order (R1's dips
        /// deepest-first, then R2's turn index if tried): fit index, 3-pt
        /// mean dQ/dx at it (/mip for the caller), per-arm scores/flags.
        /// Bounded by the local-minimum count; for the caller's debug log.
        struct Attempt { int idx; double m3, sA, sB; bool fA, fB, accepted; };
        std::vector<Attempt> attempts;
    };

    /// doc sbnd_xin/docs/pr/48 sec 6 -- the two-end residual-range
    /// back-to-back junction scan.  A single particle has at most one Bragg
    /// (stopping) end; dQ/dx rising at BOTH ends of one segment implies a
    /// junction at the interior dip (18255-51513/56211/57903/57485: the
    /// neutrino vertex mid-segment on a single unbroken fitted track, no
    /// angular kink).  Pre-gate: the length-adaptive two-end rise.
    /// Localization: route R1's candidate is the DEEPEST interior local
    /// minimum of the 3-point-mean dQ/dx at or above dip_floor (the
    /// junction dip every measured event shows; the floor vetoes
    /// instrumental dead-region dips); route R2's candidate is the interior
    /// wide-baseline PCA turn maximum, used only when R1 does not accept
    /// (the turn maximum can sit at an endpoint bend artifact, so it never
    /// overrides a dip accept).
    /// Acceptance at a candidate k: arm A = fits[0..k] hypothesized
    /// stopping at the segment START (reversed vectors) and arm B =
    /// fits[k..N-1] stopping at the END, each scored with do_track_comp
    /// against the muon/proton stopping templates over accept_range with
    /// min(muon, proton) scores under the route's cap.  Route R1 (two-end
    /// rise + absolute end floor) requires the stopping-vs-flat direction
    /// flag on BOTH arms; route R2 (turn threshold + looser rise) requires
    /// it on at least one arm -- eval_ks_ratio's margins reject a genuine
    /// but weak Bragg, and R2's primary evidence is the turn.  R1 wins when
    /// both accept.  Pure measurement: never
    /// modifies the segment; the caller owns the break (break_segment) and
    /// the protect flag.  Deterministic: index-ordered scans, first
    /// minimum/maximum wins.
    TwoEndBreakResult segment_two_end_break_scan(
        SegmentPtr seg, const Clus::ParticleDataSet::pointer& particle_data,
        const TwoEndBreakOptions& opt = {});

    /// doc sbnd_xin/docs/pr/90 sec 9.5 D3 -- route R3: the turn+activity
    /// junction scan for CHAIN-admitted candidates (the owner's "still a
    /// line, no 3-track vertex" class: muon->Michel 172832, 2-track 61681).
    /// Their junction signature is the opposite of the two-Bragg valley R1/R2
    /// look for: a BRIGHT vertex-activity spot (measured 2.50x / 3.1x MIP at
    /// the clicks) coincident with a LOCAL 10 cm-baseline turn (t10 plateau
    /// 19-23.5 deg / climb to 54 deg, vs a 5-7 deg mid-track baseline), while
    /// the production 35 cm turn averages the corner away (18.3 deg < 25) and
    /// the deepest dip is an ordinary MIP fluctuation (0.77x).  Candidate
    /// indices: arm_ok (min_arm/min_arm_pts) AND
    /// t10 = segment_wide_turn_angle(fits, k, turn_skirt, 10 cm) >= r3_turn
    /// AND max dQ/dx within +-2 cm arclength >= r3_hot x mip_dqdx_median.
    /// Winner: largest t10 (ties: smaller index) with a WELL-FORMED
    /// preference tier (both t10 windows fully achievable) ahead of the
    /// unrestricted tier -- a starved near-end t10 window reads jitter AND
    /// can sit on genuinely bright EM charge (172832's Michel end), so the
    /// activity corroboration alone cannot guard it; refined to the +-2 cm
    /// activity argmax (ties: smaller index) subject to arm_ok.
    /// Returns route3/found/break_idx/turn_deg (t10 at the winning index)
    /// and arm lengths; every other member stays at its default.  Pure
    /// measurement, deterministic.  Never call with r3_turn or r3_hot <= 0.
    TwoEndBreakResult segment_chain_turn_break_scan(
        SegmentPtr seg, const TwoEndBreakOptions& opt);

    /// Calculate track length from segment
    ///
    /// If flag == 1 and segment has fitted dx values, sum the dx values from fits.
    /// If flag == 0, calculate geometric length from wcpts.
    ///
    /// @param seg The segment to calculate length for
    /// @param flag Calculation method: 0=geometric from points, 1=from fit dx values
    /// @return Track length
    double segment_track_length(SegmentPtr seg, int flag = 0, int n1 = -1, int n2 = -1, WireCell::Vector dir_perp = WireCell::Vector(0,0,0));
    double segment_track_direct_length(SegmentPtr seg, int n1 = -1, int n2 = -1, WireCell::Vector dir_perp = WireCell::Vector(0,0,0));
    double segment_track_max_deviation(SegmentPtr seg, int n1 = -1, int n2 = -1);
    /// Calculate track length above dQ/dx threshold
    ///
    /// Extracts dQ and dx from segment's fits and calculates length above threshold.
    ///
    /// @param seg The segment containing fit data
    /// @param threshold dQ/dx threshold value
    /// @return Length of track segments above threshold
    double segment_track_length_threshold(SegmentPtr seg, double threshold = 75000./units::cm);
    /// Calculate track length from segment using geometric distance between points
    ///
    /// This is a convenience function that always uses geometric calculation
    /// regardless of available dx data.
    ///
    /// @param seg The segment to calculate length for
    /// @return Geometric track length
    double segment_geometric_length(SegmentPtr seg, int n1 = -1, int n2 = -1, WireCell::Vector dir_perp = WireCell::Vector(0,0,0));


    /// Calculate median dQ/dx for a segment
    ///
    /// Extracts dQ and dx from segment's fits and calculates median dQ/dx.
    ///
    /// @param seg The segment containing fit data
    /// @return Median dQ/dx value (0 if no valid fits)
    double segment_median_dQ_dx(SegmentPtr seg, int n1 = -1, int n2 = -1);
    double segment_rms_dQ_dx(SegmentPtr seg);

    /// doc sbnd_xin/docs/pr/40 -- shared evidence check for the F2/F3 knobs
    /// (shower_reclass_dqdx_guard, shower_topo_dqdx_guard).  Reuses the SAME
    /// median-dQ/dx thresholds segment_determine_dir_track's short-track
    /// fallback already trusts (PRSegmentFunctions.cxx, medium_dQ_dx >
    /// MIP_dQdx*1.75 => proton, < MIP_dQdx*1.2 => muon) to decide whether a
    /// segment's OWN charge profile is decisively non-electron-like, i.e. it
    /// should be spared from an unconditional track-to-electron
    /// reclassification.  A zero/absent median (no valid dQ/dx samples) never
    /// spares -- that is "no evidence", not "MIP-like evidence".
    bool segment_dqdx_spares_electron_reclass(SegmentPtr seg, double MIP_dQdx);

    /// doc sbnd_xin/docs/pr/40 round 10 -- sibling spare-test for
    /// examine_all_showers' cluster-wide reclassification
    /// (NeutrinoTrackShowerSep.cxx ~2070-2105), complementing
    /// segment_dqdx_spares_electron_reclass's flat median-dQ/dx ratio test
    /// (which only spares ratio>1.75 or ratio<1.2, leaving a gap at
    /// [1.2,1.75] where a real disconnected muon fragment near
    /// end-of-range can sit -- SBND 18255-314507 seg 17002, xMIP 1.57x,
    /// median-dQ/dx-impure but confidently muon by shape).  True iff `seg`
    /// is longer than 20 cm AND carries a real (<1.0, not the 100
    /// "unscored" sentinel) particle_score from the Bragg/dE-dx-template
    /// PID (segment_determine_dir_track / segment_do_track_pid,
    /// PRSegmentFunctions.cxx).  Above 20 cm that PID's own template
    /// competition never considers electron (~2540), so a real good score
    /// there is unambiguous muon-or-proton evidence.  Deliberately does
    /// NOT read particle_info()->pdg(): particle_score is left untouched
    /// by examine_all_showers' own reclassification (which only rewrites
    /// particle_info), so a stale-but-still-valid Bragg score survives
    /// independent of whatever pdg a caller may already have overwritten
    /// -- checking pdg here would make the test see its own target's
    /// prior verdict instead of the PID's.
    bool segment_bragg_spares_electron_reclass(SegmentPtr seg);

    /// doc sbnd_xin/docs/pr/74 round 2 P1 -- veto for examine_direction's
    /// flag_shower_in cascade (the |pdg|==13/pdg==0 electron relabel).
    /// True iff `seg` is BOTH long (track length > max_len) AND MIP-like
    /// (median dQ/dx < mip_hi * mip_dqdx_median).  A zero/absent median (no
    /// valid dQ/dx samples) never vetoes -- "no evidence" is not "MIP-like
    /// evidence", same convention as segment_dqdx_spares_electron_reclass.
    bool segment_shower_in_cascade_vetoed(SegmentPtr seg, double mip_dqdx_median,
                                          double max_len, double mip_hi);

    /// doc sbnd_xin/docs/pr/74 round 2 P2 -- total track length reachable
    /// from `start_vtx` without traversing `stem`.  Used by the F14 Michel
    /// rescue's Michel-terminal check: a genuine Michel electron is a few cm
    /// of terminal track, a shower trunk heads a large downstream tree.
    /// The sum is traversal-order independent (visited sets are only
    /// tested/inserted, never iterated) and early-exits once `cap` is
    /// exceeded -- callers only need "over cap".
    double segment_far_subtree_track_length(Graph& graph, VertexPtr start_vtx,
                                            SegmentPtr stem, double cap);

    /// doc sbnd_xin/docs/pr/74 round 4 -- the opening angle, in DEGREES,
    /// between two segments that meet at `shared_pt`, measured over the
    /// leading `dis_cut` of each with segment_cal_dir_3vector (the same 15 cm
    /// tangent idiom the daughter-shower angle tests use,
    /// NeutrinoVertexFinder.cxx:1449/1457).
    ///
    /// Both tangents point AWAY from the shared point, so a straight
    /// continuation gives ~180 deg between them; this function returns the
    /// KINK, `180 - angle`, so straight = 0 and a hard turn-back = 180.
    /// Returns -1 when either tangent is undefined (no fits, degenerate
    /// geometry) -- "unmeasurable", which callers must not read as "straight".
    double segment_pair_kink_deg(SegmentPtr a, SegmentPtr b,
                                 const WireCell::Point& shared_pt, double dis_cut);

    /// doc sbnd_xin/docs/pr/74 round 2 -- env-gated (WCT_SHOWER_TOPO_DEBUG)
    /// transition record for kShowerTopology.  Round 1 could not attribute
    /// 90055's "second demotion" because the flag can be cleared on paths
    /// that write no pdg (invisible to WCT_PID_WRITE_DEBUG); this probe is
    /// called at EVERY set/unset site.  Default off => no output, no
    /// behavior change.
    void pr74_probe_topo_flag(SegmentPtr seg, const char* what, const char* site);

    /// doc sbnd_xin/docs/pr/40 round 2 F5 -- an electron cannot father a proton.
    ///
    /// True iff `seg` emanates from `main_vertex` (graph identity at either
    /// endpoint, not a distance cut -- every measured case sits at exactly
    /// d=0.00 cm) and its FAR endpoint's out-edges include a segment that is
    /// (a) already PID'd proton (pdg 2212, i.e. has_particle_info() and
    /// pdg()==2212) and (b) independently charge-confirmed by its own
    /// median dQ/dx (> 1.75x MIP_dQdx, the same threshold
    /// segment_dqdx_spares_electron_reclass and segment_determine_dir_track's
    /// own short-track fallback already trust).  `main_vertex` may be null
    /// (e.g. before stage 4 determines one), in which case this always
    /// returns false -- there is no vertex to require the segment emanate
    /// from.  Callers own the `flag_shower` check; this only judges the
    /// proton-daughter topology.
    bool segment_has_proton_daughter(Graph& graph, SegmentPtr seg, VertexPtr main_vertex, double MIP_dQdx);

    /// doc sbnd_xin/docs/pr/43 F1 -- multi-hop generalization of the
    /// muon-candidate loop's 1-hop `n_proton` check
    /// (NeutrinoVertexFinder.cxx examine_direction) and of
    /// segment_has_proton_daughter above: a charge-confirmed proton (pdg
    /// 2212, median dQ/dx > 1.75x MIP_dQdx) may sit beyond `seg`'s
    /// immediate far vertex, reached through a bounded chain of short,
    /// non-shower, degree-2 continuation segments (pdg 13 / pdg 0 /
    /// already pdg 211).  `near_vertex` is the vertex `seg` emanates FROM
    /// (graph identity, same convention as segment_has_proton_daughter);
    /// the walk starts at seg's other endpoint and stops at the first
    /// vertex that is not a simple one-segment-in/one-segment-out
    /// continuation (a real hadronic multi-prong vertex, a shower
    /// attaching, or a dead end), so it cannot reach an unrelated proton
    /// past genuine vertex activity. `max_hops` bounds how far the chain
    /// is followed. Returns false if `near_vertex`/`seg` is null,
    /// `MIP_dQdx`/`max_hops` is non-positive, or `seg` does not emanate
    /// from `near_vertex`.
    bool segment_chain_has_proton(Graph& graph, SegmentPtr seg, VertexPtr near_vertex, double MIP_dQdx, int max_hops = 3);

    /// doc sbnd_xin/docs/pr/43 F1 -- companion to segment_chain_has_proton:
    /// collects (rather than just detects) the same short, non-shower,
    /// degree-2 continuation chain, for relabeling every segment in a
    /// disqualified muon candidate's own chain to pion, so the candidate's
    /// stub segments do not stay muon once the candidate itself is
    /// demoted. The proton that ends the chain is never included. Returns
    /// segments in walk order, nearest first; empty if `seg` does not
    /// emanate from `near_vertex` or no continuation exists.
    std::vector<SegmentPtr> segment_chain_continuation(Graph& graph, SegmentPtr seg, VertexPtr near_vertex, int max_hops = 3);

    /// doc sbnd_xin/docs/pr/40 round 4 F8 -- a muon cannot terminate in a
    /// multi-proton hadronic vertex.
    ///
    /// True iff `seg` has an endpoint vertex OTHER than `main_vertex` (graph
    /// identity, same convention as segment_has_proton_daughter) whose
    /// out-edges include at least `min_protons` OTHER segments that are (a)
    /// already PID'd proton (pdg 2212) and (b) independently charge-
    /// confirmed by their own median dQ/dx (> 1.75x MIP_dQdx, same threshold
    /// segment_has_proton_daughter uses).  The main-vertex endpoint is
    /// excluded deliberately: a muon and two protons meeting AT the
    /// neutrino vertex is the ordinary, correct numuCC topology.
    /// `main_vertex` may be null, in which case this always returns false.
    bool segment_at_multi_proton_vertex(Graph& graph, SegmentPtr seg, VertexPtr main_vertex, double MIP_dQdx, int min_protons = 2);

    /// doc sbnd_xin/docs/pr/40 round 5 F10/F11 -- shared straightness veto.
    ///
    /// True iff `seg` is long enough and straight enough that it should
    /// never be treated as (the seed of) an EM shower, regardless of what
    /// triggered the shower classification.  Same shape as the toolkit's
    /// own straightness demotion (NeutrinoVertexFinder.cxx:1432-1447, a
    /// byte-faithful port of prototype NeutrinoID_track_shower.h:2042-2054):
    /// length > min_length, and either the direct (straight-line) span is
    /// itself >= min_direct, or the direct/arc-length ratio exceeds
    /// straight_ratio (near-perfectly straight over its whole extent).  A
    /// muon or exiting track satisfies this; a genuine EM shower's bushy
    /// tail keeps its direct length well below its arc length.
    bool segment_is_straight_long_track(SegmentPtr seg, double min_length = 10*units::cm,
                                        double min_direct = 34*units::cm, double straight_ratio = 0.93);

    /// doc sbnd_xin/docs/pr/40 round 9 -- continuation-aware extension of
    /// segment_is_straight_long_track.
    ///
    /// The break-at-point paths in shower clustering (e.g.
    /// shower_clustering_with_nv_from_vertices PATH C) can hand a guard a
    /// sub-10 cm HALF of a straight track that individually fails the
    /// length floor even though the whole object is one straight muon
    /// (SBND 286906: 8.68 cm anchor seg 9002 at a 4.9 deg kink to the
    /// 126.89 cm body seg 9001).  True iff `seg` itself passes
    /// segment_is_straight_long_track, OR a SAME-cluster sibling across
    /// either endpoint vertex is collinear with it (segment_pair_kink_deg
    /// < max_kink_deg over the 15 cm tangent; -1 = unmeasurable is never
    /// collinear) and either (a) the sibling itself passes
    /// segment_is_straight_long_track or (b) the combined two-segment
    /// chain does (summed arc length > 10 cm and far-endpoint-to-far-
    /// endpoint direct span >= 34 cm or > 0.93x the summed arc length --
    /// the same three constants).  Iteration uses sorted_out_edges only
    /// (deterministic).  Pure predicate: no graph or segment mutation.
    bool segment_is_straight_long_track_or_continuation(Graph& graph, SegmentPtr seg,
                                                        double max_kink_deg = 25.0);

    /// Create and associate a DynamicPointCloud with a segment from path points
    ///
    /// @param segment The segment to associate the DynamicPointCloud with
    /// @param path_points Vector of 3D points to process
    /// @param dv Detector volume for wire plane ID determination
    /// @param cloud_name Name for the DynamicPointCloud (default: "main")
    void create_segment_point_cloud(SegmentPtr segment,
                                    const std::vector<geo_point_t>& path_points,
                                    const IDetectorVolumes::pointer& dv,
                                    const std::string& cloud_name = "main",
                                    const std::vector<size_t>& global_indices = {});

    void create_segment_fit_point_cloud(SegmentPtr segment,
                                    const IDetectorVolumes::pointer& dv,
                                    const std::string& cloud_name = "fit");

    /// doc sbnd_xin/docs/pr/85 -- carry a prong through a short connector.
    ///
    /// `stub` connects `at` <-> `anchor`; `prong` is another segment incident
    /// on `at`.  carry_prong_verify checks (read-only) that the prong can be
    /// re-anchored on `anchor` by splicing the stub's wcpts onto its `at`
    /// end: both segments' wcpt chains terminate on the vertices they claim
    /// (0.01 cm), the prong's far vertex is valid and distinct from `anchor`,
    /// and the (far, anchor) edge slot is free (setS edge aliasing --
    /// examine_structure_review.md B.7).  carry_prong_execute performs the
    /// splice (the mvga op3 re-seat mechanics), rebuilds the "main" point
    /// cloud, and rewires the prong to (anchor, far).  The SegmentPtr
    /// survives: fits/flags/particle_info are preserved (fits stale until
    /// the caller's refit).  Callers must verify ALL prongs at `at` before
    /// executing ANY (all-or-nothing: never half-mutate a vertex's arms).
    bool carry_prong_verify(Graph& graph, SegmentPtr prong, SegmentPtr stub,
                            VertexPtr at, VertexPtr anchor);
    bool carry_prong_execute(Graph& graph, SegmentPtr prong, SegmentPtr stub,
                             VertexPtr at, VertexPtr anchor,
                             const IDetectorVolumes::pointer& dv);

    std::pair<double, WireCell::Point> segment_get_closest_point(SegmentPtr seg, const WireCell::Point& point, const std::string& cloud_name = "fit", const std::string& base_cloud_name = "main");
    std::tuple<double, double, double> segment_get_closest_2d_distances(SegmentPtr seg, const WireCell::Point& point, int apa, int face, const std::string& cloud_name = "fit");
    double segment_get_closest_2d_distance(SegmentPtr seg, const WireCell::Point& point, int apa, int face, int plane, const std::string& cloud_name = "fit");


    // PID related
    bool eval_ks_ratio(double ks1, double ks2, double ratio1, double ratio2);
    // skip_stop_samples: exclude that many samples nearest the hypothesized
    // stopping end (largest L, smallest residual range) from the comparison,
    // WITHOUT moving the template anchor (end_L).  0 = legacy, byte-identical.
    // Used by the endpoint-trim retry (doc sbnd_xin/docs/pr/9 sec. 6 F1): the
    // stopping tip's dQ/dx is unreliable when the endpoint is ill-defined
    // (diluted OR piled-up), and it is compared against the template's Bragg
    // maximum, so one bad tip sample can veto an otherwise clean decision.
    // empty_abstain (doc sbnd_xin/docs/pr/31 §12, F6 was P7): when the
    // comparison window holds no samples (or a dEdx template is missing), the
    // degenerate return's element 0 becomes 0.0 ("abstain") instead of 1.0
    // ("this orientation passed the direction gate").  0.0 is the prototype's
    // degenerate answer, verified by execution: zero-bin TH1F KolmogorovTest
    // returns 0 for both templates, so eval_ks_ratio's first gate
    // (ks1-ks2 >= 0.0) fails and results[0] is false.  false = legacy 1.0 =
    // byte-identical.
    std::vector<double> do_track_comp(std::vector<double>& L , std::vector<double>& dQ_dx, double compare_range, double offset_length, const Clus::ParticleDataSet::pointer& particle_data, double MIP_dQdx = 50000/units::cm, int skip_stop_samples = 0, bool empty_abstain = false);

    // Options for segment_do_track_pid / segment_determine_dir_track.
    // Defaults reproduce legacy behavior byte-for-byte (uBooNE values, vote off).
    // - mip_dqdx: the flat-template MIP amplitude handed to do_track_comp
    //   (uBooNE legacy 50e3 e/cm; distinct from the 43e3 median-threshold
    //   scale carried by the functions' MIP_dQdx parameter).
    // - proton_dir_vote: doc sbnd_xin/docs/pr/8 improvement (beyond prototype):
    //   when the muon-vs-flat direction gate abstains in both orientations,
    //   let the proton template declare the direction, guarded by
    //   score/asymmetry/free-end conditions.  Initial thresholds are
    //   placeholders pending the pr/8 sec. 6 calibration.
    // - start_n/end_n: vertex connectivity, filled by
    //   segment_determine_dir_track (guard G4: stop end must be free).
    // - endpoint_trim_retry: doc sbnd_xin/docs/pr/9 sec. 6 F1 (beyond
    //   prototype, owner-approved 2026-07-30): when the legacy decision
    //   abstains in both orientations, retry ONCE with exactly one sample
    //   excluded at each orientation's hypothesized stopping end (template
    //   anchor unchanged).  Dynamic by construction: a well-found endpoint
    //   means the first attempt already decides and nothing is trimmed;
    //   never trims more than 1 sample; value-agnostic (robust to diluted
    //   AND piled-up tips).  Runs BEFORE the proton_dir_vote fallback.
    // - track_comp_empty_abstain: doc sbnd_xin/docs/pr/31 §12 (F6, was P7).
    //   Forwarded to do_track_comp's empty_abstain (see its comment above).
    //   false = legacy "confirmed" filler = byte-identical.
    // - dir_track_median_local: doc sbnd_xin/docs/pr/31 §12 (F4, was P8).
    //   segment_determine_dir_track takes its median dQ/dx over the SAME local
    //   dQ_dx vector it hands to the PID (prototype ProtoSegment.cxx:1574-1576,
    //   :1602-1611: nth_element over a copy of dQ_dx), instead of
    //   segment_median_dQ_dx's filtered rebuild from fits() -- which drops
    //   invalid/dx<=0/dQ<0 samples the PID vector keeps as zeros, so the two
    //   disagree about the same pathological point in opposite directions and
    //   nth_element selects a different order statistic.  Same internal-unit
    //   scale either way (charge per internal length; see the unit-convention
    //   comment in segment_determine_dir_track).  false = filtered helper =
    //   byte-identical.
    // - track_pid_persist_dqdx: doc sbnd_xin/docs/pr/40 (F1, = doc pr/7 sec 5 /
    //   pr/31 P14/F8).  segment_determine_dir_track's final store
    //   (`if (pdg_code != 0 && ((dirsign==1 && end_n==1) || (dirsign==-1 &&
    //   start_n==1)))`) gates BOTH the type/mass persistence and the 4-mom
    //   recompute on the direction pointing at a free end -- so a dQ/dx-only
    //   recovery (medium_dQ_dx > 1.75x/< 1.2x MIP with dirsign left at 0, or a
    //   confident template PID whose stop end is not topologically free) is
    //   computed and then silently discarded; the segment exits with no
    //   particle_info at all.  The prototype (ProtoSegment.cxx:1637-1639)
    //   persists type+mass unconditionally and gates ONLY cal_4mom() on
    //   get_particle_4mom(3)>0 ("an energy was already computed").  true =
    //   restore that shape: store type+mass whenever pdg_code != 0, gate only
    //   the 4-momentum recompute on the existing free-end test.  false =
    //   legacy = byte-identical.  Measured (doc pr/40 Part 0): of 9 owner-
    //   reported track-to-electron cases, this persistence gate is why 4 of
    //   them (evt 74544, 174637, 267597, 269774) had NO particle_info at all
    //   by the time a wholesale shower-reclassification site looked at them.
    // - track_pid_persist_4mom: doc sbnd_xin/docs/pr/40 round 2 (F4).  When
    //   track_pid_persist_dqdx rescues a non-free-end store, the legacy
    //   4-momentum stub is D4Vector(mass,0,0,0) -- E = mass exactly, so
    //   Aux::ParticleInfo::kinetic_energy() (= E - mass) is ZERO for that
    //   segment, byte-for-byte, forever (SBND evt 174637 seg 9050 measured
    //   as "mu- 0 MeV" in the Bee PF tree once track_pid_persist_dqdx went
    //   SBND-default-ON).  segment_cal_4mom has NO dependence on the
    //   direction being a free end -- its only direction coupling is
    //   segment_cal_dir_3vector, which already returns a zero 3-vector when
    //   dirsign()==0 -- so the free-end gate on it was purely external and
    //   unnecessarily strict.  true = call segment_cal_4mom unconditionally
    //   whenever pdg_code != 0 (correct KE, momentum direction only as good
    //   as dirsign allows).  false = legacy rest-mass-only stub =
    //   byte-identical.  Independent of track_pid_persist_dqdx.
    struct TrackPidOptions {
        double mip_dqdx{50000/units::cm};
        bool   proton_dir_vote{false};
        double proton_dir_score_max{0.25};
        double proton_dir_asym_min{1.3};
        bool   endpoint_trim_retry{false};
        int    start_n{1};
        int    end_n{1};
        bool   track_comp_empty_abstain{false};
        bool   dir_track_median_local{false};
        bool   track_pid_persist_dqdx{false};
        bool   track_pid_persist_4mom{false};
        // doc sbnd_xin/docs/pr/40 round 5 F9 -- narrows track_pid_persist_dqdx
        // (F1).  When true, F1 no longer rescues a pdg_code==11 conclusion
        // that lacks a free end (segment_determine_dir_track's free_end_dir
        // false) -- it still rescues every non-electron pdg exactly as
        // before.  See NeutrinoPatternBase.h's m_track_pid_persist_dqdx_
        // electron_guard comment for the measured case (SBND evt 84229) and
        // porting_dictionary.md for the designed-divergence status.  false =
        // today's shipped F1 behaviour = byte-identical.
        bool   track_pid_persist_dqdx_electron_guard{false};
    };

    // success, flag_dir, pdg_code, particle_score
    std::tuple<bool, int, int, double> segment_do_track_pid(SegmentPtr segment, std::vector<double>& L , std::vector<double>& dQ_dx, const Clus::ParticleDataSet::pointer& particle_data, double compare_range=35*units::cm, double offset_length = 0*units::cm, bool flag_force = false,  const TrackPidOptions& pid_opts = {});

    // direction calculation ...
    /// Returns true if the segment's direction is weakly determined.
    ///
    /// Mirrors prototype ProtoSegment::is_dir_weak(): checks particle-type-based
    /// score thresholds (muon, proton) and the static dir_weak flag.
    bool segment_is_dir_weak(SegmentPtr seg);
    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg);
    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg, WireCell::Point& p, double dis_cut);
    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg, int direction, int num_points, int start);
    void segment_determine_dir_track(SegmentPtr segment, int start_n, int end_n, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx = 43000/units::cm, bool flag_print = false, const TrackPidOptions& pid_opts = {});

    // kinemiatics calculations ...
    double segment_cal_kine_dQdx(SegmentPtr seg, const IRecombinationModel::pointer& recomb_model);
    double cal_kine_dQdx(std::vector<double>& vec_dQ, std::vector<double>& vec_dx, const IRecombinationModel::pointer& recomb_model);
    double cal_kine_range(double L, int pdg_code, const Clus::ParticleDataSet::pointer& particle_data);
    // 4-momentum: E, px, py, pz
    WireCell::D4Vector<double> segment_cal_4mom(SegmentPtr segment, int pdg_code, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx = 50000/units::cm);

    // EMshower PID
    //
    // doc pr/64 round 7: reassign_orphans (no-op unless true) re-examines
    // every point Stage C would otherwise DROP -- either because no Voronoi
    // terminal reaches it ("channel a") or because it lost the ghost-removal
    // 2D-projection contest to a segment that never claims it ("channel b1")
    // -- and hands it to whichever OTHER segment in the SAME cluster actually
    // achieves the global 2D minimum, using a duplicated copy of the Stage-C
    // acceptance rule (M10: the primary rule chain is untouched). A point
    // whose true 2D winner is in a DIFFERENT cluster is still dropped, so
    // cross-cluster ghost rejection is unaffected by construction. Default
    // false keeps this function byte-for-byte identical to today's behavior;
    // see WCT_PR64_ORPHAN_CENSUS for the log-only (no-op) diagnostic that
    // motivated this knob, and PRSegmentFunctions.cxx for the full mechanism.
    void clustering_points_segments(std::vector<SegmentPtr> segments, const IDetectorVolumes::pointer& dv, const std::string& cloud_name = "associate_points", double search_range = 1.2*units::cm, double scaling_2d = 0.7, bool reassign_orphans = false);

    /// doc sbnd_xin/docs/pr/32 §11 F2 -- knob transport for
    /// segment_is_shower_trajectory's flag semantics.
    ///
    /// The prototype's ProtoSegment::is_shower_trajectory() opens with
    /// `flag_shower_trajectory = false` (ProtoSegment.cxx:544) and sets it true
    /// only if the test fires (:608), so EVERY call re-caches the label and a
    /// segment that no longer qualifies is DEMOTED.  The toolkit's port only
    /// ever calls set_flags (PRSegmentFunctions.cxx), so kShowerTrajectory is
    /// monotone: once set it survives every later negative test.  That is why
    /// the two improve_vertex call sites recompute instead of reading the flag
    /// (doc pr/32 §10.3) -- reading it would return an OR over history.
    ///
    /// true => mirror the prototype: clear the bit on entry, set it on a
    /// positive result.  false => today's set-only behaviour, byte-identical.
    ///
    /// This is a free function reached from three files with no component
    /// configuration in scope, so the value travels through a process-wide
    /// flag written once per event by TaggerCheckNeutrino before any graph is
    /// built -- the same transport as PR::g_graph_endpoint_policy (doc pr/30
    /// §11 P8).  Read-mostly; never written concurrently with clustering.
    extern bool g_shower_traj_refresh_flag;

    // straight_guard (doc sbnd_xin/docs/pr/40 round 5 F11): after the
    // existing 5-section wiggliness test decides flag_shower_trajectory,
    // override it back to false if segment_is_straight_long_track(seg) --
    // see that function's comment and m_shower_traj_straight_guard's
    // comment in NeutrinoPatternBase.h.  false = legacy = byte-identical.
    bool segment_is_shower_trajectory(SegmentPtr seg, double step_size = 10*units::cm, double mip_dQ_dx = 50000 / units::cm, bool straight_guard = false);
    void segment_determine_shower_direction_trajectory(SegmentPtr segment, int start_n, int end_n, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx = 43000/units::cm, bool flag_print = false, const TrackPidOptions& pid_opts = {});
    
    // median_local (doc sbnd_xin/docs/pr/31 §12, F4): forwarded to the interior
    // segment_determine_dir_track call on the short-segment branch as
    // TrackPidOptions::dir_track_median_local ONLY -- that interior call
    // otherwise passes a default-constructed TrackPidOptions, and forwarding a
    // caller's full option set there would unconditionally change proton_dir_vote
    // et al. at that site.  false = legacy = byte-identical.
    bool segment_determine_shower_direction(SegmentPtr segment, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, const std::string& cloud_name = "associate_points", double MIP_dQdx = 43000/units::cm, double rms_cut= 0.4*units::cm, double mip_dqdx = 50000/units::cm, bool median_local = false);
    // sbnd_xin doc pr/25 sec 3: `demote_len` > 0 makes the existing long-track
    // guard unconditional -- any segment this function would flag
    // kShowerTopology whose geometric length (segment_track_length(seg,0),
    // the same measure the legacy >50cm guard uses) exceeds `demote_len` is
    // demoted to a track instead.  Default 0 == off == byte-identical.
    // Motivation: an owner hand-scan of every long shower-topology segment on
    // a selected nu-candidate main cluster in the 572-event valfast manifest
    // (2026-08-03, 10/10 events) found NO showers -- all tracks.  See sec 3.8.
    // reset (doc sbnd_xin/docs/pr/31 §12, F3 was P13): the prototype's
    // is_shower_topology opens with two assignments BEFORE its early returns
    // (ProtoSegment.cxx:319-321): flag_shower_topology = tmp_val (all callers
    // pass false => the flag is CLEARED on every entry) and flag_dir = 0.  The
    // toolkit's port assigns a local and never clears the segment flag, and its
    // four early returns skip dirsign() -- so a segment that qualified on an
    // earlier pass keeps a stale kShowerTopology flag and a stale direction.
    // true => clear the flag bit (unset_flags, not clear_flags -- other bits
    // survive) and zero dirsign at entry, mirroring the prototype; the tail
    // re-sets both when the answer is still yes.  false = set-only legacy =
    // byte-identical.
    // dqdx_guard (doc sbnd_xin/docs/pr/40, F3): the local vec_dQ_dx this
    // function already builds (fits[i].dQ/dx, normalized by MIP_dQ_dx) is
    // otherwise DEAD -- read only as .size() (pr/31 GOTCHA 5) -- so the
    // 5-branch geometric spread test that decides flag_shower_topology never
    // once consults the segment's own charge.  true = after that test (and
    // the existing length-based demotions) decide flag_shower_topology=true,
    // override it back to false when segment_dqdx_spares_electron_reclass
    // says the segment's median dQ/dx is decisively proton- or muon-like.
    // This does NOT touch flag_dir or the geometric test itself -- a
    // genuinely isochronous/EM-like spread pattern with ambiguous charge
    // still sets the flag.  false = legacy = byte-identical.  Measured (doc
    // pr/40 Part 0): rescues both owner-reported cases where the topology
    // test itself, not a pdg-11 write elsewhere, was the mechanism (evt
    // 256587 seg 11079 at 1.26x MIP, evt 489330 seg 4018 at 2.73x MIP -- both
    // shorter than the existing shower_topo_demote_len's reach).
    bool segment_is_shower_topology(SegmentPtr seg, bool tmp_val=false, double MIP_dQ_dx = 43000/units::cm,
                                    double demote_len = 0, bool reset = false, bool dqdx_guard = false);

    // doc sbnd_xin/docs/pr/72 round 2 -- examine_structure_3
    // (NeutrinoStructureExaminer.cxx) merges any degree-2 junction whose
    // bulk (10cm) and local (3cm) direction agreement both clear lenient
    // thresholds (18deg/27deg), with no check for whether the junction is a
    // genuine near-vertex stub meeting a shower trunk rather than one
    // particle's trajectory the tracker happened to split (18255-196649:
    // a real 6.28cm track stub, terminal at the true neutrino vertex,
    // silently absorbed into a 33cm shower trunk).  These parameters and
    // the pure predicate below are fitted from a 117-event
    // (48 nueCC + 19 NC-pi0 + 50 PR-data) merge census: at these defaults
    // exactly one event in that sample (196649 itself) is touched, zero
    // residual "suspicious but unsuppressed" merges remain -- see
    // sbnd_xin/docs/pr/72 round 2 for the grid scan and the near-miss
    // events these thresholds deliberately do NOT catch (evt235435,
    // evt423981, evt506746: same length/angle regime but a NON-terminal
    // short-arm far end, i.e. the short arm chains further into the graph
    // rather than ending at a free vertex candidate).
    struct Es3StubGuardParams {
        double stub_max{7 * units::cm};    // internal units; short-arm length ceiling
        double len_ratio{2.0};             // long/short length ratio floor
        double ang3_min{15.0};             // degrees; local-kink floor
        double ang_ratio{1.0};             // require ang3 > ang_ratio * ang10
        bool   require_terminal{true};     // far end of the short arm must be a free vertex (degree 1)
    };

    // Pure predicate, no graph/segment access -- every input is a value
    // already computed by examine_structure_3's merge branch (or its
    // WCT_ES3_MERGE_CENSUS census), so this is trivially unit-testable and
    // trivially proven inert when the caller doesn't invoke it.  Caller
    // determines len_short/len_long/deg_short by comparing len1 vs len2
    // itself (edge-index order sg1/sg2 has nothing to do with length).
    // Returns false (never suppresses) whenever either arm has fewer than 2
    // fit points -- segment_track_length returns 0.0 for such a degenerate
    // arm, which would otherwise make it "the short arm" by construction
    // and make len_long/len_short divide by (near-)zero.
    bool es3_stub_suppress(double len_short, double len_long, double ang3, double ang10,
                            int deg_short, int nfit_short, int nfit_long,
                            const Es3StubGuardParams& params = {});
}

#endif
