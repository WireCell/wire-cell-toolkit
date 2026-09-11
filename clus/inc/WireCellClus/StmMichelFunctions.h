/** Free functions behind the CheckSTM_Michel visitor (doc pdvd/48).

    Everything here operates on a PR::Graph plus scalars, so the muon-chain
    walk, the residual-range profile, the Bragg-contrast metric and the
    Michel / delta-ray / continuation classification can be unit-tested on a
    synthetic graph (clus/test/doctest_stm_michel.cxx) without a detector,
    a fitter or a particle dataset.  The visitor (clus/src/CheckSTM_Michel.cxx)
    supplies the graph the general PR chain built and reads the verdict back.

    Iteration is deterministic by construction: every out-edge walk goes
    through sorted_out_edges() and every set/map is index-ordered.  Nothing
    here iterates a pointer-keyed container.
 */
#ifndef WIRECELL_CLUS_STMMICHELFUNCTIONS
#define WIRECELL_CLUS_STMMICHELFUNCTIONS

#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellIface/IRecombinationModel.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"

#include <array>
#include <functional>
#include <vector>

namespace WireCell::Clus::PR {

    /// The vertex position the walk routes from: the fitted point when the
    /// vertex has one, else the original wcpt (same rule as
    /// find_cont_muon_segment, NeutrinoVertexFinder.cxx:1356).
    WireCell::Point stm_michel_vertex_point(const VertexPtr& vtx);

    /// Segments of the shortest (by segment_track_length) route from `from`
    /// to `to`, in walk order.  Empty when either vertex lacks a descriptor
    /// or no route exists.  Hand-rolled Dijkstra: the graph carries no weight
    /// property map and boost::out_edges on setS is pointer-ordered, so this
    /// walks sorted_out_edges and breaks ties on the edge index.
    std::vector<SegmentPtr> stm_michel_shortest_chain(Graph& g, VertexPtr from, VertexPtr to);

    /// The reachable vertex farthest (by the same segment_track_length
    /// metric) from `from`, optionally restricted to vertices `accept`
    /// admits (e.g. those stamped with the main cluster).  Ties break on
    /// the vertex index.  Null when nothing but `from` is reachable.  Used
    /// when the tagger's stop has no graph vertex (doc pdhd/03: the tagger's
    /// fit bridged into a detached fragment) -- the muon is then the longest
    /// route out of the entry, and the verdict keeps R_STOP_UNMATCHED.
    VertexPtr stm_michel_farthest_vertex(Graph& g, VertexPtr from,
                                         const std::function<bool(const VertexPtr&)>& accept = nullptr);

    /// doc pdvd/88 (doc pdvd/78 action item 7b): every vertex reachable from
    /// `from` along graph edges, `from` included, sorted by graph index.
    /// Empty when `from` is null or has no descriptor.
    std::vector<VertexPtr> stm_michel_reachable_vertices(Graph& g, VertexPtr from);

    /// doc pdvd/88: the vertex of `cands` nearest `pt` (by its wcpt, the
    /// point PatternAlgorithms::closest_cluster_vertex measures) among those
    /// `accept` admits; ties keep the earlier candidate.  {nullptr, 1e9} when
    /// none qualifies -- closest_cluster_vertex's own "none" value.
    std::pair<VertexPtr, double> stm_michel_closest_vertex_of(const std::vector<VertexPtr>& cands,
                                                              const WireCell::Point& pt,
                                                              const std::function<bool(const VertexPtr&)>& accept = nullptr);

    /// The chain's vertices in walk order: entry, every junction, the far
    /// vertex of the last segment.  Size = chain.size() + 1 (0 if chain empty
    /// or a segment is not attached to the running vertex).
    std::vector<VertexPtr> stm_michel_chain_vertices(Graph& g, const std::vector<SegmentPtr>& chain, VertexPtr entry);

    /// dQ/dx vs residual range along a chain, walked from `entry`.
    ///
    /// Each segment's fits are taken in the direction that leaves the
    /// running vertex -- decided by which of fits.front()/fits.back() is
    /// nearer that vertex, NEVER by dirsign() (unset until examine_direction
    /// runs, and re-oriented by it).  Fits with dx <= 0 or dQ < 0 are dropped
    /// (the segment_median_dQ_dx filter).  L is the cumulative point-to-point
    /// arclength from the entry, rr = L_total - L.
    ///
    /// UNITS: L, rr, pts in WCT internal length; dQdx in ELECTRONS PER CM --
    /// the STM tagger's frame, in which TaggerCheckSTM::m_mip_dqdx and the
    /// ParticleDataSet *DeDx tables are expressed
    /// (cfg/.../particle_dataset.jsonnet).  do_track_comp wants the internal
    /// frame (divide by units::cm) -- convert at that one call.
    struct StmMichelProfile {
        std::vector<double> L, dQdx, rr;
        std::vector<WireCell::Point> pts;
        std::vector<int> seg_idx;      // index into the chain vector
        double total_length{0};
        bool empty() const { return L.empty(); }
    };
    StmMichelProfile stm_michel_profile(Graph& g, const std::vector<SegmentPtr>& chain, VertexPtr entry);

    /// The profile with every point whose dQdx < min_dqdx removed (L, rr and
    /// total_length keep their geometric values, so residual range is still
    /// measured from the stop).  n_dead receives the number removed.  A fit
    /// point with (near-)zero charge is a cell the fit could not read -- a
    /// dead / unresponsive channel stretch, an APA edge, a wrapped-wire
    /// ambiguity -- and carries no particle information; left in, a 20 cm
    /// zero stretch drags the plateau median to ~0 (a FAKE Bragg contrast of
    /// 5) and hands the template PID to the electron (doc pdhd/03 sec 5).
    StmMichelProfile stm_michel_profile_live(const StmMichelProfile& prof, double min_dqdx, int& n_dead);

    /// Median of a copy of v; 0 when empty.
    double stm_michel_median(std::vector<double> v);

    /// Charge (electrons) -> energy, for a Michel piece the fitter never
    /// reached so there is no dx and hence no dQ/dx to invert (doc pdhd/15).
    ///
    /// Forked BY DUPLICATION (CLAUDE.md M10) from the one line every
    /// charge-based energy in this repo ends on,
    /// NeutrinoEnergyReco.cxx:509 (itself the port of the prototype's
    /// NeutrinoID_energy_reco.h:248):
    ///
    ///     overall / recom_factor / fudge_factor * w_value / 1e6 * units::MeV
    ///
    /// `w_value` is in eV per ion pair (23.6).  The factor pair is the
    /// caller's: KineChargeOptions offers TRACK (0.7 / 0.95) and SHOWER
    /// (0.5 / 0.8) -- doc pdhd/15 sec 6 measures the track pair against the
    /// chain's own segment_cal_kine_dQdx on 47 PDVD / 26 PDHD single-piece
    /// Michels (ratio median 1.19 / 1.01) and the shower pair at 1.98 / 1.67,
    /// so the track pair is what puts an unfitted piece on the same scale as
    /// the fitted ones.  Returns 0 for a non-positive or non-finite input.
    double stm_michel_charge_to_energy(double dQ_electrons, double recom_factor,
                                       double fudge_factor, double w_value_ev);

    /// doc pdvd/51: the capture-gamma RING test, on one companion cluster.
    ///
    /// A mu- that stops in argon is captured far more often than it decays, and
    /// the capture leaves de-excitation gammas of order an MeV.  A gamma is
    /// NEUTRAL: it leaves no track from the stop, travels a couple of Compton
    /// mean free paths and deposits a compact blob at an arbitrary angle.  So
    /// the reconstructible object is "small, detached, same-bundle, past the
    /// stop", and it is a DIFFERENT object from a Michel.
    ///
    /// The ring's inner edge is deliberately the Michel's own admission radius:
    /// inside it the charge is the Michel object's to claim, and this stage
    /// must not touch it -- widening the Michel radius instead would feed the
    /// Michel piece assembly, and PDVD 039252_15 cluster 77 is already the only
    /// object above the 52.8 MeV Michel endpoint.
    ///
    /// Boundaries, pinned here because they are the whole semantics:
    /// `d_stop` must be STRICTLY greater than `inner_cm` (the Michel's test is
    /// `<=`, so the two partitions do not overlap and do not leave a gap) and
    /// at most `outer_cm`; `len` at most `max_len_cm`.  All arguments in the
    /// same length unit.  A non-finite input admits nothing.
    bool stm_michel_stop_gamma_ring(double d_stop, double len,
                                    double inner_cm, double outer_cm, double max_len_cm);

    /// doc pdvd/51: the capture-gamma ACCEPTANCE window, a separate stage from
    /// the ring test above because a companion's energy does not exist until
    /// the fitter has run on it.  Inclusive at both ends; a non-finite energy is
    /// rejected (a NaN passes no gate and fails every one silently -- doc
    /// pdhd/15 sec 7).
    bool stm_michel_stop_gamma_energy(double ke_mev, double lo_mev, double hi_mev);

    /// doc pdvd/53: the companion ADMISSION radius -- the one number that says
    /// which clusters ever reach the fitter, and therefore which pieces can
    /// exist at all for a later stage or a hand scanner to judge.  It is a
    /// maximum over the stages that are switched on, never a stage's own test:
    /// the Michel still applies `michel_cm` at cluster level and the gamma still
    /// applies its ring, so widening this admits nothing to either object.  A
    /// non-finite or negative radius from a disabled stage contributes nothing.
    double stm_michel_admit_radius(double michel_cm,
                                   double gamma_cm, bool gamma_on,
                                   double survey_cm, bool survey_on);

    /// doc pdhd/17 sec 9: the SAME conversion, but with the survival read out
    /// of the recombination model the component is already holding instead of
    /// the hard-coded pair above.
    ///
    ///     E = dQ_electrons * (dedx * dx) / model(dedx * dx, dx)
    ///
    /// i.e. "how many MeV does one collected electron stand for, if this charge
    /// was deposited at `dedx_mev_per_cm`", asked of the model itself.  `dx`
    /// cancels, so any positive value gives the same answer; 1 cm is used.
    ///
    /// This exists because the pair (0.7 / 0.95 = 0.665 survival) and the
    /// component's dQ/dx -> dE/dx inverse are two carriers of ONE physical
    /// quantity, and doc pdhd/16 moved only the second: after it the two sit
    /// 20 % (PDVD) / 16 % (PDHD) apart at MIP, where before they were within
    /// 5 %.  Deriving the first FROM the second cannot drift again, whatever
    /// model is bound.
    ///
    /// The MIP assumption is unavoidable and it is not neutral: an unfitted
    /// cluster has no dx, so there is no dQ/dx to invert and a dE/dx must be
    /// assumed.  Because quenching rises with dE/dx, a deposit DENSER than the
    /// assumption needs MORE MeV per electron than this returns -- at the
    /// calibrated PDHD model, 4.13e-5 MeV/e at 2.1 MeV/cm against 4.96e-5 at
    /// 5 MeV/cm -- so assuming MIP UNDER-estimates a dense deposit by up to
    /// ~20 %.  Everything quoted through this function is MIP-equivalent.
    ///
    /// Returns 0 for a non-positive or non-finite input, a null model, a
    /// non-positive `dedx_mev_per_cm`, or a model that returns a non-positive
    /// charge there (the Modified Box's forward goes negative below its A < 1
    /// zero crossing, u = 1 - A: at p = 1 that is 0.205 MeV/cm on PDVD and
    /// 0.226 on PDHD, far below any MIP assumption).
    double stm_michel_charge_to_energy_model(double dQ_electrons,
                                             const IRecombinationModel::pointer& model,
                                             double dedx_mev_per_cm);

    /// doc pdvd/81: the chain's own three-plane charge combination, forked BY
    /// DUPLICATION (CLAUDE.md M10) from kine_charge_from_maps
    /// (NeutrinoEnergyReco.cxx:149-186, the port of the prototype's
    /// NeutrinoID_energy_reco.h:248 block): the weighted mean
    /// sum(w_p q_p) / sum(w) of the per-plane sums, replaced by the (median,
    /// minimum) pair's weighted mean when the two largest planes disagree by a
    /// relative asymmetry |med - max| / (med + max) above `asym_switch` (the
    /// largest plane is then treated as contaminated).  The asymmetry is only
    /// evaluated when med + max > 0; with an all-zero weight triple the result
    /// is 0.  Written for SIGNED input: a plane sum that came out negative
    /// (measured minus predicted) takes part like any other number, so the
    /// caller converts the result with a charge-to-energy that returns 0 for a
    /// non-positive charge.  `dropped_plane`, when given, receives the plane
    /// index the switch dropped, else -1.  Pure.
    double stm_michel_combine_planes(const std::array<double, 3>& sums,
                                     const std::array<double, 3>& weights,
                                     double asym_switch, int* dropped_plane = nullptr);

    /// The Bragg-contrast metric (doc pdvd/25 sec 13.9 item 3): median dQ/dx
    /// over the tail window rr in [tail_lo, tail_hi] divided by the median
    /// over the plateau window rr in [pl_lo, pl_hi].  `expected` is the SAME
    /// two medians of the injected muon table evaluated at the profile's own
    /// rr values (in cm), so the two ratios share their binning.  A chain
    /// shorter than pl_hi falls back to a plateau window of [pl_lo/2, pl_hi/2]
    /// and sets short_track.  valid is false when either window holds fewer
    /// than 3 points.
    struct StmMichelBragg {
        double tail_med{0}, plateau_med{0}, contrast{0}, expected{0};
        int n_tail{0}, n_plateau{0};
        bool short_track{false};
        bool valid{false};
    };
    StmMichelBragg stm_michel_bragg_contrast(const StmMichelProfile& prof,
                                             const std::function<double(double)>& mu_dqdx_at_rr_cm,
                                             double tail_lo, double tail_hi,
                                             double pl_lo, double pl_hi);

    /// WCP ProtoSegment::get_flag_shower(): trajectory flag OR topology flag
    /// OR |pdg| == 11.  Copy of the file-local helper in
    /// NeutrinoTaggerCosmic.cxx:76-80 (not exported there).
    bool stm_michel_seg_is_shower(const SegmentPtr& seg);

    /// Thresholds for the arm classification.  Lengths in WCT internal
    /// units, angles in degrees, mip_* as multiples of mip_dqdx_median.
    struct StmMichelArmThresholds {
        double mip_dqdx_median{43000 / units::cm};   // internal units (charge per length)
        double dir_window{15 * units::cm};            // tangent window for the kink angle
        // a collinear MIP track leaving the stop = the muon did not stop here
        double continuation_max_angle_deg{20};        // doc pdvd/42 sec 4.4: 83 % of the leftover within 20 deg
        double continuation_min_len{3 * units::cm};   // the leftover's median is 5.8 cm; a 3 cm floor kills fit jitter only
        double continuation_mip_lo{0.7};
        double continuation_mip_hi{1.3};
        // Michel electron at the stop
        double michel_max_len{25 * units::cm};       // own length + far subtree
        double michel_mip_hi{2.0};                    // upper cap: an electron is MIP-like or below
        double michel_mip_lo{0.3};                    // charge-desert guard
        double michel_min_kink_deg{30};               // required unless the arm is shower-flagged
        double michel_shower_min_kink_deg{-1};        // doc pdhd/03: a shower-flagged arm still needs this kink (measurable); -1 = any (doc pdvd/48)
        // doc pdvd/73 (P2): PDVD operating points for the Michel clause.  Each
        // is OFF at its default (-1), and with all three off the classifier
        // runs its doc pdvd/48 expression verbatim.
        double michel_mip_lo_turned{-1};              // (a) a lower charge floor, only for an arm that turns hard; -1 = off
        double michel_mip_lo_turned_kink_deg{60};     // (a) "turns hard"
        double michel_far_len_shower_max{-1};         // (b) a shower-flagged arm's far subtree is judged on its own cap, not in len + far_len; -1 = off
        double michel_kink_window{-1};                // (c) the Michel turn test may also read the kink over this shorter window; -1 = off
        // delta ray off the muon body
        double delta_max_len{8 * units::cm};
        // a long heavily-ionizing prong off the body = hadronic vertex
        double hadron_mip{1.4};
    };

    struct StmMichelArm {
        enum Kind { kOther = 0, kMichel = 1, kContinuation = 2, kDelta = 3, kHadron = 4 };
        Kind kind{kOther};
        SegmentPtr seg;
        double len{0};          // segment_track_length(arm)
        double mip{0};          // median dQ/dx / mip_dqdx_median
        double kink_deg{-1};    // segment_pair_kink_deg vs the incoming muon segment; -1 = unmeasurable
        double far_len{0};      // track length reachable beyond the arm's far vertex (capped)
        bool shower_like{false};
        bool terminal{false};   // far vertex has degree 1
        double kink_w_deg{-1};  // doc pdvd/73 (P2c): the kink over michel_kink_window; -1 = not measured (off, or unmeasurable)
    };

    /// Classify a non-muon arm at the STOP vertex.
    ///  - Continuation: not shower-like, kink < continuation_max_angle_deg
    ///    (measurable), len > continuation_min_len, mip in [lo, hi].
    ///  - Michel: len + far_len <= michel_max_len, michel_mip_lo < mip <
    ///    michel_mip_hi, and (shower_like OR kink >= michel_min_kink_deg).
    ///    Never through dQ/dx alone (see the .cxx).
    ///  - else Other.
    StmMichelArm stm_michel_classify_stop_arm(Graph& g, SegmentPtr last_muon, SegmentPtr arm,
                                              VertexPtr stop, const StmMichelArmThresholds& th);

    /// doc pdvd/73 (P2): the Michel clause of stm_michel_classify_stop_arm as a
    /// pure predicate, used when any P2 field of `th` is on (with all of them
    /// off the classifier keeps its own expression, and this function returns
    /// the same answer -- the doctest pins that on a grid).  An arm that is
    /// not a continuation is a Michel when
    ///  - reach: len + far_len <= michel_max_len; with (b) on, a shower-flagged
    ///    arm instead needs len <= michel_max_len and far_len <=
    ///    michel_far_len_shower_max;
    ///  - charge: michel_mip_lo < mip < michel_mip_hi; with (a) on, the floor
    ///    is michel_mip_lo_turned when kink_deg >= michel_mip_lo_turned_kink_deg;
    ///  - turn: shower-flagged (subject to michel_shower_min_kink_deg), or
    ///    kink_deg >= michel_min_kink_deg; with (c) on, kink_w_deg >=
    ///    michel_min_kink_deg also counts.
    /// kink_deg / kink_w_deg < 0 = unmeasurable.  Lengths in one unit.
    bool stm_michel_michel_gate(double len, double far_len, double mip, double kink_deg, double kink_w_deg,
                                bool shower_like, const StmMichelArmThresholds& th);

    /// doc pdvd/73 (P2b): the track length reachable from `far_vtx` without
    /// crossing `stem` and without stepping back into `stop_vtx` -- the
    /// segment_far_subtree_track_length walk (PRSegmentFunctions.h) with the
    /// stop vertex fenced off, so a loop back to the stop never adds the muon
    /// chain behind it.  Returns as soon as the total exceeds `cap` (the value
    /// is then a lower bound above the cap).  Deterministic: out-edges are
    /// walked in sorted order.
    double stm_michel_far_subtree_len(Graph& g, VertexPtr far_vtx, SegmentPtr stem, VertexPtr stop_vtx,
                                      double cap);

    /// Classify a non-chain arm at an INTERIOR chain vertex.
    ///  - Delta: len <= delta_max_len, far_len <= delta_max_len, terminal.
    ///  - Hadron: longer than delta_max_len and mip > hadron_mip.
    ///  - else Other.
    StmMichelArm stm_michel_classify_chain_arm(Graph& g, SegmentPtr in_seg, SegmentPtr arm,
                                               VertexPtr vtx, const StmMichelArmThresholds& th);

    /// doc pdvd/83 (doc 78 action item 4): a Michel that leaves the muon
    /// chain a few cm BEFORE its stop.  On the owner's scan the four missed
    /// Michels of doc 78 sec 3.2 all hang off the chain's penultimate vertex
    /// (2.5-6.5 cm before the stop, the last chain segment a short stub), so
    /// the stop-arm classifier never sees them and the interior-vertex one
    /// calls them kOther.  This offers every such arm the STOP-arm gate.
    ///
    /// Walks the chain's interior vertices from the stop backward
    /// (vi = chain.size()-1 ... 1; chain_vtxs[vi] joins chain[vi-1] and
    /// chain[vi]) while the along-chain distance to the stop,
    /// sum_{j >= vi} segment_track_length(chain[j]), is <= max_dist.  At each
    /// vertex, every out-edge segment that is not a chain segment and not
    /// skip(arm) is re-read with stm_michel_classify_chain_arm; kDelta and
    /// kHadron are left alone (the caller already acted on them), the rest
    /// are counted in n_examined and classified with
    /// stm_michel_classify_stop_arm(g, chain[vi-1], arm, v, th) -- the same
    /// gate, with the incoming chain segment as the kink reference.  The
    /// NEAREST vertex with at least one kMichel wins; its kMichel arms are
    /// returned longest first (graph index breaks ties).  max_dist <= 0 or a
    /// chain shorter than 2 segments returns vtx_index -1.  Deterministic:
    /// out-edges are walked in sorted order.
    struct StmMichelNearArms {
        int vtx_index{-1};                 // index into chain_vtxs of the vertex the Michel leaves; -1 = none
        double dist{-1};                   // along-chain distance stop -> that vertex (internal units)
        int n_examined{0};                 // kOther body arms offered the gate, nearest vertex out to the winner (or max_dist)
        std::vector<StmMichelArm> michel;  // kMichel arms at that vertex, longest first
    };
    StmMichelNearArms stm_michel_near_stop_arms(Graph& g, const std::vector<SegmentPtr>& chain,
                                                const std::vector<VertexPtr>& chain_vtxs, double max_dist,
                                                const StmMichelArmThresholds& th,
                                                const std::function<bool(const SegmentPtr&)>& skip);

    /// doc pdvd/57: the STOP RETREAT.  `find_first_kink`'s charge gate wants
    /// BOTH arms of a kink >= 0.6 MIP, so an asymmetric muon->Michel junction
    /// (Bragg on one side, 0.1-0.4 MIP on the other) returns the sentinel and
    /// `CheckSTM_Michel` clamps the stop to the fit's far end -- the tip of
    /// whatever the Steiner path's boundary-to-boundary search reached, not
    /// where the muon actually stopped.  `stop_extend_max` only walks the
    /// stop OUTWARD (:1298); this is the mirror, graph-only so it needs no
    /// fitter to test.
    ///
    /// Walks the chain's OWN vertices backward from its current end, up to
    /// `max_drop` times, dropping one trailing chain segment per iteration
    /// while ALL of:
    ///   - the segment(s) already dropped plus this one span <= max_drop_len
    ///     (michel_max_len_cm: the same ceiling an attached Michel arm faces);
    ///   - the LIVE points on the dropped tail (>= min_tail_pts of them) have
    ///     median dQ/dx < collapse_frac * plateau -- a collapsed tail, not a
    ///     continuation (a muon at 0.7-1.3 MIP plateau cannot pass this at
    ///     collapse_frac 0.5 without a physically implausible plateau);
    ///   - the profile that SURVIVES the drop still peaks >= peak_frac *
    ///     plateau somewhere within peak_window of its new end -- there is a
    ///     Bragg rise to retreat TO, not just charge to drop.
    /// `plateau` is the SAME median dQ/dx over [plateau_lo, plateau_hi] (with
    /// the short-track halving) that stm_michel_bragg_contrast uses, computed
    /// ONCE from the profile as first passed in -- never recomputed in the
    /// shrinking frame, so the reference does not chase the retreat.
    ///
    /// `prof` must be the profile BEFORE stm_michel_profile_live (seg_idx and
    /// L/rr must match the chain the caller will actually shrink); the
    /// function filters dead points internally via min_dqdx_live so a dead
    /// stretch cannot fake a collapse.  Returns n_drop = 0 (a default-built
    /// result) when the input is too short to judge, when no candidate tail
    /// qualifies, or when max_drop <= 0.
    struct StmMichelRetreatThresholds {
        int    max_drop{0};                   // stop_retreat_max; <= 0 = off
        double collapse_frac{0.5};            // retreat_collapse_frac
        double peak_frac{1.4};                 // retreat_peak_frac
        double peak_window{15 * units::cm};   // retreat_peak_window_cm
        double max_drop_len{25 * units::cm};  // reuses michel_max_len_cm
        double plateau_lo{20 * units::cm};    // reuses bragg_plateau_lo_cm
        double plateau_hi{40 * units::cm};    // reuses bragg_plateau_hi_cm
        double min_dqdx_live{0};              // reuses profile_min_dqdx_frac * mip_dqdx
        int    min_tail_pts{2};
        // doc pdvd/74 (P3): how the dropped tail is read.  Both false = the
        // doc pdvd/57 reading above, verbatim.
        //  tail_strict:  only rows strictly past the vertex the retreat lands
        //    on.  That vertex's row (written twice, once per segment) is the
        //    kept segment's end: on an overshoot it carries the Bragg peak and,
        //    over a 2-3 cm segment of a few rows, sets the median by itself
        //    (039252_16/32: 1.17 x plateau with it, 0.47 without).
        //  tail_sublive: the live floor is not applied to the tail.  A
        //    collapsed overshoot reads 0.1-0.2 MIP, which is exactly what the
        //    floor calls a dead cell (039253_3/61).  The peak test stays
        //    live-only.  HAZARD: a dead-channel stretch past a live Bragg rise
        //    then reads as a collapse -- nothing here knows the channel map.
        bool   tail_strict{false};
        bool   tail_sublive{false};
        // doc pdvd/82 (doc 78 action item 2): an ALTERNATIVE admission for the
        // collapsed-tail test, judged against the surviving profile's PEAK
        // instead of the plateau.  On the population doc 78 sec 2.3 named, the
        // fit rides through the Michel: after a Bragg peak of 1.5-2.9 x plateau
        // the tail falls back to 0.54-1.87 x plateau -- far below what a
        // post-Bragg muon carries, but not below the track's own plateau, so
        // collapse_frac never fires.  tail_med <= tail_peak_frac * peak admits
        // that shape.  It is much looser than the plateau test (with peak /
        // plateau in 1.4-3, tail_peak_frac 0.5 is tail <= 0.7-1.5 x plateau),
        // so it is paired with tail_peak_kink_min: the row's OWN trajectory
        // bend (stm_michel_row_kink_deg, the doc 58 discriminator) must reach
        // it.  The bend is required ONLY when this test is what admits -- a row
        // the plateau test already accepts is unaffected, which is what keeps
        // the OFF path and the doc 57/58 path byte-identical.  0 = off.
        double tail_peak_frac{0};             // stop_tail_peak_frac; <= 0 = off
        double tail_peak_kink_min{25.0};      // stop_tail_peak_kink_min_deg
        double dir_window{5 * units::cm};     // the bend's arm length; = split_dir_window_cm
    };
    struct StmMichelRetreat {
        int n_drop{0};              // chain segments to pop from the back
        double drop_len{0};         // their total segment_track_length
        double plateau{0};          // the reference plateau this was judged against
        double last_tail_med{0};    // the last-accepted drop's tail median (diagnostic)
        double last_peak{0};        // the last-accepted drop's surviving peak (diagnostic)
        bool by_tail_peak{false};   // doc pdvd/82: the peak-relative test is what admitted the last drop
        double last_kink_deg{-1};   // doc pdvd/82: the bend at the boundary row (diagnostic; -1 = not measured)
    };
    StmMichelRetreat stm_michel_stop_retreat(const StmMichelProfile& prof, int n_chain_segs,
                                             const StmMichelRetreatThresholds& th);

    /// The bend of the FITTED trajectory at profile row `i`: the angle between
    /// the incoming direction (row i back to the row `window` of arclength
    /// behind it, i.e. larger rr) and the outgoing direction (row i forward to
    /// the row `window` ahead, i.e. smaller rr).  0 deg = straight; a real
    /// second particle at the stop bends this, a through-going track's own
    /// multiple scattering mostly does not (doc 58 sec 1: median 18.5 deg on
    /// the missed-collapse population that has no chain vertex to retreat
    /// onto, vs 6.3 deg on the collapse-shaped through-going control -- the
    /// owner's kink discriminator, applied at a fit row instead of a segment
    /// pair). Returns -1 when `i` is too close to either end of `prof` to
    /// measure a `window`-long arm on both sides.
    double stm_michel_row_kink_deg(const StmMichelProfile& prof, size_t i, double window);

    /// doc pdvd/58 (T1c): the STOP SPLIT.  stm_michel_stop_retreat can only
    /// move the stop onto a vertex the chain ALREADY has; on a population of
    /// missed collapse-shaped stoppers the collapse sits INSIDE the last
    /// chain segment (no chain vertex near it at all -- doc 57 sec 1's 23-item
    /// finding), so no graph-local retreat reaches them.  This function picks
    /// a FIT ROW inside that last segment for the caller to split at
    /// (PR::break_segment -- CheckSTM_Michel already uses it in
    /// anchor_vertex()), rather than an existing vertex.
    ///
    /// Losing the vertex constraint removes the guard that kept the retreat
    /// honest (doc 57 sec 1: a vertex within 4cm existed on 24/51 missed vs
    /// 11/80 through-going), so this function requires instead that the row's
    /// OWN trajectory bend (stm_michel_row_kink_deg) reach kink_min_deg -- a
    /// real kink in the fitted path, not just a charge-shape collapse.  Doc 58
    /// sec 1 measured this discriminator BEFORE this function existed: over
    /// the population with no chain vertex, bend >= 15 deg keeps every one of
    /// 58 through-going negative-control items from qualifying while still
    /// reaching 3 of 23 missed ones (probe numbers; the real graph typically
    /// recovers fewer, as it did for stm_michel_stop_retreat's 7 -> 2).
    ///
    /// Candidate rows: seg_idx[i] == n_chain_segs - 1 (INSIDE the last chain
    /// segment only -- anywhere else already has a vertex and belongs to
    /// stm_michel_stop_retreat, never this function), min_drop <= rr[i] <=
    /// max_drop_len, kink >= kink_min_deg, PLUS the same collapsed-tail /
    /// Bragg-rise-to-retreat-to tests stm_michel_stop_retreat applies (same
    /// plateau reference, same running-median-then-window convention). Among
    /// qualifying rows, the LARGEST kink wins (doc 58 sec 1: this beat
    /// anchoring on the shape classifier's own collapse onset by costing 0
    /// rather than 2 new false positives at the same recovery).
    ///
    /// `prof` must be the profile BEFORE stm_michel_profile_live, exactly as
    /// stm_michel_stop_retreat requires. Returns ok=false (a default-built
    /// result) when max_split <= 0, prof is empty, n_chain_segs <= 0, or no
    /// candidate row qualifies.
    struct StmMichelSplitThresholds {
        int    max_split{0};                  // stop_split_max; <= 0 = off
        double kink_min_deg{15.0};             // split_kink_min_deg
        double min_drop{3 * units::cm};        // split_min_drop_cm
        double collapse_frac{0.5};             // split_collapse_frac
        double peak_frac{1.4};                 // split_peak_frac
        double peak_window{15 * units::cm};   // split_peak_window_cm
        double dir_window{5 * units::cm};     // split_dir_window_cm
        double max_drop_len{25 * units::cm};  // reuses michel_max_len_cm
        double plateau_lo{20 * units::cm};    // reuses bragg_plateau_lo_cm
        double plateau_hi{40 * units::cm};    // reuses bragg_plateau_hi_cm
        double min_dqdx_live{0};              // reuses profile_min_dqdx_frac * mip_dqdx
        int    min_tail_pts{3};
        // doc pdvd/82: the same alternative collapsed-tail admission the
        // retreat carries -- see StmMichelRetreatThresholds.  Here the row's
        // bend is already required to reach kink_min_deg, so tail_peak_kink_min
        // is the HIGHER bar this looser tail reading must additionally clear
        // (the effective requirement is max(kink_min_deg, tail_peak_kink_min)).
        double tail_peak_frac{0};             // stop_tail_peak_frac; <= 0 = off
        double tail_peak_kink_min{25.0};      // stop_tail_peak_kink_min_deg
    };
    struct StmMichelSplit {
        bool   ok{false};
        size_t index{0};            // row into prof.pts/L/rr to split the segment at
        double cut_rr{0};           // prof.rr[index]: residual range of the chosen split
        double drop_len{0};         // prof.total_length - prof.L[index]
        double kink_deg{0};         // the row's own bend, for the record
        double plateau{0};          // the reference plateau this was judged against
        double tail_med{0};         // the dropped tail's median (diagnostic)
        double peak{0};             // the surviving profile's peak (diagnostic)
        bool by_tail_peak{false};   // doc pdvd/82: the peak-relative test is what admitted this row
    };
    StmMichelSplit stm_michel_stop_split(const StmMichelProfile& prof, int n_chain_segs,
                                         const StmMichelSplitThresholds& th);

    /// Reject bits carried in the stm_michel PC / T_stm_michel.reject_bits.
    /// is_stm = (reject_bits == 0).  Michel presence is deliberately NOT a
    /// criterion (mu- capture in argon).
    enum StmMichelReject : unsigned {
        R_NO_CHAIN           = 1u << 0,   // no stm_pass/stm_fit PC, or no route entry -> stop
        R_STOP_UNMATCHED     = 1u << 1,   // no graph vertex near the tagger's stop; chain walked greedily
        R_NO_BRAGG           = 1u << 2,   // contrast < bragg_contrast_min * expected
        R_SHAPE_FLAT         = 1u << 3,   // ks_mu + ks_margin >= ks_flat
        R_NOT_MUON_PID       = 1u << 4,   // do_track_comp: proton or electron template beats muon, or gate 0
        R_CONTINUATION       = 1u << 5,   // a collinear MIP track leaves the stop
        R_STOP_NEAR_BOUNDARY = 1u << 6,   // stop outside the fiducial inset by stop_fv_margin
        R_VERTEX_HADRON      = 1u << 7,   // a long heavily-ionizing prong off the body
        R_SHORT              = 1u << 8,   // fewer than min_chain_points profile points
        R_PROFILE_SPARSE     = 1u << 9,   // too few LIVE points to judge: a Bragg window or the compare range holds < 3 (doc pdhd/03)
        R_PLATEAU_OFF_MIP    = 1u << 10,  // plateau_med / mip_dqdx outside [plateau_mip_lo, plateau_mip_hi] (doc pdhd/03)
        R_STOP_INTO_DEAD     = 1u << 11,  // the visible end walks into a dead region (FiducialUtils::check_dead_volume) (doc pdhd/03)
        R_CLUSTER_NOT_TRACK  = 1u << 12,  // too few of the cluster's points lie on the reconstructed track (doc pdhd/03)
        R_PROFILE_GEOMETRY   = 1u << 13,  // the profile is not a measurement: a coiled end (arc/span) or a long fitted segment the charge does not support (doc pdvd/66)
    };

    /// doc pdvd/70 (P1, topology_stop_evidence): the reject bits that a Michel
    /// of sufficient quality at the stop clears -- the owner's rule that a
    /// Michel at the end is by itself strong evidence of a stop, while the
    /// dQ/dx rise counts only when it genuinely matches a Bragg peak.
    /// Returns the SUBSET of `reject_bits` to clear: R_NO_BRAGG | R_SHAPE_FLAT
    /// (| R_PROFILE_SPARSE when `clears_sparse`) when `michel_found`, the
    /// object is attached (conn_type 1) or bridged (2), ke_mev >= ke_min_mev
    /// and len_cm >= len_min_cm; 0 otherwise.  Every other bit is left to the
    /// caller, so is_stm = (bits & ~clear) == 0 can only move 0 -> 1.  A
    /// non-finite energy or length never qualifies.
    unsigned stm_michel_topology_clear(unsigned reject_bits, int michel_found, int conn_type,
                                       double ke_mev, double len_cm,
                                       double ke_min_mev, double len_min_cm, bool clears_sparse);

    /// doc pdvd/71 (P4, michel_gamma_collect): the per-cluster test for an
    /// isolated gamma blob that belongs to the Michel -- the owner's three
    /// criteria: it lies along the Michel electron's direction, it is a
    /// compact dot near the stop, and it carries a gamma's energy, not an
    /// over-clustered lump's.  Returns 0 to accept, else the FIRST gate that
    /// fails, in this order: 6 a non-finite input; 1 d_stop > radius; 2 len >
    /// max_len; 3 cos_dir < cos_min (cos of the angle between the stop -> blob
    /// centroid and the stop -> Michel-object centroid directions); 4 d_body <=
    /// d_mich (the blob is at least as close to the muon body as to the Michel
    /// object and the stop); 5 ke_mev > max_ke_mev.  All boundaries inclusive
    /// on the accepting side.  Lengths in one unit, energies in MeV.
    int stm_michel_gamma_gate(double d_stop, double len, double cos_dir, double d_mich,
                              double d_body, double ke_mev, double radius, double max_len,
                              double cos_min, double max_ke_mev);

    /// doc pdvd/71 (P4): the total-energy guard -- over-clustering can hand the
    /// Michel a huge energy, so the object never grows past `total_max_mev`.
    /// Walks `ke_mev` in the caller's order (nearest first) from a running total
    /// of `core_ke_mev` and takes a blob only when the total stays <= the cap;
    /// a blob that does not fit is skipped and the walk continues.  Returns a
    /// 0/1 mask parallel to `ke_mev`.  A non-finite core takes nothing; a
    /// non-finite or negative blob energy is never taken.
    std::vector<int> stm_michel_gamma_take(double core_ke_mev, const std::vector<double>& ke_mev,
                                           double total_max_mev);

    /// doc pdvd/84 (doc 78 action item 3): which exemption, if any, spares a
    /// would-be moved-stop veto (T2c, moved_stop_michel_guard) of an attached
    /// Michel.  kKink when kink_min_deg >= 0 and kink_deg >= kink_min_deg (doc
    /// pdvd/72, checked first so its precedence and counter are unchanged);
    /// else kReach when reach_min_cm >= 0 and reach_cm >= reach_min_cm, where
    /// reach is the arm's length plus its far subtree (michel_len +
    /// michel_far_len); else kVeto.  A threshold < 0 is off; kink_deg < 0
    /// (unmeasurable) never spares by the kink; a non-finite value never
    /// spares.  Degrees and cm, plain.
    enum class StmMichelMovedStopSpare { kVeto, kKink, kReach };
    StmMichelMovedStopSpare stm_michel_moved_stop_spare(double kink_deg, double reach_cm,
                                                        double kink_min_deg, double reach_min_cm);

    /// doc pdvd/85 (doc 78 action item 5): whether a candidate's capture
    /// gammas (the doc pdvd/51 stage, role 5) are withheld -- true exactly when
    /// the knob is on and the verdict rejects the candidate (any reject bit
    /// set, i.e. is_stm 0).  A capture gamma presupposes a muon that stopped.
    bool stm_michel_stop_gamma_withhold(bool require_stm, unsigned reject_bits);

    /// doc pdvd/85: the keep-mask for an order-preserving row erase -- 0 where
    /// roles[i] == drop_role, 1 elsewhere, parallel to `roles`.
    std::vector<char> stm_michel_rows_keep(const std::vector<int>& roles, int drop_role);

}  // namespace WireCell::Clus::PR

#endif
