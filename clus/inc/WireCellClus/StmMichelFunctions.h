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

    /// Classify a non-chain arm at an INTERIOR chain vertex.
    ///  - Delta: len <= delta_max_len, far_len <= delta_max_len, terminal.
    ///  - Hadron: longer than delta_max_len and mip > hadron_mip.
    ///  - else Other.
    StmMichelArm stm_michel_classify_chain_arm(Graph& g, SegmentPtr in_seg, SegmentPtr arm,
                                               VertexPtr vtx, const StmMichelArmThresholds& th);

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
    };

}  // namespace WireCell::Clus::PR

#endif
