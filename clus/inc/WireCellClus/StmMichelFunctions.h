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

    /// Median of a copy of v; 0 when empty.
    double stm_michel_median(std::vector<double> v);

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
    };

}  // namespace WireCell::Clus::PR

#endif
