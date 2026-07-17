#ifndef WIRECELL_MATCH_QLMATCHING
#define WIRECELL_MATCH_QLMATCHING

#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IFiducial.h"
#include "WireCellIface/ITensorSetFanin.h"

#include "WireCellMatch/PhotonLibraryModel.h"
#include "WireCellMatch/SemiAnalyticalModel.h"
#include "WireCellMatch/TimingTPCBundle.h"

#include "WireCellUtil/Units.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace WireCell::Match {

    /// Charge-light matcher.
    ///
    /// Port of larwirecell/qlmatch/QLMatching to wire-cell-toolkit, without
    /// any larsoft/fhicl/art dependencies. The SBND semi-analytical optical
    /// model is loaded from a JSON file at configure time (see
    /// cfg/sbnd/semi-analytical-sbnd.json or equivalent).
    ///
    /// Input (ITensorSetFilter): pctree tensor set holding clusters, with the
    /// per-event optical flashes attached as the canonical "flash"/"light"/
    /// "flashlight" point clouds on the live root node (see
    /// Aux::FlashTensorToOpticalPCs, the same schema as root/UbooneClusterSource).
    /// Output: cluster tensor set with per-cluster t0 AND a per-cluster
    /// "flash" scalar (matched flash row index) set from the matched flash, so
    /// Clus::Facade::Cluster::get_flash() reflects the match.
    class QLMatching : public Aux::Logger,
                       public ITensorSetFanin,
                       public IConfigurable {
    public:
        QLMatching();
        virtual ~QLMatching();

        // Fanin: one input port per anode (multiplicity = #anodes). The
        // single-anode default (multiplicity 1) wires and behaves exactly as the
        // historical single-input filter; a multi-anode list matches both APAs in
        // one node and merges the results (see operator()).
        std::vector<std::string> input_types() override;
        bool operator()(const input_vector& invec, output_pointer& out) override;

        void configure(const WireCell::Configuration& cfg) override;
        WireCell::Configuration default_configuration() const override;

    private:
        std::size_t m_count{0};

        // Number of optical-detector channels (per-flash PE vector length).
        // The only bit-safe source of nchan now that flashes arrive via the
        // canonical PCs (do NOT derive from max channel ident).
        int m_nchan{312};

        // Anode face whose sensitive drift volume defines the per-TPC geometry
        // (DetectorVolumes::inner_bounds wpid). Default 0 (SBND APAs image through
        // face 0). PDHD's +x drift side (APA1/APA3 group) images through face 1, so
        // its matching node sets tpc_face=1; the -x side (APA0/APA2) keeps face 0.
        int m_tpc_face{0};
        // Optional per-input face list ("tpc_faces"), one entry per input port, for a
        // JOINT multi-side node whose inputs image through different faces (PDHD: -x
        // side face 0, +x side face 1). Empty (default) => use the single m_tpc_face
        // for every input => bit-identical for SBND / the per-side PDHD nodes.
        std::vector<int> m_tpc_faces;
        // Optional per-input SECOND face ("tpc_extra_faces") to union into the active-
        // volume box alongside tpc_face(s). PDHD/SBND anode faces span the full Y/Z and
        // differ only slightly in X, so a single face's sensitive box already covers the
        // TPC; PDVD's two faces instead split the Y range in half (adjacent, disjoint),
        // so using only one face truncates the active-volume box (and the light-
        // prediction inclusion gate) to roughly half the true Y extent. Empty (default)
        // => no second face unioned => bit-identical for PDHD/SBND.
        std::vector<int> m_tpc_extra_faces;

        // PMTs on vs off (default true => apply the OpDet type mask); see
        // ch_mask for further per-channel disables.
        bool m_pmts{true};
        // OpDet types kept by the matching mask, matched against
        // SemiAnalyticalModel::OpticalDetector::type (1 = dome PMT,
        // 0 = (X)Arapuca). Default {1} = PMTs only, reproducing the historical
        // SBND selection. The mask is derived per-channel from the injected
        // OpDet table rather than a hard-coded SBND layout, so it generalizes to
        // other detectors via config.
        std::vector<int> m_active_opdet_types{1};
        bool m_data{true};
        std::vector<int> m_ch_mask;

        // Per-event dynamic PMT auto-mask (off by default => production bit-identical).
        // Within a single event/TPC, drop a PMT that never fires (max PE over the
        // event's flashes < pe_low) while its nearest live PMTs do (>= min_contrast
        // flashes whose neighbour-median PE > pe_bright). Catches a channel that is
        // dead in THIS run but not in the static ch_mask. Folds into run.opdet_mask, so
        // prediction / chi2 / KS / ndf all inherit it. See match/docs/qlmatching-code.md.
        bool m_auto_mask{false};
        double m_auto_mask_pe_low{5.0};      // a PMT "fires" if its max event PE >= this
        int    m_auto_mask_neighbors{4};     // K nearest live PMTs for the brightness ref
        double m_auto_mask_pe_bright{50.0};  // neighbour-median PE meaning "light present"
        int    m_auto_mask_min_contrast{1};  // # bright-neighbour flashes required
        int    m_auto_mask_min_flash{3};     // skip auto-masking below this flash count
        // Restrict the auto-mask neighbour pool to OpDets of the SAME type as the
        // candidate channel. For a mixed-PD detector (PDVD: X-ARAPUCAs ~0.03 eff vs
        // PMTs ~0.12/0.036) a (Y,Z)-nearest neighbour can be a different PD type
        // whose brightness is not a meaningful reference. Default false => the
        // historical all-active-type pool (bit-identical).
        bool m_auto_mask_same_type{false};

        // ---- PDVD / vertical-drift knobs (all default OFF => other detectors
        // byte-identical). See match/docs/qlmatching-code.md. ----
        //
        // shared_flash: the detector produces ONE all-PD flash list shared by every
        // input port (PDVD: a single flash covers both drift volumes; the cathode
        // XAs at x=0 see both sides). The per-run pipeline through cull_inconsistent
        // is unchanged, but the two LASSO rounds run JOINTLY: each physical flash
        // (keyed by its input tensor row id, identical across ports) gets ONE
        // measured block and columns from EVERY port's candidate bundles, so
        // clusters from both drift sides share the explanation of one flash instead
        // of each side independently absorbing the full PE (which would bias every
        // cathode-crosser). Requires >1 input port; the merged-root opflash PC and
        // matched_flash_gid are emitted once (side-0 keyed). The empty-flash and
        // cluster rescues are NOT shared-flash-aware and are skipped with a warning.
        bool m_shared_flash{false};
        // opdet_all_volumes: keep every active OpDet in every run's mask instead of
        // splitting by cathode-plane x (PDVD: the cathode XAs sit AT x=0 and are
        // double-sided; a single flash needs all PDs on both drift-side runs).
        bool m_opdet_all_volumes{false};
        // vd_surface_flags: PD-surface-aware endpoint flags. The historical
        // flag_close_to_PMT test assumes the PDs sit behind the anode plane (SBND /
        // PDHD horizontal drift). For PDVD the PDs are at the cathode (XAs), on the
        // two y walls (membrane XAs) and behind the BOTTOM anode only (PMTs):
        //  - anode-end proximity sets flag_close_to_PMT only for an input whose
        //    anode_pd_channels list is non-empty (PDVD: bottom volume yes, top no);
        //    flag_at_x_boundary is set regardless (unchanged T0 semantics).
        //  - proximity to a y wall with a non-empty pd_wall_channels list sets
        //    flag_close_to_PMT (drift-invariant, memoized per cluster).
        //  - the cathode end additionally sets the inert flag_at_cathode.
        // Each trigger also records its surface's PD channels on the bundle
        // (TimingTPCBundle::add_relax_channels), so the chi2_relax denominator
        // widening applies only to the RELEVANT PDs, not any channel.
        bool m_vd_surface_flags{false};
        double m_pd_wall_cushion{10 * units::cm};  // wall-proximity band
        std::vector<int> m_pd_wall_channels_ylo;   // PD channels on the low-y wall (empty => wall inactive)
        std::vector<int> m_pd_wall_channels_yhi;   // PD channels on the high-y wall (empty => wall inactive)
        // Per-input-port PD channels behind that input's anode (empty list => that
        // anode hosts no PDs => no anode-end flag_close_to_PMT in vd mode).
        std::vector<std::vector<int>> m_anode_pd_channels;
        bool m_beamonly{false};
        double m_flash_minPE{500};
        double m_flash_mintime{-1.5 * units::ms};
        double m_flash_maxtime{1.5 * units::ms};
        // Channel-scoped flash admission (on top of flash_minPE): admit a flash
        // only if, over flash_sel_channels, sum(PE) >= flash_sel_minPE AND at
        // least flash_sel_min_fired channels read >= flash_sel_fired_pe.  Meant
        // for detectors where one PD family is the only reliable ruler (PDVD
        // cathode X-ARAPUCAs: full-stream readout, 100% coverage, proportional
        // response), so a flash whose evidence lives only on unreliable /
        // self-trigger-silent families never enters matching.  Empty channel
        // list (default) => no-op => bit-identical legacy admission.
        std::vector<int> m_flash_sel_channels;
        double m_flash_sel_minPE{0.0};
        int m_flash_sel_min_fired{0};
        double m_flash_sel_fired_pe{1.0};
        double m_beam_mintime{-5 * units::us};
        double m_beam_maxtime{5 * units::us};
        double m_QtoL{0.5};
        // Compute reflected/VIS light from the semi-analytical model. SBND has
        // cathode WLS-reflector foils (true => bit-identical). Detectors with no
        // reflected component (e.g. PDHD, DUNE DoReflectedLight=false) set false:
        // the model then skips loading the VISHits tables (so an empty VISHits is
        // allowed) and detectedReflectedVisibilities returns all-zero.
        bool m_doReflectedLight{true};
        // Drift speed used for the per-flash X correction (flash_x_offset =
        // sign * flash_time * drift_speed). In WCT internal units (length/time,
        // i.e. mm/ns numerically). Default is the historical hard-coded value;
        // configs should pass the common params.lar.drift_speed via jsonnet.
        double m_drift_speed{1.563 * units::mm / units::us};
        // Per-event trigger offset folded into the per-flash X correction
        // (flash_x_offset = sign * (flash_time + trigger_offset) * drift_speed).
        // Detectors that DON'T bake the readout-vs-trigger offset into the charge
        // x at imaging time (e.g. PDHD, whose BlobSampler time_offset stays 0) pass
        // the per-event offset_us here so the matching geometry lands on the same
        // trigger time base as the flash times. In WCT internal time units (ns).
        // Default 0 => detectors that DO bake it (e.g. SBND) are bit-identical.
        double m_trigger_offset{0.0};
        // Per-INPUT trigger offsets ("trigger_offsets", one per input port), for
        // detectors whose charge readout windows start at different times per
        // drift volume (PDVD: the TDE/BDE crates open up to ~32 us apart, each
        // wandering event-to-event vs the trigger-locked light window). Empty
        // (the default) => m_trigger_offset for every input, bit-identical.
        std::vector<double> m_trigger_offsets;
        double trigger_offset_for(std::size_t input_idx) const
        {
            return m_trigger_offsets.empty() ? m_trigger_offset
                                             : m_trigger_offsets.at(input_idx);
        }
        // Per-INPUT drift speeds ("drift_speeds", one per input port), for
        // detectors whose drift volumes can carry different calibrated speeds
        // (PDVD: bottom volume = port 0, top volume = port 1). Empty (the
        // default) => m_drift_speed for every input, bit-identical.
        std::vector<double> m_drift_speeds;
        double drift_speed_for(std::size_t input_idx) const
        {
            return m_drift_speeds.empty() ? m_drift_speed
                                          : m_drift_speeds.at(input_idx);
        }
        // LASSO solution threshold below which a (flash, cluster) bundle is
        // dropped after each matching round. Was hard-coded 0.05 inline; pulled
        // out so it can be widened/narrowed from the jsonnet without rebuild.
        double m_strength_cutoff{0.05};

        // Assemble the round-1/2 LASSO response matrices (P, PF) as Eigen sparse
        // and form the normal-equations X = PᵀP + PFᵀPF / y = PᵀM + PFᵀMF with
        // sparse products, instead of the dense Eigen path. P/PF are block-sparse
        // (each bundle column touches only its flash's nopdet-row block), so on a
        // busy event (e.g. PDHD run 29107 evt 1015, ~440 flashes) the dense path's
        // P/PT spike (~5 GB) and the O(nbeta²·nopdet·nflash) matrix build dominate
        // both memory and time. The sparse path removes both. Sparse and dense
        // matrix products accumulate FP sums in different orders, so X (hence the
        // diagnostic LASSO strength) can differ at the ULP level while every match
        // assignment is unchanged -> default OFF (dense, byte-identical to history);
        // PDHD qlmatching.jsonnet turns it on. See match/docs/qlmatching-perf-evt1015-pdhd.md.
        bool m_sparse_lasso{false};

        // ---- Tuning constants pulled out of the inline code (see
        // match/docs/improve_progress.md). The §B-§G defaults equal the historical
        // hard-coded literals (bit-identical if omitted); the §A cushions below
        // default to the MicroBooNE convention and intentionally differ. ----

        // §A active-volume cushions. The raw active-volume bounds (anode/cathode
        // X, Y span, Z span) now come from the IDetectorVolumes service
        // (m_dv->inner_bounds) so the code is detector-agnostic. These signed
        // cushions adjust the effective PE-inclusion and boundary-flag windows
        // relative to the true geometry, in the per-TPC anode->cathode drift
        // coordinate u (u=0 at the anode/PMT plane, u=u_cathode at the cathode
        // seam), so a single set serves both reversed-drift SBND APAs.
        //
        // Defaults follow the MicroBooNE prototype convention
        // (ToyMatching.cxx::calculate_pred_pe ~158-164): outward-positive in u.
        // They shift the inclusion window slightly vs the old SBND literals
        // (intended); overriding them from jsonnet can reproduce the old bounds
        // bit-identically. See match/docs/improve_progress.md §A.
        double m_anode_ext1{-2.0 * units::cm};    // PE-inclusion window edge, below the anode (low_x_cut_ext1)
        double m_anode_ext2{4.0 * units::cm};     // anode flag-window outer edge (low_x_cut_ext2)
        double m_cathode_ext1{1.2 * units::cm};   // PE-inclusion window edge, beyond the cathode (high_x_cut_ext1)
        double m_cathode_ext2{-2.0 * units::cm};  // cathode flag-window inner edge (high_x_cut_ext2)
        // Extra tolerance BELOW anode_ext1, subtracted to form the anode-side floor
        // used in TWO coupled places in compute_endpoint_flags: the containment gate
        // (first_u > anode_ext1 - margin) and the anode flag-window inner edge
        // (first_u > anode_ext1 - margin, which sets flag_at_x_boundary /
        // flag_close_to_PMT). The two are deliberately the SAME bound: a cluster
        // admitted by the containment slack must also acquire the at-x-boundary /
        // close-to-PMT flags, else it would be a candidate that has silently dropped
        // its T0-crosser constraint. Widening moves the floor outward (deeper past
        // the anode) only -- the flag window's outer edge stays m_anode_ext2 -- so it
        // cannot re-flag anode-SHORT ends the way a cathode_ext2 widening once did.
        // Prototype hard-codes 1.0 cm (ProtectOverClustering.cxx:296); the knob keeps
        // that as the default => byte-identical. PDVD production sets 2.0 cm (floor
        // -2 -> -4 cm; run 039252 evt298567 full-gap crosser top:22 missed by 0.5 cm).
        double m_anode_ext1_margin{1.0 * units::cm};
        double m_y_cushion{0.0 * units::cm};       // signed inward(+)/outward(-) shift of each |y| edge
        double m_z_cushion{0.0 * units::cm};       // signed inward(+)/outward(-) shift of each z edge
        // Proximity band for the diagnostic flag_two_boundary edge test: a main-PCA
        // extreme counts as "at the detector edge" if it is within this distance of
        // any of the 6 per-APA active-volume faces (anode, cathode, ±y, ±z).
        double m_two_boundary_margin{3.0 * units::cm};
        // Robust containment endpoint (default OFF => byte-identical). The cathode/anode
        // inward-walk trims in compute_endpoint_flags only break (and trim) at an
        // *isolated* deep straggle separated from the body by a >0.75 cm gap. A thin
        // off-axis overclustering tail that stays within 0.75 cm of the body is never
        // trimmed, so a 0.3%-of-points straggle drags the endpoint past the detector
        // edge and falsely fails containment. When true, a gap-independent post-pass
        // snaps each endpoint back inside the in-edge iff the outer material (points
        // beyond the edge) is sparse -- below max(frac*cluster_points, count). Point
        // mass (not blob count) is the sparsity measure, so dense genuine track-ends
        // (and at_x_boundary crossers) are untouched.
        bool m_robust_endpoint_trim{false};
        double m_robust_endpoint_frac{0.05};   // outer-straggle allowance, fraction of cluster points
        double m_robust_endpoint_count{0.0};   // outer-straggle allowance, absolute point floor
        double m_robust_endpoint_charge_frac{0.0};  // outer-straggle allowance, fraction of cluster
                                                    // charge (0 => disabled, point-count judge only).
                                                    // Trims low-charge overclustered satellites the
                                                    // point-count judge misses (dense real track-ends
                                                    // exceed the budget and are left to fail).
        double m_robust_endpoint_charge_abs{0.0};   // absolute per-point charge-DENSITY ceiling for
                                                    // trimmable outer material (0 => disabled). Size-
                                                    // independent: keeps a charge-dense real track tip
                                                    // (never trimmed) distinct from a charge-sparse
                                                    // overclustered satellite (trimmable), regardless
                                                    // of cluster size or tail length.
        // Gap-aware (detachment) anode-end trim (default OFF => byte-identical). The
        // charge_abs density gate above refuses MODERATELY-dense (~2000 q/pt)
        // overclustered junk because it overlaps the genuine track-end band
        // (~2500 q/pt). But such junk gives itself away another way: it sits in
        // gap-separated outer groups DETACHED from the contiguous track body, whereas
        // a real anode-piercing track end is continuous (adjacent drift slices < the
        // ~0.08 cm tick spacing apart). When m_robust_endpoint_gap > 0, the anode-end
        // walk trims the sub-anode material iff (a) some gap in it exceeds
        // m_robust_endpoint_gap (the detachment signal, NOT density) and (b) its charge
        // stays within m_robust_endpoint_gap_charge_frac of the cluster (the safety
        // bound: a real sub-anode stretch carrying more charge is left to fail). The
        // existing density path is untouched.
        double m_robust_endpoint_gap{0.0};            // detachment-gap threshold (0 => disabled)
        double m_robust_endpoint_gap_charge_frac{0.0};// outer-charge cap for the gap path
        // Mirror the gap path at the CATHODE end (default OFF => byte-identical, and in
        // particular PDHD -- which sets the two keys above -- is unaffected until it opts
        // in). The two knobs above were added for anode-end cases only, but nothing in
        // their justification is direction-specific: a real cathode-REACHING track end is
        // just as continuous as a real anode-piercing one, and the junk is just as
        // detached. The exposure is symmetric by construction, because the upstream cause
        // -- clustering_isolated's small->big merge (clustering_isolated.cxx: hardcoded
        // small_big_dis_cut = 80 cm, no angle/direction/gap test) -- has no directional
        // preference: it absorbs any "small" cluster into the nearest big one within
        // 80 cm, at whichever end it happens to lie. The always-on legacy trims above
        // already treat both ends as mirror images; only this knob-gated robust layer
        // broke that symmetry, as an accident of where the first cases turned up.
        // Compare PDHD evt 1007 uid-126 ("7 gap-separated groups marching from u=-111cm
        // up to the body", ~2000 q/pt, anode end) with PDVD 039252 evt298567 apa-4
        // ident 97 (6 gap-separated isolated-merged clumps, 5847 q/pt, CATHODE end):
        // same pathology, same upstream cause, opposite ends -- only the first had a
        // judge. When true the cathode-end walk reuses m_robust_endpoint_gap and
        // m_robust_endpoint_gap_charge_frac (same thresholds, same safety bound: a real
        // over-cathode stretch carrying more charge is left to fail, which is what stops
        // a wrong-T0 hypothesis being rescued). Density path untouched.
        // See pdvd/docs/qlmatch/16_pdvd-clus97-crosser-evt298567.md.
        bool m_robust_endpoint_gap_cathode{false};
        // Where the anode-end walk stops calling material "outside" (default OFF =>
        // byte-identical). The walk breaks at the first slice above anode_in
        // (= m_anode_ext1), but the containment gate it feeds is
        // first_u > anode_in - m_anode_ext1_margin -- a floor m_anode_ext1_margin cm
        // LOWER. So a cluster whose body legitimately starts in that dead band is
        // contained, yet the walk marches past the detached junk INTO the body before
        // it can break, hands the judges a slab of dense real track charge, and every
        // judge then (correctly) refuses to trim. The junk survives and containment
        // fails on it. Whether the rescue works at all then depends on cluster SIZE --
        // whether max(frac*npoints, count) happens to cover the swallowed body points
        // -- which is not physics: PDVD 039252 evt298567 apa-0 ident 34 (2649 pts,
        // allowance 26.5) is rescued while apa-4 ident 1 (1375 pts, allowance 15)
        // is not, on the same 20 swallowed points.
        // When true the walk instead breaks at the floor the gate actually uses, so it
        // never swallows material that is already inside the gate; the judges are
        // unchanged and still protect genuine track ends (a body starting BELOW the
        // floor still fails, as it must). Cathode end needs no counterpart: its gate
        // (last_u < cathode_in) and its walk share one threshold.
        // See pdvd/docs/qlmatch/15_pdvd-clus34-unmatched-evt298567.md.
        bool m_robust_endpoint_walk_to_floor{false};

        // §D pre-selection / bad-match gates.
        double m_mc_saturation_pe{5000};      // MC saturated-PMT mask trigger (total flash PE)

        // §E out-of-beam QA cuts.
        double m_outbeam_ks_max{0.2};
        double m_outbeam_chi2ndf_max{20};
        double m_outbeam_pe_frac{0.5};        // out-of-beam total-PE mismatch fraction

        // §C LASSO weights.
        double m_lasso_lambda{0.1};
        double m_delta_charge{0.01};
        double m_delta_light{0.025};
        double m_delta_shape{0.01};
        double m_bkg_weight{0.5};              // background-column weight (round 1)
        double m_pe_mismatch_knee{0.3};        // PE-mismatch weight knee (fraction of meas)
        double m_pe_mismatch_floor{0.3};       // PE-mismatch weight floor
        // Flag-aware per-column L1 down-weight (prototype ToyMatching.cxx:1437; default
        // OFF = bit-identical). When on, a boundary/near-PMT/window-truncated bundle's
        // weight is multiplied by m_lasso_boundary_weight so it is shrunk less and more
        // likely to survive the strength cutoff.
        bool   m_lasso_flag_weight{false};
        double m_lasso_boundary_weight{0.2};

        // §G flash PE-error model (forwarded to Opflash for the LASSO; the same
        // floor/frac/knee feed the bundle chi2 via BundleQualityParams).
        double m_pe_err_floor{0.3};
        double m_pe_err_frac{0.3};
        double m_pe_err_knee{1.0};
        // When true the bundle chi2 computes its per-PMT error from the PREDICTED pe
        // (not measured); the LASSO weight stays measured-based. SBND-on, default off.
        bool   m_pe_err_on_pred{false};
        // Efficiency-aware low-PE error inflation (PD detection inefficiency at low
        // light). When m_pe_err_lowpe_frac >= 0 (and m_pe_err_on_pred) the bundle-chi2
        // relative error grows as the predicted pe falls: rel = frac + (lowpe_frac-frac)
        // *exp(-pred/lowpe_knee), PE_err = sqrt((rel*pred)^2 + floor^2). <0 => disabled
        // (bit-identical). Calibrated on run-29107 hand scans (pdhd/ql_light_calib/
        // fit_lowpe.py): PDHD 2.1/4.0. Feeds the chi2 only; LASSO weight stays measured.
        double m_pe_err_lowpe_frac{-1.0};
        double m_pe_err_lowpe_knee{4.0};
        // Per-PD-family PE-error override (PDVD cathode-XA vs PMT calibration;
        // pdvd/docs/qlmatch/19_pdvd-scan-tuning-039252.md).  Parallel arrays:
        // family i covers the channels of pe_err_family_channels[i]; the
        // matching entry of pe_err_family_{floor,frac,lowpe_frac,lowpe_knee}
        // overrides the global value on those channels in BOTH error paths
        // (Opflash LASSO floor/frac; bundle-chi2 per_opdet_perr, all four).
        // A missing/negative entry falls back to the global.  Empty channel
        // list (default) => byte-identical legacy model.
        std::vector<std::vector<int>> m_pe_err_family_channels;
        std::vector<double> m_pe_err_family_floor;
        std::vector<double> m_pe_err_family_frac;
        std::vector<double> m_pe_err_family_lowpe_frac;
        std::vector<double> m_pe_err_family_lowpe_knee;
        // Resolved per-channel overrides (length m_nchan when any family is
        // configured, else empty; -1 entry => use the global value).
        std::vector<double> m_pe_err_ch_floor;
        std::vector<double> m_pe_err_ch_frac;
        std::vector<double> m_pe_err_ch_lowpe_frac;
        std::vector<double> m_pe_err_ch_lowpe_knee;
        // Optional per-channel measured-PE gain correction (length nchan), applied
        // to every flash's measured PE as it is read. Empty => identity (byte-
        // identical). PDHD uses it to scale up the gain-biased -x full-stream half
        // (optical ch 120-159) so its measured light matches the prediction.
        std::vector<double> m_measured_pe_scale;
        double m_flash_pe_threshold{0.0};      // Opflash "fired" threshold (PE)
        // Mask DAPHNE-rail-saturated channels PER FLASH (Opflash::get_sat,
        // fed by the OpHitFinder flag_saturation -> OpFlashFinder flash_sat
        // chain): the channel is dropped from that flash's opdet mask
        // (=> pred, chi2, KS via bundle_mask_ks) and its LASSO rows are
        // zeroed.  A railed channel's clipped PE can be a x2-10 underestimate
        // (pdvd/docs/qlmatch/pdvd-saturation-recovery.md) -- masking beats
        // both the veto (exact 0) and trusting the clipped value.  Default
        // OFF -> bit-identical.
        bool m_use_saturation_flag{false};

        // Whether a rail-flagged channel is also DROPPED from that flash's
        // opdet mask (=> out of chi2 and, under bundle_mask_ks, out of KS),
        // or stays in them at its measured PE.
        //
        // The clipped PE is NOT a missing measurement: it is a LOWER BOUND on
        // the true PE (11_pdvd-saturation-recovery.md §2.1 -- `clip` is a
        // "tight deterministic underestimate", -4% at depth d=1.4, -14% at
        // d=2, -40% at d=4, +-2% band), and with OpDecon saturation_repair on
        // (the PDVD production chain) it is a repaired estimate biased high by
        // ~+9% at d=2 / ~+23% at d=4 with a [1.06,2.05] band.  Either way the
        // channel carries real information, so dropping it discards a
        // constraint AND makes a wrong prediction free there.
        //
        // LASSO is deliberately NOT re-opened by this knob: a lower-bound
        // measurement in a least-squares solve biases the fitted charge DOWN.
        // The rail-flagged LASSO rows stay zeroed at every fit site, so the
        // railed channel acts as a CHECK on the charge-derived prediction
        // without constraining it.
        //
        // Pair with chi2_sat_inflate for the extra per-channel error.
        // Default TRUE = the legacy drop => bit-identical for existing cfgs.
        bool m_saturation_mask_fit{true};

        // Track readout coverage PER FLASH (Opflash::get_cov, fed by the
        // OpHitFinder emit_coverage -> OpFlashFinder flash_cov chain): a
        // self-triggered channel (PDVD membrane XA / PMT 16.4-us snippets,
        // duty ~5-30%) that carried no waveform over the flash window reads
        // measured = 0.  Turning this on makes the per-flash coverage
        // available to coverage_mask_fit and to the calib dump / ql_scan
        // "nodata" label.  Default OFF -> bit-identical; legacy archives have
        // get_cov == 1.
        bool m_use_coverage_flag{false};
        double m_coverage_min{1.0};

        // Whether an uncovered channel is also DROPPED from that flash's
        // opdet mask and LASSO rows (like a saturated one), or stays in the
        // fit at its measured PE of 0.
        //
        // Dropping was the original 2026-07-16 reading: no snippet => no
        // measurement.  The measurement says otherwise (18 evts run 039252,
        // pdvd/docs/qlmatch/scripts/analyze_coverage_deadtime.py, coverage
        // rebuilt against the chain's own flash_cov on 266071/266071 cells):
        //   - DAPHNE has no dead time -- the minimum gap between consecutive
        //     snippets on a channel is 0.336 us, so a channel is never busy
        //     when a flash arrives;
        //   - uncovered channels are QUIETER, not busier: median gap since
        //     their last snippet 93.8 us vs 54.2 us for covered ones;
        //   - the self-trigger threshold is ~1 PE (median snippet holds 3.3
        //     PE, 25% under 0.96 PE).
        // So silence is a real measurement -- "this channel saw < ~1 PE" --
        // and dropping it discards an upper limit the prediction must
        // respect, biasing the fit toward over-prediction.  The ~12% duty
        // cycle is an EFFECT of there being no light, not an independent
        // blindness.  Keep the channel: pair with pe_err_nodata so the
        // measured 0 carries the threshold band rather than a tight error.
        // Default TRUE = the legacy drop => bit-identical for existing cfgs.
        bool m_coverage_mask_fit{true};

        // PE_err floor (in PE) applied to readout-uncovered channels when
        // they stay in the fit (coverage_mask_fit false).  Their measurement
        // is the upper limit "PE < self-trigger threshold", not "PE == 0 +-
        // pe_err_floor": at the 0.3 PE default floor a 5 PE over-prediction
        // on a no-data channel costs the LASSO a pred/0.3 ~ 17 residual and
        // over-punishes the true match just as the legacy path did.
        //
        // SCOPE: this sets Opflash::PE_err, so it always reaches the LASSO
        // (pe_err = sqrt(PE + PE_err^2)).  It reaches the bundle chi2/KS ONLY
        // when pe_err_on_pred is false -- with it true, per_opdet_perr
        // (TimingTPCBundle.cxx:39-45) derives the error from the PREDICTED pe
        // and ignores the measured one entirely.
        //
        // <= 0 => disabled (bit-identical); needs the coverage chain (cov
        // empty => no-op).  PDVD leaves it unset: pe_err_floor is already 2.0
        // (with pe_err_knee at the default 1.0) so the LASSO sees +-2 PE,
        // ~2x the ~1 PE measured threshold, and pe_err_on_pred makes the chi2
        // price a no-data channel gently on its own.  See doc 14 section 12.
        double m_pe_err_nodata{-1.0};

        // §F bundle-quality thresholds (forwarded to TimingTPCBundle).
        double m_bundle_ks_merge_max{0.2};
        double m_bundle_chi2ndf_merge_max{20};
        double m_bundle_addmerge_exponent{0.8};
        double m_highconsist_ks_max{0.06};
        int    m_highconsist_min_ndf{3};
        double m_bundle_pe_ndf_knee{1.0};
        bool   m_bundle_mask_ks{false};  // apply opdet_mask to the KS shape metric too
        // Flag-aware "high-consistent" ladder (default off = single-branch, bit-identical; SBND-on).
        bool   m_highconsist_ladder{false};
        double m_hc_clean_ks{0.06};  double m_hc_clean_c2{6.0};
        double m_hc_good_ks{0.09};   double m_hc_good_c2{4.0};
        double m_hc_tb_ks{0.10};     double m_hc_tb_c2{8.0};
        double m_hc_miss_ks{0.08};   double m_hc_miss_c2{60.0};
        int    m_hc_miss_min_ndf{5};
        // Per-bundle chi2 relaxation (prototype close-to-PMT denom inflation +
        // one-inefficient-PMT subtraction); default off = bit-identical, SBND-on.
        bool   m_chi2_relax{false};
        double m_chi2_pmt_excess{350.0};
        double m_chi2_pmt_ratio{1.3};
        double m_chi2_pmt_inflate{0.5};
        // Extra chi2 denominator width on a rail-flagged channel, in the same
        // form as chi2_pmt_inflate: denom += (pe * chi2_sat_inflate)^2.  Unlike
        // the close_to_PMT widening this is gated on the sat flag ALONE, with no
        // excess/ratio test -- those fire on measured >> predicted (near-PMT
        // over-response), the OPPOSITE regime from a clipped channel, so reusing
        // them here would be a silent no-op.
        //
        // Independent of chi2_relax (a railed channel need not be near a PD).
        // Only bites when saturation_mask_fit is false (a masked channel never
        // reaches the chi2 loop).  Default 0 = no extra width => bit-identical.
        double m_chi2_sat_inflate{0.0};

        // §H raw readout-window truncation flag (T0-independent, APA-agnostic),
        // always computed. Flags a bundle whose cluster's leading/trailing time
        // slice sits within m_window_edge_ticks of the raw readout window
        // [0, m_readout_window_ticks] (consumed by the ql_scan "wtrunc" label and
        // the high-consistent ladder's miss branch).
        //
        // m_readout_window_ticks is the post-resample frame length in raw ticks.
        // Detector-dependent: the C++ default is the SBND daq.nticks; PDHD's window
        // is longer (5999), so its run script reads the real value from the SP frame
        // and passes it via -S readout_window_ticks (see run_clus_evt.sh) rather
        // than hard-coding it.
        int m_readout_window_ticks{3427};
        int m_window_edge_ticks{4};        // edge-proximity threshold (~one slice = nticks_live_slice)

        // When true, discard any (flash, cluster) bundle whose cluster is not
        // contained in the TPC drift box once the flash-derived T0 x-offset is
        // applied (the prototype flag_good_bundle gate, ToyMatching.cxx 272-275).
        // Default OFF so existing production configs stay bit-identical; SBND
        // jsonnet sets require_containment: true.
        bool m_require_containment{false};

        // Opaque-cathode "mismatched candidate" filter (PDHD). A candidate bundle that
        // pairs a cluster with a flash lit on the OPPOSITE drift side is non-physical:
        // with an opaque cathode a +x cluster cannot be the source of -x light (and
        // vice-versa), so the cluster's predicted light has ~no overlap with the flash's
        // measured light. The ONE physical exception is a genuine CATHODE-CROSSER: at the
        // flash's T0 the cluster's end reaches the cathode (flag_at_x_boundary), so its
        // far half could be the source of the other side's light (whose own-side flash
        // may be missing/dark). So: keep same-side and cross-side cathode-crossers, drop
        // every other cross-side bundle. (Brightness is irrelevant -- a bright crosser is
        // as valid as a dim one; a dim mid-drift cluster contained at some opposite
        // flash's T0 is just a coincidence.) Applied in build_bundles (fit candidate set)
        // and again in dump_calib so the hand-scan tables show the same universe. OFF by
        // default (no-op, bit-identical); PDHD jsonnet turns it on.
        bool m_cross_side_filter{false};

        // Perf only: apply the cross-side mismatch drop BEFORE the per-point visibility
        // loop in build_bundles instead of after it. The drop verdict (flash lit-side vs
        // the run's fixed cluster side, plus the pre-loop at_x_boundary flag) never uses
        // the predicted light, so hoisting it skips the SemiAnalyticalModel evaluation of
        // bundles that the post-loop test discards anyway -- same surviving candidate set,
        // bit-identical matching, lower wall time (the now-dominant cost on busy PDHD
        // events: vis_loop ~90% of QLMatching). Only meaningful with m_cross_side_filter
        // on. Default OFF => the legacy post-loop drop runs and the result is byte-for-byte
        // the old path (validated same-binary on PDHD: matching identical with skip on/off).
        bool m_crossside_skip_vis{false};

        // Approximation knob (default OFF): coarsen the per-point visibility loop of
        // build_bundles for LARGE cluster groups. When vis_sample_stride > 1 and the
        // group's total point count is >= vis_sample_min_pts, every stride-th point of
        // each blob is evaluated and weighted by the number of points it stands in for
        // (min(stride, npoints - i)), so each blob's total charge entering pred_flash
        // is conserved exactly; the approximation is that the visibility is taken
        // piecewise-constant over stride-long runs of the blob's point order. At the
        // default stride 1 the loop body performs the identical FP operations
        // (bit-identical output). Changes pred_flash on coarsened groups =>
        // result-changing => jsonnet-togglable and A/B-validated before any enable.
        int m_vis_sample_stride{1};
        int m_vis_sample_min_pts{10000};

        // Light-pattern over-prediction prefilter (the prototype's fired-fraction
        // reject, FlashTPCBundle.cxx 547-602). Drop a (flash, cluster) bundle BEFORE
        // the chi2 fit when its predicted light is much larger than the measured
        // light (a cluster emitting far more light than the flash shows cannot be the
        // match). One-directional: only over-prediction is rejected; under-prediction
        // is physically fine (a flash may be lit by several clusters). Two metrics,
        // both over the same opdet_mask the chi2 uses:
        //   R_total = sum(pred_masked) / sum(meas_masked)
        //   R_max   = pred[ch*] / max(meas[ch*], 1),  ch* = argmax(pred_masked)
        // Reject if R_total > overpred_total_ratio OR R_max > overpred_maxch_ratio.
        // Boundary/truncated bundles (close_to_PMT | window_truncated | at_x_boundary)
        // are EXEMPT (measured is an underestimate there). Default OFF (large ratios
        // = inert) so production stays bit-identical; SBND jsonnet sets the tuned
        // values (tuned on data hand-scans, validated on MC; see
        // sbnd_xin/ql_prefilter_tune.py).
        bool m_reject_overpred{false};
        double m_overpred_total_ratio{1e9};   // R_total ceiling (inert when huge)
        double m_overpred_maxch_ratio{1e9};   // R_max   ceiling (inert when huge)
        // Channel scope for the two ratios: when non-empty, R_total/R_max are
        // computed over (overpred_channels ∩ opdet_mask) instead of the full
        // masked set.  Lets the prefilter judge on one trusted PD family only
        // (PDVD: cathode XAs — the wall XAs are bimodal and the PMTs
        // self-trigger-silent ~46% of the time, so their meas~0 would fake
        // over-prediction on a good bundle).  Empty (default) => legacy full
        // masked set => bit-identical.
        std::vector<int> m_overpred_channels;
        std::vector<char> m_overpred_ch_sel;  // per-channel flag form, built at configure

        // §I empty-flash light-quality rescue (the prototype's flash-centric pick,
        // ToyMatching.cxx organize_matched_bundles). The LASSO selects by strength,
        // which can leave a flash EMPTY (no bundle above m_strength_cutoff) even when
        // a cluster is an excellent LIGHT match for it — the cluster was won by a
        // neighbouring flash on strength alone. After the fit, for each empty flash,
        // adopt its best light-quality candidate (metric = ks*(chi2/ndf)^exp, with a
        // boundary/near-PMT down-weight) from the pre-cutoff universe if the metric
        // clears m_rescue_metric_max. If that cluster is already matched elsewhere,
        // reassign it ONLY when the empty flash is a strictly better light match
        // (guard: it never removes a better match, so it cannot drop a correct pair).
        // Default OFF (huge bar = inert) so production stays bit-identical; SBND-on.
        bool   m_empty_rescue{false};
        double m_rescue_metric_max{1e9};       // light-quality bar (inert when huge)
        double m_rescue_exponent{0.8};         // chi2/ndf exponent (prototype 0.8)
        double m_rescue_boundary_weight{0.8};  // per-flag down-weight (prototype 0.8/0.64)

        // §J cluster-centric rescue (the complement of §I's flash-centric pick). The
        // empty-flash rescue only fills flashes left WHOLLY empty; a big cluster with an
        // excellent candidate (low ks, chi2/ndf~1-2, pred/meas~1) but driven to strength
        // 0 by the LASSO L1 sparsity — a rival already explains its flash — stays
        // unmatched even when that flash is non-empty. After the empty-flash rescue, for
        // each cluster STILL unmatched, adopt its best accepted candidate from the
        // pre-cutoff snapshot, attaching it even onto an already-non-empty flash
        // (many-clusters-per-flash is physical and GT-endorsed). Acceptance bar:
        //   ks < ks_max  AND  chi2/ndf < chi2ndf_max  AND  ratio_lo < pred/meas < ratio_hi
        // (pred = total predicted light, meas = flash total measured PE). Default
        // ks_max=0 => the gate ks<0 is vacuously false => no-op (doubly inert with the
        // bool guard) so production stays bit-identical; PDHD jsonnet turns it on.
        bool   m_cluster_rescue{false};
        // Shared-flash-aware rescue variants (doc 19 phase 5; PDVD).  The legacy
        // rescues above are per-run "empty flash" concepts and are SKIPPED under
        // shared_flash; these run in the JOINT rounds instead.  empty_rescue_shared:
        // a flash is empty only when NO drift side has a surviving bundle for it;
        // the best snapshot candidate ACROSS sides is adopted (same metric/
        // reassign-only-if-strictly-better/pin-locked rules as §I, same
        // rescue_metric_max bar).  cluster_rescue_shared: per-run cluster-centric
        // adoption exactly as §J (ADD-only; a shared flash legitimately holds
        // bundles of both sides), gated by the same cluster_rescue_* thresholds.
        // Both default OFF => byte-identical.
        bool   m_empty_rescue_shared{false};
        bool   m_cluster_rescue_shared{false};
        double m_cluster_rescue_ks_max{0.0};        // 0 => gate false => inert
        double m_cluster_rescue_chi2ndf_max{0.0};
        double m_cluster_rescue_ratio_lo{0.0};
        double m_cluster_rescue_ratio_hi{0.0};
        // Draw the cluster-rescue candidate pool from the PRE-cull universe
        // (run.all_bundles minus build-prefilter rejects) instead of the post-cull
        // prefit_snapshot. cull_inconsistent drops a cluster's NON-consistent rivals
        // when it has a consistent bundle; that can remove the cluster's only good
        // LIGHT match (e.g. evt983 uid50 flash127: ks 0.108, pred/meas 1.45, free flash)
        // before it ever reaches the fit OR the snapshot, so the post-cull rescue can't
        // see it. With this on, the rescue can re-adopt such cull-victims -- still gated
        // by the same PE-scale bar, so purity stays bounded. Default OFF (snapshot pool,
        // = the shipped behaviour); PDHD-on.
        bool   m_cluster_rescue_precull{false};

        // Bee-op flash gid: when a single global optical flash list is fed to
        // multiple per-side QLMatching nodes (PDHD all-PD light: both drift sides
        // read the same opflash archive), EACH node's write_opflash_pc emits the
        // FULL list keyed by its own anode ident, so every physical flash lands on
        // the merged root TWICE and the Bee op display shows each flash doubled.
        // With this set, the gid (and the cluster's matched_flash_gid) is keyed by
        // the flash's PHYSICAL side instead of the node ident, so both nodes emit
        // the same gid for one flash -> fill_bee_flashes collapses the duplicate and
        // a cross-side (xTPC) matched_flash_gid still resolves to it. Default OFF =>
        // legacy node-ident gid (byte-identical for SBND per-TPC flashes / single-
        // node configs, whose per-node flash lists are already one-per-side).
        bool   m_opflash_phys_gid{false};

        // Per-PMT non-linearity correction applied to the predicted PE total (study-grade,
        // scintillation-profile-dependent; see sbnd_xin/pmt_nonlinearity_curve.py and
        // match/docs/sbnd-opdetsim-chain.md). Maps the true predicted PE p on each PMT into
        // the saturated (observed) space so it matches the post-saturation measured op_pes:
        //   p' = p                                      for p <= knee
        //      = knee * exp(beta*L + gamma*L^2), L=ln(p/knee)   for p > knee
        // beta=1, gamma=0 (or empty arrays) => identity. Default OFF so production configs
        // stay bit-identical; the sbnd_xin standalone run config enables it.
        bool m_pmt_nonlinearity{false};
        double m_pmt_nl_knee{700.0};
        std::vector<double> m_pmt_nl_beta;   // per-OpDet; empty or 1.0 => identity
        std::vector<double> m_pmt_nl_gamma;  // per-OpDet; empty or 0.0 => identity

        // Default SBND VUV/VIS efficiency arrays, indexed by OpDet (312 entries).
        // Configuration "VUVEfficiency"/"VISEfficiency" arrays override these.
        std::vector<double> m_VUVEfficiency{
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752
        };
        std::vector<double> m_VISEfficiency{
            0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271
        };

        IAnodePlane::pointer m_anode{nullptr};
        // Anodes this node matches. The single-anode path (config "anode") yields
        // one entry and multiplicity 1, so the node has one input port and behaves
        // exactly as the historical single-input matcher. A list (config "anodes")
        // enables the joint multi-APA path: one input port per anode, each matched
        // independently in its own ApaRun, then the trees are merged.
        std::vector<IAnodePlane::pointer> m_anodes;
        std::size_t m_multiplicity{1};
        // Anodes registered on each run's Grouping for wpid->anode lookup
        // (Cluster::get_length etc.). Empty (default) => just the run's own anode,
        // bit-identical to the historical single-APA path. When one input tree
        // spans several anodes of a shared drift volume (PDHD drift-side group:
        // APA0+APA2 / APA1+APA3 merged before the all-TPC stage), list them all
        // here so blobs from every APA resolve; the run's drift geometry still
        // uses the single representative "anode".
        //
        // Stored PER INPUT (sized to m_multiplicity). Config "grouping_anodes" may be
        // a flat array (the same list for every input -- single-side node) or an array
        // of arrays (one list per input port -- JOINT multi-side node, so each run's
        // active-volume box unions only its OWN side's APAs). An input with no list
        // falls back to {run.anode}.
        std::vector<std::vector<IAnodePlane::pointer>> m_grouping_anodes;
        // Root-node local PCs concatenated across inputs when merging the per-APA
        // trees (multi-APA path only); mirrors PointTreeMerging. ['opflash'] = the
        // optical-flash display PC; everything else is dropped from non-primary roots.
        std::set<std::string> m_root_pcs_to_merge{"opflash"};
        IDetectorVolumes::pointer m_dv;

        // Optional SBND CPA structure-exclusion fiducial volume. When set (the
        // "cathode_fiducial" config tn), the cathode-end flag_at_x_boundary is
        // decided by a 3-D point-in-volume test against it; when null (default,
        // non-SBND) the cathode-end flag falls back to the original flat-cathode
        // 1-D drift-coordinate window, leaving those configs bit-identical.
        IFiducial::pointer m_cathode_fv{nullptr};

        // Per-cluster cache of the cluster's significant extreme points
        // (get_extreme_wcps, flattened). Used by the cathode-end fiducial test; the
        // geometric extremes are intrinsic to the cluster (offset-independent for
        // ranking) so they are computed once and reused across candidate flashes.
        // Cleared each event in operator().
        mutable std::unordered_map<const WireCell::Clus::Facade::Cluster*,
                                   std::vector<WireCell::Point>> m_extreme_cache;

        // Per-cluster cache of the two main-PCA-axis extreme endpoints
        // (get_extreme_wcps()[0][0] / [1][0]; high/low projection). Geometric and
        // flash-independent, so the O(npoints) scan runs once per cluster and is
        // reused across candidate flashes by compute_two_boundary_flag. Cleared
        // each event in operator() alongside m_extreme_cache.
        mutable std::unordered_map<const WireCell::Clus::Facade::Cluster*,
                                   std::pair<WireCell::Point, WireCell::Point>> m_pca_endpoints_cache;

        // Per-cluster wall-proximity verdict cache (vd_surface_flags only; a cluster
        // belongs to exactly one input port, so keying by cluster is unambiguous).
        // Cleared each event in operator() alongside m_extreme_cache.
        mutable std::unordered_map<const WireCell::Clus::Facade::Cluster*, uint8_t>
            m_wall_flag_cache;

        std::string m_inpath{"pointtrees/%d"};
        std::string m_outpath{"pointtrees/%d"};
        float m_cluster_t0{-1e12};

        // Hand-scan calibration dump. When non-empty, after all per-APA runs
        // complete, write one per-event JSON holding the FULL candidate-bundle
        // universe (run.all_bundles) with predicted/measured light, metrics, flags,
        // cluster geometry and the detector box — the input for the off-line Q/L
        // hand-scan event display (sbnd_xin/ql_scan). Empty (default) => the dump
        // method is never called, so production output is bit-identical. A "%" in
        // the path is String::format-ed with the node call count. Observation-only:
        // it reads already-populated state and never touches the matching path.
        std::string m_calib_dump{""};
        // ±80 ns flash-flash coincidence gap used to tag the dump's per-flash
        // "group" id (mirrors MultiAlgBlobClustering::store_flash_groups /
        // SBND clus.jsonnet flash_group_window), so the viewer's both-TPC
        // coincidence filter matches the Bee display. Calib-dump only.
        double m_flash_group_window{80 * units::ns};

        // Cathode-crossing TPC0/TPC1 offset diagnostic (off unless cathode_diag is
        // set). Pairs a TPC0 + TPC1 cluster sharing one T0 flash group, both reaching
        // the cathode, finds their closest pair of points in the T0-corrected frame,
        // and logs the two local Hough directions plus the connecting vector so the
        // drift-x gap (t0/velocity, degenerate) can be separated from a transverse
        // y-z position shift. Empty (default) => never invoked, production unchanged.
        // Observation-only: reads finished runs, never touches the matching path.
        std::string m_cathode_diag{""};
        // Radius for the local Hough direction (vhough_transform) at each cathode-end
        // point. Cathode-diag only.
        double m_cathode_diag_radius{15 * units::cm};

        // Cross-TPC cathode-crossing consistency confirm-stamp (off unless xtpc_flag).
        // Post-matching: for each coincident matched-main-cluster pair across the two
        // TPCs, set flag_xtpc_consistent if the two halves connect as one track across
        // the cathode -- closest approach < xtpc_dmax (scenario 1), or, when a half is
        // window-truncated, the connecting vector is collinear with both local Hough
        // directions to within xtpc_angle_max (scenario 2). Cuts tuned on the 10 SBND
        // hand-scan data events. Default off => never invoked, output bit-identical.
        // Observation-only: matching assignments unchanged, only adds the flag.
        bool   m_xtpc_flag{false};
        double m_xtpc_dmax{5 * units::cm};
        double m_xtpc_dmax2{300 * units::cm};     // scenario-2 closest-approach ceiling
        double m_xtpc_angle_max{20.0};            // degrees
        double m_xtpc_hough_radius{15 * units::cm};

        // Cross-TPC JOINT-PIN (off unless xtpc_joint_pin; requires xtpc_flag + 2 sides).
        // Strengthens the scenario-1 priority cull: instead of each half keeping all its
        // scenario-1 bundles independently (which lets a half drift to whichever flash the
        // light prefers, splitting a true crosser across two flashes), a DIRECTION-CONFIRMED
        // pair is pinned to ONE coincident flash -- both halves keep only that bundle. The
        // confirmation adds a track-axis collinearity test to scenario-1's d<xtpc_dmax: the
        // two halves' axes must agree within m_xtpc_pin_angle, taken as the OR (min angle) of
        // the LOCAL vhough direction (at the cathode-end closest point) and the GLOBAL cluster
        // PCA axis. Neither alone is reliable -- local vhough reads a spurious end-curl kink on
        // straight crossers (PDHD run29107: vhough 37/43deg but global 2/1.7deg), global PCA is
        // meaningless on bent/messy clusters; the OR passes a genuine crosser via either. The
        // connecting-vector angles are NOT used (a ~1cm inter-TPC transverse shift makes them
        // ~perpendicular even for true crossers). Pairs greedily by smallest closest-approach d
        // so each cluster is pinned once, to the flash where its halves actually meet at the
        // cathode. Default off => flag_xtpc_pin never set => bit-identical.
        bool   m_xtpc_joint_pin{false};
        double m_xtpc_pin_angle{20.0};            // degrees, folded track-axis collinearity

        // Cross-TPC CATHODE RESCUE (off unless xtpc_cathode_tol > 0; needs xtpc_flag +
        // 2 sides to ever confirm anything). Motivation: a genuine cathode crosser
        // whose two halves meet a few cm off the nominal cathode plane fails BOTH
        // admission gates at once -- the overshooting half fails containment
        // (require_containment drops the bundle at build), and the short half misses
        // the at_cathode window so it never gets flag_at_x_boundary and never enters
        // cull_cross_tpc's candidate set. The pair the xtpc joint-pin exists for is
        // then invisible to it (PDVD 039252 evt298567 top-97 + bot-139 vs flash 96:
        // ends agree to 0.4 cm but meet ~3.5 cm below the cathode face; see
        // wcp-porting-img/pdvd/docs/qlmatch/16_pdvd-clus97-crosser-evt298567.md §10).
        // Mechanism, all gated on m_xtpc_cathode_tol > 0 (flat-cathode window only,
        // never when m_cathode_fv is configured):
        //  - a bundle failing containment ONLY by cathode overshoot is kept
        //    PROVISIONALLY when its junk-tolerant cathode endpoint (the deepest slice
        //    that cannot be discarded within a charge budget of m_xtpc_cathode_qfrac
        //    of the cluster charge, walking from the deep end; qfrac<=0 => the raw
        //    endpoint) lies within m_xtpc_cathode_tol past the containment gate;
        //  - a CONTAINED bundle whose cathode end sits within m_xtpc_cathode_tol
        //    BELOW the at_cathode window gets xtpc candidate admission only
        //    (flag_xtpc_cathode_cand, NOT flag_at_x_boundary -- that flag also feeds
        //    the ladder/cross-side/LASSO-weight paths, which stay legacy);
        //  - cull_cross_tpc admits candidates by the new flag, but never confirms two
        //    provisional (uncontained) halves with each other;
        //  - after cull_cross_tpc, provisional bundles WITHOUT a scenario-1
        //    confirmation are purged (purge_unconfirmed_cathode_rescue), BEFORE
        //    cull_inconsistent/fit -- downstream then sees exactly the legacy set.
        // So the tolerance is only ever exercised for a pair the existing scenario-1 /
        // joint-pin geometry independently confirms (opposite-volume partner on the
        // SAME flash, clouds meeting within xtpc_dmax at the cathode). Default 0 =>
        // no flag is ever set, no bundle is kept or removed => byte-identical.
        double m_xtpc_cathode_tol{0.0};           // length; 0 = off
        double m_xtpc_cathode_qfrac{0.0};         // charge fraction discardable as junk

        // -------- xtpc / selection quality gates (039252 scan tuning, doc 19) ----------
        // All default to the legacy behaviour => byte-identical when unset.
        //
        // Pin strength floor: a pinned crosser half is exempt from the LASSO
        // strength-cutoff prune ("the pin ignores light by design").  The 18-evt
        // hand/AI scan showed that exemption is the pin path's phantom source:
        // phantom pins have median strength 0.00 while agreed pins sit at p10
        // 0.88.  When > 0, a pinned bundle whose solution <= this floor LOSES the
        // exemption (pruned like any bundle; its cluster can then be rescued or
        // stay unmatched).  0 (default) = pins always exempt.
        double m_xtpc_pin_min_strength{0.0};
        // Scenario-1/xtpc-consistent light gate: legacy sets the xtpc flags on
        // both halves from GEOMETRY alone, at every coincident flash whose T0
        // offset makes the halves touch -- junk bundles then monopolize
        // cull_inconsistent's xtpc tiers (scan: sc1 phantoms ks p50 0.40 /
        // chi2/ndf p50 76 vs agrees 0.17 / 1.3).  When on, a bundle acquires
        // flag_xtpc_consistent / flag_xtpc_scenario1 only if its own light
        // passes (ks <= sc1_ks_max AND chi2/ndf <= sc1_c2n_max).  The joint-pin
        // candidacy itself is untouched (see pin_min_strength above).
        bool   m_xtpc_sc1_light_gate{false};
        double m_xtpc_sc1_ks_max{0.3};
        double m_xtpc_sc1_c2n_max{50.0};
        // Cathode-rescue ks ceiling: purge_unconfirmed_cathode_rescue keeps a
        // provisional overshoot bundle on scenario-1 confirmation alone; scan
        // shows those survivors at 69% phantom (ks p50 0.435 vs agrees 0.20).
        // When > 0 a confirmed-but-dim bundle (ks > ceiling) is purged too.
        double m_xtpc_cathode_ks_max{0.0};
        // Post-fit cull of UNFLAGGED low-quality selections: a selected bundle
        // carrying no quality flag (not consistent / xtpc_consistent /
        // scenario1 / pin) survived on LASSO strength alone -- the largest
        // phantom bucket (103/147 = 70% phantom; ks p50 0.35, chi2/ndf p50 27
        // vs agrees 0.18 / 3.3).  When on, such a bundle is removed from the
        // matched output if ks > postcull_ks_max OR chi2/ndf > postcull_c2n_max
        // (cut 0.30/20 kills 72/92 flagged-phantoms at a cost of 8/40 agrees).
        bool   m_postcull_unflagged{false};
        double m_postcull_ks_max{0.30};
        double m_postcull_c2n_max{20.0};

        // Path to the JSON file holding VUVHits, VISHits, geometry and the
        // SBND OpDet array.
        std::string m_semimodel_file{"sbnd/photodet/semi-analytical-sbnd.json"};

        std::unique_ptr<SemiAnalyticalModel> m_semi_model;

        // Visibility backend selector: "semi" (default, the semi-analytical
        // model) or "library" (gridded lookup sampled offline from a detector
        // optical model, e.g. the PDVD PDFastSimANN computable graph; see
        // PhotonLibraryModel.h).  In library mode semimodel_file is still
        // loaded (it supplies the OpDet table used by masks and cathode-side
        // logic) but the per-point visibilities come from
        // photon_library_file, with reflected light 0 (the library holds
        // total photon arrival) and no same-TPC x-sign gate (the library
        // encodes cathode opacity and double-sided cathode PDs).  Default
        // "semi" leaves the existing path untouched (bit-identical).
        std::string m_light_model{"semi"};
        std::string m_photon_library_file{""};
        std::unique_ptr<PhotonLibraryModel> m_lib_model;
        // Cross-check tolerance (cm) between the library's own optional
        // per-channel chan_pos_cm and m_opdets[i].center, at configure time.
        // Only exercised when the library file actually carries chan_pos_cm
        // (older/hand-written files without it skip the check silently, so
        // this is not a new requirement on existing files) -- catches a
        // channel-order mismatch between the library export and the OpDet
        // table that the bare nchan==nopdets count check cannot.
        double m_photon_library_pos_tol{5.0};

        // Per-OpDet table (type + position) loaded from semimodel_file. Used to
        // build the SemiAnalyticalModel and, at execute() time, to derive the
        // per-channel on/off mask and the per-TPC OpDet split.
        std::vector<SemiAnalyticalModel::OpticalDetector> m_opdets;
        double m_cathode_x{0.0};  // cathode-plane x; OpDets on its low/high side belong to TPC 0/1

        // True iff this (flash, cluster) bundle should be dropped by the opaque-cathode
        // mismatched-candidate filter: the flash is lit on the opposite drift side from
        // the cluster AND the cluster is NOT a cathode-crosser at this flash's T0
        // (at_x_boundary == false), so it cannot physically be the far half of the track
        // that lit the flash. The flash's lit side is taken from where its measured PE
        // actually is (OpDet x vs the cathode plane), matching the opflash-PC tagging.
        // Always false when the filter is off (m_cross_side_filter == false), so the
        // default path is bit-identical.
        bool cross_side_mismatch_drop(const Opflash* flash, int cluster_side,
                                      bool at_x_boundary) const;

        // Significant extreme points of a cluster (get_extreme_wcps, flattened),
        // memoized in m_extreme_cache. Used to locate the cathode endpoint.
        const std::vector<WireCell::Point>&
            cluster_extreme_points(WireCell::Clus::Facade::Cluster* cluster) const;

        void remove_bundle_selection(TimingTPCBundleSelection to_be_removed,
                                     TimingTPCBundleSet& bundle_set);
        void remove_bundle_selection(
            TimingTPCBundleSelection to_be_removed,
            FlashBundlesMap& flash_bundles_map,
            ClusterBundlesMap& cluster_bundles_map,
            std::map<std::pair<Opflash*, WireCell::Clus::Facade::Cluster*>,
                     TimingTPCBundle::pointer>& flash_cluster_bundles_map);
        void organize_bundles(
            TimingTPCBundleSelection& results_bundles,
            std::map<std::pair<Opflash*, WireCell::Clus::Facade::Cluster*>,
                     TimingTPCBundle::pointer>& flash_cluster_bundles_map);

        // ---- Per-APA run state (operator() refactor) ----
        // Every container the matching pipeline threads between stages lives here,
        // not as a QLMatching member, so each APA is processed in a fresh, isolated
        // ApaRun. That isolation keeps the pointer-keyed map iteration deterministic
        // when more than one APA is processed in one node (joint path): no cross-APA
        // pointer ordering can leak in. operator() builds one ApaRun per input.
        struct ApaRun {
            // input / identity
            IAnodePlane::pointer anode;
            std::unique_ptr<WireCell::PointCloud::Tree::Points::node_t> root_live;
            WireCell::Clus::Facade::Grouping* grouping{nullptr};
            std::string inpath;
            int charge_ident{0};
            std::size_t input_idx{0};   // this run's input-port index (k), for per-input config

            // optical-detector mask / kept-channel index
            std::vector<unsigned int> opdet_mask;
            std::set<int> auto_masked;          // channels dropped by the dynamic auto-mask
            std::vector<Opflash::pointer> flashes;
            unsigned int nopdet{0};
            std::vector<int> opdet_idx_v;

            // per-TPC drift geometry (computed once per anode)
            int sign_offset{1};
            double s{1.0}, anode_x{0.0}, u_cathode{0.0};
            double y_lo{0.0}, y_hi{0.0}, z_lo{0.0}, z_hi{0.0};
            // per-TPC transverse position offset (Y,Z), read from DetectorVolumes
            // metadata "pos_offset" in compute_geometry().  Parked for the future
            // cross-TPC matching judgement; not yet consumed.  0 => inert.
            double dy{0.0}, dz{0.0};

            // cluster-group decomposition
            std::vector<std::pair<WireCell::Clus::Facade::Cluster*,
                                  std::vector<WireCell::Clus::Facade::Cluster*>>> match_groups;
            std::vector<WireCell::Clus::Facade::Cluster*> clusters;
            std::map<Opflash*, int> global_flash_idx_map;
            std::map<WireCell::Clus::Facade::Cluster*, int> global_cluster_idx_map;

            // bundles + the three lookup maps
            std::vector<TimingTPCBundle::pointer> all_bundles;
            TimingTPCBundleSet pre_bundles;
            FlashBundlesMap flash_bundles_map;
            ClusterBundlesMap cluster_bundles_map;
            std::map<std::pair<Opflash*, WireCell::Clus::Facade::Cluster*>,
                     TimingTPCBundle::pointer> flash_cluster_bundles_map;

            // Full pre-LASSO candidate universe (flash -> every candidate bundle),
            // captured at fit_round1 start before any strength prune. Only filled when
            // m_empty_rescue; consumed by rescue_empty_flashes after the round-2 prune
            // so it can reach the strength-0-but-light-good bundles the fit drops.
            FlashBundlesMap prefit_snapshot;

            BundleQualityParams qp;

            // charge bookkeeping (debug only)
            double total_charge_blob{0.0};
            double total_charge_point{0.0};
            double total_charge_blob_all{0.0};
        };

        // Run the full single-APA matching pipeline on one ApaRun.
        void run_one_apa(ApaRun& run);

        // The pipeline split at the LASSO-fit boundary, so a cross-TPC pre-fit cull can
        // run between them (m_xtpc_flag, multi-APA): _prefit builds bundles + the
        // per-TPC consistency cull; _fit runs the two LASSO rounds + output. run_one_apa
        // = _prefit then _fit (the historical single-APA order).
        void run_one_apa_prefit(ApaRun& run);
        void run_one_apa_fit(ApaRun& run);

        // Pipeline stages, extracted verbatim from the old operator().
        void build_opdet_mask(ApaRun& run);          // base OpDet on/off mask
        void read_flashes(ApaRun& run);              // canonical flash PCs -> Opflash
        void decompose_cluster_groups(ApaRun& run);  // main+associated split, idx maps
        void recompose_cluster_groups(ApaRun& run);  // undo the split before output (T0-annotation only)
        void compute_geometry(ApaRun& run);          // per-TPC drift geometry, mask cull, opdet idx
        void compute_dynamic_opdet_mask(ApaRun& run, unsigned int tpc);  // per-event dead-PMT auto-mask
        void build_bundles(ApaRun& run);             // (flash,group) bundles + predicted light  [Stage 1]

        // Fill the per-bundle boundary flags (flag_close_to_PMT,
        // flag_at_x_boundary, flag_spec_end, flag_at_cathode) for one
        // (flash, cluster) bundle. Ports the MicroBooNE prototype end-trimming
        // walk (ToyMatching.cxx ~176-290) into the per-TPC drift coordinate u,
        // using run.s / run.anode_x / run.u_cathode (computed in compute_geometry).
        // With m_vd_surface_flags, PD-surface awareness applies (see the knob doc).
        //
        // Returns the prototype's flag_good_bundle / TPC-containment verdict: true
        // iff the (post-trim) cluster endpoints fall inside the TPC drift box (the
        // 4-part in-window guard, == ToyMatching.cxx 272-275). Returns false when
        // the cluster has no 3d points at all (it cannot be contained). The
        // endpoint walk unions the cluster's slices over EVERY anode it touches (a
        // PDHD drift group spans two APAs offset in Z), so a cluster living wholly
        // in the group's second APA is no longer spuriously declared uncontained.
        // The caller discards uncontained bundles when m_require_containment.
        bool compute_endpoint_flags(TimingTPCBundle* bundle,
                                    WireCell::Clus::Facade::Cluster* cluster,
                                    double flash_x_offset,
                                    const ApaRun& run) const;

        // Wall-proximity verdict for the vd_surface_flags mode: bit 0 = the cluster
        // reaches within m_pd_wall_cushion of the low-y active edge (run.y_lo),
        // bit 1 = of the high-y edge (run.y_hi). y is drift-invariant so the verdict
        // is flash-independent; memoized per cluster (m_wall_flag_cache).
        uint8_t wall_proximity(WireCell::Clus::Facade::Cluster* cluster,
                               const ApaRun& run) const;

        // Set the diagnostic flag_two_boundary on one bundle: true iff the two
        // main-PCA extremes of the main cluster each lie within m_two_boundary_margin
        // of a per-APA active-volume face AND those are two SEPARATE faces (different
        // ones of the 6: anode/cathode/±y/±z) — the cluster enters through one
        // surface and exits through a different one. Each endpoint is assigned its
        // nearest face. flash_x_offset (drift X correction) shifts x only; y/z are
        // drift-invariant. The two endpoints are cached per cluster
        // (m_pca_endpoints_cache). Observation-only: never read by the matching path;
        // build_bundles calls it only when a calib dump is requested. (The cathode
        // face could later swap to the CPA structure-exclusion m_cathode_fv used by
        // flag_at_x_boundary; the uniform box-margin test is used for now.)
        void compute_two_boundary_flag(TimingTPCBundle* bundle,
                                       WireCell::Clus::Facade::Cluster* main_cluster,
                                       double flash_x_offset, const ApaRun& run) const;
        void build_bundle_maps(ApaRun& run);         // flash/cluster/pair maps + deterministic sort
        // xtpc cathode rescue resolution (m_xtpc_cathode_tol > 0 only): drop
        // provisional cathode-overshoot bundles that cull_cross_tpc did not confirm
        // as scenario-1 crosser halves; stamp survivors contained. Runs between
        // cull_cross_tpc and cull_inconsistent.
        void purge_unconfirmed_cathode_rescue(ApaRun& run);
        void cull_inconsistent(ApaRun& run);         // drop non-consistent rivals               [Stage 1]
        // Pin exemption from the strength-cutoff prune, with the optional
        // m_xtpc_pin_min_strength floor (doc 19); 0 => legacy always-exempt.
        bool pin_exempt(const TimingTPCBundle* b, double strength) const;
        // Post-fit cull of unflagged low-quality selections (m_postcull_unflagged).
        void cull_unflagged_lowquality(ApaRun& run);
        void fit_round1(ApaRun& run);                // LASSO, per-flash background DOF           [Stage 2]
        void fit_round2(ApaRun& run);                // LASSO + KS-shape, keep best per cluster   [Stage 3]

        // Shared-flash JOINT LASSO rounds (m_shared_flash, >1 input port; PDVD).
        // Mirror fit_round1/fit_round2 but build ONE system over the physical
        // flashes (keyed by flash row id, identical across ports) with columns from
        // every port's candidate bundles, so both drift sides share each flash's
        // explanation. Prunes route to each bundle's owning run's maps. The per-run
        // fit_roundN code paths are untouched (knob off => byte-identical).
        void fit_round1_shared(std::vector<ApaRun>& runs);
        void fit_round2_shared(std::vector<ApaRun>& runs);

        // Per-column L1 down-weight (m_lasso_boundary_weight) for a boundary / near-PMT /
        // window-truncated bundle when m_lasso_flag_weight is on (prototype's flag-based
        // weight), else 1.0. Multiplies the pe-mismatch (+KS) base in both rounds.
        double lasso_flag_factor(const TimingTPCBundle::pointer& bundle) const;

        // Empty-flash light-quality rescue (m_empty_rescue; see §I). snapshot is the
        // full pre-strength-cutoff flash->candidate map; this mutates run.flash_bundles_map
        // in place, adopting the best light-quality candidate of each emptied flash.
        void rescue_empty_flashes(ApaRun& run, const FlashBundlesMap& snapshot);

        // Cluster-centric rescue (m_cluster_rescue; see §J). After the empty-flash
        // rescue, adopt the best accepted candidate for each cluster the LASSO left
        // UNMATCHED (its flash was won by a rival under L1 sparsity), attaching it even
        // onto an already-non-empty flash. snapshot is the same pre-strength-cutoff
        // flash->candidate map; this mutates run.flash_bundles_map in place.
        void rescue_unmatched_clusters(ApaRun& run, const FlashBundlesMap& snapshot);
        // Shared-flash empty-flash rescue (m_empty_rescue_shared): joint emptiness
        // across all runs, best candidate across sides; mutates the owning run's
        // flash_bundles_map.
        void rescue_empty_flashes_shared(std::vector<ApaRun>& runs);
        void apply_matched_t0s(ApaRun& run);         // write cluster t0 / flash / matched gid
        void write_opflash_pc(ApaRun& run);          // merge-safe per-root "opflash" PC

        // Physical drift side (0 = low-x, 1 = high-x of the cathode) of a flash,
        // from where its measured light actually sits. When m_opflash_phys_gid is
        // set this seeds the Bee-op flash gid (and matched_flash_gid) instead of the
        // processing node's anode ident, so a flash that two per-side nodes both
        // dump from one shared global flash list gets ONE gid and collapses in the
        // Bee op display (see write_opflash_pc).
        int  flash_phys_side(const Opflash* flash) const;

        // Hand-scan calibration dump (m_calib_dump only). Builds one per-event JSON
        // from the finished per-APA runs (the full all_bundles candidate set across
        // both TPCs) and writes it via WireCell::Persist::dump. Never touches the
        // matching path.
        void dump_calib(const std::vector<ApaRun>& runs);

        // Cathode-crossing offset diagnostic (m_cathode_diag only). Logs, per
        // cross-TPC cathode-crossing pair, the two local directions and the
        // connecting vector of the closest point pair. Never touches the matching.
        void dump_cathode_diag(const std::vector<ApaRun>& runs);

        // Cross-TPC cathode-crossing pre-fit cull (m_xtpc_flag, multi-APA). Runs after
        // each TPC's bundles + per-TPC cull but BEFORE the LASSO: pairs candidate
        // main-cluster bundles across the two TPCs whose flashes coincide, sets
        // flag_xtpc_consistent on a geometrically cross-TPC-consistent pair (the two
        // halves connect as one track across the cathode), then drops each marked
        // cluster's non-consistent rivals so fewer bundles enter the fit. Off by default
        // => not called => bit-identical.
        void cull_cross_tpc(std::vector<ApaRun>& runs);

        // One candidate cross-TPC pair test, shared by cull_cross_tpc. m{0,1} carry the
        // two clusters with their T0 x-offset, per-TPC (y,z) pos_offset, and truncation
        // flag. Returns the scenario code: 1 = scenario 1 (closest approach < m_xtpc_dmax,
        // tight/self-vetoing), 2 = scenario 2 (a half truncated AND conn,dir0,dir1 mutually
        // collinear < m_xtpc_angle_max AND d < m_xtpc_dmax2), 0 = not consistent.
        struct XtpcMC {
            TimingTPCBundle* b;
            WireCell::Clus::Facade::Cluster* c;
            double off, dy, dz;
            bool wt;
        };
        // Raw-coordinate (T0-independent) bounding boxes of one cluster's 3D points:
        // the whole-cluster box plus per-chunk boxes over runs of XTPC_CHUNK
        // consecutive points, used by xtpc_pair_consistent to prune its brute-force
        // closest-approach loop against the running best distance. Pruning only
        // skips point pairs that cannot strictly improve the minimum, and the
        // iteration order of the surviving pairs is the legacy row-major order, so
        // the closest pair (value AND argmin/tie-break) is bit-identical.
        static constexpr int XTPC_CHUNK = 256;
        struct XtpcBoxes {
            std::array<double, 6> whole;                 // {xlo,xhi,ylo,yhi,zlo,zhi}
            std::vector<std::array<double, 6>> chunks;   // per XTPC_CHUNK-point run
        };
        static XtpcBoxes xtpc_boxes(const WireCell::Clus::Facade::Cluster* c);

        // d_out (closest approach) and pin_collinear_out (scenario 1 AND the combined
        // local/global track-axis collinearity test for xtpc_joint_pin) are optional
        // outputs; the global-axis test is computed only when m_xtpc_joint_pin (off-path
        // unchanged). The scenario return value is unchanged. bx0/bx1 (optional) enable
        // the box-pruned closest-approach path; null falls back to the plain loop.
        int xtpc_pair_consistent(const XtpcMC& m0, const XtpcMC& m1,
                                 double* d_out = nullptr,
                                 bool* pin_collinear_out = nullptr,
                                 const XtpcBoxes* bx0 = nullptr,
                                 const XtpcBoxes* bx1 = nullptr) const;

        // Deterministic iteration orders over the bundle maps (pointer-keyed maps
        // would otherwise iterate in heap-address order). Static: no this-state.
        static std::vector<Opflash*> flash_iter_order(const FlashBundlesMap& m);
        static std::vector<WireCell::Clus::Facade::Cluster*>
            cluster_iter_order(const ClusterBundlesMap& m,
                               const std::map<WireCell::Clus::Facade::Cluster*, int>& idx);
    };

} // namespace WireCell::Match

#endif
