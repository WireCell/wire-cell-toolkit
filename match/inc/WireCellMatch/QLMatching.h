#ifndef WIRECELL_MATCH_QLMATCHING
#define WIRECELL_MATCH_QLMATCHING

#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IFiducial.h"
#include "WireCellIface/ITensorSetFanin.h"

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
        bool m_beamonly{false};
        double m_flash_minPE{500};
        double m_flash_mintime{-1.5 * units::ms};
        double m_flash_maxtime{1.5 * units::ms};
        double m_beam_mintime{-5 * units::us};
        double m_beam_maxtime{5 * units::us};
        double m_QtoL{0.5};
        // Drift speed used for the per-flash X correction (flash_x_offset =
        // sign * flash_time * drift_speed). In WCT internal units (length/time,
        // i.e. mm/ns numerically). Default is the historical hard-coded value;
        // configs should pass the common params.lar.drift_speed via jsonnet.
        double m_drift_speed{1.563 * units::mm / units::us};
        // LASSO solution threshold below which a (flash, cluster) bundle is
        // dropped after each matching round. Was hard-coded 0.05 inline; pulled
        // out so it can be widened/narrowed from the jsonnet without rebuild.
        double m_strength_cutoff{0.05};

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
        double m_y_cushion{0.0 * units::cm};       // signed inward(+)/outward(-) shift of each |y| edge
        double m_z_cushion{0.0 * units::cm};       // signed inward(+)/outward(-) shift of each z edge
        // Proximity band for the diagnostic flag_two_boundary edge test: a main-PCA
        // extreme counts as "at the detector edge" if it is within this distance of
        // any of the 6 per-APA active-volume faces (anode, cathode, ±y, ±z).
        double m_two_boundary_margin{3.0 * units::cm};

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
        double m_flash_pe_threshold{0.0};      // Opflash "fired" threshold (PE)

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

        // §H raw readout-window truncation flag (T0-independent, APA-agnostic),
        // always computed. Flags a bundle whose cluster's leading/trailing time
        // slice sits within m_window_edge_ticks of the raw readout window
        // [0, m_readout_window_ticks]. Inert (no consumer yet); since nothing
        // reads it, always filling it leaves production output unchanged.
        int m_readout_window_ticks{3427};  // window end in raw ticks (SBND daq.nticks)
        int m_window_edge_ticks{4};        // edge-proximity threshold (~one slice = nticks_live_slice)

        // When true, discard any (flash, cluster) bundle whose cluster is not
        // contained in the TPC drift box once the flash-derived T0 x-offset is
        // applied (the prototype flag_good_bundle gate, ToyMatching.cxx 272-275).
        // Default OFF so existing production configs stay bit-identical; SBND
        // jsonnet sets require_containment: true.
        bool m_require_containment{false};

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

        // Path to the JSON file holding VUVHits, VISHits, geometry and the
        // SBND OpDet array.
        std::string m_semimodel_file{"sbnd/photodet/semi-analytical-sbnd.json"};

        std::unique_ptr<SemiAnalyticalModel> m_semi_model;

        // Per-OpDet table (type + position) loaded from semimodel_file. Used to
        // build the SemiAnalyticalModel and, at execute() time, to derive the
        // per-channel on/off mask and the per-TPC OpDet split.
        std::vector<SemiAnalyticalModel::OpticalDetector> m_opdets;
        double m_cathode_x{0.0};  // cathode-plane x; OpDets on its low/high side belong to TPC 0/1

        // Fill the per-bundle boundary flags (flag_close_to_PMT,
        // flag_at_x_boundary, flag_spec_end) for one (flash, cluster) bundle.
        // Ports the MicroBooNE prototype end-trimming walk (ToyMatching.cxx
        // ~176-290) into the per-TPC drift coordinate u. s/anode_x/u_cathode are
        // the geometry scalars computed once per anode in operator().
        //
        // Returns the prototype's flag_good_bundle / TPC-containment verdict: true
        // iff the (post-trim) cluster endpoints fall inside the TPC drift box (the
        // 4-part in-window guard, == ToyMatching.cxx 272-275). Returns false when
        // the cluster has no 3d points under this anode (it cannot be contained).
        // The caller discards uncontained bundles when m_require_containment.
        bool compute_endpoint_flags(TimingTPCBundle* bundle,
                                    WireCell::Clus::Facade::Cluster* cluster,
                                    double flash_x_offset,
                                    double s, double anode_x, double u_cathode,
                                    int anode_ident) const;

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
        void compute_geometry(ApaRun& run);          // per-TPC drift geometry, mask cull, opdet idx
        void compute_dynamic_opdet_mask(ApaRun& run, unsigned int tpc);  // per-event dead-PMT auto-mask
        void build_bundles(ApaRun& run);             // (flash,group) bundles + predicted light  [Stage 1]

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
        void cull_inconsistent(ApaRun& run);         // drop non-consistent rivals               [Stage 1]
        void fit_round1(ApaRun& run);                // LASSO, per-flash background DOF           [Stage 2]
        void fit_round2(ApaRun& run);                // LASSO + KS-shape, keep best per cluster   [Stage 3]

        // Per-column L1 down-weight (m_lasso_boundary_weight) for a boundary / near-PMT /
        // window-truncated bundle when m_lasso_flag_weight is on (prototype's flag-based
        // weight), else 1.0. Multiplies the pe-mismatch (+KS) base in both rounds.
        double lasso_flag_factor(const TimingTPCBundle::pointer& bundle) const;

        // Empty-flash light-quality rescue (m_empty_rescue; see §I). snapshot is the
        // full pre-strength-cutoff flash->candidate map; this mutates run.flash_bundles_map
        // in place, adopting the best light-quality candidate of each emptied flash.
        void rescue_empty_flashes(ApaRun& run, const FlashBundlesMap& snapshot);
        void apply_matched_t0s(ApaRun& run);         // write cluster t0 / flash / matched gid
        void write_opflash_pc(ApaRun& run);          // merge-safe per-root "opflash" PC

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
        int xtpc_pair_consistent(const XtpcMC& m0, const XtpcMC& m1) const;

        // Deterministic iteration orders over the bundle maps (pointer-keyed maps
        // would otherwise iterate in heap-address order). Static: no this-state.
        static std::vector<Opflash*> flash_iter_order(const FlashBundlesMap& m);
        static std::vector<WireCell::Clus::Facade::Cluster*>
            cluster_iter_order(const ClusterBundlesMap& m,
                               const std::map<WireCell::Clus::Facade::Cluster*, int>& idx);
    };

} // namespace WireCell::Match

#endif
