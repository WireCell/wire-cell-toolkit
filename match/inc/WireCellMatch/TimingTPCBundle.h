#ifndef WIRECELL_MATCH_TIMINGTPCBUNDLE
#define WIRECELL_MATCH_TIMINGTPCBUNDLE

#include "WireCellMatch/Opflash.h"
#include "WireCellClus/Facade.h"

#include <map>
#include <memory>
#include <set>
#include <vector>

namespace WireCell::Match {

    using ClusterSelection = std::vector<WireCell::Clus::Facade::Cluster*>;

    /// Bundle-quality thresholds used by examine_bundle()/add_bundle(). Defaults
    /// reproduce the historical hard-coded values; QLMatching pushes its
    /// configured values via set_quality_params() so they are tunable without a
    /// rebuild.
    struct BundleQualityParams {
        double ks_merge_max     = 0.2;   // max KS to accept a merge candidate
        double chi2ndf_merge_max = 20.0;  // max chi2/ndf to accept a merge candidate
        double addmerge_exponent = 0.8;   // exponent in the add_bundle tie-break
        double highconsist_ks_max = 0.06; // max KS for the "high-consistent" flag
        int    highconsist_min_ndf = 3;   // min ndf for the "high-consistent" flag
        double pe_ndf_knee = 1.0;         // per-opdet PE level below which ndf is not counted
        bool   mask_ks = false;           // also apply opdet_mask to the KS shape metric
                                          // (the chi2/LASSO paths always mask). Default OFF
                                          // so existing configs are bit-identical; SBND-on.
        // Per-PMT light-error model for the bundle chi2. When pe_err_on_pred is true the
        // chi2's per-opdet error is PE_err = (pred < knee ? floor : frac*pred) computed from
        // the PREDICTED pe (not the measured-based flash->get_PE_err); sigma^2 = meas + PE_err^2.
        // Default off + floor/frac/knee = 0.3/0.3/1.0 reproduce the measured-based chi2.
        double pe_err_floor = 0.3;
        double pe_err_frac  = 0.3;
        double pe_err_knee  = 1.0;
        bool   pe_err_on_pred = false;
        // Efficiency-aware low-PE error inflation (PD detection inefficiency at low
        // light: a channel with a few PE predicted often measures zero -- far more
        // than Poisson). When pe_err_lowpe_frac >= 0 (and pe_err_on_pred) the relative
        // error grows as the predicted PE falls:
        //   rel(pred) = pe_err_frac + (pe_err_lowpe_frac - pe_err_frac)*exp(-pred/lowpe_knee)
        //   PE_err    = sqrt((rel*pred)^2 + pe_err_floor^2)
        // so rel -> frac at high pred (unchanged) and -> lowpe_frac (~unity) at low pred,
        // letting "predicted a few PE, measured zero" be tolerated rather than penalized.
        // pe_err_lowpe_frac < 0 => disabled (use the floor/frac/knee branch, bit-identical).
        // Calibrated on run-29107 hand scans (pdhd/ql_light_calib/fit_lowpe.py): PDHD 2.1/4.0.
        double pe_err_lowpe_frac = -1.0;
        double pe_err_lowpe_knee = 4.0;
        // Flag-aware multi-branch "high-consistent" ladder (ported from the MicroBooNE
        // prototype FlashTPCBundle, re-tuned for the SBND post-recipe chi2 scale from the
        // 10 hand-scan data events). When highconsist_ladder is false the single-branch
        // (highconsist_ks_max/min_ndf) logic above is used -> bit-identical. SBND-on.
        // Branches (OR), c2n = chi2/ndf: B1 clean very-good; B2 general good; B3 two_boundary;
        // B4 x-boundary/close-PMT/window-truncated (ks-led, chi2 relaxed for missing charge).
        bool   highconsist_ladder = false;
        double hc_clean_ks = 0.06;  double hc_clean_c2 = 6.0;   // B1
        double hc_good_ks  = 0.09;  double hc_good_c2  = 4.0;   // B2
        double hc_tb_ks    = 0.10;  double hc_tb_c2    = 8.0;   // B3
        double hc_miss_ks  = 0.08;  double hc_miss_c2  = 60.0;  // B4
        int    hc_miss_min_ndf = 5;                             // B4 ndf floor (B1-B3 use highconsist_min_ndf)
        // Per-bundle chi2 relaxation ported from the prototype FlashTPCBundle.cxx:480-502
        // (default off -> bit-identical; SBND-on). Two behaviors, on the toolkit's
        // (meas + perr^2) Poisson base: (1) when flag_close_to_PMT and a channel shows a
        // big measured excess (pe-pred > chi2_pmt_excess && pe > chi2_pmt_ratio*pred),
        // widen its denominator by (pe*chi2_pmt_inflate)^2 (softens near-PMT over-
        // response); (2) tolerate one inefficient PMT -- if the single worst-chi2 channel
        // is (pe==0 && pred>0), subtract (max_chi2 - 1) from the total. chi2_pmt_excess is
        // absolute PE -> re-tuned for SBND (prototype 350 is a MicroBooNE light-yield value).
        bool   chi2_relax = false;
        double chi2_pmt_excess  = 350.0;   // measured-excess PE threshold (RE-TUNE for SBND)
        double chi2_pmt_ratio   = 1.3;     // measured/predicted ratio threshold
        double chi2_pmt_inflate = 0.5;     // fraction of pe added in quadrature to the denom
        // Rail-saturated channel widening: denom += (pe*chi2_sat_inflate)^2 when the
        // flash flags this channel saturated (Opflash::get_sat). Gated on the flag ALONE
        // -- the close_to_PMT excess/ratio tests fire on measured >> predicted, the
        // opposite regime from a clipped channel, so they would never fire here.
        // Independent of chi2_relax. Reaches the chi2 only when QLMatching leaves the
        // channel in the opdet mask (saturation_mask_fit false). 0 = off, bit-identical.
        double chi2_sat_inflate = 0.0;
    };

    class TimingTPCBundle {
    public:
        using pointer = std::shared_ptr<TimingTPCBundle>;
        using Cluster = WireCell::Clus::Facade::Cluster;

        TimingTPCBundle(Opflash* flash, Cluster* main_cluster,
                        int flash_index_id, int cluster_index_id);
        ~TimingTPCBundle();
        // Member-wise copy is safe: the Cluster*/Opflash* members are non-owning
        // (the destructor is =default and frees nothing) and the rest are values.
        // QLMatching deep-copies matched bundles before the discarded organize pass
        // so its in-place merge cannot corrupt the live flash_bundles_map.
        TimingTPCBundle(const TimingTPCBundle&) = default;
        TimingTPCBundle& operator=(const TimingTPCBundle&) = default;

        void set_flag_close_to_PMT(bool v) { flag_close_to_PMT = v; }
        void set_flag_at_x_boundary(bool v) { flag_at_x_boundary = v; }
        bool get_flag_close_to_PMT() const { return flag_close_to_PMT; }
        bool get_flag_at_x_boundary() const { return flag_at_x_boundary; }

        std::vector<double>& get_pred_flash() { return pred_flash; }
        void set_pred_flash(const std::vector<double>& v) { pred_flash = v; }
        // Full (per-flash rail-flag-unmasked) prediction for the calib dump /
        // display. Empty unless QLMatching stores it (use_saturation_flag on);
        // falls back to the fit vector so consumers can read unconditionally.
        const std::vector<double>& get_pred_flash_full() const
        {
            return pred_flash_full.empty() ? pred_flash : pred_flash_full;
        }
        void set_pred_flash_full(const std::vector<double>& v) { pred_flash_full = v; }

        Opflash* get_flash() const { return flash; }
        void set_flash(Opflash* f) { flash = f; }
        double get_total_pred_light() const;

        Cluster* get_main_cluster() const { return main_cluster; }
        void set_main_cluster(Cluster* c) { main_cluster = c; }
        Cluster* get_orig_cluster() const { return orig_main_cluster; }
        void set_orig_cluster(Cluster* c) { orig_main_cluster = c; }

        ClusterSelection& get_other_clusters() { return other_clusters; }
        ClusterSelection& get_more_clusters() { return more_clusters; }
        void clear_other_clusters() { other_clusters.clear(); }
        void clear_more_clusters() { more_clusters.clear(); }
        void add_other_cluster(Cluster* c) { other_clusters.push_back(c); }

        std::vector<unsigned int>& get_opdet_mask() { return opdet_mask; }
        void set_opdet_mask(const std::vector<unsigned int>& m) { opdet_mask = m; }

        void set_quality_params(const BundleQualityParams& qp) { m_qp = qp; }
        const BundleQualityParams& get_quality_params() const { return m_qp; }

        bool examine_bundle();
        bool examine_bundle(TimingTPCBundle* candidate_bundle);
        void add_bundle(TimingTPCBundle* candidate_bundle);

        double get_chi2() const { return chi2; }
        void set_chi2(double v) { chi2 = v; }
        int    get_ndf()  const { return ndf; }
        void set_ndf(int v) { ndf = v; }
        double get_ks_dis() const { return ks_dis; }
        void set_ks_dis(double v) { ks_dis = v; }

        void set_consistent_flag(bool v) { flag_high_consistent = v; }
        bool get_consistent_flag() const { return flag_high_consistent; }
        void set_spec_end_flag(bool v) { flag_spec_end = v; }
        bool get_spec_end_flag() const { return flag_spec_end; }
        void set_flag_window_truncated(bool v) { flag_window_truncated = v; }
        bool get_flag_window_truncated() const { return flag_window_truncated; }
        bool get_potential_bad_match_flag() const { return flag_potential_bad_match; }
        void set_potential_bad_match_flag(bool v) { flag_potential_bad_match = v; }
        // TPC-containment verdict from compute_endpoint_flags (calib/diagnostic
        // only; nothing in the matching path reads it). Default true.
        void set_contained(bool v) { flag_contained = v; }
        bool get_contained() const { return flag_contained; }
        // Both main-PCA extremes of the main cluster lie at a detector edge (the
        // cluster enters AND exits the per-APA active volume). Diagnostic only
        // (calib dump); nothing in the matching path reads it. Default false.
        void set_flag_two_boundary(bool v) { flag_two_boundary = v; }
        bool get_flag_two_boundary() const { return flag_two_boundary; }
        // Cathode-end proximity (set alongside the cathode-end at_x_boundary in
        // compute_endpoint_flags). Inert diagnostic for now — nothing in the
        // matching path reads it; exported in the calib dump so a PDVD
        // cathode-PD treatment can be designed from hand-scan data. Default false.
        void set_flag_at_cathode(bool v) { flag_at_cathode = v; }
        bool get_flag_at_cathode() const { return flag_at_cathode; }
        // Per-channel chi2-relax eligibility (QLMatching vd_surface_flags mode).
        // EMPTY (default) => every channel is eligible, which is the historical
        // behavior — the close-to-PMT chi2 denominator inflation may fire on any
        // channel with a big measured excess. When filled (PDVD: only the PD
        // channels of the surface the activity is actually near — that wall's
        // membrane XAs, or the bottom PMTs), the inflation is restricted to those
        // channels. add_relax_channels unions surfaces across calls.
        void add_relax_channels(const std::vector<int>& chs);
        const std::vector<uint8_t>& get_relax_channels() const { return relax_channels; }
        // Cross-TPC cathode-crossing confirmation (post-matching): the matched main
        // cluster connects geometrically, across the cathode, to the coincident other
        // TPC's matched main cluster. Confirm-stamp only; nothing in the matching path
        // reads it. Default false. See QLMatching::flag_cross_tpc_consistency.
        void set_flag_xtpc_consistent(bool v) { flag_xtpc_consistent = v; }
        bool get_flag_xtpc_consistent() const { return flag_xtpc_consistent; }
        // Subset of flag_xtpc_consistent: the cross-TPC pair passed scenario 1 (both
        // cathode ends present, closest approach < xtpc_dmax — a tight, self-vetoing
        // match). Drives the xtpc-PRIORITY cull (a scenario-1 crosser half overrides a
        // cluster's high-consistent bundle on a different flash). Default false.
        void set_flag_xtpc_scenario1(bool v) { flag_xtpc_scenario1 = v; }
        bool get_flag_xtpc_scenario1() const { return flag_xtpc_scenario1; }
        // Joint-pin winner: this bundle is the chosen (cluster,flash) of a confirmed
        // direction-collinear cross-TPC pair (scenario 1 AND track axes collinear by the
        // combined local-vhough / global-PCA test). Set only when xtpc_joint_pin is on;
        // drives the TOP-priority cull (a cluster owning a pinned bundle keeps ONLY it),
        // binding both crosser halves to one coincident flash. Default false => off-path
        // never sets it => bit-identical.
        void set_flag_xtpc_pin(bool v) { flag_xtpc_pin = v; }
        bool get_flag_xtpc_pin() const { return flag_xtpc_pin; }

        double get_strength() const { return strength; }
        void set_strength(double v) { strength = v; }

        int get_flash_index_id()   const { return flash_index_id; }
        int get_cluster_index_id() const { return cluster_index_id; }

    private:
        Opflash* flash;
        Cluster* main_cluster;
        Cluster* orig_main_cluster;
        int flash_index_id;
        int cluster_index_id;
        std::vector<unsigned int> opdet_mask;

        bool flag_close_to_PMT;
        bool flag_at_x_boundary;
        bool flag_spec_end;
        bool flag_window_truncated;
        bool flag_potential_bad_match;
        bool flag_high_consistent;
        bool flag_contained;
        bool flag_two_boundary;
        bool flag_at_cathode{false};
        bool flag_xtpc_consistent{false};
        bool flag_xtpc_scenario1{false};
        bool flag_xtpc_pin{false};

        double ks_dis;
        double chi2;
        int    ndf;
        double strength;

        int m_nchan;
        BundleQualityParams m_qp;
        // Per-channel chi2-relax eligibility (see add_relax_channels). Empty =>
        // all channels eligible (historical behavior, bit-identical).
        std::vector<uint8_t> relax_channels;
        std::vector<double> pred_flash;
        std::vector<double> pred_flash_full;
        std::vector<Cluster*> other_clusters;
        std::vector<Cluster*> more_clusters;
    };

    using TimingTPCBundleSelection = std::vector<TimingTPCBundle::pointer>;
    using TimingTPCBundleSet       = std::set<TimingTPCBundle::pointer>;
    using FlashBundlesMap          = std::map<Opflash*, TimingTPCBundleSelection>;
    using ClusterBundlesMap        =
        std::map<WireCell::Clus::Facade::Cluster*, TimingTPCBundleSelection>;

} // namespace WireCell::Match

#endif
