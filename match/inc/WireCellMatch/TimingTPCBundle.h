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
        bool flag_xtpc_consistent{false};
        bool flag_xtpc_scenario1{false};

        double ks_dis;
        double chi2;
        int    ndf;
        double strength;

        int m_nchan;
        BundleQualityParams m_qp;
        std::vector<double> pred_flash;
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
