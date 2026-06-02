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
    };

    class TimingTPCBundle {
    public:
        using pointer = std::shared_ptr<TimingTPCBundle>;
        using Cluster = WireCell::Clus::Facade::Cluster;

        TimingTPCBundle(Opflash* flash, Cluster* main_cluster,
                        int flash_index_id, int cluster_index_id);
        ~TimingTPCBundle();

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
