/** CheckSTM_Michel -- stopping-muon + Michel-electron reconstruction (doc pdvd/48).
 *
 * A light replacement for the neutrino PR tail (TaggerCheckNeutrino) on a
 * cosmic detector.  Runs AFTER tagger_check_stm and reconstructs every
 * STM-tagged main cluster as
 *
 *     entry point  --(muon body, delta rays)-->  Bragg stop  --> Michel e-
 *                                                              (+ nearby dots)
 *
 * with the particle flow ROOTED AT THE ENTRY POINT (the STM tagger's L = 0
 * end), and a REJECT verdict for tagged clusters that do not look like a
 * stopped muon (still reconstructed, flagged).
 *
 * What is reused from the general PR chain and how (all callable without
 * TaggerCheckNeutrino):
 *   - PR::PatternAlgorithms as a plain local object (TaggerCheckNeutrino.cxx:2326):
 *     find_proto_vertex -> clustering_points -> separate_track_shower ->
 *     determine_direction on the main cluster (:2876-2893) and on the
 *     nearby companion clusters (:2927-2962), then examine_direction from
 *     the caller-chosen entry vertex (:3147, here with flag_final=true).
 *   - TrackFitting owned per candidate exactly as :2253-2267 (fresh fitter,
 *     parameters copied from the member fitter that load_trackfitting_config
 *     filled from the runtime JSON, add_graph, preload_clusters).
 *   - PR::Shower for the Michel (PRShower.h), IndexedShowerSet -> set_showers.
 *   - The STM tagger's own products as the anchor: Flags::STM, the stm_pass /
 *     stm_fit cluster PCs (TaggerCheckSTM.cxx:936-1022): L = 0 is the entry,
 *     the kink_num-th fit row of the accepted pass is the stop.
 *   - Publication identical to TaggerCheckNeutrino.cxx:3588-3590 so every
 *     existing consumer (Bee track_fit/shower_track/vertices/mc layers,
 *     PdvdPrMagnifyTrackingVisitor, PrDisplayDump) renders this stage's
 *     graph unchanged: the unnamed slot = candidate 0, and "nu<i>" per
 *     candidate (the slot-name convention those writers walk; inherited, not
 *     chosen here).
 *
 * The graph-only logic (chain walk, profile, Bragg contrast, arm
 * classification) lives in StmMichelFunctions.{h,cxx} so it is doctested
 * on synthetic graphs.
 *
 * Nothing here runs unless the component is named in a pipeline, so no
 * existing detector output can change.  Shared C++ touched by this round:
 * none in clus/ (two new files), root/PdvdPrMagnifyTrackingVisitor.cxx gains
 * two trees written only when the stm_michel PC exists.
 */

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/ParticleDataSet.h"
#include "WireCellClus/FiducialUtils.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRShower.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRTrajectoryView.h"
#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/TrackFittingPresets.h"
#include "WireCellClus/StmMichelFunctions.h"
#include "WireCellAux/ParticleInfo.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IScalarFunction.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/KSTest.h"
#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/Units.h"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>

class CheckSTM_Michel;
WIRECELL_FACTORY(CheckSTM_Michel, CheckSTM_Michel,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::Clus::PR;

static auto s_log = WireCell::Log::logger("clus.CheckSTM_Michel");

class CheckSTM_Michel : public IConfigurable, public Clus::IEnsembleVisitor,
                        private Clus::NeedDV, private Clus::NeedPCTS,
                        private Clus::NeedRecombModel, private Clus::NeedParticleData,
                        private Clus::NeedFiducial {
public:
    CheckSTM_Michel() {
        m_track_fitter = std::make_shared<TrackFitting>(TrackFittingPresets::create_with_current_values());
    }
    virtual ~CheckSTM_Michel() {}

    virtual void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config);
        NeedRecombModel::configure(config);
        NeedParticleData::configure(config);
        // Same guard as TaggerCheckSTM.cxx:102: the mixin's type-only default
        // "DetectorVolumes" is not an instantiated component in the PR config.
        m_use_fiducial = !config["fiducial"].isNull();
        if (m_use_fiducial) NeedFiducial::configure(config);

        m_cfg = config;
        m_grouping_name = get<std::string>(config, "grouping", m_grouping_name);
        m_perf = get<bool>(config, "perf", m_perf);
        m_require_stm_flag = get<bool>(config, "require_stm_flag", m_require_stm_flag);
        m_publish_nu_slots = get<bool>(config, "publish_nu_slots", m_publish_nu_slots);
        m_build_michel_shower = get<bool>(config, "build_michel_shower", m_build_michel_shower);
        m_max_candidates = get<int>(config, "max_candidates", m_max_candidates);
        m_min_chain_points = get<int>(config, "min_chain_points", m_min_chain_points);

        // e/cm, the TaggerCheckNeutrino convention (converted at the PA copy).
        m_mip_dqdx = get<double>(config, "mip_dqdx", m_mip_dqdx);
        m_mip_dqdx_median = get<double>(config, "mip_dqdx_median", m_mip_dqdx_median);

        // verdict thresholds (cm / deg / ratios; see default_configuration)
        m_stop_snap_tol_cm = get<double>(config, "stop_snap_tol_cm", m_stop_snap_tol_cm);
        m_entry_snap_tol_cm = get<double>(config, "entry_snap_tol_cm", m_entry_snap_tol_cm);
        m_stop_fv_margin_cm = get<double>(config, "stop_fv_margin_cm", m_stop_fv_margin_cm);
        m_bragg_contrast_min = get<double>(config, "bragg_contrast_min", m_bragg_contrast_min);
        m_bragg_tail_lo_cm = get<double>(config, "bragg_tail_lo_cm", m_bragg_tail_lo_cm);
        m_bragg_tail_hi_cm = get<double>(config, "bragg_tail_hi_cm", m_bragg_tail_hi_cm);
        m_bragg_plateau_lo_cm = get<double>(config, "bragg_plateau_lo_cm", m_bragg_plateau_lo_cm);
        m_bragg_plateau_hi_cm = get<double>(config, "bragg_plateau_hi_cm", m_bragg_plateau_hi_cm);
        m_ks_margin = get<double>(config, "ks_margin", m_ks_margin);
        m_compare_range_cm = get<double>(config, "compare_range_cm", m_compare_range_cm);
        m_offset_length_cm = get<double>(config, "offset_length_cm", m_offset_length_cm);
        m_continuation_max_angle_deg = get<double>(config, "continuation_max_angle_deg", m_continuation_max_angle_deg);
        m_continuation_min_len_cm = get<double>(config, "continuation_min_len_cm", m_continuation_min_len_cm);
        m_continuation_mip_lo = get<double>(config, "continuation_mip_lo", m_continuation_mip_lo);
        m_continuation_mip_hi = get<double>(config, "continuation_mip_hi", m_continuation_mip_hi);
        m_michel_max_len_cm = get<double>(config, "michel_max_len_cm", m_michel_max_len_cm);
        m_michel_mip_hi = get<double>(config, "michel_mip_hi", m_michel_mip_hi);
        m_michel_mip_lo = get<double>(config, "michel_mip_lo", m_michel_mip_lo);
        m_michel_min_kink_deg = get<double>(config, "michel_min_kink_deg", m_michel_min_kink_deg);
        m_michel_dot_radius_cm = get<double>(config, "michel_dot_radius_cm", m_michel_dot_radius_cm);
        m_dot_max_len_cm = get<double>(config, "dot_max_len_cm", m_dot_max_len_cm);
        m_dot_body_exclusion_cm = get<double>(config, "dot_body_exclusion_cm", m_dot_body_exclusion_cm);
        m_delta_max_len_cm = get<double>(config, "delta_max_len_cm", m_delta_max_len_cm);
        m_vertex_hadron_mip = get<double>(config, "vertex_hadron_mip", m_vertex_hadron_mip);

        // TrackFitting parameters carried by TaggerCheckNeutrino's config
        // rather than by the runtime JSON (TaggerCheckNeutrino.cxx:2777-2791).
        m_fit_blob_coverage = get<double>(config, "fit_blob_coverage", m_fit_blob_coverage);
        m_dqdx_fit_keep_all_points = get<bool>(config, "dqdx_fit_keep_all_points", m_dqdx_fit_keep_all_points);
        m_excl_t0_frame = get<bool>(config, "excl_t0_frame", m_excl_t0_frame);

        auto tf_config = get<std::string>(config, "trackfitting_config_file", "");
        if (!tf_config.empty()) load_trackfitting_config(tf_config);
    }

    virtual WireCell::Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping_name;
        cfg["detector_volumes"] = "DetectorVolumes";
        cfg["pc_transforms"] = "PCTransformSet";
        cfg["recombination_model"] = "BoxRecombination";
        cfg["particle_dataset"] = "ParticleDataSet";
        cfg["fiducial"] = Json::Value();
        cfg["trackfitting_config_file"] = "";
        cfg["perf"] = m_perf;
        cfg["require_stm_flag"] = m_require_stm_flag;
        cfg["publish_nu_slots"] = m_publish_nu_slots;
        cfg["build_michel_shower"] = m_build_michel_shower;
        cfg["max_candidates"] = m_max_candidates;
        cfg["min_chain_points"] = m_min_chain_points;
        cfg["mip_dqdx"] = m_mip_dqdx;                 // e/cm
        cfg["mip_dqdx_median"] = m_mip_dqdx_median;   // e/cm
        cfg["stop_snap_tol_cm"] = m_stop_snap_tol_cm;
        cfg["entry_snap_tol_cm"] = m_entry_snap_tol_cm;   // looser: an existing end vertex beats a split that leaves a stub
        cfg["stop_fv_margin_cm"] = m_stop_fv_margin_cm;
        cfg["bragg_contrast_min"] = m_bragg_contrast_min;   // fraction of the tabulated rise
        cfg["bragg_tail_lo_cm"] = m_bragg_tail_lo_cm;
        cfg["bragg_tail_hi_cm"] = m_bragg_tail_hi_cm;
        cfg["bragg_plateau_lo_cm"] = m_bragg_plateau_lo_cm;
        cfg["bragg_plateau_hi_cm"] = m_bragg_plateau_hi_cm;
        cfg["ks_margin"] = m_ks_margin;
        cfg["compare_range_cm"] = m_compare_range_cm;
        cfg["offset_length_cm"] = m_offset_length_cm;
        cfg["continuation_max_angle_deg"] = m_continuation_max_angle_deg;
        cfg["continuation_min_len_cm"] = m_continuation_min_len_cm;
        cfg["continuation_mip_lo"] = m_continuation_mip_lo;
        cfg["continuation_mip_hi"] = m_continuation_mip_hi;
        cfg["michel_max_len_cm"] = m_michel_max_len_cm;
        cfg["michel_mip_hi"] = m_michel_mip_hi;
        cfg["michel_mip_lo"] = m_michel_mip_lo;
        cfg["michel_min_kink_deg"] = m_michel_min_kink_deg;
        cfg["michel_dot_radius_cm"] = m_michel_dot_radius_cm;
        cfg["dot_max_len_cm"] = m_dot_max_len_cm;
        cfg["dot_body_exclusion_cm"] = m_dot_body_exclusion_cm;
        cfg["delta_max_len_cm"] = m_delta_max_len_cm;
        cfg["vertex_hadron_mip"] = m_vertex_hadron_mip;
        cfg["fit_blob_coverage"] = m_fit_blob_coverage;
        cfg["dqdx_fit_keep_all_points"] = m_dqdx_fit_keep_all_points;
        cfg["excl_t0_frame"] = m_excl_t0_frame;
        // The PR-partition knobs (the four early stages) are read straight
        // from the config with PatternAlgorithms' own C++ defaults -- see
        // apply_pattern_knobs().  Published here so a compiled config can be
        // checked against the list; absent keys ride the C++ defaults.
        for (const auto& k : pattern_knob_keys()) cfg[k] = Json::Value();
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const;

private:
    using Clock = std::chrono::steady_clock;
    using MS = std::chrono::duration<double, std::milli>;

    // ---- config ------------------------------------------------------------
    Configuration m_cfg;
    std::string m_grouping_name{"live"};
    bool m_perf{false};
    bool m_use_fiducial{false};
    bool m_require_stm_flag{true};
    bool m_publish_nu_slots{true};
    bool m_build_michel_shower{true};
    int m_max_candidates{8};
    int m_min_chain_points{10};
    double m_mip_dqdx{50000.0};          // e/cm
    double m_mip_dqdx_median{43000.0};   // e/cm
    double m_stop_snap_tol_cm{2.0};
    double m_entry_snap_tol_cm{5.0};
    double m_stop_fv_margin_cm{5.0};
    double m_bragg_contrast_min{0.6};
    double m_bragg_tail_lo_cm{0.5}, m_bragg_tail_hi_cm{3.0};
    double m_bragg_plateau_lo_cm{20.0}, m_bragg_plateau_hi_cm{40.0};
    double m_ks_margin{0.0};
    double m_compare_range_cm{35.0}, m_offset_length_cm{0.0};
    double m_continuation_max_angle_deg{20.0}, m_continuation_min_len_cm{3.0};
    double m_continuation_mip_lo{0.7}, m_continuation_mip_hi{1.3};
    double m_michel_max_len_cm{25.0}, m_michel_mip_hi{2.0}, m_michel_mip_lo{0.3}, m_michel_min_kink_deg{30.0};
    double m_michel_dot_radius_cm{15.0}, m_dot_max_len_cm{10.0}, m_dot_body_exclusion_cm{5.0};
    double m_delta_max_len_cm{8.0};
    double m_vertex_hadron_mip{1.4};
    double m_fit_blob_coverage{-1.0};
    bool m_dqdx_fit_keep_all_points{false};
    bool m_excl_t0_frame{false};

    // The member fitter is only the parameter holder the runtime JSON lands
    // in; every candidate gets its own fitter seeded from it (see visit()).
    std::shared_ptr<TrackFitting> m_track_fitter;
    mutable std::string m_evt_tag;

    // ---- per-candidate record ---------------------------------------------
    struct Record {
        int cluster_id{-1}, gid{-1};
        double t0_us{0};
        unsigned reject_bits{0};
        int has_pass{0}, pass{-1}, kink_num{-1};
        Point entry_pt, stop_pt, tagger_stop_pt;
        int entry_vtx_id{-1}, stop_vtx_id{-1};
        double stop_dis{0};             // graph stop vertex vs tagger stop point
        int n_chain_segs{0}, n_profile_pts{0};
        double muon_len{0};
        int n_delta{0}, n_body_other{0}, n_body_hadron{0};
        double delta_len{0};
        // dQ/dx shape
        double ks_mu{0}, ks_flat{0}, ratio_mu{0}, ratio_flat{0};
        std::vector<double> comp_fwd{0, 1e9, 1e9, 1e9}, comp_bwd{0, 1e9, 1e9, 1e9};
        StmMichelBragg bragg;
        // stop arms
        int n_stop_arms{0}, n_michel_segs{0}, michel_found{0}, michel_conn_type{0};
        double michel_len{0}, michel_mip{0}, michel_kink_deg{-1}, michel_far_len{0};
        double michel_ke_dqdx{0}, michel_ke_range{0}, michel_ke_best{0};
        double cont_len{0}, cont_angle_deg{-1}, cont_mip{0};
        // dots
        int n_dots{0}, n_dot_clusters_unfit{0};
        double dots_ke_dqdx{0}, dots_charge_unfit{0};
        int in_fv{-1};
        // points for the Bee/ROOT layer
        std::vector<double> px, py, pz, pq, pL, prr;
        std::vector<int> prole, pseg;
    };

    // ---- helpers -------------------------------------------------------------
    static const std::vector<std::string>& pattern_knob_keys() {
        // The subset of TaggerCheckNeutrino's config keys that the four early
        // stages (find_proto_vertex, clustering_points, separate_track_shower,
        // determine_direction) and examine_direction read.  Same names, same
        // units (TaggerCheckNeutrino.cxx:2327-2530), so a PDVD job can forward
        // its production partition knobs unchanged.
        static const std::vector<std::string> keys = {
            "dir_weak_use_score", "proton_dir_vote", "proton_dir_score_max", "proton_dir_asym_min",
            "endpoint_trim_retry", "cathode_x", "cathode_kink_xcut",
            "cathode_wide_kink_angle", "cathode_wide_kink_skirt", "cathode_wide_kink_baseline",
            "two_end_break", "teb_min_len", "teb_min_arm", "teb_min_arm_pts", "teb_stub_max",
            "teb_accept_range", "teb_rise_r1", "teb_rise_r2", "teb_abs_end_min", "teb_dip_floor",
            "teb_score_cap_r1", "teb_score_cap_r2", "teb_turn_angle", "teb_turn_baseline",
            "teb_turn_skirt", "teb_turn_min_arm_frac", "teb_bragg_veto_turn",
            "kink_walk_dqdx_stop", "kink_break_protect", "kink_dqdx_hot_ratio",
            "shower_topo_demote_len", "shower_topo_reset", "shower_topo_dqdx_guard", "shower_topo_proto_dir",
            "shower_traj_straight_guard",
            "fit_exclusion", "oov_prototype_parity", "first_seg_local_pca",
            "other_seg_relaxed_accept", "other_seg_empty_2d_guard", "other_seg_keep_isolated",
            "other_seg_keep_isolated_min_points", "other_seg_keep_isolated_min_length",
            "other_seg_keep_isolated_len_admit", "iso_snap_min_dir_mag",
            "assoc_full_recluster", "assoc_reassign_orphans", "assoc_clear_on_merge",
            "track_comp_empty_abstain", "reclass_preserve_4mom", "reclass_never_computed_ke_floor",
            "dir_track_median_local",
            "iso_endpoint", "iso_endpoint_min_length", "iso_endpoint_max_xext", "iso_endpoint_xext_frac",
            "iso_endpoint_xext_quantile", "iso_endpoint_tube_radius", "iso_endpoint_min_aspect",
            "pr_find_other_rounds", "v3_extension_guard", "v3_extension_min_gain",
            "es3_stub_guard", "es3sg_stub_max", "es3sg_len_ratio", "es3sg_ang3_min", "es3sg_ang_ratio",
            "es3sg_require_terminal",
            "steiner_gap_penalty", "sgp_dead_alpha", "sgp_min_edge", "sgp_sample_step", "sgp_point_radius",
            "sgp_weak_scale", "sgp_weak_qref", "sgp_max_sep", "good_point_pitch_frac",
            "break_seg_orient", "graph_endpoint_tol",
        };
        return keys;
    }

    // Copy the PR-partition knobs onto a PatternAlgorithms, reading the
    // config directly with the PA's own C++ defaults as fallback.  Unit
    // conventions are TaggerCheckNeutrino.cxx:2327-2530's: cm keys scaled by
    // units::cm, dQ/dx keys divided by units::cm, degrees/ratios/bools raw.
    void apply_pattern_knobs(PatternAlgorithms& pa) const {
        auto B  = [&](bool& f, const char* k)   { f = get<bool>(m_cfg, k, f); };
        auto I  = [&](int& f, const char* k)    { f = get<int>(m_cfg, k, f); };
        auto D  = [&](double& f, const char* k) { f = get<double>(m_cfg, k, f); };
        auto CM = [&](double& f, const char* k) { if (m_cfg.isMember(k)) f = m_cfg[k].asDouble() * units::cm; };

        pa.m_perf = m_perf;
        pa.m_mip_dqdx        = m_mip_dqdx / units::cm;
        pa.m_mip_dqdx_median = m_mip_dqdx_median / units::cm;
        pa.m_sgp_dv   = m_dv;
        pa.m_sgp_pcts = m_pcts;
        pa.m_recomb_model = m_recomb_model;

        B(pa.m_dir_weak_use_score, "dir_weak_use_score");
        B(pa.m_proton_dir_vote, "proton_dir_vote");
        D(pa.m_proton_dir_score_max, "proton_dir_score_max");
        D(pa.m_proton_dir_asym_min, "proton_dir_asym_min");
        B(pa.m_endpoint_trim_retry, "endpoint_trim_retry");
        CM(pa.m_cathode_x, "cathode_x");
        CM(pa.m_cathode_kink_xcut, "cathode_kink_xcut");
        D(pa.m_cathode_wide_kink_angle, "cathode_wide_kink_angle");
        CM(pa.m_cathode_wide_kink_skirt, "cathode_wide_kink_skirt");
        CM(pa.m_cathode_wide_kink_baseline, "cathode_wide_kink_baseline");
        B(pa.m_two_end_break, "two_end_break");
        CM(pa.m_teb_min_len, "teb_min_len");
        CM(pa.m_teb_min_arm, "teb_min_arm");
        I(pa.m_teb_min_arm_pts, "teb_min_arm_pts");
        CM(pa.m_teb_stub_max, "teb_stub_max");
        CM(pa.m_teb_accept_range, "teb_accept_range");
        D(pa.m_teb_rise_r1, "teb_rise_r1");
        D(pa.m_teb_rise_r2, "teb_rise_r2");
        D(pa.m_teb_abs_end_min, "teb_abs_end_min");
        D(pa.m_teb_dip_floor, "teb_dip_floor");
        D(pa.m_teb_score_cap_r1, "teb_score_cap_r1");
        D(pa.m_teb_score_cap_r2, "teb_score_cap_r2");
        D(pa.m_teb_turn_angle, "teb_turn_angle");
        CM(pa.m_teb_turn_baseline, "teb_turn_baseline");
        CM(pa.m_teb_turn_skirt, "teb_turn_skirt");
        D(pa.m_teb_turn_min_arm_frac, "teb_turn_min_arm_frac");
        D(pa.m_teb_bragg_veto_turn, "teb_bragg_veto_turn");
        B(pa.m_kink_walk_dqdx_stop, "kink_walk_dqdx_stop");
        B(pa.m_kink_break_protect, "kink_break_protect");
        D(pa.m_kink_dqdx_hot_ratio, "kink_dqdx_hot_ratio");
        CM(pa.m_shower_topo_demote_len, "shower_topo_demote_len");
        B(pa.m_shower_topo_reset, "shower_topo_reset");
        B(pa.m_shower_topo_dqdx_guard, "shower_topo_dqdx_guard");
        B(pa.m_shower_topo_proto_dir, "shower_topo_proto_dir");
        B(pa.m_shower_traj_straight_guard, "shower_traj_straight_guard");
        B(pa.m_fit_exclusion, "fit_exclusion");
        B(pa.m_oov_prototype_parity, "oov_prototype_parity");
        B(pa.m_first_seg_local_pca, "first_seg_local_pca");
        B(pa.m_other_seg_relaxed_accept, "other_seg_relaxed_accept");
        B(pa.m_other_seg_empty_2d_guard, "other_seg_empty_2d_guard");
        B(pa.m_other_seg_keep_isolated, "other_seg_keep_isolated");
        I(pa.m_other_seg_keep_isolated_min_points, "other_seg_keep_isolated_min_points");
        CM(pa.m_other_seg_keep_isolated_min_length, "other_seg_keep_isolated_min_length");
        CM(pa.m_other_seg_keep_isolated_len_admit, "other_seg_keep_isolated_len_admit");
        CM(pa.m_iso_snap_min_dir_mag, "iso_snap_min_dir_mag");
        B(pa.m_assoc_full_recluster, "assoc_full_recluster");
        B(pa.m_assoc_reassign_orphans, "assoc_reassign_orphans");
        B(pa.m_assoc_clear_on_merge, "assoc_clear_on_merge");
        B(pa.m_track_comp_empty_abstain, "track_comp_empty_abstain");
        B(pa.m_reclass_preserve_4mom, "reclass_preserve_4mom");
        B(pa.m_reclass_never_computed_ke_floor, "reclass_never_computed_ke_floor");
        B(pa.m_dir_track_median_local, "dir_track_median_local");
        B(pa.m_iso_endpoint, "iso_endpoint");
        CM(pa.m_iso_endpoint_min_length, "iso_endpoint_min_length");
        CM(pa.m_iso_endpoint_max_xext, "iso_endpoint_max_xext");
        D(pa.m_iso_endpoint_xext_frac, "iso_endpoint_xext_frac");
        D(pa.m_iso_endpoint_xext_quantile, "iso_endpoint_xext_quantile");
        CM(pa.m_iso_endpoint_tube_radius, "iso_endpoint_tube_radius");
        D(pa.m_iso_endpoint_min_aspect, "iso_endpoint_min_aspect");
        I(pa.m_pr_find_other_rounds, "pr_find_other_rounds");
        B(pa.m_v3_extension_guard, "v3_extension_guard");
        CM(pa.m_v3_extension_min_gain, "v3_extension_min_gain");
        B(pa.m_es3_stub_guard, "es3_stub_guard");
        CM(pa.m_es3sg_stub_max, "es3sg_stub_max");
        D(pa.m_es3sg_len_ratio, "es3sg_len_ratio");
        D(pa.m_es3sg_ang3_min, "es3sg_ang3_min");
        D(pa.m_es3sg_ang_ratio, "es3sg_ang_ratio");
        B(pa.m_es3sg_require_terminal, "es3sg_require_terminal");
        D(pa.m_steiner_gap_penalty, "steiner_gap_penalty");
        D(pa.m_sgp_dead_alpha, "sgp_dead_alpha");
        CM(pa.m_sgp_min_edge, "sgp_min_edge");
        CM(pa.m_sgp_sample_step, "sgp_sample_step");
        CM(pa.m_sgp_point_radius, "sgp_point_radius");
        D(pa.m_sgp_weak_scale, "sgp_weak_scale");
        D(pa.m_sgp_weak_qref, "sgp_weak_qref");
        CM(pa.m_sgp_max_sep, "sgp_max_sep");
        D(pa.m_good_point_pitch_frac, "good_point_pitch_frac");
        B(pa.m_break_seg_orient, "break_seg_orient");
        // Process-wide state (TaggerCheckNeutrino.cxx:2505): only written
        // when the job configures it, so a job that does not name the key
        // leaves whatever the process already had.
        if (m_cfg.isMember("graph_endpoint_tol")) {
            pa.m_graph_endpoint_tol = m_cfg["graph_endpoint_tol"].asDouble() * units::cm;
            PR::g_graph_endpoint_policy.tol = pa.m_graph_endpoint_tol;
        }
    }

    // Verbatim copy of TaggerCheckSTM::load_trackfitting_config (:1024-1070):
    // the runtime TrackFitting parameter JSON, resolved on WIRECELL_PATH.
    void load_trackfitting_config(const std::string& config_file) {
        try {
            const std::string resolved = Persist::resolve(config_file);
            std::ifstream file(resolved.empty() ? config_file : resolved);
            if (!file.is_open()) {
                std::cerr << "CheckSTM_Michel: Cannot open config file: " << config_file
                          << " (not found on WIRECELL_PATH either)" << std::endl;
                return;
            }
            Json::Value root;
            Json::CharReaderBuilder builder;
            std::string errs;
            if (!Json::parseFromStream(builder, file, &root, &errs)) {
                std::cerr << "CheckSTM_Michel: Failed to parse JSON: " << errs << std::endl;
                return;
            }
            for (const auto& param_name : root.getMemberNames()) {
                if (param_name.substr(0, 1) == "_") continue;  // comments
                try {
                    double value = root[param_name].asDouble();
                    m_track_fitter->set_parameter(param_name, value);
                } catch (const std::exception& e) {
                    std::cerr << "CheckSTM_Michel: Failed to set parameter " << param_name
                              << ": " << e.what() << std::endl;
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "CheckSTM_Michel: Exception loading config: " << e.what() << std::endl;
        }
    }

    static std::string bits_string(unsigned bits) {
        if (!bits) return "STM";
        std::string s;
        auto add = [&](unsigned b, const char* n) { if (bits & b) { if (!s.empty()) s += "|"; s += n; } };
        add(R_NO_CHAIN, "no_chain"); add(R_STOP_UNMATCHED, "stop_unmatched");
        add(R_NO_BRAGG, "no_bragg"); add(R_SHAPE_FLAT, "shape_flat");
        add(R_NOT_MUON_PID, "not_muon_pid"); add(R_CONTINUATION, "continuation");
        add(R_STOP_NEAR_BOUNDARY, "stop_near_boundary"); add(R_VERTEX_HADRON, "vertex_hadron");
        add(R_SHORT, "short");
        return s;
    }

    // Read the STM tagger's accepted pass: entry (L == 0 row) and stop (row
    // kink_num, bounds-clamped as TaggerCheckSTM.cxx:2141).  False when the
    // PCs are absent (save_stm_fit off upstream) or no pass was accepted.
    bool read_stm_anchor(const Cluster& cluster, Record& rec) const {
        if (!cluster.has_pc("stm_pass") || !cluster.has_pc("stm_fit")) return false;
        const auto& pass_pc = cluster.get_pc("stm_pass");
        const auto& fit_pc = cluster.get_pc("stm_fit");
        if (pass_pc.empty() || fit_pc.empty()) return false;
        const auto v_pass = pass_pc.get("pass")->elements<int>();
        const auto v_status = pass_pc.get("status")->elements<int>();
        const auto v_kink = pass_pc.get("kink_num")->elements<int>();
        int chosen = -1;
        for (size_t i = 0; i < v_pass.size(); ++i) {
            if (v_status[i] != 0) continue;
            if (chosen < 0 || v_pass[i] < v_pass[chosen]) chosen = static_cast<int>(i);
        }
        if (chosen < 0) return false;
        rec.has_pass = 1;
        rec.pass = v_pass[chosen];
        rec.kink_num = v_kink[chosen];

        const auto fx = fit_pc.get("x")->elements<double>();
        const auto fy = fit_pc.get("y")->elements<double>();
        const auto fz = fit_pc.get("z")->elements<double>();
        const auto fL = fit_pc.get("L")->elements<double>();
        const auto fp = fit_pc.get("pass")->elements<int>();
        std::vector<size_t> rows;
        for (size_t i = 0; i < fp.size(); ++i) if (fp[i] == rec.pass) rows.push_back(i);
        if (rows.empty()) return false;
        size_t ientry = rows.front();
        for (auto r : rows) if (fL[r] < fL[ientry]) ientry = r;
        const int n = static_cast<int>(rows.size());
        const int k = (rec.kink_num >= 0 && rec.kink_num < n) ? rec.kink_num : n - 1;
        const size_t istop = rows[k];
        rec.entry_pt = Point(fx[ientry], fy[ientry], fz[ientry]);
        rec.tagger_stop_pt = Point(fx[istop], fy[istop], fz[istop]);
        rec.stop_pt = rec.tagger_stop_pt;
        return true;
    }

    // Nearest graph vertex of `cluster` to `pt`; if none within tol, split the
    // nearest segment there (PR::break_segment) so the entry/stop become
    // vertices the orientation walk and the PF tree can use.
    VertexPtr anchor_vertex(PatternAlgorithms& pa, Graph& g, Cluster& cluster, const Point& pt,
                            double tol, bool allow_split, double& out_dis) const {
        auto [vtx, dis] = pa.closest_cluster_vertex(g, cluster, pt);
        out_dis = dis;
        if (vtx && dis <= tol) return vtx;
        if (allow_split) {
            SegmentPtr best; double best_d = 1e9; Point best_p;
            for (auto& seg : pa.find_cluster_segments(g, cluster)) {
                if (!seg || seg->fits().size() < 4) continue;
                auto [d, p] = segment_get_closest_point(seg, pt, "fit", "main");
                if (d < best_d) { best_d = d; best = seg; best_p = p; }
            }
            if (best && best_d <= tol) {
                try {
                    auto [ok, pair, nvtx] = break_segment(g, best, best_p, particle_data(), m_recomb_model, m_dv,
                                                          1e9 * units::cm, get<bool>(m_cfg, "break_seg_orient", false));
                    if (ok && nvtx) {
                        // break_segment (PRSegmentFunctions.cxx:1105+) stamps the
                        // cluster on the two child segments but NOT on the new
                        // vertex.  A clusterless vertex is fatal downstream:
                        // examine_direction returns false at once
                        // (NeutrinoVertexFinder.cxx:1503) so nothing gets
                        // oriented, and fill_bee_pf_tree's main-cluster test
                        // (pf_track_main_cluster_only) then rejects every seed
                        // from the main vertex -- smoke 039252/2 cluster 86 lost
                        // its mu- node and its Michel fell back to a ROOT shower.
                        if (!nvtx->cluster()) nvtx->cluster(&cluster);
                        out_dis = best_d;
                        return nvtx;
                    }
                } catch (const std::exception& e) {
                    SPDLOG_LOGGER_WARN(s_log, "{}anchor_vertex: break_segment threw: {}", m_evt_tag, e.what());
                }
            }
        }
        return vtx;  // may be null
    }

    void set_pdg(const SegmentPtr& seg, int pdg) const {
        if (!seg) return;
        auto p4 = segment_cal_4mom(seg, pdg, particle_data(), m_recomb_model, m_mip_dqdx / units::cm);
        auto pinfo = std::make_shared<Aux::ParticleInfo>(pdg, particle_data()->get_particle_mass(pdg),
                                                         particle_data()->pdg_to_name(pdg), p4);
        seg->particle_info(pinfo);
        seg->particle_score(100.0);
    }

    void add_points(Record& rec, const SegmentPtr& seg, int role, const StmMichelProfile* prof = nullptr) const {
        const int seg_id = (seg && seg->cluster() ? seg->cluster()->get_cluster_id() : 0) * 1000
                         + (seg ? static_cast<int>(seg->get_graph_index()) : 0);
        if (prof) {
            for (size_t i = 0; i < prof->pts.size(); ++i) {
                rec.px.push_back(prof->pts[i].x()); rec.py.push_back(prof->pts[i].y()); rec.pz.push_back(prof->pts[i].z());
                rec.pq.push_back(prof->dQdx[i]); rec.pL.push_back(prof->L[i]); rec.prr.push_back(prof->rr[i]);
                rec.prole.push_back(role); rec.pseg.push_back(seg_id);
            }
            return;
        }
        if (!seg) return;
        for (const auto& f : seg->fits()) {
            if (f.dx <= 0) continue;
            rec.px.push_back(f.point.x()); rec.py.push_back(f.point.y()); rec.pz.push_back(f.point.z());
            rec.pq.push_back(f.dQ / (f.dx / units::cm)); rec.pL.push_back(-1); rec.prr.push_back(-1);
            rec.prole.push_back(role); rec.pseg.push_back(seg_id);
        }
    }

    // Units in the persisted rows (and hence in T_stm_michel): lengths and
    // coordinates in CM, energies in MeV, dQ/dx in e/cm, angles in degrees --
    // the generic PC->TTree writer has no unit knowledge, so the PC carries
    // the human units directly (unlike stm_fit, whose dedicated writer
    // divides by units::cm).
    void persist(Cluster& cluster, const Record& r) const {
        using WireCell::PointCloud::Array;
        using WireCell::PointCloud::Dataset;
        const double cm = units::cm;
        std::map<std::string, Array> a;
        auto I1 = [&](const char* k, int v)    { a.emplace(k, Array(std::vector<int>{v})); };
        auto D1 = [&](const char* k, double v) { a.emplace(k, Array(std::vector<double>{v})); };
        I1("cluster_id", r.cluster_id); I1("gid", r.gid); D1("t0_us", r.t0_us);
        I1("is_stm", r.reject_bits == 0 ? 1 : 0); I1("reject_bits", static_cast<int>(r.reject_bits));
        I1("has_pass", r.has_pass); I1("pass", r.pass); I1("kink_num", r.kink_num);
        D1("entry_x", r.entry_pt.x() / cm); D1("entry_y", r.entry_pt.y() / cm); D1("entry_z", r.entry_pt.z() / cm);
        D1("stop_x", r.stop_pt.x() / cm); D1("stop_y", r.stop_pt.y() / cm); D1("stop_z", r.stop_pt.z() / cm);
        D1("tagger_stop_x", r.tagger_stop_pt.x() / cm); D1("tagger_stop_y", r.tagger_stop_pt.y() / cm); D1("tagger_stop_z", r.tagger_stop_pt.z() / cm);
        I1("entry_vtx_id", r.entry_vtx_id); I1("stop_vtx_id", r.stop_vtx_id); D1("stop_dis", r.stop_dis / cm);
        I1("n_chain_segs", r.n_chain_segs); I1("n_profile_pts", r.n_profile_pts); D1("muon_len", r.muon_len / cm);
        I1("n_delta", r.n_delta); I1("n_body_other", r.n_body_other); I1("n_body_hadron", r.n_body_hadron); D1("delta_len", r.delta_len / cm);
        D1("ks_mu", r.ks_mu); D1("ks_flat", r.ks_flat); D1("ratio_mu", r.ratio_mu); D1("ratio_flat", r.ratio_flat);
        for (int i = 0; i < 4; ++i) {
            D1(("comp_fwd" + std::to_string(i)).c_str(), r.comp_fwd[i]);
            D1(("comp_bwd" + std::to_string(i)).c_str(), r.comp_bwd[i]);
        }
        D1("tail_med", r.bragg.tail_med); D1("plateau_med", r.bragg.plateau_med);
        D1("contrast", r.bragg.contrast); D1("contrast_expected", r.bragg.expected);
        I1("n_tail", r.bragg.n_tail); I1("n_plateau", r.bragg.n_plateau);
        I1("short_track", r.bragg.short_track ? 1 : 0); I1("bragg_valid", r.bragg.valid ? 1 : 0);
        I1("n_stop_arms", r.n_stop_arms); I1("michel_found", r.michel_found); I1("n_michel_segs", r.n_michel_segs);
        I1("michel_conn_type", r.michel_conn_type);
        D1("michel_len", r.michel_len / cm); D1("michel_mip", r.michel_mip); D1("michel_kink_deg", r.michel_kink_deg);
        D1("michel_far_len", r.michel_far_len / cm);
        D1("michel_ke_dqdx", r.michel_ke_dqdx); D1("michel_ke_range", r.michel_ke_range); D1("michel_ke_best", r.michel_ke_best);
        D1("cont_len", r.cont_len / cm); D1("cont_angle_deg", r.cont_angle_deg); D1("cont_mip", r.cont_mip);
        I1("n_dots", r.n_dots); I1("n_dot_clusters_unfit", r.n_dot_clusters_unfit);
        D1("dots_ke_dqdx", r.dots_ke_dqdx); D1("dots_charge_unfit", r.dots_charge_unfit);
        I1("in_fv", r.in_fv);
        cluster.local_pcs()["stm_michel"] = Dataset(a);

        if (!r.px.empty()) {
            std::map<std::string, Array> p;
            auto tocm = [&](const std::vector<double>& v) { std::vector<double> o(v); for (auto& x : o) if (x > -0.5) x /= cm; return o; };
            std::vector<double> x(r.px), y(r.py), z(r.pz);
            for (auto& v : x) v /= cm; for (auto& v : y) v /= cm; for (auto& v : z) v /= cm;
            p.emplace("x", Array(x)); p.emplace("y", Array(y)); p.emplace("z", Array(z));
            p.emplace("q", Array(r.pq)); p.emplace("L", Array(tocm(r.pL))); p.emplace("rr", Array(tocm(r.prr)));
            p.emplace("role", Array(r.prole)); p.emplace("seg_id", Array(r.pseg));
            cluster.local_pcs()["stm_michel_pts"] = Dataset(p);
        }
    }
};

void CheckSTM_Michel::visit(Ensemble& ensemble) const
{
    const auto t_total = Clock::now();
    m_evt_tag = ensemble.rse_valid()
        ? WireCell::String::format("evt%d ", ensemble.ident()) : std::string();

    // The member fitter only holds parameters, but it may have cached
    // geometry from a previous event through set_parameters copies -- reset
    // is the documented per-event contract (TrackFitting.h:439).
    m_track_fitter->reset_for_new_event();

    auto groupings = ensemble.with_name(m_grouping_name);
    if (groupings.empty()) return;
    auto& grouping = *groupings.at(0);

    // ---- candidates ---------------------------------------------------------
    std::vector<Cluster*> candidates;
    for (auto* cluster : grouping.children()) {
        if (!cluster->get_flag(Flags::main_cluster)) continue;
        if (cluster->get_flag(Flags::TGM)) continue;
        if (m_require_stm_flag && !cluster->get_flag(Flags::STM)) continue;
        if (!m_require_stm_flag && !cluster->has_pc("stm_pass")) continue;
        candidates.push_back(cluster);
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const Cluster* a, const Cluster* b) { return a->ident() < b->ident(); });
    if (static_cast<int>(candidates.size()) > m_max_candidates) {
        SPDLOG_LOGGER_INFO(s_log, "{}CheckSTM_Michel: {} candidates, keeping the first {} (max_candidates)",
                           m_evt_tag, candidates.size(), m_max_candidates);
        candidates.resize(m_max_candidates);
    }
    if (candidates.empty()) {
        SPDLOG_LOGGER_INFO(s_log, "{}CheckSTM_Michel: no STM-tagged main cluster; nothing to reconstruct", m_evt_tag);
        return;
    }

    auto fiducial_utils = grouping.get_fiducialutils();
    std::vector<double> fv_tol(6, -m_stop_fv_margin_cm * units::cm);
    auto mu_fn = particle_data() ? particle_data()->get_dEdx_function("muon") : nullptr;

    int n_stm = 0, n_michel = 0;
    for (size_t ci = 0; ci < candidates.size(); ++ci) {
        auto t0 = Clock::now();
        Cluster* main = candidates[ci];
        Record rec;
        rec.cluster_id = main->get_cluster_id();
        rec.gid = main->get_scalar<int>("matched_flash_gid", -1);
        rec.t0_us = main->get_cluster_t0() / units::us;

        // ---- the anchor from the STM tagger ------------------------------
        const bool anchored = read_stm_anchor(*main, rec);
        if (!anchored) {
            rec.reject_bits |= R_NO_CHAIN;
            SPDLOG_LOGGER_WARN(s_log, "{}CheckSTM_Michel: cluster {} has no accepted stm_pass / stm_fit PC (save_stm_fit off upstream?); verdict only",
                               m_evt_tag, rec.cluster_id);
            persist(*main, rec);
            continue;
        }

        // ---- companions: same flash bundle (the only clusters with a t0 and
        // hence a drift coordinate), near the stop, short ------------------
        std::vector<Cluster*> companions;
        std::vector<Cluster*> dot_clusters;
        if (rec.gid >= 0) {
            for (auto* oc : grouping.children()) {
                if (oc == main) continue;
                if (oc->get_scalar<int>("matched_flash_gid", -1) != rec.gid) continue;
                if (oc->get_length() > m_dot_max_len_cm * units::cm) continue;
                if (oc->npoints() == 0) continue;
                const auto [cp, blob] = oc->get_closest_point_blob(rec.stop_pt);
                if ((cp - rec.stop_pt).magnitude() > m_michel_dot_radius_cm * units::cm) continue;
                companions.push_back(oc);
            }
        }
        std::sort(companions.begin(), companions.end(),
                  [](const Cluster* a, const Cluster* b) { return a->ident() < b->ident(); });

        // ---- fitter + graph, as TaggerCheckNeutrino.cxx:2253-2324 ---------
        auto tf = std::make_shared<TrackFitting>(TrackFittingPresets::create_with_current_values());
        tf->set_parameters(m_track_fitter->get_parameters());
        tf->set_perf(m_perf);
        tf->set_detector_volume(m_dv);
        tf->set_pc_transforms(m_pcts);
        tf->set_parameter("fit_blob_coverage", m_fit_blob_coverage);
        tf->set_parameter("dqdx_fit_keep_all_points", m_dqdx_fit_keep_all_points ? 1.0 : 0.0);
        if (m_excl_t0_frame) tf->set_parameter("excl_t0_frame", 1.0);
        {
            std::vector<Cluster*> to_preload{main};
            for (auto* c : companions) to_preload.push_back(c);
            tf->preload_clusters(to_preload);
        }
        auto pr_graph = std::make_shared<Graph>();
        tf->add_graph(pr_graph);
        Graph& g = *pr_graph;

        PatternAlgorithms pa;
        apply_pattern_knobs(pa);

        // ---- the four PR stages on the main cluster ----------------------
        const bool ok_main = pa.find_proto_vertex(g, *main, *tf, m_dv, true, 2, true, particle_data());
        if (ok_main) {
            pa.clustering_points(g, *main, m_dv);
            pa.separate_track_shower(g, *main);
            pa.determine_direction(g, *main, particle_data(), m_recomb_model);
        }
        if (m_perf) SPDLOG_LOGGER_DEBUG(s_log, "{}CheckSTM_Michel timing: main PR took {} ms", m_evt_tag, MS(Clock::now() - t0).count());

        // ---- the same on the companions (TaggerCheckNeutrino.cxx:2927-2962) -
        for (auto* cluster : companions) {
            bool ok = false;
            if (cluster->get_length() > 6 * units::cm) {
                ok = pa.find_proto_vertex(g, *cluster, *tf, m_dv, true, 2, false);
            }
            else {
                ok = pa.find_proto_vertex(g, *cluster, *tf, m_dv, false, 1, false);
                if (!ok) { pa.init_point_segment(g, *cluster, *tf, m_dv); ok = true; }
            }
            if (!ok) continue;
            pa.clustering_points(g, *cluster, m_dv);
            pa.separate_track_shower(g, *cluster);
            pa.determine_direction(g, *cluster, particle_data(), m_recomb_model);
        }

        // ---- entry and stop vertices -------------------------------------
        VertexPtr entry_v, stop_v;
        double entry_dis = 1e9, stop_dis = 1e9;
        if (ok_main) {
            // The entry sits at the tagger's boundary end; the graph's own end
            // vertex is typically within a few cm of it, and splitting there
            // would leave a sub-cm stub beyond the entry (smoke 039252/2
            // cluster 86: a 0.8 cm "pi+" leaf).  So the entry snaps loosely;
            // the stop, whose position the Michel search keys on, snaps tightly.
            entry_v = anchor_vertex(pa, g, *main, rec.entry_pt, m_entry_snap_tol_cm * units::cm, true, entry_dis);
            stop_v  = anchor_vertex(pa, g, *main, rec.stop_pt,  m_stop_snap_tol_cm * units::cm, true, stop_dis);
            if (stop_v && stop_v == entry_v) stop_v = nullptr;
        }
        if (!entry_v) {
            rec.reject_bits |= R_NO_CHAIN;
        }
        else {
            // The PF root and the Bee/JSON "main vertex" -- the ENTRY point.
            entry_v->set_flags(VertexFlags::kNeutrinoVertex);
            tf->set_main_vertex(entry_v);
            rec.entry_vtx_id = rec.cluster_id * 1000 + static_cast<int>(entry_v->get_graph_index());
            // Orient everything outward from the entry (re-assigns dirsign on
            // every reached segment; flag_final so a prior strong direction
            // does not survive; provisional PID we override below).
            IndexedVertexSet v_lm; IndexedSegmentSet s_lm;
            pa.examine_direction(g, entry_v, entry_v, v_lm, s_lm, particle_data(), m_recomb_model, true);
        }

        // ---- the muon chain ----------------------------------------------
        std::vector<SegmentPtr> chain;
        if (entry_v && stop_v) chain = stm_michel_shortest_chain(g, entry_v, stop_v);
        if (entry_v && chain.empty()) {
            // No vertex near the tagger's stop (or unreachable): walk greedily
            // from the entry with find_cont_muon_segment, ignoring dQ/dx so the
            // Bragg segment is admitted, and stop at the vertex nearest the
            // tagger's stop point.
            rec.reject_bits |= R_STOP_UNMATCHED;
            SegmentPtr cur; double best_len = -1;
            for (auto e : sorted_out_edges(entry_v->get_descriptor(), g)) {
                auto s = g[e].segment; if (!s) continue;
                const double l = segment_track_length(s);
                if (l > best_len) { best_len = l; cur = s; }
            }
            VertexPtr vtx = entry_v;
            std::set<size_t> seen;
            while (cur && seen.insert(cur->get_graph_index()).second) {
                chain.push_back(cur);
                VertexPtr far = find_other_vertex(g, cur, vtx);
                if (!far) break;
                vtx = far;
                if ((stm_michel_vertex_point(far) - rec.stop_pt).magnitude() < m_stop_snap_tol_cm * units::cm) break;
                auto [nseg, nvtx] = pa.find_cont_muon_segment(g, cur, far, true);
                cur = nseg;
            }
            stop_v = vtx;
            if (stop_v == entry_v) stop_v = nullptr;
        }
        if (chain.empty()) rec.reject_bits |= R_NO_CHAIN;

        std::vector<VertexPtr> chain_vtxs = stm_michel_chain_vertices(g, chain, entry_v);
        if (!chain.empty() && chain_vtxs.size() != chain.size() + 1) { chain.clear(); chain_vtxs.clear(); rec.reject_bits |= R_NO_CHAIN; }
        if (stop_v) {
            rec.stop_vtx_id = rec.cluster_id * 1000 + static_cast<int>(stop_v->get_graph_index());
            rec.stop_pt = stm_michel_vertex_point(stop_v);
            rec.stop_dis = (rec.stop_pt - rec.tagger_stop_pt).magnitude();
        }
        rec.n_chain_segs = static_cast<int>(chain.size());

        StmMichelProfile prof;
        IndexedSegmentSet chain_set;
        if (!chain.empty()) {
            for (auto& s : chain) {
                chain_set.insert(s);
                s->unset_flags(SegmentFlags::kShowerTrajectory);
                s->unset_flags(SegmentFlags::kShowerTopology);
                set_pdg(s, 13);
            }
            prof = stm_michel_profile(g, chain, entry_v);
            rec.n_profile_pts = static_cast<int>(prof.L.size());
            rec.muon_len = prof.total_length;
            add_points(rec, chain.back(), 1, &prof);
            if (rec.n_profile_pts < m_min_chain_points) rec.reject_bits |= R_SHORT;
        }

        // ---- dQ/dx vs residual range: contrast, KS shape, template PID ----
        if (!prof.empty()) {
            std::function<double(double)> mu_at = nullptr;
            if (mu_fn) mu_at = [mu_fn](double rr_cm) { return mu_fn->scalar_function(rr_cm); };
            rec.bragg = stm_michel_bragg_contrast(prof, mu_at,
                                                  m_bragg_tail_lo_cm * units::cm, m_bragg_tail_hi_cm * units::cm,
                                                  m_bragg_plateau_lo_cm * units::cm, m_bragg_plateau_hi_cm * units::cm);
            if (!rec.bragg.valid || rec.bragg.expected <= 0 ||
                rec.bragg.contrast < m_bragg_contrast_min * rec.bragg.expected) {
                rec.reject_bits |= R_NO_BRAGG;
            }
            // The TaggerCheckSTM::eval_stm_core_impl recipe (:2780-2797) over
            // the last compare_range of residual range, e/cm frame.
            if (mu_fn) {
                std::vector<double> test, ref_mu, ref_flat;
                for (size_t i = 0; i < prof.rr.size(); ++i) {
                    if (prof.rr[i] > m_compare_range_cm * units::cm) continue;
                    test.push_back(prof.dQdx[i]);
                    ref_mu.push_back(mu_fn->scalar_function(prof.rr[i] / units::cm + m_offset_length_cm));
                    ref_flat.push_back(m_mip_dqdx);
                }
                if (test.size() >= 3) {
                    auto sum = [](const std::vector<double>& v) { double s = 0; for (double x : v) s += x; return s; };
                    rec.ks_mu = WireCell::kslike_compare(test, ref_mu);
                    rec.ratio_mu = sum(ref_mu) / (sum(test) + 1e-9);
                    rec.ks_flat = WireCell::kslike_compare(test, ref_flat);
                    rec.ratio_flat = sum(ref_flat) / (sum(test) + 1e-9);
                    if (rec.ks_mu + m_ks_margin >= rec.ks_flat) rec.reject_bits |= R_SHAPE_FLAT;
                }
                // do_track_comp: internal-unit arrays; forward = stop at L.back().
                std::vector<double> L_f(prof.L), q_f(prof.dQdx.size());
                for (size_t i = 0; i < q_f.size(); ++i) q_f[i] = prof.dQdx[i] / units::cm;
                std::vector<double> L_b(L_f.size()), q_b(q_f.size());
                const size_t n = L_f.size();
                for (size_t i = 0; i < n; ++i) { L_b[i] = prof.total_length - L_f[n - 1 - i]; q_b[i] = q_f[n - 1 - i]; }
                rec.comp_fwd = do_track_comp(L_f, q_f, m_compare_range_cm * units::cm, m_offset_length_cm * units::cm,
                                             particle_data(), m_mip_dqdx / units::cm);
                rec.comp_bwd = do_track_comp(L_b, q_b, m_compare_range_cm * units::cm, m_offset_length_cm * units::cm,
                                             particle_data(), m_mip_dqdx / units::cm);
                if (rec.comp_fwd.size() == 4) {
                    const bool muon_like = rec.comp_fwd[0] > 0.5 && rec.comp_fwd[1] < rec.comp_fwd[2] && rec.comp_fwd[1] < rec.comp_fwd[3];
                    if (!muon_like) rec.reject_bits |= R_NOT_MUON_PID;
                }
            }
        }

        // ---- arms: body (delta rays / hadrons) and stop (Michel / continuation)
        StmMichelArmThresholds th;
        th.mip_dqdx_median = m_mip_dqdx_median / units::cm;
        th.continuation_max_angle_deg = m_continuation_max_angle_deg;
        th.continuation_min_len = m_continuation_min_len_cm * units::cm;
        th.continuation_mip_lo = m_continuation_mip_lo;
        th.continuation_mip_hi = m_continuation_mip_hi;
        th.michel_max_len = m_michel_max_len_cm * units::cm;
        th.michel_mip_hi = m_michel_mip_hi;
        th.michel_mip_lo = m_michel_mip_lo;
        th.michel_min_kink_deg = m_michel_min_kink_deg;
        th.delta_max_len = m_delta_max_len_cm * units::cm;
        th.hadron_mip = m_vertex_hadron_mip;

        IndexedShowerSet showers;
        std::shared_ptr<Shower> michel_shower;
        std::vector<StmMichelArm> michel_arms;
        if (!chain.empty()) {
            // interior vertices
            for (size_t vi = 1; vi + 1 < chain_vtxs.size(); ++vi) {
                VertexPtr v = chain_vtxs[vi];
                SegmentPtr in_seg = chain[vi - 1], out_seg = chain[vi];
                for (auto e : sorted_out_edges(v->get_descriptor(), g)) {
                    auto arm = g[e].segment;
                    if (!arm || arm == in_seg || arm == out_seg || chain_set.count(arm)) continue;
                    auto a = stm_michel_classify_chain_arm(g, in_seg, arm, v, th);
                    if (a.kind == StmMichelArm::kDelta) {
                        ++rec.n_delta; rec.delta_len += a.len;
                        set_pdg(arm, 11);
                        add_points(rec, arm, 2);
                    }
                    else if (a.kind == StmMichelArm::kHadron) {
                        ++rec.n_body_hadron;
                        rec.reject_bits |= R_VERTEX_HADRON;
                    }
                    else {
                        ++rec.n_body_other;
                    }
                }
            }
            // the stop
            if (stop_v) {
                SegmentPtr last = chain.back();
                for (auto e : sorted_out_edges(stop_v->get_descriptor(), g)) {
                    auto arm = g[e].segment;
                    if (!arm || arm == last || chain_set.count(arm)) continue;
                    ++rec.n_stop_arms;
                    auto a = stm_michel_classify_stop_arm(g, last, arm, stop_v, th);
                    if (a.kind == StmMichelArm::kContinuation) {
                        rec.reject_bits |= R_CONTINUATION;
                        if (a.len > rec.cont_len) { rec.cont_len = a.len; rec.cont_angle_deg = a.kink_deg; rec.cont_mip = a.mip; }
                    }
                    else if (a.kind == StmMichelArm::kMichel) {
                        michel_arms.push_back(a);
                    }
                }
            }
        }

        // ---- the Michel shower ---------------------------------------------
        if (!michel_arms.empty() && stop_v) {
            std::sort(michel_arms.begin(), michel_arms.end(),
                      [](const StmMichelArm& a, const StmMichelArm& b) {
                          if (a.len != b.len) return a.len > b.len;
                          return a.seg->get_graph_index() < b.seg->get_graph_index();
                      });
            const auto& seed = michel_arms.front();
            rec.michel_found = 1;
            rec.michel_len = seed.len; rec.michel_mip = seed.mip; rec.michel_kink_deg = seed.kink_deg; rec.michel_far_len = seed.far_len;
            rec.michel_conn_type = 1;
            for (auto& a : michel_arms) { set_pdg(a.seg, 11); add_points(rec, a.seg, 3); }
            if (m_build_michel_shower) {
                michel_shower = std::make_shared<Shower>(g);
                michel_shower->set_start_vertex(stop_v, 1);
                michel_shower->set_start_segment(seed.seg, false, "fit", "associate_points");
                IndexedSegmentSet used(chain_set);
                michel_shower->complete_structure_with_start_segment(used, "fit", "associate_points", true);
                // every member is the electron's
                IndexedVertexSet mv; IndexedSegmentSet ms;
                michel_shower->fill_sets(mv, ms, false);
                for (auto& s : ms) { if (!chain_set.count(s)) set_pdg(s, 11); }
                rec.n_michel_segs = static_cast<int>(ms.size());
                michel_shower->set_particle_type(11);
                michel_shower->calculate_kinematics(particle_data(), m_recomb_model);
                rec.michel_ke_dqdx = michel_shower->get_kine_dQdx() / units::MeV;
                rec.michel_ke_range = michel_shower->get_kine_range() / units::MeV;
                // range energy is meaningless for an electron
                michel_shower->set_kine_best(michel_shower->get_kine_dQdx());
                rec.michel_ke_best = michel_shower->get_kine_best() / units::MeV;
            }
            else {
                rec.n_michel_segs = static_cast<int>(michel_arms.size());
                rec.michel_ke_dqdx = segment_cal_kine_dQdx(seed.seg, m_recomb_model) / units::MeV;
                rec.michel_ke_best = rec.michel_ke_dqdx;
            }
        }

        // ---- dots: companion segments near the stop, closer to the stop
        // than to the muon body ------------------------------------------------
        if (stop_v && !companions.empty()) {
            std::set<int> fitted_ids;
            for (auto* oc : companions) {
                auto segs = pa.find_cluster_segments(g, *oc);
                if (segs.empty()) {
                    ++rec.n_dot_clusters_unfit;
                    double q = 0; for (const auto* b : oc->children()) q += b->charge();
                    rec.dots_charge_unfit += q;
                    continue;
                }
                fitted_ids.insert(oc->get_cluster_id());
                for (auto& seg : segs) {
                    auto [d_stop, cp] = segment_get_closest_point(seg, rec.stop_pt, "fit", "main");
                    if (d_stop > m_michel_dot_radius_cm * units::cm) continue;
                    if (segment_track_length(seg) > m_dot_max_len_cm * units::cm) continue;
                    // distance to the muon body beyond its last dot_body_exclusion
                    double d_body = 1e9;
                    for (size_t i = 0; i < prof.pts.size(); ++i) {
                        if (prof.rr[i] < m_dot_body_exclusion_cm * units::cm) continue;
                        d_body = std::min(d_body, (prof.pts[i] - cp).magnitude());
                    }
                    if (d_body < d_stop) continue;   // a delta ray / body fragment, not a Michel dot
                    ++rec.n_dots;
                    rec.dots_ke_dqdx += segment_cal_kine_dQdx(seg, m_recomb_model) / units::MeV;
                    set_pdg(seg, 11);
                    add_points(rec, seg, 4);
                    if (m_build_michel_shower) {
                        if (!michel_shower) {
                            michel_shower = std::make_shared<Shower>(g);
                            michel_shower->set_start_vertex(stop_v, 2);
                            michel_shower->set_start_segment(seg, false, "fit", "associate_points");
                            michel_shower->set_particle_type(11);
                            rec.michel_conn_type = 2;
                        }
                        else {
                            michel_shower->add_segment(seg, true, "fit", "associate_points");
                        }
                    }
                }
            }
            if (michel_shower && rec.michel_conn_type == 2) {
                michel_shower->calculate_kinematics(particle_data(), m_recomb_model);
                michel_shower->set_kine_best(michel_shower->get_kine_dQdx());
                rec.michel_ke_best = michel_shower->get_kine_best() / units::MeV;
            }
        }
        if (michel_shower) showers.insert(michel_shower);
        tf->set_showers(showers);

        // ---- containment of the stop -----------------------------------------
        if (fiducial_utils) {
            rec.in_fv = fiducial_utils->inside_fiducial_volume(rec.stop_pt, fv_tol) ? 1 : 0;
            if (!rec.in_fv) rec.reject_bits |= R_STOP_NEAR_BOUNDARY;
        }

        // ---- publish (TaggerCheckNeutrino.cxx:3580-3590) --------------------
        tf->assemble_fitted_charge_2d();
        if (ci == 0) grouping.set_track_fitting(tf);
        if (m_publish_nu_slots) grouping.set_track_fitting("nu" + std::to_string(ci), tf);
        persist(*main, rec);

        if (rec.reject_bits == 0) ++n_stm;
        if (rec.michel_found) ++n_michel;
        SPDLOG_LOGGER_INFO(s_log,
            "{}CheckSTM_Michel: cluster {} gid {} verdict {} bits {} | chain {} segs {:.1f} cm ({} pts) stop_dis {:.1f} cm | "
            "contrast {:.2f}/{:.2f} ks_mu {:.3f} ks_flat {:.3f} comp_fwd {:.0f}/{:.2f}/{:.2f}/{:.2f} | "
            "delta {} hadron {} | michel {} ({} segs, {:.1f} cm, {:.2f} mip, kink {:.0f} deg, {:.1f} MeV) dots {} ({:.1f} MeV) unfit {} | cont {:.1f} cm @ {:.0f} deg | in_fv {} | {:.0f} ms",
            m_evt_tag, rec.cluster_id, rec.gid, bits_string(rec.reject_bits), rec.reject_bits,
            rec.n_chain_segs, rec.muon_len / units::cm, rec.n_profile_pts, rec.stop_dis / units::cm,
            rec.bragg.contrast, rec.bragg.expected, rec.ks_mu, rec.ks_flat,
            rec.comp_fwd[0], rec.comp_fwd[1], rec.comp_fwd[2], rec.comp_fwd[3],
            rec.n_delta, rec.n_body_hadron,
            rec.michel_found, rec.n_michel_segs, rec.michel_len / units::cm, rec.michel_mip, rec.michel_kink_deg, rec.michel_ke_best,
            rec.n_dots, rec.dots_ke_dqdx, rec.n_dot_clusters_unfit,
            rec.cont_len / units::cm, rec.cont_angle_deg, rec.in_fv,
            MS(Clock::now() - t0).count());
    }

    SPDLOG_LOGGER_INFO(s_log, "{}CheckSTM_Michel: {} candidate(s), {} pass every check, {} with a Michel; {:.0f} ms",
                       m_evt_tag, candidates.size(), n_stm, n_michel, MS(Clock::now() - t_total).count());
}
