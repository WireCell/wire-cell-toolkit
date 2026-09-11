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
#include "WireCellMcs/MuonMCS.h"
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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
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
        // doc pdvd/88 (doc 78 action item 7b): on a candidate where the
        // stop-local keep fired, snap the stop only onto what the entry can reach.
        m_stop_snap_reachable = get<bool>(config, "stop_snap_reachable", m_stop_snap_reachable);
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
        m_companion_max_len_cm = get<double>(config, "companion_max_len_cm", m_companion_max_len_cm);
        m_michel_unfit_recom = get<double>(config, "michel_unfit_recom", m_michel_unfit_recom);
        m_michel_unfit_fudge = get<double>(config, "michel_unfit_fudge", m_michel_unfit_fudge);
        m_michel_unfit_w_ev = get<double>(config, "michel_unfit_w_ev", m_michel_unfit_w_ev);
        // doc pdhd/17: derive the unfitted-charge survival FROM the bound
        // recombination model at michel_unfit_dedx instead of the hard-coded
        // michel_unfit_recom x michel_unfit_fudge pair.  C++ default false =>
        // an absent key leaves dots_ke_unfit where doc pdhd/15 put it.
        m_michel_unfit_from_model = get<bool>(config, "michel_unfit_from_model", m_michel_unfit_from_model);
        m_michel_unfit_dedx = get<double>(config, "michel_unfit_dedx", m_michel_unfit_dedx);
        m_dot_body_exclusion_cm = get<double>(config, "dot_body_exclusion_cm", m_dot_body_exclusion_cm);
        // doc pdvd/51: the muon-capture gamma at the stop.
        m_stop_gamma_enable = get<bool>(config, "stop_gamma_enable", m_stop_gamma_enable);
        m_stop_gamma_radius_cm = get<double>(config, "stop_gamma_radius_cm", m_stop_gamma_radius_cm);
        m_stop_gamma_max_len_cm = get<double>(config, "stop_gamma_max_len_cm", m_stop_gamma_max_len_cm);
        m_stop_gamma_min_ke_mev = get<double>(config, "stop_gamma_min_ke_mev", m_stop_gamma_min_ke_mev);
        m_stop_gamma_max_ke_mev = get<double>(config, "stop_gamma_max_ke_mev", m_stop_gamma_max_ke_mev);
        m_stop_gamma_max_n = get<int>(config, "stop_gamma_max_n", m_stop_gamma_max_n);
        // doc pdvd/85 (doc 78 action item 5): a capture gamma only on a stopper.
        m_stop_gamma_require_stm = get<bool>(config, "stop_gamma_require_stm", m_stop_gamma_require_stm);
        // doc pdvd/53: the survey.
        m_survey_enable = get<bool>(config, "survey_enable", m_survey_enable);
        m_survey_radius_cm = get<double>(config, "survey_radius_cm", m_survey_radius_cm);
        m_survey_max_len_cm = get<double>(config, "survey_max_len_cm", m_survey_max_len_cm);
        m_delta_max_len_cm = get<double>(config, "delta_max_len_cm", m_delta_max_len_cm);
        m_vertex_hadron_mip = get<double>(config, "vertex_hadron_mip", m_vertex_hadron_mip);
        m_profile_min_dqdx_frac = get<double>(config, "profile_min_dqdx_frac", m_profile_min_dqdx_frac);
        m_pid_mode = get<int>(config, "pid_mode", m_pid_mode);
        m_plateau_mip_lo = get<double>(config, "plateau_mip_lo", m_plateau_mip_lo);
        m_plateau_mip_hi = get<double>(config, "plateau_mip_hi", m_plateau_mip_hi);
        m_stop_extend_max = get<int>(config, "stop_extend_max", m_stop_extend_max);
        // doc pdvd/57
        m_stop_retreat_max = get<int>(config, "stop_retreat_max", m_stop_retreat_max);
        m_retreat_collapse_frac = get<double>(config, "retreat_collapse_frac", m_retreat_collapse_frac);
        m_retreat_peak_frac = get<double>(config, "retreat_peak_frac", m_retreat_peak_frac);
        m_retreat_peak_window_cm = get<double>(config, "retreat_peak_window_cm", m_retreat_peak_window_cm);
        // doc pdvd/58 (T1c)
        m_stop_split_max = get<int>(config, "stop_split_max", m_stop_split_max);
        m_split_kink_min_deg = get<double>(config, "split_kink_min_deg", m_split_kink_min_deg);
        m_split_min_drop_cm = get<double>(config, "split_min_drop_cm", m_split_min_drop_cm);
        m_split_collapse_frac = get<double>(config, "split_collapse_frac", m_split_collapse_frac);
        m_split_peak_frac = get<double>(config, "split_peak_frac", m_split_peak_frac);
        m_split_peak_window_cm = get<double>(config, "split_peak_window_cm", m_split_peak_window_cm);
        m_split_dir_window_cm = get<double>(config, "split_dir_window_cm", m_split_dir_window_cm);
        // doc pdvd/61 (T2c): veto an attached (conn_type 1) Michel whose stop was
        // moved this event (a retreat or split fired) when its assembled KE
        // cannot support the claim -- the mechanism that moves the stop is not
        // itself a Michel test, and 5 of the items it moves the stop on are
        // through-going tracks that pick up a spurious attached arm (docs 57/58/59).
        m_moved_stop_michel_guard = get<bool>(config, "moved_stop_michel_guard", m_moved_stop_michel_guard);
        m_moved_stop_michel_ke_min = get<double>(config, "moved_stop_michel_ke_min", m_moved_stop_michel_ke_min);
        // doc pdvd/62 (T3): the Michel with no trajectory.  A: keep the
        // pr54 isolated residual when it sits within stop_local_residual_cm
        // of the tagger's stop (0 = off); B: admit a DISCONNECTED same-cluster
        // piece near the stop into the Michel object the way a companion
        // piece is; C: doc 55 sec 15.1's range-energy guard on a bridged /
        // charge-only Michel (dis > michel_range_energy_dis_cm and KE <
        // michel_range_energy_ke_min cannot be an electron born at the stop).
        m_stop_local_residual_cm = get<double>(config, "stop_local_residual_cm", m_stop_local_residual_cm);
        // doc pdvd/87 (doc 78 action item 7): a size floor on A's keep --
        // Steiner terminals and fitted length (cm); 0 = no floor.  Read only
        // when stop_local_residual_cm > 0.
        m_stop_local_residual_min_points = get<int>(config, "stop_local_residual_min_points", m_stop_local_residual_min_points);
        m_stop_local_residual_min_len_cm = get<double>(config, "stop_local_residual_min_len_cm", m_stop_local_residual_min_len_cm);
        m_stop_local_michel_pieces = get<bool>(config, "stop_local_michel_pieces", m_stop_local_michel_pieces);
        m_michel_range_energy_guard = get<bool>(config, "michel_range_energy_guard", m_michel_range_energy_guard);
        m_michel_range_energy_dis_cm = get<double>(config, "michel_range_energy_dis_cm", m_michel_range_energy_dis_cm);
        m_michel_range_energy_ke_min = get<double>(config, "michel_range_energy_ke_min", m_michel_range_energy_ke_min);
        // doc pdvd/64 (T6): give every kOther arm (interior or stop) a role-7 row
        // in stm_michel_pts so display and scan can see what the chain fitted
        // and classified as neither delta, hadron, Michel nor continuation.
        m_publish_other_arms = get<bool>(config, "publish_other_arms", m_publish_other_arms);
        m_segment_census = get<bool>(config, "segment_census", m_segment_census);   // doc pdvd/80
        // doc pdvd/65 (T7): read the two shape tests from a residual-range origin
        // anchored on the profile's Bragg peak (the prototype eval_stm recipe)
        // rather than on the chain's geometric far end.
        m_bragg_peak_anchor = get<bool>(config, "bragg_peak_anchor", m_bragg_peak_anchor);
        m_bragg_peak_search_cm = get<double>(config, "bragg_peak_search_cm", m_bragg_peak_search_cm);
        // doc pdvd/70 (P1): topology-first stop evidence -- a Michel of
        // sufficient quality at the stop clears the two dQ/dx-shape bits.
        m_topology_stop_evidence = get<bool>(config, "topology_stop_evidence", m_topology_stop_evidence);
        m_topology_michel_ke_min = get<double>(config, "topology_michel_ke_min", m_topology_michel_ke_min);
        m_topology_michel_len_min_cm = get<double>(config, "topology_michel_len_min_cm", m_topology_michel_len_min_cm);
        m_topology_clears_sparse = get<bool>(config, "topology_clears_sparse", m_topology_clears_sparse);
        // doc pdvd/71 (P4): collect the Michel's isolated gamma blobs as role-4
        // members, with their own energy -- michel_ke_best is never touched.
        m_michel_gamma_collect = get<bool>(config, "michel_gamma_collect", m_michel_gamma_collect);
        m_michel_gamma_radius_cm = get<double>(config, "michel_gamma_radius_cm", m_michel_gamma_radius_cm);
        m_michel_gamma_max_len_cm = get<double>(config, "michel_gamma_max_len_cm", m_michel_gamma_max_len_cm);
        m_michel_gamma_cos_min = get<double>(config, "michel_gamma_cos_min", m_michel_gamma_cos_min);
        m_michel_gamma_max_ke_mev = get<double>(config, "michel_gamma_max_ke_mev", m_michel_gamma_max_ke_mev);
        m_michel_gamma_total_ke_max_mev = get<double>(config, "michel_gamma_total_ke_max_mev", m_michel_gamma_total_ke_max_mev);
        // doc pdvd/81: the charge-based Michel energy (2-D charge minus the muon
        // fit's prediction) and the Michel / STM 2-D cell table.  Rows and
        // branches only; michel_ke_best is never touched.
        m_michel_q2d = get<bool>(config, "michel_q2d", m_michel_q2d);
        m_michel_q2d_cells = get<bool>(config, "michel_q2d_cells", m_michel_q2d_cells);
        m_michel_q2d_dis_cm = get<double>(config, "michel_q2d_dis_cm", m_michel_q2d_dis_cm);
        m_michel_q2d_stm_window_cm = get<double>(config, "michel_q2d_stm_window_cm", m_michel_q2d_stm_window_cm);
        // doc pdvd/72 (P3b): the moved-stop veto (T2c) spares an attached
        // Michel that turns at least this hard at the stop.  -1 = off.
        m_moved_stop_michel_kink_min = get<double>(config, "moved_stop_michel_kink_min", m_moved_stop_michel_kink_min);
        // doc pdvd/84 (doc 78 item 3): ... or that reaches at least this far
        // (arm length + far subtree, cm).  -1 = off.
        m_moved_stop_michel_reach_min_cm = get<double>(config, "moved_stop_michel_reach_min_cm", m_moved_stop_michel_reach_min_cm);
        // doc pdvd/73 (P2): PDVD operating points for the attached Michel gate.
        m_michel_mip_lo_turned = get<double>(config, "michel_mip_lo_turned", m_michel_mip_lo_turned);
        m_michel_mip_lo_turned_kink_deg = get<double>(config, "michel_mip_lo_turned_kink_deg", m_michel_mip_lo_turned_kink_deg);
        m_michel_far_len_shower_max_cm = get<double>(config, "michel_far_len_shower_max_cm", m_michel_far_len_shower_max_cm);
        m_michel_kink_window_cm = get<double>(config, "michel_kink_window_cm", m_michel_kink_window_cm);
        // doc pdvd/74 (P3): how the stop retreat reads the dropped tail, and
        // whether the retreat / split may run on a Bragg-confirmed chain.
        m_retreat_tail_strict = get<bool>(config, "retreat_tail_strict", m_retreat_tail_strict);
        m_retreat_tail_sublive = get<bool>(config, "retreat_tail_sublive", m_retreat_tail_sublive);
        m_michel_collinear_split = get<bool>(config, "michel_collinear_split", m_michel_collinear_split);
        // doc pdvd/82 (doc 78 action item 2): the peak-relative collapsed-tail
        // admission shared by the retreat and the split.
        m_stop_tail_peak_frac = get<double>(config, "stop_tail_peak_frac", m_stop_tail_peak_frac);
        m_stop_tail_peak_kink_min_deg = get<double>(config, "stop_tail_peak_kink_min_deg", m_stop_tail_peak_kink_min_deg);
        // doc pdvd/83 (doc 78 action item 4): the Michel gate offered to an arm
        // leaving the chain within this many cm before the stop.
        m_michel_near_stop_arm_cm = get<double>(config, "michel_near_stop_arm_cm", m_michel_near_stop_arm_cm);
        // doc pdvd/75 (P1b): when the peak anchor (T7) rejects a profile on a
        // shape bit although its anchored peak is prominent, and the same
        // tests pass at the geometric origin, the geometric reading stands.
        m_bragg_anchor_geo_fallback = get<bool>(config, "bragg_anchor_geo_fallback", m_bragg_anchor_geo_fallback);
        m_bragg_anchor_rise_min = get<double>(config, "bragg_anchor_rise_min", m_bragg_anchor_rise_min);
        // doc pdvd/66 (T8): the profile-geometry fields are always written; the guard
        // turns them into a reject bit (R_PROFILE_GEOMETRY).
        m_profile_geometry_guard = get<bool>(config, "profile_geometry_guard", m_profile_geometry_guard);
        m_profile_arc_span_max = get<double>(config, "profile_arc_span_max", m_profile_arc_span_max);
        m_unsupported_min_len_cm = get<double>(config, "unsupported_min_len_cm", m_unsupported_min_len_cm);
        m_unsupported_frac = get<double>(config, "unsupported_frac", m_unsupported_frac);
        m_end_window_cm = get<double>(config, "end_window_cm", m_end_window_cm);
        m_dead_volume_check = get<bool>(config, "dead_volume_check", m_dead_volume_check);
        m_min_chain_coverage = get<double>(config, "min_chain_coverage", m_min_chain_coverage);
        m_michel_guards_stop = get<bool>(config, "michel_guards_stop", m_michel_guards_stop);
        m_michel_shower_min_kink_deg = get<double>(config, "michel_shower_min_kink_deg", m_michel_shower_min_kink_deg);
        m_absorb_bragg_stub = get<bool>(config, "absorb_bragg_stub", m_absorb_bragg_stub);
        m_stop_fv_use_config_tolerance = get<bool>(config, "stop_fv_use_config_tolerance", m_stop_fv_use_config_tolerance);
        m_fv_tolerance.clear();
        if (config.isMember("fv_tolerance") && config["fv_tolerance"].isArray())
            for (const auto& t : config["fv_tolerance"]) m_fv_tolerance.push_back(t.asDouble());
        m_coverage_radius_cm = get<double>(config, "coverage_radius_cm", m_coverage_radius_cm);

        // TrackFitting parameters carried by TaggerCheckNeutrino's config
        // rather than by the runtime JSON (TaggerCheckNeutrino.cxx:2777-2791).
        m_fit_blob_coverage = get<double>(config, "fit_blob_coverage", m_fit_blob_coverage);
        m_dqdx_fit_keep_all_points = get<bool>(config, "dqdx_fit_keep_all_points", m_dqdx_fit_keep_all_points);
        m_excl_t0_frame = get<bool>(config, "excl_t0_frame", m_excl_t0_frame);

        // doc pdhd/16
        m_mcs_enable = get<bool>(config, "mcs_enable", m_mcs_enable);
        m_mcs_min_len_cm = get<double>(config, "mcs_min_len_cm", m_mcs_min_len_cm);
        m_mcs_cathode_x = get<double>(config, "mcs_cathode_x", m_mcs_cathode_x);
        m_mcs_cathode_xcut = get<double>(config, "mcs_cathode_xcut", m_mcs_cathode_xcut);

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
        // doc pdvd/88 (doc pdvd/78 action item 7b): default off.  The stop
        // snap takes the main cluster's nearest vertex (or splits its nearest
        // segment) with no test that the entry can reach it.  A pr54 residual
        // the stop-local keep (stop_local_residual_cm) added is disconnected
        // by construction, and doc pdvd/87 sec 6.2 traced it capturing the
        // stop: 039349_36/63's kept residual ends 0.64 cm from the tagger's
        // stop, nearer than the chain's own end (2.29 cm) -> empty chain ->
        // stop_unmatched.  When on, and ONLY on a candidate where that keep
        // fired on the main cluster, the snap considers only vertices and
        // segments the entry reaches; every other candidate runs the legacy
        // snap.  stop_snap_skipped (written only when on) = 1 where the
        // legacy snap's nearest vertex was unreachable.
        cfg["stop_snap_reachable"] = m_stop_snap_reachable;
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
        // doc pdhd/15 sec 3: the per-PIECE length cap at the dot loop.  Default
        // raised 10 -> 25 cm (= michel_max_len_cm, the ceiling an ATTACHED Michel
        // arm already faces) so a detached Michel is judged by the same size rule
        // as an attached one.  DECLARED DEFAULT CHANGE: a job that does not set
        // this key changes behaviour.  Set it back to 10 with
        // companion_max_len_cm 10 to reproduce the doc pdhd/14 tree.
        cfg["dot_max_len_cm"] = m_dot_max_len_cm;
        // doc pdhd/15 sec 3: which same-bundle clusters enter the PR at all.
        // Was the same knob as the piece cap, so raising one loosened the other
        // (doc pdhd/13 defect D2).  PDVD 039252_15 cluster 77's Michel is a
        // 20.1 cm cluster 1.91 cm from the stop: the 10 cm cap kept it out of
        // `companions` entirely, so it could become neither an arm nor a dot.
        cfg["companion_max_len_cm"] = m_companion_max_len_cm;
        // doc pdhd/15 sec 6: charge -> energy for a companion cluster the fitter
        // never reached (no dx, so no dQ/dx to invert).  The KineChargeOptions
        // TRACK pair; measured against the chain's own segment_cal_kine_dQdx at
        // ratio median 1.19 (PDVD) / 1.01 (PDHD), against 1.98 / 1.67 for the
        // SHOWER pair (0.5 / 0.8).
        cfg["michel_unfit_recom"] = m_michel_unfit_recom;
        cfg["michel_unfit_fudge"] = m_michel_unfit_fudge;
        cfg["michel_unfit_w_ev"] = m_michel_unfit_w_ev;
        cfg["michel_unfit_from_model"] = m_michel_unfit_from_model;   // doc pdhd/17
        cfg["michel_unfit_dedx"] = m_michel_unfit_dedx;
        cfg["dot_body_exclusion_cm"] = m_dot_body_exclusion_cm;
        // doc pdvd/51.  stop_gamma_enable false reproduces the doc pdhd/17 tree.
        cfg["stop_gamma_enable"] = m_stop_gamma_enable;
        cfg["stop_gamma_radius_cm"] = m_stop_gamma_radius_cm;
        cfg["stop_gamma_max_len_cm"] = m_stop_gamma_max_len_cm;
        cfg["stop_gamma_min_ke_mev"] = m_stop_gamma_min_ke_mev;
        cfg["stop_gamma_max_ke_mev"] = m_stop_gamma_max_ke_mev;
        cfg["stop_gamma_max_n"] = m_stop_gamma_max_n;
        // doc pdvd/85.  false reproduces the doc pdvd/84 tree: the capture
        // gamma is published on every candidate, stopper or not.
        cfg["stop_gamma_require_stm"] = m_stop_gamma_require_stm;
        // doc pdvd/53.  survey_enable false reproduces the doc pdvd/51 tree,
        // stm_michel_pts included (no role-6 rows, no rej/d_stop/d_body columns).
        cfg["survey_enable"] = m_survey_enable;
        cfg["survey_radius_cm"] = m_survey_radius_cm;
        cfg["survey_max_len_cm"] = m_survey_max_len_cm;
        cfg["delta_max_len_cm"] = m_delta_max_len_cm;
        cfg["vertex_hadron_mip"] = m_vertex_hadron_mip;
        // doc pdhd/03: profile points with dQ/dx < frac * mip_dqdx are DEAD
        // (a cell the fit could not read) and are dropped from the verdict
        // metrics (Bragg medians, KS shape, template PID).  0 = keep every
        // point (the doc pdvd/48 behaviour).
        cfg["profile_min_dqdx_frac"] = m_profile_min_dqdx_frac;
        // doc pdhd/03 sec 6: what the template PID's "muon" verdict requires.
        //  0 (doc pdvd/48): do_track_comp's direction gate AND the muon
        //    template beats BOTH the proton and the electron template.
        //  1: the gate AND muon beats proton.
        //  2: muon beats proton only (the proton veto).
        // The electron table is near-flat at ~MIP, so on a real stopping muon
        // whose Bragg rise the fit smears it scores within 20 % of the muon
        // template and often below it; the gate itself fails on textbook
        // muons whose last-35-cm rise is weaker than the table's (PDHD
        // 029107/18 cluster 143: 269 cm, contrast 2.04 = the tabulated 1.95,
        // ratio_mu 1.40, gate 0).  The Bragg-contrast and KS bits already
        // carry the shape verdict.
        cfg["pid_mode"] = m_pid_mode;
        // doc pdhd/03 sec 6: plateau_med / mip_dqdx must lie in [lo, hi] or
        // the candidate gets R_PLATEAU_OFF_MIP -- a track whose plateau reads
        // 0.3 MIP has no trustworthy charge (PDHD: 8 of 15 passers at 0.27-0.45
        // once the electron test was dropped).  hi <= lo disables (default).
        cfg["plateau_mip_lo"] = m_plateau_mip_lo;
        cfg["plateau_mip_hi"] = m_plateau_mip_hi;
        // doc pdhd/03 sec 6: when the stop vertex has a collinear MIP arm
        // (a continuation), extend the muon chain along it and re-judge at
        // the new end, up to this many times.  0 = doc pdvd/48 (the tagger's
        // stop is final and the continuation only rejects).
        cfg["stop_extend_max"] = m_stop_extend_max;
        // doc pdvd/57: the mirror of stop_extend_max.  find_first_kink's
        // charge gate wants BOTH arms of a kink >= 0.6 MIP, so an asymmetric
        // muon->Michel junction returns the sentinel and the stop is clamped
        // to the fit's far end -- the Michel's own tip on 51 of 125 missed
        // stoppers (doc pdvd/56 sec 2).  stop_extend_max only walks the stop
        // OUTWARD; this walks it back onto an existing chain vertex while the
        // dropped tail is collapsed and a Bragg rise survives before it (see
        // StmMichelFunctions.h stm_michel_stop_retreat).  0 = off (doc
        // pdvd/48/56 behaviour); measured on the d53v census at 2 to recover
        // 8 of 51 missed collapse-shaped stoppers with 1 new false positive
        // at this frac (doc pdvd/57); the default frac below is the point
        // that measured 0 new false positives (7 recovered).
        cfg["stop_retreat_max"] = m_stop_retreat_max;
        // doc pdvd/57: the dropped tail's median dQ/dx must be below this
        // fraction of the plateau (a muon plateau, even at the top of PDVD's
        // [0.6, 1.6] MIP admission window, cannot fall this low without
        // actually collapsing).
        cfg["retreat_collapse_frac"] = m_retreat_collapse_frac;
        // doc pdvd/57: the profile surviving the drop must still peak at this
        // fraction of the plateau within retreat_peak_window_cm of the new
        // end -- there must be a Bragg rise to retreat TO, not just charge to
        // drop (same ratio and window stm_michel's own shape census uses).
        cfg["retreat_peak_frac"] = m_retreat_peak_frac;
        cfg["retreat_peak_window_cm"] = m_retreat_peak_window_cm;
        // doc pdvd/58 (T1c): stop_retreat_max can only move the stop onto a
        // vertex the chain ALREADY has.  On a population of missed
        // collapse-shaped stoppers the collapse sits INSIDE the fit's last
        // chain segment -- no chain vertex anywhere near it (doc pdvd/57 sec
        // 1's 23-item finding) -- so no graph-local retreat reaches them.
        // This walks the SAME collapsed-tail/Bragg-rise test as the retreat,
        // but at a FIT ROW (PR::break_segment splits the segment there,
        // exactly as anchor_vertex() already does at the tagger's own
        // entry/stop), gated additionally on the fitted trajectory's own bend
        // at that row (stm_michel_row_kink_deg) -- losing the vertex
        // constraint removes the guard that kept the retreat honest, and the
        // kink is what replaces it (the owner's kink discriminator: charge
        // shape alone cannot tell a second particle).  0 = off (doc pdvd/56/57
        // behaviour).  Skipped whenever the retreat already fired.
        cfg["stop_split_max"] = m_stop_split_max;
        // doc pdvd/58: the fitted-trajectory bend (deg) a candidate split row
        // must reach.  Measured (over the missed-collapse population with no
        // chain vertex to retreat onto vs. the collapse-shaped through-going
        // control): median bend 18.5 deg (n=23) vs 6.3 deg (n=58); 15 deg
        // costs 0 of 58 new false positives while still reaching 3 of 23.
        // Deliberately its OWN knob, not stop_retreat_max's: coupling T1c's
        // operating point to the retreat's would silently move it if either
        // is retuned later.
        cfg["split_kink_min_deg"] = m_split_kink_min_deg;
        // doc pdvd/58: minimum residual range of the split row from the
        // fit's current far end -- a candidate closer than this trims almost
        // nothing (the offline probe's "first satisfying row" rule degenerated
        // to 1-2 cm trims before this floor was added).
        cfg["split_min_drop_cm"] = m_split_min_drop_cm;
        // doc pdvd/58: same meaning as retreat_collapse_frac, judged against
        // the SAME plateau reference (profile_plateau, shared with the
        // retreat) -- kept as a separate key for the reason given above.
        cfg["split_collapse_frac"] = m_split_collapse_frac;
        cfg["split_peak_frac"] = m_split_peak_frac;
        cfg["split_peak_window_cm"] = m_split_peak_window_cm;
        // doc pdvd/58: the +/- arm (cm of arclength) stm_michel_row_kink_deg
        // measures the incoming/outgoing direction over.
        cfg["split_dir_window_cm"] = m_split_dir_window_cm;
        // doc pdvd/61 (T2c): default off.  When on, an attached Michel
        // (michel_conn_type == 1) whose stop was moved this event (n_retreat
        // > 0 or n_split > 0) is demoted (michel_conn_type reset to 0, so
        // michel_found follows) when its assembled michel_ke_best falls below
        // moved_stop_michel_ke_min -- touches michel_found only, never
        // reject_bits/is_stm.  Measured against the census: 5 named
        // through-going items acquired a spurious attached Michel exactly
        // this way across T1a/T1b/T1c (docs 57/58/59), median KE 7.97 MeV,
        // against a 21.4 MeV median for the same conn_type==1 population's
        // true Michels.
        cfg["moved_stop_michel_guard"] = m_moved_stop_michel_guard;
        // doc pdvd/61: MeV, doc 55 sec 15.1's own KE floor re-used here for a
        // narrower, moved-stop-only population (that section's own cut
        // targets conn_type==2, dis_cm>3; this targets conn_type==1, dis_cm==0
        // -- the two do not overlap).
        cfg["moved_stop_michel_ke_min"] = m_moved_stop_michel_ke_min;
        // doc pdvd/62 (T3a): cm, 0 = off.  When > 0, the PR partition's pr54
        // isolated-residual test (NeutrinoOtherSegments.cxx) keeps a residual
        // whose fitted endpoint lies within this distance of the STM tagger's
        // stop, whatever the terminal-count / length floors said -- doc
        // pdvd/54 sec 2: 363 residuals per 120 events are dropped in the
        // 2-24 terminal band, which is a Michel fragment's size, 74 of them
        // within 20 cm of a scan candidate's stop.  Applies to the main
        // cluster's partition AND the admitted companions' (same pa).
        cfg["stop_local_residual_cm"] = m_stop_local_residual_cm;
        // doc pdvd/87 (doc pdvd/78 action item 7): 0 / 0 = no floor, doc 62's
        // keep unchanged.  A residual inside stop_local_residual_cm is kept
        // only with at least this many Steiner terminals AND at least this
        // fitted length (cm).  Doc pdvd/86 sec 8.1 on PDVD production: at 5 cm
        // with 5 / 5 the keep reaches 5 judged items (the owner's two unfitted
        // Michels 039253_0/44 and 039349_30/45, 0 through-going) where the
        // unfloored 5 cm keep reaches 11 (3 through-going 1-7 cm stubs) and
        // doc 62's 20 cm keep 31 (14 through-going).  Inert when
        // stop_local_residual_cm is 0.
        cfg["stop_local_residual_min_points"] = m_stop_local_residual_min_points;
        cfg["stop_local_residual_min_len_cm"] = m_stop_local_residual_min_len_cm;
        // doc pdvd/62 (T3b): default off.  A kept residual is a DISCONNECTED
        // piece of the main cluster's graph, which neither the attached path
        // (needs an arm at stop_v) nor the dots loop (companion clusters
        // only) ever looks at; when on, such a piece within michel_dot_radius_cm
        // of the stop, no longer than dot_max_len_cm and passing the same
        // body test as a companion piece joins the Michel object exactly as a
        // companion piece does (conn_type 2).  Attached interior arms are
        // excluded by construction.
        cfg["stop_local_michel_pieces"] = m_stop_local_michel_pieces;
        // doc pdvd/62 (T3c): default off.  doc 55 sec 15.1's kinematic test as
        // a knob: a bridged (conn_type 2) or charge-only (3) Michel farther
        // than michel_range_energy_dis_cm from the stop carrying less than
        // michel_range_energy_ke_min MeV cannot be an electron born at the
        // stop (26 of 45 conn_type-2 objects on the census fail it; 0 of 113
        // attached ones).  Demotes via michel_conn_type only -- michel_found
        // follows, reject_bits / is_stm untouched.  The 5 cm default is the
        // argued distance ("cannot travel 5 cm ... under 10 MeV"); doc 55's
        // best-F1 row was 3 cm, chosen on the same record, so it is offered
        // as the owner's option rather than shipped.
        cfg["michel_range_energy_guard"] = m_michel_range_energy_guard;
        cfg["michel_range_energy_dis_cm"] = m_michel_range_energy_dis_cm;
        cfg["michel_range_energy_ke_min"] = m_michel_range_energy_ke_min;
        // doc pdvd/64 (T6): default off.  When on, every arm the classifier
        // returned kOther for -- at an interior chain vertex (> delta_max_len_cm
        // and cooler than vertex_hadron_mip, or a short arm whose far vertex
        // continues) or at the stop (a hot stub not absorbed, debris beside a
        // Michel, a demoted continuation) -- gets role-7 rows in stm_michel_pts,
        // written LAST and claiming nothing, so no verdict row changes and an
        // arm the Michel object walk absorbed keeps its role 3.  Rows only:
        // stm_michel_pts grows, every verdict branch is unchanged.  The
        // interior count (n_body_other) has always been persisted; the stop
        // count (n_stop_other) and the published count are new.
        cfg["publish_other_arms"] = m_publish_other_arms;
        // doc pdvd/80 (doc 78 action item 1): THE SEGMENT CENSUS.  Every PR
        // segment of the MAIN cluster that no stage claimed and no stage
        // published gets role-8 rows, plus the same rej / d_stop / d_body
        // columns the survey (role 6) carries: the doc 62 T3b piece gates read
        // on it (2 too far from the stop, 3 longer than the piece cap, 4 body
        // exclusion), or 13 attached to the chain and not taken, 14 passes
        // every T3b gate yet unclaimed, 15 no endpoint vertices, 9 T3b off,
        // 10 no stop vertex.  Rows only, written LAST; no verdict reads them.
        // On the production record 42 of the scanner's michel tags sit on such
        // segments (doc 78 sec 3.2) and nothing in the output names them.
        // Default false: the knob-off T_stm_michel_pts schema and every branch
        // are byte-identical.
        cfg["segment_census"] = m_segment_census;
        // doc pdvd/65 (T7): default off.  When on, the verdict-stage contrast
        // and KS tests read the live profile from rr' = end_L - L with end_L =
        // L[peak] + 0.2 cm, the peak being the row maximising the 5-point
        // running mean of dQ/dx within bragg_peak_search_cm of the geometric
        // end (TaggerCheckSTM::eval_stm_core_impl's recipe, ToyFiducial.cxx:
        // 1551); rows past the peak are dropped.  The chain walk's own Bragg
        // guards keep the geometric origin.  bragg_anchor_shift_cm records the
        // move (0 when off, or when the peak is the last row).
        cfg["bragg_peak_anchor"] = m_bragg_peak_anchor;
        cfg["bragg_peak_search_cm"] = m_bragg_peak_search_cm;
        // doc pdvd/70 (P1): default off.  The owner's rule: a Michel at the end
        // is by itself strong evidence of a stop; the dQ/dx rise counts only
        // when it genuinely matches a Bragg peak.  When on, a candidate whose
        // Michel object exists (michel_found, after the T2c / T3c vetoes),
        // attached (conn_type 1) or bridged (2), with michel_ke_best >=
        // topology_michel_ke_min MeV and michel_len >=
        // topology_michel_len_min_cm, has R_NO_BRAGG and R_SHAPE_FLAT cleared
        // (and R_PROFILE_SPARSE with topology_clears_sparse) just before the
        // verdict is persisted -- after every other bit is set, so is_stm moves
        // only 0 -> 1 and only when nothing else rejects.  michel_found is read,
        // never written.  10 MeV = the T2c / T3c floor already in the bag, 3 cm
        // = the range-energy distance.  topology_cleared_bits (the bits this
        // cleared) is written only when the knob is on, so the tree is
        // byte-identical when off.
        cfg["topology_stop_evidence"] = m_topology_stop_evidence;
        cfg["topology_michel_ke_min"] = m_topology_michel_ke_min;
        cfg["topology_michel_len_min_cm"] = m_topology_michel_len_min_cm;
        cfg["topology_clears_sparse"] = m_topology_clears_sparse;
        // doc pdvd/71 (P4): default off.  The owner's Q4 (doc pdvd/70 sec 5):
        // isolated gamma blobs near the stop -- the Michel electron's brems and
        // Compton pieces -- belong to the Michel object, and through doc 70
        // ~150 of the ~165 the owner tagged sat outside it.  When on, after the
        // capture-gamma stage, every UNCLAIMED companion cluster with a fitted
        // segment is offered to stm_michel_gamma_gate (StmMichelFunctions.h):
        // within michel_gamma_radius_cm of the FINAL stop, no longer than
        // michel_gamma_max_len_cm, inside a cone of cos >= michel_gamma_cos_min
        // about the stop -> Michel-object direction, closer to the Michel than
        // to the muon body, and at most michel_gamma_max_ke_mev.  Survivors are
        // taken nearest first while michel_ke_best + michel_ke_gamma stays
        // <= michel_gamma_total_ke_max_mev (over-clustering must not hand the
        // Michel a huge energy), once michel_found is final.  They get role-4
        // rows and nothing else: no PDG, no Shower, no PF node, and
        // michel_ke_best / michel_found / is_stm are not read back, so every
        // pre-existing branch is unchanged.  35 cm = the capture-gamma radius,
        // i.e. today's admission radius, so switching the knob on admits no new
        // companion (a larger radius widens admission, and doc pdvd/53 measured
        // what that does to the fit through preload_clusters).  The six new
        // branches are written only when the knob is on.
        cfg["michel_gamma_collect"] = m_michel_gamma_collect;
        cfg["michel_gamma_radius_cm"] = m_michel_gamma_radius_cm;
        cfg["michel_gamma_max_len_cm"] = m_michel_gamma_max_len_cm;
        cfg["michel_gamma_cos_min"] = m_michel_gamma_cos_min;
        cfg["michel_gamma_max_ke_mev"] = m_michel_gamma_max_ke_mev;
        cfg["michel_gamma_total_ke_max_mev"] = m_michel_gamma_total_ke_max_mev;
        // doc pdvd/81: the charge-based Michel energy.  The owner's estimator:
        // a Michel trajectory wiggles, so neither its range nor the fitted
        // dQ/dx integral (michel_ke_best) is trusted; instead every 2-D cell the
        // chain's own association rule (kine_charge_from_maps: nearest
        // associated point within michel_q2d_dis_cm, per plane) assigns to the
        // Michel object is summed as measured charge MINUS the charge the muon
        // chain's own multi-track fit predicts on that cell (the fit's stored
        // response, TrackFitting keep_dqdx_response; the muon rows are the
        // chain's Fit::index set, the shared stop-vertex row counted as muon),
        // the three planes are combined with the chain's rule
        // (stm_michel_combine_planes: 0.25/0.25/1.0, 4 % asymmetry switch) and
        // the charge is converted with the unfitted-piece constant
        // (michel_unfit_from_model / michel_unfit_dedx).  The taken role-4 gamma
        // blobs get the same treatment as a separate term.  Nothing reads the
        // result back: michel_ke_best, michel_found, is_stm and every reject bit
        // are unchanged, and the branches are written only when on.
        cfg["michel_q2d"] = m_michel_q2d;
        // doc pdvd/81: the per-cell table behind it (PC stm_michel_2d ->
        // T_stm_michel_2d): every Michel (role 3) and gamma (role 4) cell with
        // its measured charge, error, flag, the muon-only and the full fit
        // prediction, plus the STM's own footprint (role 1: cells the fit
        // predicts from chain rows within michel_q2d_stm_window_cm of the stop;
        // -1 = the whole chain).  Read only when michel_q2d is on.
        cfg["michel_q2d_cells"] = m_michel_q2d_cells;
        cfg["michel_q2d_dis_cm"] = m_michel_q2d_dis_cm;              // cm, = kine_charge_from_maps's 0.6 cm literal
        cfg["michel_q2d_stm_window_cm"] = m_michel_q2d_stm_window_cm;
        // doc pdvd/72 (P3b): deg, -1 = off.  When >= 0, the moved-stop veto
        // (moved_stop_michel_guard) does not demote an attached Michel whose
        // kink at the stop (michel_kink_deg) is at least this -- the turn is
        // the Michel's own evidence, which the KE floor alone does not read.
        // n_michel_veto_exempt counts the spared ones and is written only
        // when the knob is on.  is_stm is untouched: T2c fires below
        // moved_stop_michel_ke_min and topology_stop_evidence needs at least
        // topology_michel_ke_min, both 10 MeV, so a spared Michel cannot
        // reach P1 -- true only while those two stay equal.
        cfg["moved_stop_michel_kink_min"] = m_moved_stop_michel_kink_min;
        // doc pdvd/84 (doc 78 action item 3): cm, -1 = off.  When >= 0, the
        // moved-stop veto also spares an attached Michel the kink test did
        // not, if its reach -- michel_len + michel_far_len, the arm plus the
        // subtree past its far end -- is at least this.  On the PDVD record
        // the veto's instances (12 on 7 items, 28 arms) put every
        // through-going arm at <= 5.8 cm and every owner-confirmed Michel at
        // >= 7.5 cm, while the kink leaves a 0.9 deg window and the KE floor
        // does not separate at all.  n_michel_veto_reach_exempt counts the
        // spared ones and is written only when the knob is on.  The same P1
        // argument as the kink test holds: a spared Michel is under
        // moved_stop_michel_ke_min, so topology_stop_evidence cannot fire on
        // it while that and topology_michel_ke_min stay equal.  The far
        // subtree is the stop-arm classifier's walk, which can re-enter the
        // muon chain through a side loop (doc pdvd/83 sec 9.4) and read the
        // muon as the arm's reach; no T2c instance on the record does.
        cfg["moved_stop_michel_reach_min_cm"] = m_moved_stop_michel_reach_min_cm;
        // doc pdvd/73 (P2): three operating points for the attached Michel
        // gate (stm_michel_michel_gate), each off at -1; with all three off
        // the classifier runs its doc pdvd/48 expression verbatim.
        //   michel_mip_lo_turned (+ _kink_deg 60): the charge floor for an
        //     arm that turns at least _kink_deg (a diluted Michel admitted by
        //     its topology, never by low dQ/dx alone);
        //   michel_far_len_shower_max_cm: a shower-flagged arm's far subtree
        //     (its own brems) is capped on its own instead of counting in
        //     len + far_len <= michel_max_len_cm;
        //   michel_kink_window_cm: the Michel turn test also accepts the
        //     kink over this (shorter) window -- the continuation test and
        //     michel_kink_deg keep the classifier's window.
        cfg["michel_mip_lo_turned"] = m_michel_mip_lo_turned;
        cfg["michel_mip_lo_turned_kink_deg"] = m_michel_mip_lo_turned_kink_deg;
        cfg["michel_far_len_shower_max_cm"] = m_michel_far_len_shower_max_cm;
        cfg["michel_kink_window_cm"] = m_michel_kink_window_cm;
        // doc pdvd/74 (P3): the Michel the fit carried inside the muon chain.
        //   retreat_tail_strict: stop_retreat_max's collapse test reads only
        //     the rows strictly past the vertex it would retreat onto.  That
        //     vertex's own row is the kept segment's end -- on an overshoot it
        //     is the Bragg peak, and over a 2-3 cm segment it sets the tail
        //     median by itself (039252_16/32).
        //   retreat_tail_sublive: the same test also counts rows below the
        //     profile_min_dqdx_frac floor -- a collapsed overshoot reads
        //     0.1-0.2 MIP, which the floor calls dead (039253_3/61).  It
        //     cannot tell a dead-channel stretch from a collapse.
        //   michel_collinear_split: doc 70's P3 as proposed -- the retreat and
        //     the split also run on a chain whose profile already shows the
        //     Bragg rise.
        // stop_move_p3_bits (bit0: the tail reading changed the retreat's
        // answer; bit1: the stop moved on a Bragg-confirmed chain) is written
        // only when one of the three is on.  A P3 retreat is a retreat: the
        // moved-stop veto (T2c) reads it like any other.
        cfg["retreat_tail_strict"] = m_retreat_tail_strict;
        cfg["retreat_tail_sublive"] = m_retreat_tail_sublive;
        cfg["michel_collinear_split"] = m_michel_collinear_split;
        // doc pdvd/82 (doc 78 action item 2): an ALTERNATIVE admission for the
        // collapsed-tail test both stop movers apply.  Today a tail counts as
        // collapsed only below retreat_collapse_frac/split_collapse_frac x the
        // PLATEAU.  On the 14 items doc 78 sec 2.3 named, the fit runs through
        // the Michel: after a Bragg peak of 1.5-2.9 x plateau the tail reads
        // 0.54-1.87 x plateau -- far under a post-Bragg muon, never under the
        // track's own plateau.  stop_tail_peak_frac admits tail_med <= frac x
        // the surviving PEAK instead.  That reading is much looser (with peak /
        // plateau in 1.4-3, 0.5 x peak is 0.7-1.5 x plateau), so it is allowed
        // only when the row's own trajectory bend reaches
        // stop_tail_peak_kink_min_deg -- the doc 58 discriminator, which is
        // what separates a second particle from a muon that keeps going.  A row
        // the plateau test already accepts never sees the bend requirement, so
        // 0 (off) leaves the doc 57/58/74 path byte-identical.
        // Both movers still refuse to run at all on a Bragg-confirmed chain
        // unless michel_collinear_split (doc 74) is also on -- which is the
        // real reason they are silent on 14 of doc 78's 21 items.
        cfg["stop_tail_peak_frac"] = m_stop_tail_peak_frac;
        cfg["stop_tail_peak_kink_min_deg"] = m_stop_tail_peak_kink_min_deg;
        // doc pdvd/83 (doc 78 action item 4): default 0 (off).  The four
        // Michels doc 78 sec 3.2 found missing on michel_found-0 stoppers are
        // PR segments ATTACHED to the muon chain (doc 80: rej 13), each hanging
        // off the chain's penultimate vertex 2.5-6.5 cm before the stop, the
        // last chain segment a short stub.  The stop-arm classifier never sees
        // them and the interior-vertex one calls them kOther, which nothing
        // reads.  When > 0 and NO Michel was found (michel_conn_type still 0
        // after the attached and companion stages), every kOther body arm at a
        // chain vertex within this along-chain distance of the stop is offered
        // the stop-arm gate (stm_michel_classify_stop_arm, kink against the
        // incoming chain segment); the nearest vertex with a kMichel becomes
        // an ATTACHED Michel (conn 1) started at that vertex.  Offline on the
        // owner's scan: 2 of the 4 within 5-7 cm, 0 control items; at 10 cm
        // the first STM_ONLY and THRU arms enter -- the distance is what holds
        // the purity.  The block writes the Michel object and nothing else a
        // verdict reads; topology_stop_evidence is the one declared channel
        // through which is_stm can then move.
        cfg["michel_near_stop_arm_cm"] = m_michel_near_stop_arm_cm;
        // doc pdvd/75 (P1b): default off.  The 3 cm peak anchor (T7) discards
        // the rows past the running-mean maximum; on a rise that runs to the
        // fit's last row the low partial-step end row pulls that maximum 2-3
        // rows back and the anchor throws away the top of the Bragg rise, so
        // the anchored profile reads flat.  When on: if the anchor fired, the
        // anchored reading set any of no_bragg / shape_flat / plateau_off_mip
        // / profile_sparse, and the anchored peak (the winning 5-point mean)
        // is at least bragg_anchor_rise_min x the anchored plateau median,
        // the contrast, plateau and KS tests are re-read at the geometric
        // origin; if none of the four bits is set there, that reading (bragg,
        // ks_*, ratio_*) replaces the anchored one and the four bits are
        // cleared.  It can only clear bits.  bragg_anchor_fallback (1 when it
        // fired) is written only when on.  The do_track_comp and dead-volume
        // probes keep the anchored profile.
        cfg["bragg_anchor_geo_fallback"] = m_bragg_anchor_geo_fallback;
        cfg["bragg_anchor_rise_min"] = m_bragg_anchor_rise_min;
        // doc pdvd/66 (T8): default off.  The fields end_arc_span (fit arc over 3-D
        // span in the last end_window_cm of the chain profile -- doc 55 sec 17.1
        // item 5, the coiled end) and n_unsupported_segs (fitted segments of the
        // main cluster at least unsupported_min_len_cm long whose median dQ/dx is
        // below unsupported_frac of the chain's plateau median -- item 4, the fit
        // that spans ground the imaging never covered) are ALWAYS written; with
        // the guard on, either condition sets R_PROFILE_GEOMETRY ("the profile is
        // not a measurement"), a hard veto like every other bit.  The literals
        // are doc 55 sec 17.1's own class definitions (1.5, 20 cm, 0.25).
        cfg["profile_geometry_guard"] = m_profile_geometry_guard;
        cfg["profile_arc_span_max"] = m_profile_arc_span_max;
        cfg["unsupported_min_len_cm"] = m_unsupported_min_len_cm;
        cfg["unsupported_frac"] = m_unsupported_frac;
        cfg["end_window_cm"] = m_end_window_cm;
        // doc pdhd/03 sec 6: FiducialUtils::check_dead_volume from the end of
        // the live profile along the muon direction; a stop that walks into a
        // dead region is R_STOP_INTO_DEAD (three PDHD tracks end on the same
        // y = 493 cm line where a dead block begins).  Needs the fiducialutils
        // stage; false = off.
        cfg["dead_volume_check"] = m_dead_volume_check;
        // doc pdhd/03 sec 6: the fraction of the main cluster's own 3-D points
        // within coverage_radius_cm of a reconstructed point (muon chain,
        // deltas, Michel, dots).  Below min_chain_coverage the cluster is not a
        // track with attachments but a shower / blob the fit threaded a path
        // through (029107/11 cluster 18: a 65 cm chain inside a 1591-point EM
        // blob, "contrast 3.2" on the shower core).  0 = off.
        cfg["min_chain_coverage"] = m_min_chain_coverage;
        // doc pdhd/03 sec 6.3: at a stop that already shows the Bragg rise
        // (live contrast >= bragg_contrast_min x expected) AND has a
        // Michel-class arm, a collinear MIP arm no longer than michel_max_len
        // is stop debris (a delta, or the electron's other branch), not a
        // continuation: it neither extends the chain nor sets R_CONTINUATION.
        // 029107/12 cluster 112: an 11 cm arm at 10 deg beside a 108-deg,
        // 12.7 MeV Michel -- extending along it walked the stop to the anode
        // face and lost the Michel.  false = doc pdvd/48.
        cfg["michel_guards_stop"] = m_michel_guards_stop;
        // doc pdhd/03 sec 6.8: a shower-flagged stop arm with a MEASURABLE kink
        // below this is not a Michel (the muon's own Bragg stub, mis-split);
        // -1 = doc pdvd/48 (the flag alone admits).
        cfg["michel_shower_min_kink_deg"] = m_michel_shower_min_kink_deg;
        // doc pdhd/03 sec 6.8: with stop_extend_max > 0, also absorb a short
        // collinear arm HOTTER than a MIP (the muon's own Bragg stub the
        // partition split off) into the chain.  Measured on PDHD 029107/1
        // cluster 113: absorbing its 5.5 cm 1.64-MIP stub moved the tail window
        // onto the fading tip and turned a clean STM (contrast 1.76) into
        // no_bragg -- so OFF by default and OFF in the PDHD bag; the stub is
        // still not a Michel (michel_shower_min_kink_deg).
        cfg["absorb_bragg_stub"] = m_absorb_bragg_stub;
        // doc pdhd/03 sec 6.9: containment of the stop with the SAME per-wall
        // margins the cosmic taggers get (config fv_tolerance, internal units,
        // negative = inset) instead of the flat stop_fv_margin_cm inset; a
        // 746 cm through-going muon leaving PDHD at y = 597 cm passed the flat
        // 5 cm inset where TGM's 17.5 cm y-margin would have caught it.
        // false = doc pdvd/48.
        cfg["stop_fv_use_config_tolerance"] = m_stop_fv_use_config_tolerance;
        cfg["fv_tolerance"] = Json::Value(Json::arrayValue);
        cfg["coverage_radius_cm"] = m_coverage_radius_cm;
        cfg["fit_blob_coverage"] = m_fit_blob_coverage;
        cfg["dqdx_fit_keep_all_points"] = m_dqdx_fit_keep_all_points;
        cfg["excl_t0_frame"] = m_excl_t0_frame;
        // doc pdhd/16.  A stopping muon's range energy is the baseline and MCS
        // is the independent cross-check -- the one estimator that reads no
        // charge at all, so it is blind to gain, lifetime and recombination.
        // The engine (mcs/, WireCell::Mcs::MuonMCS) is already a clus
        // dependency; PR::mcs_fill_kine is NOT reused because it writes a
        // KineInfo this component does not have and its beam_window_only
        // default returns silently on cosmics.
        // mcs_cathode_xcut: half-width (cm) of the band around mcs_cathode_x
        // whose segments are dropped; both ProtoDUNEs are cathode-centred at
        // x = 0 and lose charge at the seam.
        cfg["mcs_enable"] = m_mcs_enable;
        cfg["mcs_min_len_cm"] = m_mcs_min_len_cm;
        cfg["mcs_cathode_x"] = m_mcs_cathode_x;
        cfg["mcs_cathode_xcut"] = m_mcs_cathode_xcut;
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
    double m_michel_dot_radius_cm{15.0}, m_dot_max_len_cm{25.0}, m_dot_body_exclusion_cm{5.0};
    double m_companion_max_len_cm{25.0};                                   // doc pdhd/15
    double m_michel_unfit_recom{0.7}, m_michel_unfit_fudge{0.95}, m_michel_unfit_w_ev{23.6};
    // doc pdhd/17.  false = the doc pdhd/15 pair, so an absent key is
    // byte-identical.  2.1 MeV/cm is the MIP operating point the
    // PowerBoxRecombination fit also pivots on (RecombinationModels.h).
    bool m_michel_unfit_from_model{false};
    double m_michel_unfit_dedx{2.1};
    // doc pdvd/51.  A mu- that stops in argon is CAPTURED far more often than
    // it decays, and the capture leaves de-excitation gammas.  A gamma is
    // neutral: it travels, then deposits a compact blob off the muon's axis --
    // 039252_0 cluster 77's is 6 points, 2.17 cm, 0.8 MeV, 26.5 cm from the
    // stop at 73.5 deg to the muon direction.  That is a DIFFERENT object from a
    // Michel and it must not be reached by widening michel_dot_radius_cm, which
    // also feeds the Michel piece assembly: 039252_15 cluster 77 is already the
    // only PDVD object above the 52.8 MeV Michel endpoint (76.8 of 160), and
    // handing it more distant charge makes exactly the item the owner is
    // scanning worse.  So: its own ring, its own compactness and energy caps,
    // its own branches, its own PF nodes.
    bool m_stop_gamma_enable{true};
    // 35 cm: the same-bundle stop-anchored blob density knee (1.66 -> 0.60 per
    // candidate per 1e6 cm^3 between the 20-30 and 30-40 shells, against an
    // ambient floor the foreign-bundle control measures FLAT at ~0.6), and the
    // plateau of the mu- anti-correlation, which peaks at 2.05 here and
    // collapses to 1.41 by 45 cm.  039252_0 cluster 77's gamma is at 26.5 cm.
    double m_stop_gamma_radius_cm{35.0};      // ring outer edge; inner edge is michel_dot_radius_cm
    double m_stop_gamma_max_len_cm{10.0};     // a gamma deposit is a blob, not a track
    double m_stop_gamma_min_ke_mev{0.2};
    double m_stop_gamma_max_ke_mev{20.0};
    int m_stop_gamma_max_n{8};
    // doc pdvd/85 (doc 78 action item 5): a capture gamma is the claim that
    // the muon STOPPED and was captured, so on a candidate the verdict rejects
    // it has no stop to belong to.  On the owner's merged scan record the
    // stage's clusters are 36 gamma / 0 delta-other on is_stm 1 candidates and
    // 7 gamma / 33 delta-other on is_stm 0 ones (27 of those on through-going
    // muons).  On: after the verdict is final, a rejected candidate's capture
    // gammas are withheld -- no role-5 rows, no PF showers, the segments'
    // particle info restored, stop_gamma_* at their defaults, and the count in
    // n_stop_gammas_withheld.  Off (default): published on every candidate.
    bool m_stop_gamma_require_stm{false};
    // doc pdvd/53: the SURVEY.  Scaffolding for the hand scan, not a physics
    // selection.  A companion admitted only by the survey is fitted into the
    // graph and given role-6 point rows so the display can draw it, click it
    // and let the scanner group it; it can become neither a Michel piece (the
    // cluster re-test at michel_dot_radius_cm is untouched) nor a capture gamma
    // (the ring's outer edge is stop_gamma_radius_cm).  The ONLY physics delta
    // is the preload perturbation of the candidate's own fit.  Default OFF, so
    // an absent key leaves every branch AND every stm_michel_pts row
    // byte-identical; the two ProtoDUNE pr.jsonnet turn it on.
    bool m_survey_enable{false};
    double m_survey_radius_cm{60.0};      // admission ring, outer edge, from the stop
    double m_survey_max_len_cm{25.0};     // mirrors companion_max_len_cm
    double m_delta_max_len_cm{8.0};
    double m_vertex_hadron_mip{1.4};
    double m_profile_min_dqdx_frac{0.0};
    int m_pid_mode{0};
    double m_plateau_mip_lo{0.0}, m_plateau_mip_hi{0.0};
    int m_stop_extend_max{0};
    int m_stop_retreat_max{0};                    // doc pdvd/57
    double m_retreat_collapse_frac{0.5};
    double m_retreat_peak_frac{1.4};
    double m_retreat_peak_window_cm{15.0};
    int m_stop_split_max{0};                      // doc pdvd/58 (T1c)
    double m_split_kink_min_deg{15.0};
    double m_split_min_drop_cm{3.0};
    double m_split_collapse_frac{0.5};
    double m_split_peak_frac{1.4};
    double m_split_peak_window_cm{15.0};
    double m_split_dir_window_cm{5.0};
    bool m_moved_stop_michel_guard{false};        // doc pdvd/61 (T2c)
    double m_moved_stop_michel_ke_min{10.0};      // MeV
    double m_stop_local_residual_cm{0.0};         // doc pdvd/62 (T3a): 0 = off
    int m_stop_local_residual_min_points{0};      // doc pdvd/87: Steiner terminals, 0 = no floor
    double m_stop_local_residual_min_len_cm{0.0}; // doc pdvd/87: cm, 0 = no floor
    bool m_stop_snap_reachable{false};            // doc pdvd/88: reachable-only stop snap on keep-fire candidates
    bool m_stop_local_michel_pieces{false};       // doc pdvd/62 (T3b)
    bool m_michel_range_energy_guard{false};      // doc pdvd/62 (T3c)
    double m_michel_range_energy_dis_cm{5.0};     // cm
    double m_michel_range_energy_ke_min{10.0};    // MeV
    bool m_publish_other_arms{false};             // doc pdvd/64 (T6): role-7 rows for kOther arms
    bool m_segment_census{false};                 // doc pdvd/80: role-8 rows for every unclaimed PR segment of the main cluster
    bool m_bragg_peak_anchor{false};              // doc pdvd/65 (T7): peak-anchored rr origin for the verdict shape tests
    double m_bragg_peak_search_cm{10.0};          // cm
    bool m_topology_stop_evidence{false};         // doc pdvd/70 (P1)
    double m_topology_michel_ke_min{10.0};        // MeV
    double m_topology_michel_len_min_cm{3.0};     // cm
    bool m_topology_clears_sparse{false};         // doc pdvd/70 sec 9.4: also clear R_PROFILE_SPARSE
    bool m_michel_gamma_collect{false};           // doc pdvd/71 (P4)
    double m_michel_gamma_radius_cm{35.0};        // from the FINAL stop; = the admission radius today
    double m_michel_gamma_max_len_cm{10.0};       // a dot, not a track
    double m_michel_gamma_cos_min{0.5};           // 60 deg about the stop -> Michel direction
    double m_michel_gamma_max_ke_mev{20.0};       // per blob
    double m_michel_gamma_total_ke_max_mev{60.0}; // the Michel object with its blobs: the 52.8 MeV endpoint plus resolution
    bool m_michel_q2d{false};                     // doc pdvd/81: the charge-based Michel energy
    bool m_michel_q2d_cells{false};               // doc pdvd/81: the Michel / STM 2-D cell table (stm_michel_2d)
    double m_michel_q2d_dis_cm{0.6};              // cm, the chain's association radius (kine_charge_from_maps)
    double m_michel_q2d_stm_window_cm{30.0};      // cm from the stop for the role-1 STM footprint rows; -1 = whole chain
    double m_moved_stop_michel_kink_min{-1.0};    // doc pdvd/72 (P3b): deg, -1 = off
    double m_moved_stop_michel_reach_min_cm{-1.0};// doc pdvd/84: cm, -1 = off
    double m_michel_mip_lo_turned{-1.0};          // doc pdvd/73 (P2a): -1 = off
    double m_michel_mip_lo_turned_kink_deg{60.0}; // deg
    double m_michel_far_len_shower_max_cm{-1.0};  // doc pdvd/73 (P2b): cm, -1 = off
    double m_michel_kink_window_cm{-1.0};         // doc pdvd/73 (P2c): cm, -1 = off
    bool m_retreat_tail_strict{false};            // doc pdvd/74 (P3)
    bool m_retreat_tail_sublive{false};           // doc pdvd/74 (P3)
    bool m_michel_collinear_split{false};         // doc pdvd/74 (P3, doc 70's literal)
    double m_stop_tail_peak_frac{0.0};            // doc pdvd/82: <= 0 = off
    double m_stop_tail_peak_kink_min_deg{25.0};   // doc pdvd/82: deg
    double m_michel_near_stop_arm_cm{0.0};        // doc pdvd/83: cm, <= 0 = off
    bool m_bragg_anchor_geo_fallback{false};      // doc pdvd/75 (P1b): the geometric reading may stand when the anchor rejects a prominent peak
    double m_bragg_anchor_rise_min{1.5};          // x the anchored plateau median
    bool m_profile_geometry_guard{false};         // doc pdvd/66 (T8)
    double m_profile_arc_span_max{1.5}, m_unsupported_min_len_cm{20.0}, m_unsupported_frac{0.25}, m_end_window_cm{20.0};
    bool m_dead_volume_check{false};
    double m_min_chain_coverage{0.0}, m_coverage_radius_cm{3.0};
    bool m_michel_guards_stop{false};
    double m_michel_shower_min_kink_deg{-1.0};
    bool m_absorb_bragg_stub{false};
    bool m_stop_fv_use_config_tolerance{false};
    std::vector<double> m_fv_tolerance;
    double m_fit_blob_coverage{-1.0};
    bool m_dqdx_fit_keep_all_points{false};
    bool m_excl_t0_frame{false};
    // doc pdhd/16: multiple-Coulomb-scattering momentum for the STM muon.
    // Default OFF so an absent bag leaves the compiled config and the tree
    // values exactly where doc pdhd/15 left them; both ProtoDUNE drivers set
    // mcs_enable true.
    bool m_mcs_enable{false};
    double m_mcs_min_len_cm{40.0};
    double m_mcs_cathode_x{0.0}, m_mcs_cathode_xcut{0.0};

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
        int n_live_pts{0}, n_dead_pts{0}, n_cmp_live{0};   // doc pdhd/03: live = dQ/dx >= profile_min_dqdx_frac * mip
        double dead_frac_cmp{0};                          // dead fraction within compare_range of the stop
        double muon_len{0};
        int n_delta{0}, n_body_other{0}, n_body_hadron{0};
        int n_stop_other{0}, n_other_published{0};   // doc pdvd/64 (T6): kOther arms at the stop; role-7 segments actually published
        int n_census_segs{0}, n_census_admissible{0};   // doc pdvd/80: role-8 segments written; of them, rej 14 (passes every T3b gate, unclaimed)
        double bragg_anchor_shift{0};                // doc pdvd/65 (T7): how far back of the geometric end the peak anchor put rr = 0 (0 = off / peak at the end)
        // doc pdvd/66 (T8): the profile's end geometry and the fitted segments' charge support -- reject-class inputs.  -1 = not computed.
        double end_arc_cm{-1}, end_span_cm{-1}, end_arc_span{-1}; int n_end_pts{0};
        int n_unsupported_segs{-1}; double unsupported_len_cm{-1}, unsupported_frac_min{-1}, chain_support_min{-1};
        double delta_len{0};
        // dQ/dx shape
        double ks_mu{0}, ks_flat{0}, ratio_mu{0}, ratio_flat{0};
        std::vector<double> comp_fwd{0, 1e9, 1e9, 1e9}, comp_bwd{0, 1e9, 1e9, 1e9};
        StmMichelBragg bragg;
        // stop arms
        int n_stop_arms{0}, n_michel_segs{0}, michel_found{0}, michel_conn_type{0};
        double michel_len{0}, michel_mip{0}, michel_kink_deg{-1}, michel_far_len{0};
        double michel_ke_dqdx{0}, michel_ke_range{0}, michel_ke_best{0};
        // doc pdhd/14
        double muon_ke_range{0}, muon_ke_dqdx{0}, muon_ke_best{0};
        // doc pdhd/16.  -1 means "not computed" and must survive as -1: the
        // engine returns -1 for every path it refuses (bad_path, < 20 trimmed
        // points, trimmed end < 2*14 cm from the stop, < 2 fitted segments),
        // and a 0 there would read as a measured zero energy.
        double muon_ke_mcs{-1}, muon_mcs_amb{-1}, muon_mcs_tracklen{-1}, muon_mcs_range_ke{-1};
        int muon_mcs_nsegs{0}, muon_mcs_bad_path{0};
        // doc pdhd/16 sec 9: how often the cathode-band excision actually
        // fired on THIS muon.  Both are 0 when mcs_cathode_xcut is 0, and
        // also when it is on but no 14 cm segment reaches the band -- in
        // which case the engine's angle mask stays empty and the result is
        // bit-identical to the excision being off (MuonMCS.cxx:1181-1210).
        // A cathode-caused abort needs no branch of its own: it is exactly
        // muon_ke_mcs < 0 && muon_mcs_nsegs >= 2 && muon_mcs_cathode_angles > 0.
        int muon_mcs_cathode_segs{0}, muon_mcs_cathode_angles{0};
        double muon_p_range{-1}, muon_p_dqdx{-1}, muon_p_mcs{-1};
        int michel_seg_id{-1};                     // the daughter, named the way stop_vtx_id names the shared vertex
        // doc pdhd/15: the Michel as ONE object -- core + pieces
        double michel_ke_core{0};                  // the core only: what michel_ke_dqdx meant before doc 15
        double michel_ke_charge{0};                // Shower::get_kine_charge(), for comparison -- never `best`
        int michel_n_pieces{0};                    // fitted member segments + unfitted companion clusters
        int michel_parent_vtx_id{-1};              // the muon vertex the object hangs from (= stop_vtx_id)
        double michel_dis_cm{-1};                  // stop -> object start; 0 when attached (doc pdvd/83: the along-chain distance when michel_near_arm)
        Point michel_start_pt;
        double cont_len{0}, cont_angle_deg{-1}, cont_mip{0};
        int n_ext{0}; double ext_len{0};          // doc pdhd/03: chain extensions past the tagger's stop
        int n_stub_absorb{0};                     // doc pdvd/63 (T5): of those, absorb_bragg_stub's (a hot collinear stub taken as the true end)
        int n_retreat{0}; double retreat_len{0};  // doc pdvd/57: chain segments retreated off the fit's far end
        int n_split{0}; double split_len{0}, split_kink_deg{0};  // doc pdvd/58: T1c fit-row split
        int stop_move_p3_bits{0};  // doc pdvd/74 (P3): bit0 the P3 tail reading changed the retreat's answer, bit1 the stop moved on a Bragg-confirmed chain; doc pdvd/82: bit2 the peak-relative tail is what admitted the move
        int bragg_anchor_fallback{0};  // doc pdvd/75 (P1b): 1 = the anchor's rejection was replaced by the geometric reading
        int michel_near_arm{0};        // doc pdvd/83: 1 = the Michel is an arm leaving the chain before the stop
        double near_arm_dist_cm{-1};   // doc pdvd/83: along-chain distance stop -> that arm's vertex; -1 = none
        int n_near_arms_examined{0};   // doc pdvd/83: kOther body arms near the stop offered the gate
        int n_michel_veto{0};  // doc pdvd/61: T2c fired -- an attached moved-stop Michel with too little charge was demoted
        int n_michel_veto_exempt{0};  // doc pdvd/72 (P3b): T2c would have fired, and the Michel's turn spared it
        int n_michel_veto_reach_exempt{0};  // doc pdvd/84: T2c would have fired, the turn did not spare it, its reach did
        int n_kept_near_stop_main{0}, n_kept_near_stop_comp{0};  // doc pdvd/62 (T3a): pr54 residuals kept by the stop anchor, main cluster / companions
        int n_floored_near_stop{0};  // doc pdvd/87: residuals inside the stop anchor's radius the size floor refused (main + companions)
        int stop_snap_skipped{0};    // doc pdvd/88: 1 = the legacy stop snap's nearest vertex was unreachable from the entry, and was skipped
        int n_local_pieces{0};      // doc pdvd/62 (T3b): disconnected same-cluster pieces admitted into the Michel object
        int n_michel_range_veto{0}; // doc pdvd/62 (T3c): a bridged / charge-only Michel demoted by the range-energy guard
        int topology_cleared_bits{0}; // doc pdvd/70 (P1): the reject bits topology_stop_evidence cleared (0 = none)
        // doc pdvd/71 (P4): the Michel's gamma blobs.  cand = passed every gate;
        // capped = passed and then refused by the total-energy guard.
        int n_michel_gammas{0}, n_michel_gamma_cand{0}, n_michel_gamma_capped{0};
        double michel_ke_gamma{0}, michel_ke_total{0}, michel_gamma_dis_max{-1};
        // doc pdvd/81: the charge-based Michel energy.  valid 0/1; reason 0 ok,
        // 1 no stored response, 2 a chain row's index is outside the response,
        // 3 a chain row's dQ is not the response's solution (a later refit),
        // 4 empty chain, 5 no grouping.  Per-plane sums are SIGNED
        // (measured - muon prediction) electrons; *_mu_* the subtracted muon
        // charge; *_n_* the cells.  q2d / q2d_gamma the plane-combined charge,
        // ke_* their MeV, total = the two added.
        int michel_q2d_valid{0}, michel_q2d_reason{0}, michel_q2d_dropped_plane{-1};
        double michel_q2d_u{0}, michel_q2d_v{0}, michel_q2d_w{0};
        double michel_q2d_mu_u{0}, michel_q2d_mu_v{0}, michel_q2d_mu_w{0};
        int michel_q2d_n_u{0}, michel_q2d_n_v{0}, michel_q2d_n_w{0};
        // the cross-shared cells (charge_err at the fitter's share sentinel: the
        // channel-slice also carries another, non-preloaded cluster's blob) --
        // their count per plane, and the plain sum (measured - muon) over ALL
        // Michel cells, shared included, that the headline replaces on them
        // with the fit's own non-muon prediction.
        int michel_q2d_nx_u{0}, michel_q2d_nx_v{0}, michel_q2d_nx_w{0};
        double michel_q2d_raw_u{0}, michel_q2d_raw_v{0}, michel_q2d_raw_w{0};
        double michel_q2d{0}, michel_q2d_gamma{0};
        double michel_ke_q2d{0}, michel_ke_q2d_gamma{0}, michel_ke_q2d_total{0};
        // doc pdvd/81: every segment that received a row, by role (in row
        // order; a segment can repeat) -- the estimator's fallback source of
        // the Michel (role 3) and its source of the taken gammas (role 4).
        std::map<int, std::vector<SegmentPtr>> role_segs;
        // doc pdvd/81: the admitted companion clusters the fitter produced no
        // segment for (dots_charge_unfit's owners; a conn_type 3 Michel is
        // nothing but these) -- their cells are the Michel's, by blob coverage.
        std::vector<const Cluster*> unfit_dot_clusters;
        // doc pdvd/81: the cell table (PC stm_michel_2d), parallel vectors.
        std::vector<int> c_apa, c_face, c_plane, c_wire, c_time, c_time_slice, c_channel, c_flag, c_role, c_shared, c_xshared, c_sel;
        std::vector<double> c_q, c_qerr, c_pred_mu, c_pred_all;
        int dead_ahead{-1};                        // doc pdhd/03: 1 = the live end walks into a dead region
        int n_cluster_pts{0}; double chain_coverage{-1};   // doc pdhd/03: cluster points within coverage_radius of a reconstructed point
        // dots
        int n_dots{0}, n_dot_clusters_unfit{0};
        double dots_ke_dqdx{0}, dots_charge_unfit{0};
        double dots_ke_unfit{0};                   // doc pdhd/15: dots_charge_unfit converted to MeV
        int michel_n_clusters{0};                  // doc pdvd/51: how many clusters the Michel object spans
        // doc pdvd/51: the capture gammas.  A SEPARATE object list -- never
        // folded into michel_ke_*, so the Michel spectrum keeps its meaning
        // against the free 52.8 MeV endpoint.
        int n_stop_gammas{0}, stop_gamma_n_unfit{0}, stop_gamma_seg_id{-1};
        double stop_gamma_ke_tot{0}, stop_gamma_ke_max{0}, stop_gamma_charge{0};
        double stop_gamma_dis_min{-1}, stop_gamma_dis_max{-1};
        int n_stop_gammas_withheld{0};   // doc pdvd/85
        // doc pdvd/53: the survey.  n_survey_clusters counts admitted companion
        // clusters that yielded at least one UNCLAIMED fitted segment;
        // n_survey_segs counts those segments; n_survey_unfit counts admitted
        // companions that produced no segment at all (so they have charge but
        // no dx, and the display must fall back to their image points).
        int n_survey_clusters{0}, n_survey_segs{0}, n_survey_unfit{0};
        int in_fv{-1};
        // points for the Bee/ROOT layer
        std::vector<double> px, py, pz, pq, pL, prr;
        std::vector<double> pmed;   // doc pdvd/66: the owning segment's median dQ/dx (e/cm) per point; persisted as q_sup = pmed / plateau_med
        std::vector<int> prole, pseg;
        // doc pdvd/53.  Written only when the survey is on, so the knob-off
        // stm_michel_pts schema is unchanged.  prej names the gate that dropped
        // a role-6 segment (0 for a claimed role); pdstop/pdbody are the
        // distances the SHIPPED predicates measured, in WCT units, emitted here
        // rather than recomputed offline -- doc pdvd/51 sec 6.5 is the record of
        // what an offline re-derivation costs.
        std::vector<int> prej;
        std::vector<double> pdstop, pdbody;
        std::set<int> claimed;      // seg_ids already given a role 1/2/3/5
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
            "traj_cover_probe",   // doc pdvd/62: the pr/67 fos census lines (log-only, byte-identical when off)
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
        // doc pdhd/15 sec 7.  A coincident pair of fit points gives dx == 0;
        // cal_kine_dQdx then evaluates the Box model at 0/0, and the NaN survives
        // every clamp into the sum, so ONE such point zeroes the whole object's
        // energy.  doc pdvd/45 sec 5.4 added this knob for exactly that; it is
        // off by default for bit-identicality elsewhere, and this module turns it
        // on because a Michel assembled from several segments meets a zero-dx
        // point often (6 of 280 objects on the d15 arms).  Skipping the point
        // makes it contribute 0, which is what the prototype's (dx + 1e-9) did.
        pa.m_kine_charge.dqdx_skip_zero_dx = true;
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
        // doc pdvd/62: forward the pr/67 find_other_segments census switch so
        // the "already covered" tagging of doc pdvd/54 sec 2.5's 27 lumps can
        // be read on an arm.  Log lines only (NeutrinoOtherSegments.cxx:179).
        B(pa.m_traj_cover_probe, "traj_cover_probe");
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
        add(R_SHORT, "short"); add(R_PROFILE_SPARSE, "profile_sparse");
        add(R_PLATEAU_OFF_MIP, "plateau_off_mip"); add(R_STOP_INTO_DEAD, "stop_into_dead");
        add(R_CLUSTER_NOT_TRACK, "cluster_not_track");
        add(R_PROFILE_GEOMETRY, "profile_geometry");   // doc pdvd/66
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
                        // The cluster is stamped by break_segment itself since
                        // doc sbnd_xin/docs/pr/143; `best` is chosen from
                        // find_cluster_segments(g, cluster) just above, so the
                        // factory writes this cluster.  A clusterless vertex
                        // here is fatal downstream -- examine_direction returns
                        // false at once (NeutrinoVertexFinder.cxx:1503) so
                        // nothing gets oriented, and fill_bee_pf_tree's
                        // main-cluster test (pf_track_main_cluster_only) then
                        // rejects every seed from the main vertex; 039252/2
                        // cluster 86 lost its mu- node and its Michel fell back
                        // to a ROOT shower (doc pdvd/48 sec 8.0).  The explicit
                        // stamp that used to stand here is now redundant.
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

    // doc pdvd/88 (doc pdvd/78 action item 7b): anchor_vertex restricted to
    // what `from` reaches.  Fork by duplication of anchor_vertex above
    // (CLAUDE.md M10) -- the same nearest-vertex-then-split rule, except that
    // a vertex, or a segment with an endpoint, the entry cannot reach is never
    // a candidate.  A pr54 residual kept by the stop-local keep is exactly
    // such a piece (doc pdvd/87 sec 6.2).
    VertexPtr anchor_vertex_reachable(PatternAlgorithms& pa, Graph& g, Cluster& cluster, const Point& pt,
                                      double tol, bool allow_split, double& out_dis, VertexPtr from) const {
        const auto reach_v = stm_michel_reachable_vertices(g, from);
        std::set<VertexPtr> reach(reach_v.begin(), reach_v.end());   // membership only, never iterated
        const Cluster* cl = &cluster;
        auto [vtx, dis] = stm_michel_closest_vertex_of(reach_v, pt,
                                                       [cl](const VertexPtr& v) { return v->cluster() == cl; });
        out_dis = dis;
        if (vtx && dis <= tol) return vtx;
        if (allow_split) {
            SegmentPtr best; double best_d = 1e9; Point best_p;
            for (auto& seg : pa.find_cluster_segments(g, cluster)) {
                if (!seg || seg->fits().size() < 4) continue;
                auto [va, vb] = find_vertices(g, seg);
                if (!va || !vb || !reach.count(va) || !reach.count(vb)) continue;   // doc pdvd/88
                auto [d, p] = segment_get_closest_point(seg, pt, "fit", "main");
                if (d < best_d) { best_d = d; best = seg; best_p = p; }
            }
            if (best && best_d <= tol) {
                try {
                    auto [ok, pair, nvtx] = break_segment(g, best, best_p, particle_data(), m_recomb_model, m_dv,
                                                          1e9 * units::cm, get<bool>(m_cfg, "break_seg_orient", false));
                    if (ok && nvtx) {
                        // The cluster is stamped by break_segment itself since
                        // doc sbnd_xin/docs/pr/143; `best` is chosen from
                        // find_cluster_segments(g, cluster) just above, so the
                        // factory writes this cluster.  A clusterless vertex
                        // here is fatal downstream -- examine_direction returns
                        // false at once (NeutrinoVertexFinder.cxx:1503) so
                        // nothing gets oriented, and fill_bee_pf_tree's
                        // main-cluster test (pf_track_main_cluster_only) then
                        // rejects every seed from the main vertex; 039252/2
                        // cluster 86 lost its mu- node and its Michel fell back
                        // to a ROOT shower (doc pdvd/48 sec 8.0).  The explicit
                        // stamp that used to stand here is now redundant.
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

    // doc pdvd/53: rej/d_stop/d_body ride along on every row so the three
    // arrays stay parallel to the other six.  They are pushed ONLY when the
    // survey is on, which is what keeps the knob-off PC schema identical.
    void add_survey_cols(Record& rec, size_t n, int rej, double d_stop, double d_body) const {
        if (!m_survey_enable && !m_segment_census) return;   // doc pdvd/80: the census carries the same three columns
        for (size_t i = 0; i < n; ++i) {
            rec.prej.push_back(rej); rec.pdstop.push_back(d_stop); rec.pdbody.push_back(d_body);
        }
    }

    // keep_dead: emit a row even for a fit with dx <= 0, carrying the -1 dQ/dx
    // sentinel roles 2/3/4 already use for "the fit found no charge here".  doc
    // pdvd/53: without it a SURVEYED segment whose every fit has dx <= 0 is
    // counted in n_survey_segs, emits no point row, is therefore never named in
    // stm_michel_pts, and is therefore dropped by the display's role-keyed
    // selector -- which is precisely the doc pdvd/51 defect this round exists to
    // fix, surviving on 33 of 872 PDVD and 26 of 610 PDHD surveyed segments.
    // The points are real and drawable: 039252_12 cluster 130's companion
    // segment 429055 has two rows in T_rec_charge, with nq 0 and q at the
    // dQdx_offset sentinel.  Role 6 only, so nothing else moves.
    void add_points(Record& rec, const SegmentPtr& seg, int role, const StmMichelProfile* prof = nullptr,
                    int rej = 0, double d_stop = -1, double d_body = -1,
                    bool keep_dead = false) const {
        const int seg_id = (seg && seg->cluster() ? seg->cluster()->get_cluster_id() : 0) * 1000
                         + (seg ? static_cast<int>(seg->get_graph_index()) : 0);
        if (role != 6 && role != 7 && role != 8) rec.claimed.insert(seg_id);   // 6 survey, 7 other (doc pdvd/64), 8 census (doc pdvd/80): rows, not claims
        if (seg) rec.role_segs[role].push_back(seg);   // doc pdvd/81: bookkeeping only, no row or branch depends on it
        const double seg_med = seg ? segment_median_dQ_dx(seg) * units::cm : -1.0;   // doc pdvd/66: e/cm, the offline dqdx_med's C++ twin
        const size_t n0 = rec.px.size();
        if (prof) {
            for (size_t i = 0; i < prof->pts.size(); ++i) {
                rec.px.push_back(prof->pts[i].x()); rec.py.push_back(prof->pts[i].y()); rec.pz.push_back(prof->pts[i].z());
                rec.pq.push_back(prof->dQdx[i]); rec.pL.push_back(prof->L[i]); rec.prr.push_back(prof->rr[i]);
                rec.prole.push_back(role); rec.pseg.push_back(seg_id);
                rec.pmed.push_back(seg_med);
            }
            add_survey_cols(rec, rec.px.size() - n0, rej, d_stop, d_body);
            return;
        }
        if (!seg) return;
        for (const auto& f : seg->fits()) {
            if (f.dx <= 0 && !keep_dead) continue;
            rec.px.push_back(f.point.x()); rec.py.push_back(f.point.y()); rec.pz.push_back(f.point.z());
            rec.pq.push_back(f.dx > 0 ? f.dQ / (f.dx / units::cm) : -1.0);
            rec.pL.push_back(-1); rec.prr.push_back(-1);
            rec.prole.push_back(role); rec.pseg.push_back(seg_id);
            rec.pmed.push_back(seg_med);
        }
        add_survey_cols(rec, rec.px.size() - n0, rej, d_stop, d_body);
    }

    // doc pdhd/16: p = sqrt((KE + m)^2 - m^2), the segment_cal_4mom idiom
    // (PRSegmentFunctions.cxx:2912).  -1 in, -1 out: a KE that was never
    // computed must not become a momentum of zero.
    static double mom_from_ke(double ke_mev) {
        static const double mmu = 105.658;      // MeV, mcs/src/MuonMCS.cxx:37
        if (!(ke_mev > 0) || !std::isfinite(ke_mev)) return -1;
        const double e = ke_mev + mmu;
        return std::sqrt(std::max(e * e - mmu * mmu, 0.0));
    }

    // doc pdhd/16: the MCS momentum of the muon chain.
    //
    // Mcs::MuonMCS::run() is a plain numeric library -- three vectors of
    // doubles in CENTIMETRES, no graph, no fitter, no particle data (mcs/inc/
    // WireCellMcs/MuonMCS.h:207-222).  vtx_start MUST be the high-energy end:
    // the likelihood walks the CSDA range from there, so for a stopping muon
    // it is the entry and vtx_end is the stop.  The profile's own points are
    // handed over rather than seg->fits() so that the trimmed path is the same
    // object muon_len and muon_ke_dqdx were measured on.
    void fill_mcs(Record& r, const StmMichelProfile& prof) const {
        if (!m_mcs_enable) return;
        if (r.muon_len < m_mcs_min_len_cm * units::cm) return;
        if (prof.pts.size() < 2) return;
        std::vector<std::vector<double>> points;
        points.reserve(prof.pts.size());
        for (const auto& p : prof.pts)
            points.push_back({p.x() / units::cm, p.y() / units::cm, p.z() / units::cm});
        const std::vector<double> start{r.entry_pt.x() / units::cm, r.entry_pt.y() / units::cm,
                                        r.entry_pt.z() / units::cm};
        const std::vector<double> stop{r.stop_pt.x() / units::cm, r.stop_pt.y() / units::cm,
                                       r.stop_pt.z() / units::cm};
        Mcs::McsOptions opt;                 // the five upstream-bug fixes stay ON
        opt.cathode_x = m_mcs_cathode_x;
        opt.cathode_xcut = m_mcs_cathode_xcut;
        const Mcs::McsResult res = Mcs::MuonMCS(opt).run(start, stop, points);
        r.muon_ke_mcs = res.ke_MCS;                  // MeV, -1 when not computed
        r.muon_mcs_amb = res.ambiguity_MCS;          // 1 = maximally ambiguous
        r.muon_mcs_tracklen = res.mu_tracklen;       // cm, the TRIMMED path
        r.muon_mcs_range_ke = res.ke_tracklen;       // MeV, the engine's own CSDA
        r.muon_mcs_nsegs = res.nsegs;
        r.muon_mcs_bad_path = res.bad_path ? 1 : 0;
        r.muon_mcs_cathode_segs = res.counters.cathode_seg_dropped;    // doc pdhd/16 sec 9
        r.muon_mcs_cathode_angles = res.counters.cathode_angle_masked; // incl. the bridging angle
    }

    // Units in the persisted rows (and hence in T_stm_michel): lengths and
    // coordinates in CM, energies in MeV, dQ/dx in e/cm, angles in degrees --
    // the generic PC->TTree writer has no unit knowledge, so the PC carries
    // the human units directly (unlike stm_fit, whose dedicated writer
    // divides by units::cm).
    // ---- doc pdvd/81: the charge-based Michel energy ------------------------
    // The owner's estimator (default_configuration, "michel_q2d"): the Michel's
    // 2-D charge, cell by cell, MINUS the charge the muon chain's own
    // multi-track fit predicts on those cells, combined across planes with the
    // chain's rule and converted with the unfitted-piece constant.  The muon
    // prediction is scale * R * (pos_3D masked to the chain's Fit::index rows)
    // from the response the fitter stored at its LAST multi fit of the main
    // cluster (TrackFitting::Parameters::keep_dqdx_response) -- guarded row by
    // row: every chain row's dQ must BE the stored solution, bit for bit, or the
    // estimate is declared invalid (a refit or re-index after the store).
    //
    // Cells: the union charge maps of every preloaded cluster (the chain's own
    // caches, pa.m_charge_2d_*), walked in map order; a cell is the Michel's
    // (role 3) when the nearest point of a Michel segment's associate_points /
    // fit cloud lies within michel_q2d_dis_cm in that plane (kine_charge_from_
    // maps's rule, NeutrinoEnergyReco.cxx:85-140, cloud fallback included),
    // else a taken gamma's (role 4) by the same test, else the STM's footprint
    // (role 1) when the windowed muon rows predict charge there.  `sel` bit 1 =
    // that nearest-point rule, bit 2 = fit support (a nonzero response entry
    // in a Michel-segment column of the MAIN cluster's fit; companion segments
    // live in their own fits and have none).  The headline sums use bit 1.
    //
    // The Michel segments are the Shower's members when there is one (the
    // in-cluster members get no rows in production, doc pdvd/53), else the
    // role-3 row owners; the gammas are the role-4 row owners.  Companion
    // segments' Fit::index values belong to THEIR cluster's fit, so only
    // main-cluster segments enter the fit-support mask.
    //
    // Writes rec.michel_q2d_* (+ the cell vectors); nothing upstream reads them.
    void michel_q2d_estimate(Record& rec, TrackFitting& tf, PatternAlgorithms& pa, Cluster* main,
                             const std::vector<SegmentPtr>& chain, const IndexedSegmentSet& chain_set,
                             const std::shared_ptr<Shower>& michel_shower) const
    {
        using CoordReadout = TrackFitting::CoordReadout;
        const auto t0 = Clock::now();
        rec.michel_q2d_valid = 0; rec.michel_q2d_reason = 0;
        if (chain.empty()) { rec.michel_q2d_reason = 4; return; }
        const TrackFitting::DqdxResponse* resp = tf.get_dqdx_response(main);
        if (!resp) {
            rec.michel_q2d_reason = 1;
            SPDLOG_LOGGER_DEBUG(s_log, "{}CheckSTM_Michel michel-q2d: cluster {} no stored fit response", m_evt_tag, rec.cluster_id);
            return;
        }
        auto* grouping = tf.grouping();
        if (!grouping) { rec.michel_q2d_reason = 5; return; }
        const Eigen::Index n3 = resp->pos_3D.size();

        // the muon rows (and the windowed subset for the footprint), guarded
        std::vector<char> mask_mu(n3, 0), mask_win(n3, 0), mask_michel(n3, 0);
        int n_rows = 0, n_oor = 0, n_mis = 0;
        const double win = m_michel_q2d_stm_window_cm * units::cm;
        for (const auto& seg : chain) {
            if (!seg) continue;
            for (const auto& f : seg->fits()) {
                ++n_rows;
                if (f.index < 0 || f.index >= n3) { ++n_oor; continue; }
                if (!(resp->pos_3D(f.index) == f.dQ)) { ++n_mis; continue; }
                mask_mu[f.index] = 1;
                if (m_michel_q2d_stm_window_cm < 0 || (f.point - rec.stop_pt).magnitude() <= win) mask_win[f.index] = 1;
            }
        }
        if (n_oor || n_mis) {
            rec.michel_q2d_reason = n_oor ? 2 : 3;
            SPDLOG_LOGGER_DEBUG(s_log, "{}CheckSTM_Michel michel-q2d: cluster {} response mismatch: {} chain rows, {} outside [0,{}), {} dQ != solution",
                                m_evt_tag, rec.cluster_id, n_rows, n_oor, n3, n_mis);
            return;
        }

        // the Michel and gamma segments, graph-index ordered, unique
        auto tidy = [](std::vector<SegmentPtr>& v) {
            v.erase(std::remove(v.begin(), v.end(), nullptr), v.end());
            std::sort(v.begin(), v.end(), [](const SegmentPtr& a, const SegmentPtr& b) { return a->get_graph_index() < b->get_graph_index(); });
            v.erase(std::unique(v.begin(), v.end()), v.end());
        };
        std::vector<SegmentPtr> michel_segs, gamma_segs;
        if (michel_shower) {
            IndexedVertexSet mv; IndexedSegmentSet ms;
            michel_shower->fill_sets(mv, ms, false);
            for (const auto& sg : ms) if (!chain_set.count(sg)) michel_segs.push_back(sg);
        }
        else if (auto it = rec.role_segs.find(3); it != rec.role_segs.end()) {
            michel_segs = it->second;
        }
        if (auto it = rec.role_segs.find(4); it != rec.role_segs.end()) gamma_segs = it->second;
        tidy(michel_segs); tidy(gamma_segs);
        for (const auto* lst : {&michel_segs, &gamma_segs}) {
            for (const auto& sg : *lst) {
                if (sg->cluster() != main) continue;
                for (const auto& f : sg->fits()) if (f.index >= 0 && f.index < n3) mask_michel[f.index] = 1;
            }
        }

        // the four predictions per plane, in electrons, and the row lookup
        std::array<Eigen::VectorXd, 3> pred_mu, pred_all, pred_win, pred_mich;
        std::array<std::unordered_map<CoordReadout, int, TrackFitting::CoordReadoutHash>, 3> row_of;
        {
            const std::vector<char> ones(static_cast<size_t>(n3), 1);
            for (int p = 0; p < 3; ++p) {
                std::vector<double> scale;
                scale.reserve(resp->rows[p].size());
                for (const auto& r : resp->rows[p]) scale.push_back(r.scale);
                pred_mu[p]   = TrackFitting::masked_response_prediction(resp->R[p], resp->pos_3D, mask_mu, scale);
                pred_all[p]  = TrackFitting::masked_response_prediction(resp->R[p], resp->pos_3D, ones, scale);
                pred_win[p]  = TrackFitting::masked_response_prediction(resp->R[p], resp->pos_3D, mask_win, scale);
                pred_mich[p] = TrackFitting::masked_response_prediction(resp->R[p], resp->pos_3D, mask_michel, scale);
                for (size_t i = 0; i < resp->rows[p].size(); ++i) row_of[p].emplace(resp->rows[p][i].key, static_cast<int>(i));
            }
        }

        // the union charge maps -- the chain's caches, collected once
        if (pa.m_charge_2d_u.empty()) pa.collect_charge_maps(tf);
        const PR::ChargeMap* maps[3] = {&pa.m_charge_2d_u, &pa.m_charge_2d_v, &pa.m_charge_2d_w};

        // the clouds, kine_charge_from_maps's pair-and-fallback -- plus the frame.
        // doc pdvd/45: a cell's geometric 2-D point (convert_time_wire_2Dpoint,
        // the raw t0 = 0 drift frame) and a segment cloud's points (the
        // cluster's t0-CORRECTED frame) differ in drift by the cluster's t0
        // drift and any time-origin difference between the two conversions
        // (PDVD's trigger offset) -- metres on PDVD, which is why the chain's
        // own michel_ke_charge is 0 on nearly every PDVD Michel (doc pdhd/16).
        // Measured per segment and (apa, face) from its own fit points, the
        // way the fitter's excl_t0_frame does: the point run backward(t0) into
        // the raw frame, its x minus the point's x (the drift conversions are
        // each other's inverse, so this is the round trip without the tick
        // rounding).  A cloud with no shift for a cell's (apa, face) has no
        // points there and is skipped.
        struct CloudPair {
            std::shared_ptr<const Facade::DynamicPointCloud> a, f;
            std::map<std::pair<int, int>, double> shift;   // (apa, face) -> drift shift, WCT length
        };
        auto clouds_of = [&](const std::vector<SegmentPtr>& segs) {
            std::vector<CloudPair> out;
            for (const auto& sg : segs) {
                CloudPair c;
                c.a = sg->dpcloud("associate_points");
                c.f = sg->dpcloud("fit");
                if (!c.a && !c.f) continue;
                if (!c.a) c.a = c.f;
                if (!c.f) c.f = c.a;
                auto* cl = sg->cluster();
                if (cl && m_pcts) {
                    const auto xform = m_pcts->pc_transform(cl->get_scope_transform(cl->get_default_scope()));
                    const double t0 = cl->get_cluster_t0();
                    std::map<std::pair<int, int>, std::vector<double>> acc;
                    for (const auto& f : sg->fits()) {
                        if (f.paf.first < 0 || f.paf.second < 0) continue;
                        const auto p_raw = xform->backward(geo_point_t(f.point.x(), f.point.y(), f.point.z()), t0, f.paf.second, f.paf.first);
                        acc[f.paf].push_back(p_raw.x() - f.point.x());
                    }
                    for (auto& [paf, v] : acc) {
                        std::sort(v.begin(), v.end());
                        c.shift[paf] = v[v.size() / 2];
                    }
                }
                out.push_back(std::move(c));
            }
            return out;
        };
        const auto michel_clouds = clouds_of(michel_segs);
        const auto gamma_clouds = clouds_of(gamma_segs);
        const double dis_cut = m_michel_q2d_dis_cm * units::cm;
        double dbg_dmin = 1e9; int dbg_nq = 0, dbg_noshift = 0;   // diagnostics for the DEBUG line only
        auto within = [&](const std::vector<CloudPair>& clouds, double drift, double wp, int plane, int face, int apa) {
            for (const auto& c : clouds) {
                auto sit = c.shift.find({apa, face});
                if (sit == c.shift.end()) { ++dbg_noshift; continue; }
                const double d = drift - sit->second;
                auto r1 = c.a->get_closest_2d_point_info_direct(d, wp, plane, face, apa);
                ++dbg_nq;
                if (std::get<0>(r1) >= 0 && std::get<0>(r1) < dbg_dmin) dbg_dmin = std::get<0>(r1);
                if (std::get<0>(r1) >= 0 && std::get<0>(r1) < dis_cut && std::get<1>(r1)) return true;
                if (c.f != c.a) {
                    auto r2 = c.f->get_closest_2d_point_info_direct(d, wp, plane, face, apa);
                    if (std::get<0>(r2) >= 0 && std::get<0>(r2) < dis_cut && std::get<1>(r2)) return true;
                }
            }
            return false;
        };

        // cross-shared cells: update_dQ_dx_data sets charge_err to the share
        // sentinel on a channel-slice that another, NON-preloaded cluster's blob
        // also covers (the fitter's own deweighting; the multi path never
        // restores it).  The measured charge there is not attributable -- on a
        // dense PDHD beam event nearly every Michel cell is one (028084_18/17:
        // 14.2M e measured on 103 W cells the fit predicts 0.28M on).  Rule:
        // measured - muon where the cell is the preloaded clusters' alone,
        // the fit's own non-muon prediction (pred_all - pred_mu, floored at 0)
        // where it is shared -- measured where attributable, fitted where not.
        // The plain sum is persisted beside it (michel_q2d_raw_*).
        const double share_err = tf.get_parameters().share_charge_err;
        auto nticks_map = grouping->get_nticks_per_slice();
        auto nticks_at = [&nticks_map](int apa, int face) {
            auto a = nticks_map.find(apa);
            if (a == nticks_map.end()) return 1;
            auto f = a->second.find(face);
            return f == a->second.end() ? 1 : f->second;
        };
        struct Cell { int role, plane, apa, face, wire, time, time_slice, channel, flag, shared, xshared, sel; double q, qerr, pmu, pall; };
        std::vector<Cell> cells;
        std::array<double, 3> qm{{0, 0, 0}}, mum{{0, 0, 0}}, qg{{0, 0, 0}}, mug{{0, 0, 0}}, rawm{{0, 0, 0}};
        std::array<int, 3> nm{{0, 0, 0}}, ng{{0, 0, 0}}, nxm{{0, 0, 0}};
        int n_role1 = 0, n_sel2_only = 0, n_unfit_cells = 0;
        for (int plane = 0; plane < 3; ++plane) {
            for (const auto& [key, meas] : *maps[plane]) {
                auto wit = pa.m_map_apa_ch_plane_wires.find({key.apa, key.channel});
                if (wit == pa.m_map_apa_ch_plane_wires.end()) continue;
                // every (face, wire) the readout channel maps to in this plane
                // (a channel can serve both faces / wrapped wires; the cloud
                // query is per (apa, face), so each candidate is tested and the
                // one that matches names the row -- the first otherwise)
                std::vector<std::pair<int, int>> fw;
                for (const auto& [f, pl, w] : wit->second) { if (pl == plane) fw.emplace_back(f, w); }
                if (fw.empty()) continue;
                int face = fw.front().first, wire = fw.front().second;
                double pmu = 0, pall = 0, pwin = 0, pmich = 0;
                if (auto rit = row_of[plane].find(key); rit != row_of[plane].end()) {
                    const int i = rit->second;
                    pmu = pred_mu[plane](i); pall = pred_all[plane](i); pwin = pred_win[plane](i); pmich = pred_mich[plane](i);
                }
                auto any_within = [&](const std::vector<CloudPair>& clouds) {
                    for (const auto& [f, w] : fw) {
                        const auto p2d = grouping->convert_time_wire_2Dpoint(key.time, w, key.apa, f, plane);
                        if (within(clouds, p2d.first, p2d.second, plane, f, key.apa)) { face = f; wire = w; return true; }
                    }
                    return false;
                };
                int role = 0, sel = 0;
                if (any_within(michel_clouds)) { role = 3; sel |= 1; }
                else if (any_within(gamma_clouds)) { role = 4; sel |= 1; }
                if (!role && !rec.unfit_dot_clusters.empty()) {
                    // a charge-only piece: the cell is covered by an unfitted
                    // admitted cluster's own blobs (the fitter's predicate, the
                    // fit's tick key, no tolerance)
                    for (const auto* cl : rec.unfit_dot_clusters) {
                        for (const auto& [f, w] : fw) {
                            if (tf.is_cell_covered_by_own_blobs(cl, key.apa, f, plane, w, key.time, 0, nticks_at(key.apa, f))) {
                                role = 3; sel |= 4; face = f; wire = w; ++n_unfit_cells; break;
                            }
                        }
                        if (role) break;
                    }
                }
                if (pmich > 0) { sel |= 2; if (!role) { role = 3; ++n_sel2_only; } }
                if (!role) {
                    if (pwin > 0) { role = 1; ++n_role1; }
                    else continue;
                }
                const bool xshared = meas.charge_err >= share_err;
                const bool headline = (sel & 1) || (sel & 4);
                const double contrib = xshared ? std::max(pall - pmu, 0.0) : (meas.charge - pmu);
                if (role == 3 && headline) {
                    qm[plane] += contrib; rawm[plane] += meas.charge - pmu; mum[plane] += pmu; ++nm[plane];
                    if (xshared) ++nxm[plane];
                }
                if (role == 4 && headline) { qg[plane] += contrib; mug[plane] += pmu; ++ng[plane]; }
                if (m_michel_q2d_cells)
                    cells.push_back({role, plane, key.apa, face, wire, key.time, key.time / std::max(1, nticks_at(key.apa, face)), key.channel, meas.flag,
                                     pmu > 0 ? 1 : 0, xshared ? 1 : 0, sel, meas.charge, meas.charge_err, pmu, pall});
            }
        }

        // qm already holds the per-cell contributions (measured - muon, or the
        // fit's non-muon prediction on a cross-shared cell); qg likewise
        rec.michel_q2d_u = qm[0]; rec.michel_q2d_v = qm[1]; rec.michel_q2d_w = qm[2];
        rec.michel_q2d_raw_u = rawm[0]; rec.michel_q2d_raw_v = rawm[1]; rec.michel_q2d_raw_w = rawm[2];
        rec.michel_q2d_mu_u = mum[0]; rec.michel_q2d_mu_v = mum[1]; rec.michel_q2d_mu_w = mum[2];
        rec.michel_q2d_n_u = nm[0]; rec.michel_q2d_n_v = nm[1]; rec.michel_q2d_n_w = nm[2];
        rec.michel_q2d_nx_u = nxm[0]; rec.michel_q2d_nx_v = nxm[1]; rec.michel_q2d_nx_w = nxm[2];
        // the chain's plane rule, on the planes that HAVE cells: a plane the
        // association reached nothing on (0 cells) takes weight 0 rather than
        // reading as a zero-charge plane (which would make the one populated
        // plane "the largest" and drop it); a zero weight is the rule's own
        // "ignore this plane" (KineChargeOptions::plane_weights).
        const auto& ko = pa.m_kine_charge;
        auto weights_for = [&](const std::array<int, 3>& n) {
            std::array<double, 3> w = ko.plane_weights;
            for (int p = 0; p < 3; ++p) if (n[p] == 0) w[p] = 0;
            return w;
        };
        int dropped = -1;
        rec.michel_q2d = stm_michel_combine_planes({{rec.michel_q2d_u, rec.michel_q2d_v, rec.michel_q2d_w}},
                                                   weights_for(nm), ko.plane_asym_switch, &dropped);
        rec.michel_q2d_dropped_plane = dropped;
        rec.michel_q2d_gamma = stm_michel_combine_planes({{qg[0], qg[1], qg[2]}}, weights_for(ng), ko.plane_asym_switch, nullptr);
        // the unfitted-piece conversion (doc pdhd/17 sec 9), plain MeV out
        auto to_mev = [&](double q) {
            return (m_michel_unfit_from_model
                        ? stm_michel_charge_to_energy_model(q, m_recomb_model, m_michel_unfit_dedx)
                        : stm_michel_charge_to_energy(q, m_michel_unfit_recom, m_michel_unfit_fudge, m_michel_unfit_w_ev))
                   / units::MeV;
        };
        rec.michel_ke_q2d = to_mev(rec.michel_q2d);
        rec.michel_ke_q2d_gamma = to_mev(rec.michel_q2d_gamma);
        rec.michel_ke_q2d_total = rec.michel_ke_q2d + rec.michel_ke_q2d_gamma;
        for (double* e : {&rec.michel_q2d_u, &rec.michel_q2d_v, &rec.michel_q2d_w, &rec.michel_q2d_mu_u, &rec.michel_q2d_mu_v,
                          &rec.michel_q2d_mu_w, &rec.michel_q2d_raw_u, &rec.michel_q2d_raw_v, &rec.michel_q2d_raw_w,
                          &rec.michel_q2d, &rec.michel_q2d_gamma, &rec.michel_ke_q2d,
                          &rec.michel_ke_q2d_gamma, &rec.michel_ke_q2d_total}) {
            if (!std::isfinite(*e)) {
                SPDLOG_LOGGER_WARN(s_log, "{}CheckSTM_Michel michel-q2d: cluster {} produced a non-finite value; zeroed", m_evt_tag, rec.cluster_id);
                *e = 0;
            }
        }
        rec.michel_q2d_valid = 1;

        if (m_michel_q2d_cells) {
            std::sort(cells.begin(), cells.end(), [](const Cell& a, const Cell& b) {
                return std::tie(a.role, a.plane, a.apa, a.face, a.wire, a.time) < std::tie(b.role, b.plane, b.apa, b.face, b.wire, b.time);
            });
            for (const auto& c : cells) {
                rec.c_role.push_back(c.role); rec.c_plane.push_back(c.plane); rec.c_apa.push_back(c.apa); rec.c_face.push_back(c.face);
                rec.c_wire.push_back(c.wire); rec.c_time.push_back(c.time); rec.c_time_slice.push_back(c.time_slice);
                rec.c_channel.push_back(c.channel); rec.c_flag.push_back(c.flag);
                rec.c_shared.push_back(c.shared); rec.c_xshared.push_back(c.xshared); rec.c_sel.push_back(c.sel);
                rec.c_q.push_back(c.q); rec.c_qerr.push_back(c.qerr); rec.c_pred_mu.push_back(c.pmu); rec.c_pred_all.push_back(c.pall);
            }
        }
        SPDLOG_LOGGER_DEBUG(s_log,
            "{}CheckSTM_Michel michel-q2d: cluster {} rows {} michel segs {} gamma segs {} | cells u/v/w {}/{}/{} q-mu {:.0f}/{:.0f}/{:.0f} "
            "(mu {:.0f}/{:.0f}/{:.0f}) dropped {} -> q {:.0f} e = {:.2f} MeV | gamma {}/{}/{} -> {:.2f} MeV | total {:.2f} MeV "
            "(best {:.2f}, total {:.2f}) | footprint cells {} fit-support-only {} unfit-cluster cells {} cross-shared {}/{}/{} raw q-mu {:.0f}/{:.0f}/{:.0f} | clouds {} (pts {}/{}) queries {} dmin {:.3f} cm no-shift-skips {} shift0 {:.2f} cm | {:.0f} ms",
            m_evt_tag, rec.cluster_id, n_rows, michel_segs.size(), gamma_segs.size(),
            nm[0], nm[1], nm[2], rec.michel_q2d_u, rec.michel_q2d_v, rec.michel_q2d_w, mum[0], mum[1], mum[2], dropped,
            rec.michel_q2d, rec.michel_ke_q2d, ng[0], ng[1], ng[2], rec.michel_ke_q2d_gamma, rec.michel_ke_q2d_total,
            rec.michel_ke_best, rec.michel_ke_total, n_role1, n_sel2_only, n_unfit_cells, nxm[0], nxm[1], nxm[2], rawm[0], rawm[1], rawm[2],
            michel_clouds.size(), michel_clouds.empty() ? -1 : (int)michel_clouds.front().a->npoints(),
            michel_clouds.empty() ? -1 : (int)michel_clouds.front().f->npoints(),
            dbg_nq, dbg_dmin / units::cm, dbg_noshift,
            (michel_clouds.empty() || michel_clouds.front().shift.empty()) ? 0.0 : michel_clouds.front().shift.begin()->second / units::cm,
            MS(Clock::now() - t0).count());
    }

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
        I1("n_live_pts", r.n_live_pts); I1("n_dead_pts", r.n_dead_pts); I1("n_cmp_live", r.n_cmp_live); D1("dead_frac_cmp", r.dead_frac_cmp);
        I1("n_delta", r.n_delta); I1("n_body_other", r.n_body_other); I1("n_body_hadron", r.n_body_hadron); D1("delta_len", r.delta_len / cm);
        I1("n_stop_other", r.n_stop_other); I1("n_other_published", r.n_other_published);   // doc pdvd/64
        D1("bragg_anchor_shift_cm", r.bragg_anchor_shift / cm);                                // doc pdvd/65
        D1("end_arc_cm", r.end_arc_cm); D1("end_span_cm", r.end_span_cm); D1("end_arc_span", r.end_arc_span); I1("n_end_pts", r.n_end_pts);   // doc pdvd/66
        I1("n_unsupported_segs", r.n_unsupported_segs); D1("unsupported_len_cm", r.unsupported_len_cm);
        D1("unsupported_frac_min", r.unsupported_frac_min); D1("chain_support_min", r.chain_support_min);
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
        I1("n_ext", r.n_ext); D1("ext_len", r.ext_len / cm); I1("dead_ahead", r.dead_ahead);
        I1("n_stub_absorb", r.n_stub_absorb);   // doc pdvd/63
        I1("n_retreat", r.n_retreat); D1("retreat_len", r.retreat_len / cm);
        I1("n_split", r.n_split); D1("split_len", r.split_len / cm); D1("split_kink_deg", r.split_kink_deg);
        I1("n_michel_veto", r.n_michel_veto);
        I1("n_kept_near_stop_main", r.n_kept_near_stop_main); I1("n_kept_near_stop_comp", r.n_kept_near_stop_comp);   // doc pdvd/62
        I1("n_local_pieces", r.n_local_pieces); I1("n_michel_range_veto", r.n_michel_range_veto);
        // doc pdvd/87: written only when a floor is set (the pattern below),
        // so the floor-off tree keeps its branch list byte-identical.
        if (m_stop_local_residual_min_points > 0 || m_stop_local_residual_min_len_cm > 0)
            I1("n_floored_near_stop", r.n_floored_near_stop);
        // doc pdvd/88: the same pattern -- absent when the knob is off.
        if (m_stop_snap_reachable) I1("stop_snap_skipped", r.stop_snap_skipped);
        // doc pdvd/70 (P1): written only when the knob is on (the survey's
        // pattern), so the knob-off tree keeps its branch list byte-identical.
        if (m_topology_stop_evidence) I1("topology_cleared_bits", r.topology_cleared_bits);
        // doc pdvd/71 (P4): the same pattern -- absent when the knob is off.
        if (m_michel_gamma_collect) {
            I1("n_michel_gammas", r.n_michel_gammas); I1("n_michel_gamma_cand", r.n_michel_gamma_cand);
            I1("n_michel_gamma_capped", r.n_michel_gamma_capped);
            D1("michel_ke_gamma", r.michel_ke_gamma); D1("michel_ke_total", r.michel_ke_total);
            D1("michel_gamma_dis_max", r.michel_gamma_dis_max);
        }
        // doc pdvd/81: the same pattern -- absent when the knob is off.
        if (m_michel_q2d) {
            I1("michel_q2d_valid", r.michel_q2d_valid); I1("michel_q2d_reason", r.michel_q2d_reason);
            I1("michel_q2d_dropped_plane", r.michel_q2d_dropped_plane);
            D1("michel_q2d_u", r.michel_q2d_u); D1("michel_q2d_v", r.michel_q2d_v); D1("michel_q2d_w", r.michel_q2d_w);
            D1("michel_q2d_mu_u", r.michel_q2d_mu_u); D1("michel_q2d_mu_v", r.michel_q2d_mu_v); D1("michel_q2d_mu_w", r.michel_q2d_mu_w);
            I1("michel_q2d_n_u", r.michel_q2d_n_u); I1("michel_q2d_n_v", r.michel_q2d_n_v); I1("michel_q2d_n_w", r.michel_q2d_n_w);
            I1("michel_q2d_nx_u", r.michel_q2d_nx_u); I1("michel_q2d_nx_v", r.michel_q2d_nx_v); I1("michel_q2d_nx_w", r.michel_q2d_nx_w);
            D1("michel_q2d_raw_u", r.michel_q2d_raw_u); D1("michel_q2d_raw_v", r.michel_q2d_raw_v); D1("michel_q2d_raw_w", r.michel_q2d_raw_w);
            D1("michel_q2d", r.michel_q2d); D1("michel_q2d_gamma", r.michel_q2d_gamma);
            D1("michel_ke_q2d", r.michel_ke_q2d); D1("michel_ke_q2d_gamma", r.michel_ke_q2d_gamma);
            D1("michel_ke_q2d_total", r.michel_ke_q2d_total);
        }
        // doc pdvd/72 (P3b): the same pattern.
        if (m_moved_stop_michel_kink_min >= 0) I1("n_michel_veto_exempt", r.n_michel_veto_exempt);
        if (m_moved_stop_michel_reach_min_cm >= 0) I1("n_michel_veto_reach_exempt", r.n_michel_veto_reach_exempt);
        // doc pdvd/74 (P3): the same pattern.
        if (m_retreat_tail_strict || m_retreat_tail_sublive || m_michel_collinear_split ||
            m_stop_tail_peak_frac > 0)
            I1("stop_move_p3_bits", r.stop_move_p3_bits);
        // doc pdvd/75 (P1b): the same pattern.
        if (m_bragg_anchor_geo_fallback) I1("bragg_anchor_fallback", r.bragg_anchor_fallback);
        // doc pdvd/83: the same pattern.
        if (m_michel_near_stop_arm_cm > 0) {
            I1("michel_near_arm", r.michel_near_arm); D1("near_arm_dist_cm", r.near_arm_dist_cm);
            I1("n_near_arms_examined", r.n_near_arms_examined);
        }
        I1("n_cluster_pts", r.n_cluster_pts); D1("chain_coverage", r.chain_coverage);
        I1("n_dots", r.n_dots); I1("n_dot_clusters_unfit", r.n_dot_clusters_unfit);
        D1("dots_ke_dqdx", r.dots_ke_dqdx); D1("dots_charge_unfit", r.dots_charge_unfit);
        D1("dots_ke_unfit", r.dots_ke_unfit);
        // doc pdvd/51.  michel_n_clusters is the one membership fact a consumer
        // reading T_stm_michel ALONE cannot derive: > 1 says the object is not
        // in the candidate's own cluster, which is why doc pdhd/12's PF panel
        // (selector `sub_cluster_id // 1000 == cluster_id`) could not show it.
        // The member list itself is T_stm_michel_pts (role 3/4/5 + seg_id), and
        // it joins the PF stream exactly -- both ids are
        // cluster_id * 1000 + graph_index, verified 1107/1107 on d16vnu.
        I1("michel_n_clusters", r.michel_n_clusters);
        I1("n_stop_gammas", r.n_stop_gammas); I1("stop_gamma_n_unfit", r.stop_gamma_n_unfit);
        I1("stop_gamma_seg_id", r.stop_gamma_seg_id);
        D1("stop_gamma_ke_tot", r.stop_gamma_ke_tot); D1("stop_gamma_ke_max", r.stop_gamma_ke_max);
        D1("stop_gamma_charge", r.stop_gamma_charge);
        D1("stop_gamma_dis_min", r.stop_gamma_dis_min); D1("stop_gamma_dis_max", r.stop_gamma_dis_max);
        if (m_stop_gamma_require_stm) I1("n_stop_gammas_withheld", r.n_stop_gammas_withheld);   // doc pdvd/85
        // doc pdvd/53: the survey.  Written only when the knob is on, so the
        // knob-off T_stm_michel branch list is exactly the doc pdvd/51 one.
        if (m_survey_enable) {
            I1("n_survey_clusters", r.n_survey_clusters); I1("n_survey_segs", r.n_survey_segs);
            I1("n_survey_unfit", r.n_survey_unfit);
        }
        // doc pdvd/80: the segment census.  Written only when the knob is on.
        if (m_segment_census) {
            I1("n_census_segs", r.n_census_segs); I1("n_census_admissible", r.n_census_admissible);
        }
        I1("in_fv", r.in_fv);
        // doc pdhd/14.  The muon energy was the one thing this tree never
        // carried: set_pdg (:661) builds a 4-momentum for every chain segment
        // with segment_cal_4mom(seg, 13, ...) and hangs it on the segment,
        // where nothing ever reads it again.  A hand scan cannot ask whether a
        // 40 MeV blob beside a 250 MeV muon is a Michel if the muon's energy is
        // not in the output at all.  Unconditional (owner decision 2026-09-07:
        // a new feature of a module under active development on both
        // ProtoDUNEs, not a knob) -- so T_stm_michel gains four branches and
        // this output is NOT bit-identical to the pre-doc-14 tree.
        D1("muon_ke_range", r.muon_ke_range); D1("muon_ke_dqdx", r.muon_ke_dqdx);
        D1("muon_ke_best", r.muon_ke_best); I1("michel_seg_id", r.michel_seg_id);
        // doc pdhd/16.  Three energy scales for the same muon, and their
        // momenta: range (the baseline -- charge-blind apart from where the
        // track ends), dQ/dx (calorimetric, and therefore the one that carries
        // the gain x lifetime x recombination normalization), and MCS (purely
        // geometric, blind to charge entirely).  muon_ke_best is deliberately
        // NOT changed: on a stopping muon range is the better estimator and
        // MCS is a cross-check, gated by muon_mcs_amb.
        D1("muon_ke_mcs", r.muon_ke_mcs); D1("muon_mcs_amb", r.muon_mcs_amb);
        D1("muon_mcs_tracklen", r.muon_mcs_tracklen); D1("muon_mcs_range_ke", r.muon_mcs_range_ke);
        I1("muon_mcs_nsegs", r.muon_mcs_nsegs); I1("muon_mcs_bad_path", r.muon_mcs_bad_path);
        I1("muon_mcs_cathode_segs", r.muon_mcs_cathode_segs);
        I1("muon_mcs_cathode_angles", r.muon_mcs_cathode_angles);
        D1("muon_p_range", r.muon_p_range); D1("muon_p_dqdx", r.muon_p_dqdx);
        D1("muon_p_mcs", r.muon_p_mcs);
        // doc pdhd/15.  The Michel is ONE object -- the stop arm (or the seed
        // piece when the 3-D clustering detached it), everything the shower walk
        // reaches, and every fitted segment of an admitted companion cluster.
        // michel_ke_dqdx is now that whole object; michel_ke_core is the core
        // alone, i.e. exactly what michel_ke_dqdx meant through doc pdhd/14, so
        // the old number stays in the tree.  michel_ke_charge is the chain's
        // charge-based estimate of the same object, persisted for comparison and
        // never used as `best` (its recombination factors are the SHOWER pair).
        D1("michel_ke_core", r.michel_ke_core); D1("michel_ke_charge", r.michel_ke_charge);
        I1("michel_n_pieces", r.michel_n_pieces);
        I1("michel_parent_vtx_id", r.michel_parent_vtx_id); D1("michel_dis_cm", r.michel_dis_cm);
        D1("michel_start_x", r.michel_start_pt.x() / cm); D1("michel_start_y", r.michel_start_pt.y() / cm);
        D1("michel_start_z", r.michel_start_pt.z() / cm);
        cluster.local_pcs()["stm_michel"] = Dataset(a);

        if (!r.px.empty()) {
            std::map<std::string, Array> p;
            auto tocm = [&](const std::vector<double>& v) { std::vector<double> o(v); for (auto& x : o) if (x > -0.5) x /= cm; return o; };
            std::vector<double> x(r.px), y(r.py), z(r.pz);
            for (auto& v : x) v /= cm; for (auto& v : y) v /= cm; for (auto& v : z) v /= cm;
            p.emplace("x", Array(x)); p.emplace("y", Array(y)); p.emplace("z", Array(z));
            p.emplace("q", Array(r.pq)); p.emplace("L", Array(tocm(r.pL))); p.emplace("rr", Array(tocm(r.prr)));
            p.emplace("role", Array(r.prole)); p.emplace("seg_id", Array(r.pseg));
            // doc pdvd/66 (T8): per point, the owning segment's median dQ/dx over the
            // chain's plateau median -- the offline class-H "charge_supported" ratio,
            // written by the chain.  Unconditional (every carrier has it), -1 when
            // there is no plateau.
            {
                std::vector<double> qs(r.pmed.size(), -1.0);
                if (r.bragg.plateau_med > 0) for (size_t i = 0; i < qs.size(); ++i) qs[i] = r.pmed[i] > 0 ? r.pmed[i] / r.bragg.plateau_med : -1.0;
                p.emplace("q_sup", Array(qs));
            }
            // doc pdvd/53.  Absent when the survey is off, so write_pc_tree
            // (PdvdPrMagnifyTrackingVisitor.cxx:293, which takes its column set
            // from the first carrier) reproduces the old schema exactly.
            if (m_survey_enable || m_segment_census) {   // doc pdvd/80: the census writes the same three columns
                p.emplace("rej", Array(r.prej));
                p.emplace("d_stop", Array(tocm(r.pdstop)));
                p.emplace("d_body", Array(tocm(r.pdbody)));
            }
            cluster.local_pcs()["stm_michel_pts"] = Dataset(p);
        }
        // doc pdvd/81: the Michel / STM 2-D cells, one row each (see
        // michel_q2d_estimate).  Absent unless michel_q2d_cells is on AND the
        // candidate produced a row -- write_pc_tree creates T_stm_michel_2d from
        // the first carrier, and TensorDM's as_tensors needs same-named PCs to
        // share their columns, so every column is knob-only, never per-cluster.
        if (m_michel_q2d && m_michel_q2d_cells && !r.c_q.empty()) {
            std::map<std::string, Array> c;
            c.emplace("apa", Array(r.c_apa)); c.emplace("face", Array(r.c_face)); c.emplace("plane", Array(r.c_plane));
            c.emplace("wire", Array(r.c_wire)); c.emplace("time", Array(r.c_time)); c.emplace("time_slice", Array(r.c_time_slice));
            c.emplace("channel", Array(r.c_channel));
            c.emplace("flag", Array(r.c_flag)); c.emplace("role", Array(r.c_role)); c.emplace("shared", Array(r.c_shared));
            c.emplace("xshared", Array(r.c_xshared)); c.emplace("sel", Array(r.c_sel));
            c.emplace("charge", Array(r.c_q)); c.emplace("charge_err", Array(r.c_qerr));
            c.emplace("pred_mu", Array(r.c_pred_mu)); c.emplace("pred_all", Array(r.c_pred_all));
            cluster.local_pcs()["stm_michel_2d"] = Dataset(c);
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
    if (m_stop_fv_use_config_tolerance && m_fv_tolerance.size() == 6) fv_tol = m_fv_tolerance;   // doc pdhd/03 sec 6.9
    auto mu_fn = particle_data() ? particle_data()->get_dEdx_function("muon") : nullptr;

    int n_stm = 0, n_michel = 0;
    int n_gamma_cand = 0, n_gamma = 0;            // doc pdvd/51
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
        // doc pdvd/53: the survey widens ADMISSION and nothing else.  The Michel
        // cluster re-test below still uses michel_dot_radius_cm and the gamma
        // ring's outer edge is still stop_gamma_radius_cm, so a cluster reached
        // only by this radius is fitted, given role-6 rows, and claimed by
        // neither stage.  The predicate lives in StmMichelFunctions so its
        // stage-off behaviour is doctested.
        // doc pdvd/71 (P4): the gamma-blob radius reaches only what admission
        // brought, so a radius past today's widens admission -- with the
        // preload perturbation doc pdvd/53 measured.  At its 35 cm default it
        // equals the capture-gamma radius and the max() is a no-op; with the
        // knob off the expression is the old one, untouched.
        const double admit_radius_base = stm_michel_admit_radius(
            m_michel_dot_radius_cm, m_stop_gamma_radius_cm, m_stop_gamma_enable,
            m_survey_radius_cm, m_survey_enable);
        const double admit_radius = m_michel_gamma_collect
            ? std::max(admit_radius_base, m_michel_gamma_radius_cm) : admit_radius_base;
        if (rec.gid >= 0) {
            for (auto* oc : grouping.children()) {
                if (oc == main) continue;
                if (oc->get_scalar<int>("matched_flash_gid", -1) != rec.gid) continue;
                // doc pdhd/15 sec 3: ADMISSION, a separate knob from the per-piece
                // cap below (doc pdhd/13 D2).  A cluster kept out here never
                // reaches the fitter, never gets find_proto_vertex, and can
                // therefore become neither an attached arm nor a dot -- silently,
                // no log line.  PDVD 039252_15 cluster 77's Michel is exactly
                // that: 20.1 cm, 1.91 cm from the stop, 9.12e5 e, dropped by the
                // old 10 cm cap.
                // doc pdvd/53: the survey may carry its own length cap; it
                // defaults to companion_max_len_cm, so the knob-off pool is the
                // same pool.
                const double len_cap = m_survey_enable
                    ? std::max(m_companion_max_len_cm, m_survey_max_len_cm)
                    : m_companion_max_len_cm;
                if (oc->get_length() > len_cap * units::cm) continue;
                if (oc->npoints() == 0) continue;
                const auto [cp, blob] = oc->get_closest_point_blob(rec.stop_pt);
                // doc pdvd/51: ADMISSION reaches as far as the capture-gamma
                // ring; the MICHEL's own radius is unchanged and is re-applied
                // at cluster level in the piece loop below, so nothing the
                // Michel object gathers moves.  Admitting more companions does
                // perturb the candidate's own fit through preload_clusters --
                // that is measured, not assumed (doc pdvd/51 sec 6).
                if ((cp - rec.stop_pt).magnitude() > admit_radius * units::cm) continue;
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
        if (m_michel_q2d) tf->set_parameter("keep_dqdx_response", 1.0);   // doc pdvd/81: store every multi-fit's response
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
        // doc pdvd/62 (T3a): anchor the partition's isolated-residual keep on
        // the TAGGER's stop -- the only stop known before the partition runs
        // (read_stm_anchor above; the chain's own stop does not exist yet).
        if (m_stop_local_residual_cm > 0) {
            pa.m_other_seg_keep_anchor_cm = m_stop_local_residual_cm * units::cm;
            pa.m_other_seg_keep_anchors = {rec.tagger_stop_pt};
            // doc pdvd/87: the size floor, 0 / 0 = none.
            pa.m_other_seg_keep_anchor_min_points = m_stop_local_residual_min_points;
            pa.m_other_seg_keep_anchor_min_length = m_stop_local_residual_min_len_cm * units::cm;
        }

        // ---- the four PR stages on the main cluster ----------------------
        const bool ok_main = pa.find_proto_vertex(g, *main, *tf, m_dv, true, 2, true, particle_data());
        rec.n_kept_near_stop_main = pa.m_other_seg_keep_anchor_fires;
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
        rec.n_kept_near_stop_comp = pa.m_other_seg_keep_anchor_fires - rec.n_kept_near_stop_main;   // doc pdvd/62
        rec.n_floored_near_stop = pa.m_other_seg_keep_anchor_floored;                                 // doc pdvd/87

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
            if (m_stop_snap_reachable && entry_v && rec.n_kept_near_stop_main > 0) {
                // doc pdvd/88 (item 7b): a residual the stop-local keep added
                // is disconnected by construction and must not capture the
                // stop.  The legacy nearest vertex is read (read-only) only to
                // record whether it was one the entry cannot reach.
                const auto reach0 = stm_michel_reachable_vertices(g, entry_v);
                const auto [v_any, d_any] = pa.closest_cluster_vertex(g, *main, rec.stop_pt);
                rec.stop_snap_skipped = (v_any && std::find(reach0.begin(), reach0.end(), v_any) == reach0.end()) ? 1 : 0;
                stop_v = anchor_vertex_reachable(pa, g, *main, rec.stop_pt, m_stop_snap_tol_cm * units::cm, true, stop_dis, entry_v);
                SPDLOG_LOGGER_DEBUG(s_log, "{}stop-snap-reachable: cluster {} skipped {} (nearest vertex d={:.2f} cm) -> stop d={:.2f} cm",
                                    m_evt_tag, main->get_cluster_id(), rec.stop_snap_skipped, d_any / units::cm, stop_dis / units::cm);
            }
            else {
                stop_v  = anchor_vertex(pa, g, *main, rec.stop_pt,  m_stop_snap_tol_cm * units::cm, true, stop_dis);
            }
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
            // No vertex near the tagger's stop, or it is unreachable from the
            // entry (doc pdhd/03: the tagger's single-track fit bridged into a
            // detached fragment 20-40 cm past the track's graph end, so the
            // stop it recorded is in another component).  The muon is then the
            // longest route out of the entry within the main cluster; the
            // verdict keeps R_STOP_UNMATCHED.  (doc pdvd/48 walked greedily
            // with find_cont_muon_segment instead, which on PDHD stopped at
            // the first junction: 1-segment chains on 8 of 13 such cases.)
            rec.reject_bits |= R_STOP_UNMATCHED;
            Cluster* main_ptr = main;
            VertexPtr far = stm_michel_farthest_vertex(g, entry_v,
                [main_ptr](const VertexPtr& v) { return v->cluster() == main_ptr; });
            if (far) {
                chain = stm_michel_shortest_chain(g, entry_v, far);
                stop_v = chain.empty() ? nullptr : far;
            }
        }
        if (chain.empty()) rec.reject_bits |= R_NO_CHAIN;

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
        th.michel_shower_min_kink_deg = m_michel_shower_min_kink_deg;
        // doc pdvd/73 (P2): -1 stays -1 (off); the two lengths go to internal units only when set
        th.michel_mip_lo_turned = m_michel_mip_lo_turned;
        th.michel_mip_lo_turned_kink_deg = m_michel_mip_lo_turned_kink_deg;
        th.michel_far_len_shower_max = m_michel_far_len_shower_max_cm >= 0 ? m_michel_far_len_shower_max_cm * units::cm : -1.0;
        th.michel_kink_window = m_michel_kink_window_cm > 0 ? m_michel_kink_window_cm * units::cm : -1.0;
        th.delta_max_len = m_delta_max_len_cm * units::cm;
        th.hadron_mip = m_vertex_hadron_mip;

        // ---- doc pdhd/03 sec 6: the tagger's stop is often short of the muon's
        // end (doc pdvd/42 sec 4.4: a collinear ~0.9 MIP leftover on 26 % of
        // PDVD passes).  When the stop vertex has a continuation arm, the muon
        // did not stop there: extend the chain along it and judge at the new
        // end.  Up to stop_extend_max times; 0 = doc pdvd/48 behaviour.
        std::function<double(double)> mu_at_ext = nullptr;
        if (mu_fn) mu_at_ext = [mu_fn](double rr_cm) { return mu_fn->scalar_function(rr_cm); };
        auto bragg_confirmed = [&](const std::vector<SegmentPtr>& ch) {
            // the live profile of the current chain already shows the rise
            auto pr = stm_michel_profile(g, ch, entry_v);
            int nd = 0;
            auto lv = stm_michel_profile_live(pr, m_profile_min_dqdx_frac * m_mip_dqdx, nd);
            auto b = stm_michel_bragg_contrast(lv, mu_at_ext,
                                               m_bragg_tail_lo_cm * units::cm, m_bragg_tail_hi_cm * units::cm,
                                               m_bragg_plateau_lo_cm * units::cm, m_bragg_plateau_hi_cm * units::cm);
            return b.valid && b.expected > 0 && b.contrast >= m_bragg_contrast_min * b.expected;
        };
        for (int it = 0; it < m_stop_extend_max && stop_v && !chain.empty(); ++it) {
            SegmentPtr last = chain.back();
            StmMichelArm best, stub; bool has_michel = false;
            for (auto e : sorted_out_edges(stop_v->get_descriptor(), g)) {
                auto arm = g[e].segment;
                if (!arm || arm == last) continue;
                if (std::find(chain.begin(), chain.end(), arm) != chain.end()) continue;
                auto a = stm_michel_classify_stop_arm(g, last, arm, stop_v, th);
                if (a.kind == StmMichelArm::kMichel) has_michel = true;
                if (a.kind == StmMichelArm::kContinuation && a.len > best.len) best = a;
                // doc pdhd/03 sec 6.8: a short collinear arm HOTTER than a MIP is
                // the muon's own Bragg stub that the partition split off -- it
                // belongs to the chain, whatever the shower flag says.
                if (m_absorb_bragg_stub && a.kind != StmMichelArm::kContinuation && a.kink_deg >= 0 &&
                    a.kink_deg < m_continuation_max_angle_deg && a.len <= m_delta_max_len_cm * units::cm &&
                    a.mip > m_continuation_mip_hi && a.len > stub.len) stub = a;
            }
            if (stub.seg) best = stub;
            if (!best.seg) break;
            // doc pdhd/03 sec 6.3: a Michel arm, or a Bragg rise already in the
            // profile, says the muon stopped HERE -- do not walk past it (a
            // Bragg stub is absorbed regardless: it IS the stop).
            if (m_michel_guards_stop && !stub.seg && (has_michel || bragg_confirmed(chain))) break;
            VertexPtr far = find_other_vertex(g, best.seg, stop_v);
            if (!far || far == entry_v) break;
            chain.push_back(best.seg);
            stop_v = far;
            ++rec.n_ext; rec.ext_len += best.len;
            if (stub.seg) ++rec.n_stub_absorb;   // doc pdvd/63 (T5): this extension was absorb_bragg_stub's, not a continuation's
        }
        if (rec.n_ext > 0 && stop_v) {
            rec.stop_vtx_id = rec.cluster_id * 1000 + static_cast<int>(stop_v->get_graph_index());
            rec.stop_pt = stm_michel_vertex_point(stop_v);
            rec.stop_dis = (rec.stop_pt - rec.tagger_stop_pt).magnitude();
            rec.n_chain_segs = static_cast<int>(chain.size());
        }

        // ---- doc pdvd/57: the STOP RETREAT.  Mirror of stop_extend_max above,
        // run AFTER it so a wrong extension can still be undone.  find_first_kink
        // returns its no-kink sentinel whenever the muon->Michel junction is
        // asymmetric (Bragg on one side, 0.1-0.4 MIP on the other -- doc
        // pdvd/56 sec 2.2 measured this on 74 of 131 accepted stoppers and 29
        // of 51 missed collapse-shaped ones), so the tagger's stop is the
        // Steiner path's far end, not where the muon stopped.  Graph-only:
        // it can only retreat onto a vertex the chain ALREADY has, which is
        // what keeps it from firing on a through-going track that merely
        // trails off (doc pdvd/56 sec 2.3 / doc pdvd/57 sec 3: an on-chain
        // vertex within 4 cm of the collapse exists on 24 of 51 missed
        // collapse-shaped stoppers but only 11 of 80 collapse-shaped
        // through-going ones).  R_STOP_UNMATCHED chains were never anchored
        // to the tagger's stop at all -- retreating one would mean something
        // different, so they are excluded.
        // doc pdvd/74 (P3, doc 70's literal proposal): michel_collinear_split
        // lets this retreat and the split below also run on a chain whose
        // profile already shows the Bragg rise -- the stopper is then found,
        // and what the fit carried past the peak is the Michel's.  Off, the
        // condition is the doc 57/58 one: bragg_confirmed is pure, and it is
        // evaluated exactly when it was before.
        const bool retreat_try = m_stop_retreat_max > 0 && stop_v && chain.size() > 1 &&
                                 !(rec.reject_bits & R_STOP_UNMATCHED);
        const bool retreat_bragg = retreat_try && bragg_confirmed(chain);
        if (retreat_try && (!retreat_bragg || m_michel_collinear_split)) {
            auto prof_pre = stm_michel_profile(g, chain, entry_v);
            StmMichelRetreatThresholds rth;
            rth.max_drop = m_stop_retreat_max;
            rth.collapse_frac = m_retreat_collapse_frac;
            rth.peak_frac = m_retreat_peak_frac;
            rth.peak_window = m_retreat_peak_window_cm * units::cm;
            rth.max_drop_len = m_michel_max_len_cm * units::cm;   // same ceiling an attached Michel arm faces
            rth.plateau_lo = m_bragg_plateau_lo_cm * units::cm;
            rth.plateau_hi = m_bragg_plateau_hi_cm * units::cm;
            rth.min_dqdx_live = m_profile_min_dqdx_frac * m_mip_dqdx;
            rth.tail_strict = m_retreat_tail_strict;     // doc pdvd/74 (P3): how the dropped tail is read
            rth.tail_sublive = m_retreat_tail_sublive;
            rth.tail_peak_frac = m_stop_tail_peak_frac;  // doc pdvd/82: the peak-relative admission
            rth.tail_peak_kink_min = m_stop_tail_peak_kink_min_deg;
            rth.dir_window = m_split_dir_window_cm * units::cm;   // one bend definition for both movers
            auto rr = stm_michel_stop_retreat(prof_pre, static_cast<int>(chain.size()), rth);
            // doc pdvd/74: did that reading change the answer?  The doc 57
            // reading is re-run (a pure function) only when it can differ.
            if (m_retreat_tail_strict || m_retreat_tail_sublive) {
                auto lth = rth;
                lth.tail_strict = false;
                lth.tail_sublive = false;
                const int legacy_drop = stm_michel_stop_retreat(prof_pre, static_cast<int>(chain.size()), lth).n_drop;
                if (legacy_drop != rr.n_drop) {
                    rec.stop_move_p3_bits |= 1;
                    SPDLOG_LOGGER_DEBUG(s_log, "{}CheckSTM_Michel P3: cluster {} retreat n_drop {} ({:.2f} cm) where the doc-57 tail reading gave {}",
                                        m_evt_tag, rec.cluster_id, rr.n_drop, rr.drop_len / units::cm, legacy_drop);
                }
            }
            if (rr.n_drop > 0) {
                auto vtxs = stm_michel_chain_vertices(g, chain, entry_v);
                VertexPtr new_stop = (vtxs.size() == chain.size() + 1)
                    ? vtxs[chain.size() - rr.n_drop] : nullptr;
                if (new_stop && new_stop != entry_v) {
                    chain.resize(chain.size() - rr.n_drop);
                    stop_v = new_stop;
                    rec.n_retreat = rr.n_drop;
                    rec.retreat_len = rr.drop_len;
                    if (rr.by_tail_peak) {   // doc pdvd/82
                        rec.stop_move_p3_bits |= 4;
                        SPDLOG_LOGGER_DEBUG(s_log, "{}CheckSTM_Michel doc82: cluster {} retreat admitted by the peak-relative tail, {} seg(s) {:.2f} cm, tail {:.0f} peak {:.0f} plateau {:.0f} bend {:.1f} deg",
                                            m_evt_tag, rec.cluster_id, rr.n_drop, rr.drop_len / units::cm,
                                            rr.last_tail_med, rr.last_peak, rr.plateau, rr.last_kink_deg);
                    }
                    if (retreat_bragg) {
                        rec.stop_move_p3_bits |= 2;
                        SPDLOG_LOGGER_DEBUG(s_log, "{}CheckSTM_Michel P3: cluster {} retreat on a Bragg-confirmed chain, n_drop {} ({:.2f} cm)",
                                            m_evt_tag, rec.cluster_id, rr.n_drop, rr.drop_len / units::cm);
                    }
                }
            }
        }

        // ---- doc pdvd/58 (T1c): the STOP SPLIT.  stm_michel_stop_retreat
        // above can only move the stop onto a vertex the chain already has;
        // a population of missed collapse-shaped stoppers has no such vertex
        // at all -- the collapse sits INSIDE the fit's last chain segment
        // (doc pdvd/57 sec 1's 23-item finding).  Only tried when the retreat
        // did NOT already fire (n_retreat == 0): the retreat is the cheaper,
        // safer mechanism (an existing vertex, no graph mutation) and always
        // wins when it applies. Same anchoring guards as the retreat, plus
        // the >= 4 fits anchor_vertex():901 already requires of a segment it
        // is willing to split.
        const bool split_try = m_stop_split_max > 0 && rec.n_retreat == 0 && stop_v && !chain.empty() &&
                               chain.back() && chain.back()->fits().size() >= 4 &&
                               !(rec.reject_bits & R_STOP_UNMATCHED);
        const bool split_bragg = split_try && bragg_confirmed(chain);   // doc pdvd/74 (P3): see the retreat above
        if (split_try && (!split_bragg || m_michel_collinear_split)) {
            auto prof_pre = stm_michel_profile(g, chain, entry_v);
            StmMichelSplitThresholds sth;
            sth.max_split = m_stop_split_max;
            sth.kink_min_deg = m_split_kink_min_deg;
            sth.min_drop = m_split_min_drop_cm * units::cm;
            sth.collapse_frac = m_split_collapse_frac;
            sth.peak_frac = m_split_peak_frac;
            sth.peak_window = m_split_peak_window_cm * units::cm;
            sth.dir_window = m_split_dir_window_cm * units::cm;
            sth.tail_peak_frac = m_stop_tail_peak_frac;  // doc pdvd/82: the peak-relative admission
            sth.tail_peak_kink_min = m_stop_tail_peak_kink_min_deg;
            sth.max_drop_len = m_michel_max_len_cm * units::cm;   // same ceiling the retreat uses
            sth.plateau_lo = m_bragg_plateau_lo_cm * units::cm;
            sth.plateau_hi = m_bragg_plateau_hi_cm * units::cm;
            sth.min_dqdx_live = m_profile_min_dqdx_frac * m_mip_dqdx;
            auto sp = stm_michel_stop_split(prof_pre, static_cast<int>(chain.size()), sth);
            if (sp.ok) {
                try {
                    auto [ok, pair, nvtx] = break_segment(g, chain.back(), prof_pre.pts[sp.index],
                                                          particle_data(), m_recomb_model, m_dv,
                                                          1e9 * units::cm,
                                                          get<bool>(m_cfg, "break_seg_orient", false));
                    if (ok && nvtx) {
                        // No examine_vertices*/examine_structure_final* pass
                        // runs downstream of this component, but the flag is
                        // what protects a split from the merge family in
                        // general and costs nothing to set.
                        nvtx->set_flags(VertexFlags::kProtectedBreak);
                        auto new_chain = stm_michel_shortest_chain(g, entry_v, nvtx);
                        if (!new_chain.empty()) {
                            chain = new_chain;
                            stop_v = nvtx;
                            rec.n_split = 1;
                            rec.split_len = sp.drop_len;
                            rec.split_kink_deg = sp.kink_deg;
                            if (sp.by_tail_peak) {   // doc pdvd/82
                                rec.stop_move_p3_bits |= 4;
                                SPDLOG_LOGGER_DEBUG(s_log, "{}CheckSTM_Michel doc82: cluster {} split admitted by the peak-relative tail, {:.2f} cm at {:.1f} deg, tail {:.0f} peak {:.0f} plateau {:.0f}",
                                                    m_evt_tag, rec.cluster_id, sp.drop_len / units::cm, sp.kink_deg,
                                                    sp.tail_med, sp.peak, sp.plateau);
                            }
                            if (split_bragg) {   // doc pdvd/74 (P3)
                                rec.stop_move_p3_bits |= 2;
                                SPDLOG_LOGGER_DEBUG(s_log, "{}CheckSTM_Michel P3: cluster {} split on a Bragg-confirmed chain, {:.2f} cm at {:.1f} deg",
                                                    m_evt_tag, rec.cluster_id, sp.drop_len / units::cm, sp.kink_deg);
                            }
                        }
                        // else: the split already happened to the graph (it
                        // cannot be undone), but no route from entry to the
                        // new vertex exists -- leave chain/stop_v as they
                        // were before the split and fall through unrecorded,
                        // same as any other "no chain" outcome below.
                    }
                } catch (const std::exception& e) {
                    SPDLOG_LOGGER_WARN(s_log, "{}stop_split: break_segment threw: {}", m_evt_tag, e.what());
                }
            }
        }

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
            // doc pdhd/14: the same two estimators Shower::calculate_kinematics
            // uses, and its >= 4 cm rule for choosing between them
            // (PRSegmentFunctions.cxx:2900).  Range is taken on the
            // CHAIN's total length, not per segment: cal_kine_range is not
            // additive, so summing per-segment range energies over a chain the
            // stop-extension walked would be wrong.  dQ/dx is additive and is
            // summed segment by segment, exactly as :1262 sums the dots.
            // cal_kine_range dereferences particle_data without checking it,
            // and this component tolerates a null one elsewhere (:797), so
            // guard rather than inherit a crash on a config that omits it.
            if (particle_data())
                rec.muon_ke_range = cal_kine_range(rec.muon_len, 13, particle_data()) / units::MeV;
            for (auto& s : chain) rec.muon_ke_dqdx += segment_cal_kine_dQdx(s, m_recomb_model) / units::MeV;
            rec.muon_ke_best = (rec.muon_len < 4 * units::cm) ? rec.muon_ke_dqdx : rec.muon_ke_range;
            fill_mcs(rec, prof);        // doc pdhd/16
            rec.muon_p_range = mom_from_ke(rec.muon_ke_range);
            rec.muon_p_dqdx = mom_from_ke(rec.muon_ke_dqdx);
            rec.muon_p_mcs = mom_from_ke(rec.muon_ke_mcs);
            add_points(rec, chain.back(), 1, &prof);
            if (rec.n_profile_pts < m_min_chain_points) rec.reject_bits |= R_SHORT;
        }

        // ---- dQ/dx vs residual range: contrast, KS shape, template PID ----
        // All three metrics read the LIVE profile: points with dQ/dx below
        // profile_min_dqdx_frac * mip are cells the fit could not read (dead
        // channels, APA edges, wrapped-wire ambiguity on PDHD) and carry no
        // particle information (doc pdhd/03 sec 5; frac = 0 keeps everything).
        // A window with fewer than 3 live points is "cannot judge" ->
        // R_PROFILE_SPARSE, not a shape or PID verdict.
        if (!prof.empty()) {
            StmMichelProfile live = stm_michel_profile_live(prof, m_profile_min_dqdx_frac * m_mip_dqdx, rec.n_dead_pts);
            rec.n_live_pts = static_cast<int>(live.L.size());
            {
                int n_cmp = 0, n_cmp_dead = 0;
                for (size_t i = 0; i < prof.rr.size(); ++i) {
                    if (prof.rr[i] > m_compare_range_cm * units::cm) continue;
                    ++n_cmp;
                    if (prof.dQdx[i] < m_profile_min_dqdx_frac * m_mip_dqdx) ++n_cmp_dead;
                }
                rec.n_cmp_live = n_cmp - n_cmp_dead;
                rec.dead_frac_cmp = n_cmp > 0 ? double(n_cmp_dead) / n_cmp : 0.0;
            }
            // doc pdvd/65 (T7): the peak-anchored residual-range origin -- the
            // prototype eval_stm / TaggerCheckSTM::eval_stm_core_impl recipe
            // (end_L = L[max_bin] + 0.2 cm after a 5-point running-mean peak
            // search), applied to the LIVE profile within bragg_peak_search_cm
            // of the geometric end, for the two verdict-stage shape tests
            // ONLY.  The chain walk (bragg_confirmed: extend / retreat / split
            // guards) keeps the geometric origin -- it decides WHERE the stop
            // is; this decides where the peak is read.  rr' = end_L - L, rows
            // past the peak dropped; n_live_pts / dead_frac_cmp above stay
            // geometric.  Off => `live` untouched, byte-identical.
            // doc pdvd/75 (P1b): the geometric live profile and the anchor's
            // winning 5-point mean, kept for the fallback below.  Unused when
            // the fallback is off.
            const StmMichelProfile live_geo = live;
            double anchor_peak_mean = -1;
            if (m_bragg_peak_anchor && live.L.size() >= 5) {
                const size_t n = live.L.size();
                double best = -1; size_t ibest = n - 1;
                for (size_t i = 0; i < n; ++i) {
                    if (live.rr[i] > m_bragg_peak_search_cm * units::cm) continue;
                    double s = 0; int c = 0;
                    for (int d = -2; d <= 2; ++d) {
                        const long j = static_cast<long>(i) + d;
                        if (j < 0 || j >= static_cast<long>(n)) continue;
                        s += live.dQdx[j]; ++c;
                    }
                    s /= c;
                    if (s > best) { best = s; ibest = i; }
                }
                const double end_L = live.L[ibest] + 0.2 * units::cm;
                if (end_L < live.total_length) {          // the anchor only moves the origin BACK
                    StmMichelProfile anch;
                    for (size_t i = 0; i < n; ++i) {
                        const double r = end_L - live.L[i];
                        if (r < 0) continue;
                        anch.L.push_back(live.L[i]); anch.dQdx.push_back(live.dQdx[i]); anch.rr.push_back(r);
                        anch.pts.push_back(live.pts[i]); anch.seg_idx.push_back(live.seg_idx[i]);
                    }
                    anch.total_length = end_L;
                    if (anch.L.size() >= 3) {
                        rec.bragg_anchor_shift = live.total_length - end_L;
                        live = anch;
                        anchor_peak_mean = best;   // doc pdvd/75
                    }
                }
            }
            std::function<double(double)> mu_at = nullptr;
            if (mu_fn) mu_at = [mu_fn](double rr_cm) { return mu_fn->scalar_function(rr_cm); };
            rec.bragg = stm_michel_bragg_contrast(live, mu_at,
                                                  m_bragg_tail_lo_cm * units::cm, m_bragg_tail_hi_cm * units::cm,
                                                  m_bragg_plateau_lo_cm * units::cm, m_bragg_plateau_hi_cm * units::cm);
            if (!rec.bragg.valid || rec.bragg.expected <= 0) {
                rec.reject_bits |= R_PROFILE_SPARSE;
            }
            else if (rec.bragg.contrast < m_bragg_contrast_min * rec.bragg.expected) {
                rec.reject_bits |= R_NO_BRAGG;
            }
            if (rec.bragg.valid && m_plateau_mip_hi > m_plateau_mip_lo && m_mip_dqdx > 0) {
                const double pm = rec.bragg.plateau_med / m_mip_dqdx;
                if (pm < m_plateau_mip_lo || pm > m_plateau_mip_hi) rec.reject_bits |= R_PLATEAU_OFF_MIP;
            }
            // doc pdhd/03 sec 6: does the visible end walk into a dead region?
            // Direction = the last 10 cm of the chain; probe from the last LIVE
            // point (where the detector stopped seeing) and from the stop.
            if (m_dead_volume_check && fiducial_utils && !live.pts.empty() && prof.pts.size() >= 2) {
                const Point& pend = prof.pts.back();
                size_t k = prof.pts.size() - 1;
                while (k > 0 && (prof.L.back() - prof.L[k]) < 10 * units::cm) --k;
                Vector dir = pend - prof.pts[k];
                if (dir.magnitude() > 0) {
                    const bool live_ok = fiducial_utils->check_dead_volume(*main, live.pts.back(), dir, 1 * units::cm);
                    const bool stop_ok = fiducial_utils->check_dead_volume(*main, rec.stop_pt, dir, 1 * units::cm);
                    rec.dead_ahead = (live_ok && stop_ok) ? 0 : 1;
                    if (rec.dead_ahead) rec.reject_bits |= R_STOP_INTO_DEAD;
                }
            }
            // The TaggerCheckSTM::eval_stm_core_impl recipe (:2780-2797) over
            // the last compare_range of residual range, e/cm frame.
            if (mu_fn) {
                std::vector<double> test, ref_mu, ref_flat;
                for (size_t i = 0; i < live.rr.size(); ++i) {
                    if (live.rr[i] > m_compare_range_cm * units::cm) continue;
                    test.push_back(live.dQdx[i]);
                    ref_mu.push_back(mu_fn->scalar_function(live.rr[i] / units::cm + m_offset_length_cm));
                    ref_flat.push_back(m_mip_dqdx);
                }
                if (test.size() >= 3) {
                    auto sum = [](const std::vector<double>& v) { double s = 0; for (double x : v) s += x; return s; };
                    rec.ks_mu = WireCell::kslike_compare(test, ref_mu);
                    rec.ratio_mu = sum(ref_mu) / (sum(test) + 1e-9);
                    rec.ks_flat = WireCell::kslike_compare(test, ref_flat);
                    rec.ratio_flat = sum(ref_flat) / (sum(test) + 1e-9);
                    if (rec.ks_mu + m_ks_margin >= rec.ks_flat) rec.reject_bits |= R_SHAPE_FLAT;
                    // do_track_comp: internal-unit arrays; forward = stop at L.back().
                    std::vector<double> L_f(live.L), q_f(live.dQdx.size());
                    for (size_t i = 0; i < q_f.size(); ++i) q_f[i] = live.dQdx[i] / units::cm;
                    std::vector<double> L_b(L_f.size()), q_b(q_f.size());
                    const size_t n = L_f.size();
                    for (size_t i = 0; i < n; ++i) { L_b[i] = live.total_length - L_f[n - 1 - i]; q_b[i] = q_f[n - 1 - i]; }
                    rec.comp_fwd = do_track_comp(L_f, q_f, m_compare_range_cm * units::cm, m_offset_length_cm * units::cm,
                                                 particle_data(), m_mip_dqdx / units::cm);
                    rec.comp_bwd = do_track_comp(L_b, q_b, m_compare_range_cm * units::cm, m_offset_length_cm * units::cm,
                                                 particle_data(), m_mip_dqdx / units::cm);
                    if (rec.comp_fwd.size() == 4) {
                        const bool gate = rec.comp_fwd[0] > 0.5;
                        const bool beats_p = rec.comp_fwd[1] < rec.comp_fwd[2];
                        const bool beats_e = rec.comp_fwd[1] < rec.comp_fwd[3];
                        const bool muon_like = (m_pid_mode == 2) ? beats_p
                                             : (m_pid_mode == 1) ? (gate && beats_p)
                                             : (gate && beats_p && beats_e);
                        if (!muon_like) rec.reject_bits |= R_NOT_MUON_PID;
                    }
                }
                else {
                    rec.reject_bits |= R_PROFILE_SPARSE;
                }
            }
            // ---- doc pdvd/75 (P1b): the geometric reading may stand ---------
            // The anchor above only ever moves the origin back; on a Bragg
            // rise that runs to the fit's last row it discards the top of the
            // rise (the low partial-step end row pulls the 5-point maximum
            // back), and the anchored profile reads flat.  When the anchor
            // fired, rejected on a shape bit, and its own peak is prominent
            // against the anchored plateau, the four shape tests are re-read
            // at the geometric origin with the very same expressions; a
            // reading that sets none of the four bits replaces the anchored
            // one.  Bits can only be cleared here.  do_track_comp and the
            // dead-volume probe above keep the anchored profile.  Off =>
            // nothing below runs, byte-identical.
            {
                const int shape_bits = R_NO_BRAGG | R_SHAPE_FLAT | R_PLATEAU_OFF_MIP | R_PROFILE_SPARSE;
                if (m_bragg_anchor_geo_fallback && mu_fn && rec.bragg_anchor_shift > 0 && anchor_peak_mean > 0
                    && rec.bragg.valid && rec.bragg.plateau_med > 0
                    && anchor_peak_mean >= m_bragg_anchor_rise_min * rec.bragg.plateau_med
                    && (rec.reject_bits & shape_bits)) {
                    const StmMichelBragg geo = stm_michel_bragg_contrast(live_geo, mu_at,
                                                                         m_bragg_tail_lo_cm * units::cm, m_bragg_tail_hi_cm * units::cm,
                                                                         m_bragg_plateau_lo_cm * units::cm, m_bragg_plateau_hi_cm * units::cm);
                    int geo_bits = 0;
                    if (!geo.valid || geo.expected <= 0) geo_bits |= R_PROFILE_SPARSE;
                    else if (geo.contrast < m_bragg_contrast_min * geo.expected) geo_bits |= R_NO_BRAGG;
                    if (geo.valid && m_plateau_mip_hi > m_plateau_mip_lo && m_mip_dqdx > 0) {
                        const double pm = geo.plateau_med / m_mip_dqdx;
                        if (pm < m_plateau_mip_lo || pm > m_plateau_mip_hi) geo_bits |= R_PLATEAU_OFF_MIP;
                    }
                    std::vector<double> test, ref_mu, ref_flat;
                    for (size_t i = 0; i < live_geo.rr.size(); ++i) {
                        if (live_geo.rr[i] > m_compare_range_cm * units::cm) continue;
                        test.push_back(live_geo.dQdx[i]);
                        ref_mu.push_back(mu_fn->scalar_function(live_geo.rr[i] / units::cm + m_offset_length_cm));
                        ref_flat.push_back(m_mip_dqdx);
                    }
                    double ks_mu = 0, ks_flat = 0, ratio_mu = 0, ratio_flat = 0;
                    if (test.size() >= 3) {
                        auto sum = [](const std::vector<double>& v) { double s = 0; for (double x : v) s += x; return s; };
                        ks_mu = WireCell::kslike_compare(test, ref_mu);
                        ratio_mu = sum(ref_mu) / (sum(test) + 1e-9);
                        ks_flat = WireCell::kslike_compare(test, ref_flat);
                        ratio_flat = sum(ref_flat) / (sum(test) + 1e-9);
                        if (ks_mu + m_ks_margin >= ks_flat) geo_bits |= R_SHAPE_FLAT;
                    }
                    else {
                        geo_bits |= R_PROFILE_SPARSE;
                    }
                    SPDLOG_LOGGER_DEBUG(s_log, "{}anchor_geo_fallback: cluster {} shift {:.2f} cm peak/plateau {:.2f} | anchored contrast {:.2f}/{:.2f} ks {:.3f}/{:.3f} bits {} | geometric contrast {:.2f}/{:.2f} ks {:.3f}/{:.3f} bits {} -> {}",
                                        m_evt_tag, rec.cluster_id, rec.bragg_anchor_shift / units::cm, anchor_peak_mean / rec.bragg.plateau_med,
                                        rec.bragg.contrast, rec.bragg.expected, rec.ks_mu, rec.ks_flat, rec.reject_bits & shape_bits,
                                        geo.contrast, geo.expected, ks_mu, ks_flat, geo_bits, geo_bits == 0 ? "geometric stands" : "anchored stands");
                    if (geo_bits == 0) {
                        rec.reject_bits &= ~shape_bits;
                        rec.bragg = geo;
                        rec.ks_mu = ks_mu; rec.ks_flat = ks_flat; rec.ratio_mu = ratio_mu; rec.ratio_flat = ratio_flat;
                        rec.bragg_anchor_fallback = 1;
                    }
                }
            }
        }

        // ---- doc pdvd/66 (T8): the profile's end geometry and the fit's charge support
        // Computed on the GEOMETRIC profile (every fit row, dead or live) -- the
        // offline class-G rule (census_score.py) reads the same rows -- and on every
        // fitted segment of the main cluster against the chain's plateau median
        // (class H).  Writer fields; the guard below is the only verdict path and
        // is off by default.
        if (!prof.empty()) {
            std::vector<size_t> idx;
            for (size_t i = 0; i < prof.rr.size(); ++i) if (prof.rr[i] <= m_end_window_cm * units::cm) idx.push_back(i);
            rec.n_end_pts = static_cast<int>(idx.size());
            if (idx.size() >= 2) {
                const double arc = std::abs(prof.L[idx.back()] - prof.L[idx.front()]);
                double span = 0;
                for (size_t a = 0; a < idx.size(); ++a)
                    for (size_t b = a + 1; b < idx.size(); ++b)
                        span = std::max(span, (prof.pts[idx[a]] - prof.pts[idx[b]]).magnitude());
                rec.end_arc_cm = arc / units::cm; rec.end_span_cm = span / units::cm;
                rec.end_arc_span = span > 0.01 * units::cm ? arc / span : -1.0;
            }
            if (rec.bragg.plateau_med > 0) {
                rec.n_unsupported_segs = 0; rec.unsupported_len_cm = 0;
                for (auto& seg : pa.find_cluster_segments(g, *main)) {   // ordered_edges: deterministic
                    if (!seg) continue;
                    const double med = segment_median_dQ_dx(seg) * units::cm;   // e/cm
                    if (med <= 0) continue;
                    const double frac = med / rec.bragg.plateau_med;
                    if (chain_set.count(seg) && (rec.chain_support_min < 0 || frac < rec.chain_support_min)) rec.chain_support_min = frac;
                    const double len = segment_track_length(seg);
                    if (len >= m_unsupported_min_len_cm * units::cm && frac < m_unsupported_frac) {
                        ++rec.n_unsupported_segs; rec.unsupported_len_cm += len / units::cm;
                        if (rec.unsupported_frac_min < 0 || frac < rec.unsupported_frac_min) rec.unsupported_frac_min = frac;
                    }
                }
            }
            if (m_profile_geometry_guard &&
                (rec.end_arc_span >= m_profile_arc_span_max || rec.n_unsupported_segs > 0)) {
                rec.reject_bits |= R_PROFILE_GEOMETRY;
            }
        }

        // ---- arms: body (delta rays / hadrons) and stop (Michel / continuation)
        IndexedShowerSet showers;
        std::shared_ptr<Shower> michel_shower;
        std::vector<StmMichelArm> michel_arms;
        std::vector<SegmentPtr> other_arms;      // doc pdvd/64 (T6): kOther arms, published as role 7 at the end when publish_other_arms
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
                        if (m_publish_other_arms) other_arms.push_back(arm);   // doc pdvd/64
                    }
                }
            }
            // the stop
            if (stop_v) {
                SegmentPtr last = chain.back();
                std::vector<StmMichelArm> stop_arms;
                for (auto e : sorted_out_edges(stop_v->get_descriptor(), g)) {
                    auto arm = g[e].segment;
                    if (!arm || arm == last || chain_set.count(arm)) continue;
                    stop_arms.push_back(stm_michel_classify_stop_arm(g, last, arm, stop_v, th));
                    // doc pdvd/62/63 (T5): name every stop arm with the numbers the
                    // classifier saw, so a hand-scan item's "Bragg stub past the fit
                    // end" can be matched to the C++ kink/mip rather than an offline
                    // re-derivation.  Log only.
                    // doc pdvd/73 (P2): plus the kink over 5 cm and the far
                    // subtree walked to 100 cm with the stop fenced off, so any
                    // arm sizes the (b) cap and the (c) window from the C++'s
                    // own numbers (the payload points are not the fits), and
                    // kink_w, the (c) kink this arm was classified with (-1 off).
                    const auto& sa = stop_arms.back();
                    SPDLOG_LOGGER_DEBUG(s_log,
                        "{}CheckSTM_Michel stop-arm: cluster {} seg {} kind {} len {:.2f} cm far_len {:.2f} cm mip {:.2f} kink {:.1f} deg shower {} terminal {} kink5 {:.1f} far_full {:.2f} kink_w {:.1f}",
                        m_evt_tag, rec.cluster_id,
                        (arm->cluster() ? arm->cluster()->get_cluster_id() : 0) * 1000 + static_cast<int>(arm->get_graph_index()),
                        static_cast<int>(sa.kind), sa.len / units::cm, sa.far_len / units::cm, sa.mip, sa.kink_deg,
                        sa.shower_like ? 1 : 0, sa.terminal ? 1 : 0,
                        segment_pair_kink_deg(last, arm, stm_michel_vertex_point(stop_v), 5 * units::cm),
                        sa.terminal ? 0.0 : stm_michel_far_subtree_len(g, find_other_vertex(g, arm, stop_v), arm, stop_v,
                                                                       100 * units::cm) / units::cm,
                        sa.kink_w_deg);
                }
                bool any_michel = false;
                for (const auto& a : stop_arms) any_michel = any_michel || a.kind == StmMichelArm::kMichel;
                const bool bragg_here = rec.bragg.valid && rec.bragg.expected > 0 &&
                                        rec.bragg.contrast >= m_bragg_contrast_min * rec.bragg.expected;
                for (auto a : stop_arms) {
                    ++rec.n_stop_arms;
                    if (m_michel_guards_stop && a.kind == StmMichelArm::kContinuation && any_michel && bragg_here &&
                        a.len <= m_michel_max_len_cm * units::cm) {
                        a.kind = StmMichelArm::kOther;   // stop debris beside a Michel at a Bragg-confirmed stop
                    }
                    if (a.kind == StmMichelArm::kContinuation) {
                        rec.reject_bits |= R_CONTINUATION;
                        if (a.len > rec.cont_len) { rec.cont_len = a.len; rec.cont_angle_deg = a.kink_deg; rec.cont_mip = a.mip; }
                    }
                    else if (a.kind == StmMichelArm::kMichel) {
                        michel_arms.push_back(a);
                    }
                    else {
                        // doc pdvd/64 (T6): a kOther stop arm (hot stub not absorbed,
                        // debris beside a Michel, a continuation demoted above) had
                        // no counter and no rows through doc pdvd/63.
                        ++rec.n_stop_other;
                        if (m_publish_other_arms) other_arms.push_back(a.seg);
                    }
                }
            }
        }

        // ---- the Michel: ONE object (doc pdhd/15) ---------------------------
        // A Michel is a track plus whatever the 3-D clustering broke off it.
        // The object is assembled first -- the stop arm and everything the
        // shower walk reaches inside the main cluster, then every fitted
        // segment of an admitted companion near the stop -- and energised
        // ONCE, after the last piece is in.
        //
        // Through doc pdhd/14 the attached path called calculate_kinematics
        // BEFORE the dot loop and only re-ran it for conn_type 2, so a dot the
        // loop added to the very same shower never reached the energy: PDVD
        // 039252_15 cluster 91 reported 16.2 MeV for a 29.2 MeV Michel whose
        // second piece sits 4.33 cm past the arm's tip at 20 deg, and the same
        // 16 MeV is what mc.json's mu- -> e- node shows.  Measured on the d14
        // arms: 9 PDHD / 17 PDVD candidates, missing fraction median 22.5 % /
        // 12.2 %, max 62 % / 65 %.
        std::vector<std::pair<SegmentPtr, double>> michel_pieces;   // (segment, distance to the stop)
        if (!michel_arms.empty() && stop_v) {
            std::sort(michel_arms.begin(), michel_arms.end(),
                      [](const StmMichelArm& a, const StmMichelArm& b) {
                          if (a.len != b.len) return a.len > b.len;
                          return a.seg->get_graph_index() < b.seg->get_graph_index();
                      });
            const auto& seed = michel_arms.front();
            rec.michel_len = seed.len; rec.michel_mip = seed.mip; rec.michel_kink_deg = seed.kink_deg; rec.michel_far_len = seed.far_len;
            rec.michel_conn_type = 1;
            rec.michel_dis_cm = 0;              // attached: the arm leaves the stop vertex itself
            // doc pdhd/14: the daughter's id, in the same cluster*1000 + graph
            // index encoding as stop_vtx_id and T_stm_michel_pts.seg_id.
            if (seed.seg)
                rec.michel_seg_id = rec.cluster_id * 1000 + static_cast<int>(seed.seg->get_graph_index());
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
                for (auto& sg : ms) { if (!chain_set.count(sg)) set_pdg(sg, 11); }
                // doc pdvd/53: those members are PDG 11 and belong to the Michel
                // object, but through doc pdvd/51 only `michel_arms` got point
                // rows -- so a real member such as 039253_14 cluster 49's 49004
                // carried PDG 11 and NO role, and a display grouping off roles
                // put it in "unassigned".  Behind survey_enable so the knob-off
                // stm_michel_pts is unchanged.
                if (m_survey_enable) {
                    std::vector<SegmentPtr> extra_members;
                    for (auto& sg : ms) {
                        if (!sg || chain_set.count(sg)) continue;
                        const int sid = (sg->cluster() ? sg->cluster()->get_cluster_id() : 0) * 1000
                                      + static_cast<int>(sg->get_graph_index());
                        if (rec.claimed.count(sid)) continue;
                        extra_members.push_back(sg);
                    }
                    // IndexedSegmentSet is index-ordered, but sort on the id we
                    // actually write so the row order cannot depend on it.
                    std::sort(extra_members.begin(), extra_members.end(),
                              [](const SegmentPtr& a, const SegmentPtr& b) {
                                  const int ca = a->cluster() ? a->cluster()->get_cluster_id() : 0;
                                  const int cb = b->cluster() ? b->cluster()->get_cluster_id() : 0;
                                  if (ca != cb) return ca < cb;
                                  return a->get_graph_index() < b->get_graph_index();
                              });
                    for (auto& sg : extra_members) add_points(rec, sg, 3);
                }
                michel_shower->set_particle_type(11);
                // The CORE alone, measured before any companion piece joins:
                // michel_ke_core is exactly what michel_ke_dqdx meant through
                // doc pdhd/14, and michel_ke_range is the TRACK's range energy
                // (the object call below zeroes range whenever a member is
                // graph-disconnected, PRShower.cxx:1855-1859).
                michel_shower->calculate_kinematics(particle_data(), m_recomb_model);
                rec.michel_ke_core = michel_shower->get_kine_dQdx() / units::MeV;
                rec.michel_ke_range = michel_shower->get_kine_range() / units::MeV;
            }
            else {
                rec.michel_ke_core = segment_cal_kine_dQdx(seed.seg, m_recomb_model) / units::MeV;
            }
            for (auto& a : michel_arms) if (a.seg) michel_pieces.emplace_back(a.seg, 0.0);
        }

        // ---- the companion pieces: near the stop, closer to the stop than to
        // the muon body ---------------------------------------------------------
        // doc pdhd/15 sec 3: michel_dot_radius_cm is a CLUSTER admission test
        // (applied at :871 on the cluster's closest approach), not a per-segment
        // one.  It used to be re-applied to every segment here, which clipped
        // 39.8 % of the charge of PDVD 039252_15 cluster 77's Michel -- that
        // cluster spans 1.91 -> 19.36 cm from the stop, so more than a third of
        // the object sat outside a radius its closest point had already passed.
        // A segment of an admitted cluster is bounded by radius + cluster length
        // by construction, so the cluster test is the whole bound.  The
        // body-exclusion test and the per-piece length cap stay.
        // doc pdvd/53: the survey ledger.  Every gate below that drops an
        // admitted companion records WHY here, so a role-6 row can say it.  The
        // first gate to fire wins -- note_cl/note_seg never overwrite.  Codes
        // are the doc pdvd/53 table; 9 = the survey reached it and neither stage
        // was ever offered it; 10 = there was no stop vertex, so every companion
        // was fitted and nothing examined it.
        std::map<int, int> sv_cl_rej;                       // cluster id -> gate
        std::map<int, double> sv_cl_dstop;                  // cluster id -> d_stop (cluster level)
        std::map<int, double> sv_cl_dbody;                  // cluster id -> d_body (cluster level, gate 6)
        std::map<int, int> sv_seg_rej;                      // seg_id -> gate
        std::map<int, std::pair<double, double>> sv_seg_d;  // seg_id -> (d_stop, d_body)
        // Codes 1 and 9 are PROVISIONAL: "the Michel stage was not interested"
        // and "nobody has looked yet" are both states a later stage overwrites
        // with a real verdict.  Without this the Michel loop's code 1 masks the
        // gamma stage's answer on exactly the clusters the owner asks about --
        // 039253_13 cluster 431 read `1` when the interesting answer is that the
        // gamma body exclusion rejected it by 0.03 cm.  Every other code is
        // final: the first stage to reject for a real reason owns the row.
        auto note_cl = [&](int cid, int code, double d) {
            if (!m_survey_enable) return;
            auto it = sv_cl_rej.find(cid);
            if (it != sv_cl_rej.end() && it->second != 1 && it->second != 9) return;
            sv_cl_rej[cid] = code; sv_cl_dstop[cid] = d;
        };
        auto note_seg = [&](int sid, int code, double d_stop, double d_body) {
            if (!m_survey_enable) return;
            if (sv_seg_rej.count(sid)) return;
            sv_seg_rej[sid] = code; sv_seg_d[sid] = {d_stop, d_body};
        };
        (void)note_cl; (void)note_seg;

        if (stop_v && (!companions.empty() || m_stop_local_michel_pieces)) {
            // Two passes: collect every admissible piece with its distance to the
            // stop, then assemble.  The seed of a BRIDGED object is then the
            // NEAREST piece rather than whichever the graph's edge order reached
            // first, and michel_seg_id names the piece the gap is measured to.
            struct Piece { SegmentPtr seg; double d_stop; int cluster_id, gidx; };
            std::vector<Piece> pieces;
            double d_unfit = 1e9;            // closest approach of an UNFITTED companion
            // doc pdvd/62 (T3b): a DISCONNECTED piece of the main cluster --
            // what a pr54-kept residual is (its own two vertices, no edge to
            // the chain) -- is admitted with the same three gates a companion
            // piece gets below (radius, length, body test).  "Disconnected"
            // is: neither endpoint vertex belongs to a chain segment and
            // neither has an out-edge to one; an attached interior arm fails
            // that by construction and stays T6's population.
            if (m_stop_local_michel_pieces && !chain.empty()) {
                std::set<VertexPtr> chain_vset(chain_vtxs.begin(), chain_vtxs.end());   // membership only, never iterated
                for (auto& seg : pa.find_cluster_segments(g, *main)) {   // ordered_edges: deterministic
                    if (!seg || chain_set.count(seg)) continue;
                    const int sid = main->get_cluster_id() * 1000 + static_cast<int>(seg->get_graph_index());
                    if (rec.claimed.count(sid)) continue;
                    auto [va, vb] = find_vertices(g, seg);
                    if (!va || !vb) continue;
                    bool touches_chain = chain_vset.count(va) || chain_vset.count(vb);
                    for (VertexPtr v : {va, vb}) {
                        if (touches_chain) break;
                        for (auto e : sorted_out_edges(v->get_descriptor(), g)) {
                            auto s2 = g[e].segment;
                            if (s2 && chain_set.count(s2)) { touches_chain = true; break; }
                        }
                    }
                    if (touches_chain) continue;
                    auto [d_stop, cp] = segment_get_closest_point(seg, rec.stop_pt, "fit", "main");
                    if (d_stop > m_michel_dot_radius_cm * units::cm) continue;
                    if (segment_track_length(seg) > m_dot_max_len_cm * units::cm) continue;
                    double d_body = 1e9;
                    for (size_t i = 0; i < prof.pts.size(); ++i) {
                        if (prof.rr[i] < m_dot_body_exclusion_cm * units::cm) continue;
                        d_body = std::min(d_body, (prof.pts[i] - cp).magnitude());
                    }
                    if (d_body < d_stop) continue;    // a delta ray / body fragment, not a Michel piece
                    ++rec.n_local_pieces;
                    pieces.push_back({seg, d_stop, main->get_cluster_id(), static_cast<int>(seg->get_graph_index())});
                }
            }
            for (auto* oc : companions) {
                // The radius test again, at CLUSTER level, against the stop the
                // chain actually ended on.  `companions` was selected against
                // the TAGGER's stop (:880), and the two are not the same point:
                // R_STOP_UNMATCHED walks the chain to the farthest vertex of the
                // main cluster (:922-929) and stop_extend_max walks it along a
                // continuation, so the two can be hundreds of cm apart -- PDHD
                // 029107_17 cluster 33 has stop_dis 266.3 cm.  The per-SEGMENT
                // radius test that used to stand below hid this; dropping it
                // (sec 3) exposed it, and the cluster-level test is the right
                // place for it because the object is admitted whole.
                // doc pdvd/53: a cluster the SURVEY admitted only because it is
                // longer than companion_max_len_cm is not a Michel candidate.
                // With survey_max_len_cm at its default this cannot fire; the
                // guard is here so "the survey changes no selection" is true by
                // construction rather than by a coincidence of two defaults.
                const auto [ccp, cblob] = oc->get_closest_point_blob(rec.stop_pt);
                const double d_cl = (ccp - rec.stop_pt).magnitude();
                if (oc->get_length() > m_companion_max_len_cm * units::cm) {
                    note_cl(oc->get_cluster_id(), 3, d_cl); continue;
                }
                if (d_cl > m_michel_dot_radius_cm * units::cm) {
                    note_cl(oc->get_cluster_id(), 1, d_cl); continue;
                }
                auto segs = pa.find_cluster_segments(g, *oc);   // ordered_edges: deterministic
                if (segs.empty()) {
                    ++rec.n_dot_clusters_unfit;
                    double q = 0; for (const auto* b : oc->children()) q += b->charge();
                    rec.dots_charge_unfit += q;
                    rec.unfit_dot_clusters.push_back(oc);   // doc pdvd/81: bookkeeping only
                    d_unfit = std::min(d_unfit, d_cl);
                    continue;
                }
                for (auto& seg : segs) {
                    auto [d_stop, cp] = segment_get_closest_point(seg, rec.stop_pt, "fit", "main");
                    // segment_get_closest_point returns a 1e9 sentinel when the
                    // segment carries no usable point, and the body test below
                    // cannot catch it: on a chain shorter than
                    // dot_body_exclusion_cm the body distance is the same
                    // sentinel, so `d_body < d_stop` is false and the piece is
                    // kept with michel_dis_cm = 1e8 cm (PDHD 028084_12 cluster 3,
                    // a 1.0 cm 3-point chain).  A piece of an admitted cluster
                    // can be at most radius + the cluster cap from the stop.
                    const int sid = oc->get_cluster_id() * 1000 + static_cast<int>(seg->get_graph_index());
                    if (d_stop > (m_michel_dot_radius_cm + m_companion_max_len_cm) * units::cm) {
                        note_seg(sid, 2, d_stop, -1); continue;
                    }
                    if (segment_track_length(seg) > m_dot_max_len_cm * units::cm) {
                        note_seg(sid, 3, d_stop, -1); continue;
                    }
                    // distance to the muon body beyond its last dot_body_exclusion
                    double d_body = 1e9;
                    for (size_t i = 0; i < prof.pts.size(); ++i) {
                        if (prof.rr[i] < m_dot_body_exclusion_cm * units::cm) continue;
                        d_body = std::min(d_body, (prof.pts[i] - cp).magnitude());
                    }
                    if (d_body < d_stop) {           // a delta ray / body fragment, not a Michel piece
                        note_seg(sid, 4, d_stop, d_body); continue;
                    }
                    pieces.push_back({seg, d_stop, oc->get_cluster_id(),
                                      static_cast<int>(seg->get_graph_index())});
                }
            }
            std::sort(pieces.begin(), pieces.end(), [](const Piece& a, const Piece& b) {
                if (a.d_stop != b.d_stop) return a.d_stop < b.d_stop;
                if (a.cluster_id != b.cluster_id) return a.cluster_id < b.cluster_id;
                return a.gidx < b.gidx;
            });
            for (const auto& pc : pieces) {
                ++rec.n_dots;
                rec.dots_ke_dqdx += segment_cal_kine_dQdx(pc.seg, m_recomb_model) / units::MeV;
                set_pdg(pc.seg, 11);
                // doc pdvd/51: role 3 = "a member of the Michel OBJECT", for
                // every connection type.  Through doc pdhd/17 role 3 was set on
                // the ATTACHED arms only (:1390) and every bridged piece got
                // role 4, so the hand-scan display drew every bridged Michel in
                // the `dots` colour and named it a dot -- which is exactly what
                // the owner saw on 039252_15 cluster 77.  Role 4 stays as the
                // residual bucket for a fitted piece the object does not absorb
                // (today: none -- every admitted piece joins the object).
                add_points(rec, pc.seg, 3);
                michel_pieces.emplace_back(pc.seg, pc.d_stop);
                if (rec.michel_conn_type == 0) {
                    // Nothing attached survived at the stop: the 3-D clustering
                    // detached the Michel.  The nearest piece is the seed and the
                    // object is BRIDGED, not attached.  Set outside the
                    // build_michel_shower guard -- through doc pdhd/14 this lived
                    // inside it, so with that knob off a real Michel reported
                    // conn_type 0 with n_dots > 0.
                    rec.michel_conn_type = 2;
                    rec.michel_dis_cm = pc.d_stop / units::cm;
                    rec.michel_seg_id = pc.cluster_id * 1000 + pc.gidx;
                    rec.michel_len = segment_track_length(pc.seg);
                    if (m_build_michel_shower) {
                        michel_shower = std::make_shared<Shower>(g);
                        michel_shower->set_start_vertex(stop_v, 2);
                        michel_shower->set_start_segment(pc.seg, false, "fit", "associate_points");
                        michel_shower->set_particle_type(11);
                        michel_shower->calculate_kinematics(particle_data(), m_recomb_model);
                        rec.michel_ke_core = michel_shower->get_kine_dQdx() / units::MeV;
                        rec.michel_ke_range = michel_shower->get_kine_range() / units::MeV;
                    }
                    else {
                        rec.michel_ke_core = segment_cal_kine_dQdx(pc.seg, m_recomb_model) / units::MeV;
                    }
                }
                else if (michel_shower) {
                    michel_shower->add_segment(pc.seg, true, "fit", "associate_points");
                }
            }
            // doc pdhd/15 sec 5: a companion that passed every admission test and
            // that the fitter produced no segment for is still the muon's
            // "additional activity" -- it just has no dx, so only its charge can
            // be read.  Through doc pdhd/14 it landed in n_dot_clusters_unfit and
            // nowhere else, leaving michel_found 0 beside a non-zero energy on 15
            // of 302 PDHD candidates.  conn type 3 = CHARGE ONLY: no fitted
            // segment, no shower, michel_seg_id stays -1, and the whole energy is
            // the charge conversion.
            if (rec.michel_conn_type == 0 && rec.n_dot_clusters_unfit > 0 && d_unfit < 1e8) {
                rec.michel_conn_type = 3;
                rec.michel_dis_cm = d_unfit / units::cm;
            }
        }


        // ---- doc pdvd/83 (doc 78 action item 4): a Michel that leaves the muon
        // chain a few cm BEFORE its stop.  The four missed Michels of doc 78
        // sec 3.2 (039349_64/24 s24003, 039349_9/19 s19003, 039349_64/52
        // s52018, 039349_69/56 s56006) are attached to the chain (doc 80: rej
        // 13) at its penultimate vertex, 2.5-6.5 cm before the stop: the
        // stop-arm loop above never sees them and the interior loop files them
        // as kOther, which nothing reads.  Only where NO Michel exists yet
        // (conn 0 after the attached and companion stages), so no found
        // Michel's energy or connection type can change.
        //
        // Design constraint, the same one T2c states for itself: this writes
        // the Michel object and NOTHING ELSE a verdict reads.  A kContinuation
        // answer from the classifier is ignored (never OR'd into
        // R_CONTINUATION); n_stop_arms / n_stop_other / n_body_other and the
        // michel_guards_stop demotion above are untouched.  The one declared
        // channel through which is_stm can move is topology_stop_evidence,
        // which reads michel_found like any other attached Michel.
        if (m_michel_near_stop_arm_cm > 0 && stop_v && chain.size() >= 2 && rec.michel_conn_type == 0) {
            const int mcid = main->get_cluster_id();
            auto near = stm_michel_near_stop_arms(g, chain, chain_vtxs, m_michel_near_stop_arm_cm * units::cm, th,
                [&](const SegmentPtr& sg) {
                    return chain_set.count(sg) > 0 ||
                           rec.claimed.count(mcid * 1000 + static_cast<int>(sg->get_graph_index())) > 0;
                });
            rec.n_near_arms_examined = near.n_examined;
            for (const auto& a : near.michel) {
                SPDLOG_LOGGER_DEBUG(s_log,
                    "{}CheckSTM_Michel near-arm: cluster {} seg {} vtx_index {} dist {:.2f} cm len {:.2f} cm far_len {:.2f} cm mip {:.2f} kink {:.1f} deg shower {} terminal {}",
                    m_evt_tag, rec.cluster_id, mcid * 1000 + static_cast<int>(a.seg->get_graph_index()), near.vtx_index,
                    near.dist / units::cm, a.len / units::cm, a.far_len / units::cm, a.mip, a.kink_deg,
                    a.shower_like ? 1 : 0, a.terminal ? 1 : 0);
            }
            if (near.vtx_index > 0 && !near.michel.empty()) {
                VertexPtr nv = chain_vtxs[near.vtx_index];
                const auto& seed = near.michel.front();
                rec.michel_near_arm = 1;
                rec.near_arm_dist_cm = near.dist / units::cm;
                rec.michel_len = seed.len; rec.michel_mip = seed.mip; rec.michel_kink_deg = seed.kink_deg; rec.michel_far_len = seed.far_len;
                rec.michel_conn_type = 1;                       // attached: graph-connected to the chain
                rec.michel_dis_cm = near.dist / units::cm;      // the one conn-1 Michel whose dis is not 0 (see T2c below)
                rec.michel_seg_id = mcid * 1000 + static_cast<int>(seed.seg->get_graph_index());
                for (const auto& a : near.michel) { set_pdg(a.seg, 11); add_points(rec, a.seg, 3); }
                // Duplicated from the attached path above (M10), the start
                // vertex being the arm's own chain vertex, not the stop.
                if (m_build_michel_shower) {
                    michel_shower = std::make_shared<Shower>(g);
                    michel_shower->set_start_vertex(nv, 1);
                    michel_shower->set_start_segment(seed.seg, false, "fit", "associate_points");
                    IndexedSegmentSet used(chain_set);
                    michel_shower->complete_structure_with_start_segment(used, "fit", "associate_points", true);
                    IndexedVertexSet mv; IndexedSegmentSet ms;
                    michel_shower->fill_sets(mv, ms, false);
                    for (auto& sg : ms) { if (!chain_set.count(sg)) set_pdg(sg, 11); }
                    if (m_survey_enable) {
                        std::vector<SegmentPtr> extra_members;
                        for (auto& sg : ms) {
                            if (!sg || chain_set.count(sg)) continue;
                            const int sid = (sg->cluster() ? sg->cluster()->get_cluster_id() : 0) * 1000
                                          + static_cast<int>(sg->get_graph_index());
                            if (rec.claimed.count(sid)) continue;
                            extra_members.push_back(sg);
                        }
                        std::sort(extra_members.begin(), extra_members.end(),
                                  [](const SegmentPtr& a, const SegmentPtr& b) {
                                      const int ca = a->cluster() ? a->cluster()->get_cluster_id() : 0;
                                      const int cb = b->cluster() ? b->cluster()->get_cluster_id() : 0;
                                      if (ca != cb) return ca < cb;
                                      return a->get_graph_index() < b->get_graph_index();
                                  });
                        for (auto& sg : extra_members) add_points(rec, sg, 3);
                    }
                    michel_shower->set_particle_type(11);
                    michel_shower->calculate_kinematics(particle_data(), m_recomb_model);
                    rec.michel_ke_core = michel_shower->get_kine_dQdx() / units::MeV;
                    rec.michel_ke_range = michel_shower->get_kine_range() / units::MeV;
                }
                else {
                    rec.michel_ke_core = segment_cal_kine_dQdx(seed.seg, m_recomb_model) / units::MeV;
                }
                for (const auto& a : near.michel) if (a.seg) michel_pieces.emplace_back(a.seg, near.dist);
            }
        }


        // ---- doc pdvd/51: the muon-capture gammas at the stop -------------------
        // The physics.  A mu- that ranges out in argon is captured by a nucleus
        // far more often than it decays (in LAr, capture dominates), and the
        // capture de-excitation emits gammas of order an MeV.  A gamma is
        // neutral: it leaves no track from the stop, travels a couple of
        // Compton mean free paths, and deposits a compact blob at an arbitrary
        // angle.  So the object the reconstruction can see is "a small,
        // detached, same-bundle blob past the stop", and the muon -> gamma edge
        // is a CLAIM about a neutral, not a reconstructed connection -- which is
        // precisely what the renderer's pseudo-carrier expresses
        // (MultiAlgBlobClustering.cxx:2154, selected by start_connection_type 2).
        //
        // Measured before it was built (doc pdvd/51 sec 4, d16vnu, 119 events,
        // 151 is_stm candidates): the stop-anchored density of compact
        // same-bundle blobs falls by a factor ~4000 from the 0-10 cm shell to
        // the 100-150 cm one, while the SAME predicate on foreign-bundle
        // companions at the SAME anchor is flat at ~0.6 per candidate per
        // 1e6 cm^3 -- so the population belongs to the stop rather than filling
        // the volume.  And it is ~2x commoner on candidates with NO Michel,
        // which is what the mechanism predicts (capture gives gammas and no
        // Michel; decay gives a Michel and no capture gammas).  That
        // anti-correlation owes nothing to any threshold here, and it is what
        // sec 6.5's body-exclusion defect inverted.
        //
        // NOT folded into the Michel.  One object per companion CLUSTER, its own
        // Shower, its own branches, role 5 in the point cloud.
        std::vector<std::pair<ShowerPtr, double>> gamma_showers;   // (shower, KE MeV)
        // doc pdvd/85: what the withhold step needs to undo, filled only with
        // stop_gamma_require_stm on -- each accepted segment, its distance to
        // the stop, and the particle info set_pdg is about to overwrite.
        struct GWithhold { SegmentPtr seg; double d_stop; std::shared_ptr<Aux::ParticleInfo> info; double score; };
        std::vector<GWithhold> sg_withhold;
        if (m_stop_gamma_enable && stop_v && !companions.empty()) {
            // Clusters the Michel object already owns are off limits.  The ring
            // makes this near-vacuous -- the Michel's cluster test is the ring's
            // inner edge, so a Michel cluster cannot be in the ring -- but it is
            // stated rather than relied on, because the two radii are knobs.
            std::set<int> michel_cluster_ids;
            for (const auto& [sg, d] : michel_pieces) {
                (void)d;
                if (sg && sg->cluster()) michel_cluster_ids.insert(sg->cluster()->get_cluster_id());
            }

            // The body set, snapshotted before any gamma is accepted (see the
            // body-exclusion comment below).  rec.px/py/pz hold every point
            // add_points() has taken so far: role 1 the muon chain, 2 deltas,
            // 3 the Michel object, 4 a dot the object did not absorb.
            std::vector<Point> body_pts;
            body_pts.reserve(rec.px.size());
            for (size_t bi = 0; bi < rec.px.size(); ++bi) {
                const Point bp(rec.px[bi], rec.py[bi], rec.pz[bi]);
                if ((bp - rec.stop_pt).magnitude() < m_dot_body_exclusion_cm * units::cm) continue;
                body_pts.push_back(bp);
            }

            struct GObj {
                Cluster* oc{nullptr};
                std::vector<std::pair<SegmentPtr, double>> segs;   // (segment, d to stop)
                double d_stop{1e9}, ke{0}, charge{0};
                int cluster_id{-1}, gidx{-1};
            };
            std::vector<GObj> objs;
            for (auto* oc : companions) {
                const auto [ccp, cblob] = oc->get_closest_point_blob(rec.stop_pt);
                const double d_cl = (ccp - rec.stop_pt).magnitude();
                if (michel_cluster_ids.count(oc->get_cluster_id())) {
                    note_cl(oc->get_cluster_id(), 7, d_cl); continue;
                }
                // THE RING, plus compactness.  Inside michel_dot_radius_cm the
                // charge is the Michel's to claim and this stage must not touch
                // it; outside stop_gamma_radius_cm the stop-anchored excess has
                // died.  A gamma deposit is a blob, so a 30 cm object past a
                // stop is another cosmic in the same bundle: 039252_0 cluster 77
                // has six same-bundle neighbours at 26.5 .. 168.9 cm and only
                // the first is a candidate -- "same Q-L bundle" alone is not a
                // discriminator, the ring and the cap are.  The predicate lives
                // in StmMichelFunctions so its boundaries are doctested.
                if (!stm_michel_stop_gamma_ring(d_cl, oc->get_length(),
                                                m_michel_dot_radius_cm * units::cm,
                                                m_stop_gamma_radius_cm * units::cm,
                                                m_stop_gamma_max_len_cm * units::cm)) {
                    // 9 says "the survey brought it and neither stage was ever
                    // offered it"; 5 says "the gamma stage looked and said no".
                    note_cl(oc->get_cluster_id(),
                            d_cl > m_stop_gamma_radius_cm * units::cm ? 9 : 5, d_cl);
                    continue;
                }
                // THE BODY EXCLUSION -- at CLUSTER level, and against
                // EVERYTHING THE CHAIN HAS ALREADY CLAIMED, not just the muon.
                //
                // Both halves of that sentence are measured, and each was got
                // wrong once (doc pdvd/51 sec 6.5).  With EITHER mistake -- per
                // fitted SEGMENT instead of per cluster, or against the muon
                // profile alone instead of the whole claimed object -- the mu-
                // signature, the one physics test this class has, reads
                // 0.44 on PDVD d51gv; with both fixed it reads 1.64, against
                // 2.05 for the same predicate measured offline.  The mechanism
                // is that 29 of the 30 spurious admissions land on candidates
                // that ALSO have a Michel: they are Michel SATELLITES.  A
                // Michel showers, and its outlying blobs sit past the stop, in
                // the bundle, compact, and pass every other test -- and the
                // muon profile cannot see them, because an ATTACHED Michel's
                // arms are not in it.
                //
                // `body` is therefore every point already attributed to this
                // candidate -- muon chain, deltas, and the whole Michel object
                // (roles 1-4), snapshotted BEFORE this loop so one accepted
                // gamma cannot exclude the next -- minus the
                // dot_body_exclusion_cm around the stop itself, which is the
                // one place a real daughter is allowed to start.
                {
                    double d_body_cl = 1e9;
                    for (const auto& bp0 : body_pts) {
                        const auto [bp, bblob] = oc->get_closest_point_blob(bp0);
                        d_body_cl = std::min(d_body_cl, (bp - bp0).magnitude());
                        // doc pdvd/53: the early exit is CORRECT for the verdict
                        // and WRONG for the number.  It stops at the first body
                        // point that undercuts the stop, so d_body_cl is then an
                        // upper bound on the true minimum, not the minimum --
                        // 039252_8 cluster 97's segment 466032 came out at 27.42
                        // against a true 0.41 cm, i.e. a blob sitting ON the muon
                        // read as "rejected by 0.02 cm".  With the survey on that
                        // number is shown to a scanner, so run the loop out.
                        // The verdict cannot change: finishing can only LOWER
                        // d_body_cl, and it is already below d_cl.
                        if (d_body_cl < d_cl && !m_survey_enable) break;
                    }
                    if (d_body_cl < d_cl) {
                        note_cl(oc->get_cluster_id(), 6, d_cl);
                        // Only when 6 actually took the row, so d_body cannot
                        // describe a gate that did not fire.
                        if (m_survey_enable && sv_cl_rej[oc->get_cluster_id()] == 6)
                            sv_cl_dbody[oc->get_cluster_id()] = d_body_cl;
                        continue;
                    }
                }

                GObj o; o.oc = oc; o.d_stop = d_cl; o.cluster_id = oc->get_cluster_id();
                auto segs = pa.find_cluster_segments(g, *oc);   // ordered_edges: deterministic
                if (segs.empty()) {
                    // No dx, so no dQ/dx and no segment -- and without a segment
                    // there is nothing a Shower can hold (PR::Shower is a view
                    // over graph edges; set_start_segment requires a valid
                    // descriptor).  It is counted and its charge is converted,
                    // and it produces NO particle-flow node.  PDVD measured zero
                    // of these in 579 candidates; PDHD is where this path lives.
                    double q = 0; for (const auto* b : oc->children()) q += b->charge();
                    ++rec.stop_gamma_n_unfit;
                    rec.stop_gamma_charge += q;
                    rec.stop_gamma_ke_tot +=
                        (m_michel_unfit_from_model
                             ? stm_michel_charge_to_energy_model(q, m_recomb_model, m_michel_unfit_dedx)
                             : stm_michel_charge_to_energy(q, m_michel_unfit_recom,
                                                           m_michel_unfit_fudge, m_michel_unfit_w_ev))
                        / units::MeV;
                    if (rec.stop_gamma_dis_min < 0 || d_cl / units::cm < rec.stop_gamma_dis_min)
                        rec.stop_gamma_dis_min = d_cl / units::cm;
                    if (d_cl / units::cm > rec.stop_gamma_dis_max)
                        rec.stop_gamma_dis_max = d_cl / units::cm;
                    continue;
                }
                for (auto& seg : segs) {
                    auto [d_stop, cp] = segment_get_closest_point(seg, rec.stop_pt, "fit", "main");
                    const int sid = oc->get_cluster_id() * 1000 + static_cast<int>(seg->get_graph_index());
                    if (d_stop > (m_stop_gamma_radius_cm + m_stop_gamma_max_len_cm) * units::cm) {
                        note_seg(sid, 2, d_stop, -1); continue;
                    }
                    // The body exclusion, the same test the Michel pieces face
                    // (:1520): a piece closer to the muon BODY than to the stop
                    // is a delta ray thrown along the track, not something the
                    // stop emitted.
                    double d_body = 1e9;
                    for (size_t i = 0; i < prof.pts.size(); ++i) {
                        if (prof.rr[i] < m_dot_body_exclusion_cm * units::cm) continue;
                        d_body = std::min(d_body, (prof.pts[i] - cp).magnitude());
                    }
                    if (d_body < d_stop) { note_seg(sid, 6, d_stop, d_body); continue; }
                    o.segs.emplace_back(seg, d_stop);
                    o.ke += segment_cal_kine_dQdx(seg, m_recomb_model) / units::MeV;
                }
                if (o.segs.empty()) continue;
                std::sort(o.segs.begin(), o.segs.end(),
                          [](const auto& a, const auto& b) {
                              if (a.second != b.second) return a.second < b.second;
                              return a.first->get_graph_index() < b.first->get_graph_index();
                          });
                o.gidx = static_cast<int>(o.segs.front().first->get_graph_index());
                o.d_stop = o.segs.front().second;
                objs.push_back(std::move(o));
            }
            // ACCEPTANCE is a separate stage from admission: the energy only
            // exists after the fitter has run, so it cannot gate the cluster
            // loop above.  A candidate rejected here is dropped from the object
            // list, not from the fitter -- its charge stays in the 2-D maps.
            std::sort(objs.begin(), objs.end(), [](const GObj& a, const GObj& b) {
                if (a.d_stop != b.d_stop) return a.d_stop < b.d_stop;
                return a.cluster_id < b.cluster_id;
            });
            for (const auto& o : objs) {
                if (rec.n_stop_gammas >= m_stop_gamma_max_n) {
                    // The knob-off path keeps the original `break` exactly.  With
                    // the survey on we walk the rest only to record that they were
                    // capped out; n_stop_gammas cannot fall, so the accepted set
                    // is the same either way.
                    if (!m_survey_enable) break;
                    note_cl(o.cluster_id, 8, o.d_stop);
                    continue;
                }
                if (!stm_michel_stop_gamma_energy(o.ke, m_stop_gamma_min_ke_mev,
                                                  m_stop_gamma_max_ke_mev)) {
                    note_cl(o.cluster_id, 8, o.d_stop); continue;
                }
                ++rec.n_stop_gammas;
                rec.stop_gamma_ke_tot += o.ke;
                rec.stop_gamma_ke_max = std::max(rec.stop_gamma_ke_max, o.ke);
                for (const auto* b : o.oc->children()) rec.stop_gamma_charge += b->charge();
                const double d_cm = o.d_stop / units::cm;
                if (rec.stop_gamma_dis_min < 0 || d_cm < rec.stop_gamma_dis_min) rec.stop_gamma_dis_min = d_cm;
                if (d_cm > rec.stop_gamma_dis_max) rec.stop_gamma_dis_max = d_cm;
                if (rec.stop_gamma_seg_id < 0) rec.stop_gamma_seg_id = o.cluster_id * 1000 + o.gidx;
                if (m_stop_gamma_require_stm)
                    for (const auto& [sg, d] : o.segs)
                        if (sg) sg_withhold.push_back({sg, d, sg->particle_info(), sg->particle_score()});
                for (const auto& [sg, d] : o.segs) { (void)d; set_pdg(sg, 11); add_points(rec, sg, 5); }
                if (!m_build_michel_shower) continue;
                // The particle-flow node.  set_start_vertex(stop_v, 2) is the
                // whole stitch: the gamma's segments live in a companion
                // cluster with no graph edge into the muon's component, so the
                // track BFS could never reach them at any knob setting, but
                // stop_v is BFS-reachable and the renderer hangs the shower off
                // the last muon segment through it.  particle_type stays 11 --
                // PDG 22 is never stored anywhere in this codebase
                // (get_particle_mass(22) is 0 and cal_kine_range falls back to
                // the MUON range function); the `gamma` node is synthesised by
                // the renderer from the connection type, and the e- leaf under
                // it is what was actually reconstructed: the conversion.
                auto gsh = std::make_shared<Shower>(g);
                gsh->set_start_vertex(stop_v, 2);
                gsh->set_start_segment(o.segs.front().first, false, "fit", "associate_points");
                gsh->set_particle_type(11);
                for (size_t i = 1; i < o.segs.size(); ++i)
                    gsh->add_segment(o.segs[i].first, true, "fit", "associate_points");
                gamma_showers.emplace_back(gsh, o.ke);
            }
        }

        // ---- doc pdvd/71 (P4): the Michel's isolated gamma blobs, phase 1 ------
        // The owner (2026-09-10): isolated gamma blobs near the stop -- the
        // Michel electron's brems and Compton pieces -- belong to the Michel;
        // look along the electron's direction, take the dots nearby, and be
        // careful of energy, because over-clustering can hand it a huge one.
        // Here, AFTER the capture-gamma stage (its body snapshot and its role-5
        // claims are untouched) and BEFORE the survey (which then leaves these
        // segments alone), every companion cluster that no stage claimed is
        // measured and gated; the survivors are only RESERVED.  Rows are written
        // in phase 2, once michel_found and michel_ke_best are final, because
        // the total-energy guard needs the final core.  Geometry comes from the
        // dx > 0 fit points -- exactly the points a role-4 row will carry -- so
        // an accepted blob's numbers re-derive offline from the payload.
        struct MGam { std::vector<SegmentPtr> segs; double d_stop{0}, ke{0}; int cluster_id{-1}, gidx{-1}; };
        std::vector<MGam> mgam;
        if (m_michel_gamma_collect && stop_v && !companions.empty() &&
            (rec.michel_conn_type == 1 || rec.michel_conn_type == 2)) {
            // The Michel object's rows so far (role 3: arms + pieces, plus the
            // shower walk's members when the survey is on) set the direction;
            // the muon beyond dot_body_exclusion_cm and the deltas are the body.
            std::vector<Point> mrows, brows;
            for (size_t i = 0; i < rec.px.size(); ++i) {
                const Point p(rec.px[i], rec.py[i], rec.pz[i]);
                if (rec.prole[i] == 3) mrows.push_back(p);
                else if (rec.prole[i] == 2) brows.push_back(p);
            }
            for (size_t i = 0; i < prof.pts.size(); ++i)
                if (prof.rr[i] >= m_dot_body_exclusion_cm * units::cm) brows.push_back(prof.pts[i]);
            double ux = 0, uy = 0, uz = 0;
            for (const auto& p : mrows) { ux += p.x() - rec.stop_pt.x(); uy += p.y() - rec.stop_pt.y(); uz += p.z() - rec.stop_pt.z(); }
            const double un = std::sqrt(ux * ux + uy * uy + uz * uz);
            if (un > 1e-6 * units::cm) {
                ux /= un; uy /= un; uz /= un;
                SPDLOG_LOGGER_DEBUG(s_log, "{}CheckSTM_Michel michel-gamma: cluster {} conn {} rows {} dir {:.4f} {:.4f} {:.4f}",
                                    m_evt_tag, rec.cluster_id, rec.michel_conn_type, mrows.size(), ux, uy, uz);
                // Clusters the Michel object already spans are its own, not blobs.
                std::set<int> own;   // membership only, never iterated
                for (const auto& [sg, d] : michel_pieces) {
                    (void)d;
                    if (sg && sg->cluster()) own.insert(sg->cluster()->get_cluster_id());
                }
                if (michel_shower) {
                    IndexedVertexSet mv; IndexedSegmentSet ms;
                    michel_shower->fill_sets(mv, ms, false);
                    for (const auto& sg : ms) if (sg && sg->cluster()) own.insert(sg->cluster()->get_cluster_id());
                }
                for (auto* oc : companions) {
                    const int cid = oc->get_cluster_id();
                    if (own.count(cid)) continue;
                    auto segs = pa.find_cluster_segments(g, *oc);   // ordered_edges: deterministic
                    if (segs.empty()) continue;                     // no dx, no row: PDVD fits every companion
                    bool taken = false;
                    int gidx0 = -1;
                    std::vector<Point> pts;
                    for (const auto& seg : segs) {
                        if (!seg) continue;
                        const int gi = static_cast<int>(seg->get_graph_index());
                        if (rec.claimed.count(cid * 1000 + gi)) taken = true;
                        if (gidx0 < 0 || gi < gidx0) gidx0 = gi;
                        for (const auto& f : seg->fits()) if (f.dx > 0) pts.push_back(f.point);
                    }
                    if (taken || pts.empty()) continue;   // another object's, or nothing drawable
                    double d_stop = 1e9, d_m = 1e9, d_b = 1e9, cx = 0, cy = 0, cz = 0;
                    for (const auto& p : pts) {
                        d_stop = std::min(d_stop, (p - rec.stop_pt).magnitude());
                        for (const auto& m : mrows) d_m = std::min(d_m, (p - m).magnitude());
                        for (const auto& b : brows) d_b = std::min(d_b, (p - b).magnitude());
                        cx += p.x(); cy += p.y(); cz += p.z();
                    }
                    cx = cx / pts.size() - rec.stop_pt.x(); cy = cy / pts.size() - rec.stop_pt.y(); cz = cz / pts.size() - rec.stop_pt.z();
                    const double cn = std::sqrt(cx * cx + cy * cy + cz * cz);
                    const double cosv = cn > 0 ? (cx * ux + cy * uy + cz * uz) / cn : 1.0;
                    const double d_mich = std::min(d_stop, d_m);
                    double ke = 0;
                    for (const auto& seg : segs) if (seg) ke += segment_cal_kine_dQdx(seg, m_recomb_model) / units::MeV;
                    const int gate = stm_michel_gamma_gate(
                        d_stop / units::cm, oc->get_length() / units::cm, cosv, d_mich / units::cm, d_b / units::cm, ke,
                        m_michel_gamma_radius_cm, m_michel_gamma_max_len_cm, m_michel_gamma_cos_min, m_michel_gamma_max_ke_mev);
                    SPDLOG_LOGGER_DEBUG(s_log,
                        "{}CheckSTM_Michel michel-gamma-cl: cluster {} comp {} d_stop {:.3f} len {:.3f} cos {:.4f} d_mich {:.3f} d_body {:.3f} ke {:.4f} gate {}",
                        m_evt_tag, rec.cluster_id, cid, d_stop / units::cm, oc->get_length() / units::cm, cosv,
                        d_mich / units::cm, d_b / units::cm, ke, gate);
                    if (gate != 0) continue;
                    MGam c; c.segs = segs; c.d_stop = d_stop; c.ke = ke; c.cluster_id = cid; c.gidx = gidx0;
                    mgam.push_back(std::move(c));
                }
                // Nearest first: the total-energy guard makes the taken set
                // depend on the order, so the order is stated, not inherited.
                std::sort(mgam.begin(), mgam.end(), [](const MGam& a, const MGam& b) {
                    if (a.d_stop != b.d_stop) return a.d_stop < b.d_stop;
                    if (a.cluster_id != b.cluster_id) return a.cluster_id < b.cluster_id;
                    return a.gidx < b.gidx;
                });
                for (const auto& c : mgam)
                    for (const auto& seg : c.segs)
                        if (seg) rec.claimed.insert(c.cluster_id * 1000 + static_cast<int>(seg->get_graph_index()));
                rec.n_michel_gamma_cand = static_cast<int>(mgam.size());
            }
        }

        // doc pdvd/53: THE SURVEY.  Every admitted companion segment that no
        // stage claimed gets role-6 rows plus the gate that dropped it, so the
        // hand-scan display can draw it, select it and let the scanner group it.
        // This deliberately covers the companions INSIDE stop_gamma_radius_cm as
        // well as the ones only the survey radius reached: 039253_13 cluster
        // 102's isolated gamma is cluster 431 / segment 431017, already fitted
        // at doc pdvd/51 and invisible only because nothing named it.
        if (m_survey_enable && !companions.empty()) {
            struct SvSeg { SegmentPtr seg; int cluster_id, gidx, rej; double d_stop, d_body; };
            std::vector<SvSeg> sv;
            std::set<int> sv_clusters;
            for (auto* oc : companions) {
                const int cid = oc->get_cluster_id();
                auto segs = pa.find_cluster_segments(g, *oc);
                if (segs.empty()) { ++rec.n_survey_unfit; continue; }
                for (auto& seg : segs) {
                    if (!seg) continue;
                    const int gidx = static_cast<int>(seg->get_graph_index());
                    const int sid = cid * 1000 + gidx;
                    if (rec.claimed.count(sid)) continue;
                    int rej = 0; double d_stop = -1, d_body = -1;
                    auto its = sv_seg_rej.find(sid);
                    if (its != sv_seg_rej.end()) {
                        rej = its->second;
                        d_stop = sv_seg_d[sid].first; d_body = sv_seg_d[sid].second;
                    }
                    else {
                        auto itc = sv_cl_rej.find(cid);
                        // 10: there was no stop vertex, so the companion was
                        // fitted and then examined by nothing at all.
                        rej = (itc != sv_cl_rej.end()) ? itc->second : (stop_v ? 9 : 10);
                        if (sv_cl_dstop.count(cid)) d_stop = sv_cl_dstop[cid];
                        if (sv_cl_dbody.count(cid)) d_body = sv_cl_dbody[cid];
                        if (d_stop < 0) {
                            const auto [ccp, cblob] = oc->get_closest_point_blob(rec.stop_pt);
                            d_stop = (ccp - rec.stop_pt).magnitude();
                        }
                    }
                    sv.push_back({seg, cid, gidx, rej, d_stop, d_body});
                    sv_clusters.insert(cid);
                }
            }
            std::sort(sv.begin(), sv.end(), [](const SvSeg& a, const SvSeg& b) {
                if (a.cluster_id != b.cluster_id) return a.cluster_id < b.cluster_id;
                return a.gidx < b.gidx;
            });
            for (const auto& e : sv)
                add_points(rec, e.seg, 6, nullptr, e.rej, e.d_stop, e.d_body, /*keep_dead*/ true);
            rec.n_survey_segs = static_cast<int>(sv.size());
            rec.n_survey_clusters = static_cast<int>(sv_clusters.size());
        }

        // doc pdvd/64 (T6): publish the kOther arms as role 7, LAST -- after the
        // Michel object, the capture gamma and the survey have claimed what
        // they claim -- so an interior arm the Michel shower walk pulled in
        // keeps its role-3 row, and role 7 itself claims nothing (add_points
        // exempts it, like the survey's 6).  Rows only; no verdict reads them.
        if (m_publish_other_arms && !other_arms.empty()) {
            std::set<int> seen;   // membership only, never iterated
            for (const auto& seg : other_arms) {
                if (!seg) continue;
                const int sid = (seg->cluster() ? seg->cluster()->get_cluster_id() : 0) * 1000 + static_cast<int>(seg->get_graph_index());
                if (rec.claimed.count(sid) || seen.count(sid)) continue;
                seen.insert(sid);
                add_points(rec, seg, 7);
                ++rec.n_other_published;
            }
        }

        // Unfitted companion clusters have no dx, so dQ/dx cannot be inverted;
        // the only route is charge (doc pdhd/15 sec 6).  Computed before the
        // shower is energised because the object energy is stamped back onto it.
        // doc pdhd/17: the flat pair and this component's own dQ/dx -> dE/dx
        // inverse are two carriers of ONE quantity, and doc pdhd/16 moved only
        // the second.  With michel_unfit_from_model the survival comes out of
        // the model actually bound here, so the two cannot drift apart again.
        // MIP-EQUIVALENT either way -- an unfitted cluster has no dx.
        rec.dots_ke_unfit =
            (m_michel_unfit_from_model
                 ? stm_michel_charge_to_energy_model(rec.dots_charge_unfit, m_recomb_model,
                                                     m_michel_unfit_dedx)
                 : stm_michel_charge_to_energy(rec.dots_charge_unfit, m_michel_unfit_recom,
                                               m_michel_unfit_fudge, m_michel_unfit_w_ev))
            / units::MeV;

        // ---- energise the object, once ---------------------------------------
        // doc pdhd/15 sec 5.  PatternAlgorithms::calculate_shower_kinematics
        // (NeutrinoEnergyReco.cxx:303) is the chain's own production entry
        // point: it runs Shower::calculate_kinematics AND sets kine_charge from
        // the 2-D charge maps.  Using it rather than a local recipe keeps the
        // owner's standing rule from doc pdhd/14 -- every number in this tree is
        // the chain's, not this file's.
        // doc pdvd/51: the capture gammas join the SAME shower set, so they are
        // energised by the same call and reach the same tf->set_showers()
        // publication.  When the feature is off gamma_showers is empty and this
        // block is line-for-line the doc pdhd/17 one.
        for (const auto& [gsh, gke] : gamma_showers) { (void)gke; showers.insert(gsh); }
        if (michel_shower) {
            showers.insert(michel_shower);
            michel_shower->set_particle_type(11);
        }
        if (!showers.empty()) {
            IndexedVertexSet lm_v; IndexedSegmentSet lm_s;
            pa.calculate_shower_kinematics(showers, lm_v, lm_s, g, *tf, m_dv,
                                           particle_data(), m_recomb_model);
        }
        for (const auto& [gsh, gke] : gamma_showers) {
            // Same reason as the Michel's stamp below: for a conn-2 EM shower
            // kenergy_best is 0 and get_kine_best() falls back to
            // kenergy_charge, which is 0 on this path (doc pdhd/15 sec 6), so
            // the renderer's em_ke_min prune would delete the node entirely.
            gsh->set_kine_best(gke * units::MeV);
        }
        if (michel_shower) {
            rec.michel_ke_dqdx = michel_shower->get_kine_dQdx() / units::MeV;
            rec.michel_ke_charge = michel_shower->get_kine_charge() / units::MeV;
            IndexedVertexSet mv2; IndexedSegmentSet ms2;
            michel_shower->fill_sets(mv2, ms2, false);
            rec.n_michel_segs = static_cast<int>(ms2.size());
            // doc pdvd/51: > 1 says the object is not in the candidate's own
            // cluster -- the fact a T_stm_michel-only consumer cannot derive.
            {
                std::set<int> mcl;
                for (const auto& sg : ms2)
                    if (sg && sg->cluster()) mcl.insert(sg->cluster()->get_cluster_id());
                rec.michel_n_clusters = static_cast<int>(mcl.size());
            }
            rec.michel_start_pt = michel_shower->get_start_point();
            if (auto sv = michel_shower->start_vertex())
                rec.michel_parent_vtx_id = rec.cluster_id * 1000 + static_cast<int>(sv->get_graph_index());
            // Stamp the object energy back onto the Shower.  Shower::kenergy_best
            // is 0 whenever a member is graph-disconnected (PRShower.cxx:1855) --
            // which is EVERY bridged Michel and every attached one that gathered
            // a companion piece -- and get_kine_best() then falls back to
            // kenergy_charge, 0 on this path (sec 6).  fill_bee_pf_tree prunes an
            // EM leaf whose ke is below em_ke_min (MultiAlgBlobClustering.cxx:2082),
            // so leaving it at 0 DELETES the e- node from mc.json: 039252_15
            // cluster 91 lost the daughter it had had since doc pdvd/48.  The
            // module's own object energy is the right value for it to carry.
            michel_shower->set_kine_best((rec.michel_ke_dqdx + rec.dots_ke_unfit) * units::MeV);
        }
        else if (!michel_pieces.empty()) {
            // build_michel_shower off: no Shower object exists, so the object
            // energy is the plain sum over its fitted pieces.
            rec.n_michel_segs = static_cast<int>(michel_pieces.size());
            double sum = 0;
            for (const auto& [sg, d] : michel_pieces) { (void)d; sum += segment_cal_kine_dQdx(sg, m_recomb_model) / units::MeV; }
            rec.michel_ke_dqdx = sum;
        }
        if (rec.michel_parent_vtx_id < 0 && rec.michel_conn_type > 0) rec.michel_parent_vtx_id = rec.stop_vtx_id;

        rec.michel_n_pieces = rec.n_michel_segs + rec.n_dot_clusters_unfit;

        // The owner's rule, doc pdhd/15 sec 1: "there should be range, dQ/dx for
        // the track.  For the dots etc ... either dQ/dx -> dE/dx or the charge
        // conversion."  So `best` is dQ/dx over everything fitted plus the
        // charge term for what is not -- NOT Shower::get_kine_best(), which for
        // a shower-flagged or graph-disconnected object falls back to
        // kenergy_charge computed with the SHOWER recombination pair and
        // overshoots the chain's own dQ/dx by ~1.66x (sec 6).
        rec.michel_ke_best = rec.michel_ke_dqdx + rec.dots_ke_unfit;

        // doc pdvd/61 (T2c): an ATTACHED Michel (conn_type 1 -- the arm hangs
        // off the stop vertex itself, dis_cm 0 by construction, :1742; since
        // doc pdvd/83 also an arm leaving the chain up to
        // michel_near_stop_arm_cm before the stop, dis_cm = that distance --
        // this veto reads neither dis_cm nor where the arm leaves) whose
        // stop was MOVED this event (a retreat or split fired -- the mechanism
        // that finds a stop the tagger's own fit missed, docs pdvd/57/58) is
        // not itself evidence of a Michel: the census found 5 named
        // through-going items (039252_2/79, 039252_4/55, 039349_20/41,
        // 039349_48/21, 039349_61/62) picking up exactly this spurious
        // attachment, 4.6-8.9 MeV, against a 21.4 MeV median for genuine
        // conn_type==1 Michels.  Demote by resetting michel_conn_type (which
        // michel_found reads next) -- deliberately NOT via reject_bits, so
        // is_stm is untouched (that is the whole design constraint: this is
        // a michel_found-only veto).  n_michel_veto records that it fired
        // even when nothing else about the record changes visibly.
        if (m_moved_stop_michel_guard && rec.michel_conn_type == 1 &&
            (rec.n_retreat > 0 || rec.n_split > 0) &&
            rec.michel_ke_best < m_moved_stop_michel_ke_min) {   // both plain MeV, no units:: scale (michel_ke_best is already / units::MeV)
            // doc pdvd/72 (P3b): the KE floor does not read the arm's turn.
            // On the owner's scan the two real Michels this veto demoted turn
            // 59.6 and 132.6 deg, the three through-going items 17-59 deg; a
            // hard turn is the Michel's own evidence, so it is spared.  Off at
            // -1; michel_kink_deg is -1 when unmeasurable, which never spares.
            // doc pdvd/84 (doc 78 item 3): then the arm's reach (length + far
            // subtree), for the owner-confirmed Michel at 59.6 deg that sits
            // 0.9 deg above a through-going arm; off at -1.  With it off the
            // predicate is the doc 72 expression above, verbatim.
            switch (stm_michel_moved_stop_spare(rec.michel_kink_deg,
                                                (rec.michel_len + rec.michel_far_len) / units::cm,
                                                m_moved_stop_michel_kink_min, m_moved_stop_michel_reach_min_cm)) {
            case StmMichelMovedStopSpare::kKink:
                ++rec.n_michel_veto_exempt;
                break;
            case StmMichelMovedStopSpare::kReach:
                ++rec.n_michel_veto_reach_exempt;
                break;
            case StmMichelMovedStopSpare::kVeto:
                rec.michel_conn_type = 0;
                ++rec.n_michel_veto;
                break;
            }
        }
        // doc pdvd/62 (T3c): doc 55 sec 15.1's kinematic test.  A bridged
        // (2) or charge-only (3) Michel more than michel_range_energy_dis_cm
        // from the stop with less than michel_range_energy_ke_min MeV cannot
        // be an electron born at the stop (0 of 113 attached objects fail
        // this; 26 of 45 bridged ones do, median 9.1 cm / 4.1 MeV -- a gamma
        // or debris).  Same demotion route as the veto above: michel_conn_type
        // only, so is_stm is untouched by construction; the role-3 point rows
        // already written stay, as they do for the T2c veto.
        if (m_michel_range_energy_guard && (rec.michel_conn_type == 2 || rec.michel_conn_type == 3) &&
            rec.michel_dis_cm > m_michel_range_energy_dis_cm &&
            rec.michel_ke_best < m_michel_range_energy_ke_min) {   // cm and MeV, both plain
            rec.michel_conn_type = 0;
            ++rec.n_michel_range_veto;
        }

        // doc pdhd/15 sec 4: michel_found now means "a Michel object exists" --
        // 1 attached, 2 bridged across a clustering gap, 3 charge only.  Through
        // doc pdhd/14 it was set on the attached path only, so 33 PDHD / 41 PDVD
        // reconstructed Michels reported 0, about 30 % of them (doc pdhd/13
        // defect D1).  The old meaning is exactly
        // (michel_found && michel_conn_type == 1).
        rec.michel_found = (rec.michel_conn_type > 0) ? 1 : 0;

        // A NaN passes no gate and fails every one silently (PDVD 039349_3
        // cluster 26 persisted michel_ke_best = NaN through doc pdhd/14).
        for (double* e : {&rec.michel_ke_dqdx, &rec.michel_ke_range, &rec.michel_ke_best,
                          &rec.michel_ke_core, &rec.michel_ke_charge, &rec.dots_ke_dqdx,
                          &rec.dots_ke_unfit, &rec.muon_ke_range, &rec.muon_ke_dqdx, &rec.muon_ke_best,
                          &rec.stop_gamma_ke_tot, &rec.stop_gamma_ke_max, &rec.stop_gamma_charge}) {
            if (!std::isfinite(*e)) {
                SPDLOG_LOGGER_WARN(s_log, "{}CheckSTM_Michel: cluster {} produced a non-finite energy; zeroed",
                                   m_evt_tag, rec.cluster_id);
                *e = 0;
            }
        }
        // doc pdhd/16: the MCS fields carry a -1 = "not computed" sentinel, so
        // a non-finite one goes back to -1, not to 0 (a zero there would read
        // as a measured zero energy and pass every `>= 0` gate).
        for (double* e : {&rec.muon_ke_mcs, &rec.muon_mcs_amb, &rec.muon_mcs_tracklen,
                          &rec.muon_mcs_range_ke, &rec.muon_p_range, &rec.muon_p_dqdx,
                          &rec.muon_p_mcs}) {
            if (!std::isfinite(*e)) {
                SPDLOG_LOGGER_WARN(s_log, "{}CheckSTM_Michel: cluster {} produced a non-finite MCS field; set to -1",
                                   m_evt_tag, rec.cluster_id);
                *e = -1;
            }
        }
        tf->set_showers(showers);

        // ---- doc pdhd/03 sec 6: is the cluster a track (+ attachments) at all?
        // Fraction of the main cluster's 3-D points within coverage_radius of any
        // reconstructed point (chain, deltas, Michel, dots).  Brute force: a few
        // thousand cluster points x a few hundred reconstructed points.
        if (m_min_chain_coverage > 0 && !rec.px.empty()) {
            const double r2 = std::pow(m_coverage_radius_cm * units::cm, 2);
            const int n = main->npoints();
            int nin = 0;
            for (int i = 0; i < n; ++i) {
                const auto p = main->point3d(i);
                bool in = false;
                for (size_t j = 0; j < rec.px.size() && !in; ++j) {
                    const double dx = p.x() - rec.px[j], dy = p.y() - rec.py[j], dz = p.z() - rec.pz[j];
                    in = (dx * dx + dy * dy + dz * dz) < r2;
                }
                if (in) ++nin;
            }
            rec.n_cluster_pts = n;
            rec.chain_coverage = n > 0 ? double(nin) / n : -1.0;
            if (n > 0 && rec.chain_coverage < m_min_chain_coverage) rec.reject_bits |= R_CLUSTER_NOT_TRACK;
        }

        // ---- containment of the stop -----------------------------------------
        if (fiducial_utils) {
            rec.in_fv = fiducial_utils->inside_fiducial_volume(rec.stop_pt, fv_tol) ? 1 : 0;
            if (!rec.in_fv) rec.reject_bits |= R_STOP_NEAR_BOUNDARY;
        }

        // ---- doc pdvd/70 (P1): topology-first stop evidence -------------------
        // Here and not earlier: R_CLUSTER_NOT_TRACK and R_STOP_NEAR_BOUNDARY are
        // OR'd in just above, so every bit is final and a clear that fires is
        // exactly an item whose only objections were the dQ/dx-shape tests (or
        // the sparse profile) -- topology_cleared_bits then names what moved.
        if (m_topology_stop_evidence) {
            const unsigned clr = stm_michel_topology_clear(
                rec.reject_bits, rec.michel_found, rec.michel_conn_type,
                rec.michel_ke_best, rec.michel_len / units::cm,   // plain MeV (see :2703), internal length -> cm
                m_topology_michel_ke_min, m_topology_michel_len_min_cm, m_topology_clears_sparse);
            if (clr) {
                rec.reject_bits &= ~clr;
                rec.topology_cleared_bits = static_cast<int>(clr);
                SPDLOG_LOGGER_DEBUG(s_log, "{}CheckSTM_Michel: cluster {} topology cleared {} (michel conn {} {:.1f} MeV {:.1f} cm) -> bits {}",
                                    m_evt_tag, rec.cluster_id, bits_string(clr), rec.michel_conn_type,
                                    rec.michel_ke_best, rec.michel_len / units::cm, bits_string(rec.reject_bits));
            }
        }

        // ---- doc pdvd/85 (doc 78 action item 5): a capture gamma only on a stopper
        // The capture stage (doc pdvd/51) runs on every candidate and publishes
        // before the verdict exists; here, with every reject bit final (nothing
        // below sets one, and the coverage test above has already counted these
        // rows), a REJECTED candidate's capture gammas are withheld.  Their
        // role-5 rows go (order of the rest kept), their segments get back the
        // particle info they had before set_pdg, their showers leave the
        // published set -- calculate_shower_kinematics energised every shower
        // on its own, so the Michel's numbers above are the knob-off ones -- and
        // stop_gamma_* read as if the stage had accepted nothing.  rec.claimed
        // keeps them, so the gamma collect and the census see the knob-off pool.
        // With the survey on they come back as role-6 rows, rej 16, so the scan
        // display still draws them.  An accepted candidate is untouched.
        if (stm_michel_stop_gamma_withhold(m_stop_gamma_require_stm, rec.reject_bits) &&
            (rec.n_stop_gammas > 0 || rec.stop_gamma_n_unfit > 0)) {
            const std::vector<char> keep = stm_michel_rows_keep(rec.prole, 5);
            auto squeeze = [&keep](auto& v) {
                if (v.size() != keep.size()) return;   // rej / d_stop / d_body exist only with the survey or census on
                size_t j = 0;
                for (size_t i = 0; i < v.size(); ++i)
                    if (keep[i]) v[j++] = v[i];
                v.resize(j);
            };
            squeeze(rec.px); squeeze(rec.py); squeeze(rec.pz); squeeze(rec.pq); squeeze(rec.pL); squeeze(rec.prr);
            squeeze(rec.pseg); squeeze(rec.pmed); squeeze(rec.prej); squeeze(rec.pdstop); squeeze(rec.pdbody);
            squeeze(rec.prole);
            rec.role_segs.erase(5);
            for (const auto& w : sg_withhold) {
                w.seg->particle_info(w.info);
                w.seg->particle_score(w.score);
            }
            if (!gamma_showers.empty()) {
                for (const auto& [gsh, gke] : gamma_showers) { (void)gke; showers.erase(gsh); }
                tf->set_showers(showers);
            }
            rec.n_stop_gammas_withheld = rec.n_stop_gammas + rec.stop_gamma_n_unfit;
            SPDLOG_LOGGER_DEBUG(s_log, "{}CheckSTM_Michel stop-gamma-withheld: cluster {} bits {} gammas {} unfit {} segs {} showers {}",
                                m_evt_tag, rec.cluster_id, bits_string(rec.reject_bits), rec.n_stop_gammas,
                                rec.stop_gamma_n_unfit, sg_withhold.size(), gamma_showers.size());
            rec.n_stop_gammas = 0; rec.stop_gamma_n_unfit = 0; rec.stop_gamma_seg_id = -1;
            rec.stop_gamma_ke_tot = 0; rec.stop_gamma_ke_max = 0; rec.stop_gamma_charge = 0;
            rec.stop_gamma_dis_min = -1; rec.stop_gamma_dis_max = -1;
            if (m_survey_enable)
                for (const auto& w : sg_withhold)
                    add_points(rec, w.seg, 6, nullptr, 16, w.d_stop, -1, /*keep_dead*/ true);
        }

        // ---- doc pdvd/71 (P4): the Michel's gamma blobs, phase 2 --------------
        // Last, after every verdict input is final: michel_found has been
        // through the T2c / T3c vetoes, michel_ke_best through the NaN guard,
        // and the chain_coverage test above has counted rec.px WITHOUT these
        // rows (it reads every row, so writing them earlier would move
        // R_CLUSTER_NOT_TRACK).  A reserved blob of a vetoed Michel, or one the
        // total-energy guard refuses, gets no member row; with the survey on it
        // gets its role-6 row back, rej 12 (vetoed) or 11 (energy guard).
        if (m_michel_gamma_collect) {
            const bool live = rec.michel_found && (rec.michel_conn_type == 1 || rec.michel_conn_type == 2);
            std::vector<double> kes;
            kes.reserve(mgam.size());
            for (const auto& c : mgam) kes.push_back(c.ke);
            const std::vector<int> take = live
                ? stm_michel_gamma_take(rec.michel_ke_best, kes, m_michel_gamma_total_ke_max_mev)
                : std::vector<int>(mgam.size(), 0);
            for (size_t i = 0; i < mgam.size(); ++i) {
                const auto& c = mgam[i];
                if (take[i]) {
                    ++rec.n_michel_gammas;
                    rec.michel_ke_gamma += c.ke;
                    rec.michel_gamma_dis_max = std::max(rec.michel_gamma_dis_max, c.d_stop / units::cm);
                    for (const auto& sg : c.segs) if (sg) add_points(rec, sg, 4);
                }
                else {
                    if (live) ++rec.n_michel_gamma_capped;
                    if (m_survey_enable)
                        for (const auto& sg : c.segs)
                            if (sg) add_points(rec, sg, 6, nullptr, live ? 11 : 12, c.d_stop, -1, /*keep_dead*/ true);
                }
                SPDLOG_LOGGER_DEBUG(s_log, "{}CheckSTM_Michel michel-gamma-take: cluster {} comp {} ke {:.4f} take {} live {} core {:.4f}",
                                    m_evt_tag, rec.cluster_id, c.cluster_id, c.ke, take[i], live ? 1 : 0, rec.michel_ke_best);
            }
            rec.michel_ke_total = rec.michel_ke_best + rec.michel_ke_gamma;
        }

        // doc pdvd/80 (doc 78 action item 1): THE SEGMENT CENSUS, written last
        // of all -- after the Michel object, the capture gamma, the survey, the
        // kOther arms and the gamma collect have claimed or published what they
        // do.  Every remaining PR segment of the MAIN cluster gets role-8 rows
        // with the doc 62 T3b piece gates read on it, so the hand-scan display
        // can draw it, name it and say why the chain left it: on the production
        // record 42 of the scanner's michel tags (32 items) sit on exactly such
        // segments, 21 of them touching the stop, and four missed Michels run
        // backward along the muon body inside the body exclusion (doc 78 sec
        // 3.2).  Rows only: role 8 claims nothing (add_points exempts it) and no
        // verdict reads it; the chain_coverage test above counted rec.px before
        // these rows exist.  Deterministic: find_cluster_segments walks
        // ordered_edges, and the rows are sorted by graph index.
        if (m_segment_census) {
            const std::set<int> written(rec.pseg.begin(), rec.pseg.end());          // membership only, never iterated
            const std::set<VertexPtr> chain_vset(chain_vtxs.begin(), chain_vtxs.end());   // membership only, never iterated
            struct CsSeg { SegmentPtr seg; int gidx, rej; double d_stop, d_body; };
            std::vector<CsSeg> cs;
            for (auto& seg : pa.find_cluster_segments(g, *main)) {
                if (!seg || chain_set.count(seg)) continue;
                const int gidx = static_cast<int>(seg->get_graph_index());
                const int sid = main->get_cluster_id() * 1000 + gidx;
                if (rec.claimed.count(sid) || written.count(sid)) continue;
                int rej = 10; double d_stop = -1, d_body = -1;
                if (stop_v) {
                    auto [va, vb] = find_vertices(g, seg);
                    bool touches_chain = false;
                    if (va && vb) {
                        touches_chain = chain_vset.count(va) || chain_vset.count(vb);
                        for (VertexPtr v : {va, vb}) {
                            if (touches_chain) break;
                            for (auto e : sorted_out_edges(v->get_descriptor(), g)) {
                                auto s2 = g[e].segment;
                                if (s2 && chain_set.count(s2)) { touches_chain = true; break; }
                            }
                        }
                    }
                    auto [ds, cp] = segment_get_closest_point(seg, rec.stop_pt, "fit", "main");
                    d_stop = ds;
                    double db = 1e9;
                    for (size_t i = 0; i < prof.pts.size(); ++i) {
                        if (prof.rr[i] < m_dot_body_exclusion_cm * units::cm) continue;
                        db = std::min(db, (prof.pts[i] - cp).magnitude());
                    }
                    if (db < 1e9) d_body = db;
                    if (!va || !vb) rej = 15;
                    else if (touches_chain) rej = 13;
                    else if (d_stop > m_michel_dot_radius_cm * units::cm) rej = 2;
                    else if (segment_track_length(seg) > m_dot_max_len_cm * units::cm) rej = 3;
                    else if (d_body < d_stop) rej = 4;
                    else rej = m_stop_local_michel_pieces ? 14 : 9;
                }
                cs.push_back({seg, gidx, rej, d_stop, d_body});
            }
            std::sort(cs.begin(), cs.end(), [](const CsSeg& a, const CsSeg& b) { return a.gidx < b.gidx; });
            for (const auto& e : cs) {
                add_points(rec, e.seg, 8, nullptr, e.rej, e.d_stop, e.d_body, /*keep_dead*/ true);
                if (e.rej == 14) ++rec.n_census_admissible;
            }
            rec.n_census_segs = static_cast<int>(cs.size());
            SPDLOG_LOGGER_DEBUG(s_log, "{}CheckSTM_Michel segment-census: cluster {} unclaimed PR segments {} (rej 14: {})",
                                m_evt_tag, rec.cluster_id, rec.n_census_segs, rec.n_census_admissible);
        }

        // ---- doc pdvd/81: the charge-based Michel energy and the 2-D cells ---
        // After every row and verdict is final (the gamma take above decides
        // the role-4 set); reads the fit's stored response, writes only its
        // own branches.
        if (m_michel_q2d) michel_q2d_estimate(rec, *tf, pa, main, chain, chain_set, michel_shower);

        // ---- publish (TaggerCheckNeutrino.cxx:3580-3590) --------------------
        tf->assemble_fitted_charge_2d();
        if (ci == 0) grouping.set_track_fitting(tf);
        if (m_publish_nu_slots) grouping.set_track_fitting("nu" + std::to_string(ci), tf);
        persist(*main, rec);

        if (rec.reject_bits == 0) ++n_stm;
        if (rec.michel_found) ++n_michel;
        if (rec.n_stop_gammas > 0) { ++n_gamma_cand; n_gamma += rec.n_stop_gammas; }
        SPDLOG_LOGGER_INFO(s_log,
            "{}CheckSTM_Michel: cluster {} gid {} verdict {} bits {} | chain {} segs {:.1f} cm ({} pts, {} dead) stop_dis {:.1f} cm | "
            "contrast {:.2f}/{:.2f} ks_mu {:.3f} ks_flat {:.3f} comp_fwd {:.0f}/{:.2f}/{:.2f}/{:.2f} | "
            "mu E range {:.1f} dQ/dx {:.1f} MCS {:.1f} MeV (amb {:.2f}, {} segs, {:.1f} cm, cath {}/{}) | "
            "delta {} hadron {} | michel {} conn {} ({} segs / {} pieces, {:.1f} cm, kink {:.0f} deg, gap {:.1f} cm) "
            "E {:.1f} MeV = dQ/dx {:.1f} + unfit {:.1f} (core {:.1f}, range {:.1f}, charge {:.1f}) dots {} ({:.1f} MeV) unfit_cl {} | "
            "gammas {} ({:.1f} MeV, max {:.1f}, {:.1f}-{:.1f} cm, {} unfit) | "
            "cont {:.1f} cm @ {:.0f} deg ext {} ({:.1f} cm) dead_ahead {} cov {:.2f} | in_fv {} | {:.0f} ms",
            m_evt_tag, rec.cluster_id, rec.gid, bits_string(rec.reject_bits), rec.reject_bits,
            rec.n_chain_segs, rec.muon_len / units::cm, rec.n_profile_pts, rec.n_dead_pts, rec.stop_dis / units::cm,
            rec.bragg.contrast, rec.bragg.expected, rec.ks_mu, rec.ks_flat,
            rec.comp_fwd[0], rec.comp_fwd[1], rec.comp_fwd[2], rec.comp_fwd[3],
            rec.muon_ke_range, rec.muon_ke_dqdx, rec.muon_ke_mcs, rec.muon_mcs_amb,
            rec.muon_mcs_nsegs, rec.muon_mcs_tracklen,
            rec.muon_mcs_cathode_segs, rec.muon_mcs_cathode_angles,
            rec.n_delta, rec.n_body_hadron,
            rec.michel_found, rec.michel_conn_type, rec.n_michel_segs, rec.michel_n_pieces,
            rec.michel_len / units::cm, rec.michel_kink_deg, rec.michel_dis_cm,
            rec.michel_ke_best, rec.michel_ke_dqdx, rec.dots_ke_unfit,
            rec.michel_ke_core, rec.michel_ke_range, rec.michel_ke_charge,
            rec.n_dots, rec.dots_ke_dqdx, rec.n_dot_clusters_unfit,
            rec.n_stop_gammas, rec.stop_gamma_ke_tot, rec.stop_gamma_ke_max,
            rec.stop_gamma_dis_min, rec.stop_gamma_dis_max, rec.stop_gamma_n_unfit,
            rec.cont_len / units::cm, rec.cont_angle_deg, rec.n_ext, rec.ext_len / units::cm, rec.dead_ahead, rec.chain_coverage, rec.in_fv,
            MS(Clock::now() - t0).count());
    }

    // doc pdhd/15: "with a Michel" now counts the detached ones too (michel_found
    // means "a Michel object exists"); it used to count the attached only.
    SPDLOG_LOGGER_INFO(s_log,
                       "{}CheckSTM_Michel: {} candidate(s), {} pass every check, {} with a Michel, "
                       "{} with a capture gamma ({} gammas); {:.0f} ms",
                       m_evt_tag, candidates.size(), n_stm, n_michel, n_gamma_cand, n_gamma,
                       MS(Clock::now() - t_total).count());
}
