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
        double michel_dis_cm{-1};                  // stop -> object start; 0 when attached
        Point michel_start_pt;
        double cont_len{0}, cont_angle_deg{-1}, cont_mip{0};
        int n_ext{0}; double ext_len{0};          // doc pdhd/03: chain extensions past the tagger's stop
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
        // doc pdvd/53: the survey.  n_survey_clusters counts admitted companion
        // clusters that yielded at least one UNCLAIMED fitted segment;
        // n_survey_segs counts those segments; n_survey_unfit counts admitted
        // companions that produced no segment at all (so they have charge but
        // no dx, and the display must fall back to their image points).
        int n_survey_clusters{0}, n_survey_segs{0}, n_survey_unfit{0};
        int in_fv{-1};
        // points for the Bee/ROOT layer
        std::vector<double> px, py, pz, pq, pL, prr;
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
        if (!m_survey_enable) return;
        for (size_t i = 0; i < n; ++i) {
            rec.prej.push_back(rej); rec.pdstop.push_back(d_stop); rec.pdbody.push_back(d_body);
        }
    }

    void add_points(Record& rec, const SegmentPtr& seg, int role, const StmMichelProfile* prof = nullptr,
                    int rej = 0, double d_stop = -1, double d_body = -1) const {
        const int seg_id = (seg && seg->cluster() ? seg->cluster()->get_cluster_id() : 0) * 1000
                         + (seg ? static_cast<int>(seg->get_graph_index()) : 0);
        if (role != 6) rec.claimed.insert(seg_id);
        const size_t n0 = rec.px.size();
        if (prof) {
            for (size_t i = 0; i < prof->pts.size(); ++i) {
                rec.px.push_back(prof->pts[i].x()); rec.py.push_back(prof->pts[i].y()); rec.pz.push_back(prof->pts[i].z());
                rec.pq.push_back(prof->dQdx[i]); rec.pL.push_back(prof->L[i]); rec.prr.push_back(prof->rr[i]);
                rec.prole.push_back(role); rec.pseg.push_back(seg_id);
            }
            add_survey_cols(rec, rec.px.size() - n0, rej, d_stop, d_body);
            return;
        }
        if (!seg) return;
        for (const auto& f : seg->fits()) {
            if (f.dx <= 0) continue;
            rec.px.push_back(f.point.x()); rec.py.push_back(f.point.y()); rec.pz.push_back(f.point.z());
            rec.pq.push_back(f.dQ / (f.dx / units::cm)); rec.pL.push_back(-1); rec.prr.push_back(-1);
            rec.prole.push_back(role); rec.pseg.push_back(seg_id);
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
        // doc pdvd/53: the survey.  Written only when the knob is on, so the
        // knob-off T_stm_michel branch list is exactly the doc pdvd/51 one.
        if (m_survey_enable) {
            I1("n_survey_clusters", r.n_survey_clusters); I1("n_survey_segs", r.n_survey_segs);
            I1("n_survey_unfit", r.n_survey_unfit);
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
            // doc pdvd/53.  Absent when the survey is off, so write_pc_tree
            // (PdvdPrMagnifyTrackingVisitor.cxx:293, which takes its column set
            // from the first carrier) reproduces the old schema exactly.
            if (m_survey_enable) {
                p.emplace("rej", Array(r.prej));
                p.emplace("d_stop", Array(tocm(r.pdstop)));
                p.emplace("d_body", Array(tocm(r.pdbody)));
            }
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
        const double admit_radius = stm_michel_admit_radius(
            m_michel_dot_radius_cm, m_stop_gamma_radius_cm, m_stop_gamma_enable,
            m_survey_radius_cm, m_survey_enable);
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
        }
        if (rec.n_ext > 0 && stop_v) {
            rec.stop_vtx_id = rec.cluster_id * 1000 + static_cast<int>(stop_v->get_graph_index());
            rec.stop_pt = stm_michel_vertex_point(stop_v);
            rec.stop_dis = (rec.stop_pt - rec.tagger_stop_pt).magnitude();
            rec.n_chain_segs = static_cast<int>(chain.size());
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
        }

        // ---- arms: body (delta rays / hadrons) and stop (Michel / continuation)
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
                std::vector<StmMichelArm> stop_arms;
                for (auto e : sorted_out_edges(stop_v->get_descriptor(), g)) {
                    auto arm = g[e].segment;
                    if (!arm || arm == last || chain_set.count(arm)) continue;
                    stop_arms.push_back(stm_michel_classify_stop_arm(g, last, arm, stop_v, th));
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

        if (stop_v && !companions.empty()) {
            // Two passes: collect every admissible piece with its distance to the
            // stop, then assemble.  The seed of a BRIDGED object is then the
            // NEAREST piece rather than whichever the graph's edge order reached
            // first, and michel_seg_id names the piece the gap is measured to.
            struct Piece { SegmentPtr seg; double d_stop; int cluster_id, gidx; };
            std::vector<Piece> pieces;
            double d_unfit = 1e9;            // closest approach of an UNFITTED companion
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
            for (const auto& e : sv) add_points(rec, e.seg, 6, nullptr, e.rej, e.d_stop, e.d_body);
            rec.n_survey_segs = static_cast<int>(sv.size());
            rec.n_survey_clusters = static_cast<int>(sv_clusters.size());
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
