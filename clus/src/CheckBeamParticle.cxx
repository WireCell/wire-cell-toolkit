/** CheckBeamParticle -- beam-particle pattern recognition rooted at the beam
 *  entry point (doc pdvd/120).
 *
 * On a test-beam detector (ProtoDUNE-VD) the in-time activity is a charged
 * beam particle ENTERING the detector, so its "main vertex" is a known
 * geometric fact -- the point where it crosses the active face -- and not a
 * topological choice.  This visitor runs the neutrino PR chain on the
 * beam-flash-matched bundle only, with three differences from
 * TaggerCheckNeutrino (which stays byte-for-byte untouched; this is a fork by
 * duplication, M10):
 *
 *   1. the bundle is the one matched to the BEAM flash (the brightest flash
 *      inside the per-event beam window on cluster_t0, the rule the Bee
 *      op_beam label uses), never "every bundle";
 *   2. the main cluster is the bundle member closest to the nominal beam
 *      entry point, and the main vertex is that cluster's axis end nearer
 *      the nominal entry, snapped onto the PR graph -- the DL/SCN vertex and
 *      the whole determine_overall_main_vertex ladder are not run;
 *   3. no cosmic taggers (a crossing beam particle IS a "through-going" track
 *      to them) and no BDT features: TaggerInfo stays at its defaults except
 *      match_isFC.
 *
 * Everything else -- proto-segment partition, companion clusters, deghosting,
 * orientation from the entry, EM shower clustering, particle id, kinematics
 * -- is the same PatternAlgorithms call sequence TaggerCheckNeutrino::visit
 * runs (line references in the body), and publication is identical
 * (unnamed TrackFitting slot + "nu0"), so the Bee PR layers
 * (track_fit/shower_track/vertices/mc), PdvdPrMagnifyTrackingVisitor,
 * UbooneTaggerOutputVisitor (T_kine) and PrDisplayDump render this stage's
 * output unchanged with visitor 'CheckBeamParticle:<prefix>'.
 *
 * The three selection rules are pure functions in BeamParticleFunctions.h,
 * doctested on plain rows.  Nothing here runs unless the component is named
 * in a pipeline, and beam_window_low >= beam_window_high (the C++ default)
 * makes the stage a no-op, so no existing output can change.
 */

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/ParticleDataSet.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRShower.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/NeutrinoTaggerInfo.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/TrackFittingPresets.h"
#include "WireCellClus/BeamParticleFunctions.h"
#include "WireCellClus/Facade_Flash.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Persist.h"
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

class CheckBeamParticle;
WIRECELL_FACTORY(CheckBeamParticle, CheckBeamParticle,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::Clus::PR;

static auto s_log = WireCell::Log::logger("clus.CheckBeamParticle");

class CheckBeamParticle : public IConfigurable, public Clus::IEnsembleVisitor,
                          private Clus::NeedDV, private Clus::NeedPCTS,
                          private Clus::NeedRecombModel, private Clus::NeedParticleData,
                          private Clus::NeedFiducial {
public:
    CheckBeamParticle() {
        m_track_fitter = std::make_shared<TrackFitting>(TrackFittingPresets::create_with_current_values());
    }
    virtual ~CheckBeamParticle() {}

    virtual void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config);
        NeedRecombModel::configure(config);
        NeedParticleData::configure(config);
        // Same guard as CheckSTM_Michel.cxx:111: the mixin's type-only default
        // "DetectorVolumes" is not an instantiated component in the PR config.
        m_use_fiducial = !config["fiducial"].isNull();
        if (m_use_fiducial) NeedFiducial::configure(config);

        m_cfg = config;
        m_grouping_name = get<std::string>(config, "grouping", m_grouping_name);
        m_perf = get<bool>(config, "perf", m_perf);
        m_publish_nu_slots = get<bool>(config, "publish_nu_slots", m_publish_nu_slots);

        // e/cm, the TaggerCheckNeutrino convention (converted at the PA copy).
        m_mip_dqdx = get<double>(config, "mip_dqdx", m_mip_dqdx);
        m_mip_dqdx_median = get<double>(config, "mip_dqdx_median", m_mip_dqdx_median);

        // The beam window [low, high) on cluster_t0 (= the RAW matched flash
        // time, TaggerCheckNeutrino.cxx:486), internal units.  low >= high =
        // the gate is off, and off means this stage selects NOTHING: a
        // beam-particle stage must never fall back to "every bundle".
        m_beam_window_low = get<double>(config, "beam_window_low", m_beam_window_low);
        m_beam_window_high = get<double>(config, "beam_window_high", m_beam_window_high);

        // The nominal beam entry (cm, detector frame, T0-corrected x) and the
        // direction the particles travel.  doc pdvd/120 sec 1.
        auto read3 = [&](const char* key, std::vector<double>& out) {
            if (!config.isMember(key) || !config[key].isArray() || config[key].size() != 3) return;
            for (int i = 0; i < 3; ++i) out[i] = config[key][i].asDouble();
        };
        read3("beam_entry_point_cm", m_beam_entry_point_cm);
        read3("beam_dir", m_beam_dir);
        m_beam_entry_max_dist_cm = get<double>(config, "beam_entry_max_dist_cm", m_beam_entry_max_dist_cm);
        m_beam_dir_max_angle_deg = get<double>(config, "beam_dir_max_angle_deg", m_beam_dir_max_angle_deg);
        m_entry_tie_tol_cm = get<double>(config, "entry_tie_tol_cm", m_entry_tie_tol_cm);
        m_entry_snap_tol_cm = get<double>(config, "entry_snap_tol_cm", m_entry_snap_tol_cm);
        m_min_main_length_cm = get<double>(config, "min_main_length_cm", m_min_main_length_cm);
        m_improve_entry_vertex = get<bool>(config, "improve_entry_vertex", m_improve_entry_vertex);
        m_entry_fail_fallback_geo = get<bool>(config, "entry_fail_fallback_geo", m_entry_fail_fallback_geo);
        m_fv_tolerance.clear();
        if (config.isMember("fv_tolerance") && config["fv_tolerance"].isArray())
            for (const auto& t : config["fv_tolerance"]) m_fv_tolerance.push_back(t.asDouble());

        // TrackFitting parameters carried by TaggerCheckNeutrino's config
        // rather than by the runtime JSON (TaggerCheckNeutrino.cxx:3558-3572).
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
        cfg["fv_tolerance"] = Json::Value(Json::arrayValue);
        cfg["trackfitting_config_file"] = "";
        cfg["perf"] = m_perf;
        cfg["publish_nu_slots"] = m_publish_nu_slots;
        cfg["mip_dqdx"] = m_mip_dqdx;
        cfg["mip_dqdx_median"] = m_mip_dqdx_median;
        cfg["beam_window_low"] = m_beam_window_low;
        cfg["beam_window_high"] = m_beam_window_high;
        Json::Value ep(Json::arrayValue), bd(Json::arrayValue);
        for (int i = 0; i < 3; ++i) { ep.append(m_beam_entry_point_cm[i]); bd.append(m_beam_dir[i]); }
        cfg["beam_entry_point_cm"] = ep;
        cfg["beam_dir"] = bd;
        cfg["beam_entry_max_dist_cm"] = m_beam_entry_max_dist_cm;
        cfg["beam_dir_max_angle_deg"] = m_beam_dir_max_angle_deg;
        cfg["entry_tie_tol_cm"] = m_entry_tie_tol_cm;
        cfg["entry_snap_tol_cm"] = m_entry_snap_tol_cm;
        cfg["min_main_length_cm"] = m_min_main_length_cm;
        cfg["improve_entry_vertex"] = m_improve_entry_vertex;
        cfg["entry_fail_fallback_geo"] = m_entry_fail_fallback_geo;
        cfg["fit_blob_coverage"] = m_fit_blob_coverage;
        cfg["dqdx_fit_keep_all_points"] = m_dqdx_fit_keep_all_points;
        cfg["excl_t0_frame"] = m_excl_t0_frame;
        // the PR-partition knobs: null = ride the PatternAlgorithms C++ default
        for (const auto& k : pattern_knob_keys()) cfg[k] = Json::Value();
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const;

private:
    using Clock = std::chrono::steady_clock;
    using MS = std::chrono::duration<double, std::milli>;

    Configuration m_cfg;
    std::string m_grouping_name{"live"};
    bool m_perf{false};
    bool m_publish_nu_slots{true};
    bool m_use_fiducial{false};
    double m_mip_dqdx{50000.0};
    double m_mip_dqdx_median{43000.0};
    double m_beam_window_low{0.0};
    double m_beam_window_high{0.0};
    // doc pdvd/120 sec 1: the entry the three run-39305 beam-matched tracks
    // point at (x, y from the data; z the active face), and the GDML beam-plug
    // axis direction (protodunevd_v5_ggd.gdml), travel sense.
    std::vector<double> m_beam_entry_point_cm{110.0, 159.0, 0.6};
    std::vector<double> m_beam_dir{-0.095, -0.704, 0.704};
    double m_beam_entry_max_dist_cm{50.0};
    double m_beam_dir_max_angle_deg{-1.0};   // -1 = off
    double m_entry_tie_tol_cm{5.0};
    double m_entry_snap_tol_cm{10.0};
    double m_min_main_length_cm{0.0};
    bool m_improve_entry_vertex{false};
    bool m_entry_fail_fallback_geo{false};
    std::vector<double> m_fv_tolerance;
    double m_fit_blob_coverage{-1.0};
    bool m_dqdx_fit_keep_all_points{false};
    bool m_excl_t0_frame{false};
    std::shared_ptr<TrackFitting> m_track_fitter;
    mutable std::string m_evt_tag;

    // Everything the stage decided, one row on the main cluster ("beam_particle" PC).
    struct Record {
        int cluster_id{-1}, gid{-1}, n_in_window{0}, n_gids{0}, n_companions{0};
        double t0_us{0}, flash_pe{-1}, main_len_cm{0};
        Point nominal, entry, exit, entry_vtx, geo_vtx;
        double entry_dist{1e9}, entry_cos_beam{0}, entry_snap{1e9}, geo_entry_dist{-1};
        int entry_ok{0}, entry_split{0}, entry_tie_by_dir{0}, has_geo_vtx{0}, entry_vtx_id{-1};
        int n_showers{0}, match_isFC{0};
        double kine_reco_Enu{0};
    };

    static const std::vector<std::string>& pattern_knob_keys() {
        // The subset of TaggerCheckNeutrino's config keys the PR stages read
        // (same names, same units: TaggerCheckNeutrino.cxx:3107-3300), so the
        // PDVD job forwards its production partition knobs unchanged.  Verbatim
        // the CheckSTM_Michel list (fork by duplication, M10).
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
            "traj_cover_probe",
        };
        return keys;
    }

    // Copy the PR-partition knobs onto a PatternAlgorithms, reading the
    // config directly with the PA's own C++ defaults as fallback.  Unit
    // conventions are TaggerCheckNeutrino.cxx:3107-3300's: cm keys scaled by
    // units::cm, dQ/dx keys divided by units::cm, degrees/ratios/bools raw.
    // Verbatim CheckSTM_Michel::apply_pattern_knobs MINUS its
    // dqdx_skip_zero_dx override: this stage keeps TaggerCheckNeutrino's
    // kinematics behaviour unless the key is in the bag.
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
        B(pa.m_traj_cover_probe, "traj_cover_probe");
        // Process-wide state (TaggerCheckNeutrino.cxx): only written when the
        // job configures it, so a job that does not name the key leaves
        // whatever the process already had.
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
                std::cerr << "CheckBeamParticle: Cannot open config file: " << config_file
                          << " (not found on WIRECELL_PATH either)" << std::endl;
                return;
            }
            Json::Value root;
            Json::CharReaderBuilder builder;
            std::string errs;
            if (!Json::parseFromStream(builder, file, &root, &errs)) {
                std::cerr << "CheckBeamParticle: Failed to parse JSON: " << errs << std::endl;
                return;
            }
            for (const auto& param_name : root.getMemberNames()) {
                if (param_name.substr(0, 1) == "_") continue;  // comments
                try {
                    double value = root[param_name].asDouble();
                    m_track_fitter->set_parameter(param_name, value);
                } catch (const std::exception& e) {
                    std::cerr << "CheckBeamParticle: Failed to set parameter " << param_name
                              << ": " << e.what() << std::endl;
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "CheckBeamParticle: Exception loading config: " << e.what() << std::endl;
        }
    }

    // Nearest graph vertex of `cluster` to `pt`; if none within tol, split the
    // nearest segment there (PR::break_segment) so the entry becomes a vertex
    // the orientation walk and the PF tree can use.  Duplicate of
    // CheckSTM_Michel::anchor_vertex (:1566-1604), M10.
    VertexPtr anchor_vertex(PatternAlgorithms& pa, Graph& g, Cluster& cluster, const Point& pt,
                            double tol, double& out_dis, bool& split) const {
        split = false;
        auto [vtx, dis] = pa.closest_cluster_vertex(g, cluster, pt);
        out_dis = dis;
        if (vtx && dis <= tol) return vtx;
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
                    // break_segment stamps the cluster itself (doc sbnd_xin/pr/143;
                    // see the CheckSTM_Michel note on a clusterless vertex).
                    out_dis = best_d;
                    split = true;
                    return nvtx;
                }
            } catch (const std::exception& e) {
                SPDLOG_LOGGER_WARN(s_log, "{}anchor_vertex: break_segment threw: {}", m_evt_tag, e.what());
            }
        }
        return vtx;  // may be null
    }

    static Point vertex_point(const VertexPtr& v) {
        return v->fit().valid() ? v->fit().point : v->wcpt().point;
    }

    void persist(Cluster& cluster, const Record& r) const {
        using WireCell::PointCloud::Array;
        using WireCell::PointCloud::Dataset;
        const double cm = units::cm;
        std::map<std::string, Array> a;
        auto I1 = [&](const char* k, int v)    { a.emplace(k, Array(std::vector<int>{v})); };
        auto D1 = [&](const char* k, double v) { a.emplace(k, Array(std::vector<double>{v})); };
        auto P3 = [&](const char* k, const Point& p) {
            a.emplace(std::string(k) + "_x", Array(std::vector<double>{p.x() / cm}));
            a.emplace(std::string(k) + "_y", Array(std::vector<double>{p.y() / cm}));
            a.emplace(std::string(k) + "_z", Array(std::vector<double>{p.z() / cm}));
        };
        I1("cluster_id", r.cluster_id); I1("gid", r.gid); D1("t0_us", r.t0_us); D1("flash_pe", r.flash_pe);
        I1("n_in_window", r.n_in_window); I1("n_gids", r.n_gids); I1("n_companions", r.n_companions);
        D1("main_len_cm", r.main_len_cm);
        P3("nominal", r.nominal); P3("entry", r.entry); P3("exit", r.exit);
        D1("entry_dist_cm", r.entry_dist / cm); D1("entry_cos_beam", r.entry_cos_beam);
        I1("entry_tie_by_dir", r.entry_tie_by_dir); I1("entry_ok", r.entry_ok);
        P3("entry_vtx", r.entry_vtx); I1("entry_vtx_id", r.entry_vtx_id);
        D1("entry_snap_cm", r.entry_snap / cm); I1("entry_split", r.entry_split);
        I1("has_geo_vtx", r.has_geo_vtx); P3("geo_vtx", r.geo_vtx); D1("geo_entry_dist_cm", r.geo_entry_dist / cm);
        I1("n_showers", r.n_showers); D1("kine_reco_Enu", r.kine_reco_Enu); I1("match_isFC", r.match_isFC);
        cluster.local_pcs()["beam_particle"] = Dataset(a);
    }
};

void CheckBeamParticle::visit(Ensemble& ensemble) const
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

    // ---- 1. the beam bundle ---------------------------------------------
    std::vector<BeamBundleRow> rows;
    for (auto* cluster : grouping.children()) {
        if (!cluster->get_flag(Flags::main_cluster)) continue;
        rows.push_back({cluster->get_cluster_id(), cluster->get_scalar<int>("matched_flash_gid", -1),
                        cluster->get_cluster_t0(), cluster->get_length()});
    }
    std::set<int> gids;   // int-keyed => sorted
    for (const auto& r : rows) {
        if (r.gid >= 0 && r.t0 >= m_beam_window_low && r.t0 < m_beam_window_high) gids.insert(r.gid);
    }
    std::vector<BeamFlashInfo> flashes;
    for (int gid : gids) {
        auto fl = grouping.flash_by_gid(gid);
        flashes.push_back({gid, fl ? fl.value() : -1.0, static_cast<bool>(fl)});
    }
    const auto pick = beam_particle_pick_bundle(rows, flashes, m_beam_window_low, m_beam_window_high);
    SPDLOG_LOGGER_INFO(s_log,
        "{}CheckBeamParticle: {} main(s), {} in window [{:.3f}, {:.3f}) us over {} bundle(s); picked gid {} ({})",
        m_evt_tag, rows.size(), pick.n_in_window, m_beam_window_low / units::us, m_beam_window_high / units::us,
        pick.n_gids, pick.gid, pick.why);
    if (pick.gid < 0) {
        SPDLOG_LOGGER_INFO(s_log, "{}CheckBeamParticle: no beam bundle; nothing to reconstruct", m_evt_tag);
        return;
    }
    double flash_pe = -1;
    for (const auto& f : flashes) if (f.gid == pick.gid && f.valid) flash_pe = f.pe;

    // ---- 2. the main cluster: the bundle member closest to the nominal entry
    const Point nominal(m_beam_entry_point_cm[0] * units::cm, m_beam_entry_point_cm[1] * units::cm,
                        m_beam_entry_point_cm[2] * units::cm);
    const Vector beam_dir(m_beam_dir[0], m_beam_dir[1], m_beam_dir[2]);
    std::vector<Cluster*> bundle;
    for (auto* c : grouping.children()) {
        if (c->get_scalar<int>("matched_flash_gid", -1) != pick.gid) continue;
        if (c->npoints() == 0) continue;
        bundle.push_back(c);
    }
    std::sort(bundle.begin(), bundle.end(), [](const Cluster* a, const Cluster* b) { return a->ident() < b->ident(); });
    std::vector<BeamMainCand> cands;
    for (auto* c : bundle) {
        const auto [cp, blob] = c->get_closest_point_blob(nominal);
        cands.push_back({c->get_cluster_id(), (cp - nominal).magnitude(), c->get_length()});
        SPDLOG_LOGGER_INFO(s_log, "{}CheckBeamParticle: gid {} member cluster {} L {:.1f} cm closest point ({:.1f}, {:.1f}, {:.1f}) d_entry {:.1f} cm",
                           m_evt_tag, pick.gid, c->get_cluster_id(), c->get_length() / units::cm,
                           cp.x() / units::cm, cp.y() / units::cm, cp.z() / units::cm, cands.back().dist / units::cm);
    }
    const int imain = beam_particle_pick_main(cands, m_min_main_length_cm * units::cm);
    if (imain < 0) {
        SPDLOG_LOGGER_INFO(s_log, "{}CheckBeamParticle: bundle gid {} has no usable cluster ({} members); nothing to reconstruct",
                           m_evt_tag, pick.gid, bundle.size());
        return;
    }
    Cluster* main = bundle[imain];
    std::vector<Cluster*> companions;
    for (auto* c : bundle) if (c != main) companions.push_back(c);
    for (auto* c : bundle) c->set_scalar<int>("beam_particle_gid", pick.gid);
    main->set_scalar<int>("beam_particle_main", 1);

    Record rec;
    rec.cluster_id = main->get_cluster_id();
    rec.gid = pick.gid;
    rec.n_in_window = pick.n_in_window;
    rec.n_gids = pick.n_gids;
    rec.n_companions = static_cast<int>(companions.size());
    rec.t0_us = main->get_cluster_t0() / units::us;
    rec.flash_pe = flash_pe;
    rec.main_len_cm = main->get_length() / units::cm;
    rec.nominal = nominal;
    SPDLOG_LOGGER_INFO(s_log,
        "{}CheckBeamParticle: gid {} t0 {:.3f} us pe {:.0f}: main cluster {} d_entry {:.1f} cm L {:.1f} cm; {} companion(s)",
        m_evt_tag, pick.gid, rec.t0_us, flash_pe, rec.cluster_id, cands[imain].dist / units::cm, rec.main_len_cm,
        companions.size());

    // ---- 3. fitter + graph, as CheckSTM_Michel.cxx:2716-2735 ------------
    auto t0 = Clock::now();
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

    IndexedVertexSet vertices_in_long_muon;
    IndexedSegmentSet segments_in_long_muon;
    VertexPtr main_vertex = nullptr;
    ClusterVertexMap map_cluster_main_vertices;

    // ---- 4. the main cluster's PR, TaggerCheckNeutrino.cxx:3657-3699 -----
    const bool ok_main = pa.find_proto_vertex(g, *main, *tf, m_dv, true, 2, true, particle_data());
    if (!ok_main) {
        SPDLOG_LOGGER_WARN(s_log, "{}CheckBeamParticle: find_proto_vertex failed on main cluster {}; nothing published",
                           m_evt_tag, rec.cluster_id);
        persist(*main, rec);
        return;
    }
    pa.clustering_points(g, *main, m_dv);
    pa.separate_track_shower(g, *main);
    pa.determine_direction(g, *main, particle_data(), m_recomb_model);
    pa.shower_determining_in_main_cluster(g, *main, particle_data(), m_recomb_model, m_dv);
    // determine_main_vertex is kept for its graph work (examine_structure_*,
    // examine_vertices); the vertex it chooses is recorded as "what the
    // neutrino chain would pick" and then REPLACED by the entry below.
    pa.determine_main_vertex(g, *main, main_vertex, vertices_in_long_muon, segments_in_long_muon, *tf, m_dv,
                             particle_data(), m_recomb_model);
    pa.reassociate_cluster_orphans(g, *main, m_dv);
    VertexPtr geo_v = main_vertex;
    if (main_vertex != nullptr) {
        map_cluster_main_vertices[main] = main_vertex;
        main_vertex = nullptr;
    }
    if (m_perf) SPDLOG_LOGGER_DEBUG(s_log, "{}CheckBeamParticle timing: main PR took {} ms", m_evt_tag, MS(Clock::now() - t0).count());
    t0 = Clock::now();

    // ---- 5. companions, TaggerCheckNeutrino.cxx:3707-3756 ----------------
    if (!companions.empty()) {
        for (auto* cluster : companions) {
            if (cluster->get_length() > 6 * units::cm) {
                pa.find_proto_vertex(g, *cluster, *tf, m_dv, true, 2, false);
            }
            else {
                if (!pa.find_proto_vertex(g, *cluster, *tf, m_dv, false, 1, false)) {
                    pa.init_point_segment(g, *cluster, *tf, m_dv);
                }
            }
            pa.clustering_points(g, *cluster, m_dv);
            pa.separate_track_shower(g, *cluster);
            pa.determine_direction(g, *cluster, particle_data(), m_recomb_model);
            pa.shower_determining_in_main_cluster(g, *cluster, particle_data(), m_recomb_model, m_dv);
            pa.determine_main_vertex(g, *cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, *tf, m_dv,
                                     particle_data(), m_recomb_model);
            pa.reassociate_cluster_orphans(g, *cluster, m_dv);
            if (main_vertex != nullptr) {
                map_cluster_main_vertices[cluster] = main_vertex;
                main_vertex = nullptr;
            }
        }
        std::vector<Cluster*> all_clusters;
        all_clusters.push_back(main);
        all_clusters.insert(all_clusters.end(), companions.begin(), companions.end());
        pa.deghosting(g, map_cluster_main_vertices, all_clusters, *tf, m_dv);
        if (m_perf) SPDLOG_LOGGER_DEBUG(s_log, "{}CheckBeamParticle timing: companions + deghosting took {} ms", m_evt_tag, MS(Clock::now() - t0).count());
        t0 = Clock::now();
    }

    // ---- 6. the ENTRY vertex (replaces TaggerCheckNeutrino.cxx:3765-3813) --
    const auto ends = main->get_main_axis_points();
    const auto choice = beam_particle_choose_entry(ends.first, ends.second, nominal, beam_dir,
                                                   m_entry_tie_tol_cm * units::cm);
    rec.entry = choice.entry;
    rec.exit = choice.exit;
    rec.entry_dist = choice.dist;
    rec.entry_cos_beam = choice.cos_beam;
    rec.entry_tie_by_dir = choice.tie_by_dir ? 1 : 0;
    if (geo_v) {
        rec.has_geo_vtx = 1;
        rec.geo_vtx = vertex_point(geo_v);
        rec.geo_entry_dist = (rec.geo_vtx - choice.entry).magnitude();
    }
    bool entry_ok = choice.dist <= m_beam_entry_max_dist_cm * units::cm;
    if (entry_ok && m_beam_dir_max_angle_deg >= 0) {
        entry_ok = choice.cos_beam >= std::cos(m_beam_dir_max_angle_deg * M_PI / 180.0);
    }
    VertexPtr final_main_vertex = nullptr;
    if (entry_ok) {
        bool split = false;
        final_main_vertex = anchor_vertex(pa, g, *main, choice.entry, m_entry_snap_tol_cm * units::cm,
                                          rec.entry_snap, split);
        rec.entry_split = split ? 1 : 0;
        if (!final_main_vertex) {
            SPDLOG_LOGGER_WARN(s_log, "{}CheckBeamParticle: cluster {} entry ({:.1f}, {:.1f}, {:.1f}) cm has no graph vertex within {:.1f} cm (nearest {:.1f} cm)",
                               m_evt_tag, rec.cluster_id, choice.entry.x() / units::cm, choice.entry.y() / units::cm,
                               choice.entry.z() / units::cm, m_entry_snap_tol_cm, rec.entry_snap / units::cm);
            entry_ok = false;
        }
    }
    else {
        SPDLOG_LOGGER_WARN(s_log, "{}CheckBeamParticle: cluster {} axis end ({:.1f}, {:.1f}, {:.1f}) cm is {:.1f} cm from the nominal entry (max {:.1f}), cos_beam {:.2f}: not the beam particle",
                           m_evt_tag, rec.cluster_id, choice.entry.x() / units::cm, choice.entry.y() / units::cm,
                           choice.entry.z() / units::cm, choice.dist / units::cm, m_beam_entry_max_dist_cm, choice.cos_beam);
    }
    rec.entry_ok = entry_ok ? 1 : 0;
    if (!entry_ok) {
        if (m_entry_fail_fallback_geo && geo_v) {
            final_main_vertex = geo_v;
            SPDLOG_LOGGER_INFO(s_log, "{}CheckBeamParticle: entry_fail_fallback_geo: using the neutrino-chain vertex of cluster {}",
                               m_evt_tag, rec.cluster_id);
        }
        else {
            persist(*main, rec);
            SPDLOG_LOGGER_INFO(s_log, "{}CheckBeamParticle: cluster {} gid {} entry_ok 0; nothing published | {:.0f} ms",
                               m_evt_tag, rec.cluster_id, rec.gid, MS(Clock::now() - t_total).count());
            return;
        }
    }
    final_main_vertex->set_flags(VertexFlags::kNeutrinoVertex);
    map_cluster_main_vertices[main] = final_main_vertex;

    // ---- 7. refinement, TaggerCheckNeutrino.cxx:3851-3930 (subset) -------
    if (m_improve_entry_vertex) {
        pa.improve_vertex(g, *main, final_main_vertex, vertices_in_long_muon, segments_in_long_muon,
                          *tf, m_dv, particle_data(), m_recomb_model, true, true);
        final_main_vertex->set_flags(VertexFlags::kNeutrinoVertex);   // the pointer may have moved
        map_cluster_main_vertices[main] = final_main_vertex;
    }
    pa.clustering_points(g, *main, m_dv);
    pa.reassociate_cluster_orphans(g, *main, m_dv);
    for (auto* cluster : companions) pa.reassociate_cluster_orphans(g, *cluster, m_dv);
    // examine_direction runs last and has the final word on segment
    // orientations relative to the main vertex -- here the entry.
    pa.examine_direction(g, final_main_vertex, final_main_vertex, vertices_in_long_muon, segments_in_long_muon,
                         particle_data(), m_recomb_model, true);
    rec.entry_vtx = vertex_point(final_main_vertex);
    rec.entry_vtx_id = rec.cluster_id * 1000 + static_cast<int>(final_main_vertex->get_graph_index());
    SPDLOG_LOGGER_INFO(s_log,
        "{}CheckBeamParticle: cluster {} entry ({:.1f}, {:.1f}, {:.1f}) cm d_nominal {:.1f} cm cos_beam {:.2f} -> vertex {} ({:.1f}, {:.1f}, {:.1f}) snap {:.1f} cm split {} | neutrino-chain vertex {} ({:.1f}, {:.1f}, {:.1f}) d {:.1f} cm",
        m_evt_tag, rec.cluster_id, rec.entry.x() / units::cm, rec.entry.y() / units::cm, rec.entry.z() / units::cm,
        rec.entry_dist / units::cm, rec.entry_cos_beam, rec.entry_vtx_id,
        rec.entry_vtx.x() / units::cm, rec.entry_vtx.y() / units::cm, rec.entry_vtx.z() / units::cm,
        rec.entry_snap / units::cm, rec.entry_split,
        rec.has_geo_vtx, rec.geo_vtx.x() / units::cm, rec.geo_vtx.y() / units::cm, rec.geo_vtx.z() / units::cm,
        rec.geo_entry_dist / units::cm);
    if (m_perf) SPDLOG_LOGGER_DEBUG(s_log, "{}CheckBeamParticle timing: entry + examine_direction took {} ms", m_evt_tag, MS(Clock::now() - t0).count());
    t0 = Clock::now();

    // ---- 8. EM showers, TaggerCheckNeutrino.cxx:3949-4001 ----------------
    int acc_segment_id = 0;
    IndexedShowerSet pi0_showers;
    ShowerIntMap map_shower_pio_id;
    std::map<int, std::vector<ShowerPtr>> map_pio_id_showers;
    std::map<int, std::pair<double, int>> map_pio_id_mass;
    std::map<int, std::pair<int, int>> map_pio_id_saved_pair;
    Pi0KineFeatures pio_kine{};
    ShowerVertexMap map_vertex_in_shower;
    ShowerSegmentMap map_segment_in_shower;
    VertexShowerSetMap map_vertex_to_shower;
    ClusterPtrSet used_shower_clusters;
    IndexedShowerSet showers;

    pa.demote_cross_cluster_straight_stems(g, final_main_vertex, particle_data(), m_recomb_model);
    pa.shower_clustering_with_nv(acc_segment_id, pi0_showers,
                                 map_shower_pio_id, map_pio_id_showers,
                                 map_pio_id_mass, map_pio_id_saved_pair,
                                 pio_kine,
                                 vertices_in_long_muon, segments_in_long_muon,
                                 g, final_main_vertex, showers,
                                 main, companions,
                                 map_cluster_main_vertices,
                                 map_vertex_in_shower, map_segment_in_shower,
                                 map_vertex_to_shower, used_shower_clusters,
                                 *tf, m_dv, particle_data(),
                                 m_recomb_model);
    pa.reconcile_particle_flags(g, final_main_vertex, showers,
                                map_vertex_in_shower, map_segment_in_shower,
                                map_vertex_to_shower, map_shower_pio_id,
                                particle_data(), m_recomb_model);
    rec.n_showers = static_cast<int>(showers.size());
    if (m_perf) SPDLOG_LOGGER_DEBUG(s_log, "{}CheckBeamParticle timing: shower clustering took {} ms", m_evt_tag, MS(Clock::now() - t0).count());
    t0 = Clock::now();

    // ---- 9. no taggers: TaggerInfo at its defaults + match_isFC (:4322-4337)
    TaggerInfo tagger_info;
    pa.init_tagger_info(tagger_info);
    {
        auto fc_result = Facade::cluster_fc_check(*main, m_dv, m_use_fiducial ? m_fiducial : nullptr, m_fv_tolerance);
        tagger_info.match_isFC = fc_result.is_fc ? 1.0f : 0.0f;
        rec.match_isFC = fc_result.is_fc ? 1 : 0;
    }

    // ---- 10. kinematics, TaggerCheckNeutrino.cxx:4348-4354 ---------------
    std::set<int> dropped_sat_ids;
    KineInfo kine_info = pa.fill_kine_tree(final_main_vertex, showers, pio_kine, g, *tf, m_dv,
                                           nullptr,   // no clus_geom_helper on PDVD
                                           particle_data(), m_recomb_model, pi0_showers, &dropped_sat_ids);
    kine_info.cluster_id = rec.cluster_id;
    kine_info.matched_flash_gid = rec.gid;
    kine_info.nu_index = 0;
    kine_info.has_vertex = 1;
    rec.kine_reco_Enu = kine_info.kine_reco_Enu;

    // ---- 11. publish, TaggerCheckNeutrino.cxx:4362-4448 ------------------
    tf->set_pi0_data(pi0_showers, map_shower_pio_id, map_pio_id_showers, map_pio_id_mass);
    tf->set_dropped_satellite_shower_ids(std::move(dropped_sat_ids));
    tf->set_main_vertex(final_main_vertex);
    tf->set_showers(showers);
    tf->set_kine_info(kine_info);
    tf->set_tagger_info(tagger_info);
    tf->assemble_fitted_charge_2d();
    grouping.set_track_fitting(tf);
    if (m_publish_nu_slots) grouping.set_track_fitting("nu0", tf);
    persist(*main, rec);

    SPDLOG_LOGGER_INFO(s_log,
        "{}CheckBeamParticle: cluster {} gid {} entry_ok {} | {} shower(s) Enu {:.1f} MeV vertex ({:.1f}, {:.1f}, {:.1f}) cm isFC {} | {:.0f} ms",
        m_evt_tag, rec.cluster_id, rec.gid, rec.entry_ok, rec.n_showers, rec.kine_reco_Enu,
        kine_info.kine_nu_x_corr, kine_info.kine_nu_y_corr, kine_info.kine_nu_z_corr, rec.match_isFC,
        MS(Clock::now() - t_total).count());
}
