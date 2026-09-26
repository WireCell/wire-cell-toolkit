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
        m_entry_long_muon_absorb = get<bool>(config, "entry_long_muon_absorb", m_entry_long_muon_absorb);
        m_kine_charge_t0_frame = get<bool>(config, "kine_charge_t0_frame", m_kine_charge_t0_frame);
        m_fv_tolerance.clear();
        if (config.isMember("fv_tolerance") && config["fv_tolerance"].isArray())
            for (const auto& t : config["fv_tolerance"]) m_fv_tolerance.push_back(t.asDouble());

        // TrackFitting parameters carried by TaggerCheckNeutrino's config
        // rather than by the runtime JSON (TaggerCheckNeutrino.cxx:3558-3572).
        m_fit_blob_coverage = get<double>(config, "fit_blob_coverage", m_fit_blob_coverage);
        m_dqdx_fit_keep_all_points = get<bool>(config, "dqdx_fit_keep_all_points", m_dqdx_fit_keep_all_points);
        m_excl_t0_frame = get<bool>(config, "excl_t0_frame", m_excl_t0_frame);

        read_neutrino_knobs(config);

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
        cfg["entry_long_muon_absorb"] = m_entry_long_muon_absorb;
        cfg["kine_charge_t0_frame"] = m_kine_charge_t0_frame;
        cfg["fit_blob_coverage"] = m_fit_blob_coverage;
        cfg["dqdx_fit_keep_all_points"] = m_dqdx_fit_keep_all_points;
        cfg["excl_t0_frame"] = m_excl_t0_frame;
        default_neutrino_knobs(cfg);
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
    bool m_entry_long_muon_absorb{false};   // doc pdvd/120 sec 9.3; true = neutrino-chain flood-fill
    // doc pdvd/120 sec 9.4: KineChargeOptions::t0_frame.  This stage is
    // PDVD-only and new, so ITS default is on (the PatternAlgorithms default
    // stays false: every other consumer is byte-identical).  Without it every
    // EM shower of the beam bundle has kine_charge = 0 and the Bee tree's
    // em_ke_min floor hides it.
    bool m_kine_charge_t0_frame{true};
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
        int n_showers{0}, match_isFC{0}, entry_long_muon_nseg{0};
        double kine_reco_Enu{0};
    };

    // ---- doc pdvd/120 sec 9: the neutrino PR knob set ----------------------
    // The first cut of this stage read only the PR-partition subset that
    // CheckSTM_Michel reads (79 keys) and the PDVD job forwarded only that
    // subset, so shower clustering, the cross-cluster bridges/continuations,
    // the pi0 finder and the kine accounting all ran at their uBooNE-legacy
    // C++ defaults while TaggerCheckNeutrino ran with ~250 production keys.
    // The Bee particle flow of the beam bundle then held the entry muon and
    // little else (run 39305 evt 157312: 20 companion clusters fitted, 3 in
    // the tree).  The stage now carries TaggerCheckNeutrino's members,
    // configure() reads and PatternAlgorithms copy VERBATIM (generated from
    // the TCN source by pdvd/kaon/gen_beam_knobs.py; fork by duplication,
    // M10), MINUS the families this stage never runs: DL/SCN vertex + dual
    // chain (dl_*, dual_chain_*), candidate selection (nu_*, beam_window_* is
    // its own), cosmic taggers + BDT features (cosmic_*, ssm_*, sp_*,
    // muon_dqdx_curve, tagger_ordered_segment_sets, stem_endpoint_wcpt_parity,
    // broken_muon_cluster_id_count, neutrino_type_bitmask), MCS (mcs_*), the
    // main-vertex kink/junction snaps (the entry is a geometric fact:
    // vertex_kink_snap/vks_*, vertex_junction_snap/vjs_*), the SBND cathode
    // bridge (long_muon_cathode_bridge*) and the diagnostic probes
    // (rough_path_probe, sgp_edge_probe, vertex_scoreboard, dl_vtx_harvest,
    // kine_continuation_debug).  Same key names, same units (cm/deg/MeV at
    // the config, converted at the PA copy exactly as TaggerCheckNeutrino.cxx
    // :3104-3553 does), same C++ defaults (TaggerCheckNeutrino.h), so the
    // PDVD job hands this stage the neutrino node's compiled data.
    bool                 m_dir_weak_use_score{false};
    bool                 m_proton_dir_vote{false};
    double               m_proton_dir_score_max{0.25};
    double               m_proton_dir_asym_min{1.3};
    bool                 m_endpoint_trim_retry{false};
    double               m_fit_vertex_min_seg_length{0};
    bool                 m_mvfit_robust{false};
    bool                 m_mvfit_main_only{true};
    double               m_mvfit_min_len{10.0};  // cm
    double               m_mvfit_rin_margin{2.0};  // cm
    double               m_mvfit_rout_frac{0.5};
    double               m_mvfit_rout_min{9.0};  // cm
    double               m_mvfit_rout_max{18.0};  // cm
    double               m_mvfit_angle{20.0};  // deg
    int                  m_mvfit_min_pts{5};
    double               m_mvfit_min_aniso{3.0};
    double               m_mvfit_prior_range{1.0};  // cm
    double               m_cathode_x{0};
    double               m_cathode_kink_xcut{0};
    double               m_cathode_wide_kink_angle{0};  // deg
    double               m_cathode_wide_kink_skirt{3};  // cm
    double               m_cathode_wide_kink_baseline{15};  // cm
    bool                 m_two_end_break{false};
    double               m_teb_min_len{10};  // cm
    double               m_teb_min_arm{1.8};  // cm
    int                  m_teb_min_arm_pts{4};
    double               m_teb_stub_max{4};  // cm
    double               m_teb_accept_range{15};  // cm
    double               m_teb_rise_r1{1.3};
    double               m_teb_rise_r2{1.15};
    double               m_teb_abs_end_min{1.7};  // x mip_dqdx_median
    double               m_teb_dip_floor{0.6};  // x mip_dqdx_median; dips below = instrumental
    double               m_teb_score_cap_r1{0.6};
    double               m_teb_score_cap_r2{0.9};
    double               m_teb_turn_angle{25};  // deg; <= 0 disables route R2
    double               m_teb_turn_baseline{35};  // cm
    double               m_teb_turn_skirt{3};  // cm
    double               m_teb_turn_min_arm_frac{0.0};
    double               m_teb_bragg_veto_turn{0.0};  // deg; <= 0 disables the R2 bragg veto
    bool                 m_kink_walk_dqdx_stop{false};
    double               m_kink_dqdx_hot_ratio{1.7};
    bool                 m_kink_break_protect{false};
    bool                 m_esva_ignore_empty_2d{false};
    bool                 m_main_vertex_graph_audit{false};
    double               m_mvga_radius{15.0};  // cm
    double               m_mvga_dup_tol{1.4};  // cm
    double               m_mvga_dup_frac{0.7};  // fraction
    double               m_mvga_dup_angle{20.0};  // deg; <= 0 disables the op1 near-parallel guard
    double               m_mvga_bridge_mip{0.5};  // x mip_dqdx_median (measured: 268067 bridge 0.436, real track 1.29)
    double               m_mvga_reconnect{5.0};  // cm
    double               m_mvga_stub{2.0};  // cm
    int                  m_mvga_stub_pts{4};  // valid fit points
    double               m_mvga_reseat_angle{150.0};  // deg
    double               m_mvga_satellite{0};  // cm; 0 = main-vertex-only op3 scope (round 2), byte-identical (round 3, doc pr/51)
    bool                 m_mvga_interposed{false};  // op3 interposed-stub absorb at the main-vertex anchor (doc pr/85); false = terminal-only, byte-identical
    double               m_mvga_interposed_angle{150.0};  // deg; far-end collinearity gate for the interposed absorb
    double               m_mvga_interposed_len{0.0};  // cm; interposed-splice candidate ceiling (doc pr/86); 0 = use mvga_stub, byte-identical
    double               m_mvga_sat_dup_frac{0.0};  // fraction; satellite-anchor op3 overlap threshold (doc pr/86); 0 = use mvga_dup_frac, byte-identical
    bool                 m_mvga_interposed_deg1{false};  // op3 interposed splice at degree-1 main anchors (doc pr/86); false = byte-identical
    double               m_mvga_splice_straighten{0.0};  // cm; op3 post-carry straighten reach past the junction (doc pr/86 round 2); 0 = concatenation verbatim, byte-identical
    double               m_mvga_approach_collapse{0.0};  // cm; op3.5 junction-collapse radius around the main vertex (doc pr/86 round 2); 0 = pass skipped, byte-identical
    double               m_mvga_straighten_radius{0.0};  // cm; R1/R2 straight-chain charge-veto radius (doc pr/86 round 2); 0 = prototype 0.2 cm; inert unless straighten/collapse on
    double               m_mvga_op1_radius{0.0};  // cm; op1-only scope radius (doc pr/83 r3); 0 = use mvga_radius, -1 = unscoped, byte-identical at 0
    double               m_mvga_op1_dup_frac{0.0};  // fraction; op1-only overlap threshold (doc pr/83 r3); 0 = use mvga_dup_frac, byte-identical
    bool                 m_mvga_op1_post{false};  // post-op3 duplicate-corridor pass incl. created segments (doc pr/83 r3 class A); false = byte-identical
    bool                 m_swap_orphan_dup_audit{false};  // dup-audit the abandoned main cluster inside swap_main_cluster (doc pr/83 r3 Mechanism C); false = byte-identical
    double               m_mvga_proj_dup_frac{0.0};  // 2nd-best per-view overlap threshold for the projective dup collapse (doc pr/83 r4); 0 = disabled, byte-identical
    double               m_mvga_proj_dqdx_ratio{0.4};  // stem dQ/dx asymmetry gate for the same pass (doc pr/83 r4); inert while frac == 0
    double               m_mvga_proj_angle{0.0};  // deg; op1-proj chord-angle ceiling (doc pr/83 r4b); 0 = use mvga_dup_angle, byte-identical
    double               m_mvga_ac_chord_max{0.0};  // cm; op3.5 replacement-chord length cap (doc pr/99 round 2); 0 = no cap, byte-identical
    bool                 m_mvga_ac_no_cascade{false};  // op3.5: skip candidates touching `created` products (doc pr/99 round 2); false = byte-identical
    double               m_mvga_passthru{0.0};  // cm; op0 pass-through split radius (doc pr/103, SBND 18255-405707): a prong of a junction J within this radius of the main vertex whose interior passes through the main vertex is split there; 0 = off, byte-identical
    double               m_mvga_interposed_fallback_min_angle{0.0};  // deg; op3 fallback: min measured far_angle (doc pr/103); 0 = no floor beyond 'measured'; inert while the fallback is off
    bool                 m_mvga_interposed_fallback{false};  // op3: when the interposed far-angle gate declines at the main anchor, try the per-prong charge-verified straighten splice instead (doc pr/103); false = byte-identical
    double               m_mvga_passthru_tol{1.0};  // cm; op0 miss tolerance -- max distance from the main-vertex wcpt to the passing prong (doc pr/103); inert while m_mvga_passthru == 0
    double               m_mvga_dup_starved_asym{0.0};  // pair min/max dQ/dx asymmetry gate; op1-post angle-decline starved-member override (doc pr/99 round 2); 0 = off, byte-identical
    double               m_mvga_dup_starved_mip{0.0};  // absolute cap on the loser, ratio vs mip median; same override (doc pr/99 round 2); 0 = off, byte-identical
    double               m_mvga_dup_starved_span{0.0};  // pair min/max length comparability floor; same override (doc pr/99 round 2); 0 = no span test
    double               m_shower_topo_demote_len{0};
    bool                 m_fit_blob_coverage_defer{false};
    bool                 m_fit_exclusion{false};  // P1
    double               m_graph_endpoint_tol{0.3};  // cm
    bool                 m_oov_prototype_parity{false};  // F2 (all three sites)
    bool                 m_first_seg_local_pca{true};  // P2
    bool                 m_other_seg_relaxed_accept{true};  // P4
    bool                 m_other_seg_empty_2d_guard{false};  // doc pr/45 -- -1 empty-2D-tree sentinel guard
    bool                 m_other_seg_keep_isolated{false};
    int                  m_other_seg_keep_isolated_min_points{25};
    double               m_other_seg_keep_isolated_min_length{3.0};  // cm; scaled at copy
    double               m_other_seg_keep_isolated_len_admit{0.0};  // cm; scaled at copy
    double               m_iso_snap_min_dir_mag{10.0};  // cm; scaled at copy
    bool                 m_assoc_full_recluster{false};
    bool                 m_assoc_reassign_orphans{false};
    bool                 m_assoc_clear_on_merge{false};
    bool                 m_shower_topo_proto_dir{false};  // pr/31 F2
    bool                 m_vertex_dir_use_fit_point{false};  // pr/32 F1 (was P1)
    bool                 m_shower_traj_recheck_parity{false};  // pr/32 F2 (was P3)
    bool                 m_main_vertex_require_descriptor{false};  // pr/32 F3 (was P7)
    bool                 m_main_vertex_candidate_flag{false};  // pr/32 F4 (was P12)
    bool                 m_cont_muon_dir3_30cm{false};  // pr/31 F5 (was P6)
    bool                 m_track_comp_empty_abstain{false};  // pr/31 F6 (was P7)
    bool                 m_shower_topo_reset{false};  // pr/31 F3 (was P13)
    bool                 m_reclass_preserve_4mom{false};  // pr/31 F1 (was P1+P3a+P4)
    bool                 m_reclass_never_computed_ke_floor{false};  // pr/40 round 2 F6
    bool                 m_dir_track_median_local{false};  // pr/31 F4 (was P8)
    bool                 m_examine_showers_vertex_by_index{false};  // pr/31 F7 (was P5) -- stays OFF pending pr/30 F4
    bool                 m_iso_endpoint{false};
    double               m_iso_endpoint_min_length{40};  // cm
    double               m_iso_endpoint_max_xext{25};  // cm
    double               m_iso_endpoint_xext_frac{0.35};
    double               m_iso_endpoint_xext_quantile{0.02};
    double               m_iso_endpoint_tube_radius{4};  // cm
    double               m_iso_endpoint_min_aspect{0.12};
    bool                 m_traj_cover_probe{false};
    int                  m_pr_find_other_rounds{0};
    bool                 m_v3_extension_guard{false};
    double               m_v3_extension_min_gain{-1.0};  // cm
    bool                 m_es3_stub_guard{false};
    double               m_es3sg_stub_max{7.0};  // cm; scaled at copy
    double               m_es3sg_len_ratio{2.0};
    double               m_es3sg_ang3_min{15.0};  // degrees
    double               m_es3sg_ang_ratio{1.0};
    bool                 m_es3sg_require_terminal{true};
    double               m_vertex_z_prior_scale{200};
    double               m_kine_fudge_factor{0.95};
    double               m_kine_recom_factor{0.7};
    double               m_kine_shower_fudge_factor{0.8};
    double               m_kine_shower_recom_factor{0.5};
    double               m_kine_proton_recom_factor{0.35};
    double               m_kine_plane_asym_switch{0.04};
    double               m_kine_w_value{23.6};  // eV per electron-ion pair
    bool                 m_kine_shower_pdg_live{false};
    double               m_steiner_gap_penalty{0};
    double               m_sgp_dead_alpha{0.25};  // dead-sample weight in bad_fraction
    double               m_sgp_min_edge{0.5};  // cm; shorter edges never scanned
    double               m_sgp_sample_step{0.3};  // cm; edge-interior sampling step
    double               m_sgp_point_radius{0.2};  // cm; test_good_point radius
    double               m_good_point_pitch_frac{0};
    double               m_sgp_weak_scale{0};  // weak-charge penalty scale; 0 = off
    double               m_sgp_weak_qref{2000};  // charge ref (calc_charge_wcp units, no cm conversion)
    double               m_sgp_max_sep{-1};  // cm; < 0 = off (unbounded, legacy)
    bool                 m_break_seg_orient{false};
    bool                 m_daughter_count_proto_main_vertex{false};
    bool                 m_daughter_count_proto_examine_showers{false};
    bool                 m_shower_pdg_from_start_segment{false};
    bool                 m_shower_pdg_from_shower_type{false};
    bool                 m_shower_pdg_exact_muon_test{false};
    bool                 m_pi0_id_shared_allocator{false};
    bool                 m_shower_flag_pdg_electron{false};
    bool                 m_shower_less_id_tiebreak{false};
    bool                 m_shower_endpoint_exclude_start_vertex{false};
    bool                 m_shower_endpoint_skip_orphan_vtx{false};
    bool                 m_shower_walk_visited_parity{false};
    bool                 m_track_pid_persist_dqdx{false};
    bool                 m_shower_reclass_dqdx_guard{false};
    bool                 m_shower_topo_dqdx_guard{false};
    bool                 m_track_pid_persist_4mom{false};
    bool                 m_shower_proton_daughter_pion{false};
    bool                 m_shower_proton_daughter_pion_dissolve{false};
    bool                 m_muon_multi_proton_pion{false};
    bool                 m_track_pid_persist_dqdx_electron_guard{false};  // doc pr/40 round 5 F9
    bool                 m_shower_connect_main_vertex_straight_guard{false};  // doc pr/40 round 5 F10
    bool                 m_shower_traj_straight_guard{false};  // doc pr/40 round 5 F11
    bool                 m_shower_absorb_track_guard{false};  // doc pr/40 round 6 F12
    bool                 m_shower_absorb_unreachable_main{false};  // doc pr/65 round 3
    bool                 m_michel_stem_muon_rescue{false};  // doc pr/40 round 6 F14
    bool                 m_shower_in_cascade_guard{false};  // doc pr/74 round 2 P1
    double               m_shower_in_max_len{40};  // cm; pr/74 P1 tunable
    double               m_shower_in_mip_hi{1.3};  // ratio; pr/74 P1 tunable
    bool                 m_shower_connect_from_vertices_straight_guard{false};  // doc pr/40 round 9 (round 8 Part A)
    bool                 m_shower_connect_start_seg_straight_guard{false};  // doc pr/40 round 9 (round 7 c2c, D1 re-target)
    bool                 m_examine_direction_dirsign_shower_in_guard{false};  // doc pr/40 round 9 (round 7 c2a, D2 re-scope)
    bool                 m_daughter_shower_angle_reclass_straight_guard{false};  // doc pr/40 round 9 (round 7 c2b)
    bool                 m_shower_topo_reexam_straight_guard{false};  // doc pr/40 round 9 (round 7 c1 safety net)
    double               m_sfv_kink_max{25.0};  // degrees; continuation-arm tunable
    bool                 m_shower_nv_bridge_track{false};  // doc pr/40 round 9 B2
    bool                 m_shower_nv_main_pi_init{false};  // doc pr/97 D1; false = legacy indeterminate main_pi read
    double               m_shower_nv_bridge_max_gap{1.8};  // cm; B2 gap cut (steiner-cloud closest approach)
    bool                 m_kine_drop_stray_satellites{false};  // doc pr/92 master
    double               m_kine_sat_min_energy{20.0};  // MeV; drop-candidate floor
    double               m_kine_sat_prox_max{8.0};  // cm; main-cluster proximity exemption
    double               m_kine_sat_angle_bad{60.0};  // degrees; Arm A attachment-angle cut
    double               m_kine_sat_angle_main{45.0};  // degrees; Arm B main-vertex-angle cut
    double               m_kine_sat_far_dis{90.0};  // cm; Arm B far-attachment trigger
    double               m_kine_sat_axis_dis_cut{30.0};  // cm; shower axis integration radius
    double               m_kine_sat_cont_kink{25.0};  // degrees; Arm C continuation kink
    double               m_kine_sat_track_max_nseg{3.0};  // count; round-2 track-like max segments
    double               m_kine_sat_em_far_dis{150.0};  // cm; round-2 EM-satellite far-drop distance
    double               m_kine_sat_cont_keep_deg{0.0};  // doc pr/146; deg (ang_sv), 0 = off => arm C unchanged
    bool                 m_michel_stem_michel_check{false};  // doc pr/74 round 2 P2
    double               m_michel_stem_max_far_len{40};  // cm; pr/74 P2 tunable
    bool                 m_shower_stem_backfill{false};  // doc pr/74 round 2 K4
    double               m_stem_backfill_max_len{30};  // cm; pr/74 K4 tunable
    double               m_stem_backfill_mip_lo{0.75};  // ratio; pr/74 K4 tunable
    double               m_stem_backfill_mip_hi{3.5};  // ratio; pr/74 K4 tunable
    double               m_stem_backfill_min_shower_len{40};  // cm; pr/74 K4 tunable
    bool                 m_shower_conn3_unreachable{false};  // doc pr/74 round 2 K5 (pr/65 rung 2)
    double               m_conn3_unreachable_min_len{10};  // cm; pr/74 K5 tunable
    double               m_conn3_stitch_max{0};  // cm; doc pr/84 r2 F3; 0 = off = byte-identical
    bool                 m_shower_dedup_start_seg{false};  // doc pr/84 r3 S1; false = off = byte-identical
    bool                 m_shower_traj_michel_stem{false};  // doc pr/74 round 4 K6 (18255-506746 muon+Michel)
    double               m_michel_stem_traj_min_len{15};  // cm; pr/74 K6 tunable
    double               m_michel_stem_traj_max_len{45};  // cm; pr/74 K6 tunable
    double               m_michel_stem_traj_mip_lo{1.3};  // x MIP median; pr/74 K6 tunable
    double               m_michel_stem_traj_max_far_len{40};  // cm; pr/74 K6 tunable (own ceiling, not P2's)
    double               m_michel_stem_traj_min_kink_deg{40.0};  // deg; pr/74 K6 tunable
    bool                 m_shower_long_muon_keep_type{false};  // doc pr/44
    bool                 m_shower_bragg_protect_start_segment{false};  // doc pr/40 round 10
    bool                 m_shower_reclass_case_b_dqdx_guard{false};  // doc pr/93 Cause A (55595)
    bool                 m_shower_accept_pid_guard{false};  // doc pr/93 Cause B (348471, 69314)
    double               m_shower_pid_guard_min_len{50};  // cm; shared Cause A/B floor, inert while both off
    bool                 m_shower_vote_track_pid_counts{false};  // doc pr/93 Cause C (292643)
    bool                 m_shower_cone_absorb_guard{false};  // doc pr/93 Cause D (315167)
    bool                 m_shower_detach_track_stem{false};  // doc pr/93 r4 (348471, 292643)
    bool                 m_shower_ghost_member_drop{false};  // doc pr/99 r2 (395148 projective ghost)
    double               m_shower_ghost_overlap_frac{0.7};  // 2nd-best per-view overlap gate; inert while drop off
    double               m_shower_ghost_dqdx_ratio{0.25};  // starved gate, ratio vs mip median; inert while drop off
    double               m_shower_ghost_min_len{10.0};  // cm; scaled at copy; inert while drop off
    bool                 m_kine_charge_dedup{false};  // doc pr/99 r3 C1 (168596 Enu double count)
    bool                 m_kine_charge_rebuild{false};  // doc pr/99 r3 C1b (prototype cloud-rebuild parity)
    bool                 m_kine_charge_track_ctx{false};  // doc pr/101 K1 (37112 shower<->track overlap)
    bool                 m_kine_mass_rules{false};  // doc pr/101 K2 (proton shower +938 MeV)
    bool                 m_kine_hadronic_dqdx{false};  // doc pr/101 K3 (hadronic shower KE = sum dE/dx)
    int                  m_kine_long_muon_mode{0};  // doc pr/101 K4 (0 dQdx, 1 range, 2 range w/ fallback)
    double               m_kine_long_muon_ratio_lo{0.3};  // inert unless mode 2
    double               m_kine_long_muon_ratio_hi{0.5};  // inert unless mode 2
    bool                 m_kine_dqdx_skip_zero_dx{false};  // doc pdvd/45 sec 5.4: skip dx<=0 fit points in the vector cal_kine_dQdx (NaN Enu)
    bool                 m_long_muon_range_empty_chain_fallback{false};  // doc 84 round 1 (P1): range over muon-typed members when the chain missed the shower
    bool                 m_long_muon_members_geometry{false};  // doc 84 round 2: add out-of-chain muon members to range/endpoint (313847, 281595)
    bool                 m_kine_mainvtx_used_guard{false};  // doc pr/101 K5 (main-vertex member double count)
    bool                 m_shower_hadronic_tag{false};  // doc pr/99 r3 A5 (hadronic shower labeled e-)
    double               m_shower_hadronic_min_len{10.0};  // cm; scaled at copy; inert while tag off
    double               m_shower_hadronic_scan_len{30.0};  // cm; scaled at copy; inert while tag off
    double               m_shower_hadronic_bin{3.0};  // cm; scaled at copy; inert while tag off
    double               m_shower_hadronic_r_cyl{8.0};  // cm; scaled at copy; inert while tag off
    double               m_shower_hadronic_r_core{1.2};  // cm; scaled at copy; inert while tag off
    double               m_shower_hadronic_growth_max{0.8};  // ratio; inert while tag off
    double               m_shower_hadronic_growth_bragg{1.2};  // ratio; inert while tag off
    double               m_shower_hadronic_bragg_ratio{3.0};  // ratio; inert while tag off
    double               m_shower_hadronic_stem_ratio{0.0};  // MIP units; 0 = branch off; inert while tag off
    bool                 m_kine_count_orphan_tracks{false};  // doc pr/93 r4 (315167)
    double               m_kine_orphan_track_min{50};  // cm; scaled at copy
    bool                 m_shower_pass4_best_owner{false};  // doc pr/117 r1 (48% of wrongly-held charge)
    bool                 m_shower_merge_relax{false};  // doc pr/117 r1 (20-event merge class)
    double               m_shower_merge_relax_dis{6.0};  // cm; scaled at copy; inert while off
    double               m_shower_merge_relax_angle{15.0};  // deg, no conversion; inert while off
    double               m_shower_merge_relax_min_len{5.0};  // cm; scaled at copy; fragment length floor; inert while off
    bool                 m_shower_merge_relax_continuity{false};  // doc pr/118 r1 (two-tier axis+charge merge path)
    double               m_shower_merge_relax_cont_frac{1.0};  // fraction, no conversion; inert while off
    double               m_shower_merge_relax_cont_gap{8.0};  // cm; scaled at copy; inert while off
    double               m_shower_merge_relax_cont_qmed{5000.0};  // charge units, no conversion; inert while off
    double               m_shower_merge_relax_cont_axis{7.5};  // deg, no conversion; inert while off
    double               m_shower_merge_relax_cont_dmax{120.0};  // cm; scaled at copy; inert while off
    double               m_shower_merge_relax_cont_t1_gap{1.0};  // cm; scaled at copy; inert while off
    double               m_shower_merge_relax_cont_t1_fold{30.0};  // deg, no conversion; inert while off
    bool                 m_stem_backfill_back_guard{false};  // doc pr/120 r1 (47212/281567 backward stems)
    double               m_stem_backfill_back_ang{110.0};  // deg, no conversion; inert while off
    bool                 m_shower_ex1_dedup_rehome{false};  // doc pr/121 r1 (348471 dedup orphaning)
    bool                 m_shower_pass4_prune_detached{false};  // doc pr/123 r1; false = no prune pass
    double               m_shower_pass4_prune_gap{40.0};  // doc pr/123 r1; cm, component linkage gap
    double               m_shower_pass4_prune_gap2{0.0};  // doc pr/124 A; cm, 0 = no tier-2 band prune
    double               m_shower_pass4_prune2_ang{40.0};  // doc pr/124 A; deg off kept-core centroid
    double               m_shower_pass4_prune2_mdqdx{2.5};  // doc pr/124 A; x mip_dqdx_median
    double               m_shower_pass3_cone_guard_len{0.0};  // doc pr/124 C; cm, 0 = no pass3 track-pdg decline
    bool                 m_shower_samevtx_track_absorb{false};  // doc pr/125; false = no pass
    double               m_shower_samevtx_absorb_gap{6.0};  // doc pr/125; cm, frag<->host cloud gap cap
    double               m_shower_samevtx_absorb_max_len{50.0};  // doc pr/125; cm, fragment length cap
    double               m_shower_samevtx_absorb_min_len{5.0};  // doc pr/125; cm, fragment length floor
    bool                 m_shower_satellite_absorb{false};  // doc pr/125; false = no pass
    bool                 m_shower_split{false};  // doc pr/138 B2; false = no pass
    double               m_shower_split_max_valley{0.95};  // doc pr/138 B2; sec A5.4 knee
    double               m_shower_split_min_frac{0.03};  // doc pr/138 B2; per-seed charge share floor
    int                  m_shower_split_max_parts{2};  // doc pr/138 B3; 2 = the measured-exact kernel
    double               m_shower_split_min_charge{1e6};  // doc pr/138 B1; candidate charge floor (raw Fit::dQ)
    int                  m_shower_split_min_nseg{3};  // doc pr/138 B1; candidate member-count floor
    double               m_shower_split_bundle_gap{4};  // doc pr/138 B3; cm; single-linkage bundle gap
    double               m_shower_split_snap{0.80};  // doc pr/138 B3; k>=3 bundle dominance floor
    bool                 m_shower_split_skip_shared{false};  // doc pr/139 P1.1; refuse a component holding a segment another shower also owns
    bool                 m_shower_split_shed_shared{false};  // doc pr/139 sec 15; shed an ENTIRELY co-owned refused component instead of refusing it
    int                  m_shower_split_max_seeds{4};  // doc pr/139 sec 17; angular-maxima cap (shipped 4)
    double               m_shower_split_em_type_max_len{0};  // doc pr/139 sec 25; cm; 0 = off
    double               m_shower_split_max_impact{0};  // doc pr/139 P1.2; cm; 0 = no bound
    bool                 m_shower_split_em_start{false};  // doc pr/139 P1.3; seed the daughter on its nearest EM-typed member
    bool                 m_shower_split_rehome{false};  // doc pr/139 P1.4; offer an orphan daughter to the nearest larger EM shower
    double               m_shower_split_rehome_gap{4};  // doc pr/139 P1.4; cm; max daughter->host 3-D gap
    double               m_shower_satellite_absorb_max_mev{10.0};  // doc pr/125; MeV, satellite kine cap
    double               m_shower_satellite_absorb_host_mev{20.0};  // doc pr/125; MeV, host kine floor
    double               m_shower_pass4_track_guard_len{0.0};  // doc pr/123 r1; cm, 0 = no length guard
    double               m_shower_pass4_prox_guard_len{0.0};  // doc pr/130 item 1b; cm, 0 = pass4_proximity unguarded (legacy)
    double               m_shower_pass3_backfill_guard_len{0.0};  // doc pr/130 item 1b; cm, 0 = pass3 sibling backfill ignores pr/124's decline (legacy)
    double               m_stem_backfill_back_dvtx{0.0};  // doc pr/130 item B; cm, 0 = back guard ignores the vertex distance (legacy)
    bool                 m_shower_pass4_prefilter_v1_escape{false};  // doc pr/136 r2; false = pass-4 angle_v2>30 pre-filter discards regardless of angle_v1 (legacy)
    double               m_shower_pass4_prefilter_v1_max_v2{0.0};  // doc pr/136 r2; deg, 0 = the angle_v1 escape has no angle_v2 ceiling
    double               m_shower_pass4_prefilter_v1_max_dis{0.0};  // doc pr/136 r3; cm, 0 = the angle_v1 escape has no proximity bound
    double               m_pi0_mass_offset{10.0};  // doc pr/132 K1; MeV, the finders' "+10" window offset
    double               m_pi0_assoc_angle_deg{30.0};  // doc pr/132 K2; deg, disconnected<->vertex association cut
    double               m_pi0_attached_partner_min_mev{0.0};  // doc pr/132 K3; MeV, 0 = no nueCC-fake guard
    int                  m_pi0_nv_max_prongs{2};  // doc pr/132 K5; without-vertex GATE1 prong cap
    bool                 m_pi0_readmit_retyped{false};  // doc pr/132 K7; readmit hadronic-retyped showers into pi0 pairing
    bool                 m_pi0_admit_type3{false};  // doc pr/132 K8; admit conn_type-3 showers into the with-vertex pool
    double               m_pi0_crumb_assoc_mev{0.0};  // doc pr/132 K9; MeV, 0 = crumbs keep the association-angle test
    double               m_pi0_collinear_merge_deg{0.0};  // doc pr/132 K12; deg, 0 = pairing sees each detached fragment alone
    double               m_pi0_nv_partner_min_mev{0.0};  // doc pr/132 K13; MeV, 0 = no path-2 partner floor
    double               m_shower_em_collinear_deg{0.0};  // doc pr/132 K16; deg, 0 = no build-time EM collinear merge
    double               m_shower_em_collinear_dis_cm{60.0};  // doc pr/132 K16; cm, host->fragment reach
    double               m_shower_em_collinear_host_mev{20.0};  // doc pr/132 K16; MeV, host floor
    double               m_shower_em_backext_perp_cm{0.0};  // doc pr/132 K17; cm, 0 = no start back-extension
    double               m_shower_em_backext_len_cm{40.0};  // doc pr/132 K17; cm, upstream reach
    double               m_pi0_accept_merge_dis_cm{0.0};  // doc pr/132 K18; cm, 0 = no acceptance-aware fragment merge
    double               m_pi0_bp_vertex_miss_cm{0.0};  // doc pr/132 K19; cm, 0 = no back-projection NC vertex proposer
    bool                 m_pi0_admit_muon_showers{false};  // doc pr/133 K20; admit shower-topology mu-typed objects into the pi0 pools
    double               m_pi0_nc_sig_angle_deg{0.0};  // doc pr/133 K21; deg, 0 = off; owner NC signature for the bp proposer
    double               m_pi0_nc_floor_mev{0.0};  // doc pr/133 K21 v2; MeV, 0 = legacy 20; signature-mode partner floor
    double               m_pi0_nc_pf_assoc_deg{0.0};  // doc pr/133 K21 v2.2; deg, 0 = off; post-fire PF association cone
    bool                 m_pi0_nc_frag_merge{false};  // doc pr/134 K22; NC bp pairing at merged-complex level
    double               m_pi0_pf_assoc_deg{0.0};  // doc pr/134 K23; deg, 0 = off; with-vertex post-accept PF satellite cone
    bool                 m_pi0_prefer_main_vertex{false};  // doc pr/134 K24; P1 nu-vertex preference
    double               m_pi0_nv_max_vtx_shift_cm{0.0};  // doc pr/132 K10; cm, 0 = no without-vertex decay-point shift cap
    double               m_pi0_nv_mass_window_mev{60.0};  // doc pr/132 K11; MeV, without-vertex acceptance half-window (legacy 60)
    bool                 m_kine_count_guard_freed{false};  // doc pr/123 r2; kine twin of pf_orphan_guard_freed
    double               m_kine_guard_freed_impact{0.0};  // doc pr/129; 0 = no pointing test
    double               m_kine_guard_freed_miss_deg{90.0};  // doc pr/129
    bool                 m_kine_count_near_cross_cluster{false};  // doc pr/128 (137238 class generalised)
    double               m_kine_near_gap{5};  // cm; scaled at copy
    double               m_kine_near_min_len{30};  // cm; scaled at copy
    double               m_kine_near_end_tol{10};  // cm; scaled at copy
    double               m_kine_near_kink_deg{30};  // deg
    double               m_kine_near_pointing_impact{0.0};  // doc pr/144; cm; 0 = no pointing test
    double               m_kine_near_pointing_miss_deg{90.0};  // doc pr/144; deg
    bool                 m_kine_count_conn4_near{false};  // doc pr/128 (105074 / 179048)
    double               m_kine_conn4_near_gap{20};  // cm; scaled at copy
    bool                 m_straight_cont_cross_cluster{false};  // doc pr/93 r4 (137238)
    bool                 m_sccc_bridge_body{false};  // doc pr/93 r4 second rung
    double               m_sccc_max_gap{5};  // cm; base tier
    double               m_sccc_kink_max{15.0};  // deg; base tier
    double               m_sccc_gap_aligned{12};  // cm; aligned tier
    double               m_sccc_kink_tight{7.5};  // deg; aligned tier
    bool                 m_single_muon_proton_chain_veto{false};  // doc pr/43 round 2 K1
    bool                 m_single_muon_long_muon_claim{false};  // doc pr/43 round 2 K2
    bool                 m_pid_flag_reconcile{false};  // doc pr/43 round 2 K3
    bool                 m_long_muon_stub_bridge{false};  // doc pr/46
    double               m_long_muon_stub_bridge_len{6.0};  // doc 84 round 1 (P3): stub precondition [cm]; 6.0 = legacy
    bool                 m_long_muon_angle_relax_long{false};  // doc 84 round 1 (P2): >50 cm MIP continuation angle relax
    double               m_long_muon_angle_relax_deg{16.0};  // doc 84 round 1 (P2): relaxed cap [deg]; inert unless relax on
    std::vector<double> m_kine_plane_weights{0.25, 0.25, 1.0};  // {U,V,W}

    // TaggerCheckNeutrino::configure(), the kept keys verbatim.
    void read_neutrino_knobs(const Configuration& config) {
        m_dir_weak_use_score = get(config, "dir_weak_use_score", m_dir_weak_use_score);
        m_proton_dir_vote = get(config, "proton_dir_vote", m_proton_dir_vote);
        m_proton_dir_score_max = get(config, "proton_dir_score_max", m_proton_dir_score_max);
        m_proton_dir_asym_min = get(config, "proton_dir_asym_min", m_proton_dir_asym_min);
        m_endpoint_trim_retry = get(config, "endpoint_trim_retry", m_endpoint_trim_retry);
        m_fit_vertex_min_seg_length = get(config, "fit_vertex_min_seg_length", m_fit_vertex_min_seg_length);  // cm
        m_mvfit_robust = get(config, "mvfit_robust", m_mvfit_robust);
        m_mvfit_main_only = get(config, "mvfit_main_only", m_mvfit_main_only);
        m_mvfit_min_len = get(config, "mvfit_min_len", m_mvfit_min_len);  // cm
        m_mvfit_rin_margin = get(config, "mvfit_rin_margin", m_mvfit_rin_margin);  // cm
        m_mvfit_rout_frac = get(config, "mvfit_rout_frac", m_mvfit_rout_frac);
        m_mvfit_rout_min = get(config, "mvfit_rout_min", m_mvfit_rout_min);  // cm
        m_mvfit_rout_max = get(config, "mvfit_rout_max", m_mvfit_rout_max);  // cm
        m_mvfit_angle = get(config, "mvfit_angle", m_mvfit_angle);  // deg
        m_mvfit_min_pts = get(config, "mvfit_min_pts", m_mvfit_min_pts);
        m_mvfit_min_aniso = get(config, "mvfit_min_aniso", m_mvfit_min_aniso);
        m_mvfit_prior_range = get(config, "mvfit_prior_range", m_mvfit_prior_range);  // cm
        m_cathode_x = get(config, "cathode_x", m_cathode_x);  // cm
        m_cathode_kink_xcut = get(config, "cathode_kink_xcut", m_cathode_kink_xcut);  // cm
        m_cathode_wide_kink_angle = get(config, "cathode_wide_kink_angle", m_cathode_wide_kink_angle);  // deg
        m_cathode_wide_kink_skirt = get(config, "cathode_wide_kink_skirt", m_cathode_wide_kink_skirt);  // cm
        m_cathode_wide_kink_baseline = get(config, "cathode_wide_kink_baseline", m_cathode_wide_kink_baseline);  // cm
        m_two_end_break = get(config, "two_end_break", m_two_end_break);
        m_teb_min_len = get(config, "teb_min_len", m_teb_min_len);  // cm
        m_teb_min_arm = get(config, "teb_min_arm", m_teb_min_arm);  // cm
        m_teb_min_arm_pts = get(config, "teb_min_arm_pts", m_teb_min_arm_pts);
        m_teb_stub_max = get(config, "teb_stub_max", m_teb_stub_max);  // cm
        m_teb_accept_range = get(config, "teb_accept_range", m_teb_accept_range);  // cm
        m_teb_rise_r1 = get(config, "teb_rise_r1", m_teb_rise_r1);
        m_teb_rise_r2 = get(config, "teb_rise_r2", m_teb_rise_r2);
        m_teb_abs_end_min = get(config, "teb_abs_end_min", m_teb_abs_end_min);
        m_teb_dip_floor = get(config, "teb_dip_floor", m_teb_dip_floor);
        m_teb_score_cap_r1 = get(config, "teb_score_cap_r1", m_teb_score_cap_r1);
        m_teb_score_cap_r2 = get(config, "teb_score_cap_r2", m_teb_score_cap_r2);
        m_teb_turn_angle = get(config, "teb_turn_angle", m_teb_turn_angle);  // deg
        m_teb_turn_baseline = get(config, "teb_turn_baseline", m_teb_turn_baseline);  // cm
        m_teb_turn_skirt = get(config, "teb_turn_skirt", m_teb_turn_skirt);  // cm
        m_teb_turn_min_arm_frac = get(config, "teb_turn_min_arm_frac", m_teb_turn_min_arm_frac);  // frac of baseline; doc pr/90 round 2
        m_teb_bragg_veto_turn = get(config, "teb_bragg_veto_turn", m_teb_bragg_veto_turn);  // deg; doc pr/90 round 4 (D4)
        m_kink_walk_dqdx_stop = get(config, "kink_walk_dqdx_stop", m_kink_walk_dqdx_stop);
        m_kink_break_protect = get(config, "kink_break_protect", m_kink_break_protect);
        m_kink_dqdx_hot_ratio = get(config, "kink_dqdx_hot_ratio", m_kink_dqdx_hot_ratio);
        m_esva_ignore_empty_2d = get(config, "esva_ignore_empty_2d", m_esva_ignore_empty_2d);  // docs/73 sec 12 round 3
        m_main_vertex_graph_audit = get(config, "main_vertex_graph_audit", m_main_vertex_graph_audit);
        m_mvga_radius = get(config, "mvga_radius", m_mvga_radius);  // cm
        m_mvga_dup_tol = get(config, "mvga_dup_tol", m_mvga_dup_tol);  // cm
        m_mvga_dup_frac = get(config, "mvga_dup_frac", m_mvga_dup_frac);
        m_mvga_dup_angle = get(config, "mvga_dup_angle", m_mvga_dup_angle);  // deg
        m_mvga_bridge_mip = get(config, "mvga_bridge_mip", m_mvga_bridge_mip);
        m_mvga_reconnect = get(config, "mvga_reconnect", m_mvga_reconnect);  // cm
        m_mvga_stub = get(config, "mvga_stub", m_mvga_stub);  // cm
        m_mvga_stub_pts = get(config, "mvga_stub_pts", m_mvga_stub_pts);
        m_mvga_reseat_angle = get(config, "mvga_reseat_angle", m_mvga_reseat_angle);  // deg
        m_mvga_satellite = get(config, "mvga_satellite", m_mvga_satellite);  // cm
        m_mvga_interposed = get(config, "mvga_interposed", m_mvga_interposed);  // doc pr/85
        m_mvga_interposed_angle = get(config, "mvga_interposed_angle", m_mvga_interposed_angle);  // deg
        m_mvga_interposed_len = get(config, "mvga_interposed_len", m_mvga_interposed_len);  // cm; doc pr/86
        m_mvga_sat_dup_frac = get(config, "mvga_sat_dup_frac", m_mvga_sat_dup_frac);  // fraction; doc pr/86
        m_mvga_interposed_deg1 = get(config, "mvga_interposed_deg1", m_mvga_interposed_deg1);  // doc pr/86
        m_mvga_splice_straighten = get(config, "mvga_splice_straighten", m_mvga_splice_straighten);  // cm; doc pr/86 round 2
        m_mvga_approach_collapse = get(config, "mvga_approach_collapse", m_mvga_approach_collapse);  // cm; doc pr/86 round 2
        m_mvga_straighten_radius = get(config, "mvga_straighten_radius", m_mvga_straighten_radius);  // cm; doc pr/86 round 2
        m_mvga_op1_radius = get(config, "mvga_op1_radius", m_mvga_op1_radius);  // cm; doc pr/83 r3; 0 = use mvga_radius, -1 = unscoped
        m_mvga_op1_dup_frac = get(config, "mvga_op1_dup_frac", m_mvga_op1_dup_frac);  // doc pr/83 r3; 0 = use mvga_dup_frac
        m_mvga_op1_post = get(config, "mvga_op1_post", m_mvga_op1_post);  // doc pr/83 r3 (class A)
        m_swap_orphan_dup_audit = get(config, "swap_orphan_dup_audit", m_swap_orphan_dup_audit);  // doc pr/83 r3 (Mechanism C)
        m_mvga_proj_dup_frac = get(config, "mvga_proj_dup_frac", m_mvga_proj_dup_frac);  // doc pr/83 r4; 0 = pass disabled
        m_mvga_proj_dqdx_ratio = get(config, "mvga_proj_dqdx_ratio", m_mvga_proj_dqdx_ratio);  // doc pr/83 r4; inert while frac == 0
        m_mvga_proj_angle = get(config, "mvga_proj_angle", m_mvga_proj_angle);  // deg; doc pr/83 r4b; 0 = use mvga_dup_angle
        m_mvga_ac_chord_max = get(config, "mvga_ac_chord_max", m_mvga_ac_chord_max);  // cm; doc pr/99 round 2; 0 = no cap
        m_mvga_ac_no_cascade = get(config, "mvga_ac_no_cascade", m_mvga_ac_no_cascade);  // doc pr/99 round 2
        m_mvga_passthru = get(config, "mvga_passthru", m_mvga_passthru);  // cm; doc pr/103; 0 = off
        m_mvga_passthru_tol = get(config, "mvga_passthru_tol", m_mvga_passthru_tol);  // cm; doc pr/103
        m_mvga_interposed_fallback = get(config, "mvga_interposed_fallback", m_mvga_interposed_fallback);  // doc pr/103
        m_mvga_interposed_fallback_min_angle = get(config, "mvga_interposed_fallback_min_angle", m_mvga_interposed_fallback_min_angle);  // deg; doc pr/103
        m_mvga_dup_starved_asym = get(config, "mvga_dup_starved_asym", m_mvga_dup_starved_asym);  // pair asymmetry; doc pr/99 round 2; 0 = off
        m_mvga_dup_starved_mip = get(config, "mvga_dup_starved_mip", m_mvga_dup_starved_mip);  // absolute cap on loser; doc pr/99 round 2; 0 = off
        m_mvga_dup_starved_span = get(config, "mvga_dup_starved_span", m_mvga_dup_starved_span);  // pair length comparability; doc pr/99 round 2; 0 = off
        m_shower_topo_demote_len = get(config, "shower_topo_demote_len", m_shower_topo_demote_len);  // cm
        m_fit_blob_coverage_defer = get(config, "fit_blob_coverage_defer", m_fit_blob_coverage_defer);
        m_fit_exclusion = get(config, "fit_exclusion", m_fit_exclusion);
        m_graph_endpoint_tol = get(config, "graph_endpoint_tol", m_graph_endpoint_tol);  // cm
        m_oov_prototype_parity = get(config, "oov_prototype_parity", m_oov_prototype_parity);
        m_first_seg_local_pca = get(config, "first_seg_local_pca", m_first_seg_local_pca);
        m_other_seg_relaxed_accept = get(config, "other_seg_relaxed_accept", m_other_seg_relaxed_accept);
        m_other_seg_empty_2d_guard = get(config, "other_seg_empty_2d_guard", m_other_seg_empty_2d_guard);
        m_other_seg_keep_isolated = get(config, "other_seg_keep_isolated", m_other_seg_keep_isolated);
        m_other_seg_keep_isolated_min_points = get(config, "other_seg_keep_isolated_min_points", m_other_seg_keep_isolated_min_points);
        m_other_seg_keep_isolated_min_length = get(config, "other_seg_keep_isolated_min_length", m_other_seg_keep_isolated_min_length);  // cm
        m_other_seg_keep_isolated_len_admit = get(config, "other_seg_keep_isolated_len_admit", m_other_seg_keep_isolated_len_admit);  // cm
        m_iso_snap_min_dir_mag = get(config, "iso_snap_min_dir_mag", m_iso_snap_min_dir_mag);  // cm
        m_assoc_full_recluster = get(config, "assoc_full_recluster", m_assoc_full_recluster);
        m_assoc_reassign_orphans = get(config, "assoc_reassign_orphans", m_assoc_reassign_orphans);
        m_assoc_clear_on_merge = get(config, "assoc_clear_on_merge", m_assoc_clear_on_merge);
        m_shower_topo_proto_dir = get(config, "shower_topo_proto_dir", m_shower_topo_proto_dir);
        m_vertex_dir_use_fit_point = get(config, "vertex_dir_use_fit_point", m_vertex_dir_use_fit_point);
        m_shower_traj_recheck_parity = get(config, "shower_traj_recheck_parity", m_shower_traj_recheck_parity);
        m_main_vertex_require_descriptor = get(config, "main_vertex_require_descriptor", m_main_vertex_require_descriptor);
        m_main_vertex_candidate_flag = get(config, "main_vertex_candidate_flag", m_main_vertex_candidate_flag);
        m_cont_muon_dir3_30cm = get(config, "cont_muon_dir3_30cm", m_cont_muon_dir3_30cm);
        m_track_comp_empty_abstain = get(config, "track_comp_empty_abstain", m_track_comp_empty_abstain);
        m_shower_topo_reset = get(config, "shower_topo_reset", m_shower_topo_reset);
        m_reclass_preserve_4mom = get(config, "reclass_preserve_4mom", m_reclass_preserve_4mom);
        m_reclass_never_computed_ke_floor = get(config, "reclass_never_computed_ke_floor", m_reclass_never_computed_ke_floor);
        m_dir_track_median_local = get(config, "dir_track_median_local", m_dir_track_median_local);
        m_examine_showers_vertex_by_index = get(config, "examine_showers_vertex_by_index", m_examine_showers_vertex_by_index);
        m_iso_endpoint = get(config, "iso_endpoint", m_iso_endpoint);
        m_iso_endpoint_min_length = get(config, "iso_endpoint_min_length", m_iso_endpoint_min_length);  // cm
        m_iso_endpoint_max_xext = get(config, "iso_endpoint_max_xext", m_iso_endpoint_max_xext);  // cm
        m_iso_endpoint_xext_frac = get(config, "iso_endpoint_xext_frac", m_iso_endpoint_xext_frac);
        m_iso_endpoint_xext_quantile = get(config, "iso_endpoint_xext_quantile", m_iso_endpoint_xext_quantile);
        m_iso_endpoint_tube_radius = get(config, "iso_endpoint_tube_radius", m_iso_endpoint_tube_radius);  // cm
        m_iso_endpoint_min_aspect = get(config, "iso_endpoint_min_aspect", m_iso_endpoint_min_aspect);
        m_traj_cover_probe = get(config, "traj_cover_probe", m_traj_cover_probe);
        m_pr_find_other_rounds = get(config, "pr_find_other_rounds", m_pr_find_other_rounds);
        m_v3_extension_guard = get(config, "v3_extension_guard", m_v3_extension_guard);
        m_v3_extension_min_gain = get(config, "v3_extension_min_gain", m_v3_extension_min_gain);  // cm
        m_es3_stub_guard = get(config, "es3_stub_guard", m_es3_stub_guard);
        m_es3sg_stub_max = get(config, "es3sg_stub_max", m_es3sg_stub_max);  // cm
        m_es3sg_len_ratio = get(config, "es3sg_len_ratio", m_es3sg_len_ratio);
        m_es3sg_ang3_min = get(config, "es3sg_ang3_min", m_es3sg_ang3_min);  // deg
        m_es3sg_ang_ratio = get(config, "es3sg_ang_ratio", m_es3sg_ang_ratio);
        m_es3sg_require_terminal = get(config, "es3sg_require_terminal", m_es3sg_require_terminal);
        m_vertex_z_prior_scale = get(config, "vertex_z_prior_scale", m_vertex_z_prior_scale);
        m_kine_fudge_factor = get(config, "kine_fudge_factor", m_kine_fudge_factor);
        m_kine_recom_factor = get(config, "kine_recom_factor", m_kine_recom_factor);
        m_kine_shower_fudge_factor = get(config, "kine_shower_fudge_factor", m_kine_shower_fudge_factor);
        m_kine_shower_recom_factor = get(config, "kine_shower_recom_factor", m_kine_shower_recom_factor);
        m_kine_proton_recom_factor = get(config, "kine_proton_recom_factor", m_kine_proton_recom_factor);
        m_kine_plane_asym_switch = get(config, "kine_plane_asym_switch", m_kine_plane_asym_switch);
        m_kine_w_value = get(config, "kine_w_value", m_kine_w_value);
        m_kine_shower_pdg_live = get(config, "kine_shower_pdg_live", m_kine_shower_pdg_live);
        m_steiner_gap_penalty = get(config, "steiner_gap_penalty", m_steiner_gap_penalty);
        m_sgp_dead_alpha = get(config, "sgp_dead_alpha", m_sgp_dead_alpha);
        m_sgp_min_edge = get(config, "sgp_min_edge", m_sgp_min_edge);  // cm
        m_sgp_sample_step = get(config, "sgp_sample_step", m_sgp_sample_step);  // cm
        m_sgp_point_radius = get(config, "sgp_point_radius", m_sgp_point_radius);  // cm
        m_good_point_pitch_frac = get(config, "good_point_pitch_frac", m_good_point_pitch_frac);
        m_sgp_weak_scale = get(config, "sgp_weak_scale", m_sgp_weak_scale);
        m_sgp_weak_qref = get(config, "sgp_weak_qref", m_sgp_weak_qref);  // charge units
        m_sgp_max_sep = get(config, "sgp_max_sep", m_sgp_max_sep);  // cm
        m_break_seg_orient = get(config, "break_seg_orient", m_break_seg_orient);
        m_daughter_count_proto_main_vertex = get(config, "daughter_count_proto_main_vertex", m_daughter_count_proto_main_vertex);
        m_daughter_count_proto_examine_showers = get(config, "daughter_count_proto_examine_showers", m_daughter_count_proto_examine_showers);
        m_shower_pdg_from_start_segment = get(config, "shower_pdg_from_start_segment", m_shower_pdg_from_start_segment);
        m_shower_pdg_from_shower_type = get(config, "shower_pdg_from_shower_type", m_shower_pdg_from_shower_type);
        m_shower_pdg_exact_muon_test = get(config, "shower_pdg_exact_muon_test", m_shower_pdg_exact_muon_test);
        m_pi0_id_shared_allocator = get(config, "pi0_id_shared_allocator", m_pi0_id_shared_allocator);
        m_shower_flag_pdg_electron = get(config, "shower_flag_pdg_electron", m_shower_flag_pdg_electron);
        m_shower_less_id_tiebreak = get(config, "shower_less_id_tiebreak", m_shower_less_id_tiebreak);
        m_shower_endpoint_exclude_start_vertex = get(config, "shower_endpoint_exclude_start_vertex", m_shower_endpoint_exclude_start_vertex);
        m_shower_endpoint_skip_orphan_vtx = get(config, "shower_endpoint_skip_orphan_vtx", m_shower_endpoint_skip_orphan_vtx);
        m_shower_walk_visited_parity = get(config, "shower_walk_visited_parity", m_shower_walk_visited_parity);
        m_track_pid_persist_dqdx = get(config, "track_pid_persist_dqdx", m_track_pid_persist_dqdx);
        m_shower_reclass_dqdx_guard = get(config, "shower_reclass_dqdx_guard", m_shower_reclass_dqdx_guard);
        m_shower_topo_dqdx_guard = get(config, "shower_topo_dqdx_guard", m_shower_topo_dqdx_guard);
        m_track_pid_persist_4mom = get(config, "track_pid_persist_4mom", m_track_pid_persist_4mom);
        m_shower_proton_daughter_pion = get(config, "shower_proton_daughter_pion", m_shower_proton_daughter_pion);
        m_shower_proton_daughter_pion_dissolve = get(config, "shower_proton_daughter_pion_dissolve", m_shower_proton_daughter_pion_dissolve);
        m_muon_multi_proton_pion = get(config, "muon_multi_proton_pion", m_muon_multi_proton_pion);
        m_track_pid_persist_dqdx_electron_guard = get(config, "track_pid_persist_dqdx_electron_guard", m_track_pid_persist_dqdx_electron_guard);
        m_shower_connect_main_vertex_straight_guard = get(config, "shower_connect_main_vertex_straight_guard", m_shower_connect_main_vertex_straight_guard);
        m_shower_traj_straight_guard = get(config, "shower_traj_straight_guard", m_shower_traj_straight_guard);
        m_shower_absorb_track_guard = get(config, "shower_absorb_track_guard", m_shower_absorb_track_guard);
        m_shower_absorb_unreachable_main = get(config, "shower_absorb_unreachable_main", m_shower_absorb_unreachable_main);
        m_michel_stem_muon_rescue = get(config, "michel_stem_muon_rescue", m_michel_stem_muon_rescue);
        m_shower_in_cascade_guard = get(config, "shower_in_cascade_guard", m_shower_in_cascade_guard);
        m_shower_in_max_len = get(config, "shower_in_max_len", m_shower_in_max_len);
        m_shower_in_mip_hi = get(config, "shower_in_mip_hi", m_shower_in_mip_hi);
        m_shower_connect_from_vertices_straight_guard = get(config, "shower_connect_from_vertices_straight_guard", m_shower_connect_from_vertices_straight_guard);
        m_shower_connect_start_seg_straight_guard = get(config, "shower_connect_start_seg_straight_guard", m_shower_connect_start_seg_straight_guard);
        m_examine_direction_dirsign_shower_in_guard = get(config, "examine_direction_dirsign_shower_in_guard", m_examine_direction_dirsign_shower_in_guard);
        m_daughter_shower_angle_reclass_straight_guard = get(config, "daughter_shower_angle_reclass_straight_guard", m_daughter_shower_angle_reclass_straight_guard);
        m_shower_topo_reexam_straight_guard = get(config, "shower_topo_reexam_straight_guard", m_shower_topo_reexam_straight_guard);
        m_sfv_kink_max = get(config, "sfv_kink_max", m_sfv_kink_max);
        m_shower_nv_bridge_track = get(config, "shower_nv_bridge_track", m_shower_nv_bridge_track);
        m_shower_nv_bridge_max_gap = get(config, "shower_nv_bridge_max_gap", m_shower_nv_bridge_max_gap);
        m_shower_nv_main_pi_init = get(config, "shower_nv_main_pi_init", m_shower_nv_main_pi_init);
        m_kine_drop_stray_satellites = get(config, "kine_drop_stray_satellites", m_kine_drop_stray_satellites);
        m_kine_sat_min_energy = get(config, "kine_sat_min_energy", m_kine_sat_min_energy);
        m_kine_sat_prox_max = get(config, "kine_sat_prox_max", m_kine_sat_prox_max);
        m_kine_sat_angle_bad = get(config, "kine_sat_angle_bad", m_kine_sat_angle_bad);
        m_kine_sat_angle_main = get(config, "kine_sat_angle_main", m_kine_sat_angle_main);
        m_kine_sat_far_dis = get(config, "kine_sat_far_dis", m_kine_sat_far_dis);
        m_kine_sat_axis_dis_cut = get(config, "kine_sat_axis_dis_cut", m_kine_sat_axis_dis_cut);
        m_kine_sat_cont_kink = get(config, "kine_sat_cont_kink", m_kine_sat_cont_kink);
        m_kine_sat_track_max_nseg = get(config, "kine_sat_track_max_nseg", m_kine_sat_track_max_nseg);
        m_kine_sat_em_far_dis = get(config, "kine_sat_em_far_dis", m_kine_sat_em_far_dis);
        m_kine_sat_cont_keep_deg = get(config, "kine_sat_cont_keep_deg", m_kine_sat_cont_keep_deg);
        m_michel_stem_michel_check = get(config, "michel_stem_michel_check", m_michel_stem_michel_check);
        m_michel_stem_max_far_len = get(config, "michel_stem_max_far_len", m_michel_stem_max_far_len);
        m_shower_stem_backfill = get(config, "shower_stem_backfill", m_shower_stem_backfill);
        m_stem_backfill_max_len = get(config, "stem_backfill_max_len", m_stem_backfill_max_len);
        m_stem_backfill_mip_lo = get(config, "stem_backfill_mip_lo", m_stem_backfill_mip_lo);
        m_stem_backfill_mip_hi = get(config, "stem_backfill_mip_hi", m_stem_backfill_mip_hi);
        m_stem_backfill_min_shower_len = get(config, "stem_backfill_min_shower_len", m_stem_backfill_min_shower_len);
        m_shower_conn3_unreachable = get(config, "shower_conn3_unreachable", m_shower_conn3_unreachable);
        m_conn3_unreachable_min_len = get(config, "conn3_unreachable_min_len", m_conn3_unreachable_min_len);
        m_conn3_stitch_max = get(config, "conn3_stitch_max", m_conn3_stitch_max);
        m_shower_dedup_start_seg = get(config, "shower_dedup_start_seg", m_shower_dedup_start_seg);
        m_shower_traj_michel_stem = get(config, "shower_traj_michel_stem", m_shower_traj_michel_stem);
        m_michel_stem_traj_min_len = get(config, "michel_stem_traj_min_len", m_michel_stem_traj_min_len);
        m_michel_stem_traj_max_len = get(config, "michel_stem_traj_max_len", m_michel_stem_traj_max_len);
        m_michel_stem_traj_mip_lo = get(config, "michel_stem_traj_mip_lo", m_michel_stem_traj_mip_lo);
        m_michel_stem_traj_max_far_len = get(config, "michel_stem_traj_max_far_len", m_michel_stem_traj_max_far_len);
        m_michel_stem_traj_min_kink_deg = get(config, "michel_stem_traj_min_kink_deg", m_michel_stem_traj_min_kink_deg);
        m_shower_long_muon_keep_type = get(config, "shower_long_muon_keep_type", m_shower_long_muon_keep_type);
        m_shower_bragg_protect_start_segment = get(config, "shower_bragg_protect_start_segment", m_shower_bragg_protect_start_segment);
        m_shower_reclass_case_b_dqdx_guard = get(config, "shower_reclass_case_b_dqdx_guard", m_shower_reclass_case_b_dqdx_guard);
        m_shower_accept_pid_guard = get(config, "shower_accept_pid_guard", m_shower_accept_pid_guard);
        m_shower_pid_guard_min_len = get(config, "shower_pid_guard_min_len", m_shower_pid_guard_min_len);
        m_shower_vote_track_pid_counts = get(config, "shower_vote_track_pid_counts", m_shower_vote_track_pid_counts);
        m_shower_cone_absorb_guard = get(config, "shower_cone_absorb_guard", m_shower_cone_absorb_guard);
        m_shower_detach_track_stem = get(config, "shower_detach_track_stem", m_shower_detach_track_stem);
        m_shower_ghost_member_drop = get(config, "shower_ghost_member_drop", m_shower_ghost_member_drop);
        m_shower_ghost_overlap_frac = get(config, "shower_ghost_overlap_frac", m_shower_ghost_overlap_frac);
        m_shower_ghost_dqdx_ratio = get(config, "shower_ghost_dqdx_ratio", m_shower_ghost_dqdx_ratio);
        m_shower_ghost_min_len = get(config, "shower_ghost_min_len", m_shower_ghost_min_len);
        m_kine_charge_dedup = get(config, "kine_charge_dedup", m_kine_charge_dedup);
        m_kine_charge_rebuild = get(config, "kine_charge_rebuild", m_kine_charge_rebuild);
        m_kine_charge_track_ctx = get(config, "kine_charge_track_ctx", m_kine_charge_track_ctx);
        m_kine_mass_rules = get(config, "kine_mass_rules", m_kine_mass_rules);
        m_kine_hadronic_dqdx = get(config, "kine_hadronic_dqdx", m_kine_hadronic_dqdx);
        m_kine_long_muon_mode = get(config, "kine_long_muon_mode", m_kine_long_muon_mode);
        m_kine_long_muon_ratio_lo = get(config, "kine_long_muon_ratio_lo", m_kine_long_muon_ratio_lo);
        m_kine_long_muon_ratio_hi = get(config, "kine_long_muon_ratio_hi", m_kine_long_muon_ratio_hi);
        m_kine_dqdx_skip_zero_dx = get(config, "kine_dqdx_skip_zero_dx", m_kine_dqdx_skip_zero_dx);  // doc pdvd/45 sec 5.4
        m_long_muon_range_empty_chain_fallback = get(config, "long_muon_range_empty_chain_fallback", m_long_muon_range_empty_chain_fallback);  // doc 84 round 1 (P1)
        m_long_muon_members_geometry = get(config, "long_muon_members_geometry", m_long_muon_members_geometry);  // doc 84 round 2
        m_kine_mainvtx_used_guard = get(config, "kine_mainvtx_used_guard", m_kine_mainvtx_used_guard);
        m_shower_hadronic_tag = get(config, "shower_hadronic_tag", m_shower_hadronic_tag);
        m_shower_hadronic_min_len = get(config, "shower_hadronic_min_len", m_shower_hadronic_min_len);
        m_shower_hadronic_scan_len = get(config, "shower_hadronic_scan_len", m_shower_hadronic_scan_len);
        m_shower_hadronic_bin = get(config, "shower_hadronic_bin", m_shower_hadronic_bin);
        m_shower_hadronic_r_cyl = get(config, "shower_hadronic_r_cyl", m_shower_hadronic_r_cyl);
        m_shower_hadronic_r_core = get(config, "shower_hadronic_r_core", m_shower_hadronic_r_core);
        m_shower_hadronic_growth_max = get(config, "shower_hadronic_growth_max", m_shower_hadronic_growth_max);
        m_shower_hadronic_growth_bragg = get(config, "shower_hadronic_growth_bragg", m_shower_hadronic_growth_bragg);
        m_shower_hadronic_bragg_ratio = get(config, "shower_hadronic_bragg_ratio", m_shower_hadronic_bragg_ratio);
        m_shower_hadronic_stem_ratio = get(config, "shower_hadronic_stem_ratio", m_shower_hadronic_stem_ratio);
        m_kine_count_orphan_tracks = get(config, "kine_count_orphan_tracks", m_kine_count_orphan_tracks);
        m_kine_orphan_track_min = get(config, "kine_orphan_track_min", m_kine_orphan_track_min);
        m_shower_pass4_best_owner = get(config, "shower_pass4_best_owner", m_shower_pass4_best_owner);
        m_shower_merge_relax = get(config, "shower_merge_relax", m_shower_merge_relax);
        m_shower_merge_relax_dis = get(config, "shower_merge_relax_dis", m_shower_merge_relax_dis);
        m_shower_merge_relax_angle = get(config, "shower_merge_relax_angle", m_shower_merge_relax_angle);
        m_shower_merge_relax_min_len = get(config, "shower_merge_relax_min_len", m_shower_merge_relax_min_len);
        m_shower_merge_relax_continuity = get(config, "shower_merge_relax_continuity", m_shower_merge_relax_continuity);
        m_shower_merge_relax_cont_frac = get(config, "shower_merge_relax_cont_frac", m_shower_merge_relax_cont_frac);
        m_shower_merge_relax_cont_gap = get(config, "shower_merge_relax_cont_gap", m_shower_merge_relax_cont_gap);
        m_shower_merge_relax_cont_qmed = get(config, "shower_merge_relax_cont_qmed", m_shower_merge_relax_cont_qmed);
        m_shower_merge_relax_cont_axis = get(config, "shower_merge_relax_cont_axis", m_shower_merge_relax_cont_axis);
        m_shower_merge_relax_cont_dmax = get(config, "shower_merge_relax_cont_dmax", m_shower_merge_relax_cont_dmax);
        m_shower_merge_relax_cont_t1_gap = get(config, "shower_merge_relax_cont_t1_gap", m_shower_merge_relax_cont_t1_gap);
        m_shower_merge_relax_cont_t1_fold = get(config, "shower_merge_relax_cont_t1_fold", m_shower_merge_relax_cont_t1_fold);
        m_stem_backfill_back_guard = get(config, "stem_backfill_back_guard", m_stem_backfill_back_guard);
        m_stem_backfill_back_ang = get(config, "stem_backfill_back_ang", m_stem_backfill_back_ang);
        m_shower_ex1_dedup_rehome = get(config, "shower_ex1_dedup_rehome", m_shower_ex1_dedup_rehome);
        m_shower_pass4_prune_detached = get(config, "shower_pass4_prune_detached", m_shower_pass4_prune_detached);  // doc pr/123 r1
        m_shower_pass4_prune_gap = get(config, "shower_pass4_prune_gap", m_shower_pass4_prune_gap);  // doc pr/123 r1, cm
        m_shower_pass4_prune_gap2 = get(config, "shower_pass4_prune_gap2", m_shower_pass4_prune_gap2);  // doc pr/124 A, cm
        m_shower_pass4_prune2_ang = get(config, "shower_pass4_prune2_ang", m_shower_pass4_prune2_ang);  // doc pr/124 A, deg
        m_shower_pass4_prune2_mdqdx = get(config, "shower_pass4_prune2_mdqdx", m_shower_pass4_prune2_mdqdx);  // doc pr/124 A, x MIP
        m_shower_pass3_cone_guard_len = get(config, "shower_pass3_cone_guard_len", m_shower_pass3_cone_guard_len);  // doc pr/124 C, cm
        m_shower_samevtx_track_absorb = get(config, "shower_samevtx_track_absorb", m_shower_samevtx_track_absorb);  // doc pr/125
        m_shower_samevtx_absorb_gap = get(config, "shower_samevtx_absorb_gap", m_shower_samevtx_absorb_gap);  // doc pr/125, cm
        m_shower_samevtx_absorb_max_len = get(config, "shower_samevtx_absorb_max_len", m_shower_samevtx_absorb_max_len);  // doc pr/125, cm
        m_shower_samevtx_absorb_min_len = get(config, "shower_samevtx_absorb_min_len", m_shower_samevtx_absorb_min_len);  // doc pr/125, cm
        m_shower_satellite_absorb = get(config, "shower_satellite_absorb", m_shower_satellite_absorb);  // doc pr/125
        m_shower_split = get(config, "shower_split", m_shower_split);  // doc pr/138 B2; false = no pass
        m_shower_split_max_valley = get(config, "shower_split_max_valley", m_shower_split_max_valley);  // doc pr/138 B2; sec A5.4 knee
        m_shower_split_min_frac = get(config, "shower_split_min_frac", m_shower_split_min_frac);  // doc pr/138 B2; per-seed charge share floor
        m_shower_split_max_parts = get(config, "shower_split_max_parts", m_shower_split_max_parts);  // doc pr/138 B3; 2 = the measured-exact kernel
        m_shower_split_min_charge = get(config, "shower_split_min_charge", m_shower_split_min_charge);  // doc pr/138 B1; candidate charge floor (raw Fit::dQ)
        m_shower_split_min_nseg = get(config, "shower_split_min_nseg", m_shower_split_min_nseg);  // doc pr/138 B1; candidate member-count floor
        m_shower_split_bundle_gap = get(config, "shower_split_bundle_gap", m_shower_split_bundle_gap);  // doc pr/138 B3; single-linkage bundle gap
        m_shower_split_snap = get(config, "shower_split_snap", m_shower_split_snap);  // doc pr/138 B3; k>=3 bundle dominance floor
        m_shower_split_skip_shared = get(config, "shower_split_skip_shared", m_shower_split_skip_shared);  // doc pr/139 P1.1; refuse a component holding a segment another shower also owns
        m_shower_split_shed_shared = get(config, "shower_split_shed_shared", m_shower_split_shed_shared);  // doc pr/139 sec 15; shed an ENTIRELY co-owned refused component
        m_shower_split_max_seeds = get(config, "shower_split_max_seeds", m_shower_split_max_seeds);  // doc pr/139 sec 17; angular-maxima cap (shipped 4)
        m_shower_split_em_type_max_len = get(config, "shower_split_em_type_max_len", m_shower_split_em_type_max_len);  // doc pr/139 sec 25; cm; 0 = off
        m_shower_split_max_impact = get(config, "shower_split_max_impact", m_shower_split_max_impact);  // doc pr/139 P1.2; cm; 0 = no bound
        m_shower_split_em_start = get(config, "shower_split_em_start", m_shower_split_em_start);  // doc pr/139 P1.3; seed the daughter on its nearest EM-typed member
        m_shower_split_rehome = get(config, "shower_split_rehome", m_shower_split_rehome);  // doc pr/139 P1.4; offer an orphan daughter to the nearest larger EM shower
        m_shower_split_rehome_gap = get(config, "shower_split_rehome_gap", m_shower_split_rehome_gap);  // doc pr/139 P1.4; cm; max daughter->host 3-D gap
        m_shower_satellite_absorb_max_mev = get(config, "shower_satellite_absorb_max_mev", m_shower_satellite_absorb_max_mev);  // doc pr/125, MeV
        m_shower_satellite_absorb_host_mev = get(config, "shower_satellite_absorb_host_mev", m_shower_satellite_absorb_host_mev);  // doc pr/125, MeV
        m_shower_pass4_track_guard_len = get(config, "shower_pass4_track_guard_len", m_shower_pass4_track_guard_len);  // doc pr/123 r1, cm
        m_shower_pass4_prox_guard_len = get(config, "shower_pass4_prox_guard_len", m_shower_pass4_prox_guard_len);  // doc pr/130 item 1b, cm
        m_shower_pass3_backfill_guard_len = get(config, "shower_pass3_backfill_guard_len", m_shower_pass3_backfill_guard_len);  // doc pr/130 item 1b, cm
        m_stem_backfill_back_dvtx = get(config, "stem_backfill_back_dvtx", m_stem_backfill_back_dvtx);  // doc pr/130 item B, cm
        m_shower_pass4_prefilter_v1_escape = get(config, "shower_pass4_prefilter_v1_escape", m_shower_pass4_prefilter_v1_escape);  // doc pr/136 r2
        m_shower_pass4_prefilter_v1_max_v2 = get(config, "shower_pass4_prefilter_v1_max_v2", m_shower_pass4_prefilter_v1_max_v2);  // doc pr/136 r2, deg
        m_shower_pass4_prefilter_v1_max_dis = get(config, "shower_pass4_prefilter_v1_max_dis", m_shower_pass4_prefilter_v1_max_dis);  // doc pr/136 r3, cm
        m_pi0_mass_offset = get(config, "pi0_mass_offset", m_pi0_mass_offset);  // doc pr/132 K1, MeV
        m_pi0_assoc_angle_deg = get(config, "pi0_assoc_angle_deg", m_pi0_assoc_angle_deg);  // doc pr/132 K2, deg
        m_pi0_attached_partner_min_mev = get(config, "pi0_attached_partner_min_mev", m_pi0_attached_partner_min_mev);  // doc pr/132 K3, MeV
        m_pi0_nv_max_prongs = get(config, "pi0_nv_max_prongs", m_pi0_nv_max_prongs);  // doc pr/132 K5
        m_pi0_readmit_retyped = get(config, "pi0_readmit_retyped", m_pi0_readmit_retyped);  // doc pr/132 K7
        m_pi0_admit_type3 = get(config, "pi0_admit_type3", m_pi0_admit_type3);  // doc pr/132 K8
        m_pi0_crumb_assoc_mev = get(config, "pi0_crumb_assoc_mev", m_pi0_crumb_assoc_mev);  // doc pr/132 K9, MeV
        m_pi0_collinear_merge_deg = get(config, "pi0_collinear_merge_deg", m_pi0_collinear_merge_deg);  // doc pr/132 K12, deg
        m_pi0_nv_partner_min_mev = get(config, "pi0_nv_partner_min_mev", m_pi0_nv_partner_min_mev);  // doc pr/132 K13, MeV
        m_shower_em_collinear_deg = get(config, "shower_em_collinear_deg", m_shower_em_collinear_deg);  // doc pr/132 K16, deg
        m_shower_em_collinear_dis_cm = get(config, "shower_em_collinear_dis_cm", m_shower_em_collinear_dis_cm);  // doc pr/132 K16, cm
        m_shower_em_collinear_host_mev = get(config, "shower_em_collinear_host_mev", m_shower_em_collinear_host_mev);  // doc pr/132 K16, MeV
        m_shower_em_backext_perp_cm = get(config, "shower_em_backext_perp_cm", m_shower_em_backext_perp_cm);  // doc pr/132 K17, cm
        m_shower_em_backext_len_cm = get(config, "shower_em_backext_len_cm", m_shower_em_backext_len_cm);  // doc pr/132 K17, cm
        m_pi0_accept_merge_dis_cm = get(config, "pi0_accept_merge_dis_cm", m_pi0_accept_merge_dis_cm);  // doc pr/132 K18, cm
        m_pi0_bp_vertex_miss_cm = get(config, "pi0_bp_vertex_miss_cm", m_pi0_bp_vertex_miss_cm);  // doc pr/132 K19, cm
        m_pi0_admit_muon_showers = get(config, "pi0_admit_muon_showers", m_pi0_admit_muon_showers);  // doc pr/133 K20
        m_pi0_nc_sig_angle_deg = get(config, "pi0_nc_sig_angle_deg", m_pi0_nc_sig_angle_deg);  // doc pr/133 K21, deg
        m_pi0_nc_floor_mev = get(config, "pi0_nc_floor_mev", m_pi0_nc_floor_mev);  // doc pr/133 K21 v2, MeV
        m_pi0_nc_pf_assoc_deg = get(config, "pi0_nc_pf_assoc_deg", m_pi0_nc_pf_assoc_deg);  // doc pr/133 K21 v2.2, deg
        m_pi0_nc_frag_merge = get(config, "pi0_nc_frag_merge", m_pi0_nc_frag_merge);  // doc pr/134 K22
        m_pi0_pf_assoc_deg = get(config, "pi0_pf_assoc_deg", m_pi0_pf_assoc_deg);  // doc pr/134 K23, deg
        m_pi0_prefer_main_vertex = get(config, "pi0_prefer_main_vertex", m_pi0_prefer_main_vertex);  // doc pr/134 K24
        m_pi0_nv_max_vtx_shift_cm = get(config, "pi0_nv_max_vtx_shift_cm", m_pi0_nv_max_vtx_shift_cm);  // doc pr/132 K10, cm
        m_pi0_nv_mass_window_mev = get(config, "pi0_nv_mass_window_mev", m_pi0_nv_mass_window_mev);  // doc pr/132 K11, MeV
        m_kine_count_guard_freed = get(config, "kine_count_guard_freed", m_kine_count_guard_freed);  // doc pr/123 r2
        m_kine_guard_freed_impact = get(config, "kine_guard_freed_impact", m_kine_guard_freed_impact);  // doc pr/129
        m_kine_guard_freed_miss_deg = get(config, "kine_guard_freed_miss_deg", m_kine_guard_freed_miss_deg);  // doc pr/129
        m_kine_count_near_cross_cluster = get(config, "kine_count_near_cross_cluster", m_kine_count_near_cross_cluster);  // doc pr/128
        m_kine_near_gap = get(config, "kine_near_gap", m_kine_near_gap);  // doc pr/128
        m_kine_near_min_len = get(config, "kine_near_min_len", m_kine_near_min_len);  // doc pr/128
        m_kine_near_end_tol = get(config, "kine_near_end_tol", m_kine_near_end_tol);  // doc pr/128
        m_kine_near_kink_deg = get(config, "kine_near_kink_deg", m_kine_near_kink_deg);  // doc pr/128
        m_kine_near_pointing_impact = get(config, "kine_near_pointing_impact", m_kine_near_pointing_impact);  // doc pr/144
        m_kine_near_pointing_miss_deg = get(config, "kine_near_pointing_miss_deg", m_kine_near_pointing_miss_deg);  // doc pr/144
        m_kine_count_conn4_near = get(config, "kine_count_conn4_near", m_kine_count_conn4_near);  // doc pr/128
        m_kine_conn4_near_gap = get(config, "kine_conn4_near_gap", m_kine_conn4_near_gap);  // doc pr/128
        m_straight_cont_cross_cluster = get(config, "straight_cont_cross_cluster", m_straight_cont_cross_cluster);
        m_sccc_bridge_body = get(config, "sccc_bridge_body", m_sccc_bridge_body);
        m_sccc_max_gap = get(config, "sccc_max_gap", m_sccc_max_gap);
        m_sccc_kink_max = get(config, "sccc_kink_max", m_sccc_kink_max);
        m_sccc_gap_aligned = get(config, "sccc_gap_aligned", m_sccc_gap_aligned);
        m_sccc_kink_tight = get(config, "sccc_kink_tight", m_sccc_kink_tight);
        m_single_muon_proton_chain_veto = get(config, "single_muon_proton_chain_veto", m_single_muon_proton_chain_veto);
        m_single_muon_long_muon_claim = get(config, "single_muon_long_muon_claim", m_single_muon_long_muon_claim);
        m_pid_flag_reconcile = get(config, "pid_flag_reconcile", m_pid_flag_reconcile);
        m_long_muon_stub_bridge = get(config, "long_muon_stub_bridge", m_long_muon_stub_bridge);
        m_long_muon_stub_bridge_len = get(config, "long_muon_stub_bridge_len", m_long_muon_stub_bridge_len);  // doc 84 round 1 (P3), cm
        m_long_muon_angle_relax_long = get(config, "long_muon_angle_relax_long", m_long_muon_angle_relax_long);  // doc 84 round 1 (P2)
        m_long_muon_angle_relax_deg = get(config, "long_muon_angle_relax_deg", m_long_muon_angle_relax_deg);  // doc 84 round 1 (P2), deg
        if (m_vertex_z_prior_scale <= 0) {
            SPDLOG_LOGGER_WARN(s_log, "CheckBeamParticle: vertex_z_prior_scale must be > 0; keeping 200 cm");
            m_vertex_z_prior_scale = 200;
        }
        // Per-plane weights {U,V,W} (TaggerCheckNeutrino.cxx:389-403).
        if (config.isMember("kine_plane_weights")) {
            const auto& jv = config["kine_plane_weights"];
            if (!jv.isArray() || jv.size() != 3u) {
                SPDLOG_LOGGER_WARN(s_log, "CheckBeamParticle: kine_plane_weights must be a 3-element array; keeping default");
            }
            else {
                std::vector<double> v{jv[0].asDouble(), jv[1].asDouble(), jv[2].asDouble()};
                if (v[0] + v[1] + v[2] <= 0) {
                    SPDLOG_LOGGER_WARN(s_log, "CheckBeamParticle: kine_plane_weights must sum to > 0; keeping default");
                }
                else {
                    m_kine_plane_weights = v;
                }
            }
        }
    }

    // TaggerCheckNeutrino::default_configuration(): every kept key at its
    // C++ default, so the compiled config documents the operating point.
    void default_neutrino_knobs(Configuration& cfg) const {
        cfg["dir_weak_use_score"] = m_dir_weak_use_score;
        cfg["proton_dir_vote"] = m_proton_dir_vote;
        cfg["proton_dir_score_max"] = m_proton_dir_score_max;
        cfg["proton_dir_asym_min"] = m_proton_dir_asym_min;
        cfg["endpoint_trim_retry"] = m_endpoint_trim_retry;
        cfg["fit_vertex_min_seg_length"] = m_fit_vertex_min_seg_length;
        cfg["mvfit_robust"] = m_mvfit_robust;
        cfg["mvfit_main_only"] = m_mvfit_main_only;
        cfg["mvfit_min_len"] = m_mvfit_min_len;
        cfg["mvfit_rin_margin"] = m_mvfit_rin_margin;
        cfg["mvfit_rout_frac"] = m_mvfit_rout_frac;
        cfg["mvfit_rout_min"] = m_mvfit_rout_min;
        cfg["mvfit_rout_max"] = m_mvfit_rout_max;
        cfg["mvfit_angle"] = m_mvfit_angle;
        cfg["mvfit_min_pts"] = m_mvfit_min_pts;
        cfg["mvfit_min_aniso"] = m_mvfit_min_aniso;
        cfg["mvfit_prior_range"] = m_mvfit_prior_range;
        cfg["cathode_x"] = m_cathode_x;
        cfg["cathode_kink_xcut"] = m_cathode_kink_xcut;
        cfg["cathode_wide_kink_angle"] = m_cathode_wide_kink_angle;
        cfg["cathode_wide_kink_skirt"] = m_cathode_wide_kink_skirt;
        cfg["cathode_wide_kink_baseline"] = m_cathode_wide_kink_baseline;
        cfg["two_end_break"] = m_two_end_break;
        cfg["teb_min_len"] = m_teb_min_len;
        cfg["teb_min_arm"] = m_teb_min_arm;
        cfg["teb_min_arm_pts"] = m_teb_min_arm_pts;
        cfg["teb_stub_max"] = m_teb_stub_max;
        cfg["teb_accept_range"] = m_teb_accept_range;
        cfg["teb_rise_r1"] = m_teb_rise_r1;
        cfg["teb_rise_r2"] = m_teb_rise_r2;
        cfg["teb_abs_end_min"] = m_teb_abs_end_min;
        cfg["teb_dip_floor"] = m_teb_dip_floor;
        cfg["teb_score_cap_r1"] = m_teb_score_cap_r1;
        cfg["teb_score_cap_r2"] = m_teb_score_cap_r2;
        cfg["teb_turn_angle"] = m_teb_turn_angle;
        cfg["teb_turn_baseline"] = m_teb_turn_baseline;
        cfg["teb_turn_skirt"] = m_teb_turn_skirt;
        cfg["teb_turn_min_arm_frac"] = m_teb_turn_min_arm_frac;
        cfg["teb_bragg_veto_turn"] = m_teb_bragg_veto_turn;
        cfg["kink_walk_dqdx_stop"] = m_kink_walk_dqdx_stop;
        cfg["kink_break_protect"] = m_kink_break_protect;
        cfg["kink_dqdx_hot_ratio"] = m_kink_dqdx_hot_ratio;
        cfg["esva_ignore_empty_2d"] = m_esva_ignore_empty_2d;
        cfg["main_vertex_graph_audit"] = m_main_vertex_graph_audit;
        cfg["mvga_radius"] = m_mvga_radius;
        cfg["mvga_dup_tol"] = m_mvga_dup_tol;
        cfg["mvga_dup_frac"] = m_mvga_dup_frac;
        cfg["mvga_dup_angle"] = m_mvga_dup_angle;
        cfg["mvga_bridge_mip"] = m_mvga_bridge_mip;
        cfg["mvga_reconnect"] = m_mvga_reconnect;
        cfg["mvga_stub"] = m_mvga_stub;
        cfg["mvga_stub_pts"] = m_mvga_stub_pts;
        cfg["mvga_reseat_angle"] = m_mvga_reseat_angle;
        cfg["mvga_satellite"] = m_mvga_satellite;
        cfg["mvga_interposed"] = m_mvga_interposed;
        cfg["mvga_interposed_angle"] = m_mvga_interposed_angle;
        cfg["mvga_interposed_len"] = m_mvga_interposed_len;
        cfg["mvga_sat_dup_frac"] = m_mvga_sat_dup_frac;
        cfg["mvga_interposed_deg1"] = m_mvga_interposed_deg1;
        cfg["mvga_splice_straighten"] = m_mvga_splice_straighten;
        cfg["mvga_approach_collapse"] = m_mvga_approach_collapse;
        cfg["mvga_straighten_radius"] = m_mvga_straighten_radius;
        cfg["mvga_op1_radius"] = m_mvga_op1_radius;
        cfg["mvga_op1_dup_frac"] = m_mvga_op1_dup_frac;
        cfg["mvga_op1_post"] = m_mvga_op1_post;
        cfg["swap_orphan_dup_audit"] = m_swap_orphan_dup_audit;
        cfg["mvga_proj_dup_frac"] = m_mvga_proj_dup_frac;
        cfg["mvga_proj_dqdx_ratio"] = m_mvga_proj_dqdx_ratio;
        cfg["mvga_proj_angle"] = m_mvga_proj_angle;
        cfg["mvga_ac_chord_max"] = m_mvga_ac_chord_max;
        cfg["mvga_ac_no_cascade"] = m_mvga_ac_no_cascade;
        cfg["mvga_passthru"] = m_mvga_passthru;
        cfg["mvga_passthru_tol"] = m_mvga_passthru_tol;
        cfg["mvga_interposed_fallback"] = m_mvga_interposed_fallback;
        cfg["mvga_interposed_fallback_min_angle"] = m_mvga_interposed_fallback_min_angle;
        cfg["mvga_dup_starved_asym"] = m_mvga_dup_starved_asym;
        cfg["mvga_dup_starved_mip"] = m_mvga_dup_starved_mip;
        cfg["mvga_dup_starved_span"] = m_mvga_dup_starved_span;
        cfg["shower_topo_demote_len"] = m_shower_topo_demote_len;
        cfg["fit_blob_coverage_defer"] = m_fit_blob_coverage_defer;
        cfg["fit_exclusion"] = m_fit_exclusion;
        cfg["graph_endpoint_tol"] = m_graph_endpoint_tol;
        cfg["oov_prototype_parity"] = m_oov_prototype_parity;
        cfg["first_seg_local_pca"] = m_first_seg_local_pca;
        cfg["other_seg_relaxed_accept"] = m_other_seg_relaxed_accept;
        cfg["other_seg_empty_2d_guard"] = m_other_seg_empty_2d_guard;
        cfg["other_seg_keep_isolated"] = m_other_seg_keep_isolated;
        cfg["other_seg_keep_isolated_min_points"] = m_other_seg_keep_isolated_min_points;
        cfg["other_seg_keep_isolated_min_length"] = m_other_seg_keep_isolated_min_length;
        cfg["other_seg_keep_isolated_len_admit"] = m_other_seg_keep_isolated_len_admit;
        cfg["iso_snap_min_dir_mag"] = m_iso_snap_min_dir_mag;
        cfg["assoc_full_recluster"] = m_assoc_full_recluster;
        cfg["assoc_reassign_orphans"] = m_assoc_reassign_orphans;
        cfg["assoc_clear_on_merge"] = m_assoc_clear_on_merge;
        cfg["shower_topo_proto_dir"] = m_shower_topo_proto_dir;
        cfg["vertex_dir_use_fit_point"] = m_vertex_dir_use_fit_point;
        cfg["shower_traj_recheck_parity"] = m_shower_traj_recheck_parity;
        cfg["main_vertex_require_descriptor"] = m_main_vertex_require_descriptor;
        cfg["main_vertex_candidate_flag"] = m_main_vertex_candidate_flag;
        cfg["cont_muon_dir3_30cm"] = m_cont_muon_dir3_30cm;
        cfg["track_comp_empty_abstain"] = m_track_comp_empty_abstain;
        cfg["shower_topo_reset"] = m_shower_topo_reset;
        cfg["reclass_preserve_4mom"] = m_reclass_preserve_4mom;
        cfg["reclass_never_computed_ke_floor"] = m_reclass_never_computed_ke_floor;
        cfg["dir_track_median_local"] = m_dir_track_median_local;
        cfg["examine_showers_vertex_by_index"] = m_examine_showers_vertex_by_index;
        cfg["iso_endpoint"] = m_iso_endpoint;
        cfg["iso_endpoint_min_length"] = m_iso_endpoint_min_length;
        cfg["iso_endpoint_max_xext"] = m_iso_endpoint_max_xext;
        cfg["iso_endpoint_xext_frac"] = m_iso_endpoint_xext_frac;
        cfg["iso_endpoint_xext_quantile"] = m_iso_endpoint_xext_quantile;
        cfg["iso_endpoint_tube_radius"] = m_iso_endpoint_tube_radius;
        cfg["iso_endpoint_min_aspect"] = m_iso_endpoint_min_aspect;
        cfg["traj_cover_probe"] = m_traj_cover_probe;
        cfg["pr_find_other_rounds"] = m_pr_find_other_rounds;
        cfg["v3_extension_guard"] = m_v3_extension_guard;
        cfg["v3_extension_min_gain"] = m_v3_extension_min_gain;
        cfg["es3_stub_guard"] = m_es3_stub_guard;
        cfg["es3sg_stub_max"] = m_es3sg_stub_max;
        cfg["es3sg_len_ratio"] = m_es3sg_len_ratio;
        cfg["es3sg_ang3_min"] = m_es3sg_ang3_min;
        cfg["es3sg_ang_ratio"] = m_es3sg_ang_ratio;
        cfg["es3sg_require_terminal"] = m_es3sg_require_terminal;
        cfg["vertex_z_prior_scale"] = m_vertex_z_prior_scale;
        cfg["kine_fudge_factor"] = m_kine_fudge_factor;
        cfg["kine_recom_factor"] = m_kine_recom_factor;
        cfg["kine_shower_fudge_factor"] = m_kine_shower_fudge_factor;
        cfg["kine_shower_recom_factor"] = m_kine_shower_recom_factor;
        cfg["kine_proton_recom_factor"] = m_kine_proton_recom_factor;
        cfg["kine_plane_asym_switch"] = m_kine_plane_asym_switch;
        cfg["kine_w_value"] = m_kine_w_value;
        cfg["kine_shower_pdg_live"] = m_kine_shower_pdg_live;
        cfg["steiner_gap_penalty"] = m_steiner_gap_penalty;
        cfg["sgp_dead_alpha"] = m_sgp_dead_alpha;
        cfg["sgp_min_edge"] = m_sgp_min_edge;
        cfg["sgp_sample_step"] = m_sgp_sample_step;
        cfg["sgp_point_radius"] = m_sgp_point_radius;
        cfg["good_point_pitch_frac"] = m_good_point_pitch_frac;
        cfg["sgp_weak_scale"] = m_sgp_weak_scale;
        cfg["sgp_weak_qref"] = m_sgp_weak_qref;
        cfg["sgp_max_sep"] = m_sgp_max_sep;
        cfg["break_seg_orient"] = m_break_seg_orient;
        cfg["daughter_count_proto_main_vertex"] = m_daughter_count_proto_main_vertex;
        cfg["daughter_count_proto_examine_showers"] = m_daughter_count_proto_examine_showers;
        cfg["shower_pdg_from_start_segment"] = m_shower_pdg_from_start_segment;
        cfg["shower_pdg_from_shower_type"] = m_shower_pdg_from_shower_type;
        cfg["shower_pdg_exact_muon_test"] = m_shower_pdg_exact_muon_test;
        cfg["pi0_id_shared_allocator"] = m_pi0_id_shared_allocator;
        cfg["shower_flag_pdg_electron"] = m_shower_flag_pdg_electron;
        cfg["shower_less_id_tiebreak"] = m_shower_less_id_tiebreak;
        cfg["shower_endpoint_exclude_start_vertex"] = m_shower_endpoint_exclude_start_vertex;
        cfg["shower_endpoint_skip_orphan_vtx"] = m_shower_endpoint_skip_orphan_vtx;
        cfg["shower_walk_visited_parity"] = m_shower_walk_visited_parity;
        cfg["track_pid_persist_dqdx"] = m_track_pid_persist_dqdx;
        cfg["shower_reclass_dqdx_guard"] = m_shower_reclass_dqdx_guard;
        cfg["shower_topo_dqdx_guard"] = m_shower_topo_dqdx_guard;
        cfg["track_pid_persist_4mom"] = m_track_pid_persist_4mom;
        cfg["shower_proton_daughter_pion"] = m_shower_proton_daughter_pion;
        cfg["shower_proton_daughter_pion_dissolve"] = m_shower_proton_daughter_pion_dissolve;
        cfg["muon_multi_proton_pion"] = m_muon_multi_proton_pion;
        cfg["track_pid_persist_dqdx_electron_guard"] = m_track_pid_persist_dqdx_electron_guard;
        cfg["shower_connect_main_vertex_straight_guard"] = m_shower_connect_main_vertex_straight_guard;
        cfg["shower_traj_straight_guard"] = m_shower_traj_straight_guard;
        cfg["shower_absorb_track_guard"] = m_shower_absorb_track_guard;
        cfg["shower_absorb_unreachable_main"] = m_shower_absorb_unreachable_main;
        cfg["michel_stem_muon_rescue"] = m_michel_stem_muon_rescue;
        cfg["shower_in_cascade_guard"] = m_shower_in_cascade_guard;
        cfg["shower_in_max_len"] = m_shower_in_max_len;
        cfg["shower_in_mip_hi"] = m_shower_in_mip_hi;
        cfg["shower_connect_from_vertices_straight_guard"] = m_shower_connect_from_vertices_straight_guard;
        cfg["shower_connect_start_seg_straight_guard"] = m_shower_connect_start_seg_straight_guard;
        cfg["examine_direction_dirsign_shower_in_guard"] = m_examine_direction_dirsign_shower_in_guard;
        cfg["daughter_shower_angle_reclass_straight_guard"] = m_daughter_shower_angle_reclass_straight_guard;
        cfg["shower_topo_reexam_straight_guard"] = m_shower_topo_reexam_straight_guard;
        cfg["sfv_kink_max"] = m_sfv_kink_max;
        cfg["shower_nv_bridge_track"] = m_shower_nv_bridge_track;
        cfg["shower_nv_bridge_max_gap"] = m_shower_nv_bridge_max_gap;
        cfg["shower_nv_main_pi_init"] = m_shower_nv_main_pi_init;
        cfg["kine_drop_stray_satellites"] = m_kine_drop_stray_satellites;
        cfg["kine_sat_min_energy"] = m_kine_sat_min_energy;
        cfg["kine_sat_prox_max"] = m_kine_sat_prox_max;
        cfg["kine_sat_angle_bad"] = m_kine_sat_angle_bad;
        cfg["kine_sat_angle_main"] = m_kine_sat_angle_main;
        cfg["kine_sat_far_dis"] = m_kine_sat_far_dis;
        cfg["kine_sat_axis_dis_cut"] = m_kine_sat_axis_dis_cut;
        cfg["kine_sat_cont_kink"] = m_kine_sat_cont_kink;
        cfg["kine_sat_track_max_nseg"] = m_kine_sat_track_max_nseg;
        cfg["kine_sat_em_far_dis"] = m_kine_sat_em_far_dis;
        cfg["kine_sat_cont_keep_deg"] = m_kine_sat_cont_keep_deg;
        cfg["michel_stem_michel_check"] = m_michel_stem_michel_check;
        cfg["michel_stem_max_far_len"] = m_michel_stem_max_far_len;
        cfg["shower_stem_backfill"] = m_shower_stem_backfill;
        cfg["stem_backfill_max_len"] = m_stem_backfill_max_len;
        cfg["stem_backfill_mip_lo"] = m_stem_backfill_mip_lo;
        cfg["stem_backfill_mip_hi"] = m_stem_backfill_mip_hi;
        cfg["stem_backfill_min_shower_len"] = m_stem_backfill_min_shower_len;
        cfg["shower_conn3_unreachable"] = m_shower_conn3_unreachable;
        cfg["conn3_unreachable_min_len"] = m_conn3_unreachable_min_len;
        cfg["conn3_stitch_max"] = m_conn3_stitch_max;
        cfg["shower_dedup_start_seg"] = m_shower_dedup_start_seg;
        cfg["shower_traj_michel_stem"] = m_shower_traj_michel_stem;
        cfg["michel_stem_traj_min_len"] = m_michel_stem_traj_min_len;
        cfg["michel_stem_traj_max_len"] = m_michel_stem_traj_max_len;
        cfg["michel_stem_traj_mip_lo"] = m_michel_stem_traj_mip_lo;
        cfg["michel_stem_traj_max_far_len"] = m_michel_stem_traj_max_far_len;
        cfg["michel_stem_traj_min_kink_deg"] = m_michel_stem_traj_min_kink_deg;
        cfg["shower_long_muon_keep_type"] = m_shower_long_muon_keep_type;
        cfg["shower_bragg_protect_start_segment"] = m_shower_bragg_protect_start_segment;
        cfg["shower_reclass_case_b_dqdx_guard"] = m_shower_reclass_case_b_dqdx_guard;
        cfg["shower_accept_pid_guard"] = m_shower_accept_pid_guard;
        cfg["shower_pid_guard_min_len"] = m_shower_pid_guard_min_len;
        cfg["shower_vote_track_pid_counts"] = m_shower_vote_track_pid_counts;
        cfg["shower_cone_absorb_guard"] = m_shower_cone_absorb_guard;
        cfg["shower_detach_track_stem"] = m_shower_detach_track_stem;
        cfg["shower_ghost_member_drop"] = m_shower_ghost_member_drop;
        cfg["shower_ghost_overlap_frac"] = m_shower_ghost_overlap_frac;
        cfg["shower_ghost_dqdx_ratio"] = m_shower_ghost_dqdx_ratio;
        cfg["shower_ghost_min_len"] = m_shower_ghost_min_len;
        cfg["kine_charge_dedup"] = m_kine_charge_dedup;
        cfg["kine_charge_rebuild"] = m_kine_charge_rebuild;
        cfg["kine_charge_track_ctx"] = m_kine_charge_track_ctx;
        cfg["kine_mass_rules"] = m_kine_mass_rules;
        cfg["kine_hadronic_dqdx"] = m_kine_hadronic_dqdx;
        cfg["kine_long_muon_mode"] = m_kine_long_muon_mode;
        cfg["kine_long_muon_ratio_lo"] = m_kine_long_muon_ratio_lo;
        cfg["kine_long_muon_ratio_hi"] = m_kine_long_muon_ratio_hi;
        cfg["kine_dqdx_skip_zero_dx"] = m_kine_dqdx_skip_zero_dx;
        cfg["long_muon_range_empty_chain_fallback"] = m_long_muon_range_empty_chain_fallback;
        cfg["long_muon_members_geometry"] = m_long_muon_members_geometry;
        cfg["kine_mainvtx_used_guard"] = m_kine_mainvtx_used_guard;
        cfg["shower_hadronic_tag"] = m_shower_hadronic_tag;
        cfg["shower_hadronic_min_len"] = m_shower_hadronic_min_len;
        cfg["shower_hadronic_scan_len"] = m_shower_hadronic_scan_len;
        cfg["shower_hadronic_bin"] = m_shower_hadronic_bin;
        cfg["shower_hadronic_r_cyl"] = m_shower_hadronic_r_cyl;
        cfg["shower_hadronic_r_core"] = m_shower_hadronic_r_core;
        cfg["shower_hadronic_growth_max"] = m_shower_hadronic_growth_max;
        cfg["shower_hadronic_growth_bragg"] = m_shower_hadronic_growth_bragg;
        cfg["shower_hadronic_bragg_ratio"] = m_shower_hadronic_bragg_ratio;
        cfg["shower_hadronic_stem_ratio"] = m_shower_hadronic_stem_ratio;
        cfg["kine_count_orphan_tracks"] = m_kine_count_orphan_tracks;
        cfg["kine_orphan_track_min"] = m_kine_orphan_track_min;
        cfg["shower_pass4_best_owner"] = m_shower_pass4_best_owner;
        cfg["shower_merge_relax"] = m_shower_merge_relax;
        cfg["shower_merge_relax_dis"] = m_shower_merge_relax_dis;
        cfg["shower_merge_relax_angle"] = m_shower_merge_relax_angle;
        cfg["shower_merge_relax_min_len"] = m_shower_merge_relax_min_len;
        cfg["shower_merge_relax_continuity"] = m_shower_merge_relax_continuity;
        cfg["shower_merge_relax_cont_frac"] = m_shower_merge_relax_cont_frac;
        cfg["shower_merge_relax_cont_gap"] = m_shower_merge_relax_cont_gap;
        cfg["shower_merge_relax_cont_qmed"] = m_shower_merge_relax_cont_qmed;
        cfg["shower_merge_relax_cont_axis"] = m_shower_merge_relax_cont_axis;
        cfg["shower_merge_relax_cont_dmax"] = m_shower_merge_relax_cont_dmax;
        cfg["shower_merge_relax_cont_t1_gap"] = m_shower_merge_relax_cont_t1_gap;
        cfg["shower_merge_relax_cont_t1_fold"] = m_shower_merge_relax_cont_t1_fold;
        cfg["stem_backfill_back_guard"] = m_stem_backfill_back_guard;
        cfg["stem_backfill_back_ang"] = m_stem_backfill_back_ang;
        cfg["shower_ex1_dedup_rehome"] = m_shower_ex1_dedup_rehome;
        cfg["shower_pass4_prune_detached"] = m_shower_pass4_prune_detached;
        cfg["shower_pass4_prune_gap"] = m_shower_pass4_prune_gap;
        cfg["shower_pass4_prune_gap2"] = m_shower_pass4_prune_gap2;
        cfg["shower_pass4_prune2_ang"] = m_shower_pass4_prune2_ang;
        cfg["shower_pass4_prune2_mdqdx"] = m_shower_pass4_prune2_mdqdx;
        cfg["shower_pass3_cone_guard_len"] = m_shower_pass3_cone_guard_len;
        cfg["shower_samevtx_track_absorb"] = m_shower_samevtx_track_absorb;
        cfg["shower_samevtx_absorb_gap"] = m_shower_samevtx_absorb_gap;
        cfg["shower_samevtx_absorb_max_len"] = m_shower_samevtx_absorb_max_len;
        cfg["shower_samevtx_absorb_min_len"] = m_shower_samevtx_absorb_min_len;
        cfg["shower_satellite_absorb"] = m_shower_satellite_absorb;
        cfg["shower_split"] = m_shower_split;
        cfg["shower_split_max_valley"] = m_shower_split_max_valley;
        cfg["shower_split_min_frac"] = m_shower_split_min_frac;
        cfg["shower_split_max_parts"] = m_shower_split_max_parts;
        cfg["shower_split_min_charge"] = m_shower_split_min_charge;
        cfg["shower_split_min_nseg"] = m_shower_split_min_nseg;
        cfg["shower_split_bundle_gap"] = m_shower_split_bundle_gap;
        cfg["shower_split_snap"] = m_shower_split_snap;
        cfg["shower_split_skip_shared"] = m_shower_split_skip_shared;
        cfg["shower_split_shed_shared"] = m_shower_split_shed_shared;
        cfg["shower_split_max_seeds"] = m_shower_split_max_seeds;
        cfg["shower_split_em_type_max_len"] = m_shower_split_em_type_max_len;
        cfg["shower_split_max_impact"] = m_shower_split_max_impact;
        cfg["shower_split_em_start"] = m_shower_split_em_start;
        cfg["shower_split_rehome"] = m_shower_split_rehome;
        cfg["shower_split_rehome_gap"] = m_shower_split_rehome_gap;
        cfg["shower_satellite_absorb_max_mev"] = m_shower_satellite_absorb_max_mev;
        cfg["shower_satellite_absorb_host_mev"] = m_shower_satellite_absorb_host_mev;
        cfg["shower_pass4_track_guard_len"] = m_shower_pass4_track_guard_len;
        cfg["shower_pass4_prox_guard_len"] = m_shower_pass4_prox_guard_len;
        cfg["shower_pass3_backfill_guard_len"] = m_shower_pass3_backfill_guard_len;
        cfg["stem_backfill_back_dvtx"] = m_stem_backfill_back_dvtx;
        cfg["shower_pass4_prefilter_v1_escape"] = m_shower_pass4_prefilter_v1_escape;
        cfg["shower_pass4_prefilter_v1_max_v2"] = m_shower_pass4_prefilter_v1_max_v2;
        cfg["shower_pass4_prefilter_v1_max_dis"] = m_shower_pass4_prefilter_v1_max_dis;
        cfg["pi0_mass_offset"] = m_pi0_mass_offset;
        cfg["pi0_assoc_angle_deg"] = m_pi0_assoc_angle_deg;
        cfg["pi0_attached_partner_min_mev"] = m_pi0_attached_partner_min_mev;
        cfg["pi0_nv_max_prongs"] = m_pi0_nv_max_prongs;
        cfg["pi0_readmit_retyped"] = m_pi0_readmit_retyped;
        cfg["pi0_admit_type3"] = m_pi0_admit_type3;
        cfg["pi0_crumb_assoc_mev"] = m_pi0_crumb_assoc_mev;
        cfg["pi0_collinear_merge_deg"] = m_pi0_collinear_merge_deg;
        cfg["pi0_nv_partner_min_mev"] = m_pi0_nv_partner_min_mev;
        cfg["shower_em_collinear_deg"] = m_shower_em_collinear_deg;
        cfg["shower_em_collinear_dis_cm"] = m_shower_em_collinear_dis_cm;
        cfg["shower_em_collinear_host_mev"] = m_shower_em_collinear_host_mev;
        cfg["shower_em_backext_perp_cm"] = m_shower_em_backext_perp_cm;
        cfg["shower_em_backext_len_cm"] = m_shower_em_backext_len_cm;
        cfg["pi0_accept_merge_dis_cm"] = m_pi0_accept_merge_dis_cm;
        cfg["pi0_bp_vertex_miss_cm"] = m_pi0_bp_vertex_miss_cm;
        cfg["pi0_admit_muon_showers"] = m_pi0_admit_muon_showers;
        cfg["pi0_nc_sig_angle_deg"] = m_pi0_nc_sig_angle_deg;
        cfg["pi0_nc_floor_mev"] = m_pi0_nc_floor_mev;
        cfg["pi0_nc_pf_assoc_deg"] = m_pi0_nc_pf_assoc_deg;
        cfg["pi0_nc_frag_merge"] = m_pi0_nc_frag_merge;
        cfg["pi0_pf_assoc_deg"] = m_pi0_pf_assoc_deg;
        cfg["pi0_prefer_main_vertex"] = m_pi0_prefer_main_vertex;
        cfg["pi0_nv_max_vtx_shift_cm"] = m_pi0_nv_max_vtx_shift_cm;
        cfg["pi0_nv_mass_window_mev"] = m_pi0_nv_mass_window_mev;
        cfg["kine_count_guard_freed"] = m_kine_count_guard_freed;
        cfg["kine_guard_freed_impact"] = m_kine_guard_freed_impact;
        cfg["kine_guard_freed_miss_deg"] = m_kine_guard_freed_miss_deg;
        cfg["kine_count_near_cross_cluster"] = m_kine_count_near_cross_cluster;
        cfg["kine_near_gap"] = m_kine_near_gap;
        cfg["kine_near_min_len"] = m_kine_near_min_len;
        cfg["kine_near_end_tol"] = m_kine_near_end_tol;
        cfg["kine_near_kink_deg"] = m_kine_near_kink_deg;
        cfg["kine_near_pointing_impact"] = m_kine_near_pointing_impact;
        cfg["kine_near_pointing_miss_deg"] = m_kine_near_pointing_miss_deg;
        cfg["kine_count_conn4_near"] = m_kine_count_conn4_near;
        cfg["kine_conn4_near_gap"] = m_kine_conn4_near_gap;
        cfg["straight_cont_cross_cluster"] = m_straight_cont_cross_cluster;
        cfg["sccc_bridge_body"] = m_sccc_bridge_body;
        cfg["sccc_max_gap"] = m_sccc_max_gap;
        cfg["sccc_kink_max"] = m_sccc_kink_max;
        cfg["sccc_gap_aligned"] = m_sccc_gap_aligned;
        cfg["sccc_kink_tight"] = m_sccc_kink_tight;
        cfg["single_muon_proton_chain_veto"] = m_single_muon_proton_chain_veto;
        cfg["single_muon_long_muon_claim"] = m_single_muon_long_muon_claim;
        cfg["pid_flag_reconcile"] = m_pid_flag_reconcile;
        cfg["long_muon_stub_bridge"] = m_long_muon_stub_bridge;
        cfg["long_muon_stub_bridge_len"] = m_long_muon_stub_bridge_len;
        cfg["long_muon_angle_relax_long"] = m_long_muon_angle_relax_long;
        cfg["long_muon_angle_relax_deg"] = m_long_muon_angle_relax_deg;
        Json::Value pw(Json::arrayValue);
        for (double w : m_kine_plane_weights) pw.append(w);
        cfg["kine_plane_weights"] = pw;
    }

    // TaggerCheckNeutrino.cxx:3104-3553 verbatim (the kept members): copy the
    // knobs onto the per-visit PatternAlgorithms with TCN's unit conventions.
    void apply_pattern_knobs(PatternAlgorithms& pa) const {
        pa.m_perf = m_perf;
        pa.m_mip_dqdx        = m_mip_dqdx / units::cm;         // e/cm -> internal
        pa.m_mip_dqdx_median = m_mip_dqdx_median / units::cm;  // e/cm -> internal
        pa.m_dir_weak_use_score = m_dir_weak_use_score;
        pa.m_proton_dir_vote = m_proton_dir_vote;
        pa.m_proton_dir_score_max = m_proton_dir_score_max;
        pa.m_proton_dir_asym_min = m_proton_dir_asym_min;
        pa.m_endpoint_trim_retry = m_endpoint_trim_retry;
        pa.m_fit_vertex_min_seg_length = m_fit_vertex_min_seg_length * units::cm;  // cm -> internal
        pa.m_mvfit_robust = m_mvfit_robust;
        pa.m_mvfit_main_only = m_mvfit_main_only;
        pa.m_mvfit_min_len = m_mvfit_min_len * units::cm;  // cm -> internal
        pa.m_mvfit_rin_margin = m_mvfit_rin_margin * units::cm;  // cm -> internal
        pa.m_mvfit_rout_frac = m_mvfit_rout_frac;  // unitless
        pa.m_mvfit_rout_min = m_mvfit_rout_min * units::cm;  // cm -> internal
        pa.m_mvfit_rout_max = m_mvfit_rout_max * units::cm;  // cm -> internal
        pa.m_mvfit_angle = m_mvfit_angle;  // deg, no conversion
        pa.m_mvfit_min_pts = m_mvfit_min_pts;  // count
        pa.m_mvfit_min_aniso = m_mvfit_min_aniso;  // unitless
        pa.m_mvfit_prior_range = m_mvfit_prior_range * units::cm;  // cm -> internal
        pa.m_cathode_x = m_cathode_x * units::cm;  // cm -> internal
        pa.m_cathode_kink_xcut = m_cathode_kink_xcut * units::cm;  // cm -> internal
        pa.m_cathode_wide_kink_angle = m_cathode_wide_kink_angle;  // deg, no conversion
        pa.m_cathode_wide_kink_skirt = m_cathode_wide_kink_skirt * units::cm;  // cm -> internal
        pa.m_cathode_wide_kink_baseline = m_cathode_wide_kink_baseline * units::cm;  // cm -> internal
        pa.m_two_end_break = m_two_end_break;
        pa.m_teb_min_len = m_teb_min_len * units::cm;  // cm -> internal
        pa.m_teb_min_arm = m_teb_min_arm * units::cm;  // cm -> internal
        pa.m_teb_min_arm_pts = m_teb_min_arm_pts;
        pa.m_teb_stub_max = m_teb_stub_max * units::cm;  // cm -> internal
        pa.m_teb_accept_range = m_teb_accept_range * units::cm;  // cm -> internal
        pa.m_teb_rise_r1 = m_teb_rise_r1;
        pa.m_teb_rise_r2 = m_teb_rise_r2;
        pa.m_teb_abs_end_min = m_teb_abs_end_min;
        pa.m_teb_dip_floor = m_teb_dip_floor;
        pa.m_teb_score_cap_r1 = m_teb_score_cap_r1;
        pa.m_teb_score_cap_r2 = m_teb_score_cap_r2;
        pa.m_teb_turn_angle = m_teb_turn_angle;  // deg, no conversion
        pa.m_teb_turn_baseline = m_teb_turn_baseline * units::cm;  // cm -> internal
        pa.m_teb_turn_skirt = m_teb_turn_skirt * units::cm;  // cm -> internal
        pa.m_teb_turn_min_arm_frac = m_teb_turn_min_arm_frac;  // dimensionless, no conversion
        pa.m_teb_bragg_veto_turn = m_teb_bragg_veto_turn;  // deg, no conversion
        pa.m_kink_walk_dqdx_stop = m_kink_walk_dqdx_stop;
        pa.m_kink_break_protect = m_kink_break_protect;
        pa.m_kink_dqdx_hot_ratio = m_kink_dqdx_hot_ratio;
        pa.m_esva_ignore_empty_2d = m_esva_ignore_empty_2d;
        pa.m_main_vertex_graph_audit = m_main_vertex_graph_audit;
        pa.m_mvga_radius = m_mvga_radius * units::cm;  // cm -> internal
        pa.m_mvga_dup_tol = m_mvga_dup_tol * units::cm;  // cm -> internal
        pa.m_mvga_dup_frac = m_mvga_dup_frac;  // fraction, no conversion
        pa.m_mvga_dup_angle = m_mvga_dup_angle;  // deg, no conversion
        pa.m_mvga_bridge_mip = m_mvga_bridge_mip;  // x mip median, no conversion
        pa.m_mvga_reconnect = m_mvga_reconnect * units::cm;  // cm -> internal
        pa.m_mvga_stub = m_mvga_stub * units::cm;  // cm -> internal
        pa.m_mvga_stub_pts = m_mvga_stub_pts;
        pa.m_mvga_reseat_angle = m_mvga_reseat_angle;  // deg, no conversion
        pa.m_mvga_satellite = m_mvga_satellite * units::cm;  // cm -> internal
        pa.m_mvga_interposed = m_mvga_interposed;  // doc pr/85
        pa.m_mvga_interposed_angle = m_mvga_interposed_angle;  // deg, no conversion
        pa.m_mvga_interposed_len = m_mvga_interposed_len * units::cm;  // cm -> internal (doc pr/86)
        pa.m_mvga_sat_dup_frac = m_mvga_sat_dup_frac;  // fraction, no conversion (doc pr/86)
        pa.m_mvga_interposed_deg1 = m_mvga_interposed_deg1;  // doc pr/86
        pa.m_mvga_splice_straighten = m_mvga_splice_straighten * units::cm;  // cm -> internal (doc pr/86 round 2)
        pa.m_mvga_approach_collapse = m_mvga_approach_collapse * units::cm;  // cm -> internal (doc pr/86 round 2)
        pa.m_mvga_straighten_radius = m_mvga_straighten_radius * units::cm;  // cm -> internal (doc pr/86 round 2)
        pa.m_mvga_op1_radius = m_mvga_op1_radius * units::cm;  // cm -> internal; 0 and the -1 sentinel both survive the scale (doc pr/83 r3)
        pa.m_mvga_op1_dup_frac = m_mvga_op1_dup_frac;  // fraction, no conversion (doc pr/83 r3)
        pa.m_mvga_op1_post = m_mvga_op1_post;  // doc pr/83 r3 (class A)
        pa.m_swap_orphan_dup_audit = m_swap_orphan_dup_audit;  // doc pr/83 r3 (Mechanism C)
        pa.m_mvga_proj_dup_frac = m_mvga_proj_dup_frac;  // fraction, no conversion (doc pr/83 r4)
        pa.m_mvga_proj_dqdx_ratio = m_mvga_proj_dqdx_ratio;  // ratio, no conversion (doc pr/83 r4)
        pa.m_mvga_proj_angle = m_mvga_proj_angle;  // deg, no conversion (doc pr/83 r4b)
        pa.m_mvga_ac_chord_max = m_mvga_ac_chord_max * units::cm;  // cm -> internal (doc pr/99 round 2)
        pa.m_mvga_ac_no_cascade = m_mvga_ac_no_cascade;  // doc pr/99 round 2
        pa.m_mvga_passthru = m_mvga_passthru * units::cm;  // cm -> internal (doc pr/103)
        pa.m_mvga_passthru_tol = m_mvga_passthru_tol * units::cm;  // cm -> internal (doc pr/103)
        pa.m_mvga_interposed_fallback = m_mvga_interposed_fallback;  // doc pr/103
        pa.m_mvga_interposed_fallback_min_angle = m_mvga_interposed_fallback_min_angle;  // deg (doc pr/103)
        pa.m_mvga_dup_starved_asym = m_mvga_dup_starved_asym;  // ratio, no conversion (doc pr/99 round 2)
        pa.m_mvga_dup_starved_mip = m_mvga_dup_starved_mip;  // ratio, no conversion (doc pr/99 round 2)
        pa.m_mvga_dup_starved_span = m_mvga_dup_starved_span;  // ratio, no conversion (doc pr/99 round 2)
        pa.m_steiner_gap_penalty = m_steiner_gap_penalty;
        pa.m_sgp_dead_alpha = m_sgp_dead_alpha;  // fraction, no conversion
        pa.m_sgp_min_edge = m_sgp_min_edge * units::cm;  // cm -> internal
        pa.m_sgp_sample_step = m_sgp_sample_step * units::cm;  // cm -> internal
        pa.m_sgp_point_radius = m_sgp_point_radius * units::cm;  // cm -> internal
        pa.m_good_point_pitch_frac = m_good_point_pitch_frac;  // doc pdvd/32 round 3: fraction, no conversion
        pa.m_sgp_weak_scale = m_sgp_weak_scale;
        pa.m_sgp_weak_qref = m_sgp_weak_qref;
        pa.m_sgp_max_sep = m_sgp_max_sep * units::cm;  // cm -> internal
        pa.m_break_seg_orient = m_break_seg_orient;
        pa.m_shower_topo_demote_len = m_shower_topo_demote_len * units::cm;  // cm -> internal
        pa.m_fit_exclusion = m_fit_exclusion;
        pa.m_graph_endpoint_tol = m_graph_endpoint_tol * units::cm;  // cm -> internal
        pa.m_oov_prototype_parity = m_oov_prototype_parity;
        pa.m_first_seg_local_pca = m_first_seg_local_pca;
        pa.m_other_seg_relaxed_accept = m_other_seg_relaxed_accept;
        pa.m_other_seg_empty_2d_guard = m_other_seg_empty_2d_guard;
        pa.m_other_seg_keep_isolated = m_other_seg_keep_isolated;
        pa.m_other_seg_keep_isolated_min_points = m_other_seg_keep_isolated_min_points;
        pa.m_other_seg_keep_isolated_min_length = m_other_seg_keep_isolated_min_length * units::cm;  // cm -> internal
        pa.m_other_seg_keep_isolated_len_admit = m_other_seg_keep_isolated_len_admit * units::cm;  // cm -> internal
        pa.m_iso_snap_min_dir_mag = m_iso_snap_min_dir_mag * units::cm;  // cm -> internal
        pa.m_assoc_full_recluster = m_assoc_full_recluster;
        pa.m_assoc_reassign_orphans = m_assoc_reassign_orphans;
        pa.m_assoc_clear_on_merge = m_assoc_clear_on_merge;
        pa.m_shower_topo_proto_dir = m_shower_topo_proto_dir;
        pa.m_vertex_dir_use_fit_point = m_vertex_dir_use_fit_point;
        pa.m_shower_traj_recheck_parity = m_shower_traj_recheck_parity;
        pa.m_main_vertex_require_descriptor = m_main_vertex_require_descriptor;
        pa.m_main_vertex_candidate_flag = m_main_vertex_candidate_flag;
        pa.m_cont_muon_dir3_30cm = m_cont_muon_dir3_30cm;
        pa.m_track_comp_empty_abstain = m_track_comp_empty_abstain;
        pa.m_shower_topo_reset = m_shower_topo_reset;
        pa.m_reclass_preserve_4mom = m_reclass_preserve_4mom;
        pa.m_reclass_never_computed_ke_floor = m_reclass_never_computed_ke_floor;  // doc pr/40 round 2 F6
        pa.m_dir_track_median_local = m_dir_track_median_local;
        pa.m_examine_showers_vertex_by_index = m_examine_showers_vertex_by_index;
        pa.m_iso_endpoint = m_iso_endpoint;
        pa.m_iso_endpoint_min_length = m_iso_endpoint_min_length * units::cm;  // cm -> internal
        pa.m_iso_endpoint_max_xext = m_iso_endpoint_max_xext * units::cm;  // cm -> internal
        pa.m_iso_endpoint_xext_frac = m_iso_endpoint_xext_frac;
        pa.m_iso_endpoint_xext_quantile = m_iso_endpoint_xext_quantile;
        pa.m_iso_endpoint_tube_radius = m_iso_endpoint_tube_radius * units::cm;  // cm -> internal
        pa.m_iso_endpoint_min_aspect = m_iso_endpoint_min_aspect;
        pa.m_traj_cover_probe = m_traj_cover_probe;
        pa.m_es3_stub_guard = m_es3_stub_guard;
        pa.m_es3sg_stub_max = m_es3sg_stub_max * units::cm;  // cm -> internal
        pa.m_es3sg_len_ratio = m_es3sg_len_ratio;
        pa.m_es3sg_ang3_min = m_es3sg_ang3_min;
        pa.m_es3sg_ang_ratio = m_es3sg_ang_ratio;
        pa.m_es3sg_require_terminal = m_es3sg_require_terminal;
        pa.m_pr_find_other_rounds = m_pr_find_other_rounds;
        pa.m_v3_extension_guard = m_v3_extension_guard;
        pa.m_v3_extension_min_gain = m_v3_extension_min_gain * units::cm;  // cm -> internal
        pa.m_vertex_z_prior_scale = m_vertex_z_prior_scale * units::cm;
        pa.m_kine_charge.fudge_factor = m_kine_fudge_factor;
        pa.m_kine_charge.recom_factor = m_kine_recom_factor;
        pa.m_kine_charge.shower_fudge_factor = m_kine_shower_fudge_factor;
        pa.m_kine_charge.shower_recom_factor = m_kine_shower_recom_factor;
        pa.m_kine_charge.proton_recom_factor = m_kine_proton_recom_factor;
        pa.m_kine_charge.plane_asym_switch = m_kine_plane_asym_switch;
        pa.m_kine_charge.shower_pdg_live = m_kine_shower_pdg_live;
        pa.m_kine_charge.w_value = m_kine_w_value;
        pa.m_kine_charge.dedup = m_kine_charge_dedup;  // doc pr/99 r3 C1
        pa.m_kine_charge.rebuild = m_kine_charge_rebuild;  // doc pr/99 r3 C1b
        pa.m_kine_charge.track_ctx = m_kine_charge_track_ctx;  // doc pr/101 K1
        pa.m_kine_charge.mass_rules = m_kine_mass_rules;  // doc pr/101 K2
        pa.m_kine_charge.hadronic_dqdx = m_kine_hadronic_dqdx;  // doc pr/101 K3
        pa.m_kine_charge.long_muon_mode = m_kine_long_muon_mode;  // doc pr/101 K4
        pa.m_kine_charge.long_muon_ratio_lo = m_kine_long_muon_ratio_lo;
        pa.m_kine_charge.long_muon_ratio_hi = m_kine_long_muon_ratio_hi;
        pa.m_kine_charge.dqdx_skip_zero_dx = m_kine_dqdx_skip_zero_dx;  // doc pdvd/45 sec 5.4
        pa.m_kine_charge.long_muon_range_fallback = m_long_muon_range_empty_chain_fallback;  // doc 84 round 1 (P1)
        pa.m_kine_charge.long_muon_members_geometry = m_long_muon_members_geometry;  // doc 84 round 2
        pa.m_kine_charge.mainvtx_used_guard = m_kine_mainvtx_used_guard;  // doc pr/101 K5
        pa.m_daughter_count_proto_main_vertex = m_daughter_count_proto_main_vertex;
        pa.m_daughter_count_proto_examine_showers = m_daughter_count_proto_examine_showers;
        pa.m_shower_pdg_from_start_segment = m_shower_pdg_from_start_segment;
        pa.m_shower_pdg_from_shower_type = m_shower_pdg_from_shower_type;
        pa.m_shower_pdg_exact_muon_test = m_shower_pdg_exact_muon_test;
        pa.m_pi0_id_shared_allocator = m_pi0_id_shared_allocator;
        pa.m_shower_flag_pdg_electron = m_shower_flag_pdg_electron;
        pa.m_shower_less_id_tiebreak = m_shower_less_id_tiebreak;
        pa.m_shower_endpoint_exclude_start_vertex = m_shower_endpoint_exclude_start_vertex;
        pa.m_shower_endpoint_skip_orphan_vtx = m_shower_endpoint_skip_orphan_vtx;
        pa.m_shower_walk_visited_parity = m_shower_walk_visited_parity;
        pa.m_track_pid_persist_dqdx = m_track_pid_persist_dqdx;  // F1: threaded via track_pid_options()
        pa.m_shower_reclass_dqdx_guard = m_shower_reclass_dqdx_guard;  // F2
        pa.m_shower_topo_dqdx_guard = m_shower_topo_dqdx_guard;  // F3
        pa.m_track_pid_persist_4mom = m_track_pid_persist_4mom;  // F4: threaded via track_pid_options()
        pa.m_shower_proton_daughter_pion = m_shower_proton_daughter_pion;  // F5
        pa.m_shower_proton_daughter_pion_dissolve = m_shower_proton_daughter_pion_dissolve;  // F7
        pa.m_muon_multi_proton_pion = m_muon_multi_proton_pion;  // F8
        pa.m_track_pid_persist_dqdx_electron_guard = m_track_pid_persist_dqdx_electron_guard;  // F9
        pa.m_shower_connect_main_vertex_straight_guard = m_shower_connect_main_vertex_straight_guard;  // F10
        pa.m_shower_traj_straight_guard = m_shower_traj_straight_guard;  // F11
        pa.m_shower_absorb_track_guard = m_shower_absorb_track_guard;  // F12
        pa.m_shower_absorb_unreachable_main = m_shower_absorb_unreachable_main;  // doc pr/65 round 3
        pa.m_michel_stem_muon_rescue = m_michel_stem_muon_rescue;  // F14
        pa.m_shower_in_cascade_guard = m_shower_in_cascade_guard;  // pr/74 P1
        pa.m_shower_in_max_len = m_shower_in_max_len * units::cm;  // pr/74 P1
        pa.m_shower_in_mip_hi = m_shower_in_mip_hi;  // pr/74 P1
        pa.m_shower_connect_from_vertices_straight_guard = m_shower_connect_from_vertices_straight_guard;  // pr/40 r9 (r8 Part A)
        pa.m_shower_connect_start_seg_straight_guard = m_shower_connect_start_seg_straight_guard;  // pr/40 r9 (r7 c2c)
        pa.m_examine_direction_dirsign_shower_in_guard = m_examine_direction_dirsign_shower_in_guard;  // pr/40 r9 (r7 c2a)
        pa.m_daughter_shower_angle_reclass_straight_guard = m_daughter_shower_angle_reclass_straight_guard;  // pr/40 r9 (r7 c2b)
        pa.m_shower_topo_reexam_straight_guard = m_shower_topo_reexam_straight_guard;  // pr/40 r9 (r7 c1)
        pa.m_sfv_kink_max = m_sfv_kink_max;  // pr/40 r9 (degrees)
        pa.m_shower_nv_bridge_track = m_shower_nv_bridge_track;  // pr/40 r9 B2
        pa.m_shower_nv_bridge_max_gap = m_shower_nv_bridge_max_gap * units::cm;  // pr/40 r9 B2
        pa.m_shower_nv_main_pi_init = m_shower_nv_main_pi_init;  // pr/97 D1
        pa.m_kine_drop_stray_satellites = m_kine_drop_stray_satellites;  // pr/92
        pa.m_kine_sat_min_energy = m_kine_sat_min_energy * units::MeV;  // pr/92
        pa.m_kine_sat_prox_max = m_kine_sat_prox_max * units::cm;  // pr/92
        pa.m_kine_sat_angle_bad = m_kine_sat_angle_bad;  // pr/92 (degrees)
        pa.m_kine_sat_angle_main = m_kine_sat_angle_main;  // pr/92 (degrees)
        pa.m_kine_sat_far_dis = m_kine_sat_far_dis * units::cm;  // pr/92
        pa.m_kine_sat_axis_dis_cut = m_kine_sat_axis_dis_cut * units::cm;  // pr/92
        pa.m_kine_sat_cont_kink = m_kine_sat_cont_kink;  // pr/92 (degrees)
        pa.m_kine_sat_track_max_nseg = static_cast<int>(m_kine_sat_track_max_nseg);  // pr/92 r2 (count)
        pa.m_kine_sat_em_far_dis = m_kine_sat_em_far_dis * units::cm;  // pr/92 r2
        pa.m_kine_sat_cont_keep_deg = m_kine_sat_cont_keep_deg;  // pr/146 (degrees)
        pa.m_michel_stem_michel_check = m_michel_stem_michel_check;  // pr/74 P2
        pa.m_michel_stem_max_far_len = m_michel_stem_max_far_len * units::cm;  // pr/74 P2
        pa.m_shower_stem_backfill = m_shower_stem_backfill;  // pr/74 K4
        pa.m_stem_backfill_max_len = m_stem_backfill_max_len * units::cm;  // pr/74 K4
        pa.m_stem_backfill_mip_lo = m_stem_backfill_mip_lo;  // pr/74 K4
        pa.m_stem_backfill_mip_hi = m_stem_backfill_mip_hi;  // pr/74 K4
        pa.m_stem_backfill_min_shower_len = m_stem_backfill_min_shower_len * units::cm;  // pr/74 K4
        pa.m_shower_conn3_unreachable = m_shower_conn3_unreachable;  // pr/74 K5
        pa.m_conn3_unreachable_min_len = m_conn3_unreachable_min_len * units::cm;  // pr/74 K5
        pa.m_conn3_stitch_max = m_conn3_stitch_max * units::cm;  // pr/84 r2 F3
        pa.m_shower_dedup_start_seg = m_shower_dedup_start_seg;  // pr/84 r3 S1
        pa.m_shower_traj_michel_stem = m_shower_traj_michel_stem;  // pr/74 K6
        pa.m_michel_stem_traj_min_len = m_michel_stem_traj_min_len * units::cm;  // pr/74 K6
        pa.m_michel_stem_traj_max_len = m_michel_stem_traj_max_len * units::cm;  // pr/74 K6
        pa.m_michel_stem_traj_mip_lo = m_michel_stem_traj_mip_lo;  // pr/74 K6 (dimensionless ratio)
        pa.m_michel_stem_traj_max_far_len = m_michel_stem_traj_max_far_len * units::cm;  // pr/74 K6
        pa.m_michel_stem_traj_min_kink_deg = m_michel_stem_traj_min_kink_deg;  // pr/74 K6 (degrees)
        pa.m_shower_long_muon_keep_type = m_shower_long_muon_keep_type;  // doc pr/44
        pa.m_shower_bragg_protect_start_segment = m_shower_bragg_protect_start_segment;  // doc pr/40 round 10
        pa.m_shower_reclass_case_b_dqdx_guard = m_shower_reclass_case_b_dqdx_guard;  // doc pr/93 Cause A
        pa.m_shower_accept_pid_guard = m_shower_accept_pid_guard;  // doc pr/93 Cause B
        pa.m_shower_pid_guard_min_len = m_shower_pid_guard_min_len * units::cm;  // doc pr/93 shared floor
        pa.m_shower_vote_track_pid_counts = m_shower_vote_track_pid_counts;  // doc pr/93 Cause C
        pa.m_shower_cone_absorb_guard = m_shower_cone_absorb_guard;  // doc pr/93 Cause D
        pa.m_shower_detach_track_stem = m_shower_detach_track_stem;  // doc pr/93 r4
        pa.m_shower_ghost_member_drop = m_shower_ghost_member_drop;  // doc pr/99 r2
        pa.m_shower_ghost_overlap_frac = m_shower_ghost_overlap_frac;  // fraction, no conversion (doc pr/99 r2)
        pa.m_shower_ghost_dqdx_ratio = m_shower_ghost_dqdx_ratio;  // ratio, no conversion (doc pr/99 r2)
        pa.m_shower_ghost_min_len = m_shower_ghost_min_len * units::cm;  // cm -> internal (doc pr/99 r2)
        pa.m_shower_hadronic_tag = m_shower_hadronic_tag;  // doc pr/99 r3 A5
        pa.m_shower_hadronic_min_len = m_shower_hadronic_min_len * units::cm;  // cm -> internal (doc pr/99 r3)
        pa.m_shower_hadronic_scan_len = m_shower_hadronic_scan_len * units::cm;  // cm -> internal (doc pr/99 r3)
        pa.m_shower_hadronic_bin = m_shower_hadronic_bin * units::cm;  // cm -> internal (doc pr/99 r3)
        pa.m_shower_hadronic_r_cyl = m_shower_hadronic_r_cyl * units::cm;  // cm -> internal (doc pr/99 r3)
        pa.m_shower_hadronic_r_core = m_shower_hadronic_r_core * units::cm;  // cm -> internal (doc pr/99 r3)
        pa.m_shower_hadronic_growth_max = m_shower_hadronic_growth_max;  // ratio, no conversion (doc pr/99 r3)
        pa.m_shower_hadronic_growth_bragg = m_shower_hadronic_growth_bragg;  // ratio, no conversion (doc pr/99 r3)
        pa.m_shower_hadronic_bragg_ratio = m_shower_hadronic_bragg_ratio;  // ratio, no conversion (doc pr/99 r3)
        pa.m_shower_hadronic_stem_ratio = m_shower_hadronic_stem_ratio;  // MIP units, no conversion (doc pr/99 r3)
        pa.m_kine_count_orphan_tracks = m_kine_count_orphan_tracks;  // doc pr/93 r4
        pa.m_kine_orphan_track_min = m_kine_orphan_track_min * units::cm;  // doc pr/93 r4
        pa.m_shower_pass4_best_owner = m_shower_pass4_best_owner;  // doc pr/117 r1
        pa.m_shower_merge_relax = m_shower_merge_relax;  // doc pr/117 r1
        pa.m_shower_merge_relax_dis = m_shower_merge_relax_dis * units::cm;  // cm -> internal (doc pr/117 r1)
        pa.m_shower_merge_relax_angle = m_shower_merge_relax_angle;  // deg, no conversion (doc pr/117 r1)
        pa.m_shower_merge_relax_min_len = m_shower_merge_relax_min_len * units::cm;  // cm -> internal (doc pr/117 r1)
        pa.m_shower_merge_relax_continuity = m_shower_merge_relax_continuity;  // doc pr/118 r1
        pa.m_shower_merge_relax_cont_frac = m_shower_merge_relax_cont_frac;  // fraction, no conversion (doc pr/118 r1)
        pa.m_shower_merge_relax_cont_gap = m_shower_merge_relax_cont_gap * units::cm;  // cm -> internal (doc pr/118 r1)
        pa.m_shower_merge_relax_cont_qmed = m_shower_merge_relax_cont_qmed;  // charge units, no conversion (doc pr/118 r1)
        pa.m_shower_merge_relax_cont_axis = m_shower_merge_relax_cont_axis;  // deg, no conversion (doc pr/118 r1)
        pa.m_shower_merge_relax_cont_dmax = m_shower_merge_relax_cont_dmax * units::cm;  // cm -> internal (doc pr/118 r1)
        pa.m_shower_merge_relax_cont_t1_gap = m_shower_merge_relax_cont_t1_gap * units::cm;  // cm -> internal (doc pr/118 r1)
        pa.m_shower_merge_relax_cont_t1_fold = m_shower_merge_relax_cont_t1_fold;  // deg, no conversion (doc pr/118 r1)
        pa.m_stem_backfill_back_guard = m_stem_backfill_back_guard;  // doc pr/120 r1
        pa.m_stem_backfill_back_ang = m_stem_backfill_back_ang;  // deg, no conversion (doc pr/120 r1)
        pa.m_shower_ex1_dedup_rehome = m_shower_ex1_dedup_rehome;  // doc pr/121 r1
        pa.m_shower_pass4_prune_detached = m_shower_pass4_prune_detached;  // doc pr/123 r1
        pa.m_shower_pass4_prune_gap = m_shower_pass4_prune_gap * units::cm;  // doc pr/123 r1, cm -> internal
        pa.m_shower_pass4_prune_gap2 = m_shower_pass4_prune_gap2 * units::cm;  // doc pr/124 A, cm -> internal
        pa.m_shower_pass4_prune2_ang = m_shower_pass4_prune2_ang;  // doc pr/124 A, deg
        pa.m_shower_pass4_prune2_mdqdx = m_shower_pass4_prune2_mdqdx;  // doc pr/124 A, x MIP
        pa.m_shower_pass3_cone_guard_len = m_shower_pass3_cone_guard_len * units::cm;  // doc pr/124 C, cm -> internal
        pa.m_shower_samevtx_track_absorb = m_shower_samevtx_track_absorb;  // doc pr/125
        pa.m_shower_samevtx_absorb_gap = m_shower_samevtx_absorb_gap * units::cm;  // doc pr/125, cm -> internal
        pa.m_shower_samevtx_absorb_max_len = m_shower_samevtx_absorb_max_len * units::cm;  // doc pr/125, cm -> internal
        pa.m_shower_samevtx_absorb_min_len = m_shower_samevtx_absorb_min_len * units::cm;  // doc pr/125, cm -> internal
        pa.m_shower_satellite_absorb = m_shower_satellite_absorb;  // doc pr/125
        pa.m_shower_split = m_shower_split;  // doc pr/138 B2; false = no pass
        pa.m_shower_split_max_valley = m_shower_split_max_valley;  // doc pr/138 B2; sec A5.4 knee
        pa.m_shower_split_min_frac = m_shower_split_min_frac;  // doc pr/138 B2; per-seed charge share floor
        pa.m_shower_split_max_parts = m_shower_split_max_parts;  // doc pr/138 B3; 2 = the measured-exact kernel
        pa.m_shower_split_min_charge = m_shower_split_min_charge;  // doc pr/138 B1; candidate charge floor (raw Fit::dQ)
        pa.m_shower_split_min_nseg = m_shower_split_min_nseg;  // doc pr/138 B1; candidate member-count floor
        pa.m_shower_split_bundle_gap = m_shower_split_bundle_gap * units::cm;  // doc pr/138 B3, cm -> internal
        pa.m_shower_split_snap = m_shower_split_snap;  // doc pr/138 B3; k>=3 bundle dominance floor
        pa.m_shower_split_skip_shared = m_shower_split_skip_shared;  // doc pr/139 P1.1; refuse a component holding a segment another shower also owns
        pa.m_shower_split_shed_shared = m_shower_split_shed_shared;  // doc pr/139 sec 15; shed an ENTIRELY co-owned refused component
        pa.m_shower_split_max_seeds = m_shower_split_max_seeds;  // doc pr/139 sec 17; angular-maxima cap (shipped 4)
        pa.m_shower_split_em_type_max_len = m_shower_split_em_type_max_len * units::cm;  // doc pr/139 sec 25; cm -> internal; 0 = off
        pa.m_shower_split_max_impact = m_shower_split_max_impact * units::cm;  // doc pr/139 P1.2; cm -> internal; 0 = no bound
        pa.m_shower_split_em_start = m_shower_split_em_start;  // doc pr/139 P1.3; seed the daughter on its nearest EM-typed member
        pa.m_shower_split_rehome = m_shower_split_rehome;  // doc pr/139 P1.4; offer an orphan daughter to the nearest larger EM shower
        pa.m_shower_split_rehome_gap = m_shower_split_rehome_gap * units::cm;  // doc pr/139 P1.4; cm; max daughter->host 3-D gap, cm -> internal
        pa.m_shower_satellite_absorb_max_mev = m_shower_satellite_absorb_max_mev * units::MeV;  // doc pr/125, MeV -> internal
        pa.m_shower_satellite_absorb_host_mev = m_shower_satellite_absorb_host_mev * units::MeV;  // doc pr/125, MeV -> internal
        pa.m_shower_pass4_track_guard_len = m_shower_pass4_track_guard_len * units::cm;  // doc pr/123 r1, cm -> internal
        pa.m_shower_pass4_prox_guard_len = m_shower_pass4_prox_guard_len * units::cm;  // doc pr/130 item 1b, cm -> internal
        pa.m_shower_pass3_backfill_guard_len = m_shower_pass3_backfill_guard_len * units::cm;  // doc pr/130 item 1b, cm -> internal
        pa.m_stem_backfill_back_dvtx = m_stem_backfill_back_dvtx * units::cm;  // doc pr/130 item B, cm -> internal
        pa.m_shower_pass4_prefilter_v1_escape = m_shower_pass4_prefilter_v1_escape;  // doc pr/136 r2 (bool, no scaling)
        pa.m_shower_pass4_prefilter_v1_max_v2 = m_shower_pass4_prefilter_v1_max_v2;  // doc pr/136 r2, deg (no scaling)
        pa.m_shower_pass4_prefilter_v1_max_dis = m_shower_pass4_prefilter_v1_max_dis * units::cm;  // doc pr/136 r3, cm -> internal
        pa.m_pi0_mass_offset = m_pi0_mass_offset * units::MeV;  // doc pr/132 K1, MeV -> internal
        pa.m_pi0_assoc_angle_deg = m_pi0_assoc_angle_deg;  // doc pr/132 K2, deg (no scaling)
        pa.m_pi0_attached_partner_min = m_pi0_attached_partner_min_mev * units::MeV;  // doc pr/132 K3, MeV -> internal
        pa.m_pi0_nv_max_prongs = m_pi0_nv_max_prongs;  // doc pr/132 K5
        pa.m_pi0_readmit_retyped = m_pi0_readmit_retyped;  // doc pr/132 K7
        pa.m_pi0_admit_type3 = m_pi0_admit_type3;  // doc pr/132 K8
        pa.m_pi0_crumb_assoc_max = m_pi0_crumb_assoc_mev * units::MeV;  // doc pr/132 K9, MeV -> internal
        pa.m_pi0_collinear_merge_deg = m_pi0_collinear_merge_deg;  // doc pr/132 K12, deg (no scaling)
        pa.m_pi0_nv_partner_min = m_pi0_nv_partner_min_mev * units::MeV;  // doc pr/132 K13, MeV -> internal
        pa.m_em_collinear_merge_deg = m_shower_em_collinear_deg;  // doc pr/132 K16, deg (no scaling)
        pa.m_em_collinear_merge_dis = m_shower_em_collinear_dis_cm * units::cm;  // doc pr/132 K16, cm -> internal
        pa.m_em_collinear_merge_min_host = m_shower_em_collinear_host_mev * units::MeV;  // doc pr/132 K16, MeV -> internal
        pa.m_em_backext_perp = m_shower_em_backext_perp_cm * units::cm;  // doc pr/132 K17, cm -> internal
        pa.m_em_backext_len = m_shower_em_backext_len_cm * units::cm;  // doc pr/132 K17, cm -> internal
        pa.m_pi0_am_dis = m_pi0_accept_merge_dis_cm * units::cm;  // doc pr/132 K18, cm -> internal
        pa.m_pi0_bp_miss = m_pi0_bp_vertex_miss_cm * units::cm;  // doc pr/132 K19, cm -> internal
        pa.m_pi0_admit_mu_showers = m_pi0_admit_muon_showers;  // doc pr/133 K20
        pa.m_pi0_nc_sig_angle = m_pi0_nc_sig_angle_deg;  // doc pr/133 K21, deg (unscaled, like K2)
        pa.m_pi0_nc_floor = m_pi0_nc_floor_mev * units::MeV;  // doc pr/133 K21 v2, MeV -> internal
        pa.m_pi0_nc_pf_assoc = m_pi0_nc_pf_assoc_deg;  // doc pr/133 K21 v2.2, deg (unscaled)
        pa.m_pi0_nc_frag_merge = m_pi0_nc_frag_merge;  // doc pr/134 K22
        pa.m_pi0_pf_assoc = m_pi0_pf_assoc_deg;  // doc pr/134 K23, deg (unscaled)
        pa.m_pi0_prefer_main_vertex = m_pi0_prefer_main_vertex;  // doc pr/134 K24
        pa.m_pi0_nv_max_vtx_shift = m_pi0_nv_max_vtx_shift_cm * units::cm;  // doc pr/132 K10, cm -> internal
        pa.m_pi0_nv_mass_window = m_pi0_nv_mass_window_mev * units::MeV;  // doc pr/132 K11, MeV -> internal
        pa.m_kine_count_guard_freed = m_kine_count_guard_freed;  // doc pr/123 r2
        pa.m_kine_guard_freed_impact = m_kine_guard_freed_impact * units::cm;  // doc pr/129
        pa.m_kine_guard_freed_miss_deg = m_kine_guard_freed_miss_deg;  // doc pr/129
        pa.m_kine_count_near_cross_cluster = m_kine_count_near_cross_cluster;  // doc pr/128
        pa.m_kine_near_gap = m_kine_near_gap * units::cm;  // doc pr/128
        pa.m_kine_near_min_len = m_kine_near_min_len * units::cm;  // doc pr/128
        pa.m_kine_near_end_tol = m_kine_near_end_tol * units::cm;  // doc pr/128
        pa.m_kine_near_kink_deg = m_kine_near_kink_deg;  // doc pr/128
        pa.m_kine_near_pointing_impact = m_kine_near_pointing_impact * units::cm;  // doc pr/144
        pa.m_kine_near_pointing_miss_deg = m_kine_near_pointing_miss_deg;  // doc pr/144
        pa.m_kine_count_conn4_near = m_kine_count_conn4_near;  // doc pr/128
        pa.m_kine_conn4_near_gap = m_kine_conn4_near_gap * units::cm;  // doc pr/128
        pa.m_straight_cont_cross_cluster = m_straight_cont_cross_cluster;  // doc pr/93 r4
        pa.m_sccc_bridge_body = m_sccc_bridge_body;  // doc pr/93 r4
        pa.m_sccc_max_gap = m_sccc_max_gap * units::cm;  // doc pr/93 r4
        pa.m_sccc_kink_max = m_sccc_kink_max;  // deg
        pa.m_sccc_gap_aligned = m_sccc_gap_aligned * units::cm;  // doc pr/93 r4
        pa.m_sccc_kink_tight = m_sccc_kink_tight;  // deg
        pa.m_single_muon_proton_chain_veto = m_single_muon_proton_chain_veto;  // doc pr/43 round 2 K1
        pa.m_single_muon_long_muon_claim = m_single_muon_long_muon_claim;  // doc pr/43 round 2 K2
        pa.m_pid_flag_reconcile = m_pid_flag_reconcile;  // doc pr/43 round 2 K3
        pa.m_long_muon_stub_bridge = m_long_muon_stub_bridge;  // doc pr/46
        pa.m_long_muon_stub_bridge_len_cm = m_long_muon_stub_bridge_len;  // doc 84 round 1 (P3)
        pa.m_long_muon_angle_relax_long = m_long_muon_angle_relax_long;  // doc 84 round 1 (P2)
        pa.m_long_muon_angle_relax_deg = m_long_muon_angle_relax_deg;  // doc 84 round 1 (P2)
        pa.m_kine_charge.plane_weights = {m_kine_plane_weights[0], m_kine_plane_weights[1], m_kine_plane_weights[2]};
        pa.m_kine_charge.t0_frame = m_kine_charge_t0_frame;   // doc pdvd/120 sec 9.4 (this stage's own key)
        pa.m_sgp_dv   = m_dv;
        pa.m_sgp_pcts = m_pcts;
        pa.m_recomb_model = m_recomb_model;
        // Process-wide transport, written once per visit exactly as
        // TaggerCheckNeutrino.cxx:3278-3282 / :3300 do (this stage is the
        // only PR visitor of its job).
        PR::g_shower_traj_refresh_flag = m_shower_traj_recheck_parity;
        WireCell::Clus::PR::g_graph_endpoint_policy.tol = m_graph_endpoint_tol * units::cm;
        PR::set_traj_cover_probe(m_traj_cover_probe);
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
        I1("entry_long_muon_nseg", r.entry_long_muon_nseg);
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
    tf->set_parameter("traj_cover_probe", m_traj_cover_probe ? 1.0 : 0.0);   // TaggerCheckNeutrino.cxx:3560
    if (m_excl_t0_frame) tf->set_parameter("excl_t0_frame", 1.0);
    // TaggerCheckNeutrino.cxx:3588-3594 (fit_blob_coverage_defer): the main
    // cluster's partition forms on legacy fits, deweighting restored after.
    const bool cov_defer_active = m_fit_blob_coverage_defer && m_fit_blob_coverage >= 0;
    auto cov_defer_suspend = [&]() { if (cov_defer_active) tf->set_parameter("fit_blob_coverage", -1); };
    auto cov_defer_restore = [&]() { if (cov_defer_active) tf->set_parameter("fit_blob_coverage", m_fit_blob_coverage); };
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
    cov_defer_suspend();
    const bool ok_main = pa.find_proto_vertex(g, *main, *tf, m_dv, true, 2, true, particle_data());
    cov_defer_restore();
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
    // TaggerCheckNeutrino.cxx:3867-3885, the two near-vertex graph passes the
    // neutrino chain runs after its vertex is final and BEFORE clustering_
    // points: the main-vertex graph audit (inert unless main_vertex_graph_
    // audit; may re-seat the vertex POSITION in place, never the pointer) and
    // the conn-3 stitch of disconnected main-cluster components (inert
    // unless conn3_stitch_max > 0).  The kink/junction SNAPS that precede
    // them in TaggerCheckNeutrino are deliberately not run: the entry is not
    // a topological choice.  doc pdvd/120 sec 9.
    if (pa.main_vertex_graph_audit(g, *main, final_main_vertex, *tf, m_dv)) {
        SPDLOG_LOGGER_DEBUG(s_log, "{}CheckBeamParticle: main_vertex_graph_audit changed the graph around the entry", m_evt_tag);
    }
    if (pa.stitch_disconnected_main_cluster(g, *main, final_main_vertex, *tf, m_dv)) {
        SPDLOG_LOGGER_DEBUG(s_log, "{}CheckBeamParticle: conn3 stitch bridged a disconnected main-cluster component", m_evt_tag);
    }
    pa.clustering_points(g, *main, m_dv);
    pa.reassociate_cluster_orphans(g, *main, m_dv);
    for (auto* cluster : companions) pa.reassociate_cluster_orphans(g, *cluster, m_dv);
    // examine_direction runs last and has the final word on segment
    // orientations relative to the main vertex -- here the entry.
    pa.examine_direction(g, final_main_vertex, final_main_vertex, vertices_in_long_muon, segments_in_long_muon,
                         particle_data(), m_recomb_model, true);
    // doc pdvd/120 sec 9.3 -- the entry-rooted "long muon".  examine_direction's
    // long-muon search (NeutrinoVertexFinder.cxx:1896-1955) starts at the main
    // vertex, walks find_cont_muon_segment THROUGH every junction it can
    // continue across, and shower_clustering_with_nv_in_main_cluster then seeds
    // ONE pseudo-shower on the first chain segment whose flood-fill
    // (Shower::complete_structure_with_start_segment, no stopping rule for a
    // type-13 shower) absorbs every segment connected beyond it.  For a
    // neutrino that is the exiting muon and nothing is lost; for a beam
    // particle rooted at its ENTRY it is the whole event: run 39305 evt 157312
    // rendered as "mu- 532 MeV, 11 segments" from the entry to the far end of
    // the main cluster, with the interaction vertex 77 cm in (two protons, two
    // EM showers and the outgoing muon that the neutrino chain resolves from
    // that vertex) swallowed into the one node.  Releasing the chain here --
    // the PID 13 stamps and the cleared shower flags it left on the beam
    // segments stay -- lets the beam track be ordinary track segments: the
    // shower seeder's BFS descends through them, the daughters at the
    // interaction vertex seed their own showers or stay tracks, and the Bee
    // particle flow becomes beam segment(s) -> daughters.  C++ default false
    // = release; entry_long_muon_absorb=true restores the neutrino-chain
    // behaviour for comparison.
    rec.entry_long_muon_nseg = static_cast<int>(segments_in_long_muon.size());
    if (!m_entry_long_muon_absorb && !segments_in_long_muon.empty()) {
        SPDLOG_LOGGER_INFO(s_log, "{}CheckBeamParticle: released the entry-rooted long-muon chain ({} segment(s), {} vertex(es)) to the track BFS",
                           m_evt_tag, segments_in_long_muon.size(), vertices_in_long_muon.size());
        segments_in_long_muon.clear();
        vertices_in_long_muon.clear();
    }
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
    // TaggerCheckNeutrino.cxx:3961-3973 (swap_orphan_dup_audit): one
    // duplicate-corridor audit per companion, cluster-id order, before the
    // shower maps are built.  Knob off => loop skipped.
    if (m_swap_orphan_dup_audit) {
        std::vector<Cluster*> audit_clusters(companions.begin(), companions.end());
        std::sort(audit_clusters.begin(), audit_clusters.end(),
                  [](Cluster* a, Cluster* b) { return a->get_cluster_id() < b->get_cluster_id(); });
        for (Cluster* oc : audit_clusters) {
            if (!oc || oc == main) continue;
            pa.orphan_dup_audit(g, *oc, *tf, m_dv);
        }
    }
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
    // TaggerCheckNeutrino.cxx:4014-4032 (long_muon_range_empty_chain_fallback):
    // a shower retyped |13| after its kinematics pass keeps kine_range 0;
    // clear the flag on exactly those and rerun the kinematics pass.
    if (m_long_muon_range_empty_chain_fallback) {
        int n_stale = 0;
        for (auto& shower : showers) {
            if (!shower || !shower->get_flag_kinematics()) continue;
            if (std::abs(shower->get_particle_type()) != 13) continue;
            if (shower->get_kine_range() != 0) continue;
            shower->set_flag_kinematics(false);
            ++n_stale;
        }
        if (n_stale > 0) {
            pa.calculate_shower_kinematics(showers, vertices_in_long_muon, segments_in_long_muon, g, *tf, m_dv,
                                           particle_data(), m_recomb_model);
            SPDLOG_LOGGER_DEBUG(s_log, "{}CheckBeamParticle: long_muon_range_empty_chain_fallback recomputed {} retyped muon shower(s)",
                                m_evt_tag, n_stale);
        }
    }
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
