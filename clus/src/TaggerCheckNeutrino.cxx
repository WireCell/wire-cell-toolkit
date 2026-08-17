#include "WireCellClus/TaggerCheckNeutrino.h"
#include "WireCellClus/NeutrinoPatternBase.h" // pattern recognition ...
#include "WireCellClus/PatternDebugIO.h"      // debug dump/load

#include "WireCellUtil/Persist.h"
#include <chrono>
#include <set>

#include <cstdlib>

class TaggerCheckNeutrino;
WIRECELL_FACTORY(TaggerCheckNeutrino, TaggerCheckNeutrino,
                 WireCell::INamed, WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::Clus::PR;


struct edge_base_t {
    typedef boost::edge_property_tag kind;
};

// Temporary determinism-debug helper, enabled with WCT_DET_DEBUG=1.
// Prints a content checksum of the PR graph so two runs can be diffed to find
// the first PR phase that diverges run-to-run.
namespace {
    inline uint64_t detg_fnv(uint64_t h, const void* p, size_t n) {
        const unsigned char* b = static_cast<const unsigned char*>(p);
        for (size_t i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ULL; }
        return h;
    }
    void detg_dump(const char* tag, WireCell::Clus::PR::Graph& g) {
        static const bool on = (std::getenv("WCT_DET_DEBUG") != nullptr);
        if (!on) return;
        uint64_t h = 1469598103934665603ULL;
        size_t nvtx = 0, nseg = 0;
        for (auto nd : WireCell::Clus::PR::ordered_nodes(g)) {
            auto vtx = g[nd].vertex;
            if (!vtx) continue;
            ++nvtx;
            uint64_t idx = vtx->get_graph_index();
            h = detg_fnv(h, &idx, sizeof(idx));
            const auto& p = vtx->wcpt().point;
            double xyz[3] = {p.x(), p.y(), p.z()};
            h = detg_fnv(h, xyz, sizeof(xyz));
        }
        for (auto ed : WireCell::Clus::PR::ordered_edges(g)) {
            auto seg = g[ed].segment;
            if (!seg) continue;
            ++nseg;
            uint64_t v[3] = {seg->get_graph_index(), seg->wcpts().size(), seg->fits().size()};
            h = detg_fnv(h, v, sizeof(v));
            int cid = seg->cluster() ? seg->cluster()->get_cluster_id() : -1;
            h = detg_fnv(h, &cid, sizeof(cid));
            for (const auto& w : seg->wcpts()) {
                double xyz[3] = {w.point.x(), w.point.y(), w.point.z()};
                h = detg_fnv(h, xyz, sizeof(xyz));
            }
        }
        fprintf(stderr, "WCT_DETG %-28s nvtx=%zu nseg=%zu h=%016lx\n", tag, nvtx, nseg, (unsigned long)h);
    }
}

void TaggerCheckNeutrino::configure(const WireCell::Configuration& config)
{
    m_grouping_name = get(config, "grouping_name", m_grouping_name);
    m_trackfitting_config_file = get(config, "trackfitting_config_file", m_trackfitting_config_file);
    m_perf = get(config, "perf", m_perf);
    m_dir_weak_use_score = get(config, "dir_weak_use_score", m_dir_weak_use_score);
    m_mip_dqdx            = get(config, "mip_dqdx",             m_mip_dqdx);             // e/cm
    m_mip_dqdx_median     = get(config, "mip_dqdx_median",      m_mip_dqdx_median);      // e/cm
    m_proton_dir_vote      = get(config, "proton_dir_vote",      m_proton_dir_vote);
    m_proton_dir_score_max = get(config, "proton_dir_score_max", m_proton_dir_score_max);
    m_proton_dir_asym_min  = get(config, "proton_dir_asym_min",  m_proton_dir_asym_min);
    m_endpoint_trim_retry  = get(config, "endpoint_trim_retry",  m_endpoint_trim_retry);
    m_fit_vertex_min_seg_length = get(config, "fit_vertex_min_seg_length", m_fit_vertex_min_seg_length);  // cm
    // doc sbnd_xin/docs/pr/51 round 7: robust vertex fit
    m_mvfit_robust      = get(config, "mvfit_robust",      m_mvfit_robust);
    m_mvfit_main_only   = get(config, "mvfit_main_only",   m_mvfit_main_only);
    m_mvfit_min_len     = get(config, "mvfit_min_len",     m_mvfit_min_len);      // cm
    m_mvfit_rin_margin  = get(config, "mvfit_rin_margin",  m_mvfit_rin_margin);   // cm
    m_mvfit_rout_frac   = get(config, "mvfit_rout_frac",   m_mvfit_rout_frac);
    m_mvfit_rout_min    = get(config, "mvfit_rout_min",    m_mvfit_rout_min);     // cm
    m_mvfit_rout_max    = get(config, "mvfit_rout_max",    m_mvfit_rout_max);     // cm
    m_mvfit_angle       = get(config, "mvfit_angle",       m_mvfit_angle);        // deg
    m_mvfit_min_pts     = get(config, "mvfit_min_pts",     m_mvfit_min_pts);
    m_mvfit_min_aniso   = get(config, "mvfit_min_aniso",   m_mvfit_min_aniso);
    m_mvfit_prior_range = get(config, "mvfit_prior_range", m_mvfit_prior_range);  // cm
    m_cathode_x          = get(config, "cathode_x",          m_cathode_x);           // cm
    m_cathode_kink_xcut  = get(config, "cathode_kink_xcut",  m_cathode_kink_xcut);   // cm
    m_cathode_wide_kink_angle    = get(config, "cathode_wide_kink_angle",    m_cathode_wide_kink_angle);    // deg
    m_cathode_wide_kink_skirt    = get(config, "cathode_wide_kink_skirt",    m_cathode_wide_kink_skirt);    // cm
    m_cathode_wide_kink_baseline = get(config, "cathode_wide_kink_baseline", m_cathode_wide_kink_baseline); // cm
    // doc sbnd_xin/docs/pr/48 -- back-to-back track fixes.
    m_two_end_break       = get(config, "two_end_break",       m_two_end_break);
    m_teb_min_len         = get(config, "teb_min_len",         m_teb_min_len);         // cm
    m_teb_min_arm         = get(config, "teb_min_arm",         m_teb_min_arm);         // cm
    m_teb_min_arm_pts     = get(config, "teb_min_arm_pts",     m_teb_min_arm_pts);
    m_teb_stub_max        = get(config, "teb_stub_max",        m_teb_stub_max);        // cm
    m_teb_accept_range    = get(config, "teb_accept_range",    m_teb_accept_range);    // cm
    m_teb_rise_r1         = get(config, "teb_rise_r1",         m_teb_rise_r1);
    m_teb_rise_r2         = get(config, "teb_rise_r2",         m_teb_rise_r2);
    m_teb_abs_end_min     = get(config, "teb_abs_end_min",     m_teb_abs_end_min);
    m_teb_dip_floor       = get(config, "teb_dip_floor",       m_teb_dip_floor);
    m_teb_score_cap_r1    = get(config, "teb_score_cap_r1",    m_teb_score_cap_r1);
    m_teb_score_cap_r2    = get(config, "teb_score_cap_r2",    m_teb_score_cap_r2);
    m_teb_turn_angle      = get(config, "teb_turn_angle",      m_teb_turn_angle);      // deg
    m_teb_turn_baseline   = get(config, "teb_turn_baseline",   m_teb_turn_baseline);   // cm
    m_teb_turn_skirt      = get(config, "teb_turn_skirt",      m_teb_turn_skirt);      // cm
    m_teb_turn_min_arm_frac = get(config, "teb_turn_min_arm_frac", m_teb_turn_min_arm_frac); // frac of baseline; doc pr/90 round 2
    m_teb_second_max      = get(config, "teb_second_max",      m_teb_second_max);      // cm; doc pr/90 round 2
    m_teb_chain_topology  = get(config, "teb_chain_topology",  m_teb_chain_topology);  // doc pr/90 round 4 (D1)
    m_teb_r3_turn         = get(config, "teb_r3_turn",         m_teb_r3_turn);         // deg; doc pr/90 round 4 (D3)
    m_teb_r3_hot          = get(config, "teb_r3_hot",          m_teb_r3_hot);          // x mip median; doc pr/90 round 4 (D3)
    m_teb_bragg_veto_turn = get(config, "teb_bragg_veto_turn", m_teb_bragg_veto_turn); // deg; doc pr/90 round 4 (D4)
    m_kink_walk_dqdx_stop = get(config, "kink_walk_dqdx_stop", m_kink_walk_dqdx_stop);
    m_kink_break_protect  = get(config, "kink_break_protect",  m_kink_break_protect);
    m_kink_dqdx_hot_ratio = get(config, "kink_dqdx_hot_ratio", m_kink_dqdx_hot_ratio);
    // doc sbnd_xin/docs/pr/50 -- main-vertex kink-consistency snap.
    m_vertex_kink_snap = get(config, "vertex_kink_snap", m_vertex_kink_snap);
    m_vks_radius       = get(config, "vks_radius",       m_vks_radius);       // cm
    m_vks_min_dis      = get(config, "vks_min_dis",      m_vks_min_dis);      // cm
    m_vks_angle        = get(config, "vks_angle",        m_vks_angle);        // deg
    m_vks_margin       = get(config, "vks_margin",       m_vks_margin);       // deg
    m_vks_collinear    = get(config, "vks_collinear",    m_vks_collinear);    // deg
    m_vks_skirt        = get(config, "vks_skirt",        m_vks_skirt);        // cm
    m_vks_baseline     = get(config, "vks_baseline",     m_vks_baseline);     // cm
    m_vks_min_arm      = get(config, "vks_min_arm",      m_vks_min_arm);      // cm
    m_vks_fit_miss     = get(config, "vks_fit_miss",     m_vks_fit_miss);     // cm
    m_vks_hot_ratio    = get(config, "vks_hot_ratio",    m_vks_hot_ratio);
    m_vks_carry_prong  = get(config, "vks_carry_prong",  m_vks_carry_prong);  // cm (doc pr/85)
    // doc sbnd_xin/docs/pr/51 -- main-vertex graph audit.
    m_main_vertex_graph_audit = get(config, "main_vertex_graph_audit", m_main_vertex_graph_audit);
    m_mvga_radius       = get(config, "mvga_radius",       m_mvga_radius);       // cm
    m_mvga_dup_tol      = get(config, "mvga_dup_tol",      m_mvga_dup_tol);      // cm
    m_mvga_dup_frac     = get(config, "mvga_dup_frac",     m_mvga_dup_frac);
    m_mvga_dup_angle    = get(config, "mvga_dup_angle",    m_mvga_dup_angle);    // deg
    m_mvga_bridge_mip   = get(config, "mvga_bridge_mip",   m_mvga_bridge_mip);
    m_mvga_reconnect    = get(config, "mvga_reconnect",    m_mvga_reconnect);    // cm
    m_mvga_stub         = get(config, "mvga_stub",         m_mvga_stub);         // cm
    m_mvga_stub_pts     = get(config, "mvga_stub_pts",     m_mvga_stub_pts);
    m_mvga_reseat_angle = get(config, "mvga_reseat_angle", m_mvga_reseat_angle); // deg
    m_mvga_satellite    = get(config, "mvga_satellite",    m_mvga_satellite);    // cm
    m_mvga_interposed   = get(config, "mvga_interposed",   m_mvga_interposed);   // doc pr/85
    m_mvga_interposed_angle = get(config, "mvga_interposed_angle", m_mvga_interposed_angle); // deg
    m_mvga_interposed_len = get(config, "mvga_interposed_len", m_mvga_interposed_len); // cm; doc pr/86
    m_mvga_sat_dup_frac = get(config, "mvga_sat_dup_frac", m_mvga_sat_dup_frac); // fraction; doc pr/86
    m_mvga_interposed_deg1 = get(config, "mvga_interposed_deg1", m_mvga_interposed_deg1); // doc pr/86
    m_mvga_splice_straighten = get(config, "mvga_splice_straighten", m_mvga_splice_straighten); // cm; doc pr/86 round 2
    m_mvga_approach_collapse = get(config, "mvga_approach_collapse", m_mvga_approach_collapse); // cm; doc pr/86 round 2
    m_mvga_straighten_radius = get(config, "mvga_straighten_radius", m_mvga_straighten_radius); // cm; doc pr/86 round 2
    m_shower_topo_demote_len = get(config, "shower_topo_demote_len", m_shower_topo_demote_len);  // cm
    // doc sbnd_xin/docs/pr/49 -- cross-cluster projection-ghost deweighting
    // in the trajectory fit's 2D charge association (18255-57441): live
    // cells outside the fitted cluster's own blob coverage that sit inside a
    // 3D-distant foreign cluster's keep their measurement at reduced weight
    // (deweight, not drop -- dead-channel single-view charge must stay
    // usable).  < 0 = off (legacy); >= 0 = on with value = wire/slice
    // tolerance cells.  Companion numerics (ghost_dis 15 cm, weight 0.1)
    // ride the TrackFitting::Parameters C++ defaults.
    m_fit_blob_coverage = get(config, "fit_blob_coverage", m_fit_blob_coverage);
    // doc sbnd_xin/docs/pr/50: suspend the deweighting during the
    // partition-forming find_proto_vertex stage (false = pr/49 behavior).
    m_fit_blob_coverage_defer = get(config, "fit_blob_coverage_defer", m_fit_blob_coverage_defer);
    // doc sbnd_xin/docs/pr/30 §11 port-fidelity knobs.
    m_fit_exclusion            = get(config, "fit_exclusion",            m_fit_exclusion);
    m_graph_endpoint_strict    = get(config, "graph_endpoint_strict",    m_graph_endpoint_strict);
    m_graph_endpoint_tol       = get(config, "graph_endpoint_tol",       m_graph_endpoint_tol);       // cm
    m_oov_prototype_parity     = get(config, "oov_prototype_parity",     m_oov_prototype_parity);
    m_first_seg_local_pca      = get(config, "first_seg_local_pca",      m_first_seg_local_pca);
    m_other_seg_relaxed_accept = get(config, "other_seg_relaxed_accept", m_other_seg_relaxed_accept);
    // doc sbnd_xin/docs/pr/45 -- find_other_segments empty-2D-tree sentinel guard.
    m_other_seg_empty_2d_guard = get(config, "other_seg_empty_2d_guard", m_other_seg_empty_2d_guard);
    // doc sbnd_xin/docs/pr/54 -- keep well-supported isolated residual segments.
    m_other_seg_keep_isolated            = get(config, "other_seg_keep_isolated",            m_other_seg_keep_isolated);
    m_other_seg_keep_isolated_min_points = get(config, "other_seg_keep_isolated_min_points", m_other_seg_keep_isolated_min_points);
    m_other_seg_keep_isolated_min_length = get(config, "other_seg_keep_isolated_min_length", m_other_seg_keep_isolated_min_length); // cm
    // doc sbnd_xin/docs/pr/67 round 3 -- isochronous-snap size gate.
    m_iso_snap_min_dir_mag = get(config, "iso_snap_min_dir_mag", m_iso_snap_min_dir_mag); // cm
    // doc sbnd_xin/docs/pr/59 round 2 -- per-cluster orphaned-associate_points rescue.
    m_assoc_full_recluster = get(config, "assoc_full_recluster", m_assoc_full_recluster);
    // doc sbnd_xin/docs/pr/64 round 7 -- reassign same-cluster association orphans.
    m_assoc_reassign_orphans = get(config, "assoc_reassign_orphans", m_assoc_reassign_orphans);
    // doc sbnd_xin/docs/pr/64 round 8 -- clear a merge survivor's associate_points.
    m_assoc_clear_on_merge = get(config, "assoc_clear_on_merge", m_assoc_clear_on_merge);
    // doc sbnd_xin/docs/pr/31 §11 port-fidelity knob (F2, was P2).
    m_shower_topo_proto_dir    = get(config, "shower_topo_proto_dir",    m_shower_topo_proto_dir);
    // doc sbnd_xin/docs/pr/32 §11 port-fidelity knobs (F1-F4).
    m_vertex_dir_use_fit_point       = get(config, "vertex_dir_use_fit_point",       m_vertex_dir_use_fit_point);
    m_shower_traj_recheck_parity     = get(config, "shower_traj_recheck_parity",     m_shower_traj_recheck_parity);
    m_main_vertex_require_descriptor = get(config, "main_vertex_require_descriptor", m_main_vertex_require_descriptor);
    m_main_vertex_candidate_flag     = get(config, "main_vertex_candidate_flag",     m_main_vertex_candidate_flag);
    // doc sbnd_xin/docs/pr/31 §12 port-fidelity knobs (the §10.12 round).
    m_cont_muon_dir3_30cm            = get(config, "cont_muon_dir3_30cm",            m_cont_muon_dir3_30cm);
    m_track_comp_empty_abstain       = get(config, "track_comp_empty_abstain",       m_track_comp_empty_abstain);
    m_shower_topo_reset              = get(config, "shower_topo_reset",              m_shower_topo_reset);
    m_reclass_preserve_4mom          = get(config, "reclass_preserve_4mom",          m_reclass_preserve_4mom);
    m_reclass_never_computed_ke_floor = get(config, "reclass_never_computed_ke_floor", m_reclass_never_computed_ke_floor);
    m_dir_track_median_local         = get(config, "dir_track_median_local",         m_dir_track_median_local);
    m_examine_showers_vertex_by_index = get(config, "examine_showers_vertex_by_index", m_examine_showers_vertex_by_index);
    m_iso_endpoint               = get(config, "iso_endpoint",               m_iso_endpoint);
    m_iso_endpoint_min_length    = get(config, "iso_endpoint_min_length",    m_iso_endpoint_min_length);     // cm
    m_iso_endpoint_max_xext      = get(config, "iso_endpoint_max_xext",      m_iso_endpoint_max_xext);       // cm
    m_iso_endpoint_xext_frac     = get(config, "iso_endpoint_xext_frac",     m_iso_endpoint_xext_frac);
    m_iso_endpoint_xext_quantile = get(config, "iso_endpoint_xext_quantile", m_iso_endpoint_xext_quantile);
    m_iso_endpoint_tube_radius   = get(config, "iso_endpoint_tube_radius",   m_iso_endpoint_tube_radius);     // cm
    m_iso_endpoint_min_aspect    = get(config, "iso_endpoint_min_aspect",    m_iso_endpoint_min_aspect);
    // doc sbnd_xin/docs/pr/67 (log-only probe + round-budget counterfactual).
    m_traj_cover_probe           = get(config, "traj_cover_probe",           m_traj_cover_probe);
    m_pr_find_other_rounds       = get(config, "pr_find_other_rounds",       m_pr_find_other_rounds);
    // doc sbnd_xin/docs/pr/24 §18 (round 5).
    m_v3_extension_guard         = get(config, "v3_extension_guard",         m_v3_extension_guard);
    m_v3_extension_min_gain      = get(config, "v3_extension_min_gain",      m_v3_extension_min_gain);        // cm
    // doc sbnd_xin/docs/pr/72 round 2 -- examine_structure_3 stub guard.
    m_es3_stub_guard             = get(config, "es3_stub_guard",             m_es3_stub_guard);
    m_es3sg_stub_max             = get(config, "es3sg_stub_max",             m_es3sg_stub_max);              // cm
    m_es3sg_len_ratio            = get(config, "es3sg_len_ratio",            m_es3sg_len_ratio);
    m_es3sg_ang3_min             = get(config, "es3sg_ang3_min",             m_es3sg_ang3_min);              // deg
    m_es3sg_ang_ratio            = get(config, "es3sg_ang_ratio",            m_es3sg_ang_ratio);
    m_es3sg_require_terminal     = get(config, "es3sg_require_terminal",     m_es3sg_require_terminal);
    // Detector-extent literals (docs/pr/2 sec. 2e(iv)), all cm.
    m_cosmic_y_top_main    = get(config, "cosmic_y_top_main",    m_cosmic_y_top_main);
    m_cosmic_y_top_strict  = get(config, "cosmic_y_top_strict",  m_cosmic_y_top_strict);
    m_cosmic_y_top_loose   = get(config, "cosmic_y_top_loose",   m_cosmic_y_top_loose);
    m_cosmic_y_small_piece = get(config, "cosmic_y_small_piece", m_cosmic_y_small_piece);
    m_vertex_z_prior_scale = get(config, "vertex_z_prior_scale", m_vertex_z_prior_scale);
    if (m_vertex_z_prior_scale <= 0) {
        SPDLOG_LOGGER_WARN(log, "TaggerCheckNeutrino: vertex_z_prior_scale must be > 0; keeping 200 cm");
        m_vertex_z_prior_scale = 200;
    }
    // SSM beam-line reference directions, {x,y,z} in the detector frame.
    // Anything malformed keeps the default and warns rather than falling
    // through to (0,0,0): a zero reference makes every safe_acos(dot) return
    // exactly pi/2 and silently pins all 8 ssm_*_angle_{target,absorber}
    // features to 1.5708.
    auto get_dir3 = [&](const char* key, std::vector<double> def) {
        if (!config.isMember(key)) return def;
        const auto& jv = config[key];
        if (!jv.isArray() || jv.size() != 3u) {
            SPDLOG_LOGGER_WARN(log, "TaggerCheckNeutrino: {} must be a 3-element array; keeping default", key);
            return def;
        }
        std::vector<double> v{jv[0].asDouble(), jv[1].asDouble(), jv[2].asDouble()};
        if (std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) < 1e-9) {
            SPDLOG_LOGGER_WARN(log, "TaggerCheckNeutrino: {} is the zero vector; keeping default", key);
            return def;
        }
        return v;
    };
    m_ssm_target_dir   = get_dir3("ssm_target_dir",   m_ssm_target_dir);
    m_ssm_absorber_dir = get_dir3("ssm_absorber_dir", m_ssm_absorber_dir);
    // Charge -> energy calibration constants (docs/pr/2 sec. 2e(iii)).
    m_kine_fudge_factor        = get(config, "kine_fudge_factor",        m_kine_fudge_factor);
    m_kine_recom_factor        = get(config, "kine_recom_factor",        m_kine_recom_factor);
    m_kine_shower_fudge_factor = get(config, "kine_shower_fudge_factor", m_kine_shower_fudge_factor);
    m_kine_shower_recom_factor = get(config, "kine_shower_recom_factor", m_kine_shower_recom_factor);
    m_kine_proton_recom_factor = get(config, "kine_proton_recom_factor", m_kine_proton_recom_factor);
    m_kine_plane_asym_switch   = get(config, "kine_plane_asym_switch",   m_kine_plane_asym_switch);
    m_kine_w_value             = get(config, "kine_w_value",             m_kine_w_value);
    m_kine_shower_pdg_live     = get(config, "kine_shower_pdg_live",     m_kine_shower_pdg_live);
    // Per-plane weights {U,V,W}.  Unlike the SSM directions a single zero entry
    // is legitimate ("ignore this plane"); only a malformed array or an all-zero
    // triple is rejected -- a zero sum would divide the plane average by zero.
    if (config.isMember("kine_plane_weights")) {
        const auto& jv = config["kine_plane_weights"];
        if (!jv.isArray() || jv.size() != 3u) {
            SPDLOG_LOGGER_WARN(log, "TaggerCheckNeutrino: kine_plane_weights must be a 3-element array; keeping default");
        }
        else {
            std::vector<double> v{jv[0].asDouble(), jv[1].asDouble(), jv[2].asDouble()};
            if (v[0] + v[1] + v[2] <= 0) {
                SPDLOG_LOGGER_WARN(log, "TaggerCheckNeutrino: kine_plane_weights must sum to > 0; keeping default");
            }
            else {
                m_kine_plane_weights = v;
            }
        }
    }
    // Muon median-dQ/dx-vs-length envelope {c0, c1, pivot_cm, power}
    // (docs/pr/2 sec. 2e(iv)).  A non-positive pivot or power would make the
    // envelope constant or inverted for every length; reject and keep default.
    if (config.isMember("muon_dqdx_curve")) {
        const auto& jv = config["muon_dqdx_curve"];
        if (!jv.isArray() || jv.size() != 4u) {
            SPDLOG_LOGGER_WARN(log, "TaggerCheckNeutrino: muon_dqdx_curve must be a 4-element array; keeping default");
        }
        else {
            std::vector<double> v{jv[0].asDouble(), jv[1].asDouble(), jv[2].asDouble(), jv[3].asDouble()};
            if (v[2] <= 0 || v[3] <= 0) {
                SPDLOG_LOGGER_WARN(log, "TaggerCheckNeutrino: muon_dqdx_curve pivot and power must be > 0; keeping default");
            }
            else {
                m_muon_dqdx_curve = v;
            }
        }
    }
    // Single-photon stem dE/dx conversion (docs/pr/2 sec. 2e(i)).
    m_sp_dedx_use_recomb_model = get(config, "sp_dedx_use_recomb_model", m_sp_dedx_use_recomb_model);
    m_sp_mean_dedx_cut         = get(config, "sp_mean_dedx_cut",         m_sp_mean_dedx_cut);
    auto dl_weights_raw = get(config, "dl_weights", m_dl_weights);
    if (!dl_weights_raw.empty()) {
        m_dl_weights = Persist::resolve(dl_weights_raw);
        if (m_dl_weights.empty()) {
            SPDLOG_LOGGER_WARN(log, "TaggerCheckNeutrino: dl_weights path not found: {}", dl_weights_raw);
            // doc sbnd_xin/docs/pr/75: a failed resolve leaves m_dl_weights
            // EMPTY, so at the call site doc pr/52 route 3 is indistinguishable
            // from route 1 ("DL never configured").  Remember which it was --
            // it is the only way the scoreboard can tell them apart.
            m_dl_weights_missing = true;
        }
    }
    m_dl_vtx_cut              = get(config, "dl_vtx_cut",              m_dl_vtx_cut);
    m_dQdx_scale              = get(config, "dQdx_scale",              m_dQdx_scale);
    m_dQdx_offset             = get(config, "dQdx_offset",             m_dQdx_offset);
    m_dl_vtx_rerank           = get(config, "dl_vtx_rerank",           m_dl_vtx_rerank);
    m_dl_vtx_top_k            = get(config, "dl_vtx_top_k",            m_dl_vtx_top_k);
    m_dl_vtx_min_accept_score = get(config, "dl_vtx_min_accept_score", m_dl_vtx_min_accept_score);
    m_dl_vtx_score_scale      = get(config, "dl_vtx_score_scale",      m_dl_vtx_score_scale);
    // doc sbnd_xin/docs/pr/51 (18255-506746): cross-cluster DL swap guard.
    m_dl_vtx_swap_guard       = get(config, "dl_vtx_swap_guard",       m_dl_vtx_swap_guard);
    // doc sbnd_xin/docs/pr/89 Arm C (C2): rule-1 topology term weight/center.
    m_dl_vtx_topo_weight      = get(config, "dl_vtx_topo_weight",      m_dl_vtx_topo_weight);
    m_dl_vtx_topo_center      = get(config, "dl_vtx_topo_center",      m_dl_vtx_topo_center);
    // doc sbnd_xin/docs/pr/51 round 3: traditional-path swap propagation.
    m_main_vertex_swap_apply  = get(config, "main_vertex_swap_apply",  m_main_vertex_swap_apply);
    // doc sbnd_xin/docs/pr/51 round 4: diagnostic-only rough-path probe.
    m_rough_path_probe        = get(config, "rough_path_probe",        m_rough_path_probe);
    // doc sbnd_xin/docs/pr/51 round 5: steiner gap penalty (0 = legacy).
    m_steiner_gap_penalty     = get(config, "steiner_gap_penalty",     m_steiner_gap_penalty);
    m_sgp_dead_alpha          = get(config, "sgp_dead_alpha",          m_sgp_dead_alpha);
    m_sgp_min_edge            = get(config, "sgp_min_edge",            m_sgp_min_edge);      // cm
    m_sgp_sample_step         = get(config, "sgp_sample_step",         m_sgp_sample_step);   // cm
    m_sgp_point_radius        = get(config, "sgp_point_radius",        m_sgp_point_radius);  // cm
    m_sgp_edge_probe          = get(config, "sgp_edge_probe",          m_sgp_edge_probe);
    m_vertex_scoreboard       = get(config, "vertex_scoreboard",       m_vertex_scoreboard);
    // doc sbnd_xin/docs/pr/79 §10: live-feature harvest (requires the board).
    m_dl_vtx_harvest          = get(config, "dl_vtx_harvest",          m_dl_vtx_harvest);
    if (m_dl_vtx_harvest && !m_vertex_scoreboard) {
        SPDLOG_LOGGER_WARN(log, "TaggerCheckNeutrino: dl_vtx_harvest requires vertex_scoreboard; harvest inert");
    }
    // doc sbnd_xin/docs/pr/51 round 6: weak-charge deficit term (0 = legacy round-5 flavor).
    m_sgp_weak_scale          = get(config, "sgp_weak_scale",          m_sgp_weak_scale);
    m_sgp_weak_qref           = get(config, "sgp_weak_qref",           m_sgp_weak_qref);     // charge units
    // doc sbnd_xin/docs/pr/73 round 2 F3a: route excursion cap (cm); < 0 = off.
    m_sgp_max_sep             = get(config, "sgp_max_sep",             m_sgp_max_sep);       // cm
    // doc sbnd_xin/docs/pr/83: oriented break_segment splits; false = legacy.
    m_break_seg_orient        = get(config, "break_seg_orient",        m_break_seg_orient);
    m_beam_window_low         = get(config, "beam_window_low",         m_beam_window_low);
    m_beam_window_high        = get(config, "beam_window_high",        m_beam_window_high);
    m_nu_skip_cosmic          = get(config, "nu_skip_cosmic",          m_nu_skip_cosmic);
    m_nu_skip_cosmic_bundle   = get(config, "nu_skip_cosmic_bundle",   m_nu_skip_cosmic_bundle);
    m_nu_skip_cosmic_bundle_min_length = get(config, "nu_skip_cosmic_bundle_min_length", m_nu_skip_cosmic_bundle_min_length);  // cm
    m_skip_cosmic_companions  = get(config, "skip_cosmic_companions",  m_skip_cosmic_companions);
    m_cosmic_companion_min_length = get(config, "cosmic_companion_min_length", m_cosmic_companion_min_length);  // cm
    m_sp_photon_flag          = get(config, "sp_photon_flag",          m_sp_photon_flag);

    // ---- doc sbnd_xin/docs/pr/36 §10 tagger-stage knobs ---------------------
    // F1: fiducial volume for the match_isFC recompute.  Same knob shape and
    // rationale as TaggerCheckSTM.cxx:69-116: "fiducial" absent => the
    // historical FiducialUtils fallback (sensitive-volume union, no margin),
    // which is NOT the volume TaggerCheckTGM/FC/STM test against; setting
    // "fiducial" + "fv_tolerance" to the values those taggers get restores
    // one containment definition across the stage.  NeedFiducial::configure
    // runs ONLY when the key is present: it otherwise falls back to the
    // type-only name "DetectorVolumes", which is not an instantiated
    // component in the SBND PR config (it is named dv-apa0-1) and would
    // throw.  Same guard as TaggerCheckSTM.cxx:101 / TaggerCheckFC.cxx:82.
    m_use_fiducial = !config["fiducial"].isNull();
    if (m_use_fiducial) NeedFiducial::configure(config);
    // fv_tolerance: [x_lo,x_hi,y_lo,y_hi,z_lo,z_hi] margins, negative =
    // inset (FiducialUtils::inside_fiducial_volume convention; 3 or 1
    // element also accepted).
    {
        auto tol = config["fv_tolerance"];
        if (!tol.isNull() && tol.isArray()) {
            m_fv_tolerance.clear();
            for (const auto& t : tol) m_fv_tolerance.push_back(t.asDouble());
        }
    }
    if (m_use_fiducial) {
        SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino: match_isFC uses the configured fiducial with {} tolerance value(s)",
                            m_fv_tolerance.size());
    }
    m_sp_sce_correction            = get(config, "sp_sce_correction",            m_sp_sce_correction);
    m_tagger_ordered_segment_sets  = get(config, "tagger_ordered_segment_sets",  m_tagger_ordered_segment_sets);
    m_stem_endpoint_wcpt_parity    = get(config, "stem_endpoint_wcpt_parity",    m_stem_endpoint_wcpt_parity);
    m_broken_muon_cluster_id_count = get(config, "broken_muon_cluster_id_count", m_broken_muon_cluster_id_count);
    m_neutrino_type_bitmask        = get(config, "neutrino_type_bitmask",        m_neutrino_type_bitmask);
    // doc sbnd_xin/docs/pr/33 §10 EM-shower-clustering knobs.
    m_daughter_count_proto_main_vertex     = get(config, "daughter_count_proto_main_vertex",     m_daughter_count_proto_main_vertex);
    m_daughter_count_proto_examine_showers = get(config, "daughter_count_proto_examine_showers", m_daughter_count_proto_examine_showers);
    m_shower_pdg_from_start_segment        = get(config, "shower_pdg_from_start_segment",        m_shower_pdg_from_start_segment);
    m_shower_pdg_from_shower_type          = get(config, "shower_pdg_from_shower_type",          m_shower_pdg_from_shower_type);
    m_shower_pdg_exact_muon_test           = get(config, "shower_pdg_exact_muon_test",           m_shower_pdg_exact_muon_test);
    m_pi0_id_shared_allocator              = get(config, "pi0_id_shared_allocator",              m_pi0_id_shared_allocator);
    m_shower_flag_pdg_electron             = get(config, "shower_flag_pdg_electron",             m_shower_flag_pdg_electron);
    m_shower_less_id_tiebreak              = get(config, "shower_less_id_tiebreak",              m_shower_less_id_tiebreak);
    m_shower_endpoint_exclude_start_vertex = get(config, "shower_endpoint_exclude_start_vertex", m_shower_endpoint_exclude_start_vertex);
    // doc sbnd_xin/docs/pr/40 -- track mis-identified as electron.
    m_track_pid_persist_dqdx      = get(config, "track_pid_persist_dqdx",      m_track_pid_persist_dqdx);
    m_shower_reclass_dqdx_guard   = get(config, "shower_reclass_dqdx_guard",   m_shower_reclass_dqdx_guard);
    m_shower_topo_dqdx_guard      = get(config, "shower_topo_dqdx_guard",      m_shower_topo_dqdx_guard);
    // doc sbnd_xin/docs/pr/40 round 2 -- follow-on fixes to the pr/40 round.
    m_track_pid_persist_4mom      = get(config, "track_pid_persist_4mom",      m_track_pid_persist_4mom);
    m_shower_proton_daughter_pion = get(config, "shower_proton_daughter_pion", m_shower_proton_daughter_pion);
    // doc sbnd_xin/docs/pr/40 round 4 -- two follow-on fixes to F5.
    m_shower_proton_daughter_pion_dissolve = get(config, "shower_proton_daughter_pion_dissolve", m_shower_proton_daughter_pion_dissolve);
    m_muon_multi_proton_pion               = get(config, "muon_multi_proton_pion",               m_muon_multi_proton_pion);
    m_track_pid_persist_dqdx_electron_guard     = get(config, "track_pid_persist_dqdx_electron_guard",     m_track_pid_persist_dqdx_electron_guard);
    m_shower_connect_main_vertex_straight_guard = get(config, "shower_connect_main_vertex_straight_guard", m_shower_connect_main_vertex_straight_guard);
    m_shower_traj_straight_guard                = get(config, "shower_traj_straight_guard",                m_shower_traj_straight_guard);
    m_shower_absorb_track_guard                 = get(config, "shower_absorb_track_guard",                 m_shower_absorb_track_guard);
    // doc sbnd_xin/docs/pr/65 round 3 -- reachability-relaxed absorber guards.
    m_shower_absorb_unreachable_main            = get(config, "shower_absorb_unreachable_main",            m_shower_absorb_unreachable_main);
    m_shower_connect_protected_pion_guard       = get(config, "shower_connect_protected_pion_guard",       m_shower_connect_protected_pion_guard);
    m_michel_stem_muon_rescue                   = get(config, "michel_stem_muon_rescue",                   m_michel_stem_muon_rescue);
    m_shower_in_cascade_guard                   = get(config, "shower_in_cascade_guard",                   m_shower_in_cascade_guard);
    m_shower_in_max_len                         = get(config, "shower_in_max_len",                         m_shower_in_max_len);
    m_shower_in_mip_hi                          = get(config, "shower_in_mip_hi",                          m_shower_in_mip_hi);
    m_michel_stem_michel_check                  = get(config, "michel_stem_michel_check",                  m_michel_stem_michel_check);
    m_michel_stem_max_far_len                   = get(config, "michel_stem_max_far_len",                   m_michel_stem_max_far_len);
    m_shower_stem_backfill                      = get(config, "shower_stem_backfill",                      m_shower_stem_backfill);
    m_stem_backfill_max_len                     = get(config, "stem_backfill_max_len",                     m_stem_backfill_max_len);
    m_stem_backfill_mip_lo                      = get(config, "stem_backfill_mip_lo",                      m_stem_backfill_mip_lo);
    m_stem_backfill_mip_hi                      = get(config, "stem_backfill_mip_hi",                      m_stem_backfill_mip_hi);
    m_stem_backfill_min_shower_len              = get(config, "stem_backfill_min_shower_len",              m_stem_backfill_min_shower_len);
    m_shower_conn3_unreachable                  = get(config, "shower_conn3_unreachable",                  m_shower_conn3_unreachable);
    m_conn3_unreachable_min_len                 = get(config, "conn3_unreachable_min_len",                 m_conn3_unreachable_min_len);
    m_shower_traj_michel_stem                   = get(config, "shower_traj_michel_stem",                   m_shower_traj_michel_stem);
    m_michel_stem_traj_min_len                  = get(config, "michel_stem_traj_min_len",                  m_michel_stem_traj_min_len);
    m_michel_stem_traj_max_len                  = get(config, "michel_stem_traj_max_len",                  m_michel_stem_traj_max_len);
    m_michel_stem_traj_mip_lo                   = get(config, "michel_stem_traj_mip_lo",                   m_michel_stem_traj_mip_lo);
    m_michel_stem_traj_max_far_len              = get(config, "michel_stem_traj_max_far_len",              m_michel_stem_traj_max_far_len);
    m_michel_stem_traj_min_kink_deg             = get(config, "michel_stem_traj_min_kink_deg",             m_michel_stem_traj_min_kink_deg);
    m_shower_long_muon_keep_type                = get(config, "shower_long_muon_keep_type",                m_shower_long_muon_keep_type);
    m_single_muon_proton_chain_veto             = get(config, "single_muon_proton_chain_veto",             m_single_muon_proton_chain_veto);
    m_single_muon_long_muon_claim               = get(config, "single_muon_long_muon_claim",               m_single_muon_long_muon_claim);
    m_pid_flag_reconcile                        = get(config, "pid_flag_reconcile",                        m_pid_flag_reconcile);
    m_long_muon_stub_bridge                     = get(config, "long_muon_stub_bridge",                     m_long_muon_stub_bridge);

    if (!m_trackfitting_config_file.empty()) {
        load_trackfitting_config(m_trackfitting_config_file);
    }

    NeedDV::configure(config);
    NeedPCTS::configure(config);
    NeedRecombModel::configure(config);
    NeedParticleData::configure(config);
    NeedClusGeomHelper::configure(config);
}

Configuration TaggerCheckNeutrino::default_configuration() const
{
    Configuration cfg;
    cfg["grouping"] = m_grouping_name;
    cfg["detector_volumes"] = "DetectorVolumes";
    cfg["pc_transforms"] = "PCTransformSet";
    cfg["recombination_model"] = "BoxRecombination";  
    cfg["particle_dataset"] = "ParticleDataSet"; 

    cfg["trackfitting_config_file"] = "";
    cfg["perf"] = m_perf;
    cfg["dir_weak_use_score"] = m_dir_weak_use_score;  // false = legacy raw dir_weak flag reads
    cfg["mip_dqdx"]            = m_mip_dqdx;             // e/cm; 50000 = uBooNE legacy flat-template scale
    cfg["mip_dqdx_median"]     = m_mip_dqdx_median;      // e/cm; 43000 = uBooNE legacy median-threshold scale
    cfg["proton_dir_vote"]      = m_proton_dir_vote;      // false = legacy muon-vs-flat-only direction
    cfg["proton_dir_score_max"] = m_proton_dir_score_max;
    cfg["proton_dir_asym_min"]  = m_proton_dir_asym_min;
    cfg["endpoint_trim_retry"]  = m_endpoint_trim_retry;  // false = legacy (no endpoint-trim retry on abstention)
    cfg["fit_vertex_min_seg_length"] = m_fit_vertex_min_seg_length;  // cm; 0 = legacy (all segments enter the vertex fit)
    // doc sbnd_xin/docs/pr/51 round 7: robust vertex fit; false = legacy
    // (AddSegment epilogue never runs, byte-identical)
    cfg["mvfit_robust"]      = m_mvfit_robust;
    cfg["mvfit_main_only"]   = m_mvfit_main_only;    // true = main (neutrino) vertex only
    cfg["mvfit_min_len"]     = m_mvfit_min_len;      // cm, fits-chord gate
    cfg["mvfit_rin_margin"]  = m_mvfit_rin_margin;   // cm past the re-seat radius
    cfg["mvfit_rout_frac"]   = m_mvfit_rout_frac;    // rout = clamp(frac*chord, min, max)
    cfg["mvfit_rout_min"]    = m_mvfit_rout_min;     // cm
    cfg["mvfit_rout_max"]    = m_mvfit_rout_max;     // cm
    cfg["mvfit_angle"]       = m_mvfit_angle;        // deg, folded inner-vs-outer disagreement
    cfg["mvfit_min_pts"]     = m_mvfit_min_pts;      // outer-window point floor
    cfg["mvfit_min_aniso"]   = m_mvfit_min_aniso;    // sqrt(l0/l1) floor (shower guard)
    cfg["mvfit_prior_range"] = m_mvfit_prior_range;  // cm, 2-leg substituted prior
    cfg["cathode_x"]         = m_cathode_x;          // cm, T0-corrected frame
    cfg["cathode_kink_xcut"] = m_cathode_kink_xcut;  // cm; 0 = legacy (the kink search sees every fit point)
    cfg["cathode_wide_kink_angle"]    = m_cathode_wide_kink_angle;    // deg; 0 = legacy (no wide-baseline cathode accept)
    cfg["cathode_wide_kink_skirt"]    = m_cathode_wide_kink_skirt;    // cm excluded around the crossing
    cfg["cathode_wide_kink_baseline"] = m_cathode_wide_kink_baseline; // cm PCA baseline per arm beyond the skirt
    // doc sbnd_xin/docs/pr/48 -- back-to-back track fixes.  Defaults OFF/legacy.
    cfg["two_end_break"]       = m_two_end_break;       // false = legacy (no two-end dQ/dx break pass)
    cfg["teb_min_len"]         = m_teb_min_len;         // cm
    cfg["teb_min_arm"]         = m_teb_min_arm;         // cm
    cfg["teb_min_arm_pts"]     = m_teb_min_arm_pts;
    cfg["teb_stub_max"]        = m_teb_stub_max;        // cm
    cfg["teb_accept_range"]    = m_teb_accept_range;    // cm
    cfg["teb_rise_r1"]         = m_teb_rise_r1;
    cfg["teb_rise_r2"]         = m_teb_rise_r2;
    cfg["teb_abs_end_min"]     = m_teb_abs_end_min;     // x mip_dqdx_median
    cfg["teb_dip_floor"]       = m_teb_dip_floor;       // x mip_dqdx_median
    cfg["teb_score_cap_r1"]    = m_teb_score_cap_r1;
    cfg["teb_score_cap_r2"]    = m_teb_score_cap_r2;
    cfg["teb_turn_angle"]      = m_teb_turn_angle;      // deg; <= 0 disables route R2
    cfg["teb_turn_baseline"]   = m_teb_turn_baseline;   // cm
    cfg["teb_turn_skirt"]      = m_teb_turn_skirt;      // cm
    cfg["teb_turn_min_arm_frac"] = m_teb_turn_min_arm_frac; // frac of teb_turn_baseline; 0 = legacy argmax (doc pr/90 round 2)
    cfg["teb_second_max"]      = m_teb_second_max;      // cm; 0 = legacy strict gate (doc pr/90 round 2)
    cfg["teb_chain_topology"]  = m_teb_chain_topology;  // simple-path gate admission; false = legacy gate (doc pr/90 round 4)
    cfg["teb_r3_turn"]         = m_teb_r3_turn;         // deg; <= 0 disables route R3 (doc pr/90 round 4)
    cfg["teb_r3_hot"]          = m_teb_r3_hot;          // x mip_dqdx_median; <= 0 disables route R3 (doc pr/90 round 4)
    cfg["teb_bragg_veto_turn"] = m_teb_bragg_veto_turn; // deg; <= 0 disables the R2 bragg veto (doc pr/90 round 4)
    cfg["kink_walk_dqdx_stop"] = m_kink_walk_dqdx_stop; // false = legacy (flag_search bypasses the walk gate)
    cfg["kink_break_protect"]  = m_kink_break_protect;  // false = legacy (no protected kink breaks)
    // doc sbnd_xin/docs/pr/50 -- main-vertex kink-consistency snap; false =>
    // the pass never fires => byte-identical.  Numerics cm/deg.
    cfg["vertex_kink_snap"] = m_vertex_kink_snap;
    cfg["vks_radius"]       = m_vks_radius;
    cfg["vks_min_dis"]      = m_vks_min_dis;
    cfg["vks_angle"]        = m_vks_angle;
    cfg["vks_margin"]       = m_vks_margin;
    cfg["vks_collinear"]    = m_vks_collinear;
    cfg["vks_skirt"]        = m_vks_skirt;
    cfg["vks_baseline"]     = m_vks_baseline;
    cfg["vks_min_arm"]      = m_vks_min_arm;
    cfg["vks_fit_miss"]     = m_vks_fit_miss;
    cfg["vks_hot_ratio"]    = m_vks_hot_ratio;
    cfg["vks_carry_prong"]  = m_vks_carry_prong;  // cm; 0 = off, byte-identical (doc pr/85)
    // doc sbnd_xin/docs/pr/51 -- main-vertex graph audit; false => the pass
    // never fires => byte-identical.  Numerics cm/deg.
    cfg["main_vertex_graph_audit"] = m_main_vertex_graph_audit;
    cfg["mvga_radius"]       = m_mvga_radius;
    cfg["mvga_dup_tol"]      = m_mvga_dup_tol;
    cfg["mvga_dup_frac"]     = m_mvga_dup_frac;
    cfg["mvga_dup_angle"]    = m_mvga_dup_angle;
    cfg["mvga_bridge_mip"]   = m_mvga_bridge_mip;
    cfg["mvga_reconnect"]    = m_mvga_reconnect;
    cfg["mvga_stub"]         = m_mvga_stub;
    cfg["mvga_stub_pts"]     = m_mvga_stub_pts;
    cfg["mvga_reseat_angle"] = m_mvga_reseat_angle;
    cfg["mvga_satellite"]    = m_mvga_satellite;  // 0 = main-vertex-only (round 2), byte-identical
    cfg["mvga_interposed"]   = m_mvga_interposed; // false = terminal-only op3, byte-identical (doc pr/85)
    cfg["mvga_interposed_angle"] = m_mvga_interposed_angle;  // deg; inert while mvga_interposed is false
    cfg["mvga_interposed_len"] = m_mvga_interposed_len;  // cm; 0 = use mvga_stub, byte-identical (doc pr/86)
    cfg["mvga_sat_dup_frac"] = m_mvga_sat_dup_frac;  // fraction; 0 = use mvga_dup_frac, byte-identical (doc pr/86)
    cfg["mvga_interposed_deg1"] = m_mvga_interposed_deg1;  // false = degree-1 anchors stay out of reach, byte-identical (doc pr/86)
    cfg["mvga_splice_straighten"] = m_mvga_splice_straighten;  // cm; 0 = concatenation verbatim, byte-identical (doc pr/86 round 2)
    cfg["mvga_approach_collapse"] = m_mvga_approach_collapse;  // cm; 0 = op3.5 skipped, byte-identical (doc pr/86 round 2)
    cfg["mvga_straighten_radius"] = m_mvga_straighten_radius;  // cm; 0 = prototype 0.2, inert unless straighten/collapse on (doc pr/86 round 2)
    cfg["kink_dqdx_hot_ratio"] = m_kink_dqdx_hot_ratio; // x mip_dqdx_median; inert while both above are false
    cfg["shower_topo_demote_len"] = m_shower_topo_demote_len;  // cm; 0 = legacy (long segments stay eligible for kShowerTopology)
    // doc sbnd_xin/docs/pr/49.
    cfg["fit_blob_coverage"] = m_fit_blob_coverage; // -1 = legacy (no foreign-ghost deweighting); >= 0 = tolerance cells
    // doc sbnd_xin/docs/pr/50.
    cfg["fit_blob_coverage_defer"] = m_fit_blob_coverage_defer; // false = pr/49 behavior (deweight active during find_proto_vertex)
    // doc sbnd_xin/docs/pr/30 §11.  Round-tripped so the compiled config
    // records the operating point; each default reproduces the pre-pr/30 tree.
    cfg["fit_exclusion"]            = m_fit_exclusion;             // false = legacy (all sites pass flag_exclusion=false)
    cfg["graph_endpoint_strict"]    = m_graph_endpoint_strict;     // false = legacy (WARN only, connection still made)
    cfg["graph_endpoint_tol"]       = m_graph_endpoint_tol;        // cm
    cfg["oov_prototype_parity"]     = m_oov_prototype_parity;      // false = legacy (today's three polarities)
    cfg["first_seg_local_pca"]      = m_first_seg_local_pca;       // true  = legacy (the refinement runs)
    cfg["other_seg_relaxed_accept"] = m_other_seg_relaxed_accept;  // true  = legacy (the 0.72/15cm/1.05 clause is live)
    cfg["other_seg_empty_2d_guard"] = m_other_seg_empty_2d_guard;  // false = legacy (-1 sentinel counts as covered)
    // doc sbnd_xin/docs/pr/54.
    cfg["other_seg_keep_isolated"]            = m_other_seg_keep_isolated;            // false = legacy (isolated residual discarded)
    cfg["other_seg_keep_isolated_min_points"] = m_other_seg_keep_isolated_min_points; // component-point floor when the keep is on
    cfg["other_seg_keep_isolated_min_length"] = m_other_seg_keep_isolated_min_length; // cm; fitted-length floor when the keep is on
    // doc sbnd_xin/docs/pr/67 round 3.
    cfg["iso_snap_min_dir_mag"] = m_iso_snap_min_dir_mag; // cm; 10.0 = legacy isochronous-snap size gate
    // doc sbnd_xin/docs/pr/59 round 2.
    cfg["assoc_full_recluster"] = m_assoc_full_recluster; // false = legacy (orphaned associate_points cloud stays null)
    // doc sbnd_xin/docs/pr/64 round 7.
    cfg["assoc_reassign_orphans"] = m_assoc_reassign_orphans; // false = legacy (Stage-C loss is dropped, never reassigned)
    // doc sbnd_xin/docs/pr/64 round 8.
    cfg["assoc_clear_on_merge"] = m_assoc_clear_on_merge; // false = legacy (merge survivor keeps its stale associate_points)
    // doc sbnd_xin/docs/pr/31 §11.
    cfg["shower_topo_proto_dir"]    = m_shower_topo_proto_dir;     // false = legacy (the stage-3 PCA direction call runs)
    // doc sbnd_xin/docs/pr/32 §11.
    cfg["vertex_dir_use_fit_point"]       = m_vertex_dir_use_fit_point;       // false = legacy (raw wcpt at the 11 sites)
    cfg["shower_traj_recheck_parity"]     = m_shower_traj_recheck_parity;     // false = legacy (recomputed gate, 1 cm inner)
    cfg["main_vertex_require_descriptor"] = m_main_vertex_require_descriptor; // false = legacy (unguarded argmax)
    cfg["main_vertex_candidate_flag"]     = m_main_vertex_candidate_flag;     // false = legacy (no kMainCandidate)
    // doc sbnd_xin/docs/pr/31 §12.
    cfg["cont_muon_dir3_30cm"]            = m_cont_muon_dir3_30cm;            // false = legacy (15 cm fallback for short reference segments)
    cfg["track_comp_empty_abstain"]       = m_track_comp_empty_abstain;       // false = legacy (empty window "confirms" the direction)
    cfg["shower_topo_reset"]              = m_shower_topo_reset;              // false = legacy (set-only kShowerTopology, stale dirsign)
    cfg["reclass_preserve_4mom"]          = m_reclass_preserve_4mom;          // false = legacy (unconditional 4-momentum rewrite)
    cfg["reclass_never_computed_ke_floor"] = m_reclass_never_computed_ke_floor; // false = legacy (never-computed reads KE = -mass)
    cfg["dir_track_median_local"]         = m_dir_track_median_local;         // false = legacy (filtered median helper)
    cfg["examine_showers_vertex_by_index"] = m_examine_showers_vertex_by_index; // false = legacy (proximity-ordered pair)
    cfg["iso_endpoint"]               = m_iso_endpoint;                // false = legacy wire-footprint boundary endpoints
    cfg["iso_endpoint_min_length"]    = m_iso_endpoint_min_length;     // cm
    cfg["iso_endpoint_max_xext"]      = m_iso_endpoint_max_xext;       // cm
    cfg["iso_endpoint_xext_frac"]     = m_iso_endpoint_xext_frac;
    cfg["iso_endpoint_xext_quantile"] = m_iso_endpoint_xext_quantile;
    cfg["iso_endpoint_tube_radius"]   = m_iso_endpoint_tube_radius;    // cm
    cfg["iso_endpoint_min_aspect"]    = m_iso_endpoint_min_aspect;     // trimmed transverse/axial extent ratio
    cfg["traj_cover_probe"]           = m_traj_cover_probe;            // false = no pr/67 diagnostic lines
    cfg["pr_find_other_rounds"]       = m_pr_find_other_rounds;        // 0 = keep find_proto_vertex's hardcoded budget
    cfg["v3_extension_guard"]         = m_v3_extension_guard;          // false = examine_vertices_3 unconditional accept
    cfg["v3_extension_min_gain"]      = m_v3_extension_min_gain;       // cm
    // doc sbnd_xin/docs/pr/72 round 2 -- examine_structure_3 stub guard.
    cfg["es3_stub_guard"]         = m_es3_stub_guard;         // false = examine_structure_3 unconditional merge on angle alone
    cfg["es3sg_stub_max"]         = m_es3sg_stub_max;         // cm
    cfg["es3sg_len_ratio"]        = m_es3sg_len_ratio;
    cfg["es3sg_ang3_min"]         = m_es3sg_ang3_min;         // degrees
    cfg["es3sg_ang_ratio"]        = m_es3sg_ang_ratio;
    cfg["es3sg_require_terminal"] = m_es3sg_require_terminal;
    // Detector-extent literals (docs/pr/2 sec. 2e(iv)); defaults = uBooNE prototype, cm.
    cfg["cosmic_y_top_main"]    = m_cosmic_y_top_main;     // 100 = 17 cm below the uBooNE y=+117 top
    cfg["cosmic_y_top_strict"]  = m_cosmic_y_top_strict;   // 102 = 15 cm below
    cfg["cosmic_y_top_loose"]   = m_cosmic_y_top_loose;    // 80  = 37 cm below
    cfg["cosmic_y_small_piece"] = m_cosmic_y_small_piece;  // 50  = 67 cm below
    cfg["vertex_z_prior_scale"] = m_vertex_z_prior_scale;  // cm; 200 = uBooNE (1037 cm detector)
    // SSM beam-line references, {x,y,z}; defaults = uBooNE BNB target / NuMI absorber.
    // Assign the array first: append() alone would accumulate rather than overwrite if
    // this ever ran against a reused Configuration, and a 6-element array would be
    // rejected by get_dir3() -- a key that looks set but silently is not.
    cfg["ssm_target_dir"]   = Json::arrayValue;
    cfg["ssm_absorber_dir"] = Json::arrayValue;
    for (double c : m_ssm_target_dir)   cfg["ssm_target_dir"].append(c);
    for (double c : m_ssm_absorber_dir) cfg["ssm_absorber_dir"].append(c);
    // Charge -> energy calibration (docs/pr/2 sec. 2e(iii)); defaults = uBooNE.
    cfg["kine_fudge_factor"]        = m_kine_fudge_factor;         // track-like residual scale
    cfg["kine_recom_factor"]        = m_kine_recom_factor;         // track-like recombination survival
    cfg["kine_shower_fudge_factor"] = m_kine_shower_fudge_factor;  // shower-flagged counterpart
    cfg["kine_shower_recom_factor"] = m_kine_shower_recom_factor;
    cfg["kine_proton_recom_factor"] = m_kine_proton_recom_factor;  // |pdg|==2212 (fudge stays at kine_fudge_factor)
    cfg["kine_plane_weights"] = Json::arrayValue;                  // {U,V,W} charge-average weights
    for (double w : m_kine_plane_weights) cfg["kine_plane_weights"].append(w);
    cfg["kine_plane_asym_switch"]   = m_kine_plane_asym_switch;    // (med,max) asymmetry above which the max plane is dropped
    cfg["kine_w_value"]             = m_kine_w_value;              // eV per electron-ion pair
    cfg["kine_shower_pdg_live"]     = m_kine_shower_pdg_live;      // live start-segment PDG at fill_kine_tree (doc pr/35 §10.2)
    cfg["muon_dqdx_curve"] = Json::arrayValue;                     // {c0, c1, pivot_cm, power} of the muon
    for (double c : m_muon_dqdx_curve) cfg["muon_dqdx_curve"].append(c);  // median-dQ/dx-vs-length envelope
    cfg["sp_dedx_use_recomb_model"] = m_sp_dedx_use_recomb_model;  // false = inline uBooNE-field inverse Box
    cfg["sp_mean_dedx_cut"]         = m_sp_mean_dedx_cut;          // MeV/cm; coupled to the knob above
    cfg["dl_weights"] = "";       // empty = DL vertex disabled
    cfg["dl_vtx_cut"] = 25.0;    // mm (= 2.5 cm)
    cfg["dQdx_scale"]  = 0.1;    // dQ scale factor for SCN network input
    cfg["dQdx_offset"] = -1000.0; // dQ offset for SCN network input
    cfg["dl_vtx_rerank"]           = true;    // true → use top-K + soft re-rank; false → legacy single argmax
    cfg["dl_vtx_top_k"]            = 5;       // number of top DL voxels to re-rank (only when dl_vtx_rerank==true)
    cfg["dl_vtx_min_accept_score"] = 4.0;     // min composite score to accept a re-ranked DL vertex (empirical; correct uncertain-regime picks score 8-12, failure cases 3-5)
    cfg["dl_vtx_score_scale"]      = 1000.0;  // scale factor on raw DL score in composite re-rank (1.0 = unscaled)
    cfg["dl_vtx_swap_guard"]       = m_dl_vtx_swap_guard;  // doc pr/51 (506746): false = legacy (rerank may swap the main cluster)
    cfg["dl_vtx_topo_weight"]      = m_dl_vtx_topo_weight; // doc pr/89 C2: 0 = legacy (rule-1 topology term never computed); offline C1 selected 3.0
    cfg["dl_vtx_topo_center"]      = m_dl_vtx_topo_center; // doc pr/89 C2: frac offset; frozen choice is 0.0 (the center-0.5 variant lost)
    cfg["main_vertex_swap_apply"]  = m_main_vertex_swap_apply;  // doc pr/51 round 3: false = legacy (traditional-path swap decision is computed then discarded)
    cfg["rough_path_probe"]        = m_rough_path_probe;  // doc pr/51 round 4: false = legacy (diagnostic TRACE probe never runs)
    cfg["steiner_gap_penalty"]     = m_steiner_gap_penalty;  // doc pr/51 round 5: 0 = legacy (do_rough_path stays on the unpenalized "steiner_graph")
    cfg["sgp_dead_alpha"]          = m_sgp_dead_alpha;       // doc pr/51 round 5: dead-sample weight in bad_fraction (inert at scale 0)
    cfg["sgp_min_edge"]            = m_sgp_min_edge;         // doc pr/51 round 5: cm; shorter edges never scanned (inert at scale 0)
    cfg["sgp_sample_step"]         = m_sgp_sample_step;      // doc pr/51 round 5: cm; edge-interior sampling step (inert at scale 0)
    cfg["sgp_point_radius"]        = m_sgp_point_radius;     // doc pr/51 round 5: cm; test_good_point radius (inert at scale 0)
    cfg["sgp_edge_probe"]          = m_sgp_edge_probe;       // doc pr/73: false = legacy (per-edge DEBUG sentinel never emits)
    cfg["vertex_scoreboard"]       = m_vertex_scoreboard;    // doc pr/75: false = legacy (no vertex scoreboard recorded)
    cfg["dl_vtx_harvest"]          = m_dl_vtx_harvest;       // doc pr/79 §10: false = legacy (no live-feature harvest; requires vertex_scoreboard)
    cfg["sgp_weak_scale"]          = m_sgp_weak_scale;       // doc pr/51 round 6: 0 = legacy (round-5 gap flavor verbatim)
    cfg["sgp_weak_qref"]           = m_sgp_weak_qref;        // doc pr/51 round 6: charge ref, calc_charge_wcp units (inert at weak scale 0)
    cfg["sgp_max_sep"]             = m_sgp_max_sep;          // doc pr/73 round 2 F3a: cm; NEGATIVE = legacy (no cap). 0 is a real cap, so the off-test is < 0, not <= 0
    cfg["break_seg_orient"]        = m_break_seg_orient;     // doc pr/83: false = legacy (break_segment slices by boost source/target)
    cfg["clus_geom_helper"] = ""; // empty = no SCE vertex correction
    cfg["beam_window_low"] = m_beam_window_low;   // beam window on cluster_t0; low >= high disables the
    cfg["beam_window_high"] = m_beam_window_high; // gate (uBooNE single-main selection).
    cfg["nu_skip_cosmic"] = m_nu_skip_cosmic;     // beam-gate only: skip in-window mains with flag_TGM/flag_STM/lm_flag>0
    cfg["nu_skip_cosmic_bundle"] = m_nu_skip_cosmic_bundle;  // that verdict vetoes the whole matched_flash_gid bundle
    cfg["nu_skip_cosmic_bundle_min_length"] = m_nu_skip_cosmic_bundle_min_length;  // cm; > 0 spares untagged bundle-mates at least this long (docs/pr/16 design A); 0 = veto all
    cfg["skip_cosmic_companions"] = m_skip_cosmic_companions;          // doc pr/20 I P4; drop a TGM/STM-tagged companion from other_clusters
    cfg["cosmic_companion_min_length"] = m_cosmic_companion_min_length;  // cm; a tagged companion shorter than this stays in regardless
    cfg["sp_photon_flag"] = m_sp_photon_flag;     // doc pr/26 sec. 8.2; store singlephoton_tagger()'s verdict in TaggerInfo::photon_flag (prototype NeutrinoID.cxx:271)
    // doc sbnd_xin/docs/pr/36 §10.
    cfg["fiducial"] = Json::Value();                 // null = the historical FiducialUtils containment fallback
    cfg["fv_tolerance"] = Json::Value(Json::arrayValue);  // [x_lo,x_hi,y_lo,y_hi,z_lo,z_hi] margins, negative = inset
    cfg["sp_sce_correction"]            = m_sp_sce_correction;            // false = legacy raw single-photon positions
    cfg["tagger_ordered_segment_sets"]  = m_tagger_ordered_segment_sets;  // false = legacy pointer-address iteration (M4 fix when true)
    cfg["stem_endpoint_wcpt_parity"]    = m_stem_endpoint_wcpt_parity;    // false = legacy nearest-fit-endpoint proximity rule
    cfg["broken_muon_cluster_id_count"] = m_broken_muon_cluster_id_count; // false = legacy distinct-pointer cluster count
    cfg["neutrino_type_bitmask"]        = m_neutrino_type_bitmask;        // false = legacy (no verdict bitmask, no T_tagger branch)
    // doc sbnd_xin/docs/pr/33 §10.
    cfg["daughter_count_proto_main_vertex"]     = m_daughter_count_proto_main_vertex;     // false = legacy _showers callee at the main-vertex site
    cfg["daughter_count_proto_examine_showers"] = m_daughter_count_proto_examine_showers; // false = legacy _showers callee in examine_showers
    cfg["shower_pdg_from_start_segment"]        = m_shower_pdg_from_start_segment;        // false = legacy shower->get_particle_type() at the 4 sites
    cfg["shower_pdg_from_shower_type"]          = m_shower_pdg_from_shower_type;          // false = legacy start-segment read at the inverted site
    cfg["shower_pdg_exact_muon_test"]           = m_shower_pdg_exact_muon_test;           // false = legacy abs() muon test (parity at :170 needs from_start_segment too)
    cfg["pi0_id_shared_allocator"]              = m_pi0_id_shared_allocator;              // false = legacy independent pi0-id seeds (collision possible)
    cfg["shower_flag_pdg_electron"]             = m_shower_flag_pdg_electron;             // false = legacy is_shower without the abs(pdg)==11 term
    cfg["shower_less_id_tiebreak"]              = m_shower_less_id_tiebreak;              // false = legacy pointer-address tie-break (house-rule fix when true)
    cfg["shower_endpoint_exclude_start_vertex"] = m_shower_endpoint_exclude_start_vertex; // false = legacy end_point search includes the start vertex (doc pr/39)
    // doc sbnd_xin/docs/pr/40 -- track mis-identified as electron.
    cfg["track_pid_persist_dqdx"]    = m_track_pid_persist_dqdx;    // false = legacy free-end-gated persistence
    cfg["shower_reclass_dqdx_guard"] = m_shower_reclass_dqdx_guard; // false = legacy unconditional wholesale reclassification
    cfg["shower_topo_dqdx_guard"]    = m_shower_topo_dqdx_guard;    // false = legacy (dQ/dx never consulted by the topology test)
    // doc sbnd_xin/docs/pr/40 round 2 -- two follow-on defects from the pr/40 fix round.
    cfg["track_pid_persist_4mom"]      = m_track_pid_persist_4mom;      // false = legacy rest-mass-only 4-mom stub (zero KE)
    cfg["shower_proton_daughter_pion"] = m_shower_proton_daughter_pion; // false = legacy (proton daughter never consulted)
    // doc sbnd_xin/docs/pr/40 round 4 -- two follow-on fixes to F5.
    cfg["shower_proton_daughter_pion_dissolve"] = m_shower_proton_daughter_pion_dissolve; // false = legacy (F5 relabel leaves shower flags set)
    cfg["muon_multi_proton_pion"]               = m_muon_multi_proton_pion;               // false = legacy (multi-proton muon vertex never consulted)
    cfg["track_pid_persist_dqdx_electron_guard"]     = m_track_pid_persist_dqdx_electron_guard;     // false = legacy (F1 rescues an undirected electron guess too)
    cfg["shower_connect_main_vertex_straight_guard"] = m_shower_connect_main_vertex_straight_guard; // false = legacy (straightness never consulted)
    cfg["shower_traj_straight_guard"]                = m_shower_traj_straight_guard;                // false = legacy (straightness never consulted)
    cfg["shower_absorb_track_guard"]                 = m_shower_absorb_track_guard;                 // false = legacy (flood-fill absorbs every connected segment)
    cfg["shower_absorb_unreachable_main"]            = m_shower_absorb_unreachable_main;            // false = legacy (absorbers skip ALL main-cluster segments, doc pr/65)
    cfg["shower_connect_protected_pion_guard"]       = m_shower_connect_protected_pion_guard;       // false = legacy (proton-daughter pion selectable as EM candidate)
    cfg["michel_stem_muon_rescue"]                   = m_michel_stem_muon_rescue;                   // false = legacy (Michel rescue limited to weak-dir degree-2 vertices)
    cfg["shower_in_cascade_guard"]                   = m_shower_in_cascade_guard;                   // doc pr/74 round 2 P1; false = legacy (cascade relabels unconditionally)
    cfg["shower_in_max_len"]                         = m_shower_in_max_len;                         // cm; only read when shower_in_cascade_guard
    cfg["shower_in_mip_hi"]                          = m_shower_in_mip_hi;                          // ratio; only read when shower_in_cascade_guard
    cfg["michel_stem_michel_check"]                  = m_michel_stem_michel_check;                  // doc pr/74 round 2 P2; false = legacy (any shower-like sibling passes)
    cfg["michel_stem_max_far_len"]                   = m_michel_stem_max_far_len;                   // cm; only read when michel_stem_michel_check
    cfg["shower_stem_backfill"]                      = m_shower_stem_backfill;                      // doc pr/74 round 2 K4; false = legacy (walked-past stems stay out of showers)
    cfg["stem_backfill_max_len"]                     = m_stem_backfill_max_len;                     // cm; only read when shower_stem_backfill
    cfg["stem_backfill_mip_lo"]                      = m_stem_backfill_mip_lo;                      // ratio; only read when shower_stem_backfill
    cfg["stem_backfill_mip_hi"]                      = m_stem_backfill_mip_hi;                      // ratio; only read when shower_stem_backfill
    cfg["stem_backfill_min_shower_len"]              = m_stem_backfill_min_shower_len;              // cm; only read when shower_stem_backfill
    cfg["shower_conn3_unreachable"]                  = m_shower_conn3_unreachable;                  // doc pr/74 round 2 K5 (pr/65 rung 2); false = legacy (unreachable segments stay PF-invisible)
    cfg["conn3_unreachable_min_len"]                 = m_conn3_unreachable_min_len;                 // cm; only read when shower_conn3_unreachable
    cfg["shower_traj_michel_stem"]                   = m_shower_traj_michel_stem;                   // doc pr/74 round 4 K6; false = legacy (a stopping muon + Michel stays one EM shower)
    cfg["michel_stem_traj_min_len"]                  = m_michel_stem_traj_min_len;                  // cm; only read when shower_traj_michel_stem
    cfg["michel_stem_traj_max_len"]                  = m_michel_stem_traj_max_len;                  // cm; only read when shower_traj_michel_stem
    cfg["michel_stem_traj_mip_lo"]                   = m_michel_stem_traj_mip_lo;                   // x MIP median; only read when shower_traj_michel_stem
    cfg["michel_stem_traj_max_far_len"]              = m_michel_stem_traj_max_far_len;              // cm; only read when shower_traj_michel_stem
    cfg["michel_stem_traj_min_kink_deg"]             = m_michel_stem_traj_min_kink_deg;             // deg; only read when shower_traj_michel_stem
    cfg["shower_long_muon_keep_type"]                = m_shower_long_muon_keep_type;                // false = legacy (long-muon pseudo-shower start segment majority-voted to e-)
    cfg["single_muon_proton_chain_veto"]             = m_single_muon_proton_chain_veto;             // false = legacy (1-hop proton veto only)
    cfg["single_muon_long_muon_claim"]               = m_single_muon_long_muon_claim;               // false = legacy (long-muon chain never claims the vertex muon slot)
    cfg["pid_flag_reconcile"]                        = m_pid_flag_reconcile;                        // false = legacy (no late reconciliation pass)
    cfg["long_muon_stub_bridge"]                     = m_long_muon_stub_bridge;                     // false = legacy (stub-blocked long-muon chains never form)


    return cfg;
}

void TaggerCheckNeutrino::visit(Ensemble& ensemble) const
{
    using Clock = std::chrono::steady_clock;
    using MS    = std::chrono::duration<double, std::milli>;
    auto t_total = Clock::now();
    auto t0 = Clock::now();

    // Configure the track fitter with detector volume
    m_track_fitter->set_detector_volume(m_dv);
    m_track_fitter->set_pc_transforms(m_pcts); 

    // Get the specified grouping (default: "live")
    auto groupings = ensemble.with_name(m_grouping_name);
    if (groupings.empty()) {
        return;
    }
    
    auto& grouping = *groupings.at(0);
    
    // Find clusters that have the main_cluster flag (set by clustering_recovering_bundle)
    Cluster* main_cluster = nullptr;
    std::vector<Cluster*> other_clusters;  // beam_flash clusters that are not the main cluster

    int nclusters = grouping.nchildren();
    int n_main_clusters = 0;
    int n_in_beam_clusters = 0;
    const bool beam_gate = m_beam_window_low < m_beam_window_high;
    if (!beam_gate) {
        for (auto* cluster : grouping.children()) {
            if (cluster->get_flag(Flags::main_cluster)) {
                main_cluster = cluster;
                n_main_clusters ++;
            }
            if (cluster->get_flag(Flags::beam_flash)) n_in_beam_clusters++;
        }
        for (auto* cluster : grouping.children()) {
            if (cluster != main_cluster && cluster->get_flag(Flags::beam_flash)) {
                other_clusters.push_back(cluster);
            }
        }
    }
    else {
        // Beam-bundle selection: among the (possibly many) matched main
        // clusters, neutrino PR runs on the one whose matched flash time
        // (cluster_t0) falls inside the beam window; ties broken by length.
        // Companions are the associated clusters sharing its matched_flash_gid.

        // nu_skip_cosmic_bundle: hoist the cosmic verdicts to bundle
        // granularity before selecting.  A flash bundle can hold more than one
        // main (SBND evt 18255/52195: gid 6 holds the TGM-tagged 513 cm
        // cathode-crosser 13 *and* a 5-point 1.7 cm shard 5); with the
        // per-main rule the shard becomes the nu candidate and drags the
        // bundle's untagged 400 cm associated muon through a full PR.  Keying
        // on matched_flash_gid is the same bundle definition the companion
        // loop below uses.  Empty set when the knob is off => no behavior
        // change.  int-keyed, so the iteration order is stable.
        std::set<int> cosmic_gids;
        if (m_nu_skip_cosmic && m_nu_skip_cosmic_bundle) {
            for (auto* cluster : grouping.children()) {
                if (!cluster->get_flag(Flags::main_cluster)) continue;
                const double t0c = cluster->get_cluster_t0();
                if (t0c < m_beam_window_low || t0c >= m_beam_window_high) continue;
                const int gid = cluster->get_scalar<int>("matched_flash_gid", -1);
                if (gid < 0) continue;
                if (cluster->get_flag(Flags::TGM) || cluster->get_flag(Flags::STM)
                    || cluster->get_scalar<int>("lm_flag", -1) > 0) {
                    cosmic_gids.insert(gid);
                }
            }
        }

        for (auto* cluster : grouping.children()) {
            if (!cluster->get_flag(Flags::main_cluster)) continue;
            n_main_clusters++;
            const double t0 = cluster->get_cluster_t0();
            if (t0 < m_beam_window_low || t0 >= m_beam_window_high) continue;
            n_in_beam_clusters++;
            if (m_nu_skip_cosmic) {
                // Honor upstream cosmic verdicts (TaggerCheckTGM/TaggerCheckSTM
                // flags, Q/L LM tagger lm_flag: 1=low-energy, 2=light-mismatch;
                // absent scalar = not evaluated, never skips).  Mirrors the
                // STM-skips-TGM idiom (TaggerCheckSTM).  Per-main, so a
                // cosmic-tagged longest bundle does not veto a clean runner-up.
                const bool tgm = cluster->get_flag(Flags::TGM);
                const bool stm = cluster->get_flag(Flags::STM);
                const int lm = cluster->get_scalar<int>("lm_flag", -1);
                if (tgm || stm || lm > 0) {
                    SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: in-window cluster {} (t0 {:.3f} us, L {:.1f} cm) cosmic-tagged (TGM={} STM={} lm_flag={}); skipping (nu_skip_cosmic)",
                                       cluster->get_cluster_id(), t0/units::us, cluster->get_length()/units::cm,
                                       tgm, stm, lm);
                    continue;
                }
                if (!cosmic_gids.empty()) {
                    const int gid = cluster->get_scalar<int>("matched_flash_gid", -1);
                    if (gid >= 0 && cosmic_gids.count(gid)) {
                        // min_length guard (docs/pr/16 design A): an untagged
                        // in-window main at least this long survives the
                        // bundle veto -- the taggers examined it and declined
                        // to tag it, and its cosmic-tagged bundle-mate stays
                        // out of the PR ensemble regardless (companions are
                        // gathered from associated-only clusters below).  The
                        // guard exists to keep vetoing unexamined shards
                        // (out-of-scope mains carry no verdict; SBND evt
                        // 18255/52195's 1.3 cm shard 5).  0 (default) = veto
                        // every bundle-mate, byte-identical to pre-knob.
                        const double min_len = m_nu_skip_cosmic_bundle_min_length * units::cm;
                        if (min_len > 0 && cluster->get_length() >= min_len) {
                            SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: in-window cluster {} (t0 {:.3f} us, L {:.1f} cm) shares flash bundle gid {} with a cosmic-tagged main but is untagged and >= {:.1f} cm; kept (nu_skip_cosmic_bundle_min_length)",
                                               cluster->get_cluster_id(), t0/units::us,
                                               cluster->get_length()/units::cm, gid,
                                               m_nu_skip_cosmic_bundle_min_length);
                        }
                        else {
                            SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: in-window cluster {} (t0 {:.3f} us, L {:.1f} cm) shares flash bundle gid {} with a cosmic-tagged main; skipping (nu_skip_cosmic_bundle)",
                                               cluster->get_cluster_id(), t0/units::us,
                                               cluster->get_length()/units::cm, gid);
                            continue;
                        }
                    }
                }
            }
            if (!main_cluster) {
                main_cluster = cluster;
            }
            else {
                Cluster* loser = cluster->get_length() > main_cluster->get_length() ? main_cluster : cluster;
                if (loser == main_cluster) main_cluster = cluster;
                SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: in-window bundle cluster {} (t0 {:.3f} us, L {:.1f} cm) not selected",
                                   loser->get_cluster_id(), loser->get_cluster_t0()/units::us, loser->get_length()/units::cm);
            }
        }
        if (main_cluster) {
            const int gid = main_cluster->get_scalar<int>("matched_flash_gid", -1);
            for (auto* cluster : grouping.children()) {
                if (cluster == main_cluster) continue;
                if (!cluster->get_flag(Flags::associated_cluster)) continue;
                if (gid < 0 || cluster->get_scalar<int>("matched_flash_gid", -1) == gid) {
                    // P4 (doc pr/20 Part I): a companion the cosmic taggers
                    // convicted is not neutrino charge.  Dropping it here keeps
                    // it out of the pattern-recognition input, so it cannot
                    // contribute a particle to the flow tree or its charge to
                    // kine_reco_Enu.  The length floor is the safety valve: a
                    // SHORT tagged companion stays in regardless of verdict, so
                    // a mis-tagged neutrino daughter is never silently lost and
                    // a bad tag on a fragment is bounded.  Nothing tags a
                    // companion unless the taggers ran with
                    // evaluate_demoted_mains, so this is inert without P3.
                    if (m_skip_cosmic_companions
                        && (cluster->get_flag(Flags::TGM) || cluster->get_flag(Flags::STM))
                        && cluster->get_length() >= m_cosmic_companion_min_length * units::cm) {
                        SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: companion cluster {} (L {:.1f} cm, TGM={} STM={}) dropped from other_clusters (skip_cosmic_companions, floor {:.1f} cm)",
                                           cluster->get_cluster_id(), cluster->get_length()/units::cm,
                                           cluster->get_flag(Flags::TGM), cluster->get_flag(Flags::STM),
                                           m_cosmic_companion_min_length);
                        continue;
                    }
                    other_clusters.push_back(cluster);
                }
            }
        }
    }

    if (!main_cluster) {
        SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino: no main cluster selected ({} mains, {} in-window); skipping.",
                            n_main_clusters, n_in_beam_clusters);
        return;
    }

    if (beam_gate) {
        SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: selected main cluster {} (t0 {:.3f} us, L {:.1f} cm, {} associated)",
                           main_cluster->get_cluster_id(), main_cluster->get_cluster_t0()/units::us,
                           main_cluster->get_length()/units::cm, other_clusters.size());
    }

    SPDLOG_LOGGER_TRACE(log, "Found {} clusters, {} main clusters, {} in-beam clusters, {} of blobs in main cluster id {}", nclusters, n_main_clusters, n_in_beam_clusters, main_cluster->nchildren(), main_cluster->get_cluster_id());


    // Debug dump (only when env var is set)
    if (main_cluster) {
        if (const char* dump_path = std::getenv("WCT_DUMP_INIT_FIRST_SEGMENT")) {
            DebugIO::dump_init_first_segment_inputs(
                dump_path, *main_cluster, main_cluster, true, *m_track_fitter);
        }
    }

    SPDLOG_LOGGER_TRACE(log, "Number of Main Clusters: {}", n_main_clusters);

    IndexedVertexSet vertices_in_long_muon;
    IndexedSegmentSet segments_in_long_muon;
    VertexPtr main_vertex = nullptr;
    ClusterVertexMap map_cluster_main_vertices;

    // Pre-load charge data for all beam-flash clusters once so that
    // do_multi_tracking calls throughout pattern recognition can use
    // flag_force_load_data=false and avoid redundant prepare_data() calls.
    {
        std::vector<WireCell::Clus::Facade::Cluster*> clusters_to_preload;
        clusters_to_preload.push_back(main_cluster);
        for (auto* c : other_clusters) clusters_to_preload.push_back(c);
        m_track_fitter->preload_clusters(clusters_to_preload);
    }
    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: preload_clusters took {} ms", MS(Clock::now() - t0).count());
    t0 = Clock::now();

    // Debug dump for unit-test fixture generation (only when env var is set).
    // Generate with:
    //   WCT_DUMP_TAGGER_INPUTS=./tmp/tagger_check_neutrino_input.json wire-cell ...
    // Then replay with doctest_tagger_check_neutrino end-to-end test.
    if (main_cluster) {
        if (const char* dump_path = std::getenv("WCT_DUMP_TAGGER_INPUTS")) {
            DebugIO::dump_tagger_inputs(
                dump_path, *main_cluster, other_clusters,
                /*flag_back_search=*/true, *m_track_fitter);
        }
    }

    // Create PRGraph and first segment
    auto pr_graph = std::make_shared<WireCell::Clus::PR::Graph>();
    m_track_fitter->add_graph(pr_graph);

    WireCell::Clus::PR::PatternAlgorithms pattern_algos;
    pattern_algos.m_perf = m_perf;
    pattern_algos.m_dir_weak_use_score = m_dir_weak_use_score;
    pattern_algos.m_mip_dqdx        = m_mip_dqdx / units::cm;         // e/cm -> internal
    pattern_algos.m_mip_dqdx_median = m_mip_dqdx_median / units::cm;  // e/cm -> internal
    pattern_algos.m_proton_dir_vote      = m_proton_dir_vote;
    pattern_algos.m_proton_dir_score_max = m_proton_dir_score_max;
    pattern_algos.m_proton_dir_asym_min  = m_proton_dir_asym_min;
    pattern_algos.m_endpoint_trim_retry  = m_endpoint_trim_retry;
    pattern_algos.m_fit_vertex_min_seg_length = m_fit_vertex_min_seg_length * units::cm;  // cm -> internal
    // doc sbnd_xin/docs/pr/51 round 7: robust vertex fit
    pattern_algos.m_mvfit_robust      = m_mvfit_robust;
    pattern_algos.m_mvfit_main_only   = m_mvfit_main_only;
    pattern_algos.m_mvfit_min_len     = m_mvfit_min_len * units::cm;      // cm -> internal
    pattern_algos.m_mvfit_rin_margin  = m_mvfit_rin_margin * units::cm;   // cm -> internal
    pattern_algos.m_mvfit_rout_frac   = m_mvfit_rout_frac;                // unitless
    pattern_algos.m_mvfit_rout_min    = m_mvfit_rout_min * units::cm;     // cm -> internal
    pattern_algos.m_mvfit_rout_max    = m_mvfit_rout_max * units::cm;     // cm -> internal
    pattern_algos.m_mvfit_angle       = m_mvfit_angle;                    // deg, no conversion
    pattern_algos.m_mvfit_min_pts     = m_mvfit_min_pts;                  // count
    pattern_algos.m_mvfit_min_aniso   = m_mvfit_min_aniso;                // unitless
    pattern_algos.m_mvfit_prior_range = m_mvfit_prior_range * units::cm;  // cm -> internal
    pattern_algos.m_cathode_x         = m_cathode_x * units::cm;          // cm -> internal
    pattern_algos.m_cathode_kink_xcut = m_cathode_kink_xcut * units::cm;  // cm -> internal
    pattern_algos.m_cathode_wide_kink_angle    = m_cathode_wide_kink_angle;                // deg, no conversion
    pattern_algos.m_cathode_wide_kink_skirt    = m_cathode_wide_kink_skirt * units::cm;    // cm -> internal
    pattern_algos.m_cathode_wide_kink_baseline = m_cathode_wide_kink_baseline * units::cm; // cm -> internal
    // doc sbnd_xin/docs/pr/48 -- back-to-back track fixes.
    pattern_algos.m_two_end_break       = m_two_end_break;
    pattern_algos.m_teb_min_len         = m_teb_min_len * units::cm;       // cm -> internal
    pattern_algos.m_teb_min_arm         = m_teb_min_arm * units::cm;       // cm -> internal
    pattern_algos.m_teb_min_arm_pts     = m_teb_min_arm_pts;
    pattern_algos.m_teb_stub_max        = m_teb_stub_max * units::cm;      // cm -> internal
    pattern_algos.m_teb_accept_range    = m_teb_accept_range * units::cm;  // cm -> internal
    pattern_algos.m_teb_rise_r1         = m_teb_rise_r1;
    pattern_algos.m_teb_rise_r2         = m_teb_rise_r2;
    pattern_algos.m_teb_abs_end_min     = m_teb_abs_end_min;
    pattern_algos.m_teb_dip_floor       = m_teb_dip_floor;
    pattern_algos.m_teb_score_cap_r1    = m_teb_score_cap_r1;
    pattern_algos.m_teb_score_cap_r2    = m_teb_score_cap_r2;
    pattern_algos.m_teb_turn_angle      = m_teb_turn_angle;                // deg, no conversion
    pattern_algos.m_teb_turn_baseline   = m_teb_turn_baseline * units::cm; // cm -> internal
    pattern_algos.m_teb_turn_skirt      = m_teb_turn_skirt * units::cm;    // cm -> internal
    pattern_algos.m_teb_turn_min_arm_frac = m_teb_turn_min_arm_frac;       // dimensionless, no conversion
    pattern_algos.m_teb_second_max      = m_teb_second_max * units::cm;    // cm -> internal
    pattern_algos.m_teb_chain_topology  = m_teb_chain_topology;
    pattern_algos.m_teb_r3_turn         = m_teb_r3_turn;                   // deg, no conversion
    pattern_algos.m_teb_r3_hot          = m_teb_r3_hot;                    // x mip median, no conversion
    pattern_algos.m_teb_bragg_veto_turn = m_teb_bragg_veto_turn;           // deg, no conversion
    pattern_algos.m_kink_walk_dqdx_stop = m_kink_walk_dqdx_stop;
    pattern_algos.m_kink_break_protect  = m_kink_break_protect;
    pattern_algos.m_kink_dqdx_hot_ratio = m_kink_dqdx_hot_ratio;
    // doc sbnd_xin/docs/pr/50 -- main-vertex kink-consistency snap.
    pattern_algos.m_vertex_kink_snap = m_vertex_kink_snap;
    pattern_algos.m_vks_radius       = m_vks_radius * units::cm;    // cm -> internal
    pattern_algos.m_vks_min_dis      = m_vks_min_dis * units::cm;   // cm -> internal
    pattern_algos.m_vks_angle        = m_vks_angle;                 // deg, no conversion
    pattern_algos.m_vks_margin       = m_vks_margin;                // deg, no conversion
    pattern_algos.m_vks_collinear    = m_vks_collinear;             // deg, no conversion
    pattern_algos.m_vks_skirt        = m_vks_skirt * units::cm;     // cm -> internal
    pattern_algos.m_vks_baseline     = m_vks_baseline * units::cm;  // cm -> internal
    pattern_algos.m_vks_min_arm      = m_vks_min_arm * units::cm;   // cm -> internal
    pattern_algos.m_vks_fit_miss     = m_vks_fit_miss * units::cm;  // cm -> internal
    pattern_algos.m_vks_hot_ratio    = m_vks_hot_ratio;             // x mip median, no conversion
    pattern_algos.m_vks_carry_prong  = m_vks_carry_prong * units::cm; // cm -> internal (doc pr/85)
    // doc sbnd_xin/docs/pr/51 -- main-vertex graph audit.
    pattern_algos.m_main_vertex_graph_audit = m_main_vertex_graph_audit;
    pattern_algos.m_mvga_radius       = m_mvga_radius * units::cm;    // cm -> internal
    pattern_algos.m_mvga_dup_tol      = m_mvga_dup_tol * units::cm;   // cm -> internal
    pattern_algos.m_mvga_dup_frac     = m_mvga_dup_frac;              // fraction, no conversion
    pattern_algos.m_mvga_dup_angle    = m_mvga_dup_angle;             // deg, no conversion
    pattern_algos.m_mvga_bridge_mip   = m_mvga_bridge_mip;            // x mip median, no conversion
    pattern_algos.m_mvga_reconnect    = m_mvga_reconnect * units::cm; // cm -> internal
    pattern_algos.m_mvga_stub         = m_mvga_stub * units::cm;      // cm -> internal
    pattern_algos.m_mvga_stub_pts     = m_mvga_stub_pts;
    pattern_algos.m_mvga_reseat_angle = m_mvga_reseat_angle;          // deg, no conversion
    pattern_algos.m_mvga_satellite    = m_mvga_satellite * units::cm; // cm -> internal
    pattern_algos.m_mvga_interposed   = m_mvga_interposed;            // doc pr/85
    pattern_algos.m_mvga_interposed_angle = m_mvga_interposed_angle;  // deg, no conversion
    pattern_algos.m_mvga_interposed_len = m_mvga_interposed_len * units::cm; // cm -> internal (doc pr/86)
    pattern_algos.m_mvga_sat_dup_frac = m_mvga_sat_dup_frac;          // fraction, no conversion (doc pr/86)
    pattern_algos.m_mvga_interposed_deg1 = m_mvga_interposed_deg1;    // doc pr/86
    pattern_algos.m_mvga_splice_straighten = m_mvga_splice_straighten * units::cm; // cm -> internal (doc pr/86 round 2)
    pattern_algos.m_mvga_approach_collapse = m_mvga_approach_collapse * units::cm; // cm -> internal (doc pr/86 round 2)
    pattern_algos.m_mvga_straighten_radius = m_mvga_straighten_radius * units::cm; // cm -> internal (doc pr/86 round 2)
    pattern_algos.m_rough_path_probe  = m_rough_path_probe;           // doc pr/51 round 4: diagnostic-only
    // doc pr/51 round 5: steiner gap penalty.  The two service handles are
    // unconditional copies (inert while the scale is 0).
    pattern_algos.m_steiner_gap_penalty = m_steiner_gap_penalty;
    pattern_algos.m_sgp_dead_alpha      = m_sgp_dead_alpha;                    // fraction, no conversion
    pattern_algos.m_sgp_min_edge        = m_sgp_min_edge * units::cm;          // cm -> internal
    pattern_algos.m_sgp_sample_step     = m_sgp_sample_step * units::cm;       // cm -> internal
    pattern_algos.m_sgp_point_radius    = m_sgp_point_radius * units::cm;      // cm -> internal
    pattern_algos.m_sgp_edge_probe      = m_sgp_edge_probe;                    // doc pr/73: diagnostic-only
    pattern_algos.m_vertex_scoreboard   = m_vertex_scoreboard;                 // doc pr/75: diagnostic-only
    // doc pr/79 §10: the conjunction, so fill sites may assume the board is active.
    pattern_algos.m_vtx_harvest         = m_vertex_scoreboard && m_dl_vtx_harvest;
    // doc pr/51 round 6: weak-charge deficit term (charge units, no conversion).
    pattern_algos.m_sgp_weak_scale      = m_sgp_weak_scale;
    pattern_algos.m_sgp_weak_qref       = m_sgp_weak_qref;
    // doc pr/73 round 2 F3a: a LENGTH, so it takes the cm->internal conversion.
    // -1 * units::cm stays negative, so the `< 0` off-test survives it.
    pattern_algos.m_sgp_max_sep         = m_sgp_max_sep * units::cm;           // cm -> internal
    // doc pr/83: oriented break_segment splits (bool, no conversion).
    pattern_algos.m_break_seg_orient    = m_break_seg_orient;
    pattern_algos.m_sgp_dv   = m_dv;
    pattern_algos.m_sgp_pcts = m_pcts;
    pattern_algos.m_shower_topo_demote_len = m_shower_topo_demote_len * units::cm;  // cm -> internal
    // doc sbnd_xin/docs/pr/30 §11 port-fidelity knobs.
    pattern_algos.m_fit_exclusion            = m_fit_exclusion;
    pattern_algos.m_graph_endpoint_strict    = m_graph_endpoint_strict;
    pattern_algos.m_graph_endpoint_tol       = m_graph_endpoint_tol * units::cm;   // cm -> internal
    pattern_algos.m_oov_prototype_parity     = m_oov_prototype_parity;
    pattern_algos.m_first_seg_local_pca      = m_first_seg_local_pca;
    pattern_algos.m_other_seg_relaxed_accept = m_other_seg_relaxed_accept;
    // doc sbnd_xin/docs/pr/45.
    pattern_algos.m_other_seg_empty_2d_guard = m_other_seg_empty_2d_guard;
    // doc sbnd_xin/docs/pr/54.
    pattern_algos.m_other_seg_keep_isolated            = m_other_seg_keep_isolated;
    pattern_algos.m_other_seg_keep_isolated_min_points = m_other_seg_keep_isolated_min_points;
    pattern_algos.m_other_seg_keep_isolated_min_length = m_other_seg_keep_isolated_min_length * units::cm; // cm -> internal
    // doc sbnd_xin/docs/pr/67 round 3.
    pattern_algos.m_iso_snap_min_dir_mag = m_iso_snap_min_dir_mag * units::cm; // cm -> internal
    // doc sbnd_xin/docs/pr/59 round 2.
    pattern_algos.m_assoc_full_recluster = m_assoc_full_recluster;
    // doc sbnd_xin/docs/pr/64 round 7.
    pattern_algos.m_assoc_reassign_orphans = m_assoc_reassign_orphans;
    // doc sbnd_xin/docs/pr/64 round 8.
    pattern_algos.m_assoc_clear_on_merge = m_assoc_clear_on_merge;
    // doc sbnd_xin/docs/pr/31 §11 (F2).
    pattern_algos.m_shower_topo_proto_dir    = m_shower_topo_proto_dir;
    // doc sbnd_xin/docs/pr/32 §11 (F1-F4).
    pattern_algos.m_vertex_dir_use_fit_point       = m_vertex_dir_use_fit_point;
    pattern_algos.m_shower_traj_recheck_parity     = m_shower_traj_recheck_parity;
    pattern_algos.m_main_vertex_require_descriptor = m_main_vertex_require_descriptor;
    pattern_algos.m_main_vertex_candidate_flag     = m_main_vertex_candidate_flag;
    // doc sbnd_xin/docs/pr/31 §12 (the §10.12 round: F5, F6, F3, F1, F4, F7).
    pattern_algos.m_cont_muon_dir3_30cm             = m_cont_muon_dir3_30cm;
    pattern_algos.m_track_comp_empty_abstain        = m_track_comp_empty_abstain;
    pattern_algos.m_shower_topo_reset               = m_shower_topo_reset;
    pattern_algos.m_reclass_preserve_4mom           = m_reclass_preserve_4mom;
    pattern_algos.m_reclass_never_computed_ke_floor = m_reclass_never_computed_ke_floor; // doc pr/40 round 2 F6
    pattern_algos.m_dir_track_median_local          = m_dir_track_median_local;
    pattern_algos.m_examine_showers_vertex_by_index = m_examine_showers_vertex_by_index;
    // segment_is_shower_trajectory is a free function reached from three files
    // with no component config in scope, so F2's flag-refresh half travels
    // through a process-wide flag written here, once, before any graph is
    // built -- same transport as PR::g_graph_endpoint_policy below.
    PR::g_shower_traj_refresh_flag = m_shower_traj_recheck_parity;
    // add_segment() is a free function with no component config in reach, so
    // the P8 policy travels through a process-wide struct written here, once,
    // before any graph is built (doc pr/30 §11 P8).
    WireCell::Clus::PR::g_graph_endpoint_policy.strict = m_graph_endpoint_strict;
    WireCell::Clus::PR::g_graph_endpoint_policy.tol    = m_graph_endpoint_tol * units::cm;
    pattern_algos.m_iso_endpoint               = m_iso_endpoint;
    pattern_algos.m_iso_endpoint_min_length    = m_iso_endpoint_min_length * units::cm;  // cm -> internal
    pattern_algos.m_iso_endpoint_max_xext      = m_iso_endpoint_max_xext * units::cm;    // cm -> internal
    pattern_algos.m_iso_endpoint_xext_frac     = m_iso_endpoint_xext_frac;
    pattern_algos.m_iso_endpoint_xext_quantile = m_iso_endpoint_xext_quantile;
    pattern_algos.m_iso_endpoint_tube_radius   = m_iso_endpoint_tube_radius * units::cm;  // cm -> internal
    pattern_algos.m_iso_endpoint_min_aspect    = m_iso_endpoint_min_aspect;
    pattern_algos.m_traj_cover_probe           = m_traj_cover_probe;
    // doc sbnd_xin/docs/pr/72 round 2 -- examine_structure_3 stub guard.
    pattern_algos.m_es3_stub_guard             = m_es3_stub_guard;
    pattern_algos.m_es3sg_stub_max             = m_es3sg_stub_max * units::cm;  // cm -> internal
    pattern_algos.m_es3sg_len_ratio            = m_es3sg_len_ratio;
    pattern_algos.m_es3sg_ang3_min             = m_es3sg_ang3_min;
    pattern_algos.m_es3sg_ang_ratio            = m_es3sg_ang_ratio;
    pattern_algos.m_es3sg_require_terminal     = m_es3sg_require_terminal;
    // doc pr/67 P6: remove_segment() is a free function, so the knob is mirrored
    // into a file-static in PRGraph.cxx rather than read from pattern_algos.
    PR::set_traj_cover_probe(m_traj_cover_probe);
    pattern_algos.m_pr_find_other_rounds       = m_pr_find_other_rounds;
    pattern_algos.m_v3_extension_guard         = m_v3_extension_guard;
    pattern_algos.m_v3_extension_min_gain      = m_v3_extension_min_gain * units::cm;    // cm -> internal
    // Detector-extent literals, cm -> internal (docs/pr/2 sec. 2e(iv)).
    pattern_algos.m_cosmic_y_top_main    = m_cosmic_y_top_main    * units::cm;
    pattern_algos.m_cosmic_y_top_strict  = m_cosmic_y_top_strict  * units::cm;
    pattern_algos.m_cosmic_y_top_loose   = m_cosmic_y_top_loose   * units::cm;
    pattern_algos.m_cosmic_y_small_piece = m_cosmic_y_small_piece * units::cm;
    pattern_algos.m_vertex_z_prior_scale = m_vertex_z_prior_scale * units::cm;
    // Dimensionless directions -- no unit conversion (unlike the dQ/dx scales).
    pattern_algos.m_ssm_target_dir   = WireCell::Vector(m_ssm_target_dir[0],   m_ssm_target_dir[1],   m_ssm_target_dir[2]);
    pattern_algos.m_ssm_absorber_dir = WireCell::Vector(m_ssm_absorber_dir[0], m_ssm_absorber_dir[1], m_ssm_absorber_dir[2]);
    // Charge -> energy calibration.  All dimensionless except w_value, which is
    // consumed as eV inside kine_charge_from_maps -- no unit conversion here.
    pattern_algos.m_kine_charge.fudge_factor        = m_kine_fudge_factor;
    pattern_algos.m_kine_charge.recom_factor        = m_kine_recom_factor;
    pattern_algos.m_kine_charge.shower_fudge_factor = m_kine_shower_fudge_factor;
    pattern_algos.m_kine_charge.shower_recom_factor = m_kine_shower_recom_factor;
    pattern_algos.m_kine_charge.proton_recom_factor = m_kine_proton_recom_factor;
    pattern_algos.m_kine_charge.plane_weights       = {m_kine_plane_weights[0], m_kine_plane_weights[1], m_kine_plane_weights[2]};
    pattern_algos.m_kine_charge.plane_asym_switch   = m_kine_plane_asym_switch;
    pattern_algos.m_kine_charge.shower_pdg_live     = m_kine_shower_pdg_live;
    pattern_algos.m_kine_charge.w_value             = m_kine_w_value;
    // doc sbnd_xin/docs/pr/36 §10 tagger-stage knobs (F4/F5/F6/F7).
    pattern_algos.m_tagger_ordered_segment_sets  = m_tagger_ordered_segment_sets;
    pattern_algos.m_stem_endpoint_wcpt_parity    = m_stem_endpoint_wcpt_parity;
    pattern_algos.m_broken_muon_cluster_id_count = m_broken_muon_cluster_id_count;
    pattern_algos.m_neutrino_type_bitmask        = m_neutrino_type_bitmask;
    // doc sbnd_xin/docs/pr/33 §10 EM-shower-clustering knobs.
    pattern_algos.m_daughter_count_proto_main_vertex     = m_daughter_count_proto_main_vertex;
    pattern_algos.m_daughter_count_proto_examine_showers = m_daughter_count_proto_examine_showers;
    pattern_algos.m_shower_pdg_from_start_segment        = m_shower_pdg_from_start_segment;
    pattern_algos.m_shower_pdg_from_shower_type          = m_shower_pdg_from_shower_type;
    pattern_algos.m_shower_pdg_exact_muon_test           = m_shower_pdg_exact_muon_test;
    pattern_algos.m_pi0_id_shared_allocator              = m_pi0_id_shared_allocator;
    pattern_algos.m_shower_flag_pdg_electron             = m_shower_flag_pdg_electron;
    pattern_algos.m_shower_less_id_tiebreak              = m_shower_less_id_tiebreak;
    pattern_algos.m_shower_endpoint_exclude_start_vertex = m_shower_endpoint_exclude_start_vertex;
    // doc sbnd_xin/docs/pr/40 -- track mis-identified as electron.
    pattern_algos.m_track_pid_persist_dqdx    = m_track_pid_persist_dqdx;    // F1: threaded via track_pid_options()
    pattern_algos.m_shower_reclass_dqdx_guard = m_shower_reclass_dqdx_guard; // F2
    pattern_algos.m_shower_topo_dqdx_guard    = m_shower_topo_dqdx_guard;    // F3
    // doc sbnd_xin/docs/pr/40 round 2 -- two follow-on defects from the pr/40 round.
    pattern_algos.m_track_pid_persist_4mom      = m_track_pid_persist_4mom;      // F4: threaded via track_pid_options()
    pattern_algos.m_shower_proton_daughter_pion = m_shower_proton_daughter_pion; // F5
    pattern_algos.m_shower_proton_daughter_pion_dissolve = m_shower_proton_daughter_pion_dissolve; // F7
    pattern_algos.m_muon_multi_proton_pion               = m_muon_multi_proton_pion;               // F8
    pattern_algos.m_track_pid_persist_dqdx_electron_guard     = m_track_pid_persist_dqdx_electron_guard;     // F9
    pattern_algos.m_shower_connect_main_vertex_straight_guard = m_shower_connect_main_vertex_straight_guard; // F10
    pattern_algos.m_shower_traj_straight_guard                = m_shower_traj_straight_guard;                // F11
    pattern_algos.m_shower_absorb_track_guard                 = m_shower_absorb_track_guard;                 // F12
    pattern_algos.m_shower_absorb_unreachable_main            = m_shower_absorb_unreachable_main;            // doc pr/65 round 3
    pattern_algos.m_shower_connect_protected_pion_guard       = m_shower_connect_protected_pion_guard;       // F13
    pattern_algos.m_michel_stem_muon_rescue                   = m_michel_stem_muon_rescue;                   // F14
    pattern_algos.m_shower_in_cascade_guard                   = m_shower_in_cascade_guard;                   // pr/74 P1
    pattern_algos.m_shower_in_max_len                         = m_shower_in_max_len * units::cm;             // pr/74 P1
    pattern_algos.m_shower_in_mip_hi                          = m_shower_in_mip_hi;                          // pr/74 P1
    pattern_algos.m_michel_stem_michel_check                  = m_michel_stem_michel_check;                  // pr/74 P2
    pattern_algos.m_michel_stem_max_far_len                   = m_michel_stem_max_far_len * units::cm;       // pr/74 P2
    pattern_algos.m_shower_stem_backfill                      = m_shower_stem_backfill;                      // pr/74 K4
    pattern_algos.m_stem_backfill_max_len                     = m_stem_backfill_max_len * units::cm;         // pr/74 K4
    pattern_algos.m_stem_backfill_mip_lo                      = m_stem_backfill_mip_lo;                      // pr/74 K4
    pattern_algos.m_stem_backfill_mip_hi                      = m_stem_backfill_mip_hi;                      // pr/74 K4
    pattern_algos.m_stem_backfill_min_shower_len              = m_stem_backfill_min_shower_len * units::cm;  // pr/74 K4
    pattern_algos.m_shower_conn3_unreachable                  = m_shower_conn3_unreachable;                  // pr/74 K5
    pattern_algos.m_conn3_unreachable_min_len                 = m_conn3_unreachable_min_len * units::cm;     // pr/74 K5
    pattern_algos.m_shower_traj_michel_stem                   = m_shower_traj_michel_stem;                   // pr/74 K6
    pattern_algos.m_michel_stem_traj_min_len                  = m_michel_stem_traj_min_len * units::cm;      // pr/74 K6
    pattern_algos.m_michel_stem_traj_max_len                  = m_michel_stem_traj_max_len * units::cm;      // pr/74 K6
    pattern_algos.m_michel_stem_traj_mip_lo                   = m_michel_stem_traj_mip_lo;                   // pr/74 K6 (dimensionless ratio)
    pattern_algos.m_michel_stem_traj_max_far_len              = m_michel_stem_traj_max_far_len * units::cm;  // pr/74 K6
    pattern_algos.m_michel_stem_traj_min_kink_deg             = m_michel_stem_traj_min_kink_deg;             // pr/74 K6 (degrees)
    pattern_algos.m_shower_long_muon_keep_type                = m_shower_long_muon_keep_type;                // doc pr/44
    pattern_algos.m_single_muon_proton_chain_veto             = m_single_muon_proton_chain_veto;             // doc pr/43 round 2 K1
    pattern_algos.m_single_muon_long_muon_claim               = m_single_muon_long_muon_claim;               // doc pr/43 round 2 K2
    pattern_algos.m_pid_flag_reconcile                        = m_pid_flag_reconcile;                        // doc pr/43 round 2 K3
    pattern_algos.m_long_muon_stub_bridge                     = m_long_muon_stub_bridge;                     // doc pr/46
    // Muon dQ/dx-vs-length envelope: c0/c1/power dimensionless, pivot cm -> internal.
    pattern_algos.m_muon_dqdx_curve = {m_muon_dqdx_curve[0], m_muon_dqdx_curve[1],
                                       m_muon_dqdx_curve[2] * units::cm, m_muon_dqdx_curve[3]};
    // Single-photon stem dE/dx conversion; the cut narrows to float so the
    // default compares bit-identically to the legacy 2.3f literal.
    pattern_algos.m_sp_dedx_use_recomb_model = m_sp_dedx_use_recomb_model;
    pattern_algos.m_sp_mean_dedx_cut         = static_cast<float>(m_sp_mean_dedx_cut);
    pattern_algos.m_recomb_model             = m_recomb_model;
    m_track_fitter->set_perf(m_perf);
    // doc pr/49: thread the own-blob-coverage knob into the fitter (double
    // sentinel; -1 default matches TrackFitting::Parameters, so this is a
    // no-op unless the config opts in).
    m_track_fitter->set_parameter("fit_blob_coverage", m_fit_blob_coverage);
    // doc pr/67 P3 (log-only end-trim probe); 0 = off = byte-identical.
    m_track_fitter->set_parameter("traj_cover_probe", m_traj_cover_probe ? 1.0 : 0.0);
    // doc sbnd_xin/docs/pr/50 (fit_blob_coverage_defer, default false):
    // wrap the MAIN cluster's find_proto_vertex call so its recursive
    // break partition forms on legacy (undeweighted) fits -- the partition
    // is globally sensitive to fit perturbations (172230: 200 deweight
    // firings ~90 cm away reshuffled 34->33 segments and lost the
    // true-kink main-vertex candidate; same class measured in 131357 /
    // 342199 / 360535 / 469665).  Main cluster ONLY: its later stages
    // (determine_main_vertex, improve_vertex, the final trajectory +
    // dQ/dx) all refit with the restored deweighting, so ghost protection
    // survives -- but a non-main cluster's final trajectory is essentially
    // its find_proto_vertex fit, so deferring there UN-fixes the pr/49
    // ghosts (57441 cid 20 measured 1.12 -> 1.23 cm under a global defer).
    // Local fitters spawned inside (inherit_from) copy the suspended
    // value, so the whole stage is covered.  Both lambdas are no-ops
    // unless the defer knob AND the base knob are on.
    const bool cov_defer_active = m_fit_blob_coverage_defer && m_fit_blob_coverage >= 0;
    auto cov_defer_suspend = [&]() {
        if (cov_defer_active) m_track_fitter->set_parameter("fit_blob_coverage", -1);
    };
    auto cov_defer_restore = [&]() {
        if (cov_defer_active) m_track_fitter->set_parameter("fit_blob_coverage", m_fit_blob_coverage);
    };
    if (cov_defer_active) {
        SPDLOG_LOGGER_DEBUG(log, "fit_blob_coverage_defer on: partition stage (find_proto_vertex) runs legacy fits, deweight restored for all later stages");
    }

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

    VertexPtr final_main_vertex = nullptr;
    bool flag_dl_changed = false;

    // doc pr/59: diagnostic-only, env-gated (WCT_PR59_ASSOC_CENSUS unset =>
    // no log lines, no behavior change) sentinel naming which cluster each
    // clustering_points (associate_points) pass actually ran against.
    // main_cluster can be silently repointed later (determine_overall_
    // main_vertex_DL holds it by reference and may call swap_main_cluster),
    // which would leave whichever cluster the first pass ran on NEVER
    // re-associated -- see the second call site further below.
    static const bool pr59_assoc_census = std::getenv("WCT_PR59_ASSOC_CENSUS") != nullptr;

    {
        // initial pattern recognitions
        // particle_data (doc pr/48): stopping templates for the two-end
        // break pass; inert unless two_end_break is on.
        cov_defer_suspend();
        pattern_algos.find_proto_vertex(*pr_graph, *main_cluster, *m_track_fitter, m_dv, true, 2, true, particle_data());
        cov_defer_restore();
        detg_dump("main:find_proto_vertex", *pr_graph);

        // shower related operations
        pattern_algos.clustering_points(*pr_graph, *main_cluster, m_dv);
        if (pr59_assoc_census) {
            SPDLOG_LOGGER_DEBUG(log,
                "pr59 assoc-census: first clustering_points call, main_cluster={}",
                main_cluster->get_cluster_id());
        }
        detg_dump("main:clustering_points", *pr_graph);
        pattern_algos.separate_track_shower(*pr_graph, *main_cluster);
        detg_dump("main:separate_track_shower", *pr_graph);

        // direction determination
        pattern_algos.determine_direction(*pr_graph, *main_cluster, particle_data(), m_recomb_model);
        detg_dump("main:determine_direction", *pr_graph);

        // shower clustering
        pattern_algos.shower_determining_in_main_cluster(*pr_graph, *main_cluster, particle_data(), m_recomb_model, m_dv);
        detg_dump("main:shower_determining", *pr_graph);

        // main vertex determination
        pattern_algos.determine_main_vertex(*pr_graph, *main_cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, *m_track_fitter, m_dv, particle_data(), m_recomb_model);
        detg_dump("main:determine_main_vertex", *pr_graph);

        // doc sbnd_xin/docs/pr/59 round 2 (P1): determine_main_vertex's
        // internal examine_structure_final*/examine_vertices_1 can create a
        // brand-new segment (18255-142421 seg 20; 116944-71372 segs
        // 19052/19053/136199, all confirmed via WCT_DET_DEBUG=2 backtraces).
        // Rescue it here, before determine_direction/shower_determining_in_
        // main_cluster/deghosting/shower_clustering_with_nv all consume its
        // (until now, or still, missing) associate_points and shower flags.
        // No-op unless m_assoc_full_recluster.
        pattern_algos.reassociate_cluster_orphans(*pr_graph, *main_cluster, m_dv);

        if (main_vertex !=nullptr){
            map_cluster_main_vertices[main_cluster] = main_vertex;
            main_vertex = nullptr;
        }

        std::cout << "After first round of main cluster PR" << std::endl;        pattern_algos.print_segs_info(*pr_graph, *main_cluster, 0);
    }
    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: main_cluster initial PR took {} ms", MS(Clock::now() - t0).count());
    t0 = Clock::now();

    // Loop over other (non-main) beam-flash clusters
    if (!other_clusters.empty()) {
        for (auto* cluster : other_clusters) {
            if (cluster->get_length() > 6 * units::cm) {
                // std::cout << "Long Cluster " << cluster->get_cluster_id() << " " << cluster->nchildren() << std::endl;
                // Long cluster: break tracks and do 2 rounds of other-track finding
                pattern_algos.find_proto_vertex(*pr_graph, *cluster, *m_track_fitter, m_dv, true, 2, false);
                pattern_algos.clustering_points(*pr_graph, *cluster, m_dv);
                pattern_algos.separate_track_shower(*pr_graph, *cluster);
                pattern_algos.determine_direction(*pr_graph, *cluster, particle_data(), m_recomb_model);
                pattern_algos.shower_determining_in_main_cluster(*pr_graph, *cluster, particle_data(), m_recomb_model, m_dv);
                pattern_algos.determine_main_vertex(*pr_graph, *cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, *m_track_fitter, m_dv, particle_data(), m_recomb_model);
                // doc sbnd_xin/docs/pr/59 round 2 (P1), other-cluster branch.
                pattern_algos.reassociate_cluster_orphans(*pr_graph, *cluster, m_dv);
                if (main_vertex != nullptr) {
                    map_cluster_main_vertices[cluster] = main_vertex;
                    main_vertex = nullptr;
                }
                detg_dump("other:long_cluster", *pr_graph);
            } else {
                // Short cluster: no track breaking, 1 round; fall back to init_point_segment if needed
                if (!pattern_algos.find_proto_vertex(*pr_graph, *cluster, *m_track_fitter, m_dv, false, 1, false)) {
                    // std::cout << "Point Cluster " << cluster->get_cluster_id() << " " << cluster->nchildren() <<std::endl;
                    pattern_algos.init_point_segment(*pr_graph, *cluster, *m_track_fitter, m_dv);
                }
                pattern_algos.clustering_points(*pr_graph, *cluster, m_dv);
                pattern_algos.separate_track_shower(*pr_graph, *cluster);
                pattern_algos.determine_direction(*pr_graph, *cluster, particle_data(), m_recomb_model);
                pattern_algos.shower_determining_in_main_cluster(*pr_graph, *cluster, particle_data(), m_recomb_model, m_dv);
                pattern_algos.determine_main_vertex(*pr_graph, *cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, *m_track_fitter, m_dv, particle_data(), m_recomb_model);
                // doc sbnd_xin/docs/pr/59 round 2 (P1), other-cluster branch.
                pattern_algos.reassociate_cluster_orphans(*pr_graph, *cluster, m_dv);
                if (main_vertex != nullptr) {
                    map_cluster_main_vertices[cluster] = main_vertex;
                    main_vertex = nullptr;
                }
                detg_dump("other:short_cluster", *pr_graph);
            }
        }
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: other_clusters PR took {} ms", MS(Clock::now() - t0).count());
        t0 = Clock::now();

        // Deghost across all beam-flash clusters (main + others)
        std::vector<Cluster*> all_clusters;
        all_clusters.push_back(main_cluster);
        all_clusters.insert(all_clusters.end(), other_clusters.begin(), other_clusters.end());

        pattern_algos.deghosting(*pr_graph, map_cluster_main_vertices, all_clusters, *m_track_fitter, m_dv);
        detg_dump("deghosting", *pr_graph);
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: deghosting took {} ms", MS(Clock::now() - t0).count());
        t0 = Clock::now();
    }

    // Determine the overall neutrino vertex.
    // If DL weights are configured, try DL first (matches prototype flag_dl_vtx logic).
    // Fall back to traditional algorithm if DL is disabled or does not change the vertex.
    // DL path updates map_cluster_main_vertices[main_cluster] directly (by-ref parameter).
    // Traditional path returns the chosen vertex; capture it and sync to the map.
 
    if (!m_dl_weights.empty()) {
        flag_dl_changed = pattern_algos.determine_overall_main_vertex_DL(
            *pr_graph, map_cluster_main_vertices, main_cluster, other_clusters,
            vertices_in_long_muon, segments_in_long_muon,
            *m_track_fitter, m_dv, particle_data(), m_recomb_model,
            m_dl_weights, m_dl_vtx_cut, m_dQdx_scale, m_dQdx_offset,
            m_dl_vtx_rerank, m_dl_vtx_top_k, m_dl_vtx_min_accept_score,
            m_dl_vtx_score_scale, m_dl_vtx_swap_guard,
            m_dl_vtx_topo_weight, m_dl_vtx_topo_center);
    }
    if (!flag_dl_changed) {
        // doc sbnd_xin/docs/pr/51 round 3: determine_overall_main_vertex now
        // takes the map and main_cluster by reference, so an internal cluster
        // swap (examine_main_vertices / check_switch_main_cluster[_2] ->
        // swap_main_cluster) is visible here instead of silently discarded.
        // Pass throwaway local copies -- map_copy/mc_copy end up holding the
        // function's true post-swap state, while this scope's own
        // map_cluster_main_vertices/main_cluster are untouched unless
        // m_main_vertex_swap_apply says to sync them.  Knob off => mc_copy is
        // read, compared, then discarded => byte-identical to the pre-round-3
        // by-value behaviour.
        ClusterVertexMap map_copy = map_cluster_main_vertices;
        Cluster* mc_copy = main_cluster;
        final_main_vertex = pattern_algos.determine_overall_main_vertex(
            *pr_graph, map_copy, mc_copy, other_clusters,
            vertices_in_long_muon, segments_in_long_muon,
            *m_track_fitter, m_dv, particle_data(), m_recomb_model, true);
        if (mc_copy != main_cluster) {
            SPDLOG_LOGGER_DEBUG(log,
                "mvsa: traditional path swapped main cluster {} -> {} ({})",
                main_cluster->get_cluster_id(), mc_copy->get_cluster_id(),
                m_main_vertex_swap_apply ? "applied" : "discarded");
            if (m_main_vertex_swap_apply) {
                main_cluster = mc_copy;
                map_cluster_main_vertices = map_copy;
            }
        }
        if (final_main_vertex) {
            map_cluster_main_vertices[main_cluster] = final_main_vertex;
        }
    }

    // Retrieve the chosen neutrino vertex regardless of which path ran
    {
        auto it = map_cluster_main_vertices.find(main_cluster);
        if (it != map_cluster_main_vertices.end()) {
            final_main_vertex = it->second;
        }
    }
    detg_dump("overall_main_vertex", *pr_graph);
    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: overall main vertex took {} ms", MS(Clock::now() - t0).count());
    t0 = Clock::now();

    

  

    if (final_main_vertex) {
        // doc sbnd_xin/docs/pr/50: main-vertex kink-consistency snap --
        // inert unless vertex_kink_snap.  Runs after the overall main
        // vertex is final (either DL or fallback path) and BEFORE the
        // final improve_vertex, so the local optimizer polishes a
        // corner-anchored trajectory.
        if (pattern_algos.snap_main_vertex_to_kink(*pr_graph, *main_cluster, final_main_vertex,
                                                   *m_track_fitter, m_dv, particle_data(), m_recomb_model)) {
            map_cluster_main_vertices[main_cluster] = final_main_vertex;
            detg_dump("snap_main_vertex_to_kink", *pr_graph);
        }

        pattern_algos.improve_vertex(*pr_graph, *main_cluster, final_main_vertex,
                                     vertices_in_long_muon, segments_in_long_muon,
                                     *m_track_fitter, m_dv, particle_data(), m_recomb_model,
                                     true, true);
        // improve_vertex may update final_main_vertex pointer; sync back to map
        map_cluster_main_vertices[main_cluster] = final_main_vertex;

        std::cout << "After improve vertex:" << final_main_vertex->fit().point << std::endl; pattern_algos.print_segs_info(*pr_graph, *main_cluster, 0);

        // doc sbnd_xin/docs/pr/51: main-vertex graph audit -- inert unless
        // main_vertex_graph_audit.  Runs AFTER the final improve_vertex
        // (the micro-stubs it must absorb are created there: 142421's
        // 7081/7082, 285567's 81/82/83) and BEFORE clustering_points /
        // examine_direction, which then act on the audited graph.  May
        // re-seat final_main_vertex's position in place (never the
        // pointer), so no map re-sync is needed.
        if (pattern_algos.main_vertex_graph_audit(*pr_graph, *main_cluster, final_main_vertex,
                                                  *m_track_fitter, m_dv)) {
            std::cout << "After main vertex graph audit:" << final_main_vertex->fit().point << std::endl; pattern_algos.print_segs_info(*pr_graph, *main_cluster, 0);
            detg_dump("main_vertex_graph_audit", *pr_graph);
        }

        // doc sbnd_xin/docs/pr/51 round 4: diagnostic-only rough-path probe
        // -- inert unless rough_path_probe.  Runs right after the audit
        // block (on today's production graph when mvga is off; on the
        // audited graph when it is on) so its measurements match whichever
        // near-vertex state is actually being Bee-scanned.
        pattern_algos.rough_path_probe(*pr_graph, *main_cluster, final_main_vertex, *m_track_fitter, m_dv);

        pattern_algos.clustering_points(*pr_graph, *main_cluster, m_dv);
        // doc pr/59 sentinel 3 (second call site): if this cluster_id differs
        // from the first call's, main_cluster was swapped in between and the
        // ORIGINAL main cluster's segments (any created/modified after the
        // swap) never got a second association pass.
        if (pr59_assoc_census) {
            SPDLOG_LOGGER_DEBUG(log,
                "pr59 assoc-census: second clustering_points call, main_cluster={}",
                main_cluster->get_cluster_id());
        }

        // doc sbnd_xin/docs/pr/59 round 2 (P2): a safety net for two cases P1
        // does not reach -- (a) a segment created inside improve_vertex or
        // main_vertex_graph_audit (both ran above, between the first
        // determine_main_vertex and here), and (b) main_cluster having been
        // silently repointed by determine_overall_main_vertex[_DL]'s
        // swap_main_cluster since the first clustering_points call, which
        // otherwise leaves the ORIGINAL main cluster (now in other_clusters)
        // permanently on its first-round-only association state.  No-op
        // unless m_assoc_full_recluster; a no-op per cluster (0 rescued) when
        // that cluster has no orphan, so this is cheap on the common case
        // where P1 already caught everything (measured true for both
        // 18255-142421 and 116944-71372).  Still before shower_clustering_
        // with_nv, which is the next consumer of associate_points/shower
        // flags below.
        pattern_algos.reassociate_cluster_orphans(*pr_graph, *main_cluster, m_dv);
        for (auto* cluster : other_clusters) {
            pattern_algos.reassociate_cluster_orphans(*pr_graph, *cluster, m_dv);
        }

        std::cout << "After shower clustering :" << std::endl; pattern_algos.print_segs_info(*pr_graph, *main_cluster, 0);
 
        // examine_direction runs last and has the final word on segment orientations
        // relative to the main vertex.
        pattern_algos.examine_direction(*pr_graph, final_main_vertex, final_main_vertex,
                                        vertices_in_long_muon, segments_in_long_muon,
                                        particle_data(), m_recomb_model, true);

        SPDLOG_LOGGER_TRACE(log, "Overall main vertex cluster={}", main_cluster->get_cluster_id());
        
        std::cout << "After examine direction: " << std::endl;pattern_algos.print_segs_info(*pr_graph, *main_cluster, 0);
        detg_dump("improve_vertex", *pr_graph);
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: improve_vertex + examine_direction took {} ms", MS(Clock::now() - t0).count());
        t0 = Clock::now();

        pattern_algos.shower_clustering_with_nv(acc_segment_id, pi0_showers,
                                                map_shower_pio_id, map_pio_id_showers,
                                                map_pio_id_mass, map_pio_id_saved_pair,
                                                pio_kine,
                                                vertices_in_long_muon, segments_in_long_muon,
                                                *pr_graph, final_main_vertex, showers,
                                                main_cluster, other_clusters,
                                                map_cluster_main_vertices,
                                                map_vertex_in_shower, map_segment_in_shower,
                                                map_vertex_to_shower, used_shower_clusters,
                                                *m_track_fitter, m_dv, particle_data(),
                                                m_recomb_model);

        std::cout << "After shower clustering with NV: " << std::endl; pattern_algos.print_segs_info(*pr_graph, *main_cluster, 0);
        detg_dump("shower_clustering_with_nv", *pr_graph);
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: shower_clustering_with_nv took {} ms", MS(Clock::now() - t0).count());
        t0 = Clock::now();

        // doc sbnd_xin/docs/pr/43 round 2 K3 -- late particle-info/flag
        // reconciliation, AFTER shower_clustering_with_nv and BEFORE the
        // taggers, so tagger features, kine, Bee PF tree and PR display all
        // see one consistent labeling.  No-op unless pid_flag_reconcile.
        pattern_algos.reconcile_particle_flags(*pr_graph, final_main_vertex, showers,
                                               map_vertex_in_shower, map_segment_in_shower,
                                               map_vertex_to_shower, map_shower_pio_id,
                                               particle_data(), m_recomb_model);
    }


    // Initialize tagger features to their default values unconditionally —
    // even if no vertex was found the struct must be value-initialized.
    t0 = Clock::now();
    TaggerInfo tagger_info;
    pattern_algos.init_tagger_info(tagger_info);

    // Build the full list of beam-flash clusters (main + others) once;
    // used by cosmic_tagger and potentially other taggers.
    std::vector<Cluster*> all_clusters;
    all_clusters.push_back(main_cluster);
    all_clusters.insert(all_clusters.end(), other_clusters.begin(), other_clusters.end());

    // doc pr/36 §10.3 / §10.15a (F2 = P12): the population sweep, run FIRST.
    // The prototype's get_particle_type() coerces any shower-flagged segment
    // to PDG 11 on read (ProtoSegment.cxx:10-15); the toolkit's
    // has_particle_info() gates skip such a segment.  This one pass over the
    // graph decides whether that population EXISTS on this event; the
    // per-gate counters (f2_gate_skip, incremented inside the taggers)
    // attribute it.  Zero across the manifest => F2 is dead by construction.
    {
        auto& audit = WireCell::Clus::PR::g_pr36_audit;
        for (auto ed : WireCell::Clus::PR::ordered_edges(*pr_graph)) {
            SegmentPtr sg = (*pr_graph)[ed].segment;
            if (!sg) continue;
            audit.f2_sweep_segments.fetch_add(1, std::memory_order_relaxed);
            if (!sg->has_particle_info() &&
                (sg->flags_any(PR::SegmentFlags::kShowerTrajectory) ||
                 sg->flags_any(PR::SegmentFlags::kShowerTopology)))
                audit.f2_sweep_hits.fetch_add(1, std::memory_order_relaxed);
        }
    }

    // Run cosmic and numu taggers to fill BDT input features in tagger_info.
    // Both require a valid neutrino vertex to have been found.
    if (final_main_vertex) {
        pattern_algos.cosmic_tagger(*pr_graph, final_main_vertex,
                                    showers,
                                    map_segment_in_shower,
                                    map_vertex_to_shower,
                                    segments_in_long_muon,
                                    main_cluster,
                                    all_clusters,
                                    m_dv,
                                    tagger_info);

        auto [flag_long_muon, muon_length] =
            pattern_algos.numu_tagger(*pr_graph, final_main_vertex,
                                      showers,
                                      segments_in_long_muon,
                                      main_cluster,
                                      tagger_info);
        (void)flag_long_muon;  // result stored in tagger_info.numu_cc_flag

        pattern_algos.ssm_tagger(*pr_graph, final_main_vertex,
                                 showers,
                                 map_vertex_in_shower,
                                 map_segment_in_shower,
                                 pio_kine,
                                 /*flag_ssmsp=*/-1,
                                 acc_segment_id,
                                 particle_data(),
                                 m_recomb_model,
                                 tagger_info);

        // Derive apa/face from the main-vertex position instead of assuming the
        // single-drift-volume (0,0) geometry.  nue_tagger uses these for
        // gap_identification's check_direction flags (wire_angles + drift sign
        // => flag_prolong_u/v/w, flag_parallel; NeutrinoTaggerNuE.cxx:2725-2726)
        // and for mip_quality's 2-D closest-distance queries (:1638).  In a
        // multi-drift-volume detector such as SBND a vertex in APA 1 would
        // otherwise be evaluated against APA 0's mirrored wire angles and
        // opposite drift direction.  Same derivation as
        // PatternAlgorithms::singlephoton_tagger (NeutrinoTaggerSinglePhoton.cxx).
        int nue_apa = 0, nue_face = 0;
        if (m_dv) {
            const Point nue_vtx_pt = final_main_vertex->fit().valid()
                                     ? final_main_vertex->fit().point
                                     : final_main_vertex->wcpt().point;
            const auto nue_wpid = m_dv->contained_by(nue_vtx_pt);
            // Uncontained points give apa()==-1 (and face()==1); keep the legacy
            // (0,0) rather than letting Grouping::wire_angles().at(-1) throw.
            // Same guard idiom as NeutrinoTaggerNuE.cxx:2764.
            if (nue_wpid.apa() >= 0) {
                nue_apa  = nue_wpid.apa();
                nue_face = nue_wpid.face();
            }
        }
        SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino nue_tagger volume: apa={} face={}",
                            nue_apa, nue_face);

        pattern_algos.nue_tagger(*pr_graph, main_cluster, final_main_vertex,
                                 nue_apa, nue_face,
                                 showers, map_vertex_to_shower,
                                 pi0_showers, map_shower_pio_id,
                                 map_pio_id_showers, map_pio_id_mass,
                                 m_dv, particle_data(),
                                 muon_length, tagger_info);

        // prototype (NeutrinoID.cxx lines 269-271):
        //     bool flag_sp = singlephoton_tagger(results.second);
        //     if (flag_sp){tagger_info.photon_flag = true;}
        // The port ran the tagger -- filling its ~90 shw_sp_* BDT features --
        // but dropped the verdict on the floor, leaving photon_flag at the 0
        // init_tagger_info() gives it (doc sbnd_xin/docs/pr/26 sec. 8.2).
        // C++ default false = that legacy behaviour, so the uBooNE tagger
        // ntuple's photon_flag branch is byte-identical when the knob is off.
        const bool flag_sp =
            pattern_algos.singlephoton_tagger(*pr_graph, main_cluster,
                                              final_main_vertex,
                                              showers,
                                              map_vertex_to_shower,
                                              map_shower_pio_id,
                                              map_pio_id_showers,
                                              map_pio_id_mass,
                                              m_dv,
                                              // doc pr/36 §10.4 (F3): SCE helper,
                                              // gated separately from kine's use of
                                              // clus_geom_helper.  Off (or helper
                                              // unconfigured) => nullptr => raw
                                              // positions, byte-identical legacy.
                                              m_sp_sce_correction ? m_geom_helper : nullptr,
                                              tagger_info);
        if (m_sp_photon_flag) {
            // Logged unconditionally so a knob-on run proves the branch
            // executed even on a sample where nothing is tagged.
            SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino sp_photon_flag: singlephoton_tagger returned {}",
                                flag_sp);
            if (flag_sp) tagger_info.photon_flag = 1.0f;
        }
        else {
            (void)flag_sp;  // legacy: verdict computed and discarded
        }
    }

    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: taggers took {} ms", MS(Clock::now() - t0).count());
    t0 = Clock::now();

    // Compute match_isFC: 1 if the main cluster is fully contained inside the
    // fiducial volume, 0 otherwise.  Uses the same two-round boundary check as
    // TaggerCheckSTM so the definition is consistent across both users.
    // doc pr/36 §10.2 (F1 = P1): when a "fiducial" is configured the direct
    // containment tests run against it (the TGM/FC/STM volume) instead of the
    // historical FiducialUtils fallback -- cluster_fc_check's nullptr path is
    // documented bit-for-bit (Clustering_Util.cxx:108-116), so an absent key
    // is byte-identical.  The prototype reads this verdict from the upstream
    // light-matching stage (NeutrinoID.cxx:62, branch T_eval), computed on
    // the ONE ToyFiducial volume -- so a consistent volume here is the parity
    // reading that is implementable as a knob (reading (ii), threading the
    // upstream verdict itself, is a data-flow change out of this round's
    // scope).  Both verdicts are computed when the knob is on; disagreement
    // is the §7.1 diagnostic, logged and counted.
    if (main_cluster) {
        auto fc_result = Facade::cluster_fc_check(*main_cluster, m_dv,
                                                  m_use_fiducial ? m_fiducial : nullptr,
                                                  m_fv_tolerance);
        if (m_use_fiducial) {
            auto& audit = WireCell::Clus::PR::g_pr36_audit;
            audit.f1_fc_checks.fetch_add(1, std::memory_order_relaxed);
            auto legacy = Facade::cluster_fc_check(*main_cluster, m_dv);
            if (legacy.is_fc != fc_result.is_fc) {
                audit.f1_fc_disagree.fetch_add(1, std::memory_order_relaxed);
                SPDLOG_LOGGER_INFO(log, "PR36AUDIT match_isFC disagree: legacy={} fiducial={}",
                                   legacy.is_fc, fc_result.is_fc);
            }
        }
        tagger_info.match_isFC = fc_result.is_fc ? 1.0f : 0.0f;
    }
    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: fc_check took {} ms", MS(Clock::now() - t0).count());
    t0 = Clock::now();

    // Fill reconstructed neutrino kinematics if a vertex was found.
    KineInfo kine_info{};
    if (final_main_vertex) {
        kine_info = pattern_algos.fill_kine_tree(
            final_main_vertex, showers, pio_kine,
            *pr_graph, *m_track_fitter, m_dv,
            m_geom_helper,          // nullptr when clus_geom_helper is not configured
            particle_data(), m_recomb_model);
    }
    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: fill_kine_tree took {} ms", MS(Clock::now() - t0).count());
    t0 = Clock::now(); // finalize block

    // Mark the main neutrino vertex and store neutrino results in TrackFitting
    // so that downstream consumers (e.g., Bee particle-flow output in MultiAlgBlobClustering)
    // can access them without re-running pattern recognition.
    if (final_main_vertex) {
        final_main_vertex->set_flags(PR::VertexFlags::kNeutrinoVertex);
    }
    m_track_fitter->set_pi0_data(pi0_showers, map_shower_pio_id, map_pio_id_showers, map_pio_id_mass);
    m_track_fitter->set_main_vertex(final_main_vertex);
    m_track_fitter->set_showers(showers);
    m_track_fitter->set_kine_info(kine_info);
    m_track_fitter->set_tagger_info(tagger_info);

    // doc sbnd_xin/docs/pr/75 -- stash the vertex scoreboard for PrDisplayDump.
    // Here and not earlier on purpose: this point is AFTER pr/50's
    // snap_main_vertex_to_kink and the final improve_vertex, so final_vertex_id
    // is the vertex the display actually draws, and a mismatch against
    // main_vertex in the dump means the stash moved.
    if (m_vertex_scoreboard) {
        auto& board = pattern_algos.m_vtx_board;
        board.filled = true;
        board.harvest = pattern_algos.m_vtx_harvest;  // doc pr/79 §10: serializer's emission gate
        board.weights_missing = m_dl_weights_missing;
        if (!board.dl_ran) board.route = "dl-not-run";
        if (final_main_vertex) {
            const auto pt = final_main_vertex->fit().valid()
                          ? final_main_vertex->fit().point : final_main_vertex->wcpt().point;
            const auto* cl = final_main_vertex->cluster();
            board.final_vertex_id = (cl ? cl->get_cluster_id() : 0) * 1000
                                  + static_cast<int>(final_main_vertex->get_graph_index());
            board.final_x = pt.x() / units::cm;
            board.final_y = pt.y() / units::cm;
            board.final_z = pt.z() / units::cm;
        }
        m_track_fitter->set_vertex_scoreboard(board);
    }

    // Merge every per-cluster fill_fitted_charge_2d snapshot into the flat
    // map that UbooneMagnifyTrackingVisitor::write_proj_data reads, so that
    // T_proj_data contains cells for all beam-flash clusters, not just the
    // last cluster fit by pattern recognition.
    m_track_fitter->assemble_fitted_charge_2d();

    // Store TrackFitting in the grouping for later access by bee output and tracking sink
    grouping.set_track_fitting(m_track_fitter);
    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: finalize took {} ms", MS(Clock::now() - t0).count());
    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: visit() TOTAL took {} ms", MS(Clock::now() - t_total).count());

    // doc sbnd_xin/docs/pr/30 §11 -- one machine-readable line per event.
    // The counters are process-wide and every SBND runner drives one event
    // per wire-cell process, so these cumulative totals ARE this event's.
    // Emitted in EVERY arm, knob-off included: the knob-off arm is what
    // measures how often each divergence is reachable at all.
    {
        const auto& pa = WireCell::Clus::PR::g_port_audit;
        const uint64_t moved = pa.pca_refine_moved.load();
        SPDLOG_LOGGER_INFO(log,
            "PR30AUDIT oov_iso={} oov_dead={} oov_uniq={} "
            "addseg={} addseg_reentry={} ep_mismatch={} ep_refused={} "
            "pca_calls={} pca_moved={} pca_mean_cm={:.4f} pca_max_cm={:.4f} "
            "oseg_proto={} oseg_relaxed_only={} oseg_reject={} "
            "oseg_iso_drop={} oseg_iso_keep={} "
            "knobs[fit_exclusion={} endpoint_strict={} oov_parity={} local_pca={} relaxed_accept={}]",
            pa.oov_isochronous.load(), pa.oov_dead_scan.load(), pa.oov_unique_scan.load(),
            pa.add_segment_calls.load(), pa.add_segment_reentry.load(),
            pa.endpoint_mismatch.load(), pa.endpoint_refused.load(),
            pa.pca_refine_calls.load(), moved,
            moved ? (double(pa.pca_move_um_sum.load()) / moved) * units::um / units::cm : 0.0,
            double(pa.pca_move_um_max.load()) * units::um / units::cm,
            pa.oseg_accept_proto.load(), pa.oseg_accept_relaxed.load(), pa.oseg_reject.load(),
            pa.oseg_isolated_drop.load(), pa.oseg_isolated_keep.load(),
            m_fit_exclusion, m_graph_endpoint_strict, m_oov_prototype_parity,
            m_first_seg_local_pca, m_other_seg_relaxed_accept);
    }

    // doc sbnd_xin/docs/pr/32 §11 -- same contract as PR30AUDIT above.
    {
        const auto& pb = WireCell::Clus::PR::g_pr32_audit;
        const uint64_t nfit = pb.f1_fit_valid.load();
        SPDLOG_LOGGER_INFO(log,
            "PR32AUDIT f1_reads={} f1_fit_valid={} f1_mean_cm={:.4f} f1_max_cm={:.4f} "
            "f2_gate={} f2_disagree={} f2_body={} f2_demoted={} "
            "f3_cand={} f3_dropped={} f4_flagged={} "
            "knobs[use_fit_point={} traj_parity={} require_desc={} cand_flag={}]",
            pb.f1_reads.load(), nfit,
            nfit ? (double(pb.f1_moved_um_sum.load()) / nfit) * units::um / units::cm : 0.0,
            double(pb.f1_moved_um_max.load()) * units::um / units::cm,
            pb.f2_gate_calls.load(), pb.f2_gate_disagree.load(),
            pb.f2_body_runs.load(), pb.f2_flag_cleared.load(),
            pb.f3_candidates.load(), pb.f3_dropped.load(), pb.f4_flagged.load(),
            m_vertex_dir_use_fit_point, m_shower_traj_recheck_parity,
            m_main_vertex_require_descriptor, m_main_vertex_candidate_flag);
    }

    // doc sbnd_xin/docs/pr/31 §12 -- same contract as PR30AUDIT above.  The
    // two F9 counters are the §10.10 reachability measurement: both 0 across
    // the manifest => F9 (self-loop / parallel-edge representability) closes
    // as vacuous.
    {
        const auto& pa = WireCell::Clus::PR::g_port_audit;
        SPDLOG_LOGGER_INFO(log,
            "PR31AUDIT selfloop_segment={} edge_aliased={} "
            "knobs[cont_muon_dir3_30cm={} track_comp_empty_abstain={} shower_topo_reset={} "
            "reclass_preserve_4mom={} dir_track_median_local={} examine_showers_vertex_by_index={} "
            "shower_topo_proto_dir={}]",
            pa.selfloop_segment.load(), pa.edge_aliased.load(),
            m_cont_muon_dir3_30cm, m_track_comp_empty_abstain, m_shower_topo_reset,
            m_reclass_preserve_4mom, m_dir_track_median_local, m_examine_showers_vertex_by_index,
            m_shower_topo_proto_dir);
    }

    // doc sbnd_xin/docs/pr/36 §10 -- same contract as PR30AUDIT above.  The
    // per-site vectors are fixed-width so the line is machine-parseable:
    // f2_gates = the 11 skip-gates (ids on Pr36AuditCounters::f2_gate_skip),
    // f5_* = the 18 seg_endpoint_near call sites (ids on ::f5_calls).
    {
        const auto& pc = WireCell::Clus::PR::g_pr36_audit;
        auto join = [](const std::atomic<uint64_t>* a, int n) {
            std::string s;
            for (int i = 0; i < n; ++i) {
                if (i) s += ",";
                s += std::to_string(a[i].load());
            }
            return s;
        };
        SPDLOG_LOGGER_INFO(log,
            "PR36AUDIT f1_checks={} f1_disagree={} "
            "f2_sweep={}/{} f2_gates=[{}] "
            "f5_calls=[{}] f5_disagree=[{}] f5_neither=[{}] "
            "f6_id_ptr_disagree={} "
            "knobs[use_fiducial={} sp_sce={} ordered_sets={} wcpt_parity={} "
            "cluster_id_count={} nu_type_bitmask={}]",
            pc.f1_fc_checks.load(), pc.f1_fc_disagree.load(),
            pc.f2_sweep_hits.load(), pc.f2_sweep_segments.load(),
            join(pc.f2_gate_skip, pc.f2_ngates),
            join(pc.f5_calls, pc.f5_nsites),
            join(pc.f5_disagree, pc.f5_nsites),
            join(pc.f5_neither, pc.f5_nsites),
            pc.f6_id_vs_ptr_disagree.load(),
            m_use_fiducial, m_sp_sce_correction, m_tagger_ordered_segment_sets,
            m_stem_endpoint_wcpt_parity, m_broken_muon_cluster_id_count,
            m_neutrino_type_bitmask);
    }

    // doc sbnd_xin/docs/pr/33 §10 -- same contract as PR30AUDIT above.  The
    // f2 vectors are the 6 PDG-read sites (ids on Pr33AuditCounters::f2_calls);
    // f1_* are the two daughter-count sites; f3_* count pi0s minted per
    // finder (both nonzero in one event = the collision F3 prevents);
    // f5_fallback = entries into shower_less's same-index tie-break.
    {
        const auto& pd = WireCell::Clus::PR::g_pr33_audit;
        auto join = [](const std::atomic<uint64_t>* a, int n) {
            std::string s;
            for (int i = 0; i < n; ++i) {
                if (i) s += ",";
                s += std::to_string(a[i].load());
            }
            return s;
        };
        SPDLOG_LOGGER_INFO(log,
            "PR33AUDIT f1_mv={}/{}/{} f1_ex={}/{} "
            "f2_calls=[{}] f2_disagree=[{}] "
            "f3_pi0=with:{},without:{} f4_flip={} f5_fallback={} "
            "knobs[dc_mv={} dc_ex={} pdg_startseg={} pdg_showertype={} "
            "pdg_exact={} pi0_alloc={} flag_pdg_e={} shower_less_id={}]",
            pd.f1_mv_differ.load(), pd.f1_mv_gate_flip.load(), pd.f1_mv_calls.load(),
            pd.f1_ex_differ.load(), pd.f1_ex_calls.load(),
            join(pd.f2_calls, pd.f2_nsites),
            join(pd.f2_disagree, pd.f2_nsites),
            pd.f3_pi0_with_vertex.load(), pd.f3_pi0_without_vertex.load(),
            pd.f4_flip.load(), pd.f5_fallback_hits.load(),
            m_daughter_count_proto_main_vertex, m_daughter_count_proto_examine_showers,
            m_shower_pdg_from_start_segment, m_shower_pdg_from_shower_type,
            m_shower_pdg_exact_muon_test, m_pi0_id_shared_allocator,
            m_shower_flag_pdg_electron, m_shower_less_id_tiebreak);
    }
}

void TaggerCheckNeutrino::load_trackfitting_config(const std::string& config_file)
{
    try {
        // Resolve through WIRECELL_PATH so the parameter file can be named
        // relatively (e.g. 'pgrapher/experiment/sbnd/sbnd_track_fitting.json').
        // Persist::resolve() returns an absolute path unchanged, so callers that
        // already pass one are unaffected.
        const std::string resolved = Persist::resolve(config_file);
        // Load JSON file
        std::ifstream file(resolved.empty() ? config_file : resolved);
        if (!file.is_open()) {
            std::cerr << "TaggerCheckNeutrino: Cannot open config file: " << config_file
                      << " (not found on WIRECELL_PATH either)" << std::endl;
            return;
        }
        
        Json::Value root;
        Json::CharReaderBuilder builder;
        std::string errs;
        
        if (!Json::parseFromStream(builder, file, &root, &errs)) {
            std::cerr << "TaggerCheckNeutrino: Failed to parse JSON: " << errs << std::endl;
            return;
        }
        
        // Apply each parameter from the JSON file
        for (const auto& param_name : root.getMemberNames()) {
            if (param_name.substr(0, 1) == "_") continue;  // Skip comments
            
            try {
                double value = root[param_name].asDouble();
                m_track_fitter->set_parameter(param_name, value);
                // SPDLOG_LOGGER_TRACE(log, "Set {} = {}", param_name, value);
            } catch (const std::exception& e) {
                std::cerr << "TaggerCheckNeutrino: Failed to set parameter " << param_name 
                        << ": " << e.what() << std::endl;
            }
        }
        
        SPDLOG_LOGGER_TRACE(log, "Successfully loaded TrackFitting configuration");
        
    } catch (const std::exception& e) {
        std::cerr << "TaggerCheckNeutrino: Exception loading config: " << e.what() << std::endl;
        std::cerr << "TaggerCheckNeutrino: Using default TrackFitting parameters" << std::endl;
    }
}

