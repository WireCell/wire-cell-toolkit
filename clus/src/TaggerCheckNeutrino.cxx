#include "WireCellClus/TaggerCheckNeutrino.h"
#include "WireCellClus/NeutrinoPatternBase.h" // pattern recognition ...
#include "WireCellClus/PatternDebugIO.h"      // debug dump/load

#include "WireCellUtil/Persist.h"
#include <algorithm>  // doc pr/94: stable_sort of the per-bundle candidate list
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
    // pr/83 round 2 diagnostic: where in the pipeline does a class-B
    // duplicate/near-parallel pair FIRST appear?  main_vertex_graph_audit's
    // op1 (NeutrinoGraphAudit.cxx) runs ONCE, scoped to mvga_radius of the
    // main vertex; several stages (clustering_points, examine_direction,
    // shower_clustering_with_nv, the final do_multi_tracking refit) run
    // AFTER it with no second duplicate-corridor check.  This walks every
    // pair of segments on `cluster` (main cluster, unscoped -- unlike op1
    // -- since class-B pairs are not always near-vertex), using op1's own
    // metric (point-fraction within 1.4 cm, chord angle < 20 deg), and logs
    // a TRACE line per pair clearing overlap 0.5 (looser than op1's 0.7-0.8
    // production threshold, to see a pair APPROACHING duplication, not only
    // one that already cleared it) at each `detg_dump` checkpoint. Opt-in,
    // separate from WCT_DET_DEBUG (whose checksum is unaffected either way).
    void dup_stage_census(const char* tag, WireCell::Clus::PR::Graph& g,
                          WireCell::Clus::Facade::Cluster& cluster) {
        static const bool on = (std::getenv("WCT_DUP_STAGE_DEBUG") != nullptr);
        if (!on) return;
        static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");
        std::vector<SegmentPtr> segs;
        for (auto vd : ordered_nodes(g)) {
            auto vtx = g[vd].vertex;
            if (!vtx || vtx->cluster() != &cluster) continue;
            for (auto edesc : sorted_out_edges(vd, g)) {
                SegmentPtr sg = g[edesc].segment;
                if (!sg) continue;
                if (std::find(segs.begin(), segs.end(), sg) != segs.end()) continue;
                segs.push_back(sg);
            }
        }
        SPDLOG_LOGGER_TRACE(s_log, "dupstage: {} cluster={} nseg={}",
            tag, cluster.ident(), segs.size());
        auto points_of = [](SegmentPtr sg) {
            std::vector<WireCell::Point> pts;
            if (!sg->fits().empty()) {
                for (const auto& f : sg->fits()) pts.push_back(f.point);
            } else {
                for (const auto& w : sg->wcpts()) pts.push_back(w.point);
            }
            return pts;
        };
        auto frac_within = [](const std::vector<WireCell::Point>& A,
                              const std::vector<WireCell::Point>& B, double tol) {
            if (A.empty() || B.empty()) return 0.0;
            size_t n = 0;
            for (const auto& a : A) {
                double best = 1e18;
                for (const auto& b : B) best = std::min(best, (a - b).magnitude());
                if (best < tol) ++n;
            }
            return double(n) / A.size();
        };
        for (size_t i = 0; i + 1 < segs.size(); ++i) {
            for (size_t j = i + 1; j < segs.size(); ++j) {
                double la = segment_track_length(segs[i]);
                double lb = segment_track_length(segs[j]);
                if (std::min(la, lb) < 10*units::cm) continue;
                auto sPts = points_of(la <= lb ? segs[i] : segs[j]);
                auto lPts = points_of(la <= lb ? segs[j] : segs[i]);
                if (sPts.size() < 2 || lPts.size() < 2) continue;
                double f = frac_within(sPts, lPts, 1.4*units::cm);
                if (f < 0.5) continue;
                auto u = sPts.back() - sPts.front();
                auto v = lPts.back() - lPts.front();
                double den = u.magnitude() * v.magnitude();
                double ang = den > 0
                    ? std::acos(std::clamp(std::abs(u.dot(v)) / den, 0.0, 1.0)) / 3.1415926 * 180.0
                    : -1.0;
                SPDLOG_LOGGER_TRACE(s_log,
                    "dupstage: {} cluster={} pair len {:.1f}/{:.1f}cm overlap={:.2f} angle={:.1f}deg",
                    tag, cluster.ident(), std::min(la,lb)/units::cm, std::max(la,lb)/units::cm, f, ang);
            }
        }
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
    // ---- doc 80: MCS muon momentum.  mcs_enable false => the call site is
    // skipped entirely => byte-identical legacy path.  mcs_cathode_x/_xcut
    // mirror the SBND cathode_x/cathode_kink_xcut convention (one 5 cm, one
    // justification -- doc 80 sec 7.5).
    m_mcs.enable             = get(config, "mcs_enable",             m_mcs.enable);
    m_mcs.muon_source        = get(config, "mcs_muon_source",        m_mcs.muon_source);
    m_mcs.muon_min_length_cm = get(config, "mcs_muon_min_length_cm", m_mcs.muon_min_length_cm);  // cm
    m_mcs.point_source       = get(config, "mcs_point_source",       m_mcs.point_source);
    m_mcs.beam_window_only   = get(config, "mcs_beam_window_only",   m_mcs.beam_window_only);
    m_mcs.cathode_x_cm       = get(config, "mcs_cathode_x",          m_mcs.cathode_x_cm);    // cm
    m_mcs.cathode_xcut_cm    = get(config, "mcs_cathode_xcut",       m_mcs.cathode_xcut_cm); // cm
    m_mcs.max_points         = get(config, "mcs_max_points",         m_mcs.max_points);
    m_mcs.range_comparator_chain = get(config, "mcs_range_comparator_chain", m_mcs.range_comparator_chain);  // doc 84 round 1 (P5), log-only
    m_mcs.bridged_members        = get(config, "mcs_bridged_members",        m_mcs.bridged_members);  // doc 84 round 3: fit the cathode-bridge member set too (T_kine + log only)
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
    // doc sbnd_xin/docs/pr/104 -- main-vertex junction snap.
    m_vertex_junction_snap = get(config, "vertex_junction_snap", m_vertex_junction_snap);
    m_vjs_radius       = get(config, "vjs_radius",       m_vjs_radius);       // cm
    m_vjs_min_arm      = get(config, "vjs_min_arm",      m_vjs_min_arm);      // cm
    m_vjs_min_prongs   = get(config, "vjs_min_prongs",   m_vjs_min_prongs);
    m_vjs_collinear    = get(config, "vjs_collinear",    m_vjs_collinear);    // deg
    m_vjs_fit_margin   = get(config, "vjs_fit_margin",   m_vjs_fit_margin);   // cm
    m_vjs_fit_rms      = get(config, "vjs_fit_rms",      m_vjs_fit_rms);      // cm
    m_vjs_override_kink_snap = get(config, "vjs_override_kink_snap", m_vjs_override_kink_snap);
    m_vjs_min_move     = get(config, "vjs_min_move",     m_vjs_min_move);     // cm
    // doc sbnd_xin/docs/pr/51 -- main-vertex graph audit.
    m_esva_ignore_empty_2d = get(config, "esva_ignore_empty_2d", m_esva_ignore_empty_2d);  // docs/73 sec 12 round 3
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
    m_mvga_op1_radius   = get(config, "mvga_op1_radius",   m_mvga_op1_radius);   // cm; doc pr/83 r3; 0 = use mvga_radius, -1 = unscoped
    m_mvga_op1_dup_frac = get(config, "mvga_op1_dup_frac", m_mvga_op1_dup_frac); // doc pr/83 r3; 0 = use mvga_dup_frac
    m_mvga_op1_post     = get(config, "mvga_op1_post",     m_mvga_op1_post);     // doc pr/83 r3 (class A)
    m_swap_orphan_dup_audit = get(config, "swap_orphan_dup_audit", m_swap_orphan_dup_audit); // doc pr/83 r3 (Mechanism C)
    m_mvga_proj_dup_frac  = get(config, "mvga_proj_dup_frac",  m_mvga_proj_dup_frac);  // doc pr/83 r4; 0 = pass disabled
    m_mvga_proj_dqdx_ratio = get(config, "mvga_proj_dqdx_ratio", m_mvga_proj_dqdx_ratio); // doc pr/83 r4; inert while frac == 0
    m_mvga_proj_angle = get(config, "mvga_proj_angle", m_mvga_proj_angle); // deg; doc pr/83 r4b; 0 = use mvga_dup_angle
    m_mvga_ac_veto_radius = get(config, "mvga_ac_veto_radius", m_mvga_ac_veto_radius); // cm; doc pr/99 round 2; 0 = legacy
    m_mvga_ac_chord_max   = get(config, "mvga_ac_chord_max",   m_mvga_ac_chord_max);   // cm; doc pr/99 round 2; 0 = no cap
    m_mvga_ac_no_cascade  = get(config, "mvga_ac_no_cascade",  m_mvga_ac_no_cascade);  // doc pr/99 round 2
    m_mvga_passthru       = get(config, "mvga_passthru",       m_mvga_passthru);       // cm; doc pr/103; 0 = off
    m_mvga_passthru_tol   = get(config, "mvga_passthru_tol",   m_mvga_passthru_tol);   // cm; doc pr/103
    m_mvga_interposed_fallback = get(config, "mvga_interposed_fallback", m_mvga_interposed_fallback);  // doc pr/103
    m_mvga_interposed_fallback_min_angle = get(config, "mvga_interposed_fallback_min_angle", m_mvga_interposed_fallback_min_angle);  // deg; doc pr/103
    m_mvga_dup_starved_asym = get(config, "mvga_dup_starved_asym", m_mvga_dup_starved_asym); // pair asymmetry; doc pr/99 round 2; 0 = off
    m_mvga_dup_starved_mip = get(config, "mvga_dup_starved_mip", m_mvga_dup_starved_mip); // absolute cap on loser; doc pr/99 round 2; 0 = off
    m_mvga_dup_starved_span = get(config, "mvga_dup_starved_span", m_mvga_dup_starved_span); // pair length comparability; doc pr/99 round 2; 0 = off
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
    // doc sbnd_xin/docs/pr/102 P1+P2 -- keep-isolated disjuncts + 3-D uncovered radius.
    m_other_seg_keep_isolated_min_nnf    = get(config, "other_seg_keep_isolated_min_nnf",    m_other_seg_keep_isolated_min_nnf);
    m_other_seg_keep_isolated_len_admit  = get(config, "other_seg_keep_isolated_len_admit",  m_other_seg_keep_isolated_len_admit); // cm
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
    // doc sbnd_xin/docs/pr/107: dQ/dx fit keeps every trajectory point (prototype parity).
    m_dqdx_fit_keep_all_points   = get(config, "dqdx_fit_keep_all_points",   m_dqdx_fit_keep_all_points);
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
    // doc sbnd_xin/docs/pr/106 sec 10: exclusion-free charge cloud for the DL vertex net.
    m_dl_vtx_cloud_no_exclusion = get(config, "dl_vtx_cloud_no_exclusion", m_dl_vtx_cloud_no_exclusion);
    // doc sbnd_xin/docs/pr/112 sec 11: the dual chain (exclusion-free PR pass
    // suggests the neutrino vertex).  All defaults = legacy.
    m_dl_vtx_dual_chain            = get(config, "dl_vtx_dual_chain",            m_dl_vtx_dual_chain);
    m_dual_chain_mode              = get(config, "dual_chain_mode",              m_dual_chain_mode);
    m_dual_chain_transfer          = get(config, "dual_chain_transfer",          m_dual_chain_transfer);
    m_dual_chain_transfer_max      = get(config, "dual_chain_transfer_max",      m_dual_chain_transfer_max);
    m_dual_chain_allow_cluster_swap = get(config, "dual_chain_allow_cluster_swap", m_dual_chain_allow_cluster_swap);
    m_dual_chain_vtx_weight        = get(config, "dual_chain_vtx_weight",        m_dual_chain_vtx_weight);
    if (m_dl_vtx_dual_chain && m_dual_chain_mode != "snap" && m_dual_chain_mode != "voxels" && m_dual_chain_mode != "union") {
        raise<ValueError>("TaggerCheckNeutrino: dual_chain_mode must be snap, voxels or union, got %s", m_dual_chain_mode.c_str());
    }
    // doc sbnd_xin/docs/pr/89 Arm C (C2): rule-1 topology term weight/center.
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
    m_nu_fallback_demoted_mains = get(config, "nu_fallback_demoted_mains", m_nu_fallback_demoted_mains);
    // doc pr/94 Phase 2: see the member comments in the header.
    m_nu_per_bundle = get(config, "nu_per_bundle", m_nu_per_bundle);
    m_nu_per_bundle_demoted_acts = get(config, "nu_per_bundle_demoted_acts", m_nu_per_bundle_demoted_acts);
    m_nu_per_bundle_min_length = get(config, "nu_per_bundle_min_length", m_nu_per_bundle_min_length);  // cm
    m_nu_selected_as_main = get(config, "nu_selected_as_main", m_nu_selected_as_main);
    m_nu_selected_as_main_snapshot_all = get(config, "nu_selected_as_main_snapshot_all", m_nu_selected_as_main_snapshot_all);
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
    // sbnd_xin/docs/74 G1/G2: route the SAME fiducial + fv_tolerance into
    // cosmic_tagger()'s containment tests (its inside_fv lambda and the flag-1
    // vertex test).  Requires "fiducial" to be configured; inert otherwise.
    // Default false = cosmic_tagger keeps the historical FiducialUtils
    // zero-margin sensitive-volume union, byte-identical.
    m_cosmic_consistent_fv = get(config, "cosmic_consistent_fv", m_cosmic_consistent_fv);
    if (m_cosmic_consistent_fv) {
        SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino: cosmic_consistent_fv on: cosmic_tagger containment uses the configured fiducial ({} tolerance value(s), fiducial {}configured)",
                            m_fv_tolerance.size(), m_use_fiducial ? "" : "NOT ");
    }
    // sbnd_xin/docs/75: route the SAME fiducial into the nue/single-photon
    // taggers' containment tests.  Requires "fiducial" to be configured;
    // inert otherwise.  Default false = legacy FiducialUtils zero-margin
    // path, byte-identical.
    m_nue_sp_consistent_fv = get(config, "nue_sp_consistent_fv", m_nue_sp_consistent_fv);
    if (m_nue_sp_consistent_fv) {
        SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino: nue_sp_consistent_fv on: nue/single-photon tagger containment uses the configured fiducial (fiducial {}configured)",
                            m_use_fiducial ? "" : "NOT ");
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
    m_shower_endpoint_skip_orphan_vtx = get(config, "shower_endpoint_skip_orphan_vtx", m_shower_endpoint_skip_orphan_vtx);
    m_shower_walk_visited_parity = get(config, "shower_walk_visited_parity", m_shower_walk_visited_parity);
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
    m_michel_stem_muon_rescue                   = get(config, "michel_stem_muon_rescue",                   m_michel_stem_muon_rescue);
    m_shower_in_cascade_guard                   = get(config, "shower_in_cascade_guard",                   m_shower_in_cascade_guard);
    m_shower_in_max_len                         = get(config, "shower_in_max_len",                         m_shower_in_max_len);
    m_shower_in_mip_hi                          = get(config, "shower_in_mip_hi",                          m_shower_in_mip_hi);
    // doc sbnd_xin/docs/pr/40 round 9 -- straight-track PID guard family + B2 bridge.
    m_shower_connect_from_vertices_straight_guard  = get(config, "shower_connect_from_vertices_straight_guard",  m_shower_connect_from_vertices_straight_guard);
    m_shower_connect_start_seg_straight_guard      = get(config, "shower_connect_start_seg_straight_guard",      m_shower_connect_start_seg_straight_guard);
    m_examine_direction_dirsign_shower_in_guard    = get(config, "examine_direction_dirsign_shower_in_guard",    m_examine_direction_dirsign_shower_in_guard);
    m_daughter_shower_angle_reclass_straight_guard = get(config, "daughter_shower_angle_reclass_straight_guard", m_daughter_shower_angle_reclass_straight_guard);
    m_shower_topo_reexam_straight_guard            = get(config, "shower_topo_reexam_straight_guard",            m_shower_topo_reexam_straight_guard);
    m_sfv_kink_max                                 = get(config, "sfv_kink_max",                                 m_sfv_kink_max);
    m_shower_nv_bridge_track                       = get(config, "shower_nv_bridge_track",                       m_shower_nv_bridge_track);
    m_shower_nv_bridge_max_gap                     = get(config, "shower_nv_bridge_max_gap",                     m_shower_nv_bridge_max_gap);
    m_shower_nv_main_pi_init                       = get(config, "shower_nv_main_pi_init",                       m_shower_nv_main_pi_init);
    // doc sbnd_xin/docs/pr/92 -- stray-satellite drop from kine/PF.
    m_kine_drop_stray_satellites                   = get(config, "kine_drop_stray_satellites",                   m_kine_drop_stray_satellites);
    m_kine_sat_min_energy                          = get(config, "kine_sat_min_energy",                          m_kine_sat_min_energy);
    m_kine_sat_prox_max                            = get(config, "kine_sat_prox_max",                            m_kine_sat_prox_max);
    m_kine_sat_angle_bad                           = get(config, "kine_sat_angle_bad",                           m_kine_sat_angle_bad);
    m_kine_sat_angle_main                          = get(config, "kine_sat_angle_main",                          m_kine_sat_angle_main);
    m_kine_sat_far_dis                             = get(config, "kine_sat_far_dis",                             m_kine_sat_far_dis);
    m_kine_sat_axis_dis_cut                        = get(config, "kine_sat_axis_dis_cut",                        m_kine_sat_axis_dis_cut);
    m_kine_sat_cont_kink                           = get(config, "kine_sat_cont_kink",                           m_kine_sat_cont_kink);
    m_kine_sat_track_max_nseg                      = get(config, "kine_sat_track_max_nseg",                      m_kine_sat_track_max_nseg);
    m_kine_sat_em_far_dis                          = get(config, "kine_sat_em_far_dis",                          m_kine_sat_em_far_dis);
    m_michel_stem_michel_check                  = get(config, "michel_stem_michel_check",                  m_michel_stem_michel_check);
    m_michel_stem_max_far_len                   = get(config, "michel_stem_max_far_len",                   m_michel_stem_max_far_len);
    m_shower_stem_backfill                      = get(config, "shower_stem_backfill",                      m_shower_stem_backfill);
    m_stem_backfill_max_len                     = get(config, "stem_backfill_max_len",                     m_stem_backfill_max_len);
    m_stem_backfill_mip_lo                      = get(config, "stem_backfill_mip_lo",                      m_stem_backfill_mip_lo);
    m_stem_backfill_mip_hi                      = get(config, "stem_backfill_mip_hi",                      m_stem_backfill_mip_hi);
    m_stem_backfill_min_shower_len              = get(config, "stem_backfill_min_shower_len",              m_stem_backfill_min_shower_len);
    m_shower_conn3_unreachable                  = get(config, "shower_conn3_unreachable",                  m_shower_conn3_unreachable);
    m_conn3_unreachable_min_len                 = get(config, "conn3_unreachable_min_len",                 m_conn3_unreachable_min_len);
    m_conn3_stitch_max                          = get(config, "conn3_stitch_max",                          m_conn3_stitch_max);
    m_shower_dedup_start_seg                    = get(config, "shower_dedup_start_seg",                    m_shower_dedup_start_seg);
    m_shower_traj_michel_stem                   = get(config, "shower_traj_michel_stem",                   m_shower_traj_michel_stem);
    m_michel_stem_traj_min_len                  = get(config, "michel_stem_traj_min_len",                  m_michel_stem_traj_min_len);
    m_michel_stem_traj_max_len                  = get(config, "michel_stem_traj_max_len",                  m_michel_stem_traj_max_len);
    m_michel_stem_traj_mip_lo                   = get(config, "michel_stem_traj_mip_lo",                   m_michel_stem_traj_mip_lo);
    m_michel_stem_traj_max_far_len              = get(config, "michel_stem_traj_max_far_len",              m_michel_stem_traj_max_far_len);
    m_michel_stem_traj_min_kink_deg             = get(config, "michel_stem_traj_min_kink_deg",             m_michel_stem_traj_min_kink_deg);
    m_shower_long_muon_keep_type                = get(config, "shower_long_muon_keep_type",                m_shower_long_muon_keep_type);
    m_shower_bragg_protect_start_segment        = get(config, "shower_bragg_protect_start_segment",        m_shower_bragg_protect_start_segment);
    // doc pr/93 round 3
    m_shower_reclass_case_b_dqdx_guard          = get(config, "shower_reclass_case_b_dqdx_guard",          m_shower_reclass_case_b_dqdx_guard);
    m_shower_accept_pid_guard                   = get(config, "shower_accept_pid_guard",                   m_shower_accept_pid_guard);
    m_shower_pid_guard_min_len                  = get(config, "shower_pid_guard_min_len",                  m_shower_pid_guard_min_len);
    m_shower_vote_track_pid_counts              = get(config, "shower_vote_track_pid_counts",              m_shower_vote_track_pid_counts);
    m_shower_cone_absorb_guard               = get(config, "shower_cone_absorb_guard",               m_shower_cone_absorb_guard);
    // doc pr/93 round 4
    m_shower_detach_track_stem                  = get(config, "shower_detach_track_stem",                  m_shower_detach_track_stem);
    // doc pr/99 round 2
    m_shower_ghost_member_drop                  = get(config, "shower_ghost_member_drop",                  m_shower_ghost_member_drop);
    m_shower_ghost_overlap_frac                 = get(config, "shower_ghost_overlap_frac",                 m_shower_ghost_overlap_frac);
    m_shower_ghost_dqdx_ratio                   = get(config, "shower_ghost_dqdx_ratio",                   m_shower_ghost_dqdx_ratio);
    m_shower_ghost_min_len                      = get(config, "shower_ghost_min_len",                      m_shower_ghost_min_len);
    // doc pr/99 round 3 (C1/C1b/A5)
    m_kine_charge_dedup                         = get(config, "kine_charge_dedup",                         m_kine_charge_dedup);
    m_kine_charge_rebuild                       = get(config, "kine_charge_rebuild",                       m_kine_charge_rebuild);
    // doc pr/101 Enu accounting round (K1-K5)
    m_kine_charge_track_ctx                     = get(config, "kine_charge_track_ctx",                     m_kine_charge_track_ctx);
    m_kine_mass_rules                           = get(config, "kine_mass_rules",                           m_kine_mass_rules);
    m_kine_hadronic_dqdx                        = get(config, "kine_hadronic_dqdx",                        m_kine_hadronic_dqdx);
    m_kine_long_muon_mode                       = get(config, "kine_long_muon_mode",                       m_kine_long_muon_mode);
    m_kine_long_muon_ratio_lo                   = get(config, "kine_long_muon_ratio_lo",                   m_kine_long_muon_ratio_lo);
    m_kine_long_muon_ratio_hi                   = get(config, "kine_long_muon_ratio_hi",                   m_kine_long_muon_ratio_hi);
    m_long_muon_range_empty_chain_fallback      = get(config, "long_muon_range_empty_chain_fallback",      m_long_muon_range_empty_chain_fallback);  // doc 84 round 1 (P1)
    m_long_muon_members_geometry                = get(config, "long_muon_members_geometry",                m_long_muon_members_geometry);             // doc 84 round 2
    m_long_muon_cathode_bridge                  = get(config, "long_muon_cathode_bridge",                  m_long_muon_cathode_bridge);               // doc 84 round 2
    m_long_muon_cathode_bridge_x                = get(config, "long_muon_cathode_bridge_x",                m_long_muon_cathode_bridge_x);             // doc 84 round 2, cm
    m_long_muon_cathode_bridge_xcut             = get(config, "long_muon_cathode_bridge_xcut",             m_long_muon_cathode_bridge_xcut);          // doc 84 round 2, cm
    m_long_muon_cathode_bridge_gap              = get(config, "long_muon_cathode_bridge_gap",              m_long_muon_cathode_bridge_gap);           // doc 84 round 2, cm
    m_long_muon_cathode_bridge_angle            = get(config, "long_muon_cathode_bridge_angle",            m_long_muon_cathode_bridge_angle);         // doc 84 round 2, deg
    m_long_muon_cathode_bridge_lever            = get(config, "long_muon_cathode_bridge_lever",            m_long_muon_cathode_bridge_lever);         // doc 84 round 4 G1, cm
    m_long_muon_cathode_bridge_track_partner    = get(config, "long_muon_cathode_bridge_track_partner",    m_long_muon_cathode_bridge_track_partner); // doc 84 round 4 G2
    m_long_muon_cathode_bridge_short_gap        = get(config, "long_muon_cathode_bridge_short_gap",        m_long_muon_cathode_bridge_short_gap);     // doc 84 round 4 G3, cm
    m_long_muon_cathode_bridge_short_gap_angle  = get(config, "long_muon_cathode_bridge_short_gap_angle",  m_long_muon_cathode_bridge_short_gap_angle); // doc 84 round 4 G3, deg
    m_long_muon_cathode_bridge_short_gap_len    = get(config, "long_muon_cathode_bridge_short_gap_len",    m_long_muon_cathode_bridge_short_gap_len); // doc 84 round 4 G3, cm
    m_kine_mainvtx_used_guard                   = get(config, "kine_mainvtx_used_guard",                   m_kine_mainvtx_used_guard);
    m_shower_hadronic_tag                       = get(config, "shower_hadronic_tag",                       m_shower_hadronic_tag);
    m_shower_hadronic_min_len                   = get(config, "shower_hadronic_min_len",                   m_shower_hadronic_min_len);
    m_shower_hadronic_scan_len                  = get(config, "shower_hadronic_scan_len",                  m_shower_hadronic_scan_len);
    m_shower_hadronic_bin                       = get(config, "shower_hadronic_bin",                       m_shower_hadronic_bin);
    m_shower_hadronic_r_cyl                     = get(config, "shower_hadronic_r_cyl",                     m_shower_hadronic_r_cyl);
    m_shower_hadronic_r_core                    = get(config, "shower_hadronic_r_core",                    m_shower_hadronic_r_core);
    m_shower_hadronic_growth_max                = get(config, "shower_hadronic_growth_max",                m_shower_hadronic_growth_max);
    m_shower_hadronic_growth_bragg              = get(config, "shower_hadronic_growth_bragg",              m_shower_hadronic_growth_bragg);
    m_shower_hadronic_bragg_ratio               = get(config, "shower_hadronic_bragg_ratio",               m_shower_hadronic_bragg_ratio);
    m_shower_hadronic_stem_ratio                = get(config, "shower_hadronic_stem_ratio",                m_shower_hadronic_stem_ratio);
    m_kine_count_orphan_tracks                  = get(config, "kine_count_orphan_tracks",                  m_kine_count_orphan_tracks);
    m_kine_orphan_track_min                     = get(config, "kine_orphan_track_min",                     m_kine_orphan_track_min);
    // doc sbnd_xin/docs/pr/117 round 1
    m_shower_pass4_best_owner                   = get(config, "shower_pass4_best_owner",                   m_shower_pass4_best_owner);
    m_shower_merge_relax                        = get(config, "shower_merge_relax",                        m_shower_merge_relax);
    m_shower_merge_relax_dis                    = get(config, "shower_merge_relax_dis",                    m_shower_merge_relax_dis);
    m_shower_merge_relax_angle                  = get(config, "shower_merge_relax_angle",                  m_shower_merge_relax_angle);
    m_shower_merge_relax_min_len                = get(config, "shower_merge_relax_min_len",                m_shower_merge_relax_min_len);
    m_shower_flank_absorb                       = get(config, "shower_flank_absorb",                       m_shower_flank_absorb);
    m_shower_flank_absorb_max_dis               = get(config, "shower_flank_absorb_max_dis",               m_shower_flank_absorb_max_dis);
    m_shower_flank_absorb_max_len               = get(config, "shower_flank_absorb_max_len",               m_shower_flank_absorb_max_len);
    // doc sbnd_xin/docs/pr/118 round 1
    m_shower_ex1_conn3_body_dis                 = get(config, "shower_ex1_conn3_body_dis",                 m_shower_ex1_conn3_body_dis);
    m_shower_merge_relax_continuity             = get(config, "shower_merge_relax_continuity",             m_shower_merge_relax_continuity);
    m_shower_merge_relax_cont_frac              = get(config, "shower_merge_relax_cont_frac",              m_shower_merge_relax_cont_frac);
    m_shower_merge_relax_cont_gap               = get(config, "shower_merge_relax_cont_gap",               m_shower_merge_relax_cont_gap);
    m_shower_merge_relax_cont_qmed              = get(config, "shower_merge_relax_cont_qmed",              m_shower_merge_relax_cont_qmed);
    m_shower_merge_relax_cont_axis              = get(config, "shower_merge_relax_cont_axis",              m_shower_merge_relax_cont_axis);
    m_shower_merge_relax_cont_dmax              = get(config, "shower_merge_relax_cont_dmax",              m_shower_merge_relax_cont_dmax);
    m_shower_merge_relax_cont_t1_gap            = get(config, "shower_merge_relax_cont_t1_gap",            m_shower_merge_relax_cont_t1_gap);
    m_shower_merge_relax_cont_t1_fold           = get(config, "shower_merge_relax_cont_t1_fold",           m_shower_merge_relax_cont_t1_fold);
    // doc sbnd_xin/docs/pr/120 round 1
    m_stem_backfill_back_guard                  = get(config, "stem_backfill_back_guard",                  m_stem_backfill_back_guard);
    m_stem_backfill_back_ang                    = get(config, "stem_backfill_back_ang",                    m_stem_backfill_back_ang);
    m_shower_ex1_walk_em_track_guard            = get(config, "shower_ex1_walk_em_track_guard",            m_shower_ex1_walk_em_track_guard);
    m_shower_ex1_walk_em_track_len              = get(config, "shower_ex1_walk_em_track_len",              m_shower_ex1_walk_em_track_len);
    // doc sbnd_xin/docs/pr/121 round 1
    m_shower_ex1_dedup_rehome                   = get(config, "shower_ex1_dedup_rehome",                   m_shower_ex1_dedup_rehome);
    m_shower_pass4_prune_detached               = get(config, "shower_pass4_prune_detached",               m_shower_pass4_prune_detached);            // doc pr/123 r1
    m_shower_pass4_prune_gap                    = get(config, "shower_pass4_prune_gap",                    m_shower_pass4_prune_gap);                 // doc pr/123 r1, cm
    m_shower_pass4_prune_gap2                   = get(config, "shower_pass4_prune_gap2",                   m_shower_pass4_prune_gap2);                // doc pr/124 A, cm
    m_shower_pass4_prune2_ang                   = get(config, "shower_pass4_prune2_ang",                   m_shower_pass4_prune2_ang);                // doc pr/124 A, deg
    m_shower_pass4_prune2_mdqdx                 = get(config, "shower_pass4_prune2_mdqdx",                 m_shower_pass4_prune2_mdqdx);              // doc pr/124 A, x MIP
    m_shower_pass3_cone_guard_len               = get(config, "shower_pass3_cone_guard_len",               m_shower_pass3_cone_guard_len);            // doc pr/124 C, cm
    m_shower_samevtx_track_absorb               = get(config, "shower_samevtx_track_absorb",               m_shower_samevtx_track_absorb);            // doc pr/125
    m_shower_samevtx_absorb_gap                 = get(config, "shower_samevtx_absorb_gap",                 m_shower_samevtx_absorb_gap);              // doc pr/125, cm
    m_shower_samevtx_absorb_max_len             = get(config, "shower_samevtx_absorb_max_len",             m_shower_samevtx_absorb_max_len);          // doc pr/125, cm
    m_shower_samevtx_absorb_min_len             = get(config, "shower_samevtx_absorb_min_len",             m_shower_samevtx_absorb_min_len);          // doc pr/125, cm
    m_shower_satellite_absorb                   = get(config, "shower_satellite_absorb",                   m_shower_satellite_absorb);                // doc pr/125
    m_shower_split = get(config, "shower_split", m_shower_split);   // doc pr/138 B2; false = no pass
    m_shower_split_max_valley = get(config, "shower_split_max_valley", m_shower_split_max_valley);   // doc pr/138 B2; sec A5.4 knee
    m_shower_split_min_frac = get(config, "shower_split_min_frac", m_shower_split_min_frac);   // doc pr/138 B2; per-seed charge share floor
    m_shower_split_max_parts = get(config, "shower_split_max_parts", m_shower_split_max_parts);   // doc pr/138 B3; 2 = the measured-exact kernel
    m_shower_split_min_charge = get(config, "shower_split_min_charge", m_shower_split_min_charge);   // doc pr/138 B1; candidate charge floor (raw Fit::dQ)
    m_shower_split_min_nseg = get(config, "shower_split_min_nseg", m_shower_split_min_nseg);   // doc pr/138 B1; candidate member-count floor
    m_shower_split_bundle_gap = get(config, "shower_split_bundle_gap", m_shower_split_bundle_gap);   // doc pr/138 B3; single-linkage bundle gap
    m_shower_split_snap = get(config, "shower_split_snap", m_shower_split_snap);   // doc pr/138 B3; k>=3 bundle dominance floor
    m_shower_split_skip_shared = get(config, "shower_split_skip_shared", m_shower_split_skip_shared);   // doc pr/139 P1.1; refuse a component holding a segment another shower also owns
    m_shower_split_max_impact = get(config, "shower_split_max_impact", m_shower_split_max_impact);   // doc pr/139 P1.2; cm; 0 = no bound
    m_shower_split_em_start = get(config, "shower_split_em_start", m_shower_split_em_start);   // doc pr/139 P1.3; seed the daughter on its nearest EM-typed member
    m_shower_split_rehome = get(config, "shower_split_rehome", m_shower_split_rehome);   // doc pr/139 P1.4; offer an orphan daughter to the nearest larger EM shower
    m_shower_split_rehome_gap = get(config, "shower_split_rehome_gap", m_shower_split_rehome_gap);   // doc pr/139 P1.4; cm; max daughter->host 3-D gap
    m_shower_satellite_absorb_max_mev           = get(config, "shower_satellite_absorb_max_mev",           m_shower_satellite_absorb_max_mev);        // doc pr/125, MeV
    m_shower_satellite_absorb_host_mev          = get(config, "shower_satellite_absorb_host_mev",          m_shower_satellite_absorb_host_mev);       // doc pr/125, MeV
    m_shower_pass4_track_guard_len              = get(config, "shower_pass4_track_guard_len",              m_shower_pass4_track_guard_len);           // doc pr/123 r1, cm
    m_shower_pass4_prox_guard_len               = get(config, "shower_pass4_prox_guard_len",               m_shower_pass4_prox_guard_len);            // doc pr/130 item 1b, cm
    m_shower_pass3_backfill_guard_len           = get(config, "shower_pass3_backfill_guard_len",           m_shower_pass3_backfill_guard_len);        // doc pr/130 item 1b, cm
    m_stem_backfill_back_dvtx                   = get(config, "stem_backfill_back_dvtx",                   m_stem_backfill_back_dvtx);                // doc pr/130 item B, cm
    m_shower_pass4_prefilter_v1_escape           = get(config, "shower_pass4_prefilter_v1_escape",           m_shower_pass4_prefilter_v1_escape);       // doc pr/136 r2
    m_shower_pass4_prefilter_v1_max_v2          = get(config, "shower_pass4_prefilter_v1_max_v2",          m_shower_pass4_prefilter_v1_max_v2);       // doc pr/136 r2, deg
    m_shower_pass4_prefilter_v1_max_dis         = get(config, "shower_pass4_prefilter_v1_max_dis",         m_shower_pass4_prefilter_v1_max_dis);      // doc pr/136 r3, cm
    m_pi0_mass_offset                           = get(config, "pi0_mass_offset",                           m_pi0_mass_offset);                        // doc pr/132 K1, MeV
    m_pi0_assoc_angle_deg                       = get(config, "pi0_assoc_angle_deg",                       m_pi0_assoc_angle_deg);                    // doc pr/132 K2, deg
    m_pi0_attached_partner_min_mev              = get(config, "pi0_attached_partner_min_mev",              m_pi0_attached_partner_min_mev);           // doc pr/132 K3, MeV
    m_pi0_nv_allow_type2                        = get(config, "pi0_nv_allow_type2",                        m_pi0_nv_allow_type2);                     // doc pr/132 K4
    m_pi0_nv_max_prongs                         = get(config, "pi0_nv_max_prongs",                         m_pi0_nv_max_prongs);                      // doc pr/132 K5
    m_pi0_readmit_retyped                       = get(config, "pi0_readmit_retyped",                       m_pi0_readmit_retyped);                    // doc pr/132 K7
    m_pi0_admit_type3                           = get(config, "pi0_admit_type3",                           m_pi0_admit_type3);                        // doc pr/132 K8
    m_pi0_crumb_assoc_mev                       = get(config, "pi0_crumb_assoc_mev",                       m_pi0_crumb_assoc_mev);                    // doc pr/132 K9, MeV
    m_pi0_collinear_merge_deg                   = get(config, "pi0_collinear_merge_deg",                   m_pi0_collinear_merge_deg);                // doc pr/132 K12, deg
    m_pi0_nv_partner_min_mev                    = get(config, "pi0_nv_partner_min_mev",                    m_pi0_nv_partner_min_mev);                 // doc pr/132 K13, MeV
    m_pi0_nv_retry_paired                       = get(config, "pi0_nv_retry_paired",                       m_pi0_nv_retry_paired);                    // doc pr/132 K14
    m_pi0_reseat_start_assoc                    = get(config, "pi0_reseat_start_assoc",                    m_pi0_reseat_start_assoc);                 // doc pr/132 K15
    m_shower_em_collinear_deg                   = get(config, "shower_em_collinear_deg",                   m_shower_em_collinear_deg);                // doc pr/132 K16, deg
    m_shower_em_collinear_dis_cm                = get(config, "shower_em_collinear_dis_cm",                m_shower_em_collinear_dis_cm);             // doc pr/132 K16, cm
    m_shower_em_collinear_host_mev              = get(config, "shower_em_collinear_host_mev",              m_shower_em_collinear_host_mev);           // doc pr/132 K16, MeV
    m_shower_em_backext_perp_cm                 = get(config, "shower_em_backext_perp_cm",                 m_shower_em_backext_perp_cm);              // doc pr/132 K17, cm
    m_shower_em_backext_len_cm                  = get(config, "shower_em_backext_len_cm",                  m_shower_em_backext_len_cm);               // doc pr/132 K17, cm
    m_pi0_accept_merge_dis_cm                   = get(config, "pi0_accept_merge_dis_cm",                   m_pi0_accept_merge_dis_cm);                // doc pr/132 K18, cm
    m_pi0_bp_vertex_miss_cm                     = get(config, "pi0_bp_vertex_miss_cm",                     m_pi0_bp_vertex_miss_cm);                  // doc pr/132 K19, cm
    m_pi0_admit_muon_showers                    = get(config, "pi0_admit_muon_showers",                    m_pi0_admit_muon_showers);                 // doc pr/133 K20
    m_pi0_nc_sig_angle_deg                      = get(config, "pi0_nc_sig_angle_deg",                      m_pi0_nc_sig_angle_deg);                   // doc pr/133 K21, deg
    m_pi0_nc_floor_mev                          = get(config, "pi0_nc_floor_mev",                          m_pi0_nc_floor_mev);                       // doc pr/133 K21 v2, MeV
    m_pi0_nc_pf_assoc_deg                       = get(config, "pi0_nc_pf_assoc_deg",                       m_pi0_nc_pf_assoc_deg);                    // doc pr/133 K21 v2.2, deg
    m_pi0_nc_frag_merge                         = get(config, "pi0_nc_frag_merge",                         m_pi0_nc_frag_merge);                      // doc pr/134 K22
    m_pi0_pf_assoc_deg                          = get(config, "pi0_pf_assoc_deg",                          m_pi0_pf_assoc_deg);                       // doc pr/134 K23, deg
    m_pi0_prefer_main_vertex                    = get(config, "pi0_prefer_main_vertex",                    m_pi0_prefer_main_vertex);                 // doc pr/134 K24
    m_pi0_nv_max_vtx_shift_cm                   = get(config, "pi0_nv_max_vtx_shift_cm",                   m_pi0_nv_max_vtx_shift_cm);                // doc pr/132 K10, cm
    m_pi0_nv_mass_window_mev                    = get(config, "pi0_nv_mass_window_mev",                    m_pi0_nv_mass_window_mev);                 // doc pr/132 K11, MeV
    m_kine_count_guard_freed                    = get(config, "kine_count_guard_freed",                    m_kine_count_guard_freed);                 // doc pr/123 r2
    m_kine_guard_freed_impact                   = get(config, "kine_guard_freed_impact",                   m_kine_guard_freed_impact);                // doc pr/129
    m_kine_guard_freed_miss_deg                 = get(config, "kine_guard_freed_miss_deg",                 m_kine_guard_freed_miss_deg);              // doc pr/129
    m_kine_count_near_cross_cluster             = get(config, "kine_count_near_cross_cluster",             m_kine_count_near_cross_cluster);          // doc pr/128
    m_kine_near_gap                             = get(config, "kine_near_gap",                             m_kine_near_gap);                          // doc pr/128
    m_kine_near_min_len                         = get(config, "kine_near_min_len",                         m_kine_near_min_len);                      // doc pr/128
    m_kine_near_end_tol                         = get(config, "kine_near_end_tol",                         m_kine_near_end_tol);                      // doc pr/128
    m_kine_near_kink_deg                        = get(config, "kine_near_kink_deg",                        m_kine_near_kink_deg);                     // doc pr/128
    m_kine_count_conn4_near                     = get(config, "kine_count_conn4_near",                     m_kine_count_conn4_near);                  // doc pr/128
    m_kine_conn4_near_gap                       = get(config, "kine_conn4_near_gap",                       m_kine_conn4_near_gap);                    // doc pr/128
    m_straight_cont_cross_cluster               = get(config, "straight_cont_cross_cluster",               m_straight_cont_cross_cluster);
    m_sccc_bridge_body                          = get(config, "sccc_bridge_body",                          m_sccc_bridge_body);
    m_sccc_max_gap                              = get(config, "sccc_max_gap",                              m_sccc_max_gap);
    m_sccc_kink_max                             = get(config, "sccc_kink_max",                             m_sccc_kink_max);
    m_sccc_gap_aligned                          = get(config, "sccc_gap_aligned",                          m_sccc_gap_aligned);
    m_sccc_kink_tight                           = get(config, "sccc_kink_tight",                           m_sccc_kink_tight);
    m_single_muon_proton_chain_veto             = get(config, "single_muon_proton_chain_veto",             m_single_muon_proton_chain_veto);
    m_single_muon_long_muon_claim               = get(config, "single_muon_long_muon_claim",               m_single_muon_long_muon_claim);
    m_pid_flag_reconcile                        = get(config, "pid_flag_reconcile",                        m_pid_flag_reconcile);
    m_long_muon_stub_bridge                     = get(config, "long_muon_stub_bridge",                     m_long_muon_stub_bridge);
    m_long_muon_stub_bridge_len                 = get(config, "long_muon_stub_bridge_len",                 m_long_muon_stub_bridge_len);              // doc 84 round 1 (P3), cm
    m_long_muon_angle_relax_long                = get(config, "long_muon_angle_relax_long",                m_long_muon_angle_relax_long);             // doc 84 round 1 (P2)
    m_long_muon_angle_relax_deg                 = get(config, "long_muon_angle_relax_deg",                 m_long_muon_angle_relax_deg);              // doc 84 round 1 (P2), deg

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
    // doc 80: MCS muon momentum.  All defaults = legacy no-op (enable false).
    cfg["mcs_enable"]             = m_mcs.enable;
    cfg["mcs_muon_source"]        = m_mcs.muon_source;         // pf_muon | long_muon | longest_segment
    cfg["mcs_muon_min_length_cm"] = m_mcs.muon_min_length_cm;  // cm
    cfg["mcs_point_source"]       = m_mcs.point_source;        // muon_segments | whole_event (validation)
    cfg["mcs_beam_window_only"]   = m_mcs.beam_window_only;    // correctness gate, doc 80 sec 7.4
    cfg["mcs_cathode_x"]          = m_mcs.cathode_x_cm;        // cm
    cfg["mcs_cathode_xcut"]       = m_mcs.cathode_xcut_cm;     // cm; 0 = off (SBND: 5)
    cfg["mcs_max_points"]         = m_mcs.max_points;
    cfg["mcs_range_comparator_chain"] = m_mcs.range_comparator_chain;  // doc 84 round 1 (P5), log-only
    cfg["mcs_bridged_members"]        = m_mcs.bridged_members;         // doc 84 round 3: fit the cathode-bridge member set too
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
    // doc sbnd_xin/docs/pr/104 -- main-vertex junction snap; false => the
    // pass never fires => byte-identical.  Numerics cm/deg.
    cfg["vertex_junction_snap"] = m_vertex_junction_snap;
    cfg["vjs_radius"]       = m_vjs_radius;
    cfg["vjs_min_arm"]      = m_vjs_min_arm;
    cfg["vjs_min_prongs"]   = m_vjs_min_prongs;
    cfg["vjs_collinear"]    = m_vjs_collinear;
    cfg["vjs_fit_margin"]   = m_vjs_fit_margin;
    cfg["vjs_fit_rms"]      = m_vjs_fit_rms;
    cfg["vjs_override_kink_snap"] = m_vjs_override_kink_snap;
    cfg["vjs_min_move"]     = m_vjs_min_move;
    // doc sbnd_xin/docs/pr/51 -- main-vertex graph audit; false => the pass
    // never fires => byte-identical.  Numerics cm/deg.
    cfg["esva_ignore_empty_2d"] = m_esva_ignore_empty_2d;  // docs/73 sec 12 round 3; case-5 empty-2D-index sentinel = no info, not coverage; false = legacy
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
    cfg["mvga_op1_radius"]   = m_mvga_op1_radius;   // cm; 0 = use mvga_radius, -1 = unscoped, byte-identical at 0 (doc pr/83 r3)
    cfg["mvga_op1_dup_frac"] = m_mvga_op1_dup_frac; // fraction; 0 = use mvga_dup_frac, byte-identical (doc pr/83 r3)
    cfg["mvga_op1_post"]     = m_mvga_op1_post;     // false = post-op3 dup pass skipped, byte-identical (doc pr/83 r3)
    cfg["swap_orphan_dup_audit"] = m_swap_orphan_dup_audit; // false = abandoned main cluster never audited, byte-identical (doc pr/83 r3)
    cfg["mvga_proj_dup_frac"]  = m_mvga_proj_dup_frac;  // 0 = projective dup collapse disabled, byte-identical (doc pr/83 r4)
    cfg["mvga_proj_dqdx_ratio"] = m_mvga_proj_dqdx_ratio; // stem dQ/dx asymmetry gate; inert while frac == 0 (doc pr/83 r4)
    cfg["mvga_proj_angle"] = m_mvga_proj_angle; // deg; 0 = use mvga_dup_angle, byte-identical (doc pr/83 r4b)
    cfg["mvga_ac_veto_radius"] = m_mvga_ac_veto_radius; // cm; 0 = legacy straighten_radius rule, byte-identical (doc pr/99 round 2)
    cfg["mvga_ac_chord_max"]   = m_mvga_ac_chord_max;   // cm; 0 = no cap, byte-identical (doc pr/99 round 2)
    cfg["mvga_ac_no_cascade"]  = m_mvga_ac_no_cascade;  // false = created products collapsible, byte-identical (doc pr/99 round 2)
    cfg["mvga_passthru"]       = m_mvga_passthru;       // cm; 0 = op0 pass-through split off, byte-identical (doc pr/103)
    cfg["mvga_passthru_tol"]   = m_mvga_passthru_tol;   // cm; inert while mvga_passthru == 0 (doc pr/103)
    cfg["mvga_interposed_fallback"] = m_mvga_interposed_fallback;  // false = far-angle decline stands, byte-identical (doc pr/103)
    cfg["mvga_interposed_fallback_min_angle"] = m_mvga_interposed_fallback_min_angle;  // deg; inert while the fallback is off (doc pr/103)
    cfg["mvga_dup_starved_asym"] = m_mvga_dup_starved_asym; // pair asymmetry; 0 = angle decline stands, byte-identical (doc pr/99 round 2)
    cfg["mvga_dup_starved_mip"] = m_mvga_dup_starved_mip; // absolute cap on loser; 0 = angle decline stands, byte-identical (doc pr/99 round 2)
    cfg["mvga_dup_starved_span"] = m_mvga_dup_starved_span; // pair length comparability floor; 0 = no span test (doc pr/99 round 2)
    cfg["kink_dqdx_hot_ratio"] = m_kink_dqdx_hot_ratio; // x mip_dqdx_median; inert while both above are false
    cfg["shower_topo_demote_len"] = m_shower_topo_demote_len;  // cm; 0 = legacy (long segments stay eligible for kShowerTopology)
    // doc sbnd_xin/docs/pr/49.
    cfg["fit_blob_coverage"] = m_fit_blob_coverage; // -1 = legacy (no foreign-ghost deweighting); >= 0 = tolerance cells
    // doc sbnd_xin/docs/pr/50.
    cfg["fit_blob_coverage_defer"] = m_fit_blob_coverage_defer; // false = pr/49 behavior (deweight active during find_proto_vertex)
    // doc sbnd_xin/docs/pr/30 §11.  Round-tripped so the compiled config
    // records the operating point; each default reproduces the pre-pr/30 tree.
    cfg["fit_exclusion"]            = m_fit_exclusion;             // false = legacy (all sites pass flag_exclusion=false)
    cfg["graph_endpoint_tol"]       = m_graph_endpoint_tol;        // cm
    cfg["oov_prototype_parity"]     = m_oov_prototype_parity;      // false = legacy (today's three polarities)
    cfg["first_seg_local_pca"]      = m_first_seg_local_pca;       // true  = legacy (the refinement runs)
    cfg["other_seg_relaxed_accept"] = m_other_seg_relaxed_accept;  // true  = legacy (the 0.72/15cm/1.05 clause is live)
    cfg["other_seg_empty_2d_guard"] = m_other_seg_empty_2d_guard;  // false = legacy (-1 sentinel counts as covered)
    // doc sbnd_xin/docs/pr/54.
    cfg["other_seg_keep_isolated"]            = m_other_seg_keep_isolated;            // false = legacy (isolated residual discarded)
    cfg["other_seg_keep_isolated_min_points"] = m_other_seg_keep_isolated_min_points; // component-point floor when the keep is on
    cfg["other_seg_keep_isolated_min_length"] = m_other_seg_keep_isolated_min_length; // cm; fitted-length floor when the keep is on
    // doc sbnd_xin/docs/pr/102 P1+P2.
    cfg["other_seg_keep_isolated_min_nnf"]    = m_other_seg_keep_isolated_min_nnf;    // 0 = off; nnf disjunct on the keep
    cfg["other_seg_keep_isolated_len_admit"]  = m_other_seg_keep_isolated_len_admit;  // cm; 0 = off; length disjunct on the keep
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
    cfg["dqdx_fit_keep_all_points"]   = m_dqdx_fit_keep_all_points;    // doc pr/107: false = legacy (pre-dQ/dx pass drops zero-quantity points)
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
    cfg["cosmic_consistent_fv"] = m_cosmic_consistent_fv;  // doc 74 G1/G2; false = FiducialUtils fallback
    cfg["nue_sp_consistent_fv"] = m_nue_sp_consistent_fv;  // doc 75; false = FiducialUtils fallback
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
    cfg["dl_vtx_cloud_no_exclusion"] = m_dl_vtx_cloud_no_exclusion;  // doc pr/106 sec 10: false = legacy (net sees the exclusion fit's charge)
    cfg["dl_vtx_dual_chain"]       = m_dl_vtx_dual_chain;       // doc pr/112 sec 11: false = legacy (no exclusion-free second pass)
    cfg["dual_chain_mode"]         = m_dual_chain_mode;         // doc pr/112 sec 11: snap | voxels | union
    cfg["dual_chain_transfer"]     = m_dual_chain_transfer;     // doc pr/112 sec 11: false = probe (pass runs, nothing moves)
    cfg["dual_chain_transfer_max"] = m_dual_chain_transfer_max; // doc pr/112 sec 11: cm, snap guard D (prototype dl_vtx_cut 2.0)
    cfg["dual_chain_allow_cluster_swap"] = m_dual_chain_allow_cluster_swap;  // doc pr/112 sec 5.7.8
    cfg["dual_chain_vtx_weight"]   = m_dual_chain_vtx_weight;   // doc pr/112 sec 11: union-mode proximity term, 0 = none
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
    cfg["nu_fallback_demoted_mains"] = m_nu_fallback_demoted_mains;    // docs/73 sec 12 round 3; when NO candidate survives, consider demoted mains (same gates); false = legacy
    cfg["nu_per_bundle"]             = m_nu_per_bundle;                // doc pr/94; false = legacy single event-wide candidate
    cfg["nu_per_bundle_demoted_acts"] = m_nu_per_bundle_demoted_acts;  // doc pr/94; mirror of the taggers' evaluate_demoted_mains; inert unless nu_per_bundle
    cfg["nu_per_bundle_min_length"] = m_nu_per_bundle_min_length;      // doc pr/94 Phase 5b round 2; cm; length floor for a per-bundle candidate, exempting the legacy event-wide winner; 0 = none
    cfg["nu_selected_as_main"]      = m_nu_selected_as_main;           // doc pr/94 round 3; give a demoted-main candidate the main-cluster PR treatment for the duration of its own pass; false = legacy
    cfg["nu_selected_as_main_snapshot_all"] = m_nu_selected_as_main_snapshot_all;  // doc 75; closes the DL-swap flag leak; false = legacy
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
    cfg["shower_endpoint_skip_orphan_vtx"] = m_shower_endpoint_skip_orphan_vtx; // false = legacy end_point search includes vertices no member segment touches (doc pr/91 F1)
    cfg["shower_walk_visited_parity"] = m_shower_walk_visited_parity; // false = legacy has_node()-gated flood-fill frontier (doc pr/91 round 3)
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
    cfg["michel_stem_muon_rescue"]                   = m_michel_stem_muon_rescue;                   // false = legacy (Michel rescue limited to weak-dir degree-2 vertices)
    cfg["shower_in_cascade_guard"]                   = m_shower_in_cascade_guard;                   // doc pr/74 round 2 P1; false = legacy (cascade relabels unconditionally)
    cfg["shower_in_max_len"]                         = m_shower_in_max_len;                         // cm; only read when shower_in_cascade_guard
    cfg["shower_in_mip_hi"]                          = m_shower_in_mip_hi;                          // ratio; only read when shower_in_cascade_guard
    cfg["shower_connect_from_vertices_straight_guard"]  = m_shower_connect_from_vertices_straight_guard;  // doc pr/40 round 9; false = legacy (cross-cluster anchor forced e-)
    cfg["shower_connect_start_seg_straight_guard"]      = m_shower_connect_start_seg_straight_guard;      // doc pr/40 round 9; false = legacy (accept-time set_pdg(11) unconditional)
    cfg["examine_direction_dirsign_shower_in_guard"]    = m_examine_direction_dirsign_shower_in_guard;    // doc pr/40 round 9; false = legacy (no geometry arm beside pr/74 P1)
    cfg["daughter_shower_angle_reclass_straight_guard"] = m_daughter_shower_angle_reclass_straight_guard; // doc pr/40 round 9; false = legacy (angle reclass writes e- unconditionally)
    cfg["shower_topo_reexam_straight_guard"]            = m_shower_topo_reexam_straight_guard;            // doc pr/40 round 9; false = legacy (topo re-exam escape unguarded)
    cfg["sfv_kink_max"]                                 = m_sfv_kink_max;                                 // degrees; continuation-arm tunable
    cfg["shower_nv_bridge_track"]                       = m_shower_nv_bridge_track;                       // doc pr/40 round 9 B2; false = legacy (conn-2 shower, no bridge)
    cfg["shower_nv_bridge_max_gap"]                     = m_shower_nv_bridge_max_gap;                     // cm; only read when shower_nv_bridge_track
    cfg["shower_nv_main_pi_init"]                       = m_shower_nv_main_pi_init;                       // doc pr/97 D1; false = legacy indeterminate main_pi read
    cfg["kine_drop_stray_satellites"]                   = m_kine_drop_stray_satellites;                   // doc pr/92; false = legacy (every conn-2/3 satellite summed into Enu)
    cfg["kine_sat_min_energy"]                          = m_kine_sat_min_energy;                          // MeV; only read when kine_drop_stray_satellites
    cfg["kine_sat_prox_max"]                            = m_kine_sat_prox_max;                            // cm; only read when kine_drop_stray_satellites
    cfg["kine_sat_angle_bad"]                           = m_kine_sat_angle_bad;                           // degrees; only read when kine_drop_stray_satellites
    cfg["kine_sat_angle_main"]                          = m_kine_sat_angle_main;                          // degrees; only read when kine_drop_stray_satellites
    cfg["kine_sat_far_dis"]                             = m_kine_sat_far_dis;                             // cm; only read when kine_drop_stray_satellites
    cfg["kine_sat_axis_dis_cut"]                        = m_kine_sat_axis_dis_cut;                        // cm; only read when kine_drop_stray_satellites
    cfg["kine_sat_cont_kink"]                           = m_kine_sat_cont_kink;                           // degrees; only read when kine_drop_stray_satellites
    cfg["kine_sat_track_max_nseg"]                      = m_kine_sat_track_max_nseg;                      // count; round 2 topology split
    cfg["kine_sat_em_far_dis"]                          = m_kine_sat_em_far_dis;                          // cm; round 2 EM far-drop distance
    cfg["michel_stem_michel_check"]                  = m_michel_stem_michel_check;                  // doc pr/74 round 2 P2; false = legacy (any shower-like sibling passes)
    cfg["michel_stem_max_far_len"]                   = m_michel_stem_max_far_len;                   // cm; only read when michel_stem_michel_check
    cfg["shower_stem_backfill"]                      = m_shower_stem_backfill;                      // doc pr/74 round 2 K4; false = legacy (walked-past stems stay out of showers)
    cfg["stem_backfill_max_len"]                     = m_stem_backfill_max_len;                     // cm; only read when shower_stem_backfill
    cfg["stem_backfill_mip_lo"]                      = m_stem_backfill_mip_lo;                      // ratio; only read when shower_stem_backfill
    cfg["stem_backfill_mip_hi"]                      = m_stem_backfill_mip_hi;                      // ratio; only read when shower_stem_backfill
    cfg["stem_backfill_min_shower_len"]              = m_stem_backfill_min_shower_len;              // cm; only read when shower_stem_backfill
    cfg["shower_conn3_unreachable"]                  = m_shower_conn3_unreachable;                  // doc pr/74 round 2 K5 (pr/65 rung 2); false = legacy (unreachable segments stay PF-invisible)
    cfg["conn3_unreachable_min_len"]                 = m_conn3_unreachable_min_len;                 // cm; only read when shower_conn3_unreachable
    cfg["conn3_stitch_max"]                          = m_conn3_stitch_max;                          // cm; doc pr/84 r2 F3; 0 = off = legacy = byte-identical
    cfg["shower_dedup_start_seg"]                    = m_shower_dedup_start_seg;                    // doc pr/84 r3 S1; false = off = legacy = byte-identical
    cfg["shower_traj_michel_stem"]                   = m_shower_traj_michel_stem;                   // doc pr/74 round 4 K6; false = legacy (a stopping muon + Michel stays one EM shower)
    cfg["michel_stem_traj_min_len"]                  = m_michel_stem_traj_min_len;                  // cm; only read when shower_traj_michel_stem
    cfg["michel_stem_traj_max_len"]                  = m_michel_stem_traj_max_len;                  // cm; only read when shower_traj_michel_stem
    cfg["michel_stem_traj_mip_lo"]                   = m_michel_stem_traj_mip_lo;                   // x MIP median; only read when shower_traj_michel_stem
    cfg["michel_stem_traj_max_far_len"]              = m_michel_stem_traj_max_far_len;              // cm; only read when shower_traj_michel_stem
    cfg["michel_stem_traj_min_kink_deg"]             = m_michel_stem_traj_min_kink_deg;             // deg; only read when shower_traj_michel_stem
    cfg["shower_long_muon_keep_type"]                = m_shower_long_muon_keep_type;                // false = legacy (long-muon pseudo-shower start segment majority-voted to e-)
    cfg["shower_bragg_protect_start_segment"]        = m_shower_bragg_protect_start_segment;        // false = legacy (Bragg-PID-confident muon/proton start segment majority-voted to e-)
    cfg["shower_reclass_case_b_dqdx_guard"]          = m_shower_reclass_case_b_dqdx_guard;          // doc pr/93 Cause A; false = legacy (Case B pdg-11 write unconditional)
    cfg["shower_accept_pid_guard"]                   = m_shower_accept_pid_guard;                   // doc pr/93 Cause B; false = legacy (acceptance-site set_pdg(11) unconditional)
    cfg["shower_pid_guard_min_len"]                  = m_shower_pid_guard_min_len;                  // cm; shared Cause A/B floor, inert while both off
    cfg["shower_vote_track_pid_counts"]              = m_shower_vote_track_pid_counts;              // doc pr/93 Cause C; false = legacy (only confirmed protons count as track)
    cfg["shower_cone_absorb_guard"]               = m_shower_cone_absorb_guard;               // doc pr/93 Cause D; false = legacy (pass-3 direction-cone absorber unguarded)
    cfg["shower_detach_track_stem"]                  = m_shower_detach_track_stem;                  // doc pr/93 r4; false = legacy (track-headed shower keeps its stem)
    cfg["shower_ghost_member_drop"]                  = m_shower_ghost_member_drop;                  // doc pr/99 r2; false = legacy (ghost members stay), byte-identical
    cfg["shower_ghost_overlap_frac"]                 = m_shower_ghost_overlap_frac;                 // 2nd-best per-view overlap gate; inert while drop off (doc pr/99 r2)
    cfg["shower_ghost_dqdx_ratio"]                   = m_shower_ghost_dqdx_ratio;                   // starved gate vs mip median; inert while drop off (doc pr/99 r2)
    cfg["shower_ghost_min_len"]                      = m_shower_ghost_min_len;                      // cm; inert while drop off (doc pr/99 r2)
    cfg["kine_charge_dedup"]                         = m_kine_charge_dedup;                         // doc pr/99 r3 C1; false = legacy ownership-free sum, byte-identical
    cfg["kine_charge_rebuild"]                       = m_kine_charge_rebuild;                       // doc pr/99 r3 C1b; false = legacy add-only clouds, byte-identical
    cfg["kine_charge_track_ctx"]                     = m_kine_charge_track_ctx;                     // doc pr/101 K1; false = showers-only ownership, byte-identical
    cfg["kine_mass_rules"]                           = m_kine_mass_rules;                           // doc pr/101 K2; false = legacy mass branches, byte-identical
    cfg["kine_hadronic_dqdx"]                        = m_kine_hadronic_dqdx;                        // doc pr/101 K3; false = legacy charge best, byte-identical
    cfg["kine_long_muon_mode"]                       = m_kine_long_muon_mode;                       // doc pr/101 K4; 0 = legacy dQdx, byte-identical
    cfg["kine_long_muon_ratio_lo"]                   = m_kine_long_muon_ratio_lo;                   // inert unless mode 2
    cfg["kine_long_muon_ratio_hi"]                   = m_kine_long_muon_ratio_hi;                   // inert unless mode 2
    cfg["long_muon_range_empty_chain_fallback"]      = m_long_muon_range_empty_chain_fallback;      // doc 84 round 1 (P1); false = legacy (chainless muon shower keeps range 0)
    cfg["long_muon_members_geometry"]                = m_long_muon_members_geometry;                // doc 84 round 2; false = legacy (chain-truncated range/endpoint)
    cfg["long_muon_cathode_bridge"]                  = m_long_muon_cathode_bridge;                  // doc 84 round 2; false = legacy (cathode-split muon stays split)
    cfg["long_muon_cathode_bridge_x"]                = m_long_muon_cathode_bridge_x;                // doc 84 round 2, cm; inert unless bridge on
    cfg["long_muon_cathode_bridge_xcut"]             = m_long_muon_cathode_bridge_xcut;             // doc 84 round 2, cm; inert unless bridge on
    cfg["long_muon_cathode_bridge_gap"]              = m_long_muon_cathode_bridge_gap;              // doc 84 round 2, cm; inert unless bridge on
    cfg["long_muon_cathode_bridge_angle"]            = m_long_muon_cathode_bridge_angle;            // doc 84 round 2, deg; inert unless bridge on
    cfg["long_muon_cathode_bridge_lever"]            = m_long_muon_cathode_bridge_lever;            // doc 84 round 4 G1, cm; inert unless bridge on
    cfg["long_muon_cathode_bridge_track_partner"]    = m_long_muon_cathode_bridge_track_partner;    // doc 84 round 4 G2; inert unless bridge on
    cfg["long_muon_cathode_bridge_short_gap"]        = m_long_muon_cathode_bridge_short_gap;        // doc 84 round 4 G3, cm; 0 == off
    cfg["long_muon_cathode_bridge_short_gap_angle"]  = m_long_muon_cathode_bridge_short_gap_angle;  // doc 84 round 4 G3, deg; inert while short_gap 0
    cfg["long_muon_cathode_bridge_short_gap_len"]    = m_long_muon_cathode_bridge_short_gap_len;    // doc 84 round 4 G3, cm; inert while short_gap 0
    cfg["kine_mainvtx_used_guard"]                   = m_kine_mainvtx_used_guard;                   // doc pr/101 K5; false = legacy, byte-identical
    cfg["shower_hadronic_tag"]                       = m_shower_hadronic_tag;                       // doc pr/99 r3 A5; false = legacy (label 11 stays), byte-identical
    cfg["shower_hadronic_min_len"]                   = m_shower_hadronic_min_len;                   // cm; inert while tag off (doc pr/99 r3)
    cfg["shower_hadronic_scan_len"]                  = m_shower_hadronic_scan_len;                  // cm; inert while tag off (doc pr/99 r3)
    cfg["shower_hadronic_bin"]                       = m_shower_hadronic_bin;                       // cm; inert while tag off (doc pr/99 r3)
    cfg["shower_hadronic_r_cyl"]                     = m_shower_hadronic_r_cyl;                     // cm; inert while tag off (doc pr/99 r3)
    cfg["shower_hadronic_r_core"]                    = m_shower_hadronic_r_core;                    // cm; inert while tag off (doc pr/99 r3)
    cfg["shower_hadronic_growth_max"]                = m_shower_hadronic_growth_max;                // ratio; inert while tag off (doc pr/99 r3)
    cfg["shower_hadronic_growth_bragg"]              = m_shower_hadronic_growth_bragg;              // ratio; inert while tag off (doc pr/99 r3)
    cfg["shower_hadronic_bragg_ratio"]               = m_shower_hadronic_bragg_ratio;               // ratio; inert while tag off (doc pr/99 r3)
    cfg["shower_hadronic_stem_ratio"]                = m_shower_hadronic_stem_ratio;                // MIP units; 0 = branch off (doc pr/99 r3)
    cfg["kine_count_orphan_tracks"]                  = m_kine_count_orphan_tracks;                  // doc pr/93 r4; false = legacy (graph-disconnected confident tracks absent from kine)
    cfg["kine_orphan_track_min"]                     = m_kine_orphan_track_min;                     // cm; only read when kine_count_orphan_tracks
    cfg["shower_pass4_best_owner"]                   = m_shower_pass4_best_owner;                   // doc pr/117 r1; false = legacy greedy owner, byte-identical
    cfg["shower_merge_relax"]                        = m_shower_merge_relax;                        // doc pr/117 r1; false = no pass, byte-identical
    cfg["shower_merge_relax_dis"]                    = m_shower_merge_relax_dis;                    // cm; inert while merge_relax off (doc pr/117 r1)
    cfg["shower_merge_relax_angle"]                  = m_shower_merge_relax_angle;                  // deg; inert while merge_relax off (doc pr/117 r1)
    cfg["shower_merge_relax_min_len"]                = m_shower_merge_relax_min_len;                // cm; fragment length floor; inert while merge_relax off (doc pr/117 r1)
    cfg["shower_flank_absorb"]                       = m_shower_flank_absorb;                       // doc pr/117 r1; false = no pass, byte-identical
    cfg["shower_flank_absorb_max_dis"]               = m_shower_flank_absorb_max_dis;               // cm; inert while flank_absorb off (doc pr/117 r1)
    cfg["shower_flank_absorb_max_len"]               = m_shower_flank_absorb_max_len;               // cm; inert while flank_absorb off (doc pr/117 r1)
    cfg["shower_ex1_conn3_body_dis"]                 = m_shower_ex1_conn3_body_dis;                 // doc pr/118 r1; false = start-segment gate, byte-identical
    cfg["shower_merge_relax_continuity"]             = m_shower_merge_relax_continuity;             // doc pr/118 r1; false = legacy merge_relax only, byte-identical
    cfg["shower_merge_relax_cont_frac"]              = m_shower_merge_relax_cont_frac;              // T2 charge-presence fraction; inert while continuity off (doc pr/118 r1)
    cfg["shower_merge_relax_cont_gap"]               = m_shower_merge_relax_cont_gap;               // cm, T2 stub gap; inert while continuity off (doc pr/118 r1)
    cfg["shower_merge_relax_cont_qmed"]              = m_shower_merge_relax_cont_qmed;              // T2 line-charge floor; inert while continuity off (doc pr/118 r1)
    cfg["shower_merge_relax_cont_axis"]              = m_shower_merge_relax_cont_axis;              // deg, T1+T2 axis cone; inert while continuity off (doc pr/118 r1)
    cfg["shower_merge_relax_cont_dmax"]              = m_shower_merge_relax_cont_dmax;              // cm, T2 junction reach; inert while continuity off (doc pr/118 r1)
    cfg["shower_merge_relax_cont_t1_gap"]            = m_shower_merge_relax_cont_t1_gap;            // cm, T1 touching gap; inert while continuity off (doc pr/118 r1)
    cfg["shower_merge_relax_cont_t1_fold"]           = m_shower_merge_relax_cont_t1_fold;           // deg, T1 fold; inert while continuity off (doc pr/118 r1)
    cfg["stem_backfill_back_guard"]                  = m_stem_backfill_back_guard;                  // doc pr/120 r1; false = legacy stem_backfill, byte-identical
    cfg["stem_backfill_back_ang"]                    = m_stem_backfill_back_ang;                    // deg, backward ceiling; inert while guard off (doc pr/120 r1)
    cfg["shower_ex1_walk_em_track_guard"]            = m_shower_ex1_walk_em_track_guard;            // doc pr/120 r1; false = legacy examine_shower_1 walk, byte-identical
    cfg["shower_ex1_walk_em_track_len"]              = m_shower_ex1_walk_em_track_len;              // cm, e- straight-long floor; inert while guard off (doc pr/120 r1)
    cfg["shower_ex1_dedup_rehome"]                   = m_shower_ex1_dedup_rehome;                   // doc pr/121 r1; false = legacy dedup drop, byte-identical
    cfg["shower_pass4_prune_detached"]               = m_shower_pass4_prune_detached;               // doc pr/123 r1; false = no prune pass, byte-identical
    cfg["shower_pass4_prune_gap"]                    = m_shower_pass4_prune_gap;                    // doc pr/123 r1; cm, inert while prune off
    cfg["shower_pass4_prune_gap2"]                   = m_shower_pass4_prune_gap2;                   // doc pr/124 A; cm, 0 = no tier-2 band prune
    cfg["shower_pass4_prune2_ang"]                   = m_shower_pass4_prune2_ang;                   // doc pr/124 A; deg, inert while gap2 = 0
    cfg["shower_pass4_prune2_mdqdx"]                 = m_shower_pass4_prune2_mdqdx;                 // doc pr/124 A; x MIP, inert while gap2 = 0
    cfg["shower_pass3_cone_guard_len"]               = m_shower_pass3_cone_guard_len;               // doc pr/124 C; cm, 0 = no pass3 decline
    cfg["shower_samevtx_track_absorb"]               = m_shower_samevtx_track_absorb;               // doc pr/125; false = no pass, byte-identical
    cfg["shower_samevtx_absorb_gap"]                 = m_shower_samevtx_absorb_gap;                 // doc pr/125; cm, inert while pass off
    cfg["shower_samevtx_absorb_max_len"]             = m_shower_samevtx_absorb_max_len;             // doc pr/125; cm, inert while pass off
    cfg["shower_samevtx_absorb_min_len"]             = m_shower_samevtx_absorb_min_len;             // doc pr/125; cm, inert while pass off
    cfg["shower_satellite_absorb"]                   = m_shower_satellite_absorb;                   // doc pr/125; false = no pass, byte-identical
    cfg["shower_split"] = m_shower_split;   // doc pr/138 B2; false = no pass
    cfg["shower_split_max_valley"] = m_shower_split_max_valley;   // doc pr/138 B2; sec A5.4 knee
    cfg["shower_split_min_frac"] = m_shower_split_min_frac;   // doc pr/138 B2; per-seed charge share floor
    cfg["shower_split_max_parts"] = m_shower_split_max_parts;   // doc pr/138 B3; 2 = the measured-exact kernel
    cfg["shower_split_min_charge"] = m_shower_split_min_charge;   // doc pr/138 B1; candidate charge floor (raw Fit::dQ)
    cfg["shower_split_min_nseg"] = m_shower_split_min_nseg;   // doc pr/138 B1; candidate member-count floor
    cfg["shower_split_bundle_gap"] = m_shower_split_bundle_gap;   // doc pr/138 B3; cm, inert while the pass is off
    cfg["shower_split_snap"] = m_shower_split_snap;   // doc pr/138 B3; k>=3 bundle dominance floor
    cfg["shower_split_skip_shared"] = m_shower_split_skip_shared;   // doc pr/139 P1.1; refuse a component holding a segment another shower also owns
    cfg["shower_split_max_impact"] = m_shower_split_max_impact;   // doc pr/139 P1.2; cm; 0 = no bound
    cfg["shower_split_em_start"] = m_shower_split_em_start;   // doc pr/139 P1.3; seed the daughter on its nearest EM-typed member
    cfg["shower_split_rehome"] = m_shower_split_rehome;   // doc pr/139 P1.4; offer an orphan daughter to the nearest larger EM shower
    cfg["shower_split_rehome_gap"] = m_shower_split_rehome_gap;   // doc pr/139 P1.4; cm; max daughter->host 3-D gap
    cfg["shower_satellite_absorb_max_mev"]           = m_shower_satellite_absorb_max_mev;           // doc pr/125; MeV, inert while pass off
    cfg["shower_satellite_absorb_host_mev"]          = m_shower_satellite_absorb_host_mev;          // doc pr/125; MeV, inert while pass off
    cfg["shower_pass4_track_guard_len"]              = m_shower_pass4_track_guard_len;              // doc pr/123 r1; cm, 0 = no guard, byte-identical
    cfg["shower_pass4_prox_guard_len"]               = m_shower_pass4_prox_guard_len;               // doc pr/130 item 1b; cm, 0 = unguarded, byte-identical
    cfg["shower_pass3_backfill_guard_len"]           = m_shower_pass3_backfill_guard_len;           // doc pr/130 item 1b; cm, 0 = unguarded, byte-identical
    cfg["stem_backfill_back_dvtx"]                   = m_stem_backfill_back_dvtx;                   // doc pr/130 item B; cm, 0 = off, byte-identical
    cfg["shower_pass4_prefilter_v1_escape"]          = m_shower_pass4_prefilter_v1_escape;          // doc pr/136 r2; false = legacy pre-filter, byte-identical
    cfg["shower_pass4_prefilter_v1_max_v2"]          = m_shower_pass4_prefilter_v1_max_v2;          // doc pr/136 r2; deg, 0 = no ceiling, inert while escape off
    cfg["shower_pass4_prefilter_v1_max_dis"]         = m_shower_pass4_prefilter_v1_max_dis;         // doc pr/136 r3; cm, 0 = no proximity bound, inert while escape off
    cfg["pi0_mass_offset"]                           = m_pi0_mass_offset;                           // doc pr/132 K1; MeV, 10 = legacy, byte-identical
    cfg["pi0_assoc_angle_deg"]                       = m_pi0_assoc_angle_deg;                       // doc pr/132 K2; deg, 30 = legacy, byte-identical
    cfg["pi0_attached_partner_min_mev"]              = m_pi0_attached_partner_min_mev;              // doc pr/132 K3; MeV, 0 = no guard, byte-identical
    cfg["pi0_nv_allow_type2"]                        = m_pi0_nv_allow_type2;                        // doc pr/132 K4; false = legacy pool, byte-identical
    cfg["pi0_nv_max_prongs"]                         = m_pi0_nv_max_prongs;                         // doc pr/132 K5; 2 = legacy gate, byte-identical
    cfg["pi0_readmit_retyped"]                       = m_pi0_readmit_retyped;                       // doc pr/132 K7; false = legacy exclusion, byte-identical
    cfg["pi0_admit_type3"]                           = m_pi0_admit_type3;                           // doc pr/132 K8; false = legacy pool, byte-identical
    cfg["pi0_crumb_assoc_mev"]                       = m_pi0_crumb_assoc_mev;                       // doc pr/132 K9; 0 = legacy angle test, byte-identical
    cfg["pi0_collinear_merge_deg"]                   = m_pi0_collinear_merge_deg;                   // doc pr/132 K12; 0 = legacy pairing, byte-identical
    cfg["pi0_nv_partner_min_mev"]                    = m_pi0_nv_partner_min_mev;                    // doc pr/132 K13; 0 = no floor, byte-identical
    cfg["pi0_nv_retry_paired"]                       = m_pi0_nv_retry_paired;                       // doc pr/132 K14; false = legacy early return, byte-identical
    cfg["pi0_reseat_start_assoc"]                    = m_pi0_reseat_start_assoc;                    // doc pr/132 K15; false = fit-cloud starts, byte-identical
    cfg["shower_em_collinear_deg"]                   = m_shower_em_collinear_deg;                   // doc pr/132 K16; 0 = no merge, byte-identical
    cfg["shower_em_collinear_dis_cm"]                = m_shower_em_collinear_dis_cm;                // doc pr/132 K16; inert while deg = 0
    cfg["shower_em_collinear_host_mev"]              = m_shower_em_collinear_host_mev;              // doc pr/132 K16; inert while deg = 0
    cfg["shower_em_backext_perp_cm"]                 = m_shower_em_backext_perp_cm;                 // doc pr/132 K17; 0 = no back-extension, byte-identical
    cfg["shower_em_backext_len_cm"]                  = m_shower_em_backext_len_cm;                  // doc pr/132 K17; inert while perp = 0
    cfg["pi0_accept_merge_dis_cm"]                   = m_pi0_accept_merge_dis_cm;                   // doc pr/132 K18; 0 = no acceptance merge, byte-identical
    cfg["pi0_bp_vertex_miss_cm"]                     = m_pi0_bp_vertex_miss_cm;                     // doc pr/132 K19; 0 = no NC vertex proposer, byte-identical
    cfg["pi0_admit_muon_showers"]                    = m_pi0_admit_muon_showers;                    // doc pr/133 K20; false = legacy mu exclusion, byte-identical
    cfg["pi0_nc_sig_angle_deg"]                      = m_pi0_nc_sig_angle_deg;                      // doc pr/133 K21; 0 = v3 gate, byte-identical
    cfg["pi0_nc_floor_mev"]                          = m_pi0_nc_floor_mev;                          // doc pr/133 K21 v2; 0 = legacy 20 MeV floor, byte-identical
    cfg["pi0_nc_pf_assoc_deg"]                       = m_pi0_nc_pf_assoc_deg;                       // doc pr/133 K21 v2.2; 0 = no post-fire PF update, byte-identical
    cfg["pi0_nc_frag_merge"]                         = m_pi0_nc_frag_merge;                         // doc pr/134 K22; false = fragment-level bp pairing, byte-identical
    cfg["pi0_pf_assoc_deg"]                          = m_pi0_pf_assoc_deg;                          // doc pr/134 K23; 0 = no P1 post-accept PF absorb, byte-identical
    cfg["pi0_prefer_main_vertex"]                    = m_pi0_prefer_main_vertex;                    // doc pr/134 K24; false = legacy association + ranking, byte-identical
    cfg["pi0_nv_max_vtx_shift_cm"]                   = m_pi0_nv_max_vtx_shift_cm;                   // doc pr/132 K10; 0 = no cap, byte-identical
    cfg["pi0_nv_mass_window_mev"]                    = m_pi0_nv_mass_window_mev;                    // doc pr/132 K11; 60 = legacy window, byte-identical
    cfg["kine_count_guard_freed"]                    = m_kine_count_guard_freed;                    // doc pr/123 r2; false = guard-freed tracks uncounted, byte-identical
    cfg["kine_guard_freed_impact"]                   = m_kine_guard_freed_impact;                   // doc pr/129; 0 = pointing test off, byte-identical
    cfg["kine_guard_freed_miss_deg"]                 = m_kine_guard_freed_miss_deg;                 // doc pr/129
    cfg["kine_count_near_cross_cluster"]             = m_kine_count_near_cross_cluster;             // doc pr/128; false = near cross-cluster tracks uncounted, byte-identical
    cfg["kine_near_gap"]                             = m_kine_near_gap;                             // cm; only read when kine_count_near_cross_cluster
    cfg["kine_near_min_len"]                         = m_kine_near_min_len;                         // cm; only read when kine_count_near_cross_cluster
    cfg["kine_near_end_tol"]                         = m_kine_near_end_tol;                         // cm; only read when kine_count_near_cross_cluster
    cfg["kine_near_kink_deg"]                        = m_kine_near_kink_deg;                        // deg; only read when kine_count_near_cross_cluster
    cfg["kine_count_conn4_near"]                     = m_kine_count_conn4_near;                     // doc pr/128; false = every conn-4 shower uncounted, byte-identical
    cfg["kine_conn4_near_gap"]                       = m_kine_conn4_near_gap;                       // cm; only read when kine_count_conn4_near
    cfg["straight_cont_cross_cluster"]               = m_straight_cont_cross_cluster;               // doc pr/93 r4; false = legacy (no cross-cluster continuation demotion)
    cfg["sccc_bridge_body"]                          = m_sccc_bridge_body;                          // doc pr/93 r4; false = demote-only (no bridge replay)
    cfg["sccc_max_gap"]                              = m_sccc_max_gap;                              // cm; base tier, only read when straight_cont_cross_cluster
    cfg["sccc_kink_max"]                             = m_sccc_kink_max;                             // degrees; base tier
    cfg["sccc_gap_aligned"]                          = m_sccc_gap_aligned;                          // cm; aligned tier (tight collinearity buys reach)
    cfg["sccc_kink_tight"]                           = m_sccc_kink_tight;                           // degrees; aligned tier
    cfg["single_muon_proton_chain_veto"]             = m_single_muon_proton_chain_veto;             // false = legacy (1-hop proton veto only)
    cfg["single_muon_long_muon_claim"]               = m_single_muon_long_muon_claim;               // false = legacy (long-muon chain never claims the vertex muon slot)
    cfg["pid_flag_reconcile"]                        = m_pid_flag_reconcile;                        // false = legacy (no late reconciliation pass)
    cfg["long_muon_stub_bridge"]                     = m_long_muon_stub_bridge;                     // false = legacy (stub-blocked long-muon chains never form)
    cfg["long_muon_stub_bridge_len"]                 = m_long_muon_stub_bridge_len;                 // doc 84 round 1 (P3), cm; 6.0 = legacy literal
    cfg["long_muon_angle_relax_long"]                = m_long_muon_angle_relax_long;                // doc 84 round 1 (P2); false = legacy 10 deg walk
    cfg["long_muon_angle_relax_deg"]                 = m_long_muon_angle_relax_deg;                 // doc 84 round 1 (P2), deg; inert unless relax on


    return cfg;
}

namespace {

// doc 84 round 2 (long_muon_cathode_bridge).  A muon that crosses the SBND
// cathode can reconstruct as two disconnected PR halves: the near half is a
// |13|-typed shower whose member trajectory ends at |x| < ~5 cm of the seam,
// and the far half is either a second muon-typed shower (mcp2k 53793:
// 367 + 528 MeV) or a bare not-in-any-shower segment chain (177536: 98.1 cm;
// 77978: 23.5 cm + a 9.8 cm Bragg stub typed 2212).  The halves do NOT share
// a PR cluster_id (53793: 37 vs 12), so the guards are purely geometric:
// both facing ends within `xcut` of the cathode plane, far ends on opposite
// drift sides (near ends may jitter onto one side -- 77978's far half itself
// crosses x=0), 3D gap < `gap`, and the continuation direction leaving the
// muon end within `angle` of BOTH the gap vector and the partner's tangent.
// Partner showers must themselves be |13|-typed (392901's 106 MeV EM shower
// across the seam is excluded by design); bare partners are absorbed as a
// BFS over connected not-in-any-shower segments and retyped to 13 so the
// members-geometry accumulators count them.  Only called when the knob is
// on => byte-identical off.
struct CathodeBridgeCfg {
    double x;       // cathode plane [internal units]
    double xcut;    // end-to-plane admission window
    double gap;     // max 3D gap between facing ends
    double angle;   // max continuation angle [deg]
    // doc 84 round 4 -- all default to the round-2 behaviour, so the pass is
    // byte-identical until one is moved off its legacy value.
    double lever;           // end-direction lever; 5cm == the round-2 hardcode (G1)
    bool   track_partner;   // admit a |211|-typed partner SHOWER too (G2)
    double short_gap;       // 0 == off; below this gap the gap-vector angle is waived (G3)
    double short_gap_angle; // partner-direction cap [deg] required when the waiver applies
    double short_gap_len;   // min partner length required when the waiver applies
};

struct CathodeEnd {
    ShowerPtr shower;      // nullptr for a bare segment
    SegmentPtr seg;
    WireCell::Point p;     // the end near the cathode
    WireCell::Point far_p; // the other end of the segment
    WireCell::Vector into; // unit vector from p into the segment
};

inline double cb_angle_deg(const WireCell::Vector& a, const WireCell::Vector& b)
{
    const double ma = a.magnitude(), mb = b.magnitude();
    if (ma <= 0 || mb <= 0) return 181.0;   // degenerate: fails any cut
    double c = a.dot(b) / (ma * mb);
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c) * 180.0 / 3.141592653589793;
}

// bridged_out (doc 84 round 3, may be null): collects the muon-typed
// segments this pass absorbs -- the merged partner's members in the
// shower branch, the retyped BFS chain in the bare branch -- so the MCS
// driver can optionally fit the full track (mcs_bridged_members).  Filling
// it moves no output bytes by itself.
int long_muon_cathode_bridge_pass(PatternAlgorithms& pattern_algos, Graph& graph,
    VertexPtr main_vertex, IndexedShowerSet& showers,
    ShowerSegmentMap& map_segment_in_shower,
    const Clus::ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model,
    const CathodeBridgeCfg& cfg, IndexedSegmentSet* bridged_out,
    std::shared_ptr<spdlog::logger> log)
{
    const double min_len = 5 * units::cm;   // ignore sub-5cm fragments on either side
    // doc 84 round 4 (G1): the end-direction lever, hardcoded 5cm in round 2.
    // segment_cal_dir_3vector is a CENTROID direction -- it averages every fit
    // point within `lever` of the end -- and it is evaluated exactly where the
    // cathode charge loss is worst, so a 5cm baseline can be badly bent while
    // the segment as a whole is straight (172794: both angle tests fail at
    // 40.5/38.1 deg against a corrupted reference, while the gap vector and
    // the partner direction agree with EACH OTHER to 8 deg).
    const double lever = cfg.lever;

    // Member lists per shower.  map_segment_in_shower is index-ordered
    // (SegmentIndexCmp), the grouping map is ShowerIndexCmp -- deterministic.
    std::map<ShowerPtr, std::vector<SegmentPtr>, ShowerIndexCmp> members;
    for (auto& [seg, sh] : map_segment_in_shower) {
        if (seg && sh) members[sh].push_back(seg);
    }
    auto muon_member_length = [&](const ShowerPtr& sh) {
        double tot = 0;
        auto it = members.find(sh);
        if (it == members.end()) return tot;
        for (auto& seg : it->second) {
            if (seg && seg->has_particle_info() &&
                std::abs(seg->particle_info()->pdg()) == 13)
                tot += segment_track_length(seg);
        }
        return tot;
    };

    auto collect_ends = [&](SegmentPtr seg, ShowerPtr sh, std::vector<CathodeEnd>& out) {
        if (!seg) return;
        if (segment_track_length(seg) < min_len) return;
        const auto& fits = seg->fits();
        if (fits.size() < 2) return;
        const WireCell::Point ends[2] = { fits.front().point, fits.back().point };
        for (int e = 0; e < 2; ++e) {
            WireCell::Point p = ends[e];
            if (std::abs(p.x() - cfg.x) >= cfg.xcut) continue;
            WireCell::Vector into = segment_cal_dir_3vector(seg, p, lever);
            if (into.magnitude() <= 0) continue;
            out.push_back({sh, seg, p, ends[1 - e], into});
        }
    };

    // Receiver ends: muon-typed member segments of |13|-typed showers.
    std::vector<CathodeEnd> mu_ends;
    for (auto& sh : showers) {
        if (!sh || std::abs(sh->get_particle_type()) != 13) continue;
        auto it = members.find(sh);
        if (it == members.end()) continue;
        for (auto& seg : it->second) {
            if (!seg || !seg->has_particle_info()) continue;
            if (std::abs(seg->particle_info()->pdg()) != 13) continue;
            collect_ends(seg, sh, mu_ends);
        }
    }
    if (mu_ends.empty()) return 0;

    // Partner ends: bare (not-in-any-shower) segments, or members of OTHER
    // |13|-typed showers.  ordered_edges => deterministic.
    std::vector<CathodeEnd> partner_ends;
    for (auto ed : WireCell::Clus::PR::ordered_edges(graph)) {
        SegmentPtr seg = graph[ed].segment;
        if (!seg) continue;
        auto mit = map_segment_in_shower.find(seg);
        ShowerPtr osh = (mit != map_segment_in_shower.end()) ? mit->second : nullptr;
        // doc 84 round 4 (G2): a cathode-split muon's near-seam stub is often
        // mis-PID'd -- 347890's far half is a |13| shower but its facing
        // partner is a 27cm shower typed 211 (pi+), geometrically cleaner
        // (gap 14.3cm, angles 15.4/6.2 deg) than three of the four round-2
        // successes.  track_partner admits |211| as well.  EM (11/22) is NOT
        // admitted and must not be: absorbing a genuine EM shower across the
        // seam is the failure mode this guard exists for (392901's far half is
        // a 106 MeV electron shower with comparably good geometry).
        if (osh) {
            const int optype = std::abs(osh->get_particle_type());
            const bool ok = (optype == 13) || (cfg.track_partner && optype == 211);
            if (!ok) continue;
        }
        collect_ends(seg, osh, partner_ends);
    }
    if (partner_ends.empty()) return 0;

    std::set<ShowerPtr, ShowerIndexCmp> gone;   // showers absorbed by a merge
    std::set<size_t> taken;                     // graph indices of absorbed bare segments
    int n_bridged = 0;

    for (auto& me : mu_ends) {
        if (gone.count(me.shower)) continue;
        const CathodeEnd* best = nullptr;
        double best_gap = cfg.gap;
        for (auto& pe : partner_ends) {
            if (pe.seg == me.seg) continue;
            if (pe.shower && pe.shower == me.shower) continue;
            if (pe.shower && gone.count(pe.shower)) continue;
            if (!pe.shower && taken.count(pe.seg->get_graph_index())) continue;
            if ((me.far_p.x() - cfg.x) * (pe.far_p.x() - cfg.x) > 0) continue;
            WireCell::Vector gapv = pe.p - me.p;
            const double gap = gapv.magnitude();
            if (gap >= best_gap) continue;
            const WireCell::Vector cont = -1.0 * me.into;
            const double a_gap = cb_angle_deg(cont, gapv);
            const double a_tan = cb_angle_deg(cont, pe.into);
            // doc 84 round 4 (G3): the gap vector's angular precision goes as
            // ~atan(sigma_endpoint / gap), so at a short gap it carries no
            // information -- 67026's two halves are 4 deg collinear over
            // 257.9cm of partner but sit ~4.9cm laterally apart, putting the
            // 5.6cm gap vector 78 deg off both tracks.  Waive the gap-vector
            // test there, with the partner-direction test tightened and a
            // minimum partner length standing in for the lost constraint.
            //
            // The waiver must NEVER apply to a degenerate gap: cb_angle_deg
            // returns 181 for a zero-length gap vector, and that is the ONLY
            // thing rejecting a partner that shares a vertex with the muon
            // (407798: gap exactly 0, a normal graph junction inside one
            // cluster, owner-confirmed correct reject).  Hence gap > 0 is a
            // hard precondition, not a consequence of the other two gates.
            const bool waive_gap_angle =
                cfg.short_gap > 0 && gap > 0 && gap <= cfg.short_gap &&
                a_tan <= cfg.short_gap_angle &&
                segment_track_length(pe.seg) >= cfg.short_gap_len;
            const char* reject = nullptr;
            if (!waive_gap_angle && a_gap > cfg.angle) reject = "angle_gap";
            else if (a_tan > cfg.angle) reject = "angle_tan";
            if (reject) {
                // Round 4: the pass logged accepted bridges only and was
                // silent on rejects, so a non-firing bridge could not be
                // diagnosed without re-deriving the geometry outside the
                // toolkit.  One line per geometrically-reachable candidate.
                SPDLOG_LOGGER_DEBUG(log,
                    "long_muon_cathode_bridge: reject {} sid={} seg={} partner={} "
                    "partner_kind={} partner_len={:.1f}cm gap={:.1f}cm a_gap={:.1f} "
                    "a_tan={:.1f} lever={:.1f}cm ends x=({:.1f},{:.1f})cm",
                    reject, me.shower->get_shower_id(), me.seg->get_graph_index(),
                    pe.seg->get_graph_index(), pe.shower ? "shower" : "bare",
                    segment_track_length(pe.seg) / units::cm, gap / units::cm,
                    a_gap, a_tan, cfg.lever / units::cm,
                    me.p.x() / units::cm, pe.p.x() / units::cm);
                continue;
            }
            best = &pe;
            best_gap = gap;
        }
        if (!best) continue;

        if (best->shower) {
            // Shower partner: decide the keeper -- main-vertex attachment
            // first, then the longer muon membership, then the smaller id.
            ShowerPtr a = me.shower, b = best->shower;
            // doc 84 round 4 (G2): a partner admitted only by track_partner is
            // a mis-PID'd stub, not a peer, so the muon shower must stay the
            // keeper and the stub must be retyped.  Letting the rank contest
            // run would hand 347890 to the 211-typed pi+ shower -- it is the
            // root-attached one, so it wins the first rank key -- relabelling
            // a 177.7cm muon as a pion and diverting the merged object out of
            // calculate_kinematics_long_muon.  Unreachable while
            // track_partner is off: E5 admits no non-|13| partner shower then.
            const bool track_only = std::abs(b->get_particle_type()) != 13;
            ShowerPtr keep, drop;
            if (track_only) {
                keep = a;
                drop = b;
                // Retype the absorbed stub to 13, mirroring what the
                // bare-chain branch does below, so its length reaches the
                // range/endpoint accumulators, bridged_out (MCS) and the PF
                // label as muon rather than staying a pion beside one.
                for (auto& seg : members[drop]) {
                    if (!seg->has_particle_info() ||
                        std::abs(seg->particle_info()->pdg()) != 13) {
                        auto mom = segment_cal_4mom(seg, 13, particle_data, recomb_model,
                                                    pattern_algos.m_mip_dqdx);
                        seg->particle_info(std::make_shared<Aux::ParticleInfo>(
                            13, particle_data->get_particle_mass(13),
                            particle_data->pdg_to_name(13), mom));
                    }
                }
                drop->set_particle_type(13);
            }
            else {
                auto rank = [&](const ShowerPtr& s) {
                    return std::make_tuple(s->start_vertex() == main_vertex ? 0 : 1,
                                           -muon_member_length(s),
                                           s->get_shower_id());
                };
                keep = (rank(a) <= rank(b)) ? a : b;
                drop = (keep == a) ? b : a;
            }
            SPDLOG_LOGGER_DEBUG(log,
                "long_muon_cathode_bridge: merge shower sid={} (mu_len={:.1f}cm) <- sid={} "
                "(mu_len={:.1f}cm) gap={:.1f}cm ends x=({:.1f},{:.1f})cm",
                keep->get_shower_id(), muon_member_length(keep) / units::cm,
                drop->get_shower_id(), muon_member_length(drop) / units::cm,
                best_gap / units::cm, me.p.x() / units::cm, best->p.x() / units::cm);
            if (bridged_out) {
                // The dropped partner's muon-typed members are the piece MCS
                // never saw; enumerable only BEFORE add_shower folds them in.
                for (const auto& seg : members[drop]) {
                    if (seg->has_particle_info() &&
                        std::abs(seg->particle_info()->pdg()) == 13) {
                        bridged_out->insert(seg);
                    }
                }
            }
            keep->add_shower(*drop);
            showers.erase(drop);
            gone.insert(drop);
            keep->set_flag_kinematics(false);
        }
        else {
            // Bare partner: absorb the connected not-in-any-shower chain,
            // retyping to 13 so range/endpoint accumulators count it (the
            // 77978 far half carries a 9.8 cm Bragg stub typed 2212 beyond
            // the 23.5 cm muon-typed piece).
            std::vector<SegmentPtr> chain;
            std::set<size_t> seen;
            std::vector<SegmentPtr> todo{best->seg};
            seen.insert(best->seg->get_graph_index());
            while (!todo.empty()) {
                SegmentPtr cur = todo.back();
                todo.pop_back();
                chain.push_back(cur);
                auto [va, vb] = find_vertices(graph, cur);
                for (auto* vp : {&va, &vb}) {
                    VertexPtr v = *vp;
                    if (!v || !v->descriptor_valid()) continue;
                    // Never traverse THROUGH the interaction point: other
                    // prongs at the main vertex are their own particles, not
                    // the muon's far half (77978/177536: a ~10-24 cm proton
                    // prong shares the vertex with the muon and must stay a
                    // separate PF node, not be folded into the muon).
                    if (v == main_vertex) continue;
                    for (auto ed : sorted_out_edges(v->get_descriptor(), graph)) {
                        SegmentPtr nx = graph[ed].segment;
                        if (!nx || seen.count(nx->get_graph_index())) continue;
                        if (map_segment_in_shower.count(nx)) continue;   // stays bare-only
                        if (taken.count(nx->get_graph_index())) continue;
                        seen.insert(nx->get_graph_index());
                        todo.push_back(nx);
                    }
                }
                if (chain.size() > 64) break;   // safety valve; never seen in census
            }
            double absorbed_len = 0;
            for (auto& seg : chain) {
                if (!seg->has_particle_info() ||
                    std::abs(seg->particle_info()->pdg()) != 13) {
                    auto mom = segment_cal_4mom(seg, 13, particle_data, recomb_model,
                                                pattern_algos.m_mip_dqdx);
                    seg->particle_info(std::make_shared<Aux::ParticleInfo>(
                        13, particle_data->get_particle_mass(13),
                        particle_data->pdg_to_name(13), mom));
                }
                me.shower->add_segment(seg, true);
                taken.insert(seg->get_graph_index());
                if (bridged_out) bridged_out->insert(seg);
                absorbed_len += segment_track_length(seg);
            }
            // When the absorbed chain reaches the main vertex the muon's true
            // start IS the interaction point: re-seat the shower there so the
            // kinematics start_point / init_dir and the Bee PF tree connect
            // nu -> mu instead of leaving the vertex marker dangling on the
            // far side (77978: nu at x=-7.5 with the muon drawn from +4.8).
            SegmentPtr vtx_seg = nullptr;
            for (auto& seg : chain) {
                auto [sa, sb] = find_vertices(graph, seg);
                if (sa == main_vertex || sb == main_vertex) { vtx_seg = seg; break; }
            }
            if (vtx_seg) {
                me.shower->set_start_vertex(main_vertex, 1);
                me.shower->set_start_segment(vtx_seg, true);
            }
            SPDLOG_LOGGER_DEBUG(log,
                "long_muon_cathode_bridge: absorb bare chain into sid={} nseg={} len={:.1f}cm "
                "gap={:.1f}cm ends x=({:.1f},{:.1f})cm reseat={}",
                me.shower->get_shower_id(), chain.size(), absorbed_len / units::cm,
                best_gap / units::cm, me.p.x() / units::cm, best->p.x() / units::cm,
                vtx_seg ? 1 : 0);
            me.shower->set_flag_kinematics(false);
        }
        ++n_bridged;
    }
    return n_bridged;
}

}  // namespace

void TaggerCheckNeutrino::visit(Ensemble& ensemble) const
{
    using Clock = std::chrono::steady_clock;
    using MS    = std::chrono::duration<double, std::milli>;
    auto t_total = Clock::now();
    auto t0 = Clock::now();

    // The member fitter outlives the event.  Everything it cached last event
    // points into the previous event's (destroyed) Points tree, so drop it
    // before touching anything -- once per event, never per candidate (the
    // candidate loop below deliberately reuses this fitter for candidate 0).
    // No-op on the first event, so a one-event process is unchanged.
    m_track_fitter->reset_for_new_event();

    // Configure the track fitter with detector volume
    m_track_fitter->set_detector_volume(m_dv);
    m_track_fitter->set_pc_transforms(m_pcts); 

    // Get the specified grouping (default: "live")
    auto groupings = ensemble.with_name(m_grouping_name);
    if (groupings.empty()) {
        return;
    }
    
    auto& grouping = *groupings.at(0);
    
    // ---- doc pr/94: neutrino-PR candidates -------------------------- //
    // A candidate is one main activity to reconstruct, the companion
    // (associated) clusters that share its flash bundle, and -- for the
    // per-bundle mode -- a record of EVERY main activity the cosmic taggers
    // evaluated inside that bundle together with their verdicts.  Legacy mode
    // builds at most one candidate from the single event-wide winner selected
    // by the untouched block below; per-bundle mode builds one per in-beam
    // -window bundle.  Both then run the identical PR chain, once each.
    struct ActivityRec {
        int   cluster_id{-1};
        float length_cm{0.f};
        int   is_selected{0};
        int   is_demoted{0};
        int   tgm{0};
        int   stm{0};
        int   fc{0};
        int   lm{-1};
        int   evaluated{0};
    };
    struct NuCandidate {
        Cluster* main{nullptr};
        std::vector<Cluster*> others;
        int gid{-1};
        std::vector<ActivityRec> acts;
    };
    std::vector<NuCandidate> candidates;

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
    else if (!m_nu_per_bundle) {
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
        // Round 3 (sbnd_xin/docs/73 sec 12): demoted-main fallback.  The
        // taggers can see demoted mains (evaluate_demoted_mains, doc pr/20
        // Part I P3) but this selector could not, so an event whose ONLY
        // in-window main was convicted -- e.g. a cathode-rescue join that
        // legitimately became a TGM/STM -- ended with no candidate even when
        // its bundle still held an examined, untagged former main (SBND evt
        // 65289: merged main 13 is an STM; demoted main 18 scored TGM=false
        // STM=0 by the taggers and was never a candidate).  Runs ONLY when
        // the loop above selected nothing, so any event with a surviving
        // main-cluster candidate is byte-identical.  Same gates as the
        // primary loop: beam window, nu_skip_cosmic verdicts, the bundle
        // veto with the pr/16 design-A min-length guard, longest-wins.  The
        // winner keeps its flags (still associated_cluster + demoted_main):
        // downstream PR works off this local pointer, and re-flagging would
        // change what the bundle-veto set of a later visit sees.
        if (!main_cluster && m_nu_fallback_demoted_mains) {
            int n_demoted = 0;
            for (auto* cluster : grouping.children()) {
                if (cluster->get_flag(Flags::main_cluster)) continue;
                if (!cluster->get_flag(Flags::demoted_main)) continue;
                n_demoted++;
                const double t0 = cluster->get_cluster_t0();
                if (t0 < m_beam_window_low || t0 >= m_beam_window_high) continue;
                if (m_nu_skip_cosmic) {
                    const bool tgm = cluster->get_flag(Flags::TGM);
                    const bool stm = cluster->get_flag(Flags::STM);
                    const int lm = cluster->get_scalar<int>("lm_flag", -1);
                    if (tgm || stm || lm > 0) {
                        SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: in-window demoted main {} (t0 {:.3f} us, L {:.1f} cm) cosmic-tagged (TGM={} STM={} lm_flag={}); skipping (nu_fallback_demoted_mains)",
                                           cluster->get_cluster_id(), t0/units::us, cluster->get_length()/units::cm,
                                           tgm, stm, lm);
                        continue;
                    }
                    if (!cosmic_gids.empty()) {
                        const int gid = cluster->get_scalar<int>("matched_flash_gid", -1);
                        if (gid >= 0 && cosmic_gids.count(gid)) {
                            // Same guard, same reason: an examined, untagged,
                            // long demoted main survives its convicted
                            // bundle-mate; an unexamined shard does not.
                            const double min_len = m_nu_skip_cosmic_bundle_min_length * units::cm;
                            if (!(min_len > 0 && cluster->get_length() >= min_len)) {
                                SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: in-window demoted main {} (t0 {:.3f} us, L {:.1f} cm) shares flash bundle gid {} with a cosmic-tagged main and is under the {:.1f} cm guard; skipping (nu_fallback_demoted_mains)",
                                                   cluster->get_cluster_id(), t0/units::us,
                                                   cluster->get_length()/units::cm, gid,
                                                   m_nu_skip_cosmic_bundle_min_length);
                                continue;
                            }
                        }
                    }
                }
                if (!main_cluster || cluster->get_length() > main_cluster->get_length()) {
                    main_cluster = cluster;
                }
            }
            if (main_cluster) {
                SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: no main-cluster candidate; selected demoted main {} (t0 {:.3f} us, L {:.1f} cm) of {} demoted (nu_fallback_demoted_mains)",
                                   main_cluster->get_cluster_id(), main_cluster->get_cluster_t0()/units::us,
                                   main_cluster->get_length()/units::cm, n_demoted);
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
    else {
        // ---- doc pr/94: per-bundle candidate selection --------------- //
        // One candidate per in-beam-window flash bundle instead of one winner
        // per event.  Inside a bundle the rule is today's rule -- longest
        // untagged main, then the untagged demoted-main fallback -- and the
        // PER-MAIN cosmic veto is deliberately KEPT: a TGM/STM/LM-convicted
        // activity is not a neutrino candidate, which is what leaves an
        // all-cosmic bundle with no row at all.  What is deliberately NOT
        // applied is the event-level `cosmic_gids` BUNDLE veto: it exists only
        // to stop a convicted bundle supplying the single event-wide winner,
        // and it is exactly what would discard the clean 198.9 cm demoted main
        // 21 of SBND 18255/395148 because its bundle-mate, the 533 cm STM muon
        // 10, was convicted.  Sibling conviction no longer vetoes anything;
        // every evaluated activity reports its own verdict in its act_* slot.
        //
        // Determinism: gids are int-keyed and visited in sorted order, and
        // within a gid clusters are visited in grouping.children() order with a
        // length tie-break.  No pointer-keyed container is iterated.
        // doc pr/94 Phase 5b round 2: the legacy event-wide winner, recomputed
        // here purely so the length floor below can EXEMPT it.  This is a
        // deliberate full duplication of the legacy selector at :999-1135 (M10:
        // the production branch stays textually untouched) with the logging and
        // the accumulator side effects stripped -- it must stay in step with
        // that code, and the primary-row gate is what proves it has.
        //
        // Why the exemption cannot be replaced by any length rule: mcp1k evt
        // 62583 keeps a 1.6 cm row and evt 391854 must lose a 1.7 cm row.  The
        // only thing separating them is that 62583's IS the legacy selection
        // and 391854's is not.  Exempting it makes additivity structural: the
        // row the legacy chain reports can never be floored away.
        Cluster* legacy_main = nullptr;
        {
            std::set<int> cg;   // int-keyed => stable iteration
            if (m_nu_skip_cosmic && m_nu_skip_cosmic_bundle) {
                for (auto* c : grouping.children()) {
                    if (!c->get_flag(Flags::main_cluster)) continue;
                    const double t0c = c->get_cluster_t0();
                    if (t0c < m_beam_window_low || t0c >= m_beam_window_high) continue;
                    const int gid = c->get_scalar<int>("matched_flash_gid", -1);
                    if (gid < 0) continue;
                    if (c->get_flag(Flags::TGM) || c->get_flag(Flags::STM)
                        || c->get_scalar<int>("lm_flag", -1) > 0) {
                        cg.insert(gid);
                    }
                }
            }
            const double bundle_min_len = m_nu_skip_cosmic_bundle_min_length * units::cm;
            auto survives = [&](Cluster* c) {
                if (!m_nu_skip_cosmic) return true;
                if (c->get_flag(Flags::TGM) || c->get_flag(Flags::STM)
                    || c->get_scalar<int>("lm_flag", -1) > 0) return false;
                if (!cg.empty()) {
                    const int gid = c->get_scalar<int>("matched_flash_gid", -1);
                    if (gid >= 0 && cg.count(gid)
                        && !(bundle_min_len > 0 && c->get_length() >= bundle_min_len)) {
                        return false;
                    }
                }
                return true;
            };
            for (auto* c : grouping.children()) {
                if (!c->get_flag(Flags::main_cluster)) continue;
                const double t0c = c->get_cluster_t0();
                if (t0c < m_beam_window_low || t0c >= m_beam_window_high) continue;
                if (!survives(c)) continue;
                if (!legacy_main || c->get_length() > legacy_main->get_length()) legacy_main = c;
            }
            if (!legacy_main && m_nu_fallback_demoted_mains) {
                for (auto* c : grouping.children()) {
                    if (c->get_flag(Flags::main_cluster)) continue;
                    if (!c->get_flag(Flags::demoted_main)) continue;
                    const double t0c = c->get_cluster_t0();
                    if (t0c < m_beam_window_low || t0c >= m_beam_window_high) continue;
                    if (!survives(c)) continue;
                    if (!legacy_main || c->get_length() > legacy_main->get_length()) legacy_main = c;
                }
            }
            if (legacy_main && m_nu_per_bundle_min_length > 0) {
                SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: [nu_per_bundle] legacy winner is cluster {} (L {:.1f} cm); it is exempt from the {:.1f} cm floor",
                                   legacy_main->get_cluster_id(),
                                   legacy_main->get_length()/units::cm,
                                   m_nu_per_bundle_min_length);
            }
        }

        std::map<int, std::vector<Cluster*>> gid_mains;    // in-window mains
        std::map<int, std::vector<Cluster*>> gid_demoted;  // in-window demoted mains
        for (auto* cluster : grouping.children()) {
            const bool is_main = cluster->get_flag(Flags::main_cluster);
            const bool is_dem  = !is_main && cluster->get_flag(Flags::demoted_main);
            if (!is_main && !is_dem) continue;
            if (is_main) n_main_clusters++;
            const double t0c = cluster->get_cluster_t0();
            if (t0c < m_beam_window_low || t0c >= m_beam_window_high) continue;
            if (is_main) n_in_beam_clusters++;
            const int gid = cluster->get_scalar<int>("matched_flash_gid", -1);
            if (gid < 0) continue;
            if (is_main) gid_mains[gid].push_back(cluster);
            else if (m_nu_per_bundle_demoted_acts) gid_demoted[gid].push_back(cluster);
        }
        std::set<int> gids;
        for (const auto& kv : gid_mains) gids.insert(kv.first);
        for (const auto& kv : gid_demoted) gids.insert(kv.first);

        for (int gid : gids) {
            NuCandidate cand;
            cand.gid = gid;

            // Every main activity the taggers evaluated in this bundle.  The
            // taggers set a flag only on a POSITIVE verdict, so a missing
            // TGM/STM/FC flag cannot distinguish "evaluated and exonerated"
            // from "never looked at"; membership of these two lists is the
            // reproduction of their admission gate, hence evaluated = 1.
            auto add_acts = [&](const std::vector<Cluster*>& v, int is_demoted) {
                for (auto* c : v) {
                    ActivityRec a;
                    a.cluster_id = c->get_cluster_id();
                    a.length_cm  = c->get_length() / units::cm;
                    a.is_demoted = is_demoted;
                    a.tgm = c->get_flag(Flags::TGM) ? 1 : 0;
                    a.stm = c->get_flag(Flags::STM) ? 1 : 0;
                    a.fc  = c->get_flag(Flags::FC)  ? 1 : 0;
                    a.lm  = c->get_scalar<int>("lm_flag", -1);
                    a.evaluated = 1;
                    cand.acts.push_back(a);
                }
            };
            auto mit = gid_mains.find(gid);
            auto dit = gid_demoted.find(gid);
            if (mit != gid_mains.end()) add_acts(mit->second, 0);
            if (dit != gid_demoted.end()) add_acts(dit->second, 1);

            const double per_bundle_min_len = m_nu_per_bundle_min_length * units::cm;

            auto pick = [&](const std::vector<Cluster*>& v) -> Cluster* {
                Cluster* best = nullptr;
                for (auto* c : v) {
                    if (m_nu_skip_cosmic) {
                        const bool tgm = c->get_flag(Flags::TGM);
                        const bool stm = c->get_flag(Flags::STM);
                        const int  lm  = c->get_scalar<int>("lm_flag", -1);
                        if (tgm || stm || lm > 0) {
                            SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: [nu_per_bundle] gid {} activity {} (L {:.1f} cm) cosmic-tagged (TGM={} STM={} lm_flag={}); not a candidate",
                                               gid, c->get_cluster_id(), c->get_length()/units::cm, tgm, stm, lm);
                            continue;
                        }
                    }
                    // The dot guard.  Without it, dropping the bundle veto
                    // promotes sub-cm blobs to "the neutrino" of their bundle
                    // -- they fit as a muon at rest (~107 MeV) and place the
                    // vertex on a dot (mcp1k evt 391854: 1.7 cm -> 108.9 MeV).
                    // The legacy winner is exempt, so this can only ever remove
                    // rows per-bundle mode ADDS, never a row the legacy chain
                    // reports; the two cannot be told apart by length alone
                    // (see the legacy_main derivation above).
                    if (c != legacy_main && per_bundle_min_len > 0
                        && c->get_length() < per_bundle_min_len) {
                        SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: [nu_per_bundle] gid {} activity {} (L {:.1f} cm) is under the {:.1f} cm floor and is not the legacy winner; not a candidate (nu_per_bundle_min_length)",
                                           gid, c->get_cluster_id(), c->get_length()/units::cm,
                                           m_nu_per_bundle_min_length);
                        continue;
                    }
                    if (!best || c->get_length() > best->get_length()) best = c;
                }
                return best;
            };
            if (mit != gid_mains.end()) cand.main = pick(mit->second);
            if (!cand.main && m_nu_fallback_demoted_mains && dit != gid_demoted.end()) {
                cand.main = pick(dit->second);
            }
            if (!cand.main) {
                SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: [nu_per_bundle] gid {}: no neutrino candidate among {} evaluated activit(ies)",
                                   gid, cand.acts.size());
                continue;
            }
            for (auto& a : cand.acts) {
                if (a.cluster_id == cand.main->get_cluster_id()) a.is_selected = 1;
            }

            // Companions: the same rule and the same knob as the legacy path.
            for (auto* cluster : grouping.children()) {
                if (cluster == cand.main) continue;
                if (!cluster->get_flag(Flags::associated_cluster)) continue;
                if (cluster->get_scalar<int>("matched_flash_gid", -1) != gid) continue;
                if (m_skip_cosmic_companions
                    && (cluster->get_flag(Flags::TGM) || cluster->get_flag(Flags::STM))
                    && cluster->get_length() >= m_cosmic_companion_min_length * units::cm) {
                    SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: [nu_per_bundle] gid {} companion cluster {} (L {:.1f} cm, TGM={} STM={}) dropped (skip_cosmic_companions, floor {:.1f} cm)",
                                       gid, cluster->get_cluster_id(), cluster->get_length()/units::cm,
                                       cluster->get_flag(Flags::TGM), cluster->get_flag(Flags::STM),
                                       m_cosmic_companion_min_length);
                    continue;
                }
                cand.others.push_back(cluster);
            }
            SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: [nu_per_bundle] gid {}: candidate main cluster {} (t0 {:.3f} us, L {:.1f} cm, {} associated) of {} evaluated activit(ies)",
                               gid, cand.main->get_cluster_id(), cand.main->get_cluster_t0()/units::us,
                               cand.main->get_length()/units::cm, cand.others.size(), cand.acts.size());
            candidates.push_back(std::move(cand));
        }

        // Order by the selected activity's length, longest first (ties by gid,
        // so this is deterministic).  This is NOT cosmetic: candidate 0 keeps
        // the unnamed TrackFitting slot, and every consumer not yet converted
        // to walk the "nu<i>" slots -- the three Bee PR point layers,
        // SbndPrMagnifyTrackingVisitor's T_rec_charge/T_proj_data,
        // PrDisplayDump -- still reads exactly that slot.  Enumerating by gid
        // instead put a 1.7 cm shard from the low-gid drift side in slot 0 on
        // nueCC48 evt 10550/234638/267597, which silently emptied those layers
        // (the shard reconstructs no vertex).  Longest-first is also the legacy
        // selector's own rule, so slot 0 now holds the same candidate the
        // pre-pr/94 chain would have chosen, and nu_index 0 means "primary".
        std::stable_sort(candidates.begin(), candidates.end(),
                         [](const NuCandidate& a, const NuCandidate& b) {
                             auto sel_len = [](const NuCandidate& c) {
                                 for (const auto& x : c.acts)
                                     if (x.is_selected) return x.length_cm;
                                 return 0.0f;
                             };
                             const float la = sel_len(a), lb = sel_len(b);
                             if (la != lb) return la > lb;
                             return a.gid < b.gid;
                         });
    }

    // doc pr/94 -- unify the two shapes.  The two legacy branches above select
    // at most one event-wide winner into main_cluster/other_clusters (this is
    // also the path taken with no beam gate, whatever nu_per_bundle says);
    // per-bundle mode filled `candidates` directly and left main_cluster null.
    if (candidates.empty() && main_cluster) {
        NuCandidate cand;
        cand.main   = main_cluster;
        cand.others = other_clusters;
        cand.gid    = main_cluster->get_scalar<int>("matched_flash_gid", -1);
        candidates.push_back(std::move(cand));
    }

    if (candidates.empty()) {
        SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino: no main cluster selected ({} mains, {} in-window); skipping.",
                            n_main_clusters, n_in_beam_clusters);
        return;
    }

    // ---- doc pr/94: one full PR pass per candidate ------------------- //
    // Legacy mode runs exactly ONE iteration, on the member fitter, with
    // main_cluster / other_clusters holding exactly what the untouched
    // selection above put there -- so the body below is byte-identical to the
    // pre-pr/94 single-candidate chain.  acc_segment_id is hoisted out of the
    // loop on purpose: it restarts at 0 per candidate otherwise, and shower
    // ids would then collide between bundles.
    int acc_segment_id = 0;

    // doc pr/94 round 3 -- nu_selected_as_main.  Scoped, exception-safe
    // set/restore of Flags::main_cluster on the selected candidate.  The PR
    // chain reads main-ness from that flag in three different files, so
    // threading a parameter into find_proto_vertex would only reach one of
    // them; setting the flag for exactly the candidate's own PR pass is what
    // "give it the same treatment as a main" means.  This guard's OWN
    // restore only touches the candidate's own pointer -- see
    // MainFlagSnapshotAllGuard below (sbnd_xin/docs/75) for the case where
    // that is not enough.  Disarmed (arm=false, the default knob state) or
    // already-a-main => no write at all => byte-identical.
    struct SelectedMainFlagGuard {
        Cluster* cluster{nullptr};
        bool armed{false};
        SelectedMainFlagGuard(Cluster* c, bool arm)
        {
            if (!arm || !c) return;
            if (c->get_flag(Flags::main_cluster)) return;   // already a main
            cluster = c;
            armed = true;
            cluster->set_flag(Flags::main_cluster, 1);
        }
        ~SelectedMainFlagGuard()
        {
            if (armed && cluster) cluster->set_flag(Flags::main_cluster, 0);
        }
        SelectedMainFlagGuard(const SelectedMainFlagGuard&) = delete;
        SelectedMainFlagGuard& operator=(const SelectedMainFlagGuard&) = delete;
    };

    // doc 75 -- nu_selected_as_main_snapshot_all.  Closes the gap
    // SelectedMainFlagGuard leaves open: the DL/SCN vertex path
    // (determine_overall_main_vertex_DL -> swap_main_cluster,
    // NeutrinoVertexFinder.cxx:4976, called on the REAL main_cluster/
    // other_clusters pointers at :2141-2148 below) can move
    // Flags::main_cluster onto a DIFFERENT cluster within the candidate's
    // own bundle mid-pass.  Snapshotting BEFORE SelectedMainFlagGuard runs
    // (construction order below) captures the pristine pre-pass state of
    // every cluster in {main_cluster} u other_clusters; restoring on
    // destruction -- which runs AFTER SelectedMainFlagGuard's own (reverse
    // construction order) -- is therefore authoritative regardless of how
    // many swaps happened in between, including SelectedMainFlagGuard's own
    // write.  Independent of nu_per_bundle (this swap call site is not
    // pr/94-specific); a no-op unless armed.
    struct MainFlagSnapshotAllGuard {
        std::vector<std::pair<Cluster*, bool>> snapshot;
        WireCell::Log::logptr_t glog;
        bool armed{false};
        MainFlagSnapshotAllGuard(Cluster* main, const std::vector<Cluster*>& others, bool arm,
                                  WireCell::Log::logptr_t lg)
            : glog(lg)
        {
            if (!arm) return;
            armed = true;
            if (main) snapshot.emplace_back(main, main->get_flag(Flags::main_cluster));
            for (auto* c : others) {
                if (c) snapshot.emplace_back(c, c->get_flag(Flags::main_cluster));
            }
        }
        ~MainFlagSnapshotAllGuard()
        {
            if (!armed) return;
            for (auto& cv : snapshot) {
                if (!cv.first) continue;
                const bool live = cv.first->get_flag(Flags::main_cluster);
                if (live != cv.second && glog) {
                    // A swap moved Flags::main_cluster onto/off this cluster
                    // during the pass -- the census signal doc 75 asked for
                    // (how often the narrow SelectedMainFlagGuard restore
                    // alone would have left this cluster's flag wrong).
                    SPDLOG_LOGGER_INFO(glog, "TaggerCheckNeutrino: [nu_selected_as_main_snapshot_all] "
                                       "restoring cluster {} main_cluster flag {} -> {} (DL-swap leak closed)",
                                       cv.first->get_cluster_id(), live, cv.second);
                }
                cv.first->set_flag(Flags::main_cluster, cv.second ? 1 : 0);
            }
        }
        MainFlagSnapshotAllGuard(const MainFlagSnapshotAllGuard&) = delete;
        MainFlagSnapshotAllGuard& operator=(const MainFlagSnapshotAllGuard&) = delete;
    };

    for (size_t nu_index = 0; nu_index < candidates.size(); ++nu_index) {
        main_cluster   = candidates[nu_index].main;
        other_clusters = candidates[nu_index].others;

        // Construct BEFORE SelectedMainFlagGuard so its snapshot is taken
        // pre-write (see the class comment above) -- destroyed AFTER it.
        MainFlagSnapshotAllGuard main_flag_snapshot(main_cluster, other_clusters,
                                                     m_nu_selected_as_main_snapshot_all, log);

        // doc pr/94 round 3: the selected candidate is treated as the main
        // cluster for the duration of this pass (knob nu_selected_as_main).
        SelectedMainFlagGuard selected_main_guard(main_cluster, m_nu_selected_as_main);
        if (selected_main_guard.armed) {
            SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: [nu_selected_as_main] candidate cluster {} "
                               "(L {:.1f} cm) is a demoted main; flagged main_cluster for its own PR pass",
                               main_cluster->get_cluster_id(), main_cluster->get_length()/units::cm);
        }

        // The first candidate reuses the configured member fitter (legacy,
        // byte-identical).  Later candidates need their OWN fitter: add_graph()
        // replaces the graph but sync_from_graph() ACCUMULATES clusters and
        // blobs, so reusing one fitter would leak bundle i-1's charge data into
        // bundle i.  Seeded from the member fitter's own parameters, which is
        // where load_trackfitting_config()'s set_parameter() calls landed -- a
        // freshly-preset fitter would silently drop the whole config file.
        auto track_fitter = m_track_fitter;
        if (nu_index > 0) {
            track_fitter = std::make_shared<TrackFitting>(TrackFittingPresets::create_with_current_values());
            track_fitter->set_parameters(m_track_fitter->get_parameters());
            track_fitter->set_perf(m_track_fitter->get_perf());
            track_fitter->set_detector_volume(m_dv);
            track_fitter->set_pc_transforms(m_pcts);
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
                    dump_path, *main_cluster, main_cluster, true, *track_fitter);
            }
        }

        SPDLOG_LOGGER_TRACE(log, "Number of Main Clusters: {}", n_main_clusters);

        IndexedVertexSet vertices_in_long_muon;
        IndexedSegmentSet segments_in_long_muon;
        // doc 84 round 3: cathode-bridge absorbed members, captured by
        // long_muon_cathode_bridge_pass and consumed only by the MCS driver
        // under mcs_bridged_members.  Deliberately separate from
        // segments_in_long_muon, which taggers keep reading unchanged.
        IndexedSegmentSet bridged_in_long_muon;
        VertexPtr main_vertex = nullptr;
        ClusterVertexMap map_cluster_main_vertices;

        // Pre-load charge data for all beam-flash clusters once so that
        // do_multi_tracking calls throughout pattern recognition can use
        // flag_force_load_data=false and avoid redundant prepare_data() calls.
        {
            std::vector<WireCell::Clus::Facade::Cluster*> clusters_to_preload;
            clusters_to_preload.push_back(main_cluster);
            for (auto* c : other_clusters) clusters_to_preload.push_back(c);
            track_fitter->preload_clusters(clusters_to_preload);
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
                    /*flag_back_search=*/true, *track_fitter);
            }
        }

        // Create PRGraph and first segment
        auto pr_graph = std::make_shared<WireCell::Clus::PR::Graph>();
        track_fitter->add_graph(pr_graph);

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
        // doc sbnd_xin/docs/pr/104 -- main-vertex junction snap.
        pattern_algos.m_vertex_junction_snap = m_vertex_junction_snap;
        pattern_algos.m_vjs_radius       = m_vjs_radius * units::cm;     // cm -> internal
        pattern_algos.m_vjs_min_arm      = m_vjs_min_arm * units::cm;    // cm -> internal
        pattern_algos.m_vjs_min_prongs   = m_vjs_min_prongs;
        pattern_algos.m_vjs_collinear    = m_vjs_collinear;              // deg, no conversion
        pattern_algos.m_vjs_fit_margin   = m_vjs_fit_margin * units::cm; // cm -> internal
        pattern_algos.m_vjs_fit_rms      = m_vjs_fit_rms * units::cm;    // cm -> internal
        pattern_algos.m_vjs_override_kink_snap = m_vjs_override_kink_snap;
        pattern_algos.m_vjs_min_move     = m_vjs_min_move * units::cm;   // cm -> internal
        // sbnd_xin/docs/73 sec 12 (round 3).
        pattern_algos.m_esva_ignore_empty_2d = m_esva_ignore_empty_2d;
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
        pattern_algos.m_mvga_op1_radius   = m_mvga_op1_radius * units::cm; // cm -> internal; 0 and the -1 sentinel both survive the scale (doc pr/83 r3)
        pattern_algos.m_mvga_op1_dup_frac = m_mvga_op1_dup_frac;           // fraction, no conversion (doc pr/83 r3)
        pattern_algos.m_mvga_op1_post     = m_mvga_op1_post;               // doc pr/83 r3 (class A)
        pattern_algos.m_swap_orphan_dup_audit = m_swap_orphan_dup_audit;   // doc pr/83 r3 (Mechanism C)
        pattern_algos.m_mvga_proj_dup_frac  = m_mvga_proj_dup_frac;        // fraction, no conversion (doc pr/83 r4)
        pattern_algos.m_mvga_proj_dqdx_ratio = m_mvga_proj_dqdx_ratio;     // ratio, no conversion (doc pr/83 r4)
        pattern_algos.m_mvga_proj_angle = m_mvga_proj_angle;               // deg, no conversion (doc pr/83 r4b)
        pattern_algos.m_mvga_ac_veto_radius = m_mvga_ac_veto_radius * units::cm; // cm -> internal (doc pr/99 round 2)
        pattern_algos.m_mvga_ac_chord_max   = m_mvga_ac_chord_max * units::cm;   // cm -> internal (doc pr/99 round 2)
        pattern_algos.m_mvga_ac_no_cascade  = m_mvga_ac_no_cascade;               // doc pr/99 round 2
        pattern_algos.m_mvga_passthru       = m_mvga_passthru * units::cm;        // cm -> internal (doc pr/103)
        pattern_algos.m_mvga_passthru_tol   = m_mvga_passthru_tol * units::cm;    // cm -> internal (doc pr/103)
        pattern_algos.m_mvga_interposed_fallback = m_mvga_interposed_fallback;    // doc pr/103
        pattern_algos.m_mvga_interposed_fallback_min_angle = m_mvga_interposed_fallback_min_angle;  // deg (doc pr/103)
        pattern_algos.m_mvga_dup_starved_asym = m_mvga_dup_starved_asym;          // ratio, no conversion (doc pr/99 round 2)
        pattern_algos.m_mvga_dup_starved_mip = m_mvga_dup_starved_mip;            // ratio, no conversion (doc pr/99 round 2)
        pattern_algos.m_mvga_dup_starved_span = m_mvga_dup_starved_span;          // ratio, no conversion (doc pr/99 round 2)
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
        pattern_algos.m_dl_vtx_cloud_no_exclusion = m_dl_vtx_cloud_no_exclusion;   // doc pr/106 sec 10
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
        // doc sbnd_xin/docs/pr/102 P1+P2.
        pattern_algos.m_other_seg_keep_isolated_min_nnf    = m_other_seg_keep_isolated_min_nnf;
        pattern_algos.m_other_seg_keep_isolated_len_admit  = m_other_seg_keep_isolated_len_admit * units::cm; // cm -> internal
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
        // sbnd_xin/docs/74 G1/G2: consistent-FV routing for cosmic_tagger().
        // m_fv_tolerance is already INTERNAL units (read raw in configure(),
        // same values cluster_fc_check consumes) -- no conversion.
        pattern_algos.m_cosmic_fiducial =
            (m_cosmic_consistent_fv && m_use_fiducial) ? m_fiducial : nullptr;
        pattern_algos.m_cosmic_fv_tolerance = m_fv_tolerance;
        // sbnd_xin/docs/75: consistent-FV routing for the nue/single-photon
        // taggers -- same configured fiducial, own hardcoded per-site
        // tolerances (see NeutrinoTaggerNuE.cxx / NeutrinoTaggerSinglePhoton.cxx).
        pattern_algos.m_nue_fiducial =
            (m_nue_sp_consistent_fv && m_use_fiducial) ? m_fiducial : nullptr;
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
        pattern_algos.m_kine_charge.dedup               = m_kine_charge_dedup;                          // doc pr/99 r3 C1
        pattern_algos.m_kine_charge.rebuild             = m_kine_charge_rebuild;                        // doc pr/99 r3 C1b
        pattern_algos.m_kine_charge.track_ctx           = m_kine_charge_track_ctx;                      // doc pr/101 K1
        pattern_algos.m_kine_charge.mass_rules          = m_kine_mass_rules;                            // doc pr/101 K2
        pattern_algos.m_kine_charge.hadronic_dqdx       = m_kine_hadronic_dqdx;                         // doc pr/101 K3
        pattern_algos.m_kine_charge.long_muon_mode      = m_kine_long_muon_mode;                        // doc pr/101 K4
        pattern_algos.m_kine_charge.long_muon_ratio_lo  = m_kine_long_muon_ratio_lo;
        pattern_algos.m_kine_charge.long_muon_ratio_hi  = m_kine_long_muon_ratio_hi;
        pattern_algos.m_kine_charge.long_muon_range_fallback = m_long_muon_range_empty_chain_fallback;  // doc 84 round 1 (P1)
        pattern_algos.m_kine_charge.long_muon_members_geometry = m_long_muon_members_geometry;          // doc 84 round 2
        pattern_algos.m_kine_charge.mainvtx_used_guard  = m_kine_mainvtx_used_guard;                    // doc pr/101 K5
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
        pattern_algos.m_shower_endpoint_skip_orphan_vtx = m_shower_endpoint_skip_orphan_vtx;
        pattern_algos.m_shower_walk_visited_parity = m_shower_walk_visited_parity;
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
        pattern_algos.m_michel_stem_muon_rescue                   = m_michel_stem_muon_rescue;                   // F14
        pattern_algos.m_shower_in_cascade_guard                   = m_shower_in_cascade_guard;                   // pr/74 P1
        pattern_algos.m_shower_in_max_len                         = m_shower_in_max_len * units::cm;             // pr/74 P1
        pattern_algos.m_shower_in_mip_hi                          = m_shower_in_mip_hi;                          // pr/74 P1
        pattern_algos.m_shower_connect_from_vertices_straight_guard  = m_shower_connect_from_vertices_straight_guard;  // pr/40 r9 (r8 Part A)
        pattern_algos.m_shower_connect_start_seg_straight_guard      = m_shower_connect_start_seg_straight_guard;      // pr/40 r9 (r7 c2c)
        pattern_algos.m_examine_direction_dirsign_shower_in_guard    = m_examine_direction_dirsign_shower_in_guard;    // pr/40 r9 (r7 c2a)
        pattern_algos.m_daughter_shower_angle_reclass_straight_guard = m_daughter_shower_angle_reclass_straight_guard; // pr/40 r9 (r7 c2b)
        pattern_algos.m_shower_topo_reexam_straight_guard            = m_shower_topo_reexam_straight_guard;            // pr/40 r9 (r7 c1)
        pattern_algos.m_sfv_kink_max                                 = m_sfv_kink_max;                                 // pr/40 r9 (degrees)
        pattern_algos.m_shower_nv_bridge_track                       = m_shower_nv_bridge_track;                       // pr/40 r9 B2
        pattern_algos.m_shower_nv_bridge_max_gap                     = m_shower_nv_bridge_max_gap * units::cm;         // pr/40 r9 B2
        pattern_algos.m_shower_nv_main_pi_init                       = m_shower_nv_main_pi_init;                       // pr/97 D1
        pattern_algos.m_kine_drop_stray_satellites                   = m_kine_drop_stray_satellites;                   // pr/92
        pattern_algos.m_kine_sat_min_energy                          = m_kine_sat_min_energy * units::MeV;             // pr/92
        pattern_algos.m_kine_sat_prox_max                            = m_kine_sat_prox_max * units::cm;                // pr/92
        pattern_algos.m_kine_sat_angle_bad                           = m_kine_sat_angle_bad;                           // pr/92 (degrees)
        pattern_algos.m_kine_sat_angle_main                          = m_kine_sat_angle_main;                          // pr/92 (degrees)
        pattern_algos.m_kine_sat_far_dis                             = m_kine_sat_far_dis * units::cm;                 // pr/92
        pattern_algos.m_kine_sat_axis_dis_cut                        = m_kine_sat_axis_dis_cut * units::cm;            // pr/92
        pattern_algos.m_kine_sat_cont_kink                           = m_kine_sat_cont_kink;                           // pr/92 (degrees)
        pattern_algos.m_kine_sat_track_max_nseg                      = static_cast<int>(m_kine_sat_track_max_nseg);   // pr/92 r2 (count)
        pattern_algos.m_kine_sat_em_far_dis                          = m_kine_sat_em_far_dis * units::cm;             // pr/92 r2
        pattern_algos.m_michel_stem_michel_check                  = m_michel_stem_michel_check;                  // pr/74 P2
        pattern_algos.m_michel_stem_max_far_len                   = m_michel_stem_max_far_len * units::cm;       // pr/74 P2
        pattern_algos.m_shower_stem_backfill                      = m_shower_stem_backfill;                      // pr/74 K4
        pattern_algos.m_stem_backfill_max_len                     = m_stem_backfill_max_len * units::cm;         // pr/74 K4
        pattern_algos.m_stem_backfill_mip_lo                      = m_stem_backfill_mip_lo;                      // pr/74 K4
        pattern_algos.m_stem_backfill_mip_hi                      = m_stem_backfill_mip_hi;                      // pr/74 K4
        pattern_algos.m_stem_backfill_min_shower_len              = m_stem_backfill_min_shower_len * units::cm;  // pr/74 K4
        pattern_algos.m_shower_conn3_unreachable                  = m_shower_conn3_unreachable;                  // pr/74 K5
        pattern_algos.m_conn3_unreachable_min_len                 = m_conn3_unreachable_min_len * units::cm;     // pr/74 K5
        pattern_algos.m_conn3_stitch_max                          = m_conn3_stitch_max * units::cm;              // pr/84 r2 F3
        pattern_algos.m_shower_dedup_start_seg                    = m_shower_dedup_start_seg;                    // pr/84 r3 S1
        pattern_algos.m_shower_traj_michel_stem                   = m_shower_traj_michel_stem;                   // pr/74 K6
        pattern_algos.m_michel_stem_traj_min_len                  = m_michel_stem_traj_min_len * units::cm;      // pr/74 K6
        pattern_algos.m_michel_stem_traj_max_len                  = m_michel_stem_traj_max_len * units::cm;      // pr/74 K6
        pattern_algos.m_michel_stem_traj_mip_lo                   = m_michel_stem_traj_mip_lo;                   // pr/74 K6 (dimensionless ratio)
        pattern_algos.m_michel_stem_traj_max_far_len              = m_michel_stem_traj_max_far_len * units::cm;  // pr/74 K6
        pattern_algos.m_michel_stem_traj_min_kink_deg             = m_michel_stem_traj_min_kink_deg;             // pr/74 K6 (degrees)
        pattern_algos.m_shower_long_muon_keep_type                = m_shower_long_muon_keep_type;                // doc pr/44
        pattern_algos.m_shower_bragg_protect_start_segment        = m_shower_bragg_protect_start_segment;        // doc pr/40 round 10
        pattern_algos.m_shower_reclass_case_b_dqdx_guard          = m_shower_reclass_case_b_dqdx_guard;          // doc pr/93 Cause A
        pattern_algos.m_shower_accept_pid_guard                   = m_shower_accept_pid_guard;                   // doc pr/93 Cause B
        pattern_algos.m_shower_pid_guard_min_len                  = m_shower_pid_guard_min_len * units::cm;      // doc pr/93 shared floor
        pattern_algos.m_shower_vote_track_pid_counts              = m_shower_vote_track_pid_counts;              // doc pr/93 Cause C
        pattern_algos.m_shower_cone_absorb_guard               = m_shower_cone_absorb_guard;               // doc pr/93 Cause D
        pattern_algos.m_shower_detach_track_stem                  = m_shower_detach_track_stem;                  // doc pr/93 r4
        pattern_algos.m_shower_ghost_member_drop                  = m_shower_ghost_member_drop;                  // doc pr/99 r2
        pattern_algos.m_shower_ghost_overlap_frac                 = m_shower_ghost_overlap_frac;                 // fraction, no conversion (doc pr/99 r2)
        pattern_algos.m_shower_ghost_dqdx_ratio                   = m_shower_ghost_dqdx_ratio;                   // ratio, no conversion (doc pr/99 r2)
        pattern_algos.m_shower_ghost_min_len                      = m_shower_ghost_min_len * units::cm;          // cm -> internal (doc pr/99 r2)
        pattern_algos.m_shower_hadronic_tag                       = m_shower_hadronic_tag;                       // doc pr/99 r3 A5
        pattern_algos.m_shower_hadronic_min_len                   = m_shower_hadronic_min_len * units::cm;       // cm -> internal (doc pr/99 r3)
        pattern_algos.m_shower_hadronic_scan_len                  = m_shower_hadronic_scan_len * units::cm;      // cm -> internal (doc pr/99 r3)
        pattern_algos.m_shower_hadronic_bin                       = m_shower_hadronic_bin * units::cm;           // cm -> internal (doc pr/99 r3)
        pattern_algos.m_shower_hadronic_r_cyl                     = m_shower_hadronic_r_cyl * units::cm;         // cm -> internal (doc pr/99 r3)
        pattern_algos.m_shower_hadronic_r_core                    = m_shower_hadronic_r_core * units::cm;        // cm -> internal (doc pr/99 r3)
        pattern_algos.m_shower_hadronic_growth_max                = m_shower_hadronic_growth_max;                // ratio, no conversion (doc pr/99 r3)
        pattern_algos.m_shower_hadronic_growth_bragg              = m_shower_hadronic_growth_bragg;              // ratio, no conversion (doc pr/99 r3)
        pattern_algos.m_shower_hadronic_bragg_ratio               = m_shower_hadronic_bragg_ratio;               // ratio, no conversion (doc pr/99 r3)
        pattern_algos.m_shower_hadronic_stem_ratio                = m_shower_hadronic_stem_ratio;                // MIP units, no conversion (doc pr/99 r3)
        pattern_algos.m_kine_count_orphan_tracks                  = m_kine_count_orphan_tracks;                  // doc pr/93 r4
        pattern_algos.m_kine_orphan_track_min                     = m_kine_orphan_track_min * units::cm;         // doc pr/93 r4
        pattern_algos.m_shower_pass4_best_owner                   = m_shower_pass4_best_owner;                   // doc pr/117 r1
        pattern_algos.m_shower_merge_relax                        = m_shower_merge_relax;                        // doc pr/117 r1
        pattern_algos.m_shower_merge_relax_dis                    = m_shower_merge_relax_dis * units::cm;        // cm -> internal (doc pr/117 r1)
        pattern_algos.m_shower_merge_relax_angle                  = m_shower_merge_relax_angle;                  // deg, no conversion (doc pr/117 r1)
        pattern_algos.m_shower_merge_relax_min_len                = m_shower_merge_relax_min_len * units::cm;    // cm -> internal (doc pr/117 r1)
        pattern_algos.m_shower_flank_absorb                       = m_shower_flank_absorb;                       // doc pr/117 r1
        pattern_algos.m_shower_flank_absorb_max_dis               = m_shower_flank_absorb_max_dis * units::cm;   // cm -> internal (doc pr/117 r1)
        pattern_algos.m_shower_flank_absorb_max_len               = m_shower_flank_absorb_max_len * units::cm;   // cm -> internal (doc pr/117 r1)
        pattern_algos.m_shower_ex1_conn3_body_dis                 = m_shower_ex1_conn3_body_dis;                 // doc pr/118 r1
        pattern_algos.m_shower_merge_relax_continuity             = m_shower_merge_relax_continuity;             // doc pr/118 r1
        pattern_algos.m_shower_merge_relax_cont_frac              = m_shower_merge_relax_cont_frac;              // fraction, no conversion (doc pr/118 r1)
        pattern_algos.m_shower_merge_relax_cont_gap               = m_shower_merge_relax_cont_gap * units::cm;   // cm -> internal (doc pr/118 r1)
        pattern_algos.m_shower_merge_relax_cont_qmed              = m_shower_merge_relax_cont_qmed;              // charge units, no conversion (doc pr/118 r1)
        pattern_algos.m_shower_merge_relax_cont_axis              = m_shower_merge_relax_cont_axis;              // deg, no conversion (doc pr/118 r1)
        pattern_algos.m_shower_merge_relax_cont_dmax              = m_shower_merge_relax_cont_dmax * units::cm;  // cm -> internal (doc pr/118 r1)
        pattern_algos.m_shower_merge_relax_cont_t1_gap            = m_shower_merge_relax_cont_t1_gap * units::cm; // cm -> internal (doc pr/118 r1)
        pattern_algos.m_shower_merge_relax_cont_t1_fold           = m_shower_merge_relax_cont_t1_fold;           // deg, no conversion (doc pr/118 r1)
        pattern_algos.m_stem_backfill_back_guard                  = m_stem_backfill_back_guard;                  // doc pr/120 r1
        pattern_algos.m_stem_backfill_back_ang                    = m_stem_backfill_back_ang;                    // deg, no conversion (doc pr/120 r1)
        pattern_algos.m_shower_ex1_walk_em_track_guard            = m_shower_ex1_walk_em_track_guard;            // doc pr/120 r1
        pattern_algos.m_shower_ex1_walk_em_track_len              = m_shower_ex1_walk_em_track_len * units::cm;  // cm -> internal (doc pr/120 r1)
        pattern_algos.m_shower_ex1_dedup_rehome                   = m_shower_ex1_dedup_rehome;                   // doc pr/121 r1
        pattern_algos.m_shower_pass4_prune_detached               = m_shower_pass4_prune_detached;               // doc pr/123 r1
        pattern_algos.m_shower_pass4_prune_gap                    = m_shower_pass4_prune_gap * units::cm;        // doc pr/123 r1, cm -> internal
        pattern_algos.m_shower_pass4_prune_gap2                   = m_shower_pass4_prune_gap2 * units::cm;       // doc pr/124 A, cm -> internal
        pattern_algos.m_shower_pass4_prune2_ang                   = m_shower_pass4_prune2_ang;                   // doc pr/124 A, deg
        pattern_algos.m_shower_pass4_prune2_mdqdx                 = m_shower_pass4_prune2_mdqdx;                 // doc pr/124 A, x MIP
        pattern_algos.m_shower_pass3_cone_guard_len               = m_shower_pass3_cone_guard_len * units::cm;   // doc pr/124 C, cm -> internal
        pattern_algos.m_shower_samevtx_track_absorb               = m_shower_samevtx_track_absorb;               // doc pr/125
        pattern_algos.m_shower_samevtx_absorb_gap                 = m_shower_samevtx_absorb_gap * units::cm;     // doc pr/125, cm -> internal
        pattern_algos.m_shower_samevtx_absorb_max_len             = m_shower_samevtx_absorb_max_len * units::cm; // doc pr/125, cm -> internal
        pattern_algos.m_shower_samevtx_absorb_min_len             = m_shower_samevtx_absorb_min_len * units::cm; // doc pr/125, cm -> internal
        pattern_algos.m_shower_satellite_absorb                   = m_shower_satellite_absorb;                   // doc pr/125
        pattern_algos.m_shower_split = m_shower_split;   // doc pr/138 B2; false = no pass
        pattern_algos.m_shower_split_max_valley = m_shower_split_max_valley;   // doc pr/138 B2; sec A5.4 knee
        pattern_algos.m_shower_split_min_frac = m_shower_split_min_frac;   // doc pr/138 B2; per-seed charge share floor
        pattern_algos.m_shower_split_max_parts = m_shower_split_max_parts;   // doc pr/138 B3; 2 = the measured-exact kernel
        pattern_algos.m_shower_split_min_charge = m_shower_split_min_charge;   // doc pr/138 B1; candidate charge floor (raw Fit::dQ)
        pattern_algos.m_shower_split_min_nseg = m_shower_split_min_nseg;   // doc pr/138 B1; candidate member-count floor
        pattern_algos.m_shower_split_bundle_gap = m_shower_split_bundle_gap * units::cm;   // doc pr/138 B3, cm -> internal
        pattern_algos.m_shower_split_snap = m_shower_split_snap;   // doc pr/138 B3; k>=3 bundle dominance floor
        pattern_algos.m_shower_split_skip_shared = m_shower_split_skip_shared;   // doc pr/139 P1.1; refuse a component holding a segment another shower also owns
        pattern_algos.m_shower_split_max_impact = m_shower_split_max_impact * units::cm;   // doc pr/139 P1.2; cm -> internal; 0 = no bound
        pattern_algos.m_shower_split_em_start = m_shower_split_em_start;   // doc pr/139 P1.3; seed the daughter on its nearest EM-typed member
        pattern_algos.m_shower_split_rehome = m_shower_split_rehome;   // doc pr/139 P1.4; offer an orphan daughter to the nearest larger EM shower
        pattern_algos.m_shower_split_rehome_gap = m_shower_split_rehome_gap * units::cm;   // doc pr/139 P1.4; cm; max daughter->host 3-D gap, cm -> internal
        pattern_algos.m_shower_satellite_absorb_max_mev           = m_shower_satellite_absorb_max_mev * units::MeV;  // doc pr/125, MeV -> internal
        pattern_algos.m_shower_satellite_absorb_host_mev          = m_shower_satellite_absorb_host_mev * units::MeV; // doc pr/125, MeV -> internal
        pattern_algos.m_shower_pass4_track_guard_len              = m_shower_pass4_track_guard_len * units::cm;  // doc pr/123 r1, cm -> internal
        pattern_algos.m_shower_pass4_prox_guard_len               = m_shower_pass4_prox_guard_len * units::cm;   // doc pr/130 item 1b, cm -> internal
        pattern_algos.m_shower_pass3_backfill_guard_len           = m_shower_pass3_backfill_guard_len * units::cm; // doc pr/130 item 1b, cm -> internal
        pattern_algos.m_stem_backfill_back_dvtx                   = m_stem_backfill_back_dvtx * units::cm;         // doc pr/130 item B, cm -> internal
        pattern_algos.m_shower_pass4_prefilter_v1_escape          = m_shower_pass4_prefilter_v1_escape;            // doc pr/136 r2 (bool, no scaling)
        pattern_algos.m_shower_pass4_prefilter_v1_max_v2          = m_shower_pass4_prefilter_v1_max_v2;            // doc pr/136 r2, deg (no scaling)
        pattern_algos.m_shower_pass4_prefilter_v1_max_dis         = m_shower_pass4_prefilter_v1_max_dis * units::cm; // doc pr/136 r3, cm -> internal
        pattern_algos.m_pi0_mass_offset                           = m_pi0_mass_offset * units::MeV;                // doc pr/132 K1, MeV -> internal
        pattern_algos.m_pi0_assoc_angle_deg                       = m_pi0_assoc_angle_deg;                         // doc pr/132 K2, deg (no scaling)
        pattern_algos.m_pi0_attached_partner_min                  = m_pi0_attached_partner_min_mev * units::MeV;   // doc pr/132 K3, MeV -> internal
        pattern_algos.m_pi0_nv_allow_type2                        = m_pi0_nv_allow_type2;                          // doc pr/132 K4
        pattern_algos.m_pi0_nv_max_prongs                         = m_pi0_nv_max_prongs;                           // doc pr/132 K5
        pattern_algos.m_pi0_readmit_retyped                       = m_pi0_readmit_retyped;                         // doc pr/132 K7
        pattern_algos.m_pi0_admit_type3                           = m_pi0_admit_type3;                             // doc pr/132 K8
        pattern_algos.m_pi0_crumb_assoc_max                       = m_pi0_crumb_assoc_mev * units::MeV;            // doc pr/132 K9, MeV -> internal
        pattern_algos.m_pi0_collinear_merge_deg                   = m_pi0_collinear_merge_deg;                     // doc pr/132 K12, deg (no scaling)
        pattern_algos.m_pi0_nv_partner_min                        = m_pi0_nv_partner_min_mev * units::MeV;         // doc pr/132 K13, MeV -> internal
        pattern_algos.m_pi0_nv_retry_paired                       = m_pi0_nv_retry_paired;                         // doc pr/132 K14
        pattern_algos.m_pi0_reseat_start_assoc                    = m_pi0_reseat_start_assoc;                      // doc pr/132 K15
        pattern_algos.m_em_collinear_merge_deg                    = m_shower_em_collinear_deg;                     // doc pr/132 K16, deg (no scaling)
        pattern_algos.m_em_collinear_merge_dis                    = m_shower_em_collinear_dis_cm * units::cm;      // doc pr/132 K16, cm -> internal
        pattern_algos.m_em_collinear_merge_min_host               = m_shower_em_collinear_host_mev * units::MeV;   // doc pr/132 K16, MeV -> internal
        pattern_algos.m_em_backext_perp                           = m_shower_em_backext_perp_cm * units::cm;       // doc pr/132 K17, cm -> internal
        pattern_algos.m_em_backext_len                            = m_shower_em_backext_len_cm * units::cm;        // doc pr/132 K17, cm -> internal
        pattern_algos.m_pi0_am_dis                                = m_pi0_accept_merge_dis_cm * units::cm;         // doc pr/132 K18, cm -> internal
        pattern_algos.m_pi0_bp_miss                               = m_pi0_bp_vertex_miss_cm * units::cm;           // doc pr/132 K19, cm -> internal
        pattern_algos.m_pi0_admit_mu_showers                      = m_pi0_admit_muon_showers;                      // doc pr/133 K20
        pattern_algos.m_pi0_nc_sig_angle                          = m_pi0_nc_sig_angle_deg;                        // doc pr/133 K21, deg (unscaled, like K2)
        pattern_algos.m_pi0_nc_floor                              = m_pi0_nc_floor_mev * units::MeV;               // doc pr/133 K21 v2, MeV -> internal
        pattern_algos.m_pi0_nc_pf_assoc                           = m_pi0_nc_pf_assoc_deg;                         // doc pr/133 K21 v2.2, deg (unscaled)
        pattern_algos.m_pi0_nc_frag_merge                         = m_pi0_nc_frag_merge;                           // doc pr/134 K22
        pattern_algos.m_pi0_pf_assoc                              = m_pi0_pf_assoc_deg;                            // doc pr/134 K23, deg (unscaled)
        pattern_algos.m_pi0_prefer_main_vertex                    = m_pi0_prefer_main_vertex;                      // doc pr/134 K24
        pattern_algos.m_pi0_nv_max_vtx_shift                      = m_pi0_nv_max_vtx_shift_cm * units::cm;         // doc pr/132 K10, cm -> internal
        pattern_algos.m_pi0_nv_mass_window                        = m_pi0_nv_mass_window_mev * units::MeV;         // doc pr/132 K11, MeV -> internal
        pattern_algos.m_kine_count_guard_freed                    = m_kine_count_guard_freed;                    // doc pr/123 r2
        pattern_algos.m_kine_guard_freed_impact                   = m_kine_guard_freed_impact * units::cm;       // doc pr/129
        pattern_algos.m_kine_guard_freed_miss_deg                 = m_kine_guard_freed_miss_deg;                 // doc pr/129
        pattern_algos.m_kine_count_near_cross_cluster             = m_kine_count_near_cross_cluster;             // doc pr/128
        pattern_algos.m_kine_near_gap                             = m_kine_near_gap * units::cm;                 // doc pr/128
        pattern_algos.m_kine_near_min_len                         = m_kine_near_min_len * units::cm;             // doc pr/128
        pattern_algos.m_kine_near_end_tol                         = m_kine_near_end_tol * units::cm;             // doc pr/128
        pattern_algos.m_kine_near_kink_deg                        = m_kine_near_kink_deg;                        // doc pr/128
        pattern_algos.m_kine_count_conn4_near                     = m_kine_count_conn4_near;                     // doc pr/128
        pattern_algos.m_kine_conn4_near_gap                       = m_kine_conn4_near_gap * units::cm;           // doc pr/128
        pattern_algos.m_straight_cont_cross_cluster               = m_straight_cont_cross_cluster;               // doc pr/93 r4
        pattern_algos.m_sccc_bridge_body                          = m_sccc_bridge_body;                          // doc pr/93 r4
        pattern_algos.m_sccc_max_gap                              = m_sccc_max_gap * units::cm;                  // doc pr/93 r4
        pattern_algos.m_sccc_kink_max                             = m_sccc_kink_max;                             // deg
        pattern_algos.m_sccc_gap_aligned                          = m_sccc_gap_aligned * units::cm;              // doc pr/93 r4
        pattern_algos.m_sccc_kink_tight                           = m_sccc_kink_tight;                           // deg
        pattern_algos.m_single_muon_proton_chain_veto             = m_single_muon_proton_chain_veto;             // doc pr/43 round 2 K1
        pattern_algos.m_single_muon_long_muon_claim               = m_single_muon_long_muon_claim;               // doc pr/43 round 2 K2
        pattern_algos.m_pid_flag_reconcile                        = m_pid_flag_reconcile;                        // doc pr/43 round 2 K3
        pattern_algos.m_long_muon_stub_bridge                     = m_long_muon_stub_bridge;                     // doc pr/46
        pattern_algos.m_long_muon_stub_bridge_len_cm              = m_long_muon_stub_bridge_len;                 // doc 84 round 1 (P3)
        pattern_algos.m_long_muon_angle_relax_long                = m_long_muon_angle_relax_long;                // doc 84 round 1 (P2)
        pattern_algos.m_long_muon_angle_relax_deg                 = m_long_muon_angle_relax_deg;                 // doc 84 round 1 (P2)
        // Muon dQ/dx-vs-length envelope: c0/c1/power dimensionless, pivot cm -> internal.
        pattern_algos.m_muon_dqdx_curve = {m_muon_dqdx_curve[0], m_muon_dqdx_curve[1],
                                           m_muon_dqdx_curve[2] * units::cm, m_muon_dqdx_curve[3]};
        // Single-photon stem dE/dx conversion; the cut narrows to float so the
        // default compares bit-identically to the legacy 2.3f literal.
        pattern_algos.m_sp_dedx_use_recomb_model = m_sp_dedx_use_recomb_model;
        pattern_algos.m_sp_mean_dedx_cut         = static_cast<float>(m_sp_mean_dedx_cut);
        pattern_algos.m_recomb_model             = m_recomb_model;
        track_fitter->set_perf(m_perf);
        // doc pr/49: thread the own-blob-coverage knob into the fitter (double
        // sentinel; -1 default matches TrackFitting::Parameters, so this is a
        // no-op unless the config opts in).
        track_fitter->set_parameter("fit_blob_coverage", m_fit_blob_coverage);
        // doc pr/67 P3 (log-only end-trim probe); 0 = off = byte-identical.
        track_fitter->set_parameter("traj_cover_probe", m_traj_cover_probe ? 1.0 : 0.0);
        // doc sbnd_xin/docs/pr/107: prototype-parity point retention for the
        // dQ/dx fit; 0 = legacy drop = byte-identical.  Local fitters spawned
        // via inherit_from copy m_params, so every do_multi_tracking site is
        // covered.
        track_fitter->set_parameter("dqdx_fit_keep_all_points", m_dqdx_fit_keep_all_points ? 1.0 : 0.0);
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
            if (cov_defer_active) track_fitter->set_parameter("fit_blob_coverage", -1);
        };
        auto cov_defer_restore = [&]() {
            if (cov_defer_active) track_fitter->set_parameter("fit_blob_coverage", m_fit_blob_coverage);
        };
        if (cov_defer_active) {
            SPDLOG_LOGGER_DEBUG(log, "fit_blob_coverage_defer on: partition stage (find_proto_vertex) runs legacy fits, deweight restored for all later stages");
        }

        // acc_segment_id hoisted above the candidate loop (doc pr/94).
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

        // doc sbnd_xin/docs/pr/112 sec 11 -- the dual chain's OFF pass runs
        // FIRST (it cannot inform a decision already made), on its own graph
        // and fitter, and hands production a DualChainHint.  Knob off => the
        // hint stays empty and a nullptr reaches determine_overall_main_vertex_DL.
        PR::DualChainHint dual_hint;
        if (m_dl_vtx_dual_chain) {
            dual_hint.mode              = m_dual_chain_mode;
            dual_hint.transfer          = m_dual_chain_transfer;
            dual_hint.transfer_max      = m_dual_chain_transfer_max * units::cm;
            dual_hint.allow_cluster_swap = m_dual_chain_allow_cluster_swap;
            dual_hint.vtx_weight        = m_dual_chain_vtx_weight;
            run_dual_chain_off_pass(pattern_algos, *track_fitter, main_cluster, other_clusters, cov_defer_active, dual_hint);
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: [dual-off] pass TOTAL took {} ms", MS(Clock::now() - t0).count());
            t0 = Clock::now();
        }

        {
            // initial pattern recognitions
            // particle_data (doc pr/48): stopping templates for the two-end
            // break pass; inert unless two_end_break is on.
            cov_defer_suspend();
            pattern_algos.find_proto_vertex(*pr_graph, *main_cluster, *track_fitter, m_dv, true, 2, true, particle_data());
            cov_defer_restore();
            detg_dump("main:find_proto_vertex", *pr_graph);
            dup_stage_census("main:find_proto_vertex", *pr_graph, *main_cluster);

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
            pattern_algos.determine_main_vertex(*pr_graph, *main_cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, *track_fitter, m_dv, particle_data(), m_recomb_model);
            detg_dump("main:determine_main_vertex", *pr_graph);
            dup_stage_census("main:determine_main_vertex", *pr_graph, *main_cluster);

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
                    pattern_algos.find_proto_vertex(*pr_graph, *cluster, *track_fitter, m_dv, true, 2, false);
                    pattern_algos.clustering_points(*pr_graph, *cluster, m_dv);
                    pattern_algos.separate_track_shower(*pr_graph, *cluster);
                    pattern_algos.determine_direction(*pr_graph, *cluster, particle_data(), m_recomb_model);
                    pattern_algos.shower_determining_in_main_cluster(*pr_graph, *cluster, particle_data(), m_recomb_model, m_dv);
                    pattern_algos.determine_main_vertex(*pr_graph, *cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, *track_fitter, m_dv, particle_data(), m_recomb_model);
                    // doc sbnd_xin/docs/pr/59 round 2 (P1), other-cluster branch.
                    pattern_algos.reassociate_cluster_orphans(*pr_graph, *cluster, m_dv);
                    if (main_vertex != nullptr) {
                        map_cluster_main_vertices[cluster] = main_vertex;
                        main_vertex = nullptr;
                    }
                    detg_dump("other:long_cluster", *pr_graph);
                } else {
                    // Short cluster: no track breaking, 1 round; fall back to init_point_segment if needed
                    if (!pattern_algos.find_proto_vertex(*pr_graph, *cluster, *track_fitter, m_dv, false, 1, false)) {
                        // std::cout << "Point Cluster " << cluster->get_cluster_id() << " " << cluster->nchildren() <<std::endl;
                        pattern_algos.init_point_segment(*pr_graph, *cluster, *track_fitter, m_dv);
                    }
                    pattern_algos.clustering_points(*pr_graph, *cluster, m_dv);
                    pattern_algos.separate_track_shower(*pr_graph, *cluster);
                    pattern_algos.determine_direction(*pr_graph, *cluster, particle_data(), m_recomb_model);
                    pattern_algos.shower_determining_in_main_cluster(*pr_graph, *cluster, particle_data(), m_recomb_model, m_dv);
                    pattern_algos.determine_main_vertex(*pr_graph, *cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, *track_fitter, m_dv, particle_data(), m_recomb_model);
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

            pattern_algos.deghosting(*pr_graph, map_cluster_main_vertices, all_clusters, *track_fitter, m_dv);
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
                *track_fitter, m_dv, particle_data(), m_recomb_model,
                m_dl_weights, m_dl_vtx_cut, m_dQdx_scale, m_dQdx_offset,
                m_dl_vtx_rerank, m_dl_vtx_top_k, m_dl_vtx_min_accept_score,
                m_dl_vtx_score_scale,
                m_dl_vtx_dual_chain ? &dual_hint : nullptr);   // doc pr/112 sec 11
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
                *track_fitter, m_dv, particle_data(), m_recomb_model, true);
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
        dup_stage_census("overall_main_vertex", *pr_graph, *main_cluster);
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: overall main vertex took {} ms", MS(Clock::now() - t0).count());
        t0 = Clock::now();

    

  

        if (final_main_vertex) {
            // doc sbnd_xin/docs/pr/50: main-vertex kink-consistency snap --
            // inert unless vertex_kink_snap.  Runs after the overall main
            // vertex is final (either DL or fallback path) and BEFORE the
            // final improve_vertex, so the local optimizer polishes a
            // corner-anchored trajectory.
            if (pattern_algos.snap_main_vertex_to_kink(*pr_graph, *main_cluster, final_main_vertex,
                                                       *track_fitter, m_dv, particle_data(), m_recomb_model)) {
                map_cluster_main_vertices[main_cluster] = final_main_vertex;
                detg_dump("snap_main_vertex_to_kink", *pr_graph);
                dup_stage_census("snap_main_vertex_to_kink", *pr_graph, *main_cluster);
            }

            // doc sbnd_xin/docs/pr/104: main-vertex junction snap -- inert
            // unless vertex_junction_snap.  Re-points final_main_vertex to a
            // nearby multi-prong track junction (18255-405707/65289/66712/
            // 282072/345633: the DL re-rank lands on a stub end or a
            // pass-through point 2-4 cm from the junction the prongs meet
            // at).  Runs AFTER the kink snap and BEFORE improve_vertex, so
            // the vertex fit polishes the junction with its real prongs and
            // main_vertex_graph_audit no longer has to re-anchor prongs onto
            // the wrong point.  Pointer move only: no segment is edited.
            if (pattern_algos.snap_main_vertex_to_junction(*pr_graph, *main_cluster, final_main_vertex)) {
                map_cluster_main_vertices[main_cluster] = final_main_vertex;
                detg_dump("snap_main_vertex_to_junction", *pr_graph);
                dup_stage_census("snap_main_vertex_to_junction", *pr_graph, *main_cluster);
            }

            pattern_algos.improve_vertex(*pr_graph, *main_cluster, final_main_vertex,
                                         vertices_in_long_muon, segments_in_long_muon,
                                         *track_fitter, m_dv, particle_data(), m_recomb_model,
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
                                                      *track_fitter, m_dv)) {
                std::cout << "After main vertex graph audit:" << final_main_vertex->fit().point << std::endl; pattern_algos.print_segs_info(*pr_graph, *main_cluster, 0);
                detg_dump("main_vertex_graph_audit", *pr_graph);
                dup_stage_census("main_vertex_graph_audit", *pr_graph, *main_cluster);
            }

            // doc sbnd_xin/docs/pr/84 round 2 (F3) -- inert unless
            // conn3_stitch_max > 0.  Bridges disconnected main-cluster
            // components whose closest approach to the reachable side is within
            // the radius, BEFORE clustering_points, so the piece is classified
            // conn-1 naturally instead of being promoted to a conn-3
            // "association" by shower_conn3_unreachable (which stays on as the
            // backstop for wider gaps).
            if (pattern_algos.stitch_disconnected_main_cluster(*pr_graph, *main_cluster, final_main_vertex,
                                                               *track_fitter, m_dv)) {
                detg_dump("conn3_stitch", *pr_graph);
                dup_stage_census("conn3_stitch", *pr_graph, *main_cluster);
            }

            // doc sbnd_xin/docs/pr/51 round 4: diagnostic-only rough-path probe
            // -- inert unless rough_path_probe.  Runs right after the audit
            // block (on today's production graph when mvga is off; on the
            // audited graph when it is on) so its measurements match whichever
            // near-vertex state is actually being Bee-scanned.
            pattern_algos.rough_path_probe(*pr_graph, *main_cluster, final_main_vertex, *track_fitter, m_dv);

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
            dup_stage_census("improve_vertex", *pr_graph, *main_cluster);
            if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: improve_vertex + examine_direction took {} ms", MS(Clock::now() - t0).count());
            t0 = Clock::now();

            // doc pr/93 round 4 (straight_cont_cross_cluster): demote main-vertex
            // shower-trajectory stems that are cross-cluster continuations of
            // straight long tracks (pr/57 W-gap splits).  Here on purpose: after
            // examine_direction (all clusters' segments exist in the graph, the
            // trajectory pdg-11 stamp is written, orientations are final -- the
            // pass preserves dirsign) and before shower_clustering_with_nv (the
            // seeder consumes flags/pdg; any bridge request recorded here is
            // replayed inside it, after its entry clears).  Knob off => early
            // return => byte-identical.
            pattern_algos.demote_cross_cluster_straight_stems(*pr_graph, final_main_vertex,
                                                              particle_data(), m_recomb_model);

            // doc pr/83 r3 (sec 9.5 + the 359980 follow-up): non-main clusters
            // that went through find_proto_vertex -- as a swapped-out old main
            // (Mechanism C, 350935) or as a candidate that lost the main-cluster
            // contest without any swap (359980: dup on cluster 75, main is 21)
            // -- keep their segments in the final output but never receive a
            // duplicate-corridor pass.  One unscoped audit each, BEFORE
            // shower_clustering_with_nv consumes them, so the shower maps never
            // hold a segment this pass removes.  Sorted by cluster id for
            // determinism.  Knob off (default) => loop skipped => byte-identical.
            if (m_swap_orphan_dup_audit) {
                std::vector<Cluster*> audit_clusters(other_clusters.begin(), other_clusters.end());
                std::sort(audit_clusters.begin(), audit_clusters.end(),
                          [](Cluster* a, Cluster* b) { return a->get_cluster_id() < b->get_cluster_id(); });
                audit_clusters.erase(std::unique(audit_clusters.begin(), audit_clusters.end()),
                                     audit_clusters.end());
                for (Cluster* oc : audit_clusters) {
                    if (!oc || oc == main_cluster) continue;
                    pattern_algos.orphan_dup_audit(*pr_graph, *oc, *track_fitter, m_dv);
                }
                if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: orphan_dup_audit sweep took {} ms", MS(Clock::now() - t0).count());
                t0 = Clock::now();
            }

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
                                                    *track_fitter, m_dv, particle_data(),
                                                    m_recomb_model);

            std::cout << "After shower clustering with NV: " << std::endl; pattern_algos.print_segs_info(*pr_graph, *main_cluster, 0);
            detg_dump("shower_clustering_with_nv", *pr_graph);
            dup_stage_census("shower_clustering_with_nv", *pr_graph, *main_cluster);
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

            // doc 84 round 1 (P1).  A shower (re)typed |13| AFTER its
            // kinematics pass keeps the generic multi-track verdict
            // kenergy_range = 0 / kenergy_best = 0 (PRShower.cxx generic
            // path), because retyping never clears flag_kinematics -- the
            // dominant mechanism behind the 28/242 range_zero census (doc 84
            // sec 3.2; mcp2k 497311, 332.8 cm muon reported as charge).  With
            // the labels now settled, clear the flag on exactly those showers
            // and rerun the kinematics pass: they re-dispatch through
            // calculate_kinematics_long_muon, whose empty-chain fallback
            // (same knob) supplies the muon-typed member length.  Knob off =>
            // no pass, byte-identical.
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
                    pattern_algos.calculate_shower_kinematics(showers, vertices_in_long_muon,
                                                              segments_in_long_muon, *pr_graph,
                                                              *track_fitter, m_dv,
                                                              particle_data(), m_recomb_model);
                    SPDLOG_LOGGER_DEBUG(log,
                        "long_muon_range_empty_chain_fallback: recomputed kinematics for {} retyped muon shower(s)",
                        n_stale);
                }
            }

            // doc 84 round 2 (long_muon_cathode_bridge) -- docstring at the
            // helper above.  Runs after the labels and the round-1 fallback
            // recompute settle, so the |13| typing it keys on is final.  The
            // absorbed far half joins the SHOWER only, deliberately not
            // segments_in_long_muon: tagger features (numu muon_length) and
            // the MCS chain keep their legacy inputs, so nusel stays put and
            // the change is confined to kine/PF output.  Knob off => no pass,
            // byte-identical.
            if (m_long_muon_cathode_bridge) {
                const CathodeBridgeCfg cb_cfg{
                    m_long_muon_cathode_bridge_x * units::cm,
                    m_long_muon_cathode_bridge_xcut * units::cm,
                    m_long_muon_cathode_bridge_gap * units::cm,
                    m_long_muon_cathode_bridge_angle,
                    m_long_muon_cathode_bridge_lever * units::cm,
                    m_long_muon_cathode_bridge_track_partner,
                    m_long_muon_cathode_bridge_short_gap * units::cm,
                    m_long_muon_cathode_bridge_short_gap_angle,
                    m_long_muon_cathode_bridge_short_gap_len * units::cm};
                const int n_bridged = long_muon_cathode_bridge_pass(
                    pattern_algos, *pr_graph, final_main_vertex, showers,
                    map_segment_in_shower, particle_data(), m_recomb_model,
                    cb_cfg, &bridged_in_long_muon, log);
                if (n_bridged > 0) {
                    pattern_algos.update_shower_maps(showers, map_vertex_in_shower,
                                                     map_segment_in_shower,
                                                     map_vertex_to_shower,
                                                     used_shower_clusters);
                    pattern_algos.calculate_shower_kinematics(showers, vertices_in_long_muon,
                                                              segments_in_long_muon, *pr_graph,
                                                              *track_fitter, m_dv,
                                                              particle_data(), m_recomb_model);
                    SPDLOG_LOGGER_DEBUG(log,
                        "long_muon_cathode_bridge: {} bridge(s), shower maps + kinematics recomputed",
                        n_bridged);
                }
            }
        }


        // Initialize tagger features to their default values unconditionally —
        // even if no vertex was found the struct must be value-initialized.
        t0 = Clock::now();
        TaggerInfo tagger_info;
        pattern_algos.init_tagger_info(tagger_info);

        // doc pr/94 -- per-bundle identity + per-activity cosmic block.
        // Guarded so the legacy path leaves every field at its struct default
        // (-1 / empty): the branches that expose them are booked only when
        // tagger_output's own nu_per_bundle knob is on, but TaggerInfo also
        // travels to PrDisplayDump and the Bee producer, so leaving it
        // untouched keeps those byte-identical too.
        if (m_nu_per_bundle) {
            tagger_info.cluster_id        = main_cluster->get_cluster_id();
            tagger_info.matched_flash_gid = candidates[nu_index].gid;
            tagger_info.nu_index          = static_cast<int>(nu_index);
            for (const auto& a : candidates[nu_index].acts) {
                tagger_info.act_cluster_id.push_back(a.cluster_id);
                tagger_info.act_length_cm.push_back(a.length_cm);
                tagger_info.act_is_selected.push_back(a.is_selected);
                tagger_info.act_is_demoted.push_back(a.is_demoted);
                tagger_info.act_tgm.push_back(a.tgm);
                tagger_info.act_stm.push_back(a.stm);
                tagger_info.act_fc.push_back(a.fc);
                tagger_info.act_lm.push_back(a.lm);
                tagger_info.act_evaluated.push_back(a.evaluated);
            }
        }

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
        // doc pr/92 -- ids of stray satellite showers dropped from the kine
        // tree; stashed into TrackFitting below UNCONDITIONALLY (replace
        // semantics) so no-vertex or knob-off events reset it to empty.
        std::set<int> dropped_sat_ids;
        if (final_main_vertex) {
            kine_info = pattern_algos.fill_kine_tree(
                final_main_vertex, showers, pio_kine,
                *pr_graph, *track_fitter, m_dv,
                m_geom_helper,          // nullptr when clus_geom_helper is not configured
                particle_data(), m_recomb_model,
                pi0_showers,        // pr/92: pi0-paired showers are drop-protected
                &dropped_sat_ids);
        }
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: fill_kine_tree took {} ms", MS(Clock::now() - t0).count());
        t0 = Clock::now(); // finalize block

        // Mark the main neutrino vertex and store neutrino results in TrackFitting
        // so that downstream consumers (e.g., Bee particle-flow output in MultiAlgBlobClustering)
        // can access them without re-running pattern recognition.
        if (final_main_vertex) {
            final_main_vertex->set_flags(PR::VertexFlags::kNeutrinoVertex);
        }
        track_fitter->set_pi0_data(pi0_showers, map_shower_pio_id, map_pio_id_showers, map_pio_id_mass);
        // doc pr/92 -- unconditional stash (empty when knob off / no vertex).
        track_fitter->set_dropped_satellite_shower_ids(std::move(dropped_sat_ids));
        track_fitter->set_main_vertex(final_main_vertex);
        track_fitter->set_showers(showers);
        // doc pr/94 -- T_kine[i] carries the same identity as T_tagger[i], so
        // the sync check can VERIFY the pairing instead of assuming position.
        if (m_nu_per_bundle) {
            kine_info.cluster_id        = main_cluster->get_cluster_id();
            kine_info.matched_flash_gid = candidates[nu_index].gid;
            kine_info.nu_index          = static_cast<int>(nu_index);
        }
        // doc 80 round 2: MCS muon momentum, once per bundle, default OFF.
        // Placed here so all three muon_source modes see what they need --
        // segments_in_long_muon is function-local to this scope (sec 4.2) --
        // and the result rides the existing unconditional set_kine_info store
        // below.  visit() only: the dual-chain off-pass
        // (run_dual_chain_off_pass) does not run MCS, so its rows keep -1.
        if (m_mcs.enable && final_main_vertex) {
            PR::mcs_fill_kine(kine_info, *pr_graph, segments_in_long_muon,
                              bridged_in_long_muon,
                              final_main_vertex, beam_gate, m_mcs,
                              particle_data(), m_recomb_model, log);
        }
        track_fitter->set_kine_info(kine_info);
        track_fitter->set_tagger_info(tagger_info);

        // doc pr/94 §10.1 -- the machine-readable "what was actually written"
        // sentinel, one line per emitted T_tagger/T_kine row.  It is emitted
        // HERE, not at selection time, because main_cluster can be repointed
        // between the two by swap_main_cluster (the DL and traditional overall
        // -vertex paths both may), so the selection line's cluster id is not
        // necessarily the id the row carries.  The sync check joins on this.
        if (m_nu_per_bundle) {
            SPDLOG_LOGGER_INFO(log, "TaggerCheckNeutrino: [nu_per_bundle] ROW {} gid {} cluster {} vertex ({:.4f}, {:.4f}, {:.4f}) cm Enu {:.4f} acts {}",
                               nu_index, candidates[nu_index].gid,
                               main_cluster->get_cluster_id(),
                               kine_info.kine_nu_x_corr, kine_info.kine_nu_y_corr,
                               kine_info.kine_nu_z_corr, kine_info.kine_reco_Enu,
                               tagger_info.act_cluster_id.size());
        }

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
            track_fitter->set_vertex_scoreboard(board);
        }

        // Merge every per-cluster fill_fitted_charge_2d snapshot into the flat
        // map that UbooneMagnifyTrackingVisitor::write_proj_data reads, so that
        // T_proj_data contains cells for all beam-flash clusters, not just the
        // last cluster fit by pattern recognition.
        track_fitter->assemble_fitted_charge_2d();

        // Store TrackFitting in the grouping for later access by bee output and tracking sink
        // doc pr/94 -- publish.  The UNNAMED slot keeps pointing at the first
        // candidate so every existing consumer (Bee PF, the magnify tracking
        // sink, the BDT scorers, tagger_output) is byte-identical in legacy
        // mode, where there is exactly one.  Per-bundle consumers walk the
        // "nu<i>" named slots until one comes back null.
        if (nu_index == 0) grouping.set_track_fitting(track_fitter);
        if (m_nu_per_bundle) {
            grouping.set_track_fitting("nu" + std::to_string(nu_index), track_fitter);
        }
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: finalize took {} ms", MS(Clock::now() - t0).count());
    }

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
            "addseg={} addseg_reentry={} ep_mismatch={} "
            "pca_calls={} pca_moved={} pca_mean_cm={:.4f} pca_max_cm={:.4f} "
            "oseg_proto={} oseg_relaxed_only={} oseg_reject={} "
            "oseg_iso_drop={} oseg_iso_keep={} "
            "knobs[fit_exclusion={} oov_parity={} local_pca={} relaxed_accept={}]",
            pa.oov_isochronous.load(), pa.oov_dead_scan.load(), pa.oov_unique_scan.load(),
            pa.add_segment_calls.load(), pa.add_segment_reentry.load(),
            pa.endpoint_mismatch.load(),
            pa.pca_refine_calls.load(), moved,
            moved ? (double(pa.pca_move_um_sum.load()) / moved) * units::um / units::cm : 0.0,
            double(pa.pca_move_um_max.load()) * units::um / units::cm,
            pa.oseg_accept_proto.load(), pa.oseg_accept_relaxed.load(), pa.oseg_reject.load(),
            pa.oseg_isolated_drop.load(), pa.oseg_isolated_keep.load(),
            m_fit_exclusion, m_oov_prototype_parity,
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

// doc sbnd_xin/docs/pr/112 sec 11 -- the dual chain's exclusion-free OFF pass.
//
// A DUPLICATE of visit()'s production stage sequence from the main-cluster
// initial PR through the vertex refinement block (stitch_disconnected_main_
// cluster), with fit_exclusion=false and everything non-essential removed:
// no scoreboard/harvest, no detg_dump / dup_stage_census / std::cout, no
// rough_path_probe, no env dumps.  Deliberately NOT a shared helper: the
// production block stays byte-for-byte (CLAUDE.md M10, doc sec 5.7.3).  A
// future round that inserts a stage into the production sequence must insert
// it here too, or the two chains silently stop being comparable.
//
// Mode-dependent stop: "voxels" (and "union" with vtx_weight 0) needs only
// the OFF graph's SCN top-K, so it stops after deghosting + one inference;
// "snap" runs the full vertex determination and refinement, because the
// measured +4/+5 was read off the OFF chain's SHIPPED vertex (sec 5.7.1).
//
// Facade residue: swap_main_cluster (reachable through the DL and traditional
// vertex paths) writes Flags::main_cluster.  Every cluster of the bundle is
// snapshotted before and restored after, unconditionally.
void TaggerCheckNeutrino::run_dual_chain_off_pass(const PR::PatternAlgorithms& prod_algos,
                                                  const TrackFitting& prod_fitter,
                                                  Cluster* main_cluster_in,
                                                  const std::vector<Cluster*>& other_clusters_in,
                                                  bool cov_defer_active,
                                                  PR::DualChainHint& hint) const
{
    using Clock = std::chrono::steady_clock;
    using MS = std::chrono::duration<double, std::milli>;
    auto t_total = Clock::now();
    auto t0 = Clock::now();
    const auto lap = [&](const char* what) {
        if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: [dual-off] {} took {} ms", what, MS(Clock::now() - t0).count());
        t0 = Clock::now();
    };

    hint.has_vertex = false;
    hint.voxels.clear();
    hint.n_candidates = 0;
    if (!main_cluster_in) return;

    // Facade flag snapshot/restore (see the function comment).
    std::vector<std::pair<Cluster*, bool>> flag_snapshot;
    flag_snapshot.emplace_back(main_cluster_in, main_cluster_in->get_flag(Flags::main_cluster));
    for (auto* c : other_clusters_in) if (c) flag_snapshot.emplace_back(c, c->get_flag(Flags::main_cluster));
    struct FlagRestore {
        std::vector<std::pair<Cluster*, bool>>& snap;
        ~FlagRestore() { for (auto& cv : snap) cv.first->set_flag(Flags::main_cluster, cv.second ? 1 : 0); }
    } flag_restore{flag_snapshot};

    // Local copies: the vertex paths may repoint main_cluster / edit other_clusters.
    Cluster* main_cluster = main_cluster_in;
    std::vector<Cluster*> other_clusters = other_clusters_in;

    // Own fitter, the visit() nu_index>0 recipe (parameters copied from the
    // production fitter AFTER its per-visit set_parameter calls landed).
    auto track_fitter = std::make_shared<TrackFitting>(TrackFittingPresets::create_with_current_values());
    track_fitter->set_parameters(prod_fitter.get_parameters());
    track_fitter->set_perf(prod_fitter.get_perf());
    track_fitter->set_detector_volume(m_dv);
    track_fitter->set_pc_transforms(m_pcts);
    {
        std::vector<Cluster*> clusters_to_preload;
        clusters_to_preload.push_back(main_cluster);
        for (auto* c : other_clusters) clusters_to_preload.push_back(c);
        track_fitter->preload_clusters(clusters_to_preload);
    }
    auto pr_graph = std::make_shared<PR::Graph>();
    track_fitter->add_graph(pr_graph);

    // Own PatternAlgorithms: identical configuration, exclusion off, quiet.
    PR::PatternAlgorithms pattern_algos = prod_algos;
    // WCT_DUAL_CHAIN_OFF_EXCL=1: duplicate-fidelity debug switch (sec 5.7.3) --
    // the OFF pass keeps production's exclusion setting, so its vertex must
    // equal production's event by event.  Never a cfg key.
    static const bool keep_excl = std::getenv("WCT_DUAL_CHAIN_OFF_EXCL") != nullptr;
    if (!keep_excl) pattern_algos.m_fit_exclusion = false;
    pattern_algos.m_dl_vtx_cloud_no_exclusion = false;
    pattern_algos.m_vertex_scoreboard = false;
    pattern_algos.m_vtx_harvest = false;
    pattern_algos.m_vtx_board.clear();

    IndexedVertexSet vertices_in_long_muon;
    IndexedSegmentSet segments_in_long_muon;
    VertexPtr main_vertex = nullptr;
    ClusterVertexMap map_cluster_main_vertices;
    lap("setup");

    // ---- main cluster initial PR (production :2145-2183 at b5c9f43a)
    {
        if (cov_defer_active) track_fitter->set_parameter("fit_blob_coverage", -1);
        pattern_algos.find_proto_vertex(*pr_graph, *main_cluster, *track_fitter, m_dv, true, 2, true, particle_data());
        if (cov_defer_active) track_fitter->set_parameter("fit_blob_coverage", m_fit_blob_coverage);
        pattern_algos.clustering_points(*pr_graph, *main_cluster, m_dv);
        pattern_algos.separate_track_shower(*pr_graph, *main_cluster);
        pattern_algos.determine_direction(*pr_graph, *main_cluster, particle_data(), m_recomb_model);
        pattern_algos.shower_determining_in_main_cluster(*pr_graph, *main_cluster, particle_data(), m_recomb_model, m_dv);
        pattern_algos.determine_main_vertex(*pr_graph, *main_cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, *track_fitter, m_dv, particle_data(), m_recomb_model);
        pattern_algos.reassociate_cluster_orphans(*pr_graph, *main_cluster, m_dv);
        if (main_vertex != nullptr) {
            map_cluster_main_vertices[main_cluster] = main_vertex;
            main_vertex = nullptr;
        }
    }
    lap("main_cluster initial PR");

    // ---- other clusters + deghosting (production :2196-2244)
    if (!other_clusters.empty()) {
        for (auto* cluster : other_clusters) {
            if (cluster->get_length() > 6 * units::cm) {
                pattern_algos.find_proto_vertex(*pr_graph, *cluster, *track_fitter, m_dv, true, 2, false);
            } else {
                if (!pattern_algos.find_proto_vertex(*pr_graph, *cluster, *track_fitter, m_dv, false, 1, false)) {
                    pattern_algos.init_point_segment(*pr_graph, *cluster, *track_fitter, m_dv);
                }
            }
            pattern_algos.clustering_points(*pr_graph, *cluster, m_dv);
            pattern_algos.separate_track_shower(*pr_graph, *cluster);
            pattern_algos.determine_direction(*pr_graph, *cluster, particle_data(), m_recomb_model);
            pattern_algos.shower_determining_in_main_cluster(*pr_graph, *cluster, particle_data(), m_recomb_model, m_dv);
            pattern_algos.determine_main_vertex(*pr_graph, *cluster, main_vertex, vertices_in_long_muon, segments_in_long_muon, *track_fitter, m_dv, particle_data(), m_recomb_model);
            pattern_algos.reassociate_cluster_orphans(*pr_graph, *cluster, m_dv);
            if (main_vertex != nullptr) {
                map_cluster_main_vertices[cluster] = main_vertex;
                main_vertex = nullptr;
            }
        }
        lap("other_clusters PR");
        std::vector<Cluster*> all_clusters;
        all_clusters.push_back(main_cluster);
        all_clusters.insert(all_clusters.end(), other_clusters.begin(), other_clusters.end());
        pattern_algos.deghosting(*pr_graph, map_cluster_main_vertices, all_clusters, *track_fitter, m_dv);
        lap("deghosting");
    }

    // ---- the OFF graph's own SCN top-K (voxels / union modes)
    const bool want_voxels = (hint.mode == "voxels" || hint.mode == "union");
    const bool want_vertex = (hint.mode == "snap") || (hint.mode == "union" && hint.vtx_weight != 0.0);
    if (want_voxels && !m_dl_weights.empty()) {
        hint.voxels = pattern_algos.dual_chain_scn_voxels(*pr_graph, m_dl_weights, m_dQdx_scale, m_dQdx_offset,
                                                          m_dl_vtx_top_k, hint.n_candidates);
        lap("scn voxels");
    }
    if (!want_vertex) {
        if (hint.n_candidates == 0) {
            for (const auto& nd : PR::ordered_nodes(*pr_graph)) if ((*pr_graph)[nd].vertex) ++hint.n_candidates;
        }
        hint.off_ms = MS(Clock::now() - t_total).count();
        SPDLOG_LOGGER_INFO(log, "dual_chain: OFF pass ({}) {} candidates, {} voxels, {:.0f} ms",
                           hint.mode, hint.n_candidates, hint.voxels.size() / 4, hint.off_ms);
        return;
    }

    // ---- overall main vertex (production :2250-2304)
    VertexPtr final_main_vertex = nullptr;
    bool flag_dl_changed = false;
    if (!m_dl_weights.empty()) {
        flag_dl_changed = pattern_algos.determine_overall_main_vertex_DL(
            *pr_graph, map_cluster_main_vertices, main_cluster, other_clusters,
            vertices_in_long_muon, segments_in_long_muon,
            *track_fitter, m_dv, particle_data(), m_recomb_model,
            m_dl_weights, m_dl_vtx_cut, m_dQdx_scale, m_dQdx_offset,
            m_dl_vtx_rerank, m_dl_vtx_top_k, m_dl_vtx_min_accept_score,
            m_dl_vtx_score_scale);
    }
    if (!flag_dl_changed) {
        ClusterVertexMap map_copy = map_cluster_main_vertices;
        Cluster* mc_copy = main_cluster;
        final_main_vertex = pattern_algos.determine_overall_main_vertex(
            *pr_graph, map_copy, mc_copy, other_clusters,
            vertices_in_long_muon, segments_in_long_muon,
            *track_fitter, m_dv, particle_data(), m_recomb_model, true);
        if (mc_copy != main_cluster && m_main_vertex_swap_apply) {
            main_cluster = mc_copy;
            map_cluster_main_vertices = map_copy;
        }
        if (final_main_vertex) map_cluster_main_vertices[main_cluster] = final_main_vertex;
    }
    {
        auto it = map_cluster_main_vertices.find(main_cluster);
        if (it != map_cluster_main_vertices.end()) final_main_vertex = it->second;
    }
    lap("overall main vertex");

    // ---- refinement block (production :2310-2373): part of vertex
    // determination -- the measured gain requires it (sec 5.7.1).
    if (final_main_vertex) {
        if (pattern_algos.snap_main_vertex_to_kink(*pr_graph, *main_cluster, final_main_vertex,
                                                   *track_fitter, m_dv, particle_data(), m_recomb_model)) {
            map_cluster_main_vertices[main_cluster] = final_main_vertex;
        }
        if (pattern_algos.snap_main_vertex_to_junction(*pr_graph, *main_cluster, final_main_vertex)) {
            map_cluster_main_vertices[main_cluster] = final_main_vertex;
        }
        pattern_algos.improve_vertex(*pr_graph, *main_cluster, final_main_vertex,
                                     vertices_in_long_muon, segments_in_long_muon,
                                     *track_fitter, m_dv, particle_data(), m_recomb_model,
                                     true, true);
        map_cluster_main_vertices[main_cluster] = final_main_vertex;
        pattern_algos.main_vertex_graph_audit(*pr_graph, *main_cluster, final_main_vertex, *track_fitter, m_dv);
        pattern_algos.stitch_disconnected_main_cluster(*pr_graph, *main_cluster, final_main_vertex, *track_fitter, m_dv);
    }
    lap("refinement");

    if (final_main_vertex) {
        hint.has_vertex = true;
        hint.vertex = final_main_vertex->fit().valid() ? final_main_vertex->fit().point : final_main_vertex->wcpt().point;
        hint.vertex_cluster_id = main_cluster ? main_cluster->get_cluster_id() : -1;
    }
    if (hint.n_candidates == 0) {
        for (const auto& nd : PR::ordered_nodes(*pr_graph)) if ((*pr_graph)[nd].vertex) ++hint.n_candidates;
    }
    hint.off_ms = MS(Clock::now() - t_total).count();
    SPDLOG_LOGGER_INFO(log, "dual_chain: OFF pass ({}) vertex={} ({:.2f},{:.2f},{:.2f}) cm cluster {} {} candidates, {} voxels, {:.0f} ms",
                       hint.mode, hint.has_vertex,
                       hint.vertex.x() / units::cm, hint.vertex.y() / units::cm, hint.vertex.z() / units::cm,
                       hint.vertex_cluster_id, hint.n_candidates, hint.voxels.size() / 4, hint.off_ms);
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

