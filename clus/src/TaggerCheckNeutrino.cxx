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
    m_cathode_x          = get(config, "cathode_x",          m_cathode_x);           // cm
    m_cathode_kink_xcut  = get(config, "cathode_kink_xcut",  m_cathode_kink_xcut);   // cm
    m_shower_topo_demote_len = get(config, "shower_topo_demote_len", m_shower_topo_demote_len);  // cm
    m_iso_endpoint               = get(config, "iso_endpoint",               m_iso_endpoint);
    m_iso_endpoint_min_length    = get(config, "iso_endpoint_min_length",    m_iso_endpoint_min_length);     // cm
    m_iso_endpoint_max_xext      = get(config, "iso_endpoint_max_xext",      m_iso_endpoint_max_xext);       // cm
    m_iso_endpoint_xext_frac     = get(config, "iso_endpoint_xext_frac",     m_iso_endpoint_xext_frac);
    m_iso_endpoint_xext_quantile = get(config, "iso_endpoint_xext_quantile", m_iso_endpoint_xext_quantile);
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
        }
    }
    m_dl_vtx_cut              = get(config, "dl_vtx_cut",              m_dl_vtx_cut);
    m_dQdx_scale              = get(config, "dQdx_scale",              m_dQdx_scale);
    m_dQdx_offset             = get(config, "dQdx_offset",             m_dQdx_offset);
    m_dl_vtx_rerank           = get(config, "dl_vtx_rerank",           m_dl_vtx_rerank);
    m_dl_vtx_top_k            = get(config, "dl_vtx_top_k",            m_dl_vtx_top_k);
    m_dl_vtx_min_accept_score = get(config, "dl_vtx_min_accept_score", m_dl_vtx_min_accept_score);
    m_dl_vtx_score_scale      = get(config, "dl_vtx_score_scale",      m_dl_vtx_score_scale);
    m_beam_window_low         = get(config, "beam_window_low",         m_beam_window_low);
    m_beam_window_high        = get(config, "beam_window_high",        m_beam_window_high);
    m_nu_skip_cosmic          = get(config, "nu_skip_cosmic",          m_nu_skip_cosmic);
    m_nu_skip_cosmic_bundle   = get(config, "nu_skip_cosmic_bundle",   m_nu_skip_cosmic_bundle);
    m_nu_skip_cosmic_bundle_min_length = get(config, "nu_skip_cosmic_bundle_min_length", m_nu_skip_cosmic_bundle_min_length);  // cm
    m_skip_cosmic_companions  = get(config, "skip_cosmic_companions",  m_skip_cosmic_companions);
    m_cosmic_companion_min_length = get(config, "cosmic_companion_min_length", m_cosmic_companion_min_length);  // cm

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
    cfg["cathode_x"]         = m_cathode_x;          // cm, T0-corrected frame
    cfg["cathode_kink_xcut"] = m_cathode_kink_xcut;  // cm; 0 = legacy (the kink search sees every fit point)
    cfg["shower_topo_demote_len"] = m_shower_topo_demote_len;  // cm; 0 = legacy (long segments stay eligible for kShowerTopology)
    cfg["iso_endpoint"]               = m_iso_endpoint;                // false = legacy wire-footprint boundary endpoints
    cfg["iso_endpoint_min_length"]    = m_iso_endpoint_min_length;     // cm
    cfg["iso_endpoint_max_xext"]      = m_iso_endpoint_max_xext;       // cm
    cfg["iso_endpoint_xext_frac"]     = m_iso_endpoint_xext_frac;
    cfg["iso_endpoint_xext_quantile"] = m_iso_endpoint_xext_quantile;
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
    cfg["clus_geom_helper"] = ""; // empty = no SCE vertex correction
    cfg["beam_window_low"] = m_beam_window_low;   // beam window on cluster_t0; low >= high disables the
    cfg["beam_window_high"] = m_beam_window_high; // gate (uBooNE single-main selection).
    cfg["nu_skip_cosmic"] = m_nu_skip_cosmic;     // beam-gate only: skip in-window mains with flag_TGM/flag_STM/lm_flag>0
    cfg["nu_skip_cosmic_bundle"] = m_nu_skip_cosmic_bundle;  // that verdict vetoes the whole matched_flash_gid bundle
    cfg["nu_skip_cosmic_bundle_min_length"] = m_nu_skip_cosmic_bundle_min_length;  // cm; > 0 spares untagged bundle-mates at least this long (docs/pr/16 design A); 0 = veto all
    cfg["skip_cosmic_companions"] = m_skip_cosmic_companions;          // doc pr/20 I P4; drop a TGM/STM-tagged companion from other_clusters
    cfg["cosmic_companion_min_length"] = m_cosmic_companion_min_length;  // cm; a tagged companion shorter than this stays in regardless


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
    pattern_algos.m_cathode_x         = m_cathode_x * units::cm;          // cm -> internal
    pattern_algos.m_cathode_kink_xcut = m_cathode_kink_xcut * units::cm;  // cm -> internal
    pattern_algos.m_shower_topo_demote_len = m_shower_topo_demote_len * units::cm;  // cm -> internal
    pattern_algos.m_iso_endpoint               = m_iso_endpoint;
    pattern_algos.m_iso_endpoint_min_length    = m_iso_endpoint_min_length * units::cm;  // cm -> internal
    pattern_algos.m_iso_endpoint_max_xext      = m_iso_endpoint_max_xext * units::cm;    // cm -> internal
    pattern_algos.m_iso_endpoint_xext_frac     = m_iso_endpoint_xext_frac;
    pattern_algos.m_iso_endpoint_xext_quantile = m_iso_endpoint_xext_quantile;
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
    pattern_algos.m_kine_charge.w_value             = m_kine_w_value;
    // Muon dQ/dx-vs-length envelope: c0/c1/power dimensionless, pivot cm -> internal.
    pattern_algos.m_muon_dqdx_curve = {m_muon_dqdx_curve[0], m_muon_dqdx_curve[1],
                                       m_muon_dqdx_curve[2] * units::cm, m_muon_dqdx_curve[3]};
    // Single-photon stem dE/dx conversion; the cut narrows to float so the
    // default compares bit-identically to the legacy 2.3f literal.
    pattern_algos.m_sp_dedx_use_recomb_model = m_sp_dedx_use_recomb_model;
    pattern_algos.m_sp_mean_dedx_cut         = static_cast<float>(m_sp_mean_dedx_cut);
    pattern_algos.m_recomb_model             = m_recomb_model;
    m_track_fitter->set_perf(m_perf);

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

    {
        // initial pattern recognitions
        pattern_algos.find_proto_vertex(*pr_graph, *main_cluster, *m_track_fitter, m_dv, true, 2, true);
        detg_dump("main:find_proto_vertex", *pr_graph);

        // shower related operations
        pattern_algos.clustering_points(*pr_graph, *main_cluster, m_dv);
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
            m_dl_vtx_score_scale);
    }
    if (!flag_dl_changed) {
        final_main_vertex = pattern_algos.determine_overall_main_vertex(
            *pr_graph, map_cluster_main_vertices, main_cluster, other_clusters,
            vertices_in_long_muon, segments_in_long_muon,
            *m_track_fitter, m_dv, particle_data(), m_recomb_model, true);
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
   

        pattern_algos.improve_vertex(*pr_graph, *main_cluster, final_main_vertex,
                                     vertices_in_long_muon, segments_in_long_muon,
                                     *m_track_fitter, m_dv, particle_data(), m_recomb_model,
                                     true, true);
        // improve_vertex may update final_main_vertex pointer; sync back to map
        map_cluster_main_vertices[main_cluster] = final_main_vertex;

        std::cout << "After improve vertex:" << final_main_vertex->fit().point << std::endl; pattern_algos.print_segs_info(*pr_graph, *main_cluster, 0);

        pattern_algos.clustering_points(*pr_graph, *main_cluster, m_dv);

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

        pattern_algos.singlephoton_tagger(*pr_graph, main_cluster,
                                          final_main_vertex,
                                          showers,
                                          map_vertex_to_shower,
                                          map_shower_pio_id,
                                          map_pio_id_showers,
                                          map_pio_id_mass,
                                          m_dv,
                                          tagger_info);
    }

    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: taggers took {} ms", MS(Clock::now() - t0).count());
    t0 = Clock::now();

    // Compute match_isFC: 1 if the main cluster is fully contained inside the
    // fiducial volume, 0 otherwise.  Uses the same two-round boundary check as
    // TaggerCheckSTM so the definition is consistent across both users.
    if (main_cluster) {
        auto fc_result = Facade::cluster_fc_check(*main_cluster, m_dv);
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

    // Merge every per-cluster fill_fitted_charge_2d snapshot into the flat
    // map that UbooneMagnifyTrackingVisitor::write_proj_data reads, so that
    // T_proj_data contains cells for all beam-flash clusters, not just the
    // last cluster fit by pattern recognition.
    m_track_fitter->assemble_fitted_charge_2d();

    // Store TrackFitting in the grouping for later access by bee output and tracking sink
    grouping.set_track_fitting(m_track_fitter);
    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: finalize took {} ms", MS(Clock::now() - t0).count());
    if (m_perf) SPDLOG_LOGGER_DEBUG(log, "TaggerCheckNeutrino timing: visit() TOTAL took {} ms", MS(Clock::now() - t_total).count());
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

