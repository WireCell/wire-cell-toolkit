#include "WireCellMatch/QLMatching.h"
#include "WireCellMatch/Opflash.h"

#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/TensorDMdataset.h"
#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellClus/Facade.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/ExecMon.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Ress.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Units.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <set>

WIRECELL_FACTORY(QLMatching,
                 WireCell::Match::QLMatching,
                 WireCell::INamed,
                 WireCell::ITensorSetFanin,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Match;
using namespace WireCell::Clus::Facade;

// Stride used to build a globally-unique flash id (gid = anode_ident*stride +
// per-APA flash row) so the matched flash association survives the per-APA ->
// all-APA pctree merge and stays unambiguous across APAs. Far larger than any
// realistic per-APA flash count.
namespace { constexpr int kFlashGidStride = 1000000; }

// ---- Per-phase profiling helpers (logging only; outputs bit-identical) ----
// Lightweight wall-clock split of operator() so we can attribute the per-event
// `QLMatching timing: took ...` total to its phases (prefit sub-steps, the two
// LASSO rounds, etc.). All emitted at debug level under the "QLtiming" tag, so a
// non-debug run pays only a few steady_clock reads and computes the same result.
namespace {
    using wallclock = std::chrono::steady_clock;
    inline double ms_since(const wallclock::time_point& t0)
    {
        return std::chrono::duration<double, std::milli>(wallclock::now() - t0).count();
    }
}

// ---- Per-APA tree merge (multi-APA path) ----
// Verbatim port of clus/src/PointTreeMerging.cxx so the joint (fanin) path
// reproduces, byte-for-byte, the standalone PointTreeMerging it replaces: align
// local-PC keys across all nodes of the merged tree, and concatenate the opted-in
// root PCs + take the children of each source tree into the primary (input[0]).
namespace {
    using WireCell::PointCloud::Tree::Points;

    size_t normalize_pctree_local_pcs(Points::node_t* root)
    {
        if (!root) return 0;
        std::vector<Points::node_t*> nodes;
        for (auto& noderef : root->depth()) nodes.push_back(&noderef);

        std::map<std::string, std::map<std::string, WireCell::PointCloud::Array>> templates;
        for (const auto* node : nodes) {
            for (const auto& [pcname, ds] : node->value.local_pcs()) {
                auto& pc_templates = templates[pcname];
                for (const auto& key : ds.keys()) {
                    if (pc_templates.count(key)) continue;
                    pc_templates.emplace(key, *ds.get(key));
                }
            }
        }
        size_t nadded = 0;
        for (auto* node : nodes) {
            for (auto& [pcname, ds] : node->value.local_pcs()) {
                const auto it = templates.find(pcname);
                if (it == templates.end()) continue;
                for (const auto& [key, arr_template] : it->second) {
                    if (ds.has(key)) continue;
                    ds.add(key, arr_template.zeros_like(ds.size_major()));
                    ++nadded;
                }
            }
        }
        return nadded;
    }

    void merge_pct(Points::node_t* tgt, Points::node_t* src,
                   const std::set<std::string>& root_pcs_to_merge)
    {
        if (!src) return;
        // Merge selected root-node local PCs (concatenate across inputs). NOTE:
        // local_pcs() returns a reference, so bind by reference or the merge is a
        // silent no-op. Only names in root_pcs_to_merge are merged; everything else
        // (flash/light/flashlight, per-anode ctpc_a*, ...) is dropped from the
        // source roots, exactly as the standalone PointTreeMerging did.
        auto& tgt_pc = tgt->value.local_pcs();
        for (const auto& src_pc : src->value.local_pcs()) {
            const auto& name = src_pc.first;
            if (root_pcs_to_merge.find(name) == root_pcs_to_merge.end()) continue;
            if (tgt_pc.find(name) == tgt_pc.end()) tgt_pc.emplace(name, src_pc.second);
            else tgt_pc[name].append(src_pc.second);
        }
        bool notify_value = true;
        tgt->take_children(*src, notify_value);
    }
}

QLMatching::QLMatching() : Aux::Logger("QLMatching", "match") {}
QLMatching::~QLMatching() = default;

void QLMatching::configure(const WireCell::Configuration& cfg)
{
    // Anodes: a list ("anodes") enables the joint multi-APA path (one input port
    // per anode); the historical single "anode" yields one entry => multiplicity 1
    // and the exact single-input behavior. m_anode keeps the first (back-compat).
    m_anodes.clear();
    if (cfg.isMember("anodes") && cfg["anodes"].isArray()) {
        for (const auto& a : cfg["anodes"]) m_anodes.push_back(Factory::find_tn<IAnodePlane>(a.asString()));
    }
    else {
        m_anodes.push_back(Factory::find_tn<IAnodePlane>(cfg["anode"].asString()));
    }
    m_anode = m_anodes.front();
    m_multiplicity = m_anodes.size();
    // Optional extra anodes to register on each run's Grouping (for multi-APA
    // drift-side group trees), stored per input port. "grouping_anodes" may be a flat
    // array (same list for every input -- single-side node) or an array of arrays (one
    // list per input -- JOINT multi-side node). Empty => only the run's own anode
    // (SBND-identical).
    m_grouping_anodes.assign(m_multiplicity, {});
    if (cfg.isMember("grouping_anodes") && cfg["grouping_anodes"].isArray()
        && !cfg["grouping_anodes"].empty()) {
        const auto& ga = cfg["grouping_anodes"];
        if (ga[0].isArray()) {   // per-input: one anode list per port
            for (std::size_t k = 0; k < m_multiplicity && k < ga.size(); ++k)
                for (const auto& a : ga[(Json::ArrayIndex)k])
                    m_grouping_anodes[k].push_back(Factory::find_tn<IAnodePlane>(a.asString()));
        }
        else {                   // flat: same list applied to every input (back-compat)
            std::vector<IAnodePlane::pointer> flat;
            for (const auto& a : ga) flat.push_back(Factory::find_tn<IAnodePlane>(a.asString()));
            for (std::size_t k = 0; k < m_multiplicity; ++k) m_grouping_anodes[k] = flat;
        }
    }
    if (cfg.isMember("root_pcs_to_merge") && cfg["root_pcs_to_merge"].isArray()) {
        m_root_pcs_to_merge.clear();
        for (const auto& jname : cfg["root_pcs_to_merge"]) m_root_pcs_to_merge.insert(jname.asString());
    }

    m_dv    = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());

    m_inpath        = get(cfg, "inpath", m_inpath);
    m_outpath       = get(cfg, "outpath", m_outpath);
    m_cluster_t0    = get(cfg, "cluster_t0", m_cluster_t0);
    m_semimodel_file = get(cfg, "semimodel_file", m_semimodel_file);
    m_light_model   = get(cfg, "light_model", m_light_model);
    m_photon_library_file = get(cfg, "photon_library_file", m_photon_library_file);
    m_photon_library_pos_tol = get(cfg, "photon_library_pos_tol", m_photon_library_pos_tol);
    m_calib_dump    = get(cfg, "calib_dump", m_calib_dump);
    m_flash_group_window = get(cfg, "flash_group_window", m_flash_group_window);
    m_cathode_diag  = get(cfg, "cathode_diag", m_cathode_diag);
    m_cathode_diag_radius = get(cfg, "cathode_diag_radius", m_cathode_diag_radius);
    m_xtpc_flag         = get(cfg, "xtpc_flag",         m_xtpc_flag);
    m_xtpc_dmax         = get(cfg, "xtpc_dmax",         m_xtpc_dmax);
    m_xtpc_dmax2        = get(cfg, "xtpc_dmax2",        m_xtpc_dmax2);
    m_xtpc_angle_max    = get(cfg, "xtpc_angle_max",    m_xtpc_angle_max);
    m_xtpc_hough_radius = get(cfg, "xtpc_hough_radius", m_xtpc_hough_radius);
    m_xtpc_joint_pin    = get(cfg, "xtpc_joint_pin",    m_xtpc_joint_pin);
    m_xtpc_pin_angle    = get(cfg, "xtpc_pin_angle",    m_xtpc_pin_angle);
    m_xtpc_cathode_tol   = get(cfg, "xtpc_cathode_tol",   m_xtpc_cathode_tol);
    m_xtpc_cathode_qfrac = get(cfg, "xtpc_cathode_qfrac", m_xtpc_cathode_qfrac);
    // xtpc / selection quality gates (039252 scan tuning, doc 19); legacy defaults.
    m_xtpc_pin_min_strength = get(cfg, "xtpc_pin_min_strength", m_xtpc_pin_min_strength);
    m_xtpc_sc1_light_gate   = get(cfg, "xtpc_sc1_light_gate",   m_xtpc_sc1_light_gate);
    m_xtpc_sc1_ks_max       = get(cfg, "xtpc_sc1_ks_max",       m_xtpc_sc1_ks_max);
    m_xtpc_sc1_c2n_max      = get(cfg, "xtpc_sc1_c2n_max",      m_xtpc_sc1_c2n_max);
    m_xtpc_cathode_ks_max   = get(cfg, "xtpc_cathode_ks_max",   m_xtpc_cathode_ks_max);
    m_postcull_unflagged    = get(cfg, "postcull_unflagged",    m_postcull_unflagged);
    m_postcull_ks_max       = get(cfg, "postcull_ks_max",       m_postcull_ks_max);
    m_postcull_c2n_max      = get(cfg, "postcull_c2n_max",      m_postcull_c2n_max);
    m_postcull_before_rescue = get(cfg, "postcull_before_rescue", m_postcull_before_rescue);
    m_postcull_wtrunc_overpred = get(cfg, "postcull_wtrunc_overpred", m_postcull_wtrunc_overpred);
    m_postcull_wtrunc_ratio_hi = get(cfg, "postcull_wtrunc_ratio_hi", m_postcull_wtrunc_ratio_hi);
    m_postcull_wtrunc_sat_frac = get(cfg, "postcull_wtrunc_sat_frac", m_postcull_wtrunc_sat_frac);
    m_postcull_pin_overpred    = get(cfg, "postcull_pin_overpred",    m_postcull_pin_overpred);
    m_postcull_pin_ratio_hi    = get(cfg, "postcull_pin_ratio_hi",    m_postcull_pin_ratio_hi);
    if (m_postcull_before_rescue && !m_postcull_unflagged) {
        log->warn("QLMatching: postcull_before_rescue needs postcull_unflagged "
                  "and will be inert");
    }
    if (m_xtpc_pin_min_strength > 0 || m_xtpc_sc1_light_gate ||
        m_xtpc_cathode_ks_max > 0 || m_postcull_unflagged) {
        log->debug("quality gates on: pin_min_strength={} sc1_light_gate={} "
                   "(ks<={} c2n<={}) cathode_ks_max={} postcull_unflagged={} "
                   "(ks<={} c2n<={})",
                   m_xtpc_pin_min_strength, m_xtpc_sc1_light_gate,
                   m_xtpc_sc1_ks_max, m_xtpc_sc1_c2n_max, m_xtpc_cathode_ks_max,
                   m_postcull_unflagged, m_postcull_ks_max, m_postcull_c2n_max);
    }

    m_pmts     = get(cfg, "pmts", m_pmts);
    m_data     = get(cfg, "data", m_data);
    m_beamonly = get(cfg, "beamonly", m_beamonly);

    if (cfg.isMember("active_opdet_types") && cfg["active_opdet_types"].isArray()) {
        m_active_opdet_types.clear();
        for (const auto& jt : cfg["active_opdet_types"]) m_active_opdet_types.push_back(jt.asInt());
    }

    if (cfg.isMember("ch_mask") && cfg["ch_mask"].isArray()) {
        m_ch_mask.clear();
        for (const auto& jch : cfg["ch_mask"]) m_ch_mask.push_back(jch.asInt());
    }

    m_auto_mask              = get(cfg, "auto_mask",              m_auto_mask);
    m_auto_mask_pe_low       = get(cfg, "auto_mask_pe_low",       m_auto_mask_pe_low);
    m_auto_mask_neighbors    = get(cfg, "auto_mask_neighbors",    m_auto_mask_neighbors);
    m_auto_mask_pe_bright    = get(cfg, "auto_mask_pe_bright",    m_auto_mask_pe_bright);
    m_auto_mask_min_contrast = get(cfg, "auto_mask_min_contrast", m_auto_mask_min_contrast);
    m_auto_mask_min_flash    = get(cfg, "auto_mask_min_flash",    m_auto_mask_min_flash);
    m_auto_mask_same_type    = get(cfg, "auto_mask_same_type",    m_auto_mask_same_type);

    // ---- PDVD / vertical-drift knobs (all default OFF => byte-identical). ----
    m_shared_flash      = get(cfg, "shared_flash",      m_shared_flash);
    m_opdet_all_volumes = get(cfg, "opdet_all_volumes", m_opdet_all_volumes);
    m_vd_surface_flags  = get(cfg, "vd_surface_flags",  m_vd_surface_flags);
    m_pd_wall_cushion   = get(cfg, "pd_wall_cushion",   m_pd_wall_cushion);
    if (cfg.isMember("pd_wall_channels_ylo") && cfg["pd_wall_channels_ylo"].isArray()) {
        m_pd_wall_channels_ylo.clear();
        for (const auto& jch : cfg["pd_wall_channels_ylo"]) m_pd_wall_channels_ylo.push_back(jch.asInt());
    }
    if (cfg.isMember("pd_wall_channels_yhi") && cfg["pd_wall_channels_yhi"].isArray()) {
        m_pd_wall_channels_yhi.clear();
        for (const auto& jch : cfg["pd_wall_channels_yhi"]) m_pd_wall_channels_yhi.push_back(jch.asInt());
    }
    // Per-input anode PD channels: array of arrays, one per input port (empty inner
    // list => that anode hosts no PDs). A flat array applies to every input.
    m_anode_pd_channels.assign(m_multiplicity, {});
    if (cfg.isMember("anode_pd_channels") && cfg["anode_pd_channels"].isArray()
        && !cfg["anode_pd_channels"].empty()) {
        const auto& apc = cfg["anode_pd_channels"];
        if (apc[0].isArray()) {
            for (std::size_t k = 0; k < m_multiplicity && k < apc.size(); ++k)
                for (const auto& jch : apc[(Json::ArrayIndex)k])
                    m_anode_pd_channels[k].push_back(jch.asInt());
        }
        else {
            std::vector<int> flat;
            for (const auto& jch : apc) flat.push_back(jch.asInt());
            for (std::size_t k = 0; k < m_multiplicity; ++k) m_anode_pd_channels[k] = flat;
        }
    }

    m_flash_minPE   = get(cfg, "flash_minPE",   m_flash_minPE);
    if (cfg.isMember("flash_sel_channels") && cfg["flash_sel_channels"].isArray()) {
        m_flash_sel_channels.clear();
        for (const auto& jch : cfg["flash_sel_channels"]) m_flash_sel_channels.push_back(jch.asInt());
    }
    m_flash_sel_minPE     = get(cfg, "flash_sel_minPE",     m_flash_sel_minPE);
    m_flash_sel_min_fired = get(cfg, "flash_sel_min_fired", m_flash_sel_min_fired);
    m_flash_sel_fired_pe  = get(cfg, "flash_sel_fired_pe",  m_flash_sel_fired_pe);
    m_flash_mintime = get(cfg, "flash_mintime", m_flash_mintime);
    m_flash_maxtime = get(cfg, "flash_maxtime", m_flash_maxtime);
    m_beam_mintime  = get(cfg, "beam_mintime",  m_beam_mintime);
    m_beam_maxtime  = get(cfg, "beam_maxtime",  m_beam_maxtime);
    if (m_beamonly) {
        m_flash_mintime = m_beam_mintime;
        m_flash_maxtime = m_beam_maxtime;
    }

    m_QtoL = get(cfg, "QtoL", m_QtoL);
    m_doReflectedLight = get(cfg, "doReflectedLight", m_doReflectedLight);
    m_tpc_face = get(cfg, "tpc_face", m_tpc_face);
    // Optional per-input face list for a JOINT multi-side node (PDHD -x face 0 / +x
    // face 1). Empty => m_tpc_face for every input (bit-identical).
    m_tpc_faces.clear();
    if (cfg.isMember("tpc_faces") && cfg["tpc_faces"].isArray())
        for (const auto& f : cfg["tpc_faces"]) m_tpc_faces.push_back(f.asInt());
    m_tpc_extra_faces.clear();
    if (cfg.isMember("tpc_extra_faces") && cfg["tpc_extra_faces"].isArray())
        for (const auto& f : cfg["tpc_extra_faces"]) m_tpc_extra_faces.push_back(f.asInt());
    m_strength_cutoff = get(cfg, "strength_cutoff", m_strength_cutoff);
    m_sparse_lasso = get(cfg, "sparse_lasso", m_sparse_lasso);
    m_drift_speed = get(cfg, "drift_speed", m_drift_speed);
    m_trigger_offset = get(cfg, "trigger_offset", m_trigger_offset);
    // Optional per-input offsets (see QLMatching.h). Empty => scalar everywhere.
    m_trigger_offsets.clear();
    if (cfg.isMember("trigger_offsets") && cfg["trigger_offsets"].isArray())
        for (const auto& t : cfg["trigger_offsets"]) m_trigger_offsets.push_back(t.asDouble());
    // Optional per-input drift speeds (see QLMatching.h). Empty => scalar everywhere.
    m_drift_speeds.clear();
    if (cfg.isMember("drift_speeds") && cfg["drift_speeds"].isArray())
        for (const auto& s : cfg["drift_speeds"]) m_drift_speeds.push_back(s.asDouble());
    if (!m_drift_speeds.empty() && m_drift_speeds.size() != m_multiplicity) {
        raise<ValueError>("QLMatching: drift_speeds size %d != multiplicity %d",
                          (int) m_drift_speeds.size(), (int) m_multiplicity);
    }
    m_nchan = get(cfg, "nchan", m_nchan);

    // §A active-volume cushions (see match/docs/improve_progress.md). The raw
    // bounds come from m_dv->inner_bounds; these adjust the effective windows.
    // Defaults are the MicroBooNE-convention values (NOT bit-identical to the old
    // SBND literals); override to recover the old bounds.
    m_anode_ext1   = get(cfg, "anode_ext1",   m_anode_ext1);
    m_anode_ext2   = get(cfg, "anode_ext2",   m_anode_ext2);
    m_cathode_ext1 = get(cfg, "cathode_ext1", m_cathode_ext1);
    m_cathode_ext2 = get(cfg, "cathode_ext2", m_cathode_ext2);
    m_anode_ext1_margin = get(cfg, "anode_ext1_margin", m_anode_ext1_margin);
    m_y_cushion    = get(cfg, "y_cushion",    m_y_cushion);
    m_z_cushion    = get(cfg, "z_cushion",    m_z_cushion);
    m_two_boundary_margin = get(cfg, "two_boundary_margin", m_two_boundary_margin);
    m_robust_endpoint_trim  = get(cfg, "robust_endpoint_trim",  m_robust_endpoint_trim);
    m_robust_endpoint_frac  = get(cfg, "robust_endpoint_frac",  m_robust_endpoint_frac);
    m_robust_endpoint_count = get(cfg, "robust_endpoint_count", m_robust_endpoint_count);
    m_robust_endpoint_charge_frac = get(cfg, "robust_endpoint_charge_frac", m_robust_endpoint_charge_frac);
    m_robust_endpoint_charge_abs  = get(cfg, "robust_endpoint_charge_abs",  m_robust_endpoint_charge_abs);
    m_robust_endpoint_gap         = get(cfg, "robust_endpoint_gap",         m_robust_endpoint_gap);
    m_robust_endpoint_gap_charge_frac = get(cfg, "robust_endpoint_gap_charge_frac", m_robust_endpoint_gap_charge_frac);
    m_robust_endpoint_walk_to_floor   = get(cfg, "robust_endpoint_walk_to_floor",   m_robust_endpoint_walk_to_floor);
    m_robust_endpoint_gap_cathode     = get(cfg, "robust_endpoint_gap_cathode",     m_robust_endpoint_gap_cathode);

    m_mc_saturation_pe      = get(cfg, "mc_saturation_pe",      m_mc_saturation_pe);

    m_outbeam_ks_max      = get(cfg, "outbeam_ks_max",      m_outbeam_ks_max);
    m_outbeam_chi2ndf_max = get(cfg, "outbeam_chi2ndf_max", m_outbeam_chi2ndf_max);
    m_outbeam_pe_frac     = get(cfg, "outbeam_pe_frac",     m_outbeam_pe_frac);

    m_lasso_lambda     = get(cfg, "lasso_lambda",     m_lasso_lambda);
    m_delta_charge     = get(cfg, "delta_charge",     m_delta_charge);
    m_delta_light      = get(cfg, "delta_light",      m_delta_light);
    m_delta_shape      = get(cfg, "delta_shape",      m_delta_shape);
    m_bkg_weight       = get(cfg, "bkg_weight",       m_bkg_weight);
    m_pe_mismatch_knee = get(cfg, "pe_mismatch_knee", m_pe_mismatch_knee);
    m_pe_mismatch_floor = get(cfg, "pe_mismatch_floor", m_pe_mismatch_floor);
    m_lasso_flag_weight    = get(cfg, "lasso_flag_weight",    m_lasso_flag_weight);
    m_lasso_boundary_weight = get(cfg, "lasso_boundary_weight", m_lasso_boundary_weight);

    m_pe_err_floor      = get(cfg, "pe_err_floor",      m_pe_err_floor);
    m_pe_err_frac       = get(cfg, "pe_err_frac",       m_pe_err_frac);
    m_pe_err_knee       = get(cfg, "pe_err_knee",       m_pe_err_knee);
    m_pe_err_on_pred    = get(cfg, "pe_err_on_pred",    m_pe_err_on_pred);
    m_pe_err_lowpe_frac = get(cfg, "pe_err_lowpe_frac", m_pe_err_lowpe_frac);
    m_pe_err_lowpe_knee = get(cfg, "pe_err_lowpe_knee", m_pe_err_lowpe_knee);

    // Per-PD-family PE-error override (header doc). Absent/empty => byte-identical.
    if (cfg.isMember("pe_err_family_channels") && cfg["pe_err_family_channels"].isArray()) {
        m_pe_err_family_channels.clear();
        for (const auto& fam : cfg["pe_err_family_channels"]) {
            std::vector<int> chans;
            for (const auto& ch : fam) chans.push_back(ch.asInt());
            m_pe_err_family_channels.push_back(std::move(chans));
        }
        auto grab = [&cfg](const char* key, std::vector<double>& dst) {
            if (cfg.isMember(key) && cfg[key].isArray()) {
                dst.clear();
                for (const auto& v : cfg[key]) dst.push_back(v.asDouble());
            }
        };
        grab("pe_err_family_floor",      m_pe_err_family_floor);
        grab("pe_err_family_frac",       m_pe_err_family_frac);
        grab("pe_err_family_lowpe_frac", m_pe_err_family_lowpe_frac);
        grab("pe_err_family_lowpe_knee", m_pe_err_family_lowpe_knee);
    }
    if (!m_pe_err_family_channels.empty()) {
        auto resolve = [this](const std::vector<double>& fam_vals,
                              std::vector<double>& ch_vals) {
            ch_vals.assign(m_nchan, -1.0);
            for (std::size_t i = 0; i < m_pe_err_family_channels.size(); ++i) {
                if (i >= fam_vals.size() || fam_vals[i] < 0) continue;
                for (int ch : m_pe_err_family_channels[i]) {
                    if (ch >= 0 && ch < m_nchan) ch_vals[ch] = fam_vals[i];
                }
            }
        };
        resolve(m_pe_err_family_floor,      m_pe_err_ch_floor);
        resolve(m_pe_err_family_frac,       m_pe_err_ch_frac);
        resolve(m_pe_err_family_lowpe_frac, m_pe_err_ch_lowpe_frac);
        resolve(m_pe_err_family_lowpe_knee, m_pe_err_ch_lowpe_knee);
        int nover = 0;
        for (int ch = 0; ch < m_nchan; ++ch) {
            if (m_pe_err_ch_floor[ch] >= 0 || m_pe_err_ch_frac[ch] >= 0 ||
                m_pe_err_ch_lowpe_frac[ch] >= 0 || m_pe_err_ch_lowpe_knee[ch] >= 0) ++nover;
        }
        log->debug("pe_err_family override active: {} families, {} channels overridden",
                   m_pe_err_family_channels.size(), nover);
    }

    // Optional per-channel measured-PE gain correction (length nchan). Empty =>
    // identity (byte-identical). A length mismatch is a config error.
    if (cfg["measured_pe_scale"].isArray()) {
        m_measured_pe_scale.clear();
        for (const auto& v : cfg["measured_pe_scale"]) m_measured_pe_scale.push_back(v.asDouble());
        if (!m_measured_pe_scale.empty() && (int)m_measured_pe_scale.size() != m_nchan) {
            raise<ValueError>("QLMatching: measured_pe_scale length %d != nchan %d",
                              (int)m_measured_pe_scale.size(), m_nchan);
        }
    }
    m_flash_pe_threshold = get(cfg, "flash_pe_threshold", m_flash_pe_threshold);
    m_use_saturation_flag = get(cfg, "use_saturation_flag", m_use_saturation_flag);
    m_saturation_mask_fit = get(cfg, "saturation_mask_fit", m_saturation_mask_fit);
    m_use_coverage_flag = get(cfg, "use_coverage_flag", m_use_coverage_flag);
    m_coverage_min = get(cfg, "coverage_min", m_coverage_min);
    m_coverage_mask_fit = get(cfg, "coverage_mask_fit", m_coverage_mask_fit);
    m_pe_err_nodata = get(cfg, "pe_err_nodata", m_pe_err_nodata);

    m_bundle_ks_merge_max      = get(cfg, "bundle_ks_merge_max",      m_bundle_ks_merge_max);
    m_bundle_chi2ndf_merge_max = get(cfg, "bundle_chi2ndf_merge_max", m_bundle_chi2ndf_merge_max);
    m_bundle_addmerge_exponent = get(cfg, "bundle_addmerge_exponent", m_bundle_addmerge_exponent);
    m_highconsist_ks_max       = get(cfg, "highconsist_ks_max",       m_highconsist_ks_max);
    m_highconsist_min_ndf      = get(cfg, "highconsist_min_ndf",      m_highconsist_min_ndf);
    m_bundle_pe_ndf_knee       = get(cfg, "bundle_pe_ndf_knee",       m_bundle_pe_ndf_knee);
    m_bundle_mask_ks           = get(cfg, "bundle_mask_ks",           m_bundle_mask_ks);
    m_highconsist_ladder       = get(cfg, "highconsist_ladder",       m_highconsist_ladder);
    m_hc_clean_ks = get(cfg, "hc_clean_ks", m_hc_clean_ks);  m_hc_clean_c2 = get(cfg, "hc_clean_c2", m_hc_clean_c2);
    m_hc_good_ks  = get(cfg, "hc_good_ks",  m_hc_good_ks);   m_hc_good_c2  = get(cfg, "hc_good_c2",  m_hc_good_c2);
    m_hc_tb_ks    = get(cfg, "hc_tb_ks",    m_hc_tb_ks);     m_hc_tb_c2    = get(cfg, "hc_tb_c2",    m_hc_tb_c2);
    m_hc_miss_ks  = get(cfg, "hc_miss_ks",  m_hc_miss_ks);   m_hc_miss_c2  = get(cfg, "hc_miss_c2",  m_hc_miss_c2);
    m_hc_miss_min_ndf = get(cfg, "hc_miss_min_ndf", m_hc_miss_min_ndf);
    m_chi2_relax       = get(cfg, "chi2_relax",       m_chi2_relax);
    m_chi2_pmt_excess  = get(cfg, "chi2_pmt_excess",  m_chi2_pmt_excess);
    m_chi2_pmt_ratio   = get(cfg, "chi2_pmt_ratio",   m_chi2_pmt_ratio);
    m_chi2_pmt_inflate = get(cfg, "chi2_pmt_inflate", m_chi2_pmt_inflate);
    m_chi2_sat_inflate = get(cfg, "chi2_sat_inflate", m_chi2_sat_inflate);

    m_readout_window_ticks = get(cfg, "readout_window_ticks", m_readout_window_ticks);
    m_window_edge_ticks    = get(cfg, "window_edge_ticks",    m_window_edge_ticks);

    m_require_containment  = get(cfg, "require_containment",  m_require_containment);

    m_cross_side_filter    = get(cfg, "cross_side_filter",    m_cross_side_filter);
    m_crossside_skip_vis   = get(cfg, "crossside_skip_vis",   m_crossside_skip_vis);

    m_vis_sample_stride    = get(cfg, "vis_sample_stride",    m_vis_sample_stride);
    m_vis_sample_min_pts   = get(cfg, "vis_sample_min_pts",   m_vis_sample_min_pts);

    m_reject_overpred      = get(cfg, "reject_overpred",      m_reject_overpred);
    m_overpred_total_ratio = get(cfg, "overpred_total_ratio", m_overpred_total_ratio);
    m_overpred_maxch_ratio = get(cfg, "overpred_maxch_ratio", m_overpred_maxch_ratio);
    if (cfg.isMember("overpred_channels") && cfg["overpred_channels"].isArray()) {
        m_overpred_channels.clear();
        for (const auto& jch : cfg["overpred_channels"]) m_overpred_channels.push_back(jch.asInt());
    }

    m_empty_rescue          = get(cfg, "empty_rescue",          m_empty_rescue);
    m_opflash_phys_gid      = get(cfg, "opflash_phys_gid",      m_opflash_phys_gid);
    m_rescue_metric_max     = get(cfg, "rescue_metric_max",     m_rescue_metric_max);
    m_rescue_exponent       = get(cfg, "rescue_exponent",       m_rescue_exponent);
    m_rescue_boundary_weight = get(cfg, "rescue_boundary_weight", m_rescue_boundary_weight);

    m_cluster_rescue            = get(cfg, "cluster_rescue",            m_cluster_rescue);
    m_cluster_rescue_ks_max     = get(cfg, "cluster_rescue_ks_max",     m_cluster_rescue_ks_max);
    m_cluster_rescue_chi2ndf_max = get(cfg, "cluster_rescue_chi2ndf_max", m_cluster_rescue_chi2ndf_max);
    m_cluster_rescue_ratio_lo   = get(cfg, "cluster_rescue_ratio_lo",   m_cluster_rescue_ratio_lo);
    m_cluster_rescue_ratio_hi   = get(cfg, "cluster_rescue_ratio_hi",   m_cluster_rescue_ratio_hi);
    m_cluster_rescue_precull    = get(cfg, "cluster_rescue_precull",    m_cluster_rescue_precull);
    m_cluster_rescue_precull_additive = get(cfg, "cluster_rescue_precull_additive", m_cluster_rescue_precull_additive);
    m_cluster_rescue_relaxed            = get(cfg, "cluster_rescue_relaxed",            m_cluster_rescue_relaxed);
    m_cluster_rescue_relaxed_ks_max     = get(cfg, "cluster_rescue_relaxed_ks_max",     m_cluster_rescue_relaxed_ks_max);
    m_cluster_rescue_relaxed_chi2ndf_max = get(cfg, "cluster_rescue_relaxed_chi2ndf_max", m_cluster_rescue_relaxed_chi2ndf_max);
    m_cluster_rescue_relaxed_ratio_lo   = get(cfg, "cluster_rescue_relaxed_ratio_lo",   m_cluster_rescue_relaxed_ratio_lo);
    m_cluster_rescue_relaxed_ratio_hi   = get(cfg, "cluster_rescue_relaxed_ratio_hi",   m_cluster_rescue_relaxed_ratio_hi);
    m_cluster_rescue_relaxed_min_length = get(cfg, "cluster_rescue_relaxed_min_length", m_cluster_rescue_relaxed_min_length);
    m_cluster_rescue_sat_ratio_relax = get(cfg, "cluster_rescue_sat_ratio_relax", m_cluster_rescue_sat_ratio_relax);
    m_cluster_rescue_sat_frac_min    = get(cfg, "cluster_rescue_sat_frac_min",    m_cluster_rescue_sat_frac_min);
    m_cluster_rescue_sat_ratio_mult  = get(cfg, "cluster_rescue_sat_ratio_mult",  m_cluster_rescue_sat_ratio_mult);
    if (m_shared_flash && (m_empty_rescue || m_cluster_rescue)) {
        log->warn("QLMatching: empty_rescue/cluster_rescue are not shared_flash-aware "
                  "and will be SKIPPED in the joint fit");
    }
    // Shared-flash-aware rescue variants (doc 19 phase 5); default OFF.
    // Reuse the per-run threshold knobs (rescue_metric_max/exponent/
    // boundary_weight; cluster_rescue_*).
    m_empty_rescue_shared   = get(cfg, "empty_rescue_shared",   m_empty_rescue_shared);
    m_cluster_rescue_shared = get(cfg, "cluster_rescue_shared", m_cluster_rescue_shared);
    if (!m_shared_flash && (m_empty_rescue_shared || m_cluster_rescue_shared)) {
        log->warn("QLMatching: empty_rescue_shared/cluster_rescue_shared need "
                  "shared_flash and will be inert");
    }

    // Optional CPA structure-exclusion fiducial (SBND). Empty => disabled, and the
    // cathode-end flag_at_x_boundary keeps the original flat-cathode 1-D test.
    const auto cathode_fv_tn = get<std::string>(cfg, "cathode_fiducial", std::string(""));
    if (!cathode_fv_tn.empty()) {
        m_cathode_fv = Factory::find_tn<IFiducial>(cathode_fv_tn);
    }
    if (m_xtpc_cathode_tol > 0.0 && m_cross_side_filter) {
        // Not composed: a provisional (uncontained) bundle lacks flag_at_x_boundary,
        // so the opaque-cathode cross-side drop would cull it before pairing. No
        // current detector enables both (cross_side_filter is SBND, whose CPA
        // fiducial also excludes it from the rescue's flat-window path).
        log->warn("QLMatching: xtpc_cathode_tol with cross_side_filter is untested; "
                  "cross-side-lit provisional bundles are dropped before pairing");
    }

    if (cfg["VUVEfficiency"].isArray()) {
        m_VUVEfficiency.clear();
        for (const auto& v : cfg["VUVEfficiency"]) m_VUVEfficiency.push_back(v.asDouble());
    }
    if (cfg["VISEfficiency"].isArray()) {
        m_VISEfficiency.clear();
        for (const auto& v : cfg["VISEfficiency"]) m_VISEfficiency.push_back(v.asDouble());
    }

    // Per-PMT non-linearity correction on the predicted PE (default OFF / identity).
    m_pmt_nonlinearity = get(cfg, "pmt_nonlinearity", m_pmt_nonlinearity);
    m_pmt_nl_knee      = get(cfg, "pmt_nl_knee",      m_pmt_nl_knee);
    if (cfg["pmt_nl_beta"].isArray()) {
        m_pmt_nl_beta.clear();
        for (const auto& v : cfg["pmt_nl_beta"]) m_pmt_nl_beta.push_back(v.asDouble());
    }
    if (cfg["pmt_nl_gamma"].isArray()) {
        m_pmt_nl_gamma.clear();
        for (const auto& v : cfg["pmt_nl_gamma"]) m_pmt_nl_gamma.push_back(v.asDouble());
    }

    // Load SBND semi-analytical optical model from JSON.
    auto top = Persist::load(m_semimodel_file);
    if (!top.isObject()) {
        raise<ValueError>("QLMatching: invalid semimodel_file '%s'", m_semimodel_file);
    }
    const auto vuv_cfg = top["VUVHits"];
    const auto vis_cfg = top["VISHits"];
    const auto geom_cfg = top["Geometry"];
    const auto opdets_cfg = top["OpDets"];
    // VISHits is only needed when reflected light is computed (m_doReflectedLight);
    // detectors without a reflected component (PDHD) may omit it or pass empty {}.
    if (!vuv_cfg.isObject() || !geom_cfg.isObject() || !opdets_cfg.isArray()
        || (m_doReflectedLight && !vis_cfg.isObject())) {
        raise<ValueError>("QLMatching: semimodel_file '%s' missing required sections", m_semimodel_file);
    }

    SemiAnalyticalModel::Geometry geom;
    geom.active_center_y       = get<double>(geom_cfg, "active_center_y", 0.0);
    geom.active_center_z       = get<double>(geom_cfg, "active_center_z", 0.0);
    geom.active_size_y         = get<double>(geom_cfg, "active_size_y",   0.0);
    geom.active_size_z         = get<double>(geom_cfg, "active_size_z",   0.0);
    geom.cathode_x             = get<double>(geom_cfg, "cathode_x",       0.0);
    geom.vuv_absorption_length = get<double>(geom_cfg, "vuv_absorption_length", 85.0);
    m_cathode_x = geom.cathode_x;

    m_opdets.clear();
    m_opdets.reserve(opdets_cfg.size());
    for (const auto& od : opdets_cfg) {
        SemiAnalyticalModel::OpticalDetector o;
        o.h           = get<double>(od, "h", -1.0);
        o.w           = get<double>(od, "w", -1.0);
        o.center      = WireCell::Point(od["x"].asDouble(), od["y"].asDouble(), od["z"].asDouble());
        o.type        = od.get("type", 1).asInt();
        o.orientation = od.get("orientation", 0).asInt();
        m_opdets.push_back(o);
    }

    // ch_mask entries index m_opdets/run.opdet_mask directly (build_opdet_mask,
    // below); an out-of-range value from a config typo would otherwise be a
    // silent out-of-bounds write there. All current detector configs pass
    // in-range indices, so this never fires today -- it converts a future
    // typo from memory corruption into a clear config error.
    for (int c : m_ch_mask) {
        if (c < 0 || (std::size_t) c >= m_opdets.size()) {
            raise<ValueError>("QLMatching: ch_mask entry %d out of range [0,%d)",
                              c, (int) m_opdets.size());
        }
    }

    // VUVEfficiency/VISEfficiency are indexed by the same OpDet index as
    // m_opdets (pred_flash accumulation, execute()). A too-short array would
    // otherwise throw deep in the per-point hot loop (.at(idet)) instead of
    // here at configure time; a too-long one is harmless (the extra tail is
    // simply unused), so this checks "too short", not "wrong length" --
    // every current detector config supplies an array of exactly nopdets
    // length (SBND via the 312-entry default, PDHD/PDVD via an explicit
    // nchan-length override), so this never fires today.
    if (m_VUVEfficiency.size() < m_opdets.size()) {
        raise<ValueError>("QLMatching: VUVEfficiency length %d < nopdets %d",
                          (int) m_VUVEfficiency.size(), (int) m_opdets.size());
    }
    if (m_VISEfficiency.size() < m_opdets.size()) {
        raise<ValueError>("QLMatching: VISEfficiency length %d < nopdets %d",
                          (int) m_VISEfficiency.size(), (int) m_opdets.size());
    }

    m_semi_model = std::make_unique<SemiAnalyticalModel>(
        vuv_cfg, vis_cfg, geom, m_opdets, m_doReflectedLight);

    // Optional gridded-library visibility backend (default "semi" leaves the
    // path above untouched).  The semi model and OpDet table stay loaded --
    // masks, cathode-side and boundary-flag logic use them either way.
    m_lib_model.reset();
    if (m_light_model == "library") {
        if (m_photon_library_file.empty()) {
            raise<ValueError>("QLMatching: light_model 'library' needs photon_library_file");
        }
        m_lib_model = std::make_unique<PhotonLibraryModel>(m_photon_library_file);
        if (m_lib_model->nchan() != m_opdets.size()) {
            raise<ValueError>("QLMatching: photon library nchan %d != nopdets %d",
                              (int) m_lib_model->nchan(), (int) m_opdets.size());
        }
        // The bare nchan count check above cannot catch a library whose
        // channel order does not match m_opdets' order (e.g. a future
        // re-export using a different channel convention) -- it only
        // catches a length mismatch. When the library's optional
        // chan_pos_cm is present, cross-check it against m_opdets[i].center
        // and raise on a real mismatch rather than silently mispredicting.
        // Absent in a library file (as all files shipped before this check
        // was added), has_positions() is false and this is a no-op -- adding
        // the check does not require regenerating existing library files.
        if (m_lib_model->has_positions()) {
            for (std::size_t i = 0; i < m_opdets.size(); ++i) {
                const double d = (m_lib_model->position(i) - m_opdets[i].center).magnitude();
                if (d > m_photon_library_pos_tol) {
                    raise<ValueError>(
                        "QLMatching: photon library channel %d position mismatch with OpDet "
                        "table: library=(%.1f,%.1f,%.1f) opdet=(%.1f,%.1f,%.1f) d=%.1f cm > "
                        "tol %.1f cm -- library channel order likely does not match m_opdets",
                        (int) i,
                        m_lib_model->position(i).x(), m_lib_model->position(i).y(), m_lib_model->position(i).z(),
                        m_opdets[i].center.x(), m_opdets[i].center.y(), m_opdets[i].center.z(),
                        d, m_photon_library_pos_tol);
                }
            }
        }
    }
    else if (m_light_model != "semi") {
        raise<ValueError>("QLMatching: unknown light_model '%s'", m_light_model);
    }

    log->debug("QLMatching configured: nopdets={}, semimodel_file={}, light_model={}",
               m_opdets.size(), m_semimodel_file, m_light_model);
    // Sentinel: only when the anode floor is widened off its prototype default, so a
    // smoke run can prove the knob reached the component (silent => legacy path).
    if (m_anode_ext1_margin != 1.0 * units::cm) {
        log->debug("QLMatching anode_ext1_margin={} cm => anode containment/flag floor {} cm",
                   m_anode_ext1_margin / units::cm,
                   (m_anode_ext1 - m_anode_ext1_margin) / units::cm);
    }
    // Sentinel: only when readout-uncovered channels are kept in the fit, so a
    // smoke run can prove the knob reached the component (silent => legacy drop).
    if (m_use_coverage_flag && !m_coverage_mask_fit) {
        log->debug("QLMatching coverage_mask_fit=false => uncovered channels stay in the "
                   "fit at measured 0 (pe_err_nodata={} PE, coverage_min={})",
                   m_pe_err_nodata, m_coverage_min);
    }
    // Sentinel: only when the channel-scoped flash admission is active, so a smoke
    // run can prove the knob reached the component (silent => legacy admission).
    if (!m_flash_sel_channels.empty()) {
        log->debug("QLMatching flash_sel: admission over {} channels requires sum PE >= {} "
                   "and >= {} channels fired at >= {} PE",
                   m_flash_sel_channels.size(), m_flash_sel_minPE,
                   m_flash_sel_min_fired, m_flash_sel_fired_pe);
    }
    // Per-channel flag form of overpred_channels (hoisted out of the per-bundle
    // prefilter loop; empty <=> the scope knob is off).
    m_overpred_ch_sel.clear();
    if (!m_overpred_channels.empty()) {
        m_overpred_ch_sel.assign(std::max(m_nchan, 1), 0);
        for (int ch : m_overpred_channels)
            if (ch >= 0 && ch < (int) m_overpred_ch_sel.size()) m_overpred_ch_sel[ch] = 1;
    }
    // Sentinel: only when the overpred prefilter is channel-scoped (silent =>
    // legacy full-masked-set ratios).
    if (m_reject_overpred && !m_overpred_channels.empty()) {
        log->debug("QLMatching reject_overpred scoped to {} channels "
                   "(R_total <= {}, R_max <= {})",
                   m_overpred_channels.size(), m_overpred_total_ratio, m_overpred_maxch_ratio);
    }
    // Same sentinel for the rail-flag path (silent => legacy drop from chi2/KS).
    if (m_use_saturation_flag && !m_saturation_mask_fit) {
        log->debug("QLMatching saturation_mask_fit=false => railed channels stay in the "
                   "chi2/KS at their clipped PE (chi2_sat_inflate={}); LASSO rows stay zeroed",
                   m_chi2_sat_inflate);
    }
}

WireCell::Configuration QLMatching::default_configuration() const
{
    Configuration cfg;
    cfg["inpath"]          = m_inpath;
    cfg["outpath"]         = m_outpath;
    cfg["calib_dump"]      = m_calib_dump;
    cfg["flash_group_window"] = m_flash_group_window;
    cfg["cathode_diag"]    = m_cathode_diag;
    cfg["cathode_diag_radius"] = m_cathode_diag_radius;
    cfg["xtpc_flag"]         = m_xtpc_flag;
    cfg["xtpc_dmax"]         = m_xtpc_dmax;
    cfg["xtpc_dmax2"]        = m_xtpc_dmax2;
    cfg["xtpc_angle_max"]    = m_xtpc_angle_max;
    cfg["xtpc_hough_radius"] = m_xtpc_hough_radius;
    cfg["xtpc_joint_pin"]    = m_xtpc_joint_pin;
    cfg["xtpc_pin_angle"]    = m_xtpc_pin_angle;
    cfg["xtpc_cathode_tol"]   = m_xtpc_cathode_tol;
    cfg["xtpc_cathode_qfrac"] = m_xtpc_cathode_qfrac;
    // quality gates (doc 19); legacy defaults => inert
    cfg["xtpc_pin_min_strength"] = m_xtpc_pin_min_strength;
    cfg["xtpc_sc1_light_gate"]   = m_xtpc_sc1_light_gate;
    cfg["xtpc_sc1_ks_max"]       = m_xtpc_sc1_ks_max;
    cfg["xtpc_sc1_c2n_max"]      = m_xtpc_sc1_c2n_max;
    cfg["xtpc_cathode_ks_max"]   = m_xtpc_cathode_ks_max;
    cfg["postcull_unflagged"]    = m_postcull_unflagged;
    cfg["postcull_ks_max"]       = m_postcull_ks_max;
    cfg["postcull_c2n_max"]      = m_postcull_c2n_max;
    cfg["postcull_before_rescue"] = m_postcull_before_rescue;
    cfg["postcull_wtrunc_overpred"] = m_postcull_wtrunc_overpred;
    cfg["postcull_wtrunc_ratio_hi"] = m_postcull_wtrunc_ratio_hi;
    cfg["postcull_wtrunc_sat_frac"] = m_postcull_wtrunc_sat_frac;
    cfg["postcull_pin_overpred"]    = m_postcull_pin_overpred;
    cfg["postcull_pin_ratio_hi"]    = m_postcull_pin_ratio_hi;
    cfg["nchan"]           = m_nchan;
    cfg["semimodel_file"]  = m_semimodel_file;
    cfg["light_model"]     = m_light_model;
    cfg["photon_library_file"] = m_photon_library_file;
    cfg["pmts"]            = m_pmts;
    cfg["active_opdet_types"] = Json::arrayValue;
    for (int t : m_active_opdet_types) cfg["active_opdet_types"].append(t);
    cfg["data"]            = m_data;
    cfg["beamonly"]        = m_beamonly;
    cfg["ch_mask"]         = Json::arrayValue;
    cfg["auto_mask"]              = m_auto_mask;
    cfg["auto_mask_pe_low"]       = m_auto_mask_pe_low;
    cfg["auto_mask_neighbors"]    = m_auto_mask_neighbors;
    cfg["auto_mask_pe_bright"]    = m_auto_mask_pe_bright;
    cfg["auto_mask_min_contrast"] = m_auto_mask_min_contrast;
    cfg["auto_mask_min_flash"]    = m_auto_mask_min_flash;
    cfg["auto_mask_same_type"]    = m_auto_mask_same_type;
    cfg["shared_flash"]           = m_shared_flash;
    cfg["opdet_all_volumes"]      = m_opdet_all_volumes;
    cfg["vd_surface_flags"]       = m_vd_surface_flags;
    cfg["pd_wall_cushion"]        = m_pd_wall_cushion;
    cfg["pd_wall_channels_ylo"]   = Json::arrayValue;
    cfg["pd_wall_channels_yhi"]   = Json::arrayValue;
    cfg["anode_pd_channels"]      = Json::arrayValue;
    cfg["flash_minPE"]     = m_flash_minPE;
    cfg["flash_sel_channels"]  = Json::arrayValue;
    cfg["flash_sel_minPE"]     = m_flash_sel_minPE;
    cfg["flash_sel_min_fired"] = m_flash_sel_min_fired;
    cfg["flash_sel_fired_pe"]  = m_flash_sel_fired_pe;
    cfg["flash_mintime"]   = m_flash_mintime;
    cfg["flash_maxtime"]   = m_flash_maxtime;
    cfg["beam_mintime"]    = m_beam_mintime;
    cfg["beam_maxtime"]    = m_beam_maxtime;
    cfg["QtoL"]            = m_QtoL;
    cfg["doReflectedLight"] = m_doReflectedLight;
    cfg["tpc_face"]        = m_tpc_face;
    cfg["strength_cutoff"] = m_strength_cutoff;
    cfg["sparse_lasso"]    = m_sparse_lasso;
    cfg["drift_speed"]     = m_drift_speed;
    cfg["trigger_offset"]  = m_trigger_offset;
    cfg["trigger_offsets"] = Json::Value(Json::arrayValue);
    cfg["drift_speeds"]    = Json::Value(Json::arrayValue);

    cfg["anode_ext1"]   = m_anode_ext1;
    cfg["anode_ext2"]   = m_anode_ext2;
    cfg["cathode_ext1"] = m_cathode_ext1;
    cfg["cathode_ext2"] = m_cathode_ext2;
    cfg["anode_ext1_margin"] = m_anode_ext1_margin;
    cfg["y_cushion"]    = m_y_cushion;
    cfg["z_cushion"]    = m_z_cushion;
    cfg["two_boundary_margin"] = m_two_boundary_margin;
    cfg["robust_endpoint_trim"]  = m_robust_endpoint_trim;
    cfg["robust_endpoint_frac"]  = m_robust_endpoint_frac;
    cfg["robust_endpoint_count"] = m_robust_endpoint_count;
    cfg["robust_endpoint_charge_frac"] = m_robust_endpoint_charge_frac;
    cfg["robust_endpoint_charge_abs"]  = m_robust_endpoint_charge_abs;
    cfg["robust_endpoint_gap"]         = m_robust_endpoint_gap;
    cfg["robust_endpoint_gap_charge_frac"] = m_robust_endpoint_gap_charge_frac;
    cfg["robust_endpoint_walk_to_floor"]   = m_robust_endpoint_walk_to_floor;
    cfg["robust_endpoint_gap_cathode"]     = m_robust_endpoint_gap_cathode;

    cfg["mc_saturation_pe"]      = m_mc_saturation_pe;

    cfg["outbeam_ks_max"]      = m_outbeam_ks_max;
    cfg["outbeam_chi2ndf_max"] = m_outbeam_chi2ndf_max;
    cfg["outbeam_pe_frac"]     = m_outbeam_pe_frac;

    cfg["lasso_lambda"]      = m_lasso_lambda;
    cfg["delta_charge"]      = m_delta_charge;
    cfg["delta_light"]       = m_delta_light;
    cfg["delta_shape"]       = m_delta_shape;
    cfg["bkg_weight"]        = m_bkg_weight;
    cfg["pe_mismatch_knee"]  = m_pe_mismatch_knee;
    cfg["pe_mismatch_floor"] = m_pe_mismatch_floor;
    cfg["lasso_flag_weight"]    = m_lasso_flag_weight;
    cfg["lasso_boundary_weight"] = m_lasso_boundary_weight;

    cfg["pe_err_floor"]       = m_pe_err_floor;
    cfg["pe_err_frac"]        = m_pe_err_frac;
    cfg["pe_err_knee"]        = m_pe_err_knee;
    cfg["pe_err_on_pred"]     = m_pe_err_on_pred;
    cfg["pe_err_lowpe_frac"]  = m_pe_err_lowpe_frac;
    cfg["pe_err_lowpe_knee"]  = m_pe_err_lowpe_knee;
    // empty => no per-family override (byte-identical legacy error model)
    cfg["pe_err_family_channels"]   = Json::Value(Json::arrayValue);
    cfg["pe_err_family_floor"]      = Json::Value(Json::arrayValue);
    cfg["pe_err_family_frac"]       = Json::Value(Json::arrayValue);
    cfg["pe_err_family_lowpe_frac"] = Json::Value(Json::arrayValue);
    cfg["pe_err_family_lowpe_knee"] = Json::Value(Json::arrayValue);
    cfg["measured_pe_scale"]  = Json::Value(Json::arrayValue);  // empty => identity
    cfg["flash_pe_threshold"] = m_flash_pe_threshold;
    cfg["use_saturation_flag"] = m_use_saturation_flag;
    cfg["saturation_mask_fit"] = m_saturation_mask_fit;
    cfg["use_coverage_flag"] = m_use_coverage_flag;
    cfg["coverage_min"] = m_coverage_min;
    cfg["coverage_mask_fit"] = m_coverage_mask_fit;
    cfg["pe_err_nodata"] = m_pe_err_nodata;

    cfg["bundle_ks_merge_max"]      = m_bundle_ks_merge_max;
    cfg["bundle_chi2ndf_merge_max"] = m_bundle_chi2ndf_merge_max;
    cfg["bundle_addmerge_exponent"] = m_bundle_addmerge_exponent;
    cfg["highconsist_ks_max"]       = m_highconsist_ks_max;
    cfg["highconsist_min_ndf"]      = m_highconsist_min_ndf;
    cfg["bundle_pe_ndf_knee"]       = m_bundle_pe_ndf_knee;
    cfg["bundle_mask_ks"]           = m_bundle_mask_ks;
    cfg["highconsist_ladder"]       = m_highconsist_ladder;
    cfg["hc_clean_ks"] = m_hc_clean_ks;  cfg["hc_clean_c2"] = m_hc_clean_c2;
    cfg["hc_good_ks"]  = m_hc_good_ks;   cfg["hc_good_c2"]  = m_hc_good_c2;
    cfg["hc_tb_ks"]    = m_hc_tb_ks;     cfg["hc_tb_c2"]    = m_hc_tb_c2;
    cfg["hc_miss_ks"]  = m_hc_miss_ks;   cfg["hc_miss_c2"]  = m_hc_miss_c2;
    cfg["hc_miss_min_ndf"] = m_hc_miss_min_ndf;
    cfg["chi2_relax"]       = m_chi2_relax;
    cfg["chi2_pmt_excess"]  = m_chi2_pmt_excess;
    cfg["chi2_pmt_ratio"]   = m_chi2_pmt_ratio;
    cfg["chi2_pmt_inflate"] = m_chi2_pmt_inflate;
    cfg["chi2_sat_inflate"] = m_chi2_sat_inflate;

    cfg["readout_window_ticks"] = m_readout_window_ticks;
    cfg["window_edge_ticks"]    = m_window_edge_ticks;
    cfg["require_containment"]  = m_require_containment;
    cfg["cross_side_filter"]    = m_cross_side_filter;
    cfg["crossside_skip_vis"]   = m_crossside_skip_vis;
    cfg["vis_sample_stride"]    = m_vis_sample_stride;
    cfg["vis_sample_min_pts"]   = m_vis_sample_min_pts;
    cfg["reject_overpred"]      = m_reject_overpred;
    cfg["overpred_total_ratio"] = m_overpred_total_ratio;
    cfg["overpred_channels"]    = Json::arrayValue;
    cfg["empty_rescue"]          = m_empty_rescue;
    cfg["opflash_phys_gid"]      = m_opflash_phys_gid;
    cfg["rescue_metric_max"]     = m_rescue_metric_max;
    cfg["rescue_exponent"]       = m_rescue_exponent;
    cfg["rescue_boundary_weight"] = m_rescue_boundary_weight;
    cfg["cluster_rescue"]            = m_cluster_rescue;
    cfg["empty_rescue_shared"]       = m_empty_rescue_shared;
    cfg["cluster_rescue_shared"]     = m_cluster_rescue_shared;
    cfg["cluster_rescue_ks_max"]     = m_cluster_rescue_ks_max;
    cfg["cluster_rescue_chi2ndf_max"] = m_cluster_rescue_chi2ndf_max;
    cfg["cluster_rescue_ratio_lo"]   = m_cluster_rescue_ratio_lo;
    cfg["cluster_rescue_ratio_hi"]   = m_cluster_rescue_ratio_hi;
    cfg["cluster_rescue_precull"]    = m_cluster_rescue_precull;
    cfg["cluster_rescue_precull_additive"] = m_cluster_rescue_precull_additive;
    cfg["cluster_rescue_relaxed"]            = m_cluster_rescue_relaxed;
    cfg["cluster_rescue_relaxed_ks_max"]     = m_cluster_rescue_relaxed_ks_max;
    cfg["cluster_rescue_relaxed_chi2ndf_max"] = m_cluster_rescue_relaxed_chi2ndf_max;
    cfg["cluster_rescue_relaxed_ratio_lo"]   = m_cluster_rescue_relaxed_ratio_lo;
    cfg["cluster_rescue_relaxed_ratio_hi"]   = m_cluster_rescue_relaxed_ratio_hi;
    cfg["cluster_rescue_relaxed_min_length"] = m_cluster_rescue_relaxed_min_length;
    cfg["cluster_rescue_sat_ratio_relax"] = m_cluster_rescue_sat_ratio_relax;
    cfg["cluster_rescue_sat_frac_min"]    = m_cluster_rescue_sat_frac_min;
    cfg["cluster_rescue_sat_ratio_mult"]  = m_cluster_rescue_sat_ratio_mult;
    cfg["overpred_maxch_ratio"] = m_overpred_maxch_ratio;
    cfg["cathode_fiducial"]     = "";
    cfg["pmt_nonlinearity"]     = m_pmt_nonlinearity;
    cfg["pmt_nl_knee"]          = m_pmt_nl_knee;
    cfg["pmt_nl_beta"]          = Json::arrayValue;
    cfg["pmt_nl_gamma"]         = Json::arrayValue;
    return cfg;
}

std::vector<std::string> QLMatching::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    return std::vector<std::string>(m_multiplicity, tname);
}

bool QLMatching::operator()(const input_vector& invec, output_pointer& out)
{
    out = nullptr;
    using WireCell::Clus::Facade::float_t;

    if (invec.size() != m_multiplicity) {
        raise<ValueError>("QLMatching: unexpected multiplicity got %d want %d",
                          (int)invec.size(), (int)m_multiplicity);
    }
    // EOS: succeed when ALL inputs are EOS (nullptr); a partial EOS is an error.
    std::size_t neos = 0;
    for (const auto& in : invec) if (!in) ++neos;
    if (neos == invec.size()) {
        log->debug("EOS at call {}", m_count++);
        return true;
    }
    if (neos) raise<ValueError>("QLMatching: missing %d input tensors", (int)neos);

    ExecMon em("starting QLMatching");
    // Baseline RSS at QLMatching entry (process is already loaded with the clustered
    // pctree, so this is "everything before matching"). Used below to report the
    // matcher's own incremental footprint vs the whole-process resident set.
    const auto mu0 = em.mu.current();

    m_extreme_cache.clear();  // cluster facades are per-event; drop stale endpoints
    m_pca_endpoints_cache.clear();
    m_wall_flag_cache.clear();
    m_endpoint_slice_cache.clear();

    // Match each APA's input independently in its own fresh, isolated ApaRun (one
    // per input port). The per-APA masks, geometry and flash_gid stride are
    // unchanged, and the isolation keeps the pointer-keyed map iteration
    // deterministic across APAs. With multiplicity 1 this is exactly the
    // historical single-input matcher.
    const auto t_prefit0 = wallclock::now();
    std::vector<ApaRun> runs(invec.size());
    for (std::size_t k = 0; k < invec.size(); ++k) {
        ApaRun& run = runs[k];
        run.input_idx    = k;
        run.anode        = m_anodes.at(k);
        run.charge_ident = invec[k]->ident();
        run.inpath       = m_inpath;
        if (run.inpath.find("%") != std::string::npos)
            run.inpath = String::format(run.inpath, run.charge_ident);

        const auto& charge_tens = *invec[k]->tensors();
        log->debug("charge_tens.size {}", charge_tens.size());
        auto root_live = Aux::TensorDM::as_pctree(charge_tens, run.inpath + "/live");
        if (!root_live) {
            log->error("Failed to get point cloud tree from \"{}\"", run.inpath);
            return false;
        }
        log->debug("Got live pctree with {} children", root_live->nchildren());
        log->debug(em("got live pctree"));
        run.root_live = std::move(root_live);
        run.grouping  = run.root_live->value.facade<Grouping>();

        run_one_apa_prefit(run);

        if (!run.flashes.empty()) {
            log->debug("total_charge_blob {} total_charge_point {} total_charge_blob_all {}",
                       run.total_charge_blob / run.flashes.size(),
                       run.total_charge_point / run.flashes.size(),
                       run.total_charge_blob_all);
        }
        else {
            log->debug("total_charge_blob {} total_charge_point {} total_charge_blob_all {}",
                       0, 0, run.total_charge_blob_all);
        }
    }

    // Cross-TPC cathode-crossing FLAG pass (needs both TPCs' candidate bundles, so it
    // runs between the per-APA prefit above and the per-APA cull+fit below). It marks
    // geometrically cross-TPC-consistent candidate pairs (flag_xtpc_consistent, and
    // flag_xtpc_scenario1 for the tight cathode-end case) on the FULL pre-cull bundle
    // set — BEFORE cull_inconsistent collapses each cluster onto its single
    // high-consistent flash, which would otherwise drop a truncated crosser half before
    // it can be paired. Off by default (single-APA / xtpc_flag:false) => not called.
    const double t_prefit = ms_since(t_prefit0);
    const auto t_cull0 = wallclock::now();
    if (m_xtpc_flag && runs.size() > 1) cull_cross_tpc(runs);
    // xtpc cathode rescue resolution (xtpc_cathode_tol > 0 only; otherwise no
    // provisional flag exists and this is a no-op): a provisional cathode-overshoot
    // bundle survives ONLY if cull_cross_tpc confirmed it as a scenario-1 crosser
    // half (opposite-volume partner on the SAME flash, clouds meeting within
    // xtpc_dmax). Unconfirmed provisionals are removed here, BEFORE
    // cull_inconsistent/fit, so downstream sees exactly the legacy candidate set.
    // Note if xtpc_flag is off or there is a single run, cull_cross_tpc never ran
    // and EVERY provisional is purged => legacy again.
    if (m_xtpc_cathode_tol > 0.0)
        for (auto& run : runs) purge_unconfirmed_cathode_rescue(run);
    // Per-TPC consistency cull, now AFTER the cross-TPC flag pass so it can honour the
    // xtpc flags (xtpc-priority). When no xtpc flag is set (non-SBND / xtpc_flag:false)
    // this reproduces the historical cull_inconsistent exactly => bit-identical.
    for (std::size_t k = 0; k < invec.size(); ++k) cull_inconsistent(runs[k]);
    const double t_cull = ms_since(t_cull0);

    const auto t_fit0 = wallclock::now();
    if (m_shared_flash && runs.size() > 1) {
        // Shared-flash JOINT fit (PDVD): one LASSO system over the physical flashes
        // with columns from every port, so both drift sides share each flash's
        // explanation. The per-run output stages then run unchanged per port.
        auto t = wallclock::now();
        fit_round1_shared(runs);  const double t_r1 = ms_since(t); t = wallclock::now();
        fit_round2_shared(runs);  const double t_r2 = ms_since(t); t = wallclock::now();
        for (auto& run : runs) {
            apply_matched_t0s(run);
            write_opflash_pc(run);
        }
        log->debug("QLtiming fit(shared): fit_round1 {:.1f} fit_round2 {:.1f} out {:.1f} ms",
                   t_r1, t_r2, ms_since(t));
    }
    else {
        for (std::size_t k = 0; k < invec.size(); ++k) run_one_apa_fit(runs[k]);
    }
    const double t_fit = ms_since(t_fit0);

    // Hand-scan calibration / cathode diagnostic dumps read the per-APA runs while
    // the matching decomposition is still materialized (main + associated are
    // separate Facade clusters), so run them BEFORE recompose restores the merged
    // grouping. Read-only; never perturb the matching. (Moved up from after output
    // assembly; the matched tensor output is independent of them.)
    if (!m_calib_dump.empty()) dump_calib(runs);
    if (!m_cathode_diag.empty()) dump_cathode_diag(runs);

    // Restore the stage-3 cluster grouping (undo decompose's Grouping::separate)
    // so QLMatching emits a T0/flash-annotated copy of its INPUT tree rather than
    // the matching-internal main/associated split. See recompose_cluster_groups().
    for (std::size_t k = 0; k < invec.size(); ++k) recompose_cluster_groups(runs[k]);

    // ---- Build / merge outputs ----
    const auto t_out0 = wallclock::now();
    const int out_ident = runs.front().charge_ident;
    ITensor::vector outtens;
    if (runs.size() == 1) {
        // Single-APA: the historical filter output verbatim (no merge/normalize).
        ApaRun& run = runs.front();
        auto tens_live = Aux::TensorDM::as_tensors(*run.root_live, run.inpath + "/live");
        outtens.insert(outtens.end(), tens_live.begin(), tens_live.end());

        auto root_dead = Aux::TensorDM::as_pctree(*invec[0]->tensors(), run.inpath + "/dead");
        auto tens_dead = Aux::TensorDM::as_tensors(*root_dead, run.inpath + "/dead");
        outtens.insert(outtens.end(), tens_dead.begin(), tens_dead.end());
    }
    else {
        // Multi-APA: reproduce the standalone PointTreeMerging it replaces. Primary
        // = input[0]; concatenate the opted-in root PCs (opflash), take_children
        // from each source tree, then normalize the merged-tree local-PC keys.
        auto* root_live = runs.front().root_live.get();
        auto root_dead = Aux::TensorDM::as_pctree(*invec[0]->tensors(), runs.front().inpath + "/dead");
        for (std::size_t k = 1; k < runs.size(); ++k) {
            auto src_dead = Aux::TensorDM::as_pctree(*invec[k]->tensors(), runs[k].inpath + "/dead");
            merge_pct(root_live, runs[k].root_live.get(), m_root_pcs_to_merge);
            merge_pct(root_dead.get(), src_dead.get(), m_root_pcs_to_merge);
        }
        normalize_pctree_local_pcs(root_live);
        normalize_pctree_local_pcs(root_dead.get());

        std::string outpath = m_outpath;
        if (outpath.find("%") != std::string::npos) outpath = String::format(outpath, out_ident);
        auto tens_live = Aux::TensorDM::as_tensors(*root_live, outpath + "/live");
        auto tens_dead = Aux::TensorDM::as_tensors(*root_dead, outpath + "/dead");
        outtens.insert(outtens.end(), tens_live.begin(), tens_live.end());
        outtens.insert(outtens.end(), tens_dead.begin(), tens_dead.end());
    }
    out = Aux::TensorDM::as_tensorset(outtens, out_ident);

    // (calib / cathode-diag dumps moved above, before recompose_cluster_groups,
    // so they still see the matching-internal main/associated split.)

    // Per-event timing + memory. `took` is the full operator() wall time. RSS is the
    // resident set of the WHOLE matching process at this point (dominated by the
    // upstream clustering point clouds, NOT by QLMatching); `delta` is QLMatching's
    // own incremental footprint over its entry baseline (typically small/near-zero
    // since the matcher reuses the already-resident pctree).
    // MemUsage::current() = (virtual size, resident); .second is the RSS we want.
    // Debug-level (consistent with the MABC `MABC timing:` lines). The chain is
    // single-threaded, so this whole operator() wall is QLMatching compute, not log I/O.
    const auto mu1 = em.mu.current();
    log->debug("QLMatching timing: ident {} took {} ms, proc RSS {:.1f} MB (delta {:.1f} MB)",
               out_ident, em.tk.since().total_milliseconds(),
               mu1.second / 1024.0, (mu1.second - mu0.second) / 1024.0);
    // Top-level phase split (sum ≈ `took` above). prefit = all-APA run_one_apa_prefit
    // (see `QLtiming prefit` lines), fit = all-APA run_one_apa_fit (`QLtiming fit`),
    // xtpc_cull only fires multi-APA, output = tensor/merge assembly tail.
    const double t_out = ms_since(t_out0);
    log->debug("QLtiming operator: ident {} prefit {:.1f} xtpc_cull {:.1f} fit {:.1f} output {:.1f} ms",
               out_ident, t_prefit, t_cull, t_fit, t_out);

    ++m_count;
    return true;
}

// ============================================================================
// Per-APA matching pipeline. operator() builds one ApaRun (single-APA today)
// and calls run_one_apa(); the stages below were lifted verbatim from the old
// inline operator() body, reading config from the m_* members and threading all
// per-run state through ApaRun.
// ============================================================================

// Single-APA convenience driver (not used by operator(), which runs the three
// stages explicitly so cull_cross_tpc can flag between them). Kept self-consistent
// with that order: prefit (build only) -> cull_inconsistent -> fit.
void QLMatching::run_one_apa(ApaRun& run)
{
    run_one_apa_prefit(run);
    cull_inconsistent(run);
    run_one_apa_fit(run);
}

// Pre-fit half: bundles + maps (everything up to, but not including, the
// per-TPC consistency cull and the LASSO). cull_inconsistent is now run by
// operator() AFTER cull_cross_tpc (so cross-TPC flagging sees the full pre-cull
// bundle set), not here. Split out so cull_cross_tpc can run between all APAs'
// prefit and all APAs' cull+fit.
void QLMatching::run_one_apa_prefit(ApaRun& run)
{
    auto t = wallclock::now();
    build_opdet_mask(run);          const double t_mask   = ms_since(t); t = wallclock::now();
    read_flashes(run);              const double t_flash  = ms_since(t); t = wallclock::now();
    decompose_cluster_groups(run);  const double t_decomp = ms_since(t); t = wallclock::now();
    compute_geometry(run);          const double t_geom   = ms_since(t); t = wallclock::now();
    build_bundles(run);             const double t_bundle = ms_since(t); t = wallclock::now();
    build_bundle_maps(run);         const double t_maps   = ms_since(t);
    log->debug("QLtiming prefit: ident {} mask {:.1f} read_flashes {:.1f} decompose {:.1f} "
               "geometry {:.1f} build_bundles {:.1f} build_maps {:.1f} ms",
               run.charge_ident, t_mask, t_flash, t_decomp, t_geom, t_bundle, t_maps);
}

// Fit half: the two LASSO rounds + output assembly.
void QLMatching::run_one_apa_fit(ApaRun& run)
{
    auto t = wallclock::now();
    fit_round1(run);         const double t_r1   = ms_since(t); t = wallclock::now();
    fit_round2(run);         const double t_r2   = ms_since(t); t = wallclock::now();
    apply_matched_t0s(run);  const double t_t0   = ms_since(t); t = wallclock::now();
    write_opflash_pc(run);   const double t_pc   = ms_since(t);
    log->debug("QLtiming fit: ident {} fit_round1 {:.1f} fit_round2 {:.1f} apply_t0 {:.1f} write_pc {:.1f} ms",
               run.charge_ident, t_r1, t_r2, t_t0, t_pc);
}

// Per-channel OpDet on/off mask, derived from the injected OpDet table
// (SemiAnalyticalModel::OpticalDetector::type; 1 = dome PMT, 0 = (X)Arapuca)
// rather than a hard-coded SBND layout: channel i is on iff its type is in
// m_active_opdet_types (default {1} => PMTs only, reproducing the historical
// SBND mask). m_ch_mask then disables specific channels.
void QLMatching::build_opdet_mask(ApaRun& run)
{
    run.opdet_mask.assign(m_opdets.size(), 0);
    if (m_pmts) {
        for (std::size_t i = 0; i < m_opdets.size(); ++i) {
            if (std::find(m_active_opdet_types.begin(), m_active_opdet_types.end(),
                          m_opdets[i].type) != m_active_opdet_types.end())
                run.opdet_mask[i] = 1;
        }
    }
    for (std::size_t i = 0; i < m_ch_mask.size(); ++i) run.opdet_mask[m_ch_mask[i]] = 0;
}

// Flashes come from the canonical optical point clouds on the live root node
// (written by Aux::FlashTensorToOpticalPCs) via the shared facade enumerator
// Grouping::flashes(). Each flash is wrapped in an Opflash matching-adapter.
void QLMatching::read_flashes(ApaRun& run)
{
    PEErr pe_err_model{m_pe_err_floor, m_pe_err_frac, m_pe_err_knee};
    // Per-family override arrays (empty => Opflash keeps the scalar rule).
    pe_err_model.ch_floor = m_pe_err_ch_floor;
    pe_err_model.ch_frac  = m_pe_err_ch_frac;
    const std::vector<double>* pe_scale = m_measured_pe_scale.empty() ? nullptr : &m_measured_pe_scale;
    for (const auto& ff : run.grouping->flashes()) {
        auto flash = std::make_shared<Opflash>(ff, m_flash_pe_threshold, m_nchan, pe_err_model, pe_scale);
        // A readout-uncovered channel measured "< self-trigger threshold", not
        // "0 +- pe_err_floor"; widen its error to the threshold band before the
        // flash reaches the chi2/LASSO. No-op when pe_err_nodata <= 0 (default).
        flash->inflate_nodata_err(m_pe_err_nodata, m_coverage_min);
        if (flash->get_time() < m_flash_mintime || flash->get_time() > m_flash_maxtime) continue;
        if (flash->get_total_PE() < m_flash_minPE) continue;
        // Channel-scoped admission (see the knob doc in the header): the flash must
        // show real light on the trusted PD family, not just clear the all-PD floor.
        if (!m_flash_sel_channels.empty()) {
            double sel_sum = 0.0;
            int sel_fired = 0;
            for (int ch : m_flash_sel_channels) {
                if (ch < 0 || ch >= m_nchan) continue;
                const double p = flash->get_PE(ch);
                sel_sum += p;
                if (p >= m_flash_sel_fired_pe) ++sel_fired;
            }
            if (sel_sum < m_flash_sel_minPE || sel_fired < m_flash_sel_min_fired) continue;
        }
        run.flashes.push_back(flash);
    }

    // Register the run's own anode plus any extra group anodes (PDHD drift-side
    // group spanning >1 APA) so every blob's wpid resolves in get_anode().
    const std::vector<IAnodePlane::pointer>& grp = m_grouping_anodes.at(run.input_idx);
    if (grp.empty()) {
        run.grouping->set_anodes({run.anode});
    }
    else {
        std::vector<IAnodePlane::pointer> all = grp;
        all.push_back(run.anode);
        run.grouping->set_anodes(all);
    }
    run.grouping->set_detector_volumes(m_dv);
}

// Cluster-group decomposition (MicroBooNE-style main + associated). Per-APA
// clustering (examine_bundles) merges a group into a single Facade Cluster
// carrying the "isolated"/"perblob" per-blob array: blobs tagged -1 are the main
// cluster, other ids are sub-clusters. Split each such group back into separate
// Facade Clusters via Grouping::separate(); bundle building anchors on the main.
void QLMatching::decompose_cluster_groups(ApaRun& run)
{
    auto* grouping = run.grouping;

    auto sort_by_length_then_id = [](const Cluster* a, const Cluster* b) {
        if (a->get_length() != b->get_length()) return a->get_length() > b->get_length();
        return a->get_cluster_id() < b->get_cluster_id();
    };
    // Snapshot the pre-split clusters (separate() appends new siblings; iterate a
    // copy). Sub-cluster idents are drawn from a counter seeded above every
    // existing ident, so they are unique and collision-free with the mains. The
    // counter advances in children() (tree/serialization) order, so it is build-stable.
    std::vector<Cluster*> original_clusters = grouping->children();
    int next_sub_ident = 0;
    for (auto* c : original_clusters) next_sub_ident = std::max(next_sub_ident, c->ident());
    ++next_sub_ident;
    for (auto* cluster : original_clusters) {
        if (!cluster->has_pcarray("isolated", "perblob")) {
            run.match_groups.emplace_back(cluster, std::vector<Cluster*>{});
            continue;
        }
        auto cc_span = cluster->get_pcarray("isolated", "perblob");
        std::vector<int> cc(cc_span.begin(), cc_span.end());
        const bool has_main = std::find(cc.begin(), cc.end(), -1) != cc.end();
        std::set<int> subs;
        for (int v : cc) if (v >= 0) subs.insert(v);
        if (!has_main || subs.empty()) {   // single component -> nothing to split
            run.match_groups.emplace_back(cluster, std::vector<Cluster*>{});
            continue;
        }
        // separate(): groups with id<0 stay in `cluster` (the main), each id>=0
        // becomes a new sibling cluster returned (ascending id -> deterministic).
        auto splits = grouping->separate(cluster, cc);
        cluster->set_flag("main_cluster");
        std::vector<Cluster*> others;
        for (auto& [gid, nc] : splits) {
            (void)gid;
            nc->set_ident(next_sub_ident++);
            nc->set_flag("associated_cluster");
            others.push_back(nc);
        }
        run.match_groups.emplace_back(cluster, std::move(others));
    }
    // Deterministic match-unit order: by main length (desc), ident tie-break.
    std::sort(run.match_groups.begin(), run.match_groups.end(),
              [&](const auto& a, const auto& b) { return sort_by_length_then_id(a.first, b.first); });

    // All separated clusters (mains + associated): charge bookkeeping, scalar
    // reset, and the global (debug/iteration) index. Sorted by length (desc) with
    // a stable ident tie-break (std::sort is not stable).
    run.clusters = grouping->children();
    std::sort(run.clusters.begin(), run.clusters.end(), sort_by_length_then_id);

    for (auto cluster : run.clusters) {
        for (auto blob : cluster->children()) run.total_charge_blob_all += blob->charge();
    }

    std::for_each(run.clusters.begin(), run.clusters.end(),
                  [](Cluster* c) {
                      c->set_cluster_t0(-1e12);
                      c->set_scalar<int>("flash", -1);
                      c->set_scalar<int>("matched_flash_gid", -1);
                      // Materialize the main/associated flag scalars on EVERY
                      // cluster so the cluster_scalar PC has uniform keys at
                      // as_tensors() time. get_flag() returns the current value
                      // (1 if the split set it) or 0 if unset -- no clobbering.
                      c->set_flag("main_cluster", c->get_flag("main_cluster"));
                      c->set_flag("associated_cluster", c->get_flag("associated_cluster"));
                  });

    for (std::size_t i = 0; i < run.flashes.size(); ++i) run.global_flash_idx_map[run.flashes[i].get()] = i;
    for (std::size_t i = 0; i < run.clusters.size(); ++i) run.global_cluster_idx_map[run.clusters[i]] = i;
}

// Inverse of decompose_cluster_groups(): merge each match-group's associated
// sub-clusters back into their main, restoring the stage-3 cluster grouping
// before the matched tree is serialized. decompose() split via Grouping::separate()
// ONLY to anchor each bundle's geometry on its main cluster; that split is a
// matching-internal device and must NOT leak into the emitted cluster tree, where
// it would land in the downstream all-TPC clustering (PDHD's stage 4 has no
// flash-T0 merge to re-absorb it, so a single stage-3 cluster came out split —
// see match/QLMatching split-cluster bug). After this, grouping->children() ==
// the stage-3 cluster set: each main keeps its original ident and the matched
// cluster_t0 / flash scalars written by apply_matched_t0s, so QLMatching is a pure
// T0/flash ANNOTATION step. The matching results are untouched. Groups with no
// associated sub-clusters (the common / single-component / single-anode case) are
// a no-op. Must run AFTER dump_calib / dump_cathode_diag (which read the split
// main+associated clusters) and BEFORE output assembly.
void QLMatching::recompose_cluster_groups(ApaRun& run)
{
    auto* grouping = run.grouping;
    for (auto& [main_cluster, others] : run.match_groups) {
        if (others.empty()) continue;
        // merge(): main adopts every associated cluster's blobs (take_children),
        // then the emptied associated facades are destroyed (keep=false default).
        grouping->merge(others.begin(), others.end(), main_cluster);
    }
    // run.clusters / global_cluster_idx_map / the bundles now hold dangling pointers
    // to the destroyed associated clusters; nothing reads them after this point
    // (output assembly walks the tree via root_live only).
}

// Active-volume geometry from the IDetectorVolumes service, in a per-TPC drift
// coordinate u = s*(x-anode_x) (u=0 at the anode, u=u_cathode at the cathode) so
// the prototype's single-TPC inequalities port directly and both reversed-drift
// SBND APAs share one set of cushions. Also reduces the OpDet mask to this TPC
// and caches the kept-channel index used to size the LASSO matrices.
void QLMatching::compute_geometry(ApaRun& run)
{
    const unsigned int tpc = run.anode->ident();
    run.sign_offset  = (tpc == 0) ? -1 : 1;

    // Active-volume bounds are the UNION of the sensitive boxes of every APA in this
    // drift-side group. PDHD images two APAs per drift volume, offset in Z (e.g. APA0
    // z[-1,2306] mm + APA2 z[2320,4626] mm share the -x side); using only the
    // representative anode's box would gate out all charge in the group's other APA (a
    // different Z half) from the predicted-light sum at build_bundles. With an empty
    // grouping list (SBND / single-APA) this reduces to the representative anode's box
    // => bit-identical. X is identical across same-side APAs, so the drift math
    // (anode_x / u_cathode / s) is unchanged; only the Y/Z extent widens.
    // Per-input imaging face (joint multi-side node) or the single m_tpc_face. The
    // run's group anodes are all on its own drift side, so they share this face.
    const int face = m_tpc_faces.empty() ? m_tpc_face : m_tpc_faces.at(run.input_idx);
    // Optional second face to union in (PDVD: its two faces split Y in half rather
    // than duplicating the full Y/Z like PDHD/SBND's do; see m_tpc_extra_faces).
    // -1 => none => bit-identical for PDHD/SBND.
    const int extra_face = m_tpc_extra_faces.empty() ? -1 : m_tpc_extra_faces.at(run.input_idx);
    const WirePlaneId wpid(WirePlaneLayer_t::kAllLayers, face, static_cast<int>(tpc));
    const std::vector<IAnodePlane::pointer>& grp = m_grouping_anodes.at(run.input_idx);
    const std::vector<IAnodePlane::pointer> group_anodes =
        grp.empty() ? std::vector<IAnodePlane::pointer>{run.anode} : grp;
    BoundingBox bb;
    for (const auto& a : group_anodes) {
        const WirePlaneId awpid(WirePlaneLayer_t::kAllLayers, face,
                                static_cast<int>(a->ident()));
        const BoundingBox abb = m_dv->inner_bounds(awpid);
        if (!abb.empty()) bb(abb.bounds());
        if (extra_face >= 0) {
            const WirePlaneId awpid2(WirePlaneLayer_t::kAllLayers, extra_face,
                                     static_cast<int>(a->ident()));
            const BoundingBox abb2 = m_dv->inner_bounds(awpid2);
            if (!abb2.empty()) bb(abb2.bounds());
        }
    }
    if (bb.empty()) {
        raise<ValueError>("QLMatching: empty detector-volume bounds for anode %d", tpc);
    }
    const Ray bray = bb.bounds();
    const double x_lo = bray.first.x(), x_hi = bray.second.x();
    const bool lo_is_cathode = std::abs(x_lo - m_cathode_x) < std::abs(x_hi - m_cathode_x);
    const double cathode_x = lo_is_cathode ? x_lo : x_hi;
    run.anode_x   = lo_is_cathode ? x_hi : x_lo;
    run.s         = (run.anode_x < cathode_x) ? +1.0 : -1.0;
    run.u_cathode = run.s * (cathode_x - run.anode_x);   // > 0
    // Y/Z active bounds from the same bbox, shrunk(+)/grown(-) by the cushions.
    run.y_lo = bray.first.y()  + m_y_cushion;
    run.y_hi = bray.second.y() - m_y_cushion;
    run.z_lo = bray.first.z()  + m_z_cushion;
    run.z_hi = bray.second.z() - m_z_cushion;
    // Per-TPC transverse position offset (Y,Z), read from the SAME DetectorVolumes
    // metadata the post-QLMatching T0Correction scope uses (single source of truth;
    // see match/docs/cathode-offset-correction.md).  Parked on the ApaRun for the
    // forthcoming cross-TPC matching judgement -- NOT yet consumed here (the photon
    // model and active-volume gate still read raw y,z).  Absent key => 0 => inert.
    run.dy = 0.0;
    run.dz = 0.0;
    {
        const auto md = m_dv->metadata(wpid);
        if (md.isMember("pos_offset") && md["pos_offset"].isArray() &&
            md["pos_offset"].size() >= 3) {
            run.dy = md["pos_offset"][1].asDouble();
            run.dz = md["pos_offset"][2].asDouble();
        }
    }
    log->debug("anode {} pos_offset (dy,dz) = ({:.3f},{:.3f}) cm [parked, not yet applied]",
               tpc, run.dy / units::cm, run.dz / units::cm);
    log->debug("anode {} group-bbox ({} apa) x[{:.2f},{:.2f}] y[{:.2f},{:.2f}] z[{:.2f},{:.2f}] cm; "
               "anode_x {:.2f} cathode_x {:.2f} u_cathode {:.2f} cm s {}",
               tpc, group_anodes.size(), x_lo / units::cm, x_hi / units::cm,
               bray.first.y() / units::cm, bray.second.y() / units::cm,
               bray.first.z() / units::cm, bray.second.z() / units::cm,
               run.anode_x / units::cm, cathode_x / units::cm, run.u_cathode / units::cm, run.s);

    // Bundle-quality thresholds forwarded to every TimingTPCBundle below.
    run.qp = BundleQualityParams{
        m_bundle_ks_merge_max, m_bundle_chi2ndf_merge_max, m_bundle_addmerge_exponent,
        m_highconsist_ks_max, m_highconsist_min_ndf, m_bundle_pe_ndf_knee,
        m_bundle_mask_ks, m_pe_err_floor, m_pe_err_frac, m_pe_err_knee, m_pe_err_on_pred,
        m_pe_err_lowpe_frac, m_pe_err_lowpe_knee,
        m_highconsist_ladder,
        m_hc_clean_ks, m_hc_clean_c2, m_hc_good_ks, m_hc_good_c2,
        m_hc_tb_ks, m_hc_tb_c2, m_hc_miss_ks, m_hc_miss_c2, m_hc_miss_min_ndf,
        m_chi2_relax, m_chi2_pmt_excess, m_chi2_pmt_ratio, m_chi2_pmt_inflate,
        m_chi2_sat_inflate};
    // Per-channel family overrides sit after the positional members (see
    // BundleQualityParams); empty (default) => scalar model, bit-identical.
    run.qp.pe_err_ch_floor      = m_pe_err_ch_floor;
    run.qp.pe_err_ch_frac       = m_pe_err_ch_frac;
    run.qp.pe_err_ch_lowpe_frac = m_pe_err_ch_lowpe_frac;
    run.qp.pe_err_ch_lowpe_knee = m_pe_err_ch_lowpe_knee;

    // Shared-flash mode relies on the ident-based sign_offset above encoding the
    // physical relation sign_offset == -s (true for SBND TPC0/1 and the PDVD
    // bottom/top representative anodes 0/4); guard against a future detector
    // whose idents break it.
    if (m_shared_flash && run.sign_offset != -(int)run.s) {
        raise<ValueError>("QLMatching: shared_flash sign_offset %d != -s %d for anode %d",
                          run.sign_offset, -(int)run.s, tpc);
    }

    // Reduce mask to OpDets on this TPC: an OpDet belongs to TPC 0 / TPC 1 if it
    // sits on the low / high side of the cathode plane. opdet_all_volumes (PDVD:
    // double-sided cathode XAs at x=0, one flash over both drift volumes) keeps
    // every active OpDet on every run instead.
    if (!m_opdet_all_volumes) {
        for (std::size_t idet = 0; idet < run.opdet_mask.size(); ++idet) {
            const bool low_side = m_opdets[idet].center.x() < m_cathode_x;
            if ((tpc == 0) && !low_side) run.opdet_mask[idet] = 0;
            if ((tpc == 1) &&  low_side) run.opdet_mask[idet] = 0;
        }
    }

    // Per-event dynamic dead-PMT auto-mask (optional; off => bit-identical). Folds into
    // run.opdet_mask before opdet_idx_v below, so the fit inherits it consistently.
    if (m_auto_mask) compute_dynamic_opdet_mask(run, tpc);

    // Kept-channel index (the surviving on-channels), used to size the LASSO matrices.
    run.nopdet = 0;
    run.opdet_idx_v.clear();
    for (std::size_t idet = 0; idet < run.opdet_mask.size(); ++idet) {
        if (run.opdet_mask.at(idet) == 1) {
            run.opdet_idx_v.push_back(int(idet));
            ++run.nopdet;
        }
    }
    log->debug("nopdet {} opdet_idx_v size {}", run.nopdet, run.opdet_idx_v.size());
}

// Per-event dynamic dead-PMT auto-mask. Within this single event/TPC, drop a PMT that
// never fires (max measured PE over the event's flashes < pe_low) while its nearest live
// PMTs do see light (in >= min_contrast flashes the median PE of its K nearest live PMTs
// exceeds pe_bright). This catches a channel that is dead in THIS run but absent from the
// static ch_mask, without touching one that is merely in a quiet region (no bright
// neighbours => not masked). All decisions come from a snapshot (per-channel event max +
// the live reference pool), applied after the scan, in channel order, with (distance,
// channel) tie-breaking — deterministic and Python-parity with sbnd_xin/automask_prototype.py.
void QLMatching::compute_dynamic_opdet_mask(ApaRun& run, unsigned int tpc)
{
    if (run.flashes.size() < (std::size_t)m_auto_mask_min_flash) return;
    const int nchan = m_nchan;

    // Per-channel max measured PE over this event's flashes.
    std::vector<double> maxpe(nchan, 0.0);
    for (const auto& fl : run.flashes) {
        const int nc = std::min(nchan, fl->get_num_channels());
        for (int ch = 0; ch < nc; ++ch)
            maxpe[ch] = std::max(maxpe[ch], fl->get_PE(ch));
    }

    auto is_pmt = [&](int ch) {
        return std::find(m_active_opdet_types.begin(), m_active_opdet_types.end(),
                         m_opdets[ch].type) != m_active_opdet_types.end();
    };
    auto in_tpc = [&](int ch) {
        if (m_opdet_all_volumes) return true;   // all PDs belong to every run (PDVD)
        const bool low_side = m_opdets[ch].center.x() < m_cathode_x;
        return (tpc == 0) ? low_side : !low_side;
    };

    // Live brightness reference: PMTs in this TPC that fire somewhere this event
    // (includes ones that happen to be in the static ch_mask but are alive this run).
    std::vector<int> live;
    for (int ch = 0; ch < nchan && ch < (int)m_opdets.size(); ++ch)
        if (is_pmt(ch) && in_tpc(ch) && maxpe[ch] >= m_auto_mask_pe_low)
            live.push_back(ch);

    std::vector<int> to_mask;
    for (int ch = 0; ch < nchan && ch < (int)run.opdet_mask.size(); ++ch) {
        if (run.opdet_mask[ch] != 1) continue;      // only currently-active PMTs
        if (maxpe[ch] >= m_auto_mask_pe_low) continue;  // it fires => healthy

        // K nearest live PMTs by (Y,Z) distance; tie-break by channel for determinism.
        const auto& o = m_opdets[ch].center;
        std::vector<std::pair<double,int>> d;
        d.reserve(live.size());
        for (int c : live) {
            if (c == ch) continue;
            // Same-type pool (PDVD: XA vs PMT efficiencies differ ~4x, so a
            // cross-type (Y,Z) neighbour is not a meaningful brightness reference).
            if (m_auto_mask_same_type && m_opdets[c].type != m_opdets[ch].type) continue;
            const double dy = m_opdets[c].center.y() - o.y();
            const double dz = m_opdets[c].center.z() - o.z();
            d.push_back({dy * dy + dz * dz, c});
        }
        if (d.empty()) continue;
        const int K = std::min((int)d.size(), m_auto_mask_neighbors);
        std::partial_sort(d.begin(), d.begin() + K, d.end(),
            [](const std::pair<double,int>& a, const std::pair<double,int>& b) {
                return a.first != b.first ? a.first < b.first : a.second < b.second; });

        // Count flashes where the neighbour-median PE exceeds the bright threshold.
        int contrast = 0;
        std::vector<double> vals(K);
        for (const auto& fl : run.flashes) {
            for (int k = 0; k < K; ++k) vals[k] = fl->get_PE(d[k].second);
            std::sort(vals.begin(), vals.end());
            const double med = (K % 2) ? vals[K/2] : 0.5 * (vals[K/2 - 1] + vals[K/2]);
            if (med > m_auto_mask_pe_bright) ++contrast;
        }
        if (contrast >= m_auto_mask_min_contrast) to_mask.push_back(ch);
    }

    for (int ch : to_mask) {
        run.opdet_mask[ch] = 0;
        run.auto_masked.insert(ch);
        log->debug("QLAUTOMASK tpc {} ch {} auto-masked (event max PE {:.2f} < {:.2f})",
                   tpc, ch, maxpe[ch], m_auto_mask_pe_low);
    }
    if (!to_mask.empty())
        log->info("QLAUTOMASK tpc {} auto-masked {} dead-in-event PMT(s) over {} flashes",
                  tpc, to_mask.size(), run.flashes.size());
}

// Opaque-cathode mismatched-candidate filter (see header). The flash's lit drift
// side is read from where its measured PE actually is -- OpDet x vs the cathode plane
// -- the SAME rule the merged-root opflash PC uses to tag a flash's side (the flash's
// gid encodes the processing node's anode, not its lit volume). A cross-side bundle is
// kept only when the cluster is a CATHODE-CROSSER at this flash's T0 (at_x_boundary):
// only then can its far half be the source of the opposite side's light. Every other
// cross-side bundle is dropped; same-side bundles are always kept.
bool QLMatching::cross_side_mismatch_drop(const Opflash* flash, int cluster_side,
                                          bool at_x_boundary) const
{
    if (!m_cross_side_filter) return false;
    double pe_lo = 0.0, pe_hi = 0.0;
    for (int ch = 0; ch < m_nchan; ++ch) {
        const double p = flash->get_PE(ch);
        if (ch < (int) m_opdets.size() && m_opdets[ch].center.x() >= m_cathode_x) pe_hi += p;
        else pe_lo += p;
    }
    const int flash_side = (pe_hi > pe_lo) ? 1 : 0;   // lit drift side (ties -> side 0)
    if (flash_side == cluster_side) return false;       // same-side: always keep
    return !at_x_boundary;                              // cross-side: keep only cathode-crossers
}

// [Stage 1] Build a candidate (flash, cluster-group) bundle for every pair,
// predict its light (summed over the whole group), fill the boundary flags, and
// drop the KS==1 degenerate (no measured/predicted overlap). Surviving bundles
// land in run.pre_bundles.
void QLMatching::build_bundles(ApaRun& run)
{
    // Profiling: accumulate time spent in the per-point visibility/PE loop (the
    // SemiAnalyticalModel inner kernel, the prime suspect) vs the rest of the
    // bundle assembly. Timed at per-bundle granularity to avoid perturbing the
    // millions of per-opdet calls. Logging only; no effect on results.
    double vis_ms = 0.0;
    std::size_t total_pts = 0;
    // Reuse the visibility buffers across all points (detected*Visibilities re-zero
    // them via assign), avoiding a per-point allocation pair.
    std::vector<double> direct_visibilities;
    std::vector<double> reflected_visibilities;

    // Flash-independent per-cluster point cache for the vis loop below. The
    // legacy loop re-fetched blob->points("3d") per (flash x group x blob) --
    // a string-keyed dataset lookup plus a fresh point-array copy for data that
    // only ever shifts by the per-flash constant flash_x_offset (round-0
    // profile, pdvd/docs/15_pdvd-light-ql-perf.md). Fetch each cluster's blob
    // coordinates once per ApaRun into contiguous arrays (same values, same
    // blob and point order) and record each blob's coordinate bounds for the
    // conservative AABB pre-gate in the loop. Keyed lookup only (never
    // iterated); local to this call, freed on return.
    struct BlobPts {
        double charge;
        int npts;
        std::vector<double> xs, ys, zs;
        double xlo, xhi, ylo, yhi, zlo, zhi;
    };
    std::unordered_map<const Cluster*, std::vector<BlobPts>> cluster_pts;
    for (const auto& mg : run.match_groups) {
        std::vector<Cluster*> gcs{mg.first};
        gcs.insert(gcs.end(), mg.second.begin(), mg.second.end());
        for (auto* gc : gcs) {
            auto [it, fresh] = cluster_pts.try_emplace(gc);
            if (!fresh) continue;
            auto& blobs = it->second;
            for (auto blob : gc->children()) {
                BlobPts b;
                b.charge = blob->charge();
                b.npts = blob->npoints();
                b.xlo = b.ylo = b.zlo = std::numeric_limits<double>::infinity();
                b.xhi = b.yhi = b.zhi = -std::numeric_limits<double>::infinity();
                auto points = blob->points("3d", {"x", "y", "z"});
                b.xs.reserve(b.npts);
                b.ys.reserve(b.npts);
                b.zs.reserve(b.npts);
                for (int i = 0; i < b.npts; ++i) {
                    const auto& p = points.at(i);
                    const double x = p.x(), y = p.y(), z = p.z();
                    b.xs.push_back(x);
                    b.ys.push_back(y);
                    b.zs.push_back(z);
                    if (x < b.xlo) b.xlo = x;
                    if (x > b.xhi) b.xhi = x;
                    if (y < b.ylo) b.ylo = y;
                    if (y > b.yhi) b.yhi = y;
                    if (z < b.zlo) b.zlo = z;
                    if (z > b.zhi) b.zhi = z;
                }
                blobs.push_back(std::move(b));
            }
        }
    }

    for (auto flash : run.flashes) {
        const auto flash_time = flash->get_time();
        // trigger_offset_for folds the per-event readout-vs-trigger offset into the
        // drift x correction for detectors that don't bake it into the charge x at
        // imaging time (default 0 => bit-identical); per-input when the charge
        // windows start per-crate (PDVD TDE/BDE). See QLMatching.h.
        const double flash_x_offset = run.sign_offset * (flash_time + trigger_offset_for(run.input_idx)) * drift_speed_for(run.input_idx);

        // per-flash mask (also catches simulated saturated PMTs in MC).
        std::vector<unsigned int> flash_opdet_mask = run.opdet_mask;
        // Prediction-evaluation mask: identical to flash_opdet_mask EXCEPT
        // that the data DAPHNE-rail drop is not applied, so the full
        // prediction on rail-flagged channels survives for the calib dump /
        // viewer (the flag only excludes channels from chi2/KS/LASSO; a
        // pred shown as 0 on a railed channel is confusing and wrong).
        // Diverges from flash_opdet_mask only when use_saturation_flag is
        // on => the legacy path is bit-identical.
        std::vector<unsigned int> vis_opdet_mask = run.opdet_mask;
        for (std::size_t idet = 0; idet < std::size_t(flash->get_num_channels()); ++idet) {
            auto pe_det = flash->get_PE(idet);
            if ((flash->get_total_PE() > m_mc_saturation_pe) && (pe_det == 0) && (m_data == false)) {
                flash_opdet_mask[idet] = 0;
                vis_opdet_mask[idet] = 0;
            }
            // Data DAPHNE-rail flag (flag_saturation chain).  Under
            // saturation_mask_fit the clipped channel leaves this flash's
            // pred/chi2/KS (bundle_mask_ks); otherwise it stays in them at its
            // measured PE, which is a LOWER BOUND on the true light (and a
            // repaired estimate when OpDecon saturation_repair is on) rather
            // than a missing measurement -- see QLMatching.h
            // saturation_mask_fit / chi2_sat_inflate.  The LASSO rows are
            // zeroed separately and stay zeroed either way: a lower bound must
            // not constrain a least-squares solve.  Default (mask on)
            // bit-identical.
            if (m_use_saturation_flag && m_saturation_mask_fit && flash->get_sat((int)idet))
                flash_opdet_mask[idet] = 0;
            // Readout-coverage flag (emit_coverage chain): a self-trigger
            // channel with no waveform over this flash's window reads 0 PE.
            // Under coverage_mask_fit it leaves pred/chi2/KS/LASSO like a
            // railed channel; otherwise it stays in the fit at that 0, which
            // is what the hardware reports -- the ~1 PE self-trigger
            // threshold makes the silence an upper limit, not a missing
            // measurement (QLMatching.h coverage_mask_fit; pair with
            // pe_err_nodata for the error).  Either way the full prediction
            // survives in vis_opdet_mask for the calib dump / viewer.
            // Default (mask on) bit-identical.
            if (m_use_coverage_flag && m_coverage_mask_fit
                && flash->get_cov((int)idet) < m_coverage_min)
                flash_opdet_mask[idet] = 0;
        }

        log->debug("flash time {} flash PE {} flash_x_offset {}",
                   int(flash_time) / 100.,
                   int(flash->get_total_PE() * 100) / 100.,
                   int(flash_x_offset * 100) / 100.);

        // Lit drift side of this flash from its PE pattern (the same test as
        // cross_side_mismatch_drop), hoisted out of the per-group loop so the
        // cross-side skip below can run BEFORE the per-point visibility loop. Only
        // computed when the hoisted skip is active; left at 0 (and unused) otherwise.
        int flash_side = 0;
        if (m_cross_side_filter && m_crossside_skip_vis) {
            double pe_lo = 0.0, pe_hi = 0.0;
            for (int ch = 0; ch < m_nchan; ++ch) {
                const double p = flash->get_PE(ch);
                if (ch < (int) m_opdets.size() && m_opdets[ch].center.x() >= m_cathode_x) pe_hi += p;
                else pe_lo += p;
            }
            flash_side = (pe_hi > pe_lo) ? 1 : 0;   // lit drift side (ties -> side 0)
        }

        for (std::size_t ig = 0; ig < run.match_groups.size(); ++ig) {
            Cluster* main_cluster = run.match_groups[ig].first;
            const auto& others = run.match_groups[ig].second;
            // Stable cluster index of the group anchor (the main cluster).
            const int cidx = run.global_cluster_idx_map.at(main_cluster);
            auto bundle = std::make_shared<TimingTPCBundle>(
                flash.get(), main_cluster, flash->get_flash_id(), cidx);
            bundle->set_quality_params(run.qp);
            run.all_bundles.push_back(bundle);
            bundle->set_opdet_mask(flash_opdet_mask);
            // Attach the associated sub-clusters so the bundle carries the whole
            // group (MicroBooNE-style main + others); predicted light below sums
            // over all of them.
            for (auto* oc : others) bundle->add_other_cluster(oc);

            // Fill the boundary flags (close-to-PMT / at-x-boundary / spec-end)
            // from the main cluster's drift endpoints (the group anchor). The
            // return is the prototype flag_good_bundle / TPC-containment verdict.
            const bool contained =
                compute_endpoint_flags(bundle.get(), main_cluster, flash_x_offset, run);
            bundle->set_contained(contained);   // diagnostic only (calib dump)
            // flag_two_boundary (both PCA ends at a detector edge). The calib dump and the
            // high-consistent ladder's B3 branch consume it; compute it only when one of those
            // needs it — production with the ladder off and no dump stays bit-identical.
            if (!m_calib_dump.empty() || m_highconsist_ladder)
                compute_two_boundary_flag(bundle.get(), main_cluster, flash_x_offset, run);
            // Discard bundles whose cluster is not contained in the TPC box once
            // the flash T0 x-offset is applied. Off by default. xtpc cathode rescue
            // (xtpc_cathode_tol > 0, else the flag is never set): a provisional
            // cathode-overshoot bundle is kept through the vis/eval below so
            // cull_cross_tpc can try to confirm it as a crosser half; operator()
            // purges it before cull_inconsistent/fit if unconfirmed.
            if (m_require_containment && !contained &&
                !bundle->get_flag_xtpc_cathode_provisional()) continue;

            // Opaque-cathode cross-side mismatch (hoisted form): a flash lit on the
            // opposite drift side from this cluster can only match a cathode-crosser
            // (at_x_boundary). The verdict uses only the flash lit-side, the run's fixed
            // cluster side, and the already-computed at_x_boundary flag -- it never
            // touches the predicted light -- so when crossside_skip_vis is on we drop it
            // HERE, before the expensive per-point visibility loop (exactly as
            // require_containment is hoisted above), skipping the visibility of bundles
            // the post-loop test below would discard anyway. Same final candidate set,
            // so the matching is bit-identical; the time saved is the vis_loop cost of
            // the cross-side-but-contained bundles. Default OFF => the post-loop test
            // runs instead and the result is byte-for-byte the legacy path.
            if (m_cross_side_filter && m_crossside_skip_vis
                    && flash_side != (int) run.anode->ident()
                    && !bundle->get_flag_at_x_boundary()) {
                bundle->set_potential_bad_match_flag(true);
                continue;
            }

            const std::size_t nopdets = flash->get_num_channels();
            std::vector<double> pred_flash(nopdets, 0.0);

            // Predicted light is summed over the whole group (main + associated).
            std::vector<Cluster*> group_clusters{main_cluster};
            group_clusters.insert(group_clusters.end(), others.begin(), others.end());
            std::size_t npt = 0;
            for (auto* gc : group_clusters) npt += gc->npoints();
            total_pts += npt;
            // Option D coarsening (default stride 1 = OFF): sample every stride-th
            // point of each blob and weight it by the run of points it stands in
            // for, conserving each blob's total charge. Active only on groups at
            // least vis_sample_min_pts big. At stride 1 the weight is exactly 1
            // point (qw == q, no extra FP op) and the iteration is unchanged, so
            // the default path is bit-identical.
            const int vstride =
                (m_vis_sample_stride > 1 && npt >= std::size_t(m_vis_sample_min_pts))
                    ? m_vis_sample_stride : 1;
            const auto t_vis0 = wallclock::now();
            for (auto* gc : group_clusters) {
              for (const auto& b : cluster_pts.find(gc)->second) {
                run.total_charge_blob += b.charge;
                const double q = b.charge / b.npts;

                // Conservative AABB pre-gate against the PE-inclusion box at
                // this flash's offset. Every per-coordinate op below (add a
                // constant, subtract a constant, multiply by run.s) is monotone
                // in FP, so u of every blob point lies within [ulo, uhi] built
                // from the cached x bounds, and a box strictly outside the gate
                // on any one axis implies every point fails that axis' test
                // exactly as the per-point compares would. Skipped blobs replay
                // only the total_charge_point accumulation, which never touched
                // the coordinates -- the running sums, and hence the debug
                // totals, are FP-identical to the legacy per-point loop.
                const double ua = run.s * ((b.xlo + flash_x_offset) - run.anode_x);
                const double ub = run.s * ((b.xhi + flash_x_offset) - run.anode_x);
                const double ulo = std::min(ua, ub);
                const double uhi = std::max(ua, ub);
                if (uhi < m_anode_ext1 || ulo > run.u_cathode + m_cathode_ext1 ||
                    b.yhi < run.y_lo || b.ylo > run.y_hi ||
                    b.zhi < run.z_lo || b.zlo > run.z_hi) {
                    for (int i = 0; i < b.npts; i += vstride) {
                        const double qw =
                            (vstride == 1) ? q : q * double(std::min(vstride, b.npts - i));
                        run.total_charge_point += qw;
                    }
                    continue;
                }

                for (int i = 0; i < b.npts; i += vstride) {
                    const double qw =
                        (vstride == 1) ? q : q * double(std::min(vstride, b.npts - i));
                    run.total_charge_point += qw;
                    const double x = b.xs[i] + flash_x_offset;
                    const double y = b.ys[i];
                    const double z = b.zs[i];

                    // PE-inclusion gate in the per-TPC drift coordinate u.
                    const double u = run.s * (x - run.anode_x);
                    if (u < m_anode_ext1 || u > run.u_cathode + m_cathode_ext1) continue;
                    if (y < run.y_lo || y > run.y_hi || z < run.z_lo || z > run.z_hi) continue;

                    // SemiAnalyticalModel expects positions in cm; blob points are mm.
                    const WireCell::Point xyz_cm(x / units::cm, y / units::cm, z / units::cm);
                    if (m_lib_model) {
                        // Gridded library holds total photon arrival (incl. any
                        // reflected component) with the detector's own optical
                        // shadowing; reflected term stays 0.
                        m_lib_model->visibilities(direct_visibilities, xyz_cm, &vis_opdet_mask);
                        reflected_visibilities.assign(nopdets, 0.0);
                    }
                    else {
                    // Skip masked opdets inside the model: their predicted light is
                    // discarded by the vis_opdet_mask gate below, so evaluating the
                    // per-opdet solid-angle/Gaisser-Hillas correction for them is wasted
                    // (~half of the same-TPC opdets are masked). Bit-identical result
                    // (vis_opdet_mask == flash_opdet_mask unless use_saturation_flag).
                    m_semi_model->detectedDirectVisibilities(direct_visibilities, xyz_cm, &vis_opdet_mask);
                    m_semi_model->detectedReflectedVisibilities(reflected_visibilities, xyz_cm, &vis_opdet_mask);
                    }

                    for (std::size_t idet = 0; idet < nopdets; ++idet) {
                        if (vis_opdet_mask.at(idet) == 0) continue;
                        const auto dir_vis = direct_visibilities.at(idet);
                        const auto ref_vis = reflected_visibilities.at(idet);
                        const auto dir_eff = m_VUVEfficiency.at(idet);
                        const auto ref_eff = m_VISEfficiency.at(idet);
                        pred_flash.at(idet) +=
                            qw * m_QtoL * dir_vis * dir_eff + qw * m_QtoL * ref_vis * ref_eff;
                    }
                }
              }
            }
            vis_ms += ms_since(t_vis0);

            // Per-PMT non-linearity: map the predicted PE total into the saturated
            // (observed) space so it matches the post-saturation measured PE. Acts only
            // above the knee; identity (no-op) when disabled or beta=1/gamma=0. See the
            // header and match/docs/sbnd-opdetsim-chain.md.
            if (m_pmt_nonlinearity) {
                for (std::size_t idet = 0; idet < nopdets; ++idet) {
                    const double p = pred_flash.at(idet);
                    if (p <= m_pmt_nl_knee) continue;
                    const double beta  = idet < m_pmt_nl_beta.size()  ? m_pmt_nl_beta[idet]  : 1.0;
                    const double gamma = idet < m_pmt_nl_gamma.size() ? m_pmt_nl_gamma[idet] : 0.0;
                    const double L = std::log(p / m_pmt_nl_knee);
                    pred_flash.at(idet) = m_pmt_nl_knee * std::exp(beta * L + gamma * L * L);
                }
            }

            // Full (per-flash-flag-unmasked) prediction for the calib dump /
            // viewer; the fit continues on the masked vector below. Stored
            // only when a per-flash flag knob is on (get_pred_flash_full
            // falls back to the fit vector otherwise), so the legacy path is
            // untouched.
            if (m_use_saturation_flag || m_use_coverage_flag) {
                bundle->set_pred_flash_full(pred_flash);
                for (std::size_t idet = 0; idet < nopdets; ++idet)
                    if (flash_opdet_mask.at(idet) == 0) pred_flash.at(idet) = 0.0;
            }

            // Cluster-group selection drives matching; the KS==1 guard below is
            // kept (degenerate, no measured/predicted overlap).
            bundle->set_pred_flash(pred_flash);
            bundle->examine_bundle();
            if (bundle->get_ks_dis() == 1) {
                bundle->set_potential_bad_match_flag(true);
                continue;
            }

            // Light-pattern over-prediction prefilter (prototype fired-fraction
            // reject, FlashTPCBundle.cxx 547-602). Drop a bundle whose predicted
            // light hugely exceeds the measured light over the SAME masked PMT set
            // the chi2 uses (R_total = sum pred/meas; R_max = pred/meas at the
            // brightest predicted channel). One-directional: only over-prediction is
            // cut. Boundary/truncated bundles are exempt (measured underestimates
            // there). OFF by default (huge ratios = inert); SBND sets the data-tuned
            // ceilings (sbnd_xin/ql_prefilter_tune.py).
            if (m_reject_overpred &&
                !(bundle->get_flag_close_to_PMT() || bundle->get_flag_window_truncated() ||
                  bundle->get_flag_at_x_boundary())) {
                const auto& mask = bundle->get_opdet_mask();
                // Channel scope (see the knob doc in the header): non-empty
                // overpred_channels restricts both ratios to that family (still
                // requiring the per-bundle mask); empty => legacy full masked set.
                double tot_pred = 0.0, tot_meas = 0.0, max_pred = 0.0, meas_at_max = 0.0;
                for (std::size_t j = 0; j < pred_flash.size(); ++j) {
                    if (j >= mask.size() || mask[j] == 0) continue;
                    if (!m_overpred_ch_sel.empty() &&
                        (j >= m_overpred_ch_sel.size() || !m_overpred_ch_sel[j])) continue;
                    const double p = pred_flash[j];
                    const double m = flash->get_PE(static_cast<int>(j));
                    tot_pred += p;
                    tot_meas += m;
                    if (p > max_pred) { max_pred = p; meas_at_max = m; }
                }
                const double r_total = (tot_meas > 0.0) ? tot_pred / tot_meas
                                                        : (tot_pred > 0.0 ? 1e30 : 0.0);
                const double r_max = (max_pred > 0.0)
                                         ? max_pred / (meas_at_max > 1.0 ? meas_at_max : 1.0)
                                         : 0.0;
                if (r_total > m_overpred_total_ratio || r_max > m_overpred_maxch_ratio) {
                    bundle->set_potential_bad_match_flag(true);
                    continue;
                }
            }

            log->debug("initial eval: flash {} and cluster {}, meas PE {}, pred PE {}, npts {}, "
                       "ks_dis {}, chi2/ndf {}, ndf {}",
                       flash->get_flash_id(),
                       cidx,
                       int(flash->get_total_PE() * 100) / 100.,
                       int(bundle->get_total_pred_light() * 100) / 100.,
                       npt,
                       int(bundle->get_ks_dis() * 1000) / 1000.,
                       int(bundle->get_chi2() / bundle->get_ndf() * 100) / 100.,
                       bundle->get_ndf());

            // Opaque-cathode mismatched-candidate filter (legacy post-loop form): drop a
            // bundle whose flash is lit on the opposite drift side from the cluster
            // unless the cluster is a cathode-crosser at this flash's T0 (at_x_boundary).
            // OFF by default => no-op. When crossside_skip_vis is on, the identical drop
            // already fired above (before the vis loop) so this is a no-op for those
            // bundles; same-side bundles fall through here either way. dump_calib applies
            // the same test so the hand-scan candidate tables match this fit candidate set.
            if (cross_side_mismatch_drop(flash.get(), (int)run.anode->ident(),
                                         bundle->get_flag_at_x_boundary())) {
                bundle->set_potential_bad_match_flag(true);
                continue;
            }

            run.pre_bundles.insert(bundle);
        } // cluster loop
    }     // flash loop
    log->debug("n preselected bundles: {}", run.pre_bundles.size());
    log->debug("QLtiming build_bundles: ident {} nflash {} ngroups {} nbundles {} "
               "vis_loop {:.1f} ms over {} candidate-points",
               run.charge_ident, run.flashes.size(), run.match_groups.size(),
               run.all_bundles.size(), vis_ms, total_pts);
}

// Build the flash-keyed / cluster-keyed / (flash,cluster)-keyed bundle maps and
// impose a deterministic iteration order on them (pointer-keyed std::map/std::set
// would otherwise sort by heap address and permute the LASSO matrix columns).
void QLMatching::build_bundle_maps(ApaRun& run)
{
    for (auto bundle : run.pre_bundles) {
        auto flash   = bundle->get_flash();
        auto cluster = bundle->get_main_cluster();
        run.flash_bundles_map[flash].push_back(bundle);
        run.cluster_bundles_map[cluster].push_back(bundle);
        run.flash_cluster_bundles_map[std::make_pair(flash, cluster)] = bundle;
    }

    // Inner order within a flash: bundles by cluster_index_id (the stable global
    // index). Inner order within a cluster: bundles by flash_index_id.
    auto sort_inner_by_cluster_idx = [](FlashBundlesMap& m) {
        for (auto& kv : m) {
            std::sort(kv.second.begin(), kv.second.end(),
                      [](const TimingTPCBundle::pointer& a,
                         const TimingTPCBundle::pointer& b) {
                          return a->get_cluster_index_id() < b->get_cluster_index_id();
                      });
        }
    };
    auto sort_inner_by_flash_idx = [](ClusterBundlesMap& m) {
        for (auto& kv : m) {
            std::sort(kv.second.begin(), kv.second.end(),
                      [](const TimingTPCBundle::pointer& a,
                         const TimingTPCBundle::pointer& b) {
                          return a->get_flash_index_id() < b->get_flash_index_id();
                      });
        }
    };
    sort_inner_by_cluster_idx(run.flash_bundles_map);
    sort_inner_by_flash_idx(run.cluster_bundles_map);
}

// Pin exemption from the strength-cutoff prune. Legacy (min_strength 0):
// a pinned bundle is always exempt.  With the doc-19 floor, a pinned bundle
// whose LASSO solution is at/below the floor loses the exemption (the scan
// showed phantom pins carry strength ~0 while agreed pins sit at p10 0.88).
bool QLMatching::pin_exempt(const TimingTPCBundle* b, double strength) const
{
    if (!b->get_flag_xtpc_pin()) return false;
    if (m_xtpc_pin_min_strength <= 0.0) return true;   // legacy always-exempt
    return strength > m_xtpc_pin_min_strength;
}

// Post-fit cull of unflagged low-quality selections (doc 19). A bundle still
// in flash_bundles_map after the rounds/rescues is a matched output; one that
// carries NO quality flag survived on LASSO strength alone -- the largest
// phantom bucket in the 18-evt scan. Remove it when its light disagrees
// (ks > postcull_ks_max OR chi2/ndf > postcull_c2n_max). Purely a filter,
// iteration-order independent; OFF (default) => no-op, byte-identical.
void QLMatching::cull_unflagged_lowquality(ApaRun& run)
{
    if (!m_postcull_unflagged && !m_postcull_wtrunc_overpred && !m_postcull_pin_overpred) return;
    TimingTPCBundleSelection to_be_removed;
    for (auto* flash : flash_iter_order(run.flash_bundles_map)) {
        for (const auto& b : run.flash_bundles_map.at(flash)) {
            const bool qflag = b->get_consistent_flag() || b->get_flag_xtpc_consistent() ||
                               b->get_flag_xtpc_scenario1() || b->get_flag_xtpc_pin();
            if (m_postcull_unflagged && !qflag) {
                const double c2n = b->get_ndf() > 0 ? b->get_chi2() / b->get_ndf() : 1e9;
                if (b->get_ks_dis() > m_postcull_ks_max || c2n > m_postcull_c2n_max) {
                    to_be_removed.push_back(b);
                    log->debug("QLPOSTCULL apa{} cluster {} flash {} ks={:.3f} chi2/ndf={:.1f} "
                               "(unflagged, fails {:.2f}/{:.1f})",
                               (int)run.anode->ident(), b->get_main_cluster()->ident(),
                               b->get_flash()->get_flash_id(), b->get_ks_dis(), c2n,
                               m_postcull_ks_max, m_postcull_c2n_max);
                    continue;
                }
            }
            // Window-truncated overprediction cull (doc 23 phase 2; default OFF).
            // A selection flagged window_truncated whose total pred/meas exceeds
            // the ceiling is the dominant clean phantom signature (039252: ratio
            // p90 11.3 vs 1.5 for scan-agreed wtrunc matches). xtpc pin/scenario-1
            // pairs are protected (geometric evidence overrides), and a
            // sat-dominated flash is exempt (railed meas is a lower bound, its
            // ratio is an overestimate -- same physics as the phase-1b rescue
            // extension). Applies regardless of the high-consistent flag.
            if (m_postcull_wtrunc_overpred && b->get_flag_window_truncated() &&
                !b->get_flag_xtpc_pin() && !b->get_flag_xtpc_scenario1()) {
                auto* f = b->get_flash();
                const double meas = f->get_total_PE();
                const double ratio = (meas > 0.0) ? b->get_total_pred_light() / meas : 1e9;
                if (ratio > m_postcull_wtrunc_ratio_hi) {
                    double satpe = 0;
                    const int nch = f->get_num_channels();
                    for (int j = 0; j < nch; ++j)
                        if (f->get_sat(j)) satpe += f->get_PE(j);
                    const double satfrac = (meas > 0.0) ? satpe / meas : 0.0;
                    if (satfrac <= m_postcull_wtrunc_sat_frac) {
                        to_be_removed.push_back(b);
                        log->debug("QLPOSTCULL-WTRUNC apa{} cluster {} flash {} ratio={:.2f} "
                                   "(window-truncated overpred, ceiling {:.1f}, satfrac {:.2f})",
                                   (int)run.anode->ident(), b->get_main_cluster()->ident(),
                                   f->get_flash_id(), ratio, m_postcull_wtrunc_ratio_hi, satfrac);
                        continue;
                    }
                }
            }
            // xtpc-pin overprediction cull (doc 23 phase 2): ratio-ONLY (a ks
            // gate on pins kills legitimate geometric matches); same
            // sat-dominated exemption as the wtrunc branch. The pin's partner
            // half keeps its own flag -- this removes only the overpredicting
            // half's selection.
            if (m_postcull_pin_overpred && b->get_flag_xtpc_pin()) {
                auto* f = b->get_flash();
                const double meas = f->get_total_PE();
                const double ratio = (meas > 0.0) ? b->get_total_pred_light() / meas : 1e9;
                if (ratio > m_postcull_pin_ratio_hi) {
                    double satpe = 0;
                    const int nch = f->get_num_channels();
                    for (int j = 0; j < nch; ++j)
                        if (f->get_sat(j)) satpe += f->get_PE(j);
                    const double satfrac = (meas > 0.0) ? satpe / meas : 0.0;
                    if (satfrac <= m_postcull_wtrunc_sat_frac) {
                        to_be_removed.push_back(b);
                        log->debug("QLPOSTCULL-PIN apa{} cluster {} flash {} ratio={:.2f} "
                                   "(pin overpred, ceiling {:.1f}, satfrac {:.2f})",
                                   (int)run.anode->ident(), b->get_main_cluster()->ident(),
                                   f->get_flash_id(), ratio, m_postcull_pin_ratio_hi, satfrac);
                    }
                }
            }
        }
    }
    if (to_be_removed.empty()) return;
    remove_bundle_selection(to_be_removed, run.flash_bundles_map, run.cluster_bundles_map,
                            run.flash_cluster_bundles_map);
    remove_bundle_selection(to_be_removed, run.pre_bundles);
}

// xtpc cathode rescue resolution (see the operator() call site and QLMatching.h).
// Removes provisional cathode-overshoot bundles that did not acquire
// flag_xtpc_scenario1 in cull_cross_tpc; stamps survivors contained (accepted with
// the overshoot tolerance) so dump_calib shows them. Purely a filter: result is
// independent of iteration order.
void QLMatching::purge_unconfirmed_cathode_rescue(ApaRun& run)
{
    TimingTPCBundleSelection to_be_removed;
    for (auto* flash : flash_iter_order(run.flash_bundles_map)) {
        for (const auto& b : run.flash_bundles_map.at(flash)) {
            if (!b->get_flag_xtpc_cathode_provisional()) continue;
            // Optional doc-19 ks ceiling: a scenario-1-confirmed but light-dim
            // provisional bundle is purged too (0 = off, legacy).
            if (b->get_flag_xtpc_scenario1() &&
                (m_xtpc_cathode_ks_max <= 0.0 || b->get_ks_dis() <= m_xtpc_cathode_ks_max)) {
                b->set_contained(true);
                log->debug("QLXTPC cathode-rescue KEEP apa{} cluster {} flash {} "
                           "(scenario-1 confirmed crosser half)",
                           (int)run.anode->ident(), b->get_main_cluster()->ident(),
                           b->get_flash()->get_flash_id());
            }
            else {
                to_be_removed.push_back(b);
                log->debug("QLXTPC cathode-rescue DROP apa{} cluster {} flash {} "
                           "(no cross-volume confirmation)",
                           (int)run.anode->ident(), b->get_main_cluster()->ident(),
                           b->get_flash()->get_flash_id());
            }
        }
    }
    if (to_be_removed.empty()) return;
    remove_bundle_selection(to_be_removed, run.flash_bundles_map, run.cluster_bundles_map,
                            run.flash_cluster_bundles_map);
    remove_bundle_selection(to_be_removed, run.pre_bundles);
}

// [Stage 1] Per-cluster consistency cull, cross-TPC-aware (runs AFTER cull_cross_tpc's
// flag pass). For each cluster:
//   - has a SCENARIO-1 cross-TPC bundle -> keep the xtpc bundle(s), drop ALL others
//     (incl. high-consistent bundles on other flashes). A tight cathode-crossing match
//     (closest approach < dmax, self-vetoing) overrides an accidental high-consistent
//     match on the wrong flash -- this is the steering fix.
//   - else has a high-consistent OR (scenario-2) xtpc bundle -> keep those, drop the
//     rest. Reproduces the historical cull_inconsistent + the old cull_cross_tpc weak
//     rival cull combined; scenario-2 xtpc keeps the weaker (compete-in-LASSO) behaviour.
//   - else -> keep all.
// When no bundle is xtpc-flagged (single-APA / xtpc_flag:false) only the middle branch
// fires with high-consistent bundles == today's behaviour => BIT-IDENTICAL.
void QLMatching::cull_inconsistent(ApaRun& run)
{
    auto keep_flag = [](const TimingTPCBundle::pointer& b) {
        return b->get_consistent_flag() || b->get_flag_xtpc_consistent();
    };
    TimingTPCBundleSelection to_be_removed;
    for (auto* cluster : cluster_iter_order(run.cluster_bundles_map, run.global_cluster_idx_map)) {
        auto& bundles = run.cluster_bundles_map[cluster];
        bool has_pin = false, has_sc1 = false, has_keep = false;
        for (auto& b : bundles) {
            if (b->get_flag_xtpc_pin()) has_pin = true;
            if (b->get_flag_xtpc_scenario1()) has_sc1 = true;
            if (keep_flag(b)) has_keep = true;
        }
        if (has_pin) {
            // Joint-pin (xtpc_joint_pin only): a direction-confirmed crosser half keeps ONLY
            // its pinned bundle, binding the pair to one coincident flash. Overrides even the
            // scenario-1 priority below (which keeps all sc1 bundles, allowing a split).
            for (auto& b : bundles)
                if (!b->get_flag_xtpc_pin()) {
                    to_be_removed.push_back(b);
                    log->debug("QLCULLINC apa{} cluster {} drop bundle flash {} "
                               "(cluster kept xtpc JOINT-PIN crosser)",
                               (int)run.anode->ident(), cluster->ident(),
                               b->get_flash()->get_flash_id());
                }
        }
        else if (has_sc1) {
            for (auto& b : bundles)
                if (!b->get_flag_xtpc_scenario1()) {
                    to_be_removed.push_back(b);
                    log->debug("QLCULLINC apa{} cluster {} drop bundle flash {} "
                               "(cluster kept xtpc scenario-1 crosser)",
                               (int)run.anode->ident(), cluster->ident(),
                               b->get_flash()->get_flash_id());
                }
        }
        else if (has_keep) {
            for (auto& b : bundles)
                if (!keep_flag(b)) {
                    to_be_removed.push_back(b);
                    log->debug("QLCULLINC apa{} cluster {} drop bundle flash {} "
                               "(cluster kept high/xtpc-consistent bundle)",
                               (int)run.anode->ident(), cluster->ident(),
                               b->get_flash()->get_flash_id());
                }
        }
        // else: no consistent/xtpc bundle -> keep all
    }
    remove_bundle_selection(to_be_removed, run.flash_bundles_map, run.cluster_bundles_map,
                            run.flash_cluster_bundles_map);
    remove_bundle_selection(to_be_removed, run.pre_bundles);
}

// Prototype flag-based per-column L1 down-weight (ToyMatching.cxx:1437), generalized to
// the same boundary/near-PMT/window-truncated group as the ladder's B4 branch. A
// down-weighted bundle is shrunk less and so more likely to survive the strength cutoff.
double QLMatching::lasso_flag_factor(const TimingTPCBundle::pointer& bundle) const
{
    if (m_lasso_flag_weight &&
        (bundle->get_flag_at_x_boundary() || bundle->get_flag_close_to_PMT()
         || bundle->get_flag_window_truncated()))
        return m_lasso_boundary_weight;
    return 1.0;
}

// [Stage 2] First LASSO round, with a per-flash background/light DOF column;
// post-fit prune bundles whose strength <= m_strength_cutoff.
void QLMatching::fit_round1(ApaRun& run)
{
    const auto t_build0 = wallclock::now();
    // Snapshot the full pre-LASSO candidate universe for the empty-flash (§I) and
    // cluster-centric (§J) rescues, before either round prunes by strength. Capturing
    // it is unobservable unless a rescue consumes it, so this stays byte-identical for
    // configs with both rescues off.
    if (m_empty_rescue || m_cluster_rescue) run.prefit_snapshot = run.flash_bundles_map;

    const double lambda       = m_lasso_lambda;
    const double delta_charge = m_delta_charge;
    const double delta_light  = m_delta_light;

    const unsigned int nbundle  = run.pre_bundles.size();
    const unsigned int nflash   = run.flash_bundles_map.size();
    const unsigned int ncluster = run.cluster_bundles_map.size();

    auto flashes_ordered  = flash_iter_order(run.flash_bundles_map);
    auto clusters_ordered = cluster_iter_order(run.cluster_bundles_map, run.global_cluster_idx_map);
    std::map<Opflash*, int> flash_idx_map;
    std::map<Cluster*, int> cluster_idx_map;
    int idx = 0;
    for (auto* c : clusters_ordered) cluster_idx_map[c] = idx++;
    idx = 0;
    for (auto* f : flashes_ordered)  flash_idx_map[f]   = idx++;

    Ress::vector_t M = Ress::vector_t::Zero(run.nopdet * nflash);
    Ress::vector_t MF = Ress::vector_t::Zero(ncluster + nflash);
    Ress::vector_t weights = Ress::vector_t::Zero(nbundle + nflash);
    // P (nopdet*nflash x nbundle+nflash) and PF (ncluster+nflash x nbundle+nflash)
    // are block-sparse; collect their entries as triplets and realise them sparse
    // (m_sparse_lasso) or dense (default) below.  The (row,col,value) set is the
    // same either way, so the dense realisation is byte-identical to the historical
    // dense fill.
    std::vector<Eigen::Triplet<double>> P_trip, PF_trip;

    std::vector<std::pair<Opflash*, Cluster*>> pairs;
    std::size_t i = 0, ik = 0;
    for (auto* flash : flashes_ordered) {
        auto& bundles = run.flash_bundles_map[flash];
        for (unsigned int j = 0; j < run.nopdet; ++j) {
            const int opdet_idx = run.opdet_idx_v.at(j);
            // Saturated-channel row of a flagged flash: leave the whole LASSO
            // row zero (measured, background and pred entries) -- the clipped
            // PE must not constrain the fit.  Default off, bit-identical.
            if (m_use_saturation_flag && flash->get_sat(opdet_idx)) continue;
            if (m_use_coverage_flag && m_coverage_mask_fit && flash->get_cov(opdet_idx) < m_coverage_min) continue;
            const double pe = flash->get_PE(opdet_idx);
            const double pe_err = std::sqrt(flash->get_PE(opdet_idx) + std::pow(flash->get_PE_err(opdet_idx), 2));
            M(i * run.nopdet + j) = pe / pe_err;
            P_trip.emplace_back((int)(i * run.nopdet + j), (int)(nbundle + i), pe / pe_err);
        }
        for (auto bundle : bundles) {
            const auto& pred_flash = bundle->get_pred_flash();
            for (unsigned int j = 0; j < run.nopdet; ++j) {
                const int opdet_idx = run.opdet_idx_v.at(j);
                if (m_use_saturation_flag && flash->get_sat(opdet_idx)) continue;
                if (m_use_coverage_flag && m_coverage_mask_fit && flash->get_cov(opdet_idx) < m_coverage_min) continue;
                const double pred_pe = pred_flash.at(opdet_idx);
                const double pe_err = std::sqrt(flash->get_PE(opdet_idx) + std::pow(flash->get_PE_err(opdet_idx), 2));
                P_trip.emplace_back((int)(i * run.nopdet + j), (int)pairs.size(), pred_pe / pe_err);
            }
            pairs.emplace_back(flash, bundle->get_main_cluster());
            const auto meas_pe_tot = flash->get_total_PE();
            const auto pred_pe_tot = bundle->get_total_pred_light();
            const double base = (std::abs(pred_pe_tot - meas_pe_tot) > m_pe_mismatch_knee * meas_pe_tot)
                                    ? std::abs(pred_pe_tot - meas_pe_tot) / meas_pe_tot
                                    : m_pe_mismatch_floor;
            weights(ik++) = base * lasso_flag_factor(bundle);
        }
        PF_trip.emplace_back((int)(ncluster + i), (int)(nbundle + i), 1. / delta_light);
        flash_idx_map[flash] = nbundle + i;
        ++i;
    }
    for (unsigned int k = 0; k < nflash; ++k) weights(nbundle + k) = m_bkg_weight;
    for (unsigned int k = 0; k < ncluster; ++k) MF(k) = 1. / delta_charge;
    for (std::size_t n = 0; n < pairs.size(); ++n) {
        PF_trip.emplace_back(cluster_idx_map[pairs.at(n).second], (int)n, 1. / delta_charge);
    }

    // Realise P/PF and form the normal equations.  Default: dense (byte-identical
    // to the historical Ress path).  m_sparse_lasso: sparse products -- removes the
    // dense P/PT spike and the O(nbeta²·nopdet·nflash) build on busy events (FP
    // accumulation order differs, so X may move at the ULP level; assignments do not).
    Eigen::SparseMatrix<double> P_sp((int)(run.nopdet * nflash), (int)(nbundle + nflash));
    Eigen::SparseMatrix<double> PF_sp((int)(ncluster + nflash), (int)(nbundle + nflash));
    P_sp.setFromTriplets(P_trip.begin(), P_trip.end());
    PF_sp.setFromTriplets(PF_trip.begin(), PF_trip.end());
    Ress::vector_t initial = Ress::vector_t::Zero(nbundle + nflash);
    for (std::size_t n = 0; n < pairs.size(); ++n) initial(n) = 1.0;

    Ress::Params params;
    params.model = Ress::lasso;
    params.lambda = lambda;
    log->debug("solving (round 1)");
    const int xdim = (int)(nbundle + nflash);
    Ress::vector_t solution;
    if (m_sparse_lasso) {
        // Sparse normal equations fed straight to the solver (B2): no dense X, and the
        // solver forms XᵀX / Xᵀy by sparse products -- skipping the ~98% zero work.
        Eigen::SparseMatrix<double> PT  = P_sp.transpose();
        Eigen::SparseMatrix<double> PFT = PF_sp.transpose();
        Ress::vector_t y = PT * M + PFT * MF;
        Eigen::SparseMatrix<double> Xs = PT * P_sp + PFT * PF_sp;
        const double t_build = ms_since(t_build0);
        const auto t_solve0 = wallclock::now();
        solution = Ress::solve(Xs, y, params, initial, weights);
        log->debug("QLtiming fit_round1: ident {} nbundle {} nflash {} nopdet {} matrix_build {:.1f} "
                   "lasso_solve {:.1f} ms (X {}x{})",
                   run.charge_ident, nbundle, nflash, run.nopdet, t_build, ms_since(t_solve0),
                   xdim, xdim);
    }
    else {
        Ress::matrix_t P  = P_sp.toDense();
        Ress::matrix_t PF = PF_sp.toDense();
        Ress::matrix_t PT  = P.transpose();
        Ress::matrix_t PFT = PF.transpose();
        Ress::vector_t y = PT * M + PFT * MF;
        Ress::matrix_t X = PT * P + PFT * PF;
        const double t_build = ms_since(t_build0);
        const auto t_solve0 = wallclock::now();
        solution = Ress::solve(X, y, params, initial, weights);
        log->debug("QLtiming fit_round1: ident {} nbundle {} nflash {} nopdet {} matrix_build {:.1f} "
                   "lasso_solve {:.1f} ms (X {}x{})",
                   run.charge_ident, nbundle, nflash, run.nopdet, t_build, ms_since(t_solve0),
                   (int)X.rows(), (int)X.cols());
    }

    TimingTPCBundleSelection to_be_removed;
    int n = 0;
    for (auto* flash : flashes_ordered) {
        for (auto bundle : run.flash_bundles_map[flash]) {
            // xtpc joint-pin: a pinned crosser half keeps its flash regardless of LASSO
            // strength (the pin ignores light by design), so it survives the cutoff prune
            // and stays in flash_bundles_map -> apply_matched_t0s matches it. Off-path the
            // flag is never set => identical.  pin_exempt applies the optional
            // xtpc_pin_min_strength floor (doc 19); 0 => legacy always-exempt.
            if (solution(n) <= m_strength_cutoff && !m_beamonly &&
                !pin_exempt(bundle.get(), solution(n)))
                to_be_removed.push_back(bundle);
            ++n;
        }
    }
    remove_bundle_selection(to_be_removed, run.flash_bundles_map, run.cluster_bundles_map,
                            run.flash_cluster_bundles_map);
    remove_bundle_selection(to_be_removed, run.pre_bundles);
}

// [Stage 3] Second LASSO round: background DOF dropped, KS-shape term added to
// the weights; prune by strength and keep the best flash per cluster, then merge
// and apply out-of-beam QA via organize_bundles().
void QLMatching::fit_round2(ApaRun& run)
{
    const auto t_build0 = wallclock::now();
    const double lambda       = m_lasso_lambda;
    const double delta_charge = m_delta_charge;
    const double delta_shape  = m_delta_shape;

    const unsigned int nbundle  = run.pre_bundles.size();
    const unsigned int nflash   = run.flash_bundles_map.size();
    const unsigned int ncluster = run.cluster_bundles_map.size();

    // Rebuild ordered iteration (rd-1 may have removed bundles/flashes/clusters).
    auto flashes_ordered  = flash_iter_order(run.flash_bundles_map);
    auto clusters_ordered = cluster_iter_order(run.cluster_bundles_map, run.global_cluster_idx_map);
    std::map<Cluster*, int> cluster_idx_map;
    std::map<Opflash*, int> flash_idx_map;
    int idx = 0;
    for (auto* c : clusters_ordered) cluster_idx_map[c] = idx++;
    idx = 0;
    for (auto* f : flashes_ordered)  flash_idx_map[f]   = idx++;

    Ress::vector_t M = Ress::vector_t::Zero(run.nopdet * nflash);
    Ress::vector_t MF = Ress::vector_t::Zero(ncluster);
    Ress::vector_t weights = Ress::vector_t::Zero(nbundle);
    // Block-sparse P (nopdet*nflash x nbundle) / PF (ncluster x nbundle) as triplets;
    // realised dense (default, byte-identical) or sparse (m_sparse_lasso) below.
    std::vector<Eigen::Triplet<double>> P_trip, PF_trip;
    std::vector<std::pair<Opflash*, Cluster*>> pairs;

    std::size_t i = 0, ik = 0;
    for (auto* flash : flashes_ordered) {
        auto& bundles = run.flash_bundles_map[flash];
        for (unsigned int j = 0; j < run.nopdet; ++j) {
            const int opdet_idx = run.opdet_idx_v.at(j);
            // Saturated-channel rows stay zero (see fit_round1).
            if (m_use_saturation_flag && flash->get_sat(opdet_idx)) continue;
            if (m_use_coverage_flag && m_coverage_mask_fit && flash->get_cov(opdet_idx) < m_coverage_min) continue;
            const double pe = flash->get_PE(opdet_idx);
            const double pe_err = std::sqrt(flash->get_PE(opdet_idx) + std::pow(flash->get_PE_err(opdet_idx), 2));
            M(i * run.nopdet + j) = pe / pe_err;
        }
        for (auto bundle : bundles) {
            const auto ks_dis = bundle->get_ks_dis();
            const auto& pred_flash = bundle->get_pred_flash();
            for (unsigned int j = 0; j < run.nopdet; ++j) {
                const int opdet_idx = run.opdet_idx_v.at(j);
                if (m_use_saturation_flag && flash->get_sat(opdet_idx)) continue;
                if (m_use_coverage_flag && m_coverage_mask_fit && flash->get_cov(opdet_idx) < m_coverage_min) continue;
                const double pred_pe = pred_flash.at(opdet_idx);
                const double pe_err = std::sqrt(flash->get_PE(opdet_idx) + std::pow(flash->get_PE_err(opdet_idx), 2));
                P_trip.emplace_back((int)(i * run.nopdet + j), (int)pairs.size(), pred_pe / pe_err);
            }
            pairs.emplace_back(flash, bundle->get_main_cluster());

            const auto meas_pe_tot = flash->get_total_PE();
            const auto pred_pe_tot = bundle->get_total_pred_light();
            const double base = (std::abs(pred_pe_tot - meas_pe_tot) > m_pe_mismatch_knee * meas_pe_tot)
                                    ? std::abs(pred_pe_tot - meas_pe_tot) / meas_pe_tot
                                    : m_pe_mismatch_floor;
            weights(ik++) = (base + delta_shape * run.nopdet * ks_dis / lambda) * lasso_flag_factor(bundle);
        }
        ++i;
    }
    for (unsigned int k = 0; k < ncluster; ++k) MF(k) = 1. / delta_charge;
    for (std::size_t n = 0; n < pairs.size(); ++n) {
        PF_trip.emplace_back(cluster_idx_map[pairs.at(n).second], (int)n, 1. / delta_charge);
    }

    Eigen::SparseMatrix<double> P_sp((int)(run.nopdet * nflash), (int)nbundle);
    Eigen::SparseMatrix<double> PF_sp((int)ncluster, (int)nbundle);
    P_sp.setFromTriplets(P_trip.begin(), P_trip.end());
    PF_sp.setFromTriplets(PF_trip.begin(), PF_trip.end());
    Ress::vector_t initial = Ress::vector_t::Zero(nbundle);
    for (std::size_t n = 0; n < pairs.size(); ++n) initial(n) = 1.0;

    Ress::Params params;
    params.model = Ress::lasso;
    params.lambda = lambda;
    log->debug("solving (round 2)");
    const int xdim = (int)nbundle;
    Ress::vector_t solution;
    if (m_sparse_lasso) {
        // Sparse normal equations straight to the solver (B2): no dense X.
        Eigen::SparseMatrix<double> PT  = P_sp.transpose();
        Eigen::SparseMatrix<double> PFT = PF_sp.transpose();
        Ress::vector_t y = PT * M + PFT * MF;
        Eigen::SparseMatrix<double> Xs = PT * P_sp + PFT * PF_sp;
        const double t_build = ms_since(t_build0);
        const auto t_solve0 = wallclock::now();
        solution = Ress::solve(Xs, y, params, initial, weights);
        log->debug("QLtiming fit_round2: ident {} nbundle {} nflash {} nopdet {} matrix_build {:.1f} "
                   "lasso_solve {:.1f} ms (X {}x{})",
                   run.charge_ident, nbundle, nflash, run.nopdet, t_build, ms_since(t_solve0),
                   xdim, xdim);
    }
    else {
        Ress::matrix_t P  = P_sp.toDense();
        Ress::matrix_t PF = PF_sp.toDense();
        Ress::matrix_t PT  = P.transpose();
        Ress::matrix_t PFT = PF.transpose();
        Ress::vector_t y = PT * M + PFT * MF;
        Ress::matrix_t X = PT * P + PFT * PF;
        const double t_build = ms_since(t_build0);
        const auto t_solve0 = wallclock::now();
        solution = Ress::solve(X, y, params, initial, weights);
        log->debug("QLtiming fit_round2: ident {} nbundle {} nflash {} nopdet {} matrix_build {:.1f} "
                   "lasso_solve {:.1f} ms (X {}x{})",
                   run.charge_ident, nbundle, nflash, run.nopdet, t_build, ms_since(t_solve0),
                   (int)X.rows(), (int)X.cols());
    }

    TimingTPCBundleSelection to_be_removed;
    int n = 0;
    for (auto* flash : flashes_ordered) {
        for (auto bundle : run.flash_bundles_map[flash]) {
            bundle->set_strength(solution(n));
            // xtpc joint-pin: keep a pinned crosser half regardless of strength (see round1).
            if (!(solution(n) > m_strength_cutoff || m_beamonly ||
                  pin_exempt(bundle.get(), solution(n))))
                to_be_removed.push_back(bundle);
            ++n;
        }
    }
    remove_bundle_selection(to_be_removed, run.flash_bundles_map, run.cluster_bundles_map,
                            run.flash_cluster_bundles_map);
    remove_bundle_selection(to_be_removed, run.pre_bundles);
    to_be_removed.clear();

    // Rescue blind-spot fix (doc 23 phase 1a; default OFF = bit-identical): cull the
    // postcull-doomed unflagged selections BEFORE the rescues, so a bundle destined
    // for removal cannot mark its cluster "matched" and hide it from §I/§J. The
    // legacy post-rescue call below still runs (rescue adoptions stay subject to
    // the same quality bar as before).
    if (m_postcull_before_rescue) cull_unflagged_lowquality(run);

    // Empty-flash light-quality rescue (§I; default OFF = bit-identical, SBND-on).
    // Uses the pre-LASSO snapshot captured at fit_round1 start, so it can reach the
    // strength-0 but light-good bundles both rounds pruned.
    if (m_empty_rescue) rescue_empty_flashes(run, run.prefit_snapshot);

    // Cluster-centric rescue (§J; default OFF = bit-identical, PDHD-on). After the
    // empty-flash pass, adopt the best accepted candidate for each cluster the LASSO
    // left unmatched (its flash won by a rival), attaching it even onto a non-empty
    // flash. Same pre-LASSO snapshot source as the empty-flash rescue.
    if (m_cluster_rescue) rescue_unmatched_clusters(run, run.prefit_snapshot);

    // Post-fit unflagged low-quality cull (doc 19; OFF default = bit-identical).
    cull_unflagged_lowquality(run);

    // Keep best match per cluster.
    std::map<int, std::pair<Opflash*, double>> matched_pairs;
    for (std::size_t k = 0; k < pairs.size(); ++k) {
        if (solution(k) <= m_strength_cutoff) continue;
        const int cidx = cluster_idx_map[pairs.at(k).second];
        auto flash = pairs.at(k).first;
        auto it = matched_pairs.find(cidx);
        if (it == matched_pairs.end() || solution(k) > it->second.second) {
            matched_pairs[cidx] = std::make_pair(flash, solution(k));
        }
    }
    TimingTPCBundleSelection results_bundles;
    for (auto cluster : run.clusters) {
        auto it = cluster_idx_map.find(cluster);
        if (it == cluster_idx_map.end()) continue;
        const int cidx = it->second;
        auto mit = matched_pairs.find(cidx);
        if (mit != matched_pairs.end()) {
            auto flash = mit->second.first;
            // Deep-copy the matched bundle. organize_bundles()'s result is discarded
            // (the matched output is the strength-cutoff survivors left in
            // run.flash_bundles_map), but its add_bundle() merge mutates bundles IN
            // PLACE -- and results_bundles shares the very shared_ptr objects that
            // flash_bundles_map (read later by apply_matched_t0s and the calib dump)
            // still references. Sharing them lets the merge fold a rival cluster's
            // predicted light into the surviving bundle while leaving the rival in
            // flash_bundles_map, double-counting that cluster's light on its flash.
            // Operating on copies isolates the (discarded) organize pass from the
            // live result.
            results_bundles.push_back(std::make_shared<TimingTPCBundle>(
                *run.flash_cluster_bundles_map[std::make_pair(flash, cluster)]));
        }
        else {
            auto bundle = std::make_shared<TimingTPCBundle>(nullptr, cluster, 0, cidx);
            bundle->set_quality_params(run.qp);
            bundle->set_strength(0);
            results_bundles.push_back(bundle);
        }
    }
    organize_bundles(results_bundles, run.flash_cluster_bundles_map);

    FlashBundlesMap results_flash_bundles_map;
    for (auto bundle : results_bundles) {
        results_flash_bundles_map[bundle->get_flash()].push_back(bundle);
    }

    log->debug("done with matching");
}

// ---------------------------------------------------------------------------
// Shared-flash JOINT LASSO rounds (m_shared_flash; PDVD). Every input port reads
// the SAME global flash list (identical filters), so a physical flash is keyed
// by its input tensor row id (get_flash_id) and appears as one per-run Opflash
// instance per port. The joint system has ONE measured block per physical flash
// and bundle columns from EVERY port, so clusters in both drift volumes share
// the explanation of one flash (a per-run fit would let each side independently
// absorb the full PE — biased for every cathode-crosser). The matrix recipe
// (rows, weights, background columns, sparse/dense realisation, strength prune)
// mirrors fit_round1/fit_round2 verbatim; only the column universe is joint.
// fit_round2's matched_pairs/organize_bundles tail is NOT replicated: its result
// is discarded there too (the matched output is the strength-cutoff survivors
// left in flash_bundles_map). The §I/§J rescues are per-run "empty flash"
// concepts, wrong under a shared flash — skipped (configure warns).
namespace {
    // One physical flash: per-run (run index, per-run Opflash instance) pairs,
    // in input-port order.
    struct SharedFlash {
        std::vector<std::pair<std::size_t, WireCell::Match::Opflash*>> insts;
    };
}

void QLMatching::fit_round1_shared(std::vector<ApaRun>& runs)
{
    const auto t_build0 = wallclock::now();
    const double lambda       = m_lasso_lambda;
    const double delta_charge = m_delta_charge;
    const double delta_light  = m_delta_light;

    // Pre-LASSO snapshot for the shared rescues (doc 19 phase 5), mirroring
    // fit_round1's capture; gated so the OFF path allocates nothing new.
    if (m_empty_rescue_shared || m_cluster_rescue_shared) {
        for (auto& run : runs) run.prefit_snapshot = run.flash_bundles_map;
    }

    // OpDet universe: identical across runs by construction (opdet_all_volumes +
    // deterministic auto-mask on identical flash lists); guard config mistakes.
    const ApaRun& r0 = runs.front();
    for (const auto& run : runs) {
        if (run.opdet_idx_v != r0.opdet_idx_v) {
            raise<ValueError>("QLMatching: shared_flash requires identical OpDet masks "
                              "across inputs (set opdet_all_volumes)");
        }
    }
    const unsigned int nopdet = r0.nopdet;

    // Physical flashes = union of the runs' candidate flashes, keyed by flash row
    // id (ascending => same ordering semantics as flash_iter_order).
    std::map<int, SharedFlash> phys;
    for (std::size_t k = 0; k < runs.size(); ++k) {
        for (auto* f : flash_iter_order(runs[k].flash_bundles_map)) {
            phys[f->get_flash_id()].insts.emplace_back(k, f);
        }
    }
    for (auto& [fid, pf] : phys) {
        Opflash* f0 = pf.insts.front().second;
        for (std::size_t n = 1; n < pf.insts.size(); ++n) {
            Opflash* fn = pf.insts[n].second;
            if (std::abs(fn->get_time() - f0->get_time()) > 1e-6 ||
                std::abs(fn->get_total_PE() - f0->get_total_PE()) >
                    1e-9 * std::max(1.0, std::abs(f0->get_total_PE()))) {
                raise<ValueError>("QLMatching: shared_flash flash id %d differs across "
                                  "inputs (time/PE) — ports must read one flash list", fid);
            }
        }
    }

    const unsigned int nflash = phys.size();
    unsigned int nbundle = 0;
    for (auto& [fid, pf] : phys) {
        (void)fid;
        for (auto& [k, f] : pf.insts) nbundle += runs[k].flash_bundles_map.at(f).size();
    }

    // Joint cluster index: run-major, each run's clusters in its own stable order.
    // Cluster facades are distinct objects across runs, so a flat map suffices.
    std::map<Cluster*, int> cluster_idx_map;
    {
        int idx = 0;
        for (auto& run : runs) {
            for (auto* c : cluster_iter_order(run.cluster_bundles_map, run.global_cluster_idx_map)) {
                cluster_idx_map[c] = idx++;
            }
        }
    }
    const unsigned int ncluster = cluster_idx_map.size();

    Ress::vector_t M = Ress::vector_t::Zero(nopdet * nflash);
    Ress::vector_t MF = Ress::vector_t::Zero(ncluster + nflash);
    Ress::vector_t weights = Ress::vector_t::Zero(nbundle + nflash);
    std::vector<Eigen::Triplet<double>> P_trip, PF_trip;

    // Column bookkeeping: the bundle behind each column and its owning run (for
    // the post-solve prune routing).
    std::vector<TimingTPCBundle::pointer> col_bundles;
    std::vector<std::size_t> col_run;
    col_bundles.reserve(nbundle);
    col_run.reserve(nbundle);

    std::size_t i = 0, ik = 0;
    for (auto& [fid, pf] : phys) {
        (void)fid;
        Opflash* f0 = pf.insts.front().second;
        for (unsigned int j = 0; j < nopdet; ++j) {
            const int opdet_idx = r0.opdet_idx_v.at(j);
            // Saturated-channel rows stay zero (see fit_round1).
            if (m_use_saturation_flag && f0->get_sat(opdet_idx)) continue;
            if (m_use_coverage_flag && m_coverage_mask_fit && f0->get_cov(opdet_idx) < m_coverage_min) continue;
            const double pe = f0->get_PE(opdet_idx);
            const double pe_err = std::sqrt(f0->get_PE(opdet_idx) + std::pow(f0->get_PE_err(opdet_idx), 2));
            M(i * nopdet + j) = pe / pe_err;
            P_trip.emplace_back((int)(i * nopdet + j), (int)(nbundle + i), pe / pe_err);
        }
        for (auto& [k, f] : pf.insts) {
            for (auto bundle : runs[k].flash_bundles_map.at(f)) {
                const auto& pred_flash = bundle->get_pred_flash();
                for (unsigned int j = 0; j < nopdet; ++j) {
                    const int opdet_idx = r0.opdet_idx_v.at(j);
                    if (m_use_saturation_flag && f0->get_sat(opdet_idx)) continue;
                    if (m_use_coverage_flag && m_coverage_mask_fit && f0->get_cov(opdet_idx) < m_coverage_min) continue;
                    const double pred_pe = pred_flash.at(opdet_idx);
                    const double pe_err = std::sqrt(f0->get_PE(opdet_idx) + std::pow(f0->get_PE_err(opdet_idx), 2));
                    P_trip.emplace_back((int)(i * nopdet + j), (int)col_bundles.size(), pred_pe / pe_err);
                }
                const auto meas_pe_tot = f0->get_total_PE();
                const auto pred_pe_tot = bundle->get_total_pred_light();
                const double base = (std::abs(pred_pe_tot - meas_pe_tot) > m_pe_mismatch_knee * meas_pe_tot)
                                        ? std::abs(pred_pe_tot - meas_pe_tot) / meas_pe_tot
                                        : m_pe_mismatch_floor;
                weights(ik++) = base * lasso_flag_factor(bundle);
                col_bundles.push_back(bundle);
                col_run.push_back(k);
            }
        }
        PF_trip.emplace_back((int)(ncluster + i), (int)(nbundle + i), 1. / delta_light);
        ++i;
    }
    for (unsigned int k = 0; k < nflash; ++k) weights(nbundle + k) = m_bkg_weight;
    for (unsigned int k = 0; k < ncluster; ++k) MF(k) = 1. / delta_charge;
    for (std::size_t n = 0; n < col_bundles.size(); ++n) {
        PF_trip.emplace_back(cluster_idx_map.at(col_bundles[n]->get_main_cluster()), (int)n, 1. / delta_charge);
    }

    Eigen::SparseMatrix<double> P_sp((int)(nopdet * nflash), (int)(nbundle + nflash));
    Eigen::SparseMatrix<double> PF_sp((int)(ncluster + nflash), (int)(nbundle + nflash));
    P_sp.setFromTriplets(P_trip.begin(), P_trip.end());
    PF_sp.setFromTriplets(PF_trip.begin(), PF_trip.end());
    Ress::vector_t initial = Ress::vector_t::Zero(nbundle + nflash);
    for (std::size_t n = 0; n < col_bundles.size(); ++n) initial(n) = 1.0;

    Ress::Params params;
    params.model = Ress::lasso;
    params.lambda = lambda;
    log->debug("solving (round 1, shared)");
    const int xdim = (int)(nbundle + nflash);
    Ress::vector_t solution;
    if (m_sparse_lasso) {
        Eigen::SparseMatrix<double> PT  = P_sp.transpose();
        Eigen::SparseMatrix<double> PFT = PF_sp.transpose();
        Ress::vector_t y = PT * M + PFT * MF;
        Eigen::SparseMatrix<double> Xs = PT * P_sp + PFT * PF_sp;
        const double t_build = ms_since(t_build0);
        const auto t_solve0 = wallclock::now();
        solution = Ress::solve(Xs, y, params, initial, weights);
        log->debug("QLtiming fit_round1_shared: nbundle {} nflash {} nopdet {} ncluster {} "
                   "matrix_build {:.1f} lasso_solve {:.1f} ms (X {}x{})",
                   nbundle, nflash, nopdet, ncluster, t_build, ms_since(t_solve0), xdim, xdim);
    }
    else {
        Ress::matrix_t P  = P_sp.toDense();
        Ress::matrix_t PF = PF_sp.toDense();
        Ress::matrix_t PT  = P.transpose();
        Ress::matrix_t PFT = PF.transpose();
        Ress::vector_t y = PT * M + PFT * MF;
        Ress::matrix_t X = PT * P + PFT * PF;
        const double t_build = ms_since(t_build0);
        const auto t_solve0 = wallclock::now();
        solution = Ress::solve(X, y, params, initial, weights);
        log->debug("QLtiming fit_round1_shared: nbundle {} nflash {} nopdet {} ncluster {} "
                   "matrix_build {:.1f} lasso_solve {:.1f} ms (X {}x{})",
                   nbundle, nflash, nopdet, ncluster, t_build, ms_since(t_solve0), xdim, xdim);
    }

    std::vector<TimingTPCBundleSelection> to_be_removed(runs.size());
    for (std::size_t n = 0; n < col_bundles.size(); ++n) {
        if (solution(n) <= m_strength_cutoff && !m_beamonly &&
            !pin_exempt(col_bundles[n].get(), solution(n)))
            to_be_removed[col_run[n]].push_back(col_bundles[n]);
    }
    for (std::size_t k = 0; k < runs.size(); ++k) {
        remove_bundle_selection(to_be_removed[k], runs[k].flash_bundles_map,
                                runs[k].cluster_bundles_map, runs[k].flash_cluster_bundles_map);
        remove_bundle_selection(to_be_removed[k], runs[k].pre_bundles);
    }
}

void QLMatching::fit_round2_shared(std::vector<ApaRun>& runs)
{
    const auto t_build0 = wallclock::now();
    const double lambda       = m_lasso_lambda;
    const double delta_charge = m_delta_charge;
    const double delta_shape  = m_delta_shape;

    const ApaRun& r0 = runs.front();
    const unsigned int nopdet = r0.nopdet;

    // Rebuild the physical-flash universe (round 1 may have emptied flashes).
    std::map<int, SharedFlash> phys;
    for (std::size_t k = 0; k < runs.size(); ++k) {
        for (auto* f : flash_iter_order(runs[k].flash_bundles_map)) {
            phys[f->get_flash_id()].insts.emplace_back(k, f);
        }
    }
    const unsigned int nflash = phys.size();
    unsigned int nbundle = 0;
    for (auto& [fid, pf] : phys) {
        (void)fid;
        for (auto& [k, f] : pf.insts) nbundle += runs[k].flash_bundles_map.at(f).size();
    }
    std::map<Cluster*, int> cluster_idx_map;
    {
        int idx = 0;
        for (auto& run : runs) {
            for (auto* c : cluster_iter_order(run.cluster_bundles_map, run.global_cluster_idx_map)) {
                cluster_idx_map[c] = idx++;
            }
        }
    }
    const unsigned int ncluster = cluster_idx_map.size();

    Ress::vector_t M = Ress::vector_t::Zero(nopdet * nflash);
    Ress::vector_t MF = Ress::vector_t::Zero(ncluster);
    Ress::vector_t weights = Ress::vector_t::Zero(nbundle);
    std::vector<Eigen::Triplet<double>> P_trip, PF_trip;
    std::vector<TimingTPCBundle::pointer> col_bundles;
    std::vector<std::size_t> col_run;
    col_bundles.reserve(nbundle);
    col_run.reserve(nbundle);

    std::size_t i = 0, ik = 0;
    for (auto& [fid, pf] : phys) {
        (void)fid;
        Opflash* f0 = pf.insts.front().second;
        for (unsigned int j = 0; j < nopdet; ++j) {
            const int opdet_idx = r0.opdet_idx_v.at(j);
            const double pe = f0->get_PE(opdet_idx);
            const double pe_err = std::sqrt(f0->get_PE(opdet_idx) + std::pow(f0->get_PE_err(opdet_idx), 2));
            M(i * nopdet + j) = pe / pe_err;
        }
        for (auto& [k, f] : pf.insts) {
            for (auto bundle : runs[k].flash_bundles_map.at(f)) {
                const auto ks_dis = bundle->get_ks_dis();
                const auto& pred_flash = bundle->get_pred_flash();
                for (unsigned int j = 0; j < nopdet; ++j) {
                    const int opdet_idx = r0.opdet_idx_v.at(j);
                    const double pred_pe = pred_flash.at(opdet_idx);
                    const double pe_err = std::sqrt(f0->get_PE(opdet_idx) + std::pow(f0->get_PE_err(opdet_idx), 2));
                    P_trip.emplace_back((int)(i * nopdet + j), (int)col_bundles.size(), pred_pe / pe_err);
                }
                const auto meas_pe_tot = f0->get_total_PE();
                const auto pred_pe_tot = bundle->get_total_pred_light();
                const double base = (std::abs(pred_pe_tot - meas_pe_tot) > m_pe_mismatch_knee * meas_pe_tot)
                                        ? std::abs(pred_pe_tot - meas_pe_tot) / meas_pe_tot
                                        : m_pe_mismatch_floor;
                weights(ik++) = (base + delta_shape * nopdet * ks_dis / lambda) * lasso_flag_factor(bundle);
                col_bundles.push_back(bundle);
                col_run.push_back(k);
            }
        }
        ++i;
    }
    for (unsigned int k = 0; k < ncluster; ++k) MF(k) = 1. / delta_charge;
    for (std::size_t n = 0; n < col_bundles.size(); ++n) {
        PF_trip.emplace_back(cluster_idx_map.at(col_bundles[n]->get_main_cluster()), (int)n, 1. / delta_charge);
    }

    Eigen::SparseMatrix<double> P_sp((int)(nopdet * nflash), (int)nbundle);
    Eigen::SparseMatrix<double> PF_sp((int)ncluster, (int)nbundle);
    P_sp.setFromTriplets(P_trip.begin(), P_trip.end());
    PF_sp.setFromTriplets(PF_trip.begin(), PF_trip.end());
    Ress::vector_t initial = Ress::vector_t::Zero(nbundle);
    for (std::size_t n = 0; n < col_bundles.size(); ++n) initial(n) = 1.0;

    Ress::Params params;
    params.model = Ress::lasso;
    params.lambda = lambda;
    log->debug("solving (round 2, shared)");
    const int xdim = (int)nbundle;
    Ress::vector_t solution;
    if (m_sparse_lasso) {
        Eigen::SparseMatrix<double> PT  = P_sp.transpose();
        Eigen::SparseMatrix<double> PFT = PF_sp.transpose();
        Ress::vector_t y = PT * M + PFT * MF;
        Eigen::SparseMatrix<double> Xs = PT * P_sp + PFT * PF_sp;
        const double t_build = ms_since(t_build0);
        const auto t_solve0 = wallclock::now();
        solution = Ress::solve(Xs, y, params, initial, weights);
        log->debug("QLtiming fit_round2_shared: nbundle {} nflash {} nopdet {} ncluster {} "
                   "matrix_build {:.1f} lasso_solve {:.1f} ms (X {}x{})",
                   nbundle, nflash, nopdet, ncluster, t_build, ms_since(t_solve0), xdim, xdim);
    }
    else {
        Ress::matrix_t P  = P_sp.toDense();
        Ress::matrix_t PF = PF_sp.toDense();
        Ress::matrix_t PT  = P.transpose();
        Ress::matrix_t PFT = PF.transpose();
        Ress::vector_t y = PT * M + PFT * MF;
        Ress::matrix_t X = PT * P + PFT * PF;
        const double t_build = ms_since(t_build0);
        const auto t_solve0 = wallclock::now();
        solution = Ress::solve(X, y, params, initial, weights);
        log->debug("QLtiming fit_round2_shared: nbundle {} nflash {} nopdet {} ncluster {} "
                   "matrix_build {:.1f} lasso_solve {:.1f} ms (X {}x{})",
                   nbundle, nflash, nopdet, ncluster, t_build, ms_since(t_solve0), xdim, xdim);
    }

    std::vector<TimingTPCBundleSelection> to_be_removed(runs.size());
    for (std::size_t n = 0; n < col_bundles.size(); ++n) {
        col_bundles[n]->set_strength(solution(n));
        if (!(solution(n) > m_strength_cutoff || m_beamonly ||
              pin_exempt(col_bundles[n].get(), solution(n))))
            to_be_removed[col_run[n]].push_back(col_bundles[n]);
    }
    for (std::size_t k = 0; k < runs.size(); ++k) {
        remove_bundle_selection(to_be_removed[k], runs[k].flash_bundles_map,
                                runs[k].cluster_bundles_map, runs[k].flash_cluster_bundles_map);
        remove_bundle_selection(to_be_removed[k], runs[k].pre_bundles);
    }

    // §I/§J rescues: per-run "empty flash" concepts, wrong when the flash is
    // shared across drift sides (a flash explained by the other side's clusters
    // is not empty). Skipped here; configure() warned if they were requested.
    // fit_round2's matched_pairs/organize_bundles tail is likewise not replicated
    // (its result is discarded there — the matched output is what remains in
    // flash_bundles_map).

    // Rescue blind-spot fix (doc 23 phase 1a; default OFF = bit-identical): see
    // the fit_round2 counterpart. Cull postcull-doomed unflagged selections
    // BEFORE the rescues so they cannot mark their cluster off-limits; the
    // legacy post-rescue call below still runs.
    if (m_postcull_before_rescue) {
        for (auto& run : runs) cull_unflagged_lowquality(run);
    }

    // Shared-flash-aware rescues (doc 19 phase 5; default OFF = bit-identical).
    // Empty-flash rescue uses JOINT emptiness (no side holds the flash); the
    // cluster-centric rescue is side-local by construction (ADD-only, a shared
    // flash legitimately holds bundles of both sides) so §J applies per run.
    if (m_empty_rescue_shared) rescue_empty_flashes_shared(runs);
    if (m_cluster_rescue_shared) {
        for (auto& run : runs) rescue_unmatched_clusters(run, run.prefit_snapshot);
    }

    // Post-fit unflagged low-quality cull (doc 19; OFF default = bit-identical).
    for (auto& run : runs) cull_unflagged_lowquality(run);

    log->debug("done with matching (shared)");
}

// Empty-flash light-quality rescue (§I; the prototype's flash-centric light pick).
// The LASSO selects by strength and can leave a flash with NO surviving bundle even
// when a cluster is an excellent LIGHT match for it (that cluster was won, on strength
// alone, by a neighbouring flash). For each emptied flash, adopt its best light-quality
// candidate from the pre-cutoff snapshot if it clears the m_rescue_metric_max bar. If
// that cluster is already matched elsewhere, reassign it ONLY when the empty flash is a
// strictly better light match (the guard makes the rescue non-regressive: it never
// removes a better match, so it cannot drop a currently-correct pair). Mutates
// run.flash_bundles_map in place.
void QLMatching::rescue_empty_flashes(ApaRun& run, const FlashBundlesMap& snapshot)
{
    // Light-quality metric: ks*(chi2/ndf)^exp, with a per-flag down-weight for
    // boundary / near-PMT bundles (prototype 0.8 then 0.64). Lower = better.
    auto metric = [this](const TimingTPCBundle::pointer& b) {
        const int ndf = std::max(b->get_ndf(), 1);
        double m = b->get_ks_dis() * std::pow(b->get_chi2() / ndf, m_rescue_exponent);
        if (b->get_flag_at_x_boundary()) m *= m_rescue_boundary_weight;
        if (b->get_flag_close_to_PMT()) m *= m_rescue_boundary_weight;
        return m;
    };

    // Where each currently-matched cluster lives: main_cluster -> (flash, metric).
    std::map<Cluster*, std::pair<Opflash*, double>> matched;
    for (auto& kv : run.flash_bundles_map) {
        for (auto& b : kv.second) matched[b->get_main_cluster()] = {kv.first, metric(b)};
    }
    // xtpc joint-pin: clusters bound to a pinned crosser flash must never be reassigned to
    // an empty flash by the light-quality guard below (the pin overrides light). Empty if off.
    std::set<Cluster*> pin_locked;
    for (auto& kv : run.flash_bundles_map)
        for (auto& b : kv.second)
            if (b->get_flag_xtpc_pin()) pin_locked.insert(b->get_main_cluster());

    // Empty flashes (in snapshot, absent/empty in the live map), in flash-id order.
    std::vector<Opflash*> empty_flashes;
    for (auto& kv : snapshot) {
        auto* flash = kv.first;
        auto it = run.flash_bundles_map.find(flash);
        if (it == run.flash_bundles_map.end() || it->second.empty()) {
            empty_flashes.push_back(flash);
        }
    }
    std::sort(empty_flashes.begin(), empty_flashes.end(),
              [](Opflash* a, Opflash* b) { return a->get_flash_id() < b->get_flash_id(); });

    int n_rescued = 0;
    for (auto* flash : empty_flashes) {
        // Best light-quality candidate on this flash (tie-break by cluster_index_id).
        TimingTPCBundle::pointer best;
        double best_m = 0;
        for (const auto& b : snapshot.at(flash)) {
            const double m = metric(b);
            if (!best || m < best_m ||
                (m == best_m && b->get_cluster_index_id() < best->get_cluster_index_id())) {
                best = b;
                best_m = m;
            }
        }
        if (!best || best_m > m_rescue_metric_max) continue;  // quality bar

        auto* C = best->get_main_cluster();
        if (pin_locked.count(C)) continue;  // xtpc joint-pin owns C: do not reassign
        // One flash per cluster: never leave C double-matched.
        auto mit = matched.find(C);
        if (mit == matched.end()) {
            // C is unmatched: pure addition.
            run.flash_bundles_map[flash].push_back(best);
        }
        else {
            // C is matched to X: REASSIGN (remove X,C; add F,C) only when this flash is
            // a strictly better light match (mF < mX). The guard never removes a better
            // match, so it cannot drop a correct pair on light-confident cases; the
            // tight m_rescue_metric_max bar restricts it to high-confidence rescues.
            if (!(best_m < mit->second.second)) continue;
            auto* X = mit->second.first;
            auto& xv = run.flash_bundles_map[X];
            xv.erase(std::remove_if(xv.begin(), xv.end(),
                                    [C](const TimingTPCBundle::pointer& b) {
                                        return b->get_main_cluster() == C;
                                    }),
                     xv.end());
            if (xv.empty()) run.flash_bundles_map.erase(X);
            run.flash_bundles_map[flash].push_back(best);
        }
        matched[C] = {flash, best_m};
        ++n_rescued;
        log->debug("QLrescue: flash id {} adopted cluster ident {} (metric {})",
                   flash->get_flash_id(), C->ident(), best_m);
    }

    // Re-sort inner vectors by cluster_index_id for build-stable output.
    for (auto& kv : run.flash_bundles_map) {
        std::sort(kv.second.begin(), kv.second.end(),
                  [](const TimingTPCBundle::pointer& a, const TimingTPCBundle::pointer& b) {
                      return a->get_cluster_index_id() < b->get_cluster_index_id();
                  });
    }
    log->debug("QLrescue: rescued {} empty flash(es)", n_rescued);
}

// Cluster-centric rescue (§J; the complement of the flash-centric empty-flash pick).
// The empty-flash rescue fills only flashes left WHOLLY empty; a big cluster with an
// excellent light candidate but driven to strength 0 by the LASSO L1 sparsity (a rival
// already explains its flash) stays unmatched even when that flash is non-empty. For
// each cluster still unmatched after the empty-flash pass, adopt its best ACCEPTED
// candidate from the pre-cutoff snapshot, attaching it even onto an already-non-empty
// flash (many-clusters-per-flash is physical and GT-endorsed). We only ADD (no reassign),
// for clusters with zero live bundles, one bundle each => one-flash-per-cluster holds and
// the addition set is independent of cluster iteration order. Mutates run.flash_bundles_map.
void QLMatching::rescue_unmatched_clusters(ApaRun& run, const FlashBundlesMap& snapshot)
{
    // Acceptance bar: ks < ks_max AND chi2/ndf < chi2ndf_max AND ratio_lo < pred/meas <
    // ratio_hi, with ratio = total predicted light / flash total measured PE (the same
    // pairing the LASSO PE-mismatch weight uses).
    // Saturation-aware ratio-high extension (doc 23 phase 1b; default OFF =>
    // the plain ratio < ratio_hi test). When the railed fraction of the
    // flash's measured PE exceeds sat_frac_min the measurement is
    // right-censored (railed channels report a lower bound) and pred/meas is
    // an overestimate: extend the high gate to ratio_hi * sat_ratio_mult.
    auto sat_ratio_hi_ok = [this](const TimingTPCBundle::pointer& b,
                                  double ratio, double ratio_hi) {
        if (ratio < ratio_hi) return true;
        if (!m_cluster_rescue_sat_ratio_relax) return false;
        if (!(ratio < ratio_hi * m_cluster_rescue_sat_ratio_mult)) return false;
        auto* f = b->get_flash();
        const double tot = f->get_total_PE();
        if (tot <= 0.0) return false;
        double satpe = 0;
        const int nch = f->get_num_channels();
        for (int j = 0; j < nch; ++j)
            if (f->get_sat(j)) satpe += f->get_PE(j);
        return satpe / tot > m_cluster_rescue_sat_frac_min;
    };

    auto accept = [this, &sat_ratio_hi_ok](const TimingTPCBundle::pointer& b) {
        const int ndf = std::max(b->get_ndf(), 1);
        const double c2ndf = b->get_chi2() / ndf;
        const double meas = b->get_flash()->get_total_PE();
        const double ratio = (meas > 0.0) ? b->get_total_pred_light() / meas : 1e9;
        return b->get_ks_dis() < m_cluster_rescue_ks_max
            && c2ndf < m_cluster_rescue_chi2ndf_max
            && ratio > m_cluster_rescue_ratio_lo
            && sat_ratio_hi_ok(b, ratio, m_cluster_rescue_ratio_hi);
    };
    // Light-quality score (lower = better): KS shape, chi2/ndf, and PE-scale agreement.
    auto score = [](const TimingTPCBundle::pointer& b) {
        const int ndf = std::max(b->get_ndf(), 1);
        const double meas = b->get_flash()->get_total_PE();
        const double ratio = (meas > 0.0) ? b->get_total_pred_light() / meas : 1e9;
        return b->get_ks_dis() * std::sqrt(std::max(b->get_chi2() / ndf, 1.0))
             + std::abs(std::log(ratio));
    };

    // Clusters already matched (any live bundle) are off-limits.
    std::set<Cluster*> matched;
    for (auto& kv : run.flash_bundles_map)
        for (auto& b : kv.second) matched.insert(b->get_main_cluster());

    // Group candidates by main cluster. Default pool = the post-cull snapshot (the
    // shipped behaviour). With m_cluster_rescue_precull, use the PRE-cull universe
    // instead: run.all_bundles minus the build-prefilter rejects (potential_bad_match).
    // {all_bundles : !bad} == pre_bundles BEFORE cull_inconsistent removed the
    // non-consistent rivals, so it restores cull-victim candidates (valid ks/chi2/pred,
    // since the non-bad bundles were fully examine_bundle'd) the snapshot can't reach.
    // Both pools are vectors in deterministic build order.
    //
    // m_cluster_rescue_precull_additive: keep the snapshot pool as the PRIMARY and use
    // the pre-cull universe only as a per-cluster FALLBACK. The primary decides every
    // cluster the shipped snapshot rescue already handles (so precull can only add,
    // never re-decide/mis-switch); the fallback supplies the cull-victim clusters the
    // snapshot cannot reach.
    std::map<Cluster*, std::vector<TimingTPCBundle::pointer>> by_cluster, by_cluster_fallback;
    const bool additive = m_cluster_rescue_precull && m_cluster_rescue_precull_additive;
    if (m_cluster_rescue_precull && !additive) {
        for (auto& b : run.all_bundles)
            if (!b->get_potential_bad_match_flag())
                by_cluster[b->get_main_cluster()].push_back(b);
    }
    else {
        for (auto& kv : snapshot)
            for (auto& b : kv.second) by_cluster[b->get_main_cluster()].push_back(b);
    }
    if (additive) {
        for (auto& b : run.all_bundles)
            if (!b->get_potential_bad_match_flag())
                by_cluster_fallback[b->get_main_cluster()].push_back(b);
    }

    // Iterate unmatched clusters in the canonical global_cluster_idx order (not pointer
    // order) for determinism. In additive mode the fallback pool is the superset, so
    // draw the cluster list from it; otherwise the single primary pool.
    const auto& cluster_source = additive ? by_cluster_fallback : by_cluster;
    std::vector<Cluster*> clusters;
    for (auto& kv : cluster_source)
        if (!matched.count(kv.first)) clusters.push_back(kv.first);
    std::sort(clusters.begin(), clusters.end(), [&run](Cluster* a, Cluster* b) {
        return run.global_cluster_idx_map.at(a) < run.global_cluster_idx_map.at(b);
    });

    // Best accepting bundle for one cluster (lower score wins; deterministic tie-break
    // by flash id then cluster_index_id). Returns null when no candidate accepts.
    auto pick_best = [&accept, &score](const std::vector<TimingTPCBundle::pointer>& pool) {
        TimingTPCBundle::pointer best;
        double best_s = 0;
        for (const auto& b : pool) {
            if (!accept(b)) continue;
            const double s = score(b);
            if (!best || s < best_s ||
                (s == best_s && b->get_flash()->get_flash_id() < best->get_flash()->get_flash_id()) ||
                (s == best_s && b->get_flash()->get_flash_id() == best->get_flash()->get_flash_id()
                             && b->get_cluster_index_id() < best->get_cluster_index_id())) {
                best = b;
                best_s = s;
            }
        }
        return best;
    };

    int n_rescued = 0;
    for (auto* C : clusters) {
        // Primary pool (snapshot when additive, else the selected pool); fall back to
        // the pre-cull pool only when additive and the primary rescued nothing.
        TimingTPCBundle::pointer best = pick_best(by_cluster[C]);
        if (!best && additive) best = pick_best(by_cluster_fallback[C]);
        if (!best) continue;
        run.flash_bundles_map[best->get_flash()].push_back(best);  // attach (flash may be non-empty)
        ++n_rescued;
        log->debug("QLclusrescue: flash id {} adopted cluster ident {} (score {})",
                   best->get_flash()->get_flash_id(), C->ident(), score(best));
    }

    // Relaxed SECOND-CHANCE tier (§J2, doc 21; default OFF => byte-identical).
    // Clusters the tight pass above left unmatched get one more accept() with the
    // relaxed gates, restricted to LONG clusters. Same pools, same score() and
    // deterministic tie-break, same additive-only guarantee (only ADD onto
    // flash_bundles_map, never reassign) — an existing match can never change.
    // Adoptions are stamped flag_cluster_rescue_relaxed = the low-confidence mark.
    if (m_cluster_rescue_relaxed) {
        log->debug("QLclusrescue-relaxed on: ks<{} c2ndf<{} ratio {}-{} min_len {} cm",
                   m_cluster_rescue_relaxed_ks_max, m_cluster_rescue_relaxed_chi2ndf_max,
                   m_cluster_rescue_relaxed_ratio_lo, m_cluster_rescue_relaxed_ratio_hi,
                   m_cluster_rescue_relaxed_min_length / units::cm);
        auto accept_relaxed = [this, &sat_ratio_hi_ok](const TimingTPCBundle::pointer& b) {
            const int ndf = std::max(b->get_ndf(), 1);
            const double c2ndf = b->get_chi2() / ndf;
            const double meas = b->get_flash()->get_total_PE();
            const double ratio = (meas > 0.0) ? b->get_total_pred_light() / meas : 1e9;
            return b->get_ks_dis() < m_cluster_rescue_relaxed_ks_max
                && c2ndf < m_cluster_rescue_relaxed_chi2ndf_max
                && ratio > m_cluster_rescue_relaxed_ratio_lo
                && sat_ratio_hi_ok(b, ratio, m_cluster_rescue_relaxed_ratio_hi);
        };
        auto pick_best_relaxed = [&accept_relaxed, &score](const std::vector<TimingTPCBundle::pointer>& pool) {
            TimingTPCBundle::pointer best;
            double best_s = 0;
            for (const auto& b : pool) {
                if (!accept_relaxed(b)) continue;
                const double s = score(b);
                if (!best || s < best_s ||
                    (s == best_s && b->get_flash()->get_flash_id() < best->get_flash()->get_flash_id()) ||
                    (s == best_s && b->get_flash()->get_flash_id() == best->get_flash()->get_flash_id()
                                 && b->get_cluster_index_id() < best->get_cluster_index_id())) {
                    best = b;
                    best_s = s;
                }
            }
            return best;
        };
        // Re-derive the matched set AFTER the tight pass (its adoptions are already
        // inside flash_bundles_map), so the tight loop above stays textually legacy.
        std::set<Cluster*> matched_now;   // membership only, never iterated
        for (auto& kv : run.flash_bundles_map)
            for (auto& b : kv.second) matched_now.insert(b->get_main_cluster());
        int n_relaxed = 0;
        for (auto* C : clusters) {
            if (matched_now.count(C)) continue;
            if (C->get_length() < m_cluster_rescue_relaxed_min_length) continue;
            TimingTPCBundle::pointer best = pick_best_relaxed(by_cluster[C]);
            if (!best && additive) best = pick_best_relaxed(by_cluster_fallback[C]);
            if (!best) continue;
            best->set_flag_cluster_rescue_relaxed(true);
            run.flash_bundles_map[best->get_flash()].push_back(best);
            ++n_relaxed;
            log->debug("QLclusrescue-relaxed: flash id {} adopted cluster ident {} (len {} cm, score {})",
                       best->get_flash()->get_flash_id(), C->ident(),
                       C->get_length() / units::cm, score(best));
        }
        log->debug("QLclusrescue-relaxed: rescued {} long unmatched cluster(s)", n_relaxed);
    }

    // Re-sort inner vectors by cluster_index_id for build-stable output (as §I does).
    for (auto& kv : run.flash_bundles_map) {
        std::sort(kv.second.begin(), kv.second.end(),
                  [](const TimingTPCBundle::pointer& a, const TimingTPCBundle::pointer& b) {
                      return a->get_cluster_index_id() < b->get_cluster_index_id();
                  });
    }
    log->debug("QLclusrescue: rescued {} unmatched cluster(s)", n_rescued);
}

// Shared-flash empty-flash rescue (m_empty_rescue_shared; doc 19 phase 5).
// The per-run §I notion of "empty flash" is wrong under shared_flash: a flash
// explained by the OTHER drift side's clusters is not empty.  Here a physical
// flash (keyed by flash row id, identical across ports by construction — see
// fit_round1_shared's guard) is empty only when NO run has a surviving bundle
// for its instance.  The best snapshot candidate ACROSS sides is adopted with
// the same metric / metric_max bar / reassign-only-if-strictly-better /
// pin-locked rules as §I, applied within the candidate's own run.
void QLMatching::rescue_empty_flashes_shared(std::vector<ApaRun>& runs)
{
    auto metric = [this](const TimingTPCBundle::pointer& b) {
        const int ndf = std::max(b->get_ndf(), 1);
        double m = b->get_ks_dis() * std::pow(b->get_chi2() / ndf, m_rescue_exponent);
        if (b->get_flag_at_x_boundary()) m *= m_rescue_boundary_weight;
        if (b->get_flag_close_to_PMT()) m *= m_rescue_boundary_weight;
        return m;
    };

    // Where each currently-matched cluster lives (clusters are per-run objects,
    // so one global map is collision-free): cluster -> (run idx, flash, metric).
    struct Cur { std::size_t k; Opflash* flash; double m; };
    std::map<Cluster*, Cur> matched;
    std::set<Cluster*> pin_locked;
    for (std::size_t k = 0; k < runs.size(); ++k) {
        for (auto& kv : runs[k].flash_bundles_map) {
            for (auto& b : kv.second) {
                matched[b->get_main_cluster()] = {k, kv.first, metric(b)};
                if (b->get_flag_xtpc_pin()) pin_locked.insert(b->get_main_cluster());
            }
        }
    }

    // Physical flashes present in any snapshot, by ascending flash row id;
    // live occupancy checked across ALL runs.
    struct Phys { std::vector<std::pair<std::size_t, Opflash*>> insts; };
    std::map<int, Phys> phys;
    for (std::size_t k = 0; k < runs.size(); ++k) {
        for (auto& kv : runs[k].prefit_snapshot) {
            phys[kv.first->get_flash_id()].insts.emplace_back(k, kv.first);
        }
    }

    int n_rescued = 0;
    for (auto& [fid, pf] : phys) {
        bool live = false;
        for (auto& [k, f] : pf.insts) {
            auto it = runs[k].flash_bundles_map.find(f);
            if (it != runs[k].flash_bundles_map.end() && !it->second.empty()) {
                live = true;
                break;
            }
        }
        if (live) continue;

        // Best candidate across sides; deterministic tie-break (run idx, then
        // cluster_index_id) — insts are in ascending run order by construction.
        TimingTPCBundle::pointer best;
        std::size_t best_k = 0;
        Opflash* best_f = nullptr;
        double best_m = 0;
        for (auto& [k, f] : pf.insts) {
            auto sit = runs[k].prefit_snapshot.find(f);
            if (sit == runs[k].prefit_snapshot.end()) continue;
            for (const auto& b : sit->second) {
                const double m = metric(b);
                if (!best || m < best_m ||
                    (m == best_m && b->get_cluster_index_id() < best->get_cluster_index_id())) {
                    best = b;
                    best_k = k;
                    best_f = f;
                    best_m = m;
                }
            }
        }
        if (!best || best_m > m_rescue_metric_max) continue;

        auto* C = best->get_main_cluster();
        if (pin_locked.count(C)) continue;
        auto mit = matched.find(C);
        if (mit == matched.end()) {
            runs[best_k].flash_bundles_map[best_f].push_back(best);
        }
        else {
            if (!(best_m < mit->second.m)) continue;
            auto& run_x = runs[mit->second.k];
            auto* X = mit->second.flash;
            auto& xv = run_x.flash_bundles_map[X];
            xv.erase(std::remove_if(xv.begin(), xv.end(),
                                    [C](const TimingTPCBundle::pointer& b) {
                                        return b->get_main_cluster() == C;
                                    }),
                     xv.end());
            if (xv.empty()) run_x.flash_bundles_map.erase(X);
            runs[best_k].flash_bundles_map[best_f].push_back(best);
        }
        matched[C] = {best_k, best_f, best_m};
        ++n_rescued;
        log->debug("QLrescue-shared: flash id {} adopted cluster ident {} "
                   "(run {} metric {})",
                   fid, C->ident(), best_k, best_m);
    }

    for (auto& run : runs) {
        for (auto& kv : run.flash_bundles_map) {
            std::sort(kv.second.begin(), kv.second.end(),
                      [](const TimingTPCBundle::pointer& a, const TimingTPCBundle::pointer& b) {
                          return a->get_cluster_index_id() < b->get_cluster_index_id();
                      });
        }
    }
    log->debug("QLrescue-shared: rescued {} empty flash(es)", n_rescued);
}

// Stamp the matched flash/t0 onto each cluster (and its associated sub-clusters):
// cluster_t0, the "flash" row index (Cluster::get_flash()), a globally-unique
// matched_flash_gid (survives the all-APA merge), and the per-channel flashpred.
void QLMatching::apply_matched_t0s(ApaRun& run)
{
    for (auto* flash : flash_iter_order(run.flash_bundles_map)) {
        for (auto bundle : run.flash_bundles_map[flash]) {
            auto* cluster = bundle->get_main_cluster();
            const double t0 = flash->get_time() * units::ns;
            // gid side: legacy = this node's anode ident; m_opflash_phys_gid = the
            // flash's PHYSICAL side, so a cross-side (xTPC) match on the opposite
            // node still keys the same gid the owning side emits in write_opflash_pc.
            // m_shared_flash = side 0 always: one global flash list shared by every
            // port, emitted once (input 0) in write_opflash_pc, so every run's
            // matched_flash_gid resolves to that single emission.
            const int gid_side = m_shared_flash ? 0
                                 : m_opflash_phys_gid ? flash_phys_side(flash)
                                                      : (int) run.anode->ident();
            const int flash_gid = gid_side * kFlashGidStride + run.global_flash_idx_map.at(flash);
            cluster->set_cluster_t0(t0);
            cluster->set_scalar<int>("flash", flash->get_flash_id());
            cluster->set_scalar<int>("matched_flash_gid", flash_gid);
            cluster->put_pcarray<double>(bundle->get_pred_flash(), "pe", "flashpred");
            // Propagate the group's matched flash/t0 to its associated sub-clusters.
            for (auto* oc : bundle->get_other_clusters()) {
                oc->set_cluster_t0(t0);
                oc->set_scalar<int>("flash", flash->get_flash_id());
                oc->set_scalar<int>("matched_flash_gid", flash_gid);
            }
            log->debug("flash_bundles_map: flash id {} time {} ns, cluster gidx {} "
                       "total_pred_light {} t0 {}",
                       flash->get_flash_id(), flash->get_time(),
                       run.global_cluster_idx_map[bundle->get_main_cluster()],
                       bundle->get_total_pred_light(),
                       bundle->get_main_cluster()->get_cluster_t0());
            log->debug("bundle_flags: flash id {} cluster ident {} gidx {} | "
                       "at_x_boundary {} close_to_PMT {} spec_end {} "
                       "window_truncated {} high_consistent {}",
                       flash->get_flash_id(),
                       bundle->get_main_cluster()->ident(),
                       run.global_cluster_idx_map[bundle->get_main_cluster()],
                       bundle->get_flag_at_x_boundary(),
                       bundle->get_flag_close_to_PMT(),
                       bundle->get_spec_end_flag(),
                       bundle->get_flag_window_truncated(),
                       bundle->get_consistent_flag());
        }
    }
}

// Persist the per-flash measured light into a merge-safe, self-contained per-root
// "opflash" PC (one row per (flash, channel): gid, time, ch, pe), so the all-APA
// MABC can dump the op/flash Bee display after the per-APA trees are merged.
// Physical drift side (0 = low-x / TPC0, 1 = high-x / TPC1) of a flash, from where
// its measured light sits relative to the cathode.  m_opdets[ch].center.x() is the
// OpDet position; ties (e.g. an all-zero flash) resolve to side 0.
int QLMatching::flash_phys_side(const Opflash* flash) const
{
    double pe_lo = 0.0, pe_hi = 0.0;
    for (int ch = 0; ch < m_nchan; ++ch) {
        const double p = flash->get_PE(ch);
        if (ch < (int) m_opdets.size() && m_opdets[ch].center.x() >= m_cathode_x) pe_hi += p;
        else pe_lo += p;
    }
    return (pe_hi > pe_lo) ? 1 : 0;
}

void QLMatching::write_opflash_pc(ApaRun& run)
{
    // Shared-flash mode: every port carries the SAME global flash list, so emit it
    // once (from input 0, gid side 0) — otherwise the merged root would hold every
    // physical flash N times. merge_pct simply skips roots without an opflash PC.
    if (m_shared_flash && run.input_idx != 0) return;
    // Per-side clocks: when the charge crates frame independently (per-input
    // trigger_offsets configured, e.g. PDVD BDE/TDE), also emit the flash time
    // on input-1's (top) clock as a "time1" column so display consumers need
    // not re-derive the inter-crate skew.  Key absent when trigger_offsets is
    // empty (PDHD/SBND scalar path) => those outputs stay byte-identical.
    const bool per_side_clocks = m_shared_flash && !m_trigger_offsets.empty();
    std::vector<int> op_gid, op_ch, op_apa;
    std::vector<double> op_time, op_pe, op_time1;
    for (std::size_t fi = 0; fi < run.flashes.size(); ++fi) {
        const auto& flash = run.flashes[fi];
        // Physical drift side of the flash (TPC 0 = low-x, TPC 1 = high-x of the
        // cathode), taken from where its measured light actually is rather than
        // from the node's anode ident (the per-side matcher runs against one global
        // flash list, so the node ident is NOT the lit volume).  Used for the Bee
        // "apa" side tag always, and for the gid when m_opflash_phys_gid so the two
        // per-side nodes emit ONE gid per physical flash and the duplicate collapses
        // in fill_bee_flashes.  Legacy gid = node anode ident (byte-identical when
        // each node's flash list is already one-per-side, e.g. SBND per-TPC flashes).
        const int apa = flash_phys_side(flash.get());
        const int gid_side = m_shared_flash ? 0
                             : m_opflash_phys_gid ? apa : (int) run.anode->ident();
        const int gid = gid_side * kFlashGidStride + static_cast<int>(fi);
        for (int ch = 0; ch < m_nchan; ++ch) {
            op_gid.push_back(gid);
            // Fold the per-event readout-vs-trigger offset into the displayed flash
            // time so the Bee red box lands on the (raw-x) charge it matches; the
            // matching geometry already carries the offset (see flash_x_offset).
            // 0 => bit-identical. Under shared_flash this emits for input 0 only,
            // so the per-input value is this run's own.
            op_time.push_back(flash->get_time() + trigger_offset_for(run.input_idx));
            op_ch.push_back(ch);
            op_pe.push_back(flash->get_PE(ch));
            op_apa.push_back(apa);
            if (per_side_clocks) op_time1.push_back(flash->get_time() + trigger_offset_for(1));
        }
    }
    run.grouping->put_pcarray<int>(op_gid, "gid", "opflash");
    run.grouping->put_pcarray<double>(op_time, "time", "opflash");
    run.grouping->put_pcarray<int>(op_ch, "ch", "opflash");
    run.grouping->put_pcarray<double>(op_pe, "pe", "opflash");
    run.grouping->put_pcarray<int>(op_apa, "apa", "opflash");
    if (per_side_clocks) run.grouping->put_pcarray<double>(op_time1, "time1", "opflash");
}

// ---------------------------------------------------------------------------
// Hand-scan calibration dump. Serialize the FULL candidate-bundle universe
// (run.all_bundles across every APA) plus the per-flash measured light, the
// cluster geometry and the detector box, to one per-event JSON consumed by the
// off-line Q/L hand-scan event display (sbnd_xin/ql_scan). Observation-only:
// nothing here feeds back into the matching, and it runs only when calib_dump is
// set, so production output is unchanged. See match/docs/joint-qlmatching-design.md.
//
// Conventions: all spatial quantities in cm, flash time in us, drift speed in
// cm/us, so the viewer recovers a bundle's T0 x-shift as
//   dx_cm = sign_offset * flash_time_us * drift_speed_cm_per_us
// and adds it to the (fixed-box) cluster points. The per-event readout-vs-trigger
// offset is ALREADY folded into f["time"] below (so the viewer needs no separate
// trigger term); top["trigger_offset"] is therefore 0 here.
// Cluster idents are per-APA, so a
// globally-unique cluster uid = apa*kFlashGidStride + ident disambiguates the two
// TPCs (same stride convention as the flash gid).
void QLMatching::dump_calib(const std::vector<ApaRun>& runs)
{
    const double cm   = units::cm;
    const double us   = units::us;
    const double v_cm_us = m_drift_speed / (cm / us);
    auto cluster_uid = [](int apa, int ident) { return apa * kFlashGidStride + ident; };

    Json::Value top;
    top["nchan"]        = m_nchan;
    top["drift_speed"]  = v_cm_us;             // cm/us
    // The per-event readout-vs-trigger offset is folded into f["time"] below, so the
    // viewer must NOT re-add it; keep this 0 so a folding viewer stays a no-op.
    top["trigger_offset"] = 0.0;   // us (already in f["time"])
    // Per-input offsets, for reference only (already folded into f["time"] /
    // the per-bundle geometry). Emitted only when configured => byte-identical
    // dumps for scalar-offset detectors.
    if (!m_trigger_offsets.empty()) {
        Json::Value tov(Json::arrayValue);
        for (double t : m_trigger_offsets) tov.append(t / us);
        top["trigger_offsets_us"] = tov;
    }
    // Per-input drift speeds, cm/us (already applied to the per-bundle x
    // offsets). Emitted only when configured => byte-identical dumps for
    // scalar-speed detectors; top["drift_speed"] stays the scalar.
    if (!m_drift_speeds.empty()) {
        Json::Value dsv(Json::arrayValue);
        for (double s : m_drift_speeds) dsv.append(s / (cm / us));
        top["drift_speeds"] = dsv;
    }
    top["count"]        = (Json::UInt64)m_count;
    top["charge_ident"] = runs.empty() ? 0 : runs.front().charge_ident;
    // Readout window used by the window-truncation flag (raw post-resample ticks).
    // For PDHD this is the value the run script read from the SP frame and passed
    // via readout_window_ticks; lets the viewer/validation see the window in use.
    top["readout_nticks"] = m_readout_window_ticks;

    // Quality params / uncertainty model (so the viewer shows the chi2 /
    // flag_high_consistent ingredients and the PE_err rule next to the metrics).
    Json::Value qp;
    qp["highconsist_ks_max"]  = m_highconsist_ks_max;
    qp["highconsist_min_ndf"] = m_highconsist_min_ndf;
    qp["pe_ndf_knee"]         = m_bundle_pe_ndf_knee;
    qp["mask_ks"]             = m_bundle_mask_ks;
    qp["ks_merge_max"]        = m_bundle_ks_merge_max;
    qp["chi2ndf_merge_max"]   = m_bundle_chi2ndf_merge_max;
    qp["strength_cutoff"]     = m_strength_cutoff;
    qp["pe_err_floor"]        = m_pe_err_floor;
    qp["pe_err_frac"]         = m_pe_err_frac;
    qp["pe_err_knee"]         = m_pe_err_knee;
    qp["pe_err_on_pred"]      = m_pe_err_on_pred;
    qp["pe_err_lowpe_frac"]   = m_pe_err_lowpe_frac;
    qp["pe_err_lowpe_knee"]   = m_pe_err_lowpe_knee;
    // Per-family override, resolved per channel (empty arrays => off). Dumped
    // so ql_pull_diag.py can rebuild sigma without guessing the active model.
    if (!m_pe_err_family_channels.empty()) {
        auto arr = [](const std::vector<double>& v) {
            Json::Value a(Json::arrayValue);
            for (double x : v) a.append(x);
            return a;
        };
        qp["pe_err_ch_floor"]      = arr(m_pe_err_ch_floor);
        qp["pe_err_ch_frac"]       = arr(m_pe_err_ch_frac);
        qp["pe_err_ch_lowpe_frac"] = arr(m_pe_err_ch_lowpe_frac);
        qp["pe_err_ch_lowpe_knee"] = arr(m_pe_err_ch_lowpe_knee);
    }
    {
        Json::Value mps(Json::arrayValue);
        for (double s : m_measured_pe_scale) mps.append(s);
        qp["measured_pe_scale"] = mps;  // empty => no correction applied
    }
    qp["QtoL"]                = m_QtoL;
    qp["highconsist_ladder"]  = m_highconsist_ladder;
    qp["hc_clean_ks"] = m_hc_clean_ks;  qp["hc_clean_c2"] = m_hc_clean_c2;
    qp["hc_good_ks"]  = m_hc_good_ks;   qp["hc_good_c2"]  = m_hc_good_c2;
    qp["hc_tb_ks"]    = m_hc_tb_ks;     qp["hc_tb_c2"]    = m_hc_tb_c2;
    qp["hc_miss_ks"]  = m_hc_miss_ks;   qp["hc_miss_c2"]  = m_hc_miss_c2;
    qp["hc_miss_min_ndf"] = m_hc_miss_min_ndf;
    qp["chi2_relax"]       = m_chi2_relax;
    qp["chi2_pmt_excess"]  = m_chi2_pmt_excess;
    qp["chi2_pmt_ratio"]   = m_chi2_pmt_ratio;
    qp["chi2_pmt_inflate"] = m_chi2_pmt_inflate;
    qp["chi2_sat_inflate"] = m_chi2_sat_inflate;
    qp["xtpc_pin_min_strength"] = m_xtpc_pin_min_strength;
    qp["xtpc_sc1_light_gate"]   = m_xtpc_sc1_light_gate;
    qp["xtpc_sc1_ks_max"]       = m_xtpc_sc1_ks_max;
    qp["xtpc_sc1_c2n_max"]      = m_xtpc_sc1_c2n_max;
    qp["xtpc_cathode_ks_max"]   = m_xtpc_cathode_ks_max;
    qp["postcull_unflagged"]    = m_postcull_unflagged;
    qp["postcull_ks_max"]       = m_postcull_ks_max;
    qp["postcull_c2n_max"]      = m_postcull_c2n_max;
    qp["empty_rescue_shared"]   = m_empty_rescue_shared;
    qp["cluster_rescue_shared"] = m_cluster_rescue_shared;
    qp["saturation_mask_fit"] = m_saturation_mask_fit;
    qp["auto_mask"]              = m_auto_mask;
    qp["auto_mask_pe_low"]       = m_auto_mask_pe_low;
    qp["auto_mask_neighbors"]    = m_auto_mask_neighbors;
    qp["auto_mask_pe_bright"]    = m_auto_mask_pe_bright;
    qp["auto_mask_min_contrast"] = m_auto_mask_min_contrast;
    qp["auto_mask_min_flash"]    = m_auto_mask_min_flash;
    top["quality_params"] = qp;

    // OpDet table (all channels). apa side and active flag mirror compute_geometry:
    // side by cathode-plane x; active iff the channel survives its own TPC's mask.
    std::vector<bool> active(m_nchan, false);
    std::vector<bool> automasked(m_nchan, false);
    for (const auto& run : runs) {
        for (int ch = 0; ch < m_nchan && ch < (int)run.opdet_mask.size(); ++ch)
            if (run.opdet_mask[ch] == 1) active[ch] = true;
        for (int ch : run.auto_masked)
            if (ch >= 0 && ch < m_nchan) automasked[ch] = true;
    }
    Json::Value opdets(Json::arrayValue);
    for (int ch = 0; ch < m_nchan && ch < (int)m_opdets.size(); ++ch) {
        const auto& od = m_opdets[ch];
        Json::Value o;
        o["ch"]     = ch;
        o["x"]      = od.center.x();           // OpDet table already in cm
        o["y"]      = od.center.y();
        o["z"]      = od.center.z();
        o["type"]   = od.type;
        o["apa"]    = (od.center.x() < m_cathode_x) ? 0 : 1;
        o["active"] = (bool)active[ch];
        o["auto_masked"] = (bool)automasked[ch];   // dynamic (event) mask vs static ch_mask
        opdets.append(o);
    }
    top["opdets"] = opdets;

    Json::Value geometry(Json::objectValue);
    Json::Value flashes(Json::arrayValue);
    Json::Value clusters(Json::arrayValue);
    Json::Value bundles(Json::arrayValue);

    for (const auto& run : runs) {
        const int apa = (int)run.anode->ident();

        // Detector box for this TPC. cathode_x = anode_x + s*u_cathode (compute_geometry).
        Json::Value g;
        const double anode_x_cm   = run.anode_x / cm;
        const double cathode_x_cm = (run.anode_x + run.s * run.u_cathode) / cm;
        g["anode_x"]     = anode_x_cm;
        g["cathode_x"]   = cathode_x_cm;
        g["u_cathode"]   = run.u_cathode / cm;
        g["s"]           = run.s;
        g["sign_offset"] = run.sign_offset;
        // Per-side speed the matching actually used for this volume; key
        // absent when the drift_speeds knob is off => byte-identical dumps.
        if (!m_drift_speeds.empty())
            g["drift_speed"] = drift_speed_for(run.input_idx) / (cm / us);
        g["y_lo"]        = run.y_lo / cm;
        g["y_hi"]        = run.y_hi / cm;
        g["z_lo"]        = run.z_lo / cm;
        g["z_hi"]        = run.z_hi / cm;
        geometry[std::to_string(apa)] = g;

        // Per-flash measured light + coincidence group id (computed below).
        // Shared-flash mode: the ports carry the SAME global flash list — emit it
        // once (input 0) with side-0 gids, matching write_opflash_pc/apply_matched_t0s.
        for (std::size_t fi = 0; fi < run.flashes.size(); ++fi) {
            if (m_shared_flash && run.input_idx != 0) break;
            const auto& fl = run.flashes[fi];
            Json::Value f;
            f["gid"]      = (m_shared_flash ? 0 : apa) * kFlashGidStride + (int)fi;
            f["id"]       = fl->get_flash_id();
            f["apa"]      = apa;
            f["time"]     = (fl->get_time() + trigger_offset_for(run.input_idx)) / us;   // us (trigger-folded)
            // Per-side clocks (per-input trigger_offsets, e.g. PDVD BDE/TDE):
            // "time" above is on input-0's (bottom) charge clock; also emit the
            // input-1 (top) clock so top-volume drift math needs no re-basing.
            // Key absent otherwise => PDHD/SBND dumps byte-identical.
            if (m_shared_flash && !m_trigger_offsets.empty())
                f["time1"] = (fl->get_time() + trigger_offset_for(1)) / us;
            f["total_PE"] = fl->get_total_PE();
            Json::Value pe(Json::arrayValue), pe_err(Json::arrayValue);
            for (int ch = 0; ch < m_nchan; ++ch) { pe.append(fl->get_PE(ch)); pe_err.append(fl->get_PE_err(ch)); }
            f["pe"]     = pe;
            f["pe_err"] = pe_err;
            // Per-channel DAPHNE-rail flags (only when the saturation-flag
            // masking is on; key absent otherwise => dump byte-identical).
            if (m_use_saturation_flag) {
                Json::Value sat(Json::arrayValue);
                for (int ch = 0; ch < m_nchan; ++ch) sat.append(fl->get_sat(ch) ? 1 : 0);
                f["sat"] = sat;
            }
            // Per-channel readout-coverage fractions (only when coverage
            // tracking is on; key absent otherwise => dump byte-identical).
            // The viewer labels cov < coverage_min channels "nodata".  This
            // is independent of coverage_mask_fit: the label is worth having
            // on a hand scan whether or not the channel left the fit.
            if (m_use_coverage_flag) {
                Json::Value cov(Json::arrayValue);
                for (int ch = 0; ch < m_nchan; ++ch) cov.append(fl->get_cov(ch));
                f["cov"] = cov;
            }
            f["group"]  = -1;                        // filled in the coincidence pass
            flashes.append(f);
        }

        // Cluster geometry (raw, un-shifted; the viewer applies the per-bundle dx).
        for (auto* cl : run.clusters) {
            Json::Value c;
            c["uid"]   = cluster_uid(apa, cl->ident());
            c["ident"] = cl->ident();
            c["apa"]   = apa;
            Json::Value xs(Json::arrayValue), ys(Json::arrayValue), zs(Json::arrayValue), qs(Json::arrayValue);
            std::size_t np = 0;
            for (auto* blob : cl->children()) {
                const double qpt = blob->npoints() ? blob->charge() / blob->npoints() : 0.0;
                auto pts = blob->points("3d", {"x", "y", "z"});
                for (int i = 0; i < blob->npoints(); ++i) {
                    xs.append(pts.at(i).x() / cm);
                    ys.append(pts.at(i).y() / cm);
                    zs.append(pts.at(i).z() / cm);
                    qs.append(qpt);
                    ++np;
                }
            }
            c["npoints"] = (Json::UInt64)np;
            c["x"] = xs; c["y"] = ys; c["z"] = zs; c["q"] = qs;
            clusters.append(c);
        }

        // Final auto-selected bundles = those left in flash_bundles_map after the fit.
        std::set<const TimingTPCBundle*> selected;
        for (const auto& kv : run.flash_bundles_map)
            for (const auto& b : kv.second) selected.insert(b.get());

        // The candidate universe, restricted to TPC-contained bundles: an
        // uncontained bundle (cluster outside the box after the T0 x-shift)
        // carries zeroed metrics and a zero pred_flash, so it only clutters the
        // hand-scan. None are ever auto_selected (require_containment drops them
        // before the fit), so skipping them here loses no matcher decision.
        for (const auto& b : run.all_bundles) {
            if (!b->get_contained()) continue;
            // Same opaque-cathode mismatched-candidate filter the fit uses in
            // build_bundles, so the hand-scan tables show the same bundle universe:
            // keep same-side and cross-side cathode-crossers, drop other cross-side.
            if (cross_side_mismatch_drop(b->get_flash(), apa, b->get_flag_at_x_boundary()))
                continue;
            auto* fl = b->get_flash();
            Json::Value jb;
            jb["apa"]          = apa;
            jb["flash_gid"]    = (m_shared_flash ? 0 : apa) * kFlashGidStride
                                 + run.global_flash_idx_map.at(fl);
            jb["flash_id"]     = fl->get_flash_id();
            jb["main_cluster"] = cluster_uid(apa, b->get_main_cluster()->ident());
            Json::Value others(Json::arrayValue);
            for (auto* oc : b->get_other_clusters()) others.append(cluster_uid(apa, oc->ident()));
            jb["other_clusters"]   = others;
            jb["ks_dis"]           = b->get_ks_dis();
            jb["chi2"]             = b->get_chi2();
            jb["ndf"]              = b->get_ndf();
            jb["strength"]         = b->get_strength();
            jb["total_pred_light"] = b->get_total_pred_light();
            jb["total_PE"]         = fl->get_total_PE();
            // Containment verdict from compute_endpoint_flags. Uncontained bundles
            // (dropped by require_containment) keep the ctor's zeroed metrics and a
            // zero-filled pred_flash; the flag is what tells them apart.
            // Dump the FULL prediction: rail-flagged channels keep their
            // predicted light (they are only excluded from chi2/KS/LASSO);
            // identical to get_pred_flash when use_saturation_flag is off.
            const auto& pred = b->get_pred_flash_full();
            jb["contained"]            = b->get_contained();
            jb["consistent"]           = b->get_consistent_flag();
            jb["potential_bad_match"]  = b->get_potential_bad_match_flag();
            jb["close_to_PMT"]         = b->get_flag_close_to_PMT();
            jb["at_x_boundary"]        = b->get_flag_at_x_boundary();
            jb["at_cathode"]           = b->get_flag_at_cathode();
            jb["spec_end"]             = b->get_spec_end_flag();
            jb["window_truncated"]     = b->get_flag_window_truncated();
            jb["two_boundary"]         = b->get_flag_two_boundary();
            jb["xtpc_consistent"]      = b->get_flag_xtpc_consistent();
            jb["xtpc_scenario1"]       = b->get_flag_xtpc_scenario1();
            jb["xtpc_pin"]             = b->get_flag_xtpc_pin();
            // Key emitted only when the rescue knob is on => knob-off dumps stay
            // byte-identical. True = this bundle failed raw containment by a cathode
            // overshoot within tolerance and survived only via the scenario-1
            // cross-volume confirmation (purged provisionals never reach here: they
            // keep contained=false and are skipped above).
            if (m_xtpc_cathode_tol > 0.0)
                jb["xtpc_cathode_rescued"] = b->get_flag_xtpc_cathode_provisional();
            // Key emitted only when the relaxed rescue tier is on => knob-off dumps
            // stay byte-identical. True = adopted by the RELAXED second-chance
            // gates (low-confidence match; failed the tight cluster_rescue bar).
            if (m_cluster_rescue_relaxed)
                jb["cluster_rescue_relaxed"] = b->get_flag_cluster_rescue_relaxed();
            jb["auto_selected"]        = (bool)selected.count(b.get());
            Json::Value jpred(Json::arrayValue);
            for (double v : pred) jpred.append(v);
            jb["pred_pe"] = jpred;
            bundles.append(jb);
        }
    }

    // ±80 ns flash-flash coincidence groups across BOTH TPCs (replicates
    // MultiAlgBlobClustering::store_flash_groups on the merged flash times).
    {
        std::vector<std::pair<int, double>> ft;   // (array index into flashes, time_us)
        for (Json::ArrayIndex i = 0; i < flashes.size(); ++i)
            ft.emplace_back((int)i, flashes[i]["time"].asDouble());
        std::sort(ft.begin(), ft.end(), [&](const auto& a, const auto& b) {
            if (a.second != b.second) return a.second < b.second;
            return flashes[a.first]["gid"].asInt() < flashes[b.first]["gid"].asInt();
        });
        const double window_us = m_flash_group_window / us;
        int next_group = 0;
        double prev_t = 0;
        for (std::size_t k = 0; k < ft.size(); ++k) {
            if (k == 0 || (ft[k].second - prev_t) > window_us) ++next_group;
            flashes[ft[k].first]["group"] = next_group;
            prev_t = ft[k].second;
        }
    }

    top["geometry"] = geometry;
    top["flashes"]  = flashes;
    top["clusters"] = clusters;
    top["bundles"]  = bundles;

    std::string fname = m_calib_dump;
    if (fname.find('%') != std::string::npos) fname = String::format(fname, (int)m_count);
    Persist::dump(fname, top, /*pretty=*/false);
    log->debug("calib dump: wrote {} ({} bundles, {} flashes, {} clusters across {} APAs)",
               fname, bundles.size(), flashes.size(), clusters.size(), runs.size());
}

// ---------------------------------------------------------------------------
// Cathode-crossing TPC0/TPC1 offset diagnostic. For every TPC0 + TPC1 cluster
// pair that shares one T0 flash group and whose closest T0-corrected points are
// near the central cathode, log three vectors: the two local Hough track
// directions (dir0, dir1, at the cathode-end points -- offset-invariant, so they
// give the true track direction) and the connecting vector conn = p1 - p0 (the
// quantity the miscalibration distorts). Decompose conn into along-track (a
// benign sampling slide, since the two closest points need not share an
// arc-length on an angled track) and perpendicular, and split the gap into the
// drift x (t0/velocity, DEGENERATE) and transverse y-z (position shift / SCE).
// Observation-only: runs after matching, off unless cathode_diag is set.
void QLMatching::dump_cathode_diag(const std::vector<ApaRun>& runs)
{
    const double cm = units::cm;
    const double us = units::us;
    const double R  = m_cathode_diag_radius;     // Hough radius, internal length
    const double band = 15 * cm;                 // cathode region half-width (corrected x)
    const double dmax = 10 * cm;                 // closest-pair distance cut
    const double win_us = m_flash_group_window / us;

    // Matched clusters from both drift sides with their per-cluster T0 x-offset. Side
    // (low-x volume 0 / high-x volume 1) comes from the run's anode-vs-cathode position,
    // not the raw APA ident -- consistent with cull_cross_tpc's cross-side pairing.
    struct CC { Cluster* c; int side; double off; double t_us; };
    std::vector<CC> ccs;
    for (const auto& run : runs) {
        const int side = (run.anode_x < m_cathode_x) ? 0 : 1;
        for (auto* flash : flash_iter_order(run.flash_bundles_map)) {
            const double ft  = flash->get_time();
            const double off = run.sign_offset * (ft + trigger_offset_for(run.input_idx)) * drift_speed_for(run.input_idx);
            for (const auto& bundle : run.flash_bundles_map.at(flash)) {
                ccs.push_back({bundle->get_main_cluster(), side, off, ft / us});
                for (auto* oc : bundle->get_other_clusters())
                    ccs.push_back({oc, side, off, ft / us});
            }
        }
    }

    auto deg = [](const geo_vector_t& a, const geo_vector_t& b) {
        const double m = a.magnitude() * b.magnitude();
        if (m <= 0) return 0.0;
        double c = a.dot(b) / m;
        c = std::max(-1.0, std::min(1.0, c));
        return std::acos(c) * 180.0 / M_PI;
    };

    log->info("QLCATHODE header: count side0_id side1_id t0_us t1_us d_cm "
              "p0(x,y,z)_cm p1(x,y,z)_cm dir0 dir1 conn_cm "
              "ang(d0,d1) ang(d0,conn) ang(d1,conn) ang(conn,dhat)_deg "
              "along_cm perp_cm | dX_cm[DRIFT t0/v-DEGENERATE] dY_cm dZ_cm[transverse]");

    std::set<std::pair<Cluster*, Cluster*>> seen;
    for (size_t i = 0; i < ccs.size(); ++i) {
        if (ccs[i].side != 0) continue;
        for (size_t j = 0; j < ccs.size(); ++j) {
            if (ccs[j].side != 1) continue;
            if (std::abs(ccs[i].t_us - ccs[j].t_us) > win_us) continue;
            Cluster* c0 = ccs[i].c;
            Cluster* c1 = ccs[j].c;
            if (!seen.insert({c0, c1}).second) continue;   // dedupe pairs
            const double off0 = ccs[i].off, off1 = ccs[j].off;

            // Cathode-band point subsets in the T0-corrected frame.
            std::vector<int> s0, s1;
            for (int k = 0; k < c0->npoints(); ++k)
                if (std::abs(c0->point3d(k).x() + off0) < band) s0.push_back(k);
            for (int k = 0; k < c1->npoints(); ++k)
                if (std::abs(c1->point3d(k).x() + off1) < band) s1.push_back(k);
            if (s0.empty() || s1.empty()) continue;

            // Closest pair in corrected coords (brute force over the bands).
            double best = 1e18; int bi = -1, bj = -1;
            for (int a : s0) {
                const geo_point_t ra = c0->point3d(a);
                const double ax = ra.x() + off0;
                for (int b : s1) {
                    const geo_point_t rb = c1->point3d(b);
                    const double dx = ax - (rb.x() + off1);
                    const double dy = ra.y() - rb.y();
                    const double dz = ra.z() - rb.z();
                    const double d2 = dx * dx + dy * dy + dz * dz;
                    if (d2 < best) { best = d2; bi = a; bj = b; }
                }
            }
            const double d = std::sqrt(best);
            if (d > dmax) continue;

            const geo_point_t r0 = c0->point3d(bi), r1 = c1->point3d(bj);
            const geo_point_t p0(r0.x() + off0, r0.y(), r0.z());
            const geo_point_t p1(r1.x() + off1, r1.y(), r1.z());

            geo_vector_t dir0 = c0->vhough_transform(r0, R);   // raw point: offset-invariant
            geo_vector_t dir1 = c1->vhough_transform(r1, R);
            geo_vector_t conn = p1 - p0;
            if (dir0.dot(conn) < 0) dir0 = dir0 * -1.0;        // sign-fix the Hough axes
            if (dir1.dot(conn) < 0) dir1 = dir1 * -1.0;
            geo_vector_t dhat = dir0.norm() + dir1.norm();
            if (dhat.magnitude() > 0) dhat = dhat.norm();
            const double along = conn.dot(dhat);
            const geo_vector_t perp = conn - dhat * along;

            log->info("QLCATHODE {} 0/{} 1/{} {:.3f} {:.3f} d={:.3f} "
                      "p0=({:.3f},{:.3f},{:.3f}) p1=({:.3f},{:.3f},{:.3f}) "
                      "dir0=({:.4f},{:.4f},{:.4f}) dir1=({:.4f},{:.4f},{:.4f}) "
                      "conn=({:.3f},{:.3f},{:.3f}) "
                      "a_d0d1={:.2f} a_d0c={:.2f} a_d1c={:.2f} a_cdh={:.2f} "
                      "along={:.3f} perp={:.3f} dX={:.3f} dY={:.3f} dZ={:.3f}",
                      m_count, c0->ident(), c1->ident(),
                      ccs[i].t_us, ccs[j].t_us, d / cm,
                      p0.x() / cm, p0.y() / cm, p0.z() / cm,
                      p1.x() / cm, p1.y() / cm, p1.z() / cm,
                      dir0.x(), dir0.y(), dir0.z(), dir1.x(), dir1.y(), dir1.z(),
                      conn.x() / cm, conn.y() / cm, conn.z() / cm,
                      deg(dir0, dir1), deg(dir0, conn), deg(dir1, conn), deg(conn, dhat),
                      along / cm, perp.magnitude() / cm,
                      conn.x() / cm, conn.y() / cm, conn.z() / cm);
        }
    }
}

// One candidate cross-TPC pair test (the -cathode-diag geometry on the FULL clusters,
// with the per-TPC (y,z) pos_offset applied to the closest-point search and conn).
// Returns the scenario code: 1 = scenario 1 (closest approach < dmax, cathode ends
// present — tight, self-vetoing), 2 = scenario 2 (a half window-truncated AND
// conn,dir0,dir1 mutually collinear < angle_max AND d < dmax2), 0 = not consistent. A
// wrong-flash pairing carries the wrong T0 x-offset => large d => no match (self-vetoing).
// Whole-cluster + per-256-point-chunk bounding boxes over the raw (T0-independent)
// 3D points, for the box-pruned closest-approach in xtpc_pair_consistent.
QLMatching::XtpcBoxes QLMatching::xtpc_boxes(const WireCell::Clus::Facade::Cluster* c)
{
    static constexpr std::array<double, 6> empty_box{1e300, -1e300, 1e300, -1e300, 1e300, -1e300};
    XtpcBoxes bx;
    bx.whole = empty_box;
    const auto& P = c->points();
    const std::size_t n = P[0].size();
    bx.chunks.assign((n + XTPC_CHUNK - 1) / XTPC_CHUNK, empty_box);
    for (std::size_t k = 0; k < n; ++k) {
        auto& cb = bx.chunks[k / XTPC_CHUNK];
        cb[0] = std::min(cb[0], P[0][k]); cb[1] = std::max(cb[1], P[0][k]);
        cb[2] = std::min(cb[2], P[1][k]); cb[3] = std::max(cb[3], P[1][k]);
        cb[4] = std::min(cb[4], P[2][k]); cb[5] = std::max(cb[5], P[2][k]);
    }
    for (const auto& cb : bx.chunks) {
        bx.whole[0] = std::min(bx.whole[0], cb[0]); bx.whole[1] = std::max(bx.whole[1], cb[1]);
        bx.whole[2] = std::min(bx.whole[2], cb[2]); bx.whole[3] = std::max(bx.whole[3], cb[3]);
        bx.whole[4] = std::min(bx.whole[4], cb[4]); bx.whole[5] = std::max(bx.whole[5], cb[5]);
    }
    return bx;
}

int QLMatching::xtpc_pair_consistent(const XtpcMC& m0, const XtpcMC& m1,
                                     double* d_out, bool* pin_collinear_out,
                                     const XtpcBoxes* bx0, const XtpcBoxes* bx1) const
{
    const double R = m_xtpc_hough_radius;
    const double amax = m_xtpc_angle_max;
    if (pin_collinear_out) *pin_collinear_out = false;

    auto deg = [](const geo_vector_t& a, const geo_vector_t& b) {
        const double m = a.magnitude() * b.magnitude();
        if (m <= 0) return 180.0;
        double c = a.dot(b) / m;
        c = std::max(-1.0, std::min(1.0, c));
        return std::acos(c) * 180.0 / M_PI;
    };

    // Closest approach over the full clusters in the T0-corrected frame, with the
    // per-TPC transverse (y,z) pos_offset applied (mostly shifts this vector).
    // Hoist the scoped-view resolution out of this O(n0*n1) loop: point3d() (and the
    // npoints() loop bound) each re-resolve the 3D scope (Scope hash + hashtable lookup)
    // per call, which is ~all of the cross-TPC cull cost (the inner point3d(b) is
    // re-resolved n0*n1 times; see match/docs/qlmatching-perf-event60933.md). The raw
    // coordinate arrays are T0-independent, so bind them once; the per-point T0 +
    // pos_offset shift stays a scalar add, so the result (closest indices, iteration
    // order, first-min tie-break) is bit-identical.
    const auto& P0 = m0.c->points();   // {xs,ys,zs}, resolves the 3D scope once
    const auto& P1 = m1.c->points();
    const int n0 = (int)P0[0].size();
    const int n1 = (int)P1[0].size();
    double best = 1e30; int bi = -1, bj = -1;
    if (bx0 && bx1) {
        // Box-pruned exact closest approach. Same row-major point order as the
        // plain loop below, but (i) a whole 256-point chunk of m0 is skipped when
        // its box gap to m1's whole box cannot beat the running best, and (ii) a
        // single m0 point skips each m1 chunk its point-to-box gap cannot beat.
        // A pruned block contains only pairs with d2 >= best, which the strict
        // `<` update ignores anyway, and the first pair attaining each strictly
        // smaller distance is never inside a pruned block -- so best/bi/bj
        // (including FP-tie behavior) are bit-identical to the plain loop.
        auto gap1 = [](double lo0, double hi0, double lo1, double hi1) {
            return std::max(0.0, std::max(lo1 - hi0, lo0 - hi1));
        };
        const auto& w1 = bx1->whole;
        const int nc1 = (int)bx1->chunks.size();
        for (int ca = 0; ca * XTPC_CHUNK < n0; ++ca) {
            const auto& c0 = bx0->chunks[ca];
            const double wgx = gap1(c0[0] + m0.off, c0[1] + m0.off, w1[0] + m1.off, w1[1] + m1.off);
            const double wgy = gap1(c0[2] + m0.dy,  c0[3] + m0.dy,  w1[2] + m1.dy,  w1[3] + m1.dy);
            const double wgz = gap1(c0[4] + m0.dz,  c0[5] + m0.dz,  w1[4] + m1.dz,  w1[5] + m1.dz);
            if (wgx * wgx + wgy * wgy + wgz * wgz >= best) continue;
            const int a_end = std::min(n0, (ca + 1) * XTPC_CHUNK);
            for (int a = ca * XTPC_CHUNK; a < a_end; ++a) {
                const double ax = P0[0][a] + m0.off, ay = P0[1][a] + m0.dy, az = P0[2][a] + m0.dz;
                for (int cb = 0; cb < nc1; ++cb) {
                    const auto& c1 = bx1->chunks[cb];
                    const double gx = gap1(ax, ax, c1[0] + m1.off, c1[1] + m1.off);
                    const double gy = gap1(ay, ay, c1[2] + m1.dy,  c1[3] + m1.dy);
                    const double gz = gap1(az, az, c1[4] + m1.dz,  c1[5] + m1.dz);
                    if (gx * gx + gy * gy + gz * gz >= best) continue;
                    const int b_end = std::min(n1, (cb + 1) * XTPC_CHUNK);
                    for (int b = cb * XTPC_CHUNK; b < b_end; ++b) {
                        const double dx = ax - (P1[0][b] + m1.off);
                        const double dy = ay - (P1[1][b] + m1.dy);
                        const double dz = az - (P1[2][b] + m1.dz);
                        const double d2 = dx * dx + dy * dy + dz * dz;
                        if (d2 < best) { best = d2; bi = a; bj = b; }
                    }
                }
            }
        }
    }
    else {
        for (int a = 0; a < n0; ++a) {
            const double ax = P0[0][a] + m0.off, ay = P0[1][a] + m0.dy, az = P0[2][a] + m0.dz;
            for (int b = 0; b < n1; ++b) {
                const double dx = ax - (P1[0][b] + m1.off);
                const double dy = ay - (P1[1][b] + m1.dy);
                const double dz = az - (P1[2][b] + m1.dz);
                const double d2 = dx * dx + dy * dy + dz * dz;
                if (d2 < best) { best = d2; bi = a; bj = b; }
            }
        }
    }
    if (bi < 0) return false;
    const double d = std::sqrt(best);
    if (d_out) *d_out = d;

    // Three-vector collinearity at the closest pair (scenario 2). Hough axes are
    // computed on raw points (a constant T0 x-shift does not rotate them).
    const geo_point_t r0 = m0.c->point3d(bi), r1 = m1.c->point3d(bj);
    const geo_point_t p0(r0.x() + m0.off, r0.y() + m0.dy, r0.z() + m0.dz);
    const geo_point_t p1(r1.x() + m1.off, r1.y() + m1.dy, r1.z() + m1.dz);
    geo_vector_t dir0 = m0.c->vhough_transform(r0, R);
    geo_vector_t dir1 = m1.c->vhough_transform(r1, R);
    const geo_vector_t conn = p1 - p0;
    if (dir0.dot(conn) < 0) dir0 = dir0 * -1.0;
    if (dir1.dot(conn) < 0) dir1 = dir1 * -1.0;
    const double a01 = deg(dir0, dir1), a0c = deg(dir0, conn), a1c = deg(dir1, conn);

    // Scenario 1 (cathode end present): closest approach below xtpc_dmax. Scenario 2
    // (a half window-truncated, cathode end missing): the connecting vector collinear
    // with both local directions AND the halves no farther apart than xtpc_dmax2 -- the
    // distance ceiling is load-bearing in the pre-fit candidate population, where
    // far-apart collinear pairs (d>300cm) would otherwise sneak in (a real crosser's two
    // halves cannot be arbitrarily far apart; the true sc2 pairs sit at d<=264cm).
    const bool scenario1 = d < m_xtpc_dmax;
    const bool scenario2 = (m0.wt || m1.wt) && d < m_xtpc_dmax2 &&
                           a01 < amax && a0c < amax && a1c < amax;
    const bool pass = scenario1 || scenario2;

    // Joint-pin direction confirmation (only when xtpc_joint_pin): a scenario-1 pair is a
    // genuine collinear crosser if the two TRACK AXES agree within m_xtpc_pin_angle. The
    // connecting-vector angles (a0c,a1c) are deliberately NOT used -- a ~1cm inter-TPC
    // transverse shift makes the closest-point connector ~perpendicular even for true
    // crossers. Take the OR (min) of two axis estimates: the LOCAL vhough direction (a01
    // above, at the cathode-end closest point) and the GLOBAL cluster PCA axis. Fold both
    // to [0,90] (a Hough/PCA axis has no sign). The local estimate over-reads on straight
    // crossers whose ends fan/curl (so global rescues them); the global estimate is
    // meaningless on bent/messy clusters (so local rescues them).
    double a01_global = 180.0;
    if (m_xtpc_joint_pin) {
        const geo_vector_t g0 = m0.c->get_pca().axis.at(0);
        const geo_vector_t g1 = m1.c->get_pca().axis.at(0);
        a01_global = deg(g0, g1);
    }
    auto foldd = [](double a) { return std::min(a, 180.0 - a); };
    const double a01_pin = std::min(foldd(a01), foldd(a01_global));
    if (pin_collinear_out)
        *pin_collinear_out = scenario1 && a01_pin < m_xtpc_pin_angle;

    // Diagnostic: log EVERY coincident pair we evaluate (not only passes), with the
    // cut values, so a near-miss (geometry just outside dmax/angle) is visible.
    log->debug("QLXTPC pair 0/{} 1/{} d={:.2f}cm (dmax={:.1f} dmax2={:.1f}) "
               "a01={:.1f} a0c={:.1f} a1c={:.1f} (amax={:.1f}) trunc={} sc1={} sc2={} pass={} "
               "a01g={:.1f} a01pin={:.1f} (pinmax={:.1f}) pin={}",
               m0.c->ident(), m1.c->ident(), d / units::cm,
               m_xtpc_dmax / units::cm, m_xtpc_dmax2 / units::cm,
               a01, a0c, a1c, amax, (m0.wt || m1.wt), scenario1, scenario2, pass,
               foldd(a01_global), a01_pin, m_xtpc_pin_angle,
               (scenario1 && a01_pin < m_xtpc_pin_angle));
    // Scenario 1 takes precedence (it is the tight, self-vetoing match that drives the
    // xtpc-priority cull); scenario 2 is the looser collinear-truncated fallback.
    if (scenario1) return 1;
    if (scenario2) return 2;
    return 0;
}

// Cross-TPC cathode-crossing PRE-FIT cull. A cathode-crosser lights one coincident flash
// group and is reconstructed as two halves (one per TPC). Pair every candidate
// main-cluster bundle across the two TPCs whose flashes coincide; if the two halves are
// geometrically cross-TPC consistent (xtpc_pair_consistent), flag both bundles, then drop
// each marked cluster's non-consistent rivals so fewer bundles enter the LASSO. Tuned on
// the 10 SBND hand-scan data events (MC validation); purity-first.
void QLMatching::cull_cross_tpc(std::vector<ApaRun>& runs)
{
    const double win_us = m_flash_group_window / units::us;

    // Gather candidate main-cluster bundles from both APAs, each tagged by APA side and
    // flash time (for coincidence) and carrying its T0 x-offset + transverse offset.
    // This now runs on the FULL pre-cull bundle set (cull_inconsistent has not run yet),
    // so restrict to cathode-relevant bundles -- a cathode-end present (at_x_boundary,
    // scenario 1) OR a truncated half (window_truncated, scenario 2). A clean
    // non-crossing bundle carries neither and cannot be a crosser half, so dropping it
    // here keeps the O(n0*n1*npts^2) pairing bounded without losing any real pair.
    // Tag each candidate by its DRIFT SIDE (low-x volume = 0, high-x volume = 1),
    // derived from the run's anode position vs the cathode -- NOT the raw APA ident.
    // For SBND the two APAs sit on opposite sides of the cathode, so side == ident and
    // the pairing is unchanged; for a JOINT PDHD node side 0 = {APA0,APA2} and side 1 =
    // {APA1,APA3}, so the two cathode-crossing halves (each in its own drift volume) are
    // paired across the central cathode.
    struct Cand { XtpcMC mc; int side; double t_us; };
    std::vector<Cand> cands;
    for (auto& run : runs) {
        const int side = (run.anode_x < m_cathode_x) ? 0 : 1;
        for (auto* flash : flash_iter_order(run.flash_bundles_map)) {
            const double off = run.sign_offset * (flash->get_time() + trigger_offset_for(run.input_idx)) * drift_speed_for(run.input_idx);
            for (const auto& bundle : run.flash_bundles_map.at(flash)) {
                // xtpc cathode rescue (xtpc_cathode_tol > 0, else the flag is never
                // set): admit near/past-cathode halves that miss at_x_boundary.
                if (!bundle->get_flag_at_x_boundary() && !bundle->get_flag_window_truncated()
                    && !bundle->get_flag_xtpc_cathode_cand())
                    continue;
                cands.push_back({XtpcMC{bundle.get(), bundle->get_main_cluster(), off,
                                        run.dy, run.dz, bundle->get_flag_window_truncated()},
                                 side, flash->get_time() / units::us});
            }
        }
    }
    // Diagnostic: list every cross-TPC candidate (main-cluster bundle surviving prefit).
    for (const auto& c : cands)
        log->debug("QLXTPC cand side{} t={:.3f}us mainclus={} npts={} trunc={}",
                   c.side, c.t_us, c.mc.c->ident(), c.mc.c->npoints(), c.mc.wt);

    // AABB prefilter for xtpc_pair_consistent's brute-force O(npts0 x npts1)
    // closest-approach step, two levels, both exact (bound-then-exact => flags,
    // pins and the closest pair itself are bit-identical):
    //  1. whole-box gap: the axis-aligned box gap between the two candidates'
    //     T0-shifted point clouds is a LOWER bound on their true closest
    //     approach, so a coincident pair whose gap already exceeds the scenario
    //     ceiling (xtpc_dmax, widened to xtpc_dmax2 when a window-truncated half
    //     makes scenario 2 possible) cannot pass and is skipped outright. This
    //     kills wrong-flash pairings whose T0 x-offset is wildly wrong.
    //  2. chunk boxes: pairs that survive 1 run the exact loop pruned per
    //     256-point chunk against the running best (see xtpc_pair_consistent).
    // Raw-point boxes are T0-independent -> cached per cluster; the
    // per-candidate (off, dy, dz) shift is applied to the box bounds at use.
    std::map<const WireCell::Clus::Facade::Cluster*, XtpcBoxes> box_cache;
    auto boxes_of = [&box_cache](const WireCell::Clus::Facade::Cluster* c)
        -> const XtpcBoxes& {
        auto it = box_cache.find(c);
        if (it == box_cache.end()) it = box_cache.emplace(c, xtpc_boxes(c)).first;
        return it->second;
    };
    // Gap between two shifted 1-D intervals (0 when they overlap).
    auto axis_gap = [](double lo0, double hi0, double lo1, double hi1) {
        return std::max(0.0, std::max(lo1 - hi0, lo0 - hi1));
    };

    // Profiling: the pairing loop calls xtpc_pair_consistent, whose closest-approach
    // step is a brute-force O(npts0 x npts1) double loop over both clusters' 3D points.
    // Count coincident pairs actually evaluated and the total point-pair distance ops
    // they incur, so the wall (operator's `xtpc_cull`) can be read as cost/point-pair.
    const auto t_pair0 = wallclock::now();
    long long n_pairs_eval = 0, n_pairs_coincident = 0, n_pairs_boxskip = 0;
    double point_pairs = 0;  // double: can exceed 2^31
    // Joint-pin: collect direction-confirmed scenario-1 pairs (with their closest-approach
    // d and the two bundles) for the greedy min-d pin below. Empty unless m_xtpc_joint_pin.
    struct PinRec { TimingTPCBundle* b0; TimingTPCBundle* b1;
                    WireCell::Clus::Facade::Cluster* c0; WireCell::Clus::Facade::Cluster* c1;
                    double d; };
    std::vector<PinRec> pins;
    for (size_t i = 0; i < cands.size(); ++i) {
        if (cands[i].side != 0) continue;
        for (size_t j = 0; j < cands.size(); ++j) {
            if (cands[j].side != 1) continue;
            ++n_pairs_eval;
            if (std::abs(cands[i].t_us - cands[j].t_us) > win_us) continue;  // coincident
            // xtpc cathode rescue: never confirm two PROVISIONAL (uncontained) halves
            // with each other -- at least one half of the pair must be a fully
            // contained bundle. No-op when the knob is off (flags never set).
            if (cands[i].mc.b->get_flag_xtpc_cathode_provisional() &&
                cands[j].mc.b->get_flag_xtpc_cathode_provisional())
                continue;
            ++n_pairs_coincident;
            log->debug("QLXTPC coincident 0/{} (t={:.3f}us) 1/{} (t={:.3f}us)",
                       cands[i].mc.c->ident(), cands[i].t_us,
                       cands[j].mc.c->ident(), cands[j].t_us);
            // Level 1: whole-box gap vs the widest ceiling this pair could pass.
            const XtpcBoxes& bx0 = boxes_of(cands[i].mc.c);
            const XtpcBoxes& bx1 = boxes_of(cands[j].mc.c);
            {
                const auto& b0 = bx0.whole;
                const auto& b1 = bx1.whole;
                const auto& m0 = cands[i].mc;
                const auto& m1 = cands[j].mc;
                const double gx = axis_gap(b0[0] + m0.off, b0[1] + m0.off,
                                           b1[0] + m1.off, b1[1] + m1.off);
                const double gy = axis_gap(b0[2] + m0.dy, b0[3] + m0.dy,
                                           b1[2] + m1.dy, b1[3] + m1.dy);
                const double gz = axis_gap(b0[4] + m0.dz, b0[5] + m0.dz,
                                           b1[4] + m1.dz, b1[5] + m1.dz);
                const double thr = (m0.wt || m1.wt) ? std::max(m_xtpc_dmax, m_xtpc_dmax2)
                                                    : m_xtpc_dmax;
                if (gx * gx + gy * gy + gz * gz > thr * thr) {
                    ++n_pairs_boxskip;
                    continue;
                }
            }
            point_pairs += (double)cands[i].mc.c->npoints() * cands[j].mc.c->npoints();
            double pair_d = 0; bool pin_ok = false;
            const int scenario = xtpc_pair_consistent(cands[i].mc, cands[j].mc, &pair_d, &pin_ok,
                                                      &bx0, &bx1);
            if (scenario == 0) continue;
            // Doc-19 scenario light gate: geometry alone sets the xtpc flags at
            // EVERY coincident flash whose T0 offset makes the halves touch;
            // when the gate is on, a bundle acquires the flags (and with them
            // the cull_inconsistent xtpc-tier privileges) only if its own
            // light passes.  OFF (default) => always true, bit-identical.
            // Joint-pin candidacy below is deliberately untouched (pin quality
            // is handled by the xtpc_pin_min_strength floor in the fit).
            auto sc1_light_pass = [this](const TimingTPCBundle* b) {
                if (!m_xtpc_sc1_light_gate) return true;
                const double c2n = b->get_ndf() > 0 ? b->get_chi2() / b->get_ndf() : 1e9;
                return b->get_ks_dis() <= m_xtpc_sc1_ks_max && c2n <= m_xtpc_sc1_c2n_max;
            };
            if (sc1_light_pass(cands[i].mc.b)) {
                cands[i].mc.b->set_flag_xtpc_consistent(true);
                if (scenario == 1) cands[i].mc.b->set_flag_xtpc_scenario1(true);
            }
            if (sc1_light_pass(cands[j].mc.b)) {
                cands[j].mc.b->set_flag_xtpc_consistent(true);
                if (scenario == 1) cands[j].mc.b->set_flag_xtpc_scenario1(true);
            }
            if (m_xtpc_joint_pin && pin_ok)
                pins.push_back({cands[i].mc.b, cands[j].mc.b, cands[i].mc.c, cands[j].mc.c, pair_d});
        }
    }
    log->debug("QLtiming cull_cross_tpc: ncands {} pairs_eval {} pairs_coincident {} "
               "pairs_boxskip {} point_pairs {:.3g} pairing_loop {:.1f} ms",
               cands.size(), n_pairs_eval, n_pairs_coincident, n_pairs_boxskip,
               point_pairs, ms_since(t_pair0));

    // Joint-pin: greedily bind the direction-confirmed crosser pairs, tightest (smallest
    // closest-approach d) first, so each cluster is pinned at most once -- to the coincident
    // flash where its halves actually meet at the cathode. A pinned bundle gets flag_xtpc_pin;
    // cull_inconsistent then keeps ONLY pinned bundles for a pinned cluster, binding both
    // halves to one flash (and keeping an otherwise-culled dim-flash partner alive). Empty
    // unless m_xtpc_joint_pin, so the off path is bit-identical. Deterministic: sort by
    // (d, ident0, ident1, flash_id0, flash_id1).
    if (m_xtpc_joint_pin && !pins.empty()) {
        // Per cluster-PAIR, choose the FLASH by best combined light (min ks-sum), not by
        // closest approach: every confirmed pairing already has d<dmax and collinear axes, so
        // among the coincident flashes light only picks which same-time wall -- it can never
        // pull the match off the crosser time (a non-coincident bright flash is not a pairing
        // here). This avoids forcing a same-time tie onto a worse-light wall (run29107 evt983
        // gid55 ks=1.0 vs gid56 ks=0.06). Tie-break d, then flash-ids.
        auto kssum = [](const PinRec& p) { return p.b0->get_ks_dis() + p.b1->get_ks_dis(); };
        std::map<std::pair<WireCell::Clus::Facade::Cluster*, WireCell::Clus::Facade::Cluster*>,
                 PinRec> best;
        for (const auto& p : pins) {
            auto key = std::make_pair(p.c0, p.c1);
            auto it = best.find(key);
            if (it == best.end()) { best.emplace(key, p); continue; }
            const PinRec& q = it->second;
            const double kp = kssum(p), kq = kssum(q);
            const int pf0 = p.b0->get_flash()->get_flash_id(), qf0 = q.b0->get_flash()->get_flash_id();
            bool better = (kp != kq) ? (kp < kq)
                        : (p.d != q.d) ? (p.d < q.d)
                        : (pf0 != qf0) ? (pf0 < qf0)
                        : (p.b1->get_flash()->get_flash_id() < q.b1->get_flash()->get_flash_id());
            if (better) it->second = p;
        }
        // Greedy over cluster-pairs, tightest geometry (min d) first, each cluster pinned once.
        std::vector<PinRec> reps;
        for (auto& kv : best) reps.push_back(kv.second);
        std::stable_sort(reps.begin(), reps.end(), [](const PinRec& a, const PinRec& b) {
            if (a.d != b.d) return a.d < b.d;
            if (a.c0->ident() != b.c0->ident()) return a.c0->ident() < b.c0->ident();
            return a.c1->ident() < b.c1->ident();
        });
        std::set<WireCell::Clus::Facade::Cluster*> pinned;
        for (const auto& p : reps) {
            if (pinned.count(p.c0) || pinned.count(p.c1)) continue;
            p.b0->set_flag_xtpc_pin(true);
            p.b1->set_flag_xtpc_pin(true);
            pinned.insert(p.c0);
            pinned.insert(p.c1);
            log->debug("QLXTPCPIN pair 0/{} (flash {}) 1/{} (flash {}) d={:.2f}cm kssum={:.3f}",
                       p.c0->ident(), p.b0->get_flash()->get_flash_id(),
                       p.c1->ident(), p.b1->get_flash()->get_flash_id(),
                       p.d / units::cm, kssum(p));
        }
    }
    // The actual culling (dropping a flagged cluster's rivals) is done by the unified
    // cull_inconsistent, which operator() runs per-APA right after this flag pass.
}

// Flashes ordered by get_flash_id() (the stable input tensor row index).
std::vector<Opflash*> QLMatching::flash_iter_order(const FlashBundlesMap& m)
{
    std::vector<Opflash*> v;
    v.reserve(m.size());
    for (auto& kv : m) v.push_back(kv.first);
    std::sort(v.begin(), v.end(),
              [](Opflash* a, Opflash* b) { return a->get_flash_id() < b->get_flash_id(); });
    return v;
}

// Clusters ordered by their stable global index (length-sorted vector position).
std::vector<Cluster*> QLMatching::cluster_iter_order(
    const ClusterBundlesMap& m, const std::map<Cluster*, int>& idx)
{
    std::vector<Cluster*> v;
    v.reserve(m.size());
    for (auto& kv : m) v.push_back(kv.first);
    std::sort(v.begin(), v.end(), [&](Cluster* a, Cluster* b) {
        return idx.at(a) < idx.at(b);
    });
    return v;
}

// Significant extreme points of a cluster (PCA-axis + coordinate extremes),
// flattened and memoized. The cathode endpoint is the cathode-most of these.
const std::vector<WireCell::Point>&
QLMatching::cluster_extreme_points(Cluster* cluster) const
{
    auto it = m_extreme_cache.find(cluster);
    if (it != m_extreme_cache.end()) return it->second;
    std::vector<WireCell::Point> flat;
    for (const auto& grp : cluster->get_extreme_wcps()) {
        for (const auto& p : grp) flat.push_back(p);
    }
    return m_extreme_cache.emplace(cluster, std::move(flat)).first->second;
}

// ----- boundary-flag filling (port of ToyMatching.cxx::calculate_pred_pe ~176-290) -----
bool QLMatching::compute_endpoint_flags(TimingTPCBundle* bundle,
                                        Cluster* cluster,
                                        double flash_x_offset,
                                        const ApaRun& run) const
{
    const double s = run.s, anode_x = run.anode_x, u_cathode = run.u_cathode;
    // Collect the per-time-slice representative drift coordinate u and blob
    // count. The toolkit's time_blob_map is nested apa->face->time->BlobSet. A
    // PDHD drift group spans TWO APAs offset in Z (e.g. APA1 z[0,232]cm + APA3
    // z[232,462]cm share the +x side), and the two APAs share the drift geometry
    // (anode_x / u_cathode / s) -- so union the slices over EVERY anode the
    // cluster touches. Keying on a single representative anode (as the SBND-origin
    // code did) declared any cluster living wholly in the group's second APA
    // "not contained" (no blobs under that key), hiding it from the matcher's
    // boundary flags and the calib dump. This mirrors the build_bundles
    // PE-inclusion box, which already unions grouping_anodes; for single-anode
    // runs (SBND) the cluster has one key, so the union is identical. We then walk
    // in u-order from each drift end (anode = min u, cathode = max u). Walking by u
    // rather than the prototype's strict time order is a faithful adaptation: for a
    // normally drifting track time and u are monotonic, and the boundary tests are
    // themselves in u.
    // The time_blob_map walk below is flash-independent except for the constant
    // flash_x_offset shift of the drift coordinate, yet this function runs once
    // per (flash x group) bundle (~10k calls/event on PDVD) and the per-slice
    // Blob::points() fetch (string-keyed dataset lookup + full point-array copy,
    // read only at index 0) dominated build_bundles CPU (66% on 039253 evt
    // 49846, round-0 profile of pdvd/docs/15_pdvd-light-ql-perf.md). Scan once
    // per cluster into m_endpoint_slice_cache (raw x0 and tallies, traversal
    // order) and rebuild u per flash with the identical FP ops; the per-call
    // sort below is kept so any rounding-induced tie ordering matches the
    // legacy per-call behavior exactly.
    auto cit = m_endpoint_slice_cache.find(cluster);
    if (cit == m_endpoint_slice_cache.end()) {
        EndpointSliceCache ec;
        const auto& tbm = cluster->time_blob_map();

        // Raw readout-window truncation bounds (T0-independent, APA-agnostic).
        // The blob slice indices are raw ticks (SamplingHelpers writes
        // slice_index = islice->start()/tick). No flash_x_offset enters: this is
        // a property of the raw window, identical for both reversed-drift APAs.
        // Collected over all slices with blobs, independently of the u-walk data
        // (which keeps only slices with 3d points).
        for (const auto& [anode, faces] : tbm) {
            (void)anode;
            for (const auto& [face, slices] : faces) {
                for (const auto& [t, bset] : slices) {
                    if (bset.empty()) continue;
                    const int lo = t;                               // slice_index_min (raw tick)
                    const int hi = (*bset.begin())->slice_index_max();
                    if (!ec.have_ticks) { ec.min_tick = lo; ec.max_tick = hi; ec.have_ticks = true; }
                    else { if (lo < ec.min_tick) ec.min_tick = lo; if (hi > ec.max_tick) ec.max_tick = hi; }
                }
            }
        }

        for (const auto& [anode, faces] : tbm) {
            (void)anode;
            for (const auto& [face, slices] : faces) {
                for (const auto& [t, bset] : slices) {
                    if (bset.empty()) continue;
                    const Blob* b0 = *bset.begin();
                    auto pts = b0->points("3d", {"x", "y", "z"});
                    if (pts.empty()) continue;
                    // Slice point count (charge mass of the slice) for the robust-endpoint
                    // trim; npoints() is a cached int per blob, so this stays cheap.
                    int slice_npts = 0;
                    double slice_q = 0.0;
                    for (const Blob* b : bset) { slice_npts += b->npoints(); slice_q += b->charge(); }
                    ec.slices.push_back({ pts.at(0).x(), static_cast<int>(bset.size()),
                                          slice_npts, slice_q });
                }
            }
        }
        cit = m_endpoint_slice_cache.emplace(cluster, std::move(ec)).first;
    }
    const EndpointSliceCache& ec = cit->second;

    // Window-truncation flag from the cached raw-tick bounds. The window end is
    // the post-resample frame length, m_readout_window_ticks. Its default is an
    // SBND number (3427); a detector with a longer window (PDHD: 5999) would
    // mislabel every mid-drift cluster past tick 3427, so PDHD feeds the real
    // value (read from the SP frame by its run script) via the
    // readout_window_ticks config -- see cfg/.../pdhd/qlmatching.jsonnet.
    // (An all-empty time_blob_map leaves have_ticks false: no flag, as before.)
    if (ec.have_ticks) {
        const int win = m_readout_window_ticks;
        const bool truncated =
            (ec.min_tick - 0 <= m_window_edge_ticks) ||
            (win - ec.max_tick <= m_window_edge_ticks);
        bundle->set_flag_window_truncated(truncated);
    }

    struct SliceU { double u; int nblobs; int npts; double q; };
    std::vector<SliceU> sv;
    sv.reserve(ec.slices.size());
    for (const auto& sl : ec.slices) {
        const double x = sl.x0 + flash_x_offset;
        sv.push_back({ s * (x - anode_x), sl.nblobs, sl.npts, sl.q });
    }
    if (sv.empty()) return false;  // no blobs / no 3d points: cannot be contained

    std::sort(sv.begin(), sv.end(),
              [](const SliceU& a, const SliceU& b) { return a.u < b.u; });

    const double N = static_cast<double>(cluster->nchildren());
    double first_u = sv.front().u;   // anode-end sampling u  (prototype first_pos_x - offset)
    double last_u  = sv.back().u;    // cathode-end sampling u (prototype last_pos_x  - offset)
    bool flag_spec_end = false;

    const double anode_in   = m_anode_ext1;                // low_x_cut  + low_x_cut_ext1
    const double cathode_in = u_cathode + m_cathode_ext1;  // high_x_cut + high_x_cut_ext1

    // ---- anode-end trim: walk inward (increasing u) ----
    if (first_u <= anode_in && first_u > -120 * units::cm) {
        int n_slices_out = 0, n_blobs_out = 0, n_blobs_def_out = 0;
        double prev_u = first_u, cur_u = first_u;
        for (const auto& sl : sv) {
            cur_u = sl.u;
            if (cur_u > anode_in && (cur_u - prev_u) > 0.75 * units::cm) break;
            if (n_slices_out > 60) break;
            if (cur_u < anode_in) n_blobs_def_out += sl.nblobs;
            n_slices_out += 1;
            n_blobs_out  += sl.nblobs;
            prev_u = cur_u;
        }
        if (n_slices_out <= 36 && n_blobs_out < 0.05 * N) {
            first_u = cur_u;
            if (n_slices_out > 10 && std::abs(cur_u - prev_u) < 10 * units::cm) flag_spec_end = true;
        } else if (n_slices_out <= 60 && n_blobs_out < 0.06 * N && std::abs(cur_u - prev_u) > 10 * units::cm) {
            first_u = cur_u;
        } else if (n_slices_out <= 25 && n_blobs_out < 0.12 * N && std::abs(cur_u - prev_u) > 20 * units::cm) {
            first_u = cur_u;
        }
        if (n_blobs_def_out < 0.0015 * N && n_blobs_def_out > 0) first_u = 0.0;  // snap to anode
    }

    // ---- cathode-end trim: walk inward (decreasing u) ----
    if (last_u >= cathode_in && last_u < u_cathode + 120 * units::cm) {
        int n_slices_out = 0, n_blobs_out = 0, n_blobs_def_out = 0;
        double prev_u = last_u, cur_u = last_u;
        for (auto it = sv.rbegin(); it != sv.rend(); ++it) {
            cur_u = it->u;
            if (cur_u < cathode_in && std::abs(cur_u - prev_u) > 0.75 * units::cm) break;
            if (n_slices_out > 60) break;
            if (cur_u > cathode_in) n_blobs_def_out += it->nblobs;
            n_slices_out += 1;
            n_blobs_out  += it->nblobs;
            prev_u = cur_u;
        }
        if (n_slices_out <= 36 && n_blobs_out < 0.05 * N) {
            last_u = cur_u;
            if (n_slices_out > 10 && std::abs(cur_u - prev_u) < 10 * units::cm) flag_spec_end = true;
        } else if (n_slices_out <= 60 && n_blobs_out < 0.06 * N && std::abs(cur_u - prev_u) > 10 * units::cm) {
            last_u = cur_u;
        } else if (n_slices_out <= 25 && n_blobs_out < 0.12 * N && std::abs(cur_u - prev_u) > 20 * units::cm) {
            last_u = cur_u;
        }
        if (n_blobs_def_out < 0.0015 * N && n_blobs_def_out > 0) last_u = u_cathode;  // snap to cathode
    }

    // ---- robust-endpoint trim (default OFF, byte-identical when m_robust_endpoint_trim
    // is false) ----
    // The gap-based trims above only fire for an *isolated* deep straggle separated
    // from the body by a >0.75 cm gap. A thin off-axis overclustering tail that stays
    // within 0.75 cm of the body is never trimmed, so a 0.3%-of-points straggle drags
    // the endpoint past the detector edge and falsely fails containment. This pass is
    // gap-independent: at each end it sums the POINTS in the outer region (beyond the
    // in-edge) and, if that material is sparse -- below max(frac*cluster_points, count)
    // -- snaps the endpoint back to the first slice inside the edge. Point-mass (not
    // blob count) is the sparsity measure, so a dense genuine track-end (hundreds of
    // points past the edge => grossly mis-positioned t0) exceeds the allowance and is
    // left to fail containment, while the thin tail is removed.
    if (m_robust_endpoint_trim) {
        const double allow =
            std::max(m_robust_endpoint_frac * static_cast<double>(cluster->npoints()),
                     m_robust_endpoint_count);
        // Charge budget for the outer-straggle judge: low-charge overclustered satellites
        // (~0.1% of cluster charge here) slip past the point-count allowance but stay under a
        // charge fraction, while a dense genuine track-end exceeds it. Disabled when
        // m_robust_endpoint_charge_frac == 0 (byte-identical to the point-count-only judge).
        double q_total = 0.0;
        for (const auto& sl : sv) q_total += sl.q;
        const double q_allow = m_robust_endpoint_charge_frac * q_total;
        // Absolute per-point charge-DENSITY gate (size-independent, does NOT scale with
        // cluster size or tail length): only trim outer material that is charge-sparse
        // (diffuse overclustered satellites, ~150 q/pt here), never a charge-dense real
        // track tip running into the boundary (~2500-8800 q/pt). Disabled when
        // m_robust_endpoint_charge_abs <= 0 (then the fraction/count judge stands alone).
        auto sparse = [&](int pts_out, double q_out) {
            if (pts_out <= allow) return true;   // original point-count path, untouched
            if (m_robust_endpoint_charge_frac <= 0.0) return false;
            if (q_out > q_allow) return false;   // new charge-fraction path
            // Density gate applies ONLY to the charge-fraction path: trim charge-sparse
            // (overclustered-satellite) material, never a charge-dense real track tip.
            return m_robust_endpoint_charge_abs <= 0.0 ||
                   (pts_out > 0 && q_out / pts_out < m_robust_endpoint_charge_abs);
        };
        if (last_u >= cathode_in) {                         // cathode end (deepest)
            int pts_out = 0; double q_out = 0.0; double in_u = last_u; bool reached = false;
            // gap_detached: mirror of the anode-end detachment path below, enabled by
            // m_robust_endpoint_gap_cathode (default false => this block is exactly the
            // legacy sparse()-only test, byte-identical; PDHD, which sets
            // m_robust_endpoint_gap for its anode cases, is unaffected until it opts in).
            // The over-cathode material contains a gap wider than m_robust_endpoint_gap,
            // i.e. it is overclustered junk DETACHED from the contiguous body rather than
            // a real continuous track end. Density-blind, as at the anode. Walking from
            // the deepest slice inward, consecutive u DECREASE, so the gap is
            // (prev_u - it->u).
            bool gap_detached = false; double prev_u = sv.back().u;
            for (auto it = sv.rbegin(); it != sv.rend(); ++it) {
                if (it->u < cathode_in) { in_u = it->u; reached = true; break; }
                if (m_robust_endpoint_gap_cathode && m_robust_endpoint_gap > 0.0
                    && (prev_u - it->u) > m_robust_endpoint_gap)
                    gap_detached = true;
                pts_out += it->npts; q_out += it->q; prev_u = it->u;
            }
            // Gap path: trim detached junk whose charge stays within the gap budget,
            // regardless of its per-point density (which the charge_abs gate would
            // refuse). The fraction cap is the safety bound: a real over-cathode stretch
            // carrying more charge is left to fail -- that is what stops a wrong-T0
            // hypothesis from being rescued here.
            const bool gap_trim = gap_detached
                && m_robust_endpoint_gap_charge_frac > 0.0
                && q_out <= m_robust_endpoint_gap_charge_frac * q_total;
            if (reached && (gap_trim || sparse(pts_out, q_out))) last_u = in_u;
        }
        if (first_u <= anode_in) {                          // anode end (shallowest)
            int pts_out = 0; double q_out = 0.0; double in_u = first_u; bool reached = false;
            // Where the walk stops calling material "outside". By default anode_in,
            // which is m_anode_ext1_margin cm ABOVE the floor the containment gate below
            // actually uses (first_u > anode_in - m_anode_ext1_margin) -- so material
            // already inside the gate gets counted as outside, and a body starting in
            // that dead band is fed to the judges as if it were straggle. With
            // m_robust_endpoint_walk_to_floor the walk breaks at that same floor, so it
            // never swallows what the gate would have accepted. The judges are unchanged:
            // a body starting BELOW the floor still exceeds them and still fails, which
            // is the point. Default (anode_in) => byte-identical. See QLMatching.h.
            const double walk_in = m_robust_endpoint_walk_to_floor
                                   ? anode_in - m_anode_ext1_margin
                                   : anode_in;
            // gap_detached: the sub-anode material contains a gap wider than
            // m_robust_endpoint_gap, i.e. it is overclustered junk DETACHED from the
            // contiguous body rather than a real continuous track end. Density-blind.
            bool gap_detached = false; double prev_u = sv.front().u;
            for (const auto& sl : sv) {
                if (sl.u > walk_in) { in_u = sl.u; reached = true; break; }
                if (m_robust_endpoint_gap > 0.0 && (sl.u - prev_u) > m_robust_endpoint_gap)
                    gap_detached = true;
                pts_out += sl.npts; q_out += sl.q; prev_u = sl.u;
            }
            // Gap path: trim detached junk whose charge stays within the gap budget,
            // regardless of its per-point density (which the charge_abs gate would
            // refuse). The fraction cap is the safety bound: a real sub-anode stretch
            // carrying more charge is left to fail. Disabled (=> sparse()-only,
            // byte-identical) when m_robust_endpoint_gap <= 0.
            const bool gap_trim = gap_detached
                && m_robust_endpoint_gap_charge_frac > 0.0
                && q_out <= m_robust_endpoint_gap_charge_frac * q_total;
            if (reached && (gap_trim || sparse(pts_out, q_out))) first_u = in_u;
        }
    }

    // ---- flag block (prototype 272-290), guarded by the in-window check ----
    // This guard is the prototype's flag_good_bundle / TPC-containment gate
    // (ToyMatching.cxx 272-275): the (post-trim) endpoints must lie inside the
    // TPC drift box. Its value is returned so the caller can discard uncontained
    // bundles when m_require_containment.
    const bool contained =
        first_u > anode_in - m_anode_ext1_margin &&
        last_u  > 0.0 &&
        last_u  < cathode_in &&
        first_u < u_cathode;
    if (contained) {
        bundle->set_spec_end_flag(flag_spec_end);
        // Anode end inside the flag window => close to the PMTs (which sit at the
        // anode plane) and at the x-boundary. vd_surface_flags: only an anode that
        // actually hosts PDs (non-empty anode_pd_channels for this input; PDVD:
        // the bottom volume's PMTs, the top CRP has none) sets flag_close_to_PMT,
        // and the chi2 relaxation is restricted to those channels. The at-x-boundary
        // (T0-constraining crosser) semantics are unchanged either way.
        if (first_u <= m_anode_ext2 && first_u > anode_in - m_anode_ext1_margin) {
            if (!m_vd_surface_flags) {
                bundle->set_flag_close_to_PMT(true);
            }
            else if (run.input_idx < m_anode_pd_channels.size() &&
                     !m_anode_pd_channels[run.input_idx].empty()) {
                bundle->set_flag_close_to_PMT(true);
                bundle->add_relax_channels(m_anode_pd_channels[run.input_idx]);
            }
            bundle->set_flag_at_x_boundary(true);
        }
        // vd_surface_flags: activity close to a PD-bearing side wall (PDVD membrane
        // XAs at the two y walls) => close to those PDs. y is drift-invariant so the
        // verdict is flash-independent (memoized per cluster); it does NOT imply
        // at_x_boundary (a wall graze constrains no T0).
        if (m_vd_surface_flags &&
            (!m_pd_wall_channels_ylo.empty() || !m_pd_wall_channels_yhi.empty())) {
            const uint8_t wf = wall_proximity(cluster, run);
            if ((wf & 1) && !m_pd_wall_channels_ylo.empty()) {
                bundle->set_flag_close_to_PMT(true);
                bundle->add_relax_channels(m_pd_wall_channels_ylo);
            }
            if ((wf & 2) && !m_pd_wall_channels_yhi.empty()) {
                bundle->set_flag_close_to_PMT(true);
                bundle->add_relax_channels(m_pd_wall_channels_yhi);
            }
        }
        // Cathode end => at the x-boundary (no PMTs there). With a CPA fiducial
        // configured (SBND) this is a 3-D point-in-volume test: among the cluster's
        // significant extreme points (get_extreme_wcps, cached per cluster), take the
        // cathode-most (largest drift u; the flash offset is a uniform shift so it
        // does not change which extreme is cathode-most) and test it against the
        // structure-exclusion volume, x corrected by the flash T0 offset. Without a
        // fiducial it falls back to the original flat-cathode 1-D drift window.
        bool at_cathode = false;
        if (m_cathode_fv) {
            const auto& extremes = cluster_extreme_points(cluster);
            bool have = false;
            geo_point_t cathode_pt;
            double best_u = 0.0;
            for (const auto& p : extremes) {
                const double u = s * (p.x() - anode_x);  // offset-independent ranking
                if (!have || u > best_u) { have = true; best_u = u; cathode_pt = p; }
            }
            if (have) {
                at_cathode = m_cathode_fv->contained(WireCell::Point(
                    cathode_pt.x() + flash_x_offset, cathode_pt.y(), cathode_pt.z()));
            }
        } else {
            at_cathode = (last_u >= u_cathode + m_cathode_ext2 && last_u < cathode_in);
        }
        if (at_cathode) {
            bundle->set_flag_at_x_boundary(true);
            // Inert diagnostic (nothing in the matching path reads it): lets a PDVD
            // cathode-PD treatment be designed from the calib dump. Set for every
            // detector; harmless where unused.
            bundle->set_flag_at_cathode(true);
        }
    }

    // ---- xtpc cathode rescue admission (default 0 => byte-identical; see the
    // m_xtpc_cathode_tol block in QLMatching.h and
    // wcp-porting-img/pdvd/docs/qlmatch/16_pdvd-clus97-crosser-evt298567.md §10) ----
    // Flat-cathode window only: with a CPA fiducial configured (SBND) the at_cathode
    // geometry is 3-D and this 1-D tolerance is not meaningful.
    if (m_xtpc_cathode_tol > 0.0 && !m_cathode_fv) {
        if (!contained && m_require_containment) {
            // Provisional admission: containment failed ONLY by cathode overshoot,
            // and the junk-tolerant cathode endpoint sits within the tolerance.
            // Junk tolerance: overclustered satellites merged into the cluster can
            // drag the RAW extreme tens of cm past the gate while carrying a tiny
            // charge fraction (evt298567 uid 4000097: 2.6% of charge at up to
            // u=386 vs gate 338.91, real end 342.24). Walk from the deep end
            // discarding slices while the discarded charge stays within
            // m_xtpc_cathode_qfrac of the cluster total; the endpoint is the first
            // slice that cannot be discarded. qfrac <= 0 => the raw endpoint
            // (strict). The discarded material is NOT removed from anything -- this
            // endpoint exists only for this admission test.
            const bool others_pass =
                first_u > anode_in - m_anode_ext1_margin &&
                last_u  > 0.0 &&
                first_u < u_cathode;
            if (others_pass && last_u >= cathode_in) {
                double end_u = last_u;
                if (m_xtpc_cathode_qfrac > 0.0) {
                    double q_total = 0.0;
                    for (const auto& sl : sv) q_total += sl.q;
                    const double q_budget = m_xtpc_cathode_qfrac * q_total;
                    double q_disc = 0.0;
                    for (auto it = sv.rbegin(); it != sv.rend(); ++it) {
                        if (q_disc + it->q > q_budget) { end_u = it->u; break; }
                        q_disc += it->q;
                    }
                }
                if (end_u < cathode_in + m_xtpc_cathode_tol) {
                    bundle->set_flag_xtpc_cathode_provisional(true);
                    bundle->set_flag_xtpc_cathode_cand(true);
                    bundle->set_flag_at_cathode(true);   // inert diagnostic, accurate here
                }
            }
        }
        else if (contained &&
                 last_u >= u_cathode + m_cathode_ext2 - m_xtpc_cathode_tol &&
                 last_u < cathode_in) {
            // Candidate admission only (the short side of the cathode window): the
            // partner half of a displaced crossing ends BELOW u_cathode+cathode_ext2,
            // misses at_x_boundary, and would be invisible to cull_cross_tpc. Do NOT
            // set at_x_boundary itself -- it also feeds the ladder / cross-side /
            // LASSO-weight paths, which must stay legacy.
            bundle->set_flag_xtpc_cathode_cand(true);
        }
    }
    return contained;
}

// Wall-proximity verdict for vd_surface_flags (bit 0 = low-y wall, bit 1 = high-y
// wall): does the cluster's y extent reach within m_pd_wall_cushion of the active
// volume's y edge? Uses the cached significant extreme points (which include the
// coordinate extremes); flash-independent, memoized per cluster.
uint8_t QLMatching::wall_proximity(Cluster* cluster, const ApaRun& run) const
{
    auto it = m_wall_flag_cache.find(cluster);
    if (it != m_wall_flag_cache.end()) return it->second;
    double y_min = 1e300, y_max = -1e300;
    for (const auto& p : cluster_extreme_points(cluster)) {
        if (p.y() < y_min) y_min = p.y();
        if (p.y() > y_max) y_max = p.y();
    }
    uint8_t f = 0;
    if (y_min <= y_max) {   // non-empty
        if (y_min < run.y_lo + m_pd_wall_cushion) f |= 1;
        if (y_max > run.y_hi - m_pd_wall_cushion) f |= 2;
    }
    return m_wall_flag_cache.emplace(cluster, f).first->second;
}

// ----- diagnostic flag_two_boundary (both PCA ends at a detector edge) -----
void QLMatching::compute_two_boundary_flag(TimingTPCBundle* bundle,
                                           Cluster* main_cluster,
                                           double flash_x_offset,
                                           const ApaRun& run) const
{
    bundle->set_flag_two_boundary(false);
    if (!main_cluster) return;

    // The two main-PCA-axis endpoints (get_extreme_wcps groups [0]=high, [1]=low
    // projection) are intrinsic to the cluster and flash-independent, so compute
    // once and reuse across candidate flashes. Guard degenerate clusters.
    auto it = m_pca_endpoints_cache.find(main_cluster);
    if (it == m_pca_endpoints_cache.end()) {
        const auto groups = main_cluster->get_extreme_wcps();
        if (groups.size() < 2 || groups[0].empty() || groups[1].empty()) return;
        it = m_pca_endpoints_cache.emplace(
            main_cluster, std::make_pair(groups[0][0], groups[1][0])).first;
    }
    const WireCell::Point& p_high = it->second.first;
    const WireCell::Point& p_low  = it->second.second;

    // Nearest of the 6 per-APA active-volume faces to a point, and its distance.
    // flash_x_offset shifts x only (drift); y/z are drift-invariant. The point is
    // "at" that face iff the distance <= m_two_boundary_margin.
    const double m = m_two_boundary_margin;
    auto nearest_face = [&](const WireCell::Point& p, int& face) {
        const double u = run.s * (p.x() + flash_x_offset - run.anode_x);  // drift coord
        const double d[6] = {
            u,                  // 0 anode    (u = 0)
            run.u_cathode - u,  // 1 cathode
            p.y() - run.y_lo,   // 2 bottom   (-y)
            run.y_hi - p.y(),   // 3 top      (+y)
            p.z() - run.z_lo,   // 4 upstream (-z)
            run.z_hi - p.z(),   // 5 downstream (+z)
        };
        face = 0;
        double best = d[0];
        for (int i = 1; i < 6; ++i) if (d[i] < best) { best = d[i]; face = i; }
        return best;
    };

    int face_high = 0, face_low = 0;
    const bool high_at_edge = nearest_face(p_high, face_high) <= m;
    const bool low_at_edge  = nearest_face(p_low,  face_low)  <= m;
    // At least one end must touch an x-boundary face (0 anode / 1 cathode).
    // A pair of purely transverse (y/z) edges constrains drift x — hence T0 —
    // hardly at all, so we don't count it as a two-boundary cluster.
    const bool any_x_face = (face_high <= 1) || (face_low <= 1);
    // Both PCA ends at a boundary AND at two SEPARATE faces (the cluster enters
    // through one surface and exits through a different one), with at least one
    // of those being an x-boundary.
    bundle->set_flag_two_boundary(high_at_edge && low_at_edge
                                  && face_high != face_low && any_x_face);
}

// ----- bundle map maintenance -----
void QLMatching::remove_bundle_selection(TimingTPCBundleSelection to_be_removed,
                                         TimingTPCBundleSet& bundle_set)
{
    for (auto& b : to_be_removed) bundle_set.erase(b);
}

void QLMatching::remove_bundle_selection(
    TimingTPCBundleSelection to_be_removed,
    FlashBundlesMap& flash_bundles_map,
    ClusterBundlesMap& cluster_bundles_map,
    std::map<std::pair<Opflash*, Cluster*>, TimingTPCBundle::pointer>& flash_cluster_bundles_map)
{
    for (auto& rm_bundle : to_be_removed) {
        auto rm_flash = rm_bundle->get_flash();
        auto rm_cluster = rm_bundle->get_main_cluster();
        flash_cluster_bundles_map.erase(std::make_pair(rm_flash, rm_cluster));
        {
            auto it = flash_bundles_map.find(rm_flash);
            if (it != flash_bundles_map.end()) {
                auto& v = it->second;
                auto vit = std::find(v.begin(), v.end(), rm_bundle);
                if (vit != v.end()) v.erase(vit);
                if (v.empty()) flash_bundles_map.erase(it);
            }
        }
        {
            auto it = cluster_bundles_map.find(rm_cluster);
            if (it != cluster_bundles_map.end()) {
                auto& v = it->second;
                auto vit = std::find(v.begin(), v.end(), rm_bundle);
                if (vit != v.end()) v.erase(vit);
                if (v.empty()) cluster_bundles_map.erase(it);
            }
        }
    }
}

void QLMatching::organize_bundles(
    TimingTPCBundleSelection& results_bundles,
    std::map<std::pair<Opflash*, Cluster*>, TimingTPCBundle::pointer>& /*flash_cluster_bundles_map*/)
{
    log->debug("organizing bundles");
    std::map<Opflash*, TimingTPCBundleSelection> eval_flash_bundles_map;
    for (auto bundle : results_bundles) {
        auto flash = bundle->get_flash();
        if (!flash) continue;
        eval_flash_bundles_map[flash].push_back(bundle);
    }

    // Deterministic outer order: process flashes by get_flash_id() (mirrors
    // flash_iter_order in operator()). Iterating eval_flash_bundles_map directly
    // is pointer-address order; since the per-flash best-bundle pick and the
    // add_bundle merges mutate bundle state, that order would otherwise change
    // the final bundle membership build-to-build.
    std::vector<Opflash*> eval_flash_order;
    eval_flash_order.reserve(eval_flash_bundles_map.size());
    for (auto& kv : eval_flash_bundles_map) eval_flash_order.push_back(kv.first);
    std::sort(eval_flash_order.begin(), eval_flash_order.end(),
              [](Opflash* a, Opflash* b) { return a->get_flash_id() < b->get_flash_id(); });
    for (auto* flash : eval_flash_order) {
        auto& orig_bundles = eval_flash_bundles_map[flash];
        // Inner order: by cluster_index_id, so the best-bundle tie-break and the
        // merge-candidate sequence are data-derived, not insertion/pointer order.
        std::sort(orig_bundles.begin(), orig_bundles.end(),
                  [](const TimingTPCBundle::pointer& a, const TimingTPCBundle::pointer& b) {
                      return a->get_cluster_index_id() < b->get_cluster_index_id();
                  });
        TimingTPCBundleSelection to_be_removed;
        TimingTPCBundle* best_bundle = nullptr;
        double best_strength = 0;
        for (auto b : orig_bundles) {
            if (b->get_strength() > best_strength) {
                best_strength = b->get_strength();
                best_bundle = b.get();
            }
        }
        log->debug("best bundle strength {} for flash {}", best_strength, flash->get_flash_id());
        for (auto b : orig_bundles) {
            if (b.get() == best_bundle) continue;
            if (best_bundle->examine_bundle(b.get())) {
                best_bundle->add_bundle(b.get());
                to_be_removed.push_back(b);
            }
            else {
                to_be_removed.push_back(b);
            }
        }
        for (auto& rm : to_be_removed) {
            auto it = std::find(results_bundles.begin(), results_bundles.end(), rm);
            if (it != results_bundles.end()) results_bundles.erase(it);
        }
        to_be_removed.clear();

        if (best_bundle && flash->get_time() > m_beam_mintime &&
            flash->get_time() < m_beam_maxtime) {
            log->debug("after merge, meas pe {}, pred pe {}, ks_dis {}, chi2/ndf {}",
                       int(flash->get_total_PE() * 100) / 100.,
                       int(best_bundle->get_total_pred_light() * 100) / 100.,
                       int(best_bundle->get_ks_dis() * 1000) / 1000.,
                       int(best_bundle->get_chi2() / best_bundle->get_ndf() * 100) / 100.);
        }
    }

    TimingTPCBundleSelection to_be_removed;
    for (auto bundle : results_bundles) {
        auto flash = bundle->get_flash();
        if (!flash) continue;
        if (flash->get_time() < m_beam_mintime || flash->get_time() > m_beam_maxtime) {
            if (bundle->get_ks_dis() > m_outbeam_ks_max ||
                bundle->get_chi2() / bundle->get_ndf() > m_outbeam_chi2ndf_max) {
                to_be_removed.push_back(bundle);
                continue;
            }
            if (std::abs(flash->get_total_PE() - bundle->get_total_pred_light()) >
                m_outbeam_pe_frac * flash->get_total_PE()) {
                to_be_removed.push_back(bundle);
                continue;
            }
        }
    }
    for (auto& rm : to_be_removed) {
        auto it = std::find(results_bundles.begin(), results_bundles.end(), rm);
        if (it != results_bundles.end()) results_bundles.erase(it);
    }
}
