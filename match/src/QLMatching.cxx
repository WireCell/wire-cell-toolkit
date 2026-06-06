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
    if (cfg.isMember("root_pcs_to_merge") && cfg["root_pcs_to_merge"].isArray()) {
        m_root_pcs_to_merge.clear();
        for (const auto& jname : cfg["root_pcs_to_merge"]) m_root_pcs_to_merge.insert(jname.asString());
    }

    m_dv    = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());

    m_inpath        = get(cfg, "inpath", m_inpath);
    m_outpath       = get(cfg, "outpath", m_outpath);
    m_cluster_t0    = get(cfg, "cluster_t0", m_cluster_t0);
    m_semimodel_file = get(cfg, "semimodel_file", m_semimodel_file);
    m_calib_dump    = get(cfg, "calib_dump", m_calib_dump);
    m_flash_group_window = get(cfg, "flash_group_window", m_flash_group_window);
    m_cathode_diag  = get(cfg, "cathode_diag", m_cathode_diag);
    m_cathode_diag_radius = get(cfg, "cathode_diag_radius", m_cathode_diag_radius);
    m_xtpc_flag         = get(cfg, "xtpc_flag",         m_xtpc_flag);
    m_xtpc_dmax         = get(cfg, "xtpc_dmax",         m_xtpc_dmax);
    m_xtpc_dmax2        = get(cfg, "xtpc_dmax2",        m_xtpc_dmax2);
    m_xtpc_angle_max    = get(cfg, "xtpc_angle_max",    m_xtpc_angle_max);
    m_xtpc_hough_radius = get(cfg, "xtpc_hough_radius", m_xtpc_hough_radius);

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

    m_flash_minPE   = get(cfg, "flash_minPE",   m_flash_minPE);
    m_flash_mintime = get(cfg, "flash_mintime", m_flash_mintime);
    m_flash_maxtime = get(cfg, "flash_maxtime", m_flash_maxtime);
    m_beam_mintime  = get(cfg, "beam_mintime",  m_beam_mintime);
    m_beam_maxtime  = get(cfg, "beam_maxtime",  m_beam_maxtime);
    if (m_beamonly) {
        m_flash_mintime = m_beam_mintime;
        m_flash_maxtime = m_beam_maxtime;
    }

    m_QtoL = get(cfg, "QtoL", m_QtoL);
    m_strength_cutoff = get(cfg, "strength_cutoff", m_strength_cutoff);
    m_drift_speed = get(cfg, "drift_speed", m_drift_speed);
    m_nchan = get(cfg, "nchan", m_nchan);

    // §A active-volume cushions (see match/docs/improve_progress.md). The raw
    // bounds come from m_dv->inner_bounds; these adjust the effective windows.
    // Defaults are the MicroBooNE-convention values (NOT bit-identical to the old
    // SBND literals); override to recover the old bounds.
    m_anode_ext1   = get(cfg, "anode_ext1",   m_anode_ext1);
    m_anode_ext2   = get(cfg, "anode_ext2",   m_anode_ext2);
    m_cathode_ext1 = get(cfg, "cathode_ext1", m_cathode_ext1);
    m_cathode_ext2 = get(cfg, "cathode_ext2", m_cathode_ext2);
    m_y_cushion    = get(cfg, "y_cushion",    m_y_cushion);
    m_z_cushion    = get(cfg, "z_cushion",    m_z_cushion);
    m_two_boundary_margin = get(cfg, "two_boundary_margin", m_two_boundary_margin);

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
    m_flash_pe_threshold = get(cfg, "flash_pe_threshold", m_flash_pe_threshold);

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

    m_readout_window_ticks = get(cfg, "readout_window_ticks", m_readout_window_ticks);
    m_window_edge_ticks    = get(cfg, "window_edge_ticks",    m_window_edge_ticks);

    m_require_containment  = get(cfg, "require_containment",  m_require_containment);

    m_reject_overpred      = get(cfg, "reject_overpred",      m_reject_overpred);
    m_overpred_total_ratio = get(cfg, "overpred_total_ratio", m_overpred_total_ratio);
    m_overpred_maxch_ratio = get(cfg, "overpred_maxch_ratio", m_overpred_maxch_ratio);

    m_empty_rescue          = get(cfg, "empty_rescue",          m_empty_rescue);
    m_rescue_metric_max     = get(cfg, "rescue_metric_max",     m_rescue_metric_max);
    m_rescue_exponent       = get(cfg, "rescue_exponent",       m_rescue_exponent);
    m_rescue_boundary_weight = get(cfg, "rescue_boundary_weight", m_rescue_boundary_weight);

    // Optional CPA structure-exclusion fiducial (SBND). Empty => disabled, and the
    // cathode-end flag_at_x_boundary keeps the original flat-cathode 1-D test.
    const auto cathode_fv_tn = get<std::string>(cfg, "cathode_fiducial", std::string(""));
    if (!cathode_fv_tn.empty()) {
        m_cathode_fv = Factory::find_tn<IFiducial>(cathode_fv_tn);
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
    if (!vuv_cfg.isObject() || !vis_cfg.isObject() || !geom_cfg.isObject() || !opdets_cfg.isArray()) {
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

    m_semi_model = std::make_unique<SemiAnalyticalModel>(
        vuv_cfg, vis_cfg, geom, m_opdets, /*doReflectedLight=*/true);

    log->debug("QLMatching configured: nopdets={}, semimodel_file={}",
               m_opdets.size(), m_semimodel_file);
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
    cfg["nchan"]           = m_nchan;
    cfg["semimodel_file"]  = m_semimodel_file;
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
    cfg["flash_minPE"]     = m_flash_minPE;
    cfg["flash_mintime"]   = m_flash_mintime;
    cfg["flash_maxtime"]   = m_flash_maxtime;
    cfg["beam_mintime"]    = m_beam_mintime;
    cfg["beam_maxtime"]    = m_beam_maxtime;
    cfg["QtoL"]            = m_QtoL;
    cfg["strength_cutoff"] = m_strength_cutoff;
    cfg["drift_speed"]     = m_drift_speed;

    cfg["anode_ext1"]   = m_anode_ext1;
    cfg["anode_ext2"]   = m_anode_ext2;
    cfg["cathode_ext1"] = m_cathode_ext1;
    cfg["cathode_ext2"] = m_cathode_ext2;
    cfg["y_cushion"]    = m_y_cushion;
    cfg["z_cushion"]    = m_z_cushion;
    cfg["two_boundary_margin"] = m_two_boundary_margin;

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
    cfg["flash_pe_threshold"] = m_flash_pe_threshold;

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

    cfg["readout_window_ticks"] = m_readout_window_ticks;
    cfg["window_edge_ticks"]    = m_window_edge_ticks;
    cfg["require_containment"]  = m_require_containment;
    cfg["reject_overpred"]      = m_reject_overpred;
    cfg["overpred_total_ratio"] = m_overpred_total_ratio;
    cfg["empty_rescue"]          = m_empty_rescue;
    cfg["rescue_metric_max"]     = m_rescue_metric_max;
    cfg["rescue_exponent"]       = m_rescue_exponent;
    cfg["rescue_boundary_weight"] = m_rescue_boundary_weight;
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

    // Match each APA's input independently in its own fresh, isolated ApaRun (one
    // per input port). The per-APA masks, geometry and flash_gid stride are
    // unchanged, and the isolation keeps the pointer-keyed map iteration
    // deterministic across APAs. With multiplicity 1 this is exactly the
    // historical single-input matcher.
    const auto t_prefit0 = wallclock::now();
    std::vector<ApaRun> runs(invec.size());
    for (std::size_t k = 0; k < invec.size(); ++k) {
        ApaRun& run = runs[k];
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
    // Per-TPC consistency cull, now AFTER the cross-TPC flag pass so it can honour the
    // xtpc flags (xtpc-priority). When no xtpc flag is set (non-SBND / xtpc_flag:false)
    // this reproduces the historical cull_inconsistent exactly => bit-identical.
    for (std::size_t k = 0; k < invec.size(); ++k) cull_inconsistent(runs[k]);
    const double t_cull = ms_since(t_cull0);

    const auto t_fit0 = wallclock::now();
    for (std::size_t k = 0; k < invec.size(); ++k) run_one_apa_fit(runs[k]);
    const double t_fit = ms_since(t_fit0);

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

    // Hand-scan calibration dump (off unless calib_dump is set). Reads the finished
    // per-APA runs only; never perturbs the matching above.
    if (!m_calib_dump.empty()) dump_calib(runs);
    if (!m_cathode_diag.empty()) dump_cathode_diag(runs);

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
    const PEErr pe_err_model{m_pe_err_floor, m_pe_err_frac, m_pe_err_knee};
    for (const auto& ff : run.grouping->flashes()) {
        auto flash = std::make_shared<Opflash>(ff, m_flash_pe_threshold, m_nchan, pe_err_model);
        if (flash->get_time() < m_flash_mintime || flash->get_time() > m_flash_maxtime) continue;
        if (flash->get_total_PE() < m_flash_minPE) continue;
        run.flashes.push_back(flash);
    }

    run.grouping->set_anodes({run.anode});
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

// Active-volume geometry from the IDetectorVolumes service, in a per-TPC drift
// coordinate u = s*(x-anode_x) (u=0 at the anode, u=u_cathode at the cathode) so
// the prototype's single-TPC inequalities port directly and both reversed-drift
// SBND APAs share one set of cushions. Also reduces the OpDet mask to this TPC
// and caches the kept-channel index used to size the LASSO matrices.
void QLMatching::compute_geometry(ApaRun& run)
{
    const unsigned int tpc = run.anode->ident();
    run.sign_offset  = (tpc == 0) ? -1 : 1;

    const WirePlaneId wpid(WirePlaneLayer_t::kAllLayers, 0, static_cast<int>(tpc));
    const BoundingBox bb = m_dv->inner_bounds(wpid);
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
    log->debug("anode {} bbox x[{:.2f},{:.2f}] y[{:.2f},{:.2f}] z[{:.2f},{:.2f}] cm; "
               "anode_x {:.2f} cathode_x {:.2f} u_cathode {:.2f} cm s {}",
               tpc, x_lo / units::cm, x_hi / units::cm,
               bray.first.y() / units::cm, bray.second.y() / units::cm,
               bray.first.z() / units::cm, bray.second.z() / units::cm,
               run.anode_x / units::cm, cathode_x / units::cm, run.u_cathode / units::cm, run.s);

    // Bundle-quality thresholds forwarded to every TimingTPCBundle below.
    run.qp = BundleQualityParams{
        m_bundle_ks_merge_max, m_bundle_chi2ndf_merge_max, m_bundle_addmerge_exponent,
        m_highconsist_ks_max, m_highconsist_min_ndf, m_bundle_pe_ndf_knee,
        m_bundle_mask_ks, m_pe_err_floor, m_pe_err_frac, m_pe_err_knee, m_pe_err_on_pred,
        m_highconsist_ladder,
        m_hc_clean_ks, m_hc_clean_c2, m_hc_good_ks, m_hc_good_c2,
        m_hc_tb_ks, m_hc_tb_c2, m_hc_miss_ks, m_hc_miss_c2, m_hc_miss_min_ndf,
        m_chi2_relax, m_chi2_pmt_excess, m_chi2_pmt_ratio, m_chi2_pmt_inflate};

    // Reduce mask to OpDets on this TPC: an OpDet belongs to TPC 0 / TPC 1 if it
    // sits on the low / high side of the cathode plane.
    for (std::size_t idet = 0; idet < run.opdet_mask.size(); ++idet) {
        const bool low_side = m_opdets[idet].center.x() < m_cathode_x;
        if ((tpc == 0) && !low_side) run.opdet_mask[idet] = 0;
        if ((tpc == 1) &&  low_side) run.opdet_mask[idet] = 0;
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
    for (auto flash : run.flashes) {
        const auto flash_time = flash->get_time();
        const double flash_x_offset = run.sign_offset * flash_time * m_drift_speed;

        // per-flash mask (also catches simulated saturated PMTs in MC).
        std::vector<unsigned int> flash_opdet_mask = run.opdet_mask;
        for (std::size_t idet = 0; idet < std::size_t(flash->get_num_channels()); ++idet) {
            auto pe_det = flash->get_PE(idet);
            if ((flash->get_total_PE() > m_mc_saturation_pe) && (pe_det == 0) && (m_data == false))
                flash_opdet_mask[idet] = 0;
        }

        log->debug("flash time {} flash PE {} flash_x_offset {}",
                   int(flash_time) / 100.,
                   int(flash->get_total_PE() * 100) / 100.,
                   int(flash_x_offset * 100) / 100.);

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
                compute_endpoint_flags(bundle.get(), main_cluster, flash_x_offset, run.s, run.anode_x, run.u_cathode,
                                       static_cast<int>(run.anode->ident()));
            bundle->set_contained(contained);   // diagnostic only (calib dump)
            // flag_two_boundary (both PCA ends at a detector edge). The calib dump and the
            // high-consistent ladder's B3 branch consume it; compute it only when one of those
            // needs it — production with the ladder off and no dump stays bit-identical.
            if (!m_calib_dump.empty() || m_highconsist_ladder)
                compute_two_boundary_flag(bundle.get(), main_cluster, flash_x_offset, run);
            // Discard bundles whose cluster is not contained in the TPC box once
            // the flash T0 x-offset is applied. Off by default.
            if (m_require_containment && !contained) continue;

            const std::size_t nopdets = flash->get_num_channels();
            std::vector<double> pred_flash(nopdets, 0.0);

            // Predicted light is summed over the whole group (main + associated).
            std::vector<Cluster*> group_clusters{main_cluster};
            group_clusters.insert(group_clusters.end(), others.begin(), others.end());
            std::size_t npt = 0;
            for (auto* gc : group_clusters) npt += gc->npoints();
            total_pts += npt;
            const auto t_vis0 = wallclock::now();
            for (auto* gc : group_clusters) {
              for (auto blob : gc->children()) {
                run.total_charge_blob += blob->charge();
                const double q = blob->charge() / blob->npoints();
                auto points = blob->points("3d", {"x", "y", "z"});

                for (int i = 0; i < blob->npoints(); ++i) {
                    run.total_charge_point += q;
                    const double x = points.at(i).x() + flash_x_offset;
                    const double y = points.at(i).y();
                    const double z = points.at(i).z();

                    // PE-inclusion gate in the per-TPC drift coordinate u.
                    const double u = run.s * (x - run.anode_x);
                    if (u < m_anode_ext1 || u > run.u_cathode + m_cathode_ext1) continue;
                    if (y < run.y_lo || y > run.y_hi || z < run.z_lo || z > run.z_hi) continue;

                    // SemiAnalyticalModel expects positions in cm; blob points are mm.
                    const WireCell::Point xyz_cm(x / units::cm, y / units::cm, z / units::cm);
                    // Skip masked opdets inside the model: their predicted light is
                    // discarded by the flash_opdet_mask gate below, so evaluating the
                    // per-opdet solid-angle/Gaisser-Hillas correction for them is wasted
                    // (~half of the same-TPC opdets are masked). Bit-identical result.
                    m_semi_model->detectedDirectVisibilities(direct_visibilities, xyz_cm, &flash_opdet_mask);
                    m_semi_model->detectedReflectedVisibilities(reflected_visibilities, xyz_cm, &flash_opdet_mask);

                    for (std::size_t idet = 0; idet < nopdets; ++idet) {
                        if (flash_opdet_mask.at(idet) == 0) continue;
                        const auto dir_vis = direct_visibilities.at(idet);
                        const auto ref_vis = reflected_visibilities.at(idet);
                        const auto dir_eff = m_VUVEfficiency.at(idet);
                        const auto ref_eff = m_VISEfficiency.at(idet);
                        pred_flash.at(idet) +=
                            q * m_QtoL * dir_vis * dir_eff + q * m_QtoL * ref_vis * ref_eff;
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
                double tot_pred = 0.0, tot_meas = 0.0, max_pred = 0.0, meas_at_max = 0.0;
                for (std::size_t j = 0; j < pred_flash.size(); ++j) {
                    if (j >= mask.size() || mask[j] == 0) continue;
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
        bool has_sc1 = false, has_keep = false;
        for (auto& b : bundles) {
            if (b->get_flag_xtpc_scenario1()) has_sc1 = true;
            if (keep_flag(b)) has_keep = true;
        }
        if (has_sc1) {
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
    // Snapshot the full pre-LASSO candidate universe for the empty-flash rescue,
    // before either round prunes by strength (§I; only when enabled).
    if (m_empty_rescue) run.prefit_snapshot = run.flash_bundles_map;

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
    Ress::matrix_t P = Ress::matrix_t::Zero(run.nopdet * nflash, nbundle + nflash);
    Ress::vector_t MF = Ress::vector_t::Zero(ncluster + nflash);
    Ress::matrix_t PF = Ress::matrix_t::Zero(ncluster + nflash, nbundle + nflash);
    Ress::vector_t weights = Ress::vector_t::Zero(nbundle + nflash);

    std::vector<std::pair<Opflash*, Cluster*>> pairs;
    std::size_t i = 0, ik = 0;
    for (auto* flash : flashes_ordered) {
        auto& bundles = run.flash_bundles_map[flash];
        for (unsigned int j = 0; j < run.nopdet; ++j) {
            const int opdet_idx = run.opdet_idx_v.at(j);
            const double pe = flash->get_PE(opdet_idx);
            const double pe_err = std::sqrt(flash->get_PE(opdet_idx) + std::pow(flash->get_PE_err(opdet_idx), 2));
            M(i * run.nopdet + j) = pe / pe_err;
            P(i * run.nopdet + j, nbundle + i) = pe / pe_err;
        }
        for (auto bundle : bundles) {
            const auto& pred_flash = bundle->get_pred_flash();
            for (unsigned int j = 0; j < run.nopdet; ++j) {
                const int opdet_idx = run.opdet_idx_v.at(j);
                const double pred_pe = pred_flash.at(opdet_idx);
                const double pe_err = std::sqrt(flash->get_PE(opdet_idx) + std::pow(flash->get_PE_err(opdet_idx), 2));
                P(i * run.nopdet + j, pairs.size()) = pred_pe / pe_err;
            }
            pairs.emplace_back(flash, bundle->get_main_cluster());
            const auto meas_pe_tot = flash->get_total_PE();
            const auto pred_pe_tot = bundle->get_total_pred_light();
            const double base = (std::abs(pred_pe_tot - meas_pe_tot) > m_pe_mismatch_knee * meas_pe_tot)
                                    ? std::abs(pred_pe_tot - meas_pe_tot) / meas_pe_tot
                                    : m_pe_mismatch_floor;
            weights(ik++) = base * lasso_flag_factor(bundle);
        }
        PF(ncluster + i, nbundle + i) = 1. / delta_light;
        flash_idx_map[flash] = nbundle + i;
        ++i;
    }
    for (unsigned int k = 0; k < nflash; ++k) weights(nbundle + k) = m_bkg_weight;
    for (unsigned int k = 0; k < ncluster; ++k) MF(k) = 1. / delta_charge;
    for (std::size_t n = 0; n < pairs.size(); ++n) {
        PF(cluster_idx_map[pairs.at(n).second], n) = 1. / delta_charge;
    }

    Ress::matrix_t PT  = P.transpose();
    Ress::matrix_t PFT = PF.transpose();
    Ress::vector_t y = PT * M + PFT * MF;
    Ress::matrix_t X = PT * P + PFT * PF;
    Ress::vector_t initial = Ress::vector_t::Zero(nbundle + nflash);
    for (std::size_t n = 0; n < pairs.size(); ++n) initial(n) = 1.0;

    Ress::Params params;
    params.model = Ress::lasso;
    params.lambda = lambda;
    log->debug("solving (round 1)");
    const double t_build = ms_since(t_build0);
    const auto t_solve0 = wallclock::now();
    Ress::vector_t solution = Ress::solve(X, y, params, initial, weights);
    log->debug("QLtiming fit_round1: ident {} nbundle {} nflash {} nopdet {} matrix_build {:.1f} "
               "lasso_solve {:.1f} ms (X {}x{})",
               run.charge_ident, nbundle, nflash, run.nopdet, t_build, ms_since(t_solve0),
               (int)X.rows(), (int)X.cols());

    TimingTPCBundleSelection to_be_removed;
    int n = 0;
    for (auto* flash : flashes_ordered) {
        for (auto bundle : run.flash_bundles_map[flash]) {
            if (solution(n) <= m_strength_cutoff && !m_beamonly) to_be_removed.push_back(bundle);
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
    Ress::matrix_t P = Ress::matrix_t::Zero(run.nopdet * nflash, nbundle);
    Ress::vector_t MF = Ress::vector_t::Zero(ncluster);
    Ress::matrix_t PF = Ress::matrix_t::Zero(ncluster, nbundle);
    Ress::vector_t weights = Ress::vector_t::Zero(nbundle);
    std::vector<std::pair<Opflash*, Cluster*>> pairs;

    std::size_t i = 0, ik = 0;
    for (auto* flash : flashes_ordered) {
        auto& bundles = run.flash_bundles_map[flash];
        for (unsigned int j = 0; j < run.nopdet; ++j) {
            const int opdet_idx = run.opdet_idx_v.at(j);
            const double pe = flash->get_PE(opdet_idx);
            const double pe_err = std::sqrt(flash->get_PE(opdet_idx) + std::pow(flash->get_PE_err(opdet_idx), 2));
            M(i * run.nopdet + j) = pe / pe_err;
        }
        for (auto bundle : bundles) {
            const auto ks_dis = bundle->get_ks_dis();
            const auto& pred_flash = bundle->get_pred_flash();
            for (unsigned int j = 0; j < run.nopdet; ++j) {
                const int opdet_idx = run.opdet_idx_v.at(j);
                const double pred_pe = pred_flash.at(opdet_idx);
                const double pe_err = std::sqrt(flash->get_PE(opdet_idx) + std::pow(flash->get_PE_err(opdet_idx), 2));
                P(i * run.nopdet + j, pairs.size()) = pred_pe / pe_err;
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
        PF(cluster_idx_map[pairs.at(n).second], n) = 1. / delta_charge;
    }

    Ress::matrix_t PT  = P.transpose();
    Ress::matrix_t PFT = PF.transpose();
    Ress::vector_t y = PT * M + PFT * MF;
    Ress::matrix_t X = PT * P + PFT * PF;
    Ress::vector_t initial = Ress::vector_t::Zero(nbundle);
    for (std::size_t n = 0; n < pairs.size(); ++n) initial(n) = 1.0;

    Ress::Params params;
    params.model = Ress::lasso;
    params.lambda = lambda;
    log->debug("solving (round 2)");
    const double t_build = ms_since(t_build0);
    const auto t_solve0 = wallclock::now();
    Ress::vector_t solution = Ress::solve(X, y, params, initial, weights);
    log->debug("QLtiming fit_round2: ident {} nbundle {} nflash {} nopdet {} matrix_build {:.1f} "
               "lasso_solve {:.1f} ms (X {}x{})",
               run.charge_ident, nbundle, nflash, run.nopdet, t_build, ms_since(t_solve0),
               (int)X.rows(), (int)X.cols());

    TimingTPCBundleSelection to_be_removed;
    int n = 0;
    for (auto* flash : flashes_ordered) {
        for (auto bundle : run.flash_bundles_map[flash]) {
            bundle->set_strength(solution(n));
            if (!(solution(n) > m_strength_cutoff || m_beamonly)) to_be_removed.push_back(bundle);
            ++n;
        }
    }
    remove_bundle_selection(to_be_removed, run.flash_bundles_map, run.cluster_bundles_map,
                            run.flash_cluster_bundles_map);
    remove_bundle_selection(to_be_removed, run.pre_bundles);
    to_be_removed.clear();

    // Empty-flash light-quality rescue (§I; default OFF = bit-identical, SBND-on).
    // Uses the pre-LASSO snapshot captured at fit_round1 start, so it can reach the
    // strength-0 but light-good bundles both rounds pruned.
    if (m_empty_rescue) rescue_empty_flashes(run, run.prefit_snapshot);

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

// Stamp the matched flash/t0 onto each cluster (and its associated sub-clusters):
// cluster_t0, the "flash" row index (Cluster::get_flash()), a globally-unique
// matched_flash_gid (survives the all-APA merge), and the per-channel flashpred.
void QLMatching::apply_matched_t0s(ApaRun& run)
{
    for (auto* flash : flash_iter_order(run.flash_bundles_map)) {
        for (auto bundle : run.flash_bundles_map[flash]) {
            auto* cluster = bundle->get_main_cluster();
            const double t0 = flash->get_time() * units::ns;
            const int flash_gid = run.anode->ident() * kFlashGidStride + run.global_flash_idx_map.at(flash);
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
void QLMatching::write_opflash_pc(ApaRun& run)
{
    std::vector<int> op_gid, op_ch;
    std::vector<double> op_time, op_pe;
    for (std::size_t fi = 0; fi < run.flashes.size(); ++fi) {
        const auto& flash = run.flashes[fi];
        const int gid = run.anode->ident() * kFlashGidStride + static_cast<int>(fi);
        for (int ch = 0; ch < m_nchan; ++ch) {
            op_gid.push_back(gid);
            op_time.push_back(flash->get_time());
            op_ch.push_back(ch);
            op_pe.push_back(flash->get_PE(ch));
        }
    }
    run.grouping->put_pcarray<int>(op_gid, "gid", "opflash");
    run.grouping->put_pcarray<double>(op_time, "time", "opflash");
    run.grouping->put_pcarray<int>(op_ch, "ch", "opflash");
    run.grouping->put_pcarray<double>(op_pe, "pe", "opflash");
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
// and adds it to the (fixed-box) cluster points. Cluster idents are per-APA, so a
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
    top["count"]        = (Json::UInt64)m_count;
    top["charge_ident"] = runs.empty() ? 0 : runs.front().charge_ident;

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
        g["y_lo"]        = run.y_lo / cm;
        g["y_hi"]        = run.y_hi / cm;
        g["z_lo"]        = run.z_lo / cm;
        g["z_hi"]        = run.z_hi / cm;
        geometry[std::to_string(apa)] = g;

        // Per-flash measured light + coincidence group id (computed below).
        for (std::size_t fi = 0; fi < run.flashes.size(); ++fi) {
            const auto& fl = run.flashes[fi];
            Json::Value f;
            f["gid"]      = apa * kFlashGidStride + (int)fi;
            f["id"]       = fl->get_flash_id();
            f["apa"]      = apa;
            f["time"]     = fl->get_time() / us;     // us
            f["total_PE"] = fl->get_total_PE();
            Json::Value pe(Json::arrayValue), pe_err(Json::arrayValue);
            for (int ch = 0; ch < m_nchan; ++ch) { pe.append(fl->get_PE(ch)); pe_err.append(fl->get_PE_err(ch)); }
            f["pe"]     = pe;
            f["pe_err"] = pe_err;
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
            auto* fl = b->get_flash();
            Json::Value jb;
            jb["apa"]          = apa;
            jb["flash_gid"]    = apa * kFlashGidStride + run.global_flash_idx_map.at(fl);
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
            const auto& pred = b->get_pred_flash();
            jb["contained"]            = b->get_contained();
            jb["consistent"]           = b->get_consistent_flag();
            jb["potential_bad_match"]  = b->get_potential_bad_match_flag();
            jb["close_to_PMT"]         = b->get_flag_close_to_PMT();
            jb["at_x_boundary"]        = b->get_flag_at_x_boundary();
            jb["spec_end"]             = b->get_spec_end_flag();
            jb["window_truncated"]     = b->get_flag_window_truncated();
            jb["two_boundary"]         = b->get_flag_two_boundary();
            jb["xtpc_consistent"]      = b->get_flag_xtpc_consistent();
            jb["xtpc_scenario1"]       = b->get_flag_xtpc_scenario1();
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
    const double v  = m_drift_speed;             // internal units (length/time)
    const double R  = m_cathode_diag_radius;     // Hough radius, internal length
    const double band = 15 * cm;                 // cathode region half-width (corrected x)
    const double dmax = 10 * cm;                 // closest-pair distance cut
    const double win_us = m_flash_group_window / us;

    // Matched clusters from both APAs with their per-cluster T0 x-offset.
    struct CC { Cluster* c; int apa; double off; double t_us; };
    std::vector<CC> ccs;
    for (const auto& run : runs) {
        const int apa = (int)run.anode->ident();
        for (auto* flash : flash_iter_order(run.flash_bundles_map)) {
            const double ft  = flash->get_time();
            const double off = run.sign_offset * ft * v;
            for (const auto& bundle : run.flash_bundles_map.at(flash)) {
                ccs.push_back({bundle->get_main_cluster(), apa, off, ft / us});
                for (auto* oc : bundle->get_other_clusters())
                    ccs.push_back({oc, apa, off, ft / us});
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

    log->info("QLCATHODE header: count apa0/id0 apa1/id1 t0_us t1_us d_cm "
              "p0(x,y,z)_cm p1(x,y,z)_cm dir0 dir1 conn_cm "
              "ang(d0,d1) ang(d0,conn) ang(d1,conn) ang(conn,dhat)_deg "
              "along_cm perp_cm | dX_cm[DRIFT t0/v-DEGENERATE] dY_cm dZ_cm[transverse]");

    std::set<std::pair<Cluster*, Cluster*>> seen;
    for (size_t i = 0; i < ccs.size(); ++i) {
        if (ccs[i].apa != 0) continue;
        for (size_t j = 0; j < ccs.size(); ++j) {
            if (ccs[j].apa != 1) continue;
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
int QLMatching::xtpc_pair_consistent(const XtpcMC& m0, const XtpcMC& m1) const
{
    const double R = m_xtpc_hough_radius;
    const double amax = m_xtpc_angle_max;

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
    if (bi < 0) return false;
    const double d = std::sqrt(best);

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

    // Diagnostic: log EVERY coincident pair we evaluate (not only passes), with the
    // cut values, so a near-miss (geometry just outside dmax/angle) is visible.
    log->debug("QLXTPC pair 0/{} 1/{} d={:.2f}cm (dmax={:.1f} dmax2={:.1f}) "
               "a01={:.1f} a0c={:.1f} a1c={:.1f} (amax={:.1f}) trunc={} sc1={} sc2={} pass={}",
               m0.c->ident(), m1.c->ident(), d / units::cm,
               m_xtpc_dmax / units::cm, m_xtpc_dmax2 / units::cm,
               a01, a0c, a1c, amax, (m0.wt || m1.wt), scenario1, scenario2, pass);
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
    const double v = m_drift_speed;
    const double win_us = m_flash_group_window / units::us;

    // Gather candidate main-cluster bundles from both APAs, each tagged by APA side and
    // flash time (for coincidence) and carrying its T0 x-offset + transverse offset.
    // This now runs on the FULL pre-cull bundle set (cull_inconsistent has not run yet),
    // so restrict to cathode-relevant bundles -- a cathode-end present (at_x_boundary,
    // scenario 1) OR a truncated half (window_truncated, scenario 2). A clean
    // non-crossing bundle carries neither and cannot be a crosser half, so dropping it
    // here keeps the O(n0*n1*npts^2) pairing bounded without losing any real pair.
    struct Cand { XtpcMC mc; int apa; double t_us; };
    std::vector<Cand> cands;
    for (auto& run : runs) {
        const int apa = (int)run.anode->ident();
        for (auto* flash : flash_iter_order(run.flash_bundles_map)) {
            const double off = run.sign_offset * flash->get_time() * v;
            for (const auto& bundle : run.flash_bundles_map.at(flash)) {
                if (!bundle->get_flag_at_x_boundary() && !bundle->get_flag_window_truncated())
                    continue;
                cands.push_back({XtpcMC{bundle.get(), bundle->get_main_cluster(), off,
                                        run.dy, run.dz, bundle->get_flag_window_truncated()},
                                 apa, flash->get_time() / units::us});
            }
        }
    }
    // Diagnostic: list every cross-TPC candidate (main-cluster bundle surviving prefit).
    for (const auto& c : cands)
        log->debug("QLXTPC cand apa{} t={:.3f}us mainclus={} npts={} trunc={}",
                   c.apa, c.t_us, c.mc.c->ident(), c.mc.c->npoints(), c.mc.wt);

    // Profiling: the pairing loop calls xtpc_pair_consistent, whose closest-approach
    // step is a brute-force O(npts0 x npts1) double loop over both clusters' 3D points.
    // Count coincident pairs actually evaluated and the total point-pair distance ops
    // they incur, so the wall (operator's `xtpc_cull`) can be read as cost/point-pair.
    const auto t_pair0 = wallclock::now();
    long long n_pairs_eval = 0, n_pairs_coincident = 0;
    double point_pairs = 0;  // double: can exceed 2^31
    for (size_t i = 0; i < cands.size(); ++i) {
        if (cands[i].apa != 0) continue;
        for (size_t j = 0; j < cands.size(); ++j) {
            if (cands[j].apa != 1) continue;
            ++n_pairs_eval;
            if (std::abs(cands[i].t_us - cands[j].t_us) > win_us) continue;  // coincident
            ++n_pairs_coincident;
            log->debug("QLXTPC coincident 0/{} (t={:.3f}us) 1/{} (t={:.3f}us)",
                       cands[i].mc.c->ident(), cands[i].t_us,
                       cands[j].mc.c->ident(), cands[j].t_us);
            point_pairs += (double)cands[i].mc.c->npoints() * cands[j].mc.c->npoints();
            const int scenario = xtpc_pair_consistent(cands[i].mc, cands[j].mc);
            if (scenario == 0) continue;
            cands[i].mc.b->set_flag_xtpc_consistent(true);
            cands[j].mc.b->set_flag_xtpc_consistent(true);
            if (scenario == 1) {
                cands[i].mc.b->set_flag_xtpc_scenario1(true);
                cands[j].mc.b->set_flag_xtpc_scenario1(true);
            }
        }
    }
    log->debug("QLtiming cull_cross_tpc: ncands {} pairs_eval {} pairs_coincident {} "
               "point_pairs {:.3g} pairing_loop {:.1f} ms",
               cands.size(), n_pairs_eval, n_pairs_coincident, point_pairs, ms_since(t_pair0));
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
                                        double s, double anode_x, double u_cathode,
                                        int anode_ident) const
{
    // Collect the per-time-slice representative drift coordinate u and blob
    // count. The toolkit's time_blob_map is nested apa->face->time->BlobSet; the
    // matched anode has a single active face for SBND but we iterate all faces
    // under the apa for robustness. We then walk in u-order from each drift end
    // (anode = min u, cathode = max u). Walking by u rather than the prototype's
    // strict time order is a faithful adaptation: for a normally drifting track
    // time and u are monotonic, and the boundary tests are themselves in u.
    const auto& tbm = cluster->time_blob_map();
    auto ait = tbm.find(anode_ident);
    if (ait == tbm.end()) return false;  // no blobs under this anode: cannot be contained

    // Raw readout-window truncation flag (T0-independent, APA-agnostic). The
    // blob slice indices are raw ticks (SamplingHelpers writes slice_index =
    // islice->start()/tick), so a cluster is window-truncated if its leading
    // slice sits within m_window_edge_ticks of tick 0, or its trailing slice
    // within m_window_edge_ticks of the window end. No flash_x_offset enters:
    // this is a property of the raw window, identical for both reversed-drift
    // APAs. Computed independently of the u-walk below (which skips slices with
    // no 3d points and can early-return on sv.empty()). Always filled; inert (no
    // consumer reads it yet), so filling it leaves production output unchanged.
    {
        bool have = false;
        int min_tick = 0, max_tick = 0;
        for (const auto& [face, slices] : ait->second) {
            for (const auto& [t, bset] : slices) {
                if (bset.empty()) continue;
                const int lo = t;                               // slice_index_min (raw tick)
                const int hi = (*bset.begin())->slice_index_max();
                if (!have) { min_tick = lo; max_tick = hi; have = true; }
                else { if (lo < min_tick) min_tick = lo; if (hi > max_tick) max_tick = hi; }
            }
        }
        if (have) {
            const bool truncated =
                (min_tick - 0 <= m_window_edge_ticks) ||
                (m_readout_window_ticks - max_tick <= m_window_edge_ticks);
            bundle->set_flag_window_truncated(truncated);
        }
    }

    struct SliceU { double u; int nblobs; };
    std::vector<SliceU> sv;
    for (const auto& [face, slices] : ait->second) {
        for (const auto& [t, bset] : slices) {
            if (bset.empty()) continue;
            const Blob* b0 = *bset.begin();
            auto pts = b0->points("3d", {"x", "y", "z"});
            if (pts.empty()) continue;
            const double x = pts.at(0).x() + flash_x_offset;
            sv.push_back({ s * (x - anode_x), static_cast<int>(bset.size()) });
        }
    }
    if (sv.empty()) return false;  // no 3d points here: cannot be contained

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

    // ---- flag block (prototype 272-290), guarded by the in-window check ----
    // This guard is the prototype's flag_good_bundle / TPC-containment gate
    // (ToyMatching.cxx 272-275): the (post-trim) endpoints must lie inside the
    // TPC drift box. Its value is returned so the caller can discard uncontained
    // bundles when m_require_containment.
    const bool contained =
        first_u > anode_in - 1.0 * units::cm &&
        last_u  > 0.0 &&
        last_u  < cathode_in &&
        first_u < u_cathode;
    if (contained) {
        bundle->set_spec_end_flag(flag_spec_end);
        // Anode end inside the flag window => close to the PMTs (which sit at the
        // anode plane) and at the x-boundary.
        if (first_u <= m_anode_ext2 && first_u > anode_in - 1.0 * units::cm) {
            bundle->set_flag_close_to_PMT(true);
            bundle->set_flag_at_x_boundary(true);
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
        if (at_cathode) bundle->set_flag_at_x_boundary(true);
    }
    return contained;
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
