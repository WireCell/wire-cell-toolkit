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

    m_pe_err_floor      = get(cfg, "pe_err_floor",      m_pe_err_floor);
    m_pe_err_frac       = get(cfg, "pe_err_frac",       m_pe_err_frac);
    m_pe_err_knee       = get(cfg, "pe_err_knee",       m_pe_err_knee);
    m_flash_pe_threshold = get(cfg, "flash_pe_threshold", m_flash_pe_threshold);

    m_bundle_ks_merge_max      = get(cfg, "bundle_ks_merge_max",      m_bundle_ks_merge_max);
    m_bundle_chi2ndf_merge_max = get(cfg, "bundle_chi2ndf_merge_max", m_bundle_chi2ndf_merge_max);
    m_bundle_addmerge_exponent = get(cfg, "bundle_addmerge_exponent", m_bundle_addmerge_exponent);
    m_highconsist_ks_max       = get(cfg, "highconsist_ks_max",       m_highconsist_ks_max);
    m_highconsist_min_ndf      = get(cfg, "highconsist_min_ndf",      m_highconsist_min_ndf);
    m_bundle_pe_ndf_knee       = get(cfg, "bundle_pe_ndf_knee",       m_bundle_pe_ndf_knee);
    m_bundle_mask_ks           = get(cfg, "bundle_mask_ks",           m_bundle_mask_ks);

    m_readout_window_ticks = get(cfg, "readout_window_ticks", m_readout_window_ticks);
    m_window_edge_ticks    = get(cfg, "window_edge_ticks",    m_window_edge_ticks);

    m_require_containment  = get(cfg, "require_containment",  m_require_containment);

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
    cfg["nchan"]           = m_nchan;
    cfg["semimodel_file"]  = m_semimodel_file;
    cfg["pmts"]            = m_pmts;
    cfg["active_opdet_types"] = Json::arrayValue;
    for (int t : m_active_opdet_types) cfg["active_opdet_types"].append(t);
    cfg["data"]            = m_data;
    cfg["beamonly"]        = m_beamonly;
    cfg["ch_mask"]         = Json::arrayValue;
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

    cfg["pe_err_floor"]       = m_pe_err_floor;
    cfg["pe_err_frac"]        = m_pe_err_frac;
    cfg["pe_err_knee"]        = m_pe_err_knee;
    cfg["flash_pe_threshold"] = m_flash_pe_threshold;

    cfg["bundle_ks_merge_max"]      = m_bundle_ks_merge_max;
    cfg["bundle_chi2ndf_merge_max"] = m_bundle_chi2ndf_merge_max;
    cfg["bundle_addmerge_exponent"] = m_bundle_addmerge_exponent;
    cfg["highconsist_ks_max"]       = m_highconsist_ks_max;
    cfg["highconsist_min_ndf"]      = m_highconsist_min_ndf;
    cfg["bundle_pe_ndf_knee"]       = m_bundle_pe_ndf_knee;
    cfg["bundle_mask_ks"]           = m_bundle_mask_ks;

    cfg["readout_window_ticks"] = m_readout_window_ticks;
    cfg["window_edge_ticks"]    = m_window_edge_ticks;
    cfg["require_containment"]  = m_require_containment;
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

    m_extreme_cache.clear();  // cluster facades are per-event; drop stale endpoints

    // Match each APA's input independently in its own fresh, isolated ApaRun (one
    // per input port). The per-APA masks, geometry and flash_gid stride are
    // unchanged, and the isolation keeps the pointer-keyed map iteration
    // deterministic across APAs. With multiplicity 1 this is exactly the
    // historical single-input matcher.
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

        run_one_apa(run);

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

    // ---- Build / merge outputs ----
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

    ++m_count;
    return true;
}

// ============================================================================
// Per-APA matching pipeline. operator() builds one ApaRun (single-APA today)
// and calls run_one_apa(); the stages below were lifted verbatim from the old
// inline operator() body, reading config from the m_* members and threading all
// per-run state through ApaRun.
// ============================================================================

void QLMatching::run_one_apa(ApaRun& run)
{
    build_opdet_mask(run);
    read_flashes(run);
    decompose_cluster_groups(run);
    compute_geometry(run);
    build_bundles(run);
    build_bundle_maps(run);
    cull_inconsistent(run);
    fit_round1(run);
    fit_round2(run);
    apply_matched_t0s(run);
    write_opflash_pc(run);
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
        m_bundle_mask_ks};

    // Reduce mask to OpDets on this TPC: an OpDet belongs to TPC 0 / TPC 1 if it
    // sits on the low / high side of the cathode plane.
    for (std::size_t idet = 0; idet < run.opdet_mask.size(); ++idet) {
        const bool low_side = m_opdets[idet].center.x() < m_cathode_x;
        if ((tpc == 0) && !low_side) run.opdet_mask[idet] = 0;
        if ((tpc == 1) &&  low_side) run.opdet_mask[idet] = 0;
    }

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

// [Stage 1] Build a candidate (flash, cluster-group) bundle for every pair,
// predict its light (summed over the whole group), fill the boundary flags, and
// drop the KS==1 degenerate (no measured/predicted overlap). Surviving bundles
// land in run.pre_bundles.
void QLMatching::build_bundles(ApaRun& run)
{
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
                    std::vector<double> direct_visibilities;
                    m_semi_model->detectedDirectVisibilities(direct_visibilities, xyz_cm);
                    std::vector<double> reflected_visibilities;
                    m_semi_model->detectedReflectedVisibilities(reflected_visibilities, xyz_cm);

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
}

// Build the flash-keyed / cluster-keyed / (flash,cluster)-keyed bundle maps and
// impose a deterministic iteration order on them (pointer-keyed std::map/std::set
// would otherwise sort by heap address and permute the LASSO matrix columns).
void QLMatching::build_bundle_maps(ApaRun& run)
{
    for (auto bundle : run.pre_bundles) {
        auto flash   = bundle->get_flash();
        auto cluster = bundle->get_main_cluster();
        if (bundle->get_consistent_flag()) run.consistent_bundles.push_back(bundle);
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

// [Stage 1] For each cluster that has a consistent bundle, drop its other
// non-consistent bundles (they cannot win).
void QLMatching::cull_inconsistent(ApaRun& run)
{
    TimingTPCBundleSelection to_be_removed;
    for (auto good_bundle : run.consistent_bundles) {
        auto cluster = good_bundle->get_main_cluster();
        for (auto bundle : run.cluster_bundles_map[cluster]) {
            if (bundle == good_bundle) continue;
            if (bundle->get_consistent_flag()) continue;
            to_be_removed.push_back(bundle);
        }
    }
    remove_bundle_selection(to_be_removed, run.flash_bundles_map, run.cluster_bundles_map,
                            run.flash_cluster_bundles_map);
    remove_bundle_selection(to_be_removed, run.pre_bundles);
}

// [Stage 2] First LASSO round, with a per-flash background/light DOF column;
// post-fit prune bundles whose strength <= m_strength_cutoff.
void QLMatching::fit_round1(ApaRun& run)
{
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
            weights(ik++) = (std::abs(pred_pe_tot - meas_pe_tot) > m_pe_mismatch_knee * meas_pe_tot)
                              ? std::abs(pred_pe_tot - meas_pe_tot) / meas_pe_tot
                              : m_pe_mismatch_floor;
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
    Ress::vector_t solution = Ress::solve(X, y, params, initial, weights);

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
            weights(ik++) = base + delta_shape * run.nopdet * ks_dis / lambda;
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
    Ress::vector_t solution = Ress::solve(X, y, params, initial, weights);

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
            results_bundles.push_back(run.flash_cluster_bundles_map[std::make_pair(flash, cluster)]);
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
    top["quality_params"] = qp;

    // OpDet table (all channels). apa side and active flag mirror compute_geometry:
    // side by cathode-plane x; active iff the channel survives its own TPC's mask.
    std::vector<bool> active(m_nchan, false);
    for (const auto& run : runs) {
        for (int ch = 0; ch < m_nchan && ch < (int)run.opdet_mask.size(); ++ch)
            if (run.opdet_mask[ch] == 1) active[ch] = true;
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
