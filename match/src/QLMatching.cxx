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
#include "WireCellUtil/Ress.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Units.h"

#include <algorithm>

WIRECELL_FACTORY(QLMatching,
                 WireCell::Match::QLMatching,
                 WireCell::INamed,
                 WireCell::ITensorSetFilter,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Match;
using namespace WireCell::Clus::Facade;

// Stride used to build a globally-unique flash id (gid = anode_ident*stride +
// per-APA flash row) so the matched flash association survives the per-APA ->
// all-APA pctree merge and stays unambiguous across APAs. Far larger than any
// realistic per-APA flash count.
namespace { constexpr int kFlashGidStride = 1000000; }

QLMatching::QLMatching() : Aux::Logger("QLMatching", "match") {}
QLMatching::~QLMatching() = default;

void QLMatching::configure(const WireCell::Configuration& cfg)
{
    m_anode = Factory::find_tn<IAnodePlane>(cfg["anode"].asString());
    m_dv    = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());

    m_inpath        = get(cfg, "inpath", m_inpath);
    m_outpath       = get(cfg, "outpath", m_outpath);
    m_cluster_t0    = get(cfg, "cluster_t0", m_cluster_t0);
    m_semimodel_file = get(cfg, "semimodel_file", m_semimodel_file);

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
    m_drift_out_frac        = get(cfg, "drift_out_frac",        m_drift_out_frac);
    m_min_pred_pe           = get(cfg, "min_pred_pe",           m_min_pred_pe);
    m_preselect_chi2ndf_max = get(cfg, "preselect_chi2ndf_max", m_preselect_chi2ndf_max);

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

    if (cfg["VUVEfficiency"].isArray()) {
        m_VUVEfficiency.clear();
        for (const auto& v : cfg["VUVEfficiency"]) m_VUVEfficiency.push_back(v.asDouble());
    }
    if (cfg["VISEfficiency"].isArray()) {
        m_VISEfficiency.clear();
        for (const auto& v : cfg["VISEfficiency"]) m_VISEfficiency.push_back(v.asDouble());
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
    cfg["drift_out_frac"]        = m_drift_out_frac;
    cfg["min_pred_pe"]           = m_min_pred_pe;
    cfg["preselect_chi2ndf_max"] = m_preselect_chi2ndf_max;

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
    return cfg;
}

bool QLMatching::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    using WireCell::Clus::Facade::float_t;

    if (!in) {
        log->debug("EOS at call {}", m_count++);
        return true;
    }

    ExecMon em("starting QLMatching");

    // Per-channel OpDet on/off mask, derived from the injected OpDet table
    // (SemiAnalyticalModel::OpticalDetector::type; 1 = dome PMT, 0 = (X)Arapuca)
    // rather than a hard-coded SBND layout: channel i is on iff its type is in
    // m_active_opdet_types (default {1} => PMTs only, reproducing the historical
    // SBND mask). m_ch_mask then disables specific channels.
    std::vector<unsigned int> opdet_mask(m_opdets.size(), 0);
    if (m_pmts) {
        for (std::size_t i = 0; i < m_opdets.size(); ++i) {
            if (std::find(m_active_opdet_types.begin(), m_active_opdet_types.end(),
                          m_opdets[i].type) != m_active_opdet_types.end())
                opdet_mask[i] = 1;
        }
    }
    for (std::size_t i = 0; i < m_ch_mask.size(); ++i) opdet_mask[m_ch_mask[i]] = 0;

    // ---- Read inputs ----
    const auto& charge_ts = in;
    const int charge_ident = charge_ts->ident();
    std::string inpath = m_inpath;
    if (inpath.find("%") != std::string::npos) inpath = String::format(inpath, charge_ident);

    const auto& charge_tens = *charge_ts->tensors();
    log->debug("charge_tens.size {}", charge_tens.size());
    auto root_live = Aux::TensorDM::as_pctree(charge_tens, inpath + "/live");
    if (!root_live) {
        log->error("Failed to get point cloud tree from \"{}\"", inpath);
        return false;
    }
    log->debug("Got live pctree with {} children", root_live->nchildren());
    log->debug(em("got live pctree"));

    auto grouping = root_live->value.facade<Grouping>();

    // Flashes come from the canonical optical point clouds on the live root node
    // (written by Aux::FlashTensorToOpticalPCs, the same schema as
    // root/UbooneClusterSource) via the shared facade enumerator
    // Grouping::flashes() — which owns the flashlight-join walk. Each flash is
    // wrapped in an Opflash matching-adapter (time + dense per-channel PE, with
    // the matching-specific PE_err/fired synthesized in Opflash).
    const PEErr pe_err_model{m_pe_err_floor, m_pe_err_frac, m_pe_err_knee};
    std::vector<Opflash::pointer> flashes;
    for (const auto& ff : grouping->flashes()) {
        auto flash = std::make_shared<Opflash>(ff, m_flash_pe_threshold, m_nchan, pe_err_model);
        if (flash->get_time() < m_flash_mintime || flash->get_time() > m_flash_maxtime) continue;
        if (flash->get_total_PE() < m_flash_minPE) continue;
        flashes.push_back(flash);
    }

    grouping->set_anodes({m_anode});
    grouping->set_detector_volumes(m_dv);
    std::vector<Cluster*> clusters = grouping->children();
    std::sort(clusters.begin(), clusters.end(),
              [](const Cluster* a, const Cluster* b) { return a->get_length() > b->get_length(); });

    double total_charge_blob = 0.0;
    double total_charge_point = 0.0;
    double total_charge_blob_all = 0.0;
    for (auto cluster : clusters) {
        for (auto blob : cluster->children()) total_charge_blob_all += blob->charge();
    }

    std::for_each(clusters.begin(), clusters.end(),
                  [](Cluster* c) {
                      c->set_cluster_t0(-1e12);
                      c->set_scalar<int>("flash", -1);
                      c->set_scalar<int>("matched_flash_gid", -1);
                  });

    std::map<Opflash*, int>  global_flash_idx_map;
    std::map<Cluster*, int>  global_cluster_idx_map;
    for (std::size_t i = 0; i < flashes.size(); ++i) global_flash_idx_map[flashes[i].get()] = i;
    for (std::size_t i = 0; i < clusters.size(); ++i) global_cluster_idx_map[clusters[i]] = i;

    std::vector<TimingTPCBundle::pointer> all_bundles;
    TimingTPCBundleSet pre_bundles;

    const unsigned int tpc = m_anode->ident();
    const int sign_offset  = (tpc == 0) ? -1 : 1;

    // Active-volume geometry from the IDetectorVolumes service (replaces the old
    // SBND-specific x/y/z literals). inner_bounds() is the per-face sensitive
    // bounding box. The cathode is the bbox X corner nearest the cathode seam
    // (m_cathode_x, the same x=0 reference the OpDet split uses); the anode/PMT
    // plane is the far corner. We work in a per-TPC drift coordinate
    //   u = s * (x - anode_x)   with   u=0 at the anode, u=u_cathode (>0) at the
    // cathode, so the prototype's single-TPC inequalities port directly and the
    // two reversed-drift SBND APAs share one set of cushions.
    const WirePlaneId wpid(WirePlaneLayer_t::kAllLayers, 0, static_cast<int>(tpc));
    const BoundingBox bb = m_dv->inner_bounds(wpid);
    if (bb.empty()) {
        raise<ValueError>("QLMatching: empty detector-volume bounds for anode %d", tpc);
    }
    const Ray bray = bb.bounds();
    const double x_lo = bray.first.x(), x_hi = bray.second.x();
    const bool lo_is_cathode = std::abs(x_lo - m_cathode_x) < std::abs(x_hi - m_cathode_x);
    const double cathode_x = lo_is_cathode ? x_lo : x_hi;
    const double anode_x   = lo_is_cathode ? x_hi : x_lo;
    const double s         = (anode_x < cathode_x) ? +1.0 : -1.0;
    const double u_cathode = s * (cathode_x - anode_x);   // > 0
    // Y/Z active bounds from the same bbox, shrunk(+)/grown(-) by the cushions.
    const double y_lo = bray.first.y()  + m_y_cushion;
    const double y_hi = bray.second.y() - m_y_cushion;
    const double z_lo = bray.first.z()  + m_z_cushion;
    const double z_hi = bray.second.z() - m_z_cushion;
    log->debug("anode {} bbox x[{:.2f},{:.2f}] y[{:.2f},{:.2f}] z[{:.2f},{:.2f}] cm; "
               "anode_x {:.2f} cathode_x {:.2f} u_cathode {:.2f} cm s {}",
               tpc, x_lo / units::cm, x_hi / units::cm,
               bray.first.y() / units::cm, bray.second.y() / units::cm,
               bray.first.z() / units::cm, bray.second.z() / units::cm,
               anode_x / units::cm, cathode_x / units::cm, u_cathode / units::cm, s);

    // Bundle-quality thresholds forwarded to every TimingTPCBundle below.
    const BundleQualityParams qp{
        m_bundle_ks_merge_max, m_bundle_chi2ndf_merge_max, m_bundle_addmerge_exponent,
        m_highconsist_ks_max, m_highconsist_min_ndf, m_bundle_pe_ndf_knee};

    // Reduce mask to OpDets on this TPC: an OpDet belongs to TPC 0 / TPC 1 if it
    // sits on the low / high side of the cathode plane. Position-based, matching
    // the same-TPC test the optical model itself uses (center.x sign vs cathode);
    // reproduces the historical even/odd index split for the SBND PMTs and
    // generalizes via cathode_x. (The enclosing two-TPC-around-x=0 drift bounds
    // above are still SBND-specific; see improve_progress.md §A.)
    for (std::size_t idet = 0; idet < opdet_mask.size(); ++idet) {
        const bool low_side = m_opdets[idet].center.x() < m_cathode_x;
        if ((tpc == 0) && !low_side) opdet_mask[idet] = 0;
        if ((tpc == 1) &&  low_side) opdet_mask[idet] = 0;
    }

    for (auto flash : flashes) {
        const auto flash_time = flash->get_time();
        const double flash_x_offset = sign_offset * flash_time * m_drift_speed;

        // per-flash mask (also catches simulated saturated PMTs in MC).
        std::vector<unsigned int> flash_opdet_mask = opdet_mask;
        for (std::size_t idet = 0; idet < std::size_t(flash->get_num_channels()); ++idet) {
            auto pe_det = flash->get_PE(idet);
            if ((flash->get_total_PE() > m_mc_saturation_pe) && (pe_det == 0) && (m_data == false))
                flash_opdet_mask[idet] = 0;
        }

        log->debug("flash time {} flash PE {} flash_x_offset {}",
                   int(flash_time) / 100.,
                   int(flash->get_total_PE() * 100) / 100.,
                   int(flash_x_offset * 100) / 100.);

        for (std::size_t icluster = 0; icluster < clusters.size(); ++icluster) {
            Cluster* cluster = clusters[icluster];
            auto bundle = std::make_shared<TimingTPCBundle>(
                flash.get(), cluster, flash->get_flash_id(), icluster);
            bundle->set_quality_params(qp);
            all_bundles.push_back(bundle);
            bundle->set_opdet_mask(flash_opdet_mask);

            // Fill the boundary flags (close-to-PMT / at-x-boundary / spec-end)
            // for this (flash, cluster) pair from the cluster's drift endpoints.
            compute_endpoint_flags(bundle.get(), cluster, flash_x_offset, s, anode_x, u_cathode);

            const std::size_t nopdets = flash->get_num_channels();
            std::vector<double> pred_flash(nopdets, 0.0);

            std::size_t npt = cluster->npoints();
            int npt_outside_drift = 0;
            int npt_outside_bounds = 0;
            bool drifted_outside = false;

            for (auto blob : cluster->children()) {
                total_charge_blob += blob->charge();
                const double q = blob->charge() / blob->npoints();
                auto points = blob->points("3d", {"x", "y", "z"});

                for (int i = 0; i < blob->npoints(); ++i) {
                    total_charge_point += q;
                    const double x = points.at(i).x() + flash_x_offset;
                    const double y = points.at(i).y();
                    const double z = points.at(i).z();

                    // PE-inclusion gate in the per-TPC drift coordinate u. The
                    // window runs from just below the anode to just beyond the
                    // cathode (prototype low_x_cut+ext1 .. high_x_cut+ext1).
                    const double u = s * (x - anode_x);
                    if (u < m_anode_ext1 || u > u_cathode + m_cathode_ext1) { ++npt_outside_drift; continue; }
                    if (y < y_lo || y > y_hi || z < z_lo || z > z_hi) { ++npt_outside_bounds; continue; }

                    if (npt_outside_drift > m_drift_out_frac * npt) { drifted_outside = true; break; }

                    // SemiAnalyticalModel expects positions in cm. Blob points
                    // are in WCT units (mm).
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
                if (drifted_outside) break;
            }

            if (drifted_outside) {
                bundle->set_potential_bad_match_flag(true);
                continue;
            }
            bundle->set_pred_flash(pred_flash);
            if (bundle->get_total_pred_light() < m_min_pred_pe) continue;
            bundle->examine_bundle();
            if (bundle->get_ks_dis() == 1) {
                bundle->set_potential_bad_match_flag(true);
                continue;
            }
            if (bundle->get_chi2() / bundle->get_ndf() > m_preselect_chi2ndf_max) {
                bundle->set_potential_bad_match_flag(true);
                continue;
            }

            log->debug("initial eval: flash {} and cluster {}, meas PE {}, pred PE {}, npts {}, "
                       "ks_dis {}, chi2/ndf {}, ndf {}",
                       flash->get_flash_id(),
                       global_cluster_idx_map[cluster],
                       int(flash->get_total_PE() * 100) / 100.,
                       int(bundle->get_total_pred_light() * 100) / 100.,
                       npt,
                       int(bundle->get_ks_dis() * 1000) / 1000.,
                       int(bundle->get_chi2() / bundle->get_ndf() * 100) / 100.,
                       bundle->get_ndf());

            pre_bundles.insert(bundle);
        } // cluster loop
    }     // flash loop
    log->debug("n preselected bundles: {}", pre_bundles.size());

    // ---- Build maps ----
    FlashBundlesMap flash_bundles_map;
    ClusterBundlesMap cluster_bundles_map;
    std::map<std::pair<Opflash*, Cluster*>, TimingTPCBundle::pointer> flash_cluster_bundles_map;
    std::vector<TimingTPCBundle::pointer> consistent_bundles;

    for (auto bundle : pre_bundles) {
        auto flash   = bundle->get_flash();
        auto cluster = bundle->get_main_cluster();
        if (bundle->get_consistent_flag()) consistent_bundles.push_back(bundle);
        flash_bundles_map[flash].push_back(bundle);
        cluster_bundles_map[cluster].push_back(bundle);
        flash_cluster_bundles_map[std::make_pair(flash, cluster)] = bundle;
    }

    // Deterministic iteration order over flashes/clusters/bundles.
    //
    // Without these, the LASSO matrix column / row order depends on heap
    // allocator ordering of Opflash* / Cluster* / shared_ptr addresses,
    // because std::map<Pointer*, ...> and std::set<shared_ptr> sort by
    // pointer value. Two runs with identical inputs then permute matrix
    // columns and produce slightly different solution() vectors --- enough
    // to flip bundles across the m_strength_cutoff threshold and lose
    // run-to-run reproducibility.
    //
    // Outer order: flashes by flash_id (set from the input tensor row index,
    // stable). Clusters by their global index (set from the length-sorted
    // 'clusters' vector, stable). Inner: bundles within a flash sorted by
    // cluster_index_id, again the global index from the sorted vector.
    auto sort_inner_by_cluster_idx = [](FlashBundlesMap& m) {
        for (auto& kv : m) {
            std::sort(kv.second.begin(), kv.second.end(),
                      [](const TimingTPCBundle::pointer& a,
                         const TimingTPCBundle::pointer& b) {
                          return a->get_cluster_index_id() < b->get_cluster_index_id();
                      });
        }
    };
    auto flash_iter_order = [](const FlashBundlesMap& m) {
        std::vector<Opflash*> v;
        v.reserve(m.size());
        for (auto& kv : m) v.push_back(kv.first);
        std::sort(v.begin(), v.end(),
                  [](Opflash* a, Opflash* b) { return a->get_flash_id() < b->get_flash_id(); });
        return v;
    };
    auto cluster_iter_order = [&global_cluster_idx_map](const ClusterBundlesMap& m) {
        std::vector<Cluster*> v;
        v.reserve(m.size());
        for (auto& kv : m) v.push_back(kv.first);
        std::sort(v.begin(), v.end(), [&](Cluster* a, Cluster* b) {
            return global_cluster_idx_map.at(a) < global_cluster_idx_map.at(b);
        });
        return v;
    };
    sort_inner_by_cluster_idx(flash_bundles_map);

    TimingTPCBundleSelection to_be_removed;
    for (auto good_bundle : consistent_bundles) {
        auto cluster = good_bundle->get_main_cluster();
        for (auto bundle : cluster_bundles_map[cluster]) {
            if (bundle == good_bundle) continue;
            if (bundle->get_consistent_flag()) continue;
            to_be_removed.push_back(bundle);
        }
    }
    remove_bundle_selection(to_be_removed, flash_bundles_map, cluster_bundles_map,
                            flash_cluster_bundles_map);
    remove_bundle_selection(to_be_removed, pre_bundles);
    to_be_removed.clear();

    const double lambda       = m_lasso_lambda;
    const double delta_charge = m_delta_charge;
    const double delta_light  = m_delta_light;
    const double delta_shape  = m_delta_shape;

    unsigned int nopdet = 0;
    std::vector<int> opdet_idx_v;
    for (std::size_t idet = 0; idet < opdet_mask.size(); ++idet) {
        if (opdet_mask.at(idet) == 1) {
            opdet_idx_v.push_back(int(idet));
            ++nopdet;
        }
    }
    log->debug("nopdet {} opdet_idx_v size {}", nopdet, opdet_idx_v.size());

    // ---- First matching round ----
    {
        const unsigned int nbundle  = pre_bundles.size();
        const unsigned int nflash   = flash_bundles_map.size();
        const unsigned int ncluster = cluster_bundles_map.size();

        auto flashes_ordered  = flash_iter_order(flash_bundles_map);
        auto clusters_ordered = cluster_iter_order(cluster_bundles_map);
        std::map<Opflash*, int> flash_idx_map;
        std::map<Cluster*, int> cluster_idx_map;
        int idx = 0;
        for (auto* c : clusters_ordered) cluster_idx_map[c] = idx++;
        idx = 0;
        for (auto* f : flashes_ordered)  flash_idx_map[f]   = idx++;

        Ress::vector_t M = Ress::vector_t::Zero(nopdet * nflash);
        Ress::matrix_t P = Ress::matrix_t::Zero(nopdet * nflash, nbundle + nflash);
        Ress::vector_t MF = Ress::vector_t::Zero(ncluster + nflash);
        Ress::matrix_t PF = Ress::matrix_t::Zero(ncluster + nflash, nbundle + nflash);
        Ress::vector_t weights = Ress::vector_t::Zero(nbundle + nflash);

        std::vector<std::pair<Opflash*, Cluster*>> pairs;
        std::size_t i = 0, ik = 0;
        for (auto* flash : flashes_ordered) {
            auto& bundles = flash_bundles_map[flash];
            for (unsigned int j = 0; j < nopdet; ++j) {
                const int opdet_idx = opdet_idx_v.at(j);
                const double pe = flash->get_PE(opdet_idx);
                const double pe_err = std::sqrt(flash->get_PE(opdet_idx) + std::pow(flash->get_PE_err(opdet_idx), 2));
                M(i * nopdet + j) = pe / pe_err;
                P(i * nopdet + j, nbundle + i) = pe / pe_err;
            }
            for (auto bundle : bundles) {
                const auto& pred_flash = bundle->get_pred_flash();
                for (unsigned int j = 0; j < nopdet; ++j) {
                    const int opdet_idx = opdet_idx_v.at(j);
                    const double pred_pe = pred_flash.at(opdet_idx);
                    const double pe_err = std::sqrt(flash->get_PE(opdet_idx) + std::pow(flash->get_PE_err(opdet_idx), 2));
                    P(i * nopdet + j, pairs.size()) = pred_pe / pe_err;
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

        int n = 0;
        for (auto* flash : flashes_ordered) {
            for (auto bundle : flash_bundles_map[flash]) {
                if (solution(n) <= m_strength_cutoff && !m_beamonly) to_be_removed.push_back(bundle);
                ++n;
            }
        }
        remove_bundle_selection(to_be_removed, flash_bundles_map, cluster_bundles_map,
                                flash_cluster_bundles_map);
        remove_bundle_selection(to_be_removed, pre_bundles);
        to_be_removed.clear();
    }

    // ---- Second matching round ----
    {
        const unsigned int nbundle  = pre_bundles.size();
        const unsigned int nflash   = flash_bundles_map.size();
        const unsigned int ncluster = cluster_bundles_map.size();

        // Rebuild ordered iteration (rd-1 may have removed bundles/flashes/clusters).
        auto flashes_ordered  = flash_iter_order(flash_bundles_map);
        auto clusters_ordered = cluster_iter_order(cluster_bundles_map);
        std::map<Cluster*, int> cluster_idx_map;
        std::map<Opflash*, int> flash_idx_map;
        int idx = 0;
        for (auto* c : clusters_ordered) cluster_idx_map[c] = idx++;
        idx = 0;
        for (auto* f : flashes_ordered)  flash_idx_map[f]   = idx++;

        Ress::vector_t M = Ress::vector_t::Zero(nopdet * nflash);
        Ress::matrix_t P = Ress::matrix_t::Zero(nopdet * nflash, nbundle);
        Ress::vector_t MF = Ress::vector_t::Zero(ncluster);
        Ress::matrix_t PF = Ress::matrix_t::Zero(ncluster, nbundle);
        Ress::vector_t weights = Ress::vector_t::Zero(nbundle);
        std::vector<std::pair<Opflash*, Cluster*>> pairs;

        std::size_t i = 0, ik = 0;
        for (auto* flash : flashes_ordered) {
            auto& bundles = flash_bundles_map[flash];
            for (unsigned int j = 0; j < nopdet; ++j) {
                const int opdet_idx = opdet_idx_v.at(j);
                const double pe = flash->get_PE(opdet_idx);
                const double pe_err = std::sqrt(flash->get_PE(opdet_idx) + std::pow(flash->get_PE_err(opdet_idx), 2));
                M(i * nopdet + j) = pe / pe_err;
            }
            for (auto bundle : bundles) {
                const auto ks_dis = bundle->get_ks_dis();
                const auto& pred_flash = bundle->get_pred_flash();
                for (unsigned int j = 0; j < nopdet; ++j) {
                    const int opdet_idx = opdet_idx_v.at(j);
                    const double pred_pe = pred_flash.at(opdet_idx);
                    const double pe_err = std::sqrt(flash->get_PE(opdet_idx) + std::pow(flash->get_PE_err(opdet_idx), 2));
                    P(i * nopdet + j, pairs.size()) = pred_pe / pe_err;
                }
                pairs.emplace_back(flash, bundle->get_main_cluster());

                const auto meas_pe_tot = flash->get_total_PE();
                const auto pred_pe_tot = bundle->get_total_pred_light();
                const double base = (std::abs(pred_pe_tot - meas_pe_tot) > m_pe_mismatch_knee * meas_pe_tot)
                                        ? std::abs(pred_pe_tot - meas_pe_tot) / meas_pe_tot
                                        : m_pe_mismatch_floor;
                weights(ik++) = base + delta_shape * nopdet * ks_dis / lambda;
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

        int n = 0;
        for (auto* flash : flashes_ordered) {
            for (auto bundle : flash_bundles_map[flash]) {
                bundle->set_strength(solution(n));
                if (!(solution(n) > m_strength_cutoff || m_beamonly)) to_be_removed.push_back(bundle);
                ++n;
            }
        }
        remove_bundle_selection(to_be_removed, flash_bundles_map, cluster_bundles_map,
                                flash_cluster_bundles_map);
        remove_bundle_selection(to_be_removed, pre_bundles);
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
        for (auto cluster : clusters) {
            auto it = cluster_idx_map.find(cluster);
            if (it == cluster_idx_map.end()) continue;
            const int cidx = it->second;
            auto mit = matched_pairs.find(cidx);
            if (mit != matched_pairs.end()) {
                auto flash = mit->second.first;
                results_bundles.push_back(flash_cluster_bundles_map[std::make_pair(flash, cluster)]);
            }
            else {
                auto bundle = std::make_shared<TimingTPCBundle>(nullptr, cluster, 0, cidx);
                bundle->set_quality_params(qp);
                bundle->set_strength(0);
                results_bundles.push_back(bundle);
            }
        }
        organize_bundles(results_bundles, flash_cluster_bundles_map);

        FlashBundlesMap results_flash_bundles_map;
        for (auto bundle : results_bundles) {
            results_flash_bundles_map[bundle->get_flash()].push_back(bundle);
        }

        log->debug("done with matching");
    }

    // Apply matched t0s.
    for (auto* flash : flash_iter_order(flash_bundles_map)) {
        for (auto bundle : flash_bundles_map[flash]) {
            auto* cluster = bundle->get_main_cluster();
            cluster->set_cluster_t0(flash->get_time() * units::ns);
            // Record the matched flash row index so Cluster::get_flash() reflects
            // the match (flash_id == canonical "flash" PC row index).
            cluster->set_scalar<int>("flash", flash->get_flash_id());
            // Persist the match for the downstream (all-APA) MABC op/flash Bee
            // dump: a globally-unique flash id (survives the per-APA -> all-APA
            // merge) and this bundle's predicted per-channel PE (op_pes_pred).
            // The gid is keyed on the flash's index within `flashes` (unique
            // per APA), NOT get_flash_id() (the stored ident, which can repeat).
            cluster->set_scalar<int>("matched_flash_gid",
                                     m_anode->ident() * kFlashGidStride + global_flash_idx_map.at(flash));
            cluster->put_pcarray<double>(bundle->get_pred_flash(), "pe", "flashpred");
            log->debug("flash_bundles_map: flash id {} time {} ns, cluster gidx {} "
                       "total_pred_light {} t0 {}",
                       flash->get_flash_id(), flash->get_time(),
                       global_cluster_idx_map[bundle->get_main_cluster()],
                       bundle->get_total_pred_light(),
                       bundle->get_main_cluster()->get_cluster_t0());
        }
    }

    // Persist the per-flash measured light into a merge-safe, self-contained
    // per-root "opflash" PC, so the all-APA MABC can dump the op/flash Bee
    // display after the per-APA trees are merged (the canonical flash/light
    // /flashlight join uses per-APA row indices that collide across APAs; this
    // table is keyed by the global flash id instead). One row per (flash,
    // channel): gid, time (raw ns), ch, pe. Holds ALL flashes considered for
    // matching (incl. unmatched), as the op display shows every flash.
    {
        std::vector<int> op_gid, op_ch;
        std::vector<double> op_time, op_pe;
        for (std::size_t fi = 0; fi < flashes.size(); ++fi) {
            const auto& flash = flashes[fi];
            // gid keyed on the per-APA flash index (unique), matching the
            // matched_flash_gid stamped on clusters above.
            const int gid = m_anode->ident() * kFlashGidStride + static_cast<int>(fi);
            for (int ch = 0; ch < m_nchan; ++ch) {
                op_gid.push_back(gid);
                op_time.push_back(flash->get_time());
                op_ch.push_back(ch);
                op_pe.push_back(flash->get_PE(ch));
            }
        }
        grouping->put_pcarray<int>(op_gid, "gid", "opflash");
        grouping->put_pcarray<double>(op_time, "time", "opflash");
        grouping->put_pcarray<int>(op_ch, "ch", "opflash");
        grouping->put_pcarray<double>(op_pe, "pe", "opflash");
    }

    // ---- Build outputs ----
    {
        ITensor::vector outtens;
        auto tens_live = Aux::TensorDM::as_tensors(*root_live, inpath + "/live");
        outtens.insert(outtens.end(), tens_live.begin(), tens_live.end());

        auto root_dead = Aux::TensorDM::as_pctree(charge_tens, inpath + "/dead");
        auto tens_dead = Aux::TensorDM::as_tensors(*root_dead, inpath + "/dead");
        outtens.insert(outtens.end(), tens_dead.begin(), tens_dead.end());

        out = Aux::TensorDM::as_tensorset(outtens, charge_ident);
    }

    if (!flashes.empty()) {
        log->debug("total_charge_blob {} total_charge_point {} total_charge_blob_all {}",
                   total_charge_blob / flashes.size(),
                   total_charge_point / flashes.size(),
                   total_charge_blob_all);
    }
    else {
        log->debug("total_charge_blob {} total_charge_point {} total_charge_blob_all {}",
                   0, 0, total_charge_blob_all);
    }

    ++m_count;
    return true;
}

// ----- boundary-flag filling (port of ToyMatching.cxx::calculate_pred_pe ~176-290) -----
void QLMatching::compute_endpoint_flags(TimingTPCBundle* bundle,
                                        Cluster* cluster,
                                        double flash_x_offset,
                                        double s, double anode_x, double u_cathode) const
{
    // Collect the per-time-slice representative drift coordinate u and blob
    // count. The toolkit's time_blob_map is nested apa->face->time->BlobSet; the
    // matched anode has a single active face for SBND but we iterate all faces
    // under the apa for robustness. We then walk in u-order from each drift end
    // (anode = min u, cathode = max u). Walking by u rather than the prototype's
    // strict time order is a faithful adaptation: for a normally drifting track
    // time and u are monotonic, and the boundary tests are themselves in u.
    const auto& tbm = cluster->time_blob_map();
    auto ait = tbm.find(static_cast<int>(m_anode->ident()));
    if (ait == tbm.end()) return;

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
    if (sv.empty()) return;

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
    if (first_u > anode_in - 1.0 * units::cm &&
        last_u  > 0.0 &&
        last_u  < cathode_in &&
        first_u < u_cathode) {
        bundle->set_spec_end_flag(flag_spec_end);
        // Anode end inside the flag window => close to the PMTs (which sit at the
        // anode plane) and at the x-boundary.
        if (first_u <= m_anode_ext2 && first_u > anode_in - 1.0 * units::cm) {
            bundle->set_flag_close_to_PMT(true);
            bundle->set_flag_at_x_boundary(true);
        }
        // Cathode end inside the flag window => at the x-boundary (no PMTs there).
        if (last_u >= u_cathode + m_cathode_ext2 && last_u < cathode_in) {
            bundle->set_flag_at_x_boundary(true);
        }
    }
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

    for (auto& kv : eval_flash_bundles_map) {
        auto flash = kv.first;
        auto& orig_bundles = kv.second;
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
