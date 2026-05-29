#include "WireCellMatch/QLMatching.h"
#include "WireCellMatch/Util.h"
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

WIRECELL_FACTORY(QLMatching,
                 WireCell::Match::QLMatching,
                 WireCell::INamed,
                 WireCell::ITensorSetFilter,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Match;
using namespace WireCell::Clus::Facade;

QLMatching::QLMatching() : Aux::Logger("QLMatching", "match") {}
QLMatching::~QLMatching() = default;

void QLMatching::configure(const WireCell::Configuration& cfg)
{
    m_anode = Factory::find_tn<IAnodePlane>(cfg["anode"].asString());
    m_dv    = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());

    m_inpath        = get(cfg, "inpath", m_inpath);
    m_outpath       = get(cfg, "outpath", m_outpath);
    m_bee_dir       = get(cfg, "bee_dir", m_bee_dir);
    m_cluster_t0    = get(cfg, "cluster_t0", m_cluster_t0);
    m_semimodel_file = get(cfg, "semimodel_file", m_semimodel_file);

    m_pmts     = get(cfg, "pmts", m_pmts);
    m_data     = get(cfg, "data", m_data);
    m_beamonly = get(cfg, "beamonly", m_beamonly);

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

    std::vector<SemiAnalyticalModel::OpticalDetector> opdets;
    opdets.reserve(opdets_cfg.size());
    for (const auto& od : opdets_cfg) {
        SemiAnalyticalModel::OpticalDetector o;
        o.h           = get<double>(od, "h", -1.0);
        o.w           = get<double>(od, "w", -1.0);
        o.center      = WireCell::Point(od["x"].asDouble(), od["y"].asDouble(), od["z"].asDouble());
        o.type        = od.get("type", 1).asInt();
        o.orientation = od.get("orientation", 0).asInt();
        opdets.push_back(o);
    }

    m_semi_model = std::make_unique<SemiAnalyticalModel>(
        vuv_cfg, vis_cfg, geom, opdets, /*doReflectedLight=*/true);

    log->debug("QLMatching configured: nopdets={}, semimodel_file={}",
               opdets.size(), m_semimodel_file);
}

WireCell::Configuration QLMatching::default_configuration() const
{
    Configuration cfg;
    cfg["inpath"]          = m_inpath;
    cfg["outpath"]         = m_outpath;
    cfg["bee_dir"]         = m_bee_dir;
    cfg["nchan"]           = m_nchan;
    cfg["semimodel_file"]  = m_semimodel_file;
    cfg["pmts"]            = m_pmts;
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

    // SBND OpDet on/off mask. Default pattern (PMTs only on PMT slots) is
    // hard-coded for SBND's 312-OpDet layout; m_ch_mask further disables
    // specific channels.
    std::vector<unsigned int> opdet_mask;
    if (m_pmts) {
        opdet_mask = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
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

    // Flashes come from the canonical optical point clouds on the live root
    // node (written by Aux::OpflashToFlashPCs, the same schema as
    // root/UbooneClusterSource): "flash"(time,value,ident,type,...),
    // "light"(ident,time,value,error), "flashlight"(flash,light) join.
    // Rebuild one Opflash per flash row, zero-filling the per-channel PE vector
    // from the light entries via the join.
    std::vector<Opflash::pointer> flashes;
    {
        const auto& lpcs = root_live->value.local_pcs();
        auto fit = lpcs.find("flash");
        if (fit == lpcs.end()) {
            log->warn("no \"flash\" point cloud on live root; 0 flashes");
        }
        else {
            auto ftime_arr = fit->second.get("time");
            const size_t nflash = ftime_arr ? ftime_arr->elements<double>().size() : 0;
            auto ftime = ftime_arr ? ftime_arr->elements<double>()
                                   : PointCloud::Array::span_t<double>{};

            // light + join (absent => flashes carry no per-channel PE).
            PointCloud::Array::span_t<int>    fl_flash{}, fl_light{}, l_ident{};
            PointCloud::Array::span_t<double> l_value{};
            auto jit = lpcs.find("flashlight");
            auto lit = lpcs.find("light");
            if (jit != lpcs.end() && lit != lpcs.end()) {
                if (auto a = jit->second.get("flash"))  fl_flash = a->elements<int>();
                if (auto a = jit->second.get("light"))  fl_light = a->elements<int>();
                if (auto a = lit->second.get("ident"))  l_ident = a->elements<int>();
                if (auto a = lit->second.get("value"))  l_value = a->elements<double>();
            }

            for (size_t f = 0; f < nflash; ++f) {
                std::vector<double> pe(m_nchan, 0.0);
                for (size_t ifl = 0; ifl < fl_flash.size(); ++ifl) {
                    if (fl_flash[ifl] != (int) f) continue;
                    const int li = fl_light[ifl];
                    pe[l_ident[li]] = l_value[li];
                }
                Opflash::pointer flash =
                    std::make_shared<Opflash>(ftime[f], std::move(pe), 0.0, m_nchan);
                flash->set_flash_id((int) f);
                if (flash->get_time() < m_flash_mintime || flash->get_time() > m_flash_maxtime) continue;
                if (flash->get_total_PE() < m_flash_minPE) continue;
                flashes.push_back(flash);
            }
        }
    }

    auto grouping = root_live->value.facade<Grouping>();
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
                  [](Cluster* c) { c->set_cluster_t0(-1e12); c->set_scalar<int>("flash", -1); });

    std::map<Opflash*, int>  global_flash_idx_map;
    std::map<Cluster*, int>  global_cluster_idx_map;
    for (std::size_t i = 0; i < flashes.size(); ++i) global_flash_idx_map[flashes[i].get()] = i;
    for (std::size_t i = 0; i < clusters.size(); ++i) global_cluster_idx_map[clusters[i]] = i;

    std::vector<TimingTPCBundle::pointer> all_bundles;
    TimingTPCBundleSet pre_bundles;

    const unsigned int tpc = m_anode->ident();
    const int sign_offset  = (tpc == 0) ? -1 : 1;
    const double lo_x_bound = (tpc == 0) ? -2000 : 0;
    const double hi_x_bound = (tpc == 0) ? 0 : 2000;

    // Reduce mask to OpDets on this TPC.
    for (std::size_t idet = 0; idet < opdet_mask.size(); ++idet) {
        if ((tpc == 0) && (idet % 2 == 1)) opdet_mask[idet] = 0;
        if ((tpc == 1) && (idet % 2 == 0)) opdet_mask[idet] = 0;
    }

    for (auto flash : flashes) {
        const auto flash_time = flash->get_time();
        const double flash_x_offset = sign_offset * flash_time * m_drift_speed;

        // per-flash mask (also catches simulated saturated PMTs in MC).
        std::vector<unsigned int> flash_opdet_mask = opdet_mask;
        for (std::size_t idet = 0; idet < std::size_t(flash->get_num_channels()); ++idet) {
            auto pe_det = flash->get_PE(idet);
            if ((flash->get_total_PE() > 5000) && (pe_det == 0) && (m_data == false))
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
            all_bundles.push_back(bundle);
            bundle->set_opdet_mask(flash_opdet_mask);

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

                    if (x < lo_x_bound || x > hi_x_bound) { ++npt_outside_drift; continue; }
                    if (std::abs(y) > 2000 || z < 0 || z > 5000) { ++npt_outside_bounds; continue; }

                    if (std::abs(x) && bundle->get_flag_at_x_boundary() == false)
                        bundle->set_flag_close_to_PMT(true);
                    if (std::abs(x) > 1950 && bundle->get_flag_close_to_PMT() == false)
                        bundle->set_flag_close_to_PMT(true);
                    if (npt_outside_drift > 0.25 * npt) { drifted_outside = true; break; }

                    // SemiAnalyticalModel expects positions in cm. Blob points
                    // are in WCT units (mm); units::cm == 10.
                    const WireCell::Point xyz_cm(x / 10., y / 10., z / 10.);
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
            if (bundle->get_total_pred_light() < 10) continue;
            bundle->examine_bundle();
            if (bundle->get_ks_dis() == 1) {
                bundle->set_potential_bad_match_flag(true);
                continue;
            }
            if (bundle->get_chi2() / bundle->get_ndf() > 1e4) {
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

    const double lambda       = 0.1;
    const double delta_charge = 0.01;
    const double delta_light  = 0.025;
    const double delta_shape  = 0.01;

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
                weights(ik++) = (std::abs(pred_pe_tot - meas_pe_tot) > 0.3 * meas_pe_tot)
                                  ? std::abs(pred_pe_tot - meas_pe_tot) / meas_pe_tot
                                  : 0.3;
            }
            PF(ncluster + i, nbundle + i) = 1. / delta_light;
            flash_idx_map[flash] = nbundle + i;
            ++i;
        }
        for (unsigned int k = 0; k < nflash; ++k) weights(nbundle + k) = 0.5;
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
                const double base = (std::abs(pred_pe_tot - meas_pe_tot) > 0.3 * meas_pe_tot)
                                        ? std::abs(pred_pe_tot - meas_pe_tot) / meas_pe_tot
                                        : 0.3;
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
        if (!m_bee_dir.empty()) {
            const std::string sub_dir = String::format("%s/%d", m_bee_dir, m_bee_index);
            Persist::assuredir(sub_dir);
            Match::dump_bee_3d(
                *root_live.get(),
                String::format("%s/%d-img-apa%d.json", sub_dir, m_bee_index, m_anode->ident()));
            Match::dump_light(
                flashes, flash_bundles_map, global_cluster_idx_map,
                String::format("%s/%d-op-apa%d.json", sub_dir, m_bee_index, m_anode->ident()));
            ++m_bee_index;
        }
        log->debug(em("dump bee"));
    }

    // Apply matched t0s.
    for (auto* flash : flash_iter_order(flash_bundles_map)) {
        for (auto bundle : flash_bundles_map[flash]) {
            auto* cluster = bundle->get_main_cluster();
            cluster->set_cluster_t0(flash->get_time() * units::ns);
            // Record the matched flash row index so Cluster::get_flash() reflects
            // the match (flash_id == canonical "flash" PC row index).
            cluster->set_scalar<int>("flash", flash->get_flash_id());
            log->debug("flash_bundles_map: flash id {} time {} ns, cluster gidx {} "
                       "total_pred_light {} t0 {}",
                       flash->get_flash_id(), flash->get_time(),
                       global_cluster_idx_map[bundle->get_main_cluster()],
                       bundle->get_total_pred_light(),
                       bundle->get_main_cluster()->get_cluster_t0());
        }
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
            if (bundle->get_ks_dis() > 0.2 || bundle->get_chi2() / bundle->get_ndf() > 20) {
                to_be_removed.push_back(bundle);
                continue;
            }
            if (std::abs(flash->get_total_PE() - bundle->get_total_pred_light()) >
                0.5 * flash->get_total_PE()) {
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
