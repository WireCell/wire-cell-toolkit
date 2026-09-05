#include "WireCellRoot/SbndMagnifyTrackingVisitor.h"

#include <cmath>

#include "TFile.h"
#include "TTree.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Units.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"

WIRECELL_FACTORY(SbndMagnifyTrackingVisitor, WireCell::Root::SbndMagnifyTrackingVisitor,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;

Root::SbndMagnifyTrackingVisitor::SbndMagnifyTrackingVisitor()
  : log(Log::logger("tracking"))
{
}

Root::SbndMagnifyTrackingVisitor::~SbndMagnifyTrackingVisitor() {}

void Root::SbndMagnifyTrackingVisitor::configure(const WireCell::Configuration& cfg)
{
    m_output_filename = get<std::string>(cfg, "output_filename", "tracking-stm.root");
    m_grouping_name = get<std::string>(cfg, "grouping", "live");
    m_track_fitting_name = get<std::string>(cfg, "track_fitting_name", "stm");
    m_runNo = get<int>(cfg, "runNo", 0);
    m_subRunNo = get<int>(cfg, "subRunNo", 0);
    m_eventNo = get<int>(cfg, "eventNo", 0);
    m_dQdx_scale = get<double>(cfg, "dQdx_scale", 0.1);
    m_dQdx_offset = get<double>(cfg, "dQdx_offset", -1000);
    m_nticks = get<int>(cfg, "nticks", m_nticks);

    auto anode_tns = cfg["anodes"];
    for (auto anode_tn : anode_tns) {
        auto anode = Factory::find_tn<IAnodePlane>(anode_tn.asString());
        m_anodes.push_back(anode);
    }

    m_dv = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());
}

WireCell::Configuration Root::SbndMagnifyTrackingVisitor::default_configuration() const
{
    Configuration cfg;
    cfg["output_filename"] = "tracking-stm.root";
    cfg["grouping"] = "live";
    cfg["track_fitting_name"] = "stm";
    cfg["anodes"] = Json::arrayValue;
    cfg["detector_volumes"] = "";
    cfg["runNo"] = 0;
    cfg["subRunNo"] = 0;
    cfg["eventNo"] = 0;
    cfg["dQdx_scale"] = 0.1;
    cfg["dQdx_offset"] = -1000;
    cfg["nticks"] = 3427;
    return cfg;
}

Root::SbndMagnifyTrackingVisitor::ChanScheme Root::SbndMagnifyTrackingVisitor::chan_scheme() const
{
    ChanScheme cs;
    for (const auto& anode : m_anodes) {
        for (const auto& face : anode->faces()) {
            if (!face) continue;
            auto planes = face->planes();
            for (size_t p = 0; p < planes.size() && p < 3; ++p) {
                const int n = static_cast<int>(planes[p]->channels().size());
                if (n > cs.nch[p]) cs.nch[p] = n;
            }
        }
    }
    const int napa = static_cast<int>(m_anodes.size());
    cs.base[0] = 0;
    cs.base[1] = napa * cs.nch[0];
    cs.base[2] = cs.base[1] + napa * cs.nch[1];
    return cs;
}


namespace {
    /// One file per event, and this event's RSE, when a group runs in one
    /// process.  A name with no '%' and an ensemble with no published RSE give
    /// back exactly the configured name and numbers, so the per-event job is
    /// unchanged.
    std::string event_filename(const std::string& tmpl, int ident)
    {
        if (tmpl.find('%') == std::string::npos) return tmpl;
        return WireCell::String::format(tmpl, ident);
    }
    struct Rse { int run, subrun, event; };
    Rse event_rse(const WireCell::Clus::Facade::Ensemble& ensemble,
                  int run, int subrun, int event)
    {
        if (ensemble.rse_valid()) {
            return {ensemble.runNo(), ensemble.subRunNo(), ensemble.eventNo()};
        }
        return {run, subrun, event};
    }
}

void Root::SbndMagnifyTrackingVisitor::visit(Clus::Facade::Ensemble& ensemble) const
{
    auto groupings = ensemble.with_name(m_grouping_name);
    if (groupings.empty()) {
        log->debug("SbndMagnifyTrackingVisitor: no grouping '{}'", m_grouping_name);
        return;
    }

    auto& grouping = *groupings.at(0);

    grouping.set_anodes(m_anodes);
    grouping.set_detector_volumes(m_dv);

    const auto cs = chan_scheme();
    log->debug("SbndMagnifyTrackingVisitor: channel scheme nch=[{},{},{}] base=[{},{},{}]",
               cs.nch[0], cs.nch[1], cs.nch[2], cs.base[0], cs.base[1], cs.base[2]);

    {
        const auto rse = event_rse(ensemble, m_runNo, m_subRunNo, m_eventNo);
        m_evt_runNo = rse.run; m_evt_subRunNo = rse.subrun; m_evt_eventNo = rse.event;
    }
    const std::string outname = event_filename(m_output_filename, ensemble.ident());
    TFile* output_tf = TFile::Open(outname.c_str(), "RECREATE");
    if (!output_tf || output_tf->IsZombie()) {
        log->error("SbndMagnifyTrackingVisitor: cannot open {}", outname);
        return;
    }

    write_bad_channels(output_tf, grouping, cs);
    write_trun(output_tf);
    write_proj_data(output_tf, grouping, cs);
    write_t_rec_data(output_tf, grouping, cs);
    write_stm_trees(output_tf, grouping);

    // Empty T_proj tree kept for reader compatibility.
    TTree* tree_proj = new TTree("T_proj", "T_proj");
    tree_proj->SetDirectory(output_tf);
    tree_proj->Fill();

    output_tf->Write();
    output_tf->Close();
    delete output_tf;

    log->debug("SbndMagnifyTrackingVisitor: wrote {}", outname);
}

void Root::SbndMagnifyTrackingVisitor::write_bad_channels(TFile* output_tf, Clus::Facade::Grouping& grouping,
                                                          const ChanScheme& cs) const
{
    TTree* tree = new TTree("T_bad_ch", "T_bad_ch");
    tree->SetDirectory(output_tf);

    int chid = 0;
    int plane = 0;
    int start_time = 0;
    int end_time = 0;
    int runNo = m_evt_runNo;
    int subRunNo = m_evt_subRunNo;
    int eventNo = m_evt_eventNo;

    tree->Branch("chid", &chid, "chid/I");
    tree->Branch("plane", &plane, "plane/I");
    tree->Branch("start_time", &start_time, "start_time/I");
    tree->Branch("end_time", &end_time, "end_time/I");
    tree->Branch("runNo", &runNo, "runNo/I");
    tree->Branch("subRunNo", &subRunNo, "subRunNo/I");
    tree->Branch("eventNo", &eventNo, "eventNo/I");

    auto wpids = grouping.wpids();
    std::set<std::pair<int, int>> apa_face_set;
    for (const auto& wpid : wpids) {
        apa_face_set.insert({wpid.apa(), wpid.face()});
    }

    for (const auto& [apa, face] : apa_face_set) {
        for (int pind = 0; pind < 3; ++pind) {
            try {
                auto dead_chs = grouping.get_all_dead_chs(apa, face, pind);
                plane = pind;
                for (const auto& [ch, time_range] : dead_chs) {
                    chid = cs.global(pind, apa, ch);
                    // get_all_dead_chs converts an x range to ticks, so a dead
                    // region with an unbounded x extent comes back as +-1e9.
                    // The channel is dead for the whole readout either way, and
                    // an unclamped range makes the Magnify projection pads
                    // unreadable.
                    start_time = std::max(0, std::min(m_nticks, time_range.first));
                    end_time = std::max(0, std::min(m_nticks, time_range.second));
                    if (end_time <= start_time) { start_time = 0; end_time = m_nticks; }
                    tree->Fill();
                }
            }
            catch (const std::exception& e) {
                log->warn("SbndMagnifyTrackingVisitor: failed to get dead channels for APA={}, face={}, plane={}: {}",
                          apa, face, pind, e.what());
            }
        }
    }

    log->debug("SbndMagnifyTrackingVisitor: wrote {} entries to T_bad_ch", tree->GetEntries());
}

// Per-pass and per-eval STM decision records (from TaggerCheckSTM's stm_pass /
// stm_eval cluster PCs) so downstream python (nusel_display dQ/dx panel) has
// the tagger's numbers without log parsing.
void Root::SbndMagnifyTrackingVisitor::write_stm_trees(TFile* output_tf, Clus::Facade::Grouping& grouping) const
{
    int cluster_id{0}, pass{0}, status{0}, kink_num{0}, npoints{0};
    double exit_L{0}, left_L{0}, exit_dqdx{0}, left_dqdx{0};
    TTree* t_pass = new TTree("T_stm_pass", "T_stm_pass");
    t_pass->SetDirectory(output_tf);
    t_pass->Branch("cluster_id", &cluster_id, "cluster_id/I");
    t_pass->Branch("pass", &pass, "pass/I");
    t_pass->Branch("status", &status, "status/I");
    t_pass->Branch("kink_num", &kink_num, "kink_num/I");
    t_pass->Branch("npoints", &npoints, "npoints/I");
    t_pass->Branch("exit_L", &exit_L, "exit_L/D");
    t_pass->Branch("left_L", &left_L, "left_L/D");
    t_pass->Branch("exit_dqdx", &exit_dqdx, "exit_dqdx/D");
    t_pass->Branch("left_dqdx", &left_dqdx, "left_dqdx/D");

    int e_pass{0}, e_strong{0}, e_verdict{0};
    double e_peak{0}, e_off{0}, e_com{0}, e_ks1{0}, e_ks2{0}, e_r1{0}, e_r2{0}, e_resl{0}, e_resq{0};
    TTree* t_eval = new TTree("T_stm_eval", "T_stm_eval");
    t_eval->SetDirectory(output_tf);
    t_eval->Branch("cluster_id", &cluster_id, "cluster_id/I");
    t_eval->Branch("pass", &e_pass, "pass/I");
    t_eval->Branch("strong", &e_strong, "strong/I");
    t_eval->Branch("verdict", &e_verdict, "verdict/I");
    t_eval->Branch("peak_range", &e_peak, "peak_range/D");
    t_eval->Branch("offset_length", &e_off, "offset_length/D");
    t_eval->Branch("com_range", &e_com, "com_range/D");
    t_eval->Branch("ks1", &e_ks1, "ks1/D");
    t_eval->Branch("ks2", &e_ks2, "ks2/D");
    t_eval->Branch("ratio1", &e_r1, "ratio1/D");
    t_eval->Branch("ratio2", &e_r2, "ratio2/D");
    t_eval->Branch("res_length", &e_resl, "res_length/D");
    t_eval->Branch("ave_res_dqdx", &e_resq, "ave_res_dqdx/D");

    for (const auto* cluster : grouping.children()) {
        cluster_id = cluster->get_cluster_id();

        const auto& pass_pc = cluster->get_pc("stm_pass");
        if (!pass_pc.empty()) {
            const auto& v_pass = pass_pc.get("pass")->elements<int>();
            const auto& v_status = pass_pc.get("status")->elements<int>();
            const auto& v_kink = pass_pc.get("kink_num")->elements<int>();
            const auto& v_npts = pass_pc.get("npoints")->elements<int>();
            const auto& v_exit_L = pass_pc.get("exit_L")->elements<double>();
            const auto& v_left_L = pass_pc.get("left_L")->elements<double>();
            const auto& v_exit_dqdx = pass_pc.get("exit_dqdx")->elements<double>();
            const auto& v_left_dqdx = pass_pc.get("left_dqdx")->elements<double>();
            for (size_t i = 0; i < v_pass.size(); ++i) {
                pass = v_pass[i]; status = v_status[i]; kink_num = v_kink[i]; npoints = v_npts[i];
                exit_L = v_exit_L[i] / units::cm;
                left_L = v_left_L[i] / units::cm;
                exit_dqdx = v_exit_dqdx[i];
                left_dqdx = v_left_dqdx[i];
                t_pass->Fill();
            }
        }

        const auto& eval_pc = cluster->get_pc("stm_eval");
        if (!eval_pc.empty()) {
            const auto& v_pass = eval_pc.get("pass")->elements<int>();
            const auto& v_strong = eval_pc.get("strong")->elements<int>();
            const auto& v_verdict = eval_pc.get("verdict")->elements<int>();
            const auto& v_peak = eval_pc.get("peak_range")->elements<double>();
            const auto& v_off = eval_pc.get("offset_length")->elements<double>();
            const auto& v_com = eval_pc.get("com_range")->elements<double>();
            const auto& v_ks1 = eval_pc.get("ks1")->elements<double>();
            const auto& v_ks2 = eval_pc.get("ks2")->elements<double>();
            const auto& v_r1 = eval_pc.get("ratio1")->elements<double>();
            const auto& v_r2 = eval_pc.get("ratio2")->elements<double>();
            const auto& v_resl = eval_pc.get("res_length")->elements<double>();
            const auto& v_resq = eval_pc.get("ave_res_dqdx")->elements<double>();
            for (size_t i = 0; i < v_pass.size(); ++i) {
                e_pass = v_pass[i]; e_strong = v_strong[i]; e_verdict = v_verdict[i];
                e_peak = v_peak[i] / units::cm;
                e_off = v_off[i] / units::cm;
                e_com = v_com[i] / units::cm;
                e_ks1 = v_ks1[i]; e_ks2 = v_ks2[i]; e_r1 = v_r1[i]; e_r2 = v_r2[i];
                e_resl = v_resl[i] / units::cm;
                e_resq = v_resq[i];
                t_eval->Fill();
            }
        }
    }

    log->debug("SbndMagnifyTrackingVisitor: wrote T_stm_pass {} / T_stm_eval {} entries",
               t_pass->GetEntries(), t_eval->GetEntries());
}

void Root::SbndMagnifyTrackingVisitor::write_trun(TFile* output_tf) const
{
    TTree* tree = new TTree("Trun", "Trun");
    tree->SetDirectory(output_tf);

    int eventNo = m_evt_eventNo;
    int runNo = m_evt_runNo;
    int subRunNo = m_evt_subRunNo;
    double dQdx_scale = m_dQdx_scale;
    double dQdx_offset = m_dQdx_offset;

    tree->Branch("eventNo", &eventNo, "eventNo/I");
    tree->Branch("runNo", &runNo, "runNo/I");
    tree->Branch("subRunNo", &subRunNo, "subRunNo/I");
    tree->Branch("dQdx_scale", &dQdx_scale, "dQdx_scale/D");
    tree->Branch("dQdx_offset", &dQdx_offset, "dQdx_offset/D");

    tree->Fill();
}

void Root::SbndMagnifyTrackingVisitor::write_proj_data(TFile* output_tf, Clus::Facade::Grouping& grouping,
                                                       const ChanScheme& cs) const
{
    auto tf = grouping.get_track_fitting(m_track_fitting_name);
    if (!tf) {
        log->warn("SbndMagnifyTrackingVisitor: no TrackFitting in grouping slot '{}'", m_track_fitting_name);
        return;
    }

    const auto& fitted = tf->get_fitted_charge_2d();
    const auto& per_cluster = tf->get_cluster_fitted_charge_2d();
    if (fitted.empty() && per_cluster.empty()) {
        log->warn("SbndMagnifyTrackingVisitor: fitted_charge_2d is empty");
        return;
    }

    // Which passes exist per cluster (from the stm_pass PC) so a cluster's
    // 2D charge can be duplicated under each pass's T_rec block id
    // (cluster_id*10 + pass), keeping the GUI's rec<->proj lookup working.
    std::map<int, std::vector<int>> cluster_passes;
    for (const auto* cluster : grouping.children()) {
        const auto& pass_pc = cluster->get_pc("stm_pass");
        if (pass_pc.empty()) continue;
        const auto& passes = pass_pc.get("pass")->elements<int>();
        auto& dst = cluster_passes[cluster->get_cluster_id()];
        dst.assign(passes.begin(), passes.end());
    }

    auto nticks_map = grouping.get_nticks_per_slice();

    // doc pr/109 §8: ordered per-cell accumulation, keyed by the pass block id.
    // prepare_data's dead-region filler writes one snapshot entry per TICK
    // (TrackFitting.cxx:896), so nticks_per_slice of them divide down to the
    // same time_slice; without this they would be duplicate keys in one block.
    // Charge adds, errors add in quadrature, predictions add -- identities for
    // a live cell, which occurs once per slice.
    // Live entries and dead-region fillers are kept apart: a slice can hold a
    // real readout at one tick and fillers at the others (prepare_data inserts
    // a filler only where no key exists, TrackFitting.cxx:904), and adding the
    // filler's spread charge onto a measured cell would inflate the measured
    // charge the tree reports.  Fillers are used only when the slice has no
    // live entry at all -- which is the case that carries the dead-blob charge.
    struct Cell {
        double live_charge{0}, live_err2{0};
        double dead_charge{0}, dead_err2{0};
        double pred{0};
        int nlive{0};
        double charge() const { return nlive ? live_charge : dead_charge; }
        double charge_err() const { return std::sqrt(nlive ? live_err2 : dead_err2); }
    };
    std::map<int, std::map<std::pair<int, int>, Cell>> cluster_cells;

    // doc pr/109 §8: prefer each cluster's OWN snapshot over the merged map,
    // whose pred_charge is last-writer-wins across clusters.  The STM slot this
    // visitor reads is a holder fed by merge_fitted_charge_2d() AND, since doc
    // pdvd/42, by one add_fitted_charge_2d_snapshot() per fitted pass (the
    // merged map alone kept only the LAST cluster's prediction).  The merged
    // fallback below is reached only by a holder without snapshots.
    auto emit = [&](std::map<int, std::map<std::pair<int, int>, Cell>>& local,
                    int cid, int apa, int face, int plane_idx,
                    const Clus::TrackFitting::WireTime& wt,
                    const Clus::TrackFitting::FittedCharge2D& fc,
                    int forced_pass = -1) {
        int nticks_per_slice = 1;
        auto ait = nticks_map.find(apa);
        if (ait != nticks_map.end()) {
            auto fit_ = ait->second.find(face);
            if (fit_ != ait->second.end()) nticks_per_slice = fit_->second;
        }
        const int time = wt.second / nticks_per_slice;
        const int channel = cs.global(plane_idx, apa, wt.first);
        auto pit = cluster_passes.find(cid);
        // Fall back to a single pass-0 block when the cluster carries
        // no stm_pass PC.  This DOES happen: a fitted cell can belong
        // to an associated cluster that STM itself never fit, so
        // T_proj_data can hold a block with no T_rec track (the GUI
        // shows it only under "all clusters").
        static const std::vector<int> pass0{0};
        // doc pdvd/42: a snapshot that names its STM pass goes under that one
        // block only -- the forward and backward fits of one cluster are two
        // different fits with two different predictions.
        const std::vector<int> one{forced_pass};
        const auto& passes = (forced_pass >= 0) ? one
                           : (pit == cluster_passes.end()) ? pass0 : pit->second;
        for (int pass : passes) {
            auto& cell = local[cid * 10 + pass][{channel, time}];
            if (fc.flag != 0) {
                cell.live_charge += fc.charge;
                cell.live_err2 += fc.charge_err * fc.charge_err;
                ++cell.nlive;
            }
            else {
                cell.dead_charge += fc.charge;
                cell.dead_err2 += fc.charge_err * fc.charge_err;
            }
            cell.pred += fc.pred_charge;
        }
    };

    // ACROSS snapshots landing in the same block (two live clusters sharing an
    // ident), charge and charge_err are properties of the READOUT and are taken
    // once, not summed -- only the predictions add.
    auto flush = [&](const std::map<int, std::map<std::pair<int, int>, Cell>>& local) {
        for (const auto& [block, cells] : local) {
            auto& dst = cluster_cells[block];
            for (const auto& [key, cell] : cells) {
                auto it = dst.find(key);
                if (it == dst.end()) dst.emplace(key, cell);
                else it->second.pred += cell.pred;
            }
        }
    };

    if (!per_cluster.empty()) {
        for (const auto& snap : per_cluster) {
            std::map<int, std::map<std::pair<int, int>, Cell>> local;
            for (const auto& [afp, wt_map] : snap.cells)
                for (const auto& [wt, fc] : wt_map) {
                    // doc pdvd/42: an STM pass snapshot spans the bounding box
                    // of main + associated clusters (TrackFitting::prepare_data),
                    // so most of its cells belong to OTHER clusters and carry
                    // pred 0 by construction.  Keep the fitted cluster's own
                    // cells (owner set from global_rb_map) plus owner-less
                    // dead-region fillers; the block is then "this cluster's
                    // measured 2D charge against THIS pass's prediction", which
                    // is what the GUI showed before and what a residual metric
                    // needs.  PR-stage snapshots (pass < 0) keep every cell.
                    if (snap.pass >= 0 && snap.cluster && !fc.clusters.empty()
                        && !fc.clusters.count(snap.cluster)) continue;
                    emit(local, snap.ident, std::get<0>(afp), std::get<1>(afp), std::get<2>(afp), wt, fc, snap.pass);
                }
            flush(local);
        }
    }
    else {
        std::map<int, std::map<std::pair<int, int>, Cell>> local;
        for (const auto& [afp, wt_map] : fitted)
            for (const auto& [wt, fc] : wt_map)
                for (auto* cl : fc.clusters)
                    emit(local, cl->get_cluster_id(), std::get<0>(afp), std::get<1>(afp), std::get<2>(afp), wt, fc);
        flush(local);
    }

    std::vector<int> v_cluster_id;
    std::vector<std::vector<int>> v_channel;
    std::vector<std::vector<int>> v_time_slice;
    std::vector<std::vector<int>> v_charge;
    std::vector<std::vector<int>> v_charge_err;
    std::vector<std::vector<int>> v_charge_pred;

    for (const auto& [cid, cells] : cluster_cells) {
        v_cluster_id.push_back(cid);
        std::vector<int> ch, ts, q, qe, qp;
        ch.reserve(cells.size()); ts.reserve(cells.size()); q.reserve(cells.size());
        qe.reserve(cells.size()); qp.reserve(cells.size());
        for (const auto& [key, cell] : cells) {           // std::map -> ordered
            ch.push_back(key.first);
            ts.push_back(key.second);
            q.push_back(static_cast<int>(cell.charge()));
            qe.push_back(static_cast<int>(cell.charge_err()));
            qp.push_back(static_cast<int>(cell.pred));
        }
        v_channel.push_back(std::move(ch));
        v_time_slice.push_back(std::move(ts));
        v_charge.push_back(std::move(q));
        v_charge_err.push_back(std::move(qe));
        v_charge_pred.push_back(std::move(qp));
    }

    TTree* tree = new TTree("T_proj_data", "T_proj_data");
    tree->SetDirectory(output_tf);
    tree->Branch("cluster_id", &v_cluster_id);
    tree->Branch("channel", &v_channel);
    tree->Branch("time_slice", &v_time_slice);
    tree->Branch("charge", &v_charge);
    tree->Branch("charge_err", &v_charge_err);
    tree->Branch("charge_pred", &v_charge_pred);
    tree->Fill();

    log->debug("SbndMagnifyTrackingVisitor: wrote T_proj_data with {} blocks", v_cluster_id.size());
}

void Root::SbndMagnifyTrackingVisitor::write_t_rec_data(TFile* output_tf, Clus::Facade::Grouping& grouping,
                                                        const ChanScheme& cs) const
{
    auto nticks_map = grouping.get_nticks_per_slice();

    double b_x{0}, b_y{0}, b_z{0}, b_q{0}, b_nq{0}, b_chi2{1}, b_ndf{0};
    double b_pu{0}, b_pv{0}, b_pw{0}, b_pt{0}, b_reduced_chi2{0}, b_rr{0};
    int b_flag_vertex{0}, b_flag_shower{0};
    int b_cluster_id{0}, b_real_cluster_id{0}, b_sub_cluster_id{0};
    int b_pass{0}, b_status{0};

    TTree* t_rec_charge = new TTree("T_rec_charge", "T_rec_charge");
    t_rec_charge->SetDirectory(output_tf);
    t_rec_charge->Branch("x", &b_x, "x/D");
    t_rec_charge->Branch("y", &b_y, "y/D");
    t_rec_charge->Branch("z", &b_z, "z/D");
    t_rec_charge->Branch("q", &b_q, "q/D");
    t_rec_charge->Branch("nq", &b_nq, "nq/D");
    t_rec_charge->Branch("chi2", &b_chi2, "chi2/D");
    t_rec_charge->Branch("ndf", &b_ndf, "ndf/D");
    t_rec_charge->Branch("pu", &b_pu, "pu/D");
    t_rec_charge->Branch("pv", &b_pv, "pv/D");
    t_rec_charge->Branch("pw", &b_pw, "pw/D");
    t_rec_charge->Branch("pt", &b_pt, "pt/D");
    t_rec_charge->Branch("reduced_chi2", &b_reduced_chi2, "reduced_chi2/D");
    t_rec_charge->Branch("flag_vertex", &b_flag_vertex, "flag_vertex/I");
    t_rec_charge->Branch("flag_shower", &b_flag_shower, "flag_shower/I");
    t_rec_charge->Branch("rr", &b_rr, "rr/D");
    t_rec_charge->Branch("cluster_id", &b_cluster_id, "cluster_id/I");
    t_rec_charge->Branch("real_cluster_id", &b_real_cluster_id, "real_cluster_id/I");
    t_rec_charge->Branch("sub_cluster_id", &b_sub_cluster_id, "sub_cluster_id/I");
    // STM-stage extras (absent in the uBooNE PR-stage writer).
    t_rec_charge->Branch("pass", &b_pass, "pass/I");
    t_rec_charge->Branch("status", &b_status, "status/I");

    int npoints = 0;
    for (const auto* cluster : grouping.children()) {
        const auto& fit_pc = cluster->get_pc("stm_fit");
        if (fit_pc.empty()) continue;

        const auto& x = fit_pc.get("x")->elements<double>();
        const auto& y = fit_pc.get("y")->elements<double>();
        const auto& z = fit_pc.get("z")->elements<double>();
        const auto& dQ = fit_pc.get("dQ")->elements<double>();
        const auto& dx = fit_pc.get("dx")->elements<double>();
        const auto& rr = fit_pc.get("rr")->elements<double>();
        const auto& pu = fit_pc.get("pu")->elements<double>();
        const auto& pv = fit_pc.get("pv")->elements<double>();
        const auto& pw = fit_pc.get("pw")->elements<double>();
        const auto& pt = fit_pc.get("pt")->elements<double>();
        const auto& chi2 = fit_pc.get("reduced_chi2")->elements<double>();
        const auto& apa = fit_pc.get("apa")->elements<int>();
        const auto& face = fit_pc.get("face")->elements<int>();
        const auto& pass = fit_pc.get("pass")->elements<int>();
        const auto& status = fit_pc.get("status")->elements<int>();

        const int cid = cluster->get_cluster_id();

        for (size_t i = 0; i < x.size(); ++i) {
            int nticks_per_slice = 1;
            auto ait = nticks_map.find(apa[i]);
            if (ait != nticks_map.end()) {
                auto fit_ = ait->second.find(face[i]);
                if (fit_ != ait->second.end()) nticks_per_slice = fit_->second;
            }

            b_x = x[i] / units::cm;
            b_y = y[i] / units::cm;
            b_z = z[i] / units::cm;
            b_q = dQ[i] * m_dQdx_scale + m_dQdx_offset;
            b_nq = dx[i] / units::cm;
            b_ndf = cid * 10 + pass[i];
            // Keep the fractional wire position (uBooNE writer does the same).
            b_pu = cs.base[0] + apa[i] * cs.nch[0] + pu[i];
            b_pv = cs.base[1] + apa[i] * cs.nch[1] + pv[i];
            b_pw = cs.base[2] + apa[i] * cs.nch[2] + pw[i];
            b_pt = pt[i] / nticks_per_slice;
            b_reduced_chi2 = chi2[i];
            b_rr = rr[i] / units::cm;
            b_cluster_id = cid;
            b_real_cluster_id = cid * 10 + pass[i];
            b_sub_cluster_id = cid * 10 + pass[i];
            b_pass = pass[i];
            b_status = status[i];
            t_rec_charge->Fill();
            ++npoints;
        }
    }

    log->debug("SbndMagnifyTrackingVisitor: wrote {} entries to T_rec_charge", npoints);
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
