#include "WireCellRoot/UbooneMagnifyTrackingVisitor.h"

#include "TFile.h"
#include "TTree.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade_Grouping.h"

WIRECELL_FACTORY(UbooneMagnifyTrackingVisitor, WireCell::Root::UbooneMagnifyTrackingVisitor,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;

Root::UbooneMagnifyTrackingVisitor::UbooneMagnifyTrackingVisitor()
  : log(Log::logger("tracking"))
{
}

Root::UbooneMagnifyTrackingVisitor::~UbooneMagnifyTrackingVisitor() {}

void Root::UbooneMagnifyTrackingVisitor::configure(const WireCell::Configuration& cfg)
{
    m_output_filename = get<std::string>(cfg, "output_filename", "tracking_proj.root");
    m_grouping_name = get<std::string>(cfg, "grouping", "live");
    m_runNo = get<int>(cfg, "runNo", 0);
    m_subRunNo = get<int>(cfg, "subRunNo", 0);
    m_eventNo = get<int>(cfg, "eventNo", 0);

    auto anode_tns = cfg["anodes"];
    for (auto anode_tn : anode_tns) {
        auto anode = Factory::find_tn<IAnodePlane>(anode_tn.asString());
        m_anodes.push_back(anode);
    }

    m_dv = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());
}

WireCell::Configuration Root::UbooneMagnifyTrackingVisitor::default_configuration() const
{
    Configuration cfg;
    cfg["output_filename"] = "tracking_proj.root";
    cfg["grouping"] = "live";
    cfg["anodes"] = Json::arrayValue;
    cfg["detector_volumes"] = "";
    cfg["runNo"] = 0;
    cfg["subRunNo"] = 0;
    cfg["eventNo"] = 0;
    return cfg;
}

void Root::UbooneMagnifyTrackingVisitor::visit(Clus::Facade::Ensemble& ensemble) const
{
    auto groupings = ensemble.with_name(m_grouping_name);
    if (groupings.empty()) {
        log->debug("UbooneMagnifyTrackingVisitor: no grouping '{}'", m_grouping_name);
        return;
    }

    auto& grouping = *groupings.at(0);

    // Set anodes and detector volumes on the grouping
    grouping.set_anodes(m_anodes);
    grouping.set_detector_volumes(m_dv);

    // Open ROOT file
    TFile* output_tf = TFile::Open(m_output_filename.c_str(), "RECREATE");
    if (!output_tf || output_tf->IsZombie()) {
        log->error("UbooneMagnifyTrackingVisitor: cannot open {}", m_output_filename);
        return;
    }

    write_bad_channels(output_tf, grouping);
    write_proj_data(output_tf, grouping);

    // Empty T_proj tree for now
    TTree* tree_proj = new TTree("T_proj", "T_proj");
    tree_proj->SetDirectory(output_tf);
    tree_proj->Fill();

    output_tf->Write();
    output_tf->Close();
    delete output_tf;

    log->debug("UbooneMagnifyTrackingVisitor: wrote {}", m_output_filename);
}

void Root::UbooneMagnifyTrackingVisitor::write_bad_channels(TFile* output_tf, Clus::Facade::Grouping& grouping) const
{
    TTree* tree = new TTree("T_bad_ch", "T_bad_ch");
    tree->SetDirectory(output_tf);

    int chid = 0;
    int plane = 0;
    int start_time = 0;
    int end_time = 0;
    int runNo = m_runNo;
    int subRunNo = m_subRunNo;
    int eventNo = m_eventNo;

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
                    chid = ch;
                    start_time = time_range.first;
                    end_time = time_range.second;
                    tree->Fill();
                }
            }
            catch (const std::exception& e) {
                log->warn("UbooneMagnifyTrackingVisitor: failed to get dead channels for APA={}, face={}, plane={}: {}",
                         apa, face, pind, e.what());
            }
        }
    }

    log->debug("UbooneMagnifyTrackingVisitor: wrote {} entries to T_bad_ch", tree->GetEntries());
}

void Root::UbooneMagnifyTrackingVisitor::write_proj_data(TFile* output_tf, Clus::Facade::Grouping& grouping) const
{
    auto tf = grouping.get_track_fitting();
    if (!tf) {
        log->warn("UbooneMagnifyTrackingVisitor: no TrackFitting in grouping");
        return;
    }

    const auto& fitted = tf->get_fitted_charge_2d();
    if (fitted.empty()) {
        log->warn("UbooneMagnifyTrackingVisitor: fitted_charge_2d is empty");
        return;
    }

    // Get ticks-per-slice map for time_slice conversion
    auto nticks_map = grouping.get_nticks_per_slice();

    // Reorganize fitted charge data by cluster_id
    std::map<int, std::vector<int>> cluster_channels;
    std::map<int, std::vector<int>> cluster_time_slices;
    std::map<int, std::vector<int>> cluster_charges;
    std::map<int, std::vector<int>> cluster_charge_errs;
    std::map<int, std::vector<int>> cluster_charge_preds;

    for (const auto& [afp, wt_map] : fitted) {
        int apa = std::get<0>(afp);
        int face = std::get<1>(afp);
        int plane_idx = std::get<2>(afp);
        // uBooNE channel convention: U=wire, V=2400+wire, W=4800+wire
        int ch_offset = (plane_idx == 1) ? 2400 : (plane_idx == 2) ? 4800 : 0;

        int nticks_per_slice = nticks_map.at(apa).at(face);

        for (const auto& [wt, fc] : wt_map) {
            int wire = wt.first;
            int time = wt.second / nticks_per_slice;
            int channel = ch_offset + wire;

            for (auto* cl : fc.clusters) {
                int cid = cl->get_cluster_id();
                cluster_channels[cid].push_back(channel);
                cluster_time_slices[cid].push_back(time);
                cluster_charges[cid].push_back(static_cast<int>(fc.charge));
                cluster_charge_errs[cid].push_back(static_cast<int>(fc.charge_err));
                cluster_charge_preds[cid].push_back(static_cast<int>(fc.pred_charge));
            }
        }
    }

    // Build vectors in cluster_id order
    std::vector<int> v_cluster_id;
    std::vector<std::vector<int>> v_channel;
    std::vector<std::vector<int>> v_time_slice;
    std::vector<std::vector<int>> v_charge;
    std::vector<std::vector<int>> v_charge_err;
    std::vector<std::vector<int>> v_charge_pred;

    for (const auto& [cid, chs] : cluster_channels) {
        v_cluster_id.push_back(cid);
        v_channel.push_back(chs);
        v_time_slice.push_back(cluster_time_slices[cid]);
        v_charge.push_back(cluster_charges[cid]);
        v_charge_err.push_back(cluster_charge_errs[cid]);
        v_charge_pred.push_back(cluster_charge_preds[cid]);
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

    log->debug("UbooneMagnifyTrackingVisitor: wrote T_proj_data with {} clusters", v_cluster_id.size());
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
