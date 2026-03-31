#include "WireCellRoot/UbooneMagnifyTrackingVisitor.h"

#include "TFile.h"
#include "TTree.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRSegment.h"

#include <boost/graph/graph_traits.hpp>

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
    m_dQdx_scale = get<double>(cfg, "dQdx_scale", 0.1);
    m_dQdx_offset = get<double>(cfg, "dQdx_offset", -1000);
    m_flag_skip_vertex = get<bool>(cfg, "flag_skip_vertex", false);

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
    cfg["dQdx_scale"] = 0.1;
    cfg["dQdx_offset"] = -1000;
    cfg["flag_skip_vertex"] = false;
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
    write_trun(output_tf);
    write_proj_data(output_tf, grouping);
    write_t_rec_data(output_tf, grouping);

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

void Root::UbooneMagnifyTrackingVisitor::write_trun(TFile* output_tf) const
{
    TTree* tree = new TTree("Trun", "Trun");
    tree->SetDirectory(output_tf);

    int eventNo = m_eventNo;
    int runNo = m_runNo;
    int subRunNo = m_subRunNo;
    double dQdx_scale = m_dQdx_scale;
    double dQdx_offset = m_dQdx_offset;

    tree->Branch("eventNo", &eventNo, "eventNo/I");
    tree->Branch("runNo", &runNo, "runNo/I");
    tree->Branch("subRunNo", &subRunNo, "subRunNo/I");
    tree->Branch("dQdx_scale", &dQdx_scale, "dQdx_scale/D");
    tree->Branch("dQdx_offset", &dQdx_offset, "dQdx_offset/D");

    tree->Fill();

    log->debug("UbooneMagnifyTrackingVisitor: wrote Trun with dQdx_scale={}, dQdx_offset={}", dQdx_scale, dQdx_offset);
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

void Root::UbooneMagnifyTrackingVisitor::write_t_rec_data(TFile* output_tf, Clus::Facade::Grouping& grouping) const
{
    using namespace WireCell::Clus;
    
    auto tf = grouping.get_track_fitting();
    if (!tf) {
        log->warn("UbooneMagnifyTrackingVisitor: no TrackFitting in grouping");
        return;
    }

    auto graph = tf->get_graph();
    if (!graph) {
        log->warn("UbooneMagnifyTrackingVisitor: no Graph in TrackFitting");
        return;
    }

    // Get nticks_per_slice for time scaling
    auto nticks_map = grouping.get_nticks_per_slice();
    int nticks_per_slice = 1;  // default value
    if (!nticks_map.empty()) {
        // Use the first available nticks_per_slice value (for uBooNE single APA)
        nticks_per_slice = nticks_map.begin()->second.begin()->second;
    }

    // Create the point tree structure 
    WCPointTree point_tree;
    
    // Create TTree with branches
    TTree* t_rec_charge = new TTree("T_rec_charge", "T_rec_charge");
    t_rec_charge->SetDirectory(output_tf);
    t_rec_charge->Branch("x", &point_tree.reco_x, "x/D");
    t_rec_charge->Branch("y", &point_tree.reco_y, "y/D");
    t_rec_charge->Branch("z", &point_tree.reco_z, "z/D");
    t_rec_charge->Branch("q", &point_tree.reco_dQ, "q/D");
    t_rec_charge->Branch("nq", &point_tree.reco_dx, "nq/D");
    t_rec_charge->Branch("chi2", &point_tree.reco_chi2, "chi2/D");
    t_rec_charge->Branch("ndf", &point_tree.reco_ndf, "ndf/D");
    t_rec_charge->Branch("pu", &point_tree.reco_pu, "pu/D");
    t_rec_charge->Branch("pv", &point_tree.reco_pv, "pv/D");
    t_rec_charge->Branch("pw", &point_tree.reco_pw, "pw/D");
    t_rec_charge->Branch("pt", &point_tree.reco_pt, "pt/D");
    t_rec_charge->Branch("reduced_chi2", &point_tree.reco_reduced_chi2, "reduced_chi2/D");
    t_rec_charge->Branch("flag_vertex", &point_tree.reco_flag_vertex, "flag_vertex/I");
    t_rec_charge->Branch("flag_shower", &point_tree.reco_flag_track_shower, "flag_shower/I");
    t_rec_charge->Branch("rr", &point_tree.reco_rr, "rr/D");
    t_rec_charge->Branch("cluster_id", &point_tree.reco_mother_cluster_id, "cluster_id/I");
    // Keep compatibility with legacy format where real/sub cluster ids
    // both carry the per-segment encoded proto id.
    t_rec_charge->Branch("real_cluster_id", &point_tree.reco_proto_cluster_id, "real_cluster_id/I");
    t_rec_charge->Branch("sub_cluster_id", &point_tree.reco_proto_cluster_id, "sub_cluster_id/I");
    t_rec_charge->Branch("particle_id", &point_tree.reco_particle_id, "particle_id/I");

    // Use calibration parameters from configuration
    const double dQdx_scale = m_dQdx_scale;
    const double dQdx_offset = m_dQdx_offset;
    const bool flag_skip_vertex = m_flag_skip_vertex;

    // Set default values
    point_tree.reco_chi2 = 1;

    // Collect all clusters from the graph
    std::set<Facade::Cluster*> all_clusters;
    
    // Iterate through all edges (segments) in the graph
    auto edge_range = boost::edges(*graph);
    for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
        auto seg = (*graph)[*eit].segment;
        if (seg && seg->cluster()) {
            all_clusters.insert(seg->cluster());
        }
    }

    // Iterate through all vertices in the graph
    auto vertex_range = boost::vertices(*graph);
    for (auto vit = vertex_range.first; vit != vertex_range.second; ++vit) {
        auto vtx = (*graph)[*vit].vertex;
        if (vtx && vtx->cluster()) {
            all_clusters.insert(vtx->cluster());
        }
    }

    // Find the main cluster ID
    int mother_cluster_id = -1;
    for (auto* cluster : all_clusters) {
        if (cluster && cluster->get_flag(Facade::Flags::main_cluster)) {
            mother_cluster_id = cluster->get_cluster_id();
            break;
        }
    }

    // Process each cluster
    for (auto* cluster : all_clusters) {
        if (!cluster) continue;

        point_tree.reco_mother_cluster_id = mother_cluster_id;

        // Process vertices in this cluster
        if (!flag_skip_vertex) {
            for (auto vit = vertex_range.first; vit != vertex_range.second; ++vit) {
                auto vtx = (*graph)[*vit].vertex;
                if (!vtx || vtx->cluster() != cluster) continue;

                // Fill vertex information
                point_tree.reco_cluster_id = cluster->get_cluster_id();
                point_tree.reco_proto_cluster_id = -1;
                point_tree.reco_particle_id = -1;
                point_tree.reco_ndf = cluster->get_cluster_id();
                point_tree.reco_flag_vertex = 1;
                point_tree.reco_flag_track_shower = 0;

                // Position from fit
                const auto& fit_pt = vtx->fit().point;
                point_tree.reco_x = fit_pt.x() / units::cm;
                point_tree.reco_y = fit_pt.y() / units::cm;
                point_tree.reco_z = fit_pt.z() / units::cm;
                
                // Charge and step size from fit
                point_tree.reco_dQ = vtx->fit().dQ * dQdx_scale + dQdx_offset;
                point_tree.reco_dx = vtx->fit().dx / units::cm;
                
                // Projection coordinates
                point_tree.reco_pu = vtx->fit().pu;
                point_tree.reco_pv = vtx->fit().pv;
                point_tree.reco_pw = vtx->fit().pw;
                point_tree.reco_pt = vtx->fit().pt / nticks_per_slice;
                
                point_tree.reco_reduced_chi2 = vtx->fit().reduced_chi2;
                point_tree.reco_rr = -1; // no residual range for vertices

                t_rec_charge->Fill();
            }
        }

        // Process segments in this cluster
        for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
            auto seg = (*graph)[*eit].segment;
            if (!seg || seg->cluster() != cluster) continue;

            const auto& fits = seg->fits();
            const auto& wcpts = seg->wcpts();
            
            if (fits.empty() || wcpts.empty()) continue;

            point_tree.reco_cluster_id = cluster->get_cluster_id();   // cluster id ...
            point_tree.reco_ndf = cluster->get_cluster_id();   
            point_tree.reco_proto_cluster_id = seg->cluster()->get_cluster_id() * 1000 + static_cast<int>(seg->get_graph_index());
            point_tree.reco_flag_vertex = 0;

            // Determine track/shower flag
            bool is_shower = seg->flags_any(PR::SegmentFlags::kShowerTrajectory) || 
                           seg->flags_any(PR::SegmentFlags::kShowerTopology);
            point_tree.reco_flag_track_shower = is_shower ? 1 : 0;
            
            // Prefer per-segment particle hypothesis when available.
            point_tree.reco_particle_id = seg->has_particle_info() ? seg->particle_info()->pdg()
                                                                    : (is_shower ? 1 : 4);

            // Calculate residual range vector
            std::vector<double> rr_vec(fits.size(), 0);
            {
                std::vector<double> L(fits.size(), 0);
                double acc_length = 0;
                
                for (size_t i = 0; i + 1 < fits.size(); i++) {
                    const auto& p1 = fits[i].point;
                    const auto& p2 = fits[i+1].point;
                    double step = std::sqrt(std::pow(p2.x() - p1.x(), 2) +
                                          std::pow(p2.y() - p1.y(), 2) +
                                          std::pow(p2.z() - p1.z(), 2));
                    acc_length += step;
                    L[i+1] = acc_length;
                }

                // Direction sign determines order
                int dirsign = seg->dirsign();
                if (dirsign == 1) { // forward direction
                    for (size_t i = 0; i < fits.size(); i++) {
                        rr_vec[fits.size() - 1 - i] = L.back() - L[fits.size() - 1 - i];
                    }
                } else if (dirsign == -1) { // reverse direction
                    rr_vec = L;
                } else { // unknown direction
                    rr_vec = L;
                }

                // Find vertices connected to this segment
                auto [start_vtx, end_vtx] = PR::find_vertices(*graph, seg);
                
                // Check if vertices have multiple connections
                if (start_vtx) {
                    auto start_degree = boost::out_degree(start_vtx->get_descriptor(), *graph);
                    if (start_degree > 1) rr_vec.front() = -1;
                }
                if (end_vtx) {
                    auto end_degree = boost::out_degree(end_vtx->get_descriptor(), *graph);
                    if (end_degree > 1) rr_vec.back() = -1;
                }
            }

            // Fill tree for each point in the segment
            for (size_t i = 0; i < fits.size(); i++) {
                const auto& fit = fits[i];
                
                point_tree.reco_x = fit.point.x() / units::cm;
                point_tree.reco_y = fit.point.y() / units::cm;
                point_tree.reco_z = fit.point.z() / units::cm;
                point_tree.reco_dQ = fit.dQ * dQdx_scale + dQdx_offset;
                point_tree.reco_dx = fit.dx / units::cm;
                point_tree.reco_pu = fit.pu;
                point_tree.reco_pv = fit.pv;
                point_tree.reco_pw = fit.pw;
                point_tree.reco_pt = fit.pt / nticks_per_slice;
                point_tree.reco_reduced_chi2 = fit.reduced_chi2;
                point_tree.reco_rr = rr_vec[i] / units::cm;
                
                t_rec_charge->Fill();
            }
        }
    }

    log->debug("UbooneMagnifyTrackingVisitor: wrote {} entries to T_rec_charge", t_rec_charge->GetEntries());
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
