#include "WireCellRoot/UbooneMagnifyTrackingVisitor.h"

#include <cmath>

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

// uBooNE absolute channel ID offsets per plane: U=[0,2400), V=[2400,4800), W=[4800,8256)
static constexpr int kPlaneChOffset[3] = {0, 2400, 4800};

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
                    chid = ch + kPlaneChOffset[pind];
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
    const auto& per_cluster = tf->get_cluster_fitted_charge_2d();
    if (fitted.empty() && per_cluster.empty()) {
        log->warn("UbooneMagnifyTrackingVisitor: fitted_charge_2d is empty");
        return;
    }

    // Get ticks-per-slice map for time_slice conversion
    auto nticks_map = grouping.get_nticks_per_slice();

    // Reorganize fitted charge data by cluster_id.
    // doc pr/109 §8: accumulate into an ORDERED per-cell map.  A snapshot is
    // keyed by (wire, TICK) and prepare_data's dead-region filler writes one
    // entry per tick (TrackFitting.cxx:896), so nticks_per_slice of them divide
    // down to the same time_slice; emitting them raw would put duplicate
    // (channel, time_slice) keys in one row.  Charge adds (the filler split the
    // blob charge across those ticks), errors add in quadrature, predictions
    // add -- identities for a live cell, which occurs once per slice.
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

    // doc pr/109 §8: emit one row per FITTED cluster, from that cluster's own
    // snapshot -- the prototype's semantics (wire-cell-prod-nue-port.cxx:3272).
    // The old code read the MERGED map (last-writer-wins across clusters) and
    // tagged each cell by BLOB OWNERSHIP, so a row labelled with cluster A
    // could carry cluster B's prediction.  On 5384-6528 the main cluster's
    // whole prediction was lost that way (cid 19: 2141 cells, Σpred = 0).
    using CellMap = std::map<std::pair<int, int>, Cell>;

    // Within ONE snapshot the tick entries of a slice are added up.
    auto accumulate = [&](CellMap& local, int apa, int face, int plane_idx,
                          const Clus::TrackFitting::WireTime& wt,
                          const Clus::TrackFitting::FittedCharge2D& fc) {
        int nticks_per_slice = nticks_map.at(apa).at(face);
        auto& cell = local[{kPlaneChOffset[plane_idx] + wt.first, wt.second / nticks_per_slice}];
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
    };

    // ACROSS snapshots landing in the same row (two live clusters sharing an
    // ident), charge and charge_err are properties of the READOUT and are
    // taken once, not summed -- only the predictions add.
    auto flush = [&](int cid, const CellMap& local) {
        auto& dst = cluster_cells[cid];
        for (const auto& [key, cell] : local) {
            auto it = dst.find(key);
            if (it == dst.end()) dst.emplace(key, cell);
            else it->second.pred += cell.pred;   // readout fields taken once
        }
    };

    if (!per_cluster.empty()) {
        for (const auto& snap : per_cluster) {
            CellMap local;
            for (const auto& [afp, wt_map] : snap.cells)
                for (const auto& [wt, fc] : wt_map)
                    accumulate(local, std::get<0>(afp), std::get<1>(afp), std::get<2>(afp), wt, fc);
            flush(snap.ident, local);
        }
    }
    else {
        // Fallback: a fitter that never ran with a cluster filter has no
        // snapshots.  Emit the merged map under its blob-ownership tags
        // rather than dropping the tree.
        std::map<int, CellMap> fallback;
        for (const auto& [afp, wt_map] : fitted)
            for (const auto& [wt, fc] : wt_map)
                for (auto* cl : fc.clusters)
                    accumulate(fallback[cl->get_cluster_id()],
                               std::get<0>(afp), std::get<1>(afp), std::get<2>(afp), wt, fc);
        for (const auto& [cid, local] : fallback) flush(cid, local);
    }

    // Build vectors in cluster_id order
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

    // Build per-cluster edge/vertex maps in a single pass (O(E+V) instead of O(C*(E+V)))
    using edge_desc = typename boost::graph_traits<PR::Graph>::edge_descriptor;
    using vertex_desc = typename boost::graph_traits<PR::Graph>::vertex_descriptor;
    std::map<Facade::Cluster*, std::vector<edge_desc>> cluster_edges;
    std::map<Facade::Cluster*, std::vector<vertex_desc>> cluster_vertices;

    auto edge_range = boost::edges(*graph);
    for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
        auto seg = (*graph)[*eit].segment;
        if (seg && seg->cluster()) {
            cluster_edges[seg->cluster()].push_back(*eit);
        }
    }

    auto vertex_range = boost::vertices(*graph);
    for (auto vit = vertex_range.first; vit != vertex_range.second; ++vit) {
        auto vtx = (*graph)[*vit].vertex;
        if (vtx && vtx->cluster()) {
            cluster_vertices[vtx->cluster()].push_back(*vit);
        }
    }

    // Find the main cluster ID
    int mother_cluster_id = -1;
    for (const auto& [cluster, _] : cluster_edges) {
        if (cluster && cluster->get_flag(Facade::Flags::main_cluster)) {
            mother_cluster_id = cluster->get_cluster_id();
            break;
        }
    }
    if (mother_cluster_id < 0) {
        for (const auto& [cluster, _] : cluster_vertices) {
            if (cluster && cluster->get_flag(Facade::Flags::main_cluster)) {
                mother_cluster_id = cluster->get_cluster_id();
                break;
            }
        }
    }

    // Collect unique clusters
    std::set<Facade::Cluster*> all_clusters;
    for (const auto& [c, _] : cluster_edges) all_clusters.insert(c);
    for (const auto& [c, _] : cluster_vertices) all_clusters.insert(c);

    // Process each cluster
    for (auto* cluster : all_clusters) {
        if (!cluster) continue;

        point_tree.reco_mother_cluster_id = mother_cluster_id;

        // Process vertices in this cluster
        if (!flag_skip_vertex) {
            for (auto vd : cluster_vertices[cluster]) {
                auto vtx = (*graph)[vd].vertex;
                if (!vtx) continue;

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
                
                // Projection coordinates (absolute channel IDs)
                point_tree.reco_pu = vtx->fit().pu + kPlaneChOffset[0];
                point_tree.reco_pv = vtx->fit().pv + kPlaneChOffset[1];
                point_tree.reco_pw = vtx->fit().pw + kPlaneChOffset[2];
                point_tree.reco_pt = vtx->fit().pt / nticks_per_slice;
                
                point_tree.reco_reduced_chi2 = vtx->fit().reduced_chi2;
                point_tree.reco_rr = -1; // no residual range for vertices

                t_rec_charge->Fill();
            }
        }

        // Process segments in this cluster
        for (auto ed : cluster_edges[cluster]) {
            auto seg = (*graph)[ed].segment;
            if (!seg) continue;

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
                point_tree.reco_pu = fit.pu + kPlaneChOffset[0];
                point_tree.reco_pv = fit.pv + kPlaneChOffset[1];
                point_tree.reco_pw = fit.pw + kPlaneChOffset[2];
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
