// Fork-by-duplication of UbooneMagnifyTrackingVisitor for the SBND
// neutrino-PR chain (see the header for the fork rationale and deltas).
// Structure and tree schema deliberately track the uBooNE original; the
// two-TPC channel scheme is copied from SbndMagnifyTrackingVisitor.

#include "WireCellRoot/SbndPrMagnifyTrackingVisitor.h"
#include "WireCellRoot/UbooneMagnifyTrackingVisitor.h"  // WCPointTree

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
#include "WireCellClus/PRShower.h"   // PR::ClusterPtrCmp / ClusterPtrSet

#include <boost/graph/graph_traits.hpp>

WIRECELL_FACTORY(SbndPrMagnifyTrackingVisitor, WireCell::Root::SbndPrMagnifyTrackingVisitor,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;

// doc pr/94 Phase 4b: the per-bundle neutrino candidates, in row order, or the
// single unnamed slot when the nu_per_bundle knob is off.  Never empty: the
// legacy resolution is returned as a one-element vector (possibly holding
// nullptr, which the callers already test for), so the knob-off path is
// textually the same sequence of operations it was before.
static std::vector<std::shared_ptr<Clus::TrackFitting>>
collect_nu_fitters(Clus::Facade::Grouping& grouping)
{
    std::vector<std::shared_ptr<Clus::TrackFitting>> out;
    for (int i = 0;; ++i) {
        auto tfi = grouping.get_track_fitting("nu" + std::to_string(i));
        if (!tfi) break;
        out.push_back(tfi);
    }
    if (out.empty()) out.push_back(grouping.get_track_fitting());
    return out;
}

Root::SbndPrMagnifyTrackingVisitor::SbndPrMagnifyTrackingVisitor()
  : log(Log::logger("tracking"))
{
}

Root::SbndPrMagnifyTrackingVisitor::~SbndPrMagnifyTrackingVisitor() {}

void Root::SbndPrMagnifyTrackingVisitor::configure(const WireCell::Configuration& cfg)
{
    m_output_filename = get<std::string>(cfg, "output_filename", "tracking-pr.root");
    m_grouping_name = get<std::string>(cfg, "grouping", "live");
    m_runNo = get<int>(cfg, "runNo", 0);
    m_subRunNo = get<int>(cfg, "subRunNo", 0);
    m_eventNo = get<int>(cfg, "eventNo", 0);
    m_dQdx_scale = get<double>(cfg, "dQdx_scale", 0.1);
    m_dQdx_offset = get<double>(cfg, "dQdx_offset", -1000);
    m_flag_skip_vertex = get<bool>(cfg, "flag_skip_vertex", false);
    m_nticks = get<int>(cfg, "nticks", m_nticks);

    auto anode_tns = cfg["anodes"];
    for (auto anode_tn : anode_tns) {
        auto anode = Factory::find_tn<IAnodePlane>(anode_tn.asString());
        m_anodes.push_back(anode);
    }

    m_dv = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());
}

WireCell::Configuration Root::SbndPrMagnifyTrackingVisitor::default_configuration() const
{
    Configuration cfg;
    cfg["output_filename"] = "tracking-pr.root";
    cfg["grouping"] = "live";
    cfg["anodes"] = Json::arrayValue;
    cfg["detector_volumes"] = "";
    cfg["runNo"] = 0;
    cfg["subRunNo"] = 0;
    cfg["eventNo"] = 0;
    cfg["dQdx_scale"] = 0.1;
    cfg["dQdx_offset"] = -1000;
    cfg["flag_skip_vertex"] = false;
    cfg["nticks"] = 3427;
    return cfg;
}

Root::SbndPrMagnifyTrackingVisitor::ChanScheme Root::SbndPrMagnifyTrackingVisitor::chan_scheme() const
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

void Root::SbndPrMagnifyTrackingVisitor::visit(Clus::Facade::Ensemble& ensemble) const
{
    auto groupings = ensemble.with_name(m_grouping_name);
    if (groupings.empty()) {
        log->debug("SbndPrMagnifyTrackingVisitor: no grouping '{}'", m_grouping_name);
        return;
    }

    auto& grouping = *groupings.at(0);

    // Set anodes and detector volumes on the grouping
    grouping.set_anodes(m_anodes);
    grouping.set_detector_volumes(m_dv);

    const auto cs = chan_scheme();
    log->debug("SbndPrMagnifyTrackingVisitor: channel scheme nch=[{},{},{}] base=[{},{},{}]",
               cs.nch[0], cs.nch[1], cs.nch[2], cs.base[0], cs.base[1], cs.base[2]);

    // Open ROOT file
    TFile* output_tf = TFile::Open(m_output_filename.c_str(), "RECREATE");
    if (!output_tf || output_tf->IsZombie()) {
        log->error("SbndPrMagnifyTrackingVisitor: cannot open {}", m_output_filename);
        return;
    }

    write_bad_channels(output_tf, grouping, cs);
    write_trun(output_tf);
    write_proj_data(output_tf, grouping, cs);
    write_t_rec_data(output_tf, grouping, cs);

    // Empty T_proj tree kept for reader compatibility.
    TTree* tree_proj = new TTree("T_proj", "T_proj");
    tree_proj->SetDirectory(output_tf);
    tree_proj->Fill();

    output_tf->Write();
    output_tf->Close();
    delete output_tf;

    log->debug("SbndPrMagnifyTrackingVisitor: wrote {}", m_output_filename);
}

void Root::SbndPrMagnifyTrackingVisitor::write_bad_channels(TFile* output_tf, Clus::Facade::Grouping& grouping,
                                                            const ChanScheme& cs) const
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
                log->warn("SbndPrMagnifyTrackingVisitor: failed to get dead channels for APA={}, face={}, plane={}: {}",
                          apa, face, pind, e.what());
            }
        }
    }

    log->debug("SbndPrMagnifyTrackingVisitor: wrote {} entries to T_bad_ch", tree->GetEntries());
}

void Root::SbndPrMagnifyTrackingVisitor::write_trun(TFile* output_tf) const
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

    log->debug("SbndPrMagnifyTrackingVisitor: wrote Trun with dQdx_scale={}, dQdx_offset={}", dQdx_scale, dQdx_offset);
}

void Root::SbndPrMagnifyTrackingVisitor::write_proj_data(TFile* output_tf, Clus::Facade::Grouping& grouping,
                                                         const ChanScheme& cs) const
{
    // doc pr/94 Phase 4b: accumulate over every per-bundle candidate.  T_proj_data
    // is ONE row of vector-of-vector branches keyed by cluster_id, so the merge
    // has to happen in these maps BEFORE the single Fill() -- calling this
    // function per candidate would instead leave ROOT cycles T_proj_data;1, ;2,
    // ... that uproot silently resolves to the last of.
    auto nu_fitters = collect_nu_fitters(grouping);
    if (!nu_fitters.front()) {
        log->warn("SbndPrMagnifyTrackingVisitor: no TrackFitting in grouping");
        return;
    }
    // Empty only if EVERY candidate is empty: bailing on candidate 0 alone
    // would drop the whole tree and with it a later candidate's projection.
    // One element when the knob is off, so this is the legacy test verbatim.
    bool any_fitted = false;
    for (const auto& t : nu_fitters) {
        // doc pr/109 §8: the rows come from the per-cluster snapshots now, so
        // either store being non-empty means there is something to write.
        if (t && (!t->get_fitted_charge_2d().empty() || !t->get_cluster_fitted_charge_2d().empty())) {
            any_fitted = true;
            break;
        }
    }
    if (!any_fitted) {
        log->warn("SbndPrMagnifyTrackingVisitor: fitted_charge_2d is empty");
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

    // doc pr/109 §8: emit one row per FITTED cluster, from that cluster's own
    // snapshot -- the prototype's semantics (wire-cell-prod-nue-port.cxx:3272,
    // one row straight out of each cluster's proj_data_*_map).  The old code
    // read the MERGED map (get_fitted_charge_2d(), last-writer-wins across
    // clusters) and tagged each cell by BLOB OWNERSHIP (fc.clusters), so a row
    // labelled with cluster A could carry cluster B's prediction -- and near
    // the vertex, where satellite clusters (higher ident, hence later in the
    // merge) hold the main cluster's cells in their padded bounding box and
    // predict 0 there, it usually did.  Measured on SBND 46363 before the fix:
    // main cluster ident 19 kept 44% of its own predicted charge; uBooNE
    // 5384-6528 kept 0%.
    auto emit = [&](int cid, int apa, int face, int plane_idx,
                    const Clus::TrackFitting::WireTime& wt,
                    const Clus::TrackFitting::FittedCharge2D& fc) {
        // Per-(apa,face) ticks-per-slice with a safe default (the uBooNE
        // original's .at() lookup would throw on an unmapped pair).
        int nticks_per_slice = 1;
        auto ait = nticks_map.find(apa);
        if (ait != nticks_map.end()) {
            auto fit_ = ait->second.find(face);
            if (fit_ != ait->second.end()) nticks_per_slice = fit_->second;
        }
        cluster_channels[cid].push_back(cs.global(plane_idx, apa, wt.first));
        cluster_time_slices[cid].push_back(wt.second / nticks_per_slice);
        cluster_charges[cid].push_back(static_cast<int>(fc.charge));
        cluster_charge_errs[cid].push_back(static_cast<int>(fc.charge_err));
        cluster_charge_preds[cid].push_back(static_cast<int>(fc.pred_charge));
    };

    for (const auto& tf : nu_fitters) {
        if (!tf) continue;
        const auto& per_cluster = tf->get_cluster_fitted_charge_2d();
        if (!per_cluster.empty()) {
            for (const auto& snap : per_cluster)
                for (const auto& [afp, wt_map] : snap.cells)
                    for (const auto& [wt, fc] : wt_map)
                        emit(snap.ident, std::get<0>(afp), std::get<1>(afp), std::get<2>(afp), wt, fc);
            continue;
        }
        // Fallback: a fitter that never ran with a cluster filter has no
        // snapshots (e.g. the "stm" holder fed by merge_fitted_charge_2d).
        // Emit the merged map under its blob-ownership tags rather than
        // dropping the tree.
        for (const auto& [afp, wt_map] : tf->get_fitted_charge_2d())
            for (const auto& [wt, fc] : wt_map)
                for (auto* cl : fc.clusters)
                    emit(cl->get_cluster_id(), std::get<0>(afp), std::get<1>(afp), std::get<2>(afp), wt, fc);
    }  // per-candidate

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

    log->debug("SbndPrMagnifyTrackingVisitor: wrote T_proj_data with {} clusters from {} candidate(s)",
               v_cluster_id.size(), nu_fitters.size());
}

void Root::SbndPrMagnifyTrackingVisitor::write_t_rec_data(TFile* output_tf, Clus::Facade::Grouping& grouping,
                                                          const ChanScheme& cs) const
{
    using namespace WireCell::Clus;

    // doc pr/94 Phase 4b: T_rec_charge covers EVERY per-bundle candidate, not
    // just the one in the unnamed slot.  The tree is already multi-row and
    // already joins back through cluster_id, so this needs no schema change --
    // the rows simply continue across candidates.  The loop is INSIDE this
    // function on purpose: calling the function per candidate would call
    // `new TTree("T_rec_charge")` per candidate and leave ROOT cycles
    // T_rec_charge;1, ;2, ... which uproot and every gate here silently
    // resolve to the last of, hiding the duplication (same hazard as
    // UbooneTaggerOutputVisitor's Write() without kOverwrite).
    auto nu_fitters = collect_nu_fitters(grouping);
    auto tf = nu_fitters.front();
    if (!tf) {
        log->warn("SbndPrMagnifyTrackingVisitor: no TrackFitting in grouping");
        return;
    }

    auto graph = tf->get_graph();
    if (!graph) {
        log->warn("SbndPrMagnifyTrackingVisitor: no Graph in TrackFitting");
        return;
    }

    // Per-(apa,face) ticks-per-slice with a safe default; resolved per point
    // from PR::Fit::paf below (the uBooNE original picks the single first
    // entry, correct only for one APA).
    auto nticks_map = grouping.get_nticks_per_slice();
    auto nticks_for = [&nticks_map](const std::pair<int, int>& paf) -> int {
        auto ait = nticks_map.find(paf.first);
        if (ait != nticks_map.end()) {
            auto fit_ = ait->second.find(paf.second);
            if (fit_ != ait->second.end()) return fit_->second;
        }
        return 1;
    };

    // Per-point channel/tick projection from the fit's own (apa,face); a fit
    // with no recorded paf (-1,-1) falls back to APA 0.
    // PR::Fit pu/pv/pw are FRACTIONAL wire coordinates (double) -- keep the
    // fraction, do not truncate through ChanScheme::global (int wire), or the
    // drawn track sits a systematic half-channel below the measured charge in
    // every view (evt 172230 seg 5030: median -0.5 ch, rms 0.27 = floor() of a
    // uniform fraction).  Same convention as UbooneMagnifyTrackingVisitor,
    // which writes fit.pu + kPlaneChOffset keeping the fraction.
    // The (apa,face) fallback has to be used for the TICK divisor too: an
    // unrecorded paf (-1,-1) misses nticks_map, nticks_for returns 1, and the
    // point lands at its raw tick instead of its time slice -- a factor
    // nticks_per_slice (4 on SBND) off the charge it was fitted to, while the
    // channel silently falls back to APA 0.  Counted below; 0 occurrences on
    // evt 172230/444187, so this is a guard, not a change of what is written.
    int n_paf_fallback = 0;
    WCPointTree point_tree;
    auto project_fit = [&](const PR::Fit& fit) {
        const bool have_paf = fit.paf.first >= 0;
        if (!have_paf) ++n_paf_fallback;
        const std::pair<int, int> paf = have_paf ? fit.paf : std::make_pair(0, 0);
        const int apa = paf.first;
        point_tree.reco_pu = cs.base[0] + apa * cs.nch[0] + fit.pu;
        point_tree.reco_pv = cs.base[1] + apa * cs.nch[1] + fit.pv;
        point_tree.reco_pw = cs.base[2] + apa * cs.nch[2] + fit.pw;
        point_tree.reco_pt = fit.pt / nticks_for(paf);
    };

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

    using edge_desc = typename boost::graph_traits<PR::Graph>::edge_descriptor;
    using vertex_desc = typename boost::graph_traits<PR::Graph>::vertex_descriptor;

    // Which candidate each cluster id came from, so an id claimed by two
    // candidates is reported rather than silently merged: cluster_id is the
    // join key the convert app and the GUI use, and a collision would make one
    // candidate's fit points appear inside the other's track block.  Activity
    // ids are disjoint across bundles by construction (a cluster belongs to one
    // flash bundle), so this is a guard, not an expected path.
    std::map<int, size_t> cid_owner;

    // One candidate's worth of rows.  Body is verbatim the pre-pr/94 code with
    // `graph` rebound to this candidate's own graph.
    auto fill_one = [&](const std::shared_ptr<PR::Graph>& graph, size_t cand) {
    // Build per-cluster edge/vertex maps in a single pass (O(E+V) instead of O(C*(E+V)))
    // Cluster-id ordered, NOT pointer ordered (CLAUDE.md determinism rule): these
    // maps and the set below decide the order of the T_rec_charge rows, which is
    // the order the convert app turns into track blocks and the GUI indexes as
    // "cluster index".  With std::less<Cluster*> that index moved from run to run
    // for the same cluster id (evt 172230 cluster 5 came out at block 1 in one run
    // and block 5 in the next), even under `setarch -R`.  Row CONTENT is unchanged:
    // the multiset of rows hashes the same before and after.
    std::map<Facade::Cluster*, std::vector<edge_desc>, PR::ClusterPtrCmp> cluster_edges;
    std::map<Facade::Cluster*, std::vector<vertex_desc>, PR::ClusterPtrCmp> cluster_vertices;

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

    // Order the rows of each cluster by CONTENT, not by graph descriptor.  The PR
    // graph's vertex descriptors are not reproducible run to run (evt 172230: the
    // 87 vertex rows come out in a different order every run, even under
    // `setarch -R`, while their multiset is identical), and that order is what the
    // convert app turns into blocks and the GUI shows as "cluster index".
    // Segments are keyed on their encoded proto id first so the id stays ascending.
    auto point_key = [](const PR::Fit& f) {
        return std::make_tuple(f.point.x(), f.point.y(), f.point.z());
    };
    for (auto& [cluster, vds] : cluster_vertices) {
        (void)cluster;
        std::sort(vds.begin(), vds.end(), [&](vertex_desc a, vertex_desc b) {
            auto va = (*graph)[a].vertex, vb = (*graph)[b].vertex;
            if (!va || !vb) return va < vb;
            return point_key(va->fit()) < point_key(vb->fit());
        });
    }
    for (auto& [cluster, eds] : cluster_edges) {
        (void)cluster;
        std::sort(eds.begin(), eds.end(), [&](edge_desc a, edge_desc b) {
            auto sa = (*graph)[a].segment, sb = (*graph)[b].segment;
            if (!sa || !sb) return sa < sb;
            const auto ia = sa->get_graph_index(), ib = sb->get_graph_index();
            if (ia != ib) return ia < ib;
            if (sa->fits().empty() || sb->fits().empty()) return sa->fits().size() < sb->fits().size();
            return point_key(sa->fits().front()) < point_key(sb->fits().front());
        });
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
    PR::ClusterPtrSet all_clusters;
    for (const auto& [c, _] : cluster_edges) all_clusters.insert(c);
    for (const auto& [c, _] : cluster_vertices) all_clusters.insert(c);

    // Process each cluster
    for (auto* cluster : all_clusters) {
        if (!cluster) continue;

        {
            auto [it, fresh] = cid_owner.emplace(cluster->get_cluster_id(), cand);
            if (!fresh && it->second != cand) {
                log->warn("SbndPrMagnifyTrackingVisitor: cluster {} appears in both candidate {} and {}"
                          " -- T_rec_charge rows for it are no longer separable by cluster_id",
                          cluster->get_cluster_id(), it->second, cand);
            }
        }

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

                // Projection coordinates (global channel IDs, per-point APA)
                project_fit(vtx->fit());

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
                project_fit(fit);
                point_tree.reco_reduced_chi2 = fit.reduced_chi2;
                point_tree.reco_rr = rr_vec[i] / units::cm;

                t_rec_charge->Fill();
            }
        }
    }
    };  // fill_one

    fill_one(graph, 0);
    for (size_t ci = 1; ci < nu_fitters.size(); ++ci) {
        auto g = nu_fitters[ci] ? nu_fitters[ci]->get_graph() : nullptr;
        if (!g) {
            log->warn("SbndPrMagnifyTrackingVisitor: candidate {} has no Graph; no T_rec_charge rows for it", ci);
            continue;
        }
        fill_one(g, ci);
    }

    log->debug("SbndPrMagnifyTrackingVisitor: wrote {} entries to T_rec_charge from {} candidate(s)"
               " ({} with no recorded (apa,face))",
               t_rec_charge->GetEntries(), nu_fitters.size(), n_paf_fallback);
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
