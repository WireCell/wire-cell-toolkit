// ai-helper issue 41: reco <-> truth matching trees of SbndPrMagnifyTrackingVisitor
// (knob truth_match; MC only -- nothing is written without the labeler's
// truth_track table, i.e. on data and on step-1 tars made without truth_detail).
//
//   T_truth_track  one row per G4 track with deposits in the readout window
//                  (the labeler's truth_track table): nu_idx, bee/ancestor/prong
//                  trackids, pdg, q_readout, q_in_blobs, start/end.
//   T_overlap      the sufficient statistics: one row per (reco piece, trackid)
//                  with the points and charge they share.  A reco piece is a
//                  candidate's segment (nu_index, cluster_id, segment_id,
//                  particle_id, kine_index, in_enu, in_pf) or the rest of a
//                  cluster (nu_index = segment_id = -1); trackid -1 = the
//                  piece's unlabeled points.  Points = the default-scope
//                  imaging points of every cluster, charge = the Bee per-point
//                  charge (mean over the non-dead planes).  n_pure/q_pure: the
//                  part from blobs whose dominant track has >= pure_blob_min of
//                  the blob's true charge.
//   T_match        the SPINE-style answer (TruthMatch.h) at four levels
//                  (truth objects not from a neutrino and holding no reco point
//                  are not written; obj keeps the index of the full list):
//                  0 cluster, 1 particle (T_kine row) vs prong, 2 candidate
//                  (T_kine particles, in_enu) vs interaction, 3 candidate
//                  (particle-flow tree, in_pf) vs interaction; one row per
//                  reco (side 0) or truth (side 1) object, top-k partners.
//   T_pf           the reco particle-flow tree of every candidate (the Bee "mc"
//                  reco nodes): one row per node with its parent.
//   T_blob_truth   (knob blob_truth) one row per blob: its truth label detail.
#include "WireCellRoot/SbndPrMagnifyTrackingVisitor.h"

#include "WireCellClus/TruthMatch.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRShower.h"
#include "WireCellUtil/Units.h"

#include "TFile.h"
#include "TTree.h"

#include <cstdlib>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus;

namespace {

    int seg_display_id(const PR::SegmentPtr& seg)
    {
        if (!seg) return -1;
        int sid = seg->id();
        if (sid < 0) sid = static_cast<int>(seg->get_graph_index());
        const auto* cl = seg->cluster();
        return cl ? cl->get_cluster_id() * 1000 + sid : sid;
    }

    struct TrackRow {
        int nu_idx{-1}, bee{-1}, ancestor{-1}, prong{-1}, pdg{0};
        double q_readout{0}, q_in_blobs{0};
    };

    struct PieceKey {
        int nu{-1}, cluster{-1}, segment{-1}, particle{-1}, kine{-1}, enu{0}, pf{0};
        auto tie() const { return std::tie(nu, cluster, segment, particle, kine, enu, pf); }
        bool operator<(const PieceKey& o) const { return tie() < o.tie(); }
    };

    struct Acc {
        double n{0}, q{0}, n_pure{0}, q_pure{0};
    };

    int pdg_from_pf_name(const std::string& name)
    {
        static const std::map<std::string, int> m = {
            {"e-", 11}, {"e+", -11}, {"mu-", 13}, {"mu+", -13}, {"gamma", 22},
            {"pi+", 211}, {"pi-", -211}, {"p", 2212}, {"proton", 2212},
            {"n", 2112}, {"neutron", 2112}, {"K+", 321}, {"K-", -321}, {"pi0", 111}};
        auto it = m.find(name);
        if (it != m.end()) return it->second;
        char* end = nullptr;
        const long v = std::strtol(name.c_str(), &end, 10);
        return (end && *end == '\0' && !name.empty()) ? static_cast<int>(v) : 0;
    }

    void collect_pf_ids(const Json::Value& node, std::set<int>& ids)
    {
        if (node.isMember("id")) ids.insert(node["id"].asInt());
        for (const auto& c : node["children"]) collect_pf_ids(c, ids);
    }
}

void Root::SbndPrMagnifyTrackingVisitor::write_truth_match(TFile* output_tf, Clus::Facade::Grouping& grouping,
                                                           const Clus::Facade::Ensemble& ensemble) const
{
    const auto tt = ensemble.aux_tensor("truth_track");
    if (!tt) {
        // truth_match is on, so the table was expected: warn (PR 544 review, item 4).
        // Its producer is the step-1 labeler with truth_detail (larwirecell 73c79bc).
        log->warn("SbndPrMagnifyTrackingVisitor: truth_match is on but the input has no truth_track "
                  "table (step 1 without the labeler's truth_detail?): no truth-match trees");
        return;
    }
    const int run = m_evt_runNo, subrun = m_evt_subRunNo, event = m_evt_eventNo;

    // ---- the truth tracks ----
    std::map<int, TrackRow> tracks;
    {
        const auto& md = tt->metadata();
        const auto shape = tt->shape();
        if (shape.size() == 2 && tt->element_type() == typeid(double)) {
            std::map<std::string, size_t> col;
            for (Json::ArrayIndex c = 0; c < md["columns"].size(); ++c) col[md["columns"][c].asString()] = c;
            const double* d = reinterpret_cast<const double*>(tt->data());
            const size_t nc = shape[1];
            auto get = [&](size_t r, const char* k) -> double {
                auto it = col.find(k);
                return it == col.end() ? 0.0 : d[r * nc + it->second];
            };
            for (size_t r = 0; r < shape[0]; ++r) {
                TrackRow t;
                const int tid = static_cast<int>(std::lround(get(r, "trackid")));
                t.nu_idx = static_cast<int>(std::lround(get(r, "nu_idx")));
                t.bee = static_cast<int>(std::lround(get(r, "bee_trackid")));
                t.ancestor = static_cast<int>(std::lround(get(r, "ancestor_trackid")));
                t.prong = static_cast<int>(std::lround(get(r, "prong_trackid")));
                t.pdg = static_cast<int>(std::lround(get(r, "pdg")));
                t.q_readout = get(r, "q_readout");
                t.q_in_blobs = get(r, "q_in_blobs");
                tracks[tid] = t;
            }
        }
        write_truth_table_public(output_tf, "T_truth_track", tt, run, subrun, event);
    }

    // ---- the particle-flow trees (MABC publishes them as metadata of an empty tensor) ----
    std::map<int, Json::Value> pf_by_nu;
    if (const auto pft = ensemble.aux_tensor("pf_tree:" + m_pf_tree_name)) {
        const Json::Value md = pft->metadata();   // by value: never iterate a temporary's member
        for (const auto& c : md["candidates"]) pf_by_nu[c["nu_index"].asInt()] = c["root"];
    }

    // ---- reco pieces: the candidates' segments claim their associate points ----
    std::vector<PieceKey> pieces;
    std::map<PieceKey, int> piece_index;
    auto piece_of = [&](const PieceKey& k) -> int {
        auto it = piece_index.find(k);
        if (it != piece_index.end()) return it->second;
        const int i = static_cast<int>(pieces.size());
        pieces.push_back(k);
        piece_index[k] = i;
        return i;
    };
    std::map<int, std::vector<int>> claim;   // cluster_id -> per point: piece (-1 = none)
    std::map<int, Facade::Cluster*> cluster_by_id;
    for (auto* cl : grouping.children()) {
        cluster_by_id[cl->get_cluster_id()] = cl;
        claim[cl->get_cluster_id()].assign(cl->npoints(), -1);
    }
    std::vector<std::shared_ptr<Clus::TrackFitting>> fitters;
    for (int i = 0;; ++i) {
        auto tfi = grouping.get_track_fitting("nu" + std::to_string(i));
        if (!tfi) break;
        fitters.push_back(tfi);
    }
    if (fitters.empty() && grouping.get_track_fitting()) fitters.push_back(grouping.get_track_fitting());
    size_t n_claim_fallback = 0, n_claim_dup = 0;
    for (size_t i = 0; i < fitters.size(); ++i) {
        const auto& tf = fitters[i];
        auto graph = tf->get_graph();
        if (!graph) continue;
        const int nu = tf->get_tagger_info().nu_index >= 0 ? tf->get_tagger_info().nu_index : static_cast<int>(i);
        const auto& ki = tf->get_kine_info();
        std::map<int, int> kine_by_gi;
        const size_t nlink = std::min(ki.link_seg_graph_index.size(), ki.link_seg_kine_index.size());
        for (size_t k = 0; k < nlink; ++k) kine_by_gi[ki.link_seg_graph_index[k]] = ki.link_seg_kine_index[k];
        std::map<PR::SegmentPtr, PR::ShowerPtr, PR::SegmentIndexCmp> seg_to_shower;
        for (const auto& shower : tf->get_showers()) {
            PR::IndexedVertexSet sv;
            PR::IndexedSegmentSet ss;
            shower->fill_sets(sv, ss, /*flag_exclude_start_segment=*/false);
            for (const auto& seg : ss) seg_to_shower[seg] = shower;
        }
        std::set<int> pf_ids;
        if (auto it = pf_by_nu.find(nu); it != pf_by_nu.end()) {
            for (const auto& c : it->second["children"]) collect_pf_ids(c, pf_ids);
        }
        for (auto ed : PR::ordered_edges(*graph)) {
            auto seg = (*graph)[ed].segment;
            if (!seg) continue;
            const auto* cl = seg->cluster();
            if (!cl) continue;
            const int cid = cl->get_cluster_id();
            auto cit = claim.find(cid);
            if (cit == claim.end()) continue;
            const int gi = static_cast<int>(seg->get_graph_index());
            PieceKey key;
            key.nu = nu;
            key.cluster = cid;
            key.segment = cid * 1000 + gi;
            const int own = seg_display_id(seg);
            auto sit = seg_to_shower.find(seg);
            key.particle = (sit == seg_to_shower.end()) ? own : seg_display_id(sit->second->start_segment());
            auto kit = kine_by_gi.find(gi);
            key.kine = (kit == kine_by_gi.end()) ? -1 : kit->second;
            key.enu = key.kine >= 0 ? 1 : 0;
            key.pf = (pf_ids.count(key.particle) || pf_ids.count(own)) ? 1 : 0;
            const int p = piece_of(key);
            auto& cl_claim = cit->second;
            std::vector<size_t> idx;
            if (seg->has_global_indices("associate_points")) {
                idx = seg->global_indices("associate_points");
            }
            else if (auto dpc = seg->dpcloud("associate_points")) {
                ++n_claim_fallback;
                for (size_t ip = 0; ip < dpc->npoints(); ++ip) idx.push_back(cl->get_closest_point_index(dpc->point3d(ip)));
            }
            for (size_t ix : idx) {
                if (ix >= cl_claim.size()) continue;
                if (cl_claim[ix] < 0) cl_claim[ix] = p;
                else ++n_claim_dup;
            }
        }
    }

    // ---- every imaging point: piece x trackid ----
    std::map<std::pair<int, int>, Acc> acc;   // (piece, trackid) -> sums
    for (auto& [cid, cl_claim] : claim) {
        auto* cl = cluster_by_id[cid];
        PieceKey rest;
        rest.cluster = cid;
        const int prest = piece_of(rest);
        const int n = cl->npoints();
        for (int ip = 0; ip < n; ++ip) {
            const auto* blob = cl->blob_with_point(ip);
            int tid = -1;
            double tq = 0, tqd = 0;
            if (blob) {
                tid = blob->get_scalar<int>("trackid", -1);
                tq = blob->get_scalar<double>("truth_q", 0.0);
                tqd = blob->get_scalar<double>("truth_q_dom", 0.0);
            }
            if (tid >= 0 && !tracks.count(tid)) tid = -1;   // no truth row: treat as unlabeled
            double sum = 0;
            int npl = 0;
            for (int plane = 0; plane < 3; ++plane) {
                if (!cl->is_wire_dead(ip, plane, 1e10)) {
                    sum += cl->charge_value(ip, plane);
                    ++npl;
                }
            }
            const double q = npl > 0 ? sum / npl : 0.0;
            const bool pure = tid >= 0 && tq > 0 && tqd / tq >= m_pure_blob_min;
            auto& a = acc[{cl_claim[ip] >= 0 ? cl_claim[ip] : prest, tid}];
            a.n += 1;
            a.q += q;
            if (pure) {
                a.n_pure += 1;
                a.q_pure += q;
            }
        }
    }

    // ---- T_overlap ----
    {
        TTree* t = new TTree("T_overlap", "reco piece x truth track: shared imaging points and charge");
        t->SetDirectory(output_tf);
        int b_run = run, b_sub = subrun, b_evt = event;
        int nu, cl, seg, part, kine, enu, pf, tid, nu_idx;
        float n_points, q, n_pure, q_pure;
        t->Branch("runNo", &b_run, "runNo/I");
        t->Branch("subRunNo", &b_sub, "subRunNo/I");
        t->Branch("eventNo", &b_evt, "eventNo/I");
        t->Branch("nu_index", &nu, "nu_index/I");
        t->Branch("cluster_id", &cl, "cluster_id/I");
        t->Branch("segment_id", &seg, "segment_id/I");
        t->Branch("particle_id", &part, "particle_id/I");
        t->Branch("kine_index", &kine, "kine_index/I");
        t->Branch("in_enu", &enu, "in_enu/I");
        t->Branch("in_pf", &pf, "in_pf/I");
        t->Branch("trackid", &tid, "trackid/I");
        t->Branch("nu_idx", &nu_idx, "nu_idx/I");
        t->Branch("n_points", &n_points, "n_points/F");
        t->Branch("q", &q, "q/F");
        t->Branch("n_pure", &n_pure, "n_pure/F");
        t->Branch("q_pure", &q_pure, "q_pure/F");
        for (const auto& [pk, a] : acc) {
            const auto& k = pieces[pk.first];
            nu = k.nu; cl = k.cluster; seg = k.segment; part = k.particle; kine = k.kine; enu = k.enu; pf = k.pf;
            tid = pk.second;
            nu_idx = tid < 0 ? -2 : tracks[tid].nu_idx;
            n_points = a.n; q = a.q; n_pure = a.n_pure; q_pure = a.q_pure;
            t->Fill();
        }
    }

    // ---- T_match: four levels ----
    {
        TTree* t = new TTree("T_match", "SPINE-style reco<->truth match (TruthMatch.h) at four levels");
        t->SetDirectory(output_tf);
        int b_run = run, b_sub = subrun, b_evt = event;
        int level, side, obj, r_nu, r_cluster, r_kine, t_nu_idx, t_trackid, mutual, nmatch;
        float size_n, size_q, unl_n, unl_q, t_q_readout, t_q_in_blobs;
        std::vector<int> m_obj;
        std::vector<float> m_iou_n, m_iou_q, m_sh_n, m_sh_q;
        t->Branch("runNo", &b_run, "runNo/I");
        t->Branch("subRunNo", &b_sub, "subRunNo/I");
        t->Branch("eventNo", &b_evt, "eventNo/I");
        t->Branch("level", &level, "level/I");
        t->Branch("side", &side, "side/I");
        t->Branch("obj", &obj, "obj/I");
        t->Branch("nu_index", &r_nu, "nu_index/I");
        t->Branch("cluster_id", &r_cluster, "cluster_id/I");
        t->Branch("kine_index", &r_kine, "kine_index/I");
        t->Branch("truth_nu_idx", &t_nu_idx, "truth_nu_idx/I");
        t->Branch("truth_trackid", &t_trackid, "truth_trackid/I");
        t->Branch("size_n", &size_n, "size_n/F");
        t->Branch("size_q", &size_q, "size_q/F");
        t->Branch("n_unlabeled", &unl_n, "n_unlabeled/F");
        t->Branch("q_unlabeled", &unl_q, "q_unlabeled/F");
        t->Branch("truth_q_readout", &t_q_readout, "truth_q_readout/F");
        t->Branch("truth_q_in_blobs", &t_q_in_blobs, "truth_q_in_blobs/F");
        t->Branch("nmatch", &nmatch, "nmatch/I");
        t->Branch("match_obj", &m_obj);
        t->Branch("match_iou_n", &m_iou_n);
        t->Branch("match_iou_q", &m_iou_q);
        t->Branch("match_shared_n", &m_sh_n);
        t->Branch("match_shared_q", &m_sh_q);
        t->Branch("mutual", &mutual, "mutual/I");

        // reco key: (kind, a, b): kind 0 = cluster remainder / whole cluster (a = cluster_id),
        // 1 = particle (a = nu, b = kine row), 2 = candidate (a = nu).
        using RKey = std::tuple<int, int, int>;
        // truth key: (kind, id): 0 = interaction (id = nu_idx), 1 = cosmic (id = bee trackid),
        // 2 = prong (id = prong trackid).
        using TKey = std::pair<int, int>;
        for (level = 0; level < 4; ++level) {
            auto rkey = [&](const PieceKey& k) -> RKey {
                if (level == 1 && k.kine >= 0) return {1, k.nu, k.kine};
                if (level == 2 && k.enu) return {2, k.nu, -1};
                if (level == 3 && k.pf) return {2, k.nu, -1};
                return {0, k.cluster, -1};
            };
            auto tkey = [&](const TrackRow& tr) -> TKey {
                if (level == 1) return {2, tr.prong};
                if (tr.nu_idx >= 0) return {0, tr.nu_idx};
                return {1, tr.bee};
            };
            std::map<RKey, int> rid;
            std::map<TKey, int> tidx;
            std::vector<RKey> rkeys;
            std::vector<TKey> tkeys;
            std::vector<TruthMatch::Size> rsize, tsize;
            std::vector<double> runl_n, runl_q, tq_readout, tq_in_blobs;
            std::vector<int> t_nu;   // the interaction of a truth object, from its member tracks
            std::vector<TruthMatch::Cell> cells;
            auto get_r = [&](const RKey& k) {
                auto it = rid.find(k);
                if (it != rid.end()) return it->second;
                const int i = static_cast<int>(rkeys.size());
                rid[k] = i;
                rkeys.push_back(k);
                rsize.push_back({});
                runl_n.push_back(0);
                runl_q.push_back(0);
                return i;
            };
            auto get_t = [&](const TKey& k) {
                auto it = tidx.find(k);
                if (it != tidx.end()) return it->second;
                const int i = static_cast<int>(tkeys.size());
                tidx[k] = i;
                tkeys.push_back(k);
                tsize.push_back({});
                tq_readout.push_back(0);
                tq_in_blobs.push_back(0);
                t_nu.push_back(-1);
                return i;
            };
            // every truth object of the event, also those no reco point carries
            for (const auto& [tid0, tr] : tracks) {
                const int ti = get_t(tkey(tr));
                tq_readout[ti] += tr.q_readout;
                tq_in_blobs[ti] += tr.q_in_blobs;
                if (tr.nu_idx >= 0) t_nu[ti] = tr.nu_idx;
            }
            for (const auto& [pk, a] : acc) {
                const int ri = get_r(rkey(pieces[pk.first]));
                rsize[ri].n += a.n;
                rsize[ri].q += a.q;
                if (pk.second < 0) {
                    runl_n[ri] += a.n;
                    runl_q[ri] += a.q;
                    continue;
                }
                const int ti = get_t(tkey(tracks[pk.second]));
                tsize[ti].n += a.n;
                tsize[ti].q += a.q;
                cells.push_back({ri, ti, a.n, a.q});
            }
            const auto res = TruthMatch::match(rsize, tsize, cells, m_truth_match_topk, 0.0);
            auto fill_partners = [&](const std::vector<TruthMatch::Partner>& ps) {
                m_obj.clear(); m_iou_n.clear(); m_iou_q.clear(); m_sh_n.clear(); m_sh_q.clear();
                for (const auto& p : ps) {
                    m_obj.push_back(p.id);
                    m_iou_n.push_back(p.iou_n);
                    m_iou_q.push_back(p.iou_q);
                    m_sh_n.push_back(p.shared_n);
                    m_sh_q.push_back(p.shared_q);
                }
                nmatch = static_cast<int>(ps.size());
            };
            side = 0;
            for (size_t r = 0; r < rkeys.size(); ++r) {
                const auto& [kind, a, b] = rkeys[r];
                obj = static_cast<int>(r);
                r_nu = kind == 0 ? -1 : a;
                r_cluster = kind == 0 ? a : -1;
                r_kine = kind == 1 ? b : -1;
                t_nu_idx = -1; t_trackid = -1;
                size_n = rsize[r].n; size_q = rsize[r].q;
                unl_n = runl_n[r]; unl_q = runl_q[r];
                t_q_readout = 0; t_q_in_blobs = 0;
                fill_partners(res.reco[r]);
                mutual = res.reco_mutual[r];
                t->Fill();
            }
            side = 1;
            for (size_t i = 0; i < tkeys.size(); ++i) {
                const auto& [kind, id] = tkeys[i];
                obj = static_cast<int>(i);
                // A non-neutrino truth object with no reco point can match nothing and is
                // most of the rows (cosmic secondaries); its charge is in T_truth_track.
                // Neutrino objects stay even without points (a missed interaction / prong).
                // obj keeps its index, so match_obj of the reco rows stays valid.
                if (tsize[i].n <= 0 && t_nu[i] < 0) continue;
                r_nu = -1; r_cluster = -1; r_kine = -1;
                t_nu_idx = t_nu[i];   // interaction: its nu_idx; prong: its member tracks' nu_idx; cosmic: -1
                t_trackid = kind == 0 ? -1 : id;
                size_n = tsize[i].n; size_q = tsize[i].q;
                unl_n = 0; unl_q = 0;
                t_q_readout = tq_readout[i]; t_q_in_blobs = tq_in_blobs[i];
                fill_partners(res.truth[i]);
                mutual = res.truth_mutual[i];
                t->Fill();
            }
        }
    }

    // ---- T_pf: the reco particle-flow tree of every candidate ----
    {
        TTree* t = new TTree("T_pf", "reco particle-flow tree (the Bee mc reco nodes), one row per node");
        t->SetDirectory(output_tf);
        int b_run = run, b_sub = subrun, b_evt = event;
        int nu, node_id, parent_id, depth, nchild, pdg, is_pseudo, kine;
        float ke_mev, sx, sy, sz, ex, ey, ez;
        std::string name;
        t->Branch("runNo", &b_run, "runNo/I");
        t->Branch("subRunNo", &b_sub, "subRunNo/I");
        t->Branch("eventNo", &b_evt, "eventNo/I");
        t->Branch("nu_index", &nu, "nu_index/I");
        t->Branch("node_id", &node_id, "node_id/I");
        t->Branch("parent_id", &parent_id, "parent_id/I");
        t->Branch("depth", &depth, "depth/I");
        t->Branch("n_children", &nchild, "n_children/I");
        t->Branch("name", &name);
        t->Branch("pdg", &pdg, "pdg/I");
        t->Branch("is_pseudo", &is_pseudo, "is_pseudo/I");   // 1: no reco segment / shower behind the node
        t->Branch("kine_index", &kine, "kine_index/I");
        t->Branch("ke_mev", &ke_mev, "ke_mev/F");
        t->Branch("start_x", &sx, "start_x/F");
        t->Branch("start_y", &sy, "start_y/F");
        t->Branch("start_z", &sz, "start_z/F");
        t->Branch("end_x", &ex, "end_x/F");
        t->Branch("end_y", &ey, "end_y/F");
        t->Branch("end_z", &ez, "end_z/F");
        // node id -> T_kine row, from the pieces (particle_id is the node id of a particle)
        std::map<std::pair<int, int>, int> kine_of;
        std::set<std::pair<int, int>> reco_particle;   // (nu, node id) backed by a reco segment / shower
        for (const auto& k : pieces) {
            if (k.nu < 0) continue;
            reco_particle.insert({k.nu, k.particle});
            reco_particle.insert({k.nu, k.segment});
            if (k.kine >= 0) kine_of[{k.nu, k.particle}] = k.kine;
        }
        std::function<void(const Json::Value&, int, int, int)> walk = [&](const Json::Value& n, int cand, int parent, int d) {
            nu = cand;
            node_id = n["id"].asInt();
            parent_id = parent;
            depth = d;
            nchild = static_cast<int>(n["children"].size());
            const std::string text = n["text"].asString();
            const auto sp = text.find("  ");
            name = sp == std::string::npos ? text : text.substr(0, sp);
            ke_mev = -1;
            if (sp != std::string::npos) ke_mev = std::strtof(text.c_str() + sp, nullptr);
            pdg = pdg_from_pf_name(name);
            // a pseudo-particle (gamma / neutron / pi0 carrier) is no reco segment or shower
            is_pseudo = reco_particle.count({cand, node_id}) ? 0 : 1;
            auto kit = kine_of.find({cand, node_id});
            kine = kit == kine_of.end() ? -1 : kit->second;
            const auto& dj = n["data"];
            sx = dj["start"][0].asFloat(); sy = dj["start"][1].asFloat(); sz = dj["start"][2].asFloat();
            ex = dj["end"][0].asFloat(); ey = dj["end"][1].asFloat(); ez = dj["end"][2].asFloat();
            t->Fill();
            const int me = node_id;
            for (const auto& c : n["children"]) walk(c, cand, me, d + 1);
        };
        for (const auto& [cand, root] : pf_by_nu) {
            for (const auto& c : root["children"]) walk(c, cand, -1, 0);
        }
    }

    // ---- T_blob_truth (optional) ----
    if (m_blob_truth) {
        TTree* t = new TTree("T_blob_truth", "one row per blob: its truth label");
        t->SetDirectory(output_tf);
        int b_run = run, b_sub = subrun, b_evt = event;
        int cl, npoints, tid, tid2;
        float charge, tq, tqd, tq2;
        t->Branch("runNo", &b_run, "runNo/I");
        t->Branch("subRunNo", &b_sub, "subRunNo/I");
        t->Branch("eventNo", &b_evt, "eventNo/I");
        t->Branch("cluster_id", &cl, "cluster_id/I");
        t->Branch("npoints", &npoints, "npoints/I");
        t->Branch("charge", &charge, "charge/F");
        t->Branch("trackid", &tid, "trackid/I");
        t->Branch("truth_q", &tq, "truth_q/F");
        t->Branch("truth_q_dom", &tqd, "truth_q_dom/F");
        t->Branch("trackid2", &tid2, "trackid2/I");
        t->Branch("truth_q2", &tq2, "truth_q2/F");
        for (auto* c : grouping.children()) {
            cl = c->get_cluster_id();
            for (const auto* blob : c->children()) {
                npoints = blob->npoints();
                charge = blob->charge();
                tid = blob->get_scalar<int>("trackid", -1);
                tq = blob->get_scalar<double>("truth_q", 0.0);
                tqd = blob->get_scalar<double>("truth_q_dom", 0.0);
                tid2 = blob->get_scalar<int>("trackid2", -1);
                tq2 = blob->get_scalar<double>("truth_q2", 0.0);
                t->Fill();
            }
        }
    }

    log->debug("SbndPrMagnifyTrackingVisitor: truth match: {} tracks, {} pieces, {} overlap rows, "
               "{} pf candidates; {} segments mapped by closest point, {} points claimed twice",
               tracks.size(), pieces.size(), acc.size(), pf_by_nu.size(), n_claim_fallback, n_claim_dup);
}
