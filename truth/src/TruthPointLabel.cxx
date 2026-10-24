#include "WireCellTruth/TruthPointLabel.h"
#include "WireCellTruth/Schema.h"

#include "WireCellAux/ClusterHelpers.h"
#include "WireCellIface/WirePlaneId.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/String.h"

#include <algorithm>
#include <map>
#include <utility>
#include <cmath>
#include <limits>

using namespace WireCell;
using WireCell::PointCloud::Array;
using WireCell::PointCloud::Dataset;
using node_t = WireCell::PointCloud::Tree::Points::node_t;

std::map<int, std::shared_ptr<Truth::DepoPolyline>>
Truth::build_polylines(const IDepo::vector& depos, const DepoBlobAssocParams& assoc,
                       const depo_position_f& position, bool tie_by_direction)
{
    std::map<int, std::shared_ptr<DepoPolyline>> out;
    for (const auto& d : depos) {
        if (!d) continue;
        const IDepo* key = d.get();
        if (assoc.id_from_prior && d->prior()) key = d->prior().get();
        int tid = key->id();
        double energy = key->energy();
        if (assoc.energy_map && tid >= 0 && (size_t)tid < assoc.energy_map->size()) energy = (*assoc.energy_map)[tid];
        if (assoc.id_map && tid >= 0 && (size_t)tid < assoc.id_map->size()) tid = (*assoc.id_map)[tid];
        if (assoc.fold_negative_ids && tid < 0) tid = -tid;   // as the association: larg4's "-ancestor"
        auto& pl = out[tid];
        if (!pl) pl = std::make_shared<DepoPolyline>();
        pl->add(position(*d), key->time(), energy);
    }
    for (auto& [tid, pl] : out) {
        if (tie_by_direction) pl->set_tie_by_direction(true);
        pl->build();
    }
    return out;
}

Truth::particle_polylines_t
Truth::build_particle_polylines(const IDepo::vector& depos, const DepoBlobAssocParams& assoc,
                                const depo_position_f& position, const std::vector<int>* orig_ids,
                                bool tie_by_direction)
{
    // key -> particle -> polyline; the keys as build_polylines makes them
    std::map<int, std::map<int, std::shared_ptr<DepoPolyline>>> byp;
    for (const auto& d : depos) {
        if (!d) continue;
        const IDepo* key = d.get();
        if (assoc.id_from_prior && d->prior()) key = d->prior().get();
        const int raw = key->id();
        int tid = raw;
        double energy = key->energy();
        if (assoc.energy_map && tid >= 0 && (size_t)tid < assoc.energy_map->size()) energy = (*assoc.energy_map)[tid];
        if (assoc.id_map && tid >= 0 && (size_t)tid < assoc.id_map->size()) tid = (*assoc.id_map)[tid];
        if (assoc.fold_negative_ids && tid < 0) tid = -tid;
        int part = raw;
        if (orig_ids && raw >= 0 && (size_t)raw < orig_ids->size()) part = (*orig_ids)[raw];
        auto& pl = byp[tid][part];
        if (!pl) pl = std::make_shared<DepoPolyline>();
        pl->add(position(*d), key->time(), energy);
    }
    particle_polylines_t out;
    for (auto& [tid, parts] : byp) {
        auto& v = out[tid];
        for (auto& [part, pl] : parts) {
            if (tie_by_direction) pl->set_tie_by_direction(true);
            pl->build();
            v.push_back(pl);
        }
    }
    return out;
}

namespace {

    int scalar_int(const Dataset& scalar, const std::string& key, int def = std::numeric_limits<int>::min())
    {
        auto arr = scalar.get(key);
        if (!arr || arr->size_major() == 0) return def;
        return arr->element<int>(0);
    }

    // Does the node's scalar PC describe this blob?  The same keys and
    // arithmetic as Aux::fill_scalar_blob.
    bool node_matches_blob(const node_t& node, const IBlob& blob, double tick, std::string& why)
    {
        const auto& lpcs = node.value.local_pcs();
        auto it = lpcs.find("scalar");
        if (it == lpcs.end()) { why = "no scalar PC"; return false; }
        const auto& sc = it->second;
        WirePlaneId wpid(kAllLayers, blob.face()->which(), blob.face()->anode());
        if (scalar_int(sc, "wpid") != (int)wpid.ident()) { why = "wpid"; return false; }
        const auto islice = blob.slice();
        if (scalar_int(sc, "slice_index_min") != (int)(islice->start() / tick)) { why = "slice_index_min"; return false; }
        static const std::map<int, std::string> layer_names = {{2, "u"}, {3, "v"}, {4, "w"}};
        for (const auto& strip : blob.shape().strips()) {
            auto ln = layer_names.find(strip.layer);
            if (ln == layer_names.end()) continue;
            if (scalar_int(sc, ln->second + "_wire_index_min") != (int)strip.bounds.first ||
                scalar_int(sc, ln->second + "_wire_index_max") != (int)strip.bounds.second) {
                why = ln->second + "_wire_index";
                return false;
            }
        }
        return true;
    }

    struct Columns {
        std::vector<int> tid, orig_tid, instance_valid, pdg, cls, iid, pf_node, end;
        std::vector<double> pid_frac, true_q_frac, d_traj, dir_x, dir_y, dir_z, s_t,
            v_vtx, v_end, v_kink, v_conv, dedx;
        // end_distance: the distance to the own node's start / end, -1 when there is none
        bool with_end_distance{false};
        std::vector<double> d_start, d_end;
        // heat_by_pf_node: the distance to the nearest true vertex / conversion point, -1 when there is none
        bool with_node_heat{false};
        std::vector<double> d_vtx, d_conv;
        void resize(size_t n)
        {
            if (with_end_distance) { d_start.assign(n, -1.0); d_end.assign(n, -1.0); }
            if (with_node_heat) { d_vtx.assign(n, -1.0); d_conv.assign(n, -1.0); }
            for (auto* v : {&tid, &orig_tid, &instance_valid, &pdg, &cls, &iid, &pf_node, &end}) v->assign(n, 0);
            for (auto* v : {&pid_frac, &true_q_frac, &d_traj, &dir_x, &dir_y, &dir_z, &s_t,
                            &v_vtx, &v_end, &v_kink, &v_conv, &dedx}) v->assign(n, 0.0);
        }
        Dataset dataset() const
        {
            Dataset ds;
            // the order of Schema::truth_pc_keys()
            ds.add("tid", Array(tid));
            ds.add("orig_tid", Array(orig_tid));
            ds.add("pid_frac", Array(pid_frac));
            ds.add("true_q_frac", Array(true_q_frac));
            ds.add("instance_valid", Array(instance_valid));
            ds.add("pdg", Array(pdg));
            ds.add("cls", Array(cls));
            ds.add("iid", Array(iid));
            ds.add("pf_node", Array(pf_node));
            ds.add("d_traj", Array(d_traj));
            ds.add("dir_x", Array(dir_x));
            ds.add("dir_y", Array(dir_y));
            ds.add("dir_z", Array(dir_z));
            ds.add("s_t", Array(s_t));
            ds.add("v_heat_vtx", Array(v_vtx));
            ds.add("v_heat_end", Array(v_end));
            ds.add("v_heat_kink", Array(v_kink));
            ds.add("v_heat_conv", Array(v_conv));
            ds.add("end", Array(end));
            ds.add("dedx_true", Array(dedx));
            if (with_end_distance) {
                ds.add("d_node_start", Array(d_start));
                ds.add("d_node_end", Array(d_end));
            }
            if (with_node_heat) {
                ds.add("d_vertex", Array(d_vtx));
                ds.add("d_conversion", Array(d_conv));
            }
            return ds;
        }
    };
}

Truth::PointLabelSummary
Truth::label_live_tree(node_t& root,
                       const ICluster::pointer& cluster,
                       const IDepo::vector& depos,
                       const DepoBlobAssocResult& res,
                       const TruthTable& table,
                       const PFResult& pf,
                       const std::map<int, std::shared_ptr<DepoPolyline>>& polylines,
                       const PointLabelParams& params,
                       const depo_position_f& position,
                       Log::logptr_t log,
                       const particle_polylines_t* particles,
                       const node_ends_t* node_ends,
                       const HeatKeys* heat_keys)
{
    PointLabelSummary sum;
    sum.void_roots = pf.void_roots;
    const auto& gr = cluster->graph();

    // Keypoints shared by every point: the starts of the surviving
    // primary nodes (vertex) and of conversion electrons.
    std::vector<DepoPolyline::Keypoint> vtx_kps, conv_kps;
    {
        for (const auto& n : pf.nodes) {
            if (!n.survived || !n.primary) continue;
            auto it = polylines.find(n.trackid);
            if (it == polylines.end() || it->second->nvertices() == 0) continue;
            vtx_kps.push_back({DepoPolyline::KeyKind::start, it->second->start(), 0.0, 0.0});
        }
        const int c_tid = column_index(truth_pf_columns(), "trackid");
        const int c_proc = column_index(truth_pf_columns(), "process");
        const int conv = process_code("conv");
        for (const auto& row : table.pf) {
            if ((int)row[c_proc] != conv) continue;
            auto it = polylines.find((int)row[c_tid]);
            if (it == polylines.end() || it->second->nvertices() == 0) continue;
            conv_kps.push_back({DepoPolyline::KeyKind::start, it->second->start(), 0.0, 0.0});
        }
    }
    // per-track keypoints, lazily
    std::map<int, std::vector<DepoPolyline::Keypoint>> track_kps;
    auto kps_of = [&](int tid) -> const std::vector<DepoPolyline::Keypoint>& {
        auto it = track_kps.find(tid);
        if (it != track_kps.end()) return it->second;
        std::vector<DepoPolyline::Keypoint> kps;
        if (particles) {
            auto qit = particles->find(tid);
            if (qit != particles->end()) {
                for (const auto& q : qit->second) {
                    auto qk = q->keypoints(params.kink_window, params.kink_angle);
                    kps.insert(kps.end(), qk.begin(), qk.end());
                }
            }
        }
        else {
            auto pit = polylines.find(tid);
            if (pit != polylines.end()) kps = pit->second->keypoints(params.kink_window, params.kink_angle);
        }
        return track_kps.emplace(tid, std::move(kps)).first->second;
    };
    sum.npolylines = polylines.size();

    const int c_pdg = column_index(truth_pf_columns(), "pdg");
    const int c_orig = column_index(truth_pf_columns(), "orig_trackid");
    const int c_nurow = column_index(truth_pf_columns(), "nu_row");

    auto clusters = Aux::blob_clusters(gr);
    auto cnodes = root.children();
    if (cnodes.size() != clusters.size()) {
        THROW(ValueError() << errmsg{String::format("truth point label: tree has %d cluster nodes, graph has %d clusters",
                                                     cnodes.size(), clusters.size())});
    }
    size_t ci = 0;
    for (const auto& [cluster_id, vdescs] : clusters) {
        node_t* cnode = cnodes[ci++];
        auto bnodes = cnode->children();
        size_t bi = 0;
        ++sum.nclusters;
        for (const auto& vdesc : vdescs) {
            if (gr[vdesc].code() != 'b') continue;
            const IBlob::pointer iblob = std::get<IBlob::pointer>(gr[vdesc].ptr);
            if (bi >= bnodes.size()) { ++sum.nblobs_skipped; continue; }
            std::string why;
            if (!node_matches_blob(*bnodes[bi], *iblob, params.tick, why)) {
                ++sum.nblobs_skipped;   // sample_live gave this blob no node
                continue;
            }
            node_t* node = bnodes[bi++];
            auto& lpcs = node->value.local_pcs();
            auto pc3d = lpcs.find("3d");
            if (pc3d == lpcs.end()) {
                THROW(ValueError() << errmsg{"truth point label: blob node without a 3d point cloud"});
            }
            const auto xs = pc3d->second.get("x")->elements<double>();
            const auto ys = pc3d->second.get("y")->elements<double>();
            const auto zs = pc3d->second.get("z")->elements<double>();
            const size_t npts = pc3d->second.size_major();

            // blob-level truth
            const auto& e = res.blobs.at(vdesc);
            int tid_top = tid_ghost;
            double q_top = 0.0;
            for (const auto& [tid, qe] : e.by_track) {
                if (qe.first > q_top) { q_top = qe.first; tid_top = tid; }
            }
            const bool ghost = e.q_total <= 0.0 || tid_top == tid_ghost;
            const double frac_top = ghost ? 0.0 : q_top / e.q_total;
            double frac_shared_test = frac_top;
            if (params.shared_by_pf_node && !ghost) {
                // the blob's true charge per PF node; a track without a node is its own key
                std::map<std::pair<int, int>, double> bynode;
                for (const auto& [tid, qe] : e.by_track) {
                    const int pfn = pf.info(tid).pf_node;
                    bynode[pfn >= 0 ? std::make_pair(0, pfn) : std::make_pair(1, tid)] += qe.first;
                }
                double q_node_top = 0.0;
                for (const auto& [key, q] : bynode) q_node_top = std::max(q_node_top, q);
                frac_shared_test = q_node_top / e.q_total;
            }
            const bool shared = !ghost && frac_shared_test < params.shared_frac;
            const double solved = iblob->value();
            const double true_q_frac = solved > 0.0 ? e.q_total / solved : 0.0;
            if (ghost) ++sum.nghost;
            if (shared) ++sum.nshared;

            Columns col;
            col.with_end_distance = params.end_distance;
            col.with_node_heat = params.heat_by_pf_node && heat_keys;
            col.resize(npts);
            for (size_t i = 0; i < npts; ++i) {
                const Point p(xs[i], ys[i], zs[i]);
                int tid = ghost ? tid_ghost : (shared ? tid_shared : tid_top);
                int raw_id = tid;
                if (params.refine && !ghost && !e.depos.empty()) {
                    double best = std::numeric_limits<double>::max();
                    const DepoBlobAssocResult::BlobEntry::DepoShare* bd = nullptr;
                    for (const auto& ds : e.depos) {
                        const double d2 = (position(*ds.depo) - p).magnitude2();
                        if (d2 < best) { best = d2; bd = &ds; }
                    }
                    if (bd) {
                        if (bd->tid != tid_top) ++sum.nrefined_differ;
                        tid = bd->tid;
                        const IDepo* key = bd->depo;
                        if (key->prior()) key = key->prior().get();
                        raw_id = key->id();
                    }
                }
                col.tid[i] = tid;
                double frac = frac_top;
                if (tid >= 0) {
                    auto bt = e.by_track.find(tid);
                    if (bt != e.by_track.end() && e.q_total > 0) frac = bt->second.first / e.q_total;
                }
                col.pid_frac[i] = frac;
                col.true_q_frac[i] = true_q_frac;

                int orig = tid;
                int pdg = 0, cls = (int)Cls::unknown, iid = -1, pfn = -1;
                bool valid = tid >= 0 && !shared;
                if (tid >= 0) {
                    if (!table.sed_ids.empty() && raw_id != tid) orig = table.sed_id(raw_id).second;
                    const int row = table.pf_row(tid);
                    if (row >= 0) {
                        pdg = (int)table.pf[row][c_pdg];
                        if (orig == tid) orig = (int)table.pf[row][c_orig];
                        iid = (int)table.pf[row][c_nurow];
                    }
                    const auto& gi = pf.info(tid);
                    cls = (int)gi.cls;
                    pfn = gi.pf_node;
                    if (pfn >= 0) {
                        if (!gi.instance_valid) valid = false;
                        if (iid < 0) iid = pf.nodes[pfn].nu_row;
                    }
                }
                col.orig_tid[i] = orig;
                col.pdg[i] = pdg;
                col.cls[i] = cls;
                col.iid[i] = iid;
                col.pf_node[i] = pfn;
                col.instance_valid[i] = valid ? 1 : 0;

                // trajectory
                double dtraj = params.d_cap;
                if (tid >= 0) {
                    const DepoPolyline* pl = nullptr;
                    DepoPolyline::Hit hit;
                    if (particles) {
                        // the nearest of the key's particles; a tie keeps the lower particle id
                        auto qit = particles->find(tid);
                        if (qit != particles->end()) {
                            for (const auto& q : qit->second) {
                                auto qhit = q->nearest(p);
                                if (qhit.valid && (!pl || qhit.dist < hit.dist)) {
                                    pl = q.get();
                                    hit = qhit;
                                }
                            }
                        }
                    }
                    else {
                        auto pit = polylines.find(tid);
                        if (pit != polylines.end()) {
                            pl = pit->second.get();
                            hit = pl->nearest(p);
                        }
                    }
                    if (pl && hit.valid) {
                        dtraj = std::min(hit.dist, params.d_cap);
                        col.dir_x[i] = hit.tangent.x();
                        col.dir_y[i] = hit.tangent.y();
                        col.dir_z[i] = hit.tangent.z();
                        col.s_t[i] = pl->length() > 0 ? hit.s / pl->length() : 0.0;
                        col.dedx[i] = hit.dedx * units::cm;
                        if (hit.s <= params.end_window) col.end[i] = 1;
                        else if (pl->length() - hit.s <= params.end_window) col.end[i] = 2;
                    }
                    if (params.end_distance && node_ends) {
                        auto nit = pfn >= 0 ? node_ends->find(pfn) : node_ends->end();
                        if (nit != node_ends->end()) {
                            if (nit->second.has_start) col.d_start[i] = (p - nit->second.start).magnitude();
                            if (nit->second.has_end) col.d_end[i] = (p - nit->second.end).magnitude();
                        }
                    }
                    if (params.end_by_pf_node && node_ends) {
                        // the point's own PF node: start / end within end_radius in 3-D
                        col.end[i] = 0;
                        auto nit = pfn >= 0 ? node_ends->find(pfn) : node_ends->end();
                        if (nit != node_ends->end()) {
                            const auto& ne = nit->second;
                            if (ne.has_start && (p - ne.start).magnitude() <= params.end_radius) col.end[i] = 1;
                            else if (ne.has_end && (p - ne.end).magnitude() <= params.end_radius) col.end[i] = 2;
                        }
                    }
                    if (params.heat_by_pf_node && heat_keys) {
                        // the point's own PF node: its end, the kinks of its root particle
                        col.v_end[i] = 0.0;
                        col.v_kink[i] = 0.0;
                        if (pfn >= 0 && node_ends) {
                            auto nit = node_ends->find(pfn);
                            if (nit != node_ends->end() && nit->second.has_end) {
                                const double d2 = (p - nit->second.end).magnitude2();
                                col.v_end[i] = std::exp(-0.5 * d2 / (params.heat_sigma * params.heat_sigma));
                            }
                        }
                        auto kit = pfn >= 0 ? heat_keys->kinks.find(pfn) : heat_keys->kinks.end();
                        if (kit != heat_keys->kinks.end()) {
                            for (const auto& k : kit->second) {
                                const double d2 = (p - k).magnitude2();
                                col.v_kink[i] = std::max(col.v_kink[i], std::exp(-0.5 * d2 / (params.heat_sigma * params.heat_sigma)));
                            }
                        }
                    }
                    else {
                        const auto& kps = kps_of(tid);
                        col.v_end[i] = DepoPolyline::heat(p, kps, DepoPolyline::KeyKind::end, params.heat_sigma);
                        col.v_kink[i] = DepoPolyline::heat(p, kps, DepoPolyline::KeyKind::kink, params.heat_sigma);
                    }
                }
                col.d_traj[i] = dtraj;
                if (params.heat_by_pf_node && heat_keys) {
                    // true interaction vertices and photon conversion points, every point
                    auto nearest = [&](const std::vector<Point>& keys) {
                        double best = -1.0;
                        for (const auto& k : keys) {
                            const double d = (p - k).magnitude();
                            if (best < 0 || d < best) best = d;
                        }
                        return best;
                    };
                    const double dv = nearest(heat_keys->vertices), dc = nearest(heat_keys->conversions);
                    const double s2 = params.heat_sigma * params.heat_sigma;
                    col.v_vtx[i] = dv < 0 ? 0.0 : std::exp(-0.5 * dv * dv / s2);
                    col.v_conv[i] = dc < 0 ? 0.0 : std::exp(-0.5 * dc * dc / s2);
                    col.d_vtx[i] = dv;
                    col.d_conv[i] = dc;
                }
                else {
                    col.v_vtx[i] = DepoPolyline::heat(p, vtx_kps, DepoPolyline::KeyKind::start, params.heat_sigma);
                    col.v_conv[i] = DepoPolyline::heat(p, conv_kps, DepoPolyline::KeyKind::start, params.heat_sigma);
                }
            }
            lpcs[truth_pc_name] = col.dataset();
            ++sum.nnodes;
            sum.npoints += npts;
        }
        if (bi != bnodes.size()) {
            THROW(ValueError() << errmsg{String::format("truth point label: cluster %d has %d blob nodes but only %d matched blobs",
                                                         cluster_id, bnodes.size(), bi)});
        }
    }
    if (log) {
        log->debug("truth point label: {} clusters, {} nodes, {} points, {} blobs skipped, {} ghost, {} shared, {} refined differ, {} polylines",
                   sum.nclusters, sum.nnodes, sum.npoints, sum.nblobs_skipped, sum.nghost, sum.nshared,
                   sum.nrefined_differ, sum.npolylines);
    }
    return sum;
}
