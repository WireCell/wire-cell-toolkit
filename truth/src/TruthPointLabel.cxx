#include "WireCellTruth/TruthPointLabel.h"
#include "WireCellTruth/Schema.h"

#include "WireCellAux/ClusterHelpers.h"
#include "WireCellIface/WirePlaneId.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/String.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace WireCell;
using WireCell::PointCloud::Array;
using WireCell::PointCloud::Dataset;
using node_t = WireCell::PointCloud::Tree::Points::node_t;

std::map<int, std::shared_ptr<Truth::DepoPolyline>>
Truth::build_polylines(const IDepo::vector& depos, const DepoBlobAssocParams& assoc,
                       const depo_position_f& position)
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
        auto& pl = out[tid];
        if (!pl) pl = std::make_shared<DepoPolyline>();
        pl->add(position(*d), key->time(), energy);
    }
    for (auto& [tid, pl] : out) pl->build();
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
        void resize(size_t n)
        {
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
                       Log::logptr_t log)
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
        auto pit = polylines.find(tid);
        if (pit != polylines.end()) kps = pit->second->keypoints(params.kink_window, params.kink_angle);
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
            const bool shared = !ghost && frac_top < params.shared_frac;
            const double solved = iblob->value();
            const double true_q_frac = solved > 0.0 ? e.q_total / solved : 0.0;
            if (ghost) ++sum.nghost;
            if (shared) ++sum.nshared;

            Columns col;
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
                    auto pit = polylines.find(tid);
                    if (pit != polylines.end()) {
                        const auto& pl = *pit->second;
                        auto hit = pl.nearest(p);
                        if (hit.valid) {
                            dtraj = std::min(hit.dist, params.d_cap);
                            col.dir_x[i] = hit.tangent.x();
                            col.dir_y[i] = hit.tangent.y();
                            col.dir_z[i] = hit.tangent.z();
                            col.s_t[i] = pl.length() > 0 ? hit.s / pl.length() : 0.0;
                            col.dedx[i] = hit.dedx * units::cm;
                            if (hit.s <= params.end_window) col.end[i] = 1;
                            else if (pl.length() - hit.s <= params.end_window) col.end[i] = 2;
                        }
                    }
                    const auto& kps = kps_of(tid);
                    col.v_end[i] = DepoPolyline::heat(p, kps, DepoPolyline::KeyKind::end, params.heat_sigma);
                    col.v_kink[i] = DepoPolyline::heat(p, kps, DepoPolyline::KeyKind::kink, params.heat_sigma);
                }
                col.d_traj[i] = dtraj;
                col.v_vtx[i] = DepoPolyline::heat(p, vtx_kps, DepoPolyline::KeyKind::start, params.heat_sigma);
                col.v_conv[i] = DepoPolyline::heat(p, conv_kps, DepoPolyline::KeyKind::start, params.heat_sigma);
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
