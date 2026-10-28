// ClusteringImageDeghost -- replace the blobs of EVERY live cluster by the
// cells of the deghosting cascade (pr_ml doc 28, "combination 3": traditional
// imaging, clustering and Q/L; the cascade image enters only here, before the
// PR).  Fork by duplication of ClusteringBeamDeghost (doc pdvd/125, which acts
// on the beam bundle only and stays untouched); differences:
//   * the clusters replaced are every live cluster with points and an
//     in-scope default view (select "all"), or the main clusters ("main"),
//     instead of the beam-flash bundle;
//   * the T0 correction is applied per cascade cluster with the t0 of the live
//     clusters taking its cells (image_deghost_cluster_t0): one bundle had one
//     t0, every cluster may not -- a cell whose owner's t0 differs from the
//     t0 its cascade cluster borrows is not taken (counted in the log),
//     unless t0_split (pr_ml doc 32 sec B.4) moves those cells to a cascade
//     cluster of their own, corrected with their owner's t0;
//   * per-cluster lines go to debug, one summary line to info.
// As in the beam stage: the cascade image arrives as a second grouping
// (default "deghost", built by PointTreeBuilding from the cascade's cluster
// archives); every live blob competes for a cell (beam_deghost_assign: the
// owner overlapping a cell most takes it), only cells overlapping the
// production footprint are taken, the cell's T0 correction is applied while
// it still sits in its cascade cluster, the cluster keeps its identity
// (flags, cluster_t0, matched flash), its row-parallel "perblob" PC is erased,
// and the Q/L match is not re-evaluated.  Runs after switch_scope and
// unmerge_* and before steiner.
//
// A NEW component: absent from every existing pipeline => no other output
// changes.  Deterministic: clusters keyed by tree position, cells moved in
// tree order.

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/BeamDeghostFunctions.h"
#include "WireCellClus/ImageDeghostFunctions.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Persist.h"

#include <map>
#include <set>

class ClusteringImageDeghost;
WIRECELL_FACTORY(ClusteringImageDeghost, ClusteringImageDeghost,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::Clus::PR;

static Log::logptr_t logger() {
    static Log::logptr_t l = Log::logger("clus.ImageDeghost");
    return l;
}

static DeghostBox make_box(const Blob& b, int cluster_key)
{
    DeghostBox x;
    x.wpid = b.wpid().ident();
    x.t0 = b.slice_index_min(); x.t1 = b.slice_index_max();
    x.u0 = b.u_wire_index_min(); x.u1 = b.u_wire_index_max();
    x.v0 = b.v_wire_index_min(); x.v1 = b.v_wire_index_max();
    x.w0 = b.w_wire_index_min(); x.w1 = b.w_wire_index_max();
    x.cluster_id = cluster_key;
    return x;
}

class ClusteringImageDeghost : public IConfigurable, public Clus::IEnsembleVisitor, private NeedPCTS {
public:
    ClusteringImageDeghost() {}
    virtual ~ClusteringImageDeghost() {}

    void configure(const WireCell::Configuration& config) {
        NeedPCTS::configure(config);
        m_grouping = get<std::string>(config, "grouping", m_grouping);
        m_deghost_grouping = get<std::string>(config, "deghost_grouping", m_deghost_grouping);
        m_correction_name = get<std::string>(config, "correction_name", m_correction_name);
        // "all": every live cluster with points and an in-scope default view; "main": main clusters only.
        m_select = get<std::string>(config, "select", m_select);
        if (m_select != "all" && m_select != "main") {
            raise<ValueError>("ClusteringImageDeghost: select must be 'all' or 'main', got '%s'", m_select.c_str());
        }
        // A replaced cluster with no overlapping cell is what the cascade removed completely.
        // true: the cluster is destroyed.  false: it is left as is.
        m_drop_empty = get<bool>(config, "drop_empty", m_drop_empty);
        // cluster_t0 values closer than this are one t0 (internal time units).
        m_t0_tolerance = get<double>(config, "t0_tolerance", m_t0_tolerance);
        // As ClusteringBeamDeghost (pr_ml doc 24 sec 6): keep the deghost grouping serializable.
        m_fill_cluster_t0 = get<bool>(config, "fill_cluster_t0", m_fill_cluster_t0);
        m_strip_t0cor_leftover = get<bool>(config, "strip_t0cor_leftover", m_strip_t0cor_leftover);
        // The live blobs may carry "3d" arrays persisted by an earlier job that the cells do not have (SBND
        // simulation: y_cor / z_cor written by the clustering job, not read by the PR's sim scope); TensorDM
        // then throws or drops them depending on node order.  true: after the replacement, erase from every
        // live blob each "3d" array that not all live blobs carry (logged).  Default false.
        m_harmonize_3d_keys = get<bool>(config, "harmonize_3d_keys", m_harmonize_3d_keys);
        // pr_ml doc 32 sec B.3: a diagnostic dump of the hand-over (every live blob, every cell with its owner, the
        // t0 decision of its cascade cluster and whether it was taken) to <dump_dir>/image_deghost.json.  Empty
        // (default) = no file; the dump never changes the output.
        m_dump_dir = get<std::string>(config, "dump_dir", m_dump_dir);
        // pr_ml doc 32 sec B.4: a cascade cluster joining tracks of different flash t0 borrows one t0 and leaves the
        // cells of the other owners behind (29107_3 cluster 118: 383 cells, 4.9 M e).  true: those cells move to a new
        // cascade cluster per owner t0, are corrected with that t0 and are taken.  Default false (one t0 per cascade
        // cluster, the cells left).
        m_t0_split = get<bool>(config, "t0_split", m_t0_split);
    }

    WireCell::Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping;
        cfg["deghost_grouping"] = m_deghost_grouping;
        cfg["correction_name"] = m_correction_name;
        cfg["select"] = m_select;
        cfg["drop_empty"] = m_drop_empty;
        cfg["t0_tolerance"] = m_t0_tolerance;
        cfg["fill_cluster_t0"] = m_fill_cluster_t0;
        cfg["strip_t0cor_leftover"] = m_strip_t0cor_leftover;
        cfg["harmonize_3d_keys"] = m_harmonize_3d_keys;
        cfg["dump_dir"] = m_dump_dir;
        cfg["t0_split"] = m_t0_split;
        return cfg;
    }

    void visit(Ensemble& ensemble) const {
        auto lv = ensemble.with_name(m_grouping);
        auto dv = ensemble.with_name(m_deghost_grouping);
        if (lv.empty() || dv.empty()) {
            logger()->warn("ClusteringImageDeghost: grouping '{}' or '{}' not found, nothing done",
                           m_grouping, m_deghost_grouping);
            return;
        }
        Grouping& live = *lv.at(0);
        Grouping& dg = *dv.at(0);

        // ---- the clusters to replace; every live blob competes for a cell -----
        const std::vector<Cluster*> lclusters = live.children();   // tree order; the index is the key
        std::vector<DeghostBox> owners;
        std::vector<char> selected(lclusters.size(), 0);
        std::vector<double> owner_t0(lclusters.size(), 0);
        std::vector<int> keys;
        for (size_t k = 0; k < lclusters.size(); ++k) {
            Cluster* c = lclusters[k];
            owner_t0[k] = c->get_cluster_t0();
            for (const Blob* b : c->children()) owners.push_back(make_box(*b, (int) k));
            if (m_select == "main" && !c->get_flag(Flags::main_cluster)) continue;
            if (c->npoints() == 0) continue;
            // An out-of-scope part (switch_scope: no corrected point in the
            // active volume) is not evaluated by the PR; leave it alone.
            if (!c->get_scope_filter(c->get_default_scope())) continue;
            selected[k] = 1;
            keys.push_back((int) k);
        }
        std::vector<Blob*> cells;
        std::vector<int> cell_dgc;      // index of the cell's cluster in the deghost grouping
        std::vector<DeghostBox> cboxes;
        std::vector<Cluster*> dclusters = dg.children();   // t0_split appends
        for (size_t k = 0; k < dclusters.size(); ++k) {
            for (Blob* b : dclusters[k]->children()) {
                cells.push_back(b);
                cell_dgc.push_back((int) k);
                cboxes.push_back(make_box(*b, -1));
            }
        }
        const auto owner_of = beam_deghost_assign(cboxes, owners);
        std::vector<char> take;
        auto dgt0 = image_deghost_cluster_t0(cell_dgc, owner_of, selected, owner_t0, (int) dclusters.size(),
                                                   m_t0_tolerance, take);
        size_t n_t0_conflict = 0, n_t0_left = 0;
        for (const auto& d : dgt0) n_t0_conflict += d.conflict ? 1 : 0;
        for (size_t i = 0; i < cells.size(); ++i) {
            if (owner_of[i] >= 0 && selected[owner_of[i]] && !take[i]) ++n_t0_left;
        }

        // Cells of a touched cascade cluster that overlap no production blob at
        // all: what the overlap rule cannot bring in.
        size_t n_outside = 0; double q_outside = 0;
        for (size_t i = 0; i < cells.size(); ++i) {
            if (owner_of[i] < 0 && dgt0[cell_dgc[i]].touched) { ++n_outside; q_outside += cells[i]->charge(); }
        }

        // ---- t0_split: the cells a t0 conflict leaves behind move, one new cascade cluster per owner t0 (appended
        // to the deghost grouping, cells in ascending index = child order), and are taken.
        size_t n_split = 0, n_split_cells = 0;
        Json::Value jsplit(Json::arrayValue);
        if (m_t0_split) {
            const auto splits = image_deghost_t0_split(cell_dgc, owner_of, selected, owner_t0, dgt0, m_t0_tolerance);
            for (const auto& sp : splits) {
                Cluster& nc = dg.make_child();
                const int nd = (int) dclusters.size();
                for (size_t i : sp.cells) {
                    auto node = cells[i]->cluster()->remove_child(*cells[i]);
                    nc.node()->insert(std::move(node));
                    cell_dgc[i] = nd;
                    take[i] = 1;
                }
                dclusters[sp.dgc]->invalidate_cache();
                nc.invalidate_cache();
                dclusters.push_back(&nc);
                ImageDeghostT0 x;
                x.touched = true;
                x.t0 = sp.t0;
                dgt0.push_back(x);
                ++n_split;
                n_split_cells += sp.cells.size();
                Json::Value r(Json::arrayValue);
                r.append(sp.dgc);
                r.append(nd);
                r.append(sp.t0);
                r.append((Json::UInt64) sp.cells.size());
                jsplit.append(r);
            }
        }

        // ---- T0 correction of the taken cells, per cascade cluster, while they
        // still sit there (a blob inserted without the corrected arrays enters
        // the cluster's corrected scope view short, doc pdvd/125).
        std::vector<int> cell_ok(cells.size(), 0);
        {
            std::map<int, std::vector<int>> passed_of;   // cascade cluster index -> per-child flag
            for (size_t k = 0; k < dclusters.size(); ++k) {
                if (!dgt0[k].touched) continue;
                dclusters[k]->set_cluster_t0(dgt0[k].t0);
                passed_of[(int) k] = dclusters[k]->add_corrected_points(m_pcts, m_correction_name);
            }
            std::map<int, size_t> seen;   // cells were collected in child order per cascade cluster
            for (size_t i = 0; i < cells.size(); ++i) {
                const size_t j = seen[cell_dgc[i]]++;
                auto it = passed_of.find(cell_dgc[i]);
                if (it != passed_of.end() && j < it->second.size()) cell_ok[i] = it->second[j];
            }
        }

        // pr_ml doc 32 sec B.3: the hand-over as it stands before any blob moves
        Json::Value jdump;
        const bool dumping = !m_dump_dir.empty();
        auto jbox = [](const DeghostBox& x) {
            Json::Value r(Json::arrayValue);
            for (int v : {x.wpid, x.t0, x.t1, x.u0, x.u1, x.v0, x.v1, x.w0, x.w1}) r.append(v);
            return r;
        };
        if (dumping) {
            for (size_t k = 0; k < lclusters.size(); ++k) {
                Cluster* c = lclusters[k];
                Json::Value jc;
                jc["key"] = (int) k;
                jc["cluster_id"] = c->get_cluster_id();
                jc["t0"] = owner_t0[k];
                jc["selected"] = (int) selected[k];
                jc["npoints"] = (Json::UInt64) c->npoints();
                Json::Value jb(Json::arrayValue);
                for (const Blob* b : c->children()) {
                    Json::Value r = jbox(make_box(*b, (int) k));
                    r.append(b->charge());
                    r.append(b->center_x());
                    r.append(b->center_y());
                    r.append(b->center_z());
                    jb.append(r);
                }
                jc["blobs"] = jb;
                jdump["clusters"].append(jc);
            }
            for (size_t i = 0; i < cells.size(); ++i) {
                Json::Value r = jbox(cboxes[i]);
                r.append(cell_dgc[i]);
                r.append(owner_of[i]);
                r.append((int) take[i]);
                r.append(cell_ok[i]);
                r.append(cells[i]->charge());
                r.append(cells[i]->center_x());
                r.append(cells[i]->center_y());
                r.append(cells[i]->center_z());
                jdump["cells"].append(r);
            }
            for (size_t k = 0; k < dgt0.size(); ++k) {
                Json::Value r(Json::arrayValue);
                r.append((int) dgt0[k].touched);
                r.append(dgt0[k].t0);
                r.append((int) dgt0[k].conflict);
                jdump["cascade_clusters"].append(r);
            }
            jdump["cell_columns"] = "wpid t0 t1 u0 u1 v0 v1 w0 w1 cascade_cluster owner_key take in_scope charge x y z";
            jdump["blob_columns"] = "wpid t0 t1 u0 u1 v0 v1 w0 w1 charge x y z";
            jdump["cascade_cluster_columns"] = "touched t0 conflict";
            if (m_t0_split) {
                jdump["t0_splits"] = jsplit;
                jdump["t0_split_columns"] = "from_cascade_cluster new_cascade_cluster t0 ncells";
            }
        }

        // ---- replace, one cluster at a time (ascending key) -------------------
        size_t tot_before = 0, tot_after = 0, n_dropped = 0, n_kept_asis = 0, n_replaced = 0, tot_oos = 0, tot_unc = 0;
        double totq_before = 0, totq_after = 0, totq_unc = 0;
        for (int key : keys) {
            Cluster* c = lclusters[key];
            const int cid = c->get_cluster_id();
            std::vector<size_t> mine;
            std::vector<DeghostBox> mine_boxes;
            size_t n_oos = 0;
            for (size_t i = 0; i < cells.size(); ++i) {
                if (owner_of[i] != key || !take[i]) continue;
                if (!cell_ok[i]) { ++n_oos; continue; }
                mine.push_back(i);
                mine_boxes.push_back(cboxes[i]);
                mine_boxes.back().cluster_id = 0;
            }
            tot_oos += n_oos;
            // production blobs of this cluster with no cell of its own on top
            std::vector<Blob*> old = c->children();
            std::vector<DeghostBox> old_boxes;
            for (const Blob* b : old) old_boxes.push_back(make_box(*b, -1));
            const auto covered = beam_deghost_assign(old_boxes, mine_boxes);
            size_t n_unc = 0, n_zero_after = 0; double q_before = 0, q_unc = 0, q_after = 0;
            for (size_t i = 0; i < old.size(); ++i) {
                q_before += old[i]->charge();
                if (covered[i] < 0) { ++n_unc; q_unc += old[i]->charge(); }
            }
            const size_t n_before = old.size();
            tot_before += n_before; totq_before += q_before;
            tot_unc += n_unc; totq_unc += q_unc;

            if (mine.empty()) {
                if (m_drop_empty) {
                    SPDLOG_LOGGER_DEBUG(logger(), "ClusteringImageDeghost: cluster {} blobs {} q {:.0f} -> no cell: dropped",
                                        cid, n_before, q_before);
                    live.destroy_child(c);
                    ++n_dropped;
                }
                else {
                    SPDLOG_LOGGER_DEBUG(logger(), "ClusteringImageDeghost: cluster {} blobs {} q {:.0f} -> no cell: kept as is",
                                        cid, n_before, q_before);
                    tot_after += n_before; totq_after += q_before;
                    ++n_kept_asis;
                }
                continue;
            }

            for (Blob* b : old) c->destroy_child(b);
            // The per-blob provenance arrays are row-parallel to the old blobs.
            c->local_pcs().erase("perblob");
            for (size_t i : mine) {
                Blob* cell = cells[i];
                auto node = cell->cluster()->remove_child(*cell);
                c->node()->insert(std::move(node));
            }
            c->invalidate_cache();
            ++n_replaced;
            size_t n_after = 0;
            for (const Blob* b : c->children()) {
                ++n_after; q_after += b->charge();
                if (!(b->charge() > 0)) ++n_zero_after;
            }
            tot_after += n_after; totq_after += q_after;
            SPDLOG_LOGGER_DEBUG(logger(), "ClusteringImageDeghost: cluster {} blobs {} -> {} ({} with zero charge, {} out of scope not taken)"
                                " q {:.0f} -> {:.0f} | production blobs without a cell {} q {:.0f} | points {} L {:.1f} cm",
                                cid, n_before, n_after, n_zero_after, n_oos, q_before, q_after, n_unc, q_unc,
                                c->npoints(), c->get_length() / units::cm);
        }
        logger()->info("ClusteringImageDeghost: select {} {} of {} cluster(s): {} replaced, {} dropped, {} kept as is; blobs {} -> {},"
                       " q {:.0f} -> {:.0f}; production blobs without a cell {} q {:.0f}; cells out of scope {};"
                       " cascade clusters with a t0 conflict {} (cells left {}); cells of touched cascade clusters outside"
                       " the production footprint {} q {:.0f}",
                       m_select, keys.size(), lclusters.size(), n_replaced, n_dropped, n_kept_asis, tot_before, tot_after,
                       totq_before, totq_after, tot_unc, totq_unc, tot_oos, n_t0_conflict, n_t0_left, n_outside, q_outside);
        if (m_t0_split) {
            logger()->info("ClusteringImageDeghost: t0_split: {} cascade cluster(s) of the {} with a t0 conflict split off,"
                           " {} cell(s) taken with their owner's t0", n_split, n_t0_conflict, n_split_cells);
        }
        if (dumping) Persist::dump(m_dump_dir + "/image_deghost.json", jdump);
        if (m_strip_t0cor_leftover) {
            const auto& pct = m_pcts->pc_transform(m_correction_name);
            const auto names = pct->stored_array_names();
            size_t nerased = 0, ntouched = 0;
            for (size_t k = 0; k < dclusters.size(); ++k) {
                if (!dgt0[k].touched) continue;
                ++ntouched;
                nerased += erase_blob_arrays(*dclusters[k], names);
            }
            logger()->info("ClusteringImageDeghost: strip_t0cor_leftover: {} array(s) erased from cells left in {} touched cascade cluster(s)",
                           nerased, ntouched);
        }
        if (m_harmonize_3d_keys) {
            std::vector<std::string> partial;
            const size_t nerased = harmonize_blob_arrays(live, "3d", partial);
            std::string names;
            for (const auto& name : partial) names += (names.empty() ? "" : ",") + name;
            logger()->info("ClusteringImageDeghost: harmonize_3d_keys: {} partial key(s) [{}], {} array(s) erased",
                           partial.size(), names, nerased);
        }
        if (m_fill_cluster_t0) {
            const size_t nfill = fill_cluster_t0_before_serialize(dg);
            logger()->info("ClusteringImageDeghost: fill_cluster_t0: {} cascade cluster(s) filled", nfill);
        }
    }

private:
    std::string m_grouping{"live"};
    std::string m_deghost_grouping{"deghost"};
    std::string m_correction_name{"T0Correction"};
    std::string m_select{"all"};
    bool m_drop_empty{true};
    double m_t0_tolerance{1.0 * units::ns};
    bool m_fill_cluster_t0{false};
    bool m_strip_t0cor_leftover{false};
    bool m_harmonize_3d_keys{false};
    std::string m_dump_dir{""};
    bool m_t0_split{false};
};
