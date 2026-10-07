// ClusteringBeamDeghost -- replace the blobs of the beam-flash matched clusters
// by the cells of the deghosting cascade (doc pdvd/125).
//
// The production chain (imaging, clustering, Q/L matching) is left as it is.
// The cascade (img CascadeDeghosting) is run on each whole anode as a side
// product and its image arrives here as a second grouping (default name
// "deghost"), built from the cascade's cluster archives by the same
// PointTreeBuilding the clustering job uses.  For each cluster of the beam
// bundle -- chosen with the arithmetic of CheckBeamParticle
// (beam_particle_pick_bundle: the brightest flash among the mains with
// lo <= cluster_t0 < hi) -- the cluster's blobs are destroyed and the cells
// that overlap them (BeamDeghostFunctions.h) are moved in.  Every live cluster
// competes for a cell, so a cell that belongs more to a crossing track outside
// the bundle is not taken.  The cluster keeps its identity: flags, cluster_t0,
// matched flash, default scope.  The T0 correction is applied to the cells
// BEFORE they are moved (a blob inserted without the corrected arrays would
// enter the cluster's corrected scope view short), with the bundle's flash
// time, hence this visitor runs AFTER switch_scope; it runs after
// unmerge_assoc (which needs the per-blob provenance of the old blobs) and
// BEFORE steiner (which builds its point cloud from the cluster's blobs).
//
// The Q/L match is NOT re-evaluated.  Only cells overlapping the production
// footprint are taken; what the cascade image has next to it is counted in the
// log ("outside") and left where it is.
//
// A NEW component: absent from every existing pipeline => no other output
// changes.  Deterministic: clusters are keyed by their position in the
// grouping (tree order), cells are moved in tree order.

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/BeamParticleFunctions.h"
#include "WireCellClus/BeamDeghostFunctions.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"

#include <map>
#include <set>

class ClusteringBeamDeghost;
WIRECELL_FACTORY(ClusteringBeamDeghost, ClusteringBeamDeghost,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::Clus::PR;

static Log::logptr_t logger() {
    static Log::logptr_t l = Log::logger("clus.BeamDeghost");
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

class ClusteringBeamDeghost : public IConfigurable, public Clus::IEnsembleVisitor, private NeedPCTS {
public:
    ClusteringBeamDeghost() {}
    virtual ~ClusteringBeamDeghost() {}

    void configure(const WireCell::Configuration& config) {
        NeedPCTS::configure(config);
        m_grouping = get<std::string>(config, "grouping", m_grouping);
        m_deghost_grouping = get<std::string>(config, "deghost_grouping", m_deghost_grouping);
        m_correction_name = get<std::string>(config, "correction_name", m_correction_name);
        // Same keys and convention as CheckBeamParticle: internal units, the
        // RAW flash time axis, half-open.  low >= high (the default) => off.
        m_beam_window_low = get<double>(config, "beam_window_low", m_beam_window_low);
        m_beam_window_high = get<double>(config, "beam_window_high", m_beam_window_high);
        // A bundle cluster with no overlapping cell is what the cascade removed
        // completely.  true: the cluster is destroyed.  false: it is left as is.
        m_drop_empty = get<bool>(config, "drop_empty", m_drop_empty);
    }

    WireCell::Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping;
        cfg["deghost_grouping"] = m_deghost_grouping;
        cfg["correction_name"] = m_correction_name;
        cfg["beam_window_low"] = m_beam_window_low;
        cfg["beam_window_high"] = m_beam_window_high;
        cfg["drop_empty"] = m_drop_empty;
        return cfg;
    }

    void visit(Ensemble& ensemble) const {
        auto lv = ensemble.with_name(m_grouping);
        auto dv = ensemble.with_name(m_deghost_grouping);
        if (lv.empty() || dv.empty()) {
            logger()->warn("ClusteringBeamDeghost: grouping '{}' or '{}' not found, nothing done",
                           m_grouping, m_deghost_grouping);
            return;
        }
        Grouping& live = *lv.at(0);
        Grouping& dg = *dv.at(0);

        // ---- the beam bundle, as CheckBeamParticle picks it -----------------
        std::vector<BeamBundleRow> rows;
        for (auto* cluster : live.children()) {
            if (!cluster->get_flag(Flags::main_cluster)) continue;
            rows.push_back({cluster->get_cluster_id(), cluster->get_scalar<int>("matched_flash_gid", -1),
                            cluster->get_cluster_t0(), cluster->get_length()});
        }
        std::set<int> gids;   // int-keyed => sorted
        for (const auto& r : rows) {
            if (r.gid >= 0 && r.t0 >= m_beam_window_low && r.t0 < m_beam_window_high) gids.insert(r.gid);
        }
        std::vector<BeamFlashInfo> flashes;
        for (int gid : gids) {
            auto fl = live.flash_by_gid(gid);
            flashes.push_back({gid, fl ? fl.value() : -1.0, static_cast<bool>(fl)});
        }
        const auto pick = beam_particle_pick_bundle(rows, flashes, m_beam_window_low, m_beam_window_high);
        logger()->info("ClusteringBeamDeghost: {} main(s), {} in window [{:.3f}, {:.3f}) us over {} bundle(s); picked gid {} ({})",
                       rows.size(), pick.n_in_window, m_beam_window_low / units::us, m_beam_window_high / units::us,
                       pick.n_gids, pick.gid, pick.why);
        if (pick.gid < 0) {
            logger()->info("ClusteringBeamDeghost: no beam bundle; nothing replaced");
            return;
        }

        // ---- boxes: every live blob competes, every deghost cell is a candidate
        const std::vector<Cluster*> lclusters = live.children();   // tree order; the index is the key
        std::vector<DeghostBox> owners;
        std::set<int> bundle_keys;
        for (size_t k = 0; k < lclusters.size(); ++k) {
            Cluster* c = lclusters[k];
            for (const Blob* b : c->children()) owners.push_back(make_box(*b, (int) k));
            if (c->get_scalar<int>("matched_flash_gid", -1) != pick.gid) continue;
            if (c->npoints() == 0) continue;
            // An out-of-scope part (switch_scope: no corrected point in the
            // active volume) is not evaluated by the PR; leave it alone.
            if (!c->get_scope_filter(c->get_default_scope())) continue;
            bundle_keys.insert((int) k);
        }
        std::vector<Blob*> cells;
        std::vector<int> cell_dgc;      // index of the cell's cluster in the deghost grouping
        std::vector<DeghostBox> cboxes;
        const std::vector<Cluster*> dclusters = dg.children();
        {
            for (size_t k = 0; k < dclusters.size(); ++k) {
                for (Blob* b : dclusters[k]->children()) {
                    cells.push_back(b);
                    cell_dgc.push_back((int) k);
                    cboxes.push_back(make_box(*b, -1));
                }
            }
        }
        const auto owner_of = beam_deghost_assign(cboxes, owners);

        // Cells of the cascade image connected (same cascade cluster) to a cell
        // taken by the bundle but overlapping no production blob at all: what
        // the overlap rule cannot bring in.
        std::set<int> dgc_touched;
        for (size_t i = 0; i < cells.size(); ++i) {
            if (owner_of[i] >= 0 && bundle_keys.count(owner_of[i])) dgc_touched.insert(cell_dgc[i]);
        }
        size_t n_outside = 0; double q_outside = 0;
        for (size_t i = 0; i < cells.size(); ++i) {
            if (owner_of[i] < 0 && dgc_touched.count(cell_dgc[i])) { ++n_outside; q_outside += cells[i]->charge(); }
        }

        // ---- T0 correction of the cells the bundle takes, while they still sit
        // in the cascade clusters.  One bundle = one flash = one cluster_t0, so
        // the cascade cluster borrows it.  A cell with no corrected point in the
        // active volume is not taken (switch_scope would have split it off).
        double bundle_t0 = 0;
        for (int key : bundle_keys) { bundle_t0 = lclusters[key]->get_cluster_t0(); break; }
        std::vector<int> cell_ok(cells.size(), 0);
        {
            std::map<int, std::vector<int>> passed_of;   // cascade cluster index -> per-child flag
            for (int k : dgc_touched) {
                dclusters[k]->set_cluster_t0(bundle_t0);
                passed_of[k] = dclusters[k]->add_corrected_points(m_pcts, m_correction_name);
            }
            std::map<int, size_t> seen;   // cells were collected in child order per cascade cluster
            for (size_t i = 0; i < cells.size(); ++i) {
                const size_t j = seen[cell_dgc[i]]++;
                auto it = passed_of.find(cell_dgc[i]);
                if (it != passed_of.end() && j < it->second.size()) cell_ok[i] = it->second[j];
            }
        }

        // ---- replace, one bundle cluster at a time (ascending key) ------------
        size_t tot_before = 0, tot_after = 0, n_dropped = 0, n_kept_asis = 0;
        double totq_before = 0, totq_after = 0;
        for (int key : bundle_keys) {
            Cluster* c = lclusters[key];
            const int cid = c->get_cluster_id();
            std::vector<size_t> mine;
            std::vector<DeghostBox> mine_boxes;
            size_t n_oos = 0;
            for (size_t i = 0; i < cells.size(); ++i) {
                if (owner_of[i] != key) continue;
                if (!cell_ok[i]) { ++n_oos; continue; }
                mine.push_back(i);
                mine_boxes.push_back(cboxes[i]);
                mine_boxes.back().cluster_id = 0;
            }
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

            if (mine.empty()) {
                if (m_drop_empty) {
                    logger()->info("ClusteringBeamDeghost: gid {} cluster {} blobs {} q {:.0f} -> no cell: cluster dropped",
                                   pick.gid, cid, n_before, q_before);
                    live.destroy_child(c);
                    ++n_dropped;
                }
                else {
                    logger()->info("ClusteringBeamDeghost: gid {} cluster {} blobs {} q {:.0f} -> no cell: kept as is",
                                   pick.gid, cid, n_before, q_before);
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
            size_t n_after = 0;
            for (const Blob* b : c->children()) {
                ++n_after; q_after += b->charge();
                if (!(b->charge() > 0)) ++n_zero_after;
            }
            tot_after += n_after; totq_after += q_after;
            logger()->info("ClusteringBeamDeghost: gid {} cluster {} blobs {} -> {} ({} with zero charge, {} out of scope not taken)"
                           " q {:.0f} -> {:.0f} | production blobs without a cell {} q {:.0f} | points {} L {:.1f} cm",
                           pick.gid, cid, n_before, n_after, n_zero_after, n_oos, q_before, q_after, n_unc, q_unc,
                           c->npoints(), c->get_length() / units::cm);
        }
        logger()->info("ClusteringBeamDeghost: gid {} bundle {} cluster(s), {} dropped, {} kept as is; blobs {} -> {}, q {:.0f} -> {:.0f};"
                       " cascade cells connected to the bundle outside the production footprint: {} q {:.0f}",
                       pick.gid, bundle_keys.size(), n_dropped, n_kept_asis, tot_before, tot_after, totq_before, totq_after,
                       n_outside, q_outside);
    }

private:
    std::string m_grouping{"live"};
    std::string m_deghost_grouping{"deghost"};
    std::string m_correction_name{"T0Correction"};
    double m_beam_window_low{0.0};
    double m_beam_window_high{0.0};
    bool m_drop_empty{true};
};
