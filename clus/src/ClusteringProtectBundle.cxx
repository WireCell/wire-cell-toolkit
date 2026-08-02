#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Util.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/Units.h"
#include "WireCellAux/Logger.h"

#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>

class ClusteringProtectBundle;
WIRECELL_FACTORY(ClusteringProtectBundle, ClusteringProtectBundle,
                 WireCell::INamed, WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

/**
 * PR-stage overclustering protection: the toolkit counterpart of the
 * prototype's second graph-examination round,
 *
 *     prototype (ProtectOverClustering.cxx lines 6-160, called from
 *     wire-cell-prod-nue.cxx line 1322 after Q/L matching):
 *       for every beam-window matched bundle:
 *         for every member cluster:
 *           Examine_graph(ct_point_cloud)      // rebuild graph, ct-test the
 *                                              // inter-fragment bridges
 *           -> split into surviving connected components;
 *         the main cluster's largest component keeps the main cluster id
 *         (lines 57-121), every other component becomes a new cluster in the
 *         same parent bundle (lines 104-137), which NeutrinoID then fits
 *         separately (wire-cell-prod-nue.cxx lines 1345, 1360).
 *
 * Without this stage the PR graph's uncapped MST bridges charge fragments
 * that Clustering_neutrino merged into the nu cluster (photons near the
 * vertex), and the track fitter interpolates fit points across genuinely
 * charge-free space (sbnd_xin doc pr/22 section 8, SBND evt 386948).
 *
 * Placement: after unmerge_assoc (bundle structure restored) and before the
 * steiner stage (separate() moves blob nodes; a pre-split steiner_pc would be
 * stale -- same constraint as ClusteringUnmergeBundle).
 *
 * The split/flag/ident idiom follows ClusteringUnmergeBundle mode
 * "component": connected_blobs() on the configured graph flavor, longest
 * component by get_length() keeps the retained cluster (the prototype keeps
 * the component with the most mcells, ProtectOverClustering.cxx lines 60-70;
 * length is the established toolkit idiom -- divergence recorded in doc
 * pr/23), fragments get flag_associated_cluster and a collision-free ident.
 *
 * SBND divergence, knob-gated OFF by default: the relaxed graph does not join
 * the two halves of a cathode crosser (the ~0.9 cm inactive band at x=0 shows
 * up as a ~4-5 cm charge gap with a ~1 cm transverse offset, sbnd_xin doc
 * pr/20), a geometry uboone did not have -- splitting there would undo the
 * cathode_connect/B0 work.  When cathode_rejoin_xcut > 0, component pairs
 * whose closest points both lie within cathode_rejoin_xcut of cathode_x,
 * within cathode_rejoin_dis in 3D and cathode_rejoin_dyz transversely, are
 * re-united before the split.
 */
class ClusteringProtectBundle : public IConfigurable,
                                public Clus::IEnsembleVisitor,
                                public Aux::Logger,
                                private NeedDV,
                                private NeedPCTS {
public:
    ClusteringProtectBundle()
        : Aux::Logger("ClusteringProtectBundle", "clus")
    {}
    virtual ~ClusteringProtectBundle() {}

    virtual void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config);
        m_grouping_name = get<std::string>(config, "grouping", m_grouping_name);
        m_graph_name = get<std::string>(config, "graph_name", m_graph_name);
        m_pcarray_name = get<std::string>(config, "pcarray_name", m_pcarray_name);
        m_require_in_scope = get<bool>(config, "require_in_scope", m_require_in_scope);
        // Restrict to the beam-coincident bundle(s), the prototype's scope:
        // to_be_checked there is built from beam-window flashes only
        // (wire-cell-prod-nue.cxx lines 1313-1320).  Same gate and same key
        // names as CreateSteinerGraph.  Default false = every scope-passing
        // cluster.
        m_beam_window_only = get<bool>(config, "beam_window_only", m_beam_window_only);
        m_beam_window_low = get<double>(config, "beam_window_low", m_beam_window_low);
        m_beam_window_high = get<double>(config, "beam_window_high", m_beam_window_high);
        // Cathode re-join (SBND).  cathode_rejoin_xcut <= 0 (default) disables
        // the pass entirely => prototype-faithful behavior.  Internal units.
        m_cathode_x = get<double>(config, "cathode_x", m_cathode_x);
        m_cathode_rejoin_xcut = get<double>(config, "cathode_rejoin_xcut", m_cathode_rejoin_xcut);
        m_cathode_rejoin_dyz = get<double>(config, "cathode_rejoin_dyz", m_cathode_rejoin_dyz);
        m_cathode_rejoin_dis = get<double>(config, "cathode_rejoin_dis", m_cathode_rejoin_dis);
    }

    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["detector_volumes"] = "DetectorVolumes";
        cfg["pc_transforms"] = "PCTransformSet";
        cfg["grouping"] = m_grouping_name;
        cfg["graph_name"] = m_graph_name;
        cfg["pcarray_name"] = m_pcarray_name;
        cfg["require_in_scope"] = m_require_in_scope;
        cfg["beam_window_only"] = m_beam_window_only;
        cfg["beam_window_low"] = m_beam_window_low;
        cfg["beam_window_high"] = m_beam_window_high;
        cfg["cathode_x"] = m_cathode_x;
        cfg["cathode_rejoin_xcut"] = m_cathode_rejoin_xcut;
        cfg["cathode_rejoin_dyz"] = m_cathode_rejoin_dyz;
        cfg["cathode_rejoin_dis"] = m_cathode_rejoin_dis;
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const {
        auto groupings = ensemble.with_name(m_grouping_name);
        if (groupings.empty()) {
            log->warn("no '{}' grouping found, skipping", m_grouping_name);
            return;
        }
        auto& grouping = *groupings.at(0);

        // Snapshot: separate() adds children to the grouping.
        const std::vector<Cluster*> clusters = grouping.children();

        std::set<int> taken;
        for (const auto* c : clusters) taken.insert(c->ident());

        // Beam gate (CreateSteinerGraph idiom): gids of in-window mains; a
        // bundle member is the main itself or a companion sharing its gid.
        std::unordered_set<int> beam_gids;
        const bool gate = m_beam_window_only && m_beam_window_low < m_beam_window_high;
        if (gate) {
            for (const auto* cluster : clusters) {
                if (!cluster->get_flag(Flags::main_cluster)) continue;
                const double t0 = cluster->get_cluster_t0();
                if (t0 < m_beam_window_low || t0 >= m_beam_window_high) continue;
                const int gid = cluster->get_scalar<int>("matched_flash_gid", -1);
                if (gid >= 0) beam_gids.insert(gid);
            }
        }

        size_t n_split = 0, n_parts = 0, n_rejoin = 0;
        for (Cluster* cluster : clusters) {
            if (m_require_in_scope
                && !cluster->get_scope_filter(cluster->get_default_scope())) continue;
            if (gate) {
                const int gid = cluster->get_scalar<int>("matched_flash_gid", -1);
                if (gid < 0 || beam_gids.count(gid) == 0) continue;
            }
            const size_t before = grouping.nchildren();
            if (process_cluster(grouping, cluster, taken, n_rejoin)) {
                ++n_split;
                n_parts += grouping.nchildren() - before;
            }
        }
        log->debug("split {} bundle cluster(s) into {} extra cluster(s) "
                   "({} cathode re-join(s), graph '{}')",
                   n_split, n_parts, n_rejoin, m_graph_name);
    }

private:
    std::string m_grouping_name{"live"};
    std::string m_graph_name{"relaxed"};
    std::string m_pcarray_name{"perblob"};
    bool m_require_in_scope{true};
    bool m_beam_window_only{false};
    double m_beam_window_low{0};
    double m_beam_window_high{0};
    double m_cathode_x{0};
    double m_cathode_rejoin_xcut{0};       // <= 0: pass disabled (prototype behavior)
    double m_cathode_rejoin_dyz{4 * units::cm};
    double m_cathode_rejoin_dis{8 * units::cm};

    /// Union components that meet across the cathode band.  `b2g` is edited
    /// in place; returns the number of pair unions applied.
    size_t cathode_rejoin(const Cluster* cluster, std::vector<int>& b2g) const {
        std::set<int> ids;
        for (int id : b2g) if (id >= 0) ids.insert(id);
        if (ids.size() < 2) return 0;

        // One point cloud per component, points mapped through their blob
        // (same construction as connect_graph_relaxed.cxx lines 71-82; the
        // point coordinates are in the default scope, i.e. the T0-corrected
        // frame that cathode_x is defined in).
        std::map<int, Simple3DPointCloud> clouds;
        const auto& points = cluster->points();
        const size_t npts = points[0].size();
        for (size_t i = 0; i < npts; ++i) {
            const int comp = b2g[cluster->kd3d().major_index(i)];
            if (comp < 0) continue;
            clouds[comp].add({points[0][i], points[1][i], points[2][i]});
        }

        // Union-find over component ids.
        std::map<int, int> parent;
        for (int id : ids) parent[id] = id;
        std::function<int(int)> find = [&](int a) {
            while (parent[a] != a) { parent[a] = parent[parent[a]]; a = parent[a]; }
            return a;
        };

        size_t n_union = 0;
        for (auto j = clouds.begin(); j != clouds.end(); ++j) {
            for (auto k = std::next(j); k != clouds.end(); ++k) {
                const auto [i1, i2, dis] = j->second.get_closest_points(k->second);
                if (dis >= m_cathode_rejoin_dis) continue;
                const auto p1 = j->second.point(i1);
                const auto p2 = k->second.point(i2);
                if (std::abs(p1.x() - m_cathode_x) >= m_cathode_rejoin_xcut) continue;
                if (std::abs(p2.x() - m_cathode_x) >= m_cathode_rejoin_xcut) continue;
                const double dyz = std::hypot(p1.y() - p2.y(), p1.z() - p2.z());
                if (dyz >= m_cathode_rejoin_dyz) continue;
                const int ra = find(j->first), rb = find(k->first);
                if (ra == rb) continue;
                parent[std::max(ra, rb)] = std::min(ra, rb);   // keep lowest id
                ++n_union;
                SPDLOG_LOGGER_DEBUG(log,
                    "cluster {}: cathode re-join comp {}+{} (gap {:.2f} cm, dyz {:.2f} cm, "
                    "x {:.2f}/{:.2f} cm)",
                    cluster->ident(), j->first, k->first, dis / units::cm,
                    dyz / units::cm, p1.x() / units::cm, p2.x() / units::cm);
            }
        }
        if (n_union) {
            for (int& id : b2g) if (id >= 0) id = find(id);
        }
        return n_union;
    }

    /// Split one bundle cluster at its surviving component boundaries.
    /// Returns true if a split happened.
    bool process_cluster(Grouping& grouping, Cluster*& cluster,
                         std::set<int>& taken, size_t& n_rejoin) const {
        const size_t nb = cluster->nchildren();
        if (nb < 2) return false;

        // prototype: Examine_graph(ct_point_cloud), ProtectOverClustering.cxx:56
        auto b2g = cluster->connected_blobs(m_dv, m_pcts, m_graph_name);
        if (b2g.size() != nb) return false;

        if (m_cathode_rejoin_xcut > 0) n_rejoin += cathode_rejoin(cluster, b2g);

        std::set<int> ids;
        for (int id : b2g) if (id >= 0) ids.insert(id);
        if (ids.size() < 2) return false;

        // prototype keeps the component with the most mcells
        // (ProtectOverClustering.cxx:60-70); here longest by get_length(),
        // ascending iteration so ties keep the lowest id.
        double max_length = -1;
        int longest = -1;
        for (int id : ids) {
            const double len = get_length(cluster, b2g, id);
            if (len > max_length) { max_length = len; longest = id; }
        }

        std::map<int, int> id2grp;
        int g = 0;
        for (int id : ids) if (id != longest) id2grp[id] = g++;

        std::vector<int> groups(nb, -1);
        for (size_t i = 0; i < nb; ++i) {
            if (b2g[i] < 0 || b2g[i] == longest) continue;   // unassigned stays put
            groups[i] = id2grp.at(b2g[i]);
        }

        // Snapshot the per-blob PC before the blobs move away (same landmine
        // as ClusteringUnmergeBundle: separate(remove=false) leaves the
        // retained cluster holding a full-length array over a shorter blob
        // list).
        PointCloud::Dataset perblob;
        bool have_perblob = false;
        {
            auto& lpcs = cluster->value().local_pcs();
            auto it = lpcs.find(m_pcarray_name);
            if (it != lpcs.end()) {
                if (it->second.size_major() == nb) {
                    perblob = it->second;      // copy
                    have_perblob = true;
                }
                else {
                    log->warn("cluster {}: '{}' has {} rows for {} blobs, dropping it",
                              cluster->ident(), m_pcarray_name,
                              it->second.size_major(), nb);
                }
            }
        }

        const int main_ident = cluster->ident();
        const bool was_main = cluster->get_flag(Flags::main_cluster);
        auto splits = grouping.separate(cluster, groups);   // remove=false: largest stays put

        auto carve = [&](Cluster* part, int want) {
            auto& lpcs = part->value().local_pcs();
            lpcs.erase(m_pcarray_name);
            if (!have_perblob) return;
            std::vector<size_t> rows;
            rows.reserve(nb);
            for (size_t i = 0; i < nb; ++i) if (groups[i] == want) rows.push_back(i);
            if (rows.size() != part->nchildren()) {
                log->warn("cluster {}: carved {} rows for {} blobs, dropping '{}'",
                          part->ident(), rows.size(), part->nchildren(), m_pcarray_name);
                return;
            }
            lpcs[m_pcarray_name] = perblob.subset(rows);
        };

        carve(cluster, -1);
        if (cluster->value().local_pcs().erase("steiner_pc")) {
            log->warn("cluster {}: erased a steiner_pc built before the split -- "
                      "run this visitor BEFORE the steiner stage", main_ident);
        }

        // prototype: the largest component keeps the main cluster id and role
        // (ProtectOverClustering.cxx:101-110); the retained cluster already
        // carries the flags and cluster_scalar (matched_flash_gid, cluster_t0)
        // -- fragments get them via separate()->from(), so a fragment of the
        // bundle main must have flag_main_cluster cleared or every fragment
        // would look like a bundle main to STM/TGM/FC.
        size_t nb_parts = 0;
        for (const auto& [gid, part] : splits) nb_parts += part->nchildren();
        log->debug("cluster {}{}: {} blobs -> retained {} + {} fragment(s) holding {}",
                   main_ident, was_main ? " (main)" : "", nb, cluster->nchildren(),
                   splits.size(), nb_parts);
        if (cluster->nchildren() + nb_parts != nb) {
            log->error("cluster {}: BLOB LOSS in separate(): {} + {} != {}",
                       main_ident, cluster->nchildren(), nb_parts, nb);
        }
        int sub_id = 1;
        for (auto& [gid, part] : splits) {
            part->set_flag(Flags::main_cluster, 0);
            part->set_flag(Flags::associated_cluster);
            // prototype: fresh ids via get_next_cluster_id
            // (ProtectOverClustering.cxx:113-124); here the collision-avoiding
            // ClusteringUnmergeBundle scheme.
            part->set_ident(alloc_ident(taken, main_ident * 100 + sub_id));
            ++sub_id;
            carve(part, gid);
        }
        return true;
    }

    /// `want` if free, else the lowest free integer above every taken ident.
    static int alloc_ident(std::set<int>& taken, int want) {
        if (!taken.count(want)) { taken.insert(want); return want; }
        int cand = taken.empty() ? 0 : *taken.rbegin() + 1;
        while (taken.count(cand)) ++cand;
        taken.insert(cand);
        return cand;
    }
};
