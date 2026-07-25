#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellAux/Logger.h"

#include <map>
#include <set>
#include <vector>

class ClusteringUnmergeBundle;
WIRECELL_FACTORY(ClusteringUnmergeBundle, ClusteringUnmergeBundle,
                 WireCell::INamed, WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

/**
 * Restore the uBooNE "main cluster + associated clusters" data product on a
 * post-Q/L-matching point-cloud tree.
 *
 * Motivation.  The prototype pattern-recognition chain is handed ONE main
 * cluster plus a vector of additional clusters per matched flash bundle:
 *
 *     prototype (wire-cell-prod-stm.cxx lines 815-830):
 *       temp_clusters = map_parentid_clusters[tpc_cluster_id];
 *       main_cluster  = the member whose get_cluster_id() == tpc_cluster_id;
 *       additional_clusters = all other members;
 *       main_cluster->create_steiner_graph(...);   // main only
 *       fid->check_stm(main_cluster, additional_clusters, ...);
 *
 * The main cluster there is a single connected PR3DCluster and every other
 * bundle member stays a SEPARATE cluster object.  check_stm() fits the main
 * as one track and only *counts* the companions (check_other_clusters()).
 *
 * On SBND the Q/L stage runs clustering_examine_bundles(use_flash_t0=true),
 * which MERGES every flash-time-coincident member of a bundle into one
 * cluster carrying flag_main_cluster.  The taggers therefore see a single
 * cluster made of many spatially detached pieces:
 *   - TaggerCheckSTM fits that composite as one track (garbage trajectories,
 *     e.g. showcase-stmfit-mc-evt18 track_com_18.root cluster 80),
 *   - check_other_clusters() has no companions left to count, so it is inert,
 *   - TaggerCheckNeutrino's "companions sharing matched_flash_gid" is empty.
 * TaggerCheckTGM was adapted internally instead (main_component_pairs, doc
 * 36/38), but that fix does not generalise to the fitting code.
 *
 * This visitor puts the structure back where it belongs -- in the data
 * product -- so every downstream consumer sees the prototype's layout without
 * each of them re-deriving it.  It splits each flag_main_cluster cluster into
 * the pre-merge main (retained, keeps flag_main_cluster and the cluster
 * scalars incl. matched_flash_gid) and one new cluster per other pre-merge
 * member (flag_associated_cluster).
 *
 * Two ways to identify the main:
 *
 *  - mode "real" (default, EXACT): the per-blob "real_cluster_main" /
 *    "real_cluster_id" provenance written by merge_clusters() during the
 *    flash-time merge and persisted through the pctree tarball by
 *    MultiAlgBlobClustering's save_real_cluster_id (doc 38).  Blobs with
 *    real_cluster_main != 0 come from the merge's REPRESENTATIVE member --
 *    the same member whose flags (including flag_main_cluster) and flash the
 *    merged cluster carries -- so "the main" is by construction the cluster
 *    the bundle verdicts already describe.  The remaining blobs are grouped
 *    by their pre-merge real_cluster_id, which reproduces the pre-merge
 *    cluster partition exactly.
 *
 *  - mode "component" (PROXY, opt-in only): recompute connected components on
 *    the cluster's own graph and call the longest one the main.  This is a
 *    clustering decision, not bookkeeping -- the relaxed graph does not join
 *    the two halves of a cathode crosser -- so mode "real" does NOT fall back
 *    to it: a cluster without provenance (a tarball saved without
 *    save_real_cluster_id) is left alone with a warning.
 *
 * MUST run before the steiner stage: separate() moves blob nodes into fresh
 * cluster nodes, and node-local point clouds (steiner_pc) are NOT carried
 * over -- while the RETAINED main would keep a steiner_pc built from the
 * pre-split blob set.  A stale steiner_pc is erased defensively here, but the
 * cheap and correct placement is right after switch_scope.
 */
class ClusteringUnmergeBundle : public IConfigurable,
                                public Clus::IEnsembleVisitor,
                                public Aux::Logger,
                                private NeedDV,
                                private NeedPCTS {
public:
    ClusteringUnmergeBundle()
        : Aux::Logger("ClusteringUnmergeBundle", "clus")
    {}
    virtual ~ClusteringUnmergeBundle() {}

    virtual void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config);
        m_grouping_name = get<std::string>(config, "grouping", m_grouping_name);
        m_mode = get<std::string>(config, "mode", m_mode);
        m_pcarray_name = get<std::string>(config, "pcarray_name", m_pcarray_name);
        m_graph_name = get<std::string>(config, "graph_name", m_graph_name);
        m_require_in_scope = get<bool>(config, "require_in_scope", m_require_in_scope);
        // Which provenance pair to split on.  Defaults are the flash merge's
        // (doc 38/45) so existing configs are unchanged; set both to
        // "assoc_cluster_id"/"assoc_cluster_main" for a second instance that
        // undoes the isolated GROUPING instead (doc 52 Stage 3).  The two are
        // meant to run in order -- flash first (outer), isolated second (inner)
        // -- which reproduces the prototype's main_cluster +
        // additional_clusters layout that check_stm expects
        // (wire-cell-prod-stm.cxx:816-828).
        m_id_aname = get<std::string>(config, "id_aname", m_id_aname);
        m_main_aname = get<std::string>(config, "main_aname", m_main_aname);
        if (m_mode != "real" && m_mode != "component") {
            raise<ValueError>("ClusteringUnmergeBundle: unknown mode \"%s\" (want \"real\" or \"component\")",
                              m_mode);
        }
    }

    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["detector_volumes"] = "DetectorVolumes";
        cfg["pc_transforms"] = "PCTransformSet";
        cfg["grouping"] = m_grouping_name;
        // "real" = per-blob flash-merge provenance (exact, needs a pctree
        // saved with save_real_cluster_id); clusters without it are skipped,
        // NOT proxied.  "component" = longest connected component (opt-in).
        cfg["mode"] = m_mode;
        cfg["pcarray_name"] = m_pcarray_name;
        cfg["graph_name"] = m_graph_name;
        cfg["require_in_scope"] = m_require_in_scope;
        cfg["id_aname"] = m_id_aname;
        cfg["main_aname"] = m_main_aname;
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const {
        auto groupings = ensemble.with_name(m_grouping_name);
        if (groupings.empty()) {
            log->warn("no '{}' grouping found, skipping", m_grouping_name);
            return;
        }
        auto& grouping = *groupings.at(0);

        // Snapshot the cluster list: separate() adds children to the grouping.
        const std::vector<Cluster*> clusters = grouping.children();

        // Idents already in use, so a split part never shadows an existing
        // cluster.  Only ever grows.
        std::set<int> taken;
        for (const auto* c : clusters) taken.insert(c->ident());

        size_t n_split = 0, n_parts = 0, n_real = 0, n_proxy = 0;
        for (Cluster* cluster : clusters) {
            if (!cluster->get_flag(Flags::main_cluster)) continue;
            // switch_scope leaves out-of-volume shards carrying an inherited
            // flag_main_cluster; the taggers skip them (require_in_scope) and
            // so do we.
            if (m_require_in_scope
                && !cluster->get_scope_filter(cluster->get_default_scope())) continue;

            const size_t before = grouping.nchildren();
            bool used_real = false;
            if (process_cluster(grouping, cluster, taken, used_real)) {
                ++n_split;
                n_parts += grouping.nchildren() - before;
                if (used_real) ++n_real; else ++n_proxy;
            }
        }
        log->debug("unmerged {} main cluster(s) into {} associated cluster(s) "
                   "({} exact/provenance, {} proxy/component)",
                   n_split, n_parts, n_real, n_proxy);
    }

private:
    std::string m_grouping_name{"live"};
    std::string m_mode{"real"};
    std::string m_id_aname{"real_cluster_id"};
    std::string m_main_aname{"real_cluster_main"};
    std::string m_pcarray_name{"perblob"};
    std::string m_graph_name{"relaxed"};
    bool m_require_in_scope{true};

    /// Outcome of reading the flash-merge provenance.
    enum class Prov {
        unusable,   ///< absent / malformed: the caller may fall back to the proxy
        no_split,   ///< usable and says this cluster was never flash-merged
        ok,         ///< usable and `groups` is filled
    };

    /// Fill `groups` (parallel to the cluster's blob/children order) with -1
    /// for the main's blobs and a dense non-negative id per other pre-merge
    /// member.
    Prov groups_from_provenance(const Cluster* cluster, size_t nb,
                                std::vector<int>& groups) const {
        if (!cluster->has_pcarray<int>(m_main_aname, m_pcarray_name)) return Prov::unusable;
        if (!cluster->has_pcarray<int>(m_id_aname, m_pcarray_name)) return Prov::unusable;
        auto rmain = cluster->get_pcarray<int>(m_main_aname, m_pcarray_name);
        auto rid = cluster->get_pcarray<int>(m_id_aname, m_pcarray_name);
        if (rmain.size() != nb || rid.size() != nb) {
            log->warn("cluster {}: provenance size {}/{} != {} blobs, using the proxy",
                      cluster->ident(), rmain.size(), rid.size(), nb);
            return Prov::unusable;
        }
        size_t nmain = 0;
        for (size_t i = 0; i < nb; ++i) if (rmain[i] != 0) ++nmain;
        if (nmain == 0) {
            // switch_scope can strand the representative member entirely out
            // of scope, leaving an in-scope shard with no main blob.  The
            // proxy still has a defensible answer; the provenance does not.
            log->warn("cluster {}: no blob marked real_cluster_main, using the proxy",
                      cluster->ident());
            return Prov::unusable;
        }
        // Nothing was merged in.  This is NOT a fallback case: the cluster is
        // a single pre-merge cluster and must be left alone.  Splitting it on
        // connectivity would undo the clustering chain's own deliberate
        // long-range merges (extend / regular / parallel_prolong), which the
        // prototype keeps inside one PR3DCluster -- only the flash-time
        // bookkeeping merge is ours to undo.
        if (nmain == nb) return Prov::no_split;

        // Dense group ids by ascending pre-merge ident => deterministic.
        std::set<int> others;
        for (size_t i = 0; i < nb; ++i) if (rmain[i] == 0) others.insert(rid[i]);
        std::map<int, int> id2grp;
        int g = 0;
        for (int oid : others) id2grp[oid] = g++;

        groups.assign(nb, -1);
        for (size_t i = 0; i < nb; ++i) {
            if (rmain[i] == 0) groups[i] = id2grp.at(rid[i]);
        }
        return Prov::ok;
    }


    /// Proxy: connected components of the cluster's own graph, longest = main.
    bool groups_from_components(Cluster* cluster, size_t nb,
                                std::vector<int>& groups) const {
        auto b2g = cluster->connected_blobs(m_dv, m_pcts, m_graph_name);
        if (b2g.size() != nb) return false;
        std::set<int> ids;
        for (int id : b2g) if (id >= 0) ids.insert(id);
        if (ids.size() < 2) return false;

        double max_length = -1;
        int longest = -1;
        for (int id : ids) {                 // ascending => ties keep lowest id
            const double len = get_length(cluster, b2g, id);
            if (len > max_length) { max_length = len; longest = id; }
        }

        std::map<int, int> id2grp;
        int g = 0;
        for (int id : ids) if (id != longest) id2grp[id] = g++;

        groups.assign(nb, -1);
        for (size_t i = 0; i < nb; ++i) {
            if (b2g[i] < 0 || b2g[i] == longest) continue;   // unassigned stays with the main
            groups[i] = id2grp.at(b2g[i]);
        }
        return true;
    }

    /// Split one merged main.  Returns true if a split happened.
    bool process_cluster(Grouping& grouping, Cluster*& cluster,
                         std::set<int>& taken, bool& used_real) const {
        const size_t nb = cluster->nchildren();
        if (nb < 2) return false;

        std::vector<int> groups;
        used_real = false;
        if (m_mode == "real") {
            const Prov p = groups_from_provenance(cluster, nb, groups);
            // No usable provenance => SKIP.  Deliberately not a fallback to
            // the component proxy: undoing the flash merge is bookkeeping,
            // but splitting on graph connectivity is a clustering decision
            // (it breaks cathode crossers, whose two halves the relaxed graph
            // does not join).  Reusing a pctree saved without
            // save_real_cluster_id must therefore be a no-op, not a silent
            // re-clustering.  Ask for the proxy explicitly with mode
            // "component".
            if (p != Prov::ok) {
                if (p == Prov::unusable) {
                    log->warn("cluster {}: no flash-merge provenance (save the pctree with "
                              "save_real_cluster_id, or set mode=\"component\"); not split",
                              cluster->ident());
                }
                return false;
            }
            used_real = true;
        }
        if (!used_real) {
            groups.clear();
            if (!groups_from_components(cluster, nb, groups)) return false;
        }

        int ngrp = 0;
        for (int g : groups) if (g >= 0) ngrp = std::max(ngrp, g + 1);
        if (ngrp < 1) return false;

        // Snapshot the per-blob PC before the blobs move away.  It is the only
        // node-local PC that must survive the split (Bee colors read
        // real_cluster_id; TaggerCheckTGM main_component_mode="real" reads
        // real_cluster_main), and it MUST be re-carved on the RETAINED cluster
        // too -- separate(remove=false) leaves it holding a full-length array
        // over a now-shorter blob list, which is a size-mismatch landmine for
        // every consumer that assumes it is parallel to blob order.
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
        auto splits = grouping.separate(cluster, groups);   // remove=false: main stays put

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

        // The retained main: re-carve the provenance and drop any node-local
        // PC that was derived from the pre-split blob set.  With the
        // documented placement (before the steiner stage) there is none; the
        // erase keeps a mis-ordered pipeline from fitting a stale point cloud.
        carve(cluster, -1);
        if (cluster->value().local_pcs().erase("steiner_pc")) {
            log->warn("cluster {}: erased a steiner_pc built before the split -- "
                      "run this visitor BEFORE the steiner stage", main_ident);
        }
        if (cluster->get_flag(Flags::associated_cluster)) {
            cluster->set_flag(Flags::associated_cluster, 0);
        }

        // The companions.  separate() -> from() copied the merged cluster's
        // flags and cluster_scalar (cluster_t0 / flash / matched_flash_gid),
        // which is what the taggers use to pair a companion with its main --
        // but that also hands each part flag_main_cluster, which would make
        // every fragment look like a bundle main to STM/TGM/FC.
        size_t nb_parts = 0;
        for (const auto& [gid, part] : splits) nb_parts += part->nchildren();
        log->debug("cluster {}: {} blobs -> main {} + {} associated cluster(s) holding {} ({} mode)",
                   main_ident, nb, cluster->nchildren(), splits.size(), nb_parts,
                   used_real ? "real" : "component");
        if (cluster->nchildren() + nb_parts != nb) {
            log->error("cluster {}: BLOB LOSS in separate(): {} + {} != {}",
                       main_ident, cluster->nchildren(), nb_parts, nb);
        }
        int sub_id = 1;
        for (auto& [gid, part] : splits) {
            part->set_flag(Flags::main_cluster, 0);
            part->set_flag(Flags::associated_cluster);
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
