// doc pdvd/28 (PDVD PR perf round 1, memory): release the per-cluster graphs,
// graph-algorithm caches (and optionally the ClusterCache) plus the fit-side
// scratch of the TrackFitting objects parked in the grouping slots, at a
// chosen point of the PR pipeline -- meant to sit right after
// tagger_check_neutrino, before the writers.
//
// Output-neutral by construction: nothing downstream of the neutrino stage in
// the PDVD pipeline (PdvdPrMagnifyTrackingVisitor, UbooneTaggerOutputVisitor,
// PrDisplayDump, PdvdMagnifyTrackingVisitor, the Bee dump) reads a stored
// graph or a GraphAlgorithms cache, and the fitter getters they call keep
// their data (TrackFitting::release_fit_scratch).  Local point clouds
// (steiner_pc, stm_fit, ...) are left alone: PrDisplayDump and the STM
// magnify writer read them.  A pipeline that does not name this visitor is
// untouched (nothing in the compiled config either).
//
// Config:
//   grouping        "live"  which grouping's clusters to release
//   graphs          true    take_graph + remove_graph_algorithms for every
//                           stored graph and cached GraphAlgorithms (the
//                           ClusteringProtectBundle split sequence)
//   fitter_scratch  true    release_fit_scratch() on the unnamed slot, every
//                           "nu<i>" slot and the slots listed in `slots`
//   slots           ["stm"] extra named TrackFitting slots
//   cluster_cache   false   also Cluster::invalidate_cache() (rebuilds lazily
//                           and deterministically if a later reader needs it)
#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Facade.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/MemUsage.h"

#include <set>

class ClusteringReleaseCaches;
WIRECELL_FACTORY(ClusteringReleaseCaches, ClusteringReleaseCaches,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

class ClusteringReleaseCaches : public IConfigurable, public Clus::IEnsembleVisitor {
  public:
    ClusteringReleaseCaches() {}
    virtual ~ClusteringReleaseCaches() {}

    virtual void configure(const WireCell::Configuration& cfg)
    {
        m_grouping = get(cfg, "grouping", m_grouping);
        m_graphs = get(cfg, "graphs", m_graphs);
        m_fitter_scratch = get(cfg, "fitter_scratch", m_fitter_scratch);
        m_cluster_cache = get(cfg, "cluster_cache", m_cluster_cache);
        if (cfg.isMember("slots")) {
            m_slots.clear();
            for (const auto& s : cfg["slots"]) m_slots.push_back(s.asString());
        }
    }

    virtual Configuration default_configuration() const
    {
        Configuration cfg;
        cfg["grouping"] = m_grouping;
        cfg["graphs"] = m_graphs;
        cfg["fitter_scratch"] = m_fitter_scratch;
        cfg["cluster_cache"] = m_cluster_cache;
        cfg["slots"] = Json::arrayValue;
        for (const auto& s : m_slots) cfg["slots"].append(s);
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const
    {
        auto log = Log::logger("clus");
        auto groupings = ensemble.with_name(m_grouping);
        if (groupings.empty()) {
            log->warn("ClusteringReleaseCaches: no grouping '{}'", m_grouping);
            return;
        }
        auto* grouping = groupings.front();
        const double res0 = memusage_resident();
        size_t ngraphs = 0, ngalgs = 0, nclusters = 0, nfitters = 0;
        if (m_graphs || m_cluster_cache) {
            for (auto* cluster : grouping->children()) {
                ++nclusters;
                if (m_graphs) {
                    std::vector<std::string> gnames;
                    for (const auto& [gname, gr] : cluster->graph_store()) gnames.push_back(gname);
                    for (const auto& gname : gnames) {
                        cluster->take_graph(gname);
                        cluster->remove_graph_algorithms(gname);
                        ++ngraphs;
                    }
                    for (const auto& gname : cluster->get_cached_graph_algorithms()) {
                        cluster->remove_graph_algorithms(gname);
                        ++ngalgs;
                    }
                }
                if (m_cluster_cache) cluster->invalidate_cache();
            }
        }
        if (m_fitter_scratch) {
            std::set<TrackFitting*> done;
            auto release = [&](std::shared_ptr<TrackFitting> tf) {
                if (!tf || done.count(tf.get())) return;
                tf->release_fit_scratch();
                done.insert(tf.get());
                ++nfitters;
            };
            release(grouping->get_track_fitting());
            for (int i = 0;; ++i) {
                auto tf = grouping->get_track_fitting("nu" + std::to_string(i));
                if (!tf) break;
                release(tf);
            }
            for (const auto& s : m_slots) release(grouping->get_track_fitting(s));
        }
        const double res1 = memusage_resident();
        log->debug("ClusteringReleaseCaches: grouping '{}': {} cluster(s), dropped {} graph(s) + {} cached algorithm(s), "
                   "released {} fitter(s); resident {:.0f} -> {:.0f} MB",
                   m_grouping, nclusters, ngraphs, ngalgs, nfitters, res0 / 1024.0, res1 / 1024.0);
    }

  private:
    std::string m_grouping{"live"};
    bool m_graphs{true};
    bool m_fitter_scratch{true};
    bool m_cluster_cache{false};
    std::vector<std::string> m_slots{"stm"};
};
