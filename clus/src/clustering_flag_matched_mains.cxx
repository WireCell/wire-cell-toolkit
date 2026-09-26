// ClusteringFlagMatchedMains -- flag every flash-matched cluster as a "main"
// (doc pdvd/25 M3).
//
// The cosmic taggers and the neutrino PR admit only clusters carrying
// Flags::main_cluster and gather companions by matched_flash_gid
// (TaggerCheckSTM.cxx visit).  QLMatching sets that flag ONLY on the main of a
// MicroBooNE-style main+associated decomposition, i.e. when the clustering
// left an "isolated"/"perblob" array (SBND: examine_bundles + isolated with
// save_assoc_id).  The PDVD chain never builds that decomposition, so after
// Q/L matching every PDVD cluster has cluster_t0 / flash / matched_flash_gid
// but flag_main_cluster == 0 for all, and the PR tail evaluates nothing.
//
// This visitor closes the gap without touching match/: each cluster whose
// matched_flash_gid >= 0 (and, by default, whose cluster_t0 is a real match
// time) becomes a main.  In a detector without a beam every matched cluster
// is its own candidate; the taggers then evaluate each one and treat the
// other clusters on the same flash as its companions.  A NEW component:
// absent from every existing pipeline => no other detector's output changes.
// Deterministic: iterates the grouping's children in tree order.

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/ClusteringFuncs.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"

class ClusteringFlagMatchedMains;
WIRECELL_FACTORY(ClusteringFlagMatchedMains, ClusteringFlagMatchedMains,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

static Log::logptr_t logger() {
    static Log::logptr_t l = Log::logger("clus.FlagMatchedMains");
    return l;
}

class ClusteringFlagMatchedMains : public IConfigurable, public Clus::IEnsembleVisitor {
public:
    ClusteringFlagMatchedMains() {}
    virtual ~ClusteringFlagMatchedMains() {}

    void configure(const WireCell::Configuration& config) {
        m_grouping = get<std::string>(config, "grouping", m_grouping);
        m_require_t0 = get<bool>(config, "require_t0", m_require_t0);
        m_min_length = get<double>(config, "min_length", m_min_length);   // internal units
        m_skip_flagged = get<bool>(config, "skip_flagged", m_skip_flagged);
        // doc pdvd/101: flag_unmatched (C++ default false) lets a cluster that
        // Q/L matching never saw (matched_flash_gid < 0) become a main too.  It
        // exists for LIGHT-LESS SIMULATION only -- a single simulated muon has no
        // flash, so without it the taggers and the track fit evaluate nothing.
        // An unmatched cluster keeps Cluster::get_cluster_t0()'s unstamped
        // default of 0 (QLMatching's -1e12 sentinel is only written when Q/L
        // matching runs), which is the true t0 of a t=0 simulation, so
        // require_t0 passes.  Absent key => the legacy path, byte-identical.
        m_flag_unmatched = get<bool>(config, "flag_unmatched", m_flag_unmatched);
    }

    WireCell::Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping;
        cfg["require_t0"] = m_require_t0;
        cfg["min_length"] = m_min_length;
        cfg["skip_flagged"] = m_skip_flagged;
        cfg["flag_unmatched"] = m_flag_unmatched;
        return cfg;
    }

    void visit(Ensemble& ensemble) const {
        auto vec = ensemble.with_name(m_grouping);
        if (vec.empty()) {
            logger()->warn("ClusteringFlagMatchedMains: no '{}' grouping found, skipping", m_grouping);
            return;
        }
        Grouping& grouping = *vec.at(0);
        int n_matched = 0, n_flagged = 0, n_already = 0, n_short = 0, n_not0 = 0;
        int n_unmatched = 0;                                    // doc pdvd/101 (flag_unmatched only)
        for (Cluster* cluster : grouping.children()) {          // tree order: deterministic
            const int gid = cluster->get_scalar<int>("matched_flash_gid", -1);
            if (gid < 0) {
                if (!m_flag_unmatched) continue;                // legacy: unmatched clusters are never mains
                ++n_unmatched;
            }
            else {
                ++n_matched;
            }
            if (m_require_t0 && !(cluster->get_cluster_t0() > -1e11)) { ++n_not0; continue; }
            if (m_min_length > 0 && cluster->get_length() < m_min_length) { ++n_short; continue; }
            if (cluster->get_flag(Flags::main_cluster)) { ++n_already; if (m_skip_flagged) continue; }
            cluster->set_flag(Flags::main_cluster, 1);
            ++n_flagged;
        }
        logger()->info("ClusteringFlagMatchedMains: {} cluster(s), {} matched (gid>=0), {} flagged main"
                       " ({} already, {} without t0, {} below min_length {:.1f} cm)",
                       grouping.children().size(), n_matched, n_flagged, n_already, n_not0, n_short,
                       m_min_length / units::cm);
        if (m_flag_unmatched) {
            // separate line so the legacy line above stays byte-identical when off
            logger()->info("ClusteringFlagMatchedMains: flag_unmatched admitted {} unmatched cluster(s)"
                           " (counted in the flagged/without-t0/below-min_length totals above)", n_unmatched);
        }
    }

private:
    std::string m_grouping{"live"};
    bool m_require_t0{true};
    double m_min_length{0.0};
    bool m_skip_flagged{true};
    bool m_flag_unmatched{false};   // doc pdvd/101: light-less simulation only
};

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
