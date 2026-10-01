// ClusteringBundleNoLight -- make the whole event ONE main + associated bundle
// for light-less simulation (DUNE FD-HD 1x2x6 neutrino PR; wcfm/docs/24-25).
//
// The neutrino PR (TaggerCheckNeutrino) takes one main cluster plus the
// companions of its bundle, and both come from Q/L matching: the main carries
// Flags::main_cluster, the companions Flags::associated_cluster, and all share
// matched_flash_gid.  With no light simulated, nothing sets them:
//   * ClusteringFlagMatchedMains(flag_unmatched) makes every cluster its own
//     main with gid -1, which nu_per_bundle drops (TaggerCheckNeutrino.cxx
//     gid < 0 check) and which has no companions;
//   * ClusteringProtectBundle opens only bundles with gid >= 0.
// A simulated FD event is one neutrino interaction at t0 = 0, so this visitor
// writes what a matched bundle would carry: the longest cluster is the main,
// every other cluster an associated companion, all with one gid in
// [0, 1000000) (TaggerCheckNeutrino reads gid >= 1000000 as "unmatched") and
// cluster_t0 = the configured t0.  The flags, gid and t0 are written on EVERY
// cluster, so all clusters carry the same scalar keys.  Run it BEFORE
// switch_scope, which builds the T0-corrected scope from cluster_t0.
//
// Not for data or for events with more than one interaction: it bundles
// everything.  A NEW component: absent from every existing pipeline => no
// other detector's output changes.  Deterministic: tree order, ties broken by
// npoints then cluster id (never pointer order).

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/ClusteringFuncs.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Exceptions.h"

class ClusteringBundleNoLight;
WIRECELL_FACTORY(ClusteringBundleNoLight, ClusteringBundleNoLight,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

static Log::logptr_t logger() {
    static Log::logptr_t l = Log::logger("clus.BundleNoLight");
    return l;
}

class ClusteringBundleNoLight : public IConfigurable, public Clus::IEnsembleVisitor {
public:
    ClusteringBundleNoLight() {}
    virtual ~ClusteringBundleNoLight() {}

    void configure(const WireCell::Configuration& config) {
        m_grouping = get<std::string>(config, "grouping", m_grouping);
        m_bundle_gid = get<int>(config, "bundle_gid", m_bundle_gid);
        m_cluster_t0 = get<double>(config, "cluster_t0", m_cluster_t0);     // internal units
        m_flag_beam_flash = get<bool>(config, "flag_beam_flash", m_flag_beam_flash);
        if (m_bundle_gid < 0 || m_bundle_gid >= 1000000) {
            THROW(ValueError() << errmsg{"ClusteringBundleNoLight: bundle_gid must be in [0, 1000000)"});
        }
    }

    WireCell::Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping;
        cfg["bundle_gid"] = m_bundle_gid;
        cfg["cluster_t0"] = m_cluster_t0;
        cfg["flag_beam_flash"] = m_flag_beam_flash;
        return cfg;
    }

    void visit(Ensemble& ensemble) const {
        auto vec = ensemble.with_name(m_grouping);
        if (vec.empty()) {
            logger()->warn("ClusteringBundleNoLight: no '{}' grouping found, skipping", m_grouping);
            return;
        }
        Grouping& grouping = *vec.at(0);
        const auto& clusters = grouping.children();   // tree order: deterministic
        if (clusters.empty()) {
            logger()->warn("ClusteringBundleNoLight: '{}' grouping has no clusters, skipping", m_grouping);
            return;
        }

        // The main: longest; ties -> more points -> smaller cluster id.
        Cluster* main = nullptr;
        double main_len = 0;
        int main_np = 0, main_id = 0;
        for (Cluster* cluster : clusters) {
            const double len = cluster->get_length();
            const int np = cluster->npoints();
            const int id = cluster->get_cluster_id();
            const bool better = !main || len > main_len ||
                                (len == main_len && (np > main_np || (np == main_np && id < main_id)));
            if (better) {
                main = cluster;
                main_len = len;
                main_np = np;
                main_id = id;
            }
        }

        for (Cluster* cluster : clusters) {
            const bool is_main = (cluster == main);
            cluster->set_flag(Flags::main_cluster, is_main ? 1 : 0);
            cluster->set_flag(Flags::associated_cluster, is_main ? 0 : 1);
            if (m_flag_beam_flash) cluster->set_flag(Flags::beam_flash, 1);
            cluster->set_scalar<int>("matched_flash_gid", m_bundle_gid);
            cluster->set_cluster_t0(m_cluster_t0);
        }
        logger()->info("ClusteringBundleNoLight: {} cluster(s) -> main id {} (L {:.1f} cm, {} points),"
                       " {} associated, gid {}, t0 {:.3f} us",
                       clusters.size(), main_id, main_len / units::cm, main_np, clusters.size() - 1,
                       m_bundle_gid, m_cluster_t0 / units::us);
    }

private:
    std::string m_grouping{"live"};
    int m_bundle_gid{0};
    double m_cluster_t0{0.0};
    bool m_flag_beam_flash{false};
};

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
