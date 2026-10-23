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
// everything -- unless truth_iids is set (pr_ml docs/23 sec 11.9, the SBND
// "perfect Q/L" stand-in for light-less simulation with CORSIKA cosmics):
// then only the clusters whose labelled points (the per-point "truth" cloud of
// TruthPointLabeler, key "iid" = truth_nu row) belong mostly (>= truth_min_frac)
// to one of the listed interactions form the bundle; the main is the longest
// of them, and every other cluster is left unmatched (main/associated 0,
// matched_flash_gid -1, as QLMatching leaves a cluster without a flash).
// truth_iids empty (the default) => the legacy everything-bundle, unchanged.  A NEW component: absent from every existing pipeline => no
// other detector's output changes.  Deterministic: tree order, ties broken by
// npoints then cluster id (never pointer order).

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/ClusteringFuncs.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Exceptions.h"

#include <algorithm>
#include <vector>

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
        m_truth_iids.clear();
        for (const auto& v : config["truth_iids"]) m_truth_iids.push_back(v.asInt());
        m_truth_pc = get<std::string>(config, "truth_pc", m_truth_pc);
        m_truth_min_frac = get<double>(config, "truth_min_frac", m_truth_min_frac);
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
        cfg["truth_iids"] = Json::arrayValue;
        for (int v : m_truth_iids) cfg["truth_iids"].append(v);
        cfg["truth_pc"] = m_truth_pc;
        cfg["truth_min_frac"] = m_truth_min_frac;
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

        // truth_iids: the bundle is the clusters of the listed interactions
        std::vector<char> in_bundle(clusters.size(), 1);
        if (!m_truth_iids.empty()) {
            for (size_t ic = 0; ic < clusters.size(); ++ic) {
                size_t nlab = 0, nsel = 0;
                for (Blob* blob : clusters[ic]->children()) {
                    const auto& lpcs = blob->value().local_pcs();
                    auto it = lpcs.find(m_truth_pc);
                    if (it == lpcs.end()) continue;
                    auto arr = it->second.get("iid");
                    if (!arr) continue;
                    for (int iid : arr->elements<int>()) {
                        if (iid < 0) continue;
                        ++nlab;
                        if (std::find(m_truth_iids.begin(), m_truth_iids.end(), iid) != m_truth_iids.end()) ++nsel;
                    }
                }
                in_bundle[ic] = (nlab > 0 && double(nsel) >= m_truth_min_frac * double(nlab)) ? 1 : 0;
            }
        }

        // The main: longest; ties -> more points -> smaller cluster id.
        Cluster* main = nullptr;
        double main_len = 0;
        int main_np = 0, main_id = 0;
        size_t nbundle = 0;
        for (size_t ic = 0; ic < clusters.size(); ++ic) {
            if (!in_bundle[ic]) continue;
            ++nbundle;
            Cluster* cluster = clusters[ic];
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

        for (size_t ic = 0; ic < clusters.size(); ++ic) {
            Cluster* cluster = clusters[ic];
            const bool is_main = (cluster == main);
            if (!in_bundle[ic]) {   // truth_iids only: a cluster with no flash
                cluster->set_flag(Flags::main_cluster, 0);
                cluster->set_flag(Flags::associated_cluster, 0);
                cluster->set_scalar<int>("matched_flash_gid", -1);
                cluster->set_cluster_t0(m_cluster_t0);
                continue;
            }
            cluster->set_flag(Flags::main_cluster, is_main ? 1 : 0);
            cluster->set_flag(Flags::associated_cluster, is_main ? 0 : 1);
            if (m_flag_beam_flash) cluster->set_flag(Flags::beam_flash, 1);
            cluster->set_scalar<int>("matched_flash_gid", m_bundle_gid);
            cluster->set_cluster_t0(m_cluster_t0);
        }
        if (!m_truth_iids.empty()) {
            logger()->info("ClusteringBundleNoLight: truth_iids: {} of {} cluster(s) in the bundle", nbundle,
                           clusters.size());
            if (!main) {
                logger()->info("ClusteringBundleNoLight: no cluster of the listed interactions, no main");
                return;
            }
        }
        logger()->info("ClusteringBundleNoLight: {} cluster(s) -> main id {} (L {:.1f} cm, {} points),"
                       " {} associated, gid {}, t0 {:.3f} us",
                       clusters.size(), main_id, main_len / units::cm, main_np, nbundle - 1,
                       m_bundle_gid, m_cluster_t0 / units::us);
    }

private:
    std::string m_grouping{"live"};
    int m_bundle_gid{0};
    double m_cluster_t0{0.0};
    bool m_flag_beam_flash{false};
    // pr_ml docs/23 sec 11.9: interactions (truth_nu rows) forming the bundle;
    // empty = legacy (every cluster)
    std::vector<int> m_truth_iids;
    std::string m_truth_pc{"truth"};
    double m_truth_min_frac{0.5};
};

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
