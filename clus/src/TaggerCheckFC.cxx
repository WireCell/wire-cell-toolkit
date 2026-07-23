/// TaggerCheckFC: fully-contained (FC) tagger.
///
/// Surfaces, as a first-class cluster flag, the containment verdict the
/// prototype records as event_type bit 2 / the eval-tree variable match_isFC
/// (prototype uboone_nusel_app/apps/prod-wire-cell-matching-nusel-port.cxx
/// lines 1002-1007, calling WCP2dToy::ToyFiducial::check_fully_contained,
/// 2dtoy/src/ToyFiducial.cxx:816).
///
/// The check itself is NOT re-implemented here: it is Facade::cluster_fc_check
/// (clus/src/Clustering_Util.cxx:75), the existing two-round steiner-boundary
/// port already used by TaggerCheckSTM (as an early exit whose verdict is
/// discarded) and by TaggerCheckNeutrino (to fill tagger_info.match_isFC).
/// This visitor is the missing third consumer: the one that records the
/// verdict where a per-bundle label table can see it.  See
/// wcp-porting-img/sbnd/sbnd_xin/docs/24_fc-lm-label-audit.md.
///
/// Why a separate visitor rather than promoting TaggerCheckSTM's internal
/// result: TaggerCheckSTM skips mains already flagged TGM (TaggerCheckSTM.cxx
/// :144), so its FC computation never runs for through-going bundles.  FC is
/// an orthogonal property of every main cluster, so it gets its own pass over
/// the full main-cluster population.  The cost is that non-TGM mains have
/// cluster_fc_check run twice per event.
///
/// Differences from the prototype, by design (inherited from cluster_fc_check):
///  - Endpoints come from the steiner boundary (two rounds, flag_cosmic
///    true then false), not from get_extreme_wcps().
///  - offset_x = 0 everywhere: this runs after switch_scope, so point
///    coordinates are already T0-corrected by the matched flash time.  Same
///    convention as TaggerCheckTGM.
///  - No flag=2 retry against the original pre-merge main cluster.
///  - The prototype's _fc_breakdown FV/SP/DC failure mask has no counterpart;
///    FCCheckResult carries no failure reason.
/// The verdict is therefore *an* FC definition consistent with the one STM and
/// TaggerCheckNeutrino already act on -- it is not bit-comparable to WCP's.
///
/// Sets Flags::FC ("FC"), the tagger-computed sibling of Flags::TGM/Flags::STM.
/// It is deliberately NOT the lowercase Flags::fully_contained, which is
/// reserved for uBooNE verdicts imported through ClusteringTaggerFlagTransfer.

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"

class TaggerCheckFC;
WIRECELL_FACTORY(TaggerCheckFC, TaggerCheckFC,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

static auto f_log = WireCell::Log::logger("clus.NeutrinoPattern");

class TaggerCheckFC : public IConfigurable, public Clus::IEnsembleVisitor,
                      private Clus::NeedDV {
public:
    TaggerCheckFC() {}
    virtual ~TaggerCheckFC() {}

    virtual void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        m_grouping_name = get<std::string>(config, "grouping", m_grouping_name);
        // require_in_scope (default false = same convention as the other
        // taggers): also require the cluster to pass the default-scope filter
        // set by switch_scope.  switch_scope SEPARATES out-of-volume blobs into
        // their own cluster which stays in the grouping and inherits
        // flag_main_cluster; those shards are non-physical and are excluded
        // from the label table, so evaluating them here only wastes time.
        m_require_in_scope = get(config, "require_in_scope", m_require_in_scope);
    }

    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping_name;
        cfg["detector_volumes"] = "DetectorVolumes";
        cfg["require_in_scope"] = m_require_in_scope;
        return cfg;
    }

    virtual void visit(Ensemble& ensemble) const {
        auto groupings = ensemble.with_name(m_grouping_name);
        if (groupings.empty()) return;
        auto& grouping = *groupings.at(0);

        std::vector<Cluster*> main_clusters;
        int n_out_of_scope = 0;
        for (auto* cluster : grouping.children()) {
            if (!cluster->get_flag(Flags::main_cluster)) continue;
            if (m_require_in_scope && !cluster->get_scope_filter(cluster->get_default_scope())) {
                ++n_out_of_scope;
                continue;
            }
            main_clusters.push_back(cluster);
        }
        if (n_out_of_scope) {
            SPDLOG_LOGGER_INFO(f_log, "visit: TaggerCheckFC: skipped {} out-of-scope main cluster(s)",
                               n_out_of_scope);
        }
        if (main_clusters.empty()) return;

        for (auto* main_cluster : main_clusters) {
            bool is_fc = false;
            try {
                // Conservative on failure: cluster_fc_check itself returns
                // is_fc=false when the cluster has no steiner_pc or when the
                // grouping carries no FiducialUtils (i.e. when the pipeline is
                // missing the steiner / fiducialutils stages).
                is_fc = Facade::cluster_fc_check(*main_cluster, m_dv).is_fc;
            }
            catch (const std::exception& err) {
                SPDLOG_LOGGER_WARN(f_log, "visit: TaggerCheckFC: cluster {} check failed: {}",
                                   main_cluster->ident(), err.what());
            }
            if (is_fc) main_cluster->set_flag(Flags::FC);
            SPDLOG_LOGGER_INFO(f_log, "visit: TaggerCheckFC: cluster {} → FC={}",
                               main_cluster->ident(), is_fc);
        }
    }

private:
    std::string m_grouping_name{"live"};
    bool m_require_in_scope{false};
};
