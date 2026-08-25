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
///
/// Fiducial volume.  By default (no "fiducial" key) containment is
/// FiducialUtils' -- the union of per-(apa,face) sensitive volumes with no
/// margin -- i.e. exactly what TaggerCheckSTM and TaggerCheckNeutrino see.
/// That volume is NOT the one TaggerCheckTGM uses, which made FC and TGM
/// mutually non-exclusive: TGM insets every wall by 2-3 cm via fv_tolerance
/// while FiducialUtils insets none, so a track ending 1 cm inside a wall is an
/// exiter for TGM and contained for FC; and FiducialUtils excludes the CPA slab
/// (|x| < 0.45 cm) that TGM's single cross-cathode box spans.  Setting
/// "fiducial" + "fv_tolerance" to the same values TaggerCheckTGM is given makes
/// the two verdicts comparable.  See sbnd_xin/docs/27_fc-tgm-consistent-fv.md.

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Logging.h"

class TaggerCheckFC;
WIRECELL_FACTORY(TaggerCheckFC, TaggerCheckFC,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

static auto f_log = WireCell::Log::logger("clus.NeutrinoPattern");

class TaggerCheckFC : public IConfigurable, public Clus::IEnsembleVisitor,
                      private Clus::NeedDV, private Clus::NeedFiducial {
public:
    TaggerCheckFC() {}
    virtual ~TaggerCheckFC() {}

    virtual void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        m_grouping_name = get<std::string>(config, "grouping", m_grouping_name);
        // "fiducial" absent => keep cluster_fc_check's historical FiducialUtils
        // containment.  NeedFiducial::configure is called ONLY when the key is
        // present: it falls back to looking up the type-only name
        // "DetectorVolumes", which is not an instantiated component in the SBND
        // PR config (it is named dv-apa0-1) and would throw.
        m_use_fiducial = !config["fiducial"].isNull();
        if (m_use_fiducial) NeedFiducial::configure(config);
        // fv_tolerance: boundary margins for the direct containment tests,
        // FiducialUtils convention [x_lo,x_hi,y_lo,y_hi,z_lo,z_hi], negative =
        // inset.  Set together with "fiducial" to evaluate containment exactly
        // the way TaggerCheckTGM does.
        auto tol = config["fv_tolerance"];
        if (!tol.isNull() && tol.isArray()) {
            m_fv_tolerance.clear();
            for (const auto& t : tol) m_fv_tolerance.push_back(t.asDouble());
        }
        // require_in_scope (default false = same convention as the other
        // taggers): also require the cluster to pass the default-scope filter
        // set by switch_scope.  switch_scope SEPARATES out-of-volume blobs into
        // their own cluster which stays in the grouping and inherits
        // flag_main_cluster; those shards are non-physical and are excluded
        // from the label table, so evaluating them here only wastes time.
        m_require_in_scope = get(config, "require_in_scope", m_require_in_scope);
        // evaluate_demoted_mains (default false = historical behavior): also
        // evaluate a cluster carrying Flags::demoted_main -- a split part that
        // was ITSELF a matched Q/L bundle main before the flash-time merge
        // demoted it (ClusteringUnmergeBundle restore_demoted_mains, doc pr/20
        // Part I).  Such a cluster is not a fragment; it is exactly the
        // population these cuts were tuned on, and today no cosmic tagger ever
        // looks at it.  The scope filter and beam-window gate below apply to it
        // unchanged.  Default false => no cluster is added => byte-identical.
        m_evaluate_demoted_mains = get<bool>(config, "evaluate_demoted_mains", m_evaluate_demoted_mains);

        // beam_window_only + beam_window_low/high (all C++-default off =>
        // byte-identical legacy behavior): evaluate only the beam-coincident
        // bundle, i.e. mains whose matched flash time (cluster_t0) lies in
        // [low, high).  Same gate and same window as the steiner stage and
        // TaggerCheckTGM/TaggerCheckSTM.  cluster_fc_check reads only the
        // cluster handed to it, so gating removes verdicts without changing any.
        m_beam_window_only = get(config, "beam_window_only", m_beam_window_only);
        m_beam_window_low = get(config, "beam_window_low", m_beam_window_low);
        m_beam_window_high = get(config, "beam_window_high", m_beam_window_high);
        if (m_beam_window_only && !(m_beam_window_low < m_beam_window_high)) {
            SPDLOG_LOGGER_WARN(f_log, "configure: TaggerCheckFC: beam_window_only set but the window is empty ([{}, {}) us) -- gate disabled, every main evaluated",
                               m_beam_window_low/units::us, m_beam_window_high/units::us);
        }
    }

    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping_name;
        cfg["detector_volumes"] = "DetectorVolumes";
        // Empty/absent => historical FiducialUtils containment (the union of
        // per-face sensitive volumes, no margin).  Naming an IFiducial here
        // redirects the direct containment tests to it.
        cfg["fiducial"] = Json::Value();   // null = use FiducialUtils
        cfg["fv_tolerance"] = Json::Value(Json::arrayValue);
        cfg["require_in_scope"] = m_require_in_scope;
        cfg["evaluate_demoted_mains"] = m_evaluate_demoted_mains;
        // Evaluate only mains whose cluster_t0 is in [low, high); false (or an
        // empty window) = every main, i.e. legacy behavior.
        cfg["beam_window_only"] = m_beam_window_only;
        cfg["beam_window_low"] = m_beam_window_low;
        cfg["beam_window_high"] = m_beam_window_high;
        return cfg;
    }

    // Event tag for the log lines nusel_extract.py parses.  Empty in a
    // one-event-per-process job (log text unchanged); "evt<ID> " in a group
    // job, where the log's several spdlog sinks do not arrive in time order and
    // cluster idents restart every event.  See
    // sbnd_xin/scripts/multi/slice_group_log.py.
    mutable std::string m_evt_tag;

    virtual void visit(Ensemble& ensemble) const {
        m_evt_tag = ensemble.rse_valid()
            ? WireCell::String::format("evt%d ", ensemble.ident()) : std::string();

        auto groupings = ensemble.with_name(m_grouping_name);
        if (groupings.empty()) return;
        auto& grouping = *groupings.at(0);

        const bool beam_gate = m_beam_window_only && m_beam_window_low < m_beam_window_high;

        std::vector<Cluster*> main_clusters;
        int n_out_of_scope = 0;
        int n_out_of_window = 0;
        int n_demoted = 0;
        for (auto* cluster : grouping.children()) {
            const bool demoted = m_evaluate_demoted_mains
                && !cluster->get_flag(Flags::main_cluster)
                && cluster->get_flag(Flags::demoted_main);
            if (!cluster->get_flag(Flags::main_cluster) && !demoted) continue;
            if (m_require_in_scope && !cluster->get_scope_filter(cluster->get_default_scope())) {
                ++n_out_of_scope;
                continue;
            }
            if (beam_gate) {
                const double t0 = cluster->get_cluster_t0();
                if (t0 < m_beam_window_low || t0 >= m_beam_window_high) {
                    ++n_out_of_window;
                    continue;
                }
            }
            main_clusters.push_back(cluster);
            if (demoted) ++n_demoted;
        }
        if (n_demoted) {
            SPDLOG_LOGGER_INFO(f_log, "visit: TaggerCheckFC: evaluate_demoted_mains: {} demoted main(s) added",
                               n_demoted);
        }
        if (n_out_of_scope) {
            SPDLOG_LOGGER_INFO(f_log, "visit: TaggerCheckFC: skipped {} out-of-scope main cluster(s)",
                               n_out_of_scope);
        }
        if (beam_gate) {
            SPDLOG_LOGGER_INFO(f_log, "visit: TaggerCheckFC: beam_window_only [{:.3f}, {:.3f}) us: {} main(s) evaluated, {} out of window",
                               m_beam_window_low/units::us, m_beam_window_high/units::us,
                               main_clusters.size(), n_out_of_window);
        }
        if (main_clusters.empty()) return;

        for (auto* main_cluster : main_clusters) {
            bool is_fc = false;
            try {
                // Conservative on failure: cluster_fc_check itself returns
                // is_fc=false when the cluster has no steiner_pc or when the
                // grouping carries no FiducialUtils (i.e. when the pipeline is
                // missing the steiner / fiducialutils stages).
                is_fc = Facade::cluster_fc_check(
                    *main_cluster, m_dv,
                    m_use_fiducial ? m_fiducial : nullptr,
                    m_fv_tolerance).is_fc;
            }
            catch (const std::exception& err) {
                SPDLOG_LOGGER_WARN(f_log, "visit: TaggerCheckFC: cluster {} check failed: {}",
                                   main_cluster->ident(), err.what());
            }
            if (is_fc) main_cluster->set_flag(Flags::FC);
            SPDLOG_LOGGER_INFO(f_log, "{}visit: TaggerCheckFC: cluster {} → FC={}", m_evt_tag,
                               main_cluster->ident(), is_fc);
        }
    }

private:
    std::string m_grouping_name{"live"};
    bool m_require_in_scope{false};
    bool m_evaluate_demoted_mains{false};
    // Off by default: with no "fiducial" key the containment tests stay on
    // FiducialUtils, i.e. identical to what TaggerCheckSTM/TaggerCheckNeutrino
    // see.  SBND sets both keys so FC and TGM share one fiducial (doc 27).
    bool m_use_fiducial{false};
    std::vector<double> m_fv_tolerance;
    // Beam-window gate on cluster_t0 (see configure()).  Off by default.
    bool m_beam_window_only{false};
    double m_beam_window_low{0};
    double m_beam_window_high{0};
};
