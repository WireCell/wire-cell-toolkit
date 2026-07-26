/** This is a private header local to clus/src/.

    It defines an "ensemble visitor" that will add a "steiner graph" to certain
    clusters in a grouping.

    The Steiner::Grapher class does the work for each cluster.
 */

#ifndef WIRECELLCLUS_CREATESTEINERGRAPH
#define WIRECELLCLUS_CREATESTEINERGRAPH

#include "SteinerGrapher.h"

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Ensemble.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellAux/Logger.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"



namespace WireCell::Clus::Steiner {

    class CreateSteinerGraph : public Aux::Logger, public IConfigurable, public Clus::IEnsembleVisitor,
                               private Clus::NeedDV, private Clus::NeedPCTS{
        std::string m_grouping_name{"live"};
        std::string m_graph_name{"steiner"};
        bool m_replace{true};
        bool m_perf{false};
        // When true (default, uBooNE) only clusters flagged beam_flash are
        // processed.  Post-QL-matching detectors (e.g. SBND) have no
        // beam_flash flag: set false to process every scope-passing cluster,
        // with each flagged main cluster getting the main treatment.
        bool m_require_beam_flash{true};
        // Beam-window gate on the matched flash time (cluster_t0), the same
        // [low, high) test TaggerCheckNeutrino/TaggerCheckTGM use.  When
        // m_beam_window_only is false (default) or the window is empty
        // (low >= high) every cluster selected by the rule above is processed,
        // i.e. legacy behavior.  When on, only the beam-coincident bundle is
        // retiled/graphed: its in-window main clusters plus the companions
        // sharing their matched_flash_gid.  The downstream taggers
        // (TGM/STM/FC) carry the same gate, so the clusters that lose their
        // steiner graph are exactly the ones no tagger evaluates.
        bool m_beam_window_only{false};
        double m_beam_window_low{0};
        double m_beam_window_high{0};

        Grapher::Config m_grapher_config;


    public:
        CreateSteinerGraph();
        virtual ~CreateSteinerGraph();

        // IConfigurable 
        virtual void configure(const WireCell::Configuration& cfg);
        virtual Configuration default_configuration() const;

        // IEnsembleVisitor
        /// Loops over each cluster in the chosen grouping.
        /// See SteinerCluster per-cluster operations.
        virtual void visit(Facade::Ensemble& ensemble) const;

    private:

        /* Xin,

           No per-cluster stuff here.  See SteinerGrapher.h for that

        */



    };
}

#endif
