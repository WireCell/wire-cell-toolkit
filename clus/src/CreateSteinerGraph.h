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

#include "WireCellAux/Logger.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"


namespace WireCell::Clus::Steiner {

    class CreateSteinerGraph : public Aux::Logger, public IConfigurable, public Clus::IEnsembleVisitor {
        std::string m_grouping_name{"live"};
        std::string m_graph_name{"steiner"};
        bool m_replace{true};

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
