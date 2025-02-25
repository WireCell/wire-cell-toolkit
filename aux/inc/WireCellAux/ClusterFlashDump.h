/**

   ClusterFlashDump gives a "sink" which logs a summary of a PC tree holding
   cluster and flash information provided in the form of a tensor-data-model
   ITensorSet.

   It can be used as an example for how to utilize this data.

 */
#ifndef WIRECELLAUX_CLUSTERFLASHDUMP
#define WIRECELLAUX_CLUSTERFLASHDUMP

#include "WireCellIface/ITensorSetSink.h"
#include "WireCellAux/Logger.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::Aux {

    class ClusterFlashDump : public Aux::Logger,
                             public ITensorSetSink,
                             public IConfigurable {
    public:
        ClusterFlashDump();
        virtual ~ClusterFlashDump();
        virtual bool operator()(const ITensorSet::pointer& in);

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;
        /** Configuration: "datapath" (optional)

            Give a tensor data model path for the resulting ITensorSet
            representation.

            This must be provided.  Default is empty.

            Example, UbooneClusterSource default places data at datapath:

                "pointtrees/%d/uboone";

        */
        std::string m_datapath = "";

    };

}


#endif
