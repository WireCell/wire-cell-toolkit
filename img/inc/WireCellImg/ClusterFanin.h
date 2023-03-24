#ifndef WIRECELLIMG_CLUSTERFANIN
#define WIRECELLIMG_CLUSTERFANIN

#include "WireCellIface/IClusterFanin.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell::Img {

    /** Fan N input clusters to 1 output cluster.
        
        Each input graph becomes a connected subgraph component of the
        output graph.

        The relative vertex descriptive ordering within each input
        subgraph is retained on ouput.

        Other than edges modified to match the change in vertex
        descriptors, no other input graph information is modified
        in its output representation.

     */
    class ClusterFanin :  public Aux::Logger,
                          public IClusterFanin, public IConfigurable {
      public:
        ClusterFanin();
        virtual ~ClusterFanin();

        // INode, override because we get multiplicity at run time.
        virtual std::vector<std::string> input_types();

        // IFanin
        virtual bool operator()(const input_vector& inv, output_pointer& out);

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

      private:

        /** Config: multiplicity

            The number of inputs to the fan. */
        size_t m_multiplicity{0};

        /** Config: ident_source

            The produced ICluster::ident() may be set from the graph
            input from the port number given by this configuration.
        */
        size_t m_ident_source{0};

        size_t m_count{0};
    };

}

#endif
