#ifndef WIRECELL_GEN_CLUSTERFANOUT
#define WIRECELL_GEN_CLUSTERFANOUT

#include "WireCellIface/IClusterFanout.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell::Img {

    /// Fan out clusters 1 to N
    class ClusterFanout : public Aux::Logger,
                          public IClusterFanout, public IConfigurable {
      public:
        ClusterFanout(size_t multiplicity = 0);
        virtual ~ClusterFanout();

        // INode, override because we get multiplicity at run time.
        virtual std::vector<std::string> output_types();

        // IFanout
        virtual bool operator()(const input_pointer& in, output_vector& outv);

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

      private:
        size_t m_multiplicity{0};
        size_t m_count{0};

        bool m_trivial{false};

    };
}

#endif
