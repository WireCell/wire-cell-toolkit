#ifndef WIRECELL_SPNG_TORCHTENSORSETREPLICATOR
#define WIRECELL_SPNG_TORCHTENSORSETREPLICATOR

#include "WireCellSpng/ITorchTensorSetFanout.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/WirePlaneId.h"

namespace WireCell {
    namespace SPNG {

        // Fan out 1 frame to N set at construction or configuration time.
        class TorchTensorSetReplicator
            : public Aux::Logger,
              public ITorchTensorSetFanout, public IConfigurable {
           public:
            TorchTensorSetReplicator();
            virtual ~TorchTensorSetReplicator() {};

            // INode, override because we get multiplicity at run time.
            virtual std::vector<std::string> output_types();

            // IFanout
            virtual bool operator()(const input_pointer& in, output_vector& outv);

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

           private:
            int m_multiplicity{1};

            WireCell::Configuration m_cfg;
            // Log::logptr_t log;
        };
    }  // namespace Aux
}  // namespace WireCell

#endif