#ifndef WIRECELL_SPNG_TORCHTENSORSETUNSTACKER
#define WIRECELL_SPNG_TORCHTENSORSETUNSTACKER

#include "WireCellSpng/ITorchTensorSetFanout.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace SPNG {

        // Fan out 1 frame to N set at construction or configuration time.
        class TorchTensorSetUnstacker
            : public Aux::Logger,
              public ITorchTensorSetFanout, public IConfigurable {
        public:
            TorchTensorSetUnstacker();
            virtual ~TorchTensorSetUnstacker() {};

            // INode, override because we get multiplicity at run time.
            virtual std::vector<std::string> output_types();

            // IFanout
            virtual bool operator()(const input_pointer& in, output_vector& outv);

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

        private:
            int m_multiplicity{1};
        };
    }  // namespace Aux
}  // namespace WireCell

#endif