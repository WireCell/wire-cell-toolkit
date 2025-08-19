#ifndef WIRECELL_SPNG_TORCHSETUNPACKER
#define WIRECELL_SPNG_TORCHSETUNPACKER

#include "WireCellSpng/ITorchSetUnpacker.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"

namespace WireCell {
    namespace SPNG {

        // Fan out 1 frame to N set at construction or configuration time.
        class TorchSetUnpacker : public ITorchSetUnpacker, public IConfigurable {
           public:
            TorchSetUnpacker(size_t multiplicity = 0);
            virtual ~TorchSetUnpacker();

            // INode, override because we get multiplicity at run time.
            virtual std::vector<std::string> output_types();

            // IFanout
            virtual bool operator()(const input_pointer& in, output_vector& outv);

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

           private:
            size_t m_multiplicity;
            WireCell::Configuration m_cfg;
            Log::logptr_t log;
        };
    }  // namespace Aux
}  // namespace WireCell

#endif