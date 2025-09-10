#ifndef WIRECELL_SPNG_FRAMETOTORCH
#define WIRECELL_SPNG_FRAMETOTORCH

#include "WireCellSpng/IFrameToTorch.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/TagRules.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace SPNG {

        /// Convert 1 frame to 1 TorchTensor
        /// TODO -- Flesh out description
        class FrameToTorch : public Aux::Logger,
                              public IFrameToTorch, public IConfigurable {
           public:
            FrameToTorch()
             : Aux::Logger("FrameToTorch", "spng") {};
            virtual ~FrameToTorch() {};

            // INode, override because we get multiplicity at run time.
            // virtual std::vector<std::string> output_types() {return };

            virtual bool operator()(const input_pointer& in, output_pointer& out);

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg) {};
            virtual WireCell::Configuration default_configuration() const;

           private:

        };
    }  // namespace SPNG
}  // namespace WireCell

#endif
