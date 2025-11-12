#ifndef WIRECELL_SPNG_FRAMETOTORCHSET
#define WIRECELL_SPNG_FRAMETOTORCHSET

#include "WireCellSpng/IFrameToTorchSet.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/Logger.h"


namespace WireCell {
    namespace SPNG {

        /**
           FrameToTorchSet converts an input IFrame into a TorchTensorSet.

           Configuration:

           - (none yet)

        */
        class FrameToTorchSet
            : public Aux::Logger,
              public IFrameToTorchSet, public IConfigurable {
           public:
            FrameToTorchSet();
            virtual ~FrameToTorchSet() {};

            // INode, override because we get multiplicity at run time.
            // virtual std::vector<std::string> output_types();

            // IFrameToTorchSet
            virtual bool operator()(const input_pointer& in, output_pointer& out);

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

           private:
                std::vector<std::string> m_intags;

        };
    }  // namespace SPNG
}  // namespace WireCell

#endif