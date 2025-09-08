#ifndef WIRECELL_SPNG_TORCHPACKER
#define WIRECELL_SPNG_TORCHPACKER

#include "WireCellSpng/ITorchPacker.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace SPNG {

        // Fan out 1 frame to N set at construction or configuration time.
        class TorchPacker : public Aux::Logger, public ITorchPacker, public IConfigurable {
           public:
            TorchPacker(size_t multiplicity = 0);
            virtual ~TorchPacker();

            // INode, override because we get multiplicity at run time.
            virtual std::vector<std::string> input_types();

            // IFanout
            virtual bool operator()(const input_vector& invec, output_pointer& out);

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

           private:
            /// Configuration: "multiplicity"
            ///
            /// Set the fan size.  Default is 2.
            size_t m_multiplicity{2};
            
            /// Configuration: "count"
            ///
            /// Set the initial "count" which is used to set the ident number of
            /// the output ITorchTensorSet.  Default is 0.
            size_t m_count{0};

        };
    }  // namespace SPNG
}  // namespace WireCell

#endif
