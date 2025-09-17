/// Reconstruct an IFrame from SPNG TDM compliant tensor set.

#ifndef WIRECELL_SPNG_TDMTOFRAME
#define WIRECELL_SPNG_TDMTOFRAME

#include "WireCellSpng/ITorchSetToFrame.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/Logger.h"

namespace WireCell::SPNG {

    /// A TdmToFrame will output an IFrame from an input TDM-compliant "frame"
    /// tensor collection.
    ///
    /// Caution: this class does not handled batched tensors!  Run an unbatcher
    /// node prior if you need to convert batched.
    class TdmToFrame: public ContextBase, public Logger, public ITorchSetToFrame  {
    public:

        TdmToFrame();
        virtual ~TdmToFrame();

        // IFunction
        virtual bool operator()(const input_pointer& in, output_pointer& out);
        
        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

    private:

        /// Configuration: "frame"
        ///
        /// Optional, how to find the frame in the input tensor set.
        ///
        /// - as integer :: return the i'th frame in the set.
        /// - as string :: use as regular expression to match datapath.
        ///
        /// default takes first frame in set.
        std::variant<int, std::string> m_find_frame;

    };
}


#endif

