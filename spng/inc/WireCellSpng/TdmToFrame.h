/// Reconstruct an IFrame from SPNG TDM compliant tensor set.

#ifndef WIRECELL_SPNG_TDMTOFRAME
#define WIRECELL_SPNG_TDMTOFRAME

#include "WireCellSpng/ITorchSetToFrame.h"
#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h" // may not be needed, we must safely call .to(cpu).
#include "WireCellUtil/HanaJsonCPP.h"

namespace WireCell::SPNG {


    /**
       @brief The configuration objects for TdmToFrame:

       - TdmToFrameConfig is the overall configuration.
       - TdmToFrameTraces is config for each set of tagged traces to produce.

       The configuration is in terms of "matching objects" used to select
       tensors from the set.  See functions and documentation in TdmTools.h for
       details.

       Select tensors to be combined into one set of tagged traces of the given tag.

       The traces, chids and summaries may each be a match object or array of
       match objects.  See TdmTools.h for "match objects".

       The tensors of each datatype will be concatenated along the channel
       dimension and thus the match must produced properly aligned arrays of
       matched tensors.

       If a "datatype" attribute is not provided it will suitably be added.
     */
    struct TdmToFrameTraces {
        /// The tag to apply to the tagged traces
        std::string tag;

        /// Match the traces tensors, required.
        Configuration traces;

        /// Match the chids tensors, required.
        Configuration chids;

        /// Match the summaries tensors, optional.
        Configuration summaries;
    };

    struct TdmToFrameConfig {
        /// A matching object (not array) to find the frame tensor.  If
        /// "datatype" attribute is missing it will be set to "frame".
        Configuration frame;

        /// Rules to form tagged traces.
        std::vector<TdmToFrameTraces> tagged_traces;

        /// A match object or array of them to match chmasks datatype.  If
        /// "datatype" is not provided, it will be set to "chmasks".
        Configuration chmasks;
        
    };

}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::TdmToFrameConfig, frame, tagged_traces, chmasks);
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::TdmToFrameTraces, tag, traces, chids, summaries);


namespace WireCell::SPNG {

    /// A TdmToFrame will output an IFrame from a TDM-compliant input "frame"
    /// tensor collection.
    ///
    /// This will copy input tensors to the CPU.
    ///
    /// Caution: this class does not handled batched tensors!  Run an unbatcher
    /// node prior if you need to convert batched.
    ///
    class TdmToFrame: public Logger,
                      public ContextBase,
                      public ITorchSetToFrame  {
    public:

        TdmToFrame();
        TdmToFrame(const TdmToFrameConfig& cfg);
        virtual ~TdmToFrame();

        // IFunction
        virtual bool operator()(const input_pointer& in, output_pointer& out);
        
        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        // Non-API methods

        // Common configuration processing
        void configme();

        // Helpers to pack from tensor to IFrame/ITrace forms.
        Waveform::ChannelMaskMap get_cmm(const ITorchTensor::vector& chmask_itensors) const;
        std::vector<int> get_chids(const ITorchTensor::vector& chids_itensors) const;
        ITrace::vector get_traces(const ITorchTensor::vector& traces_itensors,
                                  const std::vector<int>& chids) const;
        std::vector<double> get_summaries(const ITorchTensor::vector& summaries_itensors) const;

    private:

        TdmToFrameConfig m_cfg;
    };
}


#endif

