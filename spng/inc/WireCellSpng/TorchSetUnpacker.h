#ifndef WIRECELL_SPNG_TORCHSETUNPACKER
#define WIRECELL_SPNG_TORCHSETUNPACKER

#include "WireCellSpng/ITorchSetUnpacker.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/Logger.h"
#include "WireCellUtil/HanaJsonCPP.h"

namespace WireCell::SPNG {

    struct TorchSetUnpackerSelection {
        /// A tensor is selected either by an integer index into the tensor set or
        /// by a string providing a regular expression to match against datapaths of
        /// TDM compliant tensor sets.  At least one must be set to a valid value.
        int index=-1;           // invalid
        std::string datapath=""; // invalid

        std::string str() const;
    };
    struct TorchSetUnpackerConfig {

        /// A list of tensor selections, one for each output port.
        std::vector<TorchSetUnpackerSelection> selections;
    };
}
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::TorchSetUnpackerSelection, index, datapath);
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::TorchSetUnpackerConfig, selections);

namespace WireCell::SPNG {

    /// Select ITorchTensor instances from an ITorchTensorSet set and fan out
    /// the selection to per-tensor edges.
    ///
    /// This does not modify the ITorchTensor instances.
    class TorchSetUnpacker : public Logger, public ITorchSetUnpacker, virtual public IConfigurable {
    public:

        TorchSetUnpacker();
        TorchSetUnpacker(const TorchSetUnpackerConfig& cfg);
        virtual ~TorchSetUnpacker();

        // INode, override because we get multiplicity at run time.
        virtual std::vector<std::string> output_types();

        // IFanoutNode<ITorchTensorSet, ITorchTensor, 0>
        virtual bool operator()(const input_pointer& in, output_vector& outv);

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;


        // Non API methods
        void configme();
        size_t multiplicity() const { return m_selectors.size(); }

    private:
        TorchSetUnpackerConfig m_cfg;
        using selector_func = std::function<ITorchTensor::pointer(ITorchTensorSet::pointer)>;
        std::vector<selector_func> m_selectors;
    };

}

#endif
