#ifndef WIRECELL_SPNGTORCHTENSORSETTAGGER
#define WIRECELL_SPNGTORCHTENSORSETTAGGER

#include "WireCellAux/Logger.h"

#include "WireCellSpng/ITorchTensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITorchSpectrum.h"



namespace WireCell {
namespace SPNG {
    class TorchTensorSetTagger : public Aux::Logger,
                    public WireCell::ITorchTensorSetFilter, public WireCell::IConfigurable {
    public:
        TorchTensorSetTagger( );
        virtual ~TorchTensorSetTagger();

        virtual bool operator()(const input_pointer& in, output_pointer& out);
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const {
            Configuration cfg;
            return cfg;
        };
    private:
        //Normally, this node will crash if there already is a tag in the metadata
        bool m_allow_retagging{false};

        std::map<std::string, std::string> m_tag_list;
    };
}
}

#endif