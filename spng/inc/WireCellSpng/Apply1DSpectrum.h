#ifndef WIRECELL_SPNGAPPLY1DSPECTRUM
#define WIRECELL_SPNGAPPLY1DSPECTRUM

#include "WireCellAux/Logger.h"

#include "WireCellSpng/ITorchTensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITorchSpectrum.h"



namespace WireCell {
namespace SPNG {
    class Apply1DSpectrum : public Aux::Logger,
                    public WireCell::ITorchTensorSetFilter, public WireCell::IConfigurable {
    public:
        Apply1DSpectrum( );
        virtual ~Apply1DSpectrum();

        virtual bool operator()(const input_pointer& in, output_pointer& out);
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const {
            Configuration cfg;
            return cfg;
        };
    private:
        std::string m_base_spectrum_name{"Torch1DSpectrum"};
        std::shared_ptr<ITorchSpectrum> m_base_spectrum;
        int m_dimension{0};
        Json::Value m_target_tensor{"Default"},
                    m_output_set_tag,
                    m_output_tensor_tag{"Default"};
        Json::Value m_passthrough{Json::arrayValue};
    };
}
}

#endif