#ifndef WIRECELL_SPNGDECON
#define WIRECELL_SPNGDECON

#include "WireCellAux/Logger.h"

#include "WireCellSpng/ITorchTensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITorchSpectrum.h"



namespace WireCell {
namespace SPNG {
    class Decon : public Aux::Logger,
                    public WireCell::ITorchTensorSetFilter, public WireCell::IConfigurable {
    public:
        Decon( );
        virtual ~Decon();

        virtual bool operator()(const input_pointer& in, output_pointer& out);
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const {
            Configuration cfg;
            return cfg;
        };
    private:
        std::string m_frer_spectrum{"FRERSpectrum"};
        std::string m_wire_filter{"Torch1DSpectrum"};
        std::shared_ptr<ITorchSpectrum> base_frer_spectrum, base_wire_filter;
        int m_coarse_time_offset = 0;
        bool m_unsqueeze_input = false;
        bool m_debug_no_frer = false;
        bool m_debug_no_wire_filter = false;
        bool m_debug_no_roll = false;
        bool m_use_fft_best_length = false;
        bool m_debug_force_cpu = false;
        bool m_pad_wire_domain = false;
        
        Json::Value m_output_set_tag{"Decon2D"}, m_output_tensor_tag{"Default"};
        Json::Value m_passthrough{Json::arrayValue};
    };
}
}

#endif