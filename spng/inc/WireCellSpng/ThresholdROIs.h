#ifndef WIRECELL_SPNGDECON
#define WIRECELL_SPNGDECON

#include "WireCellAux/Logger.h"

#include "WireCellSpng/ITorchTensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITorchSpectrum.h"



namespace WireCell {
namespace SPNG {
    class ThresholdROIs : public Aux::Logger,
                    public WireCell::ITorchTensorSetFilter, public WireCell::IConfigurable {
    public:
        ThresholdROIs( );
        virtual ~ThresholdROIs();

        virtual bool operator()(const input_pointer& in, output_pointer& out);
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const {
            Configuration cfg;
            return cfg;
        };
    private:
        bool m_debug_force_cpu{false};
        bool m_unsqueeze_input{false};
        double m_threshold_rms_factor{1.};
        Json::Value m_output_set_tag{"ThresholdROIs"}, m_output_tensor_tag{"Default"};
        Json::Value m_passthrough{Json::arrayValue};
    };
}
}

#endif