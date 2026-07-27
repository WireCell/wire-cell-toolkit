#ifndef WIRECELL_SPNGAPPLYROI
#define WIRECELL_SPNGAPPLYROI

#include "WireCellAux/Logger.h"

#include "WireCellSpng/ITorchTensorSetFanin.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITorchSpectrum.h"



namespace WireCell {
namespace SPNG {
    class ApplyROI : public Aux::Logger,
                     public ITorchTensorSetFanin,
                     public IConfigurable {
    public:
        ApplyROI( );
        virtual ~ApplyROI();

        virtual bool operator()(const input_vector& inv, output_pointer& out);
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const {
            Configuration cfg;
            return cfg;
        };

        // INode, override because we get multiplicity at run time.
        virtual std::vector<std::string> input_types();

    private:
        // int m_dimension{0};
        int m_ROI_tensor_index{0},
                    m_value_tensor_index{1};
        const int m_multiplicity{2};
        Json::Value m_output_set_tag,
                    m_output_tensor_tag{"Default"};
    };
}
}

#endif
