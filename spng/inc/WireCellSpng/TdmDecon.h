/// A TDM-compliant version of the original SPNG::Decon node.

#ifndef WIRECELL_SPNG_TDMDECON
#define WIRECELL_SPNG_TDMDECON

#include "WireCellSpng/TorchFunctionNode.h"



namespace WireCell {
namespace SPNG {

    /// A node applying filtered deconvolution to input tensors via FFT method.
    ///
    /// In Fourier domain this effectively outputs
    ///
    ///   M*F/R
    ///
    /// Where M is an input measurement tensor, R is a response kernel and F is
    /// a filter.
    ///
    class TdmDecon : public TorchFunctionNode {
    public:
        TdmDecon( );
        virtual ~TdmDecon();

        virtual TensorIndex transform_tensors(TensorIndex ti) const;

        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

    private:
        std::string m_frer_spectrum{"FRERSpectrum"};
        std::string m_wire_filter{"Torch1DSpectrum"};

        /// See FunctionNode for configuration related to specifying
        /// input/output datapaths.
        ///
        /// Configuration: "response"
        ///
        /// Type/name of an ITorchSpectrum providing deconvolution (response)
        /// kernel.
        ///
        /// Configuration: "filter"
        ///
        /// Type/name of an ITorchSpectrum providing deconvolution filter.
        std::shared_ptr<ITorchSpectrum> m_response, m_filter;


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
