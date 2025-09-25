#ifndef WIRECELL_SPNGDECON
#define WIRECELL_SPNGDECON

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"

#include "WireCellSpng/ITorchTensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITorchSpectrum.h"



namespace WireCell {
namespace SPNG {
    class Decon : public ContextBase,
                  public Logger,
                  public ITorchTensorSetFilter,
                  virtual public IConfigurable {
    public:
        Decon( );
        virtual ~Decon();

        virtual bool operator()(const input_pointer& in, output_pointer& out);
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

    protected:

        // The actual deconvolution.  The input tensor provides ADC waveforms
        // and is assumed to be batched.  The sampling period must be provided.
        virtual torch::Tensor decon(torch::Tensor waveforms, double period) const;

        // Helper to get the input itensor
        ITorchTensor::pointer get_input(const ITorchTensorSet::pointer& in, size_t index=0) const;
        
        // Helper to prepare the output tensor.
        virtual ITorchTensorSet::pointer make_output(const ITorchTensorSet::pointer& intenset,
                                                     const ITorchTensor::pointer& inten,
                                                     torch::Tensor outten) const;

    private:

        std::string m_frer_spectrum{"FRERSpectrum"};
        std::string m_wire_filter{"Torch1DSpectrum"};
        std::shared_ptr<ITorchSpectrum> m_base_frer_spectrum, m_base_wire_filter;
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

        /// Configuration: tensor_index (default = 0)
        ///
        /// The index for the input tensor in the tensor set on which to operate.
        int m_tensor_index{0};
    };

}
}

#endif
