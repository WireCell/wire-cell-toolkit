/** The TorchFRERSpectrum provides a detector response as a simple 2D tensor.

    The response is the convolution of field and electronics responses provided
    by their usual WCT components.

    This component is a BaseContext and will place tensors returned according to
    the "device" configuration parameter.

 */

#ifndef WIRECELLSPNG_TORCHFRERSPECTRUM
#define WIRECELLSPNG_TORCHFRERSPECTRUM
#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/ITorchSpectrum.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Units.h"
#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IWaveform.h"

// #include <boost/compute/detail/lru_cache.hpp>


namespace WireCell {
    namespace SPNG {
        class TorchFRERSpectrum : Logger, 
                                  ContextBase,
                                  public ITorchSpectrum,
                                  virtual public IConfigurable {
        public:
            // Create directly with the JSON data file or delay that
            // for configuration.
            TorchFRERSpectrum();

            virtual ~TorchFRERSpectrum();

            // ITorchSpectrum
            virtual torch::Tensor spectrum() const;
            virtual torch::Tensor spectrum(const std::vector<int64_t> & shape) const;
            virtual std::vector<int64_t> shape() const { return m_shape; }

            // IConfigurable
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;
            
            /// Get any shifts of the response
            virtual std::vector<int64_t> shifts() const;
        private:

            // void redigitize(const std::vector<int64_t> & input_shape);
            torch::Tensor m_total_response, m_redigitized_response;
            std::string m_field_response_name{"FieldResponse"};
            std::string m_elec_response_name{"ColdElecResponse"};

            std::shared_ptr<IWaveform> m_elec_response;
            Response::Schema::FieldResponse m_field_response,
                                            m_field_response_avg;
            
            //Relevant for Field Response
            int m_plane_id = 0;
            bool m_do_average = false;
            int64_t m_fravg_nticks = 0, m_fravg_nchans = 0;
            double m_fravg_period;
            double m_gain, m_ADC_mV;
            double m_readout_period;

            //Relevant for Cold Elec Response
            float m_extra_scale = 1.;
            int m_default_nticks = 0;
            int m_default_nchans = 0;
            bool m_do_fft = false;
            
            std::vector<int64_t> m_shape;

        };

    }  // namespace spng

}  // namespace WireCell
#endif
