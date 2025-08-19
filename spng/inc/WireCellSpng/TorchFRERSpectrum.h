/** This component provides field response data as read in from a "WCT
 * field response" JSON file */

#ifndef WIRECELLSPNG_TORCHFRERSPECTRUM
#define WIRECELLSPNG_TORCHFRERSPECTRUM
#include "WireCellAux/Logger.h"
#include "WireCellSpng/ITorchSpectrum.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Units.h"
#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IWaveform.h"

namespace WireCell {
    namespace SPNG {
        class TorchFRERSpectrum : public Aux::Logger, 
                                public ITorchSpectrum,
                                public IConfigurable {
        public:
            // Create directly with the JSON data file or delay that
            // for configuration.
            TorchFRERSpectrum();

            virtual ~TorchFRERSpectrum();

            // ITorchSpectrum
            virtual torch::Tensor spectrum() const;
            virtual torch::Tensor spectrum(const std::vector<int64_t> & shape);
        //  virtual std::vector<int64_t> shape() const;

            // IConfigurable
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;
            
            /// Get any shifts of the response
            virtual std::vector<int64_t> shifts() const;
        private:

            void redigitize(const std::vector<int64_t> & input_shape);
            torch::Tensor m_total_response;
            boost::compute::detail::lru_cache<std::vector<int64_t>, torch::Tensor> m_cache;
            std::string m_field_response_name{"FieldResponse"};
            std::string m_elec_response_name{"ColdElecResponse"};

            std::shared_ptr<IWaveform> m_elec_response;
            Response::Schema::FieldResponse m_field_response,
                                            m_field_response_avg;
            
            bool m_debug_force_cpu = false;
            //Relevant for Field Response
            int m_plane_id = 0;
            bool m_do_average = false;
            int64_t m_fravg_nticks = 0, m_fravg_nchans = 0;
            double m_fravg_period;
            double m_inter_gain, m_ADC_mV;
            double m_default_period;

            //Relevant for Cold Elec Response
            float m_extra_scale = 1.;
            float m_tick_period = 0.;
            double m_gain = 0.;
            double m_shaping = 0.;
            int m_default_nticks = 0;
            int m_default_nchans = 0;
            bool m_do_fft = false;
            
            int m_anode_num = 0;

        };

    }  // namespace spng

}  // namespace WireCell
#endif