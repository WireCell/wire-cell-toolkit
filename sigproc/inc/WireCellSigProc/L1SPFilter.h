/** This component applies "compressed sensing" influenced signal
 * processing based on an L1 norm minimzation which fits both a
 * unipolar collection and a bipolar induction response to regions
 * channels in shorted regions known to have a mix.
 */
#ifndef WIRECELLSIGPROC_L1SPFILTER
#define WIRECELLSIGPROC_L1SPFILTER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDFT.h"

#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/Interpolate.h"

#include <vector>


namespace WireCell {
    namespace SigProc {

        class L1SPFilter : public Aux::Logger, public IFrameFilter, public IConfigurable {
           public:
            L1SPFilter(double gain = 14.0 * units::mV / units::fC, double shaping = 2.2 * units::microsecond,
                       double postgain = 1.2, double ADC_mV = 4096 / (2000. * units::mV),
                       double fine_time_offset = 0.0 * units::microsecond,
                       double coarse_time_offset = -8.0 * units::microsecond);
            virtual ~L1SPFilter();

            /// IFrameFilter interface.
            virtual bool operator()(const input_pointer& in, output_pointer& out);

            /// IConfigurable interface.
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

            void init_resp();

            int L1_fit(std::shared_ptr<WireCell::Aux::SimpleTrace>& newtrace,
                       std::shared_ptr<const WireCell::ITrace>& adctrace, int start_tick, int end_tick,
                       bool flag_shorted = false);

           private:
            Configuration m_cfg;
            IDFT::pointer m_dft;

            double m_gain;
            double m_shaping;
            double m_postgain;
            double m_ADC_mV;
            double m_fine_time_offset;
            double m_coarse_time_offset;
            double m_period;

            linterp<double>* lin_V;
            linterp<double>* lin_W;

            size_t m_count{0};

            // Cached L1_fit config values (avoid re-parsing JSON per call)
            double m_overall_time_offset{0};
            double m_collect_time_offset{3.0 * units::microsecond};
            double m_adc_l1_threshold{6};
            double m_adc_sum_threshold{160};
            double m_adc_sum_rescaling{90.0};
            double m_adc_sum_rescaling_limit{50.0};
            double m_adc_ratio_threshold{0.2};
            double m_l1_seg_length{120};
            double m_l1_scaling_factor{500};
            double m_l1_lambda{5};
            double m_l1_epsilon{0.05};
            double m_l1_niteration{100000};
            double m_l1_decon_limit{100};
            double m_l1_resp_scale{0.5};
            double m_l1_col_scale{1.15};
            double m_l1_ind_scale{0.5};
            double m_peak_threshold{10000};
            double m_mean_threshold{500};
            std::vector<double> m_smearing_vec;
        };
    }  // namespace SigProc
}  // namespace WireCell

#endif
