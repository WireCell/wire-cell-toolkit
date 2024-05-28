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
        };
    }  // namespace SigProc
}  // namespace WireCell

#endif
