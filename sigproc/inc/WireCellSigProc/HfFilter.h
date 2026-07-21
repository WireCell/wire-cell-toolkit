
#ifndef WIRECELLSIGPROC_HFFILTER
#define WIRECELLSIGPROC_HFFILTER

#include "WireCellIface/IFilterWaveform.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Response.h"

namespace WireCell {
    namespace SigProc {
        class HfFilter : public IFilterWaveform, public IConfigurable {
           public:
            HfFilter(double max_freq = 1 * units::megahertz, double sigma = 3.0 * units::megahertz, double power = 2,
                     bool flag = true);
            virtual ~HfFilter();

            virtual const Waveform::realseq_t filter_waveform(int nfbins) const;

            // IConfigurable
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

           private:
            double m_max_freq;
            double m_sigma;
            double m_power;
            bool m_flag;
            bool m_use_negative_freqs{true};
        };
    }  // namespace SigProc
}  // namespace WireCell

#endif
