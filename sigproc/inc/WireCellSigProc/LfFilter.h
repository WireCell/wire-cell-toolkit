
#ifndef WIRECELLSIGPROC_LFFILTER
#define WIRECELLSIGPROC_LFFILTER

#include "WireCellIface/IFilterWaveform.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Response.h"

namespace WireCell {
    namespace SigProc {
        class LfFilter : public IFilterWaveform, public IConfigurable {
           public:
            LfFilter(double max_freq = 1 * units::megahertz, double tau = 0.02 * units::megahertz);
            virtual ~LfFilter();

            virtual const Waveform::realseq_t filter_waveform(int nbins) const;

            // IConfigurable
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

           private:
            double m_max_freq;
            double m_tau;
            bool m_use_negative_freqs{true};

        };
    }  // namespace SigProc
}  // namespace WireCell

#endif
