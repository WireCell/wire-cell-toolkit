/** Resample frames to a new sampling period/frequency.

    This applies the LMN method to resample a frame under the interpretation
    that the original sampling is that of an instantaneous quantity.  That is,
    the resampling applies interpolation normalization.  See the LMN paper for
    details.

    The resampling may either target a resampled period or a resampled size
    (number of samples/ticks).  When the target is a period, it must satisfy the
    LMN rationality condition.  Whether the target is a sampling period or the
    number samples, the other quantity is determined by the method.

    CAVEAT: this will resample all traces assuming they are dense and with not
    tbin offset.  Though, traces need not be all the same size.

    Any frame tags, trace tags or trace summaries are carried forward to the
    output as-given. 
 */

#ifndef WIRECELLAUX_RESAMPLER
#define WIRECELLAUX_RESAMPLER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"
#include "WireCellIface/IDFT.h"

namespace WireCell::Aux {

    class Resampler : public Aux::Logger,
                      public IFrameFilter, public IConfigurable {
      public:
        Resampler();
        virtual ~Resampler();

        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_pointer& inframe, output_pointer& outframe);

      private:

        // Configure: period
        //
        // The target sampling period.  If zero then the nticks is the target.
        // See nticks for more.
        double m_period{0};

        // Configure: nticks
        //
        // The target number of samples.  If zero then the period is the target.
        // If both nticks and period are nonzero then the period is the target and
        // the output waveform will be truncated or zero-padded (see fixme
        // below).
        int m_nticks{0};

        // Configure: dft
        //
        // Name of the DFT component
        IDFT::pointer m_dft;

        // Fixme: a "padding" config may be added in the future to provide the
        // name of time and/or frequency domain padding.  In time we may pad
        // with zeros, extend last sample or join last and first sample values
        // linearly, half-cosine or other means.  When an upsampling is needed,
        // in frequency we may zero pad or extend last sample value or even
        // extrapolate some model fit to the spectrum tail.


        size_t m_count{0};
        
    };

}
#endif
