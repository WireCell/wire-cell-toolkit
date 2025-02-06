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
        // The target sampling period.
        double m_period{0};

        // Configure: dft
        //
        // Name of the DFT component
        IDFT::pointer m_dft;

        // Configure: time_pad
        //
        // Name for a strategy for padding in the time domain when needed to
        // reach rationality condition.
        //
        // - zero :: pad with zero values (default)
        // - first :: pad with first time sample value
        // - last :: pad with last time sample value
        // - median :: pad with median value
        // - linear :: pad with linear between first and last sample values
        std::string m_time_pad{"zero"};

        // Configure: time_sizing
        //
        // Name a strategy for truncating the final time domain waveform
        //
        // - duration :: approximately retain duration (deafult).
        // - padded :: include rationality condition padding, duration may change.
        std::string m_time_sizing{"duration"};

        size_t m_count{0};
        
    };

}
#endif
