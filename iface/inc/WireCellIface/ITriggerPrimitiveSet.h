/** Trigger primitives (TPs) and sets of them.

    A TP summarizes one run of samples over threshold on one channel, as made
    by the DUNE DAQ trigger primitive generation.  TPs are many and small, so
    they are values held by their set, not IData.

    Times are absolute, in WCT units, on the same clock as IFrame::time().
    Durations are counts of samples, as in the DAQ TP format; the sample
    period is a configuration of whatever needs it, not part of the data.

    The stream contract: a set covers the half-open data time interval
    [time_begin, time_end).  Successive sets on a port cover contiguous
    intervals and a gap means the stream restarted.  The end of the
    coverage is a watermark: no later set on the port will hold a TP starting
    earlier than time_end, less the producer's declared latency (a TP still
    open at the end of a set is emitted with a later set).  A set may be
    empty, which still advances the watermark.  See trig/docs/design.org.
 */

#ifndef WIRECELL_ITRIGGERPRIMITIVESET
#define WIRECELL_ITRIGGERPRIMITIVESET

#include "WireCellIface/IData.h"

#include <cstdint>
#include <vector>

namespace WireCell {

    struct TriggerPrimitive {
        int channel{-1};                     // WCT channel ident
        double time_start{0};                // time of the first sample over threshold
        uint16_t samples_over_threshold{0};
        uint16_t samples_to_peak{0};         // from the first sample to the peak
        uint32_t adc_integral{0};
        uint16_t adc_peak{0};
        uint8_t flag{0};
    };

    class ITriggerPrimitiveSet : public IData<ITriggerPrimitiveSet> {
       public:
        typedef std::vector<TriggerPrimitive> primitives_t;

        virtual ~ITriggerPrimitiveSet();

        /// The sequence number of this set in its stream.
        virtual int ident() const = 0;

        /// The detector unit (anode ident) the TPs come from.
        virtual int detector() const = 0;

        /// The data time interval [begin, end) this set covers.
        virtual double time_begin() const = 0;
        virtual double time_end() const = 0;

        /// The TPs, ordered by time_start then channel.
        virtual const primitives_t& primitives() const = 0;
    };

}  // namespace WireCell

#endif
