/** Conversions between the WCT trigger objects and the DAQ mirrors.

    Times convert with a DAQ clock: ticks = floor((t - epoch)/tick).  Times
    coming back are rounded down to the clock tick.  Channels convert as is.
    A WCT detector (anode ident) goes to the DAQ detid when it fits in it
    (0..253), else to invalid_detid.

    Losses: a DAQ TA has no plane, and a DAQ TC keeps only the data part of
    its TAs (no TPs).  A WCT TA with stripped inputs becomes a DAQ TA with
    "ninputs" default TPs, so that algorithms counting TPs still work.
 */

#ifndef WIRECELLTRIG_DAQCONVERT
#define WIRECELLTRIG_DAQCONVERT

#include "WireCellIface/ITriggerCandidateSet.h"
#include "WireCellTrig/DaqTypes.h"

namespace WireCell::Trig {

    /// Time conversions between WCT and the DAQ clock.
    struct DaqClock {
        double epoch{0};
        double tick{16};  // ns, as WCT units
        Daq::timestamp_t to_ticks(double time) const;
        double to_time(Daq::timestamp_t ticks) const;  // NaN for an invalid timestamp
    };

    Daq::TriggerPrimitive to_daq(const TriggerPrimitive& tp, const DaqClock& clock, double period);
    TriggerPrimitive from_daq(const Daq::TriggerPrimitive& tp, const DaqClock& clock);

    Daq::TriggerActivity to_daq(const TriggerActivity& ta, const DaqClock& clock, double period);
    TriggerActivity from_daq(const Daq::TriggerActivity& ta, const DaqClock& clock);
    /// The TA data alone, as held by a DAQ TC.  ninputs is zero.
    TriggerActivity from_daq(const Daq::TriggerActivityData& ta, const DaqClock& clock);

    /// The detectors are those of the inputs, sorted and unique.
    TriggerCandidate from_daq(const Daq::TriggerCandidate& tc, const DaqClock& clock);

}  // namespace WireCell::Trig

#endif
