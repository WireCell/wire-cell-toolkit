#include "WireCellTrig/DaqConvert.h"

#include "WireCellUtil/Exceptions.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

using namespace WireCell;
using namespace WireCell::Trig;

Daq::timestamp_t DaqClock::to_ticks(double time) const
{
    const double ticks = (time - epoch) / tick;
    if (ticks < 0) {
        THROW(ValueError() << errmsg{"trig: time " + std::to_string(time) + " is before the DAQ clock epoch"});
    }
    // A small allowance so times on the tick grid do not round down a tick.
    return (Daq::timestamp_t) std::floor(ticks + 1e-6);
}

double DaqClock::to_time(Daq::timestamp_t ticks) const
{
    if (ticks == Daq::invalid_timestamp) return std::numeric_limits<double>::quiet_NaN();
    return epoch + ticks * tick;
}

static int from_daq_channel(Daq::channel_t ch)
{
    if (ch == Daq::invalid_channel || ch == Daq::invalid_tp_channel) return -1;
    return (int) ch;
}

static Daq::channel_t to_daq_channel(int ch)
{
    return ch < 0 ? Daq::invalid_channel : (Daq::channel_t) ch;
}

static Daq::detid_t to_daq_detid(int detector)
{
    if (detector < 0 || detector >= (int) Daq::whole_detector) return Daq::invalid_detid;
    return (Daq::detid_t) detector;
}

static int from_daq_detid(Daq::detid_t detid)
{
    if (detid == Daq::invalid_detid || detid == Daq::whole_detector) return -1;
    return (int) detid;
}

Daq::TriggerPrimitive WireCell::Trig::to_daq(const TriggerPrimitive& tp, const DaqClock& clock, double)
{
    Daq::TriggerPrimitive dtp;
    dtp.flag = tp.flag;
    dtp.detid = 0;
    dtp.channel = (Daq::channel_t) tp.channel;
    dtp.samples_over_threshold = tp.samples_over_threshold;
    dtp.time_start = clock.to_ticks(tp.time_start);
    dtp.samples_to_peak = tp.samples_to_peak;
    dtp.adc_integral = tp.adc_integral;
    dtp.adc_peak = tp.adc_peak;
    return dtp;
}

TriggerPrimitive WireCell::Trig::from_daq(const Daq::TriggerPrimitive& dtp, const DaqClock& clock)
{
    TriggerPrimitive tp;
    tp.channel = from_daq_channel(dtp.channel);
    tp.time_start = clock.to_time(dtp.time_start);
    tp.samples_over_threshold = dtp.samples_over_threshold;
    tp.samples_to_peak = dtp.samples_to_peak;
    tp.adc_integral = dtp.adc_integral;
    tp.adc_peak = dtp.adc_peak;
    tp.flag = dtp.flag;
    return tp;
}

Daq::TriggerActivity WireCell::Trig::to_daq(const TriggerActivity& ta, const DaqClock& clock, double period)
{
    Daq::TriggerActivity dta;
    dta.time_start = clock.to_ticks(ta.time_start);
    dta.time_end = clock.to_ticks(ta.time_end);
    dta.time_peak = clock.to_ticks(ta.time_peak);
    dta.time_activity = clock.to_ticks(ta.time_activity);
    dta.channel_start = to_daq_channel(ta.channel_start);
    dta.channel_end = to_daq_channel(ta.channel_end);
    dta.channel_peak = to_daq_channel(ta.channel_peak);
    dta.adc_integral = ta.adc_integral;
    dta.adc_peak = ta.adc_peak;
    dta.detid = to_daq_detid(ta.detector);
    dta.type = Daq::TriggerActivityData::Type::kTPC;
    dta.algorithm = (Daq::TriggerActivityData::Algorithm) ta.algorithm;
    if (ta.inputs.empty()) {
        dta.inputs.resize(ta.ninputs);
    }
    else {
        dta.inputs.reserve(ta.inputs.size());
        for (const auto& tp : ta.inputs) dta.inputs.push_back(to_daq(tp, clock, period));
    }
    return dta;
}

TriggerActivity WireCell::Trig::from_daq(const Daq::TriggerActivityData& dta, const DaqClock& clock)
{
    TriggerActivity ta;
    ta.time_start = clock.to_time(dta.time_start);
    ta.time_end = clock.to_time(dta.time_end);
    ta.time_peak = clock.to_time(dta.time_peak);
    ta.time_activity = clock.to_time(dta.time_activity);
    ta.channel_start = from_daq_channel(dta.channel_start);
    ta.channel_end = from_daq_channel(dta.channel_end);
    ta.channel_peak = from_daq_channel(dta.channel_peak);
    ta.adc_integral = dta.adc_integral;
    ta.adc_peak = dta.adc_peak;
    ta.detector = from_daq_detid(dta.detid);
    ta.algorithm = (int) dta.algorithm;
    return ta;
}

TriggerActivity WireCell::Trig::from_daq(const Daq::TriggerActivity& dta, const DaqClock& clock)
{
    TriggerActivity ta = from_daq(static_cast<const Daq::TriggerActivityData&>(dta), clock);
    ta.ninputs = dta.inputs.size();
    ta.inputs.reserve(dta.inputs.size());
    for (const auto& dtp : dta.inputs) ta.inputs.push_back(from_daq(dtp, clock));
    return ta;
}

TriggerCandidate WireCell::Trig::from_daq(const Daq::TriggerCandidate& dtc, const DaqClock& clock)
{
    TriggerCandidate tc;
    tc.time_start = clock.to_time(dtc.time_start);
    tc.time_end = clock.to_time(dtc.time_end);
    tc.time_candidate = clock.to_time(dtc.time_candidate);
    tc.type = (int) dtc.type;
    tc.algorithm = (int) dtc.algorithm;
    tc.ninputs = dtc.inputs.size();
    tc.inputs.reserve(dtc.inputs.size());
    for (const auto& dta : dtc.inputs) {
        tc.inputs.push_back(from_daq(dta, clock));
        if (tc.inputs.back().detector >= 0) tc.detectors.push_back(tc.inputs.back().detector);
    }
    std::sort(tc.detectors.begin(), tc.detectors.end());
    tc.detectors.erase(std::unique(tc.detectors.begin(), tc.detectors.end()), tc.detectors.end());
    return tc;
}
