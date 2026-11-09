#include "WireCellTrig/SupernovaMakers.h"

#include <algorithm>
#include <cstdlib>

using namespace WireCell::Trig;
using namespace WireCell::Trig::Daq;

TAMakerSupernova::TAMakerSupernova()
  : m_cfg()
{
}

TAMakerSupernova::TAMakerSupernova(const Config& cfg)
  : m_cfg(cfg)
{
}

timestamp_t TAMakerSupernova::time_end_of(const input_type& tp) const
{
    return tp.time_start + (timestamp_t) (tp.samples_over_threshold * m_cfg.ticks_per_sample);
}

timestamp_t TAMakerSupernova::time_peak_of(const input_type& tp) const
{
    return tp.time_start + (timestamp_t) (tp.samples_to_peak * m_cfg.ticks_per_sample);
}

bool TAMakerSupernova::is_time_consistent(const input_type& tp) const
{
    const auto tol = m_cfg.time_tolerance;
    const auto ts = (timestamp_diff_t) tp.time_start;
    const auto te = (timestamp_diff_t) time_end_of(tp);
    const auto cs = (timestamp_diff_t) m_ta.time_start;
    const auto ce = (timestamp_diff_t) m_ta.time_end;
    const bool close_to_edge =
        tol > std::abs(ts - ce) || tol > std::abs(ts - cs) || tol > std::abs(te - ce) || tol > std::abs(te - cs);
    const bool in_between = (te > cs && te < ce) || (ts > cs && ts < ce);
    return in_between || close_to_edge;
}

bool TAMakerSupernova::is_channel_consistent(const input_type& tp) const
{
    const auto tol = m_cfg.channel_tolerance;
    const auto ch = (channel_diff_t) tp.channel;
    const auto cs = (channel_diff_t) m_ta.channel_start;
    const auto ce = (channel_diff_t) m_ta.channel_end;
    const bool close_to_edge = tol > std::abs(ch - ce) || tol > std::abs(ch - cs);
    const bool in_between = ch > cs && ch < ce;
    return in_between || close_to_edge;
}

void TAMakerSupernova::start(const input_type& tp)
{
    m_ta = TriggerActivity{};
    m_ta.time_start = tp.time_start;
    m_ta.time_end = time_end_of(tp);
    m_ta.time_peak = time_peak_of(tp);
    m_ta.time_activity = 0;  // never set by triggeralgs
    m_ta.channel_start = tp.channel;
    m_ta.channel_end = tp.channel;
    m_ta.channel_peak = tp.channel;
    m_ta.adc_integral = tp.adc_integral;
    m_ta.adc_peak = tp.adc_peak;
    m_ta.detid = tp.detid;
    m_ta.type = TriggerActivityData::Type::kTPC;
    m_ta.algorithm = TriggerActivityData::Algorithm::kSupernova;
    m_ta.inputs.push_back(tp);
    m_active = true;
}

void TAMakerSupernova::process(const input_type& tp, output_vector& out)
{
    if (!m_active) {
        start(tp);
        return;
    }
    const bool time_ok = is_time_consistent(tp);
    const bool channel_ok = is_channel_consistent(tp);
    const bool join = m_cfg.both ? (time_ok && channel_ok) : (time_ok || channel_ok);
    if (!join) {
        out.push_back(m_ta);
        start(tp);
        return;
    }
    const timestamp_t tend = time_end_of(tp);
    if (tp.time_start < m_ta.time_start) m_ta.time_start = tp.time_start;
    if (tend > m_ta.time_end) m_ta.time_end = tend;
    if (tp.adc_peak > m_ta.adc_peak) {
        m_ta.time_peak = time_peak_of(tp);
        m_ta.adc_peak = tp.adc_peak;
        m_ta.channel_peak = tp.channel;
    }
    if (tp.channel > m_ta.channel_end) m_ta.channel_end = tp.channel;
    if (tp.channel < m_ta.channel_start) m_ta.channel_start = tp.channel;
    m_ta.inputs.push_back(tp);
    m_ta.adc_integral += tp.adc_integral;
    m_ta.detid |= tp.detid;
}

void TAMakerSupernova::do_flush(timestamp_t until, output_vector& out)
{
    if (!m_active) return;
    if (until != invalid_timestamp) {
        // With either-or merging a TP on a nearby channel may join at any
        // later time, so only the end of stream closes the cluster.  With
        // both, no TP starting at or after time_end + tolerance can join.
        if (!m_cfg.both) return;
        if ((timestamp_diff_t) until < (timestamp_diff_t) m_ta.time_end + m_cfg.time_tolerance) return;
    }
    out.push_back(m_ta);
    m_ta = TriggerActivity{};
    m_active = false;
}

TCMakerSupernova::TCMakerSupernova()
  : m_cfg()
{
}

TCMakerSupernova::TCMakerSupernova(const Config& cfg)
  : m_cfg(cfg)
{
}

void TCMakerSupernova::process(const input_type& ta, output_vector& out)
{
    const timestamp_t time = ta.time_start;

    // forget TAs that fell out of the window
    const timestamp_diff_t how_far = (timestamp_diff_t) time - m_cfg.time_window;
    auto end = std::remove_if(m_activity.begin(), m_activity.end(), [how_far](const auto& one) {
        return (timestamp_diff_t) one.time_start < how_far;
    });
    m_activity.erase(end, m_activity.end());

    if (ta.inputs.size() > m_cfg.hit_threshold) {
        m_activity.push_back(static_cast<const TriggerActivityData&>(ta));
    }
    if (m_activity.size() <= m_cfg.threshold) {
        return;
    }
    TriggerCandidate tc;
    tc.time_start = how_far > 0 ? (timestamp_t) how_far : 0;
    tc.time_end = ta.time_end;
    tc.time_candidate = time;
    tc.detid = whole_detector;
    tc.type = common().tc_type;
    tc.algorithm = TriggerCandidateData::Algorithm::kSupernova;
    tc.inputs = m_activity;
    m_activity.clear();
    out.push_back(tc);
}
