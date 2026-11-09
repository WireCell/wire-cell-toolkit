#include "WireCellTrig/BurstMaker.h"

#include <algorithm>

using namespace WireCell;
using namespace WireCell::Trig;

bool TCMakerBurst::Selection::operator()(const TriggerActivity& ta) const
{
    if (!algorithms.empty() && std::find(algorithms.begin(), algorithms.end(), ta.algorithm) == algorithms.end()) {
        return false;
    }
    if (ta.ninputs < min_inputs || ta.ninputs > max_inputs) return false;
    if (ta.adc_integral < min_adc_integral || ta.adc_integral > max_adc_integral) return false;
    if (ta.channel_end - ta.channel_start > max_channel_span) return false;
    if (ta.time_end - ta.time_start > max_time_span) return false;
    return true;
}

TCMakerBurst::TCMakerBurst()
  : m_cfg()
{
}

TCMakerBurst::TCMakerBurst(const Config& cfg)
  : m_cfg(cfg)
{
}

void TCMakerBurst::reset()
{
    m_pending.clear();
    m_window.clear();
    m_processed = m_head = m_armed = -std::numeric_limits<double>::infinity();
}

void TCMakerBurst::operator()(const TriggerActivity& ta, output_vector& out)
{
    ++m_nseen;
    if (!m_cfg.select(ta)) return;
    ++m_nselected;
    if (ta.time_start < m_processed) {
        late(ta, out);
        return;
    }
    m_pending.emplace(ta.time_start, ta);
}

void TCMakerBurst::flush(double until, output_vector& out)
{
    auto it = m_pending.begin();
    while (it != m_pending.end() && it->first < until) {
        enter(it->second, out);
        it = m_pending.erase(it);
    }
    m_processed = std::max(m_processed, until);
}

// A TA in time order: it becomes the head of the window.
void TCMakerBurst::enter(const TriggerActivity& ta, output_vector& out)
{
    m_head = ta.time_start;
    m_window.push_back(ta);
    while (!m_window.empty() && m_window.front().time_start <= m_head - m_cfg.window) {
        m_window.pop_front();
    }
    check(m_head, out);
}

// A TA before what was processed: count it if it is in the window.
void TCMakerBurst::late(const TriggerActivity& ta, output_vector& out)
{
    ++m_nlate;
    if (ta.time_start >= m_head) {  // still in order
        enter(ta, out);
        return;
    }
    if (ta.time_start <= m_head - m_cfg.window) {
        ++m_ndropped;
        return;
    }
    auto pos = std::upper_bound(m_window.begin(), m_window.end(), ta.time_start,
                                [](double t, const TriggerActivity& one) { return t < one.time_start; });
    m_window.insert(pos, ta);
    check(m_head, out);
}

void TCMakerBurst::check(double now, output_vector& out)
{
    if (m_window.size() < m_cfg.threshold || now < m_armed) return;

    TriggerCandidate tc;
    tc.time_candidate = now;
    tc.time_start = m_window.front().time_start;
    tc.time_end = m_window.front().time_end;
    for (const auto& ta : m_window) {
        tc.time_end = std::max(tc.time_end, ta.time_end);
        if (ta.detector >= 0) tc.detectors.push_back(ta.detector);
    }
    std::sort(tc.detectors.begin(), tc.detectors.end());
    tc.detectors.erase(std::unique(tc.detectors.begin(), tc.detectors.end()), tc.detectors.end());
    tc.type = m_cfg.tc_type;
    tc.algorithm = m_cfg.algorithm;
    tc.ninputs = m_window.size();
    tc.inputs.assign(m_window.begin(), m_window.end());
    out.push_back(std::move(tc));

    m_armed = now + (m_cfg.holdoff < 0 ? m_cfg.window : m_cfg.holdoff);
}
