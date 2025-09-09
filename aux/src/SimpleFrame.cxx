#include "WireCellAux/SimpleFrame.h"

using namespace WireCell;
using namespace std;

Aux::SimpleFrame::SimpleFrame(int ident, double time, const ITrace::vector& traces, double tick,
                         const Waveform::ChannelMaskMap& cmm)
  : m_ident(ident)
  , m_time(time)
  , m_tick(tick)
  , m_traces(std::make_shared<ITrace::vector>(traces.begin(), traces.end()))
  , m_cmm(cmm)
{
}
Aux::SimpleFrame::SimpleFrame(int ident, double time, ITrace::shared_vector traces, double tick,
                         const Waveform::ChannelMaskMap& cmm)
  : m_ident(ident)
  , m_time(time)
  , m_tick(tick)
  , m_traces(traces)
  , m_cmm(cmm)
{
}
Aux::SimpleFrame::SimpleFrame(int ident, double time, double tick)
  : m_ident(ident)
  , m_time(time)
  , m_tick(tick)
  , m_traces(std::make_shared<ITrace::vector>())
{
}

Aux::SimpleFrame::~SimpleFrame() {}
int Aux::SimpleFrame::ident() const { return m_ident; }
double Aux::SimpleFrame::time() const { return m_time; }
double Aux::SimpleFrame::tick() const { return m_tick; }

ITrace::shared_vector Aux::SimpleFrame::traces() const { return m_traces; }

Waveform::ChannelMaskMap Aux::SimpleFrame::masks() const { return m_cmm; }

Aux::SimpleFrame::SimpleTraceInfo::SimpleTraceInfo()
  : indices(0)
  , summary(0)
{
}
const Aux::SimpleFrame::SimpleTraceInfo& Aux::SimpleFrame::get_trace_info(const IFrame::tag_t& tag) const
{
    static SimpleTraceInfo empty;
    auto const& it = m_trace_info.find(tag);
    if (it == m_trace_info.end()) {
        return empty;
    }
    return it->second;
}

const IFrame::tag_list_t& Aux::SimpleFrame::frame_tags() const { return m_frame_tags; }
const IFrame::tag_list_t& Aux::SimpleFrame::trace_tags() const { return m_trace_tags; }

const IFrame::trace_list_t& Aux::SimpleFrame::tagged_traces(const tag_t& tag) const { return get_trace_info(tag).indices; }

const IFrame::trace_summary_t& Aux::SimpleFrame::trace_summary(const tag_t& tag) const
{
    return get_trace_info(tag).summary;
}

void Aux::SimpleFrame::tag_frame(const tag_t& tag) { m_frame_tags.push_back(tag); }

// void Aux::SimpleFrame::tag_traces(const tag_t& tag, const IFrame::trace_list_t& indices,
//                              const IFrame::trace_summary_t& summary)
// {
//     auto& info = m_trace_info[tag];
//     info.indices = indices;
//     info.summary = summary;

//     // Kind of dumb way to update this but we want to be able to
//     // return a reference to it later.
//     m_trace_tags.clear();
//     for (auto& it : m_trace_info) {
//         m_trace_tags.push_back(it.first);
//     }
// }
