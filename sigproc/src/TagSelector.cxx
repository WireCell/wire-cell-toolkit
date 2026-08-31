#include "WireCellSigProc/TagSelector.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/NamedFactory.h"

#include <sstream>
#include <unordered_map>

WIRECELL_FACTORY(TagSelector, WireCell::SigProc::TagSelector, WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::SigProc;

TagSelector::TagSelector()
    : Aux::Logger("TagSelector", "glue")
{
}

TagSelector::~TagSelector() {}

WireCell::Configuration TagSelector::default_configuration() const
{
    Configuration cfg;

    /// Only traces with these tags will be in the output.
    cfg["tags"] = Json::arrayValue;

    return cfg;
}

void TagSelector::configure(const WireCell::Configuration& cfg)
{
    auto jtags = cfg["tags"];
    int ntags = jtags.size();
    m_tags.clear();
    m_tags.resize(ntags);
    for (int ind = 0; ind < ntags; ++ind) {
        m_tags[ind] = jtags[ind].asString();
    }
}

bool TagSelector::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("see EOS at call={}", m_count);
        ++m_count;
        return true;  // eos
    }

    if (m_tags.empty()) {
        log->warn("TagSelector configured with no tags, passing through input frame");
        out = in;
        return true;
    }

    ITrace::vector out_traces;
    std::unordered_map<const ITrace*, size_t> trace_index;
    std::vector<IFrame::trace_list_t> tagged_trace_indices;
    std::vector<IFrame::trace_summary_t> tagged_trace_summaries;

    tagged_trace_indices.reserve(m_tags.size());
    tagged_trace_summaries.reserve(m_tags.size());

    for (const auto& tag : m_tags) {
        auto traces = Aux::tagged_traces(in, tag);
        auto summary = in->trace_summary(tag);
        IFrame::trace_list_t tl;
        IFrame::trace_summary_t thl;
        tl.reserve(traces.size());
        if (summary.size()) {
            thl.reserve(traces.size());
        }
        for (size_t trind = 0; trind < traces.size(); ++trind) {
            const auto& trace = traces[trind];
            const auto* trace_ptr = trace.get();
            auto it = trace_index.find(trace_ptr);
            if (it == trace_index.end()) {
                size_t index = out_traces.size();
                out_traces.push_back(trace);
                trace_index.emplace(trace_ptr, index);
                it = trace_index.find(trace_ptr);
            }
            tl.push_back(it->second);
            if (summary.size()) {
                thl.push_back(summary[trind]);
            }
        }
        tagged_trace_indices.push_back(tl);
        tagged_trace_summaries.push_back(thl);
    }

    auto sf = new Aux::SimpleFrame(in->ident(), in->time(), out_traces, in->tick(), in->masks());

    for (size_t ind = 0; ind < m_tags.size(); ++ind) {
        const auto& tag = m_tags[ind];
        if (tagged_trace_summaries[ind].size()) {
            sf->tag_traces(tag, tagged_trace_indices[ind], tagged_trace_summaries[ind]);
        }
        else {
            sf->tag_traces(tag, tagged_trace_indices[ind]);
        }
    }

    std::vector<std::string> frame_tags = in->frame_tags();
    if (frame_tags.empty()) {
        frame_tags.push_back("");
    }
    for (const auto& ftag : frame_tags) {
        sf->tag_frame(ftag);
    }

    out = IFrame::pointer(sf);
    std::stringstream info;
    info << "input " << Aux::taginfo(in) << " output: " << Aux::taginfo(out);
    log->debug(info.str());

    return true;
}
