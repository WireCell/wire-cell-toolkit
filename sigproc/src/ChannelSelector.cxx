#include "WireCellSigProc/ChannelSelector.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/NamedFactory.h"

#include <sstream>

WIRECELL_FACTORY(ChannelSelector, WireCell::SigProc::ChannelSelector, WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::SigProc;

ChannelSelector::ChannelSelector()
    : Aux::Logger("ChannelSelector", "glue")
{
}

ChannelSelector::~ChannelSelector() {}

WireCell::Configuration ChannelSelector::default_configuration() const
{
    Configuration cfg;

    /// Only traces with channels in this array will be in the output.
    cfg["channels"] = Json::arrayValue;

    /// Only traces with these tags will be in the output.  If no tags
    /// are given then tags are not considered.
    cfg["tags"] = Json::arrayValue;

    /// Rules to govern the output tags based on input tags.
    cfg["tag_rules"] = Json::arrayValue;

    return cfg;
}

void ChannelSelector::configure(const WireCell::Configuration& cfg)
{
    // tags need some order
    auto jtags = cfg["tags"];
    int ntags = jtags.size();
    m_tags.clear();
    m_tags.resize(ntags);
    for (int ind = 0; ind < ntags; ++ind) {
        m_tags[ind] = jtags[ind].asString();
    }

    // channels are just a bag
    for (auto jchan : cfg["channels"]) {
        m_channels.insert(jchan.asInt());
    }

    auto tr = cfg["tag_rules"];
    if (tr.isNull() or tr.empty()) {
        return;
    }
    m_ft.configure(tr);
    m_use_rules = true;
}

void ChannelSelector::set_channels(const std::vector<int>& channels)
{
    m_channels.clear();
    for (int ch : channels) {
        m_channels.insert(ch);
    }
}

bool ChannelSelector::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("see EOS at call={}", m_count);
        ++m_count;
        return true;  // eos
    }

    std::vector<ITrace::vector> tracesvin;

    // size_t ntraces = 0;
    size_t ntags = m_tags.size();
    if (!ntags) {
        tracesvin.push_back(Aux::untagged_traces(in));
        // ntraces += tracesvin[0].size();
    }
    else {
        tracesvin.resize(ntags);
        for (size_t ind = 0; ind < ntags; ++ind) {
            std::string tag = m_tags[ind];
            tracesvin[ind] = Aux::tagged_traces(in, tag);
            // ntraces += tracesvin[ind].size();
        }
    }
    
    ITrace::vector out_traces;
    std::vector<IFrame::trace_list_t> tagged_trace_indices;

    for (size_t ind = 0; ind < tracesvin.size(); ++ind) {
        auto& traces = tracesvin[ind];

        IFrame::trace_list_t tl;
        for (size_t trind = 0; trind < traces.size(); ++trind) {
            auto& trace = traces[trind];
            if (m_channels.find(trace->channel()) == m_channels.end()) {
                continue;
            }
            tl.push_back(out_traces.size());
            out_traces.push_back(trace);
        }
        tagged_trace_indices.push_back(tl);
    }

    auto sf = new Aux::SimpleFrame(in->ident(), in->time(), out_traces, in->tick());
    if (ntags) {
        for (size_t ind = 0; ind < ntags; ++ind) {
            std::string tag = m_tags[ind];
            if (m_use_rules) {
                for (auto new_tag : m_ft.transform(0, "trace", {tag})) {
                    sf->tag_traces(new_tag, tagged_trace_indices[ind]);
                }
            }
            else {
                sf->tag_traces(tag, tagged_trace_indices[ind]);
            }
        }
    }

    std::vector<std::string> frame_tags = in->frame_tags();
    if (frame_tags.empty()) { frame_tags.push_back(""); }
    for (auto ftag : frame_tags) {
        if (m_use_rules) {
            for (auto new_tag : m_ft.transform(0, "frame", {ftag})) {
                sf->tag_frame(new_tag);
            }
        }
        else {
            sf->tag_frame(ftag);
        }
    }

    out = IFrame::pointer(sf);
    std::stringstream info;
    info << "input " << Aux::taginfo(in) << " output: " << Aux::taginfo(out);
    log->debug(info.str());

    return true;
}
