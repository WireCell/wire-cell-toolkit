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
    std::vector<IFrame::trace_summary_t> summariesvin; //added Ewerton 2023-10-04    

    // size_t ntraces = 0;
    size_t ntags = m_tags.size();
    if (!ntags) {
        tracesvin.push_back(Aux::untagged_traces(in));
        summariesvin.push_back({});
        log->warn("Untagged summary not supported, summary will be dropped.");
        // FIXME: need to support summary for untagged traces
        // ntraces += tracesvin[0].size();
    }
    else {
        tracesvin.resize(ntags);
        summariesvin.resize(ntags);
        for (size_t ind = 0; ind < ntags; ++ind) {
            std::string tag = m_tags[ind];
            tracesvin[ind] = Aux::tagged_traces(in, tag);
            summariesvin[ind] = in->trace_summary(tag);//added Ewerton 2023-10-04
            // std::cerr << "\nChannelSelector: tag=" << tag << "\n";//added Ewerton 2023-10-04
            // ntraces += tracesvin[ind].size();
        }
    }


    ITrace::vector out_traces;
    std::vector<IFrame::trace_list_t> tagged_trace_indices;
    std::vector<IFrame::trace_summary_t> tagged_trace_summaries;//added Ewerton 2023-10-04


    for (size_t ind = 0; ind < tracesvin.size(); ++ind) {
        //std::cerr << "\n tracesvin[ind].size()=" << tracesvin[ind].size() << "\n";//added Ewerton 2023-10-02
        auto& traces = tracesvin[ind];
        auto& summary = summariesvin[ind];//added Ewerton 2023-10-04
        IFrame::trace_list_t tl;
        IFrame::trace_summary_t thl; //added Ewerton 2023-10-04
        for (size_t trind = 0; trind < traces.size(); ++trind) {
            auto& trace = traces[trind];
            // DEBUG Ewerton 2024-04-15
            // std::cerr << "\n [ChannelSelector] summary.size()=" << summary.size() << "\n";
            // end DEBUG
           auto threshold = summary.size() ? summary[trind] : -999; // added Ewerton 2023-10-04
            if (m_channels.find(trace->channel()) == m_channels.end()) {
                continue;
            }
            tl.push_back(out_traces.size());
            if(summary.size()) thl.push_back(threshold); //added Ewerton 2023-10-04
            out_traces.push_back(trace);
            // summary[trind] => element
            // trind => element index
            // sl.push_back(out..) => index of trace in out_traces
        }
        tagged_trace_indices.push_back(tl);
        tagged_trace_summaries.push_back(thl); //empty vector if there is no summary for given tag. added Ewerton 2023-10-04
    }

    // auto sf = new Aux::SimpleFrame(in->ident(), in->time(), out_traces, in->tick()); // original
    auto sf = new Aux::SimpleFrame(in->ident(), in->time(), out_traces, in->tick(), in->masks()); // changed Ewerton 2023-10-??

    if (ntags) {
        for (size_t ind = 0; ind < ntags; ++ind) {
            std::string tag = m_tags[ind];
            if (m_use_rules) {
                for (auto new_tag : m_ft.transform(0, "trace", {tag})) {
                    sf->tag_traces(new_tag, tagged_trace_indices[ind]);
                }
            }
            else {
                //sf->tag_traces(tag, tagged_trace_indices[ind]); //original. commented Ewerton 2023-10-02
                if(tagged_trace_summaries[ind].size())
                  sf->tag_traces(tag, tagged_trace_indices[ind], tagged_trace_summaries[ind]); //added Ewerton 2023-10-04
                else
                  sf->tag_traces(tag, tagged_trace_indices[ind]); //added Ewerton 2023-10-04
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

       
/*
Waveform::ChannelMaskMap deb_masks = out->masks(); //default cmm from masks() function is empty cmm -> create local CMM?
std::cerr << "\n ChannelSelector deb_masks.size()=" << deb_masks.size() << "\n\n";    std::string cmm_tag("bad"); // use single tag for cmm: "bad"
Waveform::ChannelMasks chm = deb_masks[cmm_tag]; // use single tag for now: "bad"
int deleteme2=0;
for(auto cm : chm)
{
  deleteme2++;
  if(deleteme2<10) {
    auto ch = cm.first;
    std::cerr << "\n ChannelSelector ChannelMaskMap: ch=" << ch << "\n\n";
    Waveform::BinRangeList br = cm.second;
    for(auto r : br) std::cerr << "\n ChannelSelector rlow=" << r.first << ", rup=" << r.second;
  }
}
*/    
// ------------end debug Ewerton 2023-08-22 ------------

    return true;
}
