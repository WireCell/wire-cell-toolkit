#include "WireCellGen/Reframer.h"
#include "WireCellUtil/Units.h"

#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/NamedFactory.h"

#include <map>
#include <unordered_set>

WIRECELL_FACTORY(Reframer, WireCell::Gen::Reframer,
                 WireCell::INamed,
                 WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;
using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;

Gen::Reframer::Reframer()
    : Aux::Logger("Reframer", "gen")
{
}

Gen::Reframer::~Reframer() {}

WireCell::Configuration Gen::Reframer::default_configuration() const
{
    Configuration cfg;

    cfg["anode"] = "";

    // tags to find input traces/frames
    cfg["tags"] = Json::arrayValue;
    // tag to apply to output frame
    cfg["frame_tag"] = m_frame_tag;
    cfg["tbin"] = m_tbin;
    cfg["nticks"] = m_nticks;
    cfg["toffset"] = m_toffset;
    cfg["fill"] = m_fill;
    cfg["ignore_tags"] = m_ignore_tags;

    return cfg;
}

void Gen::Reframer::configure(const WireCell::Configuration& cfg)
{
    const std::string anode_tn = cfg["anode"].asString();
    log->debug("using anode: \"{}\"", anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(anode_tn);

    m_input_tags.clear();
    for (auto jtag : cfg["tags"]) {
        m_input_tags.push_back(jtag.asString());
    }
    m_frame_tag = get<std::string>(cfg, "frame_tag", m_frame_tag);
    m_toffset = get(cfg, "toffset", m_toffset);
    m_tbin = get(cfg, "tbin", m_tbin);
    m_fill = get(cfg, "fill", m_fill);
    m_nticks = get(cfg, "nticks", m_nticks);
    m_ignore_tags = get(cfg, "ignore_tags", m_ignore_tags);
    if ( m_ignore_tags && m_input_tags.size() ) {
        raise<ValueError>("providing and ignoring tags is not consistent");
    }
}

std::pair<ITrace::vector, IFrame::trace_summary_t> Gen::Reframer::process_one(const ITrace::vector& itraces, const IFrame::trace_summary_t& isummary) {
    if(isummary.size() !=0 && itraces.size() != isummary.size()) {
        log->error("itraces.size() != isummary.size()");
        THROW(RuntimeError() << errmsg{"itraces.size() != isummary.size()"});
    }
    // Storage for samples indexed by channel ident.
    std::map<int, std::vector<float> > waves;
    // chid -> threshold
    std::map<int, double > thresholds;

    // initialize a "rectangular" 2D array of samples
    for (auto chid : m_anode->channels()) {
        waves[chid].resize(m_nticks, m_fill);
        thresholds[chid] = 0;
    }

    // Lay down input traces over output waves
    for(size_t idx=0; idx< itraces.size(); ++idx) {
        auto trace = itraces[idx];
        const int chid = trace->channel();
        if (isummary.size() !=0) {
            thresholds[chid] = isummary[idx];
        }

        const int tbin_in = trace->tbin();
        auto in_it = trace->charge().begin();

        auto& wave = waves[chid];  // reference
        auto out_it = wave.begin();

        const int delta_tbin = m_tbin - tbin_in;
        if (delta_tbin > 0) {  // must truncate input
            in_it += delta_tbin;
        }
        else {  // must pad output
            out_it += -delta_tbin;
        }

        // Go as far as possible but don't walk of the end of either
        const auto in_end = trace->charge().end();
        const auto out_end = wave.end();
        while (in_it < in_end and out_it < out_end) {
            *out_it += *in_it;  // accumulate
            ++in_it;
            ++out_it;
        }
    }

    // Transfer waves into traces ordered by chid
    ITrace::vector otraces;
    for (auto& it : waves) {
        const int chid = it.first;
        auto& wave = it.second;
        auto out_trace = make_shared<SimpleTrace>(chid, 0, wave);
        otraces.push_back(out_trace);
    }

    // output summary should be same length too, default at 0
    IFrame::trace_summary_t osummary;
    std::transform(thresholds.begin(), thresholds.end(), std::back_inserter(osummary), [&](const auto& pair) { return pair.second; });

    return {otraces, osummary};
}

ITrace::vector Gen::Reframer::process_one(const ITrace::vector& itraces) {
    auto [traces, summary] = process_one(itraces, {});
    return traces;
}

bool Gen::Reframer::operator()(const input_pointer& inframe, output_pointer& outframe)
{
    if (!inframe) {
        outframe = nullptr;
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }

    // Get traces to consider
    ITrace::vector out_traces;
    std::unordered_map< std::string, IFrame::trace_list_t> tag_indicies;
    std::unordered_map< std::string, IFrame::trace_summary_t> tag_summary;

    if (m_input_tags.empty()) {
        auto ttags = inframe->trace_tags();
        if (! ttags.empty()) {
            std::stringstream ss;
            for (const auto& tag : ttags) {
                ss << " " << tag;
            }
            if (! m_ignore_tags) {
                log->warn("will combine traces from {} trace tags in frame:{}", ttags.size(), ss.str());
            }
        }
        out_traces = process_one(*(inframe->traces()));
    }
    else {
        for (const auto& tag : m_input_tags) {
            const auto& isummary = inframe->trace_summary(tag);
            ITrace::vector in_traces = Aux::tagged_traces(inframe, tag);
            auto [out_one, threshold] = process_one(in_traces, isummary);
            tag_summary[tag] = threshold;
            size_t nbeg = out_traces.size();
            out_traces.insert(out_traces.end(), out_one.begin(), out_one.end());
            auto& indicies = tag_indicies[tag];
            indicies.resize(out_one.size());
            std::iota(indicies.begin(), indicies.end(), nbeg);
        }
    }

    auto sframe = make_shared<SimpleFrame>(inframe->ident(), inframe->time() + m_toffset + m_tbin * inframe->tick(),
                                           out_traces, inframe->tick());
    if (! m_frame_tag.empty()) {
        sframe->tag_frame(m_frame_tag);
    }

    for (const auto& tag : m_input_tags) {
        sframe->tag_traces(tag, tag_indicies.at(tag), tag_summary[tag]);
    }

    outframe = sframe;

    log->debug("input : {}", Aux::taginfo(inframe));
    log->debug("output: {}", Aux::taginfo(outframe));

    ++m_count;
    return true;
}
