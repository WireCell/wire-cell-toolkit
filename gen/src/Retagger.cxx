#include "WireCellGen/Retagger.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/SimpleFrame.h"

WIRECELL_FACTORY(Retagger, WireCell::Gen::Retagger, WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;
using WireCell::Aux::SimpleFrame;

Gen::Retagger::Retagger() {}

Gen::Retagger::~Retagger() {}

void Gen::Retagger::configure(const WireCell::Configuration& cfg)
{
    // frame/trace/merge tagging.
    m_trctx.configure(cfg["tag_rules"]);
    
    auto jmm = cfg["maskmap"];
    for (auto name : jmm.getMemberNames()) {
        m_maskmap[name] = jmm[name].asString();
    }
}

WireCell::Configuration Gen::Retagger::default_configuration() const
{
    Configuration cfg;
    cfg["tag_rules"] = Json::arrayValue;
    cfg["maskmap"]["bad"] = "bad";
    return cfg;
}

bool Gen::Retagger::operator()(const input_pointer& inframe, output_pointer& outframe)
{
    outframe = nullptr;
    if (!inframe) {
        return true;  // eos
    }

    Waveform::ChannelMaskMap input_cmm = inframe->masks();
    Waveform::ChannelMaskMap output_cmm;
    Waveform::merge(output_cmm, input_cmm, m_maskmap);

    // Basic frame data is just shunted across.
    auto sfout = new SimpleFrame(inframe->ident(), inframe->time(), *inframe->traces(), inframe->tick(), output_cmm);

    //
    // frame
    auto fintags = inframe->frame_tags();
    // Add empty tag to allow rules to create an output tag even if none exist.
    // No equivalent should be made for trace tags.
    fintags.push_back("");
    auto fouttags = m_trctx.transform(0, "frame", fintags);
    for (auto ftag : fouttags) {
        sfout->tag_frame(ftag);
    }

    //
    // trace
    for (auto inttag : inframe->trace_tags()) {
        tagrules::tagset_t touttags = m_trctx.transform(0, "trace", inttag);
        if (touttags.empty()) {
            continue;
        }
        const auto& traces = inframe->tagged_traces(inttag);
        const auto& summary = inframe->trace_summary(inttag);
        for (auto otag : touttags) {
            sfout->tag_traces(otag, traces, summary);
        }
    }

    //
    // merge
    typedef std::tuple<tagrules::tag_t, IFrame::trace_list_t, IFrame::trace_summary_t> tag_trace_list_summary_t;
    typedef std::vector<tag_trace_list_summary_t> ttls_stash_t;
    typedef std::unordered_map<std::string, ttls_stash_t> ttls_stash_map_t;
    ttls_stash_map_t stashmap;
    for (auto inttag : inframe->trace_tags()) {
        tagrules::tagset_t touttags = m_trctx.transform(0, "merge", inttag);
        // std::cerr << "Retagger: " << inttag << " -> " << touttags.size() << " outtags\n";
        if (touttags.empty()) {
            continue;
        }
        const auto& traces = inframe->tagged_traces(inttag);
        const auto& summary = inframe->trace_summary(inttag);
        for (auto otag : touttags) {
            stashmap[otag].push_back(std::make_tuple(inttag, traces, summary));
        }
    }
    for (auto it : stashmap) {
        auto& otag = it.first;
        auto& stashv = it.second;
        size_t ntraces = 0, nsummary = 0;
        // If at least one summary is not empty we must assure output
        // summary is aligned with output trace index vector.
        // Otherwise, we assume input plays by the rules and keeps
        // size of traces and summary equal.
        for (auto& ttls : stashv) {
            ntraces += get<1>(ttls).size();
            nsummary += get<2>(ttls).size();
        }
        IFrame::trace_list_t otraces;
        otraces.reserve(ntraces);
        IFrame::trace_summary_t osummary;
        if (nsummary > 0) {
            osummary.reserve(ntraces);
        }
        for (auto& ttls : stashv) {  // one more time
            auto& traces = get<1>(ttls);
            otraces.insert(otraces.end(), traces.begin(), traces.end());

            // std::cerr << "Retagger: merge: " << get<0>(ttls) << " -> " << otag
            //          << " with " << otraces.size() << " / " << ntraces << std::endl;
            if (!nsummary) {
                continue;
            }
            auto& summary = get<2>(ttls);
            if (summary.empty()) {
                for (size_t ind = 0; ind < traces.size(); ++ind) {
                    summary.push_back(0);  // zero pad to cover missing summary and keep alignment with traces.
                }
            }
            else {
                osummary.insert(osummary.end(), summary.begin(), summary.end());
            }
        }
        // std::cerr << "Retagger: tagging trace set: " << otag << " with " << otraces.size() << " traces, "
        //           << osummary.size() << " summary\n";

        sfout->tag_traces(otag, otraces, osummary);
    }

    outframe = IFrame::pointer(sfout);
    return true;
}
