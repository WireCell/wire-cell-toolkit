#include "WireCellSigProc/FrameMerger.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellIface/ITrace.h"

#include "WireCellAux/FrameTools.h"

#include <string>

WIRECELL_FACTORY(FrameMerger,
                 WireCell::SigProc::FrameMerger,
                 WireCell::INamed, WireCell::IFrameJoiner, WireCell::IConfigurable)

using namespace WireCell;

Configuration SigProc::FrameMerger::default_configuration() const
{
    Configuration cfg;

    // The merge map specifies a list of triples of tags.  The first (second)
    // tag is compared to trace tags in the frame from port 0 (1).  Either or
    // both of the first two tags may be empty strings which match untagged
    // traces.  The third tag placed on the merged set of traces from these two
    // input in the output frame.  The merge is performed considering only
    // traces in the two input frames that respectively match their tags.  As
    // tags may have overlapping sets of traces, it is important to recognize
    // that the merge is progressive and in order of the merge map.  If the
    // merge map is empty then all traces in frame 1 are merged with all traces
    // in frame 2 irrespective of any tags.
    cfg["mergemap"] = Json::arrayValue;

    // The rule determins the algorithm employed in the merge.
    //
    // - replace :: the last trace encountered for a given channel is output.
    //   The traces from frame 1 will take precedence over any from frame 2.
    //
    // - include :: all traces encountered with a given channel are output.
    //
    // - tbd :: more may be added (eg, sum all traces on a given channel)
    cfg["rule"] = "replace";

    return cfg;
}

void SigProc::FrameMerger::configure(const Configuration& cfg) { m_cfg = cfg; }

bool SigProc::FrameMerger::operator()(const input_tuple_type& intup, output_pointer& out)
{
    out = nullptr;

    auto one = std::get<0>(intup);
    auto two = std::get<1>(intup);
    if (!one or !two) {
        log->debug("EOS at call={}", m_count++);
        return true;
    }

    auto jmergemap = m_cfg["mergemap"];
    const int nsets = jmergemap.size();

    log->debug("call={} frame1: {}", m_count, WireCell::Aux::taginfo(one));
    log->debug("call={} frame2: {}", m_count, WireCell::Aux::taginfo(two));

    // collect traces and matching trace summaries into a vector of vector
    // whether we are dealling with all traces or honoring tags.  Summaries
    // are per-trace (parallel to traces) when present; an empty summary
    // means the input did not attach one to that tag.
    std::vector<ITrace::vector> tracesv1, tracesv2;
    std::vector<IFrame::trace_summary_t> sumv1, sumv2;
    if (!nsets) {
        tracesv1.push_back(Aux::untagged_traces(one));
        tracesv2.push_back(Aux::untagged_traces(two));
        sumv1.emplace_back();
        sumv2.emplace_back();
    }
    else {
        for (int ind = 0; ind < nsets; ++ind) {
            auto jtags = jmergemap[ind];
            std::string tag1 = jtags[0].asString();
            std::string tag2 = jtags[1].asString();
            std::string tag3 = jtags[2].asString();
            tracesv1.push_back(Aux::tagged_traces(one, tag1));
            tracesv2.push_back(Aux::tagged_traces(two, tag2));
            sumv1.push_back(one->trace_summary(tag1));
            sumv2.push_back(two->trace_summary(tag2));
            log->debug("call={} tags: {}[{}] + {}[{}] -> {}",
                       m_count,
                       tag1, tracesv1.back().size(),
                       tag2, tracesv2.back().size(), tag3);
        }
    }

    ITrace::vector out_traces;
    std::vector<IFrame::trace_list_t> tagged_trace_indices;
    std::vector<IFrame::trace_summary_t> tagged_summaries;

    // apply rule, collect info for tags even if we may not be
    // configured to honor them.  Per-trace summaries follow the trace
    // they accompany: if the contributing input's matching tag has a
    // summary, use it; otherwise fall back to the other input's
    // summary at the same channel.  An output-tag summary is attached
    // only if at least one input had a non-empty summary for that tag.

    auto rule = get<std::string>(m_cfg, "rule");
    if (rule == "replace") {
        for (size_t ind = 0; ind < tracesv1.size(); ++ind) {
            auto& traces1 = tracesv1[ind];
            auto& traces2 = tracesv2[ind];
            const auto& s1 = sumv1[ind];
            const auto& s2 = sumv2[ind];
            const bool need_sum = !s1.empty() || !s2.empty();

            // ch -> (port, idx within that port's tagged-trace list)
            std::unordered_map<int, std::pair<int, size_t>> ch2src;
            for (size_t i = 0; i < traces2.size(); ++i) ch2src[traces2[i]->channel()] = {1, i};
            for (size_t i = 0; i < traces1.size(); ++i) ch2src[traces1[i]->channel()] = {0, i};

            std::unordered_map<int, ITrace::pointer> ch2tr;
            for (auto trace : traces2) ch2tr[trace->channel()] = trace;
            for (auto trace : traces1) ch2tr[trace->channel()] = trace;  // port 1 wins

            // Channel -> summary value lookup tables for fallback.
            auto build_lookup = [](const ITrace::vector& tr,
                                   const IFrame::trace_summary_t& s) {
                std::unordered_map<int, double> m;
                if (!s.empty()) {
                    const size_t n = std::min(tr.size(), s.size());
                    for (size_t i = 0; i < n; ++i) m[tr[i]->channel()] = s[i];
                }
                return m;
            };
            const auto ch2s1 = build_lookup(traces1, s1);
            const auto ch2s2 = build_lookup(traces2, s2);

            IFrame::trace_list_t tl;
            IFrame::trace_summary_t tlsum;
            for (auto chtr : ch2tr) {
                tl.push_back(out_traces.size());
                out_traces.push_back(chtr.second);
                if (need_sum) {
                    const int ch = chtr.first;
                    double val = 0.0;
                    auto srcit = ch2src.find(ch);
                    bool ok = false;
                    if (srcit != ch2src.end()) {
                        const int port = srcit->second.first;
                        const size_t idx = srcit->second.second;
                        if (port == 0 && !s1.empty() && idx < s1.size()) { val = s1[idx]; ok = true; }
                        else if (port == 1 && !s2.empty() && idx < s2.size()) { val = s2[idx]; ok = true; }
                    }
                    if (!ok) {
                        // Fall back to the other input's summary keyed by channel.
                        auto it2 = ch2s2.find(ch);
                        if (it2 != ch2s2.end()) val = it2->second;
                        else {
                            auto it1 = ch2s1.find(ch);
                            if (it1 != ch2s1.end()) val = it1->second;
                        }
                    }
                    tlsum.push_back(val);
                }
            }
            tagged_trace_indices.push_back(tl);
            tagged_summaries.push_back(tlsum);
        }
    }
    if (rule == "include") {
        for (size_t ind = 0; ind < tracesv1.size(); ++ind) {
            auto& traces1 = tracesv1[ind];
            auto& traces2 = tracesv2[ind];
            const auto& s1 = sumv1[ind];
            const auto& s2 = sumv2[ind];
            const bool need_sum = !s1.empty() || !s2.empty();

            IFrame::trace_list_t tl;
            IFrame::trace_summary_t tlsum;
            for (size_t trind = 0; trind < traces1.size(); ++trind) {
                tl.push_back(out_traces.size());
                out_traces.push_back(traces1[trind]);
                if (need_sum) tlsum.push_back(trind < s1.size() ? s1[trind] : 0.0);
            }
            for (size_t trind = 0; trind < traces2.size(); ++trind) {
                tl.push_back(out_traces.size());
                out_traces.push_back(traces2[trind]);
                if (need_sum) tlsum.push_back(trind < s2.size() ? s2[trind] : 0.0);
            }
            tagged_trace_indices.push_back(tl);
            tagged_summaries.push_back(tlsum);
        }
    }

    // Inherit channel masks from both inputs (port 1 as the base, then
    // union in port 0 with same mask names).  Frame-level metadata like
    // bad/noisy channel maps must survive the merge so downstream stages
    // (e.g. MagnifySink T_bad tree, imaging) keep seeing them.
    Waveform::ChannelMaskMap merged_masks = two->masks();
    Waveform::ChannelMaskMap one_masks = one->masks();
    std::map<std::string, std::string> name_map;
    Waveform::merge(merged_masks, one_masks, name_map);

    auto sf = new Aux::SimpleFrame(two->ident(), two->time(), out_traces, two->tick(), merged_masks);
    if (nsets) {
        for (int ind = 0; ind < nsets; ++ind) {
            std::string otag = jmergemap[ind][2].asString();
            const auto& tlsum = tagged_summaries[ind];
            if (tlsum.empty()) {
                sf->tag_traces(otag, tagged_trace_indices[ind]);
            }
            else {
                sf->tag_traces(otag, tagged_trace_indices[ind], tlsum);
            }
        }
    }
    out = IFrame::pointer(sf);
    return true;
}

SigProc::FrameMerger::FrameMerger()
    : Aux::Logger("FrameMerger", "sigproc")
{}

SigProc::FrameMerger::~FrameMerger() {}
