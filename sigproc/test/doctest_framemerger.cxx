// FrameMerger tests post commit 7d05a7cb (sigproc/FrameMerger:
// propagate trace summaries and channel masks).
//
// The merger combines traces from two input frames according to a
// configurable mergemap (tag1, tag2 -> out_tag) and rule:
//
//   - replace: per channel, port-0 trace wins over port-1
//   - include: union all traces from both ports
//
// In addition to traces, the merged output must carry per-trace
// summaries (one float per trace under each output tag) and the union
// of channel masks from both input frames.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellIface/IFrameJoiner.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"

#include <algorithm>

using namespace WireCell;
using namespace WireCell::Aux;

namespace {

ITrace::pointer make_trace(int ch, int tbin, std::vector<float> charges)
{
    return std::make_shared<SimpleTrace>(ch, tbin, charges);
}

// Build a SimpleFrame with a single tagged trace set + summary + a
// named channel mask.
IFrame::pointer make_frame(int ident,
                           const std::string& tag,
                           const ITrace::vector& traces,
                           const IFrame::trace_summary_t& summary,
                           const std::string& mask_name,
                           const std::vector<int>& mask_channels)
{
    Waveform::ChannelMaskMap cmm;
    if (!mask_channels.empty()) {
        Waveform::ChannelMasks cm;
        for (int ch : mask_channels) {
            cm[ch] = { {0, 100} };
        }
        cmm[mask_name] = cm;
    }
    auto sf = std::make_shared<SimpleFrame>(ident, 0.0, traces, 0.5, cmm);
    IFrame::trace_list_t tl;
    for (size_t i = 0; i < traces.size(); ++i) tl.push_back(i);
    sf->tag_traces(tag, tl, summary);
    return sf;
}

// Find traces with a tag and return parallel (channel, summary) pairs.
struct TaggedView {
    std::vector<int> chs;
    std::vector<double> sums;
};
TaggedView gather(IFrame::pointer frame, const std::string& tag)
{
    TaggedView v;
    auto traces = frame->traces();
    const auto& tl = frame->tagged_traces(tag);
    const auto& s = frame->trace_summary(tag);
    for (size_t i = 0; i < tl.size(); ++i) {
        v.chs.push_back(traces->at(tl[i])->channel());
    }
    v.sums.assign(s.begin(), s.end());
    return v;
}

}  // namespace

TEST_CASE("FrameMerger replace rule, summaries, channel masks") {
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellSigProc");
    pm.add("WireCellAux");

    // port 0: channels {1, 2} with summaries {10, 20}, mask "bad" on ch 1
    // port 1: channels {2, 3} with summaries {200, 300}, mask "bad" on ch 3
    // replace rule + mergemap [tagA, tagB, tagOut]:
    //   merged channels = {1, 2, 3}, summaries from port 0 where ch
    //   exists in port 0 (1 and 2), else from port 1 (3).
    //   mask "bad" must include both 1 and 3.

    auto f0 = make_frame(
        100, "tagA",
        { make_trace(1, 0, {1.0f}), make_trace(2, 0, {2.0f}) },
        {10.0, 20.0},
        "bad", {1});
    auto f1 = make_frame(
        100, "tagB",
        { make_trace(2, 0, {2.5f}), make_trace(3, 0, {3.0f}) },
        {200.0, 300.0},
        "bad", {3});

    auto icfg = Factory::lookup_tn<IConfigurable>("FrameMerger");
    auto cfg = icfg->default_configuration();
    cfg["rule"] = "replace";
    cfg["mergemap"][0][0] = "tagA";
    cfg["mergemap"][0][1] = "tagB";
    cfg["mergemap"][0][2] = "tagOut";
    icfg->configure(cfg);
    auto fj = Factory::lookup_tn<IFrameJoiner>("FrameMerger");

    IFrameJoiner::input_tuple_type in {f0, f1};
    IFrame::pointer out;
    REQUIRE(fj->operator()(in, out));
    REQUIRE(out);

    auto v = gather(out, "tagOut");
    // sort parallel by channel
    std::vector<size_t> order(v.chs.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](size_t a, size_t b) { return v.chs[a] < v.chs[b]; });
    std::vector<int> chs;
    std::vector<double> sums;
    for (auto i : order) { chs.push_back(v.chs[i]); sums.push_back(v.sums[i]); }

    CHECK(chs == std::vector<int>{1, 2, 3});
    CHECK(sums == std::vector<double>{10.0, 20.0, 300.0});  // 2 from port 0 (port 0 wins)

    // Mask union: "bad" must contain both 1 and 3.
    auto masks = out->masks();
    REQUIRE(masks.count("bad"));
    auto& cm = masks["bad"];
    CHECK(cm.count(1));
    CHECK(cm.count(3));
}

TEST_CASE("FrameMerger include rule sums to all-traces output") {
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellSigProc");
    pm.add("WireCellAux");

    auto f0 = make_frame(
        100, "tagA",
        { make_trace(1, 0, {1.0f}), make_trace(2, 0, {2.0f}) },
        {10.0, 20.0}, "", {});
    auto f1 = make_frame(
        100, "tagB",
        { make_trace(2, 0, {2.5f}), make_trace(3, 0, {3.0f}) },
        {200.0, 300.0}, "", {});

    auto icfg = Factory::lookup_tn<IConfigurable>("FrameMerger");
    auto cfg = icfg->default_configuration();
    cfg["rule"] = "include";
    cfg["mergemap"][0][0] = "tagA";
    cfg["mergemap"][0][1] = "tagB";
    cfg["mergemap"][0][2] = "tagOut";
    icfg->configure(cfg);
    auto fj = Factory::lookup_tn<IFrameJoiner>("FrameMerger");

    IFrameJoiner::input_tuple_type in {f0, f1};
    IFrame::pointer out;
    REQUIRE(fj->operator()(in, out));
    REQUIRE(out);

    auto v = gather(out, "tagOut");
    // include preserves both per-port channel-2 traces -> 4 traces total
    CHECK(v.chs.size() == 4);
    CHECK(v.sums.size() == 4);
    // sums vector holds [port0 sums..., port1 sums...]
    CHECK(v.sums[0] == 10.0);
    CHECK(v.sums[1] == 20.0);
    CHECK(v.sums[2] == 200.0);
    CHECK(v.sums[3] == 300.0);
}
