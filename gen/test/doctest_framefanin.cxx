// FrameFanin merges frames from several ports.  When the tag rules map a
// trace tag of more than one port to the same output tag (eg each anode's
// "gauss" traces -> "gauss" in a detector-wide frame), the output tag must
// cover the traces of all those ports, with their summaries, not only those
// of the last port.

#include "WireCellGen/FrameFanin.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellUtil/doctest.h"

#include <set>

using namespace WireCell;

namespace {
    // A frame with one trace per channel, all tagged with `tag`, with a
    // summary value per trace if `summary` is true.
    IFrame::pointer make_frame(const std::vector<int>& chans, const std::string& tag, bool summary)
    {
        ITrace::vector traces;
        IFrame::trace_list_t inds;
        IFrame::trace_summary_t summ;
        for (int ch : chans) {
            inds.push_back(traces.size());
            summ.push_back(10.0 * ch);
            traces.push_back(std::make_shared<Aux::SimpleTrace>(ch, 0, std::vector<float>(4, 1.0)));
        }
        auto sf = std::make_shared<Aux::SimpleFrame>(7, 0.0, traces, 0.5);
        if (summary) {
            sf->tag_traces(tag, inds, summ);
        }
        else {
            sf->tag_traces(tag, inds);
        }
        return sf;
    }

    Configuration fanin_config(int multiplicity, const std::string& from, const std::string& to)
    {
        Configuration cfg;
        cfg["multiplicity"] = multiplicity;
        for (int ind = 0; ind < multiplicity; ++ind) {
            Configuration rule;
            rule["trace"][from] = to;
            cfg["tag_rules"].append(rule);
        }
        return cfg;
    }

    std::set<int> tagged_channels(const IFrame::pointer& fr, const std::string& tag)
    {
        std::set<int> out;
        auto traces = fr->traces();
        for (auto ind : fr->tagged_traces(tag)) {
            out.insert(traces->at(ind)->channel());
        }
        return out;
    }
}

TEST_CASE("gen framefanin same output trace tag from several ports")
{
    Gen::FrameFanin fan;
    fan.configure(fanin_config(3, "gauss", "gauss"));

    Gen::FrameFanin::input_vector in{
        make_frame({1, 2}, "gauss", true),
        make_frame({3}, "gauss", true),
        make_frame({4, 5}, "gauss", true),
    };
    Gen::FrameFanin::output_pointer out;
    REQUIRE(fan(in, out));
    REQUIRE(out);
    CHECK(out->traces()->size() == 5);
    CHECK(tagged_channels(out, "gauss") == std::set<int>{1, 2, 3, 4, 5});

    // The summary stays aligned with the merged indices.
    const auto& inds = out->tagged_traces("gauss");
    const auto& summ = out->trace_summary("gauss");
    REQUIRE(summ.size() == inds.size());
    auto traces = out->traces();
    for (size_t i = 0; i < inds.size(); ++i) {
        CHECK(summ[i] == 10.0 * traces->at(inds[i])->channel());
    }
}

TEST_CASE("gen framefanin merged tag without summaries")
{
    Gen::FrameFanin fan;
    fan.configure(fanin_config(2, "gauss", "gauss"));

    Gen::FrameFanin::input_vector in{
        make_frame({1}, "gauss", false),
        make_frame({2, 3}, "gauss", false),
    };
    Gen::FrameFanin::output_pointer out;
    REQUIRE(fan(in, out));
    CHECK(tagged_channels(out, "gauss") == std::set<int>{1, 2, 3});
    CHECK(out->trace_summary("gauss").empty());
}

TEST_CASE("gen framefanin distinct output tags are unchanged")
{
    Gen::FrameFanin fan;
    Configuration cfg;
    cfg["multiplicity"] = 2;
    for (int ind = 0; ind < 2; ++ind) {
        Configuration rule;
        rule["trace"]["gauss"] = "gauss" + std::to_string(ind);
        cfg["tag_rules"].append(rule);
    }
    fan.configure(cfg);

    Gen::FrameFanin::input_vector in{
        make_frame({1, 2}, "gauss", true),
        make_frame({3}, "gauss", true),
    };
    Gen::FrameFanin::output_pointer out;
    REQUIRE(fan(in, out));
    CHECK(tagged_channels(out, "gauss0") == std::set<int>{1, 2});
    CHECK(tagged_channels(out, "gauss1") == std::set<int>{3});
    CHECK(out->trace_summary("gauss1").size() == 1);
}

TEST_CASE("gen framefanin takes time and tick from a frame with traces")
{
    // An empty frame (eg from DeposOrBust for an anode without depos) carries
    // no meaningful time or tick; it must not set those of the merged frame.
    Gen::FrameFanin fan;
    fan.configure(fanin_config(3, "gauss", "gauss"));

    auto empty = std::make_shared<Aux::SimpleFrame>(7, 0.0, ITrace::vector{}, 0.0);
    ITrace::vector traces{std::make_shared<Aux::SimpleTrace>(1, 0, std::vector<float>(4, 1.0))};
    auto full = std::make_shared<Aux::SimpleFrame>(7, -250.0, traces, 0.5);
    full->tag_traces("gauss", IFrame::trace_list_t{0});

    Gen::FrameFanin::input_vector in{empty, full, empty};
    Gen::FrameFanin::output_pointer out;
    REQUIRE(fan(in, out));
    CHECK(out->ident() == 7);
    CHECK(out->time() == -250.0);
    CHECK(out->tick() == 0.5);
    CHECK(out->traces()->size() == 1);

    // All empty: the first frame's.
    Gen::FrameFanin::input_vector in2{empty, empty, empty};
    REQUIRE(fan(in2, out));
    CHECK(out->time() == 0.0);
    CHECK(out->traces()->empty());
}
