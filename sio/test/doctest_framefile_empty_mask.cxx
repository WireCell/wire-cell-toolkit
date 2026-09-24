/** FrameFileSink -> FrameFileSource round trip of a frame whose channel mask
 *  map exists but is EMPTY (a clean simulation through OmnibusNoiseFilter
 *  produces {"bad": {}}).
 *
 *  The sink writes such a map as a (0,3) chanmask array.  pigenc returns a
 *  null data pointer for an empty array, pigenc::eigen::load() then returns
 *  false, and FrameFileSource used to throw IOError("numpy parse error of
 *  chanmask") on it -- so every clean-sim frame file with masks:true was
 *  unreadable (found 2026-09-24 on the dune10kt-1x2x6 iso-track sim,
 *  wcp-porting-img/wcfm/docs/02).  The source must read the frame back with
 *  the tag present and empty; a non-empty mask must still round-trip.
 */

#include "WireCellSio/FrameFileSink.h"
#include "WireCellSio/FrameFileSource.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellUtil/Units.h"

#include "WireCellUtil/doctest.h"

#include <cstdio>
#include <cstdlib>
#include <string>

using namespace WireCell;

namespace {
    std::string tmpname(const std::string& base)
    {
        const char* dir = std::getenv("TMPDIR");
        return std::string(dir ? dir : "/tmp") + "/" + base;
    }

    IFrame::pointer make_frame(int ident, const Waveform::ChannelMaskMap& cmm)
    {
        ITrace::vector traces;
        for (int ch = 0; ch < 4; ++ch) {
            ITrace::ChargeSequence q(10, 1.0f * (ch + 1));
            traces.push_back(std::make_shared<Aux::SimpleTrace>(ch, 0, q));
        }
        auto sf = std::make_shared<Aux::SimpleFrame>(ident, 0.0, traces, 0.5 * units::us, cmm);
        IFrame::trace_list_t idx = {0, 1, 2, 3};
        sf->tag_traces("gauss", idx);
        return sf;
    }

    IFrame::pointer round_trip(const std::string& fname, const IFrame::pointer& in)
    {
        std::remove(fname.c_str());
        {
            Sio::FrameFileSink sink;
            auto cfg = sink.default_configuration();
            cfg["outname"] = fname;
            cfg["tags"][0] = "gauss";
            cfg["digitize"] = false;
            cfg["masks"] = true;
            sink.configure(cfg);
            REQUIRE(sink(in));
            IFrame::pointer eos;
            REQUIRE(sink(eos));
        }
        Sio::FrameFileSource src;
        auto cfg = src.default_configuration();
        cfg["inname"] = fname;
        cfg["tags"][0] = "gauss";
        src.configure(cfg);
        IFrame::pointer out;
        REQUIRE(src(out));
        REQUIRE(out);
        return out;
    }
}

TEST_CASE("framefile empty channel mask round trip")
{
    Waveform::ChannelMaskMap cmm;
    cmm["bad"] = Waveform::ChannelMasks{};   // the map exists, no channel is masked
    auto in = make_frame(100, cmm);
    REQUIRE(in->masks().size() == 1);

    auto out = round_trip(tmpname("wct-framefile-empty-mask.tar"), in);
    CHECK(out->ident() == 100);
    CHECK(out->traces()->size() == 4);
    auto masks = out->masks();
    REQUIRE(masks.count("bad") == 1);
    CHECK(masks["bad"].empty());
}

TEST_CASE("framefile non-empty channel mask round trip")
{
    Waveform::ChannelMaskMap cmm;
    cmm["bad"][2].push_back({3, 7});
    cmm["bad"][3].push_back({0, 10});
    cmm["bad"][3].push_back({20, 30});
    auto in = make_frame(101, cmm);

    auto out = round_trip(tmpname("wct-framefile-nonempty-mask.tar"), in);
    auto masks = out->masks();
    REQUIRE(masks.count("bad") == 1);
    CHECK(masks["bad"].size() == 2);
    CHECK(masks["bad"][2].size() == 1);
    CHECK(masks["bad"][2][0].first == 3);
    CHECK(masks["bad"][2][0].second == 7);
    CHECK(masks["bad"][3].size() == 2);
}
