/** ISlice::start() is relative to the time of the frame the slice came from.
 *
 *  See ISlice.h and docs/components/BlobSetReframer.md.  A 2023 trial
 *  (SLICE_START_TIME_IS_RELATIVE) made SumSlicer, but not MaskSlicer, produce
 *  absolute start times.  The trial was abandoned, so the slicers now agree.
 *  BlobSetReframer used the earliest slice start as its output frame time,
 *  which loses the frame time.  It now takes the time origin from the slice's
 *  frame.
 */

#include "WireCellImg/SumSlice.h"
#include "WireCellImg/MaskSlice.h"
#include "WireCellImg/GridTiling.h"
#include "WireCellImg/BlobSetReframer.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/Testing.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IWirePlane.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <vector>

using namespace WireCell;

namespace {

    const double frame_time = 1.5 * units::ms;
    const double tick = 0.5 * units::us;
    const int tick_span = 4;

    IAnodePlane::pointer get_anode()
    {
        auto anodes = Testing::anodes("uboone");
        REQUIRE(anodes.size() > 0);
        return anodes[0];
    }

    std::string anode_tn(const IAnodePlane::pointer& anode)
    {
        return String::format("AnodePlane:%d", anode->ident());
    }

    // One trace on the anode's first channel with nonzero charge in tbins
    // [10, 15).  With tick_span=4 that spans slice bins 2 and 3.
    ITrace::pointer make_trace(const IAnodePlane::pointer& anode, float q)
    {
        const int chid = anode->channels().front();
        return std::make_shared<Aux::SimpleTrace>(chid, 10, std::vector<float>(5, q));
    }

    void check_slices_relative(const ISliceFrame::pointer& sf)
    {
        REQUIRE(sf);
        CHECK(sf->time() == doctest::Approx(frame_time));
        const auto slices = sf->slices();
        REQUIRE(slices.size() > 0);
        for (const auto& s : slices) {
            if (!s) continue;   // EOS marker
            const double rel = s->ident() * tick_span * tick;
            CHECK(s->start() == doctest::Approx(rel));
            REQUIRE(s->frame());
            CHECK(s->frame()->time() == doctest::Approx(frame_time));
        }
    }
}

TEST_CASE("SumSlicer slice start is relative to frame time")
{
    auto anode = get_anode();
    ITrace::vector traces{make_trace(anode, 10.0)};
    auto frame = std::make_shared<Aux::SimpleFrame>(7, frame_time, traces, tick);

    Img::SumSlicer slicer;
    auto cfg = slicer.default_configuration();
    cfg["anode"] = anode_tn(anode);
    cfg["tick_span"] = tick_span;
    cfg["tag"] = "";
    slicer.configure(cfg);

    ISliceFrame::pointer sf;
    REQUIRE(slicer(frame, sf));
    check_slices_relative(sf);
}

TEST_CASE("MaskSlicer slice start is relative to frame time")
{
    auto anode = get_anode();
    // Traces at indices 0, 1, 2 are tagged gauss, wiener and gauss_error.
    ITrace::vector traces{make_trace(anode, 10.0), make_trace(anode, 100.0), make_trace(anode, 1.0)};
    auto frame = std::make_shared<Aux::SimpleFrame>(7, frame_time, traces, tick);
    frame->tag_traces("gauss", IFrame::trace_list_t{0});
    frame->tag_traces("wiener", IFrame::trace_list_t{1}, IFrame::trace_summary_t{1.0});
    frame->tag_traces("gauss_error", IFrame::trace_list_t{2});

    Img::MaskSlicer slicer;
    auto cfg = slicer.default_configuration();
    cfg["anode"] = anode_tn(anode);
    cfg["tick_span"] = tick_span;
    for (int ind = 0; ind < 3; ++ind) {
        cfg["nthreshold"][ind] = 1.0;
    }
    slicer.configure(cfg);

    ISliceFrame::pointer sf;
    REQUIRE(slicer(frame, sf));
    check_slices_relative(sf);
}

TEST_CASE("BlobSetReframer output frame keeps the slice frame time")
{
    auto anode = get_anode();
    auto face = anode->faces().front();
    REQUIRE(face);

    // A slice starting 4us after its frame's time, 2us long.
    const double start = 4 * units::us;
    const double span = tick_span * tick;
    auto frame = std::make_shared<Aux::SimpleFrame>(7, frame_time, tick);
    ISlice::map_t activity;
    for (const auto& plane : face->planes()) {
        const auto& chans = plane->channels();
        const size_t mid = chans.size() / 2;
        for (size_t i = mid; i < std::min(mid + 12, chans.size()); ++i) {
            activity[chans[i]] = ISlice::value_t(1000.0, 10.0);
        }
    }
    auto slice = std::make_shared<Aux::SimpleSlice>(frame, 1, start, span, activity);

    Img::GridTiling gt;
    auto gcfg = gt.default_configuration();
    gcfg["anode"] = anode_tn(anode);
    gcfg["face"] = face->ident();
    gt.configure(gcfg);
    IBlobSet::pointer bs;
    REQUIRE(gt(slice, bs));
    REQUIRE(bs);
    REQUIRE(bs->blobs().size() > 0);

    for (const std::string source : {"blob", "activity"}) {
        CAPTURE(source);
        Img::BlobSetReframer brf;
        auto cfg = brf.default_configuration();
        cfg["tick"] = tick;
        cfg["source"] = source;
        brf.configure(cfg);

        IFrame::pointer out;
        REQUIRE(brf(bs, out));
        REQUIRE(out);
        CHECK(out->time() == doctest::Approx(frame_time));
        CHECK(out->tick() == doctest::Approx(tick));
        auto otraces = out->traces();
        REQUIRE(otraces);
        REQUIRE(otraces->size() > 0);
        for (const auto& tr : *otraces) {
            CHECK(tr->tbin() == (int) std::round(start / tick));
            CHECK(tr->charge().size() == (size_t) tick_span);
        }
    }
}
