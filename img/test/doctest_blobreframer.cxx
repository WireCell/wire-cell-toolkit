/** BlobReframer must handle slice frames that have no traces.
 *
 *  MaskSlicer gives its slices a TracelessFrame whose traces() returns
 *  nullptr.  BlobReframer dereferenced that pointer to find the number of
 *  output ticks, and so crashed on any cluster made from MaskSlicer output.
 *  It also dereferenced a null frame when the cluster had no blobs.
 */

#include "WireCellImg/BlobReframer.h"
#include "WireCellImg/GridTiling.h"
#include "WireCellImg/MaskSlice.h"
#include "WireCellAux/SimpleCluster.h"
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

    // Activity on a dozen channels in the middle of each plane of the face,
    // in tbins [tbin0, tbin0+nticks).  Traces are tagged as MaskSlicer wants.
    IFrame::pointer make_frame(const IAnodeFace::pointer& face, int tbin0, int nticks)
    {
        ITrace::vector traces;
        IFrame::trace_list_t gauss, wiener, error;
        IFrame::trace_summary_t summary;
        for (const auto& plane : face->planes()) {
            const auto& chans = plane->channels();
            const size_t mid = chans.size() / 2;
            for (size_t i = mid; i < std::min(mid + 12, chans.size()); ++i) {
                const int chid = chans[i]->ident();
                gauss.push_back(traces.size());
                traces.push_back(std::make_shared<Aux::SimpleTrace>(chid, tbin0, std::vector<float>(nticks, 1000.0)));
                wiener.push_back(traces.size());
                summary.push_back(1.0);
                traces.push_back(std::make_shared<Aux::SimpleTrace>(chid, tbin0, std::vector<float>(nticks, 1000.0)));
                error.push_back(traces.size());
                traces.push_back(std::make_shared<Aux::SimpleTrace>(chid, tbin0, std::vector<float>(nticks, 10.0)));
            }
        }
        auto frame = std::make_shared<Aux::SimpleFrame>(7, frame_time, traces, tick);
        frame->tag_traces("gauss", gauss);
        frame->tag_traces("wiener", wiener, summary);
        frame->tag_traces("gauss_error", error);
        return frame;
    }

    // Tile each slice into blobs and gather all blobs into one cluster.
    ICluster::pointer tile(const IAnodePlane::pointer& anode, const IAnodeFace::pointer& face,
                           const ISlice::vector& slices)
    {
        Img::GridTiling gt;
        auto gcfg = gt.default_configuration();
        gcfg["anode"] = anode_tn(anode);
        gcfg["face"] = face->ident();
        gt.configure(gcfg);

        cluster_graph_t g;
        size_t nblobs = 0;
        for (const auto& s : slices) {
            if (!s) continue;
            IBlobSet::pointer bs;
            REQUIRE(gt(s, bs));
            REQUIRE(bs);
            for (const auto& b : bs->blobs()) {
                boost::add_vertex(cluster_node_t(IBlob::pointer(b)), g);
                ++nblobs;
            }
        }
        REQUIRE(nblobs > 0);
        return std::make_shared<Aux::SimpleCluster>(g, 42);
    }

    IFrame::pointer reframe(const ICluster::pointer& icl)
    {
        Img::BlobReframer brf;
        auto cfg = brf.default_configuration();
        brf.configure(cfg);
        IFrame::pointer out;
        REQUIRE(brf(icl, out));
        REQUIRE(out);
        return out;
    }
}

TEST_CASE("BlobReframer handles MaskSlicer's traceless slice frames")
{
    auto anode = get_anode();
    auto face = anode->faces().front();
    REQUIRE(face);

    // Activity in tbins [8, 16), slice bins 2 and 3.
    auto frame = make_frame(face, 8, 8);

    Img::MaskSlicer slicer;
    auto scfg = slicer.default_configuration();
    scfg["anode"] = anode_tn(anode);
    scfg["tick_span"] = tick_span;
    for (int ind = 0; ind < 3; ++ind) {
        scfg["nthreshold"][ind] = 1.0;
    }
    slicer.configure(scfg);
    ISliceFrame::pointer sf;
    REQUIRE(slicer(frame, sf));
    REQUIRE(sf);

    // The slices' frame has no traces.  That is what used to crash.
    const auto slices = sf->slices();
    REQUIRE(slices.size() > 0);
    REQUIRE(slices.front()->frame());
    CHECK(slices.front()->frame()->traces() == nullptr);

    auto out = reframe(tile(anode, face, slices));
    CHECK(out->ident() == frame->ident());
    CHECK(out->time() == doctest::Approx(frame_time));
    CHECK(out->tick() == doctest::Approx(tick));

    // Without input traces the output spans up to the end of the last slice.
    auto traces = out->traces();
    REQUIRE(traces);
    REQUIRE(traces->size() > 0);
    double qtot = 0;
    for (const auto& tr : *traces) {
        CHECK(tr->tbin() == 0);
        CHECK(tr->charge().size() == 16);
        for (int ind = 0; ind < 8; ++ind) {
            CHECK(tr->charge()[ind] == 0.0);
        }
        for (auto q : tr->charge()) qtot += q;
    }
    CHECK(qtot > 0);
}

TEST_CASE("BlobReframer output length follows the slice frame's traces when it has them")
{
    auto anode = get_anode();
    auto face = anode->faces().front();
    REQUIRE(face);

    // A frame with traces spanning [0, 100) and a slice at tbins [8, 12).
    const int chid = face->planes().front()->channels().front()->ident();
    ITrace::vector ftraces{std::make_shared<Aux::SimpleTrace>(chid, 0, std::vector<float>(100, 0.0))};
    auto frame = std::make_shared<Aux::SimpleFrame>(7, frame_time, ftraces, tick);

    ISlice::map_t activity;
    for (const auto& plane : face->planes()) {
        const auto& chans = plane->channels();
        const size_t mid = chans.size() / 2;
        for (size_t i = mid; i < std::min(mid + 12, chans.size()); ++i) {
            activity[chans[i]] = ISlice::value_t(1000.0, 10.0);
        }
    }
    auto slice = std::make_shared<Aux::SimpleSlice>(frame, 2, 8 * tick, tick_span * tick, activity);

    auto out = reframe(tile(anode, face, ISlice::vector{slice}));
    CHECK(out->time() == doctest::Approx(frame_time));
    auto traces = out->traces();
    REQUIRE(traces);
    REQUIRE(traces->size() > 0);
    for (const auto& tr : *traces) {
        CHECK(tr->tbin() == 0);
        CHECK(tr->charge().size() == 100);
    }
}

TEST_CASE("BlobReframer gives an empty frame for a cluster with no blobs")
{
    cluster_graph_t g;
    auto icl = std::make_shared<Aux::SimpleCluster>(g, 42);
    auto out = reframe(icl);
    auto traces = out->traces();
    REQUIRE(traces);
    CHECK(traces->size() == 0);
}
