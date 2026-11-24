/** Aux::to_cluster() must not invent a frame time from slice start times.
 *
 *  ISlice::start() is relative to the slice frame's time.  Cluster arrays
 *  store a slice's frame ident but not the frame's time, so the loader can
 *  not know the true origin.  It used to set the invented frame's time to the
 *  earliest slice start, which only made sense if start times were absolute,
 *  and which made "frame time + start" count the earliest start twice.  It
 *  now uses 0, meaning "unknown".  Storing the true frame time is suggestion
 *  "B" in docs/components/BlobSetReframer.md.
 */

#include "WireCellAux/ClusterArrays.h"
#include "WireCellAux/SimpleChannel.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/SimpleWire.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <map>
#include <vector>

using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::Aux::ClusterArrays;

TEST_CASE("cluster arrays round trip keeps relative slice starts and does not invent a frame time")
{
    const int nchan = 3;
    cluster_graph_t g;
    std::vector<IChannel::pointer> chans;
    std::vector<cluster_vertex_t> cvs;
    for (int i = 0; i < nchan; ++i) {
        const int chid = 300 + i;
        const WirePlaneId wpid(kWlayer, 0, 0);
        Ray ray(Point(0, -1000, 3 * i), Point(0, 1000, 3 * i));
        auto iwire = std::make_shared<SimpleWire>(wpid, i, i, chid, ray, 0);
        auto ichan = std::make_shared<SimpleChannel>(chid, i, IWire::vector{iwire});
        chans.push_back(ichan);
        auto cv = boost::add_vertex(cluster_node_t(IChannel::pointer(ichan)), g);
        auto wv = boost::add_vertex(cluster_node_t(IWire::pointer(iwire)), g);
        boost::add_edge(cv, wv, g);
        cvs.push_back(cv);
    }

    const double frame_time = 1.5 * units::ms;
    const double span = 2.0 * units::us;
    std::map<int, double> starts{{2, 4.0 * units::us}, {3, 6.0 * units::us}};
    auto frame = std::make_shared<SimpleFrame>(100, frame_time);
    for (const auto& [sident, start] : starts) {
        ISlice::map_t activity;
        for (const auto& ich : chans) {
            activity[ich] = ISlice::value_t(10.0, 1.0);
        }
        auto islice = std::make_shared<SimpleSlice>(frame, sident, start, span, activity);
        auto sv = boost::add_vertex(cluster_node_t(ISlice::pointer(islice)), g);
        for (auto cv : cvs) boost::add_edge(sv, cv, g);
    }

    node_array_set_t nas;
    edge_array_set_t eas;
    to_arrays(g, nas, eas);
    auto g2 = to_cluster(nas, eas, {});

    int nslices = 0;
    for (const auto& v : boost::make_iterator_range(boost::vertices(g2))) {
        const auto& node = g2[v];
        if (node.code() != 's') continue;
        ++nslices;
        auto islice = std::get<ISlice::pointer>(node.ptr);
        REQUIRE(islice);
        REQUIRE(starts.count(islice->ident()) == 1);
        CHECK(islice->start() == doctest::Approx(starts[islice->ident()]));
        CHECK(islice->span() == doctest::Approx(span));
        REQUIRE(islice->frame());
        CHECK(islice->frame()->ident() == 100);
        // The true frame time is not stored, so the loader must say "unknown" (0).
        CHECK(islice->frame()->time() == 0.0);
    }
    CHECK(nslices == (int) starts.size());
}
