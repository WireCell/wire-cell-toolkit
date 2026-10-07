/** ClusterArrays::bodge_channel_slice_inplace() must produce exactly the graph
 *  that the by-value bodge_channel_slice() returns: the same vertices in the
 *  same order with the same node payloads, and the same edges in the same
 *  order.  to_arrays() now bodges its one working copy in place instead of
 *  calling the by-value form, which (boost::adjacency_list has no move
 *  constructor) held three copies of a whole-event graph at its return
 *  (wcp-porting-img wcfm/docs/13).
 */

#include "WireCellAux/ClusterArrays.h"
#include "WireCellAux/SimpleChannel.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/SimpleWire.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <utility>
#include <vector>

using namespace WireCell;
using namespace WireCell::Aux;

namespace {
    // Two slices sharing channels; each channel has one wire (c-w edges) and
    // each slice's activity names a subset of the channels (s-c edges).
    cluster_graph_t make_graph()
    {
        const int nchan = 5;
        cluster_graph_t g;
        std::vector<IChannel::pointer> chans;
        std::vector<cluster_vertex_t> cvs;
        for (int i = 0; i < nchan; ++i) {
            const int chid = 200 + i;
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
        auto frame = std::make_shared<SimpleFrame>(100, 0.0);
        for (int s = 0; s < 2; ++s) {
            ISlice::map_t activity;
            for (int i = s; i < nchan; i += 1 + s) {
                activity[chans[i]] = ISlice::value_t(10.0 * (i + 1) + s, 1.0);
            }
            auto islice = std::make_shared<SimpleSlice>(frame, s, 2.0 * s * units::us, 2.0 * units::us, activity);
            auto sv = boost::add_vertex(cluster_node_t(ISlice::pointer(islice)), g);
            for (int i = s; i < nchan; i += 1 + s) {
                boost::add_edge(sv, cvs[i], g);
            }
        }
        return g;
    }

    std::vector<cluster_node_t> nodes_of(const cluster_graph_t& g)
    {
        std::vector<cluster_node_t> ret;
        for (auto v : boost::make_iterator_range(boost::vertices(g))) {
            ret.push_back(g[v]);
        }
        return ret;
    }

    std::vector<std::pair<size_t, size_t>> edges_of(const cluster_graph_t& g)
    {
        std::vector<std::pair<size_t, size_t>> ret;
        for (auto e : boost::make_iterator_range(boost::edges(g))) {
            ret.emplace_back(boost::source(e, g), boost::target(e, g));
        }
        return ret;
    }
}

TEST_CASE("bodge_channel_slice_inplace equals the by-value bodge")
{
    const cluster_graph_t orig = make_graph();
    const auto orig_nodes = nodes_of(orig);
    const auto orig_edges = edges_of(orig);

    const cluster_graph_t byval = ClusterArrays::bodge_channel_slice(orig);

    cluster_graph_t inplace(orig);
    ClusterArrays::bodge_channel_slice_inplace(inplace);

    // The by-value form leaves its argument alone.
    CHECK(nodes_of(orig) == orig_nodes);
    CHECK(edges_of(orig) == orig_edges);

    // The bodge did something: new activity c-nodes and their edges.
    CHECK(boost::num_vertices(byval) > boost::num_vertices(orig));
    CHECK(boost::num_edges(byval) > boost::num_edges(orig));

    // Identical graphs, vertex and edge order included.
    CHECK(boost::num_vertices(inplace) == boost::num_vertices(byval));
    CHECK(boost::num_edges(inplace) == boost::num_edges(byval));
    CHECK(nodes_of(inplace) == nodes_of(byval));
    CHECK(edges_of(inplace) == edges_of(byval));
}

TEST_CASE("to_arrays leaves its input intact and is repeatable")
{
    const cluster_graph_t g = make_graph();
    ClusterArrays::node_array_set_t nas1, nas2;
    ClusterArrays::edge_array_set_t eas1, eas2;
    ClusterArrays::to_arrays(g, nas1, eas1);
    ClusterArrays::to_arrays(g, nas2, eas2);   // the input is const and not consumed
    REQUIRE(nas1.size() == nas2.size());
    for (const auto& [code, arr] : nas1) {
        REQUIRE(nas2.count(code));
        CHECK(arr == nas2.at(code));
    }
    REQUIRE(eas1.size() == eas2.size());
    for (const auto& [code, arr] : eas1) {
        REQUIRE(eas2.count(code));
        CHECK(arr == eas2.at(code));
    }
    // The activity array has one row per (slice, active channel): 5 + 2 here.
    REQUIRE(nas1.count('a'));
    CHECK(nas1.at('a').shape()[0] == 7);
}
