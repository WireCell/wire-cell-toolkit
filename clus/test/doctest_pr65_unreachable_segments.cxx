// Pins PR::unreachable_segments (doc sbnd_xin/docs/pr/65 round 3), the
// helper behind the shower_absorb_unreachable_main knob: the shower
// absorbers' cluster()==main_cluster guards encode "already claimed by the
// main_vertex graph walk", an invariant other_seg_keep_isolated (doc pr/54)
// broke by adding kept residual segments as DISCONNECTED components of the
// main cluster's PR graph.  The knob relaxes those guards to reachability,
// and this helper is the reachability test.
//
// Contract pinned here:
//  - connected graph, any root  => empty set;
//  - disconnected components    => exactly the segments of every component
//    not holding the root, regardless of which component holds the root;
//  - null / never-added root    => every segment in the graph.
//
// Reverting the helper to "reachable iff same cluster id" (the pre-pr/65
// reading of the guards) makes the two-component case fail: both components
// here belong to no cluster at all, yet one is unreachable.

#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRSegment.h"

#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::Clus;

namespace {

struct Chain {
    std::vector<PR::VertexPtr> vtxs;
    std::vector<PR::SegmentPtr> segs;
};

// Build a chain graph component: nvtx vertices joined by nvtx-1 segments.
Chain make_chain(PR::Graph& g, size_t nvtx)
{
    Chain c;
    for (size_t i = 0; i < nvtx; ++i) {
        c.vtxs.push_back(std::make_shared<PR::Vertex>());
    }
    for (size_t i = 0; i + 1 < nvtx; ++i) {
        auto seg = std::make_shared<PR::Segment>();
        PR::add_segment(g, seg, c.vtxs[i], c.vtxs[i + 1]);
        c.segs.push_back(seg);
    }
    return c;
}

}  // namespace

TEST_CASE("pr65 unreachable_segments: connected graph is fully reachable")
{
    PR::Graph g;
    auto chain = make_chain(g, 6);

    // From either end and from the middle: nothing is unreachable.
    CHECK(PR::unreachable_segments(g, chain.vtxs.front()).empty());
    CHECK(PR::unreachable_segments(g, chain.vtxs.back()).empty());
    CHECK(PR::unreachable_segments(g, chain.vtxs[3]).empty());
}

TEST_CASE("pr65 unreachable_segments: disconnected component is returned exactly")
{
    PR::Graph g;
    auto main_chain = make_chain(g, 5);   // 4 segments, holds the root
    auto orphan     = make_chain(g, 3);   // 2 segments, disconnected
    auto orphan2    = make_chain(g, 2);   // 1 segment, disconnected

    auto unreachable = PR::unreachable_segments(g, main_chain.vtxs[0]);
    REQUIRE(unreachable.size() == orphan.segs.size() + orphan2.segs.size());
    for (const auto& seg : orphan.segs)  CHECK(unreachable.count(seg) == 1);
    for (const auto& seg : orphan2.segs) CHECK(unreachable.count(seg) == 1);
    for (const auto& seg : main_chain.segs) CHECK(unreachable.count(seg) == 0);

    // Root in the small component: the big one is now the unreachable side.
    auto unreachable2 = PR::unreachable_segments(g, orphan.vtxs[1]);
    REQUIRE(unreachable2.size() == main_chain.segs.size() + orphan2.segs.size());
    for (const auto& seg : main_chain.segs) CHECK(unreachable2.count(seg) == 1);
    for (const auto& seg : orphan.segs)     CHECK(unreachable2.count(seg) == 0);
}

TEST_CASE("pr65 unreachable_segments: null or never-added root returns everything")
{
    PR::Graph g;
    auto chain = make_chain(g, 4);

    auto from_null = PR::unreachable_segments(g, nullptr);
    CHECK(from_null.size() == chain.segs.size());

    auto loose = std::make_shared<PR::Vertex>();  // never added to g
    auto from_loose = PR::unreachable_segments(g, loose);
    CHECK(from_loose.size() == chain.segs.size());
}
