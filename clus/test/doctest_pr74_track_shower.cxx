// doc sbnd_xin/docs/pr/74 round 2 -- the two new predicates behind the P1
// cascade guard (shower_in_cascade_guard) and the P2 Michel-terminal check
// (michel_stem_michel_check).
//
// P1 contract: segment_shower_in_cascade_vetoed(seg, mip_med, max_len,
// mip_hi) is true iff the segment is BOTH long (> max_len) AND MIP-like
// (median dQ/dx < mip_hi * mip_med).  A zero/absent median never vetoes --
// "no evidence" is not "MIP-like evidence" (the pr/40
// segment_dqdx_spares_electron_reclass convention).  The motivating case is
// SBND 18345-53361 seg 27004: 113.9 cm at 1.02x MIP relabelled e- by
// examine_direction's unconditional flag_shower_in cascade.
//
// P2 contract: segment_far_subtree_track_length(graph, far_vtx, stem, cap)
// sums the track length reachable beyond far_vtx without crossing the stem,
// early-exiting once cap is exceeded.  A genuine Michel electron is
// terminal (a few cm, nothing downstream); 90055's "Michel" sibling heads
// the 2020 MeV shower's ~155 cm tree.

#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegmentFunctions.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::Clus::PR;

namespace {

const double STEP = 0.6 * units::cm;
const double MIP_MED = 43000 / units::cm;   // production mip_dqdx_median

// Straight segment from `a` to `b` with a flat dQ/dx of `dqdx_ratio` x MIP,
// fits spaced at the production-like 0.6 cm.
SegmentPtr make_flat_track(Graph& g, VertexPtr v1, VertexPtr v2, double dqdx_ratio)
{
    auto seg = make_segment(g, v1, v2);
    const auto d = v2->wcpt().point - v1->wcpt().point;
    const double L = d.magnitude();
    const int n = static_cast<int>(L / STEP) + 1;
    for (int i = 0; i < n; i++) {
        Fit f;
        const double frac = (n > 1) ? double(i) / (n - 1) : 0.0;
        f.point = v1->wcpt().point + d * frac;
        f.index = i;
        f.dx = STEP;
        f.dQ = dqdx_ratio * MIP_MED * STEP;
        seg->fits().push_back(f);
    }
    return seg;
}

VertexPtr make_vtx(Graph& g, double x, double y, double z)
{
    auto v = make_vertex(g);
    v->wcpt().point = Point(x, y, z);
    return v;
}

}  // namespace

TEST_CASE("pr74 cascade veto: long AND MIP-like segment is vetoed")
{
    Graph g;
    auto v1 = make_vtx(g, 0, 0, 0);
    auto v2 = make_vtx(g, 0, 0, 60 * units::cm);
    auto seg = make_flat_track(g, v1, v2, 1.02);   // the 53361 shape
    CHECK(segment_shower_in_cascade_vetoed(seg, MIP_MED, 40 * units::cm, 1.3));
}

TEST_CASE("pr74 cascade veto: short MIP segment is spared (length conjunct)")
{
    Graph g;
    auto v1 = make_vtx(g, 0, 0, 0);
    auto v2 = make_vtx(g, 0, 0, 15 * units::cm);
    auto seg = make_flat_track(g, v1, v2, 0.9);    // 90055's shower-body shape
    CHECK_FALSE(segment_shower_in_cascade_vetoed(seg, MIP_MED, 40 * units::cm, 1.3));
}

TEST_CASE("pr74 cascade veto: long but charge-hot segment is spared")
{
    Graph g;
    auto v1 = make_vtx(g, 0, 0, 0);
    auto v2 = make_vtx(g, 0, 0, 60 * units::cm);
    auto seg = make_flat_track(g, v1, v2, 2.0);
    CHECK_FALSE(segment_shower_in_cascade_vetoed(seg, MIP_MED, 40 * units::cm, 1.3));
}

TEST_CASE("pr74 cascade veto: no dQ/dx samples never vetoes")
{
    Graph g;
    auto v1 = make_vtx(g, 0, 0, 0);
    auto v2 = make_vtx(g, 0, 0, 60 * units::cm);
    auto seg = make_segment(g, v1, v2);            // no fits at all
    CHECK_FALSE(segment_shower_in_cascade_vetoed(seg, MIP_MED, 40 * units::cm, 1.3));
    CHECK_FALSE(segment_shower_in_cascade_vetoed(nullptr, MIP_MED, 40 * units::cm, 1.3));
}

TEST_CASE("pr74 far subtree: terminal Michel-sized sibling stays under cap")
{
    // main -- stem -- far -- michel(4cm).  Beyond far: 4 cm total.
    Graph g;
    auto vm = make_vtx(g, 0, 0, 0);
    auto vf = make_vtx(g, 0, 0, 20 * units::cm);
    auto ve = make_vtx(g, 0, 4 * units::cm, 20 * units::cm);
    auto stem = make_flat_track(g, vm, vf, 1.8);
    make_flat_track(g, vf, ve, 0.9);
    const double far_len = segment_far_subtree_track_length(g, vf, stem, 40 * units::cm);
    CHECK(far_len < 10 * units::cm);
    CHECK(far_len > 2 * units::cm);
}

TEST_CASE("pr74 far subtree: shower-trunk sibling blows past the cap")
{
    // main -- stem -- far -- chain of 5 x 15 cm = 75 cm downstream.
    Graph g;
    auto vm = make_vtx(g, 0, 0, 0);
    auto vf = make_vtx(g, 0, 0, 20 * units::cm);
    auto stem = make_flat_track(g, vm, vf, 1.8);
    VertexPtr prev = vf;
    for (int i = 0; i < 5; i++) {
        auto nxt = make_vtx(g, 0, (i + 1) * 2.0 * units::cm, (20 + 15 * (i + 1)) * units::cm);
        make_flat_track(g, prev, nxt, 0.9);
        prev = nxt;
    }
    const double far_len = segment_far_subtree_track_length(g, vf, stem, 40 * units::cm);
    CHECK(far_len > 40 * units::cm);
}

TEST_CASE("pr74 far subtree: never walks back through the stem")
{
    // A long tail hangs off the MAIN vertex; beyond far there is nothing.
    Graph g;
    auto vm = make_vtx(g, 0, 0, 0);
    auto vf = make_vtx(g, 0, 0, 20 * units::cm);
    auto vt = make_vtx(g, 0, 50 * units::cm, 0);
    auto stem = make_flat_track(g, vm, vf, 1.8);
    make_flat_track(g, vm, vt, 0.9);               // 50 cm tail on the main side
    const double far_len = segment_far_subtree_track_length(g, vf, stem, 40 * units::cm);
    CHECK(far_len == 0.0);
}

// ---------------------------------------------------------------------------
// doc sbnd_xin/docs/pr/74 round 3 Q2 -- PR::reachable_vertices.
//
// The vertex-level complement of pr/65's unreachable_segments.  K5 anchors a
// graph-unreachable component onto "the nearest vertex", and without this
// filter the nearest vertex is the component's OWN endpoint at distance 0:
// 18255-142421 anchored seg 7013 to its far end, collapsing the pseudo-gamma
// to zero length and reversing the reconstructed electron.
//
// Deleting the `if (!reachable_vtxs.count(vtx)) continue;` guard in
// NeutrinoShowerClustering.cxx makes the second case below meaningless --
// which is exactly the round-2 bug.
// ---------------------------------------------------------------------------

#include "WireCellClus/PRGraph.h"

namespace {

// Chain component: nvtx vertices joined by nvtx-1 segments.
std::vector<WireCell::Clus::PR::VertexPtr>
pr74r3_make_chain(WireCell::Clus::PR::Graph& g, size_t nvtx)
{
    namespace PR = WireCell::Clus::PR;
    std::vector<PR::VertexPtr> vtxs;
    for (size_t i = 0; i < nvtx; ++i) vtxs.push_back(std::make_shared<PR::Vertex>());
    for (size_t i = 0; i + 1 < nvtx; ++i) {
        PR::add_segment(g, std::make_shared<PR::Segment>(), vtxs[i], vtxs[i + 1]);
    }
    return vtxs;
}

}  // namespace

TEST_CASE("pr74r3 reachable_vertices: connected graph yields every vertex")
{
    WireCell::Clus::PR::Graph g;
    auto chain = pr74r3_make_chain(g, 6);

    for (auto root : {chain.front(), chain[3], chain.back()}) {
        auto r = WireCell::Clus::PR::reachable_vertices(g, root);
        CHECK(r.size() == chain.size());
        for (const auto& v : chain) CHECK(r.count(v) == 1);
    }
}

TEST_CASE("pr74r3 reachable_vertices: a disconnected component is excluded")
{
    WireCell::Clus::PR::Graph g;
    auto main_comp = pr74r3_make_chain(g, 4);   // holds the root
    auto orphan    = pr74r3_make_chain(g, 3);   // the K5 promotion candidate

    auto r = WireCell::Clus::PR::reachable_vertices(g, main_comp.front());

    // The root's own component is reachable ...
    CHECK(r.size() == main_comp.size());
    for (const auto& v : main_comp) CHECK(r.count(v) == 1);
    // ... and NONE of the orphan component's endpoints is -- the guard that
    // stops K5 anchoring the promoted shower to its own far end.
    for (const auto& v : orphan) CHECK(r.count(v) == 0);

    // Symmetry with pr/65: the same split, seen segment-wise.
    CHECK(WireCell::Clus::PR::unreachable_segments(g, main_comp.front()).size() == 2);
}

TEST_CASE("pr74r3 reachable_vertices: null root reaches nothing")
{
    WireCell::Clus::PR::Graph g;
    pr74r3_make_chain(g, 4);
    CHECK(WireCell::Clus::PR::reachable_vertices(g, nullptr).empty());
}
