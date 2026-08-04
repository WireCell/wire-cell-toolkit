// Pins the ordering contract that the doc pr/28 determinism rounds depend on.
//
// The pattern-recognition graph stores vertices and edges with boost::setS, so
// boost::vertices()/boost::edges() iterate in DESCRIPTOR POINTER order, which
// varies between runs of an identical program.  Every decision loop in clus/src
// is expected to go through PR::ordered_nodes() / PR::ordered_edges() instead,
// which sort by NodeBundle::index / EdgeBundle::index.
//
// That only buys determinism if those indices are UNIQUE: std::sort is not
// stable, and on ties it falls back to whatever order the input vector had --
// i.e. pointer order again.  A duplicate index would therefore make the whole
// sweep silently ineffective while every call site still "looks" converted.
// These cases pin uniqueness across the operations the PR chain actually
// performs, including vertex/segment removal (indices come from a monotonic
// counter in GraphBundle that never decrements, so re-adding after a removal
// must not reuse an index).
//
// Reverting PR::add_vertex/add_segment to reuse boost::num_vertices/num_edges
// as the index makes "indices survive removal" fail.

#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRTrajectoryView.h"

#include "WireCellUtil/doctest.h"

#include <set>

using namespace WireCell;
using namespace WireCell::Clus;

// Build a small chain graph: nvtx vertices joined by nvtx-1 segments.
static std::vector<PR::VertexPtr> make_chain(PR::Graph& g, size_t nvtx)
{
    std::vector<PR::VertexPtr> vtxs;
    for (size_t i = 0; i < nvtx; ++i) {
        vtxs.push_back(std::make_shared<PR::Vertex>());
    }
    for (size_t i = 0; i + 1 < nvtx; ++i) {
        auto seg = std::make_shared<PR::Segment>();
        PR::add_segment(g, seg, vtxs[i], vtxs[i + 1]);
    }
    return vtxs;
}

TEST_CASE("pr graph node and edge indices are unique")
{
    PR::Graph g;
    auto vtxs = make_chain(g, 8);

    auto nodes = PR::ordered_nodes(g);
    REQUIRE(nodes.size() == 8);
    std::set<size_t> nidx;
    for (const auto& vd : nodes) nidx.insert(g[vd].index);
    CHECK(nidx.size() == nodes.size());   // no duplicate sort key

    auto edges = PR::ordered_edges(g);
    REQUIRE(edges.size() == 7);
    std::set<size_t> eidx;
    for (const auto& ed : edges) eidx.insert(g[ed].index);
    CHECK(eidx.size() == edges.size());
}

TEST_CASE("pr graph ordered_nodes and ordered_edges come back sorted by index")
{
    PR::Graph g;
    make_chain(g, 8);

    auto nodes = PR::ordered_nodes(g);
    for (size_t i = 1; i < nodes.size(); ++i) {
        CHECK(g[nodes[i - 1]].index < g[nodes[i]].index);
    }

    auto edges = PR::ordered_edges(g);
    for (size_t i = 1; i < edges.size(); ++i) {
        CHECK(g[edges[i - 1]].index < g[edges[i]].index);
    }
}

TEST_CASE("pr graph indices stay unique after removal and re-add")
{
    PR::Graph g;
    auto vtxs = make_chain(g, 8);

    // Drop an interior vertex (which drops its two segments with it), then
    // grow the graph again -- the PR chain does exactly this when it merges
    // or splits segments.
    const size_t before = g[boost::graph_bundle].num_node_indices;
    PR::remove_vertex(g, vtxs[3]);

    auto extra = std::make_shared<PR::Vertex>();
    auto seg = std::make_shared<PR::Segment>();
    PR::add_segment(g, seg, vtxs[0], extra);

    // The counter must have advanced past the removed vertex's index rather
    // than reusing it.
    CHECK(g[boost::graph_bundle].num_node_indices > before);

    auto nodes = PR::ordered_nodes(g);
    std::set<size_t> nidx;
    for (const auto& vd : nodes) nidx.insert(g[vd].index);
    CHECK(nidx.size() == nodes.size());

    auto edges = PR::ordered_edges(g);
    std::set<size_t> eidx;
    for (const auto& ed : edges) eidx.insert(g[ed].index);
    CHECK(eidx.size() == edges.size());
}

// The "edge already existed" path of PR::add_segment (PRGraph.cxx:86-89):
// adding a second segment between an already-connected vertex pair overwrites
// the edge bundle's segment and COPIES the existing edge index into the new
// segment.  The displaced segment keeps that index, so two live Segment objects
// then report the same get_graph_index().
//
// Two separate things must hold, and they are not the same thing:
//   * EdgeBundle::index stays unique -- otherwise ordered_edges' sort key ties
//     and silently degrades to pointer order;
//   * Segment::get_graph_index() does NOT stay unique -- which is why a set or
//     map keyed on it cannot be used for identity lookups.  See the note on
//     eliminate_short_vertex_activities in NeutrinoPatternBase.h; making that
//     swap moved kine_reco_Enu on SBND evt 239794 by 1.2 GeV.
//
// The 1.2 GeV move is a measurement and stands.  Its ATTRIBUTION to this
// inherit path does NOT: probing that container at its call site on six events,
// including evt 239794, found no two segments sharing an index and no unindexed
// member (sbnd_xin/docs/pr/28 §14.8).  A more common aliasing source is the
// SIZE_MAX default of Segment::m_graph_index (PRSegment.h:153), shared by every
// segment not yet handed to PR::add_segment.  This case pins the non-uniqueness
// itself, which is what the rule rests on -- not the attribution.
TEST_CASE("pr graph edge index stays unique when a segment index is inherited")
{
    PR::Graph g;
    auto vtxs = make_chain(g, 4);

    auto first = std::make_shared<PR::Segment>();
    PR::add_segment(g, first, vtxs[0], vtxs[2]);      // new edge
    const size_t nedges_before = PR::ordered_edges(g).size();

    auto second = std::make_shared<PR::Segment>();
    PR::add_segment(g, second, vtxs[0], vtxs[2]);      // same pair: inherit path

    // No new edge was created, and the two segments now alias on index.
    auto edges = PR::ordered_edges(g);
    CHECK(edges.size() == nedges_before);
    CHECK(first->get_graph_index() == second->get_graph_index());

    // ...but the sort key ordered_edges actually uses is still pairwise unique.
    std::set<size_t> eidx;
    for (const auto& ed : edges) eidx.insert(g[ed].index);
    CHECK(eidx.size() == edges.size());
}

TEST_CASE("pr graph ordered_nodes and graph_nodes hold the same set")
{
    // graph_nodes() is the raw pointer-ordered walk -- it is NOT a determinism
    // helper despite the name.  It must at least agree with ordered_nodes() on
    // membership, so converting a call site changes only the order.
    PR::Graph g;
    make_chain(g, 8);

    auto raw = PR::graph_nodes(g);
    auto ord = PR::ordered_nodes(g);
    REQUIRE(raw.size() == ord.size());

    std::set<PR::node_descriptor> a(raw.begin(), raw.end());
    std::set<PR::node_descriptor> b(ord.begin(), ord.end());
    CHECK(a == b);
}

// --- round 9 (doc pr/28 §15) ------------------------------------------------

// TrajectoryView::edges()/nodes() are std::unordered_set.  The edge hash is
// over the two node_descriptors, and with boost::setS vertices a
// node_descriptor is a HEAP ADDRESS -- so the bucket walk differs between runs
// of the identical binary.  Shower::get_total_length() and
// Shower::calculate_kinematics() accumulated FP over that walk, which moved
// showers[]/kine_dQdx by one ULP (1314.124434586102 -> ...103) between two
// `setarch -R` runs of SBND evt 388.
//
// Reverting ordered_edges(view, g) back to view.edges() in PRShower.cxx does
// NOT fail this case -- it pins the helper, not the call sites, because a
// single-process test cannot re-randomize the heap.  What it does pin is that
// the helper the call sites were converted to actually sorts by the stable
// index and covers the view exactly.  Deleting the std::sort from
// PR::ordered_edges(view, g) fails it.
TEST_CASE("pr trajectory view ordered_edges sorts an address-hashed set by index")
{
    PR::Graph g;
    auto vtxs = make_chain(g, 6);

    PR::TrajectoryView view(g);
    for (const auto& v : vtxs) view.add_vertex(v);
    for (const auto& ed : PR::ordered_edges(g)) view.add_segment(g[ed].segment);

    auto ord = PR::ordered_edges(view, g);
    REQUIRE(ord.size() == view.edges().size());

    // Sorted, strictly increasing on the stable EdgeBundle::index.
    for (size_t i = 0; i + 1 < ord.size(); ++i) {
        CHECK(g[ord[i]].index < g[ord[i + 1]].index);
    }

    // Same membership as the raw walk -- conversion changes order only.
    PR::edge_unordered_set raw(view.edges().begin(), view.edges().end(),
                               0, PR::EdgeDescriptorHash(g), PR::EdgeDescriptorEqual(g));
    PR::edge_unordered_set from_ord(ord.begin(), ord.end(),
                                    0, PR::EdgeDescriptorHash(g), PR::EdgeDescriptorEqual(g));
    CHECK(raw.size() == from_ord.size());
    for (const auto& ed : ord) CHECK(raw.count(ed) == 1);

    auto onodes = PR::ordered_nodes(view, g);
    REQUIRE(onodes.size() == view.nodes().size());
    for (size_t i = 0; i + 1 < onodes.size(); ++i) {
        CHECK(g[onodes[i]].index < g[onodes[i + 1]].index);
    }
}

// The SIZE_MAX hazard the round-9 warning guards.  Every PR::Segment that has
// never been handed to PR::add_segment carries PR::kUnindexed, so under
// SegmentIndexCmp they all compare EQUAL -- an IndexedSegmentSet keeps one and
// silently drops the rest.  This is the reason the guard exists; the case pins
// the collapse so that "index-keyed sets are safe" can never be assumed.
//
// Removing PR::kUnindexed / changing Segment's default m_graph_index fails it.
TEST_CASE("pr unindexed segments collapse in an index-ordered set")
{
    auto a = std::make_shared<PR::Segment>();
    auto b = std::make_shared<PR::Segment>();
    auto c = std::make_shared<PR::Segment>();
    CHECK(a->get_graph_index() == PR::kUnindexed);
    CHECK(b->get_graph_index() == PR::kUnindexed);

    PR::IndexedSegmentSet unindexed{a, b, c};
    CHECK(unindexed.size() == 1);            // three distinct objects, one survivor
    CHECK(unindexed.count(b) == 1);          // ...and find() answers for any of them

    // Once they are in a graph the indices are distinct and nothing collapses.
    PR::Graph g;
    auto v0 = std::make_shared<PR::Vertex>();
    auto v1 = std::make_shared<PR::Vertex>();
    auto v2 = std::make_shared<PR::Vertex>();
    auto v3 = std::make_shared<PR::Vertex>();
    PR::add_segment(g, a, v0, v1);
    PR::add_segment(g, b, v1, v2);
    PR::add_segment(g, c, v2, v3);
    PR::IndexedSegmentSet indexed{a, b, c};
    CHECK(indexed.size() == 3);
}

// ---------------------------------------------------------------------------
// doc sbnd_xin/docs/pr/30 §11, P8 -- the endpoint-consistency check that the
// port dropped (prototype NeutrinoID::add_proto_connection,
// NeutrinoID_proto_vertex.h:1952-1956, "Error! Vertex and Segment does not
// match").
//
// Two properties are pinned:
//   1. The check only judges a connection actually being MADE.  A re-call on a
//      segment already in the graph is a no-op and must NOT be counted -- the
//      first draft counted those and reported 108 "mismatches" on 48 events
//      that were all no-op re-entries.
//   2. Strict mode withholds the EDGE but still records the vertices, which is
//      what the prototype does (it returns before touching the connection maps
//      but any vertex it had already seen stays in proto_vertices).
//
// Revert proof: moving the check back above the descriptor_valid() early
// return makes "re-adding a segment is not judged" fail; dropping the
// g_graph_endpoint_policy.strict branch makes "strict mode withholds the
// edge" fail.
static PR::SegmentPtr make_seg_with_ends(const WireCell::Point& a, const WireCell::Point& b)
{
    auto seg = std::make_shared<PR::Segment>();
    std::vector<PR::WCPoint> wcpts(2);
    wcpts[0].point = a;
    wcpts[1].point = b;
    seg->wcpts(wcpts);
    return seg;
}
static PR::VertexPtr make_vtx_at(const WireCell::Point& p)
{
    auto v = std::make_shared<PR::Vertex>();
    v->wcpt().point = p;
    return v;
}

TEST_CASE("pr30 P8 inconsistent endpoints are counted only when a connection is made")
{
    const WireCell::Point a(0, 0, 0), b(10 * units::cm, 0, 0);
    const WireCell::Point far_off(0, 5 * units::cm, 0);   // 5 cm from either end

    PR::g_graph_endpoint_policy.strict = false;
    PR::g_graph_endpoint_policy.tol = 0.3 * units::cm;

    PR::Graph g;
    auto seg = make_seg_with_ends(a, b);
    auto v1 = make_vtx_at(a);
    auto v2 = make_vtx_at(far_off);          // NOT at either end

    const auto n0 = PR::g_port_audit.endpoint_mismatch.load();
    const auto r0 = PR::g_port_audit.add_segment_reentry.load();

    CHECK(PR::add_segment(g, seg, v1, v2) == true);
    CHECK(PR::g_port_audit.endpoint_mismatch.load() == n0 + 1);   // judged once
    CHECK(PR::g_port_audit.add_segment_reentry.load() == r0);     // first call is not a re-entry
    CHECK(boost::num_edges(g) == 1u);        // not strict => the edge is still made

    // Re-call on the SAME segment: it already has a descriptor, so this is a
    // no-op and must not be judged again.
    PR::add_segment(g, seg, v1, v2);
    CHECK(PR::g_port_audit.endpoint_mismatch.load() == n0 + 1);   // still 1, not 2
    CHECK(PR::g_port_audit.add_segment_reentry.load() == r0 + 1);
    CHECK(boost::num_edges(g) == 1u);

    // A consistent pair is never judged.
    PR::Graph g2;
    auto seg2 = make_seg_with_ends(a, b);
    CHECK(PR::add_segment(g2, seg2, make_vtx_at(a), make_vtx_at(b)) == true);
    CHECK(PR::g_port_audit.endpoint_mismatch.load() == n0 + 1);
}

TEST_CASE("pr30 P8 strict mode withholds the edge but keeps the vertices")
{
    const WireCell::Point a(0, 0, 0), b(10 * units::cm, 0, 0);
    const WireCell::Point far_off(0, 5 * units::cm, 0);

    PR::g_graph_endpoint_policy.strict = true;
    PR::g_graph_endpoint_policy.tol = 0.3 * units::cm;

    PR::Graph g;
    auto seg = make_seg_with_ends(a, b);
    const auto ref0 = PR::g_port_audit.endpoint_refused.load();

    PR::add_segment(g, seg, make_vtx_at(a), make_vtx_at(far_off));

    CHECK(PR::g_port_audit.endpoint_refused.load() == ref0 + 1);
    CHECK(boost::num_edges(g) == 0u);        // the connection is NOT made
    CHECK(boost::num_vertices(g) == 2u);     // ... but the vertices are recorded
    CHECK(seg->descriptor_valid() == false); // and the segment stays detached

    // A consistent pair still connects with strict on.
    PR::Graph g2;
    auto seg2 = make_seg_with_ends(a, b);
    PR::add_segment(g2, seg2, make_vtx_at(a), make_vtx_at(b));
    CHECK(boost::num_edges(g2) == 1u);
    CHECK(PR::g_port_audit.endpoint_refused.load() == ref0 + 1);

    PR::g_graph_endpoint_policy.strict = false;   // restore for other cases
}
