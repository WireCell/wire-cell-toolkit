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
