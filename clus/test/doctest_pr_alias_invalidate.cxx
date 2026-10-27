// doc sbnd_xin/docs/pr/31 F9 + pr_ml doc 26 sec 7.3: a second segment added between
// the same two vertices takes over the existing edge (the graph is boost::setS, no
// parallel edges).  Legacy: the first segment keeps a descriptor it believes valid
// while the edge now belongs to the second; removing the second frees the edge and
// the first segment's descriptor dangles -- PR::find_vertices on it dereferences a
// freed node bundle (S1 batch 2 PDVD cosmic1170: SIGSEGV in find_vertices via
// mcs_fill_kine, the long-muon chain held the aliased segment).
//
// Knob PR::set_alias_invalidate_old (TaggerCheckNeutrino pr_alias_invalidate_old,
// default false): at the takeover the first segment's descriptor is invalidated, so
// it reads as "not in the graph", which every caller of find_vertices already
// handles, and remove_segment on it no longer deletes the second segment's edge.

#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRSegment.h"

#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::Clus;

namespace {
    struct KnobReset {
        ~KnobReset() { PR::set_alias_invalidate_old(false); }
    };
}

TEST_CASE("pr alias: legacy (knob off) leaves the first segment holding the taken-over edge")
{
    KnobReset reset;
    PR::set_alias_invalidate_old(false);
    PR::Graph g;
    auto v1 = std::make_shared<PR::Vertex>();
    auto v2 = std::make_shared<PR::Vertex>();
    auto a = std::make_shared<PR::Segment>();
    auto b = std::make_shared<PR::Segment>();
    PR::add_segment(g, a, v1, v2);
    PR::add_segment(g, b, v1, v2);
    CHECK(boost::num_edges(g) == 1);
    CHECK(g[b->get_descriptor()].segment == b);
    CHECK(a->descriptor_valid());            // the hazard: a believes it is in the graph
    // removing a deletes b's edge (b's descriptor now dangles; not dereferenced here)
    CHECK(PR::remove_segment(g, a));
    CHECK(boost::num_edges(g) == 0);
}

TEST_CASE("pr alias: knob on invalidates the first segment at the takeover")
{
    KnobReset reset;
    PR::set_alias_invalidate_old(true);
    PR::Graph g;
    auto v1 = std::make_shared<PR::Vertex>();
    auto v2 = std::make_shared<PR::Vertex>();
    auto a = std::make_shared<PR::Segment>();
    auto b = std::make_shared<PR::Segment>();
    PR::add_segment(g, a, v1, v2);
    PR::add_segment(g, b, v1, v2);
    CHECK(boost::num_edges(g) == 1);
    CHECK(g[b->get_descriptor()].segment == b);
    CHECK_FALSE(a->descriptor_valid());
    CHECK(b->descriptor_valid());

    auto [x1, x2] = PR::find_vertices(g, a);    // "not in the graph": null pair
    CHECK(!x1);
    CHECK(!x2);
    CHECK_FALSE(PR::remove_segment(g, a));      // no longer deletes b's edge
    CHECK(boost::num_edges(g) == 1);

    // the S1 crash sequence: remove b, then ask for a's vertices
    CHECK(PR::remove_segment(g, b));
    CHECK(boost::num_edges(g) == 0);
    auto [y1, y2] = PR::find_vertices(g, a);
    CHECK(!y1);
    CHECK(!y2);
}

TEST_CASE("pr alias: knob on is a no-op when the same segment is re-added")
{
    KnobReset reset;
    PR::set_alias_invalidate_old(true);
    PR::Graph g;
    auto v1 = std::make_shared<PR::Vertex>();
    auto v2 = std::make_shared<PR::Vertex>();
    auto a = std::make_shared<PR::Segment>();
    PR::add_segment(g, a, v1, v2);
    PR::add_segment(g, a, v1, v2);
    CHECK(boost::num_edges(g) == 1);
    CHECK(a->descriptor_valid());
    CHECK(g[a->get_descriptor()].segment == a);
}
