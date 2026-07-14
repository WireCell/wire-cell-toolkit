#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegmentFunctions.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::PR;

TEST_CASE("clus pr segment basic") {
    PR::Segment seg;

    CHECK(! seg.descriptor_valid());
}

// F7 regression: find_vertices must not dereference wcpts().front() when wcpts
// is empty.  Before the fix, calling find_vertices on a segment built only from
// fits (no wcpts) caused undefined behaviour.  After the fix the code falls back
// to fits().front() if wcpts is empty.
TEST_CASE("clus pr segment find_vertices empty wcpts fallback") {
    Graph g;

    // Two vertices at known positions along X.
    auto vtx1 = make_vertex(g);
    vtx1->wcpt().point = Point(0, 0, 0);

    auto vtx2 = make_vertex(g);
    vtx2->wcpt().point = Point(10*units::cm, 0, 0);

    // Segment with fits but NO wcpts.
    auto seg = make_segment(g, vtx1, vtx2);
    REQUIRE(seg->wcpts().empty());

    // Add a fit closer to vtx1 so find_vertices can order the pair.
    Fit f;
    f.point = Point(1*units::cm, 0, 0);
    f.index = 0;
    seg->fits().push_back(f);

    // Must not crash and must return both vertices non-null.
    auto [va, vb] = find_vertices(g, seg);
    REQUIRE(va != nullptr);
    REQUIRE(vb != nullptr);

    // va should be closer to the fit point (i.e. vtx1 at x=0 cm).
    CHECK(va.get() == vtx1.get());
    CHECK(vb.get() == vtx2.get());
}

// F5 regression: clear_fit(nullptr) must not attempt to populate paf (no DV),
// so every fit's paf must remain at the default {-1,-1}.
// (Testing with a non-null DV requires a Facade::Cluster attached to the segment
// because create_segment_fit_point_cloud internally accesses the cluster's
// grouping.  That is covered by integration tests.)
TEST_CASE("clus pr segment clear_fit paf with null dv") {
    // Segment must be owned by a shared_ptr for shared_from_this().
    auto seg = std::make_shared<Segment>();

    // Populate wcpts at two known positions.
    WCPoint wp1; wp1.point = Point(0, 0, 0);
    WCPoint wp2; wp2.point = Point(5*units::cm, 0, 0);
    seg->wcpts({wp1, wp2});

    // clear_fit(nullptr) should reset fits to match wcpts but leave paf at default.
    seg->clear_fit(nullptr);

    REQUIRE(seg->fits().size() == 2);
    for (const auto& fit : seg->fits()) {
        CHECK(fit.paf.first  == -1);
        CHECK(fit.paf.second == -1);
    }

    // Points must be copied from wcpts.
    CHECK(seg->fits()[0].point == Point(0, 0, 0));
    CHECK(seg->fits()[1].point == Point(5*units::cm, 0, 0));

    // All other numeric fields must be reset to their defaults.
    for (const auto& fit : seg->fits()) {
        CHECK(fit.dQ           == doctest::Approx(-1.0));
        CHECK(fit.dx           == doctest::Approx( 0.0));
        CHECK(fit.reduced_chi2 == doctest::Approx(-1.0));
        CHECK(fit.index        == -1);
        CHECK(fit.flag_fix     == false);
    }
}

// F6 regression: both segment_track_direct_length and segment_track_max_deviation
// must interpret n2 = -1 as "last fit index" (fits.size()-1), not as 0.
// Before the fix, segment_track_direct_length clamped n2<0 to 0, so calling
// (seg, 0, -1) on a 5-fit segment measured from index 0 to 0 (length 0).
TEST_CASE("clus pr segment track functions n2 minus1 clamp") {
    // Build a segment with 5 collinear fits at x = 0..4 cm, y=z=0.
    auto seg = std::make_shared<Segment>();
    std::vector<Fit> fits(5);
    for (int i = 0; i < 5; ++i) {
        fits[i].point = Point(i * units::cm, 0, 0);
        fits[i].index = i;
    }
    seg->fits(fits);

    const double full_len = 4*units::cm;   // x=0 to x=4 cm

    // -- segment_track_direct_length --

    // Default (n1=-1, n2=-1): full segment.
    CHECK(segment_track_direct_length(seg) == doctest::Approx(full_len));

    // (0, -1) must equal (0, 4): length from first to last fit.
    CHECK(segment_track_direct_length(seg, 0, -1) == doctest::Approx(full_len));

    // (2, -1) must equal (2, 4): partial length from index 2 to end.
    CHECK(segment_track_direct_length(seg, 2, -1) == doctest::Approx(2*units::cm));

    // Explicit endpoints must agree.
    CHECK(segment_track_direct_length(seg, 0, 4) ==
          doctest::Approx(segment_track_direct_length(seg, 0, -1)));

    // -- segment_track_max_deviation --
    // Fits are collinear so max deviation is 0 in all sub-ranges.

    CHECK(segment_track_max_deviation(seg)        == doctest::Approx(0.0));
    CHECK(segment_track_max_deviation(seg, 0, -1) == doctest::Approx(0.0));
    CHECK(segment_track_max_deviation(seg, 2, -1) == doctest::Approx(0.0));

    // Explicit endpoints must agree.
    CHECK(segment_track_max_deviation(seg, 0, 4) ==
          doctest::Approx(segment_track_max_deviation(seg, 0, -1)));
}
