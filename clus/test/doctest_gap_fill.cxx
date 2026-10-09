// pdvd doc 131 -- unit tests for the pure decision core of the retiler's "gapfill" mode (GapFill.h).

#include "WireCellClus/GapFill.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell::Clus::GapFill;

TEST_CASE("gap fill: covers with margins")
{
    Box b{100, 104, {{{10, 12}, {20, 23}, {30, 31}}}};
    CHECK(covers(b, 100, {10, 20, 30}, 0, 0));
    CHECK(covers(b, 103, {11, 22, 30}, 0, 0));
    CHECK_FALSE(covers(b, 104, {11, 22, 30}, 0, 0));      // half-open in slice
    CHECK(covers(b, 104, {11, 22, 30}, 0, 4));            // one slice of margin
    CHECK_FALSE(covers(b, 100, {12, 20, 30}, 0, 0));      // half-open in wire
    CHECK(covers(b, 100, {12, 20, 30}, 1, 0));            // one wire of margin
    CHECK(covers(b, 100, {9, 23, 31}, 1, 0));
    CHECK_FALSE(covers(b, 100, {8, 20, 30}, 1, 0));
}

namespace {
    std::vector<PathPt> path(int n, int slice_step)
    {
        std::vector<PathPt> p(n);
        for (int i = 0; i < n; ++i) {
            p[i].slice = 1000 + (i / slice_step) * 4;   // 4-tick slices, slice_step samples per slice
            p[i].wire = {i, i, i};
            p[i].s = 0.3 * i;
            p[i].on_face = true;
        }
        return p;
    }
}

TEST_CASE("gap fill: find_gaps")
{
    auto p = path(40, 2);   // 20 slices, 0.3 cm per sample
    std::vector<bool> cov(40, true);
    SUBCASE("no gap")
    {
        CHECK(find_gaps(p, cov, 2, 0).empty());
    }
    SUBCASE("an interior gap of 10 samples = 5 slices, 2.7 cm")
    {
        for (int i = 10; i < 20; ++i) cov[i] = false;
        auto g = find_gaps(p, cov, 2, 0);
        REQUIRE(g.size() == 1);
        CHECK(g[0].first == 10);
        CHECK(g[0].last == 19);
        CHECK(g[0].nslices == 5);
        CHECK(g[0].length == doctest::Approx(2.7));
        CHECK(find_gaps(p, cov, 6, 0).empty());        // too few slices
        CHECK(find_gaps(p, cov, 2, 2.0).empty());      // longer than the cap
        CHECK(find_gaps(p, cov, 2, 3.0).size() == 1);
    }
    SUBCASE("a run touching the path's end is not a gap")
    {
        for (int i = 30; i < 40; ++i) cov[i] = false;
        CHECK(find_gaps(p, cov, 2, 0).empty());
        for (int i = 0; i < 6; ++i) cov[i] = false;
        CHECK(find_gaps(p, cov, 2, 0).empty());
    }
    SUBCASE("a one-slice hole is below min_slices 2")
    {
        cov[20] = false; cov[21] = false;   // both samples of one slice
        CHECK(find_gaps(p, cov, 2, 0).empty());
        CHECK(find_gaps(p, cov, 1, 0).size() == 1);
    }
    SUBCASE("samples on another face break the run")
    {
        for (int i = 10; i < 20; ++i) cov[i] = false;
        p[15].on_face = false;
        CHECK(find_gaps(p, cov, 2, 0).empty());   // neither half is bounded by covered samples on both sides
    }
    SUBCASE("two gaps")
    {
        for (int i = 4; i < 8; ++i) cov[i] = false;
        for (int i = 30; i < 36; ++i) cov[i] = false;
        auto g = find_gaps(p, cov, 2, 0);
        REQUIRE(g.size() == 2);
        CHECK(g[0].first == 4); CHECK(g[0].last == 7);
        CHECK(g[1].first == 30); CHECK(g[1].last == 35);
    }
}

TEST_CASE("gap fill: owner map, residual and the corridor rule")
{
    OwnerMap om;
    om.add(0, 1000, 10, 14, 400.0);   // 100 per wire
    om.add(0, 1000, 12, 13, 50.0);    // overlaps wire 12
    CHECK(om.is_other(0, 1000, 10));
    CHECK_FALSE(om.is_other(0, 1000, 14));
    CHECK_FALSE(om.is_other(1, 1000, 10));
    CHECK_FALSE(om.is_other(0, 1004, 10));
    CHECK(om.other_share(0, 1000, 11) == doctest::Approx(100.0));
    CHECK(om.other_share(0, 1000, 12) == doctest::Approx(150.0));
    CHECK(om.other_share(0, 1000, 20) == 0.0);
    CHECK(om.size() == 4);
    om.add(2, 1000, 5, 5, 10.0);      // empty range: nothing
    om.add(2, 1000, 5, 6, -3.0);      // negative charge counts as 0 but marks the wire
    CHECK(om.size() == 5);
    CHECK(om.is_other(2, 1000, 5));
    CHECK(om.other_share(2, 1000, 5) == 0.0);

    CHECK(residual(300.0, 100.0) == doctest::Approx(200.0));
    CHECK(residual(50.0, 100.0) == 0.0);

    CHECK(missing_planes({3, 2, 1}) == 0);
    CHECK(missing_planes({0, 2, 1}) == 1);
    CHECK(missing_planes({0, 0, 1}) == 2);
    CHECK(missing_planes({0, 0, 0}) == 3);
}
