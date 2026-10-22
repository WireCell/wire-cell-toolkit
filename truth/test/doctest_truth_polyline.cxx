/** DepoPolyline on a zig-zag: ordering by time, arc length, nearest
 *  point and tangent sign, a kink at the bend (and none on a straight
 *  line), end keypoints, heat, dedx, and a long polyline where the
 *  k-d tree candidate search must agree with brute force.
 */

#include "WireCellTruth/DepoPolyline.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <random>

using namespace WireCell;
using Truth::DepoPolyline;

TEST_CASE("polyline orders by time and measures a straight line")
{
    DepoPolyline pl;
    // added out of order; 1 mm steps along +z with energy 0.1 per vertex
    for (int i = 9; i >= 0; --i) pl.add(Point(0, 0, i * units::mm), i * 1.0 * units::us, 0.1);
    pl.build();
    CHECK(pl.nvertices() == 10);
    CHECK(pl.nsegments() == 9);
    CHECK(pl.length() == doctest::Approx(9 * units::mm));
    CHECK(pl.start().z() == doctest::Approx(0));
    CHECK(pl.end().z() == doctest::Approx(9 * units::mm));

    auto h = pl.nearest(Point(0.5 * units::mm, 0, 4.25 * units::mm));
    REQUIRE(h.valid);
    CHECK(h.dist == doctest::Approx(0.5 * units::mm));
    CHECK(h.seg == 4);
    CHECK(h.s == doctest::Approx(4.25 * units::mm));
    CHECK(h.tangent.z() == doctest::Approx(1.0));
    CHECK(h.closest.z() == doctest::Approx(4.25 * units::mm));
    CHECK(h.dedx == doctest::Approx(0.1 / units::mm));

    // beyond the end: clamped to the last vertex
    auto e = pl.nearest(Point(0, 0, 12 * units::mm));
    CHECK(e.dist == doctest::Approx(3 * units::mm));
    CHECK(e.s == doctest::Approx(pl.length()));

    auto kps = pl.keypoints(2 * units::mm, 15 * units::deg);
    REQUIRE(kps.size() == 2);   // start and end only
    CHECK(kps.front().kind == DepoPolyline::KeyKind::start);
    CHECK(kps.back().kind == DepoPolyline::KeyKind::end);
    CHECK(DepoPolyline::heat(Point(0, 0, 0), kps, DepoPolyline::KeyKind::start, 1 * units::cm) == doctest::Approx(1.0));
    CHECK(DepoPolyline::heat(Point(0, 0, 9 * units::mm), kps, DepoPolyline::KeyKind::kink, 1 * units::cm) == 0.0);
    const double d = 10 * units::mm, sig = 10 * units::mm;
    CHECK(DepoPolyline::heat(Point(d, 0, 0), kps, DepoPolyline::KeyKind::start, sig) == doctest::Approx(std::exp(-0.5)));
}

TEST_CASE("polyline finds the kink of a zig-zag and the tangent flips with it")
{
    DepoPolyline pl;
    // 5 cm along +z then 5 cm along +y, 1 mm steps
    int i = 0;
    for (; i <= 50; ++i) pl.add(Point(0, 0, i * units::mm), i, 0);
    for (int j = 1; j <= 50; ++j, ++i) pl.add(Point(0, j * units::mm, 50 * units::mm), i, 0);
    pl.build();
    CHECK(pl.length() == doctest::Approx(100 * units::mm));

    auto before = pl.nearest(Point(0, 0, 20 * units::mm));
    CHECK(before.tangent.z() == doctest::Approx(1.0));
    auto after = pl.nearest(Point(0, 20 * units::mm, 50 * units::mm));
    CHECK(after.tangent.y() == doctest::Approx(1.0));
    CHECK(after.s == doctest::Approx(70 * units::mm));

    auto kps = pl.keypoints(3 * units::cm, 15 * units::deg);
    size_t nkink = 0;
    for (const auto& kp : kps) {
        if (kp.kind != DepoPolyline::KeyKind::kink) continue;
        ++nkink;
        CHECK(kp.pos.z() == doctest::Approx(50 * units::mm));
        CHECK(kp.pos.y() == doctest::Approx(0.0));
        CHECK(kp.s == doctest::Approx(50 * units::mm));
        CHECK(kp.angle == doctest::Approx(90 * units::deg));
    }
    CHECK(nkink == 1);
    CHECK(DepoPolyline::heat(Point(0, 0, 50 * units::mm), kps, DepoPolyline::KeyKind::kink) == doctest::Approx(1.0));
}

TEST_CASE("polyline k-d candidate search agrees with brute force on a long curve")
{
    DepoPolyline pl;
    const int n = 2000;
    for (int i = 0; i < n; ++i) {
        const double s = i * 0.5 * units::mm;
        pl.add(Point(0.3 * s, 50 * units::mm * std::sin(s / (40 * units::mm)), s), i, 0);
    }
    pl.build();
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(-10 * units::mm, 1010 * units::mm);
    for (int k = 0; k < 200; ++k) {
        const Point p(u(rng) * 0.3, u(rng) * 0.1, u(rng));
        auto fast = pl.nearest(p, 8);
        auto slow = pl.nearest(p, 100000);   // every segment is a candidate
        REQUIRE(fast.valid);
        CHECK(fast.dist == doctest::Approx(slow.dist).epsilon(1e-9));
        CHECK(fast.seg == slow.seg);
    }
}

TEST_CASE("polyline degenerate inputs")
{
    DepoPolyline empty;
    empty.build();
    CHECK_FALSE(empty.nearest(Point(0, 0, 0)).valid);
    CHECK(empty.keypoints().empty());

    DepoPolyline one;
    one.add(Point(1, 2, 3), 0, 0);
    one.build();
    auto h = one.nearest(Point(1, 2, 7));
    CHECK(h.valid);
    CHECK(h.dist == doctest::Approx(4));
    CHECK(one.keypoints().size() == 1);

    DepoPolyline dup;   // coincident vertices merge, their energy adds
    dup.add(Point(0, 0, 0), 0, 1.0);
    dup.add(Point(0, 0, 0), 1, 1.0);
    dup.add(Point(0, 0, 2), 2, 0.0);
    dup.build();
    CHECK(dup.nvertices() == 2);
    CHECK(dup.nearest(Point(0, 0, 1)).dedx == doctest::Approx(0.5));   // (2.0 + 0.0) / 2 per 2 mm
}
