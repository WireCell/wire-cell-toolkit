// doc pdvd/111 round 2 -- Steiner::edge_support, the per-edge 2-D charge
// support count behind the log-only WCT_STEINER_GRAPH_DUMP.
//
// These cases pin:
//   * the sampling is gap_edge_bad_fraction's: round(len/step) intervals, both
//     endpoints included, at least one interval for a zero-length edge;
//   * the plane counting: live or dead makes a plane ok, ok3/ok2/ok1 are
//     cumulative, live1 needs a LIVE plane;
//   * a sample outside every TPC counts in n only;
//   * a non-positive step returns an all-zero record.

#include "WireCellClus/SteinerEdgeSupport.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;
using WireCell::Clus::Steiner::edge_support;

namespace {
    const double STEP = 0.3 * units::cm;
    const Point A(0, 0, 0);
    const Point B(0, 0, 3 * units::cm);   // 3 cm along z: 10 intervals, 11 samples
}

TEST_CASE("doc pdvd/111 r2 edge support: sampling count")
{
    auto all_live = [](const Point&, int planes[6]) {
        planes[0] = planes[1] = planes[2] = 1;
        return true;
    };
    auto s = edge_support(A, B, STEP, all_live);
    CHECK(s.n == 11);
    CHECK(s.ok3 == 11);
    CHECK(s.ok2 == 11);
    CHECK(s.ok1 == 11);
    CHECK(s.live1 == 11);

    // zero-length edge: one interval, two samples
    auto z = edge_support(A, A, STEP, all_live);
    CHECK(z.n == 2);

    // step longer than the edge: still one interval
    auto l = edge_support(A, B, 10 * units::cm, all_live);
    CHECK(l.n == 2);
}

TEST_CASE("doc pdvd/111 r2 edge support: plane counting")
{
    // first half (z < 1.5 cm): U live, V dead, W nothing  => 2 ok planes, live
    // second half: nothing at all                          => 0 ok planes
    auto cls = [](const Point& p, int planes[6]) {
        if (p.z() < 1.5 * units::cm) {
            planes[0] = 1;
            planes[4] = 1;
        }
        return true;
    };
    auto s = edge_support(A, B, STEP, cls);
    // samples at z = 0, 0.3, ..., 3.0 cm; z < 1.5 cm holds 0..1.2 = 5 samples
    CHECK(s.n == 11);
    CHECK(s.ok3 == 0);
    CHECK(s.ok2 == 5);
    CHECK(s.ok1 == 5);
    CHECK(s.live1 == 5);

    // dead in all three planes: ok3 but not live
    auto dead = [](const Point&, int planes[6]) {
        planes[3] = planes[4] = planes[5] = 1;
        return true;
    };
    auto d = edge_support(A, B, STEP, dead);
    CHECK(d.ok3 == 11);
    CHECK(d.live1 == 0);
}

TEST_CASE("doc pdvd/111 r2 edge support: outside the TPC and bad step")
{
    auto outside = [](const Point&, int planes[6]) {
        planes[0] = planes[1] = planes[2] = 1;   // ignored: the sample is outside
        return false;
    };
    auto s = edge_support(A, B, STEP, outside);
    CHECK(s.n == 11);
    CHECK(s.ok1 == 0);
    CHECK(s.live1 == 0);

    auto all_live = [](const Point&, int planes[6]) {
        planes[0] = planes[1] = planes[2] = 1;
        return true;
    };
    auto zero = edge_support(A, B, 0.0, all_live);
    CHECK(zero.n == 0);
    auto neg = edge_support(A, B, -1.0, all_live);
    CHECK(neg.n == 0);
}
