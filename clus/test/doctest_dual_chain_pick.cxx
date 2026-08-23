// doc sbnd_xin/docs/pr/112 sec 11: unit tests for dual_chain_pick, the pure
// snap-and-guard arithmetic of the dual chain's "snap" mode.  The full OFF
// pass needs an event; its validation record is the pr/112 sec 11 arms.

#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <vector>

using namespace WireCell;
using namespace WireCell::Clus::PR;

namespace {
    const double cm = units::cm;
    // three production candidates: two on cluster 7 (main), one on cluster 9
    std::vector<Point> cands() { return {Point(0, 0, 0), Point(5 * cm, 0, 0), Point(0, 3 * cm, 0)}; }
    std::vector<int> cids() { return {7, 7, 9}; }
}

TEST_CASE("pr112 dual_chain_pick: nearest candidate within D transfers")
{
    auto p = dual_chain_pick(cands(), cids(), 7, Point(4.5 * cm, 0, 0), 2.0 * cm, true);
    CHECK(p.index == 1);
    CHECK(p.dis == doctest::Approx(0.5 * cm));
    CHECK(p.accepted);
}

TEST_CASE("pr112 dual_chain_pick: 0.00 cm IS a transfer (sec 5.7.4), even at D=0")
{
    auto p = dual_chain_pick(cands(), cids(), 7, Point(5 * cm, 0, 0), 0.0, true);
    CHECK(p.index == 1);
    CHECK(p.dis == doctest::Approx(0.0));
    CHECK(p.accepted);
}

TEST_CASE("pr112 dual_chain_pick: beyond D keeps production (not accepted, index still reported)")
{
    auto p = dual_chain_pick(cands(), cids(), 7, Point(8 * cm, 0, 0), 2.0 * cm, true);
    CHECK(p.index == 1);
    CHECK(p.dis == doctest::Approx(3.0 * cm));
    CHECK_FALSE(p.accepted);
}

TEST_CASE("pr112 dual_chain_pick: cluster-swap gate excludes other-cluster candidates")
{
    const Point off(0, 3.2 * cm, 0);   // 0.2 cm from the cluster-9 candidate, 3.2 from (0,0,0)
    auto allow = dual_chain_pick(cands(), cids(), 7, off, 2.0 * cm, true);
    CHECK(allow.index == 2);
    CHECK(allow.accepted);
    auto deny = dual_chain_pick(cands(), cids(), 7, off, 2.0 * cm, false);
    CHECK(deny.index == 0);
    CHECK(deny.dis == doctest::Approx(3.2 * cm));
    CHECK_FALSE(deny.accepted);   // 3.2 > 2.0: production keeps its own pick
}

TEST_CASE("pr112 dual_chain_pick: empty candidate set never accepts")
{
    auto p = dual_chain_pick({}, {}, 7, Point(0, 0, 0), 2.0 * cm, true);
    CHECK(p.index == -1);
    CHECK_FALSE(p.accepted);
}
