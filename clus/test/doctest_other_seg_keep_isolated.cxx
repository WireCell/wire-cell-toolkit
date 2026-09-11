// doc sbnd_xin/docs/pr/54: unit tests for other_seg_keep_isolated_ok, the
// keep decision for isolated residual candidates in find_other_segments
// (18255-142421 separated EM shower with no fitted trajectory).  Pure
// predicate only -- the full path needs an event; its validation record is
// the pr/54 A/B arms.

#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::Clus::PR;

// The C++ defaults (NeutrinoPatternBase.h m_other_seg_keep_isolated_*).
static const int    MIN_POINTS = 25;
static const double MIN_LENGTH = 3.0 * units::cm;

TEST_CASE("pr54 keep-isolated: knob off discards unconditionally")
{
    // Legacy behaviour: no candidate is kept, however well supported.
    CHECK_FALSE(other_seg_keep_isolated_ok(false, 102, 8.0 * units::cm, MIN_POINTS, MIN_LENGTH));
    CHECK_FALSE(other_seg_keep_isolated_ok(false, 100000, 1000.0 * units::cm, MIN_POINTS, MIN_LENGTH));
    CHECK_FALSE(other_seg_keep_isolated_ok(false, 3, 0.5 * units::cm, MIN_POINTS, MIN_LENGTH));
}

TEST_CASE("pr54 keep-isolated: knob on separates well-supported from sparse")
{
    // The 18255-142421 candidates: the 102-point separated-shower piece must
    // be kept; the 3- and 4-point noise fragments must still be discarded.
    CHECK(other_seg_keep_isolated_ok(true, 102, 8.0 * units::cm, MIN_POINTS, MIN_LENGTH));
    CHECK_FALSE(other_seg_keep_isolated_ok(true, 4, 8.0 * units::cm, MIN_POINTS, MIN_LENGTH));
    CHECK_FALSE(other_seg_keep_isolated_ok(true, 3, 8.0 * units::cm, MIN_POINTS, MIN_LENGTH));

    // Both floors are AND-ed: many points but a stub-length fit still drops.
    CHECK_FALSE(other_seg_keep_isolated_ok(true, 102, 2.0 * units::cm, MIN_POINTS, MIN_LENGTH));

    // Boundary: floors are inclusive.
    CHECK(other_seg_keep_isolated_ok(true, MIN_POINTS, MIN_LENGTH, MIN_POINTS, MIN_LENGTH));
    CHECK_FALSE(other_seg_keep_isolated_ok(true, MIN_POINTS - 1, MIN_LENGTH, MIN_POINTS, MIN_LENGTH));
}

// doc sbnd_xin/docs/pr/102 P1 -- the length OR-disjunct.  doc 77 round 4
// retired P1's other disjunct (the nnf one) and took its TEST_CASE with it:
// the knob's validation FAILED at 4 and it never left 0 in any config.
TEST_CASE("pr102 keep-isolated: length disjunct admits long candidates at any terminal count")
{
    const double LEN_ADMIT = 30.0 * units::cm;
    // pr/102 sec 4 B1: 145.5 cm at 23 terminals, 67.1 cm at 16 -- admitted.
    CHECK(other_seg_keep_isolated_ok(true, 23, 145.5 * units::cm, MIN_POINTS, MIN_LENGTH, LEN_ADMIT));
    CHECK(other_seg_keep_isolated_ok(true, 16, 67.1 * units::cm, MIN_POINTS, MIN_LENGTH, LEN_ADMIT));
    // The pr/54 noise population is <= 10 cm: stays out at the same setting.
    CHECK_FALSE(other_seg_keep_isolated_ok(true, 10, 10.0 * units::cm, MIN_POINTS, MIN_LENGTH, LEN_ADMIT));
    // Disjunct inert at its 0 default.
    CHECK_FALSE(other_seg_keep_isolated_ok(true, 23, 145.5 * units::cm, MIN_POINTS, MIN_LENGTH, 0.0));
    // Boundary inclusive.
    CHECK(other_seg_keep_isolated_ok(true, 1, LEN_ADMIT, MIN_POINTS, MIN_LENGTH, LEN_ADMIT));
    // Knob-off master switch still discards everything.
    CHECK_FALSE(other_seg_keep_isolated_ok(false, 23, 145.5 * units::cm, MIN_POINTS, MIN_LENGTH, LEN_ADMIT));
}

TEST_CASE("pr102 keep-isolated: pre-pr/102 call shape unchanged (default args)")
{
    // The legacy 5-argument call must behave exactly as before.
    CHECK(other_seg_keep_isolated_ok(true, 102, 8.0 * units::cm, MIN_POINTS, MIN_LENGTH));
    CHECK_FALSE(other_seg_keep_isolated_ok(true, 4, 8.0 * units::cm, MIN_POINTS, MIN_LENGTH));
}

// doc pdvd/87: other_seg_keep_anchor_ok, doc pdvd/62's stop-local keep with
// the size floor.  The residuals are the ones doc pdvd/86 sec 8.1 read off
// PDVD production's pr54 drop lines.
TEST_CASE("pdvd87 keep near-anchor: no anchor or radius 0 never keeps")
{
    CHECK_FALSE(other_seg_keep_anchor_ok(-1.0, 5.0 * units::cm, 20, 10.6 * units::cm, 0, 0.0));
    CHECK_FALSE(other_seg_keep_anchor_ok(-1.0, 5.0 * units::cm, 20, 10.6 * units::cm, 5, 5.0 * units::cm));
    CHECK_FALSE(other_seg_keep_anchor_ok(3.1 * units::cm, 0.0, 20, 10.6 * units::cm, 0, 0.0));
}

TEST_CASE("pdvd87 keep near-anchor: 0/0 floors are doc 62's radius test")
{
    for (double d : {0.0, 1.1, 3.1, 4.99, 5.0, 5.01, 19.0}) {
        const bool radius = d <= 5.0;
        CHECK(other_seg_keep_anchor_ok(d * units::cm, 5.0 * units::cm, 2, 1.19 * units::cm, 0, 0.0) == radius);
    }
}

TEST_CASE("pdvd87 keep near-anchor: the 5/5 floor, inclusive and AND-ed")
{
    const double R = 5.0 * units::cm, L = 5.0 * units::cm;
    const int N = 5;
    CHECK(other_seg_keep_anchor_ok(3.1 * units::cm, R, 20, 10.61 * units::cm, N, L));    // 039253_0/44
    CHECK(other_seg_keep_anchor_ok(1.4 * units::cm, R, 6, 8.53 * units::cm, N, L));      // 039349_30/45, the 6-terminal piece
    CHECK_FALSE(other_seg_keep_anchor_ok(1.7 * units::cm, R, 4, 5.47 * units::cm, N, L)); // 039349_30/45, the 4-terminal piece
    CHECK_FALSE(other_seg_keep_anchor_ok(1.1 * units::cm, R, 2, 1.19 * units::cm, N, L)); // 039349_61/21, doc 62's lost TP
    CHECK_FALSE(other_seg_keep_anchor_ok(2.7 * units::cm, R, 3, 6.77 * units::cm, N, L)); // 039349_20/73 (THRU): long, too few terminals
    CHECK_FALSE(other_seg_keep_anchor_ok(1.9 * units::cm, R, 6, 0.82 * units::cm, N, L)); // 039252_17/91 (THRU): enough terminals, a stub
    CHECK(other_seg_keep_anchor_ok(R, R, N, L, N, L));                                    // boundaries inclusive
    CHECK_FALSE(other_seg_keep_anchor_ok(R, R, N - 1, L, N, L));
    CHECK_FALSE(other_seg_keep_anchor_ok(R, R, N, 0.999 * L, N, L));
    CHECK_FALSE(other_seg_keep_anchor_ok(1.01 * R, R, 20, 10.61 * units::cm, N, L));    // outside the radius, whatever the size
}
