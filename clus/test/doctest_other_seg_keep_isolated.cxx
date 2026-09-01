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
