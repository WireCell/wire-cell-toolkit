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
