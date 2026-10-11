// doc pdvd/140 sec 5: the cap on the chord find_other_segments makes when it relocates the attachment point of a
// piece shorter than 3 cm onto the tracked body.  Pure predicate only (the full path needs an event; its record
// is the doc 140 arms): off = always relocate, on = relocate only up to the cap, and the knob's default is off.
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::Clus::PR;

TEST_CASE("pdvd140 relocate cap: off relocates at any chord length")
{
    CHECK(other_seg_relocate_ok(0, 0.5 * units::cm));
    CHECK(other_seg_relocate_ok(0, 100 * units::cm));
    CHECK(other_seg_relocate_ok(-1, 100 * units::cm));
}

TEST_CASE("pdvd140 relocate cap: on refuses a chord longer than the cap")
{
    const double cap = 10 * units::cm;
    CHECK(other_seg_relocate_ok(cap, 2 * units::cm));
    CHECK(other_seg_relocate_ok(cap, 10 * units::cm));   // not longer than the cap
    CHECK_FALSE(other_seg_relocate_ok(cap, 10.1 * units::cm));
    CHECK_FALSE(other_seg_relocate_ok(cap, 80 * units::cm));
}

TEST_CASE("pdvd140 relocate cap: the PatternAlgorithms default is off")
{
    WireCell::Clus::PR::PatternAlgorithms pa;
    CHECK(pa.m_other_seg_relocate_max == 0.0);
}
