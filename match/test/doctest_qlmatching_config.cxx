#include "WireCellMatch/QLMatching.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;

TEST_CASE("qlmatching default configuration knobs")
{
    // The ctor is service-free, so default_configuration() is testable
    // without Factory-registered anodes/DetectorVolumes. Guard the
    // byte-identical-when-off contracts of the per-input vector knobs:
    // absent/empty arrays must leave the scalar members in force.
    Match::QLMatching qlm;
    auto cfg = qlm.default_configuration();

    // Historical scalar defaults, relied on by detectors that do not
    // override them.
    CHECK(cfg["drift_speed"].asDouble() ==
          doctest::Approx(1.563 * units::mm / units::us));
    CHECK(cfg["trigger_offset"].asDouble() == doctest::Approx(0.0));

    // Per-input vectors round-trip as EMPTY arrays: empty => the scalar is
    // used for every input port (bit-identical legacy path).
    REQUIRE(cfg.isMember("drift_speeds"));
    CHECK(cfg["drift_speeds"].isArray());
    CHECK(cfg["drift_speeds"].size() == 0);
    REQUIRE(cfg.isMember("trigger_offsets"));
    CHECK(cfg["trigger_offsets"].isArray());
    CHECK(cfg["trigger_offsets"].size() == 0);

    // Rescue blind-spot fix (doc 23 phase 1a): knob must round-trip and
    // default OFF (bit-identical legacy ordering when absent).
    REQUIRE(cfg.isMember("postcull_before_rescue"));
    CHECK(cfg["postcull_before_rescue"].asBool() == false);
}
