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

    // Multi-input optical-PC merge (sbnd_xin/docs/99): the C++ default stays
    // OFF, so the merge keeps exactly the historical name set (only
    // root_pcs_to_merge) and nothing re-bases its flash-row indices.  ON is NOT
    // byte-identical -- the merged tree gains the non-primary inputs'
    // flash/light rows and its per-cluster "flash" scalars shift.
    //
    // A GREEN RUN HERE DOES NOT MEAN PRODUCTION IS ON THE LEGACY PATH.  Since
    // 2026-09-03 SBND production runs this knob ON: the value is set in
    // cfg/pgrapher/experiment/sbnd/{qlmatching,wct-clus-matching-perevt}.jsonnet
    // (doc 99 sec 10, ref/prod-2026-09-05).  The C++ default is what keeps the
    // OTHER binders -- pdhd and pdvd, which have their own qlmatching.jsonnet
    // and were never gated for this -- on the pre-doc-99 behaviour.  Flipping it
    // here would ship an unvalidated archive change to both detectors, which is
    // why the fix reaches production through config and not through this line.
    REQUIRE(cfg.isMember("merge_flash_pcs"));
    CHECK(cfg["merge_flash_pcs"].asBool() == false);

    // Rescue blind-spot fix (doc 23 phase 1a): knob must round-trip and
    // default OFF (bit-identical legacy ordering when absent).
    REQUIRE(cfg.isMember("postcull_before_rescue"));
    CHECK(cfg["postcull_before_rescue"].asBool() == false);

    // Saturation-aware rescue ratio-high extension (doc 23 phase 1b):
    // default OFF with inert thresholds round-tripped.
    REQUIRE(cfg.isMember("cluster_rescue_sat_ratio_relax"));
    CHECK(cfg["cluster_rescue_sat_ratio_relax"].asBool() == false);
    CHECK(cfg["cluster_rescue_sat_frac_min"].asDouble() == doctest::Approx(0.5));
    CHECK(cfg["cluster_rescue_sat_ratio_mult"].asDouble() == doctest::Approx(2.0));

    // Window-truncated overprediction cull (doc 23 phase 2): default OFF.
    REQUIRE(cfg.isMember("postcull_wtrunc_overpred"));
    CHECK(cfg["postcull_wtrunc_overpred"].asBool() == false);
    CHECK(cfg["postcull_wtrunc_ratio_hi"].asDouble() == doctest::Approx(2.0));
    CHECK(cfg["postcull_wtrunc_sat_frac"].asDouble() == doctest::Approx(0.5));

    // xtpc-pin overprediction cull (doc 23 phase 2): default OFF.
    REQUIRE(cfg.isMember("postcull_pin_overpred"));
    CHECK(cfg["postcull_pin_overpred"].asBool() == false);
    CHECK(cfg["postcull_pin_ratio_hi"].asDouble() == doctest::Approx(2.0));
}
