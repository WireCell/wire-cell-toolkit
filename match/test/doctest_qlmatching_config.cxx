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

TEST_CASE("qlmatching sat_flag_ignore_channels round-trips empty (docs/qlmatch/32)")
{
    Match::QLMatching qlm;
    auto cfg = qlm.default_configuration();
    REQUIRE(cfg.isMember("sat_flag_ignore_channels"));
    CHECK(cfg["sat_flag_ignore_channels"].isArray());
    CHECK(cfg["sat_flag_ignore_channels"].size() == 0);
}

#include "WireCellMatch/Opflash.h"

TEST_CASE("qlmatching sat_skip_round2_shared defaults off (docs/qlmatch/33)")
{
    // Off => fit_round2_shared keeps its legacy fill (railed rows enter the
    // round-2 joint solve); only a config that sets the key changes the fit.
    Match::QLMatching qlm;
    auto cfg = qlm.default_configuration();
    REQUIRE(cfg.isMember("sat_skip_round2_shared"));
    CHECK(cfg["sat_skip_round2_shared"].isBool());
    CHECK(cfg["sat_skip_round2_shared"].asBool() == false);
}

TEST_CASE("qlmatching lasso_weight_unrailed defaults off (docs/qlmatch/33)")
{
    Match::QLMatching qlm;
    auto cfg = qlm.default_configuration();
    REQUIRE(cfg.isMember("lasso_weight_unrailed"));
    CHECK(cfg["lasso_weight_unrailed"].isBool());
    CHECK(cfg["lasso_weight_unrailed"].asBool() == false);
}

TEST_CASE("opflash clear_sat is safe on flag-free flashes and bad channels")
{
    Match::Opflash f(0.0, std::vector<double>(40, 5.0), 1.0, 40);
    f.clear_sat({-1, 4, 11, 40, 99});
    for (int ch = -1; ch <= 40; ++ch) CHECK_FALSE(f.get_sat(ch));
}

TEST_CASE("qlmatching ks_sat_tol defaults off (docs/qlmatch/34)")
{
    Match::QLMatching qlm;
    auto cfg = qlm.default_configuration();
    REQUIRE(cfg.isMember("ks_sat_tol"));
    CHECK(cfg["ks_sat_tol"].asDouble() == 0.0);
    Match::BundleQualityParams qp;
    CHECK(qp.ks_sat_tol == 0.0);
}

TEST_CASE("ks_sat_clamp: repaired rails in a shared flash (docs/qlmatch/34)")
{
    // 4 channels; ch 0 railed.  The bundle predicts 0.3 of the flash on every
    // unrailed channel (other clusters make the rest), so s = meas/pred = 1/0.3.
    const std::vector<char> rail{1, 0, 0, 0};
    const std::vector<char> fit{1, 1, 1, 1};
    const std::vector<double> pred{300, 30, 60, 90};         // x0.3 of the flash shape
    const double s = (100.0 + 200.0 + 300.0) / (30.0 + 60.0 + 90.0);

    SUBCASE("off: tol 0 changes nothing")
    {
        std::vector<double> meas{1000, 100, 200, 300};
        CHECK(Match::ks_sat_clamp(meas, pred, rail, fit, 0.0) == 0);
        CHECK(meas[0] == 1000);
    }
    SUBCASE("a rail matching the scaled shape is untouched")
    {
        std::vector<double> meas{300 * s, 100, 200, 300};
        const double before = meas[0];
        CHECK(Match::ks_sat_clamp(meas, pred, rail, fit, 0.615) == 0);
        CHECK(meas[0] == before);
    }
    SUBCASE("a rail 3x the scaled prediction moves down by at most 1+tol")
    {
        std::vector<double> meas{3 * 300 * s, 100, 200, 300};
        CHECK(Match::ks_sat_clamp(meas, pred, rail, fit, 0.615) == 1);
        CHECK(meas[0] == doctest::Approx(3 * 300 * s / 1.615));
        CHECK(meas[1] == 100);   // unrailed untouched
    }
    SUBCASE("a rail within the tolerance becomes the scaled prediction")
    {
        std::vector<double> meas{1.3 * 300 * s, 100, 200, 300};
        CHECK(Match::ks_sat_clamp(meas, pred, rail, fit, 0.615) == 1);
        CHECK(meas[0] == doctest::Approx(300 * s));
    }
    SUBCASE("no unrailed light in the fit mask => no clamp")
    {
        std::vector<double> meas{1000, 100, 200, 300};
        const std::vector<char> fit0{1, 0, 0, 0};
        CHECK(Match::ks_sat_clamp(meas, pred, rail, fit0, 0.615) == 0);
        CHECK(meas[0] == 1000);
    }
}
