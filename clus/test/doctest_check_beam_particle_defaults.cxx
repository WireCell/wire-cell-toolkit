// doc pdvd/120: pin CheckBeamParticle's default_configuration().  The
// component runs only when named in a pipeline, so no legacy output depends
// on these; the test pins the factory name the jsonnet uses, the "gate off"
// window (low >= high => the stage never runs on every bundle) and the
// nominal beam entry the doc derived.

#include "WireCellUtil/doctest.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include <string>

using namespace WireCell;

TEST_CASE("clus knob defaults: CheckBeamParticle")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("CheckBeamParticle", "doc120_defaults_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    auto D = [&](const char* k, double v) {
        const std::string msg = std::string("missing knob: ") + k;
        REQUIRE_MESSAGE(cfg.isMember(k), msg);
        CHECK(cfg[k].asDouble() == doctest::Approx(v));
    };
    auto B = [&](const char* k, bool v) {
        const std::string msg = std::string("missing knob: ") + k;
        REQUIRE_MESSAGE(cfg.isMember(k), msg);
        CHECK(cfg[k].asBool() == v);
    };
    CHECK(cfg["grouping"].asString() == "live");
    B("perf", false);
    B("publish_nu_slots", true);
    D("mip_dqdx", 50000.0);
    D("mip_dqdx_median", 43000.0);
    // the gate: low >= high means OFF, and off means "run on nothing"
    D("beam_window_low", 0.0);
    D("beam_window_high", 0.0);
    // doc pdvd/120 sec 1: the data-derived entry (cm) and the GDML beam-plug
    // axis direction
    REQUIRE(cfg["beam_entry_point_cm"].isArray());
    CHECK(cfg["beam_entry_point_cm"][0].asDouble() == doctest::Approx(110.0));
    CHECK(cfg["beam_entry_point_cm"][1].asDouble() == doctest::Approx(159.0));
    CHECK(cfg["beam_entry_point_cm"][2].asDouble() == doctest::Approx(0.6));
    REQUIRE(cfg["beam_dir"].isArray());
    CHECK(cfg["beam_dir"][0].asDouble() == doctest::Approx(-0.095));
    CHECK(cfg["beam_dir"][1].asDouble() == doctest::Approx(-0.704));
    CHECK(cfg["beam_dir"][2].asDouble() == doctest::Approx(0.704));
    D("beam_entry_max_dist_cm", 50.0);
    D("beam_dir_max_angle_deg", -1.0);   // off
    D("entry_tie_tol_cm", 5.0);
    D("entry_snap_tol_cm", 10.0);
    D("min_main_length_cm", 0.0);
    B("improve_entry_vertex", false);
    B("entry_fail_fallback_geo", false);
    B("dqdx_fit_keep_all_points", false);
    B("excl_t0_frame", false);
    D("fit_blob_coverage", -1.0);
    // the PR-partition knobs are published (null = ride the C++ default)
    CHECK(cfg.isMember("two_end_break"));
    CHECK(cfg.isMember("kink_dqdx_hot_ratio"));
    CHECK(cfg["fiducial"].isNull());
}
