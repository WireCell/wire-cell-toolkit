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
    B("entry_long_muon_absorb", false);   // doc pdvd/120 sec 9.3: release the entry-rooted chain
    B("kine_charge_t0_frame", true);      // doc pdvd/120 sec 9.4: this stage's own default is ON
    B("dqdx_fit_keep_all_points", false);
    B("excl_t0_frame", false);
    D("fit_blob_coverage", -1.0);
    // doc pdvd/120 sec 9: the stage reads TaggerCheckNeutrino's PR knob set
    // (same names, same C++ defaults, TaggerCheckNeutrino.h) -- a sample from
    // each family the first cut missed ...
    B("two_end_break", false);
    D("kink_dqdx_hot_ratio", 1.7);
    B("main_vertex_graph_audit", false);
    D("mvga_radius", 15.0);
    B("shower_nv_bridge_track", false);
    D("shower_nv_bridge_max_gap", 1.8);
    B("straight_cont_cross_cluster", false);
    B("kine_count_conn4_near", false);
    D("kine_conn4_near_gap", 20.0);
    B("shower_split", false);
    D("pi0_mass_offset", 10.0);
    D("conn3_stitch_max", 0.0);
    B("swap_orphan_dup_audit", false);
    B("long_muon_range_empty_chain_fallback", false);
    B("fit_blob_coverage_defer", false);
    D("kine_fudge_factor", 0.95);
    D("vertex_z_prior_scale", 200.0);
    REQUIRE(cfg["kine_plane_weights"].isArray());
    CHECK(cfg["kine_plane_weights"][2].asDouble() == doctest::Approx(1.0));
    // ... and the families it deliberately does not read
    CHECK(!cfg.isMember("dl_weights"));
    CHECK(!cfg.isMember("nu_per_bundle"));
    CHECK(!cfg.isMember("nu_skip_cosmic"));
    CHECK(!cfg.isMember("vertex_kink_snap"));
    CHECK(!cfg.isMember("vertex_junction_snap"));
    CHECK(!cfg.isMember("mcs_enable"));
    CHECK(!cfg.isMember("long_muon_cathode_bridge"));
    CHECK(!cfg.isMember("cosmic_y_top_main"));
    CHECK(!cfg.isMember("muon_dqdx_curve"));
    CHECK(cfg["fiducial"].isNull());
}
