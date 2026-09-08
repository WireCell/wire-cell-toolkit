// doc pdvd/48: pin CheckSTM_Michel's default_configuration().  The component
// runs only when named in a pipeline, so no legacy output depends on these;
// the test pins the documented verdict thresholds and the factory name the
// jsonnet uses.

#include "WireCellUtil/doctest.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include <string>

using namespace WireCell;

TEST_CASE("clus knob defaults: CheckSTM_Michel verdict thresholds")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("CheckSTM_Michel", "doc47_defaults_probe");
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
    B("require_stm_flag", true);
    B("publish_nu_slots", true);
    B("build_michel_shower", true);
    B("perf", false);
    B("dqdx_fit_keep_all_points", false);
    B("excl_t0_frame", false);
    D("pid_mode", 0);            // doc pdhd/03: 0 = the doc pdvd/48 criterion (gate && mu < p && mu < e)
    D("plateau_mip_lo", 0.0);    // doc pdhd/03: hi <= lo = window off
    D("plateau_mip_hi", 0.0);
    D("stop_extend_max", 0);     // doc pdhd/03: 0 = the tagger's stop is final
    B("dead_volume_check", false);
    D("min_chain_coverage", 0.0);  // doc pdhd/03: 0 = off
    D("coverage_radius_cm", 3.0);
    B("michel_guards_stop", false);   // doc pdhd/03 sec 6.3
    D("michel_shower_min_kink_deg", -1.0);   // doc pdhd/03 sec 6.8: -1 = the flag alone admits
    B("absorb_bragg_stub", false);           // doc pdhd/03 sec 6.8: measured harmful on its one test case
    B("stop_fv_use_config_tolerance", false); // doc pdhd/03 sec 6.9
    D("mip_dqdx", 50000.0);
    D("mip_dqdx_median", 43000.0);
    D("stop_snap_tol_cm", 2.0);
    D("entry_snap_tol_cm", 5.0);
    D("stop_fv_margin_cm", 5.0);
    D("bragg_contrast_min", 0.6);
    D("bragg_tail_lo_cm", 0.5);
    D("bragg_tail_hi_cm", 3.0);
    D("bragg_plateau_lo_cm", 20.0);
    D("bragg_plateau_hi_cm", 40.0);
    D("ks_margin", 0.0);
    D("compare_range_cm", 35.0);
    D("continuation_max_angle_deg", 20.0);
    D("continuation_min_len_cm", 3.0);
    D("continuation_mip_lo", 0.7);
    D("continuation_mip_hi", 1.3);
    D("michel_max_len_cm", 25.0);
    D("michel_mip_hi", 2.0);
    D("michel_mip_lo", 0.3);
    D("michel_min_kink_deg", 30.0);
    D("michel_dot_radius_cm", 15.0);
    // doc pdhd/15: DECLARED DEFAULT CHANGE, 10 -> 25 cm.  The per-piece cap now
    // matches michel_max_len_cm, the ceiling an attached Michel arm already
    // faces, so a detached Michel is judged by the same size rule.
    D("dot_max_len_cm", 25.0);
    // doc pdhd/15: companion ADMISSION, split off dot_max_len_cm (doc pdhd/13
    // defect D2).  25 = michel_max_len_cm; measured on the d14 arms to admit the
    // 11-24 cm Michel-sized neighbours and exclude every 57-292 cm cosmic.
    D("companion_max_len_cm", 25.0);
    // doc pdhd/15 sec 6: the KineChargeOptions TRACK pair for an unfitted piece
    D("michel_unfit_recom", 0.7);
    D("michel_unfit_fudge", 0.95);
    D("michel_unfit_w_ev", 23.6);
    D("dot_body_exclusion_cm", 5.0);
    D("delta_max_len_cm", 8.0);
    D("vertex_hadron_mip", 1.4);
    D("profile_min_dqdx_frac", 0.0);   // doc pdhd/03: 0 = keep every profile point
    D("fit_blob_coverage", -1.0);
    CHECK(cfg["max_candidates"].asInt() == 8);
    CHECK(cfg["min_chain_points"].asInt() == 10);
    // the PR-partition knobs are published (null = ride the C++ default)
    CHECK(cfg.isMember("two_end_break"));
    CHECK(cfg.isMember("kink_dqdx_hot_ratio"));
    CHECK(cfg["fiducial"].isNull());
}
