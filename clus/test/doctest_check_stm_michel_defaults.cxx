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
    D("dot_max_len_cm", 10.0);
    D("dot_body_exclusion_cm", 5.0);
    D("delta_max_len_cm", 8.0);
    D("vertex_hadron_mip", 1.4);
    D("fit_blob_coverage", -1.0);
    CHECK(cfg["max_candidates"].asInt() == 8);
    CHECK(cfg["min_chain_points"].asInt() == 10);
    // the PR-partition knobs are published (null = ride the C++ default)
    CHECK(cfg.isMember("two_end_break"));
    CHECK(cfg.isMember("kink_dqdx_hot_ratio"));
    CHECK(cfg["fiducial"].isNull());
}
