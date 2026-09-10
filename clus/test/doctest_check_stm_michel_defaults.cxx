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
    D("stop_retreat_max", 0);           // doc pdvd/57: 0 = off, the tagger's stop is final
    D("retreat_collapse_frac", 0.5);
    D("retreat_peak_frac", 1.4);
    D("retreat_peak_window_cm", 15.0);
    D("stop_split_max", 0);             // doc pdvd/58 (T1c): 0 = off, no fit-row split
    D("split_kink_min_deg", 15.0);
    D("split_min_drop_cm", 3.0);
    D("split_collapse_frac", 0.5);
    D("split_peak_frac", 1.4);
    D("split_peak_window_cm", 15.0);
    D("split_dir_window_cm", 5.0);
    B("moved_stop_michel_guard", false);   // doc pdvd/61 (T2c): 0/false = off, no attached-arm veto
    D("moved_stop_michel_ke_min", 10.0);
    D("stop_local_residual_cm", 0.0);      // doc pdvd/62 (T3a): 0 = off, the pr54 floors alone decide
    B("stop_local_michel_pieces", false);  // doc pdvd/62 (T3b)
    B("michel_range_energy_guard", false); // doc pdvd/62 (T3c)
    D("michel_range_energy_dis_cm", 5.0);  // the argued distance, not doc 55's best-F1 3 cm row
    D("michel_range_energy_ke_min", 10.0);
    B("publish_other_arms", false);        // doc pdvd/64 (T6): role-7 rows for kOther arms, rows only
    B("bragg_peak_anchor", false);         // doc pdvd/65 (T7): geometric rr origin unless on
    D("bragg_peak_search_cm", 10.0);
    B("profile_geometry_guard", false);    // doc pdvd/66 (T8): fields always written, the reject bit only when on
    D("profile_arc_span_max", 1.5);
    D("unsupported_min_len_cm", 20.0);
    D("unsupported_frac", 0.25);
    D("end_window_cm", 20.0);
    B("topology_stop_evidence", false);     // doc pdvd/70 (P1): the dQ/dx shape tests alone decide unless on
    D("topology_michel_ke_min", 10.0);      // the T2c / T3c floor
    D("topology_michel_len_min_cm", 3.0);   // the range-energy distance
    B("topology_clears_sparse", false);     // doc pdvd/70 sec 9.4: the owner's option
    B("michel_gamma_collect", false);       // doc pdvd/71 (P4): the Michel keeps its core only unless on
    D("michel_gamma_radius_cm", 35.0);      // = today's admission radius, so on never widens it
    D("michel_gamma_max_len_cm", 10.0);     // a dot, the capture-gamma stage's compactness cap
    D("michel_gamma_cos_min", 0.5);         // a 60 deg cone about the Michel direction
    D("michel_gamma_max_ke_mev", 20.0);     // per blob, the capture-gamma stage's cap
    D("michel_gamma_total_ke_max_mev", 60.0);   // the 52.8 MeV endpoint plus resolution
    D("moved_stop_michel_kink_min", -1.0);  // doc pdvd/72 (P3b): -1 = off, the KE floor alone decides
    D("michel_mip_lo_turned", -1.0);        // doc pdvd/73 (P2a): -1 = off, michel_mip_lo for every arm
    D("michel_mip_lo_turned_kink_deg", 60.0);
    D("michel_far_len_shower_max_cm", -1.0);   // doc pdvd/73 (P2b): -1 = off, len + far_len for every arm
    D("michel_kink_window_cm", -1.0);       // doc pdvd/73 (P2c): -1 = off, the classifier's window alone
    B("retreat_tail_strict", false);        // doc pdvd/74 (P3): the doc 57 tail reading unless on
    B("retreat_tail_sublive", false);
    B("michel_collinear_split", false);     // doc pdvd/74 (P3): a Bragg-confirmed chain is never retreated / split
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
    // doc pdvd/51: the muon-capture gamma at the stop.  A SEPARATE object class
    // from the Michel -- its ring starts where michel_dot_radius_cm ends, so
    // nothing the Michel gathers can move when these change.  stop_gamma_enable
    // false reproduces the doc pdhd/17 tree on all 97 scalar branches.
    B("stop_gamma_enable", true);
    // 35 cm: the same-bundle stop-anchored density knee, and the plateau of the
    // mu- (no-Michel) anti-correlation, which collapses by 45 cm.
    D("stop_gamma_radius_cm", 35.0);
    D("stop_gamma_max_len_cm", 10.0);   // a gamma deposit is a blob, not a track
    D("stop_gamma_min_ke_mev", 0.2);
    D("stop_gamma_max_ke_mev", 20.0);
    D("stop_gamma_max_n", 8);

    // doc pdvd/53: the SURVEY.  Default OFF is the load-bearing value -- with it
    // false an absent key leaves the whole component byte-identical, including
    // the stm_michel_pts column set (no rej/d_stop/d_body) and the
    // T_stm_michel branch list (no n_survey_*).  The two ProtoDUNE pr.jsonnet
    // turn it on; nothing else binds this component.
    B("survey_enable", false);
    // 60 cm: the survey exists to put the pieces a scanner can SEE in front of
    // them, and the display's own near-image radius is the same 60 cm.  It is
    // not a physics boundary -- the Michel keeps michel_dot_radius_cm and the
    // gamma keeps stop_gamma_radius_cm, so widening this admits nothing to
    // either object.
    D("survey_radius_cm", 60.0);
    // Mirrors companion_max_len_cm, so the knob-on companion POOL is the same
    // pool the knob-off path builds, only reaching further.
    D("survey_max_len_cm", 25.0);
    // doc pdhd/15 sec 6: the KineChargeOptions TRACK pair for an unfitted piece
    D("michel_unfit_recom", 0.7);
    D("michel_unfit_fudge", 0.95);
    D("michel_unfit_w_ev", 23.6);
    // doc pdhd/17: derive the unfitted-charge survival from the bound
    // recombination model instead of the pair above.  C++ default OFF, so an
    // absent key leaves dots_ke_unfit exactly where doc pdhd/15 put it; both
    // ProtoDUNE drivers set it true.  2.1 MeV/cm = the MIP assumption, and the
    // same pivot PowerBoxRecombination uses.
    B("michel_unfit_from_model", false);
    D("michel_unfit_dedx", 2.1);
    D("dot_body_exclusion_cm", 5.0);
    D("delta_max_len_cm", 8.0);
    D("vertex_hadron_mip", 1.4);
    D("profile_min_dqdx_frac", 0.0);   // doc pdhd/03: 0 = keep every profile point
    D("fit_blob_coverage", -1.0);
    // doc pdhd/16: MCS is default OFF in C++ so an absent bag leaves the
    // compiled config and the tree values exactly where doc pdhd/15 left them;
    // both ProtoDUNE drivers set mcs_enable true.  cathode_xcut 0 = the
    // excision is off, which is bit-for-bit upstream (MuonMCS.h:79).
    B("mcs_enable", false);
    D("mcs_min_len_cm", 40.0);
    D("mcs_cathode_x", 0.0);
    D("mcs_cathode_xcut", 0.0);
    CHECK(cfg["max_candidates"].asInt() == 8);
    CHECK(cfg["min_chain_points"].asInt() == 10);
    // the PR-partition knobs are published (null = ride the C++ default)
    CHECK(cfg.isMember("two_end_break"));
    CHECK(cfg.isMember("kink_dqdx_hot_ratio"));
    CHECK(cfg["fiducial"].isNull());
}
