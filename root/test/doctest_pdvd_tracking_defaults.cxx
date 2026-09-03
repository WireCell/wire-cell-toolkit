// doc pdvd/25 M6: pin the PDVD Magnify-tracking writers' knob defaults.
//
// PdvdMagnifyTrackingVisitor / PdvdPrMagnifyTrackingVisitor are forks by
// duplication of the SBND writers with a face-aware channel scheme (PDVD's
// two-sided CRPs carry different channel sets per face).  A green run here
// proves the C++ DEFAULT, not what the PDVD chain runs: the operating point
// is pdvd/wct-pr-perevt.jsonnet.  The readout length default is PDVD's
// production window (10000 ticks), not SBND's 3427; the output names must not
// drift -- run_pr_evt.sh and the STM/Michel scripts look for exactly them.

#include "WireCellUtil/doctest.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"

using namespace WireCell;

TEST_CASE("root knob defaults: PdvdMagnifyTrackingVisitor (STM stage)")
{
    PluginManager::instance().add("WireCellRoot");
    auto icfg = Factory::lookup<IConfigurable>("PdvdMagnifyTrackingVisitor",
                                               "doc25_knobdefaults_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE(cfg.isMember("output_filename"));
    CHECK(cfg["output_filename"].asString() == "tracking-stm.root");
    REQUIRE(cfg.isMember("track_fitting_name"));
    CHECK(cfg["track_fitting_name"].asString() == "stm");
    REQUIRE(cfg.isMember("nticks"));
    CHECK(cfg["nticks"].asInt() == 10000);
}

TEST_CASE("root knob defaults: PdvdPrMagnifyTrackingVisitor (PR stage)")
{
    PluginManager::instance().add("WireCellRoot");
    auto icfg = Factory::lookup<IConfigurable>("PdvdPrMagnifyTrackingVisitor",
                                               "doc25_knobdefaults_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE(cfg.isMember("output_filename"));
    CHECK(cfg["output_filename"].asString() == "tracking-pr.root");
    REQUIRE_MESSAGE(cfg.isMember("save_in_scope"), "missing knob: save_in_scope");
    CHECK(cfg["save_in_scope"].asBool() == false);
    // sbnd_xin/docs/99 sec 9.  flash_by_gid resolves T_cluster's flash columns
    // through matched_flash_gid against the "opflash" PC.  It is validated on
    // SBND ONLY: its precondition is that a gid names one flash, which holds
    // because SBND's gid side is the anode ident.  PDVD's gid encoding
    // (opflash_phys_gid / shared_flash, per-drift-side flash lists) has NOT been
    // checked against it, so this default is the tripwire that keeps PDVD off
    // until someone does.  Do not flip it to true without that check.
    REQUIRE_MESSAGE(cfg.isMember("flash_by_gid"), "missing knob: flash_by_gid");
    CHECK(cfg["flash_by_gid"].asBool() == false);
    REQUIRE(cfg.isMember("nticks"));
    CHECK(cfg["nticks"].asInt() == 10000);
}
