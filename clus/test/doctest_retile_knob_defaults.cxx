// doc pdhd/08 -- pin the four ImproveCluster_1 / ImproveCluster_2 knobs that
// govern the retiler's fabrication.  Every one of them is a length or a flag
// whose C++ default IS the legacy path, and every gate in this tree relies on
// that: a config that never mentions the key must reproduce the pre-knob
// behaviour.  Nothing pinned bad_blob_max_run / bad_blob_report before doc 08.
//
// Only default_configuration() is exercised: configure() goes through NeedDV
// and needs a live DetectorVolumes (same reason as
// doctest_steiner_terminal_charge_defaults.cxx).

#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Units.h"

using namespace WireCell;

TEST_CASE("pdhd doc08: the retiler's anti-ghost knobs all default to the legacy path")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("ImproveCluster_2", "doc08_retile_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();

    // doc pdvd/40 round 3
    REQUIRE_MESSAGE(cfg.isMember("bad_blob_max_run"), "missing knob: bad_blob_max_run");
    CHECK(cfg["bad_blob_max_run"].asDouble() == 0.0);       // 0 = the historical vote only
    REQUIRE_MESSAGE(cfg.isMember("bad_blob_report"), "missing knob: bad_blob_report");
    CHECK(cfg["bad_blob_report"].asBool() == false);        // log-only census, off

    // doc pdhd/08
    REQUIRE_MESSAGE(cfg.isMember("hack_max_bridge"), "missing knob: hack_max_bridge");
    CHECK(cfg["hack_max_bridge"].asDouble() == 0.0);        // 0 = uncapped = the prototype
    REQUIRE_MESSAGE(cfg.isMember("bad_blob_run_merge"), "missing knob: bad_blob_run_merge");
    CHECK(cfg["bad_blob_run_merge"].asDouble() == 0.0);     // 0 = judge each run alone
}
