// doc pdvd/25 M3: pin CreateSteinerGraph's terminal_charge_threshold default.
//
// The Steiner terminal finder's per-point charge floor was the hard-coded WCP
// prototype constant 4000 e; it is now a knob so PDVD (7.65 / 5.10 mm pitch,
// W-plane per-point median ~1400 e, only ~12 % of points above 4000 on all
// three planes) can set its own.  The default MUST stay 4000 so every existing
// SBND / uBooNE configuration is byte-identical.

#include "WireCellUtil/doctest.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"

using namespace WireCell;

TEST_CASE("clus knob defaults: CreateSteinerGraph terminal_charge_threshold is the prototype 4000")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("CreateSteinerGraph", "doc25_terminal_charge_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE_MESSAGE(cfg.isMember("terminal_charge_threshold"), "missing knob: terminal_charge_threshold");
    CHECK(cfg["terminal_charge_threshold"].asDouble() == doctest::Approx(4000.0));
}

TEST_CASE("pdvd doc31 round6: ImproveCluster_2 terminal_charge_threshold is the prototype 4000")
{
    // doc pdvd/31 sec 8.4 found the Steiner stage running TWO terminal
    // thresholds: CreateSteinerGraph took the detector's value (PDVD 500) while
    // ImproveCluster_2 -- the retiler feeding it -- default-constructed its own
    // Grapher::Config and ran the terminal finder twice at 4000, unreachable
    // from jsonnet.  Round 6 exposes it under the SAME key name, so one jsonnet
    // value can feed both and they cannot drift apart again.
    //
    // The default must stay 4000: it is the historical value, and every
    // configuration that omits the key (which is every detector but PDVD) has
    // to stay byte-identical.  Only default_configuration() is exercised --
    // configure() goes through NeedDV and requires a live DetectorVolumes
    // component, and the default is the load-bearing protection here.
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("ImproveCluster_2", "doc31r6_terminal_charge_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE_MESSAGE(cfg.isMember("terminal_charge_threshold"), "missing knob: terminal_charge_threshold");
    CHECK(cfg["terminal_charge_threshold"].asDouble() == doctest::Approx(4000.0));
}
