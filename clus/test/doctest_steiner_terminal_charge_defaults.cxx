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
