// doc pdvd/25 M3: pin TaggerCheckSTM's readout_edge_guard knob defaults.
//
// The guard vetoes an STM acceptance whose stop point's fitted arrival tick
// sits within guard_readout_edge_ticks of the readout window's edges -- a
// track cut by the window, not a Bragg stop (PDVD run 039252 evt 0 cluster
// 24 ran into time slice 0 and was accepted).  It MUST default OFF with the
// late-edge half disabled (readout_nticks 0) so every existing SBND / uBooNE
// configuration compiles and runs byte-identically.

#include "WireCellUtil/doctest.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"

using namespace WireCell;

TEST_CASE("clus knob defaults: TaggerCheckSTM readout_edge_guard is OFF")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("TaggerCheckSTM", "doc25_readout_edge_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE_MESSAGE(cfg.isMember("readout_edge_guard"), "missing knob: readout_edge_guard");
    CHECK(cfg["readout_edge_guard"].asBool() == false);
    REQUIRE(cfg.isMember("guard_readout_edge_ticks"));
    CHECK(cfg["guard_readout_edge_ticks"].asDouble() == doctest::Approx(60.0));
    REQUIRE(cfg.isMember("readout_nticks"));
    CHECK(cfg["readout_nticks"].asInt() == 0);
    // The sibling cathode guard stays where SBND left it.
    REQUIRE(cfg.isMember("cathode_guard"));
    CHECK(cfg["cathode_guard"].asBool() == false);
}
