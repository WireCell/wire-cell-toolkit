// doc pdvd/56 T1b: pin TaggerCheckSTM's find_first_kink asymmetric-kink knob
// defaults.
//
// find_first_kink's two charge gates both require BOTH arms >= 0.6 MIP (or a
// symmetric magnitude relaxation), so a genuine muon-into-Michel junction
// (Bragg on one side, 0.1-0.4 MIP on the other) falls through to the no-kink
// sentinel and CheckSTM_Michel clamps the stop to the fit's last row (doc 54
// sec 1.2: 039349_18/36, sum_fQ 1.77, sum_bQ 0.09).  kink_asym_enable adds a
// third, additive OR-clause to both gates.  It MUST default OFF so every
// existing SBND / ProtoDUNE configuration compiles and runs byte-identically
// -- mirrors doctest_stm_readout_edge_guard_defaults.cxx.

#include "WireCellUtil/doctest.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"

using namespace WireCell;

TEST_CASE("clus knob defaults: TaggerCheckSTM kink_asym is OFF")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("TaggerCheckSTM", "doc56_kink_asym_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE_MESSAGE(cfg.isMember("kink_asym_enable"), "missing knob: kink_asym_enable");
    CHECK(cfg["kink_asym_enable"].asBool() == false);
    REQUIRE(cfg.isMember("kink_asym_entry_mip"));
    CHECK(cfg["kink_asym_entry_mip"].asDouble() == doctest::Approx(1.2));
    REQUIRE(cfg.isMember("kink_asym_far_mip"));
    CHECK(cfg["kink_asym_far_mip"].asDouble() == doctest::Approx(0.5));
    // The sibling vertex-kink guard stays where SBND left it -- this knob is
    // deliberately its own set, not coupled to that guard's constants.
    REQUIRE(cfg.isMember("vertex_kink_guard"));
    CHECK(cfg["vertex_kink_guard"].asBool() == false);
}
