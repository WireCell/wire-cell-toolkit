// doc pdvd/25 sec 13.10: pin TaggerCheckNeutrino's nu_per_bundle_stm_only default.
//
// The knob restricts the per-bundle neutrino PR to bundles whose selected
// activity carries Flags::STM -- PDVD's working mode, where the stopping muon
// IS the object being reconstructed (doc 25 sec 2.3).  It is the exact
// inverse of nu_skip_cosmic, so it MUST default OFF, and nu_skip_cosmic must
// stay OFF beside it, or every existing SBND / uBooNE configuration would
// change which activity it calls the neutrino.

#include "WireCellUtil/doctest.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"

using namespace WireCell;

TEST_CASE("clus knob defaults: TaggerCheckNeutrino nu_per_bundle_stm_only is OFF")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("TaggerCheckNeutrino", "doc25_stm_only_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE_MESSAGE(cfg.isMember("nu_per_bundle_stm_only"), "missing knob: nu_per_bundle_stm_only");
    CHECK(cfg["nu_per_bundle_stm_only"].asBool() == false);
    // Its inverse, and the per-bundle mode it rides on, stay at their defaults:
    // both OFF means the knob is inert twice over in the shipped configuration.
    REQUIRE(cfg.isMember("nu_skip_cosmic"));
    CHECK(cfg["nu_skip_cosmic"].asBool() == false);
    REQUIRE(cfg.isMember("nu_per_bundle"));
    CHECK(cfg["nu_per_bundle"].asBool() == false);
    // The length floor is a different rule and keeps its own default (0 = none);
    // the STM filter deliberately does NOT inherit its legacy-winner exemption.
    REQUIRE(cfg.isMember("nu_per_bundle_min_length"));
    CHECK(cfg["nu_per_bundle_min_length"].asDouble() == doctest::Approx(0.0));
}
