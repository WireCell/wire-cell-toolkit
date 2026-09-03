// doc pdvd/25 sec 13.11: pin ClusteringProtectBundle's stm_only_bundles default.
//
// The knob is the mirror of TaggerCheckNeutrino's nu_per_bundle_stm_only one
// stage earlier: it opens only the bundles that hold an STM-tagged cluster, so
// the graph examination is not spent on bundles the per-bundle neutrino PR will
// never read.  It is NOT a no-op on the archives -- a cosmic bundle that is no
// longer split keeps its over-clustered shape -- so it MUST default OFF, and
// the two knobs it rides on (beam_window_only, open_convicted_bundles) must
// keep their own defaults beside it.

#include "WireCellUtil/doctest.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"

using namespace WireCell;

TEST_CASE("clus knob defaults: ClusteringProtectBundle stm_only_bundles is OFF")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("ClusteringProtectBundle", "doc25_pb_stm_only_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE_MESSAGE(cfg.isMember("stm_only_bundles"), "missing knob: stm_only_bundles");
    CHECK(cfg["stm_only_bundles"].asBool() == false);
    // The gate it needs to be meaningful is itself OFF by default, so the knob
    // is inert twice over in the shipped configuration.
    REQUIRE(cfg.isMember("beam_window_only"));
    CHECK(cfg["beam_window_only"].asBool() == false);
    // The pr/94 companion is a different rule and keeps its own default: an
    // STM main is convicted, so on PDVD it is open_convicted_bundles that
    // actually opens the STM bundle this knob selects.
    REQUIRE(cfg.isMember("open_convicted_bundles"));
    CHECK(cfg["open_convicted_bundles"].asBool() == false);
    REQUIRE(cfg.isMember("skip_convicted"));
    CHECK(cfg["skip_convicted"].asBool() == true);
}
