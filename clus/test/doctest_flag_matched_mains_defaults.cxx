// doc pdvd/25 M3: pin ClusteringFlagMatchedMains' knob defaults.
//
// A NEW ensemble visitor (it flags every flash-matched cluster as a main so
// the PR taggers evaluate it); it is named only in the PDVD PR pipeline, so
// no other detector's compiled config or output changes.  The defaults are
// the PDVD semantics: live grouping, a real match time required, no length
// floor, and a cluster already flagged main is left alone.

#include "WireCellUtil/doctest.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"

using namespace WireCell;

TEST_CASE("clus knob defaults: ClusteringFlagMatchedMains")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("ClusteringFlagMatchedMains", "doc25_knobdefaults_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE(cfg.isMember("grouping"));
    CHECK(cfg["grouping"].asString() == "live");
    REQUIRE(cfg.isMember("require_t0"));
    CHECK(cfg["require_t0"].asBool() == true);
    REQUIRE(cfg.isMember("min_length"));
    CHECK(cfg["min_length"].asDouble() == 0.0);
    REQUIRE(cfg.isMember("skip_flagged"));
    CHECK(cfg["skip_flagged"].asBool() == true);
}
