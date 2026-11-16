// ai-helper issue 44: the T_dlvtx_space knob of SbndPrMagnifyTrackingVisitor.
//
// dl_vtx_dump_space adds the T_dlvtx_space tree and the T_dlvtx_call
// space_cluster_id / n_space branches.  It must default off (tracking-pr.root
// stays byte-identical to every arm recorded before), and it is meaningless
// without dl_vtx_dump (the per-call rows it extends): configure() must refuse
// that, and a space_coords that does not name 3 arrays, before anything else
// is looked up.  cfg/pgrapher/experiment/sbnd/clus.jsonnet asserts the first
// at compile time too; this guards a hand-written configuration.

#include "WireCellUtil/doctest.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"

using namespace WireCell;

TEST_CASE("root knob defaults: SbndPrMagnifyTrackingVisitor dl_vtx_dump_space is OFF")
{
    PluginManager::instance().add("WireCellRoot");
    auto icfg = Factory::lookup<IConfigurable>("SbndPrMagnifyTrackingVisitor",
                                               "issue44_knobdefaults_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE_MESSAGE(cfg.isMember("dl_vtx_dump_space"), "missing knob: dl_vtx_dump_space");
    CHECK(cfg["dl_vtx_dump_space"].asBool() == false);
    REQUIRE(cfg.isMember("space_touch_cm"));
    CHECK(cfg["space_touch_cm"].asDouble() == 1.0);
    REQUIRE(cfg.isMember("space_pcname"));
    CHECK(cfg["space_pcname"].asString() == "3d");
}

TEST_CASE("SbndPrMagnifyTrackingVisitor: dl_vtx_dump_space without dl_vtx_dump is a configure error")
{
    PluginManager::instance().add("WireCellRoot");
    auto icfg = Factory::lookup<IConfigurable>("SbndPrMagnifyTrackingVisitor",
                                               "issue44_space_needs_dump_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    cfg["dl_vtx_dump_space"] = true;
    cfg["dl_vtx_dump"] = false;
    CHECK_THROWS_AS(icfg->configure(cfg), ValueError);
}

TEST_CASE("SbndPrMagnifyTrackingVisitor: dl_vtx_dump_space needs 3 space_coords")
{
    PluginManager::instance().add("WireCellRoot");
    auto icfg = Factory::lookup<IConfigurable>("SbndPrMagnifyTrackingVisitor",
                                               "issue44_space_coords_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    cfg["dl_vtx_dump_space"] = true;
    cfg["dl_vtx_dump"] = true;
    cfg["space_coords"] = Json::arrayValue;
    cfg["space_coords"].append("x_t0cor");
    cfg["space_coords"].append("y");
    CHECK_THROWS_AS(icfg->configure(cfg), ValueError);
}
