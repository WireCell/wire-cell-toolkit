// doc 87: pin the SbndPrMagnifyTrackingVisitor knob defaults.
//
// Companion to clus/test/doctest_clus_knob_defaults.cxx, and it carries the
// same caveat: a green run here proves the C++ DEFAULT, not what SBND
// production runs.  The operating point itself is guarded by
// sbnd_xin/scripts/cfg/prod_cfg_gate.py against ref/prod-<date>/.
//
// save_in_scope adds the T_cluster tree (the in-scope set plus the per-bundle
// summary).  It MUST default false: with it on, tracking-pr.root gains a tree
// and is no longer byte-identical to every arm recorded before doc 87.

#include "WireCellUtil/doctest.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"

using namespace WireCell;

TEST_CASE("root knob defaults: SbndPrMagnifyTrackingVisitor save_in_scope is OFF")
{
    PluginManager::instance().add("WireCellRoot");
    // A private instance name is virgin by construction.
    auto icfg = Factory::lookup<IConfigurable>("SbndPrMagnifyTrackingVisitor",
                                               "doc87_knobdefaults_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();

    REQUIRE_MESSAGE(cfg.isMember("save_in_scope"), "missing knob: save_in_scope");
    CHECK(cfg["save_in_scope"].asBool() == false);

    // doc 99.  flash_by_gid switches T_cluster's flash columns from the per-input
    // "flash" row index (Cluster::get_flash()) to matched_flash_gid resolved
    // against the merge-safe "opflash" PC.  With it on the three flash columns
    // change value on ~half the rows, and flash_id changes MEANING (it carries
    // the gid), so tracking-pr.root stops being comparable to every arm recorded
    // before doc 99.
    //
    // A GREEN RUN HERE DOES NOT MEAN PRODUCTION IS ON THE LEGACY PATH.  Since
    // 2026-09-03 SBND production runs this knob ON, set in
    // cfg/pgrapher/experiment/sbnd/{clus,wct-pr-perevt}.jsonnet (doc 99 sec 10,
    // ref/prod-2026-09-05).  What the C++ default still buys is the OTHER fork:
    // PdvdPrMagnifyTrackingVisitor shares this default and its gid-uniqueness
    // precondition has never been checked, so leaving it false here is what keeps
    // PDVD off until someone does (doctest_pdvd_tracking_defaults.cxx asserts it).
    REQUIRE_MESSAGE(cfg.isMember("flash_by_gid"), "missing knob: flash_by_gid");
    CHECK(cfg["flash_by_gid"].asBool() == false);

    // sbnd_xin/docs/109.  nu_provenance adds T_bundle/T_flash and Trun census
    // branches, fix_cluster_flags changes the VALUES of T_cluster's
    // tgm/stm/fc/lm/beam_flash columns, and a non-empty provenance map adds
    // Trun string branches.  All must default off, so every arm recorded
    // before doc 109 stays comparable; SBND production turns them on in
    // cfg/pgrapher/experiment/sbnd/wct-pr-perevt.jsonnet.
    REQUIRE_MESSAGE(cfg.isMember("nu_provenance"), "missing knob: nu_provenance");
    CHECK(cfg["nu_provenance"].asBool() == false);
    REQUIRE_MESSAGE(cfg.isMember("fix_cluster_flags"), "missing knob: fix_cluster_flags");
    CHECK(cfg["fix_cluster_flags"].asBool() == false);
    REQUIRE_MESSAGE(cfg.isMember("provenance"), "missing knob: provenance");
    CHECK(cfg["provenance"].isObject());
    CHECK(cfg["provenance"].size() == 0);

    // The legacy output name must not drift either -- the runner and every
    // gate script look for exactly this file.
    REQUIRE(cfg.isMember("output_filename"));
    CHECK(cfg["output_filename"].asString() == "tracking-pr.root");
}

// sbnd_xin/docs/109.  UbooneTaggerOutputVisitor is shared with the PDHD and
// PDVD PR configs, so its provenance branches must stay unbooked by default.
TEST_CASE("root knob defaults: UbooneTaggerOutputVisitor nu_provenance is OFF")
{
    PluginManager::instance().add("WireCellRoot");
    auto icfg = Factory::lookup<IConfigurable>("UbooneTaggerOutputVisitor",
                                               "doc109_knobdefaults_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE_MESSAGE(cfg.isMember("nu_provenance"), "missing knob: nu_provenance");
    CHECK(cfg["nu_provenance"].asBool() == false);
}
