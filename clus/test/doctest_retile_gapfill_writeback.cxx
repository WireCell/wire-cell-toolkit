// pr_ml doc 28 step 2 -- ImproveCluster_2 retile_deghost "writeback".
//
// The doc pdvd/131 gapfill retile tags the cluster so CreateSteinerGraph replaces its blobs by the base + corridor
// cells.  On simulation that drops the blobs' truth cloud (the R labels then score every segment at purity 0) and,
// on SBND, mixes cluster-scalar and "3d" key sets (TensorDM refuses the pctree).  "writeback": false keeps the
// cluster's blobs and lets the gapfill cloud feed the Steiner graph only.  These cases pin:
//   1. the key absent => true, the doc 130 behaviour of every existing gapfill config (PDVD -beam);
//   2. false is parsed (echoed by default_configuration once retile_deghost is configured);
//   3. without retile_deghost the echo key is absent, so a production config never sees it.
// configure() may stop in the base class (NeedDV wants a live DetectorVolumes; whether one exists depends on the
// cases run before) AFTER the retile_deghost object is parsed, so each case tolerates that throw and reads the echo.

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITensorForward.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;

namespace {
    struct WritebackFakeForward : public ITensorForward {
        virtual ~WritebackFakeForward() {}
        virtual ITensorSet::pointer forward(const ITensorSet::pointer& in) const { return in; }
    };
}  // namespace

WIRECELL_FACTORY(WritebackFakeForward, WritebackFakeForward, WireCell::ITensorForward)

namespace {
    Configuration gapfill_cfg()
    {
        Configuration rd;
        Configuration lv;
        lv["width"] = 0;
        lv["superwire"] = 16;
        lv["threshold"] = 0.0;
        lv["forward"] = "WritebackFakeForward:wb";
        rd["levels"].append(lv);
        rd["scope"] = "all";
        Configuration cfg;
        cfg["retile_deghost"] = rd;
        return cfg;
    }
    void configure_tolerant(IConfigurable::pointer icfg, const Configuration& cfg)
    {
        try {
            icfg->configure(cfg);
        }
        catch (const WireCell::Exception&) {   // the base class's NeedDV; the retile_deghost object is parsed first
        }
    }
    IConfigurable::pointer probe(const std::string& name)
    {
        PluginManager::instance().add("WireCellClus");
        make_WritebackFakeForward_factory();
        Factory::lookup_tn<ITensorForward>("WritebackFakeForward:wb");   // create the instance configure() finds
        return Factory::lookup<IConfigurable>("ImproveCluster_2", name);
    }
}  // namespace

TEST_CASE("pr_ml doc 28: gapfill writeback defaults to true (doc 130 behaviour)")
{
    auto icfg = probe("d28_wb_default");
    CHECK(!icfg->default_configuration().isMember("retile_deghost_writeback"));   // case 3
    configure_tolerant(icfg, gapfill_cfg());
    auto cfg = icfg->default_configuration();
    REQUIRE(cfg.isMember("retile_deghost_writeback"));
    CHECK(cfg["retile_deghost_writeback"].asBool() == true);                       // case 1
}

TEST_CASE("pr_ml doc 28: gapfill writeback false is parsed")
{
    auto icfg = probe("d28_wb_off");
    auto c = gapfill_cfg();
    c["retile_deghost"]["writeback"] = false;
    configure_tolerant(icfg, c);
    auto cfg = icfg->default_configuration();
    REQUIRE(cfg.isMember("retile_deghost_writeback"));
    CHECK(cfg["retile_deghost_writeback"].asBool() == false);                      // case 2
}

// pr_ml doc 32 sec C.3 / C.4 -- skip_base_subcells: absent => true (flipped 2026-10-28), false parsed (the doc 131
// behaviour), and no echo without a retile_deghost object.
TEST_CASE("pr_ml doc 32: gapfill skip_base_subcells defaults to true and parses false")
{
    auto icfg = probe("d32_sub_default");
    CHECK(!icfg->default_configuration().isMember("retile_deghost_skip_base_subcells"));
    configure_tolerant(icfg, gapfill_cfg());
    auto cfg = icfg->default_configuration();
    REQUIRE(cfg.isMember("retile_deghost_skip_base_subcells"));
    CHECK(cfg["retile_deghost_skip_base_subcells"].asBool() == true);

    auto icfg2 = probe("d32_sub_off");
    auto c = gapfill_cfg();
    c["retile_deghost"]["skip_base_subcells"] = false;
    configure_tolerant(icfg2, c);
    cfg = icfg2->default_configuration();
    REQUIRE(cfg.isMember("retile_deghost_skip_base_subcells"));
    CHECK(cfg["retile_deghost_skip_base_subcells"].asBool() == false);
}
