// ClusteringBundleNoLight (wcfm/docs/25): the whole light-less event becomes
// one bundle -- one main, every other cluster associated, one gid, one t0,
// written on every cluster.
//  * main selection without detector geometry: every get_length() is 0 here,
//    so the tie-break decides (more points, then smaller cluster id);
//  * flags, gid and t0 left by an upstream stage are overwritten;
//  * flag_beam_flash is off by default and sets beam_flash when on;
//  * an empty or missing grouping is skipped, not an error;
//  * bundle_gid outside [0, 1000000) is refused at configure time.
#include "WireCellClus/Facade.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <string>
#include <vector>

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using fa_float_t = WireCell::Clus::Facade::float_t;
using fa_int_t = WireCell::Clus::Facade::int_t;

namespace {
    // A cluster with one blob of npts points.  The blob's "scalar" PC is complete
    // with EMPTY wire/slice ranges, so get_length() is 0 without any detector
    // geometry (doctest_nu_band_veto.cxx add_geom_blob); "3d" gives npoints().
    Cluster& add_cluster(Grouping& g, int ident, int npts)
    {
        Cluster& cl = g.make_child();
        cl.set_ident(ident);
        Blob& b = cl.make_child();
        std::vector<fa_int_t> zero{0};
        Dataset scalar({
            {"charge", Array(std::vector<fa_float_t>{0})},
            {"center_x", Array(std::vector<fa_float_t>{0})},
            {"center_y", Array(std::vector<fa_float_t>{0})},
            {"center_z", Array(std::vector<fa_float_t>{0})},
            {"wpid", Array(zero)}, {"npoints", Array(std::vector<fa_int_t>{npts})},
            {"slice_index_min", Array(zero)}, {"slice_index_max", Array(zero)},
            {"u_wire_index_min", Array(zero)}, {"u_wire_index_max", Array(zero)},
            {"v_wire_index_min", Array(zero)}, {"v_wire_index_max", Array(zero)},
            {"w_wire_index_min", Array(zero)}, {"w_wire_index_max", Array(zero)},
            {"max_wire_interval", Array(zero)}, {"min_wire_interval", Array(zero)},
            {"max_wire_type", Array(zero)}, {"min_wire_type", Array(zero)},
        });
        b.value().local_pcs()["scalar"] = scalar;
        std::vector<double> xs(npts), ys(npts), zs(npts);
        for (int i = 0; i < npts; ++i) { xs[i] = i; ys[i] = 0; zs[i] = 0; }
        b.value().local_pcs()["3d"] = Dataset({{"x", Array(xs)}, {"y", Array(ys)}, {"z", Array(zs)}});
        return cl;
    }

    IEnsembleVisitor::pointer make_visitor(const std::string& name, Configuration over = Configuration())
    {
        PluginManager::instance().add("WireCellClus");
        auto icfg = Factory::lookup<IConfigurable>("ClusteringBundleNoLight", name);
        auto cfg = icfg->default_configuration();
        for (const auto& key : over.getMemberNames()) cfg[key] = over[key];
        icfg->configure(cfg);
        return Factory::find_tn<IEnsembleVisitor>("ClusteringBundleNoLight:" + name);
    }
}

TEST_CASE("bundle no light: defaults")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("ClusteringBundleNoLight", "bnl_defaults");
    auto cfg = icfg->default_configuration();
    CHECK(cfg["grouping"].asString() == "live");
    CHECK(cfg["bundle_gid"].asInt() == 0);
    CHECK(cfg["cluster_t0"].asDouble() == 0.0);
    CHECK(cfg["flag_beam_flash"].asBool() == false);
}

TEST_CASE("bundle no light: one main, the rest associated, gid and t0 on all")
{
    Points::node_t root;
    auto& ens = *root.value.facade<Ensemble>();
    auto& live = ens.make_grouping("live");
    Cluster& a = add_cluster(live, 11, 2);
    Cluster& b = add_cluster(live, 7, 5);
    Cluster& c = add_cluster(live, 3, 5);   // ties b on points, smaller id => main
    // what Q/L matching with no flash (or an upstream flag_mains) would have left
    a.set_flag(Flags::main_cluster, 1);
    a.set_scalar<int>("matched_flash_gid", -1);
    a.set_cluster_t0(-1e12);
    b.set_flag(Flags::associated_cluster, 0);

    auto vis = make_visitor("bnl_main");
    vis->visit(ens);

    CHECK(c.get_flag(Flags::main_cluster) == 1);
    CHECK(c.get_flag(Flags::associated_cluster) == 0);
    for (Cluster* cl : {&a, &b}) {
        CHECK(cl->get_flag(Flags::main_cluster) == 0);
        CHECK(cl->get_flag(Flags::associated_cluster) == 1);
    }
    for (Cluster* cl : {&a, &b, &c}) {
        CHECK(cl->get_scalar<int>("matched_flash_gid", -99) == 0);
        CHECK(cl->get_cluster_t0() == 0.0);
        CHECK(cl->get_flag(Flags::beam_flash, -1) == -1);   // not written by default
    }
}

TEST_CASE("bundle no light: more points wins over id; knobs gid, t0, beam_flash")
{
    Points::node_t root;
    auto& ens = *root.value.facade<Ensemble>();
    auto& live = ens.make_grouping("live");
    Cluster& a = add_cluster(live, 1, 3);
    Cluster& b = add_cluster(live, 9, 4);   // more points => main despite the larger id

    Configuration over;
    over["bundle_gid"] = 42;
    over["cluster_t0"] = 5.0 * units::us;
    over["flag_beam_flash"] = true;
    auto vis = make_visitor("bnl_knobs", over);
    vis->visit(ens);

    CHECK(b.get_flag(Flags::main_cluster) == 1);
    CHECK(a.get_flag(Flags::main_cluster) == 0);
    CHECK(a.get_flag(Flags::associated_cluster) == 1);
    for (Cluster* cl : {&a, &b}) {
        CHECK(cl->get_scalar<int>("matched_flash_gid", -99) == 42);
        CHECK(cl->get_cluster_t0() == doctest::Approx(5.0 * units::us));
        CHECK(cl->get_flag(Flags::beam_flash) == 1);
    }
}

TEST_CASE("bundle no light: empty or missing grouping is skipped")
{
    Points::node_t root;
    auto& ens = *root.value.facade<Ensemble>();
    auto vis = make_visitor("bnl_empty");
    CHECK_NOTHROW(vis->visit(ens));          // no "live" grouping
    ens.make_grouping("live");
    CHECK_NOTHROW(vis->visit(ens));          // "live" with no clusters
}

TEST_CASE("bundle no light: gid outside [0, 1000000) is refused")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("ClusteringBundleNoLight", "bnl_badgid");
    auto cfg = icfg->default_configuration();
    cfg["bundle_gid"] = -1;
    CHECK_THROWS(icfg->configure(cfg));
    cfg["bundle_gid"] = 1000000;
    CHECK_THROWS(icfg->configure(cfg));
}
