// SteinerTruthLabeler (wcp-porting-img pr_ml/docs/23): truth labels on the
// Steiner graph from the per-point "truth" cloud of the original blobs.
//  * a vertex takes the majority key of the cluster points within join_radius
//    (count vote), with the nearest point as the fallback;
//  * a blob without a "truth" cloud gives tid -3;
//  * attributes come from the nearest point carrying the chosen key;
//  * the edge cloud mirrors the graph (src, dst, weight, length);
//  * a cluster without a graph is left alone, a vertex/row mismatch throws.
#include "WireCellClus/Facade.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <boost/graph/adjacency_list.hpp>

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
    // One blob on a line along x; with_truth adds a "truth" cloud with the given keys.
    void add_blob(Cluster& cl, const std::vector<double>& xs, const std::vector<int>& tids, bool with_truth)
    {
        Blob& b = cl.make_child();
        const int npts = xs.size();
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
        std::vector<double> ys(npts, 0), zs(npts, 0);
        b.value().local_pcs()["3d"] = Dataset({{"x", Array(xs)}, {"y", Array(ys)}, {"z", Array(zs)}});
        if (!with_truth) return;
        std::vector<int> pdg(npts), cls(npts, 0), iid(npts, 0), pf(npts), orig(npts);
        std::vector<double> dtraj(npts), frac(npts, 1.0);
        for (int i = 0; i < npts; ++i) {
            pdg[i] = tids[i] > 0 ? 13 : 0;
            pf[i] = tids[i];
            orig[i] = tids[i] > 0 ? 100 + tids[i] : tids[i];
            dtraj[i] = 0.1 * (i + 1);
        }
        b.value().local_pcs()["truth"] = Dataset({
            {"tid", Array(tids)}, {"orig_tid", Array(orig)}, {"pid_frac", Array(frac)},
            {"pdg", Array(pdg)}, {"cls", Array(cls)}, {"iid", Array(iid)}, {"pf_node", Array(pf)},
            {"d_traj", Array(dtraj)}});
    }

    void add_steiner(Cluster& cl, const std::vector<double>& xs, size_t nvertices,
                     const std::vector<std::tuple<int, int, double>>& edges)
    {
        const size_t n = xs.size();
        std::vector<double> ys(n, 0), zs(n, 0);
        std::vector<int> term(n, 0);
        cl.local_pcs()["steiner_pc"] = Dataset({{"x", Array(xs)}, {"y", Array(ys)}, {"z", Array(zs)},
                                                {"flag_steiner_terminal", Array(term)}});
        auto& g = cl.make_graph("steiner_graph", nvertices);
        for (const auto& [a, b, w] : edges) boost::add_edge(a, b, w, g);
    }

    IEnsembleVisitor::pointer make_visitor(const std::string& name)
    {
        PluginManager::instance().add("WireCellClus");
        auto icfg = Factory::lookup<IConfigurable>("SteinerTruthLabeler", name);
        icfg->configure(icfg->default_configuration());
        return Factory::find_tn<IEnsembleVisitor>("SteinerTruthLabeler:" + name);
    }

    template <typename T>
    std::vector<T> col(const Dataset& ds, const std::string& key)
    {
        auto sp = ds.get(key)->elements<T>();
        return std::vector<T>(sp.begin(), sp.end());
    }
}

TEST_CASE("steiner truth labeler: defaults")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("SteinerTruthLabeler", "stl_defaults");
    auto cfg = icfg->default_configuration();
    CHECK(cfg["grouping"].asString() == "live");
    CHECK(cfg["graph"].asString() == "steiner_graph");
    CHECK(cfg["truth_pc"].asString() == "truth");
    CHECK(cfg["out_nodes"].asString() == "steiner_truth");
    CHECK(cfg["out_edges"].asString() == "steiner_edges");
    CHECK(cfg["join_radius"].asDouble() == doctest::Approx(0.5 * units::cm));
}

TEST_CASE("steiner truth labeler: vote, fallback, no truth, edges")
{
    Points::node_t root;
    auto& ens = *root.value.facade<Ensemble>();
    auto& live = ens.make_grouping("live");
    Cluster& cl = live.make_child();
    cl.set_ident(1);
    add_blob(cl, {0, 10, 20}, {5, 5, -1}, true);
    add_blob(cl, {100, 104, 108}, {7, 7, 9}, true);
    add_blob(cl, {200}, {0}, false);
    // vertices: near blob 1's first point; between 7,7,9 (vote 7); on the truth-less blob; 30 mm from everything
    add_steiner(cl, {1, 104, 200, 50}, 4, {{0, 3, 2.5}, {1, 2, 3.5}});
    Cluster& bare = live.make_child();   // no graph: untouched
    bare.set_ident(2);
    add_blob(bare, {500}, {3}, true);

    make_visitor("stl_main")->visit(ens);

    REQUIRE(cl.has_pc("steiner_truth"));
    const auto& st = cl.get_pc("steiner_truth");
    CHECK(col<int>(st, "tid") == std::vector<int>{5, 7, -3, -1});
    CHECK(col<int>(st, "npts_join") == std::vector<int>{1, 3, 1, 0});
    const auto frac = col<double>(st, "tid_frac");
    CHECK(frac[0] == doctest::Approx(1.0));
    CHECK(frac[1] == doctest::Approx(2.0 / 3.0));
    CHECK(frac[3] == doctest::Approx(0.0));
    const auto dn = col<double>(st, "d_near");
    CHECK(dn[0] == doctest::Approx(1.0));
    CHECK(dn[1] == doctest::Approx(0.0));
    CHECK(dn[3] == doctest::Approx(30.0));
    CHECK(col<int>(st, "tid_near") == std::vector<int>{5, 7, -3, -1});
    CHECK(col<int>(st, "orig_tid") == std::vector<int>{105, 107, -3, -1});
    CHECK(col<int>(st, "pdg") == std::vector<int>{13, 13, 0, 0});
    // the representative of vertex 1 is the nearest key-7 point (x = 104, its blob's second point)
    CHECK(col<double>(st, "d_traj")[1] == doctest::Approx(0.2));
    CHECK(col<double>(st, "x_raw") == std::vector<double>{0, 104, 200, 20});
    for (const char* k : {"dx_cor", "dy_cor", "dz_cor"}) {
        for (double v : col<double>(st, k)) CHECK(v == 0.0);
    }

    REQUIRE(cl.has_pc("steiner_edges"));
    const auto& se = cl.get_pc("steiner_edges");
    CHECK(col<int>(se, "src") == std::vector<int>{0, 1});
    CHECK(col<int>(se, "dst") == std::vector<int>{3, 2});
    CHECK(col<double>(se, "weight") == std::vector<double>{2.5, 3.5});
    const auto len = col<double>(se, "length");
    CHECK(len[0] == doctest::Approx(49.0));
    CHECK(len[1] == doctest::Approx(96.0));

    CHECK_FALSE(bare.has_pc("steiner_truth"));
    CHECK_FALSE(bare.has_pc("steiner_edges"));
}

TEST_CASE("steiner truth labeler: vertex count must equal the steiner_pc rows")
{
    Points::node_t root;
    auto& ens = *root.value.facade<Ensemble>();
    auto& live = ens.make_grouping("live");
    Cluster& cl = live.make_child();
    add_blob(cl, {0, 10}, {5, 5}, true);
    add_steiner(cl, {0, 10}, 3, {{0, 1, 1.0}});
    CHECK_THROWS_AS(make_visitor("stl_mismatch")->visit(ens), ValueError);
}
