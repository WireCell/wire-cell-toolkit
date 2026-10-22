// PRTruthExport (wcp-porting-img pr_ml/docs/23 sec 11): the PR chain's point
// assignment exported as cluster clouds.
//  * a segment's associate_points rows are joined to the cluster points by the
//    stored index when the point there equals the row (src 0), else by the
//    nearest cluster point (src 1, d_join);
//  * a segment without the cloud gets a pr_segments row and no point rows;
//  * a segment whose cluster is not in the grouping is skipped;
//  * seg_id is cluster_id*1000 + segment id (graph index when unset);
//  * no TrackFitting on the grouping => nothing written.
#include "WireCellClus/Facade.h"
#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/DynamicPointCloud.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/PointTree.h"
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
    // One blob on a line along x.
    void add_blob(Cluster& cl, const std::vector<double>& xs)
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
    }

    // An associate_points cloud with the given x positions and stored indices.
    void set_assoc(PR::SegmentPtr seg, const std::vector<double>& xs, const std::vector<size_t>& idx)
    {
        DPCBatch batch;
        for (double x : xs) {
            batch.add_point(x, 0, 0, 0, nullptr, nullptr, {0, 0, 0}, {0, 0, 0});
            for (int p = 0; p < 3; ++p) batch.end_plane();
        }
        auto dpc = std::make_shared<DynamicPointCloud>(std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>>{});
        dpc->add_points(std::move(batch));
        seg->dpcloud("associate_points", dpc);
        if (!idx.empty()) seg->set_global_indices("associate_points", idx);
    }

    IEnsembleVisitor::pointer make_visitor(const std::string& name)
    {
        PluginManager::instance().add("WireCellClus");
        auto icfg = Factory::lookup<IConfigurable>("PRTruthExport", name);
        icfg->configure(icfg->default_configuration());
        return Factory::find_tn<IEnsembleVisitor>("PRTruthExport:" + name);
    }

    template <typename T>
    std::vector<T> col(const Dataset& ds, const std::string& key)
    {
        auto sp = ds.get(key)->elements<T>();
        return std::vector<T>(sp.begin(), sp.end());
    }
}

TEST_CASE("pr truth export: defaults")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("PRTruthExport", "pte_defaults");
    auto cfg = icfg->default_configuration();
    CHECK(cfg["grouping"].asString() == "live");
    CHECK(cfg["assoc_cloud"].asString() == "associate_points");
    CHECK(cfg["out_segments"].asString() == "pr_segments");
    CHECK(cfg["out_points"].asString() == "pr_seg_points");
}

TEST_CASE("pr truth export: stored index, stale index, missing cloud, foreign cluster")
{
    Points::node_t root;
    auto& ens = *root.value.facade<Ensemble>();
    auto& live = ens.make_grouping("live");
    Cluster& cl = live.make_child();
    cl.set_ident(3);
    add_blob(cl, {0, 10, 20});        // points 0, 1, 2
    add_blob(cl, {100, 104, 108});    // points 3, 4, 5
    Cluster& bare = live.make_child();   // no segment: untouched
    bare.set_ident(4);
    add_blob(bare, {500});

    Points::node_t root2;                // a cluster of another ensemble
    auto& other = root2.value.facade<Ensemble>()->make_grouping("live");
    Cluster& foreign = other.make_child();
    foreign.set_ident(9);
    add_blob(foreign, {700});

    // the fitter syncs its geometry from the segments' clusters when a graph is
    // attached; attach the empty graph and fill it afterwards (the visitor
    // only reads the graph)
    auto g = std::make_shared<PR::Graph>();
    auto tf = std::make_shared<TrackFitting>();
    tf->add_graph(g);
    live.set_track_fitting(tf);
    std::vector<PR::VertexPtr> v;
    for (int i = 0; i < 5; ++i) v.push_back(std::make_shared<PR::Vertex>());
    std::vector<PR::SegmentPtr> s;
    for (int i = 0; i < 4; ++i) {
        auto seg = std::make_shared<PR::Segment>();
        PR::add_segment(*g, seg, v[i], v[i + 1]);
        s.push_back(seg);
    }
    s[0]->cluster(&cl);
    s[1]->cluster(&cl);
    s[2]->cluster(&cl);
    s[3]->cluster(&foreign);
    s[1]->set_id(7);
    set_assoc(s[0], {0, 10}, {0, 1});           // stored indices right
    set_assoc(s[1], {104, 108, 101}, {5, 5, 3}); // first stale, second right, third off the points
    set_assoc(s[3], {700}, {0});                 // foreign: skipped
    // s[2]: no associate_points cloud

    make_visitor("pte_main")->visit(ens);

    REQUIRE(cl.has_pc("pr_segments"));
    const auto& ss = cl.get_pc("pr_segments");
    const int id0 = 3000 + static_cast<int>(s[0]->get_graph_index());
    const int id2 = 3000 + static_cast<int>(s[2]->get_graph_index());
    CHECK(col<int>(ss, "seg_id") == std::vector<int>{id0, 3007, id2});
    CHECK(col<int>(ss, "nu_index") == std::vector<int>{0, 0, 0});
    CHECK(col<int>(ss, "n_cloud") == std::vector<int>{2, 3, 0});
    CHECK(col<int>(ss, "n_index") == std::vector<int>{2, 3, 0});
    CHECK(col<int>(ss, "n_exact") == std::vector<int>{2, 1, 0});

    REQUIRE(cl.has_pc("pr_seg_points"));
    const auto& sp = cl.get_pc("pr_seg_points");
    CHECK(col<int>(sp, "seg_id") == std::vector<int>{id0, id0, 3007, 3007, 3007});
    CHECK(col<int>(sp, "p_index") == std::vector<int>{0, 1, 4, 5, 3});
    CHECK(col<int>(sp, "src") == std::vector<int>{0, 0, 1, 0, 1});
    const auto dj = col<double>(sp, "d_join");
    CHECK(dj[2] == doctest::Approx(0.0));
    CHECK(dj[4] == doctest::Approx(1.0));
    CHECK(col<double>(sp, "x_raw") == std::vector<double>{0, 10, 104, 108, 100});

    CHECK_FALSE(bare.has_pc("pr_segments"));
    CHECK_FALSE(bare.has_pc("pr_seg_points"));
    CHECK_FALSE(foreign.has_pc("pr_segments"));
}

TEST_CASE("pr truth export: no TrackFitting, nothing written")
{
    Points::node_t root;
    auto& ens = *root.value.facade<Ensemble>();
    auto& live = ens.make_grouping("live");
    Cluster& cl = live.make_child();
    add_blob(cl, {0, 10});
    make_visitor("pte_none")->visit(ens);
    CHECK_FALSE(cl.has_pc("pr_segments"));
    CHECK_FALSE(cl.has_pc("pr_seg_points"));
}
