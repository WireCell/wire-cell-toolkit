// pr_ml doc 92 (Xin O9 / O13 / O18, 2026-11-25): the every-cluster assignment of an outside segment graph and its knobs.
//
// Pins (1) the legacy defaults of the new knobs -- TaggerCheckNeutrino pr_external_all_clusters / pr_external_seg_ids /
// pr_external_clustering_points = false, PRTruthExport seg_ext_id = false, Segment::ext_id() = -1 -- so that "knob absent =>
// byte-identical" stays true; (2) PR::assign_external_graph: every segment and clump goes to exactly one cluster, its nearest
// by the median point distance, when that median is at most max_med_dist, else to none (with its nearest cluster kept);
// (3) ties go to the earlier cluster of the list (the main, which TaggerCheckNeutrino puts first), so the main's items are
// the ones build_external_graph takes for the main (its "another cluster is nearer" test is strict).

#include "WireCellClus/Facade.h"
#include "WireCellClus/ExternalGraph.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
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
    // A cluster with one blob: npts points on the line (i, y0, 0) mm, i = 0 .. npts-1 (doctest_clustering_bundle_no_light.cxx).
    Cluster& add_cluster(Grouping& g, int ident, int npts, double y0)
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
        for (int i = 0; i < npts; ++i) { xs[i] = i * units::mm; ys[i] = y0 * units::mm; zs[i] = 0; }
        b.value().local_pcs()["3d"] = Dataset({{"x", Array(xs)}, {"y", Array(ys)}, {"z", Array(zs)}});
        return cl;
    }

    // a segment / clump of 5 points along x at height y (mm)
    std::vector<WireCell::Point> line_at(double y)
    {
        std::vector<WireCell::Point> p;
        for (int i = 2; i <= 6; ++i) p.emplace_back(i * units::mm, y * units::mm, 0);
        return p;
    }

    PR::ExternalGraph graph_of(const std::vector<double>& seg_y, const std::vector<double>& clump_y)
    {
        PR::ExternalGraph eg;
        int vid = 0;
        for (size_t i = 0; i < seg_y.size(); ++i) {
            PR::ExternalGraph::Seg s;
            s.id = 100 + (int) i;
            s.pts = line_at(seg_y[i]);
            s.v0 = vid++;
            s.v1 = vid++;
            eg.vertices.push_back({s.v0, s.pts.front()});
            eg.vertices.push_back({s.v1, s.pts.back()});
            s.type = "track";
            eg.segments.push_back(s);
        }
        for (size_t i = 0; i < clump_y.size(); ++i) {
            PR::ExternalGraph::Clump c;
            c.id = (int) i;
            c.pts = line_at(clump_y[i]);
            c.type = "shower";
            eg.clumps.push_back(c);
        }
        return eg;
    }
}

TEST_CASE("doc 92: the every-cluster knobs default off")
{
    PluginManager::instance().add("WireCellClus");
    auto tcn = Factory::lookup<IConfigurable>("TaggerCheckNeutrino", "extassign_probe");
    REQUIRE(tcn);
    auto cfg = tcn->default_configuration();
    for (const char* k : {"pr_external_all_clusters", "pr_external_seg_ids", "pr_external_clustering_points"}) {
        REQUIRE(cfg.isMember(k));
        CHECK(cfg[k].asBool() == false);
    }
    auto pte = Factory::lookup<IConfigurable>("PRTruthExport", "extassign_probe");
    REQUIRE(pte);
    auto pcfg = pte->default_configuration();
    REQUIRE(pcfg.isMember("seg_ext_id"));
    CHECK(pcfg["seg_ext_id"].asBool() == false);
    PR::Segment seg;
    CHECK(seg.ext_id() == -1);
    seg.set_ext_id(7);
    CHECK(seg.ext_id() == 7);
    CHECK(seg.id() == -1);   // ext_id is separate from the display id
}

TEST_CASE("doc 92: assign_external_graph partitions items by nearest cluster within max_med_dist")
{
    Points::node_t root;
    auto& ens = *root.value.facade<Ensemble>();
    auto& live = ens.make_grouping("live");
    Cluster& m = add_cluster(live, 1, 10, 0);     // main, y = 0
    Cluster& b = add_cluster(live, 2, 10, 100);   // y = 100
    Cluster& c = add_cluster(live, 3, 10, 200);   // y = 200
    std::vector<const Cluster*> cl{&m, &b, &c};
    // segments at y = 1 (main), 99 (b), 203 (c), 150 (nearest b or c at 50 mm: too far), 1000 (too far, nearest c)
    const auto eg = graph_of({1, 99, 203, 150, 1000}, {2, 198, 500});
    const auto as = PR::assign_external_graph(eg, cl, 6 * units::mm);
    REQUIRE(as.seg_cluster.size() == 5);
    CHECK(as.seg_cluster[0] == 0);
    CHECK(as.seg_cluster[1] == 1);
    CHECK(as.seg_cluster[2] == 2);
    CHECK(as.seg_cluster[3] == -1);
    CHECK(as.seg_cluster[4] == -1);
    CHECK(as.seg_nearest[4] == 2);
    CHECK(as.seg_md[0] == doctest::Approx(1 * units::mm));
    CHECK(as.seg_md[2] == doctest::Approx(3 * units::mm));
    CHECK(as.seg_md[4] == doctest::Approx(800 * units::mm));
    REQUIRE(as.clump_cluster.size() == 3);
    CHECK(as.clump_cluster[0] == 0);
    CHECK(as.clump_cluster[1] == 2);
    CHECK(as.clump_cluster[2] == -1);
    // every item goes to at most one cluster, and to the main only when the main is its nearest
    for (size_t i = 0; i < as.seg_cluster.size(); ++i) {
        if (as.seg_cluster[i] >= 0) CHECK(as.seg_cluster[i] == as.seg_nearest[i]);
    }
}

TEST_CASE("doc 92: assign_external_graph gives a tie to the earlier cluster (the main)")
{
    Points::node_t root;
    auto& ens = *root.value.facade<Ensemble>();
    auto& live = ens.make_grouping("live");
    Cluster& m = add_cluster(live, 1, 10, 0);    // y = 0
    Cluster& b = add_cluster(live, 2, 10, 8);    // y = 8: a segment at y = 4 is 4 mm from both
    const auto eg = graph_of({4}, {4});
    const auto as1 = PR::assign_external_graph(eg, {&m, &b}, 6 * units::mm);
    CHECK(as1.seg_cluster[0] == 0);
    CHECK(as1.clump_cluster[0] == 0);
    const auto as2 = PR::assign_external_graph(eg, {&b, &m}, 6 * units::mm);
    CHECK(as2.seg_cluster[0] == 0);   // index 0 is now b: the earlier one wins, whichever it is
    // a hair nearer b: b wins
    const auto eg3 = graph_of({4.01}, {});
    const auto as3 = PR::assign_external_graph(eg3, {&m, &b}, 6 * units::mm);
    CHECK(as3.seg_cluster[0] == 1);
}
