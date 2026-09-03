// doc pdvd/26: examine_partial_identical_segments never terminated on PDVD
// 039349/14 and 039349/53 (48 min, killed).  The loop splits the shared trunk
// of two overlapping segments at the steiner point closest to the far end of
// the overlap.  On those two events that overlap lay in a charge gap with no
// steiner points, so the closest steiner point was the VERTEX ITSELF: the
// branch cloned the vertex in place, moved both segments to the clone, and the
// next pass found the identical overlap on the clone -- one vertex and one
// zero-length segment added per iteration, forever (nv 5 -> 145 in 700 lines
// of instrumentation, doc pdvd/26 sec 3).
//
// The prototype (NeutrinoID_proto_vertex.h, get_closest_wcpoint(max_point))
// has the same logic and the same fixed point; uBooNE/SBND steiner clouds are
// dense along every track, so the closest point is always near max_point and
// the case never arose.  The port now asks partial_identical_split_point() for
// the split point and skips the vertex when there is none.

#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::PR;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;

namespace {

// A cluster whose steiner cloud is the given x positions on the x axis.
Facade::Cluster& make_cluster(Points::node_t& root, const std::vector<double>& xs)
{
    auto* g = root.value.facade<Facade::Grouping>();
    auto& cl = g->make_child();
    Dataset steiner;
    steiner.add("x", Array(xs));
    steiner.add("y", Array(std::vector<double>(xs.size(), 0.0)));
    steiner.add("z", Array(std::vector<double>(xs.size(), 0.0)));
    cl.value().local_pcs()["steiner_pc"] = steiner;
    REQUIRE(cl.has_pc("steiner_pc"));
    return cl;
}

std::vector<double> grid(double lo, double hi, double step)
{
    std::vector<double> xs;
    for (double x = lo; x <= hi + 1e-9; x += step) xs.push_back(x);
    return xs;
}

}  // namespace

TEST_CASE("partial_identical_split_point: dense steiner cloud splits at the overlap's far end")
{
    Points::node_t root;
    auto& cl = make_cluster(root, grid(0.0, 30 * units::cm, 0.5 * units::cm));
    PatternAlgorithms pa;
    const Facade::geo_point_t vtx(0, 0, 0);
    const Facade::geo_point_t max_point(8.69 * units::cm, 0, 0);   // the production overlap length

    auto pt = pa.partial_identical_split_point(cl, max_point, vtx);
    REQUIRE(pt.has_value());
    CHECK(std::abs(pt->x() - 8.5 * units::cm) < 1e-6);            // nearest grid point
    CHECK((*pt - vtx).magnitude() > 5 * units::cm);                // a real split, away from the vertex
}

TEST_CASE("partial_identical_split_point: overlap in a charge gap -> the vertex itself is not a split point (PDVD 039349/14)")
{
    // Production: vtx (2732.6,-1189.0,866.1) mm, max_point 8.69 cm out along the
    // two segments, no steiner point between them -- the kNN returns the vertex.
    Points::node_t root;
    std::vector<double> xs = {0.0};
    for (double x : grid(20 * units::cm, 30 * units::cm, 0.5 * units::cm)) xs.push_back(x);
    auto& cl = make_cluster(root, xs);
    PatternAlgorithms pa;
    const Facade::geo_point_t max_point(8.69 * units::cm, 0, 0);

    SUBCASE("vertex exactly on the steiner point") {
        const Facade::geo_point_t vtx(0, 0, 0);
        CHECK_FALSE(pa.partial_identical_split_point(cl, max_point, vtx).has_value());
    }
    SUBCASE("vertex within 0.1 cm of the steiner point") {
        const Facade::geo_point_t vtx(0.05 * units::cm, 0.02 * units::cm, 0);
        CHECK_FALSE(pa.partial_identical_split_point(cl, max_point, vtx).has_value());
    }
    SUBCASE("a split beyond the gap is still allowed") {
        const Facade::geo_point_t vtx(0, 0, 0);
        const Facade::geo_point_t far_point(21.2 * units::cm, 0, 0);
        auto pt = pa.partial_identical_split_point(cl, far_point, vtx);
        REQUIRE(pt.has_value());
        CHECK(std::abs(pt->x() - 21.0 * units::cm) < 1e-6);
    }
}

TEST_CASE("partial_identical_split_point: no or point-less steiner cloud -> no split")
{
    PatternAlgorithms pa;
    const Facade::geo_point_t vtx(0, 0, 0);
    const Facade::geo_point_t max_point(8.69 * units::cm, 0, 0);

    SUBCASE("no steiner_pc at all") {
        Points::node_t root;
        auto* g = root.value.facade<Facade::Grouping>();
        auto& cl = g->make_child();
        CHECK_FALSE(pa.partial_identical_split_point(cl, max_point, vtx).has_value());
    }
    SUBCASE("present but point-less (pr11: size() counts arrays, not points)") {
        Points::node_t root;
        auto& cl = make_cluster(root, {});
        CHECK_FALSE(pa.partial_identical_split_point(cl, max_point, vtx).has_value());
    }
}

// ---------------------------------------------------------------------------
// doc pdvd/26 round 2
// ---------------------------------------------------------------------------

#include "WireCellClus/PRGraph.h"
#include "WireCellClus/ExaminerPassBudget.h"

TEST_CASE("closest_cluster_vertex: nearest wcpt of THIS cluster, ignoring other clusters' vertices")
{
    Points::node_t root;
    auto* g = root.value.facade<Facade::Grouping>();
    auto& cl_a = g->make_child();
    auto& cl_b = g->make_child();

    Graph graph;
    auto va0 = make_vertex(graph); va0->wcpt().point = Facade::geo_point_t(0, 0, 0);           va0->cluster(&cl_a);
    auto va1 = make_vertex(graph); va1->wcpt().point = Facade::geo_point_t(100 * units::cm, 0, 0); va1->cluster(&cl_a);
    auto vb0 = make_vertex(graph); vb0->wcpt().point = Facade::geo_point_t(1 * units::cm, 0, 0);  vb0->cluster(&cl_b);

    PatternAlgorithms pa;
    const Facade::geo_point_t q(2 * units::cm, 0, 0);

    SUBCASE("cluster A: the 2 cm vertex, not cluster B's 1 cm one") {
        auto [v, d] = pa.closest_cluster_vertex(graph, cl_a, q);
        REQUIRE(v);
        CHECK(v == va0);
        CHECK(std::abs(d - 2 * units::cm) < 1e-6);
    }
    SUBCASE("cluster B: its own vertex") {
        auto [v, d] = pa.closest_cluster_vertex(graph, cl_b, q);
        REQUIRE(v);
        CHECK(v == vb0);
        CHECK(std::abs(d - 1 * units::cm) < 1e-6);
    }
    SUBCASE("a cluster with no vertex in the graph") {
        auto& cl_c = g->make_child();
        auto [v, d] = pa.closest_cluster_vertex(graph, cl_c, q);
        CHECK(!v);
        CHECK(d > 1e8);
    }
}

TEST_CASE("ExaminerPassCounter: silent up to the budget, true once past it")
{
    ExaminerPassCounter c("doctest", 7);
    for (int i = 0; i < kExaminerPassBudget; ++i) {
        CHECK_FALSE(c.exceeded());
    }
    CHECK(c.passes == kExaminerPassBudget);
    CHECK(c.exceeded());                 // pass kExaminerPassBudget + 1 is refused
    CHECK(c.passes == kExaminerPassBudget + 1);
    // The budget must stay far above any pass count seen in production (doc
    // pdvd/26 sec 7 census); the number is pinned so a change is a decision.
    CHECK(kExaminerPassBudget == 1000);
}
