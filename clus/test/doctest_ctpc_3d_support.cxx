// doc pdvd/140 sec 3: the 3D support radius of the CTPC bridge validation is a value CreateSteinerGraph
// puts on the grouping only around its own graph build.  These tests pin what the gate relies on: the default is
// OFF (0), a set value reads back through a child cluster's grouping (what connect_graph_ctpc_with_reference
// reads), and Grouping::from() does NOT carry it -- a derived grouping (another builder of the same graph flavor)
// must stay on the legacy 2D-only validation.
#include "WireCellUtil/doctest.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellUtil/Units.h"

using namespace WireCell;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;

TEST_CASE("ctpc_3d_support: default off, set/get, scoped (not carried by Grouping::from)")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    CHECK(g->ctpc_3d_support() == 0.0);

    g->set_ctpc_3d_support(2 * units::cm);
    CHECK(g->ctpc_3d_support() == doctest::Approx(2 * units::cm));

    Cluster& cl = g->make_child();
    REQUIRE(cl.grouping() == g);
    CHECK(cl.grouping()->ctpc_3d_support() == doctest::Approx(2 * units::cm));

    Points::node_t root2;
    Grouping* g2 = root2.value.facade<Grouping>();
    g2->from(*g);
    CHECK(g2->ctpc_3d_support() == 0.0);

    g->set_ctpc_3d_support(0);
    CHECK(g->ctpc_3d_support() == 0.0);
}
