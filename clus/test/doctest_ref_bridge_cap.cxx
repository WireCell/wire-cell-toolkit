// doc pdvd/138: the reference-bridge cap knob of connect_graph_with_reference is a
// grouping-wide setting (Grouping::set_ref_bridge_cap, read through
// cluster.grouping()).  These tests pin the plumbing the gate relies on: the
// default is OFF (0/0), a set value reads back, and Grouping::from() carries it
// to a derived grouping (the retile's shadow grouping) exactly as it carries
// ctpc_aniso_metric -- otherwise the base graph built on the retiled copy would
// silently run without the cap while the compiled config says it is on.
#include "WireCellUtil/doctest.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellUtil/Units.h"

using namespace WireCell;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;

TEST_CASE("ref_bridge_cap: default off, set/get, carried by Grouping::from")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    CHECK(g->ref_bridge_cap_cm() == 0.0);
    CHECK(g->ref_bridge_cap_pts() == 0);

    g->set_ref_bridge_cap(20 * units::cm, 20);
    CHECK(g->ref_bridge_cap_cm() == doctest::Approx(20 * units::cm));
    CHECK(g->ref_bridge_cap_pts() == 20);

    Points::node_t root2;
    Grouping* g2 = root2.value.facade<Grouping>();
    g2->from(*g);
    CHECK(g2->ref_bridge_cap_cm() == doctest::Approx(20 * units::cm));
    CHECK(g2->ref_bridge_cap_pts() == 20);

    // a child cluster answers through its grouping, which is what connect_graph_with_reference reads
    Cluster& cl = g->make_child();
    REQUIRE(cl.grouping() == g);
    CHECK(cl.grouping()->ref_bridge_cap_pts() == 20);

    // off again = legacy
    g->set_ref_bridge_cap(0, 0);
    CHECK(g->ref_bridge_cap_cm() == 0.0);
    CHECK(g->ref_bridge_cap_pts() == 0);
}
