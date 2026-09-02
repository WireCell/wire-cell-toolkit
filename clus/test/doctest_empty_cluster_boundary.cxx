// doc pdvd/25 M3: an EMPTY cluster must not abort the steiner retile.
//
// PDVD run 039252 evt 4: ImproveCluster_2 retiled a protect_bundle fragment
// into a cluster with zero points and Cluster::get_two_boundary_wcps() then
// called point3d(0) -> std::out_of_range -> wire-cell rc=134.  The primitive
// now answers a defined degenerate pair for npoints()==0, and the retile /
// steiner stages skip such clusters (they are logged, not fitted).  Any
// non-empty cluster takes exactly the former path.

#include "WireCellUtil/doctest.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"

using namespace WireCell;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;

TEST_CASE("empty cluster: get_two_boundary_wcps returns a defined pair instead of throwing")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    Cluster& cl = g->make_child();
    REQUIRE(cl.npoints() == 0);
    std::pair<geo_point_t, geo_point_t> pr;
    CHECK_NOTHROW(pr = cl.get_two_boundary_wcps());
    CHECK(pr.first == pr.second);
    CHECK(pr.first.x() == 0.0);
}
