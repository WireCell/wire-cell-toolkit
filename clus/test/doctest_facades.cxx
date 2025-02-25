#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Logging.h"

#include "WireCellClus/Facade_Cluster.h"

using namespace WireCell;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::PointCloud::Facade;
using spdlog::debug;

TEST_CASE("clustering facade scalar")
{
    Points::node_t node;
    Cluster* cluster = node.value.facade<Cluster>();

    int no1 = cluster->ident(-1);
    CHECK(no1 == -1);           // should not yet exist
    int no2 = cluster->ident(-2);
    CHECK(no2 == -2);           // should still not yet exist

    cluster->set_ident(42);
    int cid = cluster->ident(-1);
    CHECK(cid == 42);
    // debug("no1={} no2={} cide={}", no1, no2, cid);
}
