#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include "WireCellClus/Facade_Cluster.h"

using namespace WireCell;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;
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
TEST_CASE("clustering facade graphs")
{
    Points::node_t node;
    Cluster* cluster = node.value.facade<Cluster>();

    CHECK(cluster->graph_store().size() == 0);

    auto& gr1 = cluster->make_graph("test1");
    CHECK(cluster->graph_store().size() == 1);
    CHECK(boost::num_vertices(gr1) == 0);

    auto& gr2 = cluster->make_graph("test2", 10);
    CHECK(cluster->graph_store().size() == 2);
    CHECK(boost::num_vertices(gr2) == 10);

    auto gr = cluster->take_graph("test2");
    CHECK(cluster->graph_store().size() == 1);
    CHECK(boost::num_vertices(gr) == 10);    

    cluster->give_graph("test1", std::move(gr));
    // gr at this point is in some undefined state
    auto& gr3 = cluster->get_graph("test1");
    CHECK(cluster->graph_store().size() == 1);
    CHECK(boost::num_vertices(gr3) == 10);    

    auto& gr4 = cluster->get_graph("new one");
    CHECK(cluster->graph_store().size() == 2);
    CHECK(boost::num_vertices(gr4) == 0);    

}
