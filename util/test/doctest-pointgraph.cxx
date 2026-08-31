#include "WireCellUtil/PointGraph.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;

using graph_t = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS>;

TEST_SUITE("point graph") {

TEST_CASE("empty point graph") {
    PointGraph pg;
    CHECK(pg.nodes().size_major() == 0);
    CHECK(pg.edges().size_major() == 0);

    auto g = pg.boost_graph<graph_t>();
    CHECK(boost::num_vertices(g) == 0);
}

}  // TEST_SUITE("point graph")
