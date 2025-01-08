#include "WireCellUtil/Graph.h"

#include "WireCellUtil/doctest.h"

#include "WireCellUtil/Logging.h"

using namespace WireCell;
using spdlog::debug;

using cluster_connectivity_graph_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, int>;

TEST_CASE("cluster connectivity graph") {

    cluster_connectivity_graph_t ccg;

    for (int ind=0; ind<100; ++ind) {
        auto desc = boost::add_vertex(ind, ccg);
        REQUIRE(desc == ind);
        REQUIRE(ccg[desc] == ind);
    }
    


}
