#include "WireCellUtil/PointGraph.h"
#include "WireCellUtil/Testing.h"

using namespace WireCell;

using graph_t = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS>;

PointGraph make_empty()
{
    PointGraph pg;
    return pg;
}

void test_empty()
{
    PointGraph pg1 = make_empty();
    Assert(pg1.nodes().size_major() == 0);
    Assert(pg1.edges().size_major() == 0);
    
    auto g = pg1.boost_graph<graph_t>();
    Assert(boost::num_vertices(g) == 0);
}

int main()
{
    test_empty();
    return 0;
}
