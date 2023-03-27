#include "WireCellAux/TensorDM.h"
#include "WireCellUtil/Testing.h"
using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::Aux::TensorDM;


void test_empty()
{
    PointGraph pg;
    const std::string datapath = "blah";
    auto tens = as_tensors(pg, datapath);
    auto pg2 = as_pointgraph(tens, datapath);
}

using graph_t = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS>;
void test_something()
{
    Dataset nodes({
            {"one", Array({1  ,2  ,3  })},
            {"two", Array({1.1,2.2,3.3})},
        });
    Dataset edges({
            {"tails", Array({0, 1, 2})},
            {"heads", Array({1, 2, 0})},
        });

    PointGraph pg(nodes, edges);


    {
        auto g = pg.boost_graph<graph_t>();
        Assert(boost::num_vertices(g) == 3);
        Assert(boost::num_edges(g) == 3);
    }
    
    const std::string datapath = "blah";
    auto tens = as_tensors(pg, datapath);
    {
        Assert(tens[0]->metadata()["datapath"].asString() == datapath);
        Assert(tens.size() == 1 + 1+2 + 1+2);
    }

    auto pg2 = as_pointgraph(tens, datapath);
    {
        auto g = pg2.boost_graph<graph_t>();
        Assert(boost::num_vertices(g) == 3);
        Assert(boost::num_edges(g) == 3);
    }
    
}

int main()
{
    test_empty();
    test_something();

    return 0;
}
