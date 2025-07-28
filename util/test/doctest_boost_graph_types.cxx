#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Graph.h"
#include "WireCellUtil/Type.h"

using GraphSSU = boost::adjacency_list<boost::setS, boost::setS, boost::directedS>;
using GraphVVU = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS>;
using GraphLLU = boost::adjacency_list<boost::listS, boost::listS, boost::directedS>;
using GraphHHU = boost::adjacency_list<boost::hash_setS, boost::hash_setS, boost::directedS>;


using namespace WireCell;

TEST_CASE("boost graph types")
{
    using ss = boost::graph_traits<GraphSSU>;
    using vv = boost::graph_traits<GraphVVU>;
    using ll = boost::graph_traits<GraphLLU>;
    using hh = boost::graph_traits<GraphHHU>;

    std::cerr
        << " set  vertex=" << type<ss::vertex_descriptor>() << "\n"
        << " set  edge=  " << type<ss::edge_descriptor>() << "\n"
        << " hash vertex=" << type<hh::vertex_descriptor>() << "\n"
        << " hash edge=  " << type<hh::edge_descriptor>() << "\n"
        << " vec  vertex=" << type<vv::vertex_descriptor>() << "\n"
        << " vec  edge=  " << type<vv::edge_descriptor>() << "\n"
        << " list vertex=" << type<ll::vertex_descriptor>() << "\n"
        << " list edge=  " << type<ll::edge_descriptor>() << "\n";
        
}
