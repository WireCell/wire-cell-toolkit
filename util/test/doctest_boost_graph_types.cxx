#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Graph.h"
#include "WireCellUtil/Type.h"
#include <type_traits>

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
        << " null vertex=" << ss::null_vertex() << "\n" // (void*)0

        << " hash vertex=" << type<hh::vertex_descriptor>() << "\n"
        << " hash edge=  " << type<hh::edge_descriptor>() << "\n"
        << " null vertex=" << hh::null_vertex() << "\n" // (void*)0

        << " vec  vertex=" << type<vv::vertex_descriptor>() << "\n"
        << " vec  edge=  " << type<vv::edge_descriptor>() << "\n"
        << " null vertex=" << vv::null_vertex() << "\n" // (size_t)(2^64-1)

        << " list vertex=" << type<ll::vertex_descriptor>() << "\n"
        << " list edge=  " << type<ll::edge_descriptor>() << "\n"
        << " null vertex=" << vv::null_vertex() << "\n" // (void*)(2^64-1)
        ;
        
    // static_assert(std::is_same_v<ss::vertex_descriptor, int>, 
    //               "Vertex descriptor should be int");

}

struct NodeBundle {
    std::string node_label;
};
struct EdgeBundle {
    std::string edge_label;
};
struct GraphBundle {
    std::string graph_label;
};



// #define DT_BOOST_GRAPH_TYPES_USE_VEC 1
#ifdef DT_BOOST_GRAPH_TYPES_USE_VEC
using Graph = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::undirectedS,
    NodeBundle,
    EdgeBundle,
    // boost::property<boost::vertex_index_t, int, NodeBundle>,
    // boost::property<boost::edge_index_t, int, EdgeBundle>,
    GraphBundle>;
#else
using Graph = boost::adjacency_list<
    boost::setS, boost::setS, boost::undirectedS,
    boost::property<boost::vertex_index_t, int, NodeBundle>,
    boost::property<boost::edge_index_t, int, EdgeBundle>,
    GraphBundle>;
#endif

TEST_CASE("boost graph vertex index")
{
    Graph g;
    auto v0 = boost::add_vertex(g);
    g[v0].node_label = "v0";
    auto v1 = boost::add_vertex(g);
    g[v1].node_label = "v1";


    {
        auto index_map = get(boost::vertex_index, g);
        auto idx0 = index_map[v0];
        auto idx1 = index_map[v1];
    
        std::cerr << "v0:" << v0 << " i0:" << idx0 << "\n";
        std::cerr << "v1:" << v1 << " i1:" << idx1 << "\n";
    }

// can only set setS, not vecS
#ifndef DT_BOOST_GRAPH_TYPES_USE_VEC
    {
        auto index_map1 = get(boost::vertex_index, g);
        index_map1[v0] = 100;
        index_map1[v1] = 200;
        auto index_map = get(boost::vertex_index, g);
        auto idx0 = index_map[v0];
        auto idx1 = index_map[v1];
        std::cerr << "v0:" << v0 << " i0:" << idx0 << "\n";
        std::cerr << "v1:" << v1 << " i1:" << idx1 << "\n";
    }
#endif

    boost::remove_vertex(v0, g);
    {
        auto index_map = get(boost::vertex_index, g);
        auto idx0 = index_map[v0];
        auto idx1 = index_map[v1];

        std::cerr << "after removal of v0\n";
        std::cerr << "v0:" << v0 << " i0:" << idx0 << "\n";
        std::cerr << "v1:" << v1 << " i1:" << idx1 << "\n";
    }
}

