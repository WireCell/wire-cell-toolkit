/**

   Excercise Boost Graph Library's "filtered_graph".

   We do this in two ways:

   1. make a "predicate" struct

   2. provide predicate as lambda function

   The filtered graph will preserve the iterators (vertex descriptors)
   of the original graph.  The iterators that are visible via the
   filtered graph are sparse compared to the original set.

   By applying boost::copy_graph on the filtered graph we make a
   proper subgraph.  The filtered graph and its copy has the same
   vertices but their descriptors differ in general.

 */

// https://www.boost.org/doc/libs/1_79_0/libs/graph/doc/filtered_graph.html



#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/copy.hpp>

#include <functional>
#include <iostream>

// Simple graph with a vertex property.
struct node_t {
    int num;
};
using graph_t = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, node_t>;
using vertex_t = boost::graph_traits<graph_t>::vertex_descriptor;
using edge_t = boost::graph_traits<graph_t>::edge_descriptor;

// This predicate struct can work for both edge and vertex filtering.
// It retains all edges and vertices with a "num" larger than 1.
struct Predicate { 
    bool operator()(graph_t::edge_descriptor ed) const {
        return true;
    }
    bool operator()(graph_t::vertex_descriptor vd) const {
        return (*g)[vd].num > 1;
    }

    graph_t* g;
};
using filtered_t = boost::filtered_graph<graph_t, Predicate, Predicate>;

// A more free form filtered where we give a lambda
using Filtered = boost::filtered_graph<graph_t, boost::keep_all, std::function<bool(vertex_t)> >;

template<typename Graph>
void dump(const Graph& g, const std::string& msg="")
{
    std::cout << msg << std::endl;
    boost::write_graphviz(std::cout, g, boost::make_label_writer(boost::get(&node_t::num, g)));
}

int main()
{
    graph_t g;
    auto v1 = boost::add_vertex({1}, g);
    auto v2 = boost::add_vertex({2}, g);
    auto v3 = boost::add_vertex({3}, g);
    add_edge(v1, v2, g);
    add_edge(v2, v3, g);
    add_edge(v3, v1, g);
    dump(g, "g");

    Predicate p{&g};
    filtered_t fg(g, p, p);
    dump(fg, "fg");

    graph_t g2;
    boost::copy_graph(fg, g2);
    dump(g2, "g2");


    std::unordered_set dead = {v2};
    Filtered fg2(g, {}, [&](vertex_t vtx) {
        return dead.count(g[vtx].num) == 0; });
    graph_t g3;
    boost::copy_graph(fg2, g3);
    dump(g3, "g3");

    std::cerr << "g2: #v:" << boost::num_vertices(g2) << " #e:" << boost::num_edges(g2) << std::endl;
    std::cerr << "g3: #v:" << boost::num_vertices(g3) << " #e:" << boost::num_edges(g3) << std::endl;
    boost::copy_graph(g2,g3);
    std::cerr << "copy_graph(g2,g3):\n";
    std::cerr << "g2: #v:" << boost::num_vertices(g2) << " #e:" << boost::num_edges(g2) << std::endl;
    std::cerr << "g3: #v:" << boost::num_vertices(g3) << " #e:" << boost::num_edges(g3) << std::endl;
    dump(g2, "g2");
    dump(g3, "g3");

    return 0;
}
