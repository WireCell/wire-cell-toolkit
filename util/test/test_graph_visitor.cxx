
#include "WireCellUtil/Graph.h"




#include <exception>
#include <iostream>
#include <fstream>


// in WireCellUtil/GraphTools.h
template<typename Gr>
std::string dotify(const Gr& gr)
{
    using vertex_t = typename boost::graph_traits<Gr>::vertex_descriptor;
    std::stringstream ss;
    boost::write_graphviz(ss, gr, [&](std::ostream& out, vertex_t v) {
        const auto& dat = gr[v];
        out << "[label=\"" << dat.code << dat.id << "\"]";
    });
    return ss.str() + "\n";
}



struct Node {
    char code;
    int id;
};

struct Edge {
    size_t count;
};


using ugraph_t = boost::adjacency_list<boost::vecS, boost::vecS,
                                      boost::undirectedS,
                                      Node, Edge>;
using uvdesc_t = boost::graph_traits<ugraph_t>::vertex_descriptor;
using uedesc_t = boost::graph_traits<ugraph_t>::edge_descriptor;

using dgraph_t = boost::adjacency_list<boost::vecS, boost::vecS,
                                      boost::directedS,
                                      Node, Edge>;
using dvdesc_t = boost::graph_traits<dgraph_t>::vertex_descriptor;
using dedesc_t = boost::graph_traits<dgraph_t>::edge_descriptor;

/// Thrown when an iteration is over (eg, BFS/DFS graph visit)
class StopIteration : virtual public std::domain_error
{
    using std::domain_error::domain_error;
};



template<typename Gr>
struct DFSV : public boost::default_dfs_visitor
{
    using graph_t = Gr;
    using vertex_t = typename boost::graph_traits<Gr>::vertex_descriptor;
    using edge_t = typename boost::graph_traits<Gr>::edge_descriptor;

    void examine_vertex(vertex_t v, const graph_t& g) const {
        const auto& dat = g[v];
        std::cerr << dat.code << dat.id << std::endl;
        // if (dat.code == 'c') {
        //     throw StopIteration("channel");
        // }
    }
    void examine_edge(edge_t e, const graph_t& g) const {
        auto t = boost::source(e,g);
        auto h = boost::target(e,g);
        std::cerr << "e: " << g[t].code << g[t].id << " -> " << g[h].code << g[h].id << std::endl;
    }
    void tree_edge(edge_t e, const graph_t& g) const {
        auto t = boost::source(e,g);
        auto h = boost::target(e,g);
        std::cerr << "t: " << g[t].code << g[t].id << " -> " << g[h].code << g[h].id << std::endl;
    }
    void finish_vertex(vertex_t v, const graph_t& g) const {
        const auto& dat = g[v];
        std::cerr << dat.code << dat.id << std::endl;
    }    
};

template <typename Gr>
struct BFSV : public boost::default_bfs_visitor
{
    using graph_t = Gr;
    using vertex_t = typename boost::graph_traits<Gr>::vertex_descriptor;
    using edge_t = typename boost::graph_traits<Gr>::edge_descriptor;

    void examine_vertex(vertex_t v, const graph_t& g) const {
        const auto& dat = g[v];
        std::cerr << dat.code << dat.id << std::endl;
    }
    void examine_edge(edge_t e, const graph_t& g) const {
        auto t = boost::source(e,g);
        auto h = boost::target(e,g);
        std::cerr << "e: " << g[t].code << g[t].id << " -> " << g[h].code << g[h].id << std::endl;
    }
    void tree_edge(edge_t e, const graph_t& g) const {
        auto t = boost::source(e,g);
        auto h = boost::target(e,g);
        std::cerr << "t: " << g[t].code << g[t].id << " -> " << g[h].code << g[h].id << std::endl;
    }
    void finish_vertex(vertex_t v, const graph_t& g) const {
        const auto& dat = g[v];
        std::cerr << dat.code << dat.id << std::endl;
    }    
};

template<typename Gr>
void fill(Gr& gr)
{
    auto s0 = boost::add_vertex({'s', 0}, gr);
    auto s1 = boost::add_vertex({'s', 1}, gr);
    auto b0 = boost::add_vertex({'b', 0}, gr);
    auto w0 = boost::add_vertex({'w', 0}, gr);
    auto b1 = boost::add_vertex({'b', 1}, gr);
    auto b2 = boost::add_vertex({'b', 2}, gr);
    auto w1 = boost::add_vertex({'w', 1}, gr);
    auto w2 = boost::add_vertex({'w', 2}, gr);
    auto w3 = boost::add_vertex({'w', 3}, gr);
    // auto c0 = boost::add_vertex({'c', 0}, gr);
    // auto c1 = boost::add_vertex({'c', 1}, gr);
    // auto c3 = boost::add_vertex({'c', 3}, gr);

    boost::add_edge(s0, b0, gr);
    boost::add_edge(b0, w0, gr);
    boost::add_edge(b0, w2, gr);
    boost::add_edge(s0, b1, gr);
    boost::add_edge(b1, w1, gr);
    boost::add_edge(b0, w3, gr);
    boost::add_edge(b1, w3, gr);
    boost::add_edge(s1, b2, gr);
    boost::add_edge(b2, w0, gr);
    boost::add_edge(b2, w1, gr);
    
    // boost::add_edge(w0, c0, gr);
    // boost::add_edge(w2, c1, gr);
    // boost::add_edge(w1, c0, gr);
    // boost::add_edge(w3, c3, gr);
}

template<typename Gr>
void test(std::string dotfile)
{
    Gr gr;
    fill(gr);

    std::ofstream outf(dotfile);
    outf << dotify(gr) << std::endl;
    std::cerr << dotfile << std::endl;

    DFSV<Gr> dfsv;
    std::vector<boost::default_color_type> colors(boost::num_vertices(gr));
    boost::iterator_property_map color_map(colors.begin(), boost::get(boost::vertex_index, gr));
    std::cerr << dotfile << ": DFSV 0\n";
    boost::depth_first_search(gr, dfsv, color_map, 0);
    // std::cerr << name << " DFSV 1\n";
    // boost::depth_first_search(gr, 1, dfsv);

    BFSV<Gr> bfsv;

    std::cerr << dotfile << ": BFSV 0\n";
    boost::breadth_first_search(gr, 0, boost::visitor(bfsv));
    std::cerr << dotfile << ": BFSV 0\n";
    boost::breadth_first_search(gr, 1, boost::visitor(bfsv));
}
int main (int argc, char* argv[])
{
    std::string arg0 = argv[0];
    test<ugraph_t>(arg0 + "_ugraph.dot");
    test<dgraph_t>(arg0 + "_dgraph.dot");
    return 0;
}
