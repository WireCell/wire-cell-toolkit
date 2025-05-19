#include "WireCellClus/Graph.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell::Clus::Graph;

static void check_graph_content(const Ident::graph_type& g)
{                           // vector access
    REQUIRE(boost::num_vertices(g) == 2);
    REQUIRE(boost::num_edges(g) == 1);

    auto idents = Ident::idents(g);
    REQUIRE(idents.size() == 2);
    REQUIRE(idents[0] == 1);
    REQUIRE(idents[1] == 2);
    auto weights = Ident::weights(g);
    REQUIRE(weights.size() == 1);
    REQUIRE(weights[0] == 3.0);
}
    

TEST_CASE("clus graph") {
    Ident::graph_type g;

    auto v0 = boost::add_vertex(Ident::VertexProp{1}, g);
    auto v1 = boost::add_vertex(Ident::VertexProp{2}, g);
    auto [e,eok] = boost::add_edge(v0, v1, Ident::EdgeProp{3.0}, g);
    REQUIRE(eok);

    REQUIRE(boost::source(e, g) == v0);
    REQUIRE(boost::target(e, g) == v1);

    check_graph_content(g);
        

    // Dataset persistence
    store_t store;
    {
        bool save_ok = Ident::save("test", g, store);
        REQUIRE(save_ok);
        REQUIRE(store.size() == 2);
        REQUIRE(store.find("graph_test_vertices") != store.end());
        REQUIRE(store.find("graph_test_edges") != store.end());

    }
    {
        Ident::graph_type g2;
        bool load_ok = Ident::load("test", store, g2);
        REQUIRE(load_ok);
        check_graph_content(g2);
    }
}
