/** doc pdvd/115 -- charge-aware pricing of the Steiner BASE graph before the
    Voronoi step (WireCellClus/SteinerBaseWeight.h) and its knob defaults.

    What is pinned:
      * the C++ defaults (alpha 0, scope "tree") and the scope parser (a typo
        is refused, not mapped to a level);
      * the factor: alpha <= 0 => exactly 1; 1 + alpha * mean zero-plane count;
      * reweight_base_graph returns a COPY with the same vertices and edges,
        the base untouched; alpha 0 keeps every weight; alpha > 0 scales each
        edge by the factor of its own endpoints;
      * the routing effect at its public core (Graphs::Weighted::voronoi and
        the bridge cost create_enhanced_steiner_graph derives from it, which
        is src-private): on a 5-vertex graph with two terminals and two routes
        between them -- a short chord through a two-blank vertex and a longer
        path through three-plane vertices -- the historical graph makes the
        chord bridge cheaper and the priced graph makes the on-image path
        cheaper.  The full chain (topology, "tree" vs "tree+path" weights) is
        exercised by the knob-ON arm and its STGW / STGE dump (doc 115).
 */

#include "WireCellClus/SteinerBaseWeight.h"
#include "WireCellClus/Graphs.h"

#include "WireCellUtil/doctest.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"

#include <map>
#include <vector>

using namespace WireCell::Clus::Steiner;
using WireCell::Clus::Graphs::Weighted::graph_type;

namespace {
    double w_of(const graph_type& g, size_t a, size_t b)
    {
        auto e = boost::edge(a, b, g);
        REQUIRE(e.second);
        return boost::get(boost::edge_weight, g, e.first);
    }
}

TEST_CASE("steiner base weight: defaults, parser and factor")
{
    WireCell::PluginManager::instance().add("WireCellClus");
    auto icfg = WireCell::Factory::lookup<WireCell::IConfigurable>("CreateSteinerGraph", "doc115_base_weight_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE(cfg.isMember("base_weight_blank_alpha"));
    CHECK(cfg["base_weight_blank_alpha"].asDouble() == 0.0);
    REQUIRE(cfg.isMember("base_weight_scope"));
    CHECK(cfg["base_weight_scope"].asString() == "tree");

    BaseWeightScope s = BaseWeightScope::tree_path;
    CHECK(parse_base_weight_scope("tree", s));
    CHECK(s == BaseWeightScope::tree);
    CHECK(parse_base_weight_scope("tree+path", s));
    CHECK(s == BaseWeightScope::tree_path);
    CHECK_FALSE(parse_base_weight_scope("path", s));
    CHECK(s == BaseWeightScope::tree_path);
    CHECK_FALSE(parse_base_weight_scope("", s));

    CHECK(blank_weight_factor(0, 0, 1.0) == 1.0);
    CHECK(blank_weight_factor(3, 3, 0.0) == 1.0);
    CHECK(blank_weight_factor(3, 3, -1.0) == 1.0);
    CHECK(blank_weight_factor(1, 1, 1.0) == doctest::Approx(2.0));
    CHECK(blank_weight_factor(2, 2, 0.5) == doctest::Approx(2.0));
    CHECK(blank_weight_factor(0, 2, 1.0) == doctest::Approx(2.0));
    CHECK(blank_weight_factor(1, 0, 0.5) == doctest::Approx(1.25));
}

TEST_CASE("steiner base weight: reweight is a priced copy")
{
    // 0 -1- 1 -1- 2, 0 -3- 2 ; nz: 0->0, 1->2, 2->0
    graph_type g(3);
    boost::add_edge(0, 1, 1.0, g);
    boost::add_edge(1, 2, 1.0, g);
    boost::add_edge(0, 2, 3.0, g);
    std::map<size_t, int> nz = {{0, 0}, {1, 2}, {2, 0}};
    auto nzero = [&](size_t v) { return nz.at(v); };

    auto same = reweight_base_graph(g, nzero, 0.0);
    CHECK(boost::num_vertices(same) == 3);
    CHECK(boost::num_edges(same) == 3);
    CHECK(w_of(same, 0, 1) == 1.0);
    CHECK(w_of(same, 1, 2) == 1.0);
    CHECK(w_of(same, 0, 2) == 3.0);

    auto priced = reweight_base_graph(g, nzero, 1.0);
    CHECK(boost::num_vertices(priced) == 3);
    CHECK(boost::num_edges(priced) == 3);
    CHECK(w_of(priced, 0, 1) == doctest::Approx(2.0));   // 1 * (1 + 0.5*(0+2))
    CHECK(w_of(priced, 1, 2) == doctest::Approx(2.0));
    CHECK(w_of(priced, 0, 2) == doctest::Approx(3.0));   // both ends three-plane
    // the base is untouched
    CHECK(w_of(g, 0, 1) == 1.0);
    CHECK(w_of(g, 1, 2) == 1.0);
    CHECK(w_of(g, 0, 2) == 3.0);
}

TEST_CASE("steiner base weight: the priced graph moves the Voronoi bridge choice")
{
    // create_enhanced_steiner_graph is src-private (SteinerGrapher.h), so
    // the chain is pinned here at its public core, Graphs::Weighted::voronoi,
    // plus the bridge cost it feeds (vor.distance[s] + w + vor.distance[t]).
    // terminals 0 and 4.  chord: 0 -1- 1 -1- 4 through the blank vertex 1;
    // on-image path: 0 -1- 2 -1- 3 -1- 4 through three-plane vertices.
    using namespace WireCell::Clus::Graphs::Weighted;
    graph_type base(5);
    boost::add_edge(0, 1, 1.0, base);
    boost::add_edge(1, 4, 1.0, base);
    boost::add_edge(0, 2, 1.0, base);
    boost::add_edge(2, 3, 1.0, base);
    boost::add_edge(3, 4, 1.0, base);
    std::map<size_t, int> nz = {{0, 0}, {1, 2}, {2, 0}, {3, 0}, {4, 0}};
    auto nzero = [&](size_t v) { return nz.at(v); };
    std::vector<vertex_type> terminals = {0, 4};

    auto bridge_cost = [](const graph_type& g, const Voronoi& vor, size_t a, size_t b) {
        auto e = boost::edge(a, b, g);
        REQUIRE(e.second);
        return vor.distance[a] + boost::get(boost::edge_weight, g, e.first) + vor.distance[b];
    };

    SUBCASE("historical graph: the chord bridge is cheaper")
    {
        auto vor = voronoi(base, terminals);
        CHECK(vor.distance[1] == doctest::Approx(1.0));
        CHECK(vor.distance[2] == doctest::Approx(1.0));
        CHECK(vor.distance[3] == doctest::Approx(1.0));
        CHECK(vor.terminal[2] == 0);
        CHECK(vor.terminal[3] == 4);
        // vertex 1 is a tie; whichever region it fell in, its bridge costs 2 against the path's 3
        const size_t other = (vor.terminal[1] == 0) ? 4 : 0;
        CHECK(bridge_cost(base, vor, 1, other) == doctest::Approx(2.0));
        CHECK(bridge_cost(base, vor, 2, 3) == doctest::Approx(3.0));
    }
    SUBCASE("priced graph, alpha 1: the on-image path is cheaper and the base is unchanged")
    {
        auto routing = reweight_base_graph(base, nzero, 1.0);   // 0-1 and 1-4 become 2 each
        auto vor = voronoi(routing, terminals);
        CHECK(vor.distance[1] == doctest::Approx(2.0));
        CHECK(vor.distance[2] == doctest::Approx(1.0));
        CHECK(vor.distance[3] == doctest::Approx(1.0));
        const size_t other = (vor.terminal[1] == 0) ? 4 : 0;
        CHECK(bridge_cost(routing, vor, 1, other) == doctest::Approx(4.0));
        CHECK(bridge_cost(routing, vor, 2, 3) == doctest::Approx(3.0));
        // the geometric length the "tree" scope reads back for the reduced graph
        CHECK(w_of(base, 0, 2) == 1.0);
        CHECK(w_of(base, 0, 1) == 1.0);
    }
}
