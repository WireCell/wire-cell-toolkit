#include "WireCellUtil/IndexedGraph.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include "WireCellUtil/Graph.h"

#include <sstream>

using namespace WireCell;

char foo(size_t ind) { return "if"[ind]; }

typedef std::shared_ptr<int> iptr_t;
typedef std::shared_ptr<float> fptr_t;
typedef std::variant<iptr_t, fptr_t> if_t;

struct if_node_t {
    if_t ptr;
    if_node_t() : ptr() {}
    if_node_t(const if_t& ift) : ptr(ift) {}
    if_node_t(const iptr_t& i) : ptr(i) {}
    if_node_t(const fptr_t& f) : ptr(f) {}

    bool operator==(const if_node_t& other) const { return ptr == other.ptr; }
};

namespace std {
    template <>
    struct hash<if_node_t> {
        std::size_t operator()(const if_node_t& n) const
        {
            if (std::holds_alternative<iptr_t>(n.ptr)) {
                return (size_t)(std::get<iptr_t>(n.ptr).get());
            }
            if (std::holds_alternative<fptr_t>(n.ptr)) {
                return (size_t)(std::get<fptr_t>(n.ptr).get());
            }
            return 0;
        }
    };
}  // namespace std

TEST_SUITE("indexedgraph") {

TEST_CASE("variant basics") {
    if_t one = std::make_shared<int>(42);
    if_t two = std::make_shared<float>(6.9);
    if_t tre = std::make_shared<float>(33);

    CHECK(1 == tre.index());
    CHECK(*std::get<1>(two) < *std::get<1>(tre));
    CHECK(one < two);  // caution: compares against indices first!
    CHECK(*std::get<0>(one) == 42);
    CHECK(nullptr == std::get_if<1>(&one));
    CHECK('i' == foo(0));
    CHECK('f' == foo(1));

    if_t oneprime = one;
    CHECK(one == oneprime);
}

TEST_CASE("indexed graph operations") {
    typedef IndexedGraph<if_node_t> indexed_graph_t;
    indexed_graph_t g;
    if_t one = std::make_shared<int>(42);
    if_t two = std::make_shared<float>(6.9);
    if_t tre = std::make_shared<float>(33);

    auto v1 = g.vertex(one);
    auto v2 = g.vertex(one);
    CHECK(v1 == v2);

    g.edge(one, two);
    g.edge(one, two);
    g.edge(two, one);
    g.edge(one, tre);

    auto verts = boost::vertices(g.graph());
    CHECK(verts.second - verts.first == 3);

    if_t oneprime = one;
    CHECK(g.has(oneprime));

    indexed_graph_t g2(g.graph());

    for (auto n : g2.neighbors(two)) {
        CHECK(g2.has(n));
    }

    {
        std::vector<if_node_t> bigf;
        g.neighbors(back_inserter(bigf), one,
                    [](const if_node_t& vp) {
                        if (std::holds_alternative<fptr_t>(vp.ptr)) {
                            auto f = std::get<fptr_t>(vp.ptr);
                            return *f > 10.0;
                        }
                        return false;
                    });
        CHECK(bigf.size() == 1);
        float val = *std::get<fptr_t>(bigf[0].ptr);
        CHECK(val == 33);
    }

    boost::default_writer w;
    std::vector<std::string> names{"one", "two", "tre"};
    {
        std::unordered_map<indexed_graph_t::vdesc_t, std::string> ids;
        for (auto u : boost::make_iterator_range(vertices(g2.graph()))) {
            ids[u] = names[ids.size()];
        }
        std::ostringstream oss;
        boost::write_graphviz(oss, g2.graph(), w, w, w, boost::make_assoc_property_map(ids));
        spdlog::debug(oss.str());
    }
    {
        std::vector<if_node_t> ss{one, two};
        indexed_graph_t sg = g.induced_subgraph(ss.begin(), ss.end());
        std::unordered_map<indexed_graph_t::vdesc_t, std::string> ids;
        for (auto u : boost::make_iterator_range(vertices(sg.graph()))) {
            ids[u] = names[ids.size()];
        }
        std::ostringstream oss;
        boost::write_graphviz(oss, sg.graph(), w, w, w, boost::make_assoc_property_map(ids));
        spdlog::debug(oss.str());
    }
}

}  // TEST_SUITE("indexedgraph")
