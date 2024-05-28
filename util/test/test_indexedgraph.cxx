#include "WireCellUtil/IndexedGraph.h"
#include "WireCellUtil/Testing.h"

#include "WireCellUtil/Graph.h"


using namespace WireCell;

char foo(size_t ind) { return "if"[ind]; }

typedef std::shared_ptr<int> iptr_t;
typedef std::shared_ptr<float> fptr_t;
typedef std::variant<iptr_t, fptr_t> if_t;

struct if_node_t {
    if_t ptr;
    if_node_t()
      : ptr()
    {
    }
    if_node_t(const if_t& ift)
      : ptr(ift)
    {
    }
    if_node_t(const iptr_t& i)
      : ptr(i)
    {
    }
    if_node_t(const fptr_t& f)
      : ptr(f)
    {
    }

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

int main()
{
    std::unordered_map<if_t, int> p2i;

    typedef IndexedGraph<if_node_t> indexed_graph_t;
    indexed_graph_t g;
    if_t one = std::make_shared<int>(42);
    if_t two = std::make_shared<float>(6.9);
    if_t tre = std::make_shared<float>(33);

    p2i[one] = 1;
    p2i[two] = 2;
    p2i[tre] = 3;

    // test variant
    Assert(1 == tre.index());
    Assert(*std::get<1>(two) < *std::get<1>(tre));
    Assert(one < two);  // caution: compares against indices first!
    Assert(*std::get<0>(one) == 42);
    Assert(nullptr == std::get_if<1>(&one));
    Assert('i' == foo(0));
    Assert('f' == foo(1));

    if_t oneprime = one;
    Assert(one == oneprime);

    auto v1 = g.vertex(one);
    auto v2 = g.vertex(one);
    Assert(v1 == v2);
    g.edge(one, two);
    g.edge(one, two);
    g.edge(two, one);
    g.edge(one, tre);
    auto verts = boost::vertices(g.graph());
    Assert(verts.second - verts.first == 3);

    Assert(g.has(oneprime));

    indexed_graph_t g2(g.graph());

    for (auto n : g2.neighbors(two)) {
        Assert(g2.has(n));
    }

    {
        std::vector<if_node_t> bigf;
        g.neighbors(back_inserter(bigf), one,
                    [](const if_node_t& vp){
                        if (std::holds_alternative<fptr_t>(vp.ptr)) {
                            auto f = std::get<fptr_t>(vp.ptr);
                            return *f>10.0;
                        }
                        return false;
                    });
        Assert (bigf.size() == 1);
        float val = *std::get<fptr_t>(bigf[0].ptr);
        Assert (val == 33);
    }

    
    boost::default_writer w;
    std::vector<std::string> names{"one", "two", "tre"};
    {
        std::unordered_map<indexed_graph_t::vdesc_t, std::string> ids;
        for (auto u : boost::make_iterator_range(vertices(g2.graph()))) {
            ids[u] = names[ids.size()];
        }
        boost::write_graphviz(std::cout, g2.graph(), w, w, w, boost::make_assoc_property_map(ids));
    }
    {
        std::vector<if_node_t> ss{one, two};
        indexed_graph_t sg = g.induced_subgraph(ss.begin(), ss.end());
        std::unordered_map<indexed_graph_t::vdesc_t, std::string> ids;
        for (auto u : boost::make_iterator_range(vertices(sg.graph()))) {
            ids[u] = names[ids.size()];
        }
        boost::write_graphviz(std::cerr, sg.graph(), w, w, w, boost::make_assoc_property_map(ids));    
    }

    return 0;
}
