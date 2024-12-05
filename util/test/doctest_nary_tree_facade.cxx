
#include "WireCellUtil/NaryTreeFacade.h"
#include "WireCellUtil/NaryTesting.h"

#include "WireCellUtil/String.h"

#include "WireCellUtil/doctest.h"

#include "WireCellUtil/Logging.h"

#include <string>

using namespace WireCell;
using namespace WireCell::NaryTesting;
using spdlog::debug;

using Level = NaryTree::Facade<Introspective>;

struct L2 : public Level {  
    int level{2};           // leaf
};    

struct L1 : public NaryTree::FacadeParent<L2, Introspective> {
    int level{1};           // interior
};

struct L0 : public NaryTree::FacadeParent<L1, Introspective> {
    int level{0};           // root
};

static void dump_introspective_dfs(const Introspective& i)
{
    for (const auto& node : i.node->depth()) {
        debug(node.value.name);
    }
}


TEST_CASE("nary tree facade separate") {
    std::list<size_t> layer_sizes ={2,10};
    auto root = make_layered_tree(layer_sizes);
    debug("initial tree");
    dump_introspective_dfs(root->value);
    auto* l0 = root->value.facade<L0>();
    REQUIRE(l0 != nullptr);
    REQUIRE(l0->level == 0);
    REQUIRE(l0->node() != nullptr);
    REQUIRE(l0->nchildren() == 2);
    auto* l0c0 = l0->children()[0];
    REQUIRE(l0c0 != nullptr);

    // Separate children, skipping some and ignoring a remainder.
    auto splits = l0c0->separate({0,0,-1,1,1,-1,2,2});
    REQUIRE(splits.size() == 3);
    REQUIRE(l0->nchildren() == 2+3);
    REQUIRE(l0c0->nchildren() == 4);
    debug("after separate");
    dump_introspective_dfs(root->value);
}
