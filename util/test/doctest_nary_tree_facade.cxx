
#include "WireCellUtil/NaryTreeFacade.h"
#include "WireCellUtil/NaryTesting.h"

#include "WireCellUtil/String.h"

#include "WireCellUtil/doctest.h"

#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Type.h"

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
    int removed{0};

    virtual bool on_remove(const std::vector<node_type*>& path) {
        if (path.size() == 2) {
            ++removed;
        }
        debug("L1 on_remove n={} facade={} node={} path.size={}",
              removed, (void*)this, (void*)node(), path.size());
        return this->NaryTree::FacadeParent<L2, Introspective>::on_remove(path);
    }
    int inserted{0};
    virtual bool on_insert(const std::vector<node_type*>& path) {
        if (path.size() == 2) {
            ++inserted;
        }
        debug("L1 on_insert n={} facade={} node={} path.size={}",
              inserted, (void*)this, (void*)node(), path.size());
        return this->NaryTree::FacadeParent<L2, Introspective>::on_insert(path);
    }

};

struct L0 : public NaryTree::FacadeParent<L1, Introspective> {
    int level{0};           // root
    int removed{0};

    virtual bool on_remove(const std::vector<node_type*>& path) {
        if (path.size() == 2) {
            ++removed;
        }
        debug("L0 on_remove n={} facade={} node={} path.size={}",
              removed, (void*)this, (void*)node(), path.size());
        return this->NaryTree::FacadeParent<L1, Introspective>::on_remove(path);
    }
    int inserted{0};
    virtual bool on_insert(const std::vector<node_type*>& path) {
        if (path.size() == 2) {
            ++inserted;
        }
        debug("L0 on_insert n={} facade={} node={} path.size={}",
              inserted, (void*)this, (void*)node(), path.size());
        return this->NaryTree::FacadeParent<L1, Introspective>::on_insert(path);
    }
};

static void dump_introspective_dfs(const Introspective& i)
{
    for (const auto& node : i.node->depth()) {
        debug(node.value.name);
    }
}


void check_level_generic(Level* fac, std::list<size_t> layer_sizes)
{
    REQUIRE(fac);
    REQUIRE(fac->node());
    debug("check_level({} '{}', below={})",
          type(fac), fac->node()->value.name, layer_sizes.size());
}

template<typename L>
void check_level(L* fac, std::list<size_t> layer_sizes)
{
    check_level_generic(fac, layer_sizes);
    
    const size_t nkids = layer_sizes.front();
    layer_sizes.pop_front();
    REQUIRE(nkids == fac->nchildren()); // this counts child nodes
    auto kids = fac->children();
    REQUIRE(nkids == kids.size()); // this counts child facades
    for (auto kid : kids) {
        check_level(kid, layer_sizes);
    }
}
template<>
void check_level<L2>(L2* fac, std::list<size_t> layer_sizes)
{
    check_level_generic(fac, layer_sizes);
}

TEST_CASE("nary tree facade construct") {
    std::list<size_t> layer_sizes ={2,3};
    auto root = make_layered_tree(layer_sizes);
    auto* l0 = root->value.facade<L0>();
    check_level(l0, layer_sizes);
        
    // We pre-make the tree so facades are not alive to see any inserted
    // notifications.
    REQUIRE(l0->inserted == 0);

}

TEST_CASE("nary tree facade separate retain") {
    std::list<size_t> layer_sizes ={2,10};
    auto root = make_layered_tree(layer_sizes);
    debug("initial tree");
    dump_introspective_dfs(root->value);
    auto* l0 = root->value.facade<L0>();
    check_level(l0, layer_sizes);

    auto* l0c0 = l0->children()[0];

    {
        auto l0children = l0->children();
        debug("before separate: root facade has {} inserted for total of {} or {}",
              l0->inserted, l0->nchildren(), l0children.size());
    }
    // Separate L2 children of first L1 child of L0.  This makes three new L1 in
    // addition to the original 2.
    debug("calling separate on {} '{}'", type(l0c0), l0c0->node()->value.name);
    //auto splits = l0c0->separate({0,0,-1,1,1,-1,2,2});
    auto splits = l0->separate(l0c0, {0,0,-1,1,1,-1,2,2}, false); // no removal
    REQUIRE(l0c0 != nullptr);
    REQUIRE(splits.size() == 3);

    REQUIRE(l0->nchildren() == 2+3);

    {
        auto l0children = l0->children();
        debug("after separate: root facade has {} inserted for total of {} or {}",
              l0->inserted, l0->nchildren(), l0children.size());
        REQUIRE(l0children.size() == 2+3);
        REQUIRE(l0->inserted == 3);
        REQUIRE(l0->removed == 0);
    }
}

TEST_CASE("nary tree facade separate remove") {
    std::list<size_t> layer_sizes ={2,10};
    auto root = make_layered_tree(layer_sizes);
    debug("initial tree");
    dump_introspective_dfs(root->value);
    auto* l0 = root->value.facade<L0>();
    check_level(l0, layer_sizes);

    auto* l0c0 = l0->children()[0];

    {
        auto l0children = l0->children();
        debug("before separate: root facade has {} inserted for total of {} or {}",
              l0->inserted, l0->nchildren(), l0children.size());
    }
    // Separate L2 children of first L1 child of L0.  This makes three new L1 in
    // addition to the original 2.
    debug("calling separate on {} '{}'", type(l0c0), l0c0->node()->value.name);
    //auto splits = l0c0->separate({0,0,-1,1,1,-1,2,2});
    auto splits = l0->separate(l0c0, {0,0,-1,1,1,-1,2,2}, true); // with removal
    REQUIRE(l0c0 == nullptr);
    REQUIRE(splits.size() == 3);

    REQUIRE(l0->nchildren() == 2+3-1);

    {
        auto l0children = l0->children();
        debug("after separate: root facade has {} inserted for total of {} or {}",
              l0->inserted, l0->nchildren(), l0children.size());
        REQUIRE(l0children.size() == 2+3-1);
        REQUIRE(l0->inserted == 3);
        REQUIRE(l0->removed == 1);
    }
}

TEST_CASE("nary tree facade remove") {

    std::list<size_t> layer_sizes ={2,3};
    auto root = make_layered_tree(layer_sizes);
    debug("initial tree");
    dump_introspective_dfs(root->value);
    auto* l0 = root->value.facade<L0>();
    check_level(l0, layer_sizes);

    // The facade was not constructed prior to any "inserted" notifications.
    REQUIRE(l0->inserted == 0);

    auto l0children_before = l0->children();

    auto dead_node = l0->remove_child(*l0children_before[0]);
    REQUIRE(dead_node);
    dead_node = nullptr;

    auto l0children_after = l0->children();

    debug("root facade has {} removed", l0->removed);
    debug("root facade children.size()={} nchildren={}", l0children_after.size(), l0->nchildren());

    // The L0 facade was in place for the remove so should have gotten the
    // "removing" notification.
    REQUIRE(l0->removed == 1);

    // We've removed one child node, so only one child node should remain.
    REQUIRE(l0->nchildren() == 1);

    // We've removed one child node, so only one child facade should remain.
    REQUIRE(l0children_after.size() == 1);

}
