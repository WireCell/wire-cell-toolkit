/**
   Test bolting on "facades" to PC tree nodes.

   Here, we take a common facade type across each layer of the tree.
 */

#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/doctest.h"

#include <variant>

namespace WireCell::Clustering {

    using namespace PointCloud;

    using Level = NaryTree::Facade<Tree::Points>;

    struct L2 : public Level {  
        int level{2};           // leaf
    };    

    struct L1 : public NaryTree::FacadeParent<L2, Tree::Points> {
        int level{1};           // interior
    };

    struct L0 : public NaryTree::FacadeParent<L1, Tree::Points> {
        int level{0};           // root
    };

    template<typename L>
    L& layer(Level* l)
    {
        L* ll = dynamic_cast<L*>(l);
        if (!ll) {
            raise<TypeError>("can not dynamic cast point tree node facade");
        }
        return *ll;
    }
    template<typename L>
    L& layer(Tree::Points& value)
    {
        auto* l = value.facade();
        return layer<L>(l);
    }
    template<typename L>
    L& layer(Tree::Points::node_t& node)
    {
        return layer<L>(node.value);
    }
    template<typename L>
    L& layer(Tree::Points::node_t* node)
    {
        return layer<L>(node->value);
    }
}

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::Clustering;

TEST_CASE("point tree subclass construct")
{
    Tree::Points::node_t n0;

    // n0.value.set_facade(std::make_unique<L0>());
    auto* l0 = n0.value.facade<L0>();
    REQUIRE(l0 != nullptr);
    REQUIRE(l0->level == 0);
    REQUIRE(l0->node() == &n0);

    auto& n1 = *n0.insert();
    // n1.value.set_facade(std::make_unique<L1>());
    auto* l1 = n1.value.facade<L1>();
    REQUIRE(l1 != nullptr);
    REQUIRE(l1->level == 1);
    REQUIRE(l1->node() == &n1);

    auto& n2 = *n1.insert();
    // n2.value.set_facade(std::make_unique<L2>());
    auto* l2 = n2.value.facade<L2>();
    REQUIRE(l2 != nullptr);
    REQUIRE(l2->level == 2);
    REQUIRE(l2->node() == &n2);

    REQUIRE(l0->nchildren() == 1);
    REQUIRE(l0->children()[0] == l1);


    Tree::Points::node_t nn0;
    nn0.value.set_facade(std::make_unique<L0>());
    auto* ll0 = nn0.value.facade<L0>();
    REQUIRE(ll0->nchildren() == 0);
    ll0->take_children(*l0);
    REQUIRE(ll0->nchildren() == 1);
    REQUIRE(l0->nchildren() == 0);
    REQUIRE(ll0->children()[0] == l1);
}
