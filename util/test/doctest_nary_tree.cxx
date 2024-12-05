#include "WireCellUtil/NaryTree.h"
#include "WireCellUtil/NaryTesting.h"

#include "WireCellUtil/String.h"

#include "WireCellUtil/doctest.h"

#include "WireCellUtil/Logging.h"

#include <string>

using namespace WireCell::NaryTree;
using namespace WireCell::NaryTesting;
using spdlog::debug;

TEST_CASE("nary tree depth iter") {
    using node_type = Node<int>;

    depth_iter<node_type> a, b;
    depth_range<node_type> r;

    CHECK( !a.node );

    CHECK( a == b );
    
    CHECK( r.begin() == r.end() );

    depth_iter<node_type const> c(a);
    CHECK( a == c );
}





// return numerber of nodes in a tree with all nodes in a given having
// the same number of children as given by layer_sizes.  The root node
// is excluded from layer_sizes but added to the returned count.
size_t nodes_in_uniform_tree(const std::list<size_t>& layer_sizes)
{
    size_t sum=0, prod=1;

    for (auto s : layer_sizes) {
        prod *= s;
        sum += prod;
    }
    return sum + 1;
}



TEST_CASE("nary tree node construction move") 
{
    Introspective d("moved");
    Introspective::node_type n(std::move(d));
    CHECK(n.value.nctor == 1);
    CHECK(n.value.nmove == 1);
    CHECK(n.value.ncopy == 0);
    CHECK(n.value.ndef == 0);
    size_t nkids = n.ndescendants();
    CHECK(nkids == 0);
}
TEST_CASE("nary tree node construction copy") {


    Introspective d("copied");
    Introspective::node_type n(d);
    CHECK(n.value.nctor == 1);
    CHECK(n.value.nmove == 0);
    CHECK(n.value.ncopy == 1);
    CHECK(n.value.ndef == 0);
    size_t nkids = n.ndescendants();
    CHECK(nkids == 0);
}


bool compare_introspective(const Introspective::owned_ptr& a,
                           const Introspective::owned_ptr& b)
{
    return b->value.name < a->value.name;
}

void dump_introspective(const Introspective& i, const std::string& prefix="", const std::string title="")
{
    // fixme: spdlog::debug is not working so we use std io for now.

    if (title.size()) {
        debug(title);
    }

    debug("{} {}", prefix, i);
    for (const auto& [k,v] : i.nactions) {
        debug("{}\t {} {}", prefix, k, v);
    }
    for ( const auto* cnode : i.node->children() ) {
        dump_introspective(cnode->value, prefix + "\t");
    }
}

static void dump_introspective_dfs(const Introspective& i)
{
    for (const auto& node : i.node->depth()) {
        debug(node.value.name);
    }
}

TEST_CASE("nary tree node ordered") 
{
    std::list<size_t> layer_sizes ={10,};
    auto root = make_layered_tree(layer_sizes);
    dump_introspective(root->value);
    root->sort_children(compare_introspective);
    dump_introspective(root->value);
    CHECK(root->value.nactions["ordered"] == 1);
}


TEST_CASE("nary tree node deeply ordered") 
{
    // Consider a tree with 4 layers:
    // 0 - root node with 2 children
    // 1 - 2 nodes each with 2 children
    // 2 - 4 nodes each with 2 children
    // 3 - 8 leaf nodes with no children
    std::list<size_t> layer_sizes ={2,2,2};
    auto root = make_layered_tree(layer_sizes);
    auto root2 = make_layered_tree(layer_sizes);
    // This tree is WELL ORDERED in order of insertion.
    debug("original tree");
    dump_introspective_dfs(root->value);

    // For some reason we now want to apply a new ordering (which will be a
    // reversal) to the children of each node in layer 1.  We do this in two
    // steps. First, select all the nodes in the desired layer via DFS visit.
    // Second, sort their children.

    debug("two step ordering");
    std::vector<Introspective::node_type*> to_order;

    // We can limit DFS to layer 1 as we just need the parents.  Note, depth()
    // layer 1 refers to just the root so we need 1+1
    for (auto& node : root->depth(1+1)) {
        if (node.nparents() != 1) {
            continue;
        }
        // two parents mean the node is at layer 3
        to_order.push_back(&node);
    }
    
    for (auto* node : to_order) {
        node->sort_children(compare_introspective);
        CHECK(node->value.nactions["ordered"] == 1);
    }
    dump_introspective_dfs(root->value);

    CHECK(root->value.nactions["ordered"] == 2);

    debug("one step ordering");
    // Now try again but sort in place and see if the DFS freaks out.
    for (auto& node : root2->depth(1+1)) {
        if (node.nparents() != 1) {
            continue;
        }
        node.sort_children(compare_introspective);
    }

    dump_introspective_dfs(root2->value);
    CHECK(root2->value.nactions["ordered"] == 2);
}

TEST_CASE("nary tree separate") 
{
    // Make a root with 10 children each with 2 and separate the 10.
    std::list<size_t> layer_sizes = {10,2};
    auto root = make_layered_tree(layer_sizes);
    // Separate, explicitly skipping some and underspecifying.
    auto nurseries = root->separate({0,0,-1,1,1,-1,2,2,-1});
    REQUIRE(nurseries.size() == 3);
    REQUIRE(root->nchildren() == 4);
    REQUIRE(nurseries[0].size() == 2);
    REQUIRE(nurseries[1].size() == 2);
    REQUIRE(nurseries[2].size() == 2);
}

TEST_CASE("nary tree simple tree tests") {

    const size_t live_count = Introspective::live_count;
    CHECK(live_count == 0);

    SUBCASE("scope for tree live")
    {   // scope for tree life.
        auto root = make_simple_tree();

        SUBCASE("children make sense") {

            auto childs = root->children();
            auto chit = childs.begin();

            {
                CHECK ( (*chit)->parent == root.get() );

                auto& d = (*chit)->value;

                CHECK( d.name == "0.0" );
                CHECK( d.nctor == 1);
                CHECK( d.nmove == 0);
                CHECK( d.ncopy == 1);
                CHECK( d.ndef == 0);

                ++chit;
            }
            {
                CHECK ( (*chit)->parent == root.get() );

                auto& d = (*chit)->value;

                CHECK( d.name == "0.1" );
                CHECK( d.nctor == 1);
                CHECK( d.nmove == 1);
                CHECK( d.ncopy == 0);
                CHECK( d.ndef == 0);

                ++chit;
            }
            {
                CHECK ( (*chit)->parent == root.get() );

                auto& d = (*chit)->value;

                CHECK( d.name == "0.2" );
                CHECK( d.nctor == 1);
                CHECK( d.nmove == 1);
                CHECK( d.ncopy == 0);
                CHECK( d.ndef == 0);

                ++chit;
            }

            CHECK(chit == childs.end());

            CHECK( childs.front()->next()->value.name == "0.1" );
            CHECK( childs.back()->prev()->value.name == "0.1" );
        }


        SUBCASE("Iterate the depth firs search") 
        {
            size_t nnodes = 0;
            std::vector<Introspective> data;
            const size_t level = 0;   // 0=unlimited, the default
            for (const auto& node : root->depth(level)) 
            {
                auto& val = node.value;
                debug("depth {} {}", nnodes, val);
                ++nnodes;
                data.push_back(val);
            }
            CHECK( nnodes == 4 );

            CHECK( data[0].name == "0" );
            CHECK( data[1].name == "0.0" );
            CHECK( data[2].name == "0.1" );
            CHECK( data[3].name == "0.2" );
        }

        {
            // Const version
            const auto* rc = root.get();
            size_t nnodes = 0;
            std::vector<Introspective> data;
            for (const auto& node : rc->depth()) 
            {
                ++nnodes;
                data.push_back(node.value);
            }
            CHECK( nnodes == 4 );
        }    

    } // tree r dies

    CHECK(live_count == Introspective::live_count);
}


TEST_CASE("nary tree bigger tree lifetime")
{
    size_t live_count = Introspective::live_count;
    CHECK(live_count == 0);


    {
        std::list<size_t> layer_sizes ={2,4,8};
        size_t niut = nodes_in_uniform_tree(layer_sizes);
        auto root = make_layered_tree(layer_sizes);
        debug("starting live count: {}, ending live count: {}, diff should be {}",
              live_count, Introspective::live_count, niut);
        CHECK(Introspective::live_count - live_count == niut);
    } // root and children all die

    CHECK(live_count == Introspective::live_count);

}

TEST_CASE("nary tree remove node")
{
    size_t live_count = Introspective::live_count;
    CHECK(live_count == 0);

    {
        std::list<size_t> layer_sizes ={2,4,8};
        // size_t niut = nodes_in_uniform_tree(layer_sizes);
        auto root = make_layered_tree(layer_sizes);

        auto children = root->children();

        Introspective::node_type* doomed = children.front();
        CHECK(doomed);

        {                       // test find by itself.
            Introspective::node_type::sibling_iter sit = root->find(doomed);
            // iterator of list of unique_ptr
            CHECK(sit->get() == doomed);
        }

        const size_t nbefore = root->nchildren();
        auto dead = root->remove(doomed);
        CHECK(dead.get() == doomed);
        CHECK( ! dead->parent );
        const size_t nafter = root->nchildren();
        CHECK(nbefore == nafter + 1);

    } // root and children all die

    CHECK(live_count == Introspective::live_count);

    { // insert node that is in another tree
        auto r1 = make_simple_tree("0");
        auto r2 = make_simple_tree("1");
        Introspective::node_type* traitor = r2->children().front();
        r1->insert(traitor);
        CHECK(r1->nchildren() == 4);
        CHECK(r2->nchildren() == 2);
        CHECK(r1->children().back()->value.name == "1.0");
        CHECK(r2->children().front()->value.name == "1.1");
    }

    CHECK(live_count == Introspective::live_count);
}

TEST_CASE("nary tree child iter")
{
    auto root = make_simple_tree();

    std::vector<std::string> nstack = {"0.2","0.1","0.0"};
    for (const auto& c : root->child_values()) {
        debug("child value: {}", c);
        CHECK(c.name == nstack.back());
        nstack.pop_back();
    }

    // const iterator
    const auto* rc = root.get();
    for (const auto& cc : rc->child_values()) {
        debug("child value: {} (const)", cc);
        // cc.name += " extra";    // should not compile
    }

    for (auto& n : root->child_nodes()) {
        n.value.name += " extra";
        debug("change child node value: {}", n.value);
    }
    for (const auto& cn : rc->child_nodes()) {
        debug("child node value: {}", cn.value);
        // cn.value.name += " extra"; // should not compile
    }
}

TEST_CASE("nary tree notify")
{
    // Note, Introspective counts the actions tested here only when the
    // notification path is size one. See the ordered test for that notification
    // which counts more broadly.

    {
        Introspective::node_type node;
        auto& nactions = node.value.nactions;
        CHECK(nactions.size() == 1);
        CHECK(nactions["constructed"] == 1);
    }


    debug("nary tree notify make simple tree");
    auto root = make_simple_tree();
    {
        auto& nactions = root->value.nactions;
        CHECK(nactions.find("constructed") != nactions.end());
        CHECK(nactions["constructed"] == 1);
    }
    {
        Introspective::node_type* doomed = root->children().front();
        auto& nactions = doomed->value.nactions;
        debug("doomed: {}", doomed->value.name);
        CHECK(nactions.find("constructed") != nactions.end());
        CHECK(nactions.find("inserted") != nactions.end());
        CHECK(nactions["constructed"] == 1);
        CHECK(nactions["inserted"] == 1);
        auto dead = root->remove(doomed);
        CHECK(nactions.size() == 3);
        CHECK(nactions["removing"] == 1);        
    }
    
}
TEST_CASE("nary tree flatten")
{
    Introspective::node_type root(Introspective("r"));
    auto& rval = root.value;
    REQUIRE(rval.node == &root);
    root.insert(make_simple_tree("r.0"));
    root.insert(make_simple_tree("r.1"));
    debug("descend root node \"{}\"", rval.name);
    for (auto& node : rval.node->depth()) {
        auto& value = node.value;
        debug("\tpath to parents from node \"{}\":", value.name);
        std::string path="", slash="";
        for (auto n : node.sibling_path()) {
            path += slash + std::to_string(n);
            slash = "/";
        }
        debug("\t\tpath: {}  \t (\"{}\")", path, value.name);
    }
    size_t nkids = root.ndescendants();
    CHECK(nkids == 8);

}

static size_t count_delim(const std::string& s, const std::string& d = ".")
{
    if (s.empty()) return 0;
    return WireCell::String::split(s, d).size();
}

TEST_CASE("nary tree depth limits")
{
    // const std::list<size_t> layer_sizes = {2,4,8};
    auto root = make_layered_tree({2,4,8});

    const std::vector<size_t> want = {1,
        nodes_in_uniform_tree({2}),
        nodes_in_uniform_tree({2,4}),
        nodes_in_uniform_tree({2,4,8})};

    const size_t nwants = want.size();
    const size_t max_wants = want.back();

    for (size_t ind=0; ind<nwants; ++ind) {
        const size_t level = ind+1;
        size_t count = 0;
        for (const auto& node : root->depth(level)) {
            debug("[{}/{}] #{}/{}/{} {}", ind, nwants, count, want[ind], max_wants, node.value.name);
            ++count;
            size_t nlev = count_delim(node.value.name);
            // make sure never descend beyond level
            CHECK(nlev <= level); 
        }
        CHECK(count == want[ind]);
    }
    for (const auto* child : root->children()[0]->children()[0]->children()) {
        debug("gggc: {}", child->value.name);
    }

}

void remove_children(bool notify_child)
{
    auto root = make_layered_tree({2,4,8});
    CHECK(root->nchildren() == 2);
    auto children = root->remove_children(notify_child);
    CHECK(root->nchildren() == 0);
    CHECK(children.size() == 2);
}
TEST_CASE("nary tree remove children no notify")
{
    remove_children(false);
}
TEST_CASE("nary tree remove children with notify")
{
    remove_children(true);
}

TEST_CASE("nary tree remove adopt children")
{
    auto root = make_layered_tree({2,4,8});
    CHECK(root->nchildren() == 2);
    auto children = root->remove_children();
    CHECK(root->nchildren() == 0);
    CHECK(children.size() == 2);
    
    Introspective::node_type root2;
    CHECK(root2.nchildren() == 0);
    root2.adopt_children(children);
    CHECK(root2.nchildren() == 2);
    CHECK(children.size() == 0);
}
TEST_CASE("nary tree take children")
{
    auto root = make_layered_tree({2,4,8});
    CHECK(root->nchildren() == 2);
    
    Introspective::node_type root2;
    root2.take_children(*root);
    CHECK(root2.nchildren() == 2);
    CHECK(root->nchildren() == 0);

    for (auto* child : root2.children()) {
        CHECK(child != nullptr);
    }
}

// TEST_CASE("nary tree order")
// {
//     auto root = make_layered_tree({2,4,8});
//     for (const auto& n : root->depth()) {

//     }
// }
