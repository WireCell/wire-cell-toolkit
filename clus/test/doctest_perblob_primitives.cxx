// Tests for the "perblob" row-order invariant enforced by the Grouping
// primitives (doc 52 §13, option 2; clus/docs/perblob_invariant.md):
//
//  * Grouping::separate() carves the cluster-level N-row "perblob" Dataset
//    across the survivor and the splits, in blob (children) order;
//  * Grouping::merge() concatenates the parts' Datasets onto the target's in
//    adoption order, so any separate/merge sequence keeps row i == child i;
//  * inconsistent inputs (stale row counts, parts without a Dataset while
//    others carry one) DROP the Dataset rather than propagate misalignment.
//
// Also covers the caller-side helper realign_perblob_after_regroup(), which
// restores the invariant when a round trip was done with the RAW NaryTree
// primitives (take_children) instead of Grouping::merge().

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/ClusteringFuncs.h"

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;

// A cluster with nb bare blob children and a "perblob" Dataset whose array
// "a" holds the original blob index and "b" holds index+10.
static Cluster* make_cluster(Grouping& g, int nb)
{
    Cluster& cl = g.make_child();
    std::vector<int> a(nb), b(nb);
    for (int i = 0; i < nb; ++i) {
        cl.make_child();
        a[i] = i;
        b[i] = i + 10;
    }
    Dataset ds;
    ds.add("a", Array(a));
    ds.add("b", Array(b));
    cl.value().local_pcs()["perblob"] = ds;
    return &cl;
}

static std::vector<int> perblob_vals(Cluster* c, const std::string& key)
{
    std::vector<int> ret;
    auto& lpcs = c->value().local_pcs();
    auto it = lpcs.find("perblob");
    if (it == lpcs.end()) return ret;
    auto arr = it->second.get(key);
    if (!arr) return ret;
    auto sp = arr->elements<int>();
    ret.assign(sp.begin(), sp.end());
    return ret;
}

static bool has_perblob(Cluster* c)
{
    auto& lpcs = c->value().local_pcs();
    return lpcs.find("perblob") != lpcs.end();
}

TEST_CASE("perblob separate carves survivor and splits")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    Cluster* cl = make_cluster(*g, 6);
    const std::vector<int> cc = {-1, 0, -1, 1, 0, 1};

    auto splits = g->separate(cl, cc, false);
    REQUIRE(splits.size() == 2);

    // Survivor keeps the cc<0 rows in original relative order.
    CHECK(cl->nchildren() == 2);
    CHECK(perblob_vals(cl, "a") == std::vector<int>{0, 2});
    CHECK(perblob_vals(cl, "b") == std::vector<int>{10, 12});
    // Each split gets its own rows.
    CHECK(perblob_vals(splits[0], "a") == std::vector<int>{1, 4});
    CHECK(perblob_vals(splits[1], "a") == std::vector<int>{3, 5});
    CHECK(perblob_vals(splits[1], "b") == std::vector<int>{13, 15});
}

TEST_CASE("perblob separate/merge round trip keeps row i == child i")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    Cluster* cl = make_cluster(*g, 6);
    const std::vector<int> cc = {-1, 0, -1, 1, 0, 1};

    // Blob facade pointers are stable through take_children; record the
    // original index of each blob by pointer.
    std::map<const Blob*, int> orig_index;
    {
        auto kids = cl->children();
        for (int i = 0; i < (int) kids.size(); ++i) orig_index[kids[i]] = i;
    }

    auto splits = g->separate(cl, cc, false);
    auto cc2 = g->merge(splits, cl);  // map overload, ascending-gid adoption

    REQUIRE(cl->nchildren() == 6);
    CHECK(cc2.size() == 6);
    // Expected child order: kept rows (0,2) then gid 0 (1,4) then gid 1 (3,5).
    const std::vector<int> expect = {0, 2, 1, 4, 3, 5};
    auto kids = cl->children();
    auto a = perblob_vals(cl, "a");
    auto b = perblob_vals(cl, "b");
    REQUIRE(a.size() == 6);
    for (int i = 0; i < 6; ++i) {
        CHECK(orig_index.at(kids[i]) == expect[i]);  // the permutation happened
        CHECK(a[i] == orig_index.at(kids[i]));       // and the rows followed it
        CHECK(b[i] == orig_index.at(kids[i]) + 10);
    }
}

TEST_CASE("perblob vector-merge overload (QLMatching recompose shape)")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    Cluster* cl = make_cluster(*g, 5);
    const std::vector<int> cc = {0, -1, 1, 0, -1};

    auto splits = g->separate(cl, cc, false);
    std::vector<Cluster*> others;
    for (auto& [gid, c] : splits) others.push_back(c);  // ascending gid

    g->merge(others.begin(), others.end(), cl);

    REQUIRE(cl->nchildren() == 5);
    // kept (1,4) then gid 0 (0,3) then gid 1 (2)
    CHECK(perblob_vals(cl, "a") == std::vector<int>{1, 4, 0, 3, 2});
}

TEST_CASE("perblob total-separation round trip (neutrino shape)")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    Cluster* cl = make_cluster(*g, 4);
    const std::vector<int> cc = {1, 0, 1, 0};  // every blob leaves

    auto splits = g->separate(cl, cc, false);
    CHECK(cl->nchildren() == 0);
    // An empty survivor holds no rows; the entry is erased, not left keyless.
    CHECK(!has_perblob(cl));

    g->merge(splits, cl);
    REQUIRE(cl->nchildren() == 4);
    // gid 0 rows (1,3) then gid 1 rows (0,2)
    CHECK(perblob_vals(cl, "a") == std::vector<int>{1, 3, 0, 2});
    CHECK(perblob_vals(cl, "b") == std::vector<int>{11, 13, 10, 12});
}

TEST_CASE("perblob separate with remove=true carves the splits")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    Cluster* cl = make_cluster(*g, 4);
    const std::vector<int> cc = {0, 1, 0, 1};

    auto splits = g->separate(cl, cc, true);
    CHECK(cl == nullptr);
    CHECK(perblob_vals(splits[0], "a") == std::vector<int>{0, 2});
    CHECK(perblob_vals(splits[1], "a") == std::vector<int>{1, 3});
}

TEST_CASE("perblob stale rows are dropped, not carved")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    Cluster* cl = make_cluster(*g, 4);
    // Corrupt the dataset: 3 rows for 4 blobs.
    Dataset bad;
    bad.add("a", Array({7, 8, 9}));
    cl->value().local_pcs()["perblob"] = bad;

    const std::vector<int> cc = {-1, 0, -1, 0};
    auto splits = g->separate(cl, cc, false);
    CHECK(!has_perblob(cl));            // dropped loudly, not misaligned
    CHECK(!has_perblob(splits[0]));     // nothing truthful to give the split
}

TEST_CASE("perblob merge drops when a part carries no dataset")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    Cluster* cl = make_cluster(*g, 3);

    // A foreign cluster with children but no perblob dataset.
    Cluster& alien = g->make_child();
    alien.make_child();
    alien.make_child();

    std::vector<Cluster*> parts = {&alien};
    g->merge(parts.begin(), parts.end(), cl);

    CHECK(cl->nchildren() == 5);
    // Truthful concat impossible: rows for the alien blobs do not exist.
    CHECK(!has_perblob(cl));
}

TEST_CASE("perblob raw-primitive round trip fixed by realign helper")
{
    // A round trip done with the RAW node primitives (the pre-option-2 bug
    // shape): NaryTree separate/adopt bypasses the Grouping enforcement, so
    // the dataset keeps its original row order.  The caller-side helper
    // restores the invariant from the cc used to separate.
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    Cluster* cl = make_cluster(*g, 6);
    const std::vector<int> cc = {-1, 0, -1, 1, 0, 1};

    auto nurseries = cl->node()->separate(cc);        // raw: no carve
    for (auto& [gid, nursery] : nurseries) {
        cl->node()->adopt_children(nursery);          // raw: no concat
    }
    REQUIRE(cl->nchildren() == 6);
    // The children are permuted but the rows are not: rows are stale now.
    auto a0 = perblob_vals(cl, "a");
    CHECK(a0 == std::vector<int>{0, 1, 2, 3, 4, 5});

    CHECK(realign_perblob_after_regroup(*cl, cc));
    CHECK(perblob_vals(cl, "a") == std::vector<int>{0, 2, 1, 4, 3, 5});
    CHECK(perblob_vals(cl, "b") == std::vector<int>{10, 12, 11, 14, 13, 15});
}
