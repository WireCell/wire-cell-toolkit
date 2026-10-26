// ClusteringExamineBundles writes the per-blob "isolated" array only on the
// clusters inside its scope filter, so a live grouping can reach the pctree
// dump with "perblob" PCs whose key sets differ.  TensorDM::as_tensors appends
// the same-named local PCs of every node in child order: Dataset::append
// throws when the FIRST "perblob" PC carries a key a later one lacks
// ("missing keys in append: 1 missing: isolated", SBND MPV/MPR event mpvmpr_0,
// wcp-porting-img pr_ml/docs/26 sec 6), and in the opposite order the key is
// silently dropped for every cluster.
//
// fill_isolated_before_serialize() (MultiAlgBlobClustering knob fill_isolated,
// default off) gives every cluster that has a "perblob" PC but no "isolated"
// an all -1 array -- the array's documented "this blob is in the main
// sub-component" value, i.e. a cluster nothing split.  Clusters without a
// "perblob" PC are left alone.

#include "WireCellUtil/doctest.h"

#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellAux/TensorDMpointtree.h"

#include <vector>

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;

// A grouping of ncl clusters with nblob blob children each.  Every cluster in
// `perblob` gets a "perblob" PC holding "other" (one int per blob); those also
// in `isolated` get "isolated" there too (values 0 .. nblob-1).
static Grouping* make_grouping(Points::node_t& root, int ncl, int nblob,
                               const std::vector<int>& perblob, const std::vector<int>& isolated)
{
    Grouping* g = root.value.facade<Grouping>();
    for (int i = 0; i < ncl; ++i) {
        Cluster& cl = g->make_child();
        for (int b = 0; b < nblob; ++b) cl.make_child();
        cl.set_ident(100 + i);
    }
    auto& cls = g->children();
    for (int i : perblob) {
        cls[i]->put_pcarray(std::vector<int>(nblob, 3), "other", "perblob");
    }
    for (int i : isolated) {
        std::vector<int> v(nblob);
        for (int b = 0; b < nblob; ++b) v[b] = b;
        cls[i]->put_pcarray(v, "isolated", "perblob");
    }
    return g;
}

TEST_CASE("isolated on the first cluster only: throws, then serializes after the fill")
{
    Points::node_t root;
    Grouping* g = make_grouping(root, 3, 2, {0, 1, 2}, {0});

    // The bug: the first "perblob" PC has a key the later ones lack.
    CHECK_THROWS(Aux::TensorDM::as_tensors(root, "pointtrees/0"));

    CHECK(fill_isolated_before_serialize(*g) == 2);
    auto& cls = g->children();
    for (int i = 1; i < 3; ++i) {
        REQUIRE(cls[i]->has_pcarray<int>("isolated", "perblob"));
        auto a = cls[i]->get_pcarray<int>("isolated", "perblob");
        CHECK(a.size() == 2);
        for (auto v : a) CHECK(v == -1);
    }
    // The cluster that had it keeps its own values.
    auto a0 = cls[0]->get_pcarray<int>("isolated", "perblob");
    CHECK(a0[1] == 1);
    CHECK_NOTHROW(Aux::TensorDM::as_tensors(root, "pointtrees/0"));

    // A second call has nothing left to fill.
    CHECK(fill_isolated_before_serialize(*g) == 0);
}

TEST_CASE("isolated on a later cluster only: the fill gives every cluster the array")
{
    // The silently-dropping order: without the fill the serializer keeps only
    // the first PC's keys, so "isolated" vanishes from the whole grouping.
    Points::node_t root;
    Grouping* g = make_grouping(root, 3, 2, {0, 1, 2}, {1});
    CHECK(fill_isolated_before_serialize(*g) == 2);
    for (auto* cl : g->children()) {
        CHECK(cl->has_pcarray<int>("isolated", "perblob"));
    }
    CHECK_NOTHROW(Aux::TensorDM::as_tensors(root, "pointtrees/0"));
}

TEST_CASE("clusters without a perblob PC are left alone")
{
    Points::node_t root;
    Grouping* g = make_grouping(root, 3, 2, {0, 2}, {0});
    CHECK(fill_isolated_before_serialize(*g) == 1);   // cluster 2 only
    auto& cls = g->children();
    CHECK_FALSE(cls[1]->has_pcarray<int>("isolated", "perblob"));
    CHECK(cls[2]->has_pcarray<int>("isolated", "perblob"));
}

TEST_CASE("every cluster already has isolated: no-op")
{
    Points::node_t root;
    Grouping* g = make_grouping(root, 2, 3, {0, 1}, {0, 1});
    CHECK(fill_isolated_before_serialize(*g) == 0);
}
