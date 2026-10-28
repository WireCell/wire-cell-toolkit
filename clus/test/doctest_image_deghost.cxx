// pr_ml doc 28: the per-cascade-cluster t0 rule of ClusteringImageDeghost.
#include "WireCellClus/ImageDeghostFunctions.h"
#include "WireCellUtil/doctest.h"

#include <algorithm>

using namespace WireCell::Clus::PR;

TEST_CASE("image_deghost_cluster_t0: one t0 per cascade cluster")
{
    // cells: 0,1,2 in cascade cluster 0; 3,4 in cascade cluster 1; 5 in cascade cluster 2
    const std::vector<int> dgc{0, 0, 0, 1, 1, 2};
    const std::vector<int> own{0, 0, 1, 2, -1, 1};
    const std::vector<char> sel{1, 1, 1};
    const std::vector<double> t0{0, 0, 0};
    std::vector<char> take;
    auto r = image_deghost_cluster_t0(dgc, own, sel, t0, 3, 1.0, take);
    REQUIRE(r.size() == 3);
    CHECK(r[0].touched); CHECK(!r[0].conflict); CHECK(r[0].t0 == 0);
    CHECK(r[1].touched); CHECK(r[2].touched);
    CHECK(take == std::vector<char>{1, 1, 1, 1, 0, 1});   // the unowned cell stays
}

TEST_CASE("image_deghost_cluster_t0: unselected owners take nothing")
{
    const std::vector<int> dgc{0, 0, 1};
    const std::vector<int> own{0, 1, 1};
    const std::vector<char> sel{1, 0};
    const std::vector<double> t0{5, 5};
    std::vector<char> take;
    auto r = image_deghost_cluster_t0(dgc, own, sel, t0, 2, 1.0, take);
    CHECK(r[0].touched); CHECK(r[0].t0 == 5);
    CHECK(!r[1].touched);
    CHECK(take == std::vector<char>{1, 0, 0});
}

TEST_CASE("image_deghost_cluster_t0: conflicting t0, majority owner wins, ties to the smaller key")
{
    // cascade cluster 0: owner 0 (t0 10) takes 2 cells, owner 1 (t0 20) takes 1
    const std::vector<int> dgc{0, 0, 0, 1, 1};
    const std::vector<int> own{0, 1, 0, 1, 0};
    const std::vector<char> sel{1, 1};
    const std::vector<double> t0{10, 20};
    std::vector<char> take;
    auto r = image_deghost_cluster_t0(dgc, own, sel, t0, 2, 1.0, take);
    CHECK(r[0].conflict); CHECK(r[0].t0 == 10);
    // cascade cluster 1: one cell each -> tie -> owner 0 (smaller key)
    CHECK(r[1].conflict); CHECK(r[1].t0 == 10);
    CHECK(take == std::vector<char>{1, 0, 1, 0, 1});
    // within tolerance there is no conflict
    const std::vector<double> t0b{10, 10.5};
    auto rb = image_deghost_cluster_t0(dgc, own, sel, t0b, 2, 1.0, take);
    CHECK(!rb[0].conflict); CHECK(!rb[1].conflict);
    CHECK(take == std::vector<char>{1, 1, 1, 1, 1});
}

TEST_CASE("image_deghost_cluster_t0: empty input")
{
    std::vector<char> take{1};
    auto r = image_deghost_cluster_t0({}, {}, {}, {}, 0, 1.0, take);
    CHECK(r.empty()); CHECK(take.empty());
}

// pr_ml doc 32 sec B.4 (29107_3 cluster 118): the cells a t0 conflict leaves behind, grouped by their owner's t0.
TEST_CASE("image_deghost_t0_split: one group per owner t0 left behind")
{
    // cascade cluster 0: owners 0 (t0 10) x3, 1 (t0 20) x2, 2 (t0 20.5, within tol of 20) x1, 3 (t0 30) x1,
    //                    owner 4 unselected (t0 40), one unowned cell
    // cascade cluster 1: owner 1 only -> no conflict, nothing split
    const std::vector<int> dgc{0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    const std::vector<int> own{0, 1, 0, 3, 2, 1, 0, 4, -1, 1};
    const std::vector<char> sel{1, 1, 1, 1, 0};
    const std::vector<double> t0{10, 20, 20.5, 30, 40};
    std::vector<char> take;
    const auto r = image_deghost_cluster_t0(dgc, own, sel, t0, 2, 1.0, take);
    REQUIRE(r[0].conflict);
    CHECK(r[0].t0 == 10);
    CHECK(!r[1].conflict);
    const auto sp = image_deghost_t0_split(dgc, own, sel, t0, r, 1.0);
    REQUIRE(sp.size() == 2);
    CHECK(sp[0].dgc == 0); CHECK(sp[0].t0 == 20);
    CHECK(sp[0].cells == std::vector<size_t>{1, 4, 5});   // owners 1 and 2: one t0 within tol
    CHECK(sp[1].dgc == 0); CHECK(sp[1].t0 == 30);
    CHECK(sp[1].cells == std::vector<size_t>{3});
    // every split cell was left by the one-t0 rule; nothing taken is split
    for (const auto& g : sp) for (size_t i : g.cells) CHECK(!take[i]);
    // without a conflict there is nothing to split
    const std::vector<double> t0same{10, 10, 10, 10, 10};
    const auto r2 = image_deghost_cluster_t0(dgc, own, sel, t0same, 2, 1.0, take);
    CHECK(image_deghost_t0_split(dgc, own, sel, t0same, r2, 1.0).empty());
}

// harmonize_blob_arrays (ClusteringImageDeghost harmonize_3d_keys): after the cells replace a cluster's blobs, the live
// blobs can carry different "3d" key sets (SBND simulation: y_cor / z_cor persisted by the clustering job on the blobs
// not replaced, absent on the cells).  TensorDM::as_tensors then throws ("missing keys in append", SBND mpvmpr_1076) or,
// in the opposite node order, drops the key for every blob.
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellAux/TensorDMpointtree.h"

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;

// ncl clusters x nblob blobs; every blob has a "3d" PC with x, y, z; the blobs of the clusters in `extra` also y_cor.
static Grouping* make_grouping3d(Points::node_t& root, int ncl, int nblob, const std::vector<int>& extra)
{
    Grouping* g = root.value.facade<Grouping>();
    for (int i = 0; i < ncl; ++i) {
        Cluster& cl = g->make_child();
        for (int b = 0; b < nblob; ++b) {
            Blob& bl = cl.make_child();
            Dataset ds;
            ds.add("x", Array(std::vector<double>{1.0 * i, 2.0}));
            ds.add("y", Array(std::vector<double>{3.0, 4.0}));
            ds.add("z", Array(std::vector<double>{5.0, 6.0}));
            if (std::find(extra.begin(), extra.end(), i) != extra.end()) {
                ds.add("y_cor", Array(std::vector<double>{3.5, 4.5}));
            }
            bl.local_pcs()["3d"] = ds;
        }
        cl.set_ident(100 + i);
    }
    return g;
}

TEST_CASE("harmonize_blob_arrays: a key only the first cluster's blobs carry -- throws before, serializes after")
{
    Points::node_t root;
    Grouping* g = make_grouping3d(root, 3, 2, {0});
    CHECK_THROWS(Aux::TensorDM::as_tensors(root, "pointtrees/0"));   // the bug
    std::vector<std::string> partial;
    CHECK(harmonize_blob_arrays(*g, "3d", partial) == 2);            // y_cor of cluster 0's two blobs
    CHECK(partial == std::vector<std::string>{"y_cor"});
    for (auto* cl : g->children()) {
        for (auto* b : cl->children()) {
            const auto& ds = b->local_pcs().at("3d");
            CHECK_FALSE(ds.has("y_cor"));
            CHECK(ds.has("x")); CHECK(ds.has("y")); CHECK(ds.has("z"));
        }
    }
    CHECK_NOTHROW(Aux::TensorDM::as_tensors(root, "pointtrees/0"));
    CHECK(harmonize_blob_arrays(*g, "3d", partial) == 0);            // nothing left
    CHECK(partial.empty());
}

TEST_CASE("harmonize_blob_arrays: one key set already -- no-op")
{
    Points::node_t root;
    Grouping* g = make_grouping3d(root, 2, 2, {0, 1});
    std::vector<std::string> partial;
    CHECK(harmonize_blob_arrays(*g, "3d", partial) == 0);
    CHECK(partial.empty());
    for (auto* cl : g->children()) {
        for (auto* b : cl->children()) CHECK(b->local_pcs().at("3d").has("y_cor"));
    }
}
