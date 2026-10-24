// ClusteringBeamDeghost (doc pdvd/125) sets cluster_t0 only on the cascade
// clusters the beam bundle touches.  TensorDM::as_tensors appends the
// same-named local PCs of every node in child order, and Dataset::append
// throws when the FIRST cluster_scalar carries a key a later one lacks
// ("missing keys in append: 1 missing: cluster_t0", PDVD R3 beam0, pr_ml doc
// 24 sec 6); in the opposite order the key is silently dropped (PDVD data).
//
// fill_cluster_t0_before_serialize() must (1) make the throwing case
// serialize, (2) keep the key in the dropping case, in both with the value
// get_cluster_t0() already returns for an absent key (0), and (3) leave a
// grouping where no cluster or every cluster has the key byte-identical.

#include "WireCellUtil/doctest.h"

#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellAux/TensorDMpointtree.h"

#include <string>
#include <vector>

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;

// A grouping of ncl clusters, each with one blob child and a cluster_scalar
// PC holding "ident"; the clusters listed in `with_t0` also get cluster_t0.
static Grouping* make_grouping(Points::node_t& root, int ncl, const std::vector<int>& with_t0)
{
    Grouping* g = root.value.facade<Grouping>();
    for (int i = 0; i < ncl; ++i) {
        Cluster& cl = g->make_child();
        cl.make_child();
        cl.set_ident(100 + i);
    }
    auto& cls = g->children();
    for (int i : with_t0) {
        cls[i]->set_cluster_t0(1.5);
    }
    return g;
}

// Content of a tensor vector: per tensor, metadata JSON + dtype + shape + bytes.
static std::vector<std::string> content(const ITensor::vector& tens)
{
    std::vector<std::string> ret;
    for (const auto& t : tens) {
        std::string s = t->metadata().toStyledString() + "|" + t->dtype() + "|";
        for (auto d : t->shape()) s += std::to_string(d) + ",";
        s += "|" + std::string(reinterpret_cast<const char*>(t->data()), t->size());
        ret.push_back(s);
    }
    return ret;
}

TEST_CASE("cluster_t0 on the first cluster only: throws, then serializes after the fill")
{
    Points::node_t root;
    Grouping* g = make_grouping(root, 3, {0});

    // The bug: the first cluster_scalar has a key the later ones lack.
    CHECK_THROWS(Aux::TensorDM::as_tensors(root, "pointtrees/0"));

    auto& cls = g->children();
    const double before1 = cls[1]->get_cluster_t0();
    CHECK(fill_cluster_t0_before_serialize(*g) == 2);
    CHECK(cls[0]->get_cluster_t0() == 1.5);
    CHECK(cls[1]->get_cluster_t0() == before1);   // what a reader saw before
    CHECK(cls[2]->get_cluster_t0() == 0.0);
    CHECK_NOTHROW(Aux::TensorDM::as_tensors(root, "pointtrees/0"));

    // A second call has nothing left to fill.
    CHECK(fill_cluster_t0_before_serialize(*g) == 0);
}

TEST_CASE("cluster_t0 on a later cluster only: filled, so the key is no longer dropped")
{
    Points::node_t root;
    Grouping* g = make_grouping(root, 3, {1});
    CHECK_NOTHROW(Aux::TensorDM::as_tensors(root, "pointtrees/0"));
    CHECK(fill_cluster_t0_before_serialize(*g) == 2);
    auto& cls = g->children();
    CHECK(cls[0]->get_cluster_t0() == 0.0);
    CHECK(cls[1]->get_cluster_t0() == 1.5);
    CHECK(cls[2]->get_cluster_t0() == 0.0);
    CHECK_NOTHROW(Aux::TensorDM::as_tensors(root, "pointtrees/0"));
}

TEST_CASE("no cluster_t0 anywhere, or on every cluster: no-op, tensors identical")
{
    Points::node_t root_none, root_none_ref, root_all, root_all_ref;
    Grouping* gn = make_grouping(root_none, 3, {});
    make_grouping(root_none_ref, 3, {});
    Grouping* ga = make_grouping(root_all, 3, {0, 1, 2});
    make_grouping(root_all_ref, 3, {0, 1, 2});
    CHECK(fill_cluster_t0_before_serialize(*gn) == 0);
    CHECK(fill_cluster_t0_before_serialize(*ga) == 0);
    CHECK(content(Aux::TensorDM::as_tensors(root_none, "pointtrees/0"))
          == content(Aux::TensorDM::as_tensors(root_none_ref, "pointtrees/0")));
    CHECK(content(Aux::TensorDM::as_tensors(root_all, "pointtrees/0"))
          == content(Aux::TensorDM::as_tensors(root_all_ref, "pointtrees/0")));
}
