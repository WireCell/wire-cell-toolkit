// separate(tag_family=true) writes "sep_family" into the cluster_scalar PC of
// some clusters only.  TensorDM::as_tensors appends the same-named local PCs
// of every node in child order, and Dataset::append throws when the FIRST
// cluster_scalar carries a key a later one lacks ("missing keys in append:
// 1 missing: sep_family", wcfm event 124 / anode 2, wcp-porting-img
// wcfm/docs/13); in the opposite order the key is silently dropped.
//
// fill_sep_family_before_serialize() must (1) make the throwing case
// serialize, with 0 filled in, and (2) leave the silently-dropping case
// byte-identical.

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
// PC holding "ident"; the clusters listed in `tagged` also get sep_family.
static Grouping* make_grouping(Points::node_t& root, int ncl, const std::vector<int>& tagged)
{
    Grouping* g = root.value.facade<Grouping>();
    for (int i = 0; i < ncl; ++i) {
        Cluster& cl = g->make_child();
        cl.make_child();
        cl.set_ident(100 + i);
    }
    auto& cls = g->children();
    for (int i : tagged) {
        cls[i]->set_scalar<int>("sep_family", 7);
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

TEST_CASE("sep_family on the first cluster only: throws, then serializes after the fill")
{
    Points::node_t root;
    Grouping* g = make_grouping(root, 3, {0});

    // The bug: the first cluster_scalar has a key the later ones lack.
    CHECK_THROWS(Aux::TensorDM::as_tensors(root, "pointtrees/0"));

    CHECK(fill_sep_family_before_serialize(*g) == 2);
    auto& cls = g->children();
    CHECK(cls[0]->get_scalar<int>("sep_family", -1) == 7);
    CHECK(cls[1]->get_scalar<int>("sep_family", -1) == 0);
    CHECK(cls[2]->get_scalar<int>("sep_family", -1) == 0);
    CHECK_NOTHROW(Aux::TensorDM::as_tensors(root, "pointtrees/0"));

    // A second call has nothing left to fill.
    CHECK(fill_sep_family_before_serialize(*g) == 0);
}

TEST_CASE("sep_family on a later cluster only: the fill is a no-op, tensors identical")
{
    Points::node_t root_a, root_b;
    make_grouping(root_a, 3, {1});
    Grouping* gb = make_grouping(root_b, 3, {1});

    const auto before = content(Aux::TensorDM::as_tensors(root_a, "pointtrees/0"));
    CHECK(fill_sep_family_before_serialize(*gb) == 0);
    const auto after = content(Aux::TensorDM::as_tensors(root_b, "pointtrees/0"));
    CHECK(before == after);
}

TEST_CASE("no sep_family anywhere, or on every cluster: no-op")
{
    Points::node_t root_none, root_all;
    Grouping* gn = make_grouping(root_none, 3, {});
    Grouping* ga = make_grouping(root_all, 3, {0, 1, 2});
    CHECK(fill_sep_family_before_serialize(*gn) == 0);
    CHECK(fill_sep_family_before_serialize(*ga) == 0);
    CHECK_NOTHROW(Aux::TensorDM::as_tensors(root_none, "pointtrees/0"));
    CHECK_NOTHROW(Aux::TensorDM::as_tensors(root_all, "pointtrees/0"));
}
