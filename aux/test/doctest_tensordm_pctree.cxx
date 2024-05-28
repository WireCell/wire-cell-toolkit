#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/PointTesting.h"

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

using namespace WireCell;
using namespace WireCell::PointTesting;
using namespace WireCell::PointCloud;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Tree;

using spdlog::debug;

static
Points::node_ptr make_simple_pctree()
{
    Points::node_ptr root = std::make_unique<Points::node_t>();

    // Insert a child with a set of named points clouds with one point
    // cloud from a track.
    auto* n1 = root->insert(Points({ {"3d", make_janky_track()} }));

    const auto& pc1 = n1->value.local_pcs().at("3d");

    // Ibid from a different track
    auto* n2 = root->insert(Points({ {"3d",
                    make_janky_track(
                        Ray(Point(-1, 2, 3), Point(1, -2, -3)))}}));

    const auto& pc2 = n2->value.local_pcs().at("3d");

    REQUIRE(pc1 != pc2);
    REQUIRE_FALSE(pc1 == pc2);

    return root;
}

template<typename M>
void same_keys(const M& a, const M& b)
{
    for (auto ait : a) {
        auto bit = b.find(ait.first);
        REQUIRE(bit != b.end());
    }
}

template<typename P>
void same_pc(const P& a, const P& b)
{
    REQUIRE(a.size_major() == b.size_major());
    REQUIRE(a.keys().size() == b.keys().size());
    for (auto key : a.keys()) {
        REQUIRE(b.has(key));

        // check array level
        auto aa = a.get(key)->bytes();
        auto bb = b.get(key)->bytes();

        debug("PC array {} with {} =?= {} bytes", key, aa.size(), bb.size());
        REQUIRE(aa.size() == bb.size());
        size_t nbytes = aa.size();
        for (size_t ind=0; ind<nbytes; ++ind) {
            REQUIRE(aa[ind] == bb[ind]);
        }
    }
}


template<typename M>
void same_pcs(const M& a, const M& b)
{
    REQUIRE(a.size() == a.size());

    same_keys(a, b);
    same_keys(b, a);

    for (auto ait : a) {
        auto bit = b.find(ait.first);
        debug("PC {}", ait.first);
        same_pc(ait.second, bit->second);
    }
}


void nodes_equal(Points::node_t* a, Points::node_t* b)
{
    // PCs
    auto pcsa = a->value.local_pcs();
    auto pcsb = b->value.local_pcs();
    same_pcs(pcsa, pcsb);

    // tree structure
    REQUIRE (a->nchildren() == b->nchildren());
    auto acs = a->children();
    auto bcs = b->children();

    for (size_t ind=0; ind<a->nchildren(); ++ind) {
        nodes_equal(acs[ind], bcs[ind]);
    }
}


TEST_CASE("tensordm pctree")
{
    auto root = make_simple_pctree();
    REQUIRE(root);

    const std::string datapath = "root";
    auto tens = as_tensors(*root.get(), datapath);
    REQUIRE(tens.size() > 0);

    debug("pctree as tensor metadata:");
    debug("{:20} {}", "datatype", "datapath");
    for (auto ten : tens) {
        auto md = ten->metadata();
        debug("{:20} {}", md["datatype"].asString(), md["datapath"].asString());
    }

    auto root2 = as_pctree(tens, datapath);
    REQUIRE(root2);

    nodes_equal(root.get(), root2.get());
}

