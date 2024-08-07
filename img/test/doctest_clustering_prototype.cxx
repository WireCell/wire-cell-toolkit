#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/PointTesting.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include "WireCellImg/PointCloudFacade.h"

#include <unordered_map>

using namespace WireCell;
using namespace WireCell::PointTesting;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::PointCloud::Facade;
using fa_float_t = WireCell::PointCloud::Facade::float_t;
using fa_int_t = WireCell::PointCloud::Facade::int_t;
// WireCell::PointCloud::Tree::scoped_pointcloud_t
using spdlog::debug;

using node_ptr = std::unique_ptr<Points::node_t>;

// No more explicit DisjointDataset.  It is a PointCloud::Tree::scoped_pointcloud_t.
template <typename DisjointDataset>
void print_dds(const DisjointDataset& dds) {
    for (size_t idx=0; idx<dds.size(); ++idx) {
        const Dataset& ds = dds[idx];
        std::stringstream ss;
        ss << "ds: " << idx << std::endl;
        // const size_t len = ds.size_major();
        for (const auto& key : ds.keys()) {
            auto arr = ds.get(key)->elements<fa_float_t>();
            ss << key << ": ";
            for(auto elem : arr) {
                ss << elem << " ";
            }
            ss << std::endl;
        }
        std::cout << ss.str() << std::endl;
    }
}

static
Points::node_ptr make_simple_pctree()
{
    // empty root node
    Points::node_ptr root = std::make_unique<Points::node_t>();

    // Insert a child with a set of named points clouds with one point
    // cloud from a track.

    /// QUESTION: can only do this on construction?
    /// bv: see NaryTree.h for several insert()'s

    /// QUESTION: units?
    /// bv: units are always assumed in WCT system-of-units.  To be correct
    ///     here, we should be multiplying by some [length] unit.

    auto* n1 = root->insert(Points({
        /// QUESTION: proper Array initiation?
        {"scalar", Dataset({
            {"charge", Array({(fa_float_t)1.0})},
            {"center_x", Array({(fa_float_t)0.5})},
            {"center_y", Array({(fa_float_t)0.})},
            {"center_z", Array({(fa_float_t)0.})},
            {"npoints", Array({(fa_int_t)10})},
            {"slice_index_min", Array({(fa_int_t)0})},
            {"slice_index_max", Array({(fa_int_t)1})},
            {"u_wire_index_min", Array({(fa_int_t)0})},
            {"u_wire_index_max", Array({(fa_int_t)1})},
            {"v_wire_index_min", Array({(fa_int_t)0})},
            {"v_wire_index_max", Array({(fa_int_t)1})},
            {"w_wire_index_min", Array({(fa_int_t)0})},
            {"w_wire_index_max", Array({(fa_int_t)1})},
        })},
        {"3d", make_janky_track(Ray(Point(0, 0, 0), Point(1, 0, 0)))}
        }));

    const Dataset& pc1 = n1->value.local_pcs().at("3d");
    debug("pc1: {}", pc1.size_major());
    

    // Ibid from a different track
    auto* n2 = root->insert(Points({
        {"scalar", Dataset({
            {"charge", Array({(fa_float_t)2.0})},
            {"center_x", Array({(fa_float_t)1.5})},
            {"center_y", Array({(fa_float_t)0.})},
            {"center_z", Array({(fa_float_t)0.})},
            {"npoints", Array({(fa_int_t)10})},
            {"slice_index_min", Array({(fa_int_t)0})},
            {"slice_index_max", Array({(fa_int_t)1})},
            {"u_wire_index_min", Array({(fa_int_t)1})},
            {"u_wire_index_max", Array({(fa_int_t)2})},
            {"v_wire_index_min", Array({(fa_int_t)1})},
            {"v_wire_index_max", Array({(fa_int_t)2})},
            {"w_wire_index_min", Array({(fa_int_t)1})},
            {"w_wire_index_max", Array({(fa_int_t)2})},
        })},
        {"3d", make_janky_track(Ray(Point(1, 0, 0), Point(2, 0, 0)))}
        }));

    const Dataset& pc2 = n2->value.local_pcs().at("3d");
    debug("pc2: {}", pc2.size_major());

    REQUIRE(pc1 != pc2);
    REQUIRE_FALSE(pc1 == pc2);

    return root;
}

TEST_CASE("clustering point tree")
{
    // this test does not touch the facades so needs to Grouping root
    auto root = make_simple_pctree(); //  a cluster node.
    CHECK(root.get());

    // from WireCell::NaryTree::Node to WireCell::PointCloud::Tree::Points
    auto& rval = root->value;

    CHECK(root->children().size() == 2);
    CHECK(root.get() == rval.node()); // node raw pointer == point's node
    CHECK(rval.local_pcs().empty());
    
    {
        auto& cval = *(root->child_values().begin());
        const auto& pcs = cval.local_pcs();
        CHECK(pcs.size() > 0);
        // using pointcloud_t = Dataset;
        // key: std::string, val: Dataset
        for (const auto& [key,val] : pcs) {
            debug("child has pc named \"{}\" with {} points", key, val.size_major());
        }
        const auto& pc3d = pcs.at("3d");
        debug("got child PC named \"3d\" with {} points", pc3d.size_major());
        CHECK(pc3d.size_major() > 0);
        CHECK(pc3d.has("x"));
        CHECK(pc3d.has("y"));
        CHECK(pc3d.has("z"));
        CHECK(pc3d.has("q"));
    }

    // name, coords, [depth]
    Scope scope{ "3d", {"x","y","z"}};
    auto const& s3d = rval.scoped_view({ "3d", {"x","y","z"}});

    auto const& pc3d = s3d.pcs();
    CHECK(pc3d.size() == 2);
    print_dds(pc3d);

    auto const& pccenter = rval.scoped_view({ "center", {"x","y","z"}}).pcs();
    print_dds(pccenter);

    const auto& kd = s3d.kd();
    const auto& points = kd.points();

    /// QUESTION: how to get it -> node?
    ///
    /// bv: call "index()" on the disjoint range iterator.  This returns a
    /// pair<size_t,size_t> holding major/minor indices into the disjoint range.
    /// You can use that pair to access elements in other disjoint ranges.  See
    /// doctest-pointtree-example for details.

    std::vector<fa_float_t> some_point = {1, 0, 0};
    auto knn = kd.knn(2, some_point);
    for (const auto& [index, metric] : knn) {
        debug("knn: pt=({},{},{}) metric={}",
              points[0][index], points[1][index], points[2][index], metric);
    }
    CHECK(knn.size() == 2);

    for (const auto& [index, metric] : knn) {
        auto node_index = kd.major_index(index);
        // point-in-node index
        auto pin_index = kd.minor_index(index);
        debug("knn point {} at distance {} from query is in local point cloud {} at local point {}",
              index, metric, node_index, pin_index);
        const Dataset& pc = pc3d[node_index];
        for (const auto& name : scope.coords) {
            debug("\t{} = {}", name, pc.get(name)->element<fa_float_t>(pin_index));
        }
    }

    auto rad = kd.radius(.01, some_point);
    for (const auto& [index, metric] : rad) {
        debug("rad: pt=({},{},{}) metric={}",
              points[0][index], points[1][index], points[2][index], metric);
    }
    CHECK(rad.size() == 2);
}


TEST_CASE("clustering facade")
{
    Points::node_t root_node;
    Grouping* grouping = root_node.value.facade<Grouping>();
    REQUIRE(grouping != nullptr);
    root_node.insert(make_simple_pctree());
    Cluster* pccptr = grouping->children()[0];
    REQUIRE(pccptr != nullptr);
    REQUIRE(pccptr->grouping() == grouping);
    Cluster& pcc = *pccptr;

    CHECK(pcc.sanity());

    auto& blobs = pcc.children();

    // (0.5 * 1 + 1.5 * 2) / 3 = 1.1666666666666665
    debug("blob 0: q={}, r={}", blobs[0]->charge(), blobs[0]->center_x());
    debug("blob 1: q={}, r={}", blobs[1]->charge(), blobs[1]->center_x());
    double expect = 0;
    expect += blobs[0]->charge() * blobs[0]->center_x();
    expect += blobs[1]->charge() * blobs[1]->center_x();
    expect /= blobs[0]->charge() + blobs[1]->charge();
    debug("expect average pos {}", expect);
    auto ave_pos = pcc.calc_ave_pos({1,0,0}, 1);
    debug("ave_pos: {} | expecting (1.1666666666666665 0 0)", ave_pos);
    auto l1 = fabs(ave_pos[0] - 1.1666666666666665) + fabs(ave_pos[1]) + fabs(ave_pos[2]);
    CHECK(l1 < 1e-3);

    const auto vdir_alg0 = pcc.vhough_transform({1,0,0}, 1, Cluster::HoughParamSpace::costh_phi);
    debug("vdir_alg0: {} | expecting around {{1, 0, 0}}", vdir_alg0);
    l1 = fabs(vdir_alg0[0] - 1) + fabs(vdir_alg0[1]) + fabs(vdir_alg0[2]);
    CHECK(l1 < 1e-1);
    const auto vdir_alg1 = pcc.vhough_transform({1,0,0}, 1, Cluster::HoughParamSpace::theta_phi);
    debug("vdir_alg1: {} | expecting around {{1, 0, 0}}", vdir_alg1);
    l1 = fabs(vdir_alg1[0] - 1) + fabs(vdir_alg1[1]) + fabs(vdir_alg1[2]);
    CHECK(l1 < 1e-1);

    // sqrt(2./3.*(3*3*2*2*3)+(0.5*1.101)^2) = 8.50312003032
    const auto length = pcc.get_length();
    debug("length: {} | expecting 8.50312003032", length);
    l1 = fabs(length - 8.50312003032);
    CHECK(l1 < 1e-3);

    const auto [earliest, latest] = pcc.get_earliest_latest_points();
    debug("earliest_latest_points: {} {} | expecting (0 0 0) (1.9 0 0)", earliest, latest);
    l1 = fabs(earliest[0]) + fabs(earliest[1]) + fabs(earliest[2]);
    CHECK(l1 < 1e-3);
    l1 = fabs(latest[0] - 1.9) + fabs(latest[1]) + fabs(latest[2]);
    CHECK(l1 < 1e-3);

    const auto [num1, num2] = pcc.ndipole({0.5,0,0}, {1,0,0});
    debug("num_points: {} {} | expecting 15, 5", num1, num2);
    CHECK(num1 == 15);
    CHECK(num2 == 5);
}
