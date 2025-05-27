#include "WireCellClus/Facade.h"
#include "WireCellClus/IPCTransform.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/PointTesting.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"


#include <unordered_map>

using namespace WireCell;
using namespace WireCell::PointTesting;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;
using fa_float_t = WireCell::Clus::Facade::float_t;
using fa_int_t = WireCell::Clus::Facade::int_t;
// WireCell::PointCloud::Tree::scoped_pointcloud_t
using spdlog::debug;
using spdlog::warn;

using node_ptr = std::unique_ptr<Points::node_t>;

void print_ds(const Dataset& ds) {
    std::stringstream ss;
    ss << " size_major " << ds.size_major() << std::endl;
    for (const auto& key : ds.keys()) {
        ss << key << ": ";
        // auto arr = ds.get(key)->elements<double>();
        // for(auto elem : arr) {
        //     ss << elem << " ";
        // }
        ss << std::endl;
    }
    debug(ss.str());
}

// No more explicit DisjointDataset.  It is a PointCloud::Tree::scoped_pointcloud_t.
template <typename DisjointDataset>
void print_dds(const DisjointDataset& dds) {
    for (size_t idx=0; idx<dds.size(); ++idx) {
        const Dataset& ds = dds[idx];
        print_ds(ds);
    }
}

// Return node with two children nodes.
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
            {"wpid", Array({(fa_int_t)WirePlaneId(WirePlaneLayer_t::kAllLayers).ident()})},
            {"npoints", Array({(fa_int_t)10})},
            {"slice_index_min", Array({(fa_int_t)0})},
            {"slice_index_max", Array({(fa_int_t)1})},
            {"u_wire_index_min", Array({(fa_int_t)0})},
            {"u_wire_index_max", Array({(fa_int_t)1})},
            {"v_wire_index_min", Array({(fa_int_t)0})},
            {"v_wire_index_max", Array({(fa_int_t)1})},
            {"w_wire_index_min", Array({(fa_int_t)0})},
            {"w_wire_index_max", Array({(fa_int_t)1})},
            {"max_wire_interval", Array({(fa_int_t)1})},
            {"min_wire_interval", Array({(fa_int_t)1})},
            {"max_wire_type", Array({(fa_int_t)0})},
            {"min_wire_type", Array({(fa_int_t)0})},
        })},
        {"3d", make_janky_track(Ray(Point(0, 0, 0), Point(1, 0, 0)))}
        }));

    const Dataset& pc1 = n1->value.local_pcs().at("3d");
    debug("pc1: {}", pc1.size_major());
    
    fa_int_t wmin = 100;
    fa_int_t wmax = 101;
    // Ibid from a different track
    auto* n2 = root->insert(Points({
        {"scalar", Dataset({
            {"charge", Array({(fa_float_t)2.0})},
            {"center_x", Array({(fa_float_t)1.5})},
            {"center_y", Array({(fa_float_t)0.})},
            {"center_z", Array({(fa_float_t)0.})},
            {"wpid", Array({(fa_int_t)WirePlaneId(WirePlaneLayer_t::kAllLayers).ident()})},
            {"npoints", Array({(fa_int_t)10})},
            {"slice_index_min", Array({(fa_int_t)0})},
            {"slice_index_max", Array({(fa_int_t)1})},
            {"u_wire_index_min", Array({(fa_int_t)wmin})},
            {"u_wire_index_max", Array({(fa_int_t)wmax})},
            {"v_wire_index_min", Array({(fa_int_t)wmin})},
            {"v_wire_index_max", Array({(fa_int_t)wmax})},
            {"w_wire_index_min", Array({(fa_int_t)wmin})},
            {"w_wire_index_max", Array({(fa_int_t)wmax})},
            {"max_wire_interval", Array({(fa_int_t)1})},
            {"min_wire_interval", Array({(fa_int_t)1})},
            {"max_wire_type", Array({(fa_int_t)0})},
            {"min_wire_type", Array({(fa_int_t)0})},
        })},
        {"3d", make_janky_track(Ray(Point(1, 0, 0), Point(2, 0, 0)))}
        }));

    const Dataset& pc2 = n2->value.local_pcs().at("3d");
    debug("pc2: {}", pc2.size_major());

    REQUIRE(pc1 != pc2);
    REQUIRE_FALSE(pc1 == pc2);

    return root;
}

TEST_CASE("clustering prototype point tree")
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
    Scope scope_3d_raw{ "3d", {"x","y","z"}};
    auto const& s3d = rval.scoped_view({ "3d", {"x","y","z"}});

    auto const& pc3d = s3d.pcs();
    CHECK(pc3d.size() == 2);
    // print_dds(pc3d);

    // auto const& pccenter = rval.scoped_view({ "center", {"x","y","z"}}).pcs();
    // print_dds(pccenter);

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
        for (const auto& name : scope_3d_raw.coords) {
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


TEST_CASE("clustering prototype facade")
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
    REQUIRE(blobs[0]->center_x() == 0.5);
    REQUIRE(blobs[1]->center_x() == 1.5);
    REQUIRE(blobs[0]->charge() == 1);
    REQUIRE(blobs[1]->charge() == 2);

    double expect = 0;
    expect += blobs[0]->charge() * blobs[0]->center_x();
    expect += blobs[1]->charge() * blobs[1]->center_x();
    expect /= blobs[0]->charge() + blobs[1]->charge();
    debug("expect average pos {}", expect);
    // there is now another calc_ave_pos(const geo_point_t& origin, int N) const;
    auto ave_pos = pcc.calc_ave_pos({1,0,0}, 1.0);
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

    size_t idx11 = pcc.get_closest_point_index({1.1,0,0});
    size_t idx5 = pcc.get_closest_point_index({0.5,0,0});
    CHECK(idx5 == 5);
    CHECK(idx11 == 11);
    debug("idx5 {} idx11 {} | expecting 5, 11", idx5, idx11);
}


// static void print_MCUGraph(const MCUGraph& g) {
//     std::stringstream ss;
//     ss << "MCUGraph:" << std::endl;
//     ss << "Vertices: " << num_vertices(g) << std::endl;
//     ss << "Edges: " << num_edges(g) << std::endl;
//     ss << "Vertex Properties:" << std::endl;
//     auto vrange = boost::vertices(g);
//     for (auto vit = vrange.first; vit != vrange.second; ++vit) {
//         auto v = *vit;
//         ss << "Vertex " << v << ": Index = " << g[v].index << std::endl;
//     }
//     ss << "Edge Properties:" << std::endl;
//     auto erange = boost::edges(g);
//     auto weightMap = get(boost::edge_weight, g);
//     for (auto eit = erange.first; eit != erange.second; ++eit) {
//         auto e = *eit;
//         ss << "Edge " << e << ": Distance = " << get(weightMap, e) << std::endl;
//     }
//     debug(ss.str());
// }

TEST_CASE("clustering prototype pca")
{
    Points::node_t root_node;
    Grouping* grouping = root_node.value.facade<Grouping>();
    REQUIRE(grouping != nullptr);
    root_node.insert(make_simple_pctree());
    Cluster* pccptr = grouping->children()[0];
    REQUIRE(pccptr != nullptr);
    REQUIRE(pccptr->grouping() == grouping);
    Cluster& pcc = *pccptr;

    geo_point_t center = pcc.get_pca().center;
    debug("center: {} {} {}", center.x(), center.y(), center.z());
    for (size_t ind=0; ind<3; ++ind) {
        auto axis = pcc.get_pca().axis.at(ind);
        auto val = pcc.get_pca().values.at(ind);
        debug("pca{}: {} {} {} {}", ind, axis.x(), axis.y(), axis.z(), val);
    }
}

TEST_CASE("clustering prototype quickhull")
{
    Points::node_t root_node;
    Grouping* grouping = root_node.value.facade<Grouping>();
    REQUIRE(grouping != nullptr);
    root_node.insert(make_simple_pctree());
    Cluster* pccptr = grouping->children()[0];
    auto bpoints = pccptr->get_hull();
    for (const auto &bpoint : bpoints)
    {
        debug("boundary_point: {}", bpoint);
    }
}


TEST_CASE("clustering prototype create cluster graph")
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

    warn("not actually creating graph as test needs updating to avoid 'Anode is null' exception");
    // pcc.Create_graph();
    // print_MCUGraph(*pcc.get_graph());

    // pcc.dijkstra_shortest_paths(0, true);
}

TEST_CASE("clustering prototype Simple3DPointCloud")
{
    Simple3DPointCloud s3dpc;
    for (size_t ind=0; ind<5; ++ind) {
        s3dpc.add({0.1*ind, 0, 0});
    }
    {
        std::stringstream ss;
        ss << s3dpc;
        debug("s3dpc: {}", ss.str());
        geo_point_t p_test1(-1, 0, 0);
        geo_point_t dir1(1, 0, 0);
        double test_dis = 5;
        double dis_step = 0.1;
        double angle_cut = 10;
        double dis_cut = 0.2;
        auto [ind1, dis1] = s3dpc.get_closest_point_along_vec(p_test1, dir1, test_dis, dis_step, angle_cut, dis_cut);
        CHECK(ind1 == 0);
        CHECK(dis1 == 1);
        debug("ind1={} dis1={}", ind1, dis1);
        geo_point_t p_test2(2, 0, 0);
        geo_point_t dir2(-1, 0, 0);
        auto [ind2, dis2] = s3dpc.get_closest_point_along_vec(p_test2, dir2, test_dis, dis_step, angle_cut, dis_cut);
        CHECK(ind2 == 4);
        CHECK(dis2 == 1.6);
        debug("ind2={} dis2={}", ind2, dis2);
    }

    {
        std::stringstream ss;
        ss << s3dpc;
        debug("s3dpc: {}", ss.str());
        Simple3DPointCloud s3dpc2;
        for (size_t ind=0; ind<5; ++ind) {
            s3dpc2.add({1+0.1*ind, 0, 0});
        }
        auto [ind1, ind2, dis] = s3dpc.get_closest_points(s3dpc2);
        CHECK(ind1 == 4);
        CHECK(ind2 == 0);
        CHECK(dis == 0.6);
        debug("ind1={} ind2={} dis={}", ind1, ind2, dis);
    }
}


// static IPCTransformSet::pointer get_pcts()
// {
//     PluginManager& pm = PluginManager::instance();
//     pm.add("WireCellClus");
    
//     {
//         auto icfg = Factory::lookup<IConfigurable>("DetectorVolumes");
//         auto cfg = icfg->default_configuration();
//         icfg->configure(cfg);
//     }
//     {
//         auto icfg = Factory::lookup<IConfigurable>("PCTransformSet");
//         auto cfg = icfg->default_configuration();
//         icfg->configure(cfg);
//     }
//     {
//         auto icfg = Factory::lookup<IConfigurable>("PCTransformSet");
//         auto cfg = icfg->default_configuration();
//         icfg->configure(cfg);
//     }

//     return Factory::find_tn<IPCTransformSet>("PCTransformSet");
// }

// TEST_CASE("clustering prototype dijkstra_shortest_paths")
// {
//     auto pcts = get_pcts();

//     Points::node_t root_node;
//     Grouping* grouping = root_node.value.facade<Grouping>();
//     REQUIRE(grouping != nullptr);
//     root_node.insert(make_simple_pctree());
//     Cluster* pccptr = grouping->children()[0];
//     REQUIRE(pccptr != nullptr);
//     REQUIRE(pccptr->grouping() == grouping);
//     Cluster& pcc = *pccptr;
//     pcc.Create_graph(pcts, false);
//     print_MCUGraph(*pcc.get_graph());
//     pcc.dijkstra_shortest_paths(pcts, 5, false);
// }


TEST_CASE("clustering prototype Facade separate")
{
    Points::node_t root_node;
    Grouping* grouping = root_node.value.facade<Grouping>();
    REQUIRE(grouping != nullptr);
    root_node.insert(make_simple_pctree());
    REQUIRE(grouping->nchildren() == 1);

    Cluster* pccptr = grouping->children()[0];
    REQUIRE(pccptr != nullptr);
    REQUIRE(pccptr->grouping() == grouping);
    REQUIRE(pccptr->nchildren() == 2);

    std::vector<int> groups = {42, 128};
    // auto id2cluster = pcc.separate<Cluster, Grouping>(groups);
    auto id2cluster = grouping->separate(pccptr, groups, true);
    REQUIRE(pccptr == nullptr);
    debug("separate into {} clusters", id2cluster.size());
    REQUIRE(grouping->nchildren() == 2);
    CHECK(id2cluster.size() == 2);
    for (auto [id, cluster] : id2cluster) {
        REQUIRE(cluster != nullptr);
        debug("cluster {} has {} children {} points", id, cluster->nchildren(), cluster->npoints());
        CHECK(cluster->nchildren() == 1);
        CHECK(cluster->npoints() == 10);
    }

    debug("before removal, grouping has {} children", grouping->nchildren());
    // clusters[1]->node()->parent->remove(clusters[1]->node());
    auto child_node = grouping->remove_child(*id2cluster[128]);
    REQUIRE(child_node != nullptr);
    REQUIRE(child_node.get() != nullptr);
    debug("after removal, grouping has {} children", grouping->nchildren());
    REQUIRE(grouping->nchildren() == 1);
}



/**
   This gives an example of how to do the following:
   - add xc,yc,zc arrays representing "corrected coordinates" to blob-local "3d" PC.
   - create a "filtered scoped view" on xc,yc,zc.

   - create an equivalent scoped view on x,y,z

   See also the test "point tree filtered scoped view" in util/test/doctest_pointtree.cxx

 */
TEST_CASE("clustering prototype corrected coordinates")
{
    Points::node_t root_node;
    root_node.insert(make_simple_pctree()); // cluster 1
    root_node.insert(make_simple_pctree()); // cluster 2

    // We first get all blobs in a scoped view.
    //
    // There is more than one way to do this.  Here, we rely on the coincidence
    // that only blob nodes have a local PC named "3d".
    Scope all_blobs_scope{"3d",{"x","y","z"}};
    auto& all_sv = root_node.value.scoped_view(all_blobs_scope);

    // Now add xc,yc,zc by making a "correction".  The actual "correction" here
    // bogus and just serves as an example.
    for (auto* node : all_sv.nodes()) {
        Dataset& pc3d = node->value.local_pc("3d");
        auto npoints = pc3d.size_major();
        REQUIRE(npoints);

        // Make copy as we will mutate.
        auto xc = Array(*pc3d.get("x"));
        auto yc = Array(*pc3d.get("y"));
        auto zc = Array(*pc3d.get("z"));
        auto xcv = xc.elements<double>();
        auto ycv = yc.elements<double>();
        auto zcv = zc.elements<double>();

        // Here we do a totally bogus "correction" just to make some different
        // arrays.
        for (size_t ind=0; ind<npoints; ++ind) {
            xcv[ind] *= -1;
            ycv[ind] *= -1;
            zcv[ind] *= -1;
        }

        // Now we add these new arrays to the "3d" PC
        pc3d.add("xc", xc);
        pc3d.add("yc", yc);
        pc3d.add("zc", zc);
    }

    // Next we get a scoped view on the "corrected" charge and simulate some
    // blob-level fiducial or other selection.
    struct EveryOther {
        int count{0};
        bool operator()(const Points::node_t& node) {
            // We visit EVERY node in the DFS but can only decide based on nodes
            // that are otherwise "in scope".
            
            //Again, rely on the arrangement that blob nodes are the only with
            //"3d" PCs.  Skip those.
            if (! node.value.has_pc("3d")) return false;

            // At this point a "real" selector might use local_pc("3d") to get
            // the "3d" PC and check its points and make some sophisticated
            // decision.  But, here we just keep half of the blob nodes.
            return (count++)%2; // every other
        }
    };

    Scope everyother_corr_scope{"3d",{"xc","yc","zc"},0, "everyother_corr"};
    auto& eoc_sv = root_node.value.scoped_view(everyother_corr_scope, EveryOther{});

    debug("after keeping every other we have {} nodes in scope", eoc_sv.nodes().size()); // 74 for 2 janky tracks??
    

    // We next make a scoped view to get the same nodes but with the x,y,z
    // coordinates.
    std::set<const Points::node_t*> eo_nodes(eoc_sv.nodes().begin(), eoc_sv.nodes().end());
    Scope everyother_orig_scope{"3d",{"x","y","z"},0, "everyother_orig"};
    auto& eoo_sv = root_node.value.scoped_view(everyother_orig_scope, 
                                              [&](const Points::node_t& node) {
                                                  return eo_nodes.find(&node) != eo_nodes.end();
                                              });
    REQUIRE(eoo_sv.nodes().size() == eoc_sv.nodes().size());
    
    // We can now do a "k-d tree query" on the "corrected" coordinate PC and
    // then use the returned point indices to refer to points in the "original"
    // PC. 

    // Get k-d trees for each SV.  We will query kdc "corr" and index into kdo
    // "orig" points.

    const auto& kdo = eoo_sv.kd();     // orig
    const auto& kdc = eoc_sv.kd();     // corr
    const auto& kdo_pts = kdo.points();
    const auto& kdc_pts = kdc.points();
    CHECK(kdo_pts.size() == kdc_pts.size());

    // Do some random k-d tree query.
    const std::vector<double> origin = {0,0,0};
    auto kdc_nn = kdc.knn(1, origin);
    CHECK(kdc_nn.size() == 1);

    for (const auto& [index, metric] : kdc_nn) {
        debug("query=({},{},{}) orig=({},{},{}) corr=({},{},{}) metric={}",
              origin[0], origin[1], origin[2],
              kdo_pts[0][index], kdo_pts[1][index], kdo_pts[2][index],
              kdc_pts[0][index], kdc_pts[1][index], kdc_pts[2][index],
              metric);

        // We next show some ways that this data is all interrelated.  

        // We can get the node that provided the point at the index:
        auto* node = eoo_sv.node_with_point(index); 

        // That is literally indexing the nodes by the "major index" of the
        // point.  Remember, there is the pair of (major,minor) indices
        // corresponding to the "point" index.  The "major" is simply the index
        // of the node in the SV and the "minor" is the point in the local PC of
        // that node.
        REQUIRE(node == eoo_sv.nodes().at(kdo.major_index(index)));

        // Again, the minor index is the index of the point in the node's local
        // PC.  We can use it in the local PC of the node.
        auto minor_index = kdo.minor_index(index);
        const auto& pc3d = node->value.local_pc("3d");
        double xc = pc3d.get("xc")->element<double>(minor_index);
        REQUIRE(xc == kdc_pts[0][index]); // same!

    }
    
}


TEST_CASE("haiwang")
{
    Points::node_t root_node;
    root_node.insert(make_simple_pctree()); // cluster 1
    // root_node.insert(make_simple_pctree()); // cluster 2

    Scope all_scope{"3d",{"x","y","z"}};
    auto& all_sv = root_node.value.scoped_view(all_scope);
    debug("all_sv has {} nodes", all_sv.nodes().size());
    print_dds(all_sv.pcs());

    Scope smallx_scope{"3d",{"x","y","z"}, 0, "smallx"};
    auto& smallx_sv = root_node.value.scoped_view(smallx_scope,
        [&](const Points::node_t& node) {
            debug("filtering node");
            const auto& lpcs = node.value.local_pcs();
            debug("filtering node with {} local pcs", lpcs.size());
            const auto& it = lpcs.find("3d");
            if (it == lpcs.end()) {
                return false;
            }
            const auto& pc = it->second;
            const auto& x = pc.get("x");
            // const auto xv = x->elements<double>();
            const auto& wpid = pc.get("wpid");
            const auto wpidv = wpid->elements<int>();
            // for (auto val : xv) {
            //     debug("filtering x={}", val);
            //     if (val < 1.) {
            //         return true;
            //     }
            // }
            for (auto val : wpidv) {
                debug("filtering wpid={}", val);
                if (val < 11) {
                    debug("passing wpid={}", val);
                    return true;
                }
            }
            return false;
        }
    );
    debug("smallx_sv has {} nodes", smallx_sv.nodes().size());
    debug("print_dds(smallx_sv.pcs());");
    print_dds(smallx_sv.pcs());
    auto smallx_fp = smallx_sv.flat_pc("3d");
    debug("print_ds(smallx_fp);");
    print_ds(smallx_fp);
    auto smallx_fc = smallx_sv.flat_coords();
    debug("print_ds(smallx_fc);");
    print_ds(smallx_fc);

    {
        WireCell::WirePlaneId wpid{-1};
        debug("wpid: wpid.ident() {} wpid.name() {} ok? {} valid? {}", wpid.ident(), wpid.name(), wpid.valid()? true : false, wpid.valid());
    }

    {
        WireCell::WirePlaneId wpid{0};
        debug("wpid: wpid.ident() {} wpid.name() {} ok? {} valid? {}", wpid.ident(), wpid.name(), wpid.valid()? true : false, wpid.valid());
    }
}
