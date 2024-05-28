#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/PointTesting.h"

#include "WireCellUtil/doctest.h"

#include "WireCellUtil/Logging.h"

#include <unordered_map>

using namespace WireCell;
using namespace WireCell::PointTesting;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using spdlog::debug;

using node_ptr = std::unique_ptr<Points::node_t>;


TEST_CASE("point tree scope")
{
    Scope s;
    CHECK(s.pcname == "");
    CHECK(s.coords.empty());
    CHECK(s.depth == 0);

    Scope s0{ "pcname", {"x","y","z"}, 0};
    Scope s1{ "pcname", {"x","y","z"}, 1};
    Scope s2{ "PCNAME", {"x","y","z"}, 0};
    Scope s3{ "pcname", {"X","Y","Z"}, 0};
    Scope sc{ "pcname", {"x","y","z"}, 0};

    CHECK(s0 == s0);
    CHECK(s0 == sc);
    CHECK(s0 != s1);
    CHECK(s0 != s2);
    CHECK(s0 != s3);

    CHECK(s0.hash() == sc.hash());
    CHECK(s0.hash() != s1.hash());
    CHECK(s0.hash() != s2.hash());
    CHECK(s0.hash() != s3.hash());

    std::unordered_map<Scope, size_t> m;
    m[s1] = 1;
    m[s2] = 2;
    CHECK(m[s0] == 0);          // spans default
    CHECK(m[s1] == 1);
    CHECK(m[s2] == 2);

}

TEST_CASE("point tree no points")
{
    Points p;
    CHECK(p.node() == nullptr);
    auto& lpcs = p.local_pcs();
    CHECK(lpcs.empty());

    Scope s;
    auto const& sv = p.scoped_view(s);
    auto const& pcr = sv.pcs();
    CHECK(pcr.size() == 0);

    const auto& kd = sv.kd();
    CHECK(kd.points().size() == 0);
}

static
Points::node_ptr make_simple_pctree()
{
    debug("make simple pctree root node");
    Points::node_ptr root = std::make_unique<Points::node_t>();

    // Insert a child with a set of named points clouds with one point
    // cloud from a track.
    debug("node 1 insert");
    auto* n1 = root->insert(Points({ {"3d", make_janky_track()} }));

    debug("node 1 get 3d");
    const Dataset& pc1 = n1->value.local_pcs().at("3d");

    // Ibid from a different track
    debug("node 2 insert");
    auto* n2 = root->insert(Points({ {"3d", make_janky_track(
                        Ray(Point(-1, 2, 3), Point(1, -2, -3)))} }));

    debug("node 2 get 3d");
    const Dataset& pc2 = n2->value.local_pcs().at("3d");

    REQUIRE(pc1 != pc2);
    REQUIRE_FALSE(pc1 == pc2);

    debug("return simple pctree");
    return root;
}

TEST_CASE("point tree with points")
{
    debug("make simple tree");
    auto root = make_simple_pctree();
    debug("got simple tree");
    CHECK(root.get());

    Points& rval = root->value;

    CHECK(root->nchildren() == 2);
    CHECK(root.get() == rval.node());
    CHECK(rval.local_pcs().empty());
    
    debug("check children");
    {
        for (auto& cval : root->child_values()) {
            const auto& pcs = cval.local_pcs();
            CHECK(pcs.size() > 0);
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
            auto sel = pc3d.selection({"x","y","z","q"});
            for (size_t ind=0; ind< pc3d.size_major(); ++ind) {
                debug("\t{} pt=({},{},{}) q={}", ind,
                      sel[0]->element<double>(ind),
                      sel[1]->element<double>(ind),
                      sel[2]->element<double>(ind),
                      sel[3]->element<double>(ind));
            }
        }
    }

    Scope scope{ "3d", {"x","y","z"}};

    debug("request scoped PC at scope = {}", scope);
    ScopedView<double>& sview1 = rval.scoped_view(scope);
    auto& pc3d = sview1.pcs();
    debug("got scoped PC at scope = {} at {} with {} PCs",
          scope, (void*)&pc3d, pc3d.size());
    CHECK(pc3d.size() == 2);
    CHECK(sview1.nodes().size() == 2);
    CHECK(sview1.npoints() > 100);

    // reget sview to also check if that causes any monkey business 
    debug("request k-d tree at scope = {}", scope);
    ScopedView<double>& sview = rval.scoped_view(scope);
    const auto& kd = sview.kd();
    debug("got scoped k-d tree at scope = {} at {} with {} points and {} dimensions",
          scope, (void*)&kd, kd.npoints(), kd.ndim());

    CHECK(kd.ndim() == 3);      // 
    CHECK(kd.nblocks() == 2);   // 2 blocks of ~75 points
    CHECK(kd.npoints() > 100);  // 
    CHECK(kd.npoints() == sview1.npoints());

    const auto& pts = kd.points();
    {
        const size_t ndim = pts.size();
        debug("kd has {} dimensions", ndim);

        CHECK(ndim == 3);

        for (const auto& dim : pts) {
            CHECK(dim.size() == kd.npoints());
        }
    }

    const std::vector<double> origin = {0,0,0};
    auto knn = kd.knn(6, origin);
    for (const auto& [index, metric] : knn) {
        debug("knn: pt=({},{},{}) dist={}",
              pts[0][index], pts[1][index], pts[2][index], metric);
    }
    CHECK(knn.size() == 6);


    auto rad = kd.radius(.001, origin);
    for (const auto& [index, metric] : rad) {
        debug("rad: pt=({},{},{}) dist={}",
              pts[0][index], pts[1][index], pts[2][index], metric);
    }
    CHECK(rad.size() == 6);

}

TEST_CASE("point tree remove node")
{
    auto root = make_simple_pctree();
    Scope scope{ "3d", {"x","y","z"}};
    auto& rval = root->value;

    // Note, we are about to invalidate the scoped while keeping alive the nodes
    // and since the nodes ultimately hold the dataset, these references should
    // remain valid.  
    const auto& pc3d_orig = rval.scoped_view(scope).pcs();
    const Dataset& pc3d_one = pc3d_orig.at(0);
    const Dataset& pc3d_two = pc3d_orig.at(1);
    
    SUBCASE("remove child one") {
        const size_t nleft = pc3d_two.size_major();
        auto cptr = root->children()[0];
        // We get back as unique_ptr so node "dead" stays alive for the context.
        auto dead = root->remove(cptr);
        CHECK(dead);
        CHECK(root->nchildren() == 1);

        const auto& pc3d = rval.scoped_view(scope).pcs();
        CHECK(pc3d.size() == 1);
        CHECK(pc3d[0].get() == pc3d_two);

        const auto& kd = rval.scoped_view(scope).kd();
        CHECK(kd.npoints() == nleft);
    }
    SUBCASE("remove child two") {
        const size_t nleft = pc3d_one.size_major();
        auto cptr = root->children()[1];
        // We get back as unique_ptr so node "dead" stays alive for the context.
        auto dead = root->remove(cptr);
        CHECK(dead);
        CHECK(root->nchildren() == 1);

        const auto& pc3d = rval.scoped_view(scope).pcs();
        CHECK(pc3d.size() == 1);
        CHECK(pc3d[0].get() == pc3d_one);

        const auto& kd = rval.scoped_view(scope).kd();
        CHECK(kd.npoints() == nleft);
    }
    SUBCASE("shared cached point cloud") {
        Scope sxy{ "3d", {"x","y"}};
        Scope syz{ "3d", {"y","z"}};
        auto& pcxy = rval.scoped_view(sxy).pcs();
        auto& pcyz = rval.scoped_view(syz).pcs();
        CHECK(pcxy.size() == pcyz.size());
    }
}

TEST_CASE("point tree merge trees")
{
    auto r1 = make_simple_pctree();
    auto r2 = make_simple_pctree();

    Points::node_t root;
    root.insert(std::move(r1));
    root.insert(std::move(r2));
    CHECK(!r1);
    CHECK(!r2);


    Scope scope{ "3d", {"x","y","z"}};
    auto& rval = root.value;
    const auto& pc3d = rval.scoped_view(scope).pcs();
    CHECK(pc3d.size() == 4);
}

