// pr_ml doc 28: the per-cascade-cluster t0 rule of ClusteringImageDeghost.
#include "WireCellClus/ImageDeghostFunctions.h"
#include "WireCellUtil/doctest.h"

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
