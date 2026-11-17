// DV2-MAX space points (a bundle of format dv2max-bundle-v2, e.g. aurora_space):
// how TaggerCheckNeutrino reads the bundle's space section from the DV2_Load
// report, and the touch rule that picks a call's clusters
// (bee_touching_clusters, shared with the T_dlvtx_space writer).  The points
// themselves come from bee_cluster_points on a real grouping; their equality
// with T_dlvtx_space is checked in-chain (sbnd-vertexer aurora/wct/).

#include "WireCellClus/BeeClusterPoints.h"
#include "WireCellClus/TaggerCheckNeutrino.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <string>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus;

TEST_CASE("dv2 space: the load report's space section gives the touch distance")
{
    std::string err;
    CHECK(TaggerCheckNeutrino::dv2_space_touch_cm(R"({"bundle_sha256": "ab", "space": null})", err) == 0.0);
    CHECK(err.empty());                                                       // cloud-only bundle
    CHECK(TaggerCheckNeutrino::dv2_space_touch_cm(R"({"bundle_sha256": "ab"})", err) == 0.0);
    CHECK(err.empty());                                                       // report without the field
    CHECK(TaggerCheckNeutrino::dv2_space_touch_cm(
              R"({"space": {"touch_cm": 1.0, "match_cm": 1.0, "q_ratio": 0.2474, "link_cm": 3.0}})", err) == 1.0);
    CHECK(err.empty());
    CHECK(TaggerCheckNeutrino::dv2_space_touch_cm(R"({"space": {"touch_cm": 2.5}})", err) == 2.5);
}

TEST_CASE("dv2 space: a bad report or space section is an error")
{
    std::string err;
    for (const std::string& rep : {std::string("not json"), std::string("[1, 2]"),
                                   std::string(R"({"space": {"match_cm": 1.0}})"),
                                   std::string(R"({"space": {"touch_cm": "1"}})"),
                                   std::string(R"({"space": {"touch_cm": 0}})"),
                                   std::string(R"({"space": {"touch_cm": -1}})"),
                                   std::string(R"({"space": 3})")}) {
        CAPTURE(rep);
        CHECK(TaggerCheckNeutrino::dv2_space_touch_cm(rep, err) == 0.0);
        CHECK_FALSE(err.empty());
    }
}

TEST_CASE("dv2 space: a cluster is taken whole when any point is within touch_cm (inclusive)")
{
    // cloud: two points on the x axis
    const std::vector<float> cx{0.f, 10.f}, cy{0.f, 0.f}, cz{0.f, 0.f};
    std::vector<std::vector<float>> clusters{
        {50.f, 50.f, 50.f, 0.f, 1.0f, 0.f},     // 0: second point exactly 1 cm from cloud point 0 -> taken
        {0.f, 1.01f, 0.f},                      // 1: 1.01 cm -> not taken
        {10.6f, 0.6f, 0.6f, 80.f, 0.f, 0.f},    // 2: |(0.6, 0.6, 0.6)| = 1.039 cm -> not taken
        {10.5f, 0.5f, 0.5f},                    // 3: 0.866 cm from cloud point 1 -> taken
        {},                                     // 4: no points -> not taken
    };
    const auto t = bee_touching_clusters(clusters, cx, cy, cz, 1.0);
    CHECK(t == std::vector<size_t>{0, 3});                                   // increasing order
    CHECK(bee_touching_clusters(clusters, cx, cy, cz, 2.0) == std::vector<size_t>{0, 1, 2, 3});
}

TEST_CASE("dv2 space: no cloud or no touch distance takes nothing")
{
    const std::vector<std::vector<float>> clusters{{0.f, 0.f, 0.f}};
    const std::vector<float> none, zero{0.f};
    CHECK(bee_touching_clusters(clusters, none, none, none, 1.0).empty());
    CHECK(bee_touching_clusters(clusters, zero, zero, zero, 0.0).empty());
    CHECK(bee_touching_clusters(clusters, zero, zero, zero, 1.0) == std::vector<size_t>{0});
}

TEST_CASE("dv2 space: the touch test uses cells of touch_cm, so neighbours across cell edges count")
{
    // cloud point just below a cell edge, cluster point just above it, 0.2 cm apart
    const std::vector<float> cx{0.99f}, cy{-0.01f}, cz{2.99f};
    const std::vector<std::vector<float>> clusters{{1.09f, 0.09f, 3.09f}, {-5.f, -5.f, -5.f}};
    CHECK(bee_touching_clusters(clusters, cx, cy, cz, 1.0) == std::vector<size_t>{0});
}
