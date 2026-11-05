// ai-helper issue 41: the SPINE-style reco<->truth match of TruthMatch.h.
#include "WireCellUtil/doctest.h"
#include "WireCellClus/TruthMatch.h"

using namespace WireCell::Clus::TruthMatch;

TEST_CASE("truth match: IoU both ways, unlabeled points count in the reco size")
{
    // reco 0: 10 points (8 of truth 0, 2 unlabeled); reco 1: 5 points of truth 0, 5 of truth 1
    std::vector<Size> reco{{10, 100}, {10, 100}};
    std::vector<Size> truth{{13, 130}, {5, 50}};
    std::vector<Cell> cells{{0, 0, 8, 80}, {1, 0, 5, 50}, {1, 1, 5, 50}};
    auto r = match(reco, truth, cells);
    REQUIRE(r.reco[0].size() == 1);
    CHECK(r.reco[0][0].id == 0);
    CHECK(r.reco[0][0].iou_n == doctest::Approx(8.0 / (10 + 13 - 8)));
    // reco 1: truth 1 IoU 5/10 beats truth 0 IoU 5/18
    REQUIRE(r.reco[1].size() == 2);
    CHECK(r.reco[1][0].id == 1);
    CHECK(r.reco[1][0].iou_n == doctest::Approx(0.5));
    CHECK(r.reco[1][1].id == 0);
    // truth 0 prefers reco 0
    CHECK(r.truth[0][0].id == 0);
    CHECK(r.reco_mutual[0] == 1);
    CHECK(r.reco_mutual[1] == 1);
    CHECK(r.truth_mutual[0] == 1);
    CHECK(r.truth_mutual[1] == 1);
}

TEST_CASE("truth match: duplicate cells are summed, top-k, ties, min_overlap, empty")
{
    std::vector<Size> reco{{4, 4}};
    std::vector<Size> truth{{2, 2}, {2, 2}, {2, 2}};
    std::vector<Cell> cells{{0, 2, 1, 1}, {0, 2, 1, 1}, {0, 1, 2, 2}, {0, 0, 1, 1}};
    auto r = match(reco, truth, cells, 2);
    REQUIRE(r.reco[0].size() == 2);
    // truth 1 and 2 tie at 2/4: lower index first
    CHECK(r.reco[0][0].id == 1);
    CHECK(r.reco[0][1].id == 2);
    CHECK(r.reco[0][0].shared_n == doctest::Approx(2));
    auto r2 = match(reco, truth, cells, 0, 0.45);
    CHECK(r2.reco[0].size() == 2);   // 0.5, 0.5 kept; 0.2 dropped
    auto e = match({}, {}, {});
    CHECK(e.reco.empty());
    CHECK(e.truth.empty());
}

TEST_CASE("truth match: a reco object with only unlabeled points matches nothing")
{
    std::vector<Size> reco{{3, 3}};
    std::vector<Size> truth{{5, 5}};
    auto r = match(reco, truth, {});
    CHECK(r.reco[0].empty());
    CHECK(r.truth[0].empty());
    CHECK(r.reco_mutual[0] == 0);
}
