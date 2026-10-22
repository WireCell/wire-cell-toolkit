/** Matching on hand-built overlaps: purity/efficiency/IoU, the
 *  Pandora cut, best-for-reco, one-to-one assignment, PQ and ARI.
 */

#include "WireCellTruth/Matching.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell::Truth;

TEST_CASE("intersection metrics in one pass")
{
    Intersection I;
    // reco 0 = 80 of truth A(1) + 20 of truth B(2); reco 1 = 30 of B; reco 2 ghost 10
    I.add(0, 1, 80);
    I.add(0, 2, 20);
    I.add(1, 2, 30);
    I.add_reco_only(2, 10);
    I.add_truth_only(3, 50);   // lost truth

    CHECK(I.nreco() == 3);
    CHECK(I.ntruth() == 3);
    auto p = I.pair(0, 1);
    CHECK(p.shared == 80);
    CHECK(p.purity == doctest::Approx(0.8));
    CHECK(p.efficiency == doctest::Approx(1.0));
    CHECK(p.iou == doctest::Approx(0.8));
    auto q = I.pair(0, 2);
    CHECK(q.purity == doctest::Approx(0.2));
    CHECK(q.efficiency == doctest::Approx(0.4));
    CHECK(q.iou == doctest::Approx(20.0 / 130.0));
    CHECK(I.pair(2, 1).shared == 0);
    CHECK(I.pair(2, 1).iou == 0);

    CHECK(I.best_for_reco(0).truth == 1);
    CHECK(I.best_for_reco(2).truth == -1);
    CHECK(I.best_for_truth(2).reco == 1);

    CHECK(Intersection::pandora_good(p));
    CHECK_FALSE(Intersection::pandora_good(q));          // purity 0.2
    CHECK_FALSE(Intersection::pandora_good(p, 100));     // shared below the minimum

    auto m = I.one_to_one(0.5);
    REQUIRE(m.size() == 2);
    CHECK(m[0].reco == 0);
    CHECK(m[0].truth == 1);
    CHECK(m[1].reco == 1);
    CHECK(m[1].truth == 2);
    CHECK(m[1].iou == doctest::Approx(30.0 / 50.0));

    auto pq = I.panoptic(0.5);
    CHECK(pq.tp == 2);
    CHECK(pq.fp == 1);
    CHECK(pq.fn == 1);
    CHECK(pq.sq == doctest::Approx(0.5 * (0.8 + 0.6)));
    CHECK(pq.rq == doctest::Approx(2.0 / 3.0));
    CHECK(pq.pq == doctest::Approx(pq.sq * pq.rq));

    CHECK(I.pairs().size() == 3);
}

TEST_CASE("one-to-one is greedy by IoU, not by first seen")
{
    Intersection I;
    I.add(0, 1, 60);   // reco 0: iou with 1 = 60/100
    I.add(0, 2, 40);   // reco 0: iou with 2 = 40/(100+90-40)=0.267
    I.add(1, 2, 50);   // reco 1: iou with 2 = 50/90 = 0.556
    auto m = I.one_to_one(0.25);
    REQUIRE(m.size() == 2);
    CHECK(m[0].reco == 0); CHECK(m[0].truth == 1);
    CHECK(m[1].reco == 1); CHECK(m[1].truth == 2);
}

TEST_CASE("adjusted rand index")
{
    std::vector<int> a = {0, 0, 0, 1, 1, 1, 2, 2, 2};
    CHECK(Intersection::ari(a, a) == doctest::Approx(1.0));
    std::vector<int> b = {5, 5, 5, 7, 7, 7, 9, 9, 9};
    CHECK(Intersection::ari(a, b) == doctest::Approx(1.0));          // relabelled
    std::vector<int> c = {0, 0, 1, 1, 2, 2, 0, 1, 2};
    CHECK(Intersection::ari(a, c) < 0.2);
    std::vector<int> mask = {1, 1, 1, 1, 1, 1, 0, 0, 0};
    std::vector<int> d = {0, 0, 0, 1, 1, 1, 7, 8, 9};
    CHECK(Intersection::ari(a, d, mask) == doctest::Approx(1.0));     // the tail is masked out
    CHECK(Intersection::ari(a, d) < 1.0);
    CHECK(Intersection::ari({}, {}) == doctest::Approx(1.0));
}
