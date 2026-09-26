// doc pdvd/81 -- the charge-based Michel energy's two pure pieces.
//
//  * stm_michel_combine_planes: the chain's three-plane combination
//    (kine_charge_from_maps, NeutrinoEnergyReco.cxx:149-186) forked for a
//    SIGNED input -- symmetric planes give the weighted mean, an asymmetric
//    pair drops the largest plane, a non-positive (med + max) never trips the
//    switch, an all-zero weight triple gives 0.
//  * TrackFitting::masked_response_prediction: scale * R * (pos masked), the
//    muon-only prediction read off a stored multi-fit response; short masks and
//    short scales read as zero beyond their end.
//  * TrackFitting::Parameters::keep_dqdx_response defaults to 0 (nothing is
//    stored unless a caller asks) and round-trips through set/get_parameter.

#include "WireCellUtil/doctest.h"
#include "WireCellClus/StmMichelFunctions.h"
#include "WireCellClus/TrackFitting.h"

#include <Eigen/Sparse>
#include <array>
#include <vector>

using namespace WireCell;
using WireCell::Clus::TrackFitting;
using WireCell::Clus::PR::stm_michel_combine_planes;
using WireCell::Clus::PR::stm_michel_pick_wire;

// doc pdhd/28: the wrapped-channel wire lookup for the region reading.
TEST_CASE("stm_michel_pick_wire: a single incarnation is itself")
{
    int calls = 0;
    auto p = stm_michel_pick_wire({{3.0, 40.0}}, 10.0, true, true, [&](size_t) { ++calls; return 1; });
    CHECK(p.index == 0); CHECK(p.own == 1); CHECK(p.in_radius); CHECK(calls == 1);
    auto q = stm_michel_pick_wire({{30.0, 40.0}}, 10.0, true, true, [&](size_t) { ++calls; return 1; });
    CHECK(q.index == 0); CHECK(q.own == 0); CHECK_FALSE(q.in_radius); CHECK(calls == 1);   // not in radius: no own test
}

TEST_CASE("stm_michel_pick_wire: the far face listed first loses to the covered near one")
{
    std::vector<size_t> asked;
    auto p = stm_michel_pick_wire({{150.0, 160.0}, {3.0, 40.0}}, 10.0, true, true,
                                  [&](size_t i) { asked.push_back(i); return i == 1 ? 4 : 0; });
    CHECK(p.index == 1); CHECK(p.own == 4); CHECK(p.in_radius);
    CHECK(asked == std::vector<size_t>{1});                  // the out-of-radius one is never tested
}

TEST_CASE("stm_michel_pick_wire: coverage decides among in-radius incarnations, then distance")
{
    // both in radius, the nearer one NOT covered: the covered one wins
    auto p = stm_michel_pick_wire({{2.0, -1.0}, {8.0, -1.0}}, 10.0, true, false, [](size_t i) { return i == 1 ? 1 : 0; });
    CHECK(p.index == 1); CHECK(p.own == 1);
    // both covered: the first in order, and the second is never asked
    int calls = 0;
    auto q = stm_michel_pick_wire({{8.0, -1.0}, {2.0, -1.0}}, 10.0, true, false, [&](size_t) { ++calls; return 1; });
    CHECK(q.index == 0); CHECK(calls == 1);
    // none covered: the nearest in-radius one, own 0
    auto r = stm_michel_pick_wire({{8.0, -1.0}, {2.0, -1.0}, {50.0, -1.0}}, 10.0, true, false, [](size_t) { return 0; });
    CHECK(r.index == 1); CHECK(r.own == 0); CHECK(r.in_radius);
}

TEST_CASE("stm_michel_pick_wire: disabled centres, unprojectable distances, ties")
{
    // the stop distance is ignored when only the control is on
    auto p = stm_michel_pick_wire({{1.0, 90.0}, {80.0, 5.0}}, 10.0, false, true, [](size_t) { return 0; });
    CHECK(p.index == 1); CHECK(p.in_radius);
    // nothing in radius: the nearest overall
    auto q = stm_michel_pick_wire({{300.0, 280.0}, {40.0, 60.0}}, 10.0, true, true, [](size_t) { return 1; });
    CHECK(q.index == 1); CHECK_FALSE(q.in_radius); CHECK(q.own == 0);
    // nothing projectable: index 0
    auto r = stm_michel_pick_wire({{-1.0, -1.0}, {-1.0, -1.0}}, 10.0, true, true, [](size_t) { return 1; });
    CHECK(r.index == 0); CHECK_FALSE(r.in_radius);
    // a tie keeps the earlier index
    auto s = stm_michel_pick_wire({{4.0, -1.0}, {4.0, -1.0}}, 10.0, true, true, [](size_t) { return 0; });
    CHECK(s.index == 0);
    // an empty own_of is allowed: coverage is never decided, distance picks
    auto t = stm_michel_pick_wire({{9.0, -1.0}, {1.0, -1.0}}, 10.0, true, true, nullptr);
    CHECK(t.index == 1); CHECK(t.own == 0);
}

TEST_CASE("stm_michel_combine_planes: symmetric input is the weighted mean")
{
    const std::array<double, 3> w{{0.25, 0.25, 1.0}};
    int dropped = 99;
    // equal planes: the mean is the common value, no plane dropped
    CHECK(stm_michel_combine_planes({{100.0, 100.0, 100.0}}, w, 0.04, &dropped) == doctest::Approx(100.0));
    CHECK(dropped == -1);
    // within the 4 % window: weighted mean (0.25*100 + 0.25*102 + 1.0*101) / 1.5
    CHECK(stm_michel_combine_planes({{100.0, 102.0, 101.0}}, w, 0.04, &dropped) ==
          doctest::Approx((0.25 * 100 + 0.25 * 102 + 1.0 * 101) / 1.5));
    CHECK(dropped == -1);
}

TEST_CASE("stm_michel_combine_planes: an asymmetric pair drops the largest plane")
{
    const std::array<double, 3> w{{0.25, 0.25, 1.0}};
    int dropped = -5;
    // W is 50 % above the others: |med - max| / (med + max) = 50/250 = 0.2 > 0.04
    // -> the (median, minimum) pair = (U, V) with weights 0.25/0.25
    const double r = stm_michel_combine_planes({{100.0, 100.0, 150.0}}, w, 0.04, &dropped);
    CHECK(r == doctest::Approx(100.0));
    CHECK(dropped == 2);
    // U hot: dropped plane is 0, the pair is (V, W)
    dropped = -5;
    const double r2 = stm_michel_combine_planes({{200.0, 100.0, 110.0}}, w, 0.04, &dropped);
    CHECK(dropped == 0);
    CHECK(r2 == doctest::Approx((0.25 * 100 + 1.0 * 110) / 1.25));
    // the same input with the switch off (1.0) keeps all three
    dropped = -5;
    CHECK(stm_michel_combine_planes({{200.0, 100.0, 110.0}}, w, 1.0, &dropped) ==
          doctest::Approx((0.25 * 200 + 0.25 * 100 + 1.0 * 110) / 1.5));
    CHECK(dropped == -1);
    // nullptr for dropped_plane is allowed
    CHECK(stm_michel_combine_planes({{100.0, 100.0, 150.0}}, w, 0.04, nullptr) == doctest::Approx(100.0));
}

TEST_CASE("stm_michel_combine_planes: signed input and degenerate weights")
{
    const std::array<double, 3> w{{0.25, 0.25, 1.0}};
    int dropped = -5;
    // med + max <= 0: the asymmetry is never evaluated, the 3-plane mean stands
    const double r = stm_michel_combine_planes({{-10.0, -20.0, -5.0}}, w, 0.04, &dropped);
    CHECK(r == doctest::Approx((0.25 * -10 + 0.25 * -20 + 1.0 * -5) / 1.5));
    CHECK(dropped == -1);
    // a negative plane with two positive ones: switch fires on the positive pair
    dropped = -5;
    const double r2 = stm_michel_combine_planes({{-10.0, 100.0, 150.0}}, w, 0.04, &dropped);
    CHECK(dropped == 2);
    CHECK(r2 == doctest::Approx((0.25 * 100 + 0.25 * -10) / 0.5));
    // all-equal planes take the (0,1,2) fallback ordering and never drop
    dropped = -5;
    CHECK(stm_michel_combine_planes({{7.0, 7.0, 7.0}}, w, 0.0, &dropped) == doctest::Approx(7.0));
    CHECK(dropped == -1);
    // all-zero weights: 0, nothing dropped
    dropped = -5;
    CHECK(stm_michel_combine_planes({{100.0, 100.0, 150.0}}, {{0.0, 0.0, 0.0}}, 0.04, &dropped) == 0.0);
    CHECK(dropped == -1);
}

TEST_CASE("masked_response_prediction: scale * R * (pos masked), short vectors read as zero")
{
    // R (3 rows x 2 cols): r0 = [1, 0], r1 = [0.5, 2], r2 = [0, 1]
    Eigen::SparseMatrix<double> R(3, 2);
    R.insert(0, 0) = 1.0;
    R.insert(1, 0) = 0.5;
    R.insert(1, 1) = 2.0;
    R.insert(2, 1) = 1.0;
    R.makeCompressed();
    Eigen::VectorXd pos(2);
    pos << 10.0, 100.0;
    const std::vector<double> scale{2.0, 1.0, 1.0};

    // column 0 only: R * [10, 0] = [10, 5, 0], scaled [20, 5, 0]
    auto p = TrackFitting::masked_response_prediction(R, pos, {1, 0}, scale);
    REQUIRE(p.size() == 3);
    CHECK(p(0) == doctest::Approx(20.0));
    CHECK(p(1) == doctest::Approx(5.0));
    CHECK(p(2) == doctest::Approx(0.0));

    // both columns: R * [10, 100] = [10, 205, 100], scaled [20, 205, 100]
    auto q = TrackFitting::masked_response_prediction(R, pos, {1, 1}, scale);
    CHECK(q(0) == doctest::Approx(20.0));
    CHECK(q(1) == doctest::Approx(205.0));
    CHECK(q(2) == doctest::Approx(100.0));

    // a mask shorter than pos reads as 0 beyond its end == {1, 0}
    auto s = TrackFitting::masked_response_prediction(R, pos, {1}, scale);
    CHECK(s(1) == doctest::Approx(5.0));
    CHECK(s(2) == doctest::Approx(0.0));

    // a scale shorter than the rows zeroes the rows it does not cover
    auto t = TrackFitting::masked_response_prediction(R, pos, {1, 1}, {2.0});
    CHECK(t(0) == doctest::Approx(20.0));
    CHECK(t(1) == doctest::Approx(0.0));
    CHECK(t(2) == doctest::Approx(0.0));

    // an empty mask predicts nothing
    auto z = TrackFitting::masked_response_prediction(R, pos, {}, scale);
    CHECK(z.norm() == 0.0);

    // a pos of the wrong length predicts nothing rather than reading past R
    Eigen::VectorXd bad(3);
    bad << 1, 2, 3;
    auto b = TrackFitting::masked_response_prediction(R, bad, {1, 1, 1}, scale);
    REQUIRE(b.size() == 3);
    CHECK(b.norm() == 0.0);
}

TEST_CASE("TrackFitting keep_dqdx_response defaults to 0 and round-trips")
{
    TrackFitting tf;
    CHECK(tf.get_parameter("keep_dqdx_response") == 0.0);
    CHECK(tf.get_dqdx_response(nullptr) == nullptr);
    tf.set_parameter("keep_dqdx_response", 1.0);
    CHECK(tf.get_parameter("keep_dqdx_response") == 1.0);
    tf.clear_dqdx_responses();
    CHECK(tf.get_dqdx_response(nullptr) == nullptr);
}
