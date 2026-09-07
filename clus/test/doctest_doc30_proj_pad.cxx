// doc 30 round 3 -- TrackFitting::Parameters::proj_pad_wire / proj_pad_time,
// the "predicted cells and their neighbourhood" filter on the fitted 2-D charge
// DISPLAY product.
//
// fill_fitted_charge_2d() stores one FittedCharge2D per cell of the fit's whole
// charge map, whether the fit predicts anything there or not.  Measured on PDHD
// 029107/18 that is ~362 k cells per call, 103 calls per event, ~2 % of them
// carrying a prediction.  The knob keeps only cells within proj_pad_wire wires
// and proj_pad_time slices of a predicted cell.
//
// Four contracts pinned here:
//   1. OFF (the shipped default, -1) stores every cell -- the legacy product.
//   2. pad 0 stores exactly the predicted cells.
//   3. pad N stores the wire band around them, and the pad is per (apa, face,
//      plane): a cell in another plane is never kept by a U-plane seed.
//   4. The SEED is the RAW prediction (R*pos_3D), NOT pred_charge.  A dead or
//      below-threshold channel crossed by the track has pred_charge == 0 by
//      construction; seeding on pred_charge would punch holes in the display
//      exactly where the fit is most interesting.
//
// The time pad is expressed in SLICES and converted with
// Grouping::get_nticks_per_slice(); with no grouping wired (as here) the
// fallback is 1 tick per slice, which case 5 pins.

#include "WireCellUtil/doctest.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"

#include <Eigen/Dense>
#include <map>
#include <set>

using namespace WireCell;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;
using WireCell::Clus::TrackFitting;

namespace {

    using ChargeMap = std::map<TrackFitting::CoordReadout,
                               std::pair<TrackFitting::ChargeMeasurement,
                                         std::set<TrackFitting::Coord2D>>>;

    // One U-plane cell at (apa 0, face 0, wire, tick).  flag 1 = live.
    void add_cell(ChargeMap& m, int wire, int tick, double q, int flag = 1)
    {
        const int apa = 0, face = 0, channel = wire;
        TrackFitting::CoordReadout key(apa, tick, channel);
        TrackFitting::ChargeMeasurement meas(q, 100.0, flag);
        std::set<TrackFitting::Coord2D> coords;
        coords.insert(TrackFitting::Coord2D(apa, face, tick, wire, channel, kUlayer));
        m[key] = std::make_pair(meas, coords);
    }

    constexpr double REL = 0.075, ADD = 0.0;   // induction

    // Run the map flavour with map_U only.  The CoordReadout order is
    // (apa, time, channel), so pred_u's index i is the i-th cell in that order.
    void fill(TrackFitting& tf, const ChargeMap& mu, const Eigen::VectorXd& pred_u)
    {
        ChargeMap empty;
        Eigen::VectorXd none(0);
        tf.fill_fitted_charge_2d(mu, empty, empty, pred_u, none, none, REL, 0.05, ADD, 300.0);
    }

    size_t ncells(const TrackFitting& tf)
    {
        size_t n = 0;
        for (const auto& [afp, wt_map] : tf.get_fitted_charge_2d()) { (void)afp; n += wt_map.size(); }
        return n;
    }

    bool has(const TrackFitting& tf, int wire, int tick)
    {
        const auto& m = tf.get_fitted_charge_2d();
        auto it = m.find(TrackFitting::APAFacePlane{0, 0, 0});
        if (it == m.end()) return false;
        return it->second.count(TrackFitting::WireTime{wire, tick}) != 0;
    }

}  // namespace

TEST_CASE("doc30 r3: proj_pad_wire defaults OFF and stores every cell")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();
    (void)g;

    ChargeMap mu;
    for (int w = 10; w <= 20; ++w) add_cell(mu, w, 100, 1000.0);
    Eigen::VectorXd pred_u = Eigen::VectorXd::Zero(11);
    pred_u(5) = 0.5;                     // only wire 15 is predicted

    TrackFitting tf;
    CHECK(tf.get_parameter("proj_pad_wire") == doctest::Approx(-1.0));
    CHECK(tf.get_parameter("proj_pad_time") == doctest::Approx(0.0));

    fill(tf, mu, pred_u);
    CHECK(ncells(tf) == 11);             // the legacy product: everything
    CHECK(has(tf, 10, 100));
    CHECK(has(tf, 20, 100));
}

TEST_CASE("doc30 r3: pad 0 keeps exactly the predicted cells, pad N the wire band")
{
    ChargeMap mu;
    for (int w = 10; w <= 20; ++w) add_cell(mu, w, 100, 1000.0);
    Eigen::VectorXd pred_u = Eigen::VectorXd::Zero(11);
    pred_u(5) = 0.5;                     // wire 15

    TrackFitting tf;

    tf.set_parameter("proj_pad_wire", 0.0);
    fill(tf, mu, pred_u);
    CHECK(ncells(tf) == 1);
    CHECK(has(tf, 15, 100));
    CHECK_FALSE(has(tf, 14, 100));

    tf.set_parameter("proj_pad_wire", 2.0);
    fill(tf, mu, pred_u);
    CHECK(ncells(tf) == 5);              // wires 13..17
    CHECK(has(tf, 13, 100));
    CHECK(has(tf, 17, 100));
    CHECK_FALSE(has(tf, 12, 100));
    CHECK_FALSE(has(tf, 18, 100));

    // Two separated seeds give two bands, not their convex hull.
    pred_u(5) = 0.5; pred_u(10) = 0.5;   // wires 15 and 20
    tf.set_parameter("proj_pad_wire", 1.0);
    fill(tf, mu, pred_u);
    CHECK(ncells(tf) == 5);              // 14..16 and 19..20 (21 is not in the map)
    CHECK_FALSE(has(tf, 17, 100));
    CHECK(has(tf, 19, 100));
}

TEST_CASE("doc30 r3: the seed is the RAW prediction, so a dead channel on the track survives")
{
    ChargeMap mu;
    for (int w = 10; w <= 20; ++w) add_cell(mu, w, 100, 1000.0);
    // Wire 15 is DEAD: pred_charge is forced to 0 by the charge>0 && flag!=0
    // gate, but the fit's response column there is nonzero.
    add_cell(mu, 15, 100, 0.0, 0);
    Eigen::VectorXd pred_u = Eigen::VectorXd::Zero(11);
    pred_u(5) = 0.5;

    TrackFitting tf;
    tf.set_parameter("proj_pad_wire", 0.0);
    fill(tf, mu, pred_u);

    REQUIRE(ncells(tf) == 1);
    CHECK(has(tf, 15, 100));             // kept, even though pred_charge == 0
    const auto& m = tf.get_fitted_charge_2d();
    const auto& fc = m.at(TrackFitting::APAFacePlane{0, 0, 0}).at(TrackFitting::WireTime{15, 100});
    CHECK(fc.pred_charge == doctest::Approx(0.0));
    CHECK(fc.flag == 0);
}

TEST_CASE("doc30 r3: the pad is per plane -- a U seed never keeps a V cell")
{
    // Same wire/tick in two planes.  map_V's own prediction is zero, so nothing
    // in V may survive however wide the U pad is.
    ChargeMap mu, mv;
    for (int w = 10; w <= 20; ++w) add_cell(mu, w, 100, 1000.0);
    for (int w = 10; w <= 20; ++w) add_cell(mv, w, 100, 1000.0);
    Eigen::VectorXd pred_u = Eigen::VectorXd::Zero(11);
    pred_u(5) = 0.5;
    Eigen::VectorXd pred_v = Eigen::VectorXd::Zero(11);
    Eigen::VectorXd none(0);
    ChargeMap empty;

    TrackFitting tf;
    tf.set_parameter("proj_pad_wire", 3.0);
    tf.fill_fitted_charge_2d(mu, mv, empty, pred_u, pred_v, none, REL, 0.05, ADD, 300.0);

    const auto& m = tf.get_fitted_charge_2d();
    CHECK(m.count(TrackFitting::APAFacePlane{0, 0, 0}) == 1);   // U band survives
    CHECK(m.count(TrackFitting::APAFacePlane{0, 0, 1}) == 0);   // V is empty
    CHECK(ncells(tf) == 7);                                      // wires 12..18
}

TEST_CASE("doc30 r3: the time pad steps by the slice width (1 tick with no grouping)")
{
    ChargeMap mu;
    for (int t = 98; t <= 102; ++t)
        for (int w = 14; w <= 16; ++w) add_cell(mu, w, t, 1000.0);
    // CoordReadout order is (apa, time, channel): 3 cells per tick, ticks
    // 98,99,100,101,102 -> the wire-15 cell of tick 100 is index 2*3 + 1 = 7.
    Eigen::VectorXd pred_u = Eigen::VectorXd::Zero(15);
    pred_u(7) = 0.5;

    TrackFitting tf;
    tf.set_parameter("proj_pad_wire", 0.0);
    tf.set_parameter("proj_pad_time", 0.0);
    fill(tf, mu, pred_u);
    CHECK(ncells(tf) == 1);

    tf.set_parameter("proj_pad_time", 1.0);
    fill(tf, mu, pred_u);
    CHECK(ncells(tf) == 3);              // ticks 99,100,101 at wire 15
    CHECK(has(tf, 15, 99));
    CHECK(has(tf, 15, 101));
    CHECK_FALSE(has(tf, 15, 98));
    CHECK_FALSE(has(tf, 14, 100));       // wire pad is still 0

    tf.set_parameter("proj_pad_wire", 1.0);
    fill(tf, mu, pred_u);
    CHECK(ncells(tf) == 9);              // 3 wires x 3 ticks
}
