// doc pr/109 §8 -- the predicted 2-D charge a cluster's OWN dQ/dx fit produced
// must survive to the writers.
//
// Symptom this pins: T_proj_data rows tagged with the main cluster carried
// pred_charge = 0 on cells sitting right on that cluster's own trajectory,
// concentrated near the neutrino vertex (uBooNE 5384-6528: cid 19, 2141 cells,
// Σpred = 0; SBND 46363: cid 19 kept 44 % of its own predicted charge).
//
// Two mechanisms, one per TEST_CASE below:
//
//   * The MERGE.  assemble_fitted_charge_2d() flattens every per-cluster
//     snapshot into one cell map, last-writer-wins in ascending ident order.
//     A satellite cluster holds a main-cluster cell inside its padded bounding
//     box and predicts 0 there; having the higher ident it writes last and
//     wins.  The merged map is therefore NOT "this cluster's prediction" --
//     the writers must read get_cluster_fitted_charge_2d() instead.
//
//   * The IDENT COLLISION.  The snapshot store used to be a
//     std::map<Facade::Cluster*, ..., PR::ClusterPtrCmp>, i.e. keyed by ident,
//     so two live clusters sharing an ident compared EQUAL and the earlier
//     snapshot was discarded whole.  Reverting the store to that map makes the
//     second case fail (one entry instead of two, and the first cluster's
//     cells gone).
//
// Reverting the writers to the merged map does not break these cases directly
// -- it breaks the invariant case 1 states -- so read case 1 as the contract
// the writers rely on.

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

    // One live U-plane cell at (apa 0, face 0, wire, tick), charge q.
    void add_cell(ChargeMap& m, int wire, int tick, double q)
    {
        const int apa = 0, face = 0, channel = wire;
        TrackFitting::CoordReadout key(apa, tick, channel);
        TrackFitting::ChargeMeasurement meas(q, 100.0, 1);
        std::set<TrackFitting::Coord2D> coords;
        coords.insert(TrackFitting::Coord2D(apa, face, tick, wire, channel, kUlayer));
        m[key] = std::make_pair(meas, coords);
    }

    // Un-whitened prediction that fill_fitted_charge_2d() will recover:
    // it multiplies by sqrt(charge_err^2 + (charge*rel)^2 + add^2).
    double whiten(double pred, double q, double qe, double rel, double add)
    {
        return pred / std::sqrt(qe * qe + (q * rel) * (q * rel) + add * add);
    }

    constexpr double REL = 0.075, ADD = 0.0;   // induction

    void fit_one(TrackFitting& tf, Cluster* cl, const ChargeMap& mu, const Eigen::VectorXd& pred_u)
    {
        ChargeMap empty;
        Eigen::VectorXd none(0);
        tf.set_cluster_filter(cl);
        tf.fill_fitted_charge_2d(mu, empty, empty, pred_u, none, none, REL, 0.05, ADD, 300.0);
    }

    double merged_pred(const TrackFitting& tf, int wire, int tick)
    {
        const auto& m = tf.get_fitted_charge_2d();
        auto it = m.find(TrackFitting::APAFacePlane{0, 0, 0});
        REQUIRE(it != m.end());
        auto jt = it->second.find(TrackFitting::WireTime{wire, tick});
        REQUIRE(jt != it->second.end());
        return jt->second.pred_charge;
    }

}  // namespace

TEST_CASE("pr109: the merged map is not the cluster's own prediction; the snapshots are")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();
    Cluster& main_cl = g.make_child();   main_cl.set_cluster_id(19);
    Cluster& sat_cl  = g.make_child();   sat_cl.set_cluster_id(75);

    // The shared cell (wire 100) plus one cell each side that only one of the
    // two clusters carries.
    const double q = 20000.0, qe = 100.0;
    ChargeMap main_map, sat_map;
    add_cell(main_map, 99, 40, q);
    add_cell(main_map, 100, 40, q);      // shared
    add_cell(sat_map, 100, 40, q);       // shared
    add_cell(sat_map, 101, 40, q);

    // The main cluster's fit explains the shared cell; the satellite's does
    // not reach it (its trajectory is elsewhere), so it predicts 0 there.
    Eigen::VectorXd main_pred(2), sat_pred(2);
    main_pred << whiten(18000.0, q, qe, REL, ADD), whiten(19000.0, q, qe, REL, ADD);
    sat_pred  << 0.0, whiten(17000.0, q, qe, REL, ADD);

    TrackFitting tf;
    fit_one(tf, &main_cl, main_map, main_pred);
    fit_one(tf, &sat_cl, sat_map, sat_pred);

    // Both snapshots are held, in capture order, each with its OWN prediction.
    const auto& snaps = tf.get_cluster_fitted_charge_2d();
    REQUIRE(snaps.size() == 2);
    CHECK(snaps[0].ident == 19);
    CHECK(snaps[1].ident == 75);
    CHECK(snaps[0].cells.at({0, 0, 0}).at({100, 40}).pred_charge == doctest::Approx(19000.0));
    CHECK(snaps[1].cells.at({0, 0, 0}).at({100, 40}).pred_charge == doctest::Approx(0.0));
    // ... and the cells only one of them carries stay with that one.
    CHECK(snaps[0].cells.at({0, 0, 0}).count({101, 40}) == 0);
    CHECK(snaps[1].cells.at({0, 0, 0}).count({99, 40}) == 0);

    // The merge loses the main cluster's answer on the shared cell: the
    // higher ident writes last.  This is the defect the writers must not read.
    tf.assemble_fitted_charge_2d();
    CHECK(merged_pred(tf, 100, 40) == doctest::Approx(0.0));
    CHECK(merged_pred(tf, 99, 40) == doctest::Approx(18000.0));   // unshared: intact
    CHECK(merged_pred(tf, 101, 40) == doctest::Approx(17000.0));
}

TEST_CASE("pr109: two live clusters sharing an ident both keep their snapshot")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();
    Cluster& a = g.make_child();  a.set_cluster_id(19);
    Cluster& b = g.make_child();  b.set_cluster_id(19);   // the 6528 collision

    const double q = 20000.0, qe = 100.0;
    ChargeMap ma, mb;
    add_cell(ma, 10, 40, q);
    add_cell(mb, 20, 40, q);

    Eigen::VectorXd pa(1), pb(1);
    pa << whiten(11000.0, q, qe, REL, ADD);
    pb << whiten(22000.0, q, qe, REL, ADD);

    TrackFitting tf;
    fit_one(tf, &a, ma, pa);
    fit_one(tf, &b, mb, pb);

    // The ident-keyed map kept ONE of these and threw the other away.
    REQUIRE(tf.get_cluster_fitted_charge_2d().size() == 2);

    // Re-fitting the SAME cluster still replaces in place -- latest fit wins,
    // no third entry.
    Eigen::VectorXd pa2(1);
    pa2 << whiten(12000.0, q, qe, REL, ADD);
    fit_one(tf, &a, ma, pa2);
    REQUIRE(tf.get_cluster_fitted_charge_2d().size() == 2);
    CHECK(tf.get_cluster_fitted_charge_2d()[0].cells.at({0, 0, 0}).at({10, 40}).pred_charge
          == doctest::Approx(12000.0));

    tf.assemble_fitted_charge_2d();
    CHECK(merged_pred(tf, 10, 40) == doctest::Approx(12000.0));
    CHECK(merged_pred(tf, 20, 40) == doctest::Approx(22000.0));
}
