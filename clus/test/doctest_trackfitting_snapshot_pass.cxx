// doc pdvd/42 -- TrackFitting::add_fitted_charge_2d_snapshot() and the
// ClusterFitted2D::pass field.
//
// The STM tagger fits a main cluster up to twice (forward, backward), each
// pass with its own prediction.  Before this round the STM hand-off holder
// carried only the MERGED 2D map, which is last-writer-wins on cells: a later
// cluster's bounding-box cells (pred 0 off its trajectory) overwrote every
// earlier cluster's prediction, so tracking-stm.root's T_proj_data had
// charge_pred > 0 for the last fitted cluster only (PDVD 039252/2: 6923 of
// 177857 cells).  The fix hands the holder one snapshot per pass.  These cases
// pin the additive API the writers rely on:
//   * a snapshot recorded by a real fit keeps pass = -1 (PR-stage semantics,
//     doc pr/109 sec 8) -- reverting the default breaks case 1;
//   * added snapshots keep capture order, ident and pass, and two passes of
//     one cluster stay two entries with their own cells -- reverting to a
//     merge breaks case 2;
//   * adding snapshots leaves the merged map alone -- case 3.

#include "WireCellUtil/doctest.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"

#include <map>

using namespace WireCell;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus;

namespace {
    using Cells = std::map<TrackFitting::APAFacePlane, std::map<TrackFitting::WireTime, TrackFitting::FittedCharge2D>>;
    Cells one_cell(int wire, int time, double q, double pred)
    {
        Cells c;
        TrackFitting::FittedCharge2D fc;
        fc.charge = q; fc.charge_err = 0.1 * q; fc.pred_charge = pred; fc.flag = 1;
        c[{0, 0, 2}][{wire, time}] = fc;
        return c;
    }
}

TEST_CASE("doc pdvd/42 snapshot pass: default is the PR-stage value -1")
{
    TrackFitting::ClusterFitted2D snap;
    CHECK(snap.pass == -1);
    CHECK(snap.ident == -1);
    CHECK(snap.cluster == nullptr);
}

TEST_CASE("doc pdvd/42 snapshot pass: two passes of one cluster stay two entries")
{
    Points::node_t root;
    auto& grouping = *root.value.facade<Facade::Grouping>();
    auto& cluster = grouping.make_child();

    TrackFitting holder;
    CHECK(holder.get_cluster_fitted_charge_2d().empty());

    holder.add_fitted_charge_2d_snapshot(&cluster, 7, 0, one_cell(10, 100, 1000., 900.));
    holder.add_fitted_charge_2d_snapshot(&cluster, 7, 1, one_cell(10, 100, 1000., 0.));

    const auto& snaps = holder.get_cluster_fitted_charge_2d();
    REQUIRE(snaps.size() == 2);
    CHECK(snaps[0].cluster == &cluster);
    CHECK(snaps[0].ident == 7);
    CHECK(snaps[0].pass == 0);
    CHECK(snaps[1].pass == 1);
    // Each pass keeps ITS prediction on the shared cell.
    const auto& c0 = snaps[0].cells.at({0, 0, 2}).at({10, 100});
    const auto& c1 = snaps[1].cells.at({0, 0, 2}).at({10, 100});
    CHECK(c0.pred_charge == doctest::Approx(900.));
    CHECK(c1.pred_charge == doctest::Approx(0.));
    CHECK(c0.charge == doctest::Approx(1000.));
    CHECK(c1.charge == doctest::Approx(1000.));
}

TEST_CASE("doc pdvd/42 snapshot pass: adding snapshots does not touch the merged map")
{
    Points::node_t root;
    auto& grouping = *root.value.facade<Facade::Grouping>();
    auto& cluster = grouping.make_child();

    TrackFitting holder;
    holder.merge_fitted_charge_2d(one_cell(3, 30, 500., 0.));   // the legacy accumulator
    holder.add_fitted_charge_2d_snapshot(&cluster, 3, 0, one_cell(3, 30, 500., 480.));

    const auto& merged = holder.get_fitted_charge_2d();
    REQUIRE(merged.size() == 1);
    CHECK(merged.at({0, 0, 2}).at({3, 30}).pred_charge == doctest::Approx(0.));   // unchanged
    REQUIRE(holder.get_cluster_fitted_charge_2d().size() == 1);
    CHECK(holder.get_cluster_fitted_charge_2d()[0].cells.at({0, 0, 2}).at({3, 30}).pred_charge == doctest::Approx(480.));
}
