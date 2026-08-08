// Regression pins for the doc pr/49 fit_blob_coverage knob's coverage
// predicate, TrackFitting::is_cell_covered_by_own_blobs (18255-57441
// V-plane projection ghost: a physically unrelated cluster's real charge
// aliases with the fitted track in one projection view only and detours the
// fit; the knob drops live candidate cells not covered by the fitted
// cluster's OWN blobs).
//
// What these cases pin, and the plausible wrong fixes they catch:
//
//   * Wire ranges are half-open [min, max) (CLAUDE.md M7).  With tolerance
//     the accepted band is [min - tol, max + tol) -- NOT <= max + tol.  The
//     "two wires out with tol=1" cases fail against the inclusive mistake.
//
//   * The time test is an INTERVAL search over time_blob_map()'s
//     slice_index_min keys, not an exact-key lookup.  The candidate cells of
//     form_point_association's Steiner and fallback branches carry
//     floor-quantized ticks that need not equal a blob's own key
//     (doc pr/49 "time-key alignment" blocker), so a mid-interval tick MUST
//     still be found.  The "tick inside the slice but off the key" case
//     fails against a naive tbm.find(time)/at(time).
//
//   * The 57441 contamination sits ONE cell outside the fitted cluster's
//     own footprint, so tol=0 must exclude the adjacent cell and tol=1 must
//     re-admit it -- the "one cell out" cases pin both directions.
//
// Reverting the knob (or the predicate) must make this file fail to compile
// or fail these cases.

#include "WireCellUtil/doctest.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"

#include <vector>

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;

namespace {

    using fa_float_t = double;
    using fa_int_t = int;

    // Every blob here lives in (apa 0, face 0) so the predicate's map walk
    // and the probe share a key path.
    int test_wpid() { return WirePlaneId(kUlayer, 0, 0).ident(); }

    // A blob holding exactly one point, with the scalar PC Blob::fill_cache
    // needs.  Wire ranges are given as the HALF-OPEN [min, max) the facade
    // stores.  Deliberately duplicated from doctest_steiner_terminal_filter
    // rather than factored into a shared fixture (CLAUDE.md M10).
    void add_blob(Cluster& cl,
                  int u_lo, int u_hi, int v_lo, int v_hi, int w_lo, int w_hi,
                  int slice_lo, int slice_hi)
    {
        Blob& bl = cl.make_child();

        Dataset ds;
        ds.add("x", Array(std::vector<double>{0.0}));
        ds.add("y", Array(std::vector<double>{0.0}));
        ds.add("z", Array(std::vector<double>{0.0}));
        ds.add("uwire_index", Array(std::vector<int>{u_lo}));
        ds.add("vwire_index", Array(std::vector<int>{v_lo}));
        ds.add("wwire_index", Array(std::vector<int>{w_lo}));
        bl.value().local_pcs()["3d"] = ds;

        bl.value().local_pcs()["scalar"] = Dataset({
            {"charge", Array({(fa_float_t) 1.0})},
            {"center_x", Array({(fa_float_t) 0.0})},
            {"center_y", Array({(fa_float_t) 0.0})},
            {"center_z", Array({(fa_float_t) 0.0})},
            {"wpid", Array({(fa_int_t) test_wpid()})},
            {"npoints", Array({(fa_int_t) 1})},
            {"slice_index_min", Array({(fa_int_t) slice_lo})},
            {"slice_index_max", Array({(fa_int_t) slice_hi})},
            {"u_wire_index_min", Array({(fa_int_t) u_lo})},
            {"u_wire_index_max", Array({(fa_int_t) u_hi})},
            {"v_wire_index_min", Array({(fa_int_t) v_lo})},
            {"v_wire_index_max", Array({(fa_int_t) v_hi})},
            {"w_wire_index_min", Array({(fa_int_t) w_lo})},
            {"w_wire_index_max", Array({(fa_int_t) w_hi})},
            {"max_wire_interval", Array({(fa_int_t) 1})},
            {"min_wire_interval", Array({(fa_int_t) 1})},
            {"max_wire_type", Array({(fa_int_t) 0})},
            {"min_wire_type", Array({(fa_int_t) 0})},
        });
    }

    constexpr int NT = 4;  // nticks per slice (SBND)

}  // namespace

TEST_CASE("pr49 coverage predicate: half-open wire band and tolerance, per plane")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();
    Cluster& cl = g.make_child();
    // u wires 10..11, v wires 20..24, w wires 30..39; one slice [0, 4).
    add_blob(cl, 10, 12, 20, 25, 30, 40, 0, NT);

    Clus::TrackFitting tf;

    // U plane, tol=0: exactly the half-open band.
    CHECK(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 10, 0, 0, NT));
    CHECK(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 11, 0, 0, NT));
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 12, 0, 0, NT));  // one past max
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 9, 0, 0, NT));   // one below min

    // U plane, tol=1: one cell out is re-admitted, two is not, both sides.
    CHECK(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 12, 0, 1, NT));
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 13, 0, 1, NT));  // <= max + tol mistake admits this
    CHECK(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 9, 0, 1, NT));
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 8, 0, 1, NT));

    // Plane dispatch: the same wire is judged against its OWN plane's band.
    CHECK(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 1, 24, 0, 0, NT));        // v: in [20,25)
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 1, 25, 0, 0, NT));
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 24, 0, 0, NT));  // u band does not cover 24
    CHECK(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 2, 39, 0, 0, NT));        // w: in [30,40)
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 2, 40, 0, 0, NT));
}

TEST_CASE("pr49 coverage predicate: time is an interval search, not an exact key")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();
    Cluster& cl = g.make_child();
    // One slice [8, 12): time_blob_map key is 8.
    add_blob(cl, 10, 12, 20, 25, 30, 40, 8, 8 + NT);

    Clus::TrackFitting tf;

    // Ticks inside the slice interval are covered even when they are not the
    // map key -- this is the Steiner/fallback misaligned-tick case.
    CHECK(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 10, 8, 0, NT));   // the key itself
    CHECK(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 10, 10, 0, NT));  // mid-slice tick, off the key
    CHECK(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 10, 11, 0, NT));  // last covered tick
    // Half-open in time as well.
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 10, 12, 0, NT));
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 10, 7, 0, NT));

    // tol=1 widens the time window by one slice (NT ticks) each way.
    CHECK(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 10, 12, 1, NT));
    CHECK(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 10, 15, 1, NT));
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 10, 16, 1, NT));
    CHECK(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 10, 4, 1, NT));
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 10, 3, 1, NT));
}

TEST_CASE("pr49 coverage predicate: only the cluster's own (apa,face) blobs count")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();
    Cluster& cl = g.make_child();
    add_blob(cl, 10, 12, 20, 25, 30, 40, 0, NT);

    Clus::TrackFitting tf;

    // Same cell, wrong apa or face: no coverage.
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&cl, 1, 0, 0, 10, 0, 0, NT));
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&cl, 0, 1, 0, 10, 0, 0, NT));

    // A cluster with no blobs covers nothing (empty time_blob_map).
    Cluster& empty = g.make_child();
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&empty, 0, 0, 0, 10, 0, 0, NT));

    // The 57441 shape in miniature: a second cluster's blob adjacent in time
    // (next slice) and wire does NOT make the first cluster cover its cells.
    Cluster& other = g.make_child();
    add_blob(other, 12, 14, 25, 30, 40, 50, NT, 2 * NT);
    CHECK(tf.is_cell_covered_by_own_blobs(&other, 0, 0, 0, 12, NT, 0, NT));      // other covers its own cell
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&cl, 0, 0, 0, 12, NT, 0, NT));   // cl does not
}

TEST_CASE("pr49 coverage predicate: foreign coverage decides the deweight, nobody's cells stay")
{
    // The knob's deweight condition is !own && foreign(3D-far).  Cell classes:
    //   (a) own territory        -> full weight (own test true, never reaches
    //       foreign);
    //   (b) 3D-far foreign territory -> DEWEIGHTED (the 57441 ghost cells:
    //       cluster 13's own tiled footprint, claimant 163 cm away in 3D);
    //   (c) covered by NOBODY    -> full weight (own-track charge spilling
    //       just past the tiled envelope, or dead-channel single-view charge
    //       with no 3D image -- deweighting these moved 47/48 nueCC events
    //       in the strict-coverage first cut of this round);
    //   (d) 3D-NEAR foreign territory -> full weight (a genuinely touching
    //       or crossing neighbor shares the projection legitimately).
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();
    Cluster& mine = g.make_child();
    add_blob(mine, 10, 12, 20, 25, 30, 40, 0, NT);
    Cluster& foreign = g.make_child();
    add_blob(foreign, 14, 16, 40, 45, 50, 60, 0, NT);

    Clus::TrackFitting tf;
    const WireCell::Point p0(0, 0, 0);   // fixture points all sit at the origin
    const double NOGATE = -1;            // ghost_dis <= 0 disables the 3D gate

    // (b): the foreign cluster's own cell is foreign-covered from mine's view.
    CHECK(tf.is_cell_covered_by_foreign_blobs(&g, &mine, p0, NOGATE, 0, 0, 0, 14, 0, 0, NT));
    // ...and symmetric.
    CHECK(tf.is_cell_covered_by_foreign_blobs(&g, &foreign, p0, NOGATE, 0, 0, 0, 10, 0, 0, NT));

    // (a): a cluster is never "foreign" to itself -- mine's own cell is not
    // foreign-covered (only mine covers it, and mine is skipped).
    CHECK_FALSE(tf.is_cell_covered_by_foreign_blobs(&g, &mine, p0, NOGATE, 0, 0, 0, 10, 0, 0, NT));

    // (c): the gap wire 13 (between mine's [10,12) and foreign's [14,16)) is
    // covered by nobody -- not own, not foreign.  The knob keeps it.
    CHECK_FALSE(tf.is_cell_covered_by_own_blobs(&mine, 0, 0, 0, 13, 0, 0, NT));
    CHECK_FALSE(tf.is_cell_covered_by_foreign_blobs(&g, &mine, p0, NOGATE, 0, 0, 0, 13, 0, 0, NT));

    // (d) the 3D gate: with the fit point ON the foreign cluster's points
    // (distance 0 <= ghost_dis) the claim is a nearby neighbor's -- exempt;
    // with the fit point far away (distance > ghost_dis) it is a genuine
    // projection ghost -- deweighted.
    const double GDIS = 150;  // 15 cm in WCT internal units
    CHECK_FALSE(tf.is_cell_covered_by_foreign_blobs(&g, &mine, p0, GDIS, 0, 0, 0, 14, 0, 0, NT));
    const WireCell::Point pfar(1000, 0, 0);
    CHECK(tf.is_cell_covered_by_foreign_blobs(&g, &mine, pfar, GDIS, 0, 0, 0, 14, 0, 0, NT));

    // Null grouping: no foreign coverage, never a crash.
    CHECK_FALSE(tf.is_cell_covered_by_foreign_blobs(nullptr, &mine, p0, NOGATE, 0, 0, 0, 14, 0, 0, NT));
}
