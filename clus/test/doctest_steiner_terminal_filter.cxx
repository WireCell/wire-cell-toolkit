// Regression pins for the Steiner terminal-filter fidelity fixes, doc pr/29
// D1 (wire tolerance) and D12 (adjacent-slice step).
//
// Both fixes are new PARAMETERS whose defaults reproduce the historical
// behaviour exactly, so the interesting content of these cases is the
// arithmetic at the boundary -- which is exactly where a plausible-looking
// wrong fix lands and where a diff shows nothing.
//
//   D1.  Toolkit blob wire ranges are half-open [min, max); the WCP prototype's
//        are inclusive [low, high] with high == max - 1 (CLAUDE.md M7).  The
//        prototype's one wire of slack is `index <= high + 1`, i.e.
//        `index < max + 1`.  Writing `index <= max + wire_tol` instead admits
//        one extra wire on the high side only.  The "two wires out" cases below
//        fail against that mistake and pass against the correct band.
//
//   D12. The time_blob_map key is blob->slice_index_min(), documented
//        "unit: tick".  Slice starts are one nticks-per-slice apart, so the
//        historical adjacent-slice step of 1 names no real slice whenever a
//        slice is more than one tick wide (SBND: 4) and the whole
//        flag_nearby_timeslice branch is dead.  The cases below pin BOTH that
//        the legacy step still finds nothing and that the tick-span step finds
//        the neighbour.
//
// Reverting either fix must make these fail:
//   * drop the wire_tol arithmetic in Cluster::check_wire_ranges_match
//     -> "one wire out ... with tol=1" cases fail;
//   * write it as `<= max + wire_tol` -> "two wires above" case fails;
//   * restore the literal {-1, 1} in
//     Cluster::is_point_spatially_related_to_time_blobs
//     -> "adjacent slice is found with the tick-span step" fails.

#include "WireCellUtil/doctest.h"
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

    // Every blob here lives in one (apa, face) so that the reference
    // time_blob_map and the probe point share a key path.
    int test_wpid() { return WirePlaneId(kUlayer, 0, 0).ident(); }

    // A blob holding exactly one point, with the per-point wire-index arrays
    // Cluster::wire_index() reads and the scalar PC Blob::fill_cache needs.
    //
    // Wire ranges are given as the HALF-OPEN [min, max) the facade stores, so a
    // caller writing u_hi = 12 is saying "wires 10 and 11", not "through 12".
    // Deliberately duplicated from the other clus doctests rather than factored
    // into a shared fixture (CLAUDE.md M10).
    void add_blob(Cluster& cl,
                  int u_wire, int v_wire, int w_wire,      // the point's wire indices
                  int u_lo, int u_hi, int v_lo, int v_hi, int w_lo, int w_hi,
                  int slice_lo, int slice_hi)
    {
        Blob& bl = cl.make_child();

        Dataset ds;
        ds.add("x", Array(std::vector<double>{0.0}));
        ds.add("y", Array(std::vector<double>{0.0}));
        ds.add("z", Array(std::vector<double>{0.0}));
        ds.add("uwire_index", Array(std::vector<int>{u_wire}));
        ds.add("vwire_index", Array(std::vector<int>{v_wire}));
        ds.add("wwire_index", Array(std::vector<int>{w_wire}));
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

    // The reference blob every wire case is tested against: u wires 10..11
    // (half-open [10, 12)), v and w wide enough never to be the deciding plane,
    // one slice starting at tick 0.
    const int REF_U_LO = 10, REF_U_HI = 12;
    const int WIDE_LO = 0, WIDE_HI = 1000;

}  // namespace

// ---------------------------------------------------------------------------
// D1 -- the wire tolerance, and the half-open arithmetic at its boundary
// ---------------------------------------------------------------------------

// check_wire_ranges_match is private, so every wire case below goes through the
// public entry point with the reference blob in the point's OWN slice and
// flag_nearby_timeslice=false.  That leaves the wire test as the only thing
// that can decide the answer, which is what these cases are about.

TEST_CASE("pr29 D1 wire tolerance defaults to zero (legacy path unchanged)")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();

    Cluster& ref = g.make_child();
    add_blob(ref, 10, 100, 100, REF_U_LO, REF_U_HI, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 0, 4);
    const auto& tbm = ref.time_blob_map();

    Cluster& probe = g.make_child();
    // Point sits on wire 12 -- one PAST the last reference wire (11).
    add_blob(probe, 12, 100, 100, 12, 13, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 0, 4);

    // Default arguments: exact containment, which is what get_extreme_wcps
    // wants and what this code has always done.
    CHECK_FALSE(probe.is_point_spatially_related_to_time_blobs(0, tbm, false));
    CHECK_FALSE(probe.is_point_spatially_related_to_time_blobs(0, tbm, false, 0, 1));

    // Sanity: a point squarely inside is accepted with the same defaults, so
    // the CHECK_FALSEs above are about the boundary and not about the fixture
    // being broken.
    Cluster& inside = g.make_child();
    add_blob(inside, 10, 100, 100, 10, 11, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 0, 4);
    CHECK(inside.is_point_spatially_related_to_time_blobs(0, tbm, false));
}

TEST_CASE("pr29 D1 one wire outside matches with tol=1, two wires does not")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();

    Cluster& ref = g.make_child();
    add_blob(ref, 10, 100, 100, REF_U_LO, REF_U_HI, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 0, 4);
    const auto& tbm = ref.time_blob_map();

    // One wire ABOVE the last reference wire (11).  The prototype's
    // `index <= high + 1` accepts it; under [min,max) that is `< max + 1`.
    {
        Cluster& probe = g.make_child();
        add_blob(probe, 12, 100, 100, 12, 13, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 0, 4);
        CHECK(probe.is_point_spatially_related_to_time_blobs(0, tbm, false, 1, 1));
    }

    // TWO wires above.  This is the case that separates the correct
    // `< max + tol` from the plausible-but-wrong `<= max + tol`: the latter
    // would accept wire 13 and this CHECK_FALSE would fail.
    {
        Cluster& probe = g.make_child();
        add_blob(probe, 13, 100, 100, 13, 14, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 0, 4);
        CHECK_FALSE(probe.is_point_spatially_related_to_time_blobs(0, tbm, false, 1, 1));
    }

    // One wire BELOW the first reference wire (10) -- the symmetric case.
    {
        Cluster& probe = g.make_child();
        add_blob(probe, 9, 100, 100, 9, 10, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 0, 4);
        CHECK(probe.is_point_spatially_related_to_time_blobs(0, tbm, false, 1, 1));
    }

    // Two wires below.
    {
        Cluster& probe = g.make_child();
        add_blob(probe, 8, 100, 100, 8, 9, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 0, 4);
        CHECK_FALSE(probe.is_point_spatially_related_to_time_blobs(0, tbm, false, 1, 1));
    }
}

TEST_CASE("pr29 D1 the tolerance is symmetric, not high-side only")
{
    // A band that is one wire loose above and exact below would still pass both
    // "one wire" cases above.  Pinning the accepted SET catches it: with tol=1
    // the accepted wires are exactly {9, 10, 11, 12}.
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();

    Cluster& ref = g.make_child();
    add_blob(ref, 10, 100, 100, REF_U_LO, REF_U_HI, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 0, 4);
    const auto& tbm = ref.time_blob_map();

    std::vector<int> accepted;
    for (int u = 5; u <= 17; ++u) {
        Cluster& probe = g.make_child();
        add_blob(probe, u, 100, 100, u, u + 1, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 0, 4);
        if (probe.is_point_spatially_related_to_time_blobs(0, tbm, false, 1, 1)) {
            accepted.push_back(u);
        }
    }
    CHECK(accepted == std::vector<int>{9, 10, 11, 12});
}

// ---------------------------------------------------------------------------
// D12 -- the adjacent-slice step is in ticks
// ---------------------------------------------------------------------------

TEST_CASE("pr29 D12 adjacent slice is unreachable with the legacy step of 1")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();

    // Reference blob one SLICE later: 4 ticks per slice, so it starts at tick 4.
    Cluster& ref = g.make_child();
    add_blob(ref, 10, 100, 100, REF_U_LO, REF_U_HI, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 4, 8);

    // Probe point in the slice starting at tick 0, wire-contained in the
    // reference blob so ONLY the time lookup can decide the answer.
    Cluster& probe = g.make_child();
    add_blob(probe, 10, 100, 100, 10, 11, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 0, 4);

    const auto& tbm = ref.time_blob_map();

    // No same-slice reference blob, and the legacy step looks for a slice
    // starting at tick 1 or -1.  Nothing starts there -- this is the dead
    // branch doc pr/29 D12 reports.
    CHECK_FALSE(probe.is_point_spatially_related_to_time_blobs(0, tbm, true));
    CHECK_FALSE(probe.is_point_spatially_related_to_time_blobs(0, tbm, true, 0, 1));

    // And with the fallback switched off entirely, same answer -- i.e. the
    // legacy "true" is indistinguishable from "false", which is the defect.
    CHECK_FALSE(probe.is_point_spatially_related_to_time_blobs(0, tbm, false));
}

TEST_CASE("pr29 D12 adjacent slice is found with the tick-span step")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();

    Cluster& ref = g.make_child();
    add_blob(ref, 10, 100, 100, REF_U_LO, REF_U_HI, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 4, 8);

    Cluster& probe = g.make_child();
    add_blob(probe, 10, 100, 100, 10, 11, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 0, 4);

    const auto& tbm = ref.time_blob_map();

    // Stepping by the tick span reaches it.
    CHECK(probe.is_point_spatially_related_to_time_blobs(0, tbm, true, 0, 4));

    // Still gated on flag_nearby_timeslice: a caller that does not want the
    // fallback (get_extreme_wcps) does not get it, whatever the step.
    CHECK_FALSE(probe.is_point_spatially_related_to_time_blobs(0, tbm, false, 0, 4));
}

TEST_CASE("pr29 D12 the previous slice is reached too, and only the neighbours")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();

    // Reference blob in the slice starting at tick 0.
    Cluster& ref = g.make_child();
    add_blob(ref, 10, 100, 100, REF_U_LO, REF_U_HI, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 0, 4);
    const auto& tbm = ref.time_blob_map();

    // Probe one slice LATER: the -stride branch must find it.
    {
        Cluster& probe = g.make_child();
        add_blob(probe, 10, 100, 100, 10, 11, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 4, 8);
        CHECK(probe.is_point_spatially_related_to_time_blobs(0, tbm, true, 0, 4));
    }

    // Probe TWO slices later: only the immediate neighbours are searched, so
    // this must stay false -- the fix widens the search by one slice, not by
    // an unbounded amount.
    {
        Cluster& probe = g.make_child();
        add_blob(probe, 10, 100, 100, 10, 11, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 8, 12);
        CHECK_FALSE(probe.is_point_spatially_related_to_time_blobs(0, tbm, true, 0, 4));
    }
}

TEST_CASE("pr29 D1 and D12 compose: slack applies in the adjacent slice too")
{
    // The prototype applies its wire slack in all three of t, t+1, t-1.  A fix
    // that threaded wire_tol only into the same-slice call would pass every
    // case above and fail here.
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();

    Cluster& ref = g.make_child();
    add_blob(ref, 10, 100, 100, REF_U_LO, REF_U_HI, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 4, 8);
    const auto& tbm = ref.time_blob_map();

    // One wire past the reference band AND one slice away: needs both.
    Cluster& probe = g.make_child();
    add_blob(probe, 12, 100, 100, 12, 13, WIDE_LO, WIDE_HI, WIDE_LO, WIDE_HI, 0, 4);

    CHECK_FALSE(probe.is_point_spatially_related_to_time_blobs(0, tbm, true, 0, 4));  // wire misses
    CHECK_FALSE(probe.is_point_spatially_related_to_time_blobs(0, tbm, true, 1, 1));  // slice misses
    CHECK(probe.is_point_spatially_related_to_time_blobs(0, tbm, true, 1, 4));        // both
}
