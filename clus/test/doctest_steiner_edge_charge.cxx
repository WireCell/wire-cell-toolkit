// Regression pins for the Steiner edge-weight charge fidelity fix, doc pr/29 D2.
//
// create_steiner_tree is called with disable_dead_mix_cell = false and forwards
// it to find_steiner_terminals, but historically DROPPED it at the
// create_enhanced_steiner_graph call, so the Qs/Qt entering the edge weight
// were computed with that function's `= true` default instead.  The fix threads
// the caller's value through behind Config::edge_charge_forward_dead_mix.
//
// The knob itself is a one-line ternary; what needs pinning is the claim that
// justifies it -- that the two branches of Cluster::calc_charge_wcp are NOT
// interchangeable.  They select planes by independent predicates:
//
//   disable_dead_mix_cell = true   all three planes enter the RMS, then DEAD
//                                  ones (charge_uncertainty > 1e10) are removed
//   disable_dead_mix_cell = false  only planes with a NONZERO charge value enter
//
// So they agree exactly when "dead" and "zero-valued" pick out the same planes,
// and the two cases below are the two ways that fails.  If someone ever
// "simplifies" one branch into the other, these fail.
//
// Deliberately a separate file with its own fixture rather than an addition to
// doctest_steiner_terminal_filter.cxx: that file's blobs carry no charge arrays
// and D2 is a different function on a different axis.

#include "WireCellUtil/doctest.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"

#include <cmath>
#include <vector>

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;

namespace {

    using fa_float_t = double;
    using fa_int_t = int;

    int test_wpid() { return WirePlaneId(kUlayer, 0, 0).ident(); }

    // A cluster holding exactly one point, carrying the per-plane charge value
    // and uncertainty arrays calc_charge_wcp reads.  An uncertainty above the
    // hard-coded 1e10 dead_threshold makes that plane dead.
    // Duplicated rather than shared with the other clus doctests (CLAUDE.md M10).
    void add_point(Cluster& cl,
                   double qu, double qv, double qw,
                   double uu, double uv, double uw)
    {
        Blob& bl = cl.make_child();

        Dataset ds;
        ds.add("x", Array(std::vector<double>{0.0}));
        ds.add("y", Array(std::vector<double>{0.0}));
        ds.add("z", Array(std::vector<double>{0.0}));
        ds.add("uwire_index", Array(std::vector<int>{10}));
        ds.add("vwire_index", Array(std::vector<int>{10}));
        ds.add("wwire_index", Array(std::vector<int>{10}));
        ds.add("ucharge_val", Array(std::vector<double>{qu}));
        ds.add("vcharge_val", Array(std::vector<double>{qv}));
        ds.add("wcharge_val", Array(std::vector<double>{qw}));
        ds.add("ucharge_unc", Array(std::vector<double>{uu}));
        ds.add("vcharge_unc", Array(std::vector<double>{uv}));
        ds.add("wcharge_unc", Array(std::vector<double>{uw}));
        bl.value().local_pcs()["3d"] = ds;

        bl.value().local_pcs()["scalar"] = Dataset({
            {"charge", Array({(fa_float_t) 1.0})},
            {"center_x", Array({(fa_float_t) 0.0})},
            {"center_y", Array({(fa_float_t) 0.0})},
            {"center_z", Array({(fa_float_t) 0.0})},
            {"wpid", Array({(fa_int_t) test_wpid()})},
            {"npoints", Array({(fa_int_t) 1})},
            {"slice_index_min", Array({(fa_int_t) 0})},
            {"slice_index_max", Array({(fa_int_t) 4})},
            {"u_wire_index_min", Array({(fa_int_t) 10})},
            {"u_wire_index_max", Array({(fa_int_t) 11})},
            {"v_wire_index_min", Array({(fa_int_t) 10})},
            {"v_wire_index_max", Array({(fa_int_t) 11})},
            {"w_wire_index_min", Array({(fa_int_t) 10})},
            {"w_wire_index_max", Array({(fa_int_t) 11})},
            {"max_wire_interval", Array({(fa_int_t) 1})},
            {"min_wire_interval", Array({(fa_int_t) 1})},
            {"max_wire_type", Array({(fa_int_t) 0})},
            {"min_wire_type", Array({(fa_int_t) 0})},
        });
    }

    const double LIVE = 100.0;      // ordinary uncertainty
    const double DEAD = 1e11;       // above the 1e10 dead_threshold
    const double CUT  = 4000.0;     // the charge_cut every caller in this chain uses

}  // namespace

// ---------------------------------------------------------------------------
// The two ways the predicates disagree
// ---------------------------------------------------------------------------

TEST_CASE("pr29 D2 branches differ on a LIVE plane reading zero charge")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();

    Cluster& cl = g.make_child();
    // U is perfectly live but measured no charge; V and W are ordinary.
    add_point(cl, 0.0, 5000.0, 5000.0, LIVE, LIVE, LIVE);

    auto with_true  = cl.calc_charge_wcp(0, CUT, true);
    auto with_false = cl.calc_charge_wcp(0, CUT, false);

    // true: U is not dead, so its zero stays in the sum over all three planes.
    CHECK(with_true.second == doctest::Approx(std::sqrt(50e6 / 3.0)));
    // false: the zero plane is simply not summed.
    CHECK(with_false.second == doctest::Approx(5000.0));
    CHECK(with_true.second != doctest::Approx(with_false.second));

    // The boolean half of the return disagrees too, and it is the one the
    // terminal selection reads: `false` treats a zero-charge plane as
    // satisfied, `true` leaves it failing the charge_cut.
    CHECK_FALSE(with_true.first);
    CHECK(with_false.first);
}

TEST_CASE("pr29 D2 branches differ on a DEAD plane reading nonzero charge")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();

    Cluster& cl = g.make_child();
    // U is dead but still carries a value; the three values are distinct so
    // dropping one actually moves the RMS.
    add_point(cl, 3000.0, 5000.0, 7000.0, DEAD, LIVE, LIVE);

    auto with_true  = cl.calc_charge_wcp(0, CUT, true);
    auto with_false = cl.calc_charge_wcp(0, CUT, false);

    // true: all three summed (83e6), then the dead U subtracted -> 74e6 over 2.
    CHECK(with_true.second == doctest::Approx(std::sqrt(74e6 / 2.0)));
    // false: nothing is zero, so all three stay -> 83e6 over 3.
    CHECK(with_false.second == doctest::Approx(std::sqrt(83e6 / 3.0)));
    CHECK(with_true.second != doctest::Approx(with_false.second));
}

// ---------------------------------------------------------------------------
// ... and the case where they agree, so the two above are not read as "this
// knob changes every point".
// ---------------------------------------------------------------------------

TEST_CASE("pr29 D2 branches agree when every plane is live and nonzero")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();

    Cluster& cl = g.make_child();
    add_point(cl, 3000.0, 5000.0, 7000.0, LIVE, LIVE, LIVE);

    auto with_true  = cl.calc_charge_wcp(0, CUT, true);
    auto with_false = cl.calc_charge_wcp(0, CUT, false);

    CHECK(with_true.second == doctest::Approx(with_false.second));
    CHECK(with_true.first == with_false.first);
}

TEST_CASE("pr29 D2 branches agree when a dead plane also reads zero")
{
    Points::node_t root;
    Grouping& g = *root.value.facade<Grouping>();

    Cluster& cl = g.make_child();
    // The configuration the two predicates were presumably assumed to share.
    add_point(cl, 0.0, 5000.0, 7000.0, DEAD, LIVE, LIVE);

    auto with_true  = cl.calc_charge_wcp(0, CUT, true);
    auto with_false = cl.calc_charge_wcp(0, CUT, false);

    CHECK(with_true.second == doctest::Approx(with_false.second));
    CHECK(with_true.first == with_false.first);
}
