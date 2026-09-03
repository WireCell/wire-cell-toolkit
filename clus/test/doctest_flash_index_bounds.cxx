// Bounds on the matched-flash index (sbnd_xin/docs/99).
//
// Grouping::flash_at() resolves an integer as a ROW INDEX into the grouping's
// canonical "flash" point cloud.  Its declared contract (Facade_Grouping.h) is
// that an index which is negative or out of range yields an INVALID Flash, but
// the implementation used to skip the range half: get_element() forwards to
// Array::element(), an unchecked pointer offset, so a stale index returned raw
// memory while operator bool() still reported "true".  A cluster's "flash"
// scalar is exactly such a stored index, and in SBND production 758 of 50699
// matched clusters addressed a row that does not exist.
//
// The bound is deliberately taken on the "time" array -- the same array
// flashes() sizes its enumeration loop with -- so that flashes(), and with it
// the QLMatching read path, provably cannot move.  The last two cases pin that
// invariant: a MISSING companion array must default its own field and must not
// invalidate the flash or shorten the enumeration.  (A merely *short* companion
// array is not reachable: Dataset::add() refuses an array whose major axis
// disagrees with the ones already in the Dataset, so all four arrays of a
// well-formed "flash" PC have the same row count.)

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Flash.h"

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;

// Install a canonical "flash" PC of n rows.  Row i: time=100+i, value=1000*(i+1),
// ident=10+i, type=0.  with_ident=false omits the "ident" array entirely, the
// way a legacy archive can.
static void put_flash_pc(Grouping& g, size_t n, bool with_ident = true)
{
    std::vector<double> time(n), value(n);
    std::vector<int> type(n, 0), ident(n);
    for (size_t i = 0; i < n; ++i) {
        time[i] = 100.0 + (double) i;
        value[i] = 1000.0 * (double) (i + 1);
        ident[i] = 10 + (int) i;
    }

    Dataset ds;
    ds.add("time", Array(time));
    ds.add("value", Array(value));
    if (with_ident) ds.add("ident", Array(ident));
    ds.add("type", Array(type));
    g.local_pcs()["flash"] = ds;
}

TEST_CASE("flash_at yields an invalid Flash when nothing addresses a row")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();

    // No "flash" PC at all.
    CHECK(!g->flash_at(0));

    put_flash_pc(*g, 4);

    CHECK(!g->flash_at(-1));    // negative
    CHECK(!g->flash_at(4));     // one past the end -- the defect
    CHECK(!g->flash_at(13));    // far past the end (the shape seen in production)
}

TEST_CASE("flash_at returns the addressed row when the index is in range")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    put_flash_pc(*g, 4);

    for (int i = 0; i < 4; ++i) {
        auto f = g->flash_at(i);
        REQUIRE(bool(f));
        CHECK(f.time() == doctest::Approx(100.0 + i));
        CHECK(f.value() == doctest::Approx(1000.0 * (i + 1)));
        CHECK(f.ident() == 10 + i);
        CHECK(f.type() == 0);
    }
}

TEST_CASE("flash_at bound is causal: shrinking the PC invalidates the same index")
{
    // Negative control for the guard -- corrupt exactly the property it keys on
    // (the row count of the "flash" PC) and watch a previously-good index stop
    // resolving, with no other input changed.
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();

    put_flash_pc(*g, 4);
    auto before = g->flash_at(3);
    REQUIRE(bool(before));
    CHECK(before.ident() == 13);

    put_flash_pc(*g, 3);        // same index, one fewer row
    auto after = g->flash_at(3);
    CHECK(!after);
    CHECK(after.ident() == -1); // the documented invalid default, not raw memory
}

TEST_CASE("flashes() enumeration is unchanged by the bound")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();

    CHECK(g->flashes().empty());          // no PC

    put_flash_pc(*g, 5);
    auto all = g->flashes();
    REQUIRE(all.size() == 5);
    for (size_t i = 0; i < all.size(); ++i) {
        CHECK(bool(all[i]));              // every enumerated flash stays valid
        CHECK(all[i].ident() == 10 + (int) i);
    }
}

TEST_CASE("a missing companion array defaults its own field only")
{
    // No "ident" array at all.  Validity is keyed on "time", so every flash must
    // still enumerate and stay valid; only ident falls back to the documented
    // default.  This is what the old get_element(..., def) call gave, and it is
    // the reason the range guard is NOT taken across all four arrays.
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    put_flash_pc(*g, 4, /*with_ident=*/false);

    auto all = g->flashes();
    REQUIRE(all.size() == 4);
    for (size_t i = 0; i < all.size(); ++i) {
        CHECK(bool(all[i]));
        CHECK(all[i].ident() == -1);
    }
    CHECK(all[3].time() == doctest::Approx(103.0));
    CHECK(all[3].value() == doctest::Approx(4000.0));
}

TEST_CASE("Cluster::get_flash() on a stale index is invalid, not garbage")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    put_flash_pc(*g, 4);

    Cluster& good = g->make_child();
    good.set_scalar<int>("flash", 2);
    auto fg = good.get_flash();
    REQUIRE(bool(fg));
    CHECK(fg.ident() == 12);
    CHECK(fg.time() == doctest::Approx(102.0));

    // The production shape: the scalar is a row id from a flash list that is no
    // longer the one in this grouping, so it points past the end.
    Cluster& stale = g->make_child();
    stale.set_scalar<int>("flash", 13);
    auto fs = stale.get_flash();
    CHECK(!fs);
    CHECK(fs.ident() == -1);
    CHECK(fs.time() == 0.0);
    CHECK(fs.value() == 0.0);

    // A cluster that never matched keeps the -1 sentinel.
    Cluster& none = g->make_child();
    none.set_scalar<int>("flash", -1);
    CHECK(!none.get_flash());
}
