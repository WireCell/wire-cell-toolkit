// Resolving a cluster's matched flash by GID (sbnd_xin/docs/99).
//
// THE DEFECT.  QLMatching stamps two keys on every matched cluster: "flash", the
// row index of the matched flash within THAT INPUT's flash list, and
// "matched_flash_gid", a globally-unique id.  The multi-input merge concatenates
// only the names in root_pcs_to_merge ("opflash") and DROPS every other root PC
// of the non-primary inputs, so the merged grouping keeps only the primary
// input's canonical "flash" PC.  A cluster matched on any other input then has a
// row index into a list that is not in the archive: out of range it yields an
// invalid Flash (doctest_flash_index_bounds), but IN range it silently yields a
// different, real flash.  Measured on SBND production: 50.6% correct.
//
// THE FIX UNDER TEST.  Grouping::flash_by_gid() resolves against the merge-safe
// "opflash" PC (one row per (flash, channel): gid/time/ch/pe), whose gid is the
// same number the cluster carries; Cluster::get_matched_flash() delegates to it.
//
// Cases below, in order: the invalid inputs; a gid carrying the
// gid_side*1000000 offset of a non-primary input; the channel-order sum (with a
// control proving the test can tell the orders apart); a causal negative control
// that corrupts exactly the "gid" column; the gid-collision refusal; and the
// end-to-end production shape where get_flash() resolves the WRONG real flash on
// the very cluster get_matched_flash() gets right.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Flash.h"

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;

static const int kStride = 1000000;      // QLMatching's kFlashGidStride
static const int kNchan = 4;

// One flash: gid, time, and kNchan per-channel PE, emitted in `order` (a
// permutation of the channel ids) so a caller can make the row order disagree
// with the channel order the way the merge does.
struct FlashRows {
    int gid;
    double time;
    std::vector<double> pe_by_ch;        // indexed by channel
    std::vector<int> order;              // row emission order, channel ids
};

static void put_opflash_pc(Grouping& g, const std::vector<FlashRows>& flashes)
{
    std::vector<int> gid, ch;
    std::vector<double> time, pe;
    for (const auto& f : flashes) {
        for (int c : f.order) {
            gid.push_back(f.gid);
            time.push_back(f.time);
            ch.push_back(c);
            pe.push_back(f.pe_by_ch.at(c));
        }
    }
    Dataset ds;
    ds.add("gid", Array(gid));
    ds.add("time", Array(time));
    ds.add("ch", Array(ch));
    ds.add("pe", Array(pe));
    g.local_pcs()["opflash"] = ds;
}

// A flash whose channel c holds PE 10*(c+1), rows in ascending channel order.
static FlashRows simple(int gid, double time)
{
    FlashRows f{gid, time, {}, {}};
    for (int c = 0; c < kNchan; ++c) {
        f.pe_by_ch.push_back(10.0 * (c + 1));
        f.order.push_back(c);
    }
    return f;
}

TEST_CASE("flash_by_gid yields an invalid Flash when nothing resolves")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();

    CHECK(!g->flash_by_gid(0));          // no "opflash" PC at all
    CHECK(!g->flash_by_gid(-1));

    put_opflash_pc(*g, {simple(0, 100.0), simple(kStride + 2, 200.0)});

    CHECK(!g->flash_by_gid(-1));         // the unmatched sentinel
    CHECK(!g->flash_by_gid(1));          // a gid no row carries
    CHECK(!g->flash_by_gid(kStride));    // right side, wrong index
    CHECK(g->flash_by_gid(1).ident() == -1);   // and the invalid default, not garbage
}

TEST_CASE("flash_by_gid resolves both the primary and the offset input")
{
    // gid = gid_side*1000000 + index-in-that-input's-flash-list.  The whole
    // point is that the SECOND input's flashes resolve too -- those are the ones
    // the merge drops the flash PC for.
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    put_opflash_pc(*g, {simple(0, 100.0), simple(3, 130.0),
                        simple(kStride + 0, 200.0), simple(kStride + 5, 250.0)});

    const double want_pe = 10.0 + 20.0 + 30.0 + 40.0;

    auto a = g->flash_by_gid(3);
    REQUIRE(bool(a));
    CHECK(a.time() == doctest::Approx(130.0));
    CHECK(a.value() == doctest::Approx(want_pe));
    CHECK(a.ident() == 3);               // the ident IS the gid on this path
    CHECK(a.idents().size() == (size_t) kNchan);

    auto b = g->flash_by_gid(kStride + 5);
    REQUIRE(bool(b));
    CHECK(b.time() == doctest::Approx(250.0));
    CHECK(b.value() == doctest::Approx(want_pe));
    CHECK(b.ident() == kStride + 5);

    // pes() is the dense per-channel view consumers use.
    auto dense = b.pes(kNchan);
    REQUIRE(dense.size() == (size_t) kNchan);
    for (int c = 0; c < kNchan; ++c) CHECK(dense[c] == doctest::Approx(10.0 * (c + 1)));
}

TEST_CASE("flash_by_gid sums PE in channel order, not row order")
{
    // The canonical flash PC's "value" is FlashTensorToOpticalPCs' `sum += pe`
    // over channels 0..nchan-1.  Reproducing it BIT-exactly therefore requires
    // summing in channel order; the merge does not preserve that in row order.
    // Values chosen so the two orders give genuinely different doubles -- the
    // REQUIRE below is the control that makes this case able to fail.
    // One bright channel and three faint ones: added last, the three faint ones
    // vanish into the bright one's ulp (sum 1e16); added first they accumulate to
    // 3 and survive the rounding (sum 1e16+4).
    const std::vector<double> pe_by_ch{1e16, 1.0, 1.0, 1.0};
    const std::vector<int> scrambled{1, 2, 3, 0};

    double asc = 0.0;
    for (int c = 0; c < kNchan; ++c) asc += pe_by_ch[c];
    double as_emitted = 0.0;
    for (int c : scrambled) as_emitted += pe_by_ch[c];
    REQUIRE(asc != as_emitted);          // otherwise this test proves nothing

    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    put_opflash_pc(*g, {FlashRows{7, 300.0, pe_by_ch, scrambled}});

    auto f = g->flash_by_gid(7);
    REQUIRE(bool(f));
    CHECK(f.value() == asc);             // exact, not Approx
    CHECK(f.value() != as_emitted);
    // The per-channel vectors come back in channel order too.
    auto ids = f.idents();
    REQUIRE(ids.size() == (size_t) kNchan);
    for (int c = 0; c < kNchan; ++c) CHECK(ids[c] == c);
}

TEST_CASE("flash_by_gid is causal: corrupting only the gid column stops resolution")
{
    // Negative control.  Everything else about the PC is untouched -- same rows,
    // same times, same PE -- only the key the lookup is built on changes.
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();

    put_opflash_pc(*g, {simple(0, 100.0), simple(kStride + 5, 250.0)});
    auto before = g->flash_by_gid(kStride + 5);
    REQUIRE(bool(before));
    CHECK(before.time() == doctest::Approx(250.0));

    // Re-emit with that flash's gid changed and nothing else.
    put_opflash_pc(*g, {simple(0, 100.0), simple(kStride + 6, 250.0)});
    auto after = g->flash_by_gid(kStride + 5);
    CHECK(!after);
    CHECK(after.time() == 0.0);
    CHECK(after.value() == 0.0);
    // ... and the flash is still there under its new key, so nothing else broke.
    CHECK(bool(g->flash_by_gid(kStride + 6)));
}

TEST_CASE("flash_by_gid refuses a gid that names more than one flash")
{
    // The gid-uniqueness precondition (Facade_Grouping.h): a config using
    // opflash_phys_gid can have two inputs emit the same gid for different
    // flashes.  A channel appearing twice under one gid is that condition, and
    // must yield an invalid Flash rather than a doubled PE sum.
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    put_opflash_pc(*g, {simple(4, 100.0), simple(4, 900.0)});

    auto f = g->flash_by_gid(4);
    CHECK(!f);
    CHECK(f.value() == 0.0);             // NOT 2 x the single-flash total
}

TEST_CASE("get_matched_flash() is right where get_flash() is wrong")
{
    // The production shape.  The merged grouping holds the PRIMARY input's flash
    // PC (2 rows) and the merge-safe opflash PC (all inputs).  A cluster matched
    // on the second input carries flash=1 -- an in-range row of the WRONG list --
    // and matched_flash_gid=kStride+1.  cluster_t0 is the bit-exact truth: it was
    // set from the flash the cluster actually matched.
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();

    // Primary input's canonical flash PC: two flashes at t=100, 110.
    {
        std::vector<double> time{100.0, 110.0}, value{40.0, 40.0};
        std::vector<int> ident{0, 1}, type{0, 0};
        Dataset ds;
        ds.add("time", Array(time));
        ds.add("value", Array(value));
        ds.add("ident", Array(ident));
        ds.add("type", Array(type));
        g->local_pcs()["flash"] = ds;
    }
    // Merge-safe optical PC: both inputs.
    put_opflash_pc(*g, {simple(0, 100.0), simple(1, 110.0),
                        simple(kStride + 0, 700.0), simple(kStride + 1, 750.0)});

    Cluster& c = g->make_child();
    c.set_scalar<int>("flash", 1);                    // row 1 of the WRONG list
    c.set_scalar<int>("matched_flash_gid", kStride + 1);
    c.set_cluster_t0(750.0);

    auto wrong = c.get_flash();
    REQUIRE(bool(wrong));                             // a real flash ...
    CHECK(wrong.time() == doctest::Approx(110.0));    // ... and not this cluster's
    CHECK(wrong.time() != c.get_cluster_t0());

    auto right = c.get_matched_flash();
    REQUIRE(bool(right));
    CHECK(right.time() == c.get_cluster_t0());        // the archive-local truth test
    CHECK(right.ident() == kStride + 1);
    CHECK(right.value() == doctest::Approx(100.0));

    // A cluster with no matched flash resolves to nothing on either path.
    Cluster& none = g->make_child();
    none.set_scalar<int>("flash", -1);
    none.set_scalar<int>("matched_flash_gid", -1);
    CHECK(!none.get_flash());
    CHECK(!none.get_matched_flash());

    // A cluster whose gid is stale (its flash is not in this grouping) gets the
    // sentinel, not a neighbouring flash.
    Cluster& stale = g->make_child();
    stale.set_scalar<int>("flash", 0);
    stale.set_scalar<int>("matched_flash_gid", kStride + 9);
    CHECK(!stale.get_matched_flash());
    CHECK(stale.get_matched_flash().ident() == -1);
}
