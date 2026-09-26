// sbnd_xin/docs/109: the flash grouping behind NuBundleCensus::flash_group.
//
// Two flashes on DIFFERENT TPCs closer in time than dt are one physical flash
// seen by both TPCs' light detectors.  The group id is the smallest gid in
// the group, so it does not depend on the order the flashes are listed in.

#include "WireCellUtil/doctest.h"
#include "WireCellClus/NuBundleCensus.h"

#include <cmath>

using WireCell::Clus::PR::group_flashes;

TEST_CASE("nu census: a different-TPC pair inside dt is one group")
{
    auto g = group_flashes({1000003, 5}, {1, 0}, {1.230, 1.260}, 0.05);
    REQUIRE(g.size() == 2);
    CHECK(g[0] == 5);
    CHECK(g[1] == 5);
}

TEST_CASE("nu census: same TPC, outside dt, or dt <= 0 never groups")
{
    auto same = group_flashes({1, 2}, {0, 0}, {1.0, 1.0}, 0.05);
    CHECK(same[0] == 1);
    CHECK(same[1] == 2);
    auto far = group_flashes({1, 1000002}, {0, 1}, {1.0, 1.05}, 0.05);   // |dt| == dt is outside
    CHECK(far[0] == 1);
    CHECK(far[1] == 1000002);
    auto off = group_flashes({1, 1000002}, {0, 1}, {1.0, 1.0}, 0.0);
    CHECK(off[0] == 1);
    CHECK(off[1] == 1000002);
}

TEST_CASE("nu census: an unresolved (NaN) flash time never groups")
{
    auto g = group_flashes({1, 1000002}, {0, 1}, {1.0, std::nan("")}, 0.05);
    CHECK(g[0] == 1);
    CHECK(g[1] == 1000002);
}

TEST_CASE("nu census: grouping is transitive and order independent")
{
    // 7 (tpc 0) - 1000001 (tpc 1) - 3 (tpc 0): a chain of different-TPC pairs.
    auto a = group_flashes({7, 1000001, 3}, {0, 1, 0}, {1.00, 1.03, 1.06}, 0.05);
    auto b = group_flashes({3, 7, 1000001}, {0, 0, 1}, {1.06, 1.00, 1.03}, 0.05);
    CHECK(a[0] == 3);
    CHECK(a[1] == 3);
    CHECK(a[2] == 3);
    CHECK(b[0] == 3);
    CHECK(b[1] == 3);
    CHECK(b[2] == 3);
}

// ---------------------------------------------------------------------- //
// sbnd_xin/docs/109 rev 3: dedup_flash_groups -- the rule nu_dedup_flash_group
// applies on top of the grouping above.  The candidate list arrives in the
// selection's own order (longest selected activity first, ties by gid), so
// "keep the first of each group" is "keep the longest".

#include <map>

using WireCell::Clus::PR::dedup_flash_groups;
using WireCell::Clus::PR::cathode_pair_candidates;
using WireCell::Clus::PR::bundle_contact;
using WireCell::Clus::PR::merge_bundles;

TEST_CASE("nu dedup: two candidates from one physical flash collapse to the first")
{
    // gids 5 (tpc 0) and 1000006 (tpc 1) 6 ns apart -> one group, id 5.
    // This is sbnd_xin r472 s36 e40 (the colleague's case 5): the W-side row
    // came first (longer), the E-side row is the true neutrino's.  The rule
    // keeps the LONGER one -- which candidate that is, is the selection's
    // ranking, not this function's business.
    const auto g = group_flashes({1000006, 5}, {1, 0}, {1.583, 1.577}, 0.05);
    std::map<int, int> gid_group{{1000006, g[0]}, {5, g[1]}};
    REQUIRE(gid_group[1000006] == 5);
    REQUIRE(gid_group[5] == 5);

    const auto keep = dedup_flash_groups({1000006, 5}, gid_group);
    REQUIRE(keep.size() == 1);
    CHECK(keep[0] == 0);          // the first (longest) candidate survives
}

TEST_CASE("nu dedup: candidates in different groups are all kept")
{
    std::map<int, int> gid_group{{2, 2}, {9, 9}, {1000004, 1000004}};
    const auto keep = dedup_flash_groups({2, 9, 1000004}, gid_group);
    REQUIRE(keep.size() == 3);
    CHECK(keep[0] == 0);
    CHECK(keep[1] == 1);
    CHECK(keep[2] == 2);
}

TEST_CASE("nu dedup: a gid with no known group is never a duplicate")
{
    // An empty map: every candidate falls back to its own gid as its group.
    const auto keep = dedup_flash_groups({4, 1000001, 12}, {});
    CHECK(keep.size() == 3);
    // ... and a partial map only collapses the pair it knows about.
    std::map<int, int> partial{{4, 4}, {1000001, 4}};
    const auto keep2 = dedup_flash_groups({4, 1000001, 12}, partial);
    REQUIRE(keep2.size() == 2);
    CHECK(keep2[0] == 0);
    CHECK(keep2[1] == 2);
}

TEST_CASE("nu dedup: a three-way group keeps exactly one, and one candidate is a no-op")
{
    std::map<int, int> gid_group{{3, 3}, {1000001, 3}, {7, 3}};
    const auto keep = dedup_flash_groups({3, 1000001, 7}, gid_group);
    REQUIRE(keep.size() == 1);
    CHECK(keep[0] == 0);

    const auto single = dedup_flash_groups({5}, {{5, 5}});
    REQUIRE(single.size() == 1);
    CHECK(single[0] == 0);
    CHECK(dedup_flash_groups({}, {}).empty());
}

// sbnd_xin/docs/109 rev 4 -- nu_bundle_flash_group.  One physical flash seen by
// both drift volumes can be ONE interaction split across the two (merge) or
// TWO interactions, one per volume (keep both).  Whether the two bundles'
// charge touches is what tells them apart.

TEST_CASE("nu bundle: only same-group different-TPC pairs are tested")
{
    // r472 s36 e40: gids 5 (tpc 0) and 1000006 (tpc 1), 6 ns apart -> one group.
    std::map<int, int> group{{5, 5}, {1000006, 5}};
    std::map<int, int> tpc{{5, 0}, {1000006, 1}};
    const auto p = cathode_pair_candidates({1000006, 5}, group, tpc);
    REQUIRE(p.size() == 1);
    CHECK(p[0].first == 5);            // ascending (min, max), whatever the input order
    CHECK(p[0].second == 1000006);

    // r474 s72 e31: gids 2 (tpc 0) and 1000011 (tpc 1) are 709 ns apart -> two
    // groups -> never tested, so two objects can never be merged by accident.
    std::map<int, int> group2{{2, 2}, {1000011, 1000011}};
    std::map<int, int> tpc2{{2, 0}, {1000011, 1}};
    CHECK(cathode_pair_candidates({2, 1000011}, group2, tpc2).empty());

    // Same group but the same TPC (a 3-way group's same-side members): not a
    // cathode pair.  An unknown TPC is never paired either.
    std::map<int, int> group3{{3, 3}, {4, 3}, {1000001, 3}};
    std::map<int, int> tpc3{{3, 0}, {4, 0}, {1000001, 1}};
    const auto p3 = cathode_pair_candidates({3, 4, 1000001}, group3, tpc3);
    REQUIRE(p3.size() == 2);
    CHECK(p3[0] == std::make_pair(3, 1000001));
    CHECK(p3[1] == std::make_pair(4, 1000001));
    std::map<int, int> tpc4{{3, 0}};
    CHECK(cathode_pair_candidates({3, 1000001}, group3, tpc4).empty());
    // An unmapped gid is its own group.
    CHECK(cathode_pair_candidates({3, 1000001}, {}, tpc3).empty());
}

TEST_CASE("nu bundle: contact is a distance; the cathode window only when asked")
{
    const double cx = 0.0, gap = 20.0;
    // Default: xcut 0 = distance alone.  r472 s36 e40's two bundles touch at
    // 0.4 cm at x = -10.8 / -11.1 (the VERTEX, not the seam): merged.
    CHECK(bundle_contact(0.4, -10.8, -11.1, cx, 0.0, gap));
    CHECK_FALSE(bundle_contact(30.0, -2.1, 3.0, cx, 0.0, gap));          // too far apart: two objects
    CHECK_FALSE(bundle_contact(20.0, 0.0, 0.0, cx, 0.0, gap));           // boundary is exclusive
    CHECK_FALSE(bundle_contact(std::nan(""), 0.0, 0.0, cx, 0.0, gap));   // NaN never contacts
    // xcut > 0: the optional cathode window, both points within it.
    const double xcut = 6.0;
    CHECK(bundle_contact(7.8, -2.1, 3.0, cx, xcut, gap));                // halves meeting at the seam
    CHECK_FALSE(bundle_contact(0.4, -10.8, -11.1, cx, xcut, gap));       // r472 would be REFUSED by the window
    CHECK_FALSE(bundle_contact(7.8, -2.1, 9.0, cx, xcut, gap));          // one side off the seam
    CHECK(bundle_contact(7.8, 100.0 - 2.1, 103.0, 100.0, xcut, gap));    // the plane is a parameter
    CHECK_FALSE(bundle_contact(1.0, std::nan(""), 0.0, cx, xcut, gap));
    CHECK_FALSE(bundle_contact(1.0, 6.0, 0.0, cx, xcut, gap));           // window boundary is exclusive
}

TEST_CASE("nu bundle: merged bundles take the smallest gid as root, the rest stay themselves")
{
    // The r472 pair in contact: one bundle, root 5.
    auto r = merge_bundles({5, 1000006}, {{5, 1000006}});
    CHECK(r[5] == 5);
    CHECK(r[1000006] == 5);
    // The same pair NOT in contact (the two-neutrino reading): two bundles.
    auto k = merge_bundles({5, 1000006}, {});
    CHECK(k[5] == 5);
    CHECK(k[1000006] == 1000006);
    // A three-way group where only one pair touches: the third keeps its row.
    auto t = merge_bundles({3, 4, 1000001}, {{4, 1000001}});
    CHECK(t[3] == 3);
    CHECK(t[4] == 4);
    CHECK(t[1000001] == 4);
    // Transitive contact, and a contact naming an unknown gid is ignored.
    auto u = merge_bundles({3, 4, 1000001}, {{3, 1000001}, {4, 1000001}, {7, 3}});
    CHECK(u[3] == 3);
    CHECK(u[4] == 3);
    CHECK(u[1000001] == 3);
    CHECK(u.count(7) == 0);
    CHECK(merge_bundles({}, {}).empty());
}
