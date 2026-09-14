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
