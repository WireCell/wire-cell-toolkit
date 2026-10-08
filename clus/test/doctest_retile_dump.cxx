// doc pdvd/129 phase 1: the cell provenance of the retiler's study dump.
#include "WireCellUtil/doctest.h"
#include "../src/RetileDump.h"

using namespace WireCell::Clus;

static RetileDump::activity_map_t make_map()
{
    RetileDump::activity_map_t m;
    // 5 layers: 2 dummy, then U, V, W with 4 wires each
    std::vector<WireCell::RayGrid::measure_t> layers = {{1}, {1}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    m[{100, 104}] = layers;
    m[{104, 108}] = layers;
    return m;
}

TEST_CASE("retile dump: sentinel cells are found by plane, tick and wire")
{
    auto m = make_map();
    m[{100, 104}][2][1] = 1e-3;     // U wire 1
    m[{100, 104}][3][2] = 5000.0;   // measured: not a sentinel
    m[{104, 108}][4][3] = 1e-3;     // W wire 3
    auto s = RetileDump::sentinel_cells(m);
    CHECK(s.size() == 2);
    CHECK(s.count({0, 100, 1}) == 1);
    CHECK(s.count({2, 104, 3}) == 1);
    CHECK(s.count({1, 100, 2}) == 0);
}

TEST_CASE("retile dump: a sentinel present before the paint is class 1, a new one is class 2")
{
    auto m = make_map();
    m[{100, 104}][2][1] = 1e-3;
    const auto before = RetileDump::sentinel_cells(m);
    m[{104, 108}][3][0] = 1e-3;     // painted afterwards
    m[{100, 104}][4][2] = 7000.0;   // measured charge is never listed
    auto rows = RetileDump::classify(before, m);
    REQUIRE(rows.size() == 2);
    CHECK(rows[0] == std::array<int, 4>{0, 100, 1, 1});
    CHECK(rows[1] == std::array<int, 4>{1, 104, 0, 2});
}

TEST_CASE("retile dump: zero and measured values are not sentinels")
{
    CHECK(RetileDump::is_sentinel(1e-3));
    CHECK_FALSE(RetileDump::is_sentinel(0.0));
    CHECK_FALSE(RetileDump::is_sentinel(1.0e-3 + 1e-9));
}
