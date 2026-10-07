// sbnd_xin/docs/128-129: the segment -> T_kine row map behind T_segment's
// kine_index / kine_n_rows (owner ruling 3, doc 128 sec 7).
//
// Usage:
//   wcdoctest-clus -tc="kine segment links*"

#include "WireCellClus/KineSegmentLinks.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell::Clus::PR;

TEST_CASE("kine segment links: track rows map one segment each")
{
    std::vector<KineRowMembers> rows{{false, {7}}, {false, {3}}};
    auto links = build_kine_segment_links(rows);
    REQUIRE(links.size() == 2);
    // sorted by graph index, not by row
    CHECK(links[0].graph_index == 3);
    CHECK(links[0].kine_index == 1);
    CHECK(links[0].n_rows == 1);
    CHECK(links[1].graph_index == 7);
    CHECK(links[1].kine_index == 0);
    CHECK(links[1].n_rows == 1);
}

TEST_CASE("kine segment links: a shower row links every member")
{
    std::vector<KineRowMembers> rows{{false, {1}}, {true, {4, 5, 6}}};
    auto links = build_kine_segment_links(rows);
    REQUIRE(links.size() == 4);
    for (std::size_t i = 1; i < links.size(); ++i) {
        CHECK(links[i].kine_index == 1);
        CHECK(links[i].n_rows == 1);
    }
}

TEST_CASE("kine segment links: the shower row wins over a track row (guard-off double count)")
{
    // Row 0: a shower member at the main vertex pushed as its own track row
    // (kine_mainvtx_used_guard off); row 1: its shower.
    std::vector<KineRowMembers> rows{{false, {5}}, {true, {5, 8}}};
    auto links = build_kine_segment_links(rows);
    REQUIRE(links.size() == 2);
    CHECK(links[0].graph_index == 5);
    CHECK(links[0].kine_index == 1);
    CHECK(links[0].n_rows == 2);
    CHECK(links[1].graph_index == 8);
    CHECK(links[1].kine_index == 1);
    CHECK(links[1].n_rows == 1);
}

TEST_CASE("kine segment links: two showers sharing members take the lowest shower row")
{
    std::vector<KineRowMembers> rows{{false, {2}}, {true, {9, 10}}, {true, {9, 11}}};
    auto links = build_kine_segment_links(rows);
    REQUIRE(links.size() == 4);
    CHECK(links[1].graph_index == 9);
    CHECK(links[1].kine_index == 1);
    CHECK(links[1].n_rows == 2);
    CHECK(links[3].graph_index == 11);
    CHECK(links[3].kine_index == 2);
}

TEST_CASE("kine segment links: a repeated member within one row counts once; empty input")
{
    std::vector<KineRowMembers> rows{{true, {4, 4}}};
    auto links = build_kine_segment_links(rows);
    REQUIRE(links.size() == 1);
    CHECK(links[0].n_rows == 1);
    CHECK(build_kine_segment_links({}).empty());
}

TEST_CASE("kine overlap extra copies: kine_charge_dedup counts charge-priced rows once")
{
    // Two charge-priced shower rows share a segment.
    CHECK(kine_overlap_extra_copies(2, 0, true) == 0);
    CHECK(kine_overlap_extra_copies(2, 0, false) == 1);
    CHECK(kine_overlap_extra_copies(3, 0, false) == 2);
    // A charge row and a range-priced row: the range row is not deduplicated.
    CHECK(kine_overlap_extra_copies(1, 1, true) == 1);
    CHECK(kine_overlap_extra_copies(1, 1, false) == 1);
    CHECK(kine_overlap_extra_copies(2, 1, true) == 1);
    CHECK(kine_overlap_extra_copies(2, 1, false) == 2);
    // Two non-charge rows count it twice either way.
    CHECK(kine_overlap_extra_copies(0, 2, true) == 1);
    // One row, or none: nothing extra.
    CHECK(kine_overlap_extra_copies(1, 0, true) == 0);
    CHECK(kine_overlap_extra_copies(0, 1, false) == 0);
    CHECK(kine_overlap_extra_copies(0, 0, true) == 0);
}
