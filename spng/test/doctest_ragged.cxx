#include "WireCellSpng/Testing.h"
#include "WireCellSpng/Ragged.h"

using namespace WireCell::SPNG;

TEST_CASE("Ragged::range_lengths") {
    // Example 1: Standard ranges
    auto ranges1 = torch::tensor({{1, 3}, {10, 13}, {5, 5}}, torch::kInt64);
    auto expected1 = torch::tensor({2, 3, 0}, torch::kInt64);
    auto result1 = Ragged::range_lengths(ranges1);
    CHECK(torch::equal(result1, expected1));

    // Example 2: Empty input
    auto empty_ranges = torch::empty({0, 2}, torch::kInt64);
    auto empty_expected = torch::empty({0}, torch::kInt64);
    auto empty_result = Ragged::range_lengths(empty_ranges);
    CHECK(torch::equal(empty_result, empty_expected));

    // Example 3: Negative indices
    auto ranges3 = torch::tensor({{-5, -2}, {0, 1}}, torch::kInt64);
    auto expected3 = torch::tensor({3, 1}, torch::kInt64);
    auto result3 = Ragged::range_lengths(ranges3);
    CHECK(torch::equal(result3, expected3));
}

TEST_CASE("Ragged::range_index_expansion") {
    // Example 1: Standard ranges
    auto ranges1 = torch::tensor({{1, 3}, {10, 13}, {5, 5}}, torch::kInt64);
    // Lengths: 2, 3, 0. Total length 5. Indices: 0, 0, 1, 1, 1
    auto expected1 = torch::tensor({0, 0, 1, 1, 1}, torch::kInt64);
    auto result1 = Ragged::range_index_expansion(ranges1);
    CHECK(torch::equal(result1, expected1));

    // Example 2: Empty input
    auto empty_ranges = torch::empty({0, 2}, torch::kInt64);
    auto empty_expected = torch::empty({0}, torch::kInt64);
    auto empty_result = Ragged::range_index_expansion(empty_ranges);
    CHECK(torch::equal(empty_result, empty_expected));

    // Example 3: Single range
    auto ranges3 = torch::tensor({{100, 105}}, torch::kInt64);
    auto expected3 = torch::tensor({0, 0, 0, 0, 0}, torch::kInt64);
    auto result3 = Ragged::range_index_expansion(ranges3);
    CHECK(torch::equal(result3, expected3));
}

TEST_CASE("Ragged::range_value_expansion") {
    // Example 1: Standard ranges
    auto ranges1 = torch::tensor({{1, 3}, {10, 13}, {5, 5}}, torch::kInt64);
    // Values: [1, 2], [10, 11, 12], []
    auto expected1 = torch::tensor({1, 2, 10, 11, 12}, torch::kInt64);
    auto result1 = Ragged::range_value_expansion(ranges1);
    CHECK(torch::equal(result1, expected1));

    // Example 2: Empty input
    auto empty_ranges = torch::empty({0, 2}, torch::kInt64);
    auto empty_expected = torch::empty({0}, torch::kInt64);
    auto empty_result = Ragged::range_value_expansion(empty_ranges);
    CHECK(torch::equal(empty_result, empty_expected));

    // Example 3: Single range
    auto ranges3 = torch::tensor({{5, 8}}, torch::kInt64);
    auto expected3 = torch::tensor({5, 6, 7}, torch::kInt64);
    auto result3 = Ragged::range_value_expansion(ranges3);
    CHECK(torch::equal(result3, expected3));

    // Example 4: Ranges with negative indices
    auto ranges4 = torch::tensor({{-5, -2}, {-1, 1}}, torch::kInt64);
    // Range 1: [-5, -4, -3]
    // Range 2: [-1, 0]
    auto expected4 = torch::tensor({-5, -4, -3, -1, 0}, torch::kInt64);
    auto result4 = Ragged::range_value_expansion(ranges4);
    CHECK(torch::equal(result4, expected4));
}
