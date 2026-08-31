#include "WireCellUtil/RayTools.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

using namespace WireCell::RayGrid;

TEST_SUITE("raytools") {

TEST_CASE("relative distance") {
    SUBCASE("full pitch") {
        double d = relative_distance({5,60}, {60,90});
        spdlog::debug("relative_distance({5,60}, {60,90}) = {}", d);
        CHECK(std::abs(d - 1.0) < 1e-8);
    }
    SUBCASE("zero distance") {
        double d = relative_distance({5,15}, {9,11});
        spdlog::debug("relative_distance({5,15}, {9,11}) = {}", d);
        CHECK(std::abs(d) < 1e-8);
    }
    SUBCASE("three pitches") {
        double d = relative_distance({1,2}, {4,5});
        spdlog::debug("relative_distance({1,2}, {4,5}) = {}", d);
        CHECK(std::abs(d - 3.0) < 1e-8);
    }
}

}  // TEST_SUITE("raytools")
