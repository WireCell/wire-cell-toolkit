#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Math.h"

using namespace WireCell;

TEST_SUITE("coprime") {

TEST_CASE("coprime") {
    const double fraction = 0.02;
    for (size_t capacity = 3; capacity < 10000; ++capacity) {
        const size_t target = (1-fraction)*capacity;
        auto got = nearest_coprime(capacity, target);
        REQUIRE(got != 0);
        int err = got - target;
        spdlog::debug("capacity={} target={} got={} error={}", capacity, target, got, err);
    }
}

}  // TEST_SUITE("coprime")
