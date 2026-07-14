#include "WireCellUtil/Measurement.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell::Measurement;

TEST_SUITE("measurement") {

TEST_CASE("default construction and basic operations") {
    float32 f1;
    float32 f2(42);
    float64 d1;
    float64 d2(42);
    float64 d3(42, 6.9);

    CHECK(f1 == 0);
    CHECK(!f1);
    CHECK(d1 == 0);
    CHECK(!d1);
    CHECK(f2.uncertainty() == 0);
    CHECK(d2.uncertainty() == 0);

    f1 = 11;
    d1 = 12;
    d2 = d3;

    CHECK(float(f1) == 11);
    CHECK(double(d1) == 12);
}

TEST_CASE("int precision float32") {
    // not-too-large ints should be held exact
    float  x1 [[maybe_unused]] = float32(12345678);
    double x2 [[maybe_unused]] = float64(12345678);
    CHECK(x1 == x2);
}

TEST_CASE("int precision float64") {
    // even larger for doubles
    double y1 [[maybe_unused]] = float64(1234567890);
    double y2 [[maybe_unused]] = float64(1234567890);
    CHECK(y1 == y2);
}

TEST_CASE("basic arithmetic") {
    float32 a(10, 0.1);
    float32 b(20, 0.1);
    spdlog::debug("a={} b={} +:{} -:{} *:{} /:{}", a, b, a+b, a-b, a*b, a/b);
    CHECK(float(a+b) == doctest::Approx(30.0f));
    CHECK(float(a-b) == doctest::Approx(-10.0f));
    CHECK(float(a*b) == doctest::Approx(200.0f));
    CHECK(float(a/b) == doctest::Approx(0.5f));
}

}  // TEST_SUITE("measurement")
