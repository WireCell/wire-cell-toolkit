#include "WireCellUtil/Array.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include <sstream>

using namespace std;
using namespace WireCell;
using namespace WireCell::Array;

static array_xxf dup(const array_xxf& arr) { return arr; }

TEST_SUITE("array") {

TEST_CASE("copy by value") {
    const int nrows = 300;
    const int ncols = 1000;
    array_xxf arr = Eigen::ArrayXXf::Random(nrows, ncols);
    for (int ind = 0; ind < 10; ++ind) {
        array_xxf tmp = dup(arr);
        CHECK(tmp.rows() == nrows);
        CHECK(tmp.cols() == ncols);
    }
}

TEST_CASE("return by value") {
    const int nrows = 300;
    const int ncols = 1000;
    auto arr = Eigen::ArrayXXf::Random(nrows, ncols);
    REQUIRE(arr.rows() == nrows);
    REQUIRE(arr.cols() == ncols);
}

TEST_CASE("float array division nan inf zeroing") {
    array_xxf arr1(3, 2), arr2(3, 2), arr3(3, 2);
    arr1 << 0.0, 1.0, 2.0, 3.0, 4.0, 5.0;
    arr2 << 0.0, 0.5, 0.0, 2.0, 0.0, -5.0;
    arr3 = arr1 / arr2;

    {
        std::ostringstream oss;
        oss << arr3;
        spdlog::debug("arr3 before NaN zeroing\n{}", oss.str());
    }

    for (int irow = 0; irow < arr3.rows(); ++irow) {
        for (int icol = 0; icol < arr3.cols(); ++icol) {
            float val = arr3(irow, icol);
            if (std::isnan(val)) {
                arr3(irow, icol) = -0.0;
            }
            if (std::isinf(val)) {
                arr3(irow, icol) = 0.0;
            }
        }
    }

    {
        std::ostringstream oss;
        oss << arr3;
        spdlog::debug("arr3 after NaN zeroing\n{}", oss.str());
    }

    for (int irow = 0; irow < arr3.rows(); ++irow) {
        for (int icol = 0; icol < arr3.cols(); ++icol) {
            CHECK(!std::isnan(arr3(irow, icol)));
            CHECK(!std::isinf(arr3(irow, icol)));
        }
    }
    // Check specific finite values preserved
    CHECK(arr3(0, 1) == doctest::Approx(2.0f));   // 1.0/0.5
    CHECK(arr3(1, 1) == doctest::Approx(1.5f));   // 3.0/2.0
    CHECK(arr3(2, 1) == doctest::Approx(-1.0f));  // 5.0/-5.0
}

TEST_CASE("complex array division nan inf zeroing") {
    array_xxc arr1(3, 2), arr2(3, 2), arr3(3, 2);
    arr1 << 0.0, 1.0, 2.0, 3.0, 4.0, 5.0;
    arr2 << 0.0, 0.5, 0.0, 2.0, 0.0, -5.0;
    arr3 = arr1 / arr2;

    {
        std::ostringstream oss;
        oss << arr3;
        spdlog::debug("arr3 before NaN zeroing\n{}", oss.str());
    }

    for (int irow = 0; irow < arr3.rows(); ++irow) {
        for (int icol = 0; icol < arr3.cols(); ++icol) {
            float val = abs(arr3(irow, icol));
            if (std::isnan(val)) {
                arr3(irow, icol) = -0.0;
            }
            if (std::isinf(val)) {
                arr3(irow, icol) = 0.0;
            }
        }
    }

    {
        std::ostringstream oss;
        oss << arr3;
        spdlog::debug("arr3 after NaN zeroing\n{}", oss.str());
    }

    for (int irow = 0; irow < arr3.rows(); ++irow) {
        for (int icol = 0; icol < arr3.cols(); ++icol) {
            float mag = abs(arr3(irow, icol));
            CHECK(!std::isnan(mag));
            CHECK(!std::isinf(mag));
        }
    }
    // Check specific finite values preserved
    CHECK(abs(arr3(0, 1)) == doctest::Approx(2.0f));   // 1.0/0.5
    CHECK(abs(arr3(1, 1)) == doctest::Approx(1.5f));   // 3.0/2.0
    CHECK(abs(arr3(2, 1)) == doctest::Approx(1.0f));   // 5.0/-5.0
}

}  // TEST_SUITE("array")
