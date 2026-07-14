#include "WireCellUtil/Eigen.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

#include <vector>
#include <complex>

using real_t = float;
using RV = std::vector<real_t>;
using complex_t = std::complex<real_t>;
using CV = std::vector<complex_t>;
const RV real_vector{0,1,2,3,4,5};
const RV real_vector_cw{0,3,1,4,2,5};
const CV complex_vector{{0,0},{1,10},{2,20},{3,30},{4,40},{5,50}};

using RA = Eigen::Array<real_t, Eigen::Dynamic, Eigen::Dynamic>;
using RARM = Eigen::Array<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using CA = Eigen::Array<complex_t, Eigen::Dynamic, Eigen::Dynamic>;
using CARM = Eigen::Array<complex_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

TEST_SUITE("eigen cast") {

TEST_CASE("cast real to complex") {
    RA ra = Eigen::Map<RARM>((real_t*)real_vector.data(), 2, 3);
    CA ra2c = ra.cast<complex_t>();

    for (int irow = 0; irow < 2; ++irow) {
        for (int icol = 0; icol < 3; ++icol) {
            int ind = irow*3 + icol;
            complex_t c = ra2c(irow, icol);
            real_t rwant = real_vector[ind];
            spdlog::debug("{}: c={} rwant={}", ind, c.real(), rwant);
            CHECK(c.imag() == 0.0f);
            CHECK(c.real() == rwant);
        }
    }
}

TEST_CASE("get real part of complex") {
    CA ca = Eigen::Map<CARM>((complex_t*)complex_vector.data(), 2, 3);
    RA ca2r = ca.real();

    for (int irow = 0; irow < 2; ++irow) {
        for (int icol = 0; icol < 3; ++icol) {
            int ind = irow*3 + icol;
            real_t r2 = ca2r(irow, icol);
            real_t rwant = real_vector[ind];
            spdlog::debug("{}: r2={} rwant={}", ind, r2, rwant);
            CHECK(r2 == rwant);
        }
    }
}

}  // TEST_SUITE("eigen cast")
