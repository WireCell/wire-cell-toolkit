// Regression pin for the in-place 2D-FFT overloads added in commit
// 4e37cd83 (aux/DftTools+sigproc/OmnibusSigProc: skip redundant 2D-array
// copies in FFT path).
//
// The new fwd_inplace()/inv_inplace() overloads must produce bit-identical
// output to the existing copy-returning fwd()/inv() overloads, for both
// axes and both directions.  They share the same FFTW plan signature
// (same dims, direction, axis, in-place flag) so the cached plan is the
// same and butterflies execute on the same memory in the same order.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellAux/DftTools.h"
#include "WireCellIface/IDFT.h"

#include <random>

using namespace WireCell;
using namespace WireCell::Aux;

static DftTools::complex_array_t random_array(int nrows, int ncols, unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    DftTools::complex_array_t arr(nrows, ncols);
    for (int r = 0; r < nrows; ++r) {
        for (int c = 0; c < ncols; ++c) {
            arr(r, c) = DftTools::complex_t(dist(rng), dist(rng));
        }
    }
    return arr;
}

static bool bit_equal(const DftTools::complex_array_t& a,
                      const DftTools::complex_array_t& b)
{
    if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
    const auto* pa = a.data();
    const auto* pb = b.data();
    const ptrdiff_t n = a.rows() * a.cols();
    for (ptrdiff_t i = 0; i < n; ++i) {
        if (pa[i].real() != pb[i].real()) return false;
        if (pa[i].imag() != pb[i].imag()) return false;
    }
    return true;
}

TEST_CASE("DftTools fwd/inv in-place equals allocating overload") {
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellAux");
    auto idft = Factory::lookup_tn<IDFT>("FftwDFT");

    const int nrows = 64;
    const int ncols = 128;

    SUBCASE("fwd axis=0") {
        auto a = random_array(nrows, ncols, 1u);
        auto b = a;
        auto a_out = DftTools::fwd(idft, a, 0);
        DftTools::fwd_inplace(idft, b, 0);
        CHECK(bit_equal(a_out, b));
    }
    SUBCASE("fwd axis=1") {
        auto a = random_array(nrows, ncols, 2u);
        auto b = a;
        auto a_out = DftTools::fwd(idft, a, 1);
        DftTools::fwd_inplace(idft, b, 1);
        CHECK(bit_equal(a_out, b));
    }
    SUBCASE("inv axis=0") {
        auto a = random_array(nrows, ncols, 3u);
        auto b = a;
        auto a_out = DftTools::inv(idft, a, 0);
        DftTools::inv_inplace(idft, b, 0);
        CHECK(bit_equal(a_out, b));
    }
    SUBCASE("inv axis=1") {
        auto a = random_array(nrows, ncols, 4u);
        auto b = a;
        auto a_out = DftTools::inv(idft, a, 1);
        DftTools::inv_inplace(idft, b, 1);
        CHECK(bit_equal(a_out, b));
    }
    SUBCASE("fwd then inv recovers original (in-place)") {
        auto a = random_array(nrows, ncols, 5u);
        auto orig = a;
        DftTools::fwd_inplace(idft, a, 1);
        DftTools::inv_inplace(idft, a, 1);
        // 1/N normalization on inv => recover orig within float epsilon
        const auto* pa = a.data();
        const auto* po = orig.data();
        const ptrdiff_t n = a.rows() * a.cols();
        for (ptrdiff_t i = 0; i < n; ++i) {
            CHECK(std::abs(pa[i].real() - po[i].real()) < 1e-4);
            CHECK(std::abs(pa[i].imag() - po[i].imag()) < 1e-4);
        }
    }
}
