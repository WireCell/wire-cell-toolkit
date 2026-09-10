// Regression tests for the high level 1D helpers
// DftTools::convolve() and DftTools::replace().
//
// Issue #531: DftTools::replace() divided by the untransformed, unpadded
// time-domain response samples (reading past their ends), had the
// new/old ratio inverted and never applied the inverse transform.  This
// corrupted every channel run through Gen::PerChannelVariation (the
// MicroBooNE overlay simulation) and the selected channels run through
// Gen::Misconfigure.  DftTools::convolve() had the same missing inverse
// transform.
//
// These tests compare against direct O(N*M) linear convolution so they
// fail on the pre-fix implementation and pass on the fixed one.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellAux/DftTools.h"
#include "WireCellIface/IDFT.h"

#include <algorithm>
#include <cmath>
#include <random>

using namespace WireCell;
using namespace WireCell::Aux;

using RV = DftTools::real_vector_t;

static IDFT::pointer get_dft()
{
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellAux");
    return Factory::lookup_tn<IDFT>("FftwDFT");
}

// Direct linear convolution, size a.size()+b.size()-1.
static RV linear_convolve(const RV& a, const RV& b)
{
    if (a.empty() || b.empty()) return RV{};
    RV out(a.size() + b.size() - 1, 0.0f);
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            out[i + j] += a[i] * b[j];
        }
    }
    return out;
}

static RV random_wave(size_t n, unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    RV v(n);
    for (auto& x : v) x = dist(rng);
    return v;
}

// A causal, exponentially decaying "response" of n samples.  Its DFT
// is a geometric sum that never vanishes so deconvolution by it is
// well-conditioned.
static RV decay_response(size_t n, float tau, float scale)
{
    RV v(n);
    for (size_t i = 0; i < n; ++i) v[i] = scale * std::exp(-float(i) / tau);
    return v;
}

static float max_abs(const RV& v)
{
    float m = 0;
    for (auto x : v) m = std::max(m, std::abs(x));
    return m;
}

static bool all_finite(const RV& v)
{
    return std::all_of(v.begin(), v.end(), [](float x) { return std::isfinite(x); });
}

// Check that got matches want (over want's length) to within a
// tolerance relative to the largest expected value, and that any
// samples beyond want's length are compatible with zero.
static void check_close(const RV& got, const RV& want, float rel_tol)
{
    REQUIRE(got.size() >= want.size());
    REQUIRE(all_finite(got));
    const float tol = rel_tol * max_abs(want);
    for (size_t i = 0; i < want.size(); ++i) {
        INFO("index " << i);
        CHECK(std::abs(got[i] - want[i]) <= tol);
    }
    for (size_t i = want.size(); i < got.size(); ++i) {
        INFO("tail index " << i);
        CHECK(std::abs(got[i]) <= tol);
    }
}

TEST_SUITE("DftTools high level 1D helpers") {

    TEST_CASE("convolve matches direct linear convolution") {
        auto dft = get_dft();
        const RV a = random_wave(37, 1u);
        const RV b = random_wave(11, 2u);

        const RV got = DftTools::convolve(dft, a, b);
        const RV want = linear_convolve(a, b);

        CHECK(got.size() == a.size() + b.size() - 1);
        check_close(got, want, 1e-4f);
    }

    TEST_CASE("convolve is commutative") {
        auto dft = get_dft();
        const RV a = random_wave(20, 3u);
        const RV b = random_wave(33, 4u);
        const RV ab = DftTools::convolve(dft, a, b);
        const RV ba = DftTools::convolve(dft, b, a);
        check_close(ab, ba, 1e-5f);
    }

    TEST_CASE("replace swaps old response for new response") {
        auto dft = get_dft();

        // Underlying "true" signal: a few impulses of differing sign.
        RV signal(100, 0.0f);
        signal[10] = 1.0f;
        signal[40] = -2.5f;
        signal[41] = -2.5f;
        signal[77] = 0.7f;

        // Responses of differing length and shape, in the caller's
        // (Gen::Misconfigure / Gen::PerChannelVariation) convention:
        // the measurement was made with res_old, we want res_new.
        const RV res_old = decay_response(20, 5.0f, 1.0f);
        const RV res_new = decay_response(25, 3.0f, 2.0f);

        const RV meas = linear_convolve(signal, res_old);
        const RV want = linear_convolve(signal, res_new);

        const RV got = DftTools::replace(dft, meas, res_new, res_old);

        // Size rule shared with legacy Waveform::replace_convolve().
        const size_t sizes[3] = {meas.size(), res_new.size(), res_old.size()};
        const size_t expected_size = sizes[0] + sizes[1] + sizes[2]
            - *std::min_element(sizes, sizes + 3) - 1;
        CHECK(got.size() == expected_size);

        // Deconvolution amplifies round-off; 1e-3 relative is still
        // orders of magnitude below the pre-fix garbage.
        check_close(got, want, 1e-3f);
    }

    TEST_CASE("replace with identical responses is the identity") {
        auto dft = get_dft();
        const RV meas = random_wave(64, 5u);
        const RV res = decay_response(16, 4.0f, 1.5f);

        const RV got = DftTools::replace(dft, meas, res, res);
        check_close(got, meas, 1e-4f);
    }

    TEST_CASE("replace with a truncated measurement, as callers use it") {
        // Callers do wave.resize(charge.size()) after replace().  The
        // leading charge.size() samples must match the expectation to
        // within the truncation error of the tail.
        auto dft = get_dft();
        RV signal(200, 0.0f);
        signal[50] = 1.0f;
        const RV res_old = decay_response(30, 6.0f, 1.0f);
        const RV res_new = decay_response(30, 2.0f, 3.0f);

        RV meas = linear_convolve(signal, res_old);
        meas.resize(signal.size());          // as a fixed-length trace
        RV want = linear_convolve(signal, res_new);
        want.resize(signal.size());

        RV got = DftTools::replace(dft, meas, res_new, res_old);
        REQUIRE(got.size() >= signal.size());
        got.resize(signal.size());
        check_close(got, want, 1e-3f);
    }

    TEST_CASE("replace with an all-zero old response is finite and unscaled") {
        // Policy: bins where FFT(res_old) is exactly zero are left
        // unscaled.  An all-zero response gives an exactly-zero spectrum
        // in every bin, so the output must equal the (padded) input and
        // contain no inf/NaN.
        auto dft = get_dft();
        const RV meas = random_wave(48, 6u);
        const RV res_new = decay_response(8, 2.0f, 1.0f);
        const RV res_old(8, 0.0f);

        const RV got = DftTools::replace(dft, meas, res_new, res_old);
        check_close(got, meas, 1e-5f);
    }

    TEST_CASE("replace with an all-zero new response gives zero") {
        auto dft = get_dft();
        const RV meas = random_wave(48, 7u);
        const RV res_new(8, 0.0f);
        const RV res_old = decay_response(8, 2.0f, 1.0f);

        const RV got = DftTools::replace(dft, meas, res_new, res_old);
        REQUIRE(all_finite(got));
        CHECK(max_abs(got) <= 1e-5f * max_abs(meas));
    }
}
