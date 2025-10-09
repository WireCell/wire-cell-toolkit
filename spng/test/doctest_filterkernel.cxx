#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Binning.h"       // For calculating frequencies
#include "WireCellUtil/HanaJsonCPP.h"  // For converting Configuration to/from struct
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Exceptions.h"

#include "WireCellSpng/FilterKernel.h"
#include "WireCellSpng/Testing.h"


#include <cmath>      // For std::exp, std::pow, std::abs
#include <functional> // For std::function
#include <vector>
#include <limits>     // For std::numeric_limits

// Helper for float/double comparison with tolerance
// torch::kFloat corresponds to float, so compare floats.
bool almost_equal_float(float a, float b, float epsilon = 1e-6) {
    return std::abs(a - b) < epsilon;
}

// Helper to define the analytical filter functions
std::function<float(float, double, double)> get_lowpass_func =
    [](float freq, double scale, double power) {
    return std::exp(-0.5 * std::pow(freq / scale, power));
};

std::function<float(float, double, double)> get_highpass_func =
    [](float freq, double scale, double power) {
    return 1.0f - std::exp(-std::pow(freq / scale, power));
};

// Needed to find aux components.  I thought this library was linked in by
// default but apparently not.  Or, not for wcdoctest-*?
static int init_this_crap()
{
    WireCell::PluginManager::instance().add("WireCellAux");
    return 0;
}
static int crap_initialized = init_this_crap();


TEST_SUITE("spng filter kernel") {


    using namespace WireCell::SPNG;
    using WireCell::Configuration;
    using WireCell::ValueError;
    using WireCell::Binning;
    

    // Test default constructor and configure method
    TEST_CASE("spng filter kernel default ctor and configure") {

        FilterKernel fk;

        FilterKernelConfig cfg {{ {"lowpass", 0.5, 0.2, 2.0, true} }};
        REQUIRE(cfg.axis.size() > 0);

        Configuration wc_cfg = WireCell::HanaJsonCPP::to_json(cfg);
        std::cerr << wc_cfg << "\n";
        fk.configure(wc_cfg);

        // Verify some properties by getting a spectrum
        WireCell::ITorchSpectrum::shape_t shape = {100};
        torch::Tensor spec = fk.spectrum(shape);
        
        REQUIRE(spec.numel() == 100);
        REQUIRE(spec.dtype() == torch::kFloat);

        // Check the zero frequency bin with ignore_baseline=true
        CHECK(almost_equal_float(spec[0].item<float>(), 0.0f));

        // Check a mid-range frequency
        const auto& axis = cfg.axis[0];
        const double period = axis.period;
        const double scale = axis.scale;
        const double power = axis.power;
        const size_t nbins = shape[0];
        Binning bins(nbins, 0, 1.0/period); // Frequencies from 0 to fs

        size_t ind = nbins / 4; // Mid-point in frequency range
        float freq = bins.edge(ind);
        float expected_val = get_lowpass_func(freq, scale, power);
        CHECK(almost_equal_float(spec[ind].item<float>(), expected_val));
    }

// Test constructor with FilterKernelConfig struct
    TEST_CASE("spng filter kernel config ctor") {
        FilterKernelConfig cfg {{ {"highpass", 0.01, 10.0, 3.0, false} }};
        FilterKernel fk(cfg);

        // Verify some properties by getting a spectrum
        WireCell::ITorchSpectrum::shape_t shape = {50};
        torch::Tensor spec = fk.spectrum(shape);

        REQUIRE(spec.numel() == 50);
        REQUIRE(spec.dtype() == torch::kFloat);

        // Check the zero frequency bin with ignore_baseline=false
        const double period = cfg.axis[0].period;
        const double scale = cfg.axis[0].scale;
        const double power = cfg.axis[0].power;
        float expected_dc_val = get_highpass_func(0.0f, scale, power); // Should be 0.0 for highpass
        CHECK(almost_equal_float(spec[0].item<float>(), expected_dc_val));

        // Check a mid-range frequency
        const size_t nbins = shape[0];
        Binning bins(nbins, 0, 1.0/period); // Frequencies from 0 to fs

        size_t ind = nbins / 3; // Another mid-point
        float freq = bins.edge(ind);
        float expected_val = get_highpass_func(freq, scale, power);
        CHECK(almost_equal_float(spec[ind].item<float>(), expected_val));
    }

// Test configure method with invalid values
    TEST_CASE("spng filter kernel invalid configurations") {
        FilterKernelConfig base_cfg = {{ {"lowpass", 1.0, 1.0, 2.0} }};

        // Invalid kind
        FilterKernelConfig invalid_kind_cfg = base_cfg;
        invalid_kind_cfg.axis[0].kind = "unknown";
        Configuration wc_cfg_kind = WireCell::HanaJsonCPP::to_json(invalid_kind_cfg);
        FilterKernel fk_kind;
        CHECK_THROWS_AS(fk_kind.configure(wc_cfg_kind), ValueError);

        // Invalid period (non-positive)
        FilterKernelConfig invalid_period_cfg = base_cfg;
        invalid_period_cfg.axis[0].period = 0.0;
        Configuration wc_cfg_period = WireCell::HanaJsonCPP::to_json(invalid_period_cfg);
        FilterKernel fk_period;
        CHECK_THROWS_AS(fk_period.configure(wc_cfg_period), ValueError);

        invalid_period_cfg.axis[0].period = -1.0;
        wc_cfg_period = WireCell::HanaJsonCPP::to_json(invalid_period_cfg);
        CHECK_THROWS_AS(fk_period.configure(wc_cfg_period), ValueError);


        // Invalid scale (non-positive)
        FilterKernelConfig invalid_scale_cfg = base_cfg;
        invalid_scale_cfg.axis[0].scale = 0.0;
        Configuration wc_cfg_scale = WireCell::HanaJsonCPP::to_json(invalid_scale_cfg);
        FilterKernel fk_scale;
        CHECK_THROWS_AS(fk_scale.configure(wc_cfg_scale), ValueError);

        invalid_scale_cfg.axis[0].scale = -1.0;
        wc_cfg_scale = WireCell::HanaJsonCPP::to_json(invalid_scale_cfg);
        CHECK_THROWS_AS(fk_scale.configure(wc_cfg_scale), ValueError);
    }

// Test default_configuration method
    TEST_CASE("spng filter kernel default configuration values") {
        FilterKernel fk; // Default constructed

        Configuration default_cfg = fk.default_configuration();

        FilterKernelConfig recovered_cfg;
        WireCell::HanaJsonCPP::from_json(recovered_cfg, default_cfg);

        // Actual default config is empty.
        REQUIRE(recovered_cfg.axis.empty());
        // fake having one axis 
        recovered_cfg.axis.resize(1);

        // Check if recovered config matches the struct's default values
        // Note: initial FilterKernelConfig defaults are:
        // kind{""}, period=-1.0, scale=-1.0, power=2.0, ignore_baseline=true
        CHECK(recovered_cfg.axis[0].kind == "");
        CHECK(almost_equal_float(recovered_cfg.axis[0].period, -1.0f));
        CHECK(almost_equal_float(recovered_cfg.axis[0].scale, -1.0f));
        CHECK(almost_equal_float(recovered_cfg.axis[0].power, 2.0f));
        CHECK(recovered_cfg.axis[0].ignore_baseline == true);
    }

// Test spectrum generation for lowpass filter
    TEST_CASE("spng filter kernel spectrum lowpass filter details") {
        // fs = 2.0, Nyquist = 1.0
        FilterKernelConfig cfg = {{ {"lowpass", 0.5, 0.2, 2.0, true} }};

        FilterKernel fk(cfg);
        WireCell::ITorchSpectrum::shape_t shape = {100}; // Even number of bins
        torch::Tensor spec = fk.spectrum(shape);

        REQUIRE(spec.numel() == 100);
        REQUIRE(spec.dtype() == torch::kFloat);

        const double period = cfg.axis[0].period;
        const double scale = cfg.axis[0].scale;
        const double power = cfg.axis[0].power;
        const size_t nbins = shape[0];
        Binning bins(nbins, 0, 1.0/period); // Frequencies from 0 to fs

        // Check DC (zero frequency) with ignore_baseline=true
        CHECK(almost_equal_float(spec[0].item<float>(), 0.0f));

        // Check symmetry and values at various points
        for (size_t i = 1; i < nbins / 2; ++i) {
            float freq = bins.edge(i);
            float expected_val = get_lowpass_func(freq, scale, power);
            CHECK(almost_equal_float(spec[i].item<float>(), expected_val));
            CHECK(almost_equal_float(spec[nbins - i].item<float>(), expected_val));
        }

    }

// Test spectrum generation for highpass filter
    TEST_CASE("spng filter kernel spectrum highpass filter details") {
        // fs = 100.0, Nyquist = 50.0
        FilterKernelConfig cfg = {{ {"highpass", 0.01, 10.0, 3.0, true} }};

        FilterKernel fk(cfg);
        WireCell::ITorchSpectrum::shape_t shape = {51}; // Odd number of bins
        torch::Tensor spec = fk.spectrum(shape);

        REQUIRE(spec.numel() == 51);
        REQUIRE(spec.dtype() == torch::kFloat);

        const double period = cfg.axis[0].period;
        const double scale = cfg.axis[0].scale;
        const double power = cfg.axis[0].power;
        const size_t nbins = shape[0];
        Binning bins(nbins, 0, 1.0/period); // Frequencies from 0 to fs

        // Check DC (zero frequency) with ignore_baseline=true
        CHECK(almost_equal_float(spec[0].item<float>(), 0.0f));

        // Check symmetry and values at various points
        for (size_t i = 1; i <= nbins / 2; ++i) { // Loop up to half for odd nbins
            float freq = bins.edge(i);
            float expected_val = get_highpass_func(freq, scale, power);
            CHECK(almost_equal_float(spec[i].item<float>(), expected_val));
            CHECK(almost_equal_float(spec[nbins - i].item<float>(), expected_val));
        }
    }

// Test ignore_baseline option
    TEST_CASE("spng filter kernel spectrum ignore baseline option") {
        FilterKernelConfig cfg = {{ {"lowpass", 1.0, 0.5, 2.0} }};
        
        WireCell::ITorchSpectrum::shape_t shape = {20};
        const double scale = cfg.axis[0].scale;
        const double power = cfg.axis[0].power;

        // Case 1: ignore_baseline = true (default)
        cfg.axis[0].ignore_baseline = true;
        FilterKernel fk_true(cfg);
        torch::Tensor spec_true = fk_true.spectrum(shape);
        CHECK(almost_equal_float(spec_true[0].item<float>(), 0.0f));

        // Case 2: ignore_baseline = false
        cfg.axis[0].ignore_baseline = false;
        FilterKernel fk_false(cfg);
        torch::Tensor spec_false = fk_false.spectrum(shape);
        float expected_dc_lowpass = get_lowpass_func(0.0f, scale, power); // Should be 1.0 for lowpass
        CHECK(almost_equal_float(spec_false[0].item<float>(), expected_dc_lowpass));

        // For highpass, DC is always 0.0, so ignore_baseline shouldn't change s[0]
        cfg.axis[0].kind = "highpass";
        float expected_dc_highpass = get_highpass_func(0.0f, scale, power); // Should be 0.0
        
        cfg.axis[0].ignore_baseline = true;
        FilterKernel fk_hp_true(cfg);
        torch::Tensor spec_hp_true = fk_hp_true.spectrum(shape);
        CHECK(almost_equal_float(spec_hp_true[0].item<float>(), 0.0f)); // Always 0 for highpass DC

        cfg.axis[0].ignore_baseline = false;
        FilterKernel fk_hp_false(cfg);
        torch::Tensor spec_hp_false = fk_hp_false.spectrum(shape);
        CHECK(almost_equal_float(spec_hp_false[0].item<float>(), expected_dc_highpass)); // Still 0
    }

// Test that only 1D shapes are supported
    TEST_CASE("spng filter kernel spectrum shape constraint") {
        FilterKernelConfig cfg = {{ {"lowpass", 1.0, 1.0} }};
        
        FilterKernel fk(cfg);

        // 0D shape
        WireCell::ITorchSpectrum::shape_t shape0 = {};
        CHECK_THROWS_AS(fk.spectrum(shape0), ValueError);

        // 2D shape
        WireCell::ITorchSpectrum::shape_t shape2 = {10, 20};
        CHECK_THROWS_AS(fk.spectrum(shape2), ValueError);

        // Valid 1D shape (should not throw)
        WireCell::ITorchSpectrum::shape_t shape1 = {50};
        CHECK_NOTHROW(fk.spectrum(shape1));
    }

}
