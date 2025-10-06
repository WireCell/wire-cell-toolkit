#include "WireCellSpng/Testing.h"

#include "WireCellSpng/TorchLMN.h"
#include "WireCellUtil/Exceptions.h"

#include <cmath>
#include <complex>
#include <vector>
#include <iostream>

using namespace WireCell;
using namespace WireCell::SPNG::LMN;

namespace {
    // Helper function to create a 1D complex tensor from a std::vector
    static
    torch::Tensor create_complex_tensor(const std::vector<std::complex<float>>& data) {
        // LibTorch requires complex tensors to be created from float/double data
        // of shape [N, 2] if initialized from vector, or directly from complex vector data
        // For simplicity, we initialize from a vector of complex numbers.
        torch::Tensor in = torch::from_blob(
            (void*)data.data(), 
            {(long)data.size()}, 
            torch::kComplexFloat
        ).clone();
        return in;
    }

    // Helper function to check if two tensors are close, printing debug info if not.
    static
    bool check_tensor_close(const torch::Tensor& actual, const torch::Tensor& expected,
                            double rtol = 1e-5, double atol = 1e-8) {
        if (!torch::allclose(actual, expected, rtol, atol)) {
            std::cout << "--- TENSOR MISMATCH ---" << std::endl;
            std::cout << "Expected shape: " << expected.sizes() << std::endl;
            std::cout << "Actual shape: " << actual.sizes() << std::endl;
            std::cout << "Expected:\n" << expected << std::endl;
            std::cout << "Actual:\n" << actual << std::endl;
            
            // Check element-wise differences
            auto diff = (actual - expected).abs();
            std::cout << "Diff:\n" << diff << std::endl;
            std::cout << "Max absolute difference: " << diff.max().item<double>() << std::endl;
            //CHECK(false); // Fail the test explicitly
            return false;        // make doctest report fail site
        }
        //CHECK(true);
        return true;
    }
}


TEST_SUITE("LMN Math Helpers") {
    
    TEST_CASE("gcd") {
        // Exact cases
        CHECK(std::abs(gcd(10.0, 5.0) - 5.0) < 1e-6);
        CHECK(std::abs(gcd(12.0, 18.0) - 6.0) < 1e-6);
        
        // Rational cases
        CHECK(std::abs(gcd(0.2, 0.3) - 0.1) < 1e-6);
        CHECK(std::abs(gcd(1.2, 0.2) - 0.2) < 1e-6);

        // Zero/Near Zero input
        CHECK(std::abs(gcd(1.0, 0.0) - 1.0) < 1e-6);
        CHECK(std::abs(gcd(1.0, 1e-10, 1e-6) - 1.0) < 1e-6); // Treated as zero
    }
    
    TEST_CASE("rational") {
        // Ts = 1.0, Tr = 1.2. Ratio 5:6. dT = 0.2. gcd(1.2, 0.2) = 0.2. Ns = 1 * 1.2 / 0.2 = 6.
        CHECK(rational(1.0, 1.2) == 6);
        
        // Ts = 3.0, Tr = 2.0. Ratio 3:2. dT = 1.0. gcd(2.0, 1.0) = 1.0. Ns = 1 * 2.0 / 1.0 = 2.
        CHECK(rational(3.0, 2.0) == 2); 

        // No rational solution within tolerance (should fail or return 0 if error is caught)
        // Since the implementation uses raise<ValueError>, we check the catch block.
        CHECK_THROWS_AS(rational(M_PI, 1.0, 1e-8), ValueError);

        CHECK(rational(1.0, 2.0) == 2);
        CHECK(rational(1.0, 0.5) == 1);
    }
    
    TEST_CASE("nhalf") {
        // Odd size N=5 (DC, P1, P2, N2, N1) -> (5-1)/2 = 2
        CHECK(nhalf(5) == 2); 
        
        // Even size N=4 (DC, P1, NQ, N1) -> (4-2)/2 = 1
        CHECK(nhalf(4) == 1); 

        // Size N=2 (DC, N1) -> (2-2)/2 = 0
        CHECK(nhalf(2) == 0);
    }

    TEST_CASE("nbigger") {
        // N is already divisible by Nrat
        CHECK(nbigger(10, 5) == 10);
        
        // N needs to be increased
        CHECK(nbigger(11, 5) == 15);
        CHECK(nbigger(1, 10) == 10);

        CHECK(nbigger(20, 1) == 20);
    }
}

TEST_SUITE("LMN Resize (Time Domain)") {

    // --- Torch Tensor Implementation Tests (Real Float) ---

    TEST_CASE("1D Tensor resize upsample") {
        torch::Tensor in = torch::tensor({1.0f, 2.0f, 3.0f}, torch::kFloat);
        torch::Tensor expected = torch::tensor({1.0f, 2.0f, 3.0f, 0.0f, 0.0f}, torch::kFloat);
        torch::Tensor actual = resize(in, 5, 0);
        CHECK(check_tensor_close(actual, expected));
    }

    TEST_CASE("2D Tensor resize upsample rows (axis 0)") {
        torch::Tensor in = torch::tensor({{1.0f, 2.0f}, {3.0f, 4.0f}}, torch::kFloat);
        torch::Tensor expected = torch::tensor({{1.0f, 2.0f}, {3.0f, 4.0f}, {0.0f, 0.0f}}, torch::kFloat);
        torch::Tensor actual = resize(in, 3, 0);
        CHECK(check_tensor_close(actual, expected));
    }

    TEST_CASE("2D Tensor resize downsample columns (axis 1)") {
        torch::Tensor in = torch::tensor({{1.0f, 2.0f, 3.0f}, {4.0f, 5.0f, 6.0f}}, torch::kFloat);
        torch::Tensor expected = torch::tensor({{1.0f, 2.0f}, {4.0f, 5.0f}}, torch::kFloat);
        torch::Tensor actual = resize(in, 2, 1);
        CHECK(check_tensor_close(actual, expected));
    }
}

TEST_SUITE("LMN Resample (Frequency Domain)") {
    
    const std::complex<float> DC = {10.0f, 0.0f};
    const std::complex<float> P1 = {1.0f, 1.0f};
    const std::complex<float> NQ = {5.0f, 0.0f}; // Nyquist (only present for even N)
    const std::complex<float> N1 = {1.0f, -1.0f};
    
    // --- 1D Tensor Resampling ---

    TEST_CASE("1D Resampling Ns=4 to Nr=2 (edge case)") {
        // Ns=4: [DC, P1, NQ, N1].
        std::vector<std::complex<float>> in_data = {DC, P1, NQ, N1};
        torch::Tensor in = create_complex_tensor(in_data);
        
        // Nr=2: N_half(2)=0. Pos size = 1 (DC). Neg size = 0.
        // Result: [DC, 1]
        std::vector<std::complex<float>> expected_data = {
            DC, {1.0, 0.0f}
        };
        torch::Tensor expected = create_complex_tensor(expected_data);
        
        torch::Tensor actual = resample(in, 2, 0);
        CHECK(check_tensor_close(actual, expected));
    }
    
    // --- 2D Tensor Resampling ---

    TEST_CASE("2D UpSampling Ns=4 Nr=8 (Axis 1 - Columns)") {
        // Input: 2 rows, 4 columns (Ns=4)
        std::vector<std::complex<float>> in_row1 = {DC, P1, NQ, N1};
        std::vector<std::complex<float>> in_row2 = {DC, N1, NQ, P1}; // Different pattern for row 2
        
        torch::Tensor in = torch::empty({2, 4}, torch::kComplexFloat);
        in.narrow(0, 0, 1).copy_(create_complex_tensor(in_row1));
        in.narrow(0, 1, 1).copy_(create_complex_tensor(in_row2));

        // Nr=8.
        // Result [DC, P1, NQ/2, 0, 0, 0, NQ/2, N1] (Row 1 pattern)
        // Result [DC, N1, NQ/2, 0, 0, 0, NQ/2, P1] (Row 2 pattern)
        std::vector<std::complex<float>> exp_row1 = {
            DC, P1, {2.5f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {2.5f, 0.0f}, N1
        };
        std::vector<std::complex<float>> exp_row2 = {
            DC, N1, {2.5f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {2.5f, 0.0f}, P1
        };
        
        torch::Tensor expected = torch::empty({2, 8}, torch::kComplexFloat);
        expected.narrow(0, 0, 1).copy_(create_complex_tensor(exp_row1));
        expected.narrow(0, 1, 1).copy_(create_complex_tensor(exp_row2));
        
        torch::Tensor actual = resample(in, 8, 1);
        CHECK(check_tensor_close(actual, expected));
    }
    
    TEST_CASE("2D DownSampling Ns=8 Nr=4 (Axis 0 - Rows)") {
        // Input: 8 rows, 1 column (Ns=8)
        const std::complex<float> P2 = {2.0f, 2.0f};
        const std::complex<float> P3 = {3.0f, 3.0f};
        const std::complex<float> N2 = {2.0f, -2.0f};
        const std::complex<float> N3 = {3.0f, -3.0f};
        
        std::vector<std::complex<float>> in_data = {
            DC, P1, P2, P3, NQ, N3, N2, N1
        };
        
        // Convert to (8, 1) tensor
        torch::Tensor in = create_complex_tensor(in_data).reshape({8, 1});

        // Nr=4. N_half=1. Copy DC, P1, N1.
        std::vector<std::complex<float>> expected_data = {
            DC, P1, {2.0f, 0.0f}, N1
        };
        torch::Tensor expected = create_complex_tensor(expected_data).reshape({4, 1});
        
        torch::Tensor actual = resample(in, 4, 0);
        CHECK(check_tensor_close(actual, expected));
    }

    TEST_CASE("1D UpSampling with Nyquist Split (Ns=4 to Nr=8)") {
        // Ns=4: [DC, P1, NQ, N1]. N_half(4)=1. I_NQ_in = 2.
        // NQ is defined as {5.0f, 0.0f}. NQ/2 = {2.5f, 0.0f}
        
        std::vector<std::complex<float>> in_data = {DC, P1, NQ, N1};
        torch::Tensor in = create_complex_tensor(in_data);
        
        // Nr=8. I_NQ_in=2. I_NQ_neg_rs = 8 - 2 = 6.
        // P1 is at index 1.
        // Expected components:
        // Index 0 (DC): {10, 0}
        // Index 1 (P1): {1, 1}
        // Index 2 (P2/NQ split positive): NQ/2 = {2.5, 0}
        // Index 3, 4, 5: {0, 0}
        // Index 6 (N2/NQ split negative): NQ/2 = {2.5, 0}
        // Index 7 (N1): {1, -1}
        
        std::vector<std::complex<float>> expected_data = {
            DC, P1, {2.5f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {2.5f, 0.0f}, N1
        };
        torch::Tensor expected = create_complex_tensor(expected_data);
        
        torch::Tensor actual = resample(in, 8, 0);
        CHECK(check_tensor_close(actual, expected));
    }
    
    TEST_CASE("1D DownSampling with Nyquist Combine (Ns=8 to Nr=4)") {
        // Ns=8: [DC, P1, P2, P3, NQ, N3, N2, N1]. Nr=4. I_NQ_out = 2.
        // We use P2={2, 2} and N2={2, -2} as the folded frequencies.
        // I_P_fold = 2 (P2). I_N_fold = 8 - 2 = 6 (N2).
        
        const std::complex<float> P2 = {2.0f, 2.0f};
        const std::complex<float> N2 = {2.0f, -2.0f};
        
        std::vector<std::complex<float>> in_data = {
            DC, P1, P2, {3.0f, 3.0f}, NQ, {3.0f, -3.0f}, N2, N1
        };
        torch::Tensor in = create_complex_tensor(in_data);

        // NQ_combined = P2 + N2 = {4.0f, 0.0f}.
        // Nyquist result must be real: {4.0f, 0.0f}.
        
        // Expected: [DC, P1, P2.real, N1]
        std::vector<std::complex<float>> expected_data = {
            DC, P1, {2.0f, 0.0f}, N1
        };
        torch::Tensor expected = create_complex_tensor(expected_data);
        
        torch::Tensor actual = resample(in, 4, 0);
        CHECK(check_tensor_close(actual, expected));
    }
    
    TEST_CASE("2D Resample Downsample Ns=6 to Nr=4 Nyquist Combine Axis 1 (Cols)") {
        // Ns=6 columns. Nr=4 columns. I_NQ_out = 2.
        // Folded P freq at in[..., 2]. Folded N freq at in[..., 6 - 2 = 4].
        
        const std::complex<float> P2 = {2.0f, 1.0f};
        const std::complex<float> N2 = {2.0f, -1.0f}; 
        
        // Row 1: DC, P1, P2, P3, N2, N1 (Ns=6)
        std::vector<std::complex<float>> in_row1 = {DC, P1, P2, {3.0f, 3.0f}, N2, N1};
        
        torch::Tensor in = create_complex_tensor(in_row1).reshape({1, 6});

        // Nr=4. Expected: [DC, P1, {2.0f, 0.0f}, N1]
        std::vector<std::complex<float>> expected_data = {DC, P1, {2.0f, 0.0f}, N1};
        torch::Tensor expected = create_complex_tensor(expected_data).reshape({1, 4});
        
        torch::Tensor actual = resample(in, 4, 1);
        CHECK(check_tensor_close(actual, expected));
    }
}

// Helper to create a simple real-valued signal
static
torch::Tensor create_sine_wave(int64_t N, float frequency_cycles, torch::Dtype dtype = torch::kFloat) {
    auto t = torch::arange(N, dtype);
    auto phase = (2 * M_PI * frequency_cycles / N) * t;
    return torch::sin(phase).to(dtype);
}

TEST_SUITE("LMN Resample Interval (Composite)") {

    // Helper to verify the output size and check if the complex components were handled correctly
    bool check_resample_output(const torch::Tensor& input, const torch::Tensor& output, 
                               int64_t expected_size, int64_t axis) {
        // double bounce to return bool and get call site 
        CHECK(output.dim() == input.dim());
        if (output.dim() != input.dim()) return false;
        CHECK(output.size(axis) == expected_size);
        if (output.size(axis) != expected_size) return false;
        CHECK(output.is_floating_point());
        if (! output.is_floating_point()) return false;
        return true;
    }

    TEST_CASE("Upsampling: Ts=1.0 to Tr=0.5 (2x Up)") {
        double Ts = 1.0;
        double Tr = 0.5; // 2x upsampling (Nr = 2 * Ns)

        // 1. Create a signal: 10 samples, more than 1 full cycle.  Note, have
        // to be careful with later checks due to aliasing in the initial
        // sampling that fails to hit sample peak/valley
        int64_t Ns = 10;
        torch::Tensor input = create_sine_wave(Ns, 1.1f); // [10]
        
        // Expected intermediate rational factors:
        // rational(1.0, 0.5) = 1.0. R_size_factor = 1.
        // N_rat = nbigger(10, 1) = 10.
        // Nr = 10 * (1.0 / 0.5) = 20.

        int64_t Nr_expected = 20;
        
        torch::Tensor actual = resample_interval(input, Ts, Tr, 0);
        CHECK(check_resample_output(input, actual, Nr_expected, 0));
        
        // Check numerical closeness by comparing sampled values
        // A simple check: the upsampled signal should retain the same peak/trough values
        CHECK(torch::max(actual).item<float>() == doctest::Approx(torch::max(input).item<float>()));
        CHECK(torch::min(actual).item<float>() == doctest::Approx(torch::min(input).item<float>()));

        // Also check edge point: the 0th sample should remain 0 (sine wave start)
        CHECK(actual[0].item<float>() == doctest::Approx(0.0f).epsilon(1e-4));
    }

    TEST_CASE("Downsampling: Ts=1.0 to Tr=2.0 (2x Down)") {
        double Ts = 1.0;
        double Tr = 2.0; // 2x downsampling (Nr = 0.5 * Ns)

        // 1. Create a signal: 20 samples, at least 2 full cycles
        int64_t Ns = 20;
        float num_cycles = 2.0f;
        torch::Tensor input = create_sine_wave(Ns, num_cycles); // [20]

        // Rational factors:
        // rational(1.0, 2.0) = 1.0. R_size_factor = 1.
        // N_rat = nbigger(20, 1) = 20.
        // Nr = 20 * (1.0 / 2.0) = 10.
        
        int64_t Nr_expected = 10;
        
        torch::Tensor actual = resample_interval(input, Ts, Tr, 0);
        CHECK(check_resample_output(input, actual, Nr_expected, 0));

        // Downsampling might cause aliasing if not carefully chosen, 
        // but since we used 2 cycles in 20 steps, the resulting 10 steps 
        // should capture exactly 1 cycle (perfect downsampling).
        
        // Verify output frequency (1 cycle in 10 steps)
        torch::Tensor expected_1cycle = create_sine_wave(10, num_cycles);
        CHECK(check_tensor_close(actual, expected_1cycle, 1e-2, 1e-5));
    }
    
    TEST_CASE("Non-Integer Ratio Upsampling: Ts=3.0 to Tr=2.0 (1.5x)") {
        double Ts = 3.0;
        double Tr = 2.0; // 1.5x upsampling (Nr = 1.5 * Ns)

        // Ns = 10.
        int64_t Ns = 10;
        float num_cycles = 1.0f;
        torch::Tensor input = create_sine_wave(Ns, num_cycles);

        // Rational factors:
        // rational(3.0, 2.0) = 2. R_size_factor = 2.
        // N_rat = nbigger(10, 2) = 10.
        // Nr = 10 * (3.0 / 2.0) = 15.
        
        int64_t Nr_expected = 15;
        
        torch::Tensor actual = resample_interval(input, Ts, Tr, 0);
        CHECK(check_resample_output(input, actual, Nr_expected, 0));

        // std::cout << "Input:\n" << input << "\n";
        // std::cout << "Resampled:\n" << actual << "\n";

        // Peak will not be sampled in both cases so the comparison has to be sloppy.
        CHECK(torch::max(actual).item<float>() == doctest::Approx(torch::max(input).item<float>()).epsilon(1e-1));
    }

    TEST_CASE("Need for Rational Padding: Ns=9, Ts=3.0 to Tr=2.0") {
        double Ts = 3.0;
        double Tr = 2.0; 

        // Ns = 9.
        int64_t Ns = 9;
        float num_cycles = 1.0f;
        torch::Tensor input = create_sine_wave(Ns, num_cycles);

        // Rational factors:
        // R_size_factor = 2.
        // N_rat = nbigger(9, 2) = 10. (Padding required)
        // Nr = 10 * (3.0 / 2.0) = 15.
        
        //int64_t N_rat_expected = 10;
        int64_t Nr_expected = 15;
        
        // We cannot easily test N_rat directly, but we can verify the final size.
        torch::Tensor actual = resample_interval(input, Ts, Tr, 0);
        CHECK(check_resample_output(input, actual, Nr_expected, 0));
        
        // std::cout << "Input:\n" << input << "\n";
        // std::cout << "Resampled:\n" << actual << "\n";

        // Peak will not be sampled in both cases so the comparison has to be sloppy.
        CHECK(torch::max(actual).item<float>() == doctest::Approx(torch::max(input).item<float>()).epsilon(1e-1));
    }
}

