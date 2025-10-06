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
    void check_tensor_close(const torch::Tensor& actual, const torch::Tensor& expected, double atol = 1e-6) {
        if (!torch::allclose(actual, expected, atol)) {
            std::cout << "--- TENSOR MISMATCH ---" << std::endl;
            std::cout << "Expected shape: " << expected.sizes() << std::endl;
            std::cout << "Actual shape: " << actual.sizes() << std::endl;
            std::cout << "Expected:\n" << expected << std::endl;
            std::cout << "Actual:\n" << actual << std::endl;
            
            // Check element-wise differences
            auto diff = (actual - expected).abs();
            std::cout << "Max absolute difference: " << diff.max().item<double>() << std::endl;
            CHECK(false); // Fail the test explicitly
        } else {
            CHECK(true);
        }
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
    }
}

TEST_SUITE("LMN Resize (Time Domain)") {

    // --- Vector Implementation Tests ---
    TEST_CASE("vector resize upsample") {
        std::vector<float> in = {1.0f, 2.0f, 3.0f};
        std::vector<float> expected = {1.0f, 2.0f, 3.0f, 0.0f, 0.0f};
        CHECK(resize(in, 5) == expected);
    }

    TEST_CASE("vector resize downsample") {
        std::vector<float> in = {1.0f, 2.0f, 3.0f, 4.0f};
        std::vector<float> expected = {1.0f, 2.0f};
        CHECK(resize(in, 2) == expected);
    }

    // --- Torch Tensor Implementation Tests (Real Float) ---

    TEST_CASE("1D Tensor resize upsample") {
        torch::Tensor in = torch::tensor({1.0f, 2.0f, 3.0f}, torch::kFloat);
        torch::Tensor expected = torch::tensor({1.0f, 2.0f, 3.0f, 0.0f, 0.0f}, torch::kFloat);
        torch::Tensor actual = resize(in, 5, 0);
        check_tensor_close(actual, expected);
    }

    TEST_CASE("2D Tensor resize upsample rows (axis 0)") {
        torch::Tensor in = torch::tensor({{1.0f, 2.0f}, {3.0f, 4.0f}}, torch::kFloat);
        torch::Tensor expected = torch::tensor({{1.0f, 2.0f}, {3.0f, 4.0f}, {0.0f, 0.0f}}, torch::kFloat);
        torch::Tensor actual = resize(in, 3, 0);
        check_tensor_close(actual, expected);
    }

    TEST_CASE("2D Tensor resize downsample columns (axis 1)") {
        torch::Tensor in = torch::tensor({{1.0f, 2.0f, 3.0f}, {4.0f, 5.0f, 6.0f}}, torch::kFloat);
        torch::Tensor expected = torch::tensor({{1.0f, 2.0f}, {4.0f, 5.0f}}, torch::kFloat);
        torch::Tensor actual = resize(in, 2, 1);
        check_tensor_close(actual, expected);
    }
}

TEST_SUITE("LMN Resample (Frequency Domain)") {
    
    const std::complex<float> DC = {10.0f, 0.0f};
    const std::complex<float> P1 = {1.0f, 1.0f};
    const std::complex<float> NQ = {5.0f, 0.0f}; // Nyquist (only present for even N)
    const std::complex<float> N1 = {1.0f, -1.0f};
    
    // --- 1D Tensor Resampling ---

    TEST_CASE("1D UpSampling (Ns=4 to Nr=8)") {
        // Ns=4: [DC, P1, NQ, N1]. N_half(4)=1.
        std::vector<std::complex<float>> in_data = {DC, P1, NQ, N1};
        torch::Tensor in = create_complex_tensor(in_data);

        // Nr=8: N_half(4)=1. Copy DC (1), P1 (1). Copy N1 (1). Total 3 components preserved.
        // Result: [DC, P1, 0, 0, 0, 0, 0, N1]
        std::vector<std::complex<float>> expected_data = {
            DC, P1, {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, N1
        };
        torch::Tensor expected = create_complex_tensor(expected_data);
        
        torch::Tensor actual = resample(in, 8, 0);
        check_tensor_close(actual, expected);
    }

    TEST_CASE("1D DownSampling (Ns=8 to Nr=4)") {
        // Ns=8: N_half(8)=3. We use [DC, P1, P2, P3, NQ, N_3, N_2, N_1]
        // Define P2 and P3 for input
        const std::complex<float> P2 = {2.0f, 2.0f};
        const std::complex<float> P3 = {3.0f, 3.0f};
        const std::complex<float> N2 = {2.0f, -2.0f};
        const std::complex<float> N3 = {3.0f, -3.0f};
        
        std::vector<std::complex<float>> in_data = {
            DC, P1, P2, P3, NQ, N3, N2, N1
        };
        torch::Tensor in = create_complex_tensor(in_data);

        // Nr=4: N_half(4)=1. Copy DC (1), P1 (1). Copy N1 (1). Total 3 components preserved.
        // Result: [DC, P1, 0, N1] -- Wait, NQ is index 2 for N=4.
        // Result layout for Nr=4: [DC, P1, NQ_new, N1]
        // N_half=1. Pos size = 2 (DC, P1). Neg size = 1 (N1). 
        // Index 2 is padded zero because we skip the original Nyquist/high freqs.
        
        std::vector<std::complex<float>> expected_data = {
            DC, P1, {0.0f, 0.0f}, N1 // P1 is copied, P2, P3 dropped. N3, N2 dropped. N1 copied.
        };
        torch::Tensor expected = create_complex_tensor(expected_data);
        
        torch::Tensor actual = resample(in, 4, 0);
        check_tensor_close(actual, expected);
    }
    
    TEST_CASE("1D Resampling to Nr=2 (edge case)") {
        // Ns=4: [DC, P1, NQ, N1].
        std::vector<std::complex<float>> in_data = {DC, P1, NQ, N1};
        torch::Tensor in = create_complex_tensor(in_data);
        
        // Nr=2: N_half(2)=0. Pos size = 1 (DC). Neg size = 0.
        // Result: [DC, 0]
        std::vector<std::complex<float>> expected_data = {
            DC, {0.0f, 0.0f}
        };
        torch::Tensor expected = create_complex_tensor(expected_data);
        
        torch::Tensor actual = resample(in, 2, 0);
        check_tensor_close(actual, expected);
    }
    
    // --- 2D Tensor Resampling ---

    TEST_CASE("2D UpSampling (Axis 1 - Columns)") {
        // Input: 2 rows, 4 columns (Ns=4)
        std::vector<std::complex<float>> in_row1 = {DC, P1, NQ, N1};
        std::vector<std::complex<float>> in_row2 = {DC, N1, NQ, P1}; // Different pattern for row 2
        
        torch::Tensor in = torch::empty({2, 4}, torch::kComplexFloat);
        in.narrow(0, 0, 1).copy_(create_complex_tensor(in_row1));
        in.narrow(0, 1, 1).copy_(create_complex_tensor(in_row2));

        // Nr=8. Result [DC, P1, 0, 0, 0, 0, 0, N1] (Row 1 pattern)
        // Result [DC, N1, 0, 0, 0, 0, 0, P1] (Row 2 pattern)
        std::vector<std::complex<float>> exp_row1 = {
            DC, P1, {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, N1
        };
        std::vector<std::complex<float>> exp_row2 = {
            DC, N1, {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, P1
        };
        
        torch::Tensor expected = torch::empty({2, 8}, torch::kComplexFloat);
        expected.narrow(0, 0, 1).copy_(create_complex_tensor(exp_row1));
        expected.narrow(0, 1, 1).copy_(create_complex_tensor(exp_row2));
        
        torch::Tensor actual = resample(in, 8, 1);
        check_tensor_close(actual, expected);
    }
    
    TEST_CASE("2D DownSampling (Axis 0 - Rows)") {
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
            DC, P1, {0.0f, 0.0f}, N1
        };
        torch::Tensor expected = create_complex_tensor(expected_data).reshape({4, 1});
        
        torch::Tensor actual = resample(in, 4, 0);
        check_tensor_close(actual, expected);
    }

}
