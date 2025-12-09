#include "WireCellSpng/Rebaseline.h"
#include "WireCellSpng/Testing.h"
#include "WireCellSpng/Util.h"

#include <iostream>

using namespace WireCell::SPNG;


// Helper function to check if two tensors are close within a tolerance
static
bool tensors_are_close(const torch::Tensor& a, const torch::Tensor& b, double tolerance = 1e-4) {
    if (a.sizes().vec() != b.sizes().vec()) {
        spdlog::error("Tensor size mismatch: {} vs {}",
                      to_string(a), to_string(b));
        return false;
    }
    const bool ok = a.allclose(b, tolerance, tolerance);
    if (!ok) {
        spdlog::error("tensor mismatch:\n{}\n{}",
                      to_string(a, true), to_string(b, true));
    }
    return ok;
}

// --- Test Suite ---
TEST_SUITE("spng rebaseline") {

    TEST_CASE("1D_Simple_Fragment") {
        spdlog::debug("--- Testing 1D Simple Fragment ---");

        // [0, 0, 5, 10, 5, 0, 0]
        // Fragment: [5, 10, 5] (indices 2 to 4)
        // Baseline model connects (index 2, value 5) to (index 4, value 5)
        // Line: y = 5 (A=0, B=5)
        torch::Tensor input = torch::tensor({0.0f, 0.0f, 5.0f, 10.0f, 5.0f, 0.0f, 0.0f}, torch::kFloat32);
        float threshold = 2.0f; // Threshold for fragment detection

        torch::Tensor expected_output = torch::tensor({0.0f, 0.0f, 0.0f, 5.0f, 0.0f, 0.0f, 0.0f}, torch::kFloat32);
        // Calculation: 
        // Index 2: 5 - 5 = 0
        // Index 3: 10 - 5 = 5
        // Index 4: 5 - 5 = 0
        
        torch::Tensor actual_output = rebaseline(input, -1, threshold);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Expected: {}", to_string(expected_output));
        spdlog::debug("Actual: {}", to_string(actual_output));
        
        CHECK(tensors_are_close(actual_output, expected_output));
    }

    TEST_CASE("1D_Sloped_Baseline") {
        spdlog::debug("--- Testing 1D Sloped Baseline ---");

        // [1, 2, 8, 10, 6, 0]
        // Fragment: [8, 10, 6] (indices 2 to 4). Threshold 5.0
        // Baseline model connects (index 2, value 8) to (index 4, value 6)
        // Length = 4 - 2 = 2
        // Slope A = (6 - 8) / 2 = -1.0
        // Intercept B = 8 - (-1.0) * 2 = 10.0
        // Line: y = -1.0*x + 10.0
        // Expected Baseline values:
        // x=2: -2 + 10 = 8.0
        // x=3: -3 + 10 = 7.0
        // x=4: -4 + 10 = 6.0
        
        torch::Tensor input = torch::tensor({1.0f, 2.0f, 8.0f, 10.0f, 6.0f, 0.0f}, torch::kFloat32);
        float threshold = 5.0f; 

        // Expected output:
        // Index 2: 8.0 - 8.0 = 0.0
        // Index 3: 10.0 - 7.0 = 3.0
        // Index 4: 6.0 - 6.0 = 0.0
        torch::Tensor expected_output = torch::tensor({1.0f, 2.0f, 0.0f, 3.0f, 0.0f, 0.0f}, torch::kFloat32);
        
        torch::Tensor actual_output = rebaseline(input, -1, threshold);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Expected: {}", to_string(expected_output));
        spdlog::debug("Actual: {}", to_string(actual_output));

        CHECK(tensors_are_close(actual_output, expected_output));
    }
    
    TEST_CASE("1D_Multiple_Fragments") {
        spdlog::debug("--- Testing 1D Multiple Fragments ---");

        // Input: [10, 5, 0, 10, 5] (Threshold 2.0)
        // Frag 1: [10, 5] (0 to 1). Baseline (0, 10) to (1, 5). y = -5x + 10.
        //   x=0: 10, x=1: 5. Subtracted: [0, 0]
        // Frag 2: [10, 5] (3 to 4). Baseline (3, 10) to (4, 5). y = -5x + 25.
        //   x=3: 10, x=4: 5. Subtracted: [0, 0]
        
        torch::Tensor input = torch::tensor({10.0f, 5.0f, 0.0f, 10.0f, 5.0f}, torch::kFloat32);
        float threshold = 2.0f;
        
        torch::Tensor expected_output = torch::tensor({0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, torch::kFloat32);
        
        torch::Tensor actual_output = rebaseline(input, -1, threshold);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Expected: {}", to_string(expected_output));
        spdlog::debug("Actual: {}", to_string(actual_output));

        CHECK(tensors_are_close(actual_output, expected_output));
    }

    TEST_CASE("2D_Tensor_Dim_0") {
        spdlog::debug("--- Testing 2D Tensor on Dim 0 ---");

        // Input (3x3):
        // Col 0: [10, 0, 0] (Frag: 0 to 0). Length 0, should do nothing.
        // Col 1: [10, 5, 0] (Frag: 0 to 1). Baseline (0,10) to (1,5). Subtracted: [0, 0]
        // Col 2: [10, 10, 10] (Frag: 0 to 2). Baseline (0,10) to (2,10). y=10. Subtracted: [0, 0, 0]
        torch::Tensor input = torch::tensor({
            {10.0f, 10.0f, 10.0f},
            {0.0f, 5.0f, 10.0f},
            {0.0f, 0.0f, 10.0f}
        }, torch::kFloat32);
        float threshold = 2.0f; 
        int64_t dim = 0; // Operate along the column (Dim 0)

        torch::Tensor expected_output = torch::tensor({
            {10.0f, 0.0f, 0.0f},
            {0.0f, 0.0f, 0.0f},
            {0.0f, 0.0f, 0.0f}
        }, torch::kFloat32);
        
        torch::Tensor actual_output = rebaseline(input, dim, threshold);

        spdlog::debug("Input:\n{}", to_string(input));
        spdlog::debug("Expected:\n{}", to_string(expected_output));
        spdlog::debug("Actual:\n{}", to_string(actual_output));

        CHECK(tensors_are_close(actual_output, expected_output));
    }
    
    TEST_CASE("2D_Tensor_Dim_1") {
        spdlog::debug("--- Testing 2D Tensor on Dim 1 ---");

        // Input (2x4):
        // Row 0: [1, 5, 10, 5] (Frag: 1 to 3). Baseline (1, 5) to (3, 5). y=5. Subtracted: [0, 5, 0]
        // Row 1: [0, 10, 10, 0] (Frag: 1 to 2). Baseline (1, 10) to (2, 10). y=10. Subtracted: [0, 0]
        torch::Tensor input = torch::tensor({
            {1.0f, 5.0f, 10.0f, 5.0f},
            {0.0f, 10.0f, 10.0f, 0.0f}
        }, torch::kFloat32);
        float threshold = 2.0f; 
        int64_t dim = 1; // Operate along the row (Dim 1)

        torch::Tensor expected_output = torch::tensor({
            {1.0f, 0.0f, 5.0f, 0.0f}, // Row 0: 5-5=0, 10-5=5, 5-5=0
            {0.0f, 0.0f, 0.0f, 0.0f}  // Row 1: 10-10=0, 10-10=0
        }, torch::kFloat32);
        
        torch::Tensor actual_output = rebaseline(input, dim, threshold);

        spdlog::debug("Input:\n{}", to_string(input));
        spdlog::debug("Expected:\n{}", to_string(expected_output));
        spdlog::debug("Actual:\n{}", to_string(actual_output));

        CHECK(tensors_are_close(actual_output, expected_output));
    }

    TEST_CASE("Edge_Case_No_Fragments") {
        spdlog::debug("--- Testing Edge Case: No Fragments ---");

        // All values are below the threshold
        torch::Tensor input = torch::tensor({0.5f, 1.0f, -2.0f, 0.0f}, torch::kFloat32);
        float threshold = 2.0f; 

        // Output should be identical to input
        torch::Tensor expected_output = input.clone();
        torch::Tensor actual_output = rebaseline(input, -1, threshold);

        CHECK(tensors_are_close(actual_output, expected_output));
    }
}
