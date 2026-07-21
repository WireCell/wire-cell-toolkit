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
        
        torch::Tensor actual_output = rebaseline(input, -1, threshold, 1, 0, false, false);

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
        
        torch::Tensor actual_output = rebaseline(input, -1, threshold, 1, 0, false, false);

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
        
        torch::Tensor actual_output = rebaseline(input, -1, threshold, 1, 0, false, false);

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
        
        torch::Tensor actual_output = rebaseline(input, dim, threshold, 1, 0, false, false);

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
        
        torch::Tensor actual_output = rebaseline(input, dim, threshold, 1, 0, false, false);

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
        torch::Tensor actual_output = rebaseline(input, -1, threshold, 1, 0, false, false);

        CHECK(tensors_are_close(actual_output, expected_output));
    }

    // The following four tests share input [12, 3, 12, 0, 5, 0] with threshold=2
    // and min_roi_size=1.  Two fragments are present:
    //   Fragment 1: indices 0-2 (size 3, not small)
    //   Fragment 2: index  4   (size 1, small under min_roi_size=1)
    // Each test changes one (or all three) of the nominal values
    // expand_size=0, remove_small=false, remove_negative=false.

    TEST_CASE("1D_ExpandSize_One") {
        spdlog::debug("--- Testing 1D with expand_size=1 ---");

        // Fragment 1 expanded to [0,3]: start_val=12, end_val=0, slope=-4.
        //   baseline=[12,8,4,0], wave[0:4]=[12,3,12,0]-[12,8,4,0]=[0,-5,8,0]
        // Fragment 2: size 1 (small), remove_small=false -> skipped, stays 5.
        // Expected: [0, -5, 8, 0, 5, 0]

        torch::Tensor input = torch::tensor({12.0f, 3.0f, 12.0f, 0.0f, 5.0f, 0.0f}, torch::kFloat32);
        float threshold = 2.0f;

        torch::Tensor expected = torch::tensor({0.0f, -5.0f, 8.0f, 0.0f, 5.0f, 0.0f}, torch::kFloat32);
        torch::Tensor actual = rebaseline(input, -1, threshold, 1, 1, false, false);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Expected: {}", to_string(expected));
        spdlog::debug("Actual: {}", to_string(actual));

        CHECK(tensors_are_close(actual, expected));
    }

    TEST_CASE("1D_RemoveSmall_True") {
        spdlog::debug("--- Testing 1D with remove_small=true ---");

        // Fragment 1: no expand, start_val=12, end_val=12, slope=0.
        //   baseline=[12,12,12], wave[0:3]=[12,3,12]-[12,12,12]=[0,-9,0]
        // Fragment 2: size 1 (small), remove_small=true -> zeroed.
        // Expected: [0, -9, 0, 0, 0, 0]

        torch::Tensor input = torch::tensor({12.0f, 3.0f, 12.0f, 0.0f, 5.0f, 0.0f}, torch::kFloat32);
        float threshold = 2.0f;

        torch::Tensor expected = torch::tensor({0.0f, -9.0f, 0.0f, 0.0f, 0.0f, 0.0f}, torch::kFloat32);
        torch::Tensor actual = rebaseline(input, -1, threshold, 1, 0, true, false);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Expected: {}", to_string(expected));
        spdlog::debug("Actual: {}", to_string(actual));

        CHECK(tensors_are_close(actual, expected));
    }

    TEST_CASE("1D_RemoveNegative_True") {
        spdlog::debug("--- Testing 1D with remove_negative=true ---");

        // Fragment 1: no expand, baseline=[12,12,12], yields [0,-9,0].
        // Fragment 2: size 1 (small), remove_small=false -> skipped, stays 5.
        // Before clamp: [0,-9,0,0,5,0].  After clamp_min(0): [0,0,0,0,5,0].
        // Expected: [0, 0, 0, 0, 5, 0]

        torch::Tensor input = torch::tensor({12.0f, 3.0f, 12.0f, 0.0f, 5.0f, 0.0f}, torch::kFloat32);
        float threshold = 2.0f;

        torch::Tensor expected = torch::tensor({0.0f, 0.0f, 0.0f, 0.0f, 5.0f, 0.0f}, torch::kFloat32);
        torch::Tensor actual = rebaseline(input, -1, threshold, 1, 0, false, true);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Expected: {}", to_string(expected));
        spdlog::debug("Actual: {}", to_string(actual));

        CHECK(tensors_are_close(actual, expected));
    }

    TEST_CASE("1D_All_Three_Options") {
        spdlog::debug("--- Testing 1D with expand_size=1, remove_small=true, remove_negative=true ---");

        // Fragment 1 expanded to [0,3]: start_val=12, end_val=0, slope=-4.
        //   baseline=[12,8,4,0], wave[0:4]=[12,3,12,0]-[12,8,4,0]=[0,-5,8,0]
        // Fragment 2: size 1 (small), remove_small=true -> zeroed.
        // Before clamp: [0,-5,8,0,0,0].  After clamp_min(0): [0,0,8,0,0,0].
        // Expected: [0, 0, 8, 0, 0, 0]

        torch::Tensor input = torch::tensor({12.0f, 3.0f, 12.0f, 0.0f, 5.0f, 0.0f}, torch::kFloat32);
        float threshold = 2.0f;

        torch::Tensor expected = torch::tensor({0.0f, 0.0f, 8.0f, 0.0f, 0.0f, 0.0f}, torch::kFloat32);
        torch::Tensor actual = rebaseline(input, -1, threshold, 1, 1, true, true);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Expected: {}", to_string(expected));
        spdlog::debug("Actual: {}", to_string(actual));

        CHECK(tensors_are_close(actual, expected));
    }
}

// --- Test Suite for rebaseline_zero ---
TEST_SUITE("spng rebaseline_zero") {

    TEST_CASE("1D_Simple_ROI_bounded_by_zeros") {
        // [0,0,5,10,5,0,0]: positions 0-1 and 5-6 are 2-zero separators.
        // ROI: indices 2-4 = [5,10,5].
        // Baseline: wave[2]=5, wave[4]=5. Slope=0. Baseline=[5,5,5].
        // After: [0,5,0]. Output: [0,0,0,5,0,0,0].
        torch::Tensor input = torch::tensor({0.0f, 0.0f, 5.0f, 10.0f, 5.0f, 0.0f, 0.0f}, torch::kFloat32);
        torch::Tensor expected = torch::tensor({0.0f, 0.0f, 0.0f, 5.0f, 0.0f, 0.0f, 0.0f}, torch::kFloat32);
        torch::Tensor actual = rebaseline_zero(input);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Expected: {}", to_string(expected));
        spdlog::debug("Actual: {}", to_string(actual));

        CHECK(tensors_are_close(actual, expected));
    }

    TEST_CASE("1D_ROI_at_tensor_edge") {
        // [5,10,5,0,0]: tensor left-edge acts as boundary; positions 3-4 are separator.
        // ROI: indices 0-2 = [5,10,5].
        // Baseline: wave[0]=5, wave[2]=5. Slope=0. Baseline=[5,5,5].
        // After: [0,5,0]. Output: [0,5,0,0,0].
        torch::Tensor input = torch::tensor({5.0f, 10.0f, 5.0f, 0.0f, 0.0f}, torch::kFloat32);
        torch::Tensor expected = torch::tensor({0.0f, 5.0f, 0.0f, 0.0f, 0.0f}, torch::kFloat32);
        torch::Tensor actual = rebaseline_zero(input);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Expected: {}", to_string(expected));
        spdlog::debug("Actual: {}", to_string(actual));

        CHECK(tensors_are_close(actual, expected));
    }

    TEST_CASE("1D_ROI_with_zero_crossing") {
        // [5,0,8,0,0]: position 1 is a single zero (<2), so NOT a separator — it is a zero
        // crossing inside the ROI. Positions 3-4 are separator.
        // ROI: indices 0-2 = [5,0,8].
        // Baseline: wave[0]=5, wave[2]=8. Slope=(8-5)/2=1.5. Baseline=[5,6.5,8].
        // After subtraction: [0,-6.5,0]. Output: [0,-6.5,0,0,0].
        torch::Tensor input = torch::tensor({5.0f, 0.0f, 8.0f, 0.0f, 0.0f}, torch::kFloat32);
        torch::Tensor expected = torch::tensor({0.0f, -6.5f, 0.0f, 0.0f, 0.0f}, torch::kFloat32);
        torch::Tensor actual = rebaseline_zero(input);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Expected: {}", to_string(expected));
        spdlog::debug("Actual: {}", to_string(actual));

        CHECK(tensors_are_close(actual, expected));
    }

    TEST_CASE("1D_consecutive_zeros_distinguishes_boundary_from_crossing") {
        // Demonstrate that consecutive_zeros=2 treats a single zero as a crossing
        // but consecutive_zeros=1 treats it as a boundary.
        //
        // Input: [5,0,8], consecutive_zeros=2:
        //   Single zero at index 1 is NOT a separator -> one ROI [0,2].
        //   Baseline: wave[0]=5, wave[2]=8. Slope=1.5. After: [0,-6.5,0].
        //
        // Same input, consecutive_zeros=1:
        //   Single zero at index 1 IS a separator -> ROIs [0,0]=[5] and [2,2]=[8].
        //   Both size 1 (small, min_roi_size=1 default), so skipped. Output unchanged.
        torch::Tensor input = torch::tensor({5.0f, 0.0f, 8.0f}, torch::kFloat32);

        torch::Tensor actual_2 = rebaseline_zero(input, -1, /*consequtive_zeros=*/2);
        torch::Tensor expected_2 = torch::tensor({0.0f, -6.5f, 0.0f}, torch::kFloat32);
        CHECK(tensors_are_close(actual_2, expected_2));

        torch::Tensor actual_1 = rebaseline_zero(input, -1, /*consequtive_zeros=*/1);
        torch::Tensor expected_1 = input.clone(); // unchanged (both ROIs are size 1 = small)
        CHECK(tensors_are_close(actual_1, expected_1));
    }

    TEST_CASE("1D_multiple_ROIs") {
        // [0,0,5,10,5,0,0,8,12,8,0,0]: two ROIs separated by 2-zero runs.
        // ROI 1: indices 2-4 = [5,10,5].  wave[2]=5, wave[4]=5.  Slope=0.  After: [0,5,0].
        // ROI 2: indices 7-9 = [8,12,8].  wave[7]=8, wave[9]=8.  Slope=0.  After: [0,4,0].
        // Output: [0,0,0,5,0,0,0,0,4,0,0,0].
        torch::Tensor input = torch::tensor(
            {0.0f, 0.0f, 5.0f, 10.0f, 5.0f, 0.0f, 0.0f, 8.0f, 12.0f, 8.0f, 0.0f, 0.0f},
            torch::kFloat32);
        torch::Tensor expected = torch::tensor(
            {0.0f, 0.0f, 0.0f, 5.0f, 0.0f, 0.0f, 0.0f, 0.0f, 4.0f, 0.0f, 0.0f, 0.0f},
            torch::kFloat32);
        torch::Tensor actual = rebaseline_zero(input);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Expected: {}", to_string(expected));
        spdlog::debug("Actual: {}", to_string(actual));

        CHECK(tensors_are_close(actual, expected));
    }

    TEST_CASE("1D_whole_tensor_is_one_ROI") {
        // No zeros -> entire tensor is one ROI.
        // [5,10,15]: wave[0]=5, wave[2]=15. Slope=5. Baseline=[5,10,15]. After: [0,0,0].
        torch::Tensor input = torch::tensor({5.0f, 10.0f, 15.0f}, torch::kFloat32);
        torch::Tensor expected = torch::zeros({3}, torch::kFloat32);
        torch::Tensor actual = rebaseline_zero(input);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Expected: {}", to_string(expected));
        spdlog::debug("Actual: {}", to_string(actual));

        CHECK(tensors_are_close(actual, expected));
    }

    TEST_CASE("1D_shrink_size") {
        // [0,0,5,10,15,10,5,0,0]: ROI indices 2-6 = [5,10,15,10,5]. shrink_size=1.
        // After shrink: indices 3-5 = [10,15,10]. Size=3 (not small, min_roi_size=1 default).
        // wave[3]=10, wave[5]=10. Slope=0. Baseline=[10,10,10]. After: [0,5,0].
        // Positions 2 and 6 (the shrunk-off ends, values 5) are unchanged.
        // Output: [0,0,5,0,5,0,5,0,0].
        torch::Tensor input = torch::tensor(
            {0.0f, 0.0f, 5.0f, 10.0f, 15.0f, 10.0f, 5.0f, 0.0f, 0.0f}, torch::kFloat32);
        torch::Tensor expected = torch::tensor(
            {0.0f, 0.0f, 5.0f, 0.0f, 5.0f, 0.0f, 5.0f, 0.0f, 0.0f}, torch::kFloat32);
        torch::Tensor actual = rebaseline_zero(input,/*dim=*/-1, /*consequtive_zeros=*/2, 
                                               /*min_roi_size=*/1, /*shrink_size=*/1);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Expected: {}", to_string(expected));
        spdlog::debug("Actual: {}", to_string(actual));

        CHECK(tensors_are_close(actual, expected));
    }

    TEST_CASE("1D_shrink_collapses_small_roi") {
        // [0,0,5,10,5,0,0]: ROI [2,4], shrink_size=1 -> shrunken ROI [3,3], size=1 (small).
        // remove_small=false: skip without modifying.  Output unchanged from input.
        torch::Tensor input = torch::tensor({0.0f, 0.0f, 5.0f, 10.0f, 5.0f, 0.0f, 0.0f}, torch::kFloat32);
        torch::Tensor expected = input.clone();
        torch::Tensor actual = rebaseline_zero(input,  /*dim=*/-1, /*consequtive_zeros=*/2,
                                               /*min_roi_size=*/1, /*shrink_size=*/1);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Actual: {}", to_string(actual));

        CHECK(tensors_are_close(actual, expected));
    }

    TEST_CASE("1D_remove_small") {
        // [0,0,5,0,0]: ROI [2,2] = [5], size=1 (small, min_roi_size=1 default).
        // remove_small=false: skip, output unchanged.
        // remove_small=true:  zero the ROI, output all zeros.
        torch::Tensor input = torch::tensor({0.0f, 0.0f, 5.0f, 0.0f, 0.0f}, torch::kFloat32);

        torch::Tensor actual_keep = rebaseline_zero(input, -1, 2, 1, 0, /*remove_small=*/false);
        CHECK(tensors_are_close(actual_keep, input));

        torch::Tensor actual_zero = rebaseline_zero(input, -1, 2, 1, 0, /*remove_small=*/true);
        CHECK(tensors_are_close(actual_zero, torch::zeros_like(input)));
    }

    TEST_CASE("1D_remove_negative") {
        // [0,0,5,3,10,0,0]: ROI [2,4] = [5,3,10].
        // wave[2]=5, wave[4]=10. Slope=2.5. Baseline=[5,7.5,10].
        // After: [0,-4.5,0]. Full output before clamp: [0,0,0,-4.5,0,0,0].
        // With remove_negative=true: clamp -> [0,0,0,0,0,0,0].
        torch::Tensor input = torch::tensor({0.0f, 0.0f, 5.0f, 3.0f, 10.0f, 0.0f, 0.0f}, torch::kFloat32);
        torch::Tensor expected = torch::zeros({7}, torch::kFloat32);
        torch::Tensor actual = rebaseline_zero(input, -1, 2, 1, 0, false, /*remove_negative=*/true);

        spdlog::debug("Input: {}", to_string(input));
        spdlog::debug("Expected: {}", to_string(expected));
        spdlog::debug("Actual: {}", to_string(actual));

        CHECK(tensors_are_close(actual, expected));
    }

    TEST_CASE("2D_Tensor") {
        // Input (2x7), operating along dim=1 (rows).
        // Row 0: [0,0,5,10,5,0,0] -> ROI [2,4], baseline flat at 5. After: [0,0,0,5,0,0,0].
        // Row 1: [0,0,8,12,8,0,0] -> ROI [2,4], baseline flat at 8. After: [0,0,0,4,0,0,0].
        torch::Tensor input = torch::tensor({
            {0.0f, 0.0f, 5.0f, 10.0f, 5.0f, 0.0f, 0.0f},
            {0.0f, 0.0f, 8.0f, 12.0f, 8.0f, 0.0f, 0.0f}
        }, torch::kFloat32);
        torch::Tensor expected = torch::tensor({
            {0.0f, 0.0f, 0.0f, 5.0f, 0.0f, 0.0f, 0.0f},
            {0.0f, 0.0f, 0.0f, 4.0f, 0.0f, 0.0f, 0.0f}
        }, torch::kFloat32);
        torch::Tensor actual = rebaseline_zero(input, /*dim=*/1, 2);

        spdlog::debug("Input:\n{}", to_string(input));
        spdlog::debug("Expected:\n{}", to_string(expected));
        spdlog::debug("Actual:\n{}", to_string(actual));

        CHECK(tensors_are_close(actual, expected));
    }

    TEST_CASE("Edge_Case_all_zeros") {
        // All zeros: no non-separator positions, so no ROIs. Output equals input.
        torch::Tensor input = torch::zeros({5}, torch::kFloat32);
        torch::Tensor actual = rebaseline_zero(input);
        CHECK(tensors_are_close(actual, input));
    }
}
