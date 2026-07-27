#include "WireCellSpng/Testing.h"
#include "WireCellSpng/Util.h"


#include <cmath>
#include <complex>
#include <vector>
#include <iostream>

namespace {

    // Helper function to check if two tensors are close, printing debug info if not.
    inline
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

using namespace WireCell::SPNG;

TEST_SUITE("spng util various") {

    TEST_CASE("spng util vshape")
    {
        std::vector<int64_t> store = {2,3,4};
        torch::IntArrayRef tshape(store);
        auto vec = vshape(tshape);
        REQUIRE(vec.size() == 3);
        CHECK(2 == vec[0]);
        CHECK(3 == vec[1]);
        CHECK(4 == vec[2]);
    }

    TEST_CASE("nhalf") {
        // Odd size N=5 (DC, P1, P2, N2, N1) -> (5-1)/2 = 2
        CHECK(nhalf(5) == 2); 
        
        // Even size N=4 (DC, P1, NQ, N1) -> (4-2)/2 = 1
        CHECK(nhalf(4) == 1); 

        // Size N=2 (DC, N1) -> (2-2)/2 = 0
        CHECK(nhalf(2) == 0);
    }
}


TEST_SUITE("spng util resize older") {

    // --- Torch Tensor Implementation Tests (Real Float) ---

    TEST_CASE("1D Tensor resize upsample") {
        torch::Tensor in = torch::tensor({1.0f, 2.0f, 3.0f}, torch::kFloat);
        torch::Tensor expected = torch::tensor({1.0f, 2.0f, 3.0f, 0.0f, 0.0f}, torch::kFloat);
        torch::Tensor actual = resize_tensor_tail(in, 0, 5);
        CHECK(check_tensor_close(actual, expected));
    }

    TEST_CASE("2D Tensor resize upsample rows (axis 0)") {
        torch::Tensor in = torch::tensor({{1.0f, 2.0f}, {3.0f, 4.0f}}, torch::kFloat);
        torch::Tensor expected = torch::tensor({{1.0f, 2.0f}, {3.0f, 4.0f}, {0.0f, 0.0f}}, torch::kFloat);
        torch::Tensor actual = resize_tensor_tail(in, 0, 3);
        CHECK(check_tensor_close(actual, expected));
    }

    TEST_CASE("2D Tensor resize downsample columns (axis 1)") {
        torch::Tensor in = torch::tensor({{1.0f, 2.0f, 3.0f}, {4.0f, 5.0f, 6.0f}}, torch::kFloat);
        torch::Tensor expected = torch::tensor({{1.0f, 2.0f}, {4.0f, 5.0f}}, torch::kFloat);
        torch::Tensor actual = resize_tensor_tail(in, 1, 2);
        CHECK(check_tensor_close(actual, expected));
    }
}

TEST_SUITE("spng util resize tensor") {

    // Helper to check if two tensors are close (ignoring small float errors)
    auto check_equal = [](const torch::Tensor& a, const torch::Tensor& b, const std::string& message) {
        REQUIRE_MESSAGE(a.sizes().vec() == b.sizes().vec(), (message + ": Sizes mismatch."));
        REQUIRE_MESSAGE(a.dtype() == b.dtype(), (message + ": Dtype mismatch."));
        REQUIRE_MESSAGE(torch::allclose(a, b), (message + ": Tensor values mismatch."));
    };

    // Base 1D tensor: [0, 1, 2, 3, 4] (Ns=5)
    const torch::Tensor t5 = torch::arange(5, torch::kFloat32); // [0.0, 1.0, 2.0, 3.0, 4.0]

    // --- 1. General Padding Tests (Ns=5 -> Nr=8) ---
    TEST_CASE("General Padding (Ns < Nr)") {
        int64_t axis = 0;
        int64_t size = 8;
        
        // Case 1.1: Head Padding (index = 0)
        torch::Tensor r_head = resize_tensor(t5, axis, 0, size);
        torch::Tensor e_head = torch::tensor({0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0});
        check_equal(r_head, e_head, "Head Padding (index 0)");
        
        // Case 1.2: Tail Padding (index = Ns=5)
        torch::Tensor r_tail = resize_tensor(t5, axis, 5, size);
        torch::Tensor e_tail = torch::tensor({0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0});
        check_equal(r_tail, e_tail, "Tail Padding (index Ns)");
        
        // Case 1.3: Middle Padding (index = 2)
        // Input: [0, 1 | 2, 3, 4]. Padding (3 zeros) inserted at index 2.
        // Output: [0, 1, 0, 0, 0, 2, 3, 4]
        torch::Tensor r_mid = resize_tensor(t5, axis, 2, size);
        torch::Tensor e_mid = torch::tensor({0.0, 1.0, 0.0, 0.0, 0.0, 2.0, 3.0, 4.0});
        check_equal(r_mid, e_mid, "Middle Padding (index 2)");
    }

    // --- 2. General Truncation Tests (Ns=5 -> Nr=3) ---
    TEST_CASE("General Truncation (Ns > Nr)") {
        int64_t axis = 0;
        int64_t size = 3; // nloss = 2
        
        // Case 2.1: Head Truncation (index = 0)
        // Remove [0, 2). Keep [2, 5). Expected: [2, 3, 4]
        torch::Tensor r_head = resize_tensor(t5, axis, 0, size);
        torch::Tensor e_head = torch::tensor({2.0, 3.0, 4.0});
        check_equal(r_head, e_head, "Head Truncation (index 0)");

        // Case 2.2: Middle Truncation (index = 2, non-wrapping)
        // Remove [2, 4). Keep [0, 2) and [4, 5). Expected: [0, 1, 4]
        // beg=2, end=mod(2+2, 5)=4. beg < end.
        torch::Tensor r_mid = resize_tensor(t5, axis, 2, size);
        torch::Tensor e_mid = torch::tensor({0.0, 1.0, 4.0});
        check_equal(r_mid, e_mid, "Middle Truncation (index 2)");
        
        // Case 2.3: Wrap-around Truncation (Tail Truncation equivalent)
        // index = 3. nloss = 2. Remove 2 elements starting at 3: [3, 4] AND [0]. No, wait.
        // beg=3, end=mod(3+2, 5)=0. end < beg. Keep [end, beg) = [0, 3). Expected: [0, 1, 2]
        torch::Tensor r_wrap = resize_tensor(t5, axis, 3, size);
        torch::Tensor e_wrap = torch::tensor({0.0, 1.0, 2.0});
        check_equal(r_wrap, e_wrap, "Wrap Truncation (index 3)");

        // Case 2.4: Extreme Wrap-around Truncation
        // Ns=5 -> Nr=2. nloss=3. index=2.
        // beg=2, end=mod(2+3, 5)=0. end < beg. Keep [end, beg) = [0, 2). Expected: [0, 1]
        torch::Tensor r_extreme_wrap = resize_tensor(t5, axis, 2, 2);
        torch::Tensor e_extreme_wrap = torch::tensor({0.0, 1.0});
        check_equal(r_extreme_wrap, e_extreme_wrap, "Extreme Wrap Truncation (index 2, Nr=2)");
    }

    // --- 3. Special Case Helpers ---
    TEST_CASE("Helper Functions") {
        int64_t axis = 0;
        
        // Case 3.1: Head Padding (resize_tensor_head)
        torch::Tensor r_head_pad = resize_tensor_head(t5, axis, 8);
        torch::Tensor e_head_pad = torch::tensor({0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0});
        check_equal(r_head_pad, e_head_pad, "Head Padding via helper");

        // Case 3.2: Head Truncation (resize_tensor_head)
        torch::Tensor r_head_trunc = resize_tensor_head(t5, axis, 3);
        torch::Tensor e_head_trunc = torch::tensor({2.0, 3.0, 4.0});
        check_equal(r_head_trunc, e_head_trunc, "Head Truncation via helper");

        // Case 3.3: Tail Padding (resize_tensor_tail)
        torch::Tensor r_tail_pad = resize_tensor_tail(t5, axis, 8);
        torch::Tensor e_tail_pad = torch::tensor({0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0});
        check_equal(r_tail_pad, e_tail_pad, "Tail Padding via helper");

        // Case 3.4: Tail Truncation (resize_tensor_tail)
        // Ns=5 -> Nr=3. Should keep [0, 1, 2].
        torch::Tensor r_tail_trunc = resize_tensor_tail(t5, axis, 3);
        torch::Tensor e_tail_trunc = torch::tensor({0.0, 1.0, 2.0});
        check_equal(r_tail_trunc, e_tail_trunc, "Tail Truncation via helper");
    }
    
    // --- 4. Edge Cases and Multi-Dimensional Tests ---
    TEST_CASE("Multi-Dimensional and Edge Cases") {
        // 2D tensor (2, 5): 
        // [[0, 1, 2, 3, 4],
        //  [5, 6, 7, 8, 9]]
        torch::Tensor t2x5 = torch::arange(10, torch::kFloat32).reshape({2, 5});

        // Case 4.1: No change
        torch::Tensor r_no_change = resize_tensor(t5, 0, 0, 5);
        check_equal(r_no_change, t5, "No change (Ns=Nr)");
        
        // Case 4.2: Negative Index
        // resize_tensor(t5, 0, -2, 3). Index -2 maps to 3. (Wrap Truncate).
        // Expected: [0, 1, 2]
        torch::Tensor r_neg_idx = resize_tensor(t5, 0, -2, 3);
        torch::Tensor e_neg_idx = torch::tensor({0.0, 1.0, 2.0});
        check_equal(r_neg_idx, e_neg_idx, "Negative Index (index -2)");

        // Case 4.3: Resizing Axis 1 of 2D Tensor (Ns=5 -> Nr=3, Tail Truncation)
        // Input: [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]]
        // Expected: [[0, 1, 2], [5, 6, 7]]
        torch::Tensor r_2d_tail_trunc = resize_tensor_tail(t2x5, 1, 3);
        torch::Tensor e_2d_tail_trunc = torch::tensor({{0.0, 1.0, 2.0}, {5.0, 6.0, 7.0}});
        check_equal(r_2d_tail_trunc, e_2d_tail_trunc, "2D: Axis 1 Tail Truncation");
        
        // Case 4.4: Resizing Axis 0 of 2D Tensor (Ns=2 -> Nr=4, Head Padding)
        // Input: [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]]
        // Expected: [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 1, 2, 3, 4], [5, 6, 7, 8, 9]]
        torch::Tensor r_2d_head_pad = resize_tensor_head(t2x5, 0, 4);
        torch::Tensor e_2d_head_pad = torch::zeros({4, 5}, torch::kFloat32);
        e_2d_head_pad.narrow(0, 2, 2).copy_(t2x5);
        check_equal(r_2d_head_pad, e_2d_head_pad, "2D: Axis 0 Head Padding");

        // Case 4.5: Negative Axis Index
        // resize_tensor(t2x5, -1, 0, 3) -> Axis -1 is axis 1. Head Truncation.
        // Expected: [[2, 3, 4], [7, 8, 9]]
        torch::Tensor r_neg_axis = resize_tensor(t2x5, -1, 0, 3);
        torch::Tensor e_neg_axis = torch::tensor({{2.0, 3.0, 4.0}, {7.0, 8.0, 9.0}});
        check_equal(r_neg_axis, e_neg_axis, "Negative Axis Index");
    }
}

TEST_SUITE("spng util resize middle") {
    
    // --- Helper Definitions ---
    // Define unique values for easier visual tracking of truncation/padding
    const torch::Tensor IN_ODD = torch::tensor({1.0f, 2.0f, 3.0f, 4.0f, 5.0f}, torch::kFloat); // Ns=5
    // Halves based on Ns=5: [1.0, 2.0, 3.0] + [4.0, 5.0] -> P_size=3, L_size=2
    
    const torch::Tensor IN_EVEN = torch::tensor({10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f}, torch::kFloat); // Ns=6
    // Halves based on Ns=6: [10, 11, 12] + [13, 14, 15] -> P_size=3, L_size=3
    
    const float ZERO = 0.0f;

    // --- Upsampling Tests ---

    TEST_CASE("Upsample Odd (Ns=5 to Nr=7, Odd)") {
        // N_min = 5. P_size=3, L_size=2.
        // Copy [1, 2, 3] | Insert 2 zeros | Copy [4, 5]
        torch::Tensor expected = torch::tensor({1.0f, 2.0f, 3.0f, ZERO, ZERO, 4.0f, 5.0f}, torch::kFloat);
        torch::Tensor actual = resize_tensor_middle(IN_ODD, 0, 7);
        check_tensor_close(actual, expected);
    }
    
    TEST_CASE("Upsample Even (Ns=6 to Nr=8, Even)") {
        // N_min = 6. P_size=3, L_size=3.
        // Copy [10, 11, 12] | Insert 2 zeros | Copy [13, 14, 15]
        torch::Tensor expected = torch::tensor({10.0f, 11.0f, 12.0f, ZERO, ZERO, 13.0f, 14.0f, 15.0f}, torch::kFloat);
        torch::Tensor actual = resize_tensor_middle(IN_EVEN, 0, 8);
        check_tensor_close(actual, expected);
    }
    
    TEST_CASE("Upsample Odd to Even (Ns=5 to Nr=6)") {
        // N_min = 5. P_size=3, L_size=2.
        // Copy [1, 2, 3] | Insert 1 zero | Copy [4, 5]
        torch::Tensor expected = torch::tensor({1.0f, 2.0f, 3.0f, ZERO, 4.0f, 5.0f}, torch::kFloat);
        torch::Tensor actual = resize_tensor_middle(IN_ODD, 0, 6);
        check_tensor_close(actual, expected);
    }

    // --- Downsampling Tests ---
    
    TEST_CASE("Downsample Odd (Ns=5 to Nr=3, Odd)") {
        // N_min = 3. 
        // P_size = (3+1)/2 = 2. Copy [1, 2].
        // L_size = (3-1)/2 = 1. Copy [5].
        // Result: [1, 2, 5]. Samples 3, 4 truncated from the middle.
        torch::Tensor expected = torch::tensor({1.0f, 2.0f, 5.0f}, torch::kFloat);
        torch::Tensor actual = resize_tensor_middle(IN_ODD, 0, 3);
        check_tensor_close(actual, expected);
    }

    TEST_CASE("Downsample Even (Ns=6 to Nr=4, Even)") {
        // N_min = 4. 
        // P_size = 4/2 = 2. Copy [10, 11].
        // L_size = 4/2 = 2. Copy [14, 15].
        // Result: [10, 11, 14, 15]. Samples 12, 13 truncated from the middle.
        torch::Tensor expected = torch::tensor({10.0f, 11.0f, 14.0f, 15.0f}, torch::kFloat);
        torch::Tensor actual = resize_tensor_middle(IN_EVEN, 0, 4);
        check_tensor_close(actual, expected);
    }

    TEST_CASE("Downsample Even to Odd (Ns=6 to Nr=5)") {
        // N_min = 5.
        // P_size = 3. Copy [10, 11, 12].
        // L_size = 2. Copy [14, 15].
        // Result: [10, 11, 12, 14, 15]. Sample 13 truncated from the middle.
        torch::Tensor expected = torch::tensor({10.0f, 11.0f, 12.0f, 14.0f, 15.0f}, torch::kFloat);
        torch::Tensor actual = resize_tensor_middle(IN_EVEN, 0, 5);
        check_tensor_close(actual, expected);
    }
    
}
