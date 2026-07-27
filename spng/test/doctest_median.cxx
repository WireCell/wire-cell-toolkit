#include "WireCellSpng/Testing.h"
#include "WireCellSpng/Torch.h"

// Function to demonstrate the operation
static
torch::Tensor center_rows_by_median(const torch::Tensor& input_tensor) {
    // 1. Calculate the median across the columns (dimension 2)
    // The input shape is (nbatch, nrows, ncols).
    // We want the median for each row of each batch, so we calculate
    // the median over dimension 2 (the columns).
    // The result 'medians_per_row' will have shape (nbatch, nrows).
    torch::Tensor medians_per_row = std::get<0>(torch::median(input_tensor, /*dim=*/2, /*keepdim=*/false));
    REQUIRE(medians_per_row.size(0) == input_tensor.size(0));
    REQUIRE(medians_per_row.size(1) == input_tensor.size(1));


    // 2. Adjust the shape of the medians for broadcasting
    // To subtract the (nbatch, nrows) tensor from the (nbatch, nrows, ncols)
    // input tensor, we need to insert a new dimension of size 1 at position 2.
    // The adjusted shape will be (nbatch, nrows, 1).
    torch::Tensor medians_broadcastable = torch::unsqueeze(medians_per_row, /*dim=*/2);

    // 3. Subtract the broadcastable medians from the original tensor
    // LibTorch will automatically broadcast the (nbatch, nrows, 1) tensor
    // across the column dimension (dim 2) of the (nbatch, nrows, ncols) tensor.
    torch::Tensor result_tensor = input_tensor - medians_broadcastable;

    // Output the shapes for verification
    std::cout << "Input Tensor shape: " << input_tensor.sizes() << std::endl;
    std::cout << "Medians Per Row shape (pre-unsqueeze): " << medians_per_row.sizes() << std::endl;
    std::cout << "Medians Broadcastable shape: " << medians_broadcastable.sizes() << std::endl;
    std::cout << "Result Tensor shape: " << result_tensor.sizes() << std::endl;

    // Optional: Print a small part of the results for visual check
    if (input_tensor.numel() > 0) {
        std::cout << "\nOriginal Input (First Batch, First Row):\n" << input_tensor.index({0, 0, "..."}) << std::endl;
        std::cout << "Median to subtract (First Batch, First Row):\n" << medians_broadcastable.index({0, 0, "..."}) << std::endl;
        std::cout << "Result (First Batch, First Row):\n" << result_tensor.index({0, 0, "..."}) << std::endl;
    }
    return result_tensor;
}

TEST_SUITE("spng median") {

    TEST_CASE("spng median random") {
        
        torch::Tensor x = torch::rand({2, 3, 4});
        
        std::cout << "--- Running Row Median Centering ---" << std::endl;
        auto res = center_rows_by_median(x);
        std::cout << "------------------------------------" << std::endl;
        CHECK(res.size(0) == 2);
        CHECK(res.size(1) == 3);
        CHECK(res.size(2) == 4);

    }

    TEST_CASE("spng median specific") {

        // Another example with specific values
        torch::Tensor y = torch::tensor({
                {{1.0, 2.0, 3.0, 10.0},
                 {5.0, 6.0, 7.0, 8.0},
                 {0.0, 0.0, 0.0, 0.0}},
                {{10.0, 10.0, 10.0, 10.0},
                 {1.0, 1.0, 2.0, 1.0},
                 {10.0, 20.0, 30.0, 40.0}}
            });

        std::cout << "\n--- Running Specific Value Test ---" << std::endl;
        center_rows_by_median(y);
        std::cout << "------------------------------------" << std::endl;
    }
}

