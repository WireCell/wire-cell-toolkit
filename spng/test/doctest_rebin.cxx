#include "WireCellSpng/Rebin.h"
#include "WireCellSpng/Testing.h"

#include <iostream>

using namespace WireCell::SPNG;



TEST_CASE("spng torch rebin") {

    // Create a 1D tensor
    torch::Tensor input_tensor = torch::arange(0, 10, torch::kFloat32);
    std::cout << "Original Tensor:\n" << input_tensor << std::endl;
    int64_t dim = 0;
    int64_t factor = 2; // Test factor

    // 1. Downsampling Tests (factor = 2)
    std::cout << "\n--- Downsampling (Factor " << factor << ") ---" << std::endl;
    
    // Integral (Sum: [0+1, 2+3, 4+5, 6+7, 8+9] = [1, 5, 9, 13, 17])
    torch::Tensor down_integral = Rebin::downsample(input_tensor, factor, dim,
                                                    Rebin::Normalization::kIntegral);
    std::cout << "Integral (Sum):\n" << down_integral << std::endl;

    // Interpolation (Average: [0.5, 2.5, 4.5, 6.5, 8.5])
    torch::Tensor down_interp = Rebin::downsample(input_tensor, factor, dim,
                                                  Rebin::Normalization::kInterpolation);
    std::cout << "Interpolation (Average):\n" << down_interp << std::endl;
    
    // Maximum (Max: [1, 3, 5, 7, 9])
    torch::Tensor down_max = Rebin::downsample(input_tensor, factor, dim,
                                               Rebin::Normalization::kMaximum);
    std::cout << "Maximum (Max):\n" << down_max << std::endl;


    // 2. Upsampling Tests (factor = 2)
    std::cout << "\n--- Upsampling (Factor " << factor << ") ---" << std::endl;
    
    // Input for upsampling is a smaller tensor: [10, 20, 30]
    torch::Tensor up_input = torch::tensor({10.0f, 20.0f, 30.0f});
    std::cout << "Upsample Input:\n" << up_input << std::endl;

    // Integral (Value / factor, repeated: [5, 5, 10, 10, 15, 15])
    torch::Tensor up_integral = Rebin::upsample(up_input, factor, dim,
                                                Rebin::Normalization::kIntegral);
    std::cout << "Integral (Value/Factor):\n" << up_integral << std::endl;

    // Interpolation (Value repeated: [10, 10, 20, 20, 30, 30])
    torch::Tensor up_interp = Rebin::upsample(up_input, factor, dim,
                                              Rebin::Normalization::kInterpolation);
    std::cout << "Interpolation (Value Repeated):\n" << up_interp << std::endl;
    
    // Maximum (Same as Interpolation: [10, 10, 20, 20, 30, 30])
    torch::Tensor up_max = Rebin::upsample(up_input, factor, dim,
                                           Rebin::Normalization::kMaximum);
    std::cout << "Maximum (Value Repeated):\n" << up_max << std::endl;
}
