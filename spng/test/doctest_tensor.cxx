#include "WireCellSpng/Testing.h"
#include "WireCellSpng/Torch.h"

//using namespace WireCell::SPNG;

TEST_CASE("spng torch tensor metadata") {
    auto ten = torch::zeros({2,3,4}, torch::kFloat);
    std::cerr << "tensor dtype streamed : " << ten.dtype() << "\n";
    std::cerr << "tensor dtype c10::toSTring " << c10::toString(ten.dtype()) << "\n";
}
TEST_CASE("spng torch tensor metadata gemini mr smartypants") {
    torch::Tensor tensor_float = torch::ones({2, 2}, torch::kFloat32);
    std::string dtype_str_float = c10::toString(tensor_float.dtype());
    std::cout << "Float Tensor dtype: " << dtype_str_float << std::endl; // Output: "Float"

    torch::Tensor tensor_long = torch::zeros({3, 3}, torch::kInt64);
    std::string dtype_str_long = c10::toString(tensor_long.dtype());
    std::cout << "Long Tensor dtype: " << dtype_str_long << std::endl; // Output: "Long"

    // If you specifically want lowercase "float", "long", etc.:
    std::transform(dtype_str_float.begin(), dtype_str_float.end(), dtype_str_float.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    std::cout << "Lowercase Float Tensor dtype: " << dtype_str_float << std::endl; // Output: "float"


}
