#include "WireCellPytorch/Torch.h"
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
    // When run from waf we run in the build directory
    std::string ts = "../pytorch/test/extsmod.ts";
    // interactive users must give a model to load
    if (argc > 1) {
        ts = argv[1];
    }

    torch::Tensor tensor = torch::rand({2, 3});
    std::cout << tensor << std::endl;

    torch::jit::script::Module mod;

    try {
        // assume we run from top wtc source area
        mod = torch::jit::load(ts.c_str());
    }
    catch (const c10::Error& e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }
    try {
        mod.to(at::kCUDA);
    }
    catch (const c10::Error& e) {
        std::cerr << "Warning: not CUDA available:\n" << e.what() << std::endl;
    }
    try {
        mod.to(at::kCPU);
    }
    catch (const c10::Error& e) {
        std::cerr << "Error: failed to move model to CPU:\n" << e.what() << std::endl;
        return -1;
    }

    return 0;
}
