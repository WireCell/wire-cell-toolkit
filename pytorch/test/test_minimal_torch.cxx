#include "WireCellPytorch/Torch.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Locate the test TorchScript model.  The unit-test runner's working
// directory has varied (build/ vs build/pytorch/), so try the candidate
// locations rather than assuming a single one.
static std::string find_extsmod()
{
    const std::vector<std::string> candidates = {
        "../pytorch/test/extsmod.ts",     // cwd = build/
        "../../pytorch/test/extsmod.ts",  // cwd = build/pytorch/
        "pytorch/test/extsmod.ts",        // cwd = source top
    };
    for (const auto& path : candidates) {
        if (std::ifstream(path).good()) {
            return path;
        }
    }
    return candidates[0];
}

int main(int argc, char* argv[])
{
    std::string ts = find_extsmod();
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
