#include "WireCellPytorch/Torch.h"  // One-stop header.

int main(int argc, const char* argv[])
{
    std::cout << "WireCell::pytorch : test loading TorchScript Model\n";

    const std::string mname = "model.ts";

    torch::jit::script::Module module;
    try {
        // Deserialize the ScriptModule from a file using torch::jit::load().
        module = torch::jit::load(mname);
        std::cout << "model " << mname << " loaded\n";
    }
    catch (const torch::Error& e) {
        std::cout << "error loading the model\n";
        std::cout << "example models: https://github.com/HaiwangYu/Pytorch-UNet/tree/master/ts-model\n";
    }

    return 0;
}