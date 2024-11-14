#include <torch/script.h>  // One-stop header.
#include <chrono>

/**
 * need to configure with torch_cpu on gpvm
 * e.g.--with-libtorch="$LIBTORCH_FQ_DIR/" --with-libtorch-libs torch,torch_cpu,c10
 * model for testing can be founc here:
 * https://www.phy.bnl.gov/~hyu/dunefd/dnn-roi-pdvd/Pytorch-UNet/ts-model/
*/

int main(int argc, const char* argv[])
{
    std::cout << "WireCell::pytorch : test loading TorchScript Model\n";

    const std::string mname = "model.ts";
    auto dtype = torch::kFloat16;

    torch::jit::script::Module module;
    // Deserialize the ScriptModule from a file using torch::jit::load().
    auto start = std::chrono::high_resolution_clock::now();
    module = torch::jit::load(mname);
    module.to(at::kCPU, dtype);
    torch::TensorOptions options = torch::TensorOptions().dtype(dtype);
    torch::Tensor iten = torch::rand({1, 3, 800, 600}, options);
    // torch::Tensor iten = torch::zeros({1, 3, 800, 600}, options);
    std::vector<torch::IValue> itens {iten};
    auto otens = module.forward(itens).toTensor();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "timing: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

    return 0;
}