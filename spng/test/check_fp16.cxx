#include "WireCellSpng/Torch.h"  // One-stop header.
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/TorchService.h"
#include "WireCellUtil/Exceptions.h"
#include <chrono>
#include <torch/torch.h>

/**
 * need to configure with torch_cpu on gpvm
 * e.g.--with-libtorch="$LIBTORCH_FQ_DIR/" --with-libtorch-libs torch,torch_cpu,c10
 * model for testing can be found here:
 * https://www.phy.bnl.gov/~hyu/dunefd/dnn-roi-pdvd/Pytorch-UNet/ts-model/
*/

using namespace WireCell;
using namespace WireCell::SPNG;


ITorchTensorSet::pointer to_itensor(const std::vector<torch::IValue>& inputs) {
    auto itv = std::make_shared<ITorchTensor::vector>();

    for (size_t i = 0; i < inputs.size(); ++i) {
        try {
            const auto& ivalue = inputs[i];
            if (!ivalue.isTensor()) {
                std::cerr << "Error: Expected torch::IValue at index " << i << " to be a Tensor\n";
                continue;
            }
            torch::Tensor ten = ivalue.toTensor();


            if (ten.dim() != 4) {
                std::cerr << "Error: Tensor at index " << i << " must be 4D, got " << ten.dim() << std::endl;
                continue;
            }

            if(ten.scalar_type() != torch::kFloat32) {
                ten = ten.to(torch::kFloat32);
                std::cout << "Converted tensor " << i << " to float32." << std::endl;
            }

            std::cout << "Tensor " << i << ": shape=" << ten.sizes() 
            << ", dtype=" << ten.dtype() 
            << ", device=" << ten.device() << std::endl;

            // No forced dtype or device conversion here to preserve input tensors as-is
            auto stp = std::make_shared<SimpleTorchTensor>(ten);
            itv->emplace_back(stp);

        } catch (const std::exception& e) {
            std::cerr << "Exception caught while processing tensor " << i << ": " << e.what() << std::endl;
        } catch (...) {
            std::cerr << "Unknown exception caught while processing tensor " << i << std::endl;
        }
    }
    return std::make_shared<SimpleTorchTensorSet>(0, Json::nullValue, itv);
}

int main(int argc, const char* argv[])
{
    std::cout<<"Test with the GPU information on the machine "<<std::endl;
    if (torch::cuda::is_available()) {
        std::cout << "CUDA is available!" << std::endl;
        int device_count = torch::cuda::device_count();
        std::cout << "Number of GPUs: " << device_count << std::endl;
        for (int i = 0; i < device_count; ++i) {
            torch::Device device(torch::kCUDA, i);
            std::cout << "GPU " << i << ": " << device.str() << std::endl;
        }
    } else {
        std::cout << "CUDA is not available. Using CPU." << std::endl;
    }

    std::cout << "WireCell::SPNG : test loading TorchScript Model\n";
    auto torch_service = std::make_shared<SPNG::TorchService>();
    // Configure the TorchService
    WireCell::Configuration cfg;
    const std::string model_path = argv[1];
    //const std::string model_path="/nfs/data/1/abashyal/spng/spng_dev_050525/Pytorch-UNet/ts-model-2.3/unet-l23-cosmic500-e50.ts";
    cfg["model"] = model_path;
    cfg["device"] = "gpu0"; // Use GPU 0
    torch_service->configure(cfg);
    const std::string mname = model_path;
    auto dtype = torch::kFloat16;

    torch::jit::script::Module module;
    // Deserialize the ScriptModule from a file using torch::jit::load().
    auto start = std::chrono::high_resolution_clock::now();
    //torch::Device device = torch::Device(torch::kCUDA,0);
    //module = 
    // (mname, device);
    //module.to(at::kCUDA, dtype);
    torch::TensorOptions options = torch::TensorOptions().dtype(dtype);
    torch::Tensor iten = torch::rand({1, 3, 800, 600}, options);
    std::vector<torch::IValue> itens {iten};
    //auto otens = module.forward(itens).toTensor();
    std::vector<torch::IValue> inputs = {iten};
    ITorchTensorSet::pointer iitens = to_itensor(inputs);
    auto oitens = torch_service->forward(iitens);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "timing: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

    return 0;
}