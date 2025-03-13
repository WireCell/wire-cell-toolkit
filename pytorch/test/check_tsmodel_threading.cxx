#include "WireCellPytorch/Torch.h"  // One-stop header.

#include <iostream>
#include <memory>
#include <thread>
#include <chrono>

void hello(int& threadid)
{
    std::cout << "Hello World! Thread ID, " << threadid << std::endl;
    return;
}

namespace TorchScript {
    void forward(int threadid, torch::jit::script::Module& module)
    {
        // Create a vector of inputs.
        std::vector<torch::jit::IValue> inputs;
        inputs.push_back(torch::ones({1, 3, 800, 600}, torch::dtype(torch::kFloat32).device(torch::kCUDA, 0)));
        // inputs.push_back(torch::rand({1, 3, 800, 600}, torch::dtype(torch::kFloat32).device(torch::kCPU, 0)));

        float total_wait_time = 0;
        int unit_wait_time = 500;  // ms
        // Execute the model and turn its output into a tensor.
        for (int i = 0; i < 10; ++i) {
            bool success = false;
            while (!success) {
                try {
                    at::Tensor output = module.forward(inputs).toTensor().cpu();
                    // std::cout << output[0][0][0][0] << '\n';
                    // std::cout << output.sizes() << '\n';
                    success = true;
                }
                catch (...) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(unit_wait_time));
                    total_wait_time += unit_wait_time;
                }
            }
        }

        std::cout << "th " << threadid << " waited: " << total_wait_time / 1000 << " sec." << std::endl;
        return;
    }
}  // namespace TorchScript

int main(int argc, const char* argv[])
{
    // if (argc != 3) {
    //   std::cout << "usage: threads <path-to-exported-script-module> nthreads\n";
    //   return -1;
    // }
    // const std::string mname = argv[1];
    // const int nthreads = atoi(argv[2]);

    std::cout << "WireCell::pytorch : test multi-threading of TorchScript Model\n";

    const std::string mname = "model.ts";
    const int nthreads = 4;

    torch::jit::script::Module module;
    try {
        // Deserialize the ScriptModule from a file using torch::jit::load().
        module = torch::jit::load(mname);
        std::cout << "model " << mname << " loaded\n";
    }
    catch (const torch::Error& e) {
        std::cout << "error loading the model\n";
        std::cout << "example models: https://github.com/HaiwangYu/Pytorch-UNet/tree/master/ts-model\n";
        return 1;
    }

    std::vector<std::thread> threads;
    for (int i = 0; i < nthreads; i++) {
        std::cout << "main() : creating thread, " << i << std::endl;
        // threads.push_back(std::thread(hello, std::ref(i)));
        threads.push_back(std::thread(TorchScript::forward, i, std::ref(module)));
    }
    for (int i = 0; i < nthreads; i++) {
        threads[i].join();
    }

    return 0;
}