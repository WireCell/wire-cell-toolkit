#include <torch/script.h>  // One-stop header.
#include "WireCellUtil/Eigen.h"
#include <chrono>
#include "WireCellUtil/ExecMon.h"
#include "WireCellUtil/Logging.h"
#include "WireCellPytorch/Util.h"
using namespace WireCell;
using spdlog::info;
/**
 * to use single CPU thread, export OMP_NUM_THREADS=1 when running the program
*/
int main(int argc, const char* argv[])
{
    if(argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <nchunks>\n";
        return 1;
    }
    const size_t nchunks = atoi(argv[1]);
    size_t nsecsleep = 0;
    if (argc > 2) {
        nsecsleep = atoi(argv[2]);
    }

    // Disable gradient computation, use like a mutex lock
    torch::NoGradGuard no_grad;
    // if (false) {
    //     torch::Tensor a = torch::rand({1, 2, 3, 3});
    //     std::cout << "Tensor shape: " << a.sizes() << std::endl;
    //     std::cout << a << std::endl;
    //     // Extract the desired sub-tensor using the slice method
    //     torch::Tensor sub_tensor = a.index({torch::indexing::Slice(), torch::indexing::Slice(1, 2), torch::indexing::Slice(), torch::indexing::Slice()});
    //     // Print the shape of the extracted tensor
    //     std::cout << "Sub-tensor shape: " << sub_tensor.sizes() << std::endl;
    //     std::cout << sub_tensor << std::endl;
    // }
    {
        Eigen::MatrixXf ch_eigen = Eigen::MatrixXf::Identity(4, 4);
        std::cout << "matrix: \n" << ch_eigen << std::endl;
        std::vector<torch::Tensor> ch;
        for (unsigned int i = 0; i < 2; ++i) {
            ch.push_back(torch::from_blob(ch_eigen.data(), {ch_eigen.cols(), ch_eigen.rows()}));
        }
        auto img = torch::stack(ch, 0);
        auto batch = torch::stack({torch::transpose(img, 1, 2)}, 0);
        std::cout << "batch: \n" << batch << std::endl;
        auto chunks = batch.chunk(4, 2);
        std::vector<torch::Tensor> outputs;
        for (auto chunk : chunks) {
            std::cout << "ichunk size: " << chunk.sizes() << std::endl;
            std::cout << "ichunk: \n" << chunk << std::endl;
            std::vector<torch::IValue> itens {chunk};
            auto iitens = Pytorch::to_itensor(itens);
            auto oitens = iitens;
            torch::Tensor ochunk = Pytorch::from_itensor({oitens}).front().toTensor().cpu();
            std::cout << "ochunk size: " << ochunk.sizes() << std::endl;
            std::cout << "ochunk: \n" << ochunk << std::endl;
            outputs.push_back(ochunk.clone());
        }
        for (auto output : outputs) {
            std::cout << "output size: " << output.sizes() << std::endl;
            std::cout << "output: \n" << output << std::endl;
        }
        auto output = torch::cat(outputs, 2);
        std::cout << "output size: " << output.sizes() << std::endl;
        std::cout << "output: \n" << output << std::endl;
        Eigen::Map<Eigen::ArrayXXf> out_e(output[0][0].data<float>(), output.size(3), output.size(2));
        std::cout << "out_e: \n" << out_e << std::endl;
    }
    return 0;

    std::cout << "WireCell::pytorch : test loading TorchScript Model\n";

    const std::string mname = "model.ts";
    auto dtype = torch::kFloat32;

    ExecMon em("initialize");
    torch::jit::script::Module module;
    // Deserialize the ScriptModule from a file using torch::jit::load().
    auto start = std::chrono::high_resolution_clock::now();
    module = torch::jit::load(mname);
    std::cout << em("load model") << std::endl;
    module.to(at::kCPU, dtype);
    info(em("to device"));
    std::cout << em("to device") << std::endl;
    torch::TensorOptions options = torch::TensorOptions().dtype(dtype);
    torch::Tensor iten = torch::rand({1, 3, 2000, 340}, options);
    std::cout << em("make input tensor") << std::endl;
    auto chunks = iten.chunk(nchunks, 2);
    std::cout << em("chunking") << std::endl;
    std::vector<torch::Tensor> outputs;
    // use these sleeps to identify the places where the program is at
    std::this_thread::sleep_for(std::chrono::seconds(nsecsleep));
    for (auto chunk : chunks) {
        std::cout << "chunk size: " << chunk.sizes() << std::endl;
        std::vector<torch::IValue> itens {chunk};
        auto otens = module.forward(itens).toTensor().cpu();
        std::this_thread::sleep_for(std::chrono::seconds(nsecsleep));
        std::cout << "otens size: " << chunk.sizes() << std::endl;
        outputs.push_back(otens);
    }
    std::this_thread::sleep_for(std::chrono::seconds(nsecsleep));
    torch::Tensor combined = torch::cat(outputs, 2);
    std::cout << "Combined size: " << combined.sizes() << std::endl;
    std::cout << em("inference") << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "timing: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

    return 0;
}