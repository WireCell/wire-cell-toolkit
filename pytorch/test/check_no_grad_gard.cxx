#include <torch/script.h>  // One-stop header.
#include <chrono>
#include "WireCellUtil/ExecMon.h"
#include "WireCellUtil/Logging.h"
using namespace WireCell;
using spdlog::info;
/**
 * to use single CPU thread, export OMP_NUM_THREADS=1 when running the program
 */
int main(int argc, const char* argv[])
{
    auto x = torch::tensor({1.}, torch::requires_grad());
    {
        auto doubler = [](torch::Tensor x) {
            return x * 2;
        };
        auto z = doubler(x);
        std::cout << z.requires_grad() << std::endl;  // prints `true`
    }
    torch::NoGradGuard no_grad1;
    {
        auto doubler = [](torch::Tensor x) {
            return x * 2;
        };
        auto z = doubler(x);
        std::cout << z.requires_grad() << std::endl;  // prints `false`
    }
    {
        torch::NoGradGuard no_grad2;
        auto y = x * 2;
        std::cout << y.requires_grad() << std::endl;  // prints `false`
    }
    return 0;
}