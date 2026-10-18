/** Time rebaseline_masked() against rebaseline_zero() on a workflow-sized
 *  tensor (one pdhd APA: 2560 channels x 6000 ticks) on CPU and, when
 *  available, CUDA, and check that anchor_pad=0 reproduces rebaseline_zero().
 *
 *  Prints a small table; throws if the pad=0 result differs from the legacy
 *  result (max |diff| > 1e-3 on values up to ~800, i.e. float32 rounding).
 */
#include "WireCellSpng/Rebaseline.h"
#include "WireCellSpng/Torch.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace WireCell::SPNG;

namespace {

    // Positive-valued ROIs with gaps >= 2 and no exact zeros inside ROIs so the
    // legacy zero-run ROI finder and the mask agree.
    torch::Tensor make_workflow_like(int64_t nchans, int64_t nticks, int64_t seed)
    {
        torch::manual_seed(seed);
        auto tensor = torch::zeros({nchans, nticks}, torch::kFloat32);
        for (int64_t ch = 0; ch < nchans; ++ch) {
            const int64_t nrois = 3 + (ch % 5);
            for (int64_t roi = 0; roi < nrois; ++roi) {
                const int64_t base = 10 + roi * (nticks / (nrois + 1));
                const int64_t start = std::min<int64_t>(base + ((ch * 17 + roi * 31) % 200), nticks - 20);
                const int64_t len = 20 + ((ch * 13 + roi * 29) % 180);
                const int64_t stop = std::min<int64_t>(start + len, nticks - 3);
                auto vals = torch::randint(1, 800, {stop - start + 1}, torch::kInt32).to(torch::kFloat32)
                    + torch::rand({stop - start + 1}, torch::kFloat32) * 0.25;
                tensor.index_put_({ch, torch::indexing::Slice(start, stop + 1)}, vals);
            }
        }
        return tensor;
    }

    template<typename F>
    double timed(F&& f, const torch::Device& dev, int reps = 3)
    {
        for (int i = 0; i < 2; ++i) { auto t = f(); if (dev.is_cuda()) torch::cuda::synchronize(); }
        auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < reps; ++i) { auto t = f(); if (dev.is_cuda()) torch::cuda::synchronize(); }
        auto t1 = std::chrono::steady_clock::now();
        return std::chrono::duration<double>(t1 - t0).count() / reps;
    }

    void run(const torch::Device& dev)
    {
        auto wave = make_workflow_like(2560, 6000, 67890).to(dev);
        auto mask = wave != 0;
        auto masked = wave * mask;

        MaskRebaselineOptions legacy;   // pad 0, width 1
        MaskRebaselineOptions padded;
        padded.anchor_pad = 4;

        MaskRebaselineOptions padded_vec = padded;
        padded_vec.force_vectorized = true;

        const double t_zero = timed([&] { return rebaseline_zero(masked); }, dev);
        const double t_pad0 = timed([&] { return rebaseline_masked(wave, mask, legacy); }, dev);
        const double t_pad4 = timed([&] { return rebaseline_masked(wave, mask, padded); }, dev);
        const double t_vec = timed([&] { return rebaseline_masked(wave, mask, padded_vec); }, dev);

        auto ref = rebaseline_zero(masked);
        auto got = rebaseline_masked(wave, mask, legacy);
        const double maxdiff = (ref - got).abs().max().item<double>();

        std::cout << std::left << std::setw(8) << dev.str()
                  << " rebaseline_zero " << std::fixed << std::setprecision(1) << t_zero * 1e3 << " ms"
                  << " | masked pad=0 " << t_pad0 * 1e3 << " ms"
                  << " | masked pad=4 " << t_pad4 * 1e3 << " ms"
                  << " | masked pad=4 torch-vectorized " << t_vec * 1e3 << " ms"
                  << " | pad=0 max|diff| vs legacy " << std::setprecision(5) << maxdiff << "\n";
        if (maxdiff > 1e-3) {
            throw std::runtime_error("rebaseline_masked(pad=0) differs from rebaseline_zero on " + dev.str());
        }
    }
}

int main()
{
    std::cout << "2560 x 6000 tensor, mean of 3 calls after 2 warm-ups\n";
    run(torch::kCPU);
    if (torch::cuda::is_available()) {
        run(torch::Device(torch::kCUDA));
    }
    else {
        std::cout << "CUDA not available, skipping\n";
    }
    return 0;
}
