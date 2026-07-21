#include <torch/torch.h>

#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Timing {
    torch::Tensor output;
    double seconds{0.0};
};

void sync_if_cuda(const torch::Tensor& tensor)
{
    if (tensor.defined() && tensor.is_cuda()) {
        torch::cuda::synchronize();
    }
}

template<typename Function>
Timing time_call(Function&& function)
{
    auto t0 = std::chrono::steady_clock::now();
    auto output = function();
    sync_if_cuda(output);
    auto t1 = std::chrono::steady_clock::now();
    return {output, std::chrono::duration<double>(t1 - t0).count()};
}

int64_t modulo(int64_t value, int64_t modulus)
{
    const int64_t ret = value % modulus;
    return ret < 0 ? ret + modulus : ret;
}

std::vector<int64_t> make_permutation(int64_t ndim, int64_t dim)
{
    std::vector<int64_t> permutation;
    permutation.reserve(ndim);
    for (int64_t d = 0; d < ndim; ++d) {
        if (d != dim) {
            permutation.push_back(d);
        }
    }
    permutation.push_back(dim);
    return permutation;
}

torch::Tensor restore_shape(const torch::Tensor& batch_view,
                            const torch::Tensor& permuted_output,
                            const std::vector<int64_t>& permutation)
{
    std::vector<int64_t> reverse_permutation(permutation.size());
    for (int64_t k = 0; k < static_cast<int64_t>(permutation.size()); ++k) {
        reverse_permutation[permutation[k]] = k;
    }
    return batch_view.view(permuted_output.sizes()).permute(reverse_permutation).contiguous();
}

torch::Tensor original_rebaseline_zero(const torch::Tensor& tensor,
                                       int64_t dim = -1,
                                       int64_t consequtive_zeros = 2,
                                       int64_t min_roi_size = 1,
                                       int64_t shrink_size = 0,
                                       bool remove_small = false,
                                       bool remove_negative = false)
{
    torch::Tensor output = tensor.clone();
    dim = modulo(dim, output.dim());

    std::vector<int64_t> permutation;
    for (int64_t d = 0; d < output.dim(); ++d) {
        if (d != dim) permutation.push_back(d);
    }
    permutation.push_back(dim);
    torch::Tensor permuted_output = output.permute(permutation).contiguous();

    const int64_t nticks = output.size(dim);
    const int64_t nbatches = permuted_output.numel() / nticks;
    torch::Tensor batch_view = permuted_output.view({nbatches, nticks});

    for (int64_t i = 0; i < nbatches; ++i) {
        torch::Tensor wave = batch_view[i];

        std::vector<bool> is_sep(nticks, false);
        {
            int64_t k = 0;
            while (k < nticks) {
                if (wave[k].item<float>() == 0.0f) { // If the wave tensor lives in the device, there is a element transfer between device and host which is expensive due to synchronizations.
                    int64_t run_start = k;
                    while (k < nticks && wave[k].item<float>() == 0.0f) { //Same as above.
                        ++k;
                    }
                    if (k - run_start >= consequtive_zeros) {
                        for (int64_t m = run_start; m < k; ++m) {
                            is_sep[m] = true;
                        }
                    }
                } else {
                    ++k;
                }
            }
        }

        std::vector<std::pair<int64_t,int64_t>> rois;
        bool in_roi = false;
        int64_t roi_start = -1;
        for (int64_t k = 0; k < nticks; ++k) {
            if (!is_sep[k]) {
                if (!in_roi) { roi_start = k; in_roi = true; }
            } else {
                if (in_roi) { rois.push_back({roi_start, k-1}); in_roi = false; }
            }
        }
        if (in_roi) { rois.push_back({roi_start, nticks-1}); }

        for (auto [start_idx, end_idx] : rois) {
            start_idx += shrink_size;
            end_idx   -= shrink_size;

            if (start_idx > end_idx) continue;

            if (end_idx - start_idx <= min_roi_size - 1) {
                if (remove_small) {
                    const auto roi = torch::indexing::Slice(start_idx, end_idx + 1);
                    wave.index_put_({roi}, 0.0);
                }
                continue;
            }

            float start_val = wave[start_idx].item<float>();
            float end_val   = wave[end_idx].item<float>();
            const float length = static_cast<float>(end_idx - start_idx);
            const float slope  = (end_val - start_val) / length;
            torch::Tensor baseline = start_val + slope * torch::arange(0, end_idx - start_idx + 1,
                                                                        wave.options());
            const auto roi = torch::indexing::Slice(start_idx, end_idx + 1);
            wave.index_put_({roi}, wave.index({roi}) - baseline);
        }
    }

    if (remove_negative) {
        output.clamp_min_(0.0);
    }

    std::vector<int64_t> reverse_permutation(output.dim());
    for (int64_t k = 0; k < output.dim(); ++k) {
        reverse_permutation[permutation[k]] = k;
    }
    return batch_view.view(permuted_output.sizes()).permute(reverse_permutation).contiguous();
}

template<typename Wave>
void apply_rebaseline_to_roi_like_original(Wave&& wave,
                                           int64_t start_idx,
                                           int64_t end_idx,
                                           int64_t min_roi_size,
                                           int64_t shrink_size,
                                           bool remove_small)
{
    start_idx += shrink_size;
    end_idx -= shrink_size;

    if (start_idx > end_idx) {
        return;
    }

    const auto roi = torch::indexing::Slice(start_idx, end_idx + 1);
    if (end_idx - start_idx <= min_roi_size - 1) {
        if (remove_small) {
            wave.index_put_({roi}, 0.0);
        }
        return;
    }

    float start_val = wave[start_idx].template item<float>();
    float end_val = wave[end_idx].template item<float>();
    const float length = static_cast<float>(end_idx - start_idx);
    const float slope = (end_val - start_val) / length;
    auto baseline = start_val + slope * torch::arange(0, end_idx - start_idx + 1,
                                                      wave.options());
    wave.index_put_({roi}, wave.index({roi}) - baseline);
}

void new_rebaseline_zero_cuda_kernel(torch::Tensor& batch_view,
                                     int64_t nticks,
                                     int64_t nbatches,
                                     int64_t consequtive_zeros,
                                     int64_t min_roi_size,
                                     int64_t shrink_size,
                                     bool remove_small)
{
    auto zero_mask = batch_view.eq(0).to(torch::kUInt8).to(torch::kCPU).contiguous(); //masking in bulk and bulk copy of the tensor to the CPU
    const auto* zeros = zero_mask.data_ptr<uint8_t>();

    for (int64_t i = 0; i < nbatches; ++i) { // do the expensive element by element comparison in the CPU. 
        const auto* row_zeros = zeros + i * nticks;
        auto wave = batch_view[i];

        bool in_roi = false;
        int64_t roi_start = -1;

        int64_t k = 0;
        while (k < nticks) {
            if (row_zeros[k]) {
                const int64_t run_start = k;
                while (k < nticks && row_zeros[k]) {
                    ++k;
                }
                if (k - run_start >= consequtive_zeros) {
                    if (in_roi) {
                        apply_rebaseline_to_roi_like_original(wave, roi_start, run_start - 1,
                                                              min_roi_size, shrink_size, remove_small);
                        in_roi = false;
                    }
                }
                else if (!in_roi) {
                    roi_start = run_start;
                    in_roi = true;
                }
                continue;
            }

            if (!in_roi) {
                roi_start = k;
                in_roi = true;
            }
            ++k;
        }

        if (in_roi) {
            apply_rebaseline_to_roi_like_original(wave, roi_start, nticks - 1,
                                                  min_roi_size, shrink_size, remove_small);
        }
    }
}

torch::Tensor new_rebaseline_zero(const torch::Tensor& tensor,
                                  int64_t dim = -1,
                                  int64_t consequtive_zeros = 2,
                                  int64_t min_roi_size = 1,
                                  int64_t shrink_size = 0,
                                  bool remove_small = false,
                                  bool remove_negative = false)
{
    if (!tensor.is_cuda()) {
        return original_rebaseline_zero(tensor, dim, consequtive_zeros,
                                        min_roi_size, shrink_size,
                                        remove_small, remove_negative);
    }

    if (!tensor.numel()) {
        return tensor.clone();
    }

    dim = modulo(dim, tensor.dim());
    auto permutation = make_permutation(tensor.dim(), dim);

    torch::Tensor output = tensor.clone();
    torch::Tensor permuted_output = output.permute(permutation).contiguous();

    const int64_t nticks = output.size(dim);
    if (nticks == 0) {
        return restore_shape(permuted_output, permuted_output, permutation);
    }

    const int64_t nbatches = permuted_output.numel() / nticks;
    torch::Tensor batch_view = permuted_output.view({nbatches, nticks});

    new_rebaseline_zero_cuda_kernel(batch_view, nticks, nbatches, consequtive_zeros,
                                    min_roi_size, shrink_size, remove_small);

    if (remove_negative) {
        output.clamp_min_(0.0);
    }

    return restore_shape(batch_view, permuted_output, permutation);
}

void print_tensor_or_summary(const std::string& label,
                             const torch::Tensor& tensor,
                             bool verbose)
{
    auto cpu = tensor.to(torch::kCPU);
    if (verbose && cpu.numel() <= 200) {
        std::cout << label << ":\n" << cpu << "\n";
        return;
    }

    auto flat = cpu.reshape({-1});
    std::cout << label
              << ": shape=" << cpu.sizes()
              << " dtype=" << cpu.dtype()
              << " min=" << flat.min().item<float>()
              << " max=" << flat.max().item<float>()
              << " mean=" << flat.mean().item<float>()
              << " nonzero=" << flat.ne(0).sum().item<int64_t>()
              << "\n";
}

double max_abs_diff(const torch::Tensor& lhs, const torch::Tensor& rhs)
{
    auto diff = (lhs.to(torch::kCPU) - rhs.to(torch::kCPU)).abs();
    return diff.numel() ? diff.max().item<float>() : 0.0;
}

void compare_and_print(const std::string& label,
                       const Timing& original_cpu,
                       const Timing& new_cpu,
                       const Timing& new_device,
                       const Timing* original_device,
                       bool verbose)
{
    auto a = original_cpu.output.to(torch::kCPU);
    auto b = new_cpu.output.to(torch::kCPU);
    auto c = new_device.output.to(torch::kCPU);

    std::cout << "\n=== " << label << " ===\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "time original CPU: " << original_cpu.seconds << " s\n";
    std::cout << "time new CPU:      " << new_cpu.seconds << " s\n";
    std::cout << "time new device:   " << new_device.seconds << " s\n";
    if (original_device) {
        std::cout << "time original device: " << original_device->seconds << " s\n";
    }

    print_tensor_or_summary("original CPU", a, verbose);
    print_tensor_or_summary("new CPU", b, verbose);
    print_tensor_or_summary("new device", c, verbose);
    if (original_device) {
        print_tensor_or_summary("original device", original_device->output, verbose);
    }

    std::cout << "max_abs(original CPU - new CPU):    "
              << max_abs_diff(a, b) << "\n";
    std::cout << "max_abs(original CPU - new device): "
              << max_abs_diff(a, c) << "\n";
    if (original_device) {
        std::cout << "max_abs(original CPU - original device): "
                  << max_abs_diff(a, original_device->output) << "\n";
    }

    if (!torch::allclose(a, b, 1e-6, 1e-6, true)) {
        throw std::runtime_error(label + ": new CPU differs from original CPU");
    }
}

void run_case(const std::string& label,
              const torch::Tensor& input_cpu,
              int64_t dim = -1,
              int64_t consequtive_zeros = 2,
              int64_t min_roi_size = 1,
              int64_t shrink_size = 0,
              bool remove_small = false,
              bool remove_negative = false,
              bool verbose = true)
{
    auto original_cpu = time_call([&] {
        return original_rebaseline_zero(input_cpu, dim, consequtive_zeros,
                                        min_roi_size, shrink_size,
                                        remove_small, remove_negative);
    });
    auto new_cpu = time_call([&] {
        return new_rebaseline_zero(input_cpu, dim, consequtive_zeros,
                                   min_roi_size, shrink_size,
                                   remove_small, remove_negative);
    });

    Timing new_device;
    Timing original_device;
    Timing* original_device_ptr = nullptr;
    if (torch::cuda::is_available()) {
        auto input_cuda = input_cpu.to(torch::kCUDA);
        sync_if_cuda(input_cuda);
        new_device = time_call([&] {
            return new_rebaseline_zero(input_cuda, dim, consequtive_zeros,
                                       min_roi_size, shrink_size,
                                       remove_small, remove_negative);
        });
        original_device = time_call([&] {
            return original_rebaseline_zero(input_cuda, dim, consequtive_zeros,
                                            min_roi_size, shrink_size,
                                            remove_small, remove_negative);
        });
        original_device_ptr = &original_device;
    }
    else {
        std::cout << "\nCUDA unavailable; using CPU as the device path placeholder for "
                  << label << "\n";
        new_device = new_cpu;
    }

    compare_and_print(label, original_cpu, new_cpu, new_device, original_device_ptr, verbose);
}

torch::Tensor make_workflow_like_tensor(int64_t nchans, int64_t nticks, int64_t seed)
{
    torch::manual_seed(seed);
    auto tensor = torch::zeros({nchans, nticks}, torch::kFloat32);

    for (int64_t ch = 0; ch < nchans; ++ch) {
        const int64_t nrois = 3 + (ch % 5);
        for (int64_t roi_index = 0; roi_index < nrois; ++roi_index) {
            const int64_t base = 10 + roi_index * (nticks / (nrois + 1));
            const int64_t start = std::min<int64_t>(base + ((ch * 17 + roi_index * 31) % 200),
                                                    nticks - 20);
            const int64_t len = 20 + ((ch * 13 + roi_index * 29) % 180);
            const int64_t stop = std::min<int64_t>(start + len, nticks - 3);
            const int64_t width = stop - start + 1;

            auto integer_part = torch::randint(1, 800, {width}, torch::kInt32).to(torch::kFloat32);
            auto fractional_part = torch::rand({width}, torch::kFloat32) * 0.25;
            auto roi_values = integer_part + fractional_part;
            tensor.index_put_({ch, torch::indexing::Slice(start, stop + 1)}, roi_values);

            // Single exact zero inside some ROIs: this is a crossing, not a separator
            // when consequtive_zeros=2.
            if ((ch + roi_index) % 4 == 0 && width > 8) {
                tensor.index_put_({ch, start + width / 2}, 0.0f);
            }
        }
    }

    return tensor;
}

} // namespace

int main()
{
    try {
        const auto f32 = torch::TensorOptions().dtype(torch::kFloat32);
        run_case("1D bounded ROI",
                 torch::tensor({0, 0, 5, 10, 5, 0, 0}, f32));

        run_case("1D single zero crossing",
                 torch::tensor({5, 0, 8, 0, 0}, f32));

        run_case("1D shrink",
                 torch::tensor({0, 0, 5, 10, 15, 10, 5, 0, 0}, f32),
                 -1, 2, 1, 1);

        run_case("1D remove_small",
                 torch::tensor({0, 0, 5, 0, 0}, f32),
                 -1, 2, 1, 0, true);

        run_case("1D remove_negative",
                 torch::tensor({0, 0, 5, 3, 10, 0, 0}, f32),
                 -1, 2, 1, 0, false, true);

        run_case("2D dim=1",
                 torch::tensor({{0, 0, 5, 10, 5, 0, 0},
                                {0, 0, 8, 12, 8, 0, 0}}, f32),
                 1, 2);

        run_case("2D dim=0",
                 torch::tensor({{0, 5, 0},
                                {0, 10, 0},
                                {0, 5, 0},
                                {0, 0, 0},
                                {0, 0, 0}}, f32),
                 0, 2);

        run_case("workflow-like small 256x1500 dim=1",
                 make_workflow_like_tensor(256, 1500, 12345),
                 1, 2, 1, 0, false, false,
                 false);

        run_case("workflow-like full 2560x6000 dim=1",
                 make_workflow_like_tensor(2560, 6000, 67890),
                 1, 2, 1, 0, false, false,
                 false);

        std::cout<< " \nNote: Here original means the implementation of Rebaseline.cxx from: \n" ;
      std::cout<<" \nLink: https://github.com/WireCell/wire-cell-toolkit/commit/95da99676f9ad579da54adbea6b40ccf6f3b35d1#diff-f336204ad6d7172eda9af20f611637a11b4fa2a0cdf5a67e65f09f913a8aa587 \n";
      
        std::cout<<" \nNote: Here new means the implementation of Rebaseline.cxx from: \n" ;
      std::cout<<" \nLink: https://github.com/WireCell/wire-cell-toolkit/commit/9939e37ccfd27ba8f81fe331a27350da36347584 \n";
        std::cout << "\nNote: original implementation is timed on CUDA when CUDA is available.\n";
        std::cout << "The full workflow-like case may be slow because it intentionally exercises the old scalar .item() GPU path.\n";
    }
    catch (const std::exception& err) {
        std::cerr << "ERROR: " << err.what() << "\n";
        return 1;
    }
    return 0;
}
