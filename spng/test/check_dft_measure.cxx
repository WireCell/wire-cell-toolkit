#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <iomanip> // For std::fixed, std::setprecision
#include <functional> // For std::function

#include <torch/torch.h>
#include <torch/fft.h>

// for cli
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/CLI11.hpp"


int verbose = 0;                // please forgive my sin

// Helper function to create input tensors based on requirements for each FFT variant.
// - dim: The base dimension size (e.g., N for a 1D array, or N for an NxN 2D array).
// - device: The torch::Device (CPU or CUDA) where the tensor should be created.
// - is_1d: True if the tensor is for a 1D FFT, false for 2D.
// - requires_complex_input: True if the function expects a complex input tensor (kComplexFloat).
// - is_irfft_style_input_shape: True for irfft/irfft2, where the input tensor shape is
//                                 typically derived from the output of rfft/rfft2 (e.g., [dim/2 + 1]).
torch::Tensor create_input_tensor(
    int64_t dim,
    torch::Device device,
    bool is_1d,
    bool requires_complex_input,
    bool is_irfft_style_input_shape
) {
    torch::Dtype dtype = requires_complex_input ? torch::kComplexFloat : torch::kFloat;
    torch::TensorOptions options = torch::TensorOptions().dtype(dtype).device(device);

    int64_t half = dim/2+1;
    if (half < 1) half = 1;

    if (is_1d) {
        if (is_irfft_style_input_shape) {
            // Input for irfft: complex, shape typically [dim/2 + 1]
            // We use std::max to ensure dim/2 + 1 is at least 1 for small 'dim'
            return torch::zeros(half, options);
        } else {
            // Input for fft, rfft, ifft: shape [dim]
            return torch::zeros({dim}, options);
        }
    } else { // 2D
        if (is_irfft_style_input_shape) {
            // Input for irfft2: complex, shape typically [dim, dim/2 + 1]
            return torch::zeros({dim, half}, options);
        } else {
            // Input for fft2, rfft2, ifft2: shape [dim, dim]
            return torch::zeros({dim, dim}, options);
        }
    }
}

// Function to dynamically calculate the number of repetitions for a benchmark.
// It aims for a minimum total measurement time to ensure reasonable accuracy,
// especially for very fast operations.
long long calculate_repetitions(const std::function<void()>& func_to_bench, torch::Device device) {
    // Warm-up call (not timed) to ensure resources are ready and avoid startup overhead.
    func_to_bench();
    if (device.is_cuda()) {
        torch::cuda::synchronize();
    }

    // Measure a single call to estimate its duration.
    auto start_trial = std::chrono::high_resolution_clock::now();
    func_to_bench();
    if (device.is_cuda()) {
        torch::cuda::synchronize();
    }
    auto end_trial = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(end_trial - start_trial).count();

    const double TARGET_TOTAL_MS = 100.0; // Target a total measurement time of at least 100ms
    long long repetitions = 1;

    if (elapsed_ms > 0 && elapsed_ms < TARGET_TOTAL_MS) {
        repetitions = static_cast<long long>(std::ceil(TARGET_TOTAL_MS / elapsed_ms));
    } else if (elapsed_ms == 0) {
        // Fallback for extremely fast operations where elapsed_ms is 0.
        // This ensures a large enough number of repetitions for very small tensors.
        repetitions = 10000;
    }
    return std::max(1LL, repetitions); // Ensure at least 1 repetition.
}

// Struct to hold information for each FFT function to be benchmarked.
struct FFTFunctionInfo {
    std::string name;
    bool is_1d; // True for 1D functions (fft, rfft, ifft, irfft)
    bool requires_complex_input; // True if the input tensor should be kComplexFloat
    bool is_irfft_style_input_shape; // True if input shape is like [dim/2 + 1] or [dim, dim/2 + 1]
    
    // The actual function call, defined as a lambda. It accepts the pre-created
    // input tensor and the original 'dim' to correctly set 'n' or 's' parameters.
    std::function<void(const torch::Tensor& input, int64_t dim_size)> call;
};

std::vector<int64_t> getPrimeFactors(int64_t n) {
    std::vector<int64_t> factors;

    // Handle the factor 2 separately
    while (n % 2 == 0) {
        factors.push_back(2);
        n /= 2;
    }

    // Iterate for odd numbers starting from 3 up to the square root of n
    for (int64_t i = 3; i <= std::sqrt(n); i += 2) {
        while (n % i == 0) {
            factors.push_back(i);
            n /= i;
        }
    }

    // If n is a prime number greater than 2 (after all divisions)
    if (n > 2) {
        factors.push_back(n);
    }

    return factors;
}
int64_t maxPrimeFactor(int64_t n) {
    int64_t max_prime=0;

    // Handle the factor 2 separately
    while (n % 2 == 0) {
        max_prime = 2;
        n /= 2;
    }

    // Iterate for odd numbers starting from 3 up to the square root of n
    for (int64_t i = 3; i <= std::sqrt(n); i += 2) {
        while (n % i == 0) {
            if (i > max_prime) {
                max_prime = i;
            }
            n /= i;
        }
    }

    // If n is a prime number greater than 2 (after all divisions)
    if (n > 2) {
        if (n > max_prime) {
            max_prime = n;
        }
    }

    return max_prime;
}


std::vector<int64_t> sizes_with_small_primes(int64_t min_size, int64_t max_size, int64_t max_prime=7)
{
    std::vector<int64_t> sizes;
    for (int64_t dim=min_size; dim <= max_size; ++dim) {
        if (maxPrimeFactor(dim) <= max_prime) {
            sizes.push_back(dim);
        }
    }
    return sizes;
}

std::vector<int64_t> sizes_in_range(int64_t min_size, int64_t max_size)
{
    std::vector<int64_t> sizes(max_size - min_size);
    std::iota(sizes.begin(), sizes.end(), min_size);
    return sizes;
}

struct consumer_func {
    
    virtual void operator()(int64_t size,
                            const std::vector<double>& cpu_times,
                            const std::vector<double>& gpu_times) = 0;
};

struct csv_file : public consumer_func {
    std::string output_filename;
    std::ofstream out;
    std::vector<std::string> fnames;
    int64_t nlines{0};

    csv_file(const std::string& output_filename,
             const std::vector<std::string>& function_names)
        : output_filename(output_filename)
        , out(output_filename)
        , fnames(function_names)
    {
        if (!out.is_open()) {
            throw std::runtime_error("can not open file for writing: " + output_filename);
        }
    }

    ~csv_file() {
        out.close();
        if (verbose > 0 ) {
            std::cerr << "Results written to " << output_filename << std::endl;
        }

    }

    void operator()(int64_t size,
                    const std::vector<double>& cpu_times,
                    const std::vector<double>& gpu_times) {
        
        if (!nlines) {
            // header line
            out << "size";
            if (cpu_times.size()) {
                for (const auto& fname : fnames) {
                    out << "," << fname << "_cpu_ms";
                }
            }
            if (gpu_times.size()) {
                for (const auto& fname : fnames) {
                    out << "," << fname << "_gpu_ms";
                }
            }
            out << "\n";
        }
        ++nlines;

        out << size;
        if (cpu_times.size()) {
            for (double time_ms : cpu_times) {
                out << "," << std::fixed << std::setprecision(6) << time_ms;
            }
        }
        if (gpu_times.size()) {
            for (double time_ms : gpu_times) {
                out << "," << std::fixed << std::setprecision(6) << time_ms;
            }
        }
        out << "\n";
        out.flush();
    }
};



std::vector<double> do_test_size_device(const std::vector<FFTFunctionInfo>& functions_to_test,
                                        int64_t dim,
                                        torch::Device device)
{
    std::vector<double> times;
    for (const auto& func_info : functions_to_test) {
        // Create input tensor for the current function and device.
        torch::Tensor input_tensor = create_input_tensor(
            dim, device, func_info.is_1d, func_info.requires_complex_input, func_info.is_irfft_style_input_shape);

        // Dynamically determine repetitions for this specific benchmark.
        long long repetitions = calculate_repetitions([&]() {
            func_info.call(input_tensor, dim);
        }, device);

        if (verbose > 1) {
            std::cerr << repetitions << "x: " << func_info.name << " on " << device << std::endl;
        }

        // Perform the timed benchmark.
        auto start_time = std::chrono::high_resolution_clock::now();
        for (long long i = 0; i < repetitions; ++i) {
            func_info.call(input_tensor, dim);
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        double elapsed_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
        times.push_back(elapsed_ms / repetitions);
    }
    return times;

}

int do_tests(const std::vector<FFTFunctionInfo>& functions_to_test,
             const std::vector<int64_t>& sizes_to_test,
             const std::vector<torch::Device>& devices,
             consumer_func& output)
{
    // 6. Main loop for dimension sizes.
    for (int64_t dim : sizes_to_test) {
        if (verbose > 0) {
            std::cerr << "Benchmarking size: " << dim << std::endl;
        }

        std::vector<double> cpu_times;
        std::vector<double> gpu_times;

        for (const auto dev : devices) {
            if (dev == torch::kCPU) {
                cpu_times = do_test_size_device(functions_to_test, dim, dev);
            }
            else {
                gpu_times = do_test_size_device(functions_to_test, dim, dev);                
            }
        }

        // // --- CPU Benchmarks ---
        // torch::Device cpu_device(torch::kCPU);
        // for (const auto& func_info : functions_to_test) {
        //     // Create input tensor for the current function and device.
        //     torch::Tensor input_tensor = create_input_tensor(
        //         dim, cpu_device, func_info.is_1d, func_info.requires_complex_input, func_info.is_irfft_style_input_shape);

        //     // Dynamically determine repetitions for this specific benchmark.
        //     long long repetitions = calculate_repetitions([&]() {
        //         func_info.call(input_tensor, dim);
        //     }, cpu_device);

        //     if (verbose > 1) {
        //         std::cerr << "\tcpu: " << repetitions << "x: " << func_info.name << std::endl;
        //     }

        //     // Perform the timed benchmark.
        //     auto start_time = std::chrono::high_resolution_clock::now();
        //     for (long long i = 0; i < repetitions; ++i) {
        //         func_info.call(input_tensor, dim);
        //     }
        //     auto end_time = std::chrono::high_resolution_clock::now();
        //     double elapsed_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
        //     cpu_times.push_back(elapsed_ms / repetitions);
        // }

        // // --- GPU Benchmarks (if CUDA is available) ---
        // if (cuda_is_available) {
        //     torch::Device gpu_device(torch::kCUDA);
        //     for (const auto& func_info : functions_to_test) {
        //         // Create input tensor for the current function and device.
        //         torch::Tensor input_tensor = create_input_tensor(
        //             dim, gpu_device, func_info.is_1d, func_info.requires_complex_input, func_info.is_irfft_style_input_shape);
                
        //         // Dynamically determine repetitions for this specific benchmark.
        //         long long repetitions = calculate_repetitions([&]() {
        //             func_info.call(input_tensor, dim);
        //         }, gpu_device);

        //         if (verbose > 1) {
        //             std::cerr << "\tgpu: " << repetitions << "x: " << func_info.name << std::endl;
        //         }

        //         // Perform the timed benchmark.
        //         auto start_time = std::chrono::high_resolution_clock::now();
        //         for (long long i = 0; i < repetitions; ++i) {
        //             func_info.call(input_tensor, dim);
        //         }
        //         // Synchronize CUDA operations to ensure all GPU work is complete before stopping timer.
        //         torch::cuda::synchronize();
        //         auto end_time = std::chrono::high_resolution_clock::now();
        //         double elapsed_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
        //         gpu_times.push_back(elapsed_ms / repetitions);
        //     }
        // }
 
        output(dim, cpu_times, gpu_times);
    }

    return 0;
}

struct json_file : public consumer_func {
    std::string output_filename;
    std::vector<std::string> fnames;
    Json::Value obj = Json::objectValue;

    json_file(const std::string& output_filename,
              std::vector<std::string> fnames)
        : output_filename(output_filename)
        , fnames(fnames)
    {
    }

    ~json_file() {
        WireCell::Persist::dump(output_filename, obj, true);
        if (verbose > 0) {
            std::cerr << "Results written to " << output_filename << std::endl;
            std::cerr << "Place in WIRECELL_PATH location to use\n";
        }
    }

    void operator()(int64_t size,
                    const std::vector<double>& cpu_times,
                    const std::vector<double>& gpu_times) {
        
        obj["sizes"].append(size);

        const size_t nfuncs = fnames.size();
        if (cpu_times.size()) {
            for (size_t ind=0; ind<nfuncs; ++ind) {
                obj[fnames[ind] + "_cpu_ms"].append(cpu_times[ind]);
            }
        }
        if (gpu_times.size()) {
            for (size_t ind=0; ind<nfuncs; ++ind) {
                obj[fnames[ind] + "_gpu_ms"].append(gpu_times[ind]);
            }
        }
    }

};

std::vector<std::string> function_names(const std::vector<FFTFunctionInfo> funcs)
{
    std::vector<std::string> fnames;
    for (const auto& func : funcs) {
        fnames.push_back(func.name);
    }
    return fnames;
}
    

std::vector<FFTFunctionInfo> benchmark_functions() {
    std::vector<FFTFunctionInfo> functions_to_test = {
        // 1D FFTs
        {"fft",    true,  true,  false, [](const torch::Tensor& input, int64_t dim) { 
            // n parameter is 'dim' (length of transformed axis of the output).
            auto output = torch::fft::fft(input, dim, -1, c10::nullopt); 
            // Avoid compiler optimizing out the call
            (void)output; 
        }},
    };
    return functions_to_test;
}

std::vector<FFTFunctionInfo> exhaustive_functions() {
    // 4. Define the list of FFT functions to benchmark.
    // Each entry specifies the function's properties and its calling lambda.
    std::vector<FFTFunctionInfo> functions_to_test = {
        // 1D FFTs
        {"fft",    true,  true,  false, [](const torch::Tensor& input, int64_t dim) { 
            // n parameter is 'dim' (length of transformed axis of the output).
            auto output = torch::fft::fft(input, dim, -1, c10::nullopt); 
            // Avoid compiler optimizing out the call
            (void)output; 
        }},
        {"rfft",   true,  false, false, [](const torch::Tensor& input, int64_t dim) { 
            // n parameter is 'dim' (length of transformed axis of the output).
            auto output = torch::fft::rfft(input, dim, -1, c10::nullopt); 
            (void)output; 
        }},
        {"ifft",   true,  true,  false, [](const torch::Tensor& input, int64_t dim) { 
            // n parameter is 'dim' (length of transformed axis of the output).
            auto output = torch::fft::ifft(input, dim, -1, c10::nullopt); 
            (void)output; 
        }},
        {"irfft",  true,  true,  true,  [](const torch::Tensor& input, int64_t dim) { 
            // n parameter is 'dim' (length of transformed axis of the real output).
            auto output = torch::fft::irfft(input, dim, -1, c10::nullopt); 
            (void)output; 
        }},

        // 2D FFTs
        {"fft2",   false, true,  false, [](const torch::Tensor& input, int64_t dim) { 
            // s parameter is {dim, dim} (shape of last 2 dimensions of the output).
            auto output = torch::fft::fft2(input);
            (void)output; 
        }},
        {"rfft2",  false, false, false, [](const torch::Tensor& input, int64_t dim) { 
            // s parameter is {dim, dim} (shape of last 2 dimensions of the output).
            auto output = torch::fft::rfft2(input);
            (void)output; 
        }},
        {"ifft2",  false, true,  false, [](const torch::Tensor& input, int64_t dim) { 
            // s parameter is {dim, dim} (shape of last 2 dimensions of the output).
            auto output = torch::fft::ifft2(input);
            (void)output; 
        }},
        {"irfft2", false, true,  true,  [](const torch::Tensor& input, int64_t dim) { 
            // s parameter is {dim, dim} (shape of last 2 dimensions of the real output).
            auto output = torch::fft::irfft2(input);
            (void)output; 
        }}
    };
    return functions_to_test;
}

int main(int argc, char* argv[]) {
    CLI::App app{"Measure DFT times.\nWrite .csv to make plots with check_dft_measure.py.\nWrite faster_dft_size.json in WIRECELL_JSON for use by WCT/SPNG faster_dft_size() method.\nSet OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 to assure torch uses 1 CPU thread.\n"};

    std::string output_filename="/dev/stdout";
    std::string size_set = "primes";
    std::string func_set = "benchmark";
    std::string devices = "cpu,gpu";
    /// These defaults assume primes and benchmark.
    /// With range and/or exhaustive, this can take a long time to run.
    /// Adjust to taste or have patience.
    int64_t min_size=2, max_size=1<<14, max_prime=7;


    app.add_option("-s,--sizes", size_set,
                   "How to select sizes (def='primes', also 'range' and best set min/max)"
        )->type_size(1)->allow_extra_args(false);
    // fixme: add custom list of function names
    app.add_option("-f,--funcs", func_set,
                   "How to select functions (def='benchmark', also 'exhaustive' and best set min/max)"
        )->type_size(1)->allow_extra_args(false);
    app.add_option("-m,--min-size", min_size,
                   "Minimum size to test (def=2)"
        )->type_size(1)->allow_extra_args(false);
    app.add_option("-d,--devices", devices,
                   "Devices to consider (def: 'cpu,gpu')"
        )->type_size(1)->allow_extra_args(false);
    app.add_option("-M,--max-size", max_size,
                   "Minimum size to test (def=2^14)"
        )->type_size(1)->allow_extra_args(false);
    app.add_option("-p,--max-prime", max_prime,
                   "Minimum prime (def=7)"
        )->type_size(1)->allow_extra_args(false);
    app.add_option("-v,--verbose", verbose,
                   "Set verbosity level (def=0, set to, 1 or 2 for more)"
        )->type_size(1)->allow_extra_args(false);
    app.add_option("output", output_filename,
                   "Write results in .csv (def) or .json file (required, use '-' for stdout)"
        )->required();


    CLI11_PARSE(app, argc, argv);

    std::vector<FFTFunctionInfo> functions_to_test;
    if (func_set == "benchmark") {
        functions_to_test = benchmark_functions();
    }
    else if (func_set == "exhaustive") {
        functions_to_test = exhaustive_functions();
    }
    else {
        std::cerr << "unsupported function selector: " << func_set << "\n";
        return 1;
    }
    auto fnames = function_names(functions_to_test);

    std::vector<int64_t> sizes_to_test;
    if (size_set == "primes") {
        sizes_to_test = sizes_with_small_primes(min_size, max_size, max_prime);
    }
    else if (size_set == "range") {
        sizes_to_test = sizes_in_range(min_size, max_size);
    }
    else {
        std::cerr << "unsupported size selector: " << size_set << "\n";
        return 1;
    }

    if (output_filename == "-") {
        output_filename = "/dev/stdout";
    }

    std::unique_ptr<consumer_func> consumer;
    if (WireCell::String::endswith(output_filename, ".json")) {
        consumer = std::make_unique<json_file>(output_filename, fnames);
    }
    else { // assume csv
        consumer = std::make_unique<csv_file>(output_filename, fnames);        
    }
    
    std::vector<torch::Device> torch_devices;
    for (const auto& one : WireCell::String::split(devices, ",")) {
        if (one == "cpu") {
            torch_devices.push_back(torch::kCPU);
            continue;
        }
        if (one == "gpu" || one == "cuda") {
            torch_devices.push_back(torch::kCUDA);
            continue;
        }
        std::cerr << "unsupported torch device: " << one << "\n";
        return 1;
    }

    int rc = do_tests(functions_to_test, sizes_to_test, torch_devices, *consumer);
    consumer.reset();
    return rc;
}
