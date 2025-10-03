// This is some test from gemini to try to estimate fft cost from first
// principles and then measure it.  Conclusion: estimation does not match measurement..

#include <vector>
#include <algorithm> // For std::max
#include <limits>    // For std::numeric_limits
#include <cmath>     // For std::log2, std::pow
#include <map>       // For std::map
#include <chrono>    // For high_resolution_clock
#include <iostream>  // For std::cout, std::cerr
#include <iomanip>   // For std::setw, std::fixed, std::setprecision

// LibTorch headers
#include <torch/torch.h>

// Helper list of small primes that FFT libraries typically optimize for.
const std::vector<int> default_small_fft_primes = {2, 3, 5, 7, 11, 13, 17, 19, 23};

// --- Forward Declarations ---
bool is_fft_smooth(int n, const std::vector<int>& small_primes_list);
int get_optimal_fft_size(int N, const std::vector<int>& small_primes_list);
double get_smoothness_penalty_multiplier(int n, const std::vector<int>& small_primes_list);
double estimate_relative_dft_cost(int N, int M, const std::vector<int>& small_primes_list);
double measure_relative_dft_cost(int N, int M, const std::vector<int>& small_primes_list, int num_repetitions_override = 0);
void analyze_fft_sizes_in_range(int Nmax, const std::vector<int>& target_small_primes_list);


// --- Existing Functions from previous responses (updated signatures) ---

// Checks if a number's prime factors are all within the 'small_primes_list'.
bool is_fft_smooth(int n, const std::vector<int>& small_primes_list) {
    if (n <= 0) return false;
    if (n == 1) return true;

    int temp_n = n;
    for (int p : small_primes_list) {
        while (temp_n > 0 && temp_n % p == 0) {
            temp_n /= p;
        }
        if (temp_n == 1) {
            return true;
        }
        if (p > temp_n / p) { // Optimization
            break;
        }
    }
    return temp_n == 1;
}

// Given an input size N, suggests a size M >= N that is likely to result
// in a faster DFT computation by having "smoother" prime factors (according to small_primes_list).
int get_optimal_fft_size(int N, const std::vector<int>& small_primes_list) {
    if (N <= 0) return 1;

    int max_search_limit = std::max(N + 256, N + (N / 4)); // Max 25% padding or 256, whichever is larger.
    if (max_search_limit < N) max_search_limit = std::numeric_limits<int>::max(); // Overflow check, be safe.

    for (int M = N; M <= max_search_limit; ++M) {
        if (is_fft_smooth(M, small_primes_list)) {
            return M;
        }
    }

    // Fallback: If no "smooth" number (according to small_primes_list) is found
    // within the padding limit, the safest fallback is to pad to the next power of 2.
    // This is because powers of 2 are generally the fastest, even if '2' wasn't
    // explicitly in the 'small_primes_list' provided for the *search*.
    int next_power_of_2 = 1;
    while (next_power_of_2 < N) {
        if (next_power_of_2 > std::numeric_limits<int>::max() / 2) { // Prevent overflow
            return N; // Cannot find a larger power of 2
        }
        next_power_of_2 <<= 1;
    }
    return next_power_of_2;
}

// Estimates a "smoothness penalty multiplier" for a given size.
// A higher multiplier means less efficient computation.
double get_smoothness_penalty_multiplier(int n, const std::vector<int>& small_primes_list) {
    if (n <= 1) return 1.0;

    double total_penalty_multiplier = 1.0;
    int temp_n = n;

    // Build a specific penalty map for the provided small_primes_list
    // Primes *not* in small_primes_list will be treated as "large" for this function.
    std::map<int, double> current_prime_penalties;
    // Standard penalties for common optimized primes
    std::map<int, double> standard_penalties = {
        {2, 1.0}, {3, 1.1}, {5, 1.2}, {7, 1.3}, {11, 1.5},
        {13, 1.6}, {17, 1.7}, {19, 1.8}, {23, 1.9}
    };

    // Populate current_prime_penalties based on the provided small_primes_list
    for (int p : small_primes_list) {
        if (standard_penalties.count(p)) {
            current_prime_penalties[p] = standard_penalties[p];
        } else {
            current_prime_penalties[p] = 2.5; // Default penalty for an unrecognized "small" prime
        }
    }

    for (int p : small_primes_list) { // Iterate through the *provided* small_primes_list
        if (p * p > temp_n && temp_n > 1) break;

        int count = 0;
        while (temp_n > 0 && temp_n % p == 0) {
            temp_n /= p;
            count++;
        }
        if (count > 0) {
            total_penalty_multiplier *= std::pow(
                current_prime_penalties.at(p), // Use .at() for safety, though p should be in map
                count
            );
        }
    }

    // If 'temp_n' is still greater than 1, it means N has a prime factor
    // that was NOT in the 'small_primes_list' provided to this function.
    // This implies using less efficient algorithms like Bluestein's or Rader's.
    if (temp_n > 1) {
        total_penalty_multiplier *= 7.0; // Strong heuristic for large prime factor overhead
    }
    return total_penalty_multiplier;
}

// Estimates the relative computational cost of a DFT of size N compared to M.
double estimate_relative_dft_cost(int N, int M, const std::vector<int>& small_primes_list) {
    if (N <= 0 || M <= 0) return std::numeric_limits<double>::quiet_NaN();

    double cost_N_log = (N == 1) ? 1.0 : static_cast<double>(N) * std::log2(static_cast<double>(N));
    double cost_M_log = (M == 1) ? 1.0 : static_cast<double>(M) * std::log2(static_cast<double>(M));

    double penalty_N = get_smoothness_penalty_multiplier(N, small_primes_list);
    double penalty_M = get_smoothness_penalty_multiplier(M, small_primes_list);

    double total_cost_N = cost_N_log * penalty_N;
    double total_cost_M = cost_M_log * penalty_M;

    if (total_cost_M == 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    return total_cost_N / total_cost_M;
}


// --- New Function: Measure Relative DFT Cost using LibTorch ---
double measure_relative_dft_cost(int N, int M, const std::vector<int>& small_primes_list, int num_repetitions_override) {
    if (N <= 0 || M <= 0) {
        // std::cerr << "Error: Invalid input size for DFT measurement." << std::endl; // Suppress for table output
        return std::numeric_limits<double>::quiet_NaN();
    }

    torch::Device device = torch::cuda::is_available() ? torch::kCUDA : torch::kCPU;

    auto measure_single_dft_time = [&](int size) {
        if (size == 1) return 1e-9; // Assign a very small non-zero time for trivial DFT

        at::Tensor input_tensor = torch::rand({size}, torch::TensorOptions().dtype(torch::kComplexFloat).device(device));

        int num_repetitions;
        if (num_repetitions_override > 0) {
            num_repetitions = num_repetitions_override;
        } else {
            // Original dynamic logic for individual measurements
            if (size < 50) num_repetitions = 5000;
            else if (size < 500) num_repetitions = 1000;
            else if (size < 5000) num_repetitions = 100;
            else num_repetitions = 10;
        }
        num_repetitions = std::max(5, num_repetitions); // Ensure at least 5 repetitions

        // Warm-up run to ensure everything is loaded into cache/initialized
        torch::fft::fft(input_tensor);
        if (device.is_cuda()) {
            torch::cuda::synchronize();
        }

        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < num_repetitions; ++i) {
            torch::fft::fft(input_tensor);
        }
        if (device.is_cuda()) {
            torch::cuda::synchronize(); // Ensure all GPU operations are complete
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        return duration.count() / num_repetitions; // Average time per FFT
    };

    double time_N = measure_single_dft_time(N);
    double time_M = measure_single_dft_time(M);

    if (time_M == 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    return time_N / time_M;
}

// --- New Function: Analyze FFT sizes in a range ---
void analyze_fft_sizes_in_range(int Nmax, const std::vector<int>& target_small_primes_list) {
    if (Nmax <= 0) {
        std::cerr << "Error: Nmax must be greater than 0." << std::endl;
        return;
    }
    
    // Determine device once for the entire analysis
    torch::Device device_for_analysis = torch::cuda::is_available() ? torch::kCUDA : torch::kCPU;

    std::cout << "\n--- FFT Analysis for N=1 to " << Nmax << " ---" << std::endl;
    std::cout << "  (Optimizing M with factors: {";
    bool first = true;
    for (int p : target_small_primes_list) {
        if (!first) std::cout << ", ";
        std::cout << p;
        first = false;
    }
    std::cout << "})" << std::endl;
    std::cout << "  (Measurements on device: " << (device_for_analysis.is_cuda() ? "CUDA" : "CPU") << ")" << std::endl;
    std::cout << "  (Measured FFTs repeated " << 50 << " times for speed)" << std::endl;


    std::cout << std::fixed << std::setprecision(3);
    std::cout << std::setw(6) << "N"
              << std::setw(6) << "M"
              << std::setw(18) << "Est. Ratio (N/M)"
              << std::setw(18) << "Meas. Ratio (N/M)"
              << std::endl;
    std::cout << std::string(48, '-') << std::endl; // Separator

    int range_analysis_repetitions = 50; // Use a fixed smaller number for range analysis

    for (int N = 1; N <= Nmax; ++N) {
        int M = get_optimal_fft_size(N, target_small_primes_list);
        double estimated_ratio = estimate_relative_dft_cost(N, M, target_small_primes_list);
        double measured_ratio = measure_relative_dft_cost(N, M, target_small_primes_list, range_analysis_repetitions);

        std::cout << std::setw(6) << N
                  << std::setw(6) << M
                  << std::setw(18) << estimated_ratio
                  << std::setw(18) << measured_ratio
                  << std::endl;

        if (N % (Nmax / 10 + 1) == 0 && Nmax > 10) { // Progress indicator
            std::cout << "  (Progress: " << N << "/" << Nmax << ")" << std::endl;
        }
    }
    std::cout << std::string(48, '-') << std::endl;
}

// --- Main function to demonstrate and compare ---
int main() {
    std::cout << "--- FFT Size Optimization Examples (Individual Cases) ---" << std::endl;

    // Example 1: N=11 (prime) vs M=12 (2^2 * 3)
    int N1 = 11;
    int M1 = get_optimal_fft_size(N1, default_small_fft_primes); // Use default full list for general optimal
    double estimated_cost1 = estimate_relative_dft_cost(N1, M1, default_small_fft_primes);
    double measured_cost1 = measure_relative_dft_cost(N1, M1, default_small_fft_primes);
    std::cout << "--- Case 1 (N=" << N1 << " vs M=" << M1 << ") ---" << std::endl;
    std::cout << "  Estimated Cost(N)/Cost(M): " << std::fixed << std::setprecision(3) << estimated_cost1 << std::endl;
    std::cout << "  Measured Cost(N)/Cost(M):  " << std::fixed << std::setprecision(3) << measured_cost1 << std::endl;
    std::cout << "  Comment: N=11 is prime (high penalty). M=12 is very smooth. Expect N1 to be much slower than M1." << std::endl;

    // Example 2: N=97 (prime) vs M=98 (2 * 7^2)
    int N2 = 97;
    int M2 = get_optimal_fft_size(N2, default_small_fft_primes); // Use default full list
    double estimated_cost2 = estimate_relative_dft_cost(N2, M2, default_small_fft_primes);
    double measured_cost2 = measure_relative_dft_cost(N2, M2, default_small_fft_primes);
    std::cout << "\n--- Case 2 (N=" << N2 << " vs M=" << M2 << ") ---" << std::endl;
    std::cout << "  Estimated Cost(N)/Cost(M): " << std::fixed << std::setprecision(3) << estimated_cost2 << std::endl;
    std::cout << "  Measured Cost(N)/Cost(M):  " << std::fixed << std::setprecision(3) << measured_cost2 << std::endl;
    std::cout << "  Comment: N=97 is prime. M=98 is smooth. Expect N2 to be much slower than M2." << std::endl;
    
    // Example 3: N=1023 (3 * 11 * 31) vs M=1024 (2^10)
    int N3 = 1023;
    int M3 = get_optimal_fft_size(N3, default_small_fft_primes); // Use default full list
    double estimated_cost3 = estimate_relative_dft_cost(N3, M3, default_small_fft_primes);
    double measured_cost3 = measure_relative_dft_cost(N3, M3, default_small_fft_primes);
    std::cout << "\n--- Case 3 (N=" << N3 << " vs M=" << M3 << ") ---" << std::endl;
    std::cout << "  Estimated Cost(N)/Cost(M): " << std::fixed << std::setprecision(3) << estimated_cost3 << std::endl;
    std::cout << "  Measured Cost(N)/Cost(M):  " << std::fixed << std::setprecision(3) << measured_cost3 << std::endl;
    std::cout << "  Comment: N=1023 has a large prime factor 31 (not in default_small_fft_primes). M=1024 is a power of 2. Expect N3 to be much slower than M3, despite M3 being slightly larger." << std::endl;

    // Example 4: N=100 (2^2 * 5^2) vs M=96 (2^5 * 3) - Both smooth
    int N4 = 100;
    int M4 = 96; // M=96 is specifically chosen here for comparison.
    double estimated_cost4 = estimate_relative_dft_cost(N4, M4, default_small_fft_primes);
    double measured_cost4 = measure_relative_dft_cost(N4, M4, default_small_fft_primes);
    std::cout << "\n--- Case 4 (N=" << N4 << " vs M=" << M4 << ") ---" << std::endl;
    std::cout << "  Estimated Cost(N)/Cost(M): " << std::fixed << std::setprecision(3) << estimated_cost4 << std::endl;
    std::cout << "  Measured Cost(N)/Cost(M):  " << std::fixed << std::setprecision(3) << measured_cost4 << std::endl;
    std::cout << "  Comment: Both are smooth. M4 is slightly smaller and more biased towards factors of 2. Expect ratio to be around 1.0 or slightly less/more depending on specific library optimizations." << std::endl;
    
    // Example 5: N=5003 (a prime number) vs M=5005 (5 * 7 * 11 * 13)
    int N5 = 5003;
    int M5 = get_optimal_fft_size(N5, default_small_fft_primes); // Should find 5005
    double estimated_cost5 = estimate_relative_dft_cost(N5, M5, default_small_fft_primes);
    double measured_cost5 = measure_relative_dft_cost(N5, M5, default_small_fft_primes);
    std::cout << "\n--- Case 5 (N=" << N5 << " vs M=" << M5 << ") ---" << std::endl;
    std::cout << "  Estimated Cost(N)/Cost(M): " << std::fixed << std::setprecision(3) << estimated_cost5 << std::endl;
    std::cout << "  Measured Cost(N)/Cost(M):  " << std::fixed << std::setprecision(3) << measured_cost5 << std::endl;
    std::cout << "  Comment: N=5003 is prime (very high penalty). M=5005 is very smooth. Expect N5 to be significantly slower than M5." << std::endl;
    

    // --- Range Analysis ---
    int Nmax_for_analysis = 5000; // Set a small Nmax for demonstration; larger values will take time.
    analyze_fft_sizes_in_range(Nmax_for_analysis, {3, 5, 7}); // Analyze with specific small primes for M
    
    return 0;
}
