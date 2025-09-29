#include "WireCellSpng/DFT.h"
#include "WireCellUtil/Persist.h"

namespace WireCell::SPNG {


    DftSpeedData load_dft_data(torch::Device device,
                               const std::string& filename)
    {
        std::string fname = "fft_";
        if (device == torch::kCPU) {
            fname += "cpu_ms";
        }
        else {
            fname += "gpu_ms";
        }

        auto dat = Persist::load(filename);
        return DftSpeedData {
            convert<std::vector<double>>(dat[fname]),
            convert<std::vector<int>>(data["sizes"])
        };
    }

    DftSpeedData bound_dft_data(const DftSpeedData& data, int64_t min_size, int64_t max_size)
    {
        DftSpeedData ret;
        const size_t nmeas = data.durations.size();
        for (size_t ind=0; ind<nmeas; ++ind) {

            const int64_t this_size = data.sizes[ind];
            if (this_size < min_size || this_size > max_size) {
                continue;
            }
            ret.sizes.push_back(this_size);
            ret.durations.push_back(data.durations[ind]);
        }
        return ret;        
    }


    int64_t faster_dft_measured(const DftSpeedData& data, int64_t orig_size, double orig_duration)
    {
        double best_duration = 0;
        int64_t best_size = 0;

        const size_t nmeas = durations.size();
        for (size_t ind=0; ind<nmeas; ++ind) {
            if (sizes[ind] < orig_size) {
                continue;
            }
            if (sizes [ind] == orig_size && orig_duration < 0) {
                orig_duration = durations[ind];
                continue;
            }
            
            if (!best_size || best_duration > durations[ind]) {
                best_size = sizes[ind];
                best_duration = durations[ind];
            }
        }
        if (orig_duration < 0) {
            raise<ValueError>("faster_dft_measured requires an original duration for size %d", orig_size);
        }

        if (orig_duration > best_duration) {
            return best_size;
        }

        return orig_size;
    }


    int64_t largest_prime_factor(int64_t n) {
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

    int64_t faster_dft_primes(int64_t orig_size, int64_t max_prime)
    {
        while (true) {          
            int64_t lpf = largest_prime_factor(orig_size);
            if (lpf <= max_prime) {
                // guaranteed to eventually break, right?.... 
                break;
            }
            ++orig_size;
        }
        return orig_size;
    }


} // namespace


