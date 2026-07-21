#include "WireCellSpng/DFT.h"
#include "WireCellSpng/Util.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Persist.h"       // for load()
#include "WireCellUtil/Configuration.h" // for convert()

#include <boost/container_hash/hash.hpp>

namespace WireCell::SPNG {


    // Main user touching class

    FasterDftSize::FasterDftSize(torch::Device device)
        : m_fds(make_faster_dft_size(device)) { }

    int64_t FasterDftSize::operator()(int64_t orig_size) const
    {
        if (!m_fds) {
            raise<ValueError>("no FasterDftSize has been set, fix the code.");
        }
        if (orig_size == m_last_orig) { return m_last_faster; }
        m_last_orig = orig_size;
        m_last_faster = (*m_fds)(orig_size);
        return m_last_faster;
    }

    void FasterDftSize::reset(FasterDftSizeBase::pointer fds) { m_fds = fds; }


    // best effort factory
    FasterDftSizeBase::pointer make_faster_dft_size(torch::Device device)
    {
        FasterDftSizeBase::pointer ret = make_faster_dft_size_measured(device);
        if (ret) return ret;
        return make_faster_dft_size_primes();
    }


    // Measured version
    struct FasterDftSizeMeasured : public FasterDftSizeBase {

        FasterDftSizeMeasured() = default;
        virtual ~FasterDftSizeMeasured() = default;

        std::vector<int64_t> sizes;
        std::vector<float> durations;
        
        virtual int64_t operator()(int64_t orig_size) const {
            /// Remember, not such a cheap call....
            return faster_dft_size(sizes, durations, orig_size); // measured
        }
    };
    FasterDftSizeBase::pointer make_faster_dft_size_measured(
        torch::Device device, const std::string& datafile)
    {
        // not-thread-safe cache by device and data file name.
        static std::unordered_map<size_t, std::shared_ptr<const FasterDftSizeMeasured>> cache;

        std::size_t key = 0;
        boost::hash_combine(key, SPNG::to_string(device));
        boost::hash_combine(key, datafile);
        
        auto it = cache.find(key);
        if (it != cache.end()) { // hit
            return it->second;
        }

        Configuration data;
        try {
            data = Persist::load(datafile);
        }
        catch (const WireCell::IOError&) {
             cache[key] = nullptr; // save a load() next time we are called
            return nullptr;
        }
        
        // Note, this will share the same data between different callers but it
        // does duplicate the sizes vector when one caller asks for cpu and
        // another for gpu from the same datafile.
        auto res = std::make_shared<FasterDftSizeMeasured>();
        res->sizes = convert<std::vector<int64_t>>(data["sizes"]);
        if (device == torch::kCPU) {
            res->durations = convert<std::vector<float>>(data["fft_cpu_ms"]);
        }
        else {
            res->durations = convert<std::vector<float>>(data["fft_gpu_ms"]);
        }
        cache[key] = res;
        return res;
    }

    // The small-primes heuristic.
    struct FasterDftSizePrimes : public FasterDftSizeBase
    {
        FasterDftSizePrimes(int64_t max_prime) : max_prime(max_prime) {}
        virtual ~FasterDftSizePrimes() = default;

        const int64_t max_prime;
        
        virtual int64_t operator()(int64_t orig_size) const {
            return faster_dft_size(orig_size, max_prime);
        }
    };
    FasterDftSizeBase::pointer make_faster_dft_size_primes(int64_t max_prime)
    {
        // We just intern the max_prime config.  This is pretty small and simple
        // so don't bother with any caching.
        return std::make_shared<FasterDftSizePrimes>(max_prime);
    }



    int64_t faster_dft_size(const std::vector<int64_t>& sizes,
                            const std::vector<float>& durations,
                            int64_t orig_size)
    {
        if (2 == largest_prime_factor(orig_size)) {
            // Safety net to avoid unnecessary inflation of size.
            // We can not always count on this heuristic for 3 or larger.
            return orig_size;
        }

        float orig_duration = 0;
        float best_duration = 0;
        int64_t best_size = 0;

        const size_t nmeas = sizes.size();
        for (size_t ind=0; ind<nmeas; ++ind) {

            // scan past smaller sizes
            if (sizes[ind] < orig_size) {
                continue;
            }

            // grab smallest that is not smaller than orig_size as proxy
            if (orig_duration == 0) { 
                orig_size = sizes[ind];
                orig_duration = durations[ind];
                continue;
            }

            if (!best_size || best_duration > durations[ind]) {
                best_size = sizes[ind];
                best_duration = durations[ind];
            }
        }
        if (orig_duration < 0) {
            raise<ValueError>("faster_dft_size requires an original duration for size %d", orig_size);
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

    int64_t faster_dft_size(int64_t orig_size, int64_t max_prime)
    {
        while (true) {          
            int64_t lpf = largest_prime_factor(orig_size);
            if (lpf <= max_prime) {
                // guaranteed to eventually break out, right?.... 
                break;
            }
            ++orig_size;
        }
        return orig_size;
    }


} // namespace


