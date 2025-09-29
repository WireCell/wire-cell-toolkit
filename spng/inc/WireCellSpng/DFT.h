/** Functions related to applying Discrete Fourier Transform functions from torch.

    The faster_dft_*() functions take an original tensor dimension size and
    return a size that is no smaller and which may have a DFT that is not slower
    than a DFT of original size.

    This "faster dft" business is rife with caveats.

    - The function tested matters (fft vs rfft vs iffy vs fft2, etc).  The rest
      here considers only fft().

    - CPU beats GPU for sizes smaller than about 512.

    - Starting with 2048, CPU speeds for powers of 2 are not faster than numbers
      with larger small prime factors.

    - GPU speed is relatively flat as a function of sizes comprised of small
      prime factors except for isolated cases.  Isolated cases are a few times
      slower than neighbors.  Density of slower cases increases starting around
      8192, about 50-50 at 10k and 100% about 13k.

    - Single-thread CPU generally shows N log N bands which are not driven just
      by largest prime factor.  It may be a function of the largest small prime
      optimization but torch shows 3 bands.  Bands represent 2-6x slow down
      between nearby fast and slow sizes.  If torch is allowed multi-threading
      CPU, these bands become scattered and N log N is broken at about 4k
      elements.


 */

#ifndef WIRECELL_SPNG_DFT
#define WIRECELL_SPNG_DFT

#include "WireCellSpng/Torch.h"

#include <vector>

namespace WireCell::SPNG {

    /// Group a duration for a DFT with its size.
    struct DftSpeedData {
        const std::vector<double> durations,
        const std::vector<int64_t> sizes,
    };


    DftSpeedData 


    /// Load a duration data file from WIRECELL_PATH.
    ///
    /// This file can be produced by:
    ///
    ///     OMP_NUM_THREADS=1 MKL_NUM_THREADS=1  \
    ///     ./build/spng/check_dft_measure torch_dft_duration.json
    ///
    ///
    /// Caveats:
    /// - It is possible to make such JSON files with only CPU, only GPU or both.
    ///   Getting it wrong gives you an exception.
    /// - The results to not record exact device.  Your CPU/GPU may differ from what was used.
    DftSpeedData load_dft_data(torch::Device device = torch::kCPU,
                               const std::string& filename = "torch_dft_duration.json");

    /// Return new data that has sizes in the inclusive bounds
    DftSpeedData bound_dft_data(const DftSpeedData& data, int64_t min_size, int64_t max_size);

    /// Return a faster size for fft() from measurements.
    ///
    /// This will return the value in sizes greater than or equal to orig_size
    /// that gives a lesser duration.  The orig_size must be in the data or its
    /// original duration given or else a ValueError is thrown.
    int64_t faster_dft_measured(const DftSpeedData& data, int64_t orig_size, double orig_duration=-1);


    /// Return the largest prime factor of n.
    int64_t largest_prime_factor(int64_t n);

    /// Return the closest size greater or equal to original with a largest
    /// prime no larger than max_prime.  This is a pretty good heuristic that
    /// should return a size not "too" far from the original that is faster,
    /// especially for original sizes that have large prime factors.  However,
    /// some yet larger sizes with small primes can be faster still.  Also, it
    /// is not always true that the fastest in the neighborhood has the smallest
    /// largest prime factor.  Eg, sometimes a somewhat larger "seven" is faster
    /// than a smaller "five".
    int64_t faster_dft_primes(int64_t orig_size, int64_t max_prime=7);
    
    
}

#endif

