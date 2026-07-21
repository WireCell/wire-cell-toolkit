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


      Most developers of INode should have code like:

      @code{cpp}
      // Define a class member variable;
      FasterDftSize m_fds;

      // Use it:
      size_t better_size = m_fds(original_size);
      @endcode

      By default this will give you a cached version of the best-effort size
      resolver from faster_dft() for CPU.  See the class for options.

 */

#ifndef WIRECELL_SPNG_DFT
#define WIRECELL_SPNG_DFT

#include "WireCellSpng/Torch.h"

#include <vector>

namespace WireCell::SPNG {

    /// Interface with method that returns a possibly larger and possibly faster
    /// DFT size.
    ///
    /// Implementations are safe to call in a MT context but generally must be
    /// created in a ST context.
    ///
    /// However, implementations may take longer to calculate faster DFT size
    /// than it would take to simply use your original DFT size.  Use the
    /// caching FasterDftSize to ameliorate this cost when the original size
    /// does not change from call to call.
    struct FasterDftSizeBase {

        
        FasterDftSizeBase() = default;
        virtual ~FasterDftSizeBase() = default;

        using pointer = std::shared_ptr<const FasterDftSizeBase>;

        virtual int64_t operator()(int64_t orig_size) const = 0;
    };

    /// A faster DFT size resolver with caching based on calling another one.
    class FasterDftSize : public FasterDftSizeBase {
    public:

        /// This will initialize with a default resolver
        /// If you want to use another resolver, use reset().
        FasterDftSize(torch::Device device = torch::kCPU);
            
        virtual ~FasterDftSize() = default;

        /// Call to get a possibly faster, possibly larger size.
        virtual int64_t operator()(int64_t orig_size) const;

        /// Use a non-default resolver.
        void reset(FasterDftSizeBase::pointer fds);

    private:

        // For now, we have a simple, last-used cache.
        FasterDftSizeBase::pointer m_fds;
        mutable int64_t m_last_orig{0}, m_last_faster{0};
    };

    /// Return a dft faster size resolver by making a best effort to find the
    /// best resolver.
    ///
    /// Most users should call this or use FasterDftSize.
    ///
    /// This will first try to find the default file of measurements and fall
    /// back to using the small primes heuristic.
    ///
    /// See also make_faster_dft_size_measured() and make_faster_dft_size_primes().
    FasterDftSizeBase::pointer make_faster_dft_size(torch::Device device = torch::kCPU);


    /// The default filename holding dft size measurements.
    inline const std::string faster_dft_size_default_measurement_file = "torch_dft_measured.json";

    /// Return a FasterDftSize based on a data file holding measurements.
    ///
    /// See spng/test/check_dft_measured.cxx for one way to produce this file.
    ///
    /// This function is not thread safe, but the returned interface is.
    ///
    /// It will return nullptr if data file is not located.
    FasterDftSizeBase::pointer make_faster_dft_size_measured(
        torch::Device device = torch::kCPU,
        const std::string& filename = faster_dft_size_default_measurement_file);

    /// Return a FasterDftSize that implements the "next small prime factor
    /// size" heuristic.
    FasterDftSizeBase::pointer make_faster_dft_size_primes(int64_t max_prime=7);
    

    /// Return a faster size for fft() from measurements.
    ///
    /// This will return the value in sizes greater than or equal to orig_size
    /// that gives a lesser duration.
    ///
    /// If orig_size is not in the data, the duration of the first larger size
    /// is used as a proxy.  This can lead to overly large or even slower sizes
    /// in the case that the original size is best.  To protect the user a
    /// little bit from this, if orig_size is a power of two then that size is
    /// simply returned as that always beat the neighboring sizes.  This is not
    /// often but not always true when smallest prime factor of 3 or larger.
    ///
    /// This is thread safe.
    ///
    /// Caution: this function may scan the entire pair of vectors.  Caller
    /// should pre-truncate vectors to reduce scan time and/or to assure a
    /// result is bounded to some maximum size.
    int64_t faster_dft_size(const std::vector<int64_t>& sizes,
                            const std::vector<float>& durations,
                            int64_t orig_size);


    /// Return the largest prime factor of n.
    int64_t largest_prime_factor(int64_t n);

    /// Return a faster size for fft() from measurements with "small primes"
    /// heuristic.
    ///
    /// This will return the first size found that is greater or equal to the
    /// orig_size and which has no prime factor greater than max_prime.
    ///
    /// This is a reasonably good heuristic for the case that orig_size has
    /// larger primes.  In such cases, the resulting small-prime size will
    /// generally be rather faster.  However, sometimes yet more speed up is
    /// possible by yet larger small-prime sizes.  This is especially true when
    /// the DFT function will run in a single thread.  Measured data is required
    /// to achieve this level of optimization.
    int64_t faster_dft_size(int64_t orig_size, int64_t max_prime=7);
    
    
}

#endif

