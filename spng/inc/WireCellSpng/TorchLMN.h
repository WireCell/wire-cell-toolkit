// TorchLMN.h converted from LMN.h with help from Gemini 
#ifndef WIRECELL_SPNG_LMN_TORCHLMN
#define WIRECELL_SPNG_LMN_TORCHLMN

#include <torch/torch.h>
#include <cmath>
#include <vector>
#include <complex>

namespace WireCell::SPNG::LMN {

    // Helper functions operating on standard types (size_t, double, vector) are preserved.

    /// Greated common divisor of two floating  point values
    ///
    /// The returned value divides both a and b to give an integer value within
    /// given error. The ratio of numbers a and b must be rational.
    double gcd(double a, double b, double eps=1e-6);

    /// Return a minimum sampling size that allows a rational resampling from
    /// sampling period Ts to sampling period Tr within error eps.  Zero is
    /// returned if error is exceeded.
    size_t rational(double Ts, double Tr, double eps=1e-6);

    /// Return the "half size" number of samples in a spectrum of full size N.
    /// This excludes the "zero frequency" sample and the "Nyquist bin", if one
    /// exists.
    size_t nhalf(size_t N) {
        if (N%2) {
            return (N-1)/2;     // odd
        }
        return (N-2)/2;         // even
    }

    // Return number divisible by Nrat that is equal or minimally larger than N.
    size_t nbigger(size_t N, size_t Nrat) {
        if (N == 0) return Nrat;
        if (! (N%Nrat)) return N;
        return Nrat * ( N/Nrat + 1);
    }

    /// Return a new tensor with size Nr on axis.
    /// This implements simple truncation/zero-padding in the time/spatial domain.
    torch::Tensor resize(const torch::Tensor& in, int64_t Nr, int64_t axis=1);
    std::vector<float> resize(const std::vector<float>& in, size_t Nr);

    void fill_constant(std::vector<float>::iterator begin,
                       std::vector<float>::iterator end,
                       float value = 0);
    void fill_linear(std::vector<float>::iterator begin,
                     std::vector<float>::iterator end,
                     float first, float last);


    /// Resample a complex tensor interpreted as a frequency-domain spectrum
    /// along the given axis so that it has Nr samples.
    /// This uses Torch tensor manipulation (narrow/pad).
    torch::Tensor resample(const torch::Tensor& in, int64_t Nr, int64_t axis = 1);
    std::vector<std::complex<float>> resample(const std::vector<std::complex<float>>& in, size_t Nr);


}
#endif
