/** TorchLMN.h

    This provides the LMN "rational" size determination function and FOURIER
    space resampling functions.

    This was converted from WCT's WireCellUtil/LMN.h with help from Gemini to
    switch from Eigen arrays to Torch tensors.  The std::vector and some
    vestigial functions are ignored.

    Note: this also fixes Nyquist bin handling which is currently still needed
    in the original LMN.h.
*/
#ifndef WIRECELL_SPNG_LMN_TORCHLMN
#define WIRECELL_SPNG_LMN_TORCHLMN

#include <torch/torch.h>
#include <cmath>
#include <vector>
#include <complex>

namespace WireCell::SPNG::LMN {

    /// Greated common divisor of two floating  point values
    ///
    /// The returned value divides both a and b to give an integer value within
    /// given error. The ratio of numbers a and b must be rational.
    double gcd(double a, double b, double eps=1e-6);


    /// Return the "half size" number of samples in a spectrum of full size N.
    /// This excludes the "zero frequency" sample and the "Nyquist bin", if one
    /// exists.
    inline
    size_t nhalf(size_t N) {
        if (N%2) {
            return (N-1)/2;     // odd
        }
        return (N-2)/2;         // even
    }

    // Return number divisible by Nrat that is equal or minimally larger than N.
    inline
    size_t nbigger(size_t N, size_t Nrat) {
        if (N%Nrat == 0) return N;
        return Nrat * ( N/Nrat + 1);
    }

    /// Return a minimum sampling size that allows a rational resampling from
    /// sampling period Ts to sampling period Tr within error eps.  Zero is
    /// returned if error is encountered.  A tensor must be padded to a size
    /// that is a multiple of the rational size.
    size_t rational(double Ts, double Tr, double eps=1e-6);

    /// Return a new tensor with size Nr on axis.
    /// This implements simple truncation/zero-padding in the time/spatial domain.
    torch::Tensor resize(const torch::Tensor& in, int64_t Nr, int64_t axis=1);


    /// Resample a complex tensor interpreted as a frequency-domain spectrum
    /// along the given axis so that it has Nr samples.
    /// This uses Torch tensor manipulation (narrow/pad).
    torch::Tensor resample(const torch::Tensor& fourier_spectrum, int64_t Nr, int64_t axis = 1);

    /// One must assume an interpretation of the original sampling to properly
    /// normalize the resampled value.
    enum class Normalization {
        kUnknown=0,
        kInterpolation,         // for sampling of instantaneous function
        kIntegral,              // for integral over sample period
        kEnergy,                // for samples of probability amplitude
    };

    /// Return a new interval-space tensor that is a resampling of input
    /// interval-space tensor from sampling period Ts to sampling period Tr
    /// along given axis.
    ///
    /// This function wraps up the above functions.
    ///
    /// The size of the resampled dimension will be rational_size*(Tr/Ts) where
    /// rational_size comes from applying padding implied by rational(Ts, Tr).
    ///
    torch::Tensor resample_interval(const torch::Tensor& interval, double Ts, double Tr,
                                    int64_t axis = 1,
                                    Normalization ni = Normalization::kInterpolation,
                                    double eps=1e-6);


}
#endif
