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

    /// Return a minimum sampling size that allows a rational resampling from
    /// sampling period Ts to sampling period Tr within error eps.  Zero is
    /// returned if error is exceeded.
    size_t rational(double Ts, double Tr, double eps=1e-6);

    /// Return a new tensor with size Nr on axis.
    /// This implements simple truncation/zero-padding in the time/spatial domain.
    torch::Tensor resize(const torch::Tensor& in, int64_t Nr, int64_t axis=1);


    /// Resample a complex tensor interpreted as a frequency-domain spectrum
    /// along the given axis so that it has Nr samples.
    /// This uses Torch tensor manipulation (narrow/pad).
    torch::Tensor resample(const torch::Tensor& in, int64_t Nr, int64_t axis = 1);


}
#endif
