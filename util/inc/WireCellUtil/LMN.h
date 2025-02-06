/// Functions related to LMN resampling.
///
/// See LMN paper for details

#ifndef WIRECELLUTIL_LMN
#define WIRECELLUTIL_LMN

#include "WireCellUtil/Array.h"
#include <cmath>

namespace WireCell::LMN {

    /// Greated common divisor of two floating  point values
    ///
    /// The returned value divides both a and b to give an integer value within
    /// given error.  The ratio of numbers a and b must be rational.
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
        return (N-2)/2;         // event
    }

    // Return number divisible by Nrat that is equal or minimally larger than N.
    size_t nbigger(size_t N, size_t Nrat) {
        if (! N%Nrat) return N;
        return Nrat * ( N/Nrat + 1);
    }

    /// Return a new array with size Nr on axis.
    Array::array_xxf resize(const Array::array_xxf& in, size_t Nr,
                            size_t axis=1);
    std::vector<float> resize(const std::vector<float>& in, size_t Nr);

    void fill_constant(std::vector<float>::iterator begin,
                       std::vector<float>::iterator end,
                       float value = 0);
    void fill_linear(std::vector<float>::iterator begin,
                     std::vector<float>::iterator end,
                     float first, float last);


    /// Resample a complex array interpreted as a frequency-domain spectrum
    /// along the given axis so that it has Nr samples.  The layout is assumed
    /// to by nyquist-centered, ie the zero'th element represents zero
    /// frequency.  If the resampling is an upsampling, additional frequency
    /// samples of zero amplitude are inserted above the initial Nyquist
    /// frequency.  If a downsampling, samples around the Nyquist frequency are
    /// removedd.
    Array::array_xxc resample(const Array::array_xxc& in, size_t Nr, size_t axis = 1);
    std::vector<std::complex<float>> resample(const std::vector<std::complex<float>>& in, size_t Nr);


}
#endif
