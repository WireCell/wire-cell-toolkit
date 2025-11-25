/** TorchLMN.h

    This provides the LMN "rational" size determination function and FOURIER
    space resampling functions.

    This was converted from WCT's WireCellUtil/LMN.h with help from Gemini to
    switch from Eigen arrays to Torch tensors.  The std::vector and some
    vestigial functions are ignored.

    Note: this also fixes Nyquist bin handling which is currently still needed
    in the original LMN.h.

    For a node that uses this, see Resample.
*/
#ifndef WIRECELL_SPNG_LMN_TORCHLMN
#define WIRECELL_SPNG_LMN_TORCHLMN

#include <torch/torch.h>
#include <cmath>
#include <vector>
#include <complex>

namespace WireCell::SPNG::LMN {

    /**
    * @brief Calculates the Greatest Common Divisor (GCD) of two floating-point numbers.
    *
    * This function determines the largest floating-point value $d$ such that both
    * $a/d$ and $b/d$ are near-integer values, within the specified error tolerance ($\epsilon$).
    *
    * The method is an adaptation of the Euclidean algorithm for floating-point arithmetic.
    * Note that this calculation requires the ratio of the two input numbers, $a/b$, to be rational
    * for a mathematically meaningful result.
    *
    * @param a The first floating-point value.
    * @param b The second floating-point value.
    * @param eps The maximum allowable error tolerance for near-integer quotients. Defaults to 1e-6.
    * @return The greatest common divisor (GCD) of $a$ and $b$.
    */
    double gcd(double a, double b, double eps=1e-6);



    // Return number divisible by Nrat that is equal or minimally larger than N.
    inline size_t nbigger(size_t N, size_t Nrat)
    {
        if (N%Nrat == 0) return N;
        return Nrat * ( N/Nrat + 1);
    }

    /**
    * @brief Calculates the minimum sampling size required for rational resampling between two periods.
    *
    * This function determines the smallest integral size (N) such that resampling from the
    * source sampling period (\p Ts) to the target sampling period (\p Tr) can be
    * approximated by a rational factor (P/Q) with a relative error less than or equal to \p eps.
    *
    * This size is crucial because any input data structure (e.g., a tensor or array)
    * that is subject to this resampling must be padded or sized to a multiple of the
    * returned value to ensure the rational transformation is applied correctly and without
    * unacceptable truncation error.
    *
    * @param Ts The source sampling period (time per sample). Must be positive.
    * @param Tr The target sampling period (time per resampled unit). Must be positive.
    * @param eps The maximum allowable relative error tolerance when approximating the rational ratio $Tr/Ts$. Defaults to $10^{-6}$.
    * @return The minimum rational sampling size required for accurate resampling. Returns
    *         \c 0 if a suitable rational approximation cannot be found within the specified
    *         tolerance or if internal computation errors are encountered.
    */
    size_t rational(double Ts, double Tr, double eps=1e-6);

    /**
     * @brief Resamples a dimension of a COMPLEX, FREQUENCY-DOMAIN spectrum
     * tensor to a target number of samples.
     *
     * This function does NOT ACCEPT INTERVAL DOMAIN tensors.
     *
     * This function adjusts the spectral resolution of the input by
     * manipulating its Fourier components.  If the target size (Nr) is smaller
     * than the current size along the axis, the spectrum is narrowed (high
     * frequencies are truncated).  If Nr is larger, the spectrum is zero-padded
     * (which corresponds to sinc interpolation in the time domain).  The
     * implementation uses standard Torch tensor operations (narrow and pad).
     *
     * @param fourier_spectrum The complex tensor representing the frequency-domain spectrum.
     * @param Nr The target number of samples (size) along the specified axis.
     * @param axis The dimension along which resampling is performed. Defaults to last dimension.
     * @return The resampled complex tensor of size Nr along the specified axis.
     */
    torch::Tensor resample_fourier(const torch::Tensor& fourier_spectrum, int64_t Nr, int64_t axis = -1);

    /// One must assume an interpretation of the original sampling to properly
    /// normalize the resampled value.
    enum class Normalization {
        kUnknown=0,
        kInterpolation,         // for sampling of instantaneous function
        kIntegral,              // for integral over sample period
        kEnergy,                // for samples of probability amplitude

    };
    /// Convert a string describing a Normalization to a Normalization.
    Normalization norm(std::string name);

    /**
     * @brief Combine the primitive operations into a full LMN resampling.
     *
     * @param interval An interval-space tensor.
     * @param Ts The sampling period of the tensor dimension to resample.
     * @param Tr The sampling period target for the resampled dimension.
     * @param axis The dimension of the tensor to resample.  Defaults to last dimension.
     * @param ni The interpretation of the samples to assume when normalizing.
     * @param fourier_space If true, return the Fourier space representation, default will apply ifft() and return the interval-space representation.
     * @param eps A small number passed when determining the rational size.
     * @return The tensor with the given axis resampled.
     *
     * The size of the resampled dimension of the returned tensor may be larger
     * than interval.size(axis)*(Tr/Ts) due to padding to achieve a size that is
     * a multiple of the factor returned by rational().
     *
     */
    torch::Tensor resample_interval(const torch::Tensor& interval, double Ts, double Tr,
                                    int64_t axis = -1,
                                    Normalization ni = Normalization::kInterpolation,
                                    bool fourier_space = false,
                                    double eps=1e-6);


}
#endif
