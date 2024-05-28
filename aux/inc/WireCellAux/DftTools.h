/**
   DftTools provides support for 1D std::vector (called "vector" here)
   and 2D Eigen::Array (called "array" there typed interface to calls
   made on an IDFT.  Besides providing concrete typing and necessary
   data transforms, it includes higher level functions expressed in
   terms of IDFT primitives.  

   IDFT rules apply, in particular, forward transforms do not accrue
   any normalization and reverse transforms accrue a normalization of
   1/N where N is the number of elements transformed.

   The API uses variable names to hint at the semantic meaning of
   arguments.  We write a time/interval-space vector or array is
   "wave" or if complex as "cwave".  A frequency-space vector or array
   is written as "spec".  Where Hermitian-symmetry is anticipated we
   write "hspec" (h=Hermitian or "half" but the spectrum still spans
   the full 2*Nyquist frequency range).  No spectra in this API are
   real-valued.

   For 2D Eigen::Arrays, some functions take an "axis".  The axis
   identifies the logical array "dimension" over which the transform
   is applied.  The word "dimension" here means as in the
   "2-dimensions" (row = dimension 0, col = dimension 1) and here it
   does not mean "array size".

   For example, axis=1 means the transform is applied along the column
   dimension.  Or, in other words, conceptually it means that each row
   is transformed individually (the row elements being columns).
   Note: this is the same convention as held by numpy.fft.

   For 1D vectors, there is of course only one axis.
   
   Implementation note: the "axis" is interpreted in the "logical"
   sense of how Eigen arrays are indexed as array(irow, icol).  Ie,
   the dimension traversing rows is axis 0 and the dimension
   traversing columns is axis 1.  The internal storage order of an
   Eigen array may differ from the logical order and indeed that of
   the array template type order.  Neither is pertinent when deciding
   the value to give for "axis".

   Some function take 2D Eigen::Arrays and do not take an "axis"
   argument.  These perform transforms on both axes.

 */

#ifndef WIRECELL_AUX_DFTTOOLS
#define WIRECELL_AUX_DFTTOOLS

#include "WireCellIface/IDFT.h"
#include <vector>
#include "WireCellUtil/Eigen.h"

namespace WireCell::Aux::DftTools {

    /// Suported data types.
    using complex_t = IDFT::complex_t;
    using real_vector_t = std::vector<float>;
    using real_array_t = Eigen::ArrayXXf;
    using complex_vector_t = std::vector<complex_t>;
    using complex_array_t = Eigen::ArrayXXcf;

    // Perform forward DFT, returning a complex spectrum given a
    // complex waveform.  The 2D array version lacking an "axis"
    // performs foward DFT in both directions.
    complex_vector_t fwd(const IDFT::pointer& dft, const complex_vector_t& cwave);
    complex_array_t fwd(const IDFT::pointer& dft, const complex_array_t& cwave);
    complex_array_t fwd(const IDFT::pointer& dft, const complex_array_t& cwave, int axis);

    // Perform forward DFT, returning a complex spectrum given a real
    // waveform.  The spectrum will have Hermitian symmetry along the
    // axis of transform but only up to round-off errors accrued
    // during the transform.
    complex_vector_t fwd_r2c(const IDFT::pointer& dft, const real_vector_t& wave);
    complex_array_t fwd_r2c(const IDFT::pointer& dft, const real_array_t& wave, int axis);

    // Perform the inverse or reverse DFT, returning a complex
    // waveform given a complex spectrum.
    complex_vector_t inv(const IDFT::pointer& dft, const complex_vector_t& spec);
    complex_array_t inv(const IDFT::pointer& dft, const complex_array_t& spec);
    complex_array_t inv(const IDFT::pointer& dft, const complex_array_t& spec, int axis);

    // Perform inverse or reverse DFT, returning a real waveform given
    // a complex spectrum.  Prior to the DFT, the spectrum is forced
    // to have Hermitian symmetry and thus any input values above the
    // Nyquist frequency are ignored.
    real_vector_t inv_c2r(const IDFT::pointer& dft, const complex_vector_t& spec);
    real_array_t inv_c2r(const IDFT::pointer& dft, const complex_array_t& spec, int axis);


    /// Convolve in1 and in2 via DFT.  Returned vecgtor has size sum
    /// of sizes of in1 and in2 less one element in order to assure no
    /// periodic aliasing.  Caller need not (should not) pad either
    /// input.  Caller is free to truncate result as required.
    real_vector_t convolve(const IDFT::pointer& dft,
                           const real_vector_t& wave1,
                           const real_vector_t& wave2);


    /// Replace response res1 in meas with response res2.
    ///
    /// This will compute the FFT of all three, in frequency space will form:
    ///
    ///     meas * resp2 / resp1
    ///
    /// apply the inverse FFT and return its real part.
    ///
    /// The output vector is long enough to assure no periodic
    /// aliasing.  In general, caller should NOT pre-pad any input.
    /// Any subsequent truncation of result is up to caller.
    real_vector_t replace(const IDFT::pointer& dft,
                          const real_vector_t& meas,
                          const real_vector_t& res1,
                          const real_vector_t& res2);



    // Fixme: possible additions
    // - superposition of 2 reals for 2x speedup
}

#endif
