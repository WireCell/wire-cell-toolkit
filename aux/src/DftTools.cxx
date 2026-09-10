#include "WireCellAux/DftTools.h"
#include "WireCellUtil/Spectrum.h"
#include <algorithm>

#include <iostream>             // debugging


using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::Spectrum; // hermitian_mirror

/*** helpers, both vector/array types ***/


/*** vector ***/

DftTools::complex_vector_t DftTools::fwd(const IDFT::pointer& dft, const DftTools::complex_vector_t& seq)
{
    complex_vector_t ret(seq.size());
    dft->fwd1d(seq.data(), ret.data(), ret.size());
    return ret;
}

DftTools::complex_vector_t DftTools::fwd_r2c(const IDFT::pointer& dft, const DftTools::real_vector_t& vec)
{
    complex_vector_t cvec(vec.size());
    std::transform(vec.begin(), vec.end(), cvec.begin(),
                   [](float re) { return DftTools::complex_t(re,0.0); } );
    return fwd(dft, cvec);
}

DftTools::complex_vector_t DftTools::inv(const IDFT::pointer& dft, const DftTools::complex_vector_t& spec)
{
    complex_vector_t ret(spec.size());
    dft->inv1d(spec.data(), ret.data(), ret.size());
    return ret;
}

DftTools::real_vector_t DftTools::inv_c2r(const IDFT::pointer& dft, const DftTools::complex_vector_t& spec)
{
    // fixme: in future, expand IDFT to have c2r 
    complex_vector_t symspec(spec.size());
    hermitian_mirror(spec.begin(), spec.end(), symspec.begin());

    auto cvec = inv(dft, symspec);
    real_vector_t rvec(cvec.size());
    std::transform(cvec.begin(), cvec.end(), rvec.begin(),
                   [](const DftTools::complex_t& c) { return std::real(c); });
    return rvec;
}

DftTools::complex_vector_t DftTools::fwd_r2c_real(const IDFT::pointer& dft, const DftTools::real_vector_t& vec)
{
    complex_vector_t ret(vec.size());
    dft->fwd_r2c_1d(vec.data(), ret.data(), vec.size());
    return ret;
}

DftTools::real_vector_t DftTools::inv_c2r_real(const IDFT::pointer& dft, const DftTools::complex_vector_t& spec)
{
    real_vector_t ret(spec.size());
    dft->inv_c2r_1d(spec.data(), ret.data(), spec.size());
    return ret;
}

/*** array ***/

// Implementation notes for fwd()/inv():
//
// - We make an initial copy to get rid of any potential IsRowMajor
//   optimization/confusion over storage order.  This suffers a copy
//   but we need to allocate return anyways.
//
// - We then have column-wise storage order but IDFT assumes row-wise
// - so we reverse (nrows, ncols) and meaning of axis.

DftTools::complex_array_t DftTools::fwd(const IDFT::pointer& dft,
                              const DftTools::complex_array_t& arr,
                              int axis)
{
    DftTools::complex_array_t ret = arr;
    dft->fwd1b(ret.data(), ret.data(), ret.cols(), ret.rows(), !axis);
    return ret;
}

DftTools::complex_array_t DftTools::inv(const IDFT::pointer& dft,
                              const DftTools::complex_array_t& arr,
                              int axis)
{
    DftTools::complex_array_t ret = arr;
    dft->inv1b(ret.data(), ret.data(), ret.cols(), ret.rows(), !axis);
    return ret;
}

void DftTools::fwd_inplace(const IDFT::pointer& dft,
                           DftTools::complex_array_t& arr,
                           int axis)
{
    dft->fwd1b(arr.data(), arr.data(), arr.cols(), arr.rows(), !axis);
}

void DftTools::inv_inplace(const IDFT::pointer& dft,
                           DftTools::complex_array_t& arr,
                           int axis)
{
    dft->inv1b(arr.data(), arr.data(), arr.cols(), arr.rows(), !axis);
}


/*
  Big fat warning to future me: Passing by reference means the input
  array may carry the .IsRowMajor optimization for implementing
  transpose().  An extra copy would remove that complication but this
  interface tries to keep it.
 */
using ROWM = Eigen::Array<DftTools::complex_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using COLM = Eigen::Array<DftTools::complex_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

template<typename trans>
DftTools::complex_array_t doit(const DftTools::complex_array_t& arr, trans func)
{
    const DftTools::complex_t* in_data = arr.data();
    DftTools::complex_vector_t out_vec(arr.rows()*arr.cols());

    // std::cerr << "dft::doit: (" << arr.rows() << "," << arr.cols() << ") IsRowMajor:" << arr.IsRowMajor << std::endl;

    if (arr.IsRowMajor) {
        func(in_data, out_vec.data(), arr.cols(), arr.rows());
        return Eigen::Map<ROWM>(out_vec.data(), arr.rows(), arr.cols());
    }

    func(in_data, out_vec.data(), arr.rows(), arr.cols());
    return Eigen::Map<COLM>(out_vec.data(), arr.rows(), arr.cols());
}

DftTools::complex_array_t DftTools::fwd(const IDFT::pointer& dft, const DftTools::complex_array_t& arr)
{
    return doit(arr, [&](const complex_t* in_data,
                         complex_t* out_data,
                         int nrows, int ncols) {
        dft->fwd2d(in_data, out_data, nrows, ncols);
    });
}

DftTools::complex_array_t DftTools::inv(const IDFT::pointer& dft, const DftTools::complex_array_t& arr)
{
    return doit(arr, [&](const complex_t* in_data,
                         complex_t* out_data,
                         int nrows, int ncols) {
        dft->inv2d(in_data, out_data, nrows, ncols);
    });
}


DftTools::complex_array_t DftTools::fwd_r2c(const IDFT::pointer& dft, const DftTools::real_array_t& wave, int axis)
{
    complex_array_t cwave = wave.cast<complex_t>();
    fwd_inplace(dft, cwave, axis);
    return cwave;
}

DftTools::real_array_t DftTools::inv_c2r(const IDFT::pointer& dft, const DftTools::complex_array_t& spec, int axis)
{
    complex_array_t symspec = hermitian_mirror(spec, axis);
    inv_inplace(dft, symspec, axis);
    // Drops the small imaginary that is accrued due to round-off errors.
    return symspec.real();
}

// As with fwd_inplace()/inv_inplace(): Eigen arrays are column-major
// but IDFT assumes row-major, so reverse (nrows, ncols) and the
// meaning of axis.

DftTools::complex_array_t DftTools::fwd_r2c_real(const IDFT::pointer& dft, const DftTools::real_array_t& wave, int axis)
{
    complex_array_t ret(wave.rows(), wave.cols());
    dft->fwd_r2c_1b(wave.data(), ret.data(), wave.cols(), wave.rows(), !axis);
    return ret;
}

DftTools::real_array_t DftTools::inv_c2r_real(const IDFT::pointer& dft, const DftTools::complex_array_t& spec, int axis)
{
    real_array_t ret(spec.rows(), spec.cols());
    dft->inv_c2r_1b(spec.data(), ret.data(), spec.cols(), spec.rows(), !axis);
    return ret;
}



/*** high level functions ***/

// Zero-pad a real waveform to the given size and forward transform it.
static DftTools::complex_vector_t padded_spectrum(const IDFT::pointer& dft,
                                                  const DftTools::real_vector_t& wave,
                                                  size_t size)
{
    DftTools::real_vector_t padded(size, 0);
    std::copy(wave.begin(), wave.end(), padded.begin());
    return DftTools::fwd_r2c(dft, padded);
}

DftTools::real_vector_t DftTools::convolve(const IDFT::pointer& dft,
                                 const DftTools::real_vector_t& in1,
                                 const DftTools::real_vector_t& in2)
{
    const size_t size = in1.size() + in2.size() - 1;

    auto spec1 = padded_spectrum(dft, in1, size);
    const auto spec2 = padded_spectrum(dft, in2, size);

    for (size_t ind=0; ind<size; ++ind) {
        spec1[ind] *= spec2[ind];
    }

    // Inverse transform back to the time domain.  (Prior to the fix
    // for issue #531 this step was missing and the real part of the
    // spectrum was returned.)
    return DftTools::inv_c2r(dft, spec1);
}

DftTools::real_vector_t DftTools::replace(const IDFT::pointer& dft,
                                const DftTools::real_vector_t& meas,
                                const DftTools::real_vector_t& res_new,
                                const DftTools::real_vector_t& res_old)
{
    // Pad to a size large enough that the linear convolution of meas
    // with the (new/old) response ratio suffers no periodic aliasing.
    // This is the same size rule as the legacy Waveform::replace_convolve().
    const size_t sizes[3] = {meas.size(), res_new.size(), res_old.size()};
    const size_t size = sizes[0] + sizes[1] + sizes[2] - *std::min_element(sizes, sizes + 3) - 1;

    auto smeas = padded_spectrum(dft, meas, size);
    const auto snew = padded_spectrum(dft, res_new, size);
    const auto sold = padded_spectrum(dft, res_old, size);

    // Form meas * new / old in frequency space.
    //
    // Issue #531: the previous implementation divided by the untransformed,
    // unpadded time-domain response samples (reading past their ends),
    // had the new/old ratio inverted and never applied the inverse
    // transform.  See the doctest in aux/test/doctest_dfttools_replace.cxx.
    const complex_t zero{0, 0};
    for (size_t ind=0; ind<size; ++ind) {
        const complex_t den = sold[ind];
        if (den == zero) {
            // Old response has no power at this frequency so the
            // measurement should not either.  Leave the bin unscaled
            // rather than emit inf/NaN.
            continue;
        }
        smeas[ind] *= snew[ind] / den;
    }

    return DftTools::inv_c2r(dft, smeas);
}

