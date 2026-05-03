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



/*** high level functions ***/

DftTools::real_vector_t DftTools::convolve(const IDFT::pointer& dft,
                                 const DftTools::real_vector_t& in1,
                                 const DftTools::real_vector_t& in2)
{
    size_t size = in1.size() + in2.size() - 1;
    DftTools::complex_vector_t cin1(size,0), cin2(size,0);

    std::transform(in1.begin(), in1.end(), cin1.begin(),
                   [](float re) { return DftTools::complex_t(re,0.0); } );
    std::transform(in2.begin(), in2.end(), cin2.begin(),
                   [](float re) { return DftTools::complex_t(re,0.0); } );

    dft->fwd1d(cin1.data(), cin1.data(), size);
    dft->fwd1d(cin2.data(), cin2.data(), size);

    for (size_t ind=0; ind<size; ++ind) {
        cin1[ind] *= cin2[ind];
    }

    DftTools::real_vector_t ret(size);
    std::transform(cin1.begin(), cin1.end(), ret.begin(),
                   [](const complex_t& c) { return std::real(c); });
    return ret;
}

DftTools::real_vector_t DftTools::replace(const IDFT::pointer& dft,
                                const DftTools::real_vector_t& meas,
                                const DftTools::real_vector_t& res1,
                                const DftTools::real_vector_t& res2)
{
    size_t sizes[3] = {meas.size(), res1.size(), res2.size()};
    size_t size = sizes[0] + sizes[1] + sizes[2] - *std::min_element(sizes, sizes + 3) - 1;

    DftTools::complex_vector_t cmeas(size,0), cres1(size,0), cres2(size,0);
    std::transform(meas.begin(), meas.end(), cmeas.begin(),
                   [](float re) { return DftTools::complex_t(re,0.0); } );
    std::transform(res1.begin(), res1.end(), cres1.begin(),
                   [](float re) { return DftTools::complex_t(re,0.0); } );
    std::transform(res2.begin(), res2.end(), cres2.begin(),
                   [](float re) { return DftTools::complex_t(re,0.0); } );

    dft->fwd1d(cmeas.data(), cmeas.data(), size);
    dft->fwd1d(cres1.data(), cres1.data(), size);
    dft->fwd1d(cres2.data(), cres2.data(), size);

    for (size_t ind=0; ind<size; ++ind) {
        cmeas[ind] *= res2[ind]/res1[ind];
    }
    DftTools::real_vector_t ret(size);
    std::transform(cmeas.begin(), cmeas.end(), ret.begin(),
                   [](const complex_t& c) { return std::real(c); });

    return ret;
}

