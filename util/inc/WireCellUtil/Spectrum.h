#ifndef WIRECELL_SPECTRUM
#define WIRECELL_SPECTRUM

#include "WireCellUtil/Interpolate.h"
#include "WireCellUtil/Exceptions.h"
#include <complex>
#include <algorithm>

#include "Eigen/Core"           // for version test

namespace WireCell::Spectrum {

    /** This namespace holds types and functions for primitive
     * operation on frequency space spectra.  See WireCell::Waveform
     * for similar support for time series.  And higher level
     * operations are found in WireCell::Aux::NoiseTools */


    /** The Rayleigh PDF.  Real world noise spectra are approximately
     * this shape.  The sigma is the Rayleigh parameter (mode).
     */
    template<typename Real>
    Real rayleigh(Real freq, Real sigma) {
        const Real s2 = sigma*sigma;
        return (freq/s2)*exp(-0.5*freq*freq/s2);
    }
    template<typename FreqIt, typename PdfIt>
    void rayleigh(FreqIt ibeg, FreqIt iend, PdfIt out,
                  typename FreqIt::value_type sigma)
    {
        std::transform(ibeg, iend, out,
                       [=](typename FreqIt::value_type freq) {
                           return rayleigh(freq, sigma);});
    }

    /** Given a number u distributed uniformaly in (0,1) return a
     * number distributed according to a Rayleigh distribution of
     * parameter sigma. */
    template<typename Real>
    Real rayleigh_uniform(const Real& u, const Real& sigma)
    {
        return sigma * sqrt(-2*log(u));
    }

    template<typename T>
    std::complex<T> conjugate(const std::complex<T>& a) {
        return std::conj(a);
    }
    template<typename T> T conjugate(const T& a) { return a; }
    
    /// Transform input spectrum into one that exhibits Hermitian
    /// mirror symmetry around the zero and Nyquist frequencies (and
    /// also periodic symmetry).  Two function forms are given, first
    /// is in-place and second places result in a second sequence.
    template<typename It>
    void hermitian_mirror(It beg, It end) {

        // Zero frequency bin.
        *beg = std::real(*beg);

        // if len is even: mid is Nyquist bin.
        //           else: mid is last "postive" freq bin.
        const auto len = std::distance(beg, end);
        const It mid = beg + len/2;

        if (len % 2) {          // odd, no Nyquist bin
            std::reverse_copy(beg+1, mid+1, mid+1);
        }
        else {                  // even, mid is Nyquist bin
            *mid = std::real(*mid);
            std::reverse_copy(beg+1, mid, mid+1);
        }

        // "negative" frequencies must be conjugate
        std::transform(mid+1, end, mid+1,
                       [](typename It::value_type a) {
                           return conjugate(a);
                       });
    }
    template<typename InputIt, typename OutputIt>
    OutputIt hermitian_mirror(InputIt ibeg, InputIt iend, OutputIt obeg) {
        const auto len = std::distance(ibeg, iend);
        // Get index of bin above Nyquist frequency (regardless of len
        // being even or odd).
        InputIt ihalf = ibeg + len/2 + 1;
        std::copy(ibeg, ihalf, obeg);
        hermitian_mirror(obeg, obeg + len);
        return obeg + len;
    }

    // In-place Hermitian mirror along the given axis.
    template<typename Array>
    void hermitian_mirror_inplace(Array& arr, int axis) {
        if (axis == 0) {        // along rows means column-wise
#if EIGEN_VERSION_AT_LEAST(3,4,0)
            for (auto col : arr.colwise()) {
                hermitian_mirror(col.begin(),
                                 col.end(),
                                 col.begin());
            }
#else
            for (long icol=0; icol<arr.cols(); ++icol) {
                std::vector<typename Array::Scalar> tmp(arr.rows(), 0);
                for (long ind=0; ind<arr.rows(); ++ind) {
                    tmp[ind] = arr(ind, icol);
                }
                hermitian_mirror(tmp.begin(), tmp.end(), tmp.begin());
                for (long ind=0; ind<arr.rows(); ++ind) {
                    arr(ind, icol) = tmp[ind];
                }
            }                
#endif
            return;
        }
        if (axis == 1) {        // along columns mean row-wise
#if EIGEN_VERSION_AT_LEAST(3,4,0)
            for (auto row : arr.rowwise()) {
                hermitian_mirror(row.begin(),
                                 row.end(),
                                 row.begin());
            }
#else
            for (long irow=0; irow<arr.rows(); ++irow) {
                std::vector<typename Array::Scalar> tmp(arr.cols(), 0);
                for (long ind=0; ind<arr.cols(); ++ind) {
                    tmp[ind] = arr(irow, ind);
                }
                hermitian_mirror(tmp.begin(), tmp.end(), tmp.begin());
                for (long ind=0; ind<arr.cols(); ++ind) {
                    arr(irow, ind) = tmp[ind];
                }
            }                
#endif
            return;
        }
        THROW(ValueError() << errmsg{"illegal axis"});
    }
    // Return new array.
    template<typename Array>
    Array hermitian_mirror(const Array& arr, int axis) {
        Array ret = arr;
        hermitian_mirror_inplace(ret, axis);
        return ret;
    }



    /** Return interpolation of a spectrum.
        
        Result has size "newsize" and is normalized by
        sqrt(newsize/oldsize).  A wave generated from the interpolated
        spectrum will have same RMS as one generated from the original
        (though, will of course be longer so more energy).

        The sample period and the Nyquist frequency of the returned
        spectrum should be considered to be held constant.

        Note, interpolation by newsize < oldsize will not apply
        aliasing.  See the alias() function to produce a correct
        downsampling in time.
     */
    template<typename InputIt, typename OutputIt>
    void interp(InputIt ibeg, InputIt iend, OutputIt obeg, OutputIt oend)
    {
        linterpolate(ibeg, iend, obeg, oend);
        const double oldn = std::abs(std::distance(ibeg,iend));
        const double newn = std::abs(std::distance(obeg,oend));
        const double norm = sqrt(newn/oldn);
        std::transform(obeg, oend, obeg,
                       [=](const typename OutputIt::value_type& y) {
                           return y * norm; });
    }


    /** Return an extrapolation of a spectrum

        This is equivalent to an interpolation in time.

        This extrapolation adds more high-frequency bins to the
        spectrum with an amplitude that of the constant value.  If
        constant is negative, the amplitude at the Nyquist frequency
        is used to fill in the new bins.

        The total waveform time (N*T) is unchanged while sampling and
        Nyquist frequencies increase and thus the period decreases
        each by the ratio of new and old sizes.

        The new spectrum is normalized by sqrt(newsize/oldsize) in
        order that energy be conserved.  This will cause a decrease in
        RMS if constant is zero.
     */
    template<typename InputIt, typename OutputIt>
    void extrap(InputIt ibeg, InputIt iend, OutputIt obeg, OutputIt oend,
                typename OutputIt::value_type constant = 0)
    {
        const size_t ilen = std::distance(ibeg, iend);
        const size_t olen = std::distance(obeg, oend);
        if (ilen == olen) {
            std::copy(ibeg, iend, obeg); // trivial
            return;
        }
        if (ilen > olen) {
            THROW(ValueError() << errmsg{"extrap requires larger output than input"});
        }

        // Nyquist bin or last "positive" freq bin.
        const size_t half = ilen/2;
        if (constant < 0) {
            constant = *(ibeg+half);
        }
        // If we have Nyquist bin, we want to avoid copying it
        size_t lo = half;
        const size_t hi = half+1;
        size_t extra = 1;
        if (ilen % 2) {
            // If odd, we copy everything.
            ++lo;
            extra = 0;          // no Nyquist bin to remove
        }
        // Copy low half below Nyquist frequency
        std::copy(ibeg, ibeg+lo, obeg);

        // Copy high half above Nyquist frequency
        std::copy_backward(ibeg+hi, iend, oend);

        // And fill in the center
        const size_t nmid = olen-ilen + extra;
        std::fill(obeg+lo, obeg+lo+nmid, constant);

    }



    /** Return an aliasing of a spectrum

        This is equivalent to a downsample in time.

        If newsize is not a factor of oldsize, the smallest factor of
        oldsize larger than newsize will be used.  For more precise
        aliasing consider to interpolate prior to aliasing.

        Keeps Rayleigh frequency constant while size decreases so
        Nyquist frequency decreases and period increases.

        This normalizes to sqrt(newsize/oldsize) so that energy is
        conserved in the case that the original spectrum has had an
        anti-alias filter applied such that the aliased bins are
        contain zero amplitude.
    */
    template<typename InputIt, typename OutputIt>
    void alias(InputIt ibeg, InputIt iend, OutputIt obeg, OutputIt oend)
    {
        const size_t ilen = std::distance(ibeg, iend);
        const size_t olen = std::distance(obeg, oend);
        if (ilen == olen) {
            std::copy(ibeg, iend, obeg); // trivial
            return;
        }
        if (ilen < olen) {
            THROW(ValueError() << errmsg{"alias requires smaller output than input"});
        }
        std::fill(obeg, oend, 0);

        const size_t half = ilen/2;
        const size_t L = ceil((double)ilen/(double)olen);
        const size_t M = olen/2; // calculate only half spectrum
        for (size_t m=0; m <= M; ++m) {
            for (size_t l=0; l<L; ++l) {
                size_t oldind = m + l*M;
                if (oldind > half) {
                    // Effectively zero-extrapolated
                    break;
                }
                *(obeg+m) += *(ibeg+oldind);
            }
        }
        hermitian_mirror(obeg, oend);
        const double norm = sqrt((double)olen/(double)ilen);
        std::transform(obeg, oend, obeg, [=](const auto& a){return norm*a;});
    }


    /** Return spectrum with new size and new period.

        Note, the period must be expressed relative to input period ie,

        relperiod = newperiod / oldperiod.

        This performs an interpolation followed by an extrapolation or
        an alias depending on if newperiod is smaller or larger than
        oldperiod.
    */
    template<typename InputIt, typename OutputIt>
    void resample(InputIt ibeg, InputIt iend, OutputIt obeg, OutputIt oend, double relperiod)
    {
        //const size_t ilen = std::distance(ibeg, iend);
        const size_t olen = std::distance(obeg, oend);

        const size_t tmpsiz = ceil(olen*relperiod);
        std::vector<typename InputIt::value_type> temp(tmpsiz, 0);
        interp(ibeg, iend, temp.begin(), temp.end());

        if (relperiod > 1) {
            alias(temp.begin(), temp.end(), obeg, oend);
        }
        else {
            extrap(temp.begin(), temp.end(), obeg, oend);
        }
    }


}
#endif
