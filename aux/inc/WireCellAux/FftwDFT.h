#ifndef WIRECELLAUX_FFTWDFT
#define WIRECELLAUX_FFTWDFT

#include "WireCellIface/IDFT.h"

namespace WireCell::Aux {

    /** 
        The FftwDFT component provides IDFT based on FFTW3.

        All instances share a common thread-safe plan cache.  There is
        no benefit to using more than one instance in a process.

        See IDFT.h for important comments.
    */
    class FftwDFT : public IDFT {
      public:
        
        FftwDFT();
        virtual ~FftwDFT();

        // 1d 

        virtual 
        void fwd1d(const complex_t* in, complex_t* out,
                   int size) const;

        virtual
        void inv1d(const complex_t* in, complex_t* out,
                   int size) const;

        // Native FFTW real transforms (r2c/c2r plans).
        virtual
        void fwd_r2c_1d(const scalar_t* in, complex_t* out,
                        int size) const;

        virtual
        void inv_c2r_1d(const complex_t* in, scalar_t* out,
                        int size) const;

        virtual
        void fwd1b(const complex_t* in, complex_t* out,
                   int nrows, int ncols, int axis) const;

        virtual 
        void inv1b(const complex_t* in, complex_t* out,
                   int nrows, int ncols, int axis) const;

        virtual 
        void fwd2d(const complex_t* in, complex_t* out,
                   int nrows, int ncols) const;
        virtual 
        void inv2d(const complex_t* in, complex_t* out,
                   int nrows, int ncols) const;

        virtual
        void transpose(const scalar_t* in, scalar_t* out,
                       int nrows, int ncols) const;
        virtual
        void transpose(const complex_t* in, complex_t* out,
                       int nrows, int ncols) const;

    };
}

#endif
