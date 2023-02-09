#ifndef WIRECELLCUDA_CUFFTDFT
#define WIRECELLCUDA_CUFFTDFT

#include <cuda.h>
#include <cufft.h>

#include "WireCellIface/IDFT.h"

namespace WireCell::Cuda
{
    /* Used by the `(pre/post)transformTask` functions to
     * move around pointers to the input and output arrays
     * both on the CPU and GPU */
    template<class T>
    struct memArgs
    { const T* in; T* out; };

    // The `cuFftDFT` class provieds `IDFT` based on `cuFFT`
    class cuFftDFT : public IDFT
    {
      private:
        /*These two functions prepare/clean up the `in` and `out` arrays on both the
         * the GPU and CPU */
        memArgs<cufftComplex> preTransformTasks(memArgs<complex_t> hostMem, int size) const;
        void postTransformTasks( memArgs<complex_t> hmem, memArgs<cufftComplex> dmem, int size ) const;
                
        /* The `gen*` functions handle both the forward and inverse transforms since, 
         * unlike FFTW3, cuFFT uses the same plan and execution function for both 
         * the forward and inverse transforms. */
        void gen1d(const complex_t* in, complex_t* out, int size, int dir) const;
        void gen1b(const complex_t* in, complex_t* out, int nrows, int ncols, int axis, int dir) const;
        void gen2d(const complex_t* in, complex_t* out, int nrows, int ncols, int dir) const;
                                
        // std::shared_mutex mutex;
      public:
        // See `IDFT.h` for more information about these.
        cuFftDFT() = default;
        virtual ~cuFftDFT() = default;

        virtual void fwd1d(const complex_t* in, complex_t* out, int size) const;
        virtual void inv1d(const complex_t* in, complex_t* out, int size) const;

        virtual void fwd1b(const complex_t* in, complex_t* out, int nrows, int ncols, int axis) const;
        virtual void inv1b(const complex_t* in, complex_t* out, int nrows, int ncols, int axis) const;

        virtual void fwd2d(const complex_t* in, complex_t* out, int nrows, int ncols) const;
        virtual void inv2d(const complex_t* in, complex_t* out, int nrows, int ncols) const;

    };
}

#endif
