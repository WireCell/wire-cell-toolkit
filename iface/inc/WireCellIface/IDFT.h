#ifndef WIRECELL_IDFT
#define WIRECELL_IDFT

#include "WireCellUtil/IComponent.h"
#include <complex>

namespace WireCell {

    /** 
        Interface to perform discrete Fourier transforms on arrays of
        signal precision, complex floating point values.

        There are 6 DFT methods which are formed as the outer product
        of two lists:

        - fwd, inv
        - 1d, 1b, 2d

        The "fwd" methods provide forward transform, no normalization.
        The "inv" methods provide reverse/inverse transform normalized
        by 1/size.

        The 1d transforms take rank=1 / 1D arrays and perform a single
        transform.

        The 2d transforms take rank=2 / 2D arrays and perform nrows of
        transforms along rows and ncols of transforms along columns.
        The order over which each dimension is transformed is
        implementation-defined (and imaterial).

        The 1b transforms take rank=1 / 2D arrays and perform
        transforms along a single dimension as determined by the value
        of the "axis" parameter.  An axis=1 means to perform nrows
        transforms along rows.  Note, this is the same convention
        followed by numpy.fft functions.

        There is also a special rank=0 DFT on rank=2 arrays which is
        more commonly known as a "matrix transpose".

        Requirements on implementations:

        - Forward transforms SHALL NOT apply normalization.

        - Reverse transforms SHALL apply 1/n normalization.

        - The arrays SHALL be assumed to follow C-ordering aka
          row-major storage order.
      
        - Transform methods SHALL allow the input and output array
          pointers to be identical.

        - The IDFT interface provides 1b methods implemented in terms
          of 1d calls and a implementation MAY override these (for
          example, if implementation can exploit batch optimization).

        - Implementation SHALL allow safe concurrent calls to methods
          by different threads of execution.

        Requirement on callers.

        - Input and output arrays SHALL be pre-allocated and be sized
          at least as large as indicated by accompanying size arguments.

        - Input and output arrays MUST either be non-overlapping in
          memory or MUST be identical.

        Notes: 

        - All arrays are of type single precision complex floating
          point.  Functions and methods to easily convert between the
          two exist.

        - Eigen arrays are column-wise by default and so their
          arr.data() method can not directly supply input to this
          interface.  Likewise, use of arr.transpose().data() may run
          afowl of Eigen's IsRowMajor optimization flag.  Copy your
          default array in a Eigen::RowMajor array first or use IDFT
          via Aux::DftTools functions.
    */
    class IDFT  : public IComponent<IDFT> {
      public:
        virtual ~IDFT();

        /// The type for the signal in each bin.
        using scalar_t = float;

        /// The type for the spectrum in each bin.
        using complex_t = std::complex<scalar_t>;

        // 1d 

        virtual 
        void fwd1d(const complex_t* in, complex_t* out, int size) const = 0;

        virtual
        void inv1d(const complex_t* in, complex_t* out, int size) const = 0;

        // Real 1d transforms.  fwd_r2c_1d transforms a real waveform
        // of `size` samples to its full-size complex spectrum
        // (Hermitian by construction).  inv_c2r_1d transforms a
        // complex spectrum to a real waveform with 1/size
        // normalization; the spectrum is assumed Hermitian and values
        // above the Nyquist frequency are ignored.  Default
        // implementations widen to the complex transforms; an
        // implementation MAY override with true real-optimized
        // transforms, whose results agree with the complex path only
        // up to round-off.
        virtual
        void fwd_r2c_1d(const scalar_t* in, complex_t* out, int size) const;

        virtual
        void inv_c2r_1d(const complex_t* in, scalar_t* out, int size) const;

        // 1b

        virtual
        void fwd1b(const complex_t* in, complex_t* out,
                   int nrows, int ncols, int axis) const;

        virtual
        void inv1b(const complex_t* in, complex_t* out,
                   int nrows, int ncols, int axis) const;

        // Real 1b transforms, batched analogues of fwd_r2c_1d /
        // inv_c2r_1d following the fwd1b/inv1b array and axis
        // conventions.  fwd_r2c_1b fills "out" with the full-size
        // spectrum, Hermitian along the transform axis by
        // construction.  inv_c2r_1b assumes the spectrum is Hermitian
        // along the transform axis (values above the Nyquist
        // frequency are ignored) and applies 1/n normalization with n
        // the length along that axis.  Default implementations widen
        // to the complex 1b transforms; an implementation MAY
        // override with true real-optimized batched transforms, whose
        // results agree with the complex path only up to round-off.
        virtual
        void fwd_r2c_1b(const scalar_t* in, complex_t* out,
                        int nrows, int ncols, int axis) const;

        virtual
        void inv_c2r_1b(const complex_t* in, scalar_t* out,
                        int nrows, int ncols, int axis) const;

        // 2d

        virtual 
        void fwd2d(const complex_t* in, complex_t* out,
                   int nrows, int ncols) const = 0;
        virtual 
        void inv2d(const complex_t* in, complex_t* out,
                   int nrows, int ncols) const = 0;


        // Fill "out" with the transpose of "in", may be in-place.
        // The nrows/ncols refers to the shape of the input.
        virtual
        void transpose(const scalar_t* in, scalar_t* out,
                       int nrows, int ncols) const;
        virtual
        void transpose(const complex_t* in, complex_t* out,
                       int nrows, int ncols) const;
        
     };
}


#endif
