#include "WireCellIface/IDFT.h"

#include <vector>
#include <utility>              // std::swap since c++11

using namespace WireCell;

IDFT::~IDFT() {}

// Trivial default real transforms via the complex ones.  Override
// with native real-optimized transforms for the win (about half the
// flops and memory).

void IDFT::fwd_r2c_1d(const scalar_t* in, complex_t* out, int size) const
{
    std::vector<complex_t> cin(size);
    for (int ind = 0; ind < size; ++ind) {
        cin[ind] = complex_t(in[ind], 0);
    }
    fwd1d(cin.data(), out, size);
}

void IDFT::inv_c2r_1d(const complex_t* in, scalar_t* out, int size) const
{
    // Enforce the Hermitian assumption (values above Nyquist ignored)
    // so the result is exactly real, then take the real part.
    std::vector<complex_t> sym(in, in + size);
    for (int ind = size / 2 + 1; ind < size; ++ind) {
        sym[ind] = std::conj(sym[size - ind]);
    }
    std::vector<complex_t> cout(size);
    inv1d(sym.data(), cout.data(), size);
    for (int ind = 0; ind < size; ++ind) {
        out[ind] = cout[ind].real();
    }
}

void IDFT::fwd_r2c_1b(const scalar_t* in, complex_t* out,
                      int nrows, int ncols, int axis) const
{
    const int ntot = nrows * ncols;
    std::vector<complex_t> cin(ntot);
    for (int ind = 0; ind < ntot; ++ind) {
        cin[ind] = complex_t(in[ind], 0);
    }
    fwd1b(cin.data(), out, nrows, ncols, axis);
}

void IDFT::inv_c2r_1b(const complex_t* in, scalar_t* out,
                      int nrows, int ncols, int axis) const
{
    // Enforce the Hermitian assumption along the transform axis
    // (values above Nyquist ignored) so the result is exactly real,
    // then take the real part.
    const int ntot = nrows * ncols;
    std::vector<complex_t> sym(in, in + ntot);
    if (axis) {
        for (int irow = 0; irow < nrows; ++irow) {
            complex_t* row = sym.data() + irow * ncols;
            for (int icol = ncols / 2 + 1; icol < ncols; ++icol) {
                row[icol] = std::conj(row[ncols - icol]);
            }
        }
    }
    else {
        for (int irow = nrows / 2 + 1; irow < nrows; ++irow) {
            for (int icol = 0; icol < ncols; ++icol) {
                sym[irow * ncols + icol] = std::conj(sym[(nrows - irow) * ncols + icol]);
            }
        }
    }
    std::vector<complex_t> cout(ntot);
    inv1b(sym.data(), cout.data(), nrows, ncols, axis);
    for (int ind = 0; ind < ntot; ++ind) {
        out[ind] = cout[ind].real();
    }
}

// Trivial default "batched" implementations.  If your concrete
// implementation provides some kind of "batch optimization", such as
// with FFTW3's advanced interface or with some GPU FFT library,
// override these dumb methods for the win.

void IDFT::fwd1b(const complex_t* in, complex_t* out,
                 int nrows, int ncols, int axis) const
{
    if (axis) { 
        for (int irow=0; irow<nrows; ++irow) {
            fwd1d(in+irow*ncols, out+irow*ncols, ncols);
        }
    }
    else {
        this->transpose(in, out, nrows, ncols);
        this->fwd1b(out, out, ncols, nrows, 1);
        this->transpose(out, out, ncols, nrows);
    }
}

void IDFT::inv1b(const complex_t* in, complex_t* out,
                 int nrows, int ncols, int axis) const
{
    if (axis) { 
        for (int irow=0; irow<nrows; ++irow) {
            inv1d(in+irow*ncols, out+irow*ncols, ncols);
        }
    }
    else {
        this->transpose(in, out, nrows, ncols);
        this->inv1b(out, out, ncols, nrows, 1);
        this->transpose(out, out, ncols, nrows);
    }
}

// Trivial default transpose.  Implementations, please override if you
// can offer something faster.

template<typename ValueType>
void transpose_type(const ValueType* in, ValueType* out,
                    int nrows, int ncols) 
{
    if (in != out) {
        for (int irow=0; irow<nrows; ++irow) {
            for (int icol=0; icol<ncols; ++icol) {
                out[icol*nrows + irow] = in[irow*ncols + icol];
            }
        }
        return;
    }
    
    // inplace adapated from https://stackoverflow.com/a/9320349 which
    // comes from
    // https://en.wikipedia.org/wiki/In-place_matrix_transposition#Non-square_matrices:_Following_the_cycles

    const int n = nrows;
    const int size = nrows*ncols;
    const int mn1 = (size - 1);
    std::vector<bool> visited(size);
    ValueType* first = out;
    const ValueType* last = first + size;
    ValueType* cycle = out;
    while (++cycle != last) {
        if (visited[cycle - first])
            continue;
        int a = cycle - first;
        do  {
            a = a == mn1 ? mn1 : (n * a) % mn1;
            std::swap(*(first + a), *cycle);
            visited[a] = true;
        } while ((first + a) != cycle);
    }

}


void IDFT::transpose(const IDFT::scalar_t* in, IDFT::scalar_t* out,
                     int nrows, int ncols) const
{
    transpose_type(in, out, nrows, ncols);
}
void IDFT::transpose(const IDFT::complex_t* in, IDFT::complex_t* out,
                     int nrows, int ncols) const
{
    transpose_type(in, out, nrows, ncols);
}
