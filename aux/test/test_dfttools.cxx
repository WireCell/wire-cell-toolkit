#include "aux_test_dft_helpers.h"

#include "WireCellAux/DftTools.h"
#include "WireCellAux/FftwDFT.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Spectrum.h"

#include <iostream>
#include <memory>

#include <complex>
using namespace std::literals;

using namespace WireCell;
using namespace WireCell::Spectrum;
using namespace WireCell::Aux::Test;
using namespace WireCell::Aux::DftTools;

using real_t = float;
using RV = std::vector<real_t>;
using complex_t = std::complex<real_t>;
using CV = std::vector<complex_t>;

void test_1d_impulse(IDFT::pointer dft, int size = 64)
{
    RV rimp(size, 0);
    rimp[0] = 1.0;

    auto cimp = fwd(dft, Waveform::complex(rimp));
    assert_flat_value(cimp.data(), cimp.size());

    RV rimp2 = Waveform::real(inv(dft, cimp));
    assert_impulse_at_index(rimp2.data(), rimp2.size());
}

using FA = Eigen::Array<real_t, Eigen::Dynamic, Eigen::Dynamic>;
using CA = Eigen::Array<complex_t, Eigen::Dynamic, Eigen::Dynamic>;
using FARM = Eigen::Array<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using CARM = Eigen::Array<complex_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

void test_2d_impulse(IDFT::pointer dft, int nrows=16, int ncols=8)
{
    const size_t size = nrows*ncols;
    FA r = FA::Zero(nrows, ncols);
    r(0,0) = 1.0;
    dump("r", r);
    assert_impulse_at_index(r.data(), size);

    CA rc = r.cast<complex_t>();
    dump("rc", rc);
    assert_impulse_at_index(rc.data(), size);

    CA c = fwd(dft, rc);
    dump("c", c);
    assert_flat_value(c.data(), size);

    FA r2 = inv(dft, c).real();
    dump("r2", r2);
    assert_impulse_at_index(r2.data(), size);
}

void test_2d_eigen_transpose(IDFT::pointer dft)
{
    const int nrows=16;
    const int ncols=8;

    // where the impulse lives (off axis)
    const int imp_row = 1;
    const int imp_col = 10;

    FA r = FA::Zero(nrows, ncols); // shape:(16,8)
    dump("r", r);

    // do not remove the auto in this next line
    auto rt = r.transpose();    // shape:(8,16)
    dump("rt", rt);
    rt(imp_row, imp_col) = 1.0;

    auto c = fwd(dft, rt.cast<complex_t>());
    dump("c", c);

    auto r2 = inv(dft, c).real();
    dump("r2",r2);

    // transpose access
    const int nrowst = r2.rows();
    const int ncolst = r2.cols();

    for (int irow=0; irow<nrowst; ++irow) {
        for (int icol=0; icol<ncolst; ++icol) {
            float val = rt(irow, icol);
            float val2 = r2(irow, icol); // access with transposed indices
            std::cerr << "(" << irow << ","<< icol << "):" << val << " ? " << val2 << "\n";
            if (irow==imp_row and icol==imp_col) {
                assert(std::abs(val-1.0) < 1e-6);
                continue;
            }
            assert(std::abs(val) < 1e-6);    
        }
        std::cerr << "\n";
    }
}

void test_1b(IDFT::pointer dft, int axis, int nrows=8, int ncols=4)
{
    FA r = FA::Zero(nrows, ncols);
    r(0,0) = 1.0;
    dump("impulse", r);
    CA c = fwd(dft, r.cast<complex_t>(), axis);

    dump("spectra", c);
    std::cerr << c << std::endl;

    if (axis) {                 // transform along rows
        CA ct = c.transpose();      // convert to along columns (native Eigen storage order)
        c = ct;
        std::swap(nrows, ncols);
        dump("transpose", c);
        std::cerr << c << std::endl;
    }

    // first column has flat abs value of 1.0.
    assert_flat_value(c.data(), nrows, complex_t(1,0)); 
    // rest should be flat, zero value
    assert_flat_value(c.data()+nrows, nrows*ncols - nrows, complex_t(0,0)); 

}

void test_hs1d()
{
    std::cerr << "1d hermitian symmetry\n";
    //                            0      1      2      3      4     5
    const complex_vector_t even{{1,11},{2,22},{3,33},{4,44},{5,55},{6,66}};
    const complex_vector_t odd(even.begin(), even.end()-1);

    complex_vector_t evens(even.size());
    complex_vector_t odds(odd.size());
    hermitian_mirror(even.begin(), even.end(), evens.begin());
    hermitian_mirror(odd.begin(), odd.end(), odds.begin());

    // The Nyquist bin must be real.
    { size_t ind=3;
        std::cerr << " even["<< ind <<"]=" << even[ind]
                  << " evens["<<ind<<"]=" << evens[ind] << "\n";
        assert(std::abs(evens[ind].imag()) < 1e-6);
    }

    assert(std::real(even[0]) == evens[0]);
    assert(std::real(odd[0]) == odds[0]);

    for (size_t ind=1; ind<3; ++ind) {
        assert(even[ind] == evens[ind]);
        assert(odd[ind] == odds[ind]);
    }
    for (size_t ind=1; ind<3; ++ind) {
        size_t even_other = even.size() - ind;
        size_t odd_other = odd.size() - ind;
        std::cerr << " even["<< ind <<"]=" << even[ind]
                  << " evens["<<even_other<<"]=" << evens[even_other]
                  << " odd[" << ind << "]=" << odd[ind]
                  << " odds[" << odd_other << "]=" << odds[odd_other]
                  << "\n";
        assert(even[ind] == std::conj(evens[even_other]));
        assert(odd[ind] == std::conj(odds[odd_other]));
    }
}

void test_hs2d_axis1()
{
    using c = std::complex<float>;

    std::cerr << "2d hermitian symmetry\n";
    CA even(1,6), odd(1,5);
    even << c(1,11), c(2,22), c(3,33), c(4,44), c(5,55), c(6,66);
    odd  << c(1,11), c(2,22), c(3,33), c(4,44), c(5,55);
    dump("even", even);
    std::cerr << even << std::endl;
    dump("odd", odd);
    std::cerr << odd << std::endl;

    auto evens = hermitian_mirror(even, 1);
    auto odds = hermitian_mirror(odd, 1);

    // The Nyquist bin must be real.
    { size_t ind=3;
        std::cerr << " even["<< ind <<"]=" << even(0,ind)
                  << " evens["<<ind<<"]=" << evens(0,ind) << "\n";
        assert(std::abs(evens(0,ind).imag()) < 1e-6);
    }

    dump("evens", evens);
    std::cerr << evens << std::endl;
    dump("odds", odds);
    std::cerr << odds << std::endl;

    assert(std::real(even(0,0)) == evens(0,0));
    assert(std::real(odd(0,0)) == odds(0,0));

    for (size_t ind=1; ind<3; ++ind) {
        assert(even(0,ind) == evens(0,ind));
        assert(odd(0,ind) == odds(0,ind));
    }
    for (size_t ind=1; ind<3; ++ind) {
        size_t even_other = even.cols() - ind;
        size_t odd_other = odd.cols() - ind;
        std::cerr << " even["<< ind <<"]=" << even(0,ind)
                  << " evens["<<even_other<<"]=" << evens(0,even_other)
                  << " odd[" << ind << "]=" << odd(0,ind)
                  << " odds[" << odd_other << "]=" << odds(0,odd_other)
                  << "\n";
        assert(even(0,ind) == std::conj(evens(0,even_other)));
        assert(odd(0,ind) == std::conj(odds(0,odd_other)));
    }
}

                     
void test_hs2d_axis0()
{
    using c = std::complex<float>;

    std::cerr << "2d hermitian symmetry\n";
    CA even(6,1), odd(5,1);
    even << c(1,11), c(2,22), c(3,33), c(4,44), c(5,55), c(6,66);
    odd  << c(1,11), c(2,22), c(3,33), c(4,44), c(5,55);
    dump("even", even);
    std::cerr << even << std::endl;
    dump("odd", odd);
    std::cerr << odd << std::endl;

    auto evens = hermitian_mirror(even, 0);
    auto odds = hermitian_mirror(odd, 0);

    // The Nyquist bin must be real.
    { size_t ind=3;
        std::cerr << " even["<< ind <<"]=" << even(ind,0)
                  << " evens["<<ind<<"]=" << evens(ind,0) << "\n";
        assert(std::abs(evens(ind,0).imag()) < 1e-6);
    }

    dump("evens", evens);
    std::cerr << evens << std::endl;
    dump("odds", odds);
    std::cerr << odds << std::endl;

    assert(std::real(even(0,0)) == evens(0,0));
    assert(std::real(odd(0,0)) == odds(0,0));

    for (size_t ind=1; ind<3; ++ind) {
        assert(even(ind,0) == evens(ind,0));
        assert(odd(ind,0) == odds(ind,0));
    }
    for (size_t ind=1; ind<3; ++ind) {
        size_t even_other = even.rows() - ind;
        size_t odd_other = odd.rows() - ind;
        std::cerr << " even["<< ind <<"]=" << even(ind,0)
                  << " evens["<<even_other<<"]=" << evens(even_other,0)
                  << " odd[" << ind << "]=" << odd(ind,0)
                  << " odds[" << odd_other << "]=" << odds(odd_other,0)
                  << "\n";
        assert(even(ind,0) == std::conj(evens(even_other,0)));
        assert(odd(ind,0) == std::conj(odds(odd_other,0)));
    }
}

void test_1d_c2r_impulse(IDFT::pointer dft, size_t size = 64)
{
    std::cerr << "1d c2r impulse size=" << size << "\n";

    RV rimp(size, 0);
    rimp[0] = 1.0;

    auto spec = fwd_r2c(dft, rimp);
    spec[3*size/4] = 42;
    auto wave = inv_c2r(dft, spec);
    for (size_t ind=0; ind<size; ++ind) {
        // std::cerr << ind << " " << rimp[ind] << " " << wave[ind] << " " << spec[ind] << "\n";
        if (ind) {
            assert(wave[ind] == 0.0);
        }
        else {
            assert(wave[ind] == 1.0);
        }
    }        
}
void test_2d_c2r_impulse(IDFT::pointer dft, int axis, int nrows=8, int ncols=4)
{
    std::cerr << "2d c2r impulse axis=" << axis << " shape=(" << nrows << ", " << ncols <<")\n";

    const size_t size = nrows*ncols;
    FA r = FA::Zero(nrows, ncols);
    r(0,0) = 1.0;
    dump("rimp2d", r);
    std::cerr << r << std::endl;
    assert_impulse_at_index(r.data(), size);

    auto spec = fwd_r2c(dft, r, axis);
    dump("spec2d", spec);
    std::cerr << spec << std::endl;

    auto wave = inv_c2r(dft, spec, axis);
    dump("wave2d", wave);
    std::cerr << wave << std::endl;
    assert_impulse_at_index(wave.data(), size);    
}

int main(int argc, char* argv[])
{
    test_hs1d();
    test_hs2d_axis1();
    test_hs2d_axis0();
    return 0;

    DftArgs args;
    int rc = make_dft_args(args, argc, argv);
    if (rc) { return rc; }
    auto idft = make_dft(args.tn, args.pi, args.cfg);

    test_1d_c2r_impulse(idft);
    test_2d_c2r_impulse(idft, 0);
    test_2d_c2r_impulse(idft, 1);

    test_1d_impulse(idft);
    test_2d_impulse(idft);
    test_2d_eigen_transpose(idft);
    test_1b(idft, 0);
    test_1b(idft, 1);

    return 0;
}
