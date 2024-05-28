#include "WireCellUtil/Spectrum.h"


#include "WireCellUtil/Eigen.h"

#include <numeric>
#include <cassert>
#include <iostream>

using namespace WireCell::Spectrum;

void test_rayleigh()
{
    const double sigma = 10.0;
    std::vector<double> freqs = {0,1,2,3,4,5,6,7};
    std::vector<double> rpdf(freqs.size());
    rayleigh(freqs.begin(), freqs.end(), rpdf.begin(), sigma);
    for (size_t ind=0; ind<freqs.size(); ++ind) {
        assert(rpdf[ind] == rayleigh(freqs[ind], sigma));
    }
}


template<typename T>
void dump(const std::vector<T>& vec, const std::string& name="")
{
    std::cerr << name;
    for (const auto& a : vec) {
        std::cerr << " " << a;;
    }
    std::cerr << "\n";
}

template<typename T>
std::vector<T> odd() {
    return std::vector<T>{-1,-2,-3,-4,-5,-6,-7};
}
template<>
std::vector<std::complex<float>> odd<std::complex<float>>()
{
    return std::vector<std::complex<float>>{{-1,0},{-2,2},{-3,3},{-4,4},{-5,5},{-6,6},{-7,7}};
}
template<typename T>
std::vector<T> even() {
    return std::vector<T>{-1,-2,-3,-4,-5,-6,-7,-8};
}
template<>
std::vector<std::complex<float>> even<std::complex<float>>()
{
    return std::vector<std::complex<float>>{{-1,0},{-2,2},{-3,3},{-4,-4},{-5,5},{-6,6},{-7,7},{-8,8}};
}

template<typename T>
void test_hermitian()
{
    {
        std::vector<T> fo = odd<T>();
        std::vector<T> fog(fo.size());
        hermitian_mirror(fo.begin(), fo.end(), fog.begin());

        // odd
        dump(fo, "fo");
        dump(fog, "fog");
        assert(fog[0] == std::abs(fo[0]));
        assert(fog[1] == fo[1]);
        assert(fog[2] == fo[2]);
        assert(fog[3] == fo[3]);
        // -->nyquist frequency<--
        assert(fog[4] == conjugate(fo[3]));
        assert(fog[5] == conjugate(fo[2]));
        assert(fog[6] == conjugate(fo[1]));

    }
    {
        std::vector<T> fe = even<T>();
        std::vector<T> feg(fe.size());
        hermitian_mirror(fe.begin(), fe.end(), feg.begin());

        // even
        dump(fe, "fe");
        dump(feg, "feg");
        assert(feg[0] == std::abs(fe[0]));
        assert(feg[1] == fe[1]);
        assert(feg[2] == fe[2]);
        assert(feg[3] == fe[3]);
        assert(feg[4] == std::abs(fe[4])); // Nyquist bin
        assert(feg[5] == conjugate(fe[3]));
        assert(feg[6] == conjugate(fe[2]));
        assert(feg[7] == conjugate(fe[1]));
    }
}
void test_hermitian_eigen_rows()
{
    const int axis = 1;

    std::cerr << "axis="<<axis<<" along columns aka row-wise\n";
    auto e = even<float>();
    size_t nrows=3, ncols=e.size();
    Eigen::ArrayXXf arr(nrows, ncols);
    for (size_t irow=0; irow<nrows; ++irow) {
        for (size_t icol=0; icol<ncols; ++icol) {
            arr(irow, icol) = e[icol];
        }
    }
    std::cerr << arr << "\n";
    auto harr = hermitian_mirror(arr, axis);

    std::cerr << arr << "\n";

}
void test_hermitian_eigen_cols()
{
    const int axis = 0;

    std::cerr << "axis="<<axis<<" along rows aka column-wise\n";
    auto e = even<float>();
    size_t ncols=3, nrows=e.size();
    Eigen::ArrayXXf arr(nrows, ncols);
    for (size_t irow=0; irow<nrows; ++irow) {
        for (size_t icol=0; icol<ncols; ++icol) {
            arr(irow, icol) = e[irow];
        }
    }
    std::cerr << arr << "\n";
    auto harr = hermitian_mirror(arr, axis);
    std::cerr << harr << "\n";
}

template<typename V>
typename V::value_type energy(const V& amp)
{
    typename V::value_type sqr = 0;
    
    for (const auto& a : amp) {
        sqr += a*a;
    }
    return sqr/amp.size();

    // if (amp.empty()) { return start; };
    // typename V::value_type sqr = 
    //     std::accumulate(amp.begin(), amp.end(), start,
    //                     [](auto tot, auto a) {
    //                         return tot += a*conjugate(a);
    //                     });
    // return sqr/amp.size();
}

template<typename Real>
void test_interp()
{
    std::vector<Real> freqs(100), amp(100), amp2(121);
    std::iota(freqs.begin(), freqs.end(), 0);
    const Real sigma = 10.0;
    rayleigh(freqs.begin(), freqs.end(), amp.begin(), sigma);
    hermitian_mirror(amp.begin(), amp.end());
    interp(amp.begin(), amp.end(), amp2.begin(), amp2.end());
    auto rms1 = sqrt(energy(amp)/amp.size());
    auto rms2 = sqrt(energy(amp2)/amp2.size());
    auto err = std::abs(rms1-rms2)/(rms1+rms2);
    std::cerr << "rms1=" << rms1 << ", rms2=" << rms2 << ", relative error = " << err << "\n";
    assert(err < 1e-3);         // linear interpolation is not that great....
}


template<typename Real>
void test_extrap(size_t siz, size_t siz2)
{
    std::cerr << "extrap " << siz << " ->  " << siz2 << "\n";
    std::vector<Real> amp(siz), amp2(siz2, 42.0);
    for (size_t ind=0; ind<siz; ++ind) {
        amp[ind] = ind-(Real)siz;
    }
    dump(amp, "orig");
    hermitian_mirror(amp.begin(), amp.end());
    dump(amp, "herm");
    extrap(amp.begin(), amp.end(), amp2.begin(), amp2.end());
    dump(amp2, "extrap");
}

template<typename Real>
void test_alias(size_t siz, size_t siz2)
{
    std::cerr << "alias " << siz << " ->  " << siz2 << "\n";
    std::vector<Real> amp(siz), amp2(siz2, 42.0);
    for (size_t ind=0; ind<siz; ++ind) {
        amp[ind] = ind-(Real)siz;
    }
    dump(amp, "orig");
    hermitian_mirror(amp.begin(), amp.end());
    dump(amp, "herm");
    alias(amp.begin(), amp.end(), amp2.begin(), amp2.end());
    dump(amp2, "alias");

}

template<typename Real>
void test_resample(size_t siz, size_t siz2, double period)
{
    std::cerr << "period " << siz << " ->  " << siz2 << " period=" << period << "\n";
    std::vector<Real> amp(siz), amp2(siz2, 42.0);
    for (size_t ind=0; ind<siz; ++ind) {
        amp[ind] = ind-(Real)siz;
    }
    dump(amp, "orig");
    hermitian_mirror(amp.begin(), amp.end());
    dump(amp, "herm");
    resample(amp.begin(), amp.end(), amp2.begin(), amp2.end(), period);
    dump(amp2, "resample");

}

int main()
{
    // test_rayleigh();
    // test_hermitian<float>();
    // test_hermitian<double>();
    // test_hermitian<std::complex<float>>();
    test_hermitian_eigen_rows();
    test_hermitian_eigen_cols();
    return 0;
    // test_interp<double>();
    test_interp<float>();
    test_extrap<float>(10, 13);
    test_extrap<float>(10, 14);
    test_extrap<float>(9, 13);
    test_extrap<float>(9, 14);

    test_alias<float>(13, 10);
    test_alias<float>(14, 10);
    test_alias<float>(13, 9);
    test_alias<float>(14, 9);

    std::vector<double> periods {0.7, 1.3};
    for (const auto& period : periods) {

        test_resample<float>(10, 13, period);
        test_resample<float>(10, 14, period);
        test_resample<float>(9, 13, period);
        test_resample<float>(9, 14, period);

        test_resample<float>(13, 10, period);
        test_resample<float>(14, 10, period);
        test_resample<float>(13, 9, period);
        test_resample<float>(14, 9, period);

    }    
    return 0;


}
