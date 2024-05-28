#include <complex>
#include <vector>
#include <iostream>

int main()
{
    // note: this compiles but doesn't do what you may expect.
    // complex<double> numbers are a 2-array of doubles: [r,i] so the
    // reinterpret_cast from complex to double gives an "interleaved"
    // array of [r0,i0,r1,i1].  Likewise from double to complex gives
    // a "complex" number of [r0 ,r1].

    using complex_t = std::complex<double>;
    using cvec = std::vector<complex_t>;
    using dvec = std::vector<double>;

    cvec c1{{0,0}, {1,1}};
    dvec d1={0,1};

    complex_t* c2 = reinterpret_cast<complex_t*>(d1.data());
    // cvec c3(c2, c2+2); OOB!!
    cvec c3(c2, c2+1);

    double* d2 = reinterpret_cast<double*>(c1.data());
    dvec d3(d2, d2+2);

    for (auto c : c3) {
        std::cerr << c << " ";
    }
    std::cerr << "\n";
    for (auto d : d3) {
        std::cerr << d << " ";
    }
    std::cerr << "\n";
    
    return 0;
}
