#include "WireCellUtil/Interpolate.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#pragma GCC diagnostic pop

#include <iomanip>
#include <iostream>

// https://www.boost.org/doc/libs/1_75_0/libs/math/doc/html/math_toolkit/cardinal_cubic_b.html

using boost::math::interpolators::cardinal_cubic_b_spline;

using namespace WireCell;

template<typename Terp>
void do_terp(Terp& terp)
{
    assert(terp(0) == 0);
    assert(terp(9) == 9);
    assert(terp(-1) == 0);
    assert(terp(10) == 9);

    std::vector<typename Terp::ytype> terped;
    const size_t num = 100 + 1;
    terp(std::back_inserter(terped), num, 0, 0.1);
    // for (size_t ind=0; ind<num; ++ind) {
    //     std::cerr << "[" << ind << "] = " << terped[ind] << "\n";
    // }
    assert(terped.size() == num);
    assert(terped[0] == 0);
    assert(terped.back() == 9);
    std::cerr << "terped[90]-9=" << terped[90]-9 << "\n";
    assert(std::abs(terped[90]-9) < 2e-6);
}

template<typename X, typename Y=X>
void test_linterp()
{
    std::vector<Y> ydata{0,1,2,3,4,5,6,7,8,9};
    linterp<X,Y> terp(ydata.begin(), ydata.end(), 0, 1);
    do_terp(terp);
}
template<typename X, typename Y=X>
void test_irrterp()
{
    std::map<X,Y> data{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{6,6},{7,7},{8,8},{9,9}};
    irrterp<X,Y> terp(data.begin(), data.end());
    do_terp(terp);
}

void test_boost()
{
    std::vector<double> f{0.01, -0.02, 0.3, 0.8, 1.9, -8.78, -22.6};
    const double xstep = 0.01;
    const double x0 = xstep;
    cardinal_cubic_b_spline<double> spline(f.begin(), f.end(), x0, xstep);

    linterp<double> lin(f.begin(), f.end(), x0, xstep);

    for (double x = 0; x < x0 + xstep * 10; x += 0.1 * xstep) {
        std::cout << std::setprecision(3) << std::fixed << "x=" << x << "\tlin(x)=" << lin(x)
                  << "\tspline(x)=" << spline(x) << "\n";
    }
}

int main()
{
    std::cerr << "testing irrterp:\n";
    test_irrterp<float>();
    test_irrterp<double>();
    test_irrterp<float,  float>();
    test_irrterp<double, float>();
    test_irrterp<float,  double>();
    test_irrterp<double, double>();
    std::cerr << "testing linterp:\n";
    test_linterp<float>();
    test_linterp<double>();
    test_linterp<float,  float>();
    test_linterp<double, float>();
    test_linterp<float,  double>();
    test_linterp<double, double>();
    return 0;
}
