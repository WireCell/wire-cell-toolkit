// https://forum.kde.org/viewtopic.php?f=74&t=98871
// https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=736985
#include <vector>
#include <iostream>

#include "WireCellUtil/Eigen.h"
#include <unsupported/Eigen/Splines>

using namespace Eigen;

double uvalue(double x, double low, double high) { return (x - low) / (high - low); }

VectorXd uvalues(VectorXd xvals)
{
    const double low = xvals.minCoeff();
    const double high = xvals.maxCoeff();
    for (int i = 0; i < xvals.size(); ++i) {
        xvals(i) = uvalue(xvals(i), low, high);
    }
    return xvals;
}

int main(int argc, char* argv[])
{
    typedef Spline<double, 1> Spline2d;

    const VectorXd xvals = (VectorXd(5) << 1, 2, 3, 4, 6).finished();
    const VectorXd yvals = xvals.array().square();
    const Spline2d spline = SplineFitting<Spline2d>::Interpolate(yvals.transpose(), 3, uvalues(xvals).transpose());

    const double step = 0.1;
    for (double x = 1; x < 6 + 0.5 * step; x += step) {
        const double uv = uvalue(x, xvals.minCoeff(), xvals.maxCoeff());
        std::cout << "(" << x << "," << spline(uv).transpose() << "," << spline.derivatives<1>(uv) << ")\n";
    }
    std::cout << std::endl;
}
