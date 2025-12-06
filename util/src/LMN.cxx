#include "WireCellUtil/LMN.h"
#include "WireCellUtil/Exceptions.h"
#include <cmath>                // fmod
#include <algorithm>            // std::copy

using Eigen::seq;
using Eigen::seqN;
using Eigen::lastN;
// This is gives deprecated warning.
// using Eigen::placeholders::all;
using Eigen::all;

using namespace WireCell;
using namespace WireCell::Array;

double LMN::gcd(double a, double b, double eps)
{
    if (a < eps) {
        return b;
    }
    return LMN::gcd(fmod(b, a), a);
}

size_t LMN::rational(double Ts, double Tr, double eps)
{
    const double dT = std::abs(Ts - Tr);
    const double gcd = LMN::gcd(Tr, dT, eps);
    if (gcd == 0) return 0;
    const double n = dT/gcd;
    const size_t rn = round(n);
    const double err1 = std::abs(n - rn);
    if (err1 > eps) {
        raise<ValueError>("gcd error one too big %f > %f, rn=%f, n=%f", err1, eps, rn, n);
        return 0;
    }

    double Ns = rn * Tr / dT;
    size_t rNs = round(Ns);
    const double err2 = std::abs(rNs - Ns);
    if (err2 > eps) {
        raise<ValueError>("gcd error two too big %f > %f", err1, eps);
        return 0;
    }

    return rNs;    
}

Array::array_xxf LMN::resize(const Array::array_xxf& in, size_t Nr,
                             size_t axis)
{
    size_t Nr_rows, Nr_cols;
    size_t Ns = 0;
    if (axis == 0) {            // resampled along rows
        Ns = in.rows();
        Nr_rows = Nr;
        Nr_cols = in.cols();
    }
    else {                      // resample along columns
        Ns = in.cols();
        Nr_rows = in.rows();
        Nr_cols = Nr;
    }

    array_xxf rs = array_xxf::Zero(Nr_rows, Nr_cols);

    size_t N_min = std::min(Nr,Ns);
    auto safe = seqN(0,N_min);

    if (axis == 0) {
        rs(safe, all) = in(safe, all);
        rs(safe, all) = in(safe, all);
    }
    else {
        rs(all, safe) = in(all, safe);
        rs(all, safe) = in(all, safe);
    }
    return rs;
}

void LMN::fill_constant(std::vector<float>::iterator begin,
                        std::vector<float>::iterator end,
                        float value)
{
    while (begin != end) {
        *begin = value;
        ++begin;
    }
}

void LMN::fill_linear(std::vector<float>::iterator begin,
                      std::vector<float>::iterator end,
                      float first, float last)
{
    size_t N = std::distance(begin, end);
    for (size_t ind=0; ind<N; ++ind) {
        *(begin+ind) = first + ind*((last-first)/N);
    }
}


std::vector<float> LMN::resize(const std::vector<float>& in, size_t Nr)
{
    size_t Ns = in.size();
    if (Ns == Nr) return in;
    size_t N_min = std::min(Nr,Ns);
    float pad = 0;


    std::vector<float> rs(Nr, pad);
    std::copy(in.begin(), in.begin()+N_min, rs.begin());

    return rs;
}

Array::array_xxc LMN::resample(const Array::array_xxc& in, size_t Nr, size_t axis)
{
    size_t Nr_rows, Nr_cols;

    size_t Ns = 0;
    if (axis == 0) {            // resampled along rows
        Ns = in.rows();
        Nr_rows = Nr;
        Nr_cols = in.cols();
    }
    else {                      // resample along columns
        Ns = in.cols();
        Nr_rows = in.rows();
        Nr_cols = Nr;
    }

    size_t N_half = 0;
    if (Nr > Ns) {              // upsample
        N_half = LMN::nhalf(Ns);
    }
    else {                      // downsample
        N_half = LMN::nhalf(Nr);
    }

    array_xxc rs = array_xxc::Zero(Nr_rows, Nr_cols);

    // The regions to copy into the resampled array.
    auto pos_half = seqN(0,N_half+1);
    auto neg_half = lastN(N_half);

    if (axis == 0) {
        rs(pos_half, all) = in(pos_half, all);
        rs(neg_half, all) = in(neg_half, all);
    }
    else {
        rs(all, pos_half) = in(all, pos_half);
        rs(all, neg_half) = in(all, neg_half);
    }

    // FIXME: deal with Nyquist bin.
    return rs;
}

    

std::vector<std::complex<float>>
LMN::resample(const std::vector<std::complex<float>>& in, size_t Nr)
{
    size_t Ns = in.size();

    size_t N_half = 0;
    if (Nr > Ns) {              // upsample
        N_half = LMN::nhalf(Ns);
    }
    else {                      // downsample
        N_half = LMN::nhalf(Nr);
    }

    std::vector<std::complex<float>> rs(Nr);

    std::copy(in.begin(), in.begin()+N_half+1, rs.begin());
    std::copy(in.rbegin(), in.rbegin()+N_half, rs.rbegin());

    // FIXME: deal with Nyquist bin.
    return rs;
}

    

