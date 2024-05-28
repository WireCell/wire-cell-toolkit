#include "WireCellUtil/Testing.h"
#include "WireCellUtil/ExecMon.h"


#include "WireCellUtil/EigenFFT.h"

#include <cmath>
#include <vector>
#include <memory>
#include <iostream>

using namespace std;

using namespace Eigen;

// template<typename Val>
// using Array = Eigen::Array<Val, Eigen::Dynamic, 1>;
// /// Convert a std::vector to a Sequence
// template<typename Val>			      // not const correct?
// Array<Val> std2eig(/*const*/ std::vector<Val>& vec) {
//     typedef Array<Val> cont_t;
//     return Eigen::Map<cont_t>(vec.data(), vec.size());
// }

// /// Convert a Sequence to a std::vector
// template<typename Val>
// std::vector<Val> eig2std(/*const*/ Array<Val>& vec) {
//     return std::vector<Val>(vec.data(), vec.data() + vec.size());
// }

Eigen::ArrayXf vec2arr(const std::vector<float>& v)
{
    Eigen::ArrayXf ret(v.size());
    for (size_t ind = 0; ind < v.size(); ++ind) {
        ret(ind) = v[ind];
    }
    return ret;
}

// manual says care is needed for passing by value.  Unneeded
// temporaries can be created or worse.
Eigen::ArrayXf filter_array(const Eigen::ArrayXf& arr)
{
    cerr << "filter.arr = " << arr << endl;
    auto ret = arr + 1;
    cerr << "filter.ret = " << ret << endl;
    return ret;
}

Eigen::ArrayXf select_row(const Eigen::ArrayXXf& arr, int ind, WireCell::ExecMon& em)
{
    auto tmp = arr.row(ind);  // no copy
    em("after assignment to auto type");
    Eigen::ArrayXf ret = arr.row(ind);  // this does a copy
    em("after assignment to explicit type");
    return tmp;  // this does a copy
}

template <typename Derived>
using shared_dense = std::shared_ptr<Eigen::DenseBase<Derived> >;
template <typename Derived>
using const_shared_dense = std::shared_ptr<const Eigen::DenseBase<Derived> >;

typedef Eigen::ArrayXXf array_xxf;
typedef shared_dense<array_xxf> shared_array_xxf;
typedef const_shared_dense<array_xxf> const_shared_array_xxf;

typedef Eigen::ArrayXXcf array_xxc;
typedef shared_dense<array_xxc> shared_array_xxc;

template <typename Derived>
Eigen::Block<const Derived> return_block(WireCell::ExecMon& em, const_shared_dense<Derived> dense, int i, int j, int p,
                                         int q)
{
    // Eigen::Block<const Derived> b = dense->block(i,j,p,q);
    auto b = dense->block(i, j, p, q);
    cerr << em("made block") << endl;
    cerr << " " << b.rows() << " X " << b.cols() << endl;
    return b;
}

void do_fft(WireCell::ExecMon& em, const array_xxf& arr)
{
    const int nrows = arr.rows();
    const int ncols = arr.cols();

    em("fft: start");

    Eigen::MatrixXf in = arr.matrix();
    em("fft: convert to matrix");
    Eigen::MatrixXcf matc(nrows, ncols);
    em("fft: made temp complex matrix");

    Eigen::FFT<float> fft;
    em("fft: made fft object");

    for (int irow = 0; irow < nrows; ++irow) {
        Eigen::VectorXcf fspec(ncols);  // frequency spectrum
        fft.fwd(fspec, in.row(irow));
        matc.row(irow) = fspec;
    }
    em("fft: first dimension");

    for (int icol = 0; icol < ncols; ++icol) {
        Eigen::VectorXcf pspec(nrows);  // periodicity spectrum
        fft.fwd(pspec, matc.col(icol));
        matc.col(icol) = pspec;
    }
    em("fft: second dimension");

    shared_array_xxc ret = std::make_shared<array_xxc>(nrows, ncols);
    em("fft: make shared for return");
    (*ret) = matc;
    em("fft: set shared for return");
    ret = nullptr;
    em("fft: nullify return");
}

const int nbig_rows = 3000;
const int nbig_cols = 10000;

void take_pointer(WireCell::ExecMon& em, const_shared_array_xxf ba)
{
    do_fft(em, *ba);
    em("fft: done");

    cerr << "shared array is " << ba->rows() << " X " << ba->cols() << endl;
    auto b = return_block(em, ba, 1, 1, nbig_rows / 2, nbig_cols / 2);
    cerr << "block: " << b.rows() << " X " << b.cols() << endl;
    em("got block");
}

void test_bigass(WireCell::ExecMon& em)
{
    // not really *that* big....
    ArrayXXf bigass = ArrayXXf::Random(nbig_rows, nbig_cols);
    em("made big array");
    auto part = select_row(bigass, 0, em);
    em("got part");
    cerr << part.rows() << " X " << part.cols() << "\n";
    auto part2 = part * 10;
    em("used part");
    auto part3 = part2 * 0;
    em("zeroed");
    auto part4 = select_row(bigass, 0, em);
    em("select row again");
    int nzero = 0;
    for (int ind = 0; ind < part4.rows(); ++ind) {
        if (part4(ind) == 0.00001) {
            ++nzero;
        }
    }
    cerr << "got zero " << nzero << " times out of " << part3.rows() << endl;
    Assert(0 == nzero);
    em("checked");

    auto shared_bigass = std::make_shared<ArrayXXf>(nbig_rows, nbig_cols);
    em("make_shared");
    (*shared_bigass) = bigass;
    em("copy shared");

    take_pointer(em, shared_bigass);  // cast to const
    em("passed as const shared pointer");

    shared_bigass = nullptr;
    em("nullified shared");
}

int main()
{
    WireCell::ExecMon em;

    std::vector<float> v{1.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 3.0};
    ArrayXf ar1 = vec2arr(v);  // copy okay

    /// You must specify storage size at construction
    // ArrayXf ar2;
    ArrayXf ar2(v.size());

    /// this doesn't work:
    // ar2 << v;
    /// but literal comma list does:
    ar2 << 1.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 3.0;

    /// or, map the data
    ArrayXf ar3 = Map<ArrayXf>(v.data(), v.size());

    ArrayXXf table(ar1.size(), 3);
    table.col(0) = ar1;
    table.col(1) = ar2;
    table.col(2) = ar3;

    ArrayXf tmp = filter_array(ar3);
    cerr << "Tmp col:\n" << tmp << ".\n";

    cerr << "Table:\n" << table << ".\n";

    ArrayXf one_row = table.row(0);
    cerr << "One row:\n" << one_row << ".\n";

    ArrayXf one_col = table.col(0);
    cerr << "One col:\n" << one_col << ".\n";

    VectorXf v1 = ar1.matrix();

    for (size_t ind = 0; ind < v.size(); ++ind) {
        Assert(v[ind] == ar1(ind));
        Assert(v[ind] == ar2(ind));
        Assert(v[ind] == ar3(ind));
        Assert(v[ind] == v1(ind));
    }

    cerr << ar1.size() << " " << ar1.sum() << " " << ar1.prod() << " " << v1.norm() << " " << v1.squaredNorm() << endl;
    int n = v1.size();
    float sigma = sqrt(v1.squaredNorm() / n - ar1.mean() * ar1.mean());
    cerr << ar1.mean() << " +/- " << sigma << endl;

    ArrayXf::Index maxI = -1, minI = -1;
    float minV = ar1.minCoeff(&minI);
    float maxV = ar1.maxCoeff(&maxI);

    Assert(minI == 0);
    Assert(maxI == 4);

    cerr << minV << "@" << minI << " < " << maxV << "@" << maxI << endl;

    em("testing bigass");

    test_bigass(em);

    em("the end");
    cerr << em.summary() << endl;

    return 0;
}
