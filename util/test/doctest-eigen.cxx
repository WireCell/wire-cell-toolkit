#include "WireCellUtil/doctest.h"
#include "WireCellUtil/ExecMon.h"
#include "WireCellUtil/EigenFFT.h"
#include "WireCellUtil/Logging.h"

#include <cmath>
#include <vector>
#include <memory>
#include <sstream>

using namespace std;
using namespace Eigen;

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
    { std::ostringstream oss; oss << arr; spdlog::debug("filter.arr = {}", oss.str()); }
    auto ret = arr + 1;
    { std::ostringstream oss; oss << ret; spdlog::debug("filter.ret = {}", oss.str()); }
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
    auto b = dense->block(i, j, p, q);
    spdlog::debug("{}", em("made block"));
    spdlog::debug(" {} X {}", b.rows(), b.cols());
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
        Eigen::VectorXcf fspec(ncols);
        fft.fwd(fspec, in.row(irow));
        matc.row(irow) = fspec;
    }
    em("fft: first dimension");

    for (int icol = 0; icol < ncols; ++icol) {
        Eigen::VectorXcf pspec(nrows);
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

    spdlog::debug("shared array is {} X {}", ba->rows(), ba->cols());
    auto b = return_block(em, ba, 1, 1, nbig_rows / 2, nbig_cols / 2);
    spdlog::debug("block: {} X {}", b.rows(), b.cols());
    em("got block");
}

void test_bigass(WireCell::ExecMon& em)
{
    ArrayXXf bigass = ArrayXXf::Random(nbig_rows, nbig_cols);
    em("made big array");
    auto part = select_row(bigass, 0, em);
    em("got part");
    spdlog::debug("{} X {}", part.rows(), part.cols());
    auto part2 = part * 10;
    em("used part");
    auto part3 = part2 * 0;
    em("zeroed");
    auto part4 = select_row(bigass, 0, em);
    em("select row again");
    int nzero = 0;
    for (int ind = 0; ind < part4.rows(); ++ind) {
        if (std::fabs(part4(ind) - 0.00001f) < 1e-7f) {
            ++nzero;
        }
    }
    spdlog::debug("got zero {} times out of {}", nzero, part3.rows());
    REQUIRE(0 == nzero);
    em("checked");

    auto shared_bigass = std::make_shared<ArrayXXf>(nbig_rows, nbig_cols);
    em("make_shared");
    (*shared_bigass) = bigass;
    em("copy shared");

    take_pointer(em, shared_bigass);
    em("passed as const shared pointer");

    shared_bigass = nullptr;
    em("nullified shared");
}

TEST_SUITE("eigen") {

TEST_CASE("array construction and element access") {
    std::vector<float> v{1.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 3.0};
    ArrayXf ar1 = vec2arr(v);

    ArrayXf ar2(v.size());
    ar2 << 1.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 3.0;

    ArrayXf ar3 = Map<ArrayXf>(v.data(), v.size());

    ArrayXXf table(ar1.size(), 3);
    table.col(0) = ar1;
    table.col(1) = ar2;
    table.col(2) = ar3;

    ArrayXf tmp = filter_array(ar3);
    { std::ostringstream oss; oss << tmp; spdlog::debug("Tmp col:\n{}.", oss.str()); }
    { std::ostringstream oss; oss << table; spdlog::debug("Table:\n{}.", oss.str()); }

    ArrayXf one_row = table.row(0);
    { std::ostringstream oss; oss << one_row; spdlog::debug("One row:\n{}.", oss.str()); }

    ArrayXf one_col = table.col(0);
    { std::ostringstream oss; oss << one_col; spdlog::debug("One col:\n{}.", oss.str()); }

    VectorXf v1 = ar1.matrix();

    for (size_t ind = 0; ind < v.size(); ++ind) {
        CHECK(v[ind] == ar1(ind));
        CHECK(v[ind] == ar2(ind));
        CHECK(v[ind] == ar3(ind));
        CHECK(v[ind] == v1(ind));
    }

    int n = v1.size();
    float sigma = sqrt(v1.squaredNorm() / n - ar1.mean() * ar1.mean());
    spdlog::debug("{} {} {} {} {}", ar1.size(), ar1.sum(), ar1.prod(), v1.norm(), v1.squaredNorm());
    spdlog::debug("{} +/- {}", ar1.mean(), sigma);
}

TEST_CASE("min max coefficients") {
    std::vector<float> v{1.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 3.0};
    ArrayXf ar1 = vec2arr(v);

    ArrayXf::Index maxI = -1, minI = -1;
    float minV = ar1.minCoeff(&minI);
    float maxV = ar1.maxCoeff(&maxI);

    CHECK(minI == 0);
    CHECK(maxI == 4);

    spdlog::debug("{}@{} < {}@{}", minV, minI, maxV, maxI);
}

TEST_CASE("big array operations") {
    WireCell::ExecMon em;
    em("testing bigass");
    test_bigass(em);
    em("the end");
    spdlog::debug("{}", em.summary());
}

}  // TEST_SUITE("eigen")
