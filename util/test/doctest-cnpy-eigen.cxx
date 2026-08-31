/**
   Eigen defaults column-major ordering (FORTRAN).

   Cnpy defaults to row-major ordering (C).

   We thus must make a transpose before saving and after loading.

 */

#include "WireCellUtil/cnpy.h"
#include "WireCellUtil/Eigen.h"
#include "WireCellUtil/doctest.h"

using ntype = short;

// default Eigen ordering
using ArrayXXsColM = Eigen::Array<ntype, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
// non-default Eigen ordering (default for numpy)
using ArrayXXsRowM = Eigen::Array<ntype, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

const int Nrows = 3;
const int Ncols = 4;

TEST_SUITE("cnpy eigen") {

TEST_CASE("cnpy eigen") {
    const std::string name = "/tmp/wct-doctest-cnpy-eigen.npz";

    ArrayXXsColM colm = ArrayXXsColM::Zero(Nrows, Ncols);
    for (int irow = 0; irow < Nrows; ++irow) {
        for (int icol = 0; icol < Ncols; ++icol) {
            ntype val = icol + irow * Ncols;
            colm(irow, icol) = val;
        }
    }

    // Convert to Numpy's row-major from Eigen's column-major.
    ArrayXXsRowM rowm = colm;
    const ntype* data = rowm.data();
    cnpy::npz_save<ntype>(name.c_str(), "a", data, {Nrows, Ncols}, "w");

    // Load back in to a non-default row-major Eigen array.
    cnpy::NpyArray np = cnpy::npz_load(name, "a");
    REQUIRE(np.shape[0] == 3);
    REQUIRE(np.shape[1] == 4);
    auto eig = Eigen::Map<ArrayXXsRowM>(np.data<ntype>(), Nrows, Ncols);

    // Convert to Eigen default of column-major and check everything is as expected.
    ArrayXXsColM eig2 = eig;
    for (int irow = 0; irow < Nrows; ++irow) {
        for (int icol = 0; icol < Ncols; ++icol) {
            ntype val = icol + irow * Ncols;
            CHECK(eig(irow, icol) == val);
            CHECK(eig2(irow, icol) == val);
        }
    }
}

}  // TEST_SUITE("cnpy eigen")
