/**
   Eigen defaults column-major ordering (FORTRAN).

   Cnpy defaults to row-major ordering (C).

   We thus must make a transpose before saving and after loading.

   This is like doctest-cnpy-eigen.cxx but uses the numpy helpers.

 */

#include "WireCellUtil/NumpyHelper.h"
#include "WireCellUtil/Array.h"
#include "WireCellUtil/doctest.h"

using ntype = short;

using ArrayType = Eigen::Array<ntype, Eigen::Dynamic, Eigen::Dynamic>;

const int Nrows = 3;
const int Ncols = 4;

TEST_SUITE("cnpy eigen2") {

TEST_CASE("save and load 2d array") {
    const std::string fname = "/tmp/wct-doctest-cnpy-eigen2.npz";
    const std::string aname = "a";

    ArrayType arr = ArrayType::Zero(Nrows, Ncols);
    for (int irow = 0; irow < Nrows; ++irow) {
        for (int icol = 0; icol < Ncols; ++icol) {
            ntype val = icol + irow * Ncols;
            arr(irow, icol) = val;
        }
    }

    WireCell::Numpy::save2d(arr, aname, fname);

    /* With python, confirm:

      python -c'import numpy; a=numpy.load("/tmp/wct-doctest-cnpy-eigen2.npz")["a"]; print(f"{a.shape}\n{a}")'
      (3, 4)
      [[ 0  1  2  3]
       [ 4  5  6  7]
       [ 8  9 10 11]]
    */

    ArrayType arr2;
    WireCell::Numpy::load2d(arr2, aname, fname);

    for (int irow = 0; irow < Nrows; ++irow) {
        for (int icol = 0; icol < Ncols; ++icol) {
            ntype val = icol + irow * Ncols;
            CHECK(arr2(irow, icol) == val);
        }
    }
}

}  // TEST_SUITE("cnpy eigen2")
