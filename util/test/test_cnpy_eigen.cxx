/**
   Eigen defaults column-major ordering (FORTRAN).

   Cnpy defaults to row-major ordering (C).

   We thus must make a transpose before saving and after loading.

 */

#include "WireCellUtil/cnpy.h"

#include "WireCellUtil/Eigen.h"

#include <cassert>

using ntype = short;

// default Eigen ordering
using ArrayXXsColM = Eigen::Array<ntype, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
// non-default Eigen ordering (default for numpy)
using ArrayXXsRowM = Eigen::Array<ntype, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

const int Nrows = 3;
const int Ncols = 4;

int main(int argc, char* argv[])
{
    std::string name = argv[0];
    name += ".npz";

    // ArrayXXs a = ArrayXXs::Random(Nchanc,Ntick);
    ArrayXXsColM colm = ArrayXXsColM::Zero(Nrows, Ncols);
    for (int irow = 0; irow < Nrows; ++irow) {
        for (int icol=0; icol < Ncols; ++icol) {
            ntype val = icol + irow*Ncols;
            colm(irow, icol) = val;
        }
    }

    ///
    /// Convert to Numpy's row-major from Eigen's column-major.
    ///
    ArrayXXsRowM rowm = colm;
    const ntype* data = rowm.data();
    cnpy::npz_save<ntype>(name.c_str(), "a", data, {Nrows, Ncols}, "w");


    /* With python, confirm::

      python -c'import numpy; a=numpy.load("./build/util/test_cnpy_eigen2.npz")["a"]; print(f"{a.shape}\n{a}")'
      (3, 4)
      [[ 0  1  2  3]
       [ 4  5  6  7]
       [ 8  9 10 11]]
    */

    //
    // Load back in to a non-default row-major Eigen array.
    //
    cnpy::NpyArray np = cnpy::npz_load(name, "a");
    assert(np.shape[0] == 3);
    assert(np.shape[1] == 4);
    auto eig = Eigen::Map<ArrayXXsRowM>(np.data<ntype>(),
                                        Nrows, Ncols);

    //
    // Convert to Eigen defaujlt of column-major and check everythign
    // is as expected.
    // 
    ArrayXXsColM eig2 = eig;
    for (int irow = 0; irow < Nrows; ++irow) {
        for (int icol=0; icol < Ncols; ++icol) {
            ntype val [[maybe_unused]] = icol + irow*Ncols;
            assert(eig(irow, icol) == val);
            assert(eig2(irow, icol) == val);
        }
    }

    return 0;
}
