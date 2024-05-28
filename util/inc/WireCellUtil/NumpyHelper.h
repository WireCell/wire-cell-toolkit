/** NumpyHelper

    Some helper functions for using with Numpy files.
    
    Eigen default 2D array order is column-major.

    Numpy is row-major.

    This uses cnpy.  For alternative that uses custard/pigenc, see
    sio/FrameFile{Sink,Source}.

 */

#ifndef WIRECELL_UTIL_NUMPYHELPER
#define WIRECELL_UTIL_NUMPYHELPER

#include "WireCellUtil/cnpy.h"
#include "WireCellUtil/Persist.h"

#include "WireCellUtil/Eigen.h"

#include <string>
#include <vector>

namespace WireCell::Numpy {


    template <typename ARRAY>
    void save2d(ARRAY& array, std::string aname, std::string fname, std::string mode = "w") {
        using ROWM = Eigen::Array<typename ARRAY::Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        using Scalar = typename ARRAY::Scalar;
        ROWM rowm = array;     // assure we have row major
        const Scalar* data = rowm.data();
        std::vector<size_t> shape(2);
        shape[0] = static_cast<size_t>(array.rows());
        shape[1] = static_cast<size_t>(array.cols());
        WireCell::Persist::assuredir(fname);
        cnpy::npz_save<Scalar>(fname.c_str(), aname, data, 
                               shape, mode);
    }
    
    template <typename ARRAY>
    void load2d(ARRAY& array, std::string aname, std::string fname) {
        using Scalar = typename ARRAY::Scalar;

        cnpy::NpyArray np = cnpy::npz_load(fname, aname);

        if (np.fortran_order) {
            using COLM = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
            COLM temp = Eigen::Map<COLM>(np.data<Scalar>(),
                                         np.shape[0], np.shape[1]);
            array = temp;
        }
        else{
            using ROWM = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
            ROWM temp = Eigen::Map<ROWM>(np.data<Scalar>(),
                                         np.shape[0], np.shape[1]);
            array = temp;
        }
    }

}
#endif
