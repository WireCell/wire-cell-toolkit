/**
   Wire Cell uses Eigen3 arrays for holding large block data like the
   waveforms from one plane of one readout frame.  This header
   provides a shim between Eigen3 and the rest of Wire Cell.

   There are a few important rules:

   - Eigen3 Arrays have element-wise arithmetic and Matrices have
     matrix-like arithmetic otherwise are similar.

   - They have .array() and .matrix() methods to produce one from the other.

   - Arrays are indexed by (row,col) order.

   - An Eigen3 Vector is a 1D Matrix of shape (N,1), again, (row,col).

   - A row, column or block from an array references the original
     array so can not live beyond it.

   - In Wire Cell large arrays are accessed via const shared pointer.

   Usage examples are given below.
 */

#ifndef WIRECELL_ARRAY
#define WIRECELL_ARRAY

#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Eigen.h"
#include "WireCellUtil/Spdlog.h"

#include <memory>
#include <vector>

namespace WireCell {

    namespace Array {

        typedef Eigen::ArrayXf array_xf;
        typedef Eigen::ArrayXcf array_xc;

        /// A 16 bit short integer 2D array.
        typedef Eigen::Array<short, Eigen::Dynamic, Eigen::Dynamic> array_xxs;

        /// Integer
        typedef Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic> array_xxi;

        /// Integer
        typedef Eigen::Array<long, Eigen::Dynamic, Eigen::Dynamic> array_xxl;

        /// A real, 2D array
        typedef Eigen::ArrayXXf array_xxf;

        /// A complex, 2D array
        typedef Eigen::ArrayXXcf array_xxc;

        /** downsample a 2D array along one axis by k
         *  simple average of all numbers in a bin
         *  e.g: MxN -> Mxfloor(N/k)
         *  extra rows/cols are ignored
         */
        array_xxf downsample(const array_xxf& in, const unsigned int k, const int dim = 0);

        /** upsample a 2D array along one axis by k
         *  e.g: MxN -> MxN*k
         *  all numbers in a new bin are assigned with same value
         */
        array_xxf upsample(const array_xxf& in, const unsigned int k, const int dim = 0);

        /** put a mask on in
         *  in and mask need to have same shape
         *  values > th in mask are considered pass
         */
        array_xxf mask(const array_xxf& in, const array_xxf& mask, const float th = 0.5);

        /** linear baseline subtraction along row direction
         */
        array_xxf baseline_subtraction(const array_xxf& in);

        /** Principle component analysis
           

           Perform PCA on nrows number of vectors each vector of
           length ncols holding features (typically coordinates).

           The return value has .eigenvalue() and .eigenvectors()
           methods to get the results.  Eigenvalues are sorted in
           ascending order and the eigenvectors matrix has the vectors
           as columns that correspond to order of eigenvalues.
        */
        template<typename MatrixType=Eigen::MatrixXd>
        Eigen::SelfAdjointEigenSolver<MatrixType> pca(const MatrixType& mat) {
            MatrixType centered = mat.rowwise() - mat.colwise().mean();
            // Note: many online PCA articles incorrectly omit the 1/(n-1) normalization or use 1/n.
            MatrixType cov = ( centered.adjoint() * centered ) / ((double)mat.rows() - 1.0);
            return Eigen::SelfAdjointEigenSolver<MatrixType>(cov);
        }
        

    }  // namespace Array
}  // namespace WireCell

template <> struct fmt::formatter<Eigen::Matrix<double,-1,-1>> : fmt::ostream_formatter {};
template <> struct fmt::formatter<Eigen::Matrix<double,-1, 1>> : fmt::ostream_formatter {};

#endif
