/** This provides a single interface to Lasso and ElasticNet models.

    Given m = R*s, measurement vector m, response matrix R it attempts
    to solve for vector s with chi2 with a bias term linear in s.  The
    vector s is called also "source" in this code.

    References:

  - Regularization Paths for Generalized Linear Models via Coordinate
    Descent, Jerome Friedman, Trevor Hastie and Rob Tibshirani,
    Journal of Statistical Software, Volume 33, Issue 1, (2010).

  - Pathwise Corrdinate Optimization, Jerome Friedman, Trevor Hastie
    and Rob Tibshirani, https://arxiv.org/abs/0708.1485

 */
#ifndef WIRECELL_RESS_HEADER_SEEN
#define WIRECELL_RESS_HEADER_SEEN

#include "WireCellUtil/Eigen.h"

namespace WireCell {

    namespace Ress {

        typedef Eigen::VectorXd vector_t;
        typedef Eigen::MatrixXd matrix_t;

        enum Model {
            unknown = 0,
            lasso,  // Lasso model
            elnet   // elastic net
        };

        struct Params {
            Model model = elnet;
            double lambda = 1.0;
            int max_iter = 100000;
            double tolerance = 1e-3;
            bool non_negative = true;
	  bool set_init = false;
            double alpha = 1.0;
        };

        // Solve m = R*s for s, return s.
        vector_t solve(
            // matrix R (passed by const ref: the Gram can be GBs on busy QL events,
            // and solve() only forwards it to the model -- no need to copy it here)
            const matrix_t& response,
            // measured vector m
            vector_t measured,
            // params
            const Params& params = Params(),
            // optional initial source s
            vector_t source = Eigen::VectorXd(),
            // optional initial measurement weights
            vector_t weights = Eigen::VectorXd());

        // Sparse-response overload of the above: same LASSO solve, but the response
        // matrix R is sparse (a block-sparse Gram).  The model then forms X^T X / X^T y
        // by sparse products instead of the dense Gram-of-Gram, a large time+memory win
        // on busy QLMatching events.  Only the lasso model is supported here.
        vector_t solve(
            const Eigen::SparseMatrix<double>& response,
            vector_t measured,
            const Params& params = Params(),
            vector_t source = Eigen::VectorXd(),
            vector_t weights = Eigen::VectorXd());

        // These function provide values derived from a solution
        // ("solved"/"source") and the input response and measured
        // vectors.

        // Return a prediction for a measure vector.
        inline vector_t predict(matrix_t response, vector_t source)
        {
            return response * source;
        }

        // Return the unbiased part of the chi2.
        inline double chi2_base(vector_t measured, vector_t predicted)
        {
            return (measured - predicted).squaredNorm();
        }

        // Return the linear bias term
        inline double chi2_l1(vector_t measured, vector_t solved, double lambda=1.0)
        {
            return 2 * lambda * solved.lpNorm<1>() * measured.size();
        }

        // Return the average residual.
        inline double mean_residual(vector_t measured, vector_t predicted)
        {
            return (measured - predicted).norm() / measured.size();
        }


    }  // namespace Ress

}  // namespace WireCell

#endif /* WIRECELL_RESS */
