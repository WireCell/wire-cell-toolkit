/** LassoModel::Fit assembles the Gram X^T X by two code paths: a dense-X path (triplets +
 *  setFromTriplets) and, for a sparse response matrix such as the imaging blob-measure
 *  matrix, a direct compressed-column assembly over the support-overlapping column pairs
 *  (wcp-porting-img/wcfm/docs/08 sec 5.4).  A third path takes the response as an Eigen
 *  sparse matrix and forms the Gram by Eigen's sparse product.  All three must agree on the
 *  fit.  The direct assembly is bit-identical to the triplet build (the imaging archives are
 *  gated on that); against the sparse-product path the accumulation order differs, so the
 *  comparison here is to a tight relative tolerance.
 */

#include "WireCellUtil/LassoModel.h"
#include "WireCellUtil/Eigen.h"
#include "WireCellUtil/doctest.h"

#include <random>

using namespace WireCell;

namespace {
    // An imaging-like response: nmeas rows, nbeta columns, each column with `per_col`
    // unit entries in distinct rows (a blob in three planes' measures), true charge
    // sparse and non-negative.
    Eigen::MatrixXd make_response(int nmeas, int nbeta, int per_col, unsigned seed)
    {
        std::mt19937 rng(seed);
        std::uniform_int_distribution<int> row(0, nmeas - 1);
        Eigen::MatrixXd X = Eigen::MatrixXd::Zero(nmeas, nbeta);
        for (int c = 0; c < nbeta; ++c) {
            int placed = 0;
            while (placed < per_col) {
                const int r = row(rng);
                if (X(r, c) == 0.0) { X(r, c) = 1.0; ++placed; }
            }
        }
        return X;
    }

    Eigen::VectorXd fit_dense(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, double lambda)
    {
        LassoModel m(lambda, 100000, 1e-6, true);
        m.SetData(X, y);
        m.Fit();
        return m.Getbeta();
    }
    Eigen::VectorXd fit_sparse(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, double lambda)
    {
        LassoModel m(lambda, 100000, 1e-6, true);
        Eigen::SparseMatrix<double> Xs = X.sparseView();
        m.SetXsparse(Xs);
        m.Sety(y);
        m.Fit();
        return m.Getbeta();
    }
}

TEST_CASE("lassomodel gram assembly paths agree on a sparse imaging-like response")
{
    const int nmeas = 12, nbeta = 300;
    const auto X = make_response(nmeas, nbeta, 3, 42);
    // sparse-X branch: nnz*2 = 1800 < nbeta^2 = 90000
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> q(100.0, 1000.0);
    Eigen::VectorXd ctrue = Eigen::VectorXd::Zero(nbeta);
    for (int c = 0; c < nbeta; c += 17) ctrue(c) = q(rng);
    const Eigen::VectorXd y = X * ctrue;

    for (double lambda : {1e-3, 1e-1, 3.0}) {
        const auto bd = fit_dense(X, y, lambda);
        const auto bs = fit_sparse(X, y, lambda);
        REQUIRE(bd.size() == nbeta);
        REQUIRE(bs.size() == nbeta);
        const double scale = std::max(1.0, bd.cwiseAbs().maxCoeff());
        CHECK((bd - bs).cwiseAbs().maxCoeff() < 1e-8 * scale);
        CHECK(bd.minCoeff() >= 0.0);                       // non-negative fit
        CHECK((X * bd - y).norm() < 0.5 * y.norm());       // it did fit something
    }
}

TEST_CASE("lassomodel gram assembly paths agree on a dense response")
{
    // dense-X branch: every entry non-zero
    const int nmeas = 6, nbeta = 4;
    Eigen::MatrixXd X(nmeas, nbeta);
    X << 1, 2, 0.5, 1,
         2, 1, 1, 0.5,
         1, 1, 2, 1,
         0.5, 2, 1, 2,
         1, 0.5, 2, 1,
         2, 1, 0.5, 2;
    Eigen::VectorXd ctrue(nbeta);
    ctrue << 300, 0, 500, 0;
    const Eigen::VectorXd y = X * ctrue;
    const auto bd = fit_dense(X, y, 1e-2);
    const auto bs = fit_sparse(X, y, 1e-2);
    CHECK((bd - bs).cwiseAbs().maxCoeff() < 1e-8 * std::max(1.0, bd.cwiseAbs().maxCoeff()));
    CHECK((bd - ctrue).cwiseAbs().maxCoeff() < 1.0);       // small lambda: near the truth
}
