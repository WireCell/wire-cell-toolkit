#include "WireCellUtil/LassoModel.h"

#include "WireCellUtil/Eigen.h"

using namespace Eigen;

#include <algorithm>
#include <iostream>
using namespace std;

/* Minimize the following problem:
 * 1/(2) * ||Y - beta * X||_2^2 + N * lambda * ||beta||_1
 */

WireCell::LassoModel::LassoModel(double lambda, int max_iter, double TOL, bool non_negtive)
  : ElasticNetModel(lambda, 1., max_iter, TOL, non_negtive)
{
    name = "Lasso";
}

WireCell::LassoModel::~LassoModel() {}

void WireCell::LassoModel::Set_init_values(std::vector<double> values)
{
    const size_t nvals = values.size();
    _beta = VectorXd::Zero(nvals);
    for (size_t i = 0; i != nvals; ++i) {
        _beta(i) = values.at(i);
    }
}

void WireCell::LassoModel::SetXsparse(const Eigen::SparseMatrix<double>& X)
{
    // Mirror ElasticNetModel::SetX (reset beta + unit lambda weights) but store the
    // response sparsely; the caller may override the weights afterwards.
    _Xsp = X;
    _use_sparse = true;
    _beta = Eigen::VectorXd::Zero(X.cols());
    SetLambdaWeight(Eigen::VectorXd::Constant(X.cols(), 1.));
}

std::vector<size_t> WireCell::LassoModel::Fit()
{
    std::vector<size_t> below_threshold;

    // initialize solution to zero unless user set beta already
    Eigen::VectorXd beta = _beta;
    if (0 == beta.size()) {
        beta = VectorXd::Zero(_use_sparse ? _Xsp.cols() : _X.cols());
    }

    // initialize active_beta to true
    int nbeta = beta.size();
    _active_beta = vector<bool>(nbeta, true);

    Eigen::VectorXd y = Gety();

    // Build the per-column norm, the X^T y vector (ydX) and the X^T X Gram (XdX) the
    // coordinate descent needs.  Two paths, identical up to FP accumulation order: the
    // historical DENSE path over _X (kept verbatim below), and a SPARSE path over _Xsp
    // that forms ydX/XdX by sparse products -- skipping the ~98% zero work on a block-
    // sparse response (QLMatching sparse_lasso) and never densifying X.
    double tol2 = TOL * TOL * nbeta;
    VectorXd norm(nbeta);
    Eigen::VectorXd ydX(nbeta);
    Eigen::SparseMatrix<double> XdX(nbeta, nbeta);

    if (_use_sparse) {
        const Eigen::SparseMatrix<double>& X = _Xsp;
        for (int j = 0; j < nbeta; j++) {
            norm(j) = X.col(j).squaredNorm();
            if (norm(j) < 1e-6) { below_threshold.push_back(j); norm(j) = 1; }
        }
        ydX = X.transpose() * y;   // X^T y  (== y.dot(X.col(i)) per row, sparsely)
        XdX = X.transpose() * X;   // X^T X  (skips the ~98% zero dot-products)
    }
    else {
    const Eigen::MatrixXd& X = GetX();   // reference, not a copy of the (GB-scale) Gram

    // cooridate decsent
    // int N = y.size();
    for (int j = 0; j < nbeta; j++) {
        norm(j) = X.col(j).squaredNorm();
        if (norm(j) < 1e-6) {
            //cerr << "warning: the " << j << "th variable is not used, please consider removing it." << endl;
            below_threshold.push_back(j);
            norm(j) = 1;
        }
    }

    // const double inner_product_complexity = nbeta * nbeta * X.rows();

    // calculate the inner product
    // std::cout << "Lasso begin inner product nbeta "  << nbeta << " X.rows " << X.rows() << " X.cols " << X.cols() << std::endl;
    // double sum_non_zeros = 0;

    // orig
    // std::cout << " orig method " << std::endl;
    // for (int i = 0; i != nbeta; i++) {
    //     ydX(i) = y.dot(X.col(i));
    //     // beta(i) = ydX(i) / norm(i); // first time result saved here ...
    //     for (int j = 0; j != nbeta; j++) {
    //         double value = X.col(i).dot(X.col(j));
    //         if (value != 0) XdX.insert(i, j) = value;
    //         sum_non_zeros += value;
    //     }
    //     if ( nbeta == 5582) {
    //         std::cout << " inner_product nbeta == 5582 " << i << std::endl;
    //     }
    // }

    // triplet
    {
        // Right-size the XdX / triplet reservations to the actual sparsity of X.
        // XdX = X^T X.  For a DENSE X (e.g. the imaging CSGraph response) it is ~dense,
        // so keep the historical nbeta-scale reservations (no growth reallocs at the
        // hd-max sizes the author tuned for).  For a block-sparse X (the QL normal-
        // equations Gram, ~2% full) the nbeta^2 triplet reserve and the nbeta/2-per-
        // column XdX reserve over-allocate by >10x (multiple GB on run 29107 evt 1015);
        // detect that case and reserve a small seed, letting both grow (amortized;
        // transient peak << nbeta^2).  Capacity only -> the assembled XdX, and the fit,
        // are byte-identical either way.  size_t: int nbeta*nbeta overflows beyond 46340.
        size_t xnnz = 0;
        std::vector<std::vector<int>> col_rows(nbeta);
        for (int c = 0; c < nbeta; ++c)
            for (int r = 0; r < X.rows(); ++r)
                if (X(r, c) != 0.0) { ++xnnz; col_rows[c].push_back(r); }
        const size_t nb2 = size_t(nbeta) * size_t(nbeta);
        const bool dense_X = (xnnz * 2 > nb2);
        // Row -> columns adjacency of the nonzero pattern, for the support-overlap
        // pair enumeration below (sparse X only).
        std::vector<std::vector<int>> row_cols;
        if (!dense_X) {
            row_cols.assign(X.rows(), {});
            for (int c = 0; c < nbeta; ++c)
                for (int r : col_rows[c]) row_cols[r].push_back(c);
        }
        const long col_reserve =
            dense_X ? (long)(nbeta / 2)
                    : std::max(8L, std::min((long)nbeta, (long)(4 * xnnz / std::max(nbeta, 1))));
        XdX.reserve(Eigen::VectorXi::Constant(nbeta, (int)col_reserve));
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        // Worst case (dense X) is nbeta diagonal + 2 mirrored per off-diagonal pair =
        // nbeta^2 triplets; for sparse X size to ~2 per XdX nonzero (mirrored).
        tripletList.reserve(dense_X ? nb2 : (size_t)(2L * nbeta * col_reserve + nbeta));
        // XdX is symmetric and dot() is element-wise commutative (identical
        // multiply/accumulate sequence under operand swap), so compute each
        // off-diagonal dot once and mirror it -- bit-identical to the full
        // i,j loop at half the cost.  No duplicate triplets are produced, so
        // setFromTriplets() yields the same matrix regardless of list order.
        // For a sparse X, only column pairs with overlapping row support can have a
        // nonzero dot: a non-overlapping pair's dense dot sums only exact-zero
        // products (+-0.0), which the legacy `value != 0` guard already dropped.  So
        // enumerate the overlapping j > i via the row->cols adjacency and run the
        // UNCHANGED dense dot on just those pairs -- identical triplets, identical
        // XdX, at a fraction of the nbeta^2 * nrows cost (the dominant hd-max
        // imaging term).  Dense X keeps the plain double loop.
        std::vector<int> jstamp(dense_X ? 0 : nbeta, -1);
        std::vector<int> jlist;
        for (int i = 0; i < nbeta; i++) {
            ydX(i) = y.dot(X.col(i));
            {
                double value = X.col(i).dot(X.col(i));
                if (value != 0)
                    tripletList.push_back(T(i,i,value));
            }
            if (dense_X) {
                for (int j = i + 1; j < nbeta; j++) {
                    double value = X.col(i).dot(X.col(j));
                    if (value != 0) {
                        tripletList.push_back(T(i,j,value));
                        tripletList.push_back(T(j,i,value));
                    }
                }
                continue;
            }
            jlist.clear();
            for (int r : col_rows[i])
                for (int j : row_cols[r])
                    if (j > i && jstamp[j] != i) { jstamp[j] = i; jlist.push_back(j); }
            std::sort(jlist.begin(), jlist.end());
            for (int j : jlist) {
                double value = X.col(i).dot(X.col(j));
                if (value != 0) {
                    tripletList.push_back(T(i,j,value));
                    tripletList.push_back(T(j,i,value));
                }
            }
        }
        XdX.setFromTriplets(tripletList.begin(), tripletList.end());
    }
    }  // end dense Gram build (the sparse path is handled above)
    // std::cout << "Lasso end inner product. sum_non_zeros " << sum_non_zeros << " XdX.nonZeros() " << XdX.nonZeros()
    //           << " XdX.rows() " << XdX.rows() << " XdX.cols() " << XdX.cols() << std::endl;

    // start interation ...
    int double_check = 0;
    for (int i = 0; i < max_iter; i++) {
        // std::cout << "Lasso begin iter " << i << std::endl;
        VectorXd betalast = beta;

        // loop through sparse matrix ...
        for (int j = 0; j != nbeta; j++) {
            if (!_active_beta[j]) {
                continue;
            }
            beta(j) = ydX(j);
            for (SparseMatrix<double>::InnerIterator it(XdX, j); it; ++it) {
                // Skip exact-zero betas: their term is (value * +-0.0) = +-0.0 and
                // "x -= +-0.0" is a bitwise no-op for every x that is not itself -0.0
                // (zero betas here are always +0.0: VectorXd::Zero init or the literal-0
                // soft-threshold return).  After the first sweeps most betas are exactly
                // 0 (LASSO support is sparse), so this skips the bulk of the
                // multiply-subtract work in the dominant coordinate-descent loop with a
                // bit-identical result.  (This guard was present in the original
                // prototype as a comment.)
                const double br = beta(it.row());
                if (it.row() != j && br != 0.0) beta(j) -= it.value() * br;
            }
            beta(j) = _soft_thresholding(beta(j) / norm(j), lambda * lambda_weight(j));

            if (fabs(beta(j)) < 1e-6) {
                _active_beta[j] = false;
            }
            // VectorXd X_j = X.col(j);
            // VectorXd beta_tmp = betalast;
            // beta_tmp(j) = 0;
            // VectorXd r_j = (y - X * beta_tmp);
            // double delta_j = X_j.dot(r_j);
            // std::cout << i << " " << j << " " << beta(j) << " " << betalast(j) << std::endl;
        }
        double_check++;
        // cout << endl;
        VectorXd diff = beta - betalast;

        // std::cout << "Lasso iter " << i << " " << diff.squaredNorm() << " " << tol2 << std::endl;

        if (diff.squaredNorm() < tol2) {
            if (double_check != 1) {
                double_check = 0;
                for (int k = 0; k < nbeta; k++) {
                    _active_beta[k] = true;
                }
            }
            else {
                // cout << "found minimum at iteration: " << i << " " << flag_initial_values << endl;
                break;
            }
        }
    }

    // save results in the model
    Setbeta(beta);
    return below_threshold;
}

double WireCell::LassoModel::chi2_l1() { return 2 * lambda * Getbeta().lpNorm<1>() * Gety().size(); }
