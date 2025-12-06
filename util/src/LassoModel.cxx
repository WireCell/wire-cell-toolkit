#include "WireCellUtil/LassoModel.h"

#include "WireCellUtil/Eigen.h"

using namespace Eigen;

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

std::vector<size_t> WireCell::LassoModel::Fit()
{
    std::vector<size_t> below_threshold;

    // initialize solution to zero unless user set beta already
    Eigen::VectorXd beta = _beta;
    if (0 == beta.size()) {
        beta = VectorXd::Zero(_X.cols());
    }

    // initialize active_beta to true
    int nbeta = beta.size();
    _active_beta = vector<bool>(nbeta, true);

    Eigen::VectorXd y = Gety();
    Eigen::MatrixXd X = GetX();

    // cooridate decsent
    // int N = y.size();
    VectorXd norm(nbeta);
    for (int j = 0; j < nbeta; j++) {
        norm(j) = X.col(j).squaredNorm();
        if (norm(j) < 1e-6) {
            //cerr << "warning: the " << j << "th variable is not used, please consider removing it." << endl;
            below_threshold.push_back(j);
            norm(j) = 1;
        }
    }
    double tol2 = TOL * TOL * nbeta;

    // const double inner_product_complexity = nbeta * nbeta * X.rows();

    // calculate the inner product
    Eigen::VectorXd ydX(nbeta);
    Eigen::SparseMatrix<double> XdX(nbeta, nbeta);
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
        // std::cout << " triplet method " << std::endl;
        const size_t estimated_nonzeros = nbeta * (nbeta / 2);
        XdX.reserve(Eigen::VectorXi::Constant(nbeta, nbeta/2));
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        tripletList.reserve(estimated_nonzeros);
        for (int i = 0; i < nbeta; i++) {
            ydX(i) = y.dot(X.col(i));
            for (int j = 0; j < nbeta; j++) {
                double value = X.col(i).dot(X.col(j));
                if (value != 0)
                    tripletList.push_back(T(i,j,value));
                // sum_non_zeros += value;
            }
            // if ( nbeta == 5582) {
            //     std::cout << " inner_product nbeta == 5582 " << i << std::endl;
            // }
        }
        XdX.setFromTriplets(tripletList.begin(), tripletList.end());
    }
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
                // if (it.row()!=j && beta(it.row())!=0)
                if (it.row() != j) beta(j) -= it.value() * beta(it.row());
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
