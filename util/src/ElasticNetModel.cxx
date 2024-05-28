#include "WireCellUtil/ElasticNetModel.h"

#include "WireCellUtil/Eigen.h"
using namespace Eigen;

#include <iostream>
using namespace std;

/* Minimize the following problem:
 * 1/(2) * ||Y - beta * X||_2^2 + N * lambda * (
 *   alpha * ||beta||_1 + 0.5 * (1-alpha) * ||beta||_2^2
 * )
 * To control L1 and L2 separately, this is equivaletnt to a * L1 + b * L2,
 * where lambda = a + b and alpha = a / (a + b)
 */

WireCell::ElasticNetModel::ElasticNetModel(double lambda, double alpha, int max_iter, double TOL, bool non_negtive)
  : lambda(lambda)
  , alpha(alpha)
  , max_iter(max_iter)
  , TOL(TOL)
  , non_negtive(non_negtive)
{
    name = "Elastic net";
}

WireCell::ElasticNetModel::~ElasticNetModel() {}

std::vector<size_t> WireCell::ElasticNetModel::Fit()
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

    // use alias for easy notation
    Eigen::VectorXd y = Gety();
    Eigen::MatrixXd X = GetX();

    // cooridate decsent

    // int N = y.size();
    VectorXd norm(nbeta);
    for (int j = 0; j < nbeta; j++) {
        norm(j) = X.col(j).squaredNorm();
        if (norm(j) < 1e-6) {
            // cerr << "warning: the " << j << "th variable is not used, please consider removing it." << endl;
            below_threshold.push_back(j);
            norm(j) = 1;
        }
    }
    double tol2 = TOL * TOL * nbeta;

    int double_check = 0;
    for (int i = 0; i < max_iter; i++) {
        VectorXd betalast = beta;
        for (int j = 0; j < nbeta; j++) {
            if (!_active_beta[j]) {
                continue;
            }
            VectorXd X_j = X.col(j);
            VectorXd beta_tmp = beta;
            beta_tmp(j) = 0;
            VectorXd r_j = (y - X * beta_tmp);
            double delta_j = X_j.dot(r_j);
            //            beta(j) = _soft_thresholding(delta_j, N*lambda*alpha*lambda_weight(j)) / (1+lambda*(1-alpha))
            //            / norm(j);
            beta(j) =
                _soft_thresholding(delta_j / norm(j), lambda * alpha * lambda_weight(j)) / (1 + lambda * (1 - alpha));

            // cout << i << " " << j << " " << beta(j) << std::endl;
            if (fabs(beta(j)) < 1e-6) {
                _active_beta[j] = false;
            }
            // else { cout << beta(j) << endl;}
            // beta(j) = _soft_thresholding(delta_j, N*lambda*alpha, j) / (1+lambda*(1-alpha)) / norm(j);
            // if (j==0) cout << beta(j) << ", " << arg1 << endl;
        }
        double_check++;
        // cout << endl;
        VectorXd diff = beta - betalast;

        // std::cout << i << " " << diff.squaredNorm() << " " << tol2 << std::endl;
        if (diff.squaredNorm() < tol2) {
            if (double_check != 1) {
                double_check = 0;
                for (int k = 0; k < nbeta; k++) {
                    _active_beta[k] = true;
                }
            }
            else {
                //                cout << "found minimum at iteration: " << i << endl;
                break;
            }
        }
    }

    // save results in the model
    Setbeta(beta);
    return below_threshold;
}

double WireCell::ElasticNetModel::_soft_thresholding(double delta, double lambda_)
{
    if (delta > lambda_) {
        return delta - lambda_;
    }
    else {
        if (non_negtive) {
            return 0;
        }
        else {
            if (delta < -lambda_) {
                return delta + lambda_;
            }
            else {
                return 0;
            }
        }
    }
}
