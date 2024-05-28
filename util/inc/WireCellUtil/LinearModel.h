#ifndef WIRECELLUTIL_LINEARMODEL_H
#define WIRECELLUTIL_LINEARMODEL_H

#include "WireCellUtil/Eigen.h"
#include <string>
#include <vector>

namespace WireCell {

    class LinearModel {
       public:
        LinearModel();
        virtual ~LinearModel();

        Eigen::VectorXd& Gety() { return _y; }
        Eigen::MatrixXd& GetX() { return _X; }
        Eigen::VectorXd& Getbeta() { return _beta; }

        virtual void SetData(Eigen::MatrixXd X, Eigen::VectorXd y)
        {
            SetX(X);
            Sety(y);
        }
        virtual void Sety(Eigen::VectorXd y) { _y = y; }
        virtual void SetX(Eigen::MatrixXd X)
        {
            _X = X;
            _beta = Eigen::VectorXd::Zero(X.cols());
        }
        virtual void Setbeta(Eigen::VectorXd beta) { _beta = beta; }

        virtual std::vector<size_t> Fit(){ return std::vector<size_t>(); };
        Eigen::VectorXd Predict();
        double chi2_base();
        double MeanResidual();
        std::string name;

       protected:
        // Fit: y = X * beta
        // convention: lowercase: vector, uppercase matrix.
        Eigen::VectorXd _y;
        Eigen::MatrixXd _X;
        Eigen::VectorXd _beta;
    };

}  // namespace WireCell

#endif
