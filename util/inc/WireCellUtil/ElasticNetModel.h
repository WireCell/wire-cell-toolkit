#ifndef WIRECELLUTIL_ELASTICNETMODEL_H
#define WIRECELLUTIL_ELASTICNETMODEL_H

#include "LinearModel.h"
#include <vector>
#include "WireCellUtil/Eigen.h"

namespace WireCell {

    class ElasticNetModel : public LinearModel {
       public:
        ElasticNetModel(double lambda = 1., double alpha = 1., int max_iter = 100000, double TOL = 1e-3,
                        bool non_negtive = true);
        virtual ~ElasticNetModel();

        double lambda;  // regularization parameter
        double alpha;   // L1 ratio (L2 ratio = 1 - alpha)
        int max_iter;   // maximum iteration
        double TOL;
        bool non_negtive;
        Eigen::VectorXd lambda_weight;  // each beta can have a weight on its regularization, default weight is 1;

        void SetLambdaWeight(Eigen::VectorXd w) { lambda_weight = w; }
        void SetLambdaWeight(int i, double weight) { lambda_weight(i) = weight; }
        void SetX(Eigen::MatrixXd X)
        {
            LinearModel::SetX(X);
            SetLambdaWeight(Eigen::VectorXd::Zero(X.cols()) + Eigen::VectorXd::Constant(X.cols(), 1.));
        }
        // Perform the fit and return indices of variables below threshold.
        // These can be ignored or the fit may be retried with these variables removed.
        virtual std::vector<size_t> Fit();

       protected:
        double _soft_thresholding(double x, double lambda_);
        std::vector<bool> _active_beta;
    };

}  // namespace WireCell

#endif
