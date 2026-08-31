#ifndef WIRECELLUTIL_LASSOMODEL_H
#define WIRECELLUTIL_LASSOMODEL_H

#include "ElasticNetModel.h"

namespace WireCell {

    class LassoModel : public ElasticNetModel {
       public:
        LassoModel(double lambda = 1., int max_iter = 100000, double TOL = 1e-3, bool non_negtive = true);
        virtual ~LassoModel();

        // Perform the fit and return indices of variables below threshold.
        // These can be ignored or the fit may be retried with these variables removed.
        virtual std::vector<size_t> Fit();
        void Set_init_values(std::vector<double> values);

        // Sparse-response path (used by the Ress::solve sparse overload).  When set,
        // Fit() builds the X^T X Gram and X^T y it needs from this sparse response
        // matrix instead of the dense _X -- skipping the dense Gram-of-Gram on a
        // block-sparse X (e.g. the QLMatching sparse_lasso normal equations).  The
        // dense path (imaging, and any caller of the dense solve) is untouched.
        // See match/docs/qlmatching-perf-evt1015-pdhd.md.
        void SetXsparse(const Eigen::SparseMatrix<double>& X);

        double chi2_l1();

       private:
        Eigen::SparseMatrix<double> _Xsp;   // response matrix when _use_sparse
        bool _use_sparse{false};
    };

}  // namespace WireCell

#endif
