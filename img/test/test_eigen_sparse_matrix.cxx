
#include "WireCellUtil/Eigen.h"


#include <iostream>


using dense_mat_t = Eigen::MatrixXf;
using sparse_mat_t = Eigen::SparseMatrix<float>;
using triplet_t = Eigen::Triplet<float>;
using triplet_vec_t = std::vector<triplet_t>;

sparse_mat_t mask (const sparse_mat_t& sm, const float th = 0) {
    sparse_mat_t smm = sm;
    for (int k=0; k<smm.outerSize(); ++k) {
        for (sparse_mat_t::InnerIterator it(smm,k); it; ++it)
        {
            if (it.value() > th) {
                it.valueRef() = 1;
            }
        }
    }
    return smm;
}

int main () {

    dense_mat_t dm1(3,3);
    dm1 << 11, 12, 13,
           0, 0, 0,
           0, 0, 0;
    dense_mat_t dm2(3,3);
    dm2 << 11, 12, 0,
           21, 0, 0,
           31, 0, 0;

    auto sm1 = dm1.sparseView();
    std::cout << sm1.sum() << std::endl;
    auto sm2 = dm2.sparseView();
    std::cout << sm2.sum() << std::endl;

    sparse_mat_t sm1m = mask(sm1);
    auto sm = sm2.cwiseProduct(sm1m);

    std::cout << sm.sum() << std::endl;


    triplet_vec_t coeff;
    coeff.push_back({1, 2, 42});
    coeff.push_back({1, 2, 43});
    auto c2sm = sparse_mat_t(2,3);
    c2sm.setFromTriplets(coeff.begin(), coeff.end());
    std::cout << c2sm << std::endl;
    return 0;
}