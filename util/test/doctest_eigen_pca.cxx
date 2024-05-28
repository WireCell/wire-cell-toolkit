#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Array.h"
#include "WireCellUtil/doctest.h"

using namespace Eigen;

void dump(std::string msg,
          const MatrixXd& mat, const MatrixXd& centered,
          const MatrixXd& cov, const SelfAdjointEigenSolver<MatrixXd>& eig)
{
    spdlog::debug(msg);
    spdlog::debug("Input matrix:\n{}", mat);
    spdlog::debug("Centered matrix:\n{}", centered);
    spdlog::debug("Covariance matrix:\n{}", cov);
    spdlog::debug("Eigenvalues:\n{}", eig.eigenvalues());
    spdlog::debug("Eigenvectors:\n{}", eig.eigenvectors());
}

static
void nominal(const MatrixXd& mat, const std::string& name)
{
    {
        MatrixXd centered = mat.rowwise() - mat.colwise().mean();
        MatrixXd cov = ( centered.adjoint() * centered ) / ((double)mat.rows() - 1.0);
        SelfAdjointEigenSolver<MatrixXd> eig(cov);
        dump("nominal " + name, mat, centered, cov, eig);
        CHECK(eig.eigenvalues().size() == mat.cols());
        double last_ev = 0;
        for (const auto& ev : eig.eigenvalues()) {
            CHECK(last_ev < ev);
            last_ev = ev;
        }
    }
    {
        auto eig = WireCell::Array::pca(mat);
        CHECK(eig.eigenvalues().size() == mat.cols());
        double last_ev = 0;
        for (const auto& ev : eig.eigenvalues()) {
            CHECK(last_ev < ev);
            last_ev = ev;
        }
    }
}


TEST_CASE("test some pca methods") {

    SUBCASE("6x4") {
        Eigen::MatrixXd mat(6, 4);
        // http://www.cs.cornell.edu/courses/cs4786/2015sp/code/lec03PCAExample.cpp
        mat <<
            5.05595969, 7.43836919, 4.71878778, 6.88171439,
            4.70296456, 7.18632055, 5.27801363, 7.01708935,
            5.11956978, 7.92647724, 5.35804499, 6.94302973,
            7.45143999, 5.24763548, 7.54816947, 4.81724412,
            6.99920644, 4.49667605, 6.82754446, 4.84366695,
            7.30397357, 5.17392094, 6.97290216, 5.45113576;
        SUBCASE("nominal") {
            nominal(mat, "6x4");
        }
    }
    SUBCASE("4x3") {
        MatrixXd mat(4, 3);
        // https://blog.clairvoyantsoft.com/eigen-decomposition-and-pca-c50f4ca15501 
        mat <<
            20, 50, 4.0,
            25, 60, 5.0,
            30, 65, 5.5,
            40, 75, 6.0;
        SUBCASE("nominal") { nominal(mat, "4x3"); }
    }
    SUBCASE("5x3") {
        MatrixXd mat(5,3);
        // https://blog.demofox.org/2022/07/10/programming-pca-from-scratch-in-c/
        mat <<
            90.000000, 60.000000, 90.000000,
            90.000000, 90.000000, 30.000000,
            60.000000, 60.000000, 60.000000,
            60.000000, 60.000000, 90.000000,
            30.000000, 30.000000, 30.000000;
        SUBCASE("nominal") { nominal(mat, "5x3"); }
    }

}
