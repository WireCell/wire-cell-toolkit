#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Array.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include <sstream>
using spdlog::debug;

using namespace Eigen;

TEST_CASE("eigen pca simple")
{
    // Define and initialize the 3x3 matrix
    Eigen::Matrix3d cov_matrix;
    cov_matrix << 100, 0, 0, 0, 10, 0, 0, 0, 1;

    // Compute the eigenvalues and eigenvectors
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(cov_matrix);

    CHECK(eigenSolver.info() == Eigen::Success);

    // Output the eigenvalues
    {
        std::stringstream ss;
        ss << eigenSolver.eigenvalues();
        debug("Eigenvalues: {}", ss.str());
    }
    {
        std::stringstream ss;
        ss << eigenSolver.eigenvectors();
        // Output the eigenvectors
        debug("Eigenvectors: {}", ss.str());
    }

    const auto eigen = WireCell::Array::pca(cov_matrix);
    auto eigen_values = eigen.eigenvalues();
    auto eigen_vectors = eigen.eigenvectors();

    for (int i = 0; i != 3; i++) {
        double norm = sqrt(eigen_vectors(0, i) * eigen_vectors(0, i) + eigen_vectors(1, i) * eigen_vectors(1, i) +
                           eigen_vectors(2, i) * eigen_vectors(2, i));
        debug("WireCell {} eigen value={} vector=({},{},{})",
              i, eigen_values(i), 
              eigen_vectors(0, i) / norm,
              eigen_vectors(1, i) / norm,
              eigen_vectors(2, i) / norm);
    }
}
