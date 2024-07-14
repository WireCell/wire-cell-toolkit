#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Array.h"
#include "WireCellUtil/doctest.h"
#include <iostream>

using namespace Eigen;

TEST_CASE("simple pca")
{
    // Define and initialize the 3x3 matrix
    Eigen::Matrix3d cov_matrix;
    cov_matrix << 100, 0, 0, 0, 10, 0, 0, 0, 1;

    // Compute the eigenvalues and eigenvectors
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(cov_matrix);

    if (eigenSolver.info() != Eigen::Success) {
        std::cerr << "Failed to compute eigenvalues and eigenvectors." << std::endl;
    }

    // Output the eigenvalues
    std::cout << "Eigenvalues:\n" << eigenSolver.eigenvalues() << std::endl;

    // Output the eigenvectors
    std::cout << "Eigenvectors:\n" << eigenSolver.eigenvectors() << std::endl;

    const auto eigen = WireCell::Array::pca(cov_matrix);
    auto eigen_values = eigen.eigenvalues();
    auto eigen_vectors = eigen.eigenvectors();

    for (int i = 0; i != 3; i++) {
        double norm = sqrt(eigen_vectors(0, i) * eigen_vectors(0, i) + eigen_vectors(1, i) * eigen_vectors(1, i) +
                           eigen_vectors(2, i) * eigen_vectors(2, i));
        std::cout << "WireCell " << i << " E Value " << eigen_values(i) << " E Vector " << eigen_vectors(0, i) / norm << " "
                  << eigen_vectors(1, i) / norm << " " << eigen_vectors(2, i) / norm << std::endl;
    }
}
