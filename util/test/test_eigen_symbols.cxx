#include <iostream>
#include <Eigen/Core>

// First try with placeholders
#ifdef EIGEN_VERSION_AT_LEAST(3,4,90)
  using Eigen::placeholders::lastN;
  using Eigen::placeholders::all;
#else
  // Then fall back to main namespace
  using Eigen::lastN;
  using Eigen::all;
#endif

int main() {
    // Print Eigen version
    std::cout << "Eigen version: " 
              << EIGEN_WORLD_VERSION << "."
              << EIGEN_MAJOR_VERSION << "."
              << EIGEN_MINOR_VERSION << std::endl;
    
    // Check if EIGEN_HAS_CXX11 is defined and its value
    #ifdef EIGEN_HAS_CXX11
        std::cout << "EIGEN_HAS_CXX11 is defined with value: " << EIGEN_HAS_CXX11 << std::endl;
    #else
        std::cout << "EIGEN_HAS_CXX11 is not defined" << std::endl;
    #endif
    
    // Print C++ standard version
    std::cout << "C++ standard version: " << __cplusplus << std::endl;
    
    return 0;
}