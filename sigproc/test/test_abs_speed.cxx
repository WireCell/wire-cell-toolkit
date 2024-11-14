#include <complex>
#include <chrono>
#include <iostream>
#include <cmath>

int main() {
    std::complex<double> x(3.0, 4.0); // Example complex number

    int numIterations = 1e9;

    // Measure the time taken by std::abs(x)
    auto start = std::chrono::high_resolution_clock::now();
    double absSum = 0;
    for (int i = 0; i < numIterations; i++) {
        absSum += std::abs(x);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto absDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    // Measure the time taken by std::sqrt(std::norm(x))
    start = std::chrono::high_resolution_clock::now();
    double sqrtSum = 0;
    for (int i = 0; i < numIterations; i++) {
        sqrtSum += std::sqrt(std::norm(x));
    }
    end = std::chrono::high_resolution_clock::now();
    auto sqrtDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    // Print the results
    std::cout << "std::abs(x) duration: " << absDuration << " milliseconds, sum = " << absSum << std::endl; 
    std::cout << "std::sqrt(std::norm(x)) duration: " << sqrtDuration << " milliseconds, sum = " << sqrtSum << std::endl;

    return 0;
}
