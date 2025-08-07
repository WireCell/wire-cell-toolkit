#include "WireCellUtil/Logging.h"
#include "WireCellUtil/KSTest.h"
#include "WireCellUtil/doctest.h"

#include <random> // For generating random samples
#include <cmath>  // For std::erf, std::sqrt, M_SQRT1_2

// Define the CDF for a Gaussian distribution (standard normal)
// CDF(x) = 0.5 * (1 + erf((x - mean) / (std_dev * sqrt(2))))
static double gaussian_cdf(double x, double mean, double std_dev) {
    return 0.5 * (1.0 + std::erf((x - mean) / (std_dev * M_SQRT1_2))); // M_SQRT1_2 is 1/sqrt(2)
}

// A Gaussian PDF for demonstration
static double gaussian_pdf(double x, double mean, double std_dev) {
    const double xms = (x-mean)/std_dev;
    return (1.0 / std::sqrt(2.0 * std_dev * M_PI)) * std::exp(-0.5 * xms * xms);
}
static double normal_pdf(double x) {
    return (1.0 / std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * x * x);
}

using spdlog::debug;

using namespace WireCell;

TEST_CASE("kstest one sample")
{
    // --- Setup for the analytical reference CDF (Standard Normal: mean=0, std_dev=1) ---
    auto standard_normal_cdf = [](double x) { return gaussian_cdf(x, 0.0, 1.0); };

    // Create the KS test object with the reference CDF
    KSTest1Sample ks_test(standard_normal_cdf);

    // --- Example 1: Samples drawn from the *same* distribution (expect high p-value) ---
    debug("--- Example 1: Samples from the same distribution (Standard Normal) ---");
    std::mt19937 gen1(std::random_device{}());
    std::normal_distribution<> dist1(0.0, 1.0); // Standard Normal

    std::set<double> samples1;
    for (int i = 0; i < 20; ++i) { // 100 samples
        samples1.insert(dist1(gen1));
    }

    double d_stat1 = ks_test.d_statistic(samples1);
    double p_val1 = ks_test.pvalue(d_stat1, samples1.size());

    debug("Number of samples: {}", samples1.size());
    debug("D-statistic: {}", d_stat1);
    debug("P-value: {}", p_val1);
    debug("Interpretation (alpha=0.05): {}",
          (p_val1 > 0.05 ? "Fail to reject null hypothesis (samples likely from Standard Normal)" :
           "Reject null hypothesis (samples likely NOT from Standard Normal)"));

    // --- Example 2: Samples drawn from a *different* distribution (expect low p-value) ---
    debug("--- Example 2: Samples from a different distribution (Normal, mean=2.0) ---");
    std::mt19937 gen2(std::random_device{}());
    std::normal_distribution<> dist2(2.0, 1.0); // Different mean

    std::set<double> samples2;
    for (int i = 0; i < 100; ++i) { // 100 samples
        samples2.insert(dist2(gen2));
    }

    double d_stat2 = ks_test.d_statistic(samples2);
    double p_val2 = ks_test.pvalue(d_stat2, samples2.size());


    debug("Number of samples: {}", samples2.size());
    debug("D-statistic: {}", d_stat2);
    debug("P-value: {}", p_val2);
    debug("Interpretation (alpha=0.05): {}",
          (p_val2 > 0.05 ?
           "Fail to reject null hypothesis (samples likely from Standard Normal)" :
           "Reject null hypothesis (samples likely NOT from Standard Normal)"));

    
    for (double known_d_stat = 1.0; known_d_stat >= 0.001; known_d_stat /= 10.0) {

        debug("--- Example 3: Manual check (D={}, N=100) ---", known_d_stat);
        size_t known_n = 100;
        double known_p_val = ks_test.pvalue(known_d_stat, known_n);
        debug("D-statistic: {} for N={}", known_d_stat, known_n);
        debug("P-value: {}", known_p_val);
    }

    debug("--- Example 4: Empty samples ---");
    try {
        std::set<double> empty_samples;
        ks_test.d_statistic(empty_samples);
    } catch (const std::invalid_argument& e) {
        debug("Caught expected exception: {}", e.what());
    }

    debug("--- Example 5: Zero samples for pvalue ---");
    try {
        ks_test.pvalue(0.1, 0);
    } catch (const std::invalid_argument& e) {
        debug("Caught expected exception: {}", e.what());
    }
}

TEST_CASE("kstest two sample")
{
    // --- Example 1: Two samples from the *same* distribution (expect high p-value) ---
    debug("--- Example 1: Samples from the same distribution ---");
    std::mt19937 gen1(std::random_device{}());
    std::normal_distribution<> dist1(0.0, 1.0); // Standard Normal

    std::set<double> ref_samples1;
    for (int i = 0; i < 100; ++i) {
        ref_samples1.insert(dist1(gen1));
    }

    std::set<double> test_samples1;
    for (int i = 0; i < 120; ++i) { // Different size, same distribution
        test_samples1.insert(dist1(gen1));
    }

    KSTest2Sample ks_test1(ref_samples1);
    double d_stat1 = ks_test1.measure(test_samples1);
    double p_val1 = ks_test1.pvalue(d_stat1, test_samples1.size());

    debug("Reference samples N: {}", ref_samples1.size());
    debug("Test samples N: {}", test_samples1.size());
    debug("D-statistic: {}", d_stat1);
    debug("P-value: {}", p_val1);
    debug("Interpretation (alpha=0.05): {}",
          (p_val1 > 0.05 ?
           "Fail to reject null hypothesis (samples likely from same distribution)" :
           "Reject null hypothesis (samples likely NOT from same distribution)"));

    // --- Example 2: Two samples from *different* distributions (expect low p-value) ---
    debug("--- Example 2: Samples from different distributions ---");
    std::mt19937 gen2(std::random_device{}());
    std::normal_distribution<> dist_ref(0.0, 1.0); // Standard Normal
    std::normal_distribution<> dist_test(0.5, 1.2); // Different mean and std dev

    std::set<double> ref_samples2;
    for (int i = 0; i < 100; ++i) {
        ref_samples2.insert(dist_ref(gen2));
    }

    std::set<double> test_samples2;
    for (int i = 0; i < 100; ++i) {
        test_samples2.insert(dist_test(gen2));
    }

    KSTest2Sample ks_test2(ref_samples2);
    double d_stat2 = ks_test2.measure(test_samples2);
    double p_val2 = ks_test2.pvalue(d_stat2, test_samples2.size());

    debug("Reference samples N: {}", ref_samples2.size());
    debug("Test samples N: {}", test_samples2.size());
    debug("D-statistic: {}", d_stat2);
    debug("P-value: {}", p_val2);
    debug("Interpretation (alpha=0.05): {}",
          (p_val2 > 0.05 ?
           "Fail to reject null hypothesis" :
           "Reject null hypothesis"));


    // --- Example 3: Edge case - very small D-statistic (expect p-value ~ 1.0) ---
    debug("--- Example 3: Very small D-statistic ---");
    std::set<double> ref_small_d;
    for (int i = 0; i < 50; ++i) ref_small_d.insert(dist1(gen1));
    std::set<double> test_small_d;
    for (int i = 0; i < 50; ++i) test_small_d.insert(dist1(gen1));

    KSTest2Sample ks_test_small_d(ref_small_d);
    double d_stat_small_d = ks_test_small_d.measure(test_small_d);
    double p_val_small_d = ks_test_small_d.pvalue(d_stat_small_d, test_small_d.size());

    debug("Reference samples N: {}", ref_small_d.size());
    debug("Test samples N: {}", test_small_d.size());
    debug("D-statistic: {}", d_stat_small_d);
    debug("P-value: {}", p_val_small_d);
    debug("Interpretation: Small D-stat, large p-value. (Correct)");

}


TEST_CASE("kstest discrete cdf")
{
    auto my_pdf = normal_pdf;
    double start_coord = -5.0;
    double step_size = 0.01;
    size_t num_points = 1000;

    // Calculate the discrete CDF
    std::vector<double> cdf = discrete_cdf(my_pdf, start_coord, step_size, num_points);

    // Print some values to demonstrate
    debug("Discrete CDF calculated with {} points", num_points);
    debug("Start coordinate: {}, Step size: {}", start_coord, step_size);
    debug("CDF at x={} is {}", start_coord + 1 * step_size, cdf[0]);
    debug("CDF at x={} is {}", start_coord + 500 * step_size, cdf[499]);
    debug("CDF at x={} is {}", start_coord + num_points * step_size, cdf.back());

}

TEST_CASE("kstest like")
{
    const size_t nsamples = 100;
    std::vector<double> pdf1(nsamples), pdf2(nsamples), pdf3(nsamples);
    const double step = 10.0/nsamples;
    const double x0 = -10.0/(0.5*nsamples);
    double mean1=0, mean2=0.001, mean3=2.0;
    for (size_t ind=0; ind<nsamples; ++ind) {
        const double x = x0 + ind*step;
        pdf1[ind] = gaussian_pdf(x, mean1, 1.0);
        pdf2[ind] = gaussian_pdf(x, mean2, 1.0);
        pdf3[ind] = gaussian_pdf(x, mean3, 1.0);
    }
    const double d12 = kslike_compare(pdf1, pdf2);
    const double p12 = ks_pvalue(d12, nsamples, nsamples);
    const double d13 = kslike_compare(pdf1, pdf3);
    const double p13 = ks_pvalue(d13, nsamples, nsamples);

    debug("d12={}, p12={}", d12, p12);
    debug("d13={}, p13={}", d13, p13);
}
