#include "WireCellUtil/KSTest.h"

#include <stdexcept> // For std::invalid_argument
#include <cmath>     // For std::abs, std::max, std::sqrt

#include <iostream>             // debug

namespace WireCell {


    KSTest1Sample::KSTest1Sample(KSTest1Sample::cumulative_density_func ref_cdf)
        : m_ref_cdf(std::move(ref_cdf))
    {
    }
    
    double KSTest1Sample::d_statistic(const std::set<double>& samples) const
    {
        if (samples.empty()) {
            throw std::invalid_argument("Samples set cannot be empty for KS test.");
        }

        const size_t n_samples = samples.size();

        double max_diff = 0.0;

        size_t ind=0;
        for (double x_i : samples) {

            // Evaluate the theoretical CDF at x_i using the provided std::function
            double theoretical_cdf_at_xi = m_ref_cdf(x_i);

            // Calculate ECDF values at x_i (F_n(x_i)) and just before x_i (F_n(x_i^-))
            // ECDF is a step function. At x_i, it jumps from i/N to (i+1)/N.
            double ecdf_at_xi = static_cast<double>(ind + 1) / n_samples;
            double ecdf_before_xi = static_cast<double>(ind) / n_samples;

            // Calculate the two potential differences at x_i
            // The D statistic is the supremum of the absolute differences
            // between the ECDF and the theoretical CDF. We must check both
            // the jump point and the point just before the jump.
            double diff_plus = std::abs(ecdf_at_xi - theoretical_cdf_at_xi);
            double diff_minus = std::abs(ecdf_before_xi - theoretical_cdf_at_xi);

            // Update max_diff with the larger of the two differences
            max_diff = std::max({max_diff, diff_plus, diff_minus});

            ++ind;
        }

        return max_diff;
    }
    
    double KSTest1Sample::pvalue(double d_statistic, size_t n_samples) const
    {
        if (n_samples == 0) {
            throw std::invalid_argument("Number of samples (n_samples) cannot be zero for p-value calculation.");
        }


        // The scaled D-statistic, z = d * sqrt(N).
        double z = d_statistic * std::sqrt(static_cast<double>(n_samples));
        
        // very small z (very good KS test) will lead to oscillations in the
        // expansion.  We thus draw a line in the sand.  As a touch stone, for
        // D=0.01 and N=100 the series produces p=1.0.  Yet even smaller D gives
        // *smaller* p.
        if (z < 0.1) {
            return 1.0;
        }

        // Asymptotic formula for the p-value:
        // P(D > d) = 2 * (e^-2z^2 - e^-8z^2 + e^-18z^2 - ...)
        // The formula is P(D > d) = 2 * Sum_{k=1..inf} (-1)^{k-1} * exp(-2*k^2*z^2)
        double pval = 0.0;
        const int max_terms = 100; // Generous loop limit
        const double tolerance = 1e-15; // Convergence tolerance

        for (int k = 1; k <= max_terms; ++k) {
            double exponent = -2.0 * k * k * z * z;

            // Optimization: If the exponent is very large and negative, the term
            // will underflow to zero, so we can stop the loop early.
            if (exponent < -700.0) { // Approx. limit for a double
                std::cerr << "KSTest1Sample::pvalue: tiny exponent, I quit early\n";
                break;
            }
            
            double term = std::exp(exponent);
            if ((k % 2) == 0) { // If k is even
                pval -= term;
            } else { // If k is odd
                pval += term;
            }
            
            // Check for convergence based on the current term's magnitude
            if (std::abs(term) < tolerance) {
                std::cerr << "KSTest1Sample::pvalue: hit tolerance at k="<<k<<": term="<<term<<" tolerance=" << tolerance<<"\n";
                break;
            }
        }
        
        pval *= 2.0;

        // Clamp the p-value to [0, 1] to handle any small numerical inaccuracies.
        return std::min(1.0, std::max(0.0, pval));


    }

      
    KSTest2Sample::KSTest2Sample(const std::set<double>& ref_samples)
        : ref_samples_(ref_samples.begin(), ref_samples.end()) { // Copy to vector for efficiency
        if (ref_samples_.empty()) {
            throw std::invalid_argument("Reference samples set cannot be empty for KS test.");
        }
        // Ensure reference samples are sorted (std::set already provides this, but vector copy is explicit)
        std::sort(ref_samples_.begin(), ref_samples_.end());
    }

    double KSTest2Sample::measure(const std::set<double>& test_samples) const
    {
        if (test_samples.empty()) {
            throw std::invalid_argument("Test samples set cannot be empty for KS test.");
        }

        std::vector<double> test_samples_vec(test_samples.begin(), test_samples.end());
        std::sort(test_samples_vec.begin(), test_samples_vec.end()); // Ensure sorted

        double max_diff = 0.0;
        double n1 = static_cast<double>(ref_samples_.size());
        double n2 = static_cast<double>(test_samples_vec.size());

        // Combine all unique values from both samples into a single sorted list.
        // These points will be used to evaluate the ECDFs.
        std::vector<double> all_unique_values;
        all_unique_values.reserve(n1 + n2);
        std::merge(ref_samples_.begin(), ref_samples_.end(),
                   test_samples_vec.begin(), test_samples_vec.end(),
                   std::back_inserter(all_unique_values));
        all_unique_values.erase(std::unique(all_unique_values.begin(), all_unique_values.end()), all_unique_values.end());

        // Iterate through all unique values to find the maximum difference between the ECDFs.
        // ECDF(x) = (count of samples <= x) / total_samples
        for (double val : all_unique_values) {
            // ECDF for reference samples
            // std::upper_bound points to the first element *greater than* val.
            // Its distance from begin() gives the count of elements <= val.
            auto it_ref = std::upper_bound(ref_samples_.begin(), ref_samples_.end(), val);
            double ecdf_ref = static_cast<double>(std::distance(ref_samples_.begin(), it_ref)) / n1;

            // ECDF for test samples
            auto it_test = std::upper_bound(test_samples_vec.begin(), test_samples_vec.end(), val);
            double ecdf_test = static_cast<double>(std::distance(test_samples_vec.begin(), it_test)) / n2;

            max_diff = std::max(max_diff, std::abs(ecdf_test - ecdf_ref));
        }

        return max_diff;
    }



    double KSTest2Sample::pvalue(double d_statistic, size_t n_samples) const
    {
        if (n_samples == 0 || ref_samples_.empty()) {
            throw std::invalid_argument("Both sample sizes must be greater than zero for p-value calculation.");
        }

        double n1 = static_cast<double>(ref_samples_.size());
        double n2 = static_cast<double>(n_samples);

        // Calculate the effective sample size for the two-sample KS test.
        // This is used to scale the D statistic for comparison with the Kolmogorov distribution.
        double n_eff = (n1 * n2) / (n1 + n2);

        // The scaled D-statistic, z = d * sqrt(N_eff).
        double z = d_statistic * std::sqrt(n_eff);
        
        // Handle the edge case for very small z, where the asymptotic formula is unstable.
        // For these values, the p-value is essentially 1.
        if (z < 0.2) {
            return 1.0;
        }

        // Asymptotic formula for the p-value:
        // P(D > d) = 2 * Sum_{k=1..inf} (-1)^{k-1} * exp(-2*k^2*z^2)
        double pval = 0.0;
        const int max_terms = 100; // Generous loop limit
        const double tolerance = 1e-15; // Convergence tolerance

        for (int k = 1; k <= max_terms; ++k) {
            double exponent = -2.0 * k * k * z * z;

            if (exponent < -700.0) { // Approx. limit for a double to avoid underflow
                break;
            }
            
            double term = std::exp(exponent);
            if ((k % 2) == 0) { // If k is even
                pval -= term;
            } else { // If k is odd
                pval += term;
            }
            
            // Check for convergence based on the current term's magnitude
            if (std::abs(term) < tolerance) {
                break;
            }
        }
        
        pval *= 2.0;

        // Clamp the p-value to [0, 1] to handle any small numerical inaccuracies.
        return std::min(1.0, std::max(0.0, pval));
    }


    std::vector<double> discrete_cdf(std::function<double(double)> pdf,
                                     double start, double step, size_t n_samples) {
        if (n_samples == 0) {
            return {};
        }
        if (step <= 0) {
            throw std::invalid_argument("Step size must be positive.");
        }

        std::vector<double> cdf_values(n_samples);
        double current_cdf = 0.0;
    
        // We need the PDF value at the start of the first bin.
        double prev_coord = start;
        double prev_pdf_val = pdf(prev_coord);

        for (size_t i = 0; i < n_samples; ++i) {
            double current_coord = start + (i + 1) * step;
            double current_pdf_val = pdf(current_coord);

            // Use the trapezoidal rule to integrate over the current step.
            double area = 0.5 * (current_pdf_val + prev_pdf_val) * step;
            current_cdf += area;
            cdf_values[i] = current_cdf;

            // Prepare for the next iteration.
            prev_pdf_val = current_pdf_val;
        }

        // Normalize the CDF so its maximum value is 1.0.
        double total_area = cdf_values.back();
        if (total_area > 0.0) {
            for (size_t i = 0; i < n_samples; ++i) {
                cdf_values[i] /= total_area;
            }
        }
    
        return cdf_values;
    }

}
