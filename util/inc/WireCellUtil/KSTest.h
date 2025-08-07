#ifndef WIRECELL_UTIL_KSTEST
#define WIRECELL_UTIL_KSTEST

#include <vector>
#include <set>       // For std::set<double> input samples
#include <functional> // For std::function

namespace WireCell {

    /**
     * @brief Numerically integrates a PDF to produce a discrete CDF.
     *
     * This function calculates the CDF by taking the cumulative sum of the areas
     * under the PDF curve, using the trapezoidal rule for each step. The final
     * CDF is normalized to have a maximum value of 1.
     *
     * Consider wrapping the result in a WireCell::linterp to pass to KSTest1Sample.
     *
     * @param pdf The probability density function as a std::function<double(double)>.
     * @param start The starting x-coordinate of the integration.
     * @param step The size of each integration step.
     * @param n_samples The number of discrete CDF values to produce. The output
     * vector will have this size. The i-th value in the vector will correspond
     * to the CDF at coordinate `start + (i+1)*step`.
     * @return A vector of doubles representing the discrete CDF values.
     * @throws std::invalid_argument if n_samples is 0 or if step is non-positive.
     */
    std::vector<double> discrete_cdf(std::function<double(double)> pdf,
                                     double start, double step, size_t n_samples);


    /** The one-sample Kolmogorov-Smirnov (KS) test for comparing a set of raw,
     * independent, identically distributed data samples against a theoretical
     * cumulative probability distribution (CDF).
     *
     * This class calculates the KS "D" statistic, which is the maximum absolute
     * difference between the Empirical Cumulative Distribution Function (ECDF)
     * of the provided samples and the Cumulative Distribution Function
     * (CDF). It also provides a method to calculate the p-value using
     * Boost.Math, which is statistically valid for this scenario.
     */
    class KSTest1Sample {
    public:

        using cumulative_density_func = std::function<double(double)>;

        /**
         * @brief Constructor for the KS_OneSample_RawSamples_Test class.
         *
         * @param ref_cdf A std::function representing the reference CDF.
         */
        KSTest1Sample(cumulative_density_func ref_cdf);                    

        /**
         * @brief Return the KS test metric aka "D" statistic for the set of
         * independent, identically distributed samples.
         *
         * @param samples A std::set of raw data samples. std::set ensures uniqueness and sorting.
         * @return The KS "D" statistic.
         * @throws std::invalid_argument if samples is empty.
         */
        double d_statistic(const std::set<double>& samples) const;

        /**
         * @brief Converts a given KS "D" statistic and sample size into a p-value.
         *
         * The p-value indicates the probability of observing a D statistic as extreme as,
         * or more extreme than, the one calculated, assuming the null hypothesis (that the
         * samples come from the specified analytical distribution) is true.
         *
         * @param d_statistic The calculated KS "D" statistic.
         * @param n_samples The number of raw data samples used to calculate the D statistic for the one-sample case.
         * @return The p-value of the KS test.
         * @throws std::invalid_argument if n_samples is 0.
         */
        double pvalue(double d_statistic, size_t n_samples) const;

    private:

        std::function<double(double)> m_ref_cdf;
    };


    
    class KSTest2Sample {
    public:

        /**
         * @brief Constructor for the KS_TwoSample_RawSamples_Test class.
         *
         * Initializes the reference set of samples for comparison.
         *
         * @param ref_samples A const reference to a std::set of raw data samples
         * that will serve as the reference distribution.
         * @throws std::invalid_argument if ref_samples is empty.
         */
        KSTest2Sample(const std::set<double>& ref_samples);

        /**
         * @brief Calculates the Kolmogorov-Smirnov "D" statistic for two samples.
         *
         * This statistic is the maximum absolute difference between the Empirical Cumulative
         * Distribution Functions (ECDFs) of the internal reference samples and the provided
         * test samples.
         *
         * @param test_samples A std::set of raw data samples to be compared against
         * the reference samples. std::set ensures uniqueness and sorting.
         * @return The KS "D" statistic.
         * @throws std::invalid_argument if test_samples is empty.
         */
        double measure(const std::set<double>& test_samples) const;

        /**
         * @brief Converts a given KS "D" statistic and sample sizes into a p-value
         * using the asymptotic formula for the Kolmogorov distribution, scaled for
         * the two-sample case.
         *
         * @param d_statistic The calculated KS "D" statistic.
         * @param n_samples The number of samples in the *test* set (N2).
         * The *reference* sample size (N1) is taken from the constructor.
         * @return The p-value of the KS test.
         * @throws std::invalid_argument if n_samples is 0 or if the internal
         * reference sample size is 0.
         */
        double pvalue(double d_statistic, size_t n_samples) const;

    private:
        std::vector<double> ref_samples_; // Stored as vector for efficient ECDF calculation
    };


}


#endif
