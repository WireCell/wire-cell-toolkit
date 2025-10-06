// TorchLMN.cxx
#include "WireCellSpng/TorchLMN.h"
#include "WireCellUtil/Exceptions.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <cstdio> // for snprintf used in exception formatting


namespace WireCell::SPNG::LMN {

    double gcd(double a, double b, double eps)
    {
        // Ensure inputs are non-negative for standard GCD algorithm interpretation
        a = std::abs(a);
        b = std::abs(b);

        if (a < eps) {
            return b;
        }
        // Recursive Euclidean algorithm
        return LMN::gcd(fmod(b, a), a, eps);
    }

    size_t rational(double Ts, double Tr, double eps)
    {
        const double dT = std::abs(Ts - Tr);
        const double common_divisor = LMN::gcd(Tr, dT, eps);
    
        if (common_divisor < eps) {
            // If dT is near zero (Ts ~ Tr), common_divisor ~ Tr. 
            // If Tr and dT are orthogonal within eps, common_divisor might be near zero.
            if (common_divisor == 0) return 0;
        }

        const double n = dT / common_divisor;
        const size_t rn = std::round(n);
        const double err1 = std::abs(n - rn);
    
        if (err1 > eps) {
            raise<ValueError>("gcd error one too big %f > %f, rn=%zu, n=%f", err1, eps, rn, n);
            return 0;
        }

        // Calculate Ns = rn * Tr / dT
        double Ns_float;
        if (dT < eps) {
            // Handle Ts == Tr case explicitly if floating point error resulted in small but non-zero dT
            Ns_float = 1.0; 
        } else {
            Ns_float = (double)rn * Tr / dT;
        }

        size_t rNs = std::round(Ns_float);
        const double err2 = std::abs(rNs - Ns_float);
        if (err2 > eps) {
            raise<ValueError>("gcd error two too big %f > %f", err2, eps);
            return 0;
        }

        return rNs;    
    }

// Helper function to get size along axis
    int64_t get_size_at_axis(const torch::Tensor& in, int64_t axis) {
        if (axis < 0 || axis >= in.dim()) {
            throw std::out_of_range("Axis index out of bounds");
        }
        return in.size(axis);
    }

// ----------------------------------------------------------------------
// Tensor implementations (Time/Spatial Domain Resize)
// ----------------------------------------------------------------------

    torch::Tensor resize(const torch::Tensor& in, int64_t Nr, int64_t axis)
    {
        TORCH_CHECK(in.dim() <= 2, "LMN::resize expected 1D or 2D tensor.");
    
        if (in.dim() == 0) return in.clone(); // Scalar case

        const int64_t Ns = get_size_at_axis(in, axis);

        if (Ns == Nr) {
            return in.clone();
        }
    
        const int64_t N_min = std::min(Ns, Nr);
    
        // 1. Determine output shape and initialize zero tensor
        std::vector<int64_t> out_size = in.sizes().vec();
        if (in.dim() > 0) {
            out_size[axis] = Nr;
        }
    
        torch::Tensor rs = torch::zeros(out_size, in.options());
    
        // 2. Copy existing elements (narrowing the input and output)
        if (in.dim() == 1) {
            // 1D tensor
            rs.narrow(0, 0, N_min).copy_(in.narrow(0, 0, N_min));
        } else {
            // 2D tensor (or higher, though constrained by TORCH_CHECK)
            if (axis == 0) { // Rows
                // Copy [0:N_min, all]
                rs.narrow(0, 0, N_min).copy_(in.narrow(0, 0, N_min));
            }
            else { // Columns (axis == 1)
                // Copy [all, 0:N_min]
                rs.narrow(1, 0, N_min).copy_(in.narrow(1, 0, N_min));
            }
        }
    
        return rs;
    }

// ----------------------------------------------------------------------
// Standard vector helpers
// ----------------------------------------------------------------------

    std::vector<float> resize(const std::vector<float>& in, size_t Nr)
    {
        size_t Ns = in.size();
        if (Ns == Nr) return in;
        size_t N_min = std::min(Nr,Ns);
    
        // Initialize with size Nr and default value 0.0
        std::vector<float> rs(Nr, 0.0f); 
        std::copy(in.begin(), in.begin()+N_min, rs.begin());

        return rs;
    }

    void fill_constant(std::vector<float>::iterator begin,
                       std::vector<float>::iterator end,
                       float value)
    {
        std::fill(begin, end, value);
    }

    void fill_linear(std::vector<float>::iterator begin,
                     std::vector<float>::iterator end,
                     float first, float last)
    {
        size_t N = std::distance(begin, end);
        if (N == 0) return;
    
        float step = (last - first) / (float)N; 
    
        for (size_t ind = 0; ind < N; ++ind) {
            *(begin + ind) = first + ind * step;
        }
    }


// ----------------------------------------------------------------------
// Frequency Domain Resample (Tensor)
// ----------------------------------------------------------------------

    torch::Tensor resample(const torch::Tensor& in, int64_t Nr, int64_t axis)
    {
        TORCH_CHECK(in.is_complex(), "Input tensor must be complex for frequency domain resampling.");
        TORCH_CHECK(in.dim() <= 2, "LMN::resample expected 1D or 2D tensor.");
    
        const int64_t Ns = get_size_at_axis(in, axis);

        if (Ns == Nr) {
            return in.clone();
        }
    
        // N_half is calculated based on the dimension of the smaller frequency window
        int64_t N_half_limit = (Nr > Ns) ? Ns : Nr;
        size_t N_half = nhalf(N_half_limit);
    
        // 1. Determine output shape and initialize zero tensor
        std::vector<int64_t> out_size = in.sizes().vec();
        out_size[axis] = Nr;
        torch::Tensor rs = torch::zeros(out_size, in.options()); 

        // Define slice indices:
        // Positive half size: N_half + 1 (DC component included)
        int64_t pos_size = N_half + 1; 
        int64_t pos_start = 0;

        // Negative half size: N_half 
        int64_t neg_size = N_half;
        int64_t neg_start_in = Ns - N_half; // Start index in source (in)
        int64_t neg_start_rs = Nr - N_half; // Start index in target (rs)
    
    
        // Lambda to safely extract a slice along the specified axis
        auto slice = [&](const torch::Tensor& t, int64_t start, int64_t length) -> torch::Tensor {
            if (t.dim() == 1) {
                return t.narrow(0, start, length);
            } else if (axis == 0) {
                // Slicing rows
                return t.narrow(0, start, length);
            } else {
                // Slicing columns (axis == 1)
                return t.narrow(1, start, length);
            }
        };
    
        // 2. Copy Positive Half Frequencies (DC up to N_half)
        slice(rs, pos_start, pos_size).copy_(slice(in, pos_start, pos_size));

        // 3. Copy Negative Half Frequencies (highest negative frequencies, wrapped)
        // We only copy if N_half > 0
        if (neg_size > 0) {
            slice(rs, neg_start_rs, neg_size).copy_(slice(in, neg_start_in, neg_size));
        }
    
        // Note: Consistent with the original Eigen implementation structure,
        // the Nyquist bin handling (if Ns is even) is omitted/skipped 
        // based on the definition of nhalf, and requires external clarification if needed.

        return rs;
    }
    

// ----------------------------------------------------------------------
// Frequency Domain Resample (Standard Vector)
// ----------------------------------------------------------------------

    std::vector<std::complex<float>>
    resample(const std::vector<std::complex<float>>& in, size_t Nr)
    {
        size_t Ns = in.size();

        size_t N_half_limit = (Nr > Ns) ? Ns : Nr;
        size_t N_half = nhalf(N_half_limit);

        std::vector<std::complex<float>> rs(Nr);

        // Copy Positive Half (DC to N_half, total size N_half + 1)
        std::copy(in.begin(), in.begin() + N_half + 1, rs.begin());
    
        // Copy Negative Half (last N_half elements)
        // Source: using reverse iterators to pick up the tail
        if (N_half > 0) {
            std::copy(in.rbegin(), in.rbegin() + N_half, rs.rbegin());
        }

        // FIXME: deal with Nyquist bin. (Note maintained from original code.)
        return rs;
    }
}
