// TorchLMN.cxx
#include "WireCellSpng/TorchLMN.h"
#include "WireCellUtil/Exceptions.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <cstdio> // for snprintf used in exception formatting


namespace WireCell::SPNG::LMN {

    // Helper function to get size along axis
    static
    int64_t get_size_at_axis(const torch::Tensor& in, int64_t axis) {
        if (axis < 0 || axis >= in.dim()) {
            throw std::out_of_range("Axis index out of bounds");
        }
        return in.size(axis);
    }

    /// Return the "half size" number of samples in a spectrum of full size N.
    /// This excludes the "zero frequency" sample and the "Nyquist bin", if one
    /// exists.
    static
    size_t nhalf(size_t N) {
        if (N%2) {
            return (N-1)/2;     // odd
        }
        return (N-2)/2;         // even
    }


    double gcd(double a, double b, double eps)
    {
        // Ensure inputs are non-negative for standard GCD algorithm interpretation
        a = std::abs(a);
        b = std::abs(b);

        if (a < eps) {
            return b;
        }
        // Recursive Euclidean algorithm
        return gcd(fmod(b, a), a, eps);
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


    torch::Tensor resample(const torch::Tensor& in, int64_t Nr, int64_t axis)
    {
        TORCH_CHECK(in.is_complex(), "Input tensor must be complex for frequency domain resampling.");
        TORCH_CHECK(in.dim() <= 2, "LMN::resample expected 1D or 2D tensor.");
    
        const int64_t Ns = get_size_at_axis(in, axis);

        if (Ns == Nr) {
            return in.clone();
        }
    
        // N_half is based on the size of the common spectrum (excluding Nyquist bins)
        const int64_t N_half_limit = std::min(Ns, Nr);
        const size_t H = LMN::nhalf(N_half_limit);
        const int64_t P_size = H + 1; // DC + H positives
        const int64_t L_size = H;     // H negatives
    
        // 1. Determine output shape and initialize zero tensor
        std::vector<int64_t> out_size = in.sizes().vec();
        out_size[axis] = Nr;
        torch::Tensor rs = torch::zeros(out_size, in.options()); 

        // Lambda to safely extract/assign a slice along the specified axis
        auto slice = [&](const torch::Tensor& t, int64_t start, int64_t length) -> torch::Tensor {
            if (t.dim() == 1) {
                return t.narrow(0, start, length);
            } else if (axis == 0) {
                return t.narrow(0, start, length); // Slicing rows
            } else {
                return t.narrow(1, start, length); // Slicing columns (axis == 1)
            }
        };
    
        // 2. Standard Copy (DC and primary halves)
    
        // Copy Positives: [0, P_size)
        slice(rs, 0, P_size).copy_(slice(in, 0, P_size));

        // Copy Negatives: [N_s - L_size, N_s) -> [N_r - L_size, N_r)
        if (L_size > 0) {
            slice(rs, Nr - L_size, L_size).copy_(slice(in, Ns - L_size, L_size));
        }
    
        // 3. Nyquist Correction
        if (Ns % 2 == 0 && Nr > Ns) {
            // Case B1: Upsampling, Input Nyquist NQ_in exists (Ns even).
            int64_t I_NQ_in = Ns / 2;
            torch::Tensor NQ_val = slice(in, I_NQ_in, 1);
        
            // Split NQ_in into highest positive bin of new spectrum (index Ns/2) 
            // and precursor of lowest negative bin (index Nr - Ns/2)
        
            // 1. Highest positive bin in rs: index Ns/2 
            slice(rs, I_NQ_in, 1).add_(NQ_val / 2.0f);
        
            // 2. Lowest negative bin precursor in rs: index Nr - (Ns - Ns/2) = Nr - Ns/2
            int64_t I_NQ_neg_rs = Nr - I_NQ_in;
            slice(rs, I_NQ_neg_rs, 1).add_(NQ_val / 2.0f);
        }
    
        else if (Nr % 2 == 0 && Ns > Nr) {
            // Case B2: Downsampling, Output Nyquist NQ_out required (Nr even).
            int64_t I_NQ_out = Nr / 2; // Index in rs
        
            // Combine contributions from the folded frequencies: 
            // highest positive (in[Nr/2]) and lowest negative (in[Ns - Nr/2])
            int64_t I_P_fold = Nr / 2; 
            int64_t I_N_fold = Ns - Nr / 2;

            if (I_P_fold < Ns && I_N_fold < Ns) { 
                torch::Tensor P_fold = slice(in, I_P_fold, 1);
                torch::Tensor N_fold = slice(in, I_N_fold, 1);
            
                // Use average to protect against small numeric deviations from symmetry.
                torch::Tensor NQ_combined = 0.5*(P_fold + N_fold);
            
                // The Nyquist bin must be real-valued. Extract real part and zero imaginary part.
                torch::Tensor NQ_real_part = torch::real(NQ_combined);
                torch::Tensor NQ_imag_part = torch::zeros_like(NQ_real_part);
                torch::Tensor NQ_final = torch::complex(NQ_real_part, NQ_imag_part);

                slice(rs, I_NQ_out, 1).copy_(NQ_final);
            }
        }
        return rs;
    }
    

}
