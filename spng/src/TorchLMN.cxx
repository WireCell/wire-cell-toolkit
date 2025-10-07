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


    torch::Tensor resize(const torch::Tensor& in, int64_t Nr, int64_t axis)
    {
        if (in.dim() == 0) return in.clone(); // Scalar case

        // Validate axis constraint
        TORCH_CHECK(axis >= 0 && axis < in.dim(), "Axis index out of bounds in LMN::resize");

        const int64_t Ns = in.size(axis);

        if (Ns == Nr) {
            return in.clone();
        }

        const int64_t N_min = std::min(Ns, Nr);

        // 1. Determine output shape
        std::vector<int64_t> out_size = in.sizes().vec();
        out_size[axis] = Nr;

        // 2. Initialize zero tensor
        torch::Tensor rs = torch::zeros(out_size, in.options());

        // 3. Copy existing elements using narrow along the specified axis
        // rs.narrow(axis, start=0, length=N_min) copies from in.narrow(axis, start=0, length=N_min)
        rs.narrow(axis, 0, N_min).copy_(in.narrow(axis, 0, N_min));

        return rs;
    }


    torch::Tensor slice_tensor(const torch::Tensor& t, int64_t start, int64_t length, int64_t axis) {
        if (t.dim() == 1) {
            return t.narrow(0, start, length);
        } else if (axis == 0) {
            return t.narrow(0, start, length); // Slicing rows
        } else {
            return t.narrow(1, start, length); // Slicing columns
        }
    };


    torch::Tensor resize_middle(const torch::Tensor& in, int64_t Nr, int64_t axis)
    {
        TORCH_CHECK(in.dim() <= 2, "LMN::resize_middle expected 1D or 2D tensor.");

        if (in.dim() == 0) return in.clone(); 

        const int64_t Ns = get_size_at_axis(in, axis);

        if (Ns == Nr) {
            return in.clone();
        }

        // Determine the size of the positive and negative halves to preserve.
        const int64_t N_min = std::min(Ns, Nr);

        // H_limit is the half size excluding index 0, but including a potential center bin for even N.
        // Since the definition of "center" is asymmetric around N/2, we calculate the 
        // split point based on standard convention where index 0 is the start of the "positive" half.

        // We want to keep the samples closest to index 0.
        // If N is the size we are limiting by (N_min), we keep roughly N/2 samples at the start
        // and N/2 samples at the end.

        int64_t P_size; // Size of the first chunk (from index 0)
        int64_t L_size; // Size of the last chunk (the tail)

        if (N_min % 2 == 1) {
            // Odd N_min: (N_min - 1)/2 bins on either side of index 0.
            // P_size = 1 (index 0) + (N_min - 1)/2 = (N_min + 1) / 2
            // L_size = (N_min - 1) / 2
            P_size = (N_min + 1) / 2;
            L_size = (N_min - 1) / 2;
        } else {
            // Even N_min: (N_min - 2)/2 bins on either side of a central bin/gap.
            // If we split symmetrically around index 0, we take N_min/2 at the front
            // and N_min/2 at the back. 
            // Example N=4: [0, 1, 2, 3]. Keep [0, 1] and [2, 3].
            P_size = N_min / 2;
            L_size = N_min / 2;
        }

        // 1. Determine output shape and initialize zero tensor
        std::vector<int64_t> out_size = in.sizes().vec();
        if (in.dim() > 0) {
            out_size[axis] = Nr;
        }

        torch::Tensor rs = torch::zeros(out_size, in.options());

        // 2. Copy the first chunk (P_size elements starting at 0)
        // Source range: [0, P_size)
        // Target range: [0, P_size)
        slice_tensor(rs, 0, P_size, axis).copy_(slice_tensor(in, 0, P_size, axis));

        // 3. Copy the last chunk (L_size elements)
        if (L_size > 0) {
            // Source range: [Ns - L_size, Ns)
            int64_t Ns_start = Ns - L_size;

            // Target range: [Nr - L_size, Nr)
            int64_t Nr_start = Nr - L_size;

            slice_tensor(rs, Nr_start, L_size, axis).copy_(slice_tensor(in, Ns_start, L_size, axis));
        }

        // The gap in the middle [P_size, Nr - L_size) is handled by zero initialization 
        // (if Nr > Ns) or truncation (if Nr < Ns, where the gap is simply omitted).

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
        // auto slice = [&](const torch::Tensor& t, int64_t start, int64_t length) -> torch::Tensor {
        //     if (t.dim() == 1) {
        //         return t.narrow(0, start, length);
        //     } else if (axis == 0) {
        //         return t.narrow(0, start, length); // Slicing rows
        //     } else {
        //         return t.narrow(1, start, length); // Slicing columns (axis == 1)
        //     }
        // };
        // Intern axis and call the external function
        auto slice = [&](const torch::Tensor& t, int64_t start, int64_t length) -> torch::Tensor {
            return slice_tensor(t, start, length, axis);
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
    
    torch::Tensor resample_interval(const torch::Tensor& interval, 
                                    double Ts, double Tr,
                                    int64_t axis,
                                    Normalization ni,
                                    bool fourier_space,
                                    double eps)
    {
        TORCH_CHECK(interval.is_floating_point(), "Input tensor must be real-valued (floating point).");
        TORCH_CHECK(interval.dim() <= 2, "Input tensor must be 1D or 2D.");
        TORCH_CHECK(axis >= 0 && axis < interval.dim(), "Axis index out of bounds.");

        const int64_t Ns = interval.size(axis);

        if (std::abs(Ts - Tr) < eps) {
            // No resampling needed if periods are identical within tolerance
            return interval.clone();
        }
    
        // 1. Calculate rational factor (Ns_rat)
        // Ns_rat is the minimum size that allows rational resampling based on periods Ts and Tr.
        // NOTE: The LMN rational function returns the target size Ns_rat * Rr / Rt where Rt/Rr = Ts/Tr.
        // The rational function definition in the original code is slightly confusingly named 'Ns' (rNs).
        // Let's call the result R_size_factor, which is the necessary size N needed 
        // such that N * Tr / Ts is an integer ratio.
        size_t R_size_factor = LMN::rational(Ts, Tr, eps);

        if (R_size_factor == 0) {
            raise<ValueError>("LMN::rational failed to find a rational factor for Ts=%f, Tr=%f", Ts, Tr);
        }
    
        // 2. Determine target rational size (N_rat)
        // Find the smallest size >= Ns that is divisible by R_size_factor

        size_t N_rat = LMN::nbigger(Ns, R_size_factor);

        // 3. Calculate final output size (Nr)
        // New size Nr = N_rat * Ts / Tr. This must be an integer if R_size_factor was correct.
        double Nr_float = (double)N_rat * Ts / Tr;
        int64_t Nr = (int64_t)std::round(Nr_float);
        if (std::abs(Nr_float - Nr) > eps) {
            // Sanity check, should not happen if rational() and nbigger() are correct
            raise<ValueError>("Rational size calculation failed: Nr=%f is not an integer", Nr_float);
        }

        // 4. Pad input tensor to rational size N_rat
        torch::Tensor padded = LMN::resize(interval, N_rat, axis);

        // 5. Apply FFT
        torch::Tensor freq_domain = torch::fft::fft(padded, (int64_t)N_rat, axis);

        // 6. Resample in Fourier space
        // This performs the zero-padding/truncation and Nyquist handling.
        torch::Tensor resampled_freq = LMN::resample(freq_domain, Nr, axis);

        // 7. Apply inverse FFT and return real part
        torch::Tensor time_domain_complex = torch::fft::ifft(resampled_freq, Nr, axis);
    
        // Return only the real part, ensuring the output type matches the input precision
        // See LMN paper.
        auto raw = torch::real(time_domain_complex).to(interval.dtype());
        if (ni == Normalization::kInterpolation) {
            return raw * ((double)Nr) / ((double)Ns);
        }
        if (ni == Normalization::kIntegral) {
            return raw;
        }
        if (ni == Normalization::kEnergy) {
            return raw * sqrt(((double)Nr) / ((double)Ns));
        }
        return raw;
    }

}
