// Util functions to rebin a tensor.  See Rebinner for a DFP node.
#pragma once

#include "WireCellSpng/Torch.h"

namespace WireCell::SPNG::Rebin {

    // How to interpret sample normalization.  This enum is "morally equivalent"
    // to a similar one in the LMN namespace.
    enum class Normalization {
        kIntegral=1,
        kInterpolation,
        kMaximum,
    };

    Normalization parse_mode(const std::string& mode_str) {
        if (mode_str == "integral") return Normalization::kIntegral;
        if (mode_str == "interpolation") return Normalization::kInterpolation;
        if (mode_str == "maximum") return Normalization::kMaximum;
        throw std::invalid_argument("Invalid sampling mode: " + mode_str + 
                                    ". Must be 'integral', 'interpolation', or 'maximum'.");
    }

    /**
     * @brief Downsamples a tensor along a given dimension by an integer factor.
     *
     * Each output sample is derived from a run of 'factor' input samples.
     *
     * @param input The input tensor.
     * @param factor The downsampling factor (must be > 1).
     * @param dim The dimension along which to downsample.
     * @param mode The normalization interpretation.
     * @return The downsampled tensor.
     */
    torch::Tensor downsample(
        const torch::Tensor& input, 
        int64_t factor, 
        int64_t dim=-1, 
        Normalization mode = Normalization::kIntegral)
    {
        if (factor <= 1) {
            throw std::invalid_argument("Downsampling factor must be > 1.");
        }

        // Normalize the dimension index
        int64_t normalized_dim = (dim < 0) ? input.dim() + dim : dim;
        if (normalized_dim < 0 || normalized_dim >= input.dim()) {
            throw std::invalid_argument("Dimension index out of range.");
        }
    
        int64_t original_size = input.size(normalized_dim);
        if (original_size % factor != 0) {
            throw std::invalid_argument("Input size along dimension " + std::to_string(dim) + 
                                        " (" + std::to_string(original_size) + 
                                        ") must be divisible by the factor " + std::to_string(factor) + ".");
        }

        // 1. Reshape to split the dimension into (new_size, factor)
        std::vector<int64_t> new_shape = input.sizes().vec();
        new_shape[normalized_dim] = original_size / factor; // The new size
        new_shape.insert(new_shape.begin() + normalized_dim + 1, factor); // The factor dimension

        torch::Tensor reshaped = input.contiguous().view(new_shape);

        // 2. Aggregate along the 'factor' dimension (normalized_dim + 1)
        torch::Tensor result;
        int64_t factor_dim = normalized_dim + 1;

        switch (mode) {
        case Normalization::kIntegral: // sum over samples
            result = reshaped.sum(factor_dim);
            break;
        case Normalization::kInterpolation: // average over samples
            result = reshaped.mean(factor_dim);
            break;
        case Normalization::kMaximum: // maximum of samples
            // .values() gets the tensor of max values
            result = std::get<0>(reshaped.max(factor_dim)); 
            break;
        }

        return result;
    }

    /**
     * @brief Upsamples a tensor along a given dimension by an integer factor.
     *
     * Each input sample is replicated 'factor' times.
     *
     * @param input The input tensor.
     * @param factor The upsampling factor (must be > 1).
     * @param dim The dimension along which to upsample.
     * @param mode The normalization interpretation.
     * @return The upsampled tensor.
     */
    torch::Tensor upsample(
        const torch::Tensor& input, 
        int64_t factor, 
        int64_t dim, 
        Normalization mode = Normalization::kIntegral)
    {
        if (factor <= 1) {
            throw std::invalid_argument("Upsampling factor must be > 1.");
        }

        torch::Tensor working_tensor = input;
    
        // Handle the "integral" mode requirement: original sample divided by the factor
        if (mode == Normalization::kIntegral) {
            // Divide by the factor *before* replication
            working_tensor = working_tensor.div(static_cast<float>(factor));
        }
    
        // For "interpolation" and "maximum" modes, the value is just the original sample.
        // The operation is the same: replicate the sample 'factor' times.
        torch::Tensor result = torch::repeat_interleave(
            working_tensor,
            factor,
            dim
            );

        return result;
    }

}
