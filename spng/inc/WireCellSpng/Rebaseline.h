// Util functions to "rebaseline" a tensor consisting or non-zero ROIs.  See Rebaseliner for a DFP node.
#pragma once

#include "WireCellSpng/Torch.h"

namespace WireCell::SPNG {

    /**
     * @brief Performs local baseline subtraction on 1D fragments > threshold.
     *
     * The function identifies contiguous 1D fragments along the specified dimension (dim) 
     * where values are greater than the threshold. It then subtracts a linear baseline 
     * model connecting the fragment's first and last sample values.
     *
     * @param tensor The input tensor.
     * @param dim The dimension along which to operate (default: last dimension, -1).
     * @param threshold The value above which a fragment is considered active.
     * @return The rebaselined tensor.
     */
    torch::Tensor rebaseline(
        const torch::Tensor& tensor, 
        int64_t dim = -1, 
        float threshold = 0.0f) 
    {
        // Make a copy to modify
        torch::Tensor output = tensor.clone();
    
        // Normalize dimension index
        int64_t normalized_dim = (dim < 0) ? output.dim() + dim : dim;
        if (normalized_dim < 0 || normalized_dim >= output.dim()) {
            throw std::invalid_argument("Dimension index out of range.");
        }

        // 1. Prepare for Batch Processing
    
        // Create permutation order: (Non-dim dims...) then (dim)
        std::vector<int64_t> permutation;
        for (int64_t d = 0; d < output.dim(); ++d) {
            if (d != normalized_dim) {
                permutation.push_back(d);
            }
        }
        permutation.push_back(normalized_dim);
    
        // Permute and make contiguous for a reliable view
        torch::Tensor permuted_output = output.permute(permutation).contiguous();
    
        // Determine dimensions for the (Batch Size, Sequence Length) view
        int64_t sequence_length = output.size(normalized_dim);
        int64_t total_batch_size = permuted_output.numel() / sequence_length;

        // View as (Batch Size, Sequence Length)
        torch::Tensor batch_view = permuted_output.view({total_batch_size, sequence_length});

        // 2. Iterate Over Slices (the 'Batch')
        for (int64_t i = 0; i < total_batch_size; ++i) {
            // Get the current 1D slice (editable view)
            torch::Tensor slice_output = batch_view[i];
        
            // Find the condition mask for the slice (size N)
            torch::Tensor condition = slice_output.gt(threshold); 

            // Pad the condition with false (0) on both ends (size N+2)
            // P = [0, C0, C1, ..., C_{N-1}, 0]
            torch::Tensor padded_condition = torch::cat({
                    torch::zeros(1, condition.options()),
                    condition.to(torch::kFloat),
                    torch::zeros(1, condition.options())
                });

            // Slice Tensors (all sliced to size N = sequence_length)
            // Previous: P[0..N] -> [0, C0, C1, ..., C_{N-1}]
            torch::Tensor previous = padded_condition.slice(0, 0, sequence_length);

            // Current: P[1..N+1] -> [C0, C1, ..., C_{N-1}, 0]
            torch::Tensor current = padded_condition.slice(0, 1, sequence_length + 1);

            // Next: P[2..N+2] -> [C1, C2, ..., C_{N-1}, 0] 
            // Note: The next state is used to detect the end of the run. We need the full N elements.
            torch::Tensor next = padded_condition.slice(0, 2, sequence_length + 2);


            // --- Start Mask Detection (Size N) ---
            // A fragment starts at index i when Previous[i]=0 AND Current[i]=1
            // This is equivalent to (Previous < Current)
            torch::Tensor starts_mask = (previous < current); 

            // --- End Mask Detection (Size N) ---
            // A fragment ends at index i when Current[i]=1 AND Next[i]=0
            // This is equivalent to (Current > Next)
            torch::Tensor ends_mask = (current > next);
        
            // Get the actual indices of starts and ends (relative to the original slice 0..N-1)
            torch::Tensor start_indices = torch::where(starts_mask)[0];
            torch::Tensor end_indices = torch::where(ends_mask)[0]; // This is the index of the last element > threshold

            // Ensure we have a matching number of starts and ends
            if (start_indices.numel() != end_indices.numel()) {
                // Should not happen with this corrected logic unless sequence wraps.
                continue; 
            }

            // 3. Apply Baseline Subtraction for all fragments in the slice
            for (int64_t j = 0; j < start_indices.numel(); ++j) {
                int64_t start_idx = start_indices[j].item<int64_t>();
                int64_t end_idx = end_indices[j].item<int64_t>();
            
                // Fragment must have size > 1 to define a unique linear line segment (Length > 0)
                if (end_idx <= start_idx) { 
                    continue; 
                } 

                // Get the first and last values of the fragment
                float start_val = slice_output[start_idx].item<float>();
                float end_val = slice_output[end_idx].item<float>();
            
                // Calculate the linear model (y = A*x + B)
                float length = static_cast<float>(end_idx - start_idx);
                float A = (end_val - start_val) / length;
                // Intercept B: B = start_val - A * start_idx
                float B = start_val - A * static_cast<float>(start_idx);

                // Create the index tensor (x) for the fragment region [start_idx, end_idx]
                torch::Tensor x = torch::arange(start_idx, end_idx + 1, slice_output.options());
            
                // Calculate the baseline vector (A*x + B) using vectorized operations
                torch::Tensor baseline = A * x + B;

                // Subtract the baseline from the original fragment region (vectorized index_put_)
                slice_output.index_put_({torch::indexing::Slice(start_idx, end_idx + 1)}, 
                                        slice_output.index({torch::indexing::Slice(start_idx, end_idx + 1)}) - baseline);
            }
        }
    
        // 4. Restore Original Shape and Permutation
        std::vector<int64_t> reverse_permutation(output.dim());
        for(int64_t k = 0; k < output.dim(); ++k) {
            reverse_permutation[permutation[k]] = k;
        }
    
        // Reshape back to the permuted shape, and then permute back to the original shape
        return batch_view.view(permuted_output.sizes()).permute(reverse_permutation).contiguous();
    }
}

