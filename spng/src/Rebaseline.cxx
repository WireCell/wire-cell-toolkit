#include "WireCellSpng/Rebaseline.h"
#include "WireCellSpng/Util.h"  // for modulo

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
    torch::Tensor rebaseline( const torch::Tensor& tensor, 
                              int64_t min_roi_size,
                              int64_t expand_size,
                              bool remove_small,
                              bool remove_negative,
                              int64_t dim, float threshold) 
    {
        // Make a copy to modify to allow modify in place
        torch::Tensor output = tensor.clone();
    
        // Make the dimension to be non-negative.
        dim = modulo(dim, output.dim());

        // Move the dim so following can assume a layout
        std::vector<int64_t> permutation;
        for (int64_t d = 0; d < output.dim(); ++d) {
            if (d != dim) {
                permutation.push_back(d);
            }
        }
        permutation.push_back(dim);
        torch::Tensor permuted_output = output.permute(permutation).contiguous();
    
        // Get size of dimension along which ROIs are defined and size of all the rest.
        const int64_t nticks = output.size(dim);
        const int64_t nbatches = permuted_output.numel() / nticks;

        // View as (nbatch, nticks)
        torch::Tensor batch_view = permuted_output.view({nbatches, nticks});

        // 2. Iterate Over Slices (the 'Batch')
        for (int64_t i = 0; i < nbatches; ++i) {

            torch::Tensor wave = batch_view[i];
        
            // Next we do the pad-shift trick to find threshold crossings.  This
            // adds a zero on either side and then we compare less/greater than
            // for the first N and last N to the original.  Then find where
            // those comparisons are true and that gives us indices of
            // postive-going and negative-going threshold crossings.

            // Find the condition mask for the slice (size N)
            torch::Tensor condition = wave.gt(threshold); 

            // Pad the condition with false (0) on both ends (size N+2)
            // P = [0, C0, C1, ..., C_{N-1}, 0]
            torch::Tensor padded_condition = torch::cat({
                    torch::zeros(1, condition.options()),
                    condition.to(torch::kFloat),
                    torch::zeros(1, condition.options())
                });

            // Slice Tensors (all sliced to size N = nticks)
            // Previous: P[0..N] -> [0, C0, C1, ..., C_{N-1}]
            torch::Tensor previous = padded_condition.slice(0, 0, nticks);

            // Current: P[1..N+1] -> [C0, C1, ..., C_{N-1}, 0]
            torch::Tensor current = padded_condition.slice(0, 1, nticks + 1);

            // Next: P[2..N+2] -> [C1, C2, ..., C_{N-1}, 0] 
            // Note: The next state is used to detect the end of the run. We need the full N elements.
            torch::Tensor next = padded_condition.slice(0, 2, nticks + 2);


            torch::Tensor starts_mask = (previous < current); 
            torch::Tensor ends_mask = (current > next);
        
            // These are inclusive, closed intervals.
            torch::Tensor start_indices = torch::where(starts_mask)[0];
            torch::Tensor end_indices = torch::where(ends_mask)[0];

            // Ensure we have a matching number of starts and ends
            if (start_indices.numel() != end_indices.numel()) {
                // Should not happen with this corrected logic unless sequence wraps.
                continue; 
            }

            // Now a painful inner loop checking each ROI.
            for (int64_t j = 0; j < start_indices.numel(); ++j) {

                // These are the indices of the first/last above threshold
                int64_t start_idx = start_indices[j].item<int64_t>();
                int64_t end_idx = end_indices[j].item<int64_t>();
            
                // Handle "small" ROIs.  If equal, they represent an ROI of size
                // 1 so subtract one from the provided min size.
                if (end_idx - start_idx < min_roi_size-1) { 
                    if (remove_small) {
                        const auto roi = torch::indexing::Slice(start_idx, end_idx + 1);
                        wave.index_put_({roi}, 0.0);
                    }
                    continue;
                } 

                // Expand to include the first sample below threshold.  These
                // will form the anchors for the rebaseline.
                // Note: may want to make this expansion size be configurable.
                start_idx = std::max(start_idx-1, 0L);
                end_idx = std::min(end_idx+1, nticks-1);

                // Get the first and last values of the fragment
                float start_val = wave[start_idx].item<float>();
                float end_val = wave[end_idx].item<float>();
            
                const float length = end_idx - start_idx;
                const float slope = (end_val - start_val)/length;
                torch::Tensor baseline = start_val + slope * torch::arange(0, end_idx-start_idx+1,
                                                                           wave.options());

                const auto roi = torch::indexing::Slice(start_idx, end_idx + 1);
                // Subtract the baseline from the original fragment region (vectorized index_put_)
                wave.index_put_({roi}, wave.index({roi}) - baseline);
            }
        }
    
        if (remove_negative) {
            output.clamp_min_(0.0);
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

