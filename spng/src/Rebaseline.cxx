#include "WireCellSpng/Rebaseline.h"
#include "WireCellSpng/Util.h"  // for modulo, to_string
#include "WireCellUtil/Exceptions.h"

#include <ATen/Parallel.h>

#include <algorithm>

namespace WireCell::SPNG {

    namespace {

        std::vector<int64_t> rebaseline_permutation(int64_t ndim, int64_t dim)
        {
            std::vector<int64_t> permutation;
            permutation.reserve(ndim);
            for (int64_t d = 0; d < ndim; ++d) {
                if (d != dim) {
                    permutation.push_back(d);
                }
            }
            permutation.push_back(dim);
            return permutation;
        }

        torch::Tensor restore_rebaseline_shape(const torch::Tensor& batch_view,
                                               const torch::Tensor& permuted_output,
                                               const std::vector<int64_t>& permutation)
        {
            std::vector<int64_t> reverse_permutation(permutation.size());
            for (int64_t k = 0; k < static_cast<int64_t>(permutation.size()); ++k) {
                reverse_permutation[permutation[k]] = k;
            }
            return batch_view.view(permuted_output.sizes()).permute(reverse_permutation).contiguous();
        }

        template<typename Wave>
        void apply_rebaseline_to_roi_like_original(Wave&& wave,
                                                   int64_t start_idx,
                                                   int64_t end_idx,
                                                   int64_t min_roi_size,
                                                   int64_t shrink_size,
                                                   bool remove_small)
        {
            start_idx += shrink_size;
            end_idx -= shrink_size;

            if (start_idx > end_idx) {
                return;
            }

            const auto roi = torch::indexing::Slice(start_idx, end_idx + 1);
            if (end_idx - start_idx <= min_roi_size - 1) {
                if (remove_small) {
                    wave.index_put_({roi}, 0.0);
                }
                return;
            }

            float start_val = wave[start_idx].template item<float>();
            float end_val = wave[end_idx].template item<float>();
            const float length = static_cast<float>(end_idx - start_idx);
            const float slope = (end_val - start_val) / length;
            auto baseline = start_val + slope * torch::arange(0, end_idx - start_idx + 1,
                                                              wave.options());
            wave.index_put_({roi}, wave.index({roi}) - baseline);
        }

        void rebaseline_zero_cuda(torch::Tensor& batch_view,
                                  int64_t nticks,
                                  int64_t nbatches,
                                  int64_t consequtive_zeros,
                                  int64_t min_roi_size,
                                  int64_t shrink_size,
                                  bool remove_small)
        {
            auto zero_mask = batch_view.eq(0).to(torch::kUInt8).to(torch::kCPU).contiguous();
            const auto* zeros = zero_mask.data_ptr<uint8_t>();

            for (int64_t i = 0; i < nbatches; ++i) {
                const auto* row_zeros = zeros + i * nticks;
                auto wave = batch_view[i];

                bool in_roi = false;
                int64_t roi_start = -1;

                int64_t k = 0;
                while (k < nticks) {
                    if (row_zeros[k]) {
                        const int64_t run_start = k;
                        while (k < nticks && row_zeros[k]) {
                            ++k;
                        }
                        if (k - run_start >= consequtive_zeros) {
                            if (in_roi) {
                                apply_rebaseline_to_roi_like_original(wave, roi_start, run_start - 1,
                                                                      min_roi_size, shrink_size, remove_small);
                                in_roi = false;
                            }
                        }
                        else if (!in_roi) {
                            roi_start = run_start;
                            in_roi = true;
                        }
                        continue;
                    }

                    if (!in_roi) {
                        roi_start = k;
                        in_roi = true;
                    }
                    ++k;
                }

                if (in_roi) {
                    apply_rebaseline_to_roi_like_original(wave, roi_start, nticks - 1,
                                                          min_roi_size, shrink_size, remove_small);
                }
            }

        }

        torch::Tensor rebaseline_zero_original_cpu(const torch::Tensor& tensor,
                                                   int64_t dim,
                                                   int64_t consequtive_zeros,
                                                   int64_t min_roi_size,
                                                   int64_t shrink_size,
                                                   bool remove_small,
                                                   bool remove_negative)
        {
            torch::Tensor output = tensor.clone();
            dim = modulo(dim, output.dim());

            std::vector<int64_t> permutation;
            for (int64_t d = 0; d < output.dim(); ++d) {
                if (d != dim) permutation.push_back(d);
            }
            permutation.push_back(dim);
            torch::Tensor permuted_output = output.permute(permutation).contiguous();

            const int64_t nticks = output.size(dim);
            const int64_t nbatches = permuted_output.numel() / nticks;
            torch::Tensor batch_view = permuted_output.view({nbatches, nticks});

            for (int64_t i = 0; i < nbatches; ++i) {
                torch::Tensor wave = batch_view[i];

                // Mark separator positions: runs of exactly-zero values of length >= consequtive_zeros.
                // These form the boundaries between ROIs (not part of any ROI).
                std::vector<bool> is_sep(nticks, false);
                {
                    int64_t k = 0;
                    while (k < nticks) {
                        if (wave[k].item<float>() == 0.0f) {
                            int64_t run_start = k;
                            while (k < nticks && wave[k].item<float>() == 0.0f) {
                                ++k;
                            }
                            if (k - run_start >= consequtive_zeros) {
                                for (int64_t m = run_start; m < k; ++m) {
                                    is_sep[m] = true;
                                }
                            }
                        } else {
                            ++k;
                        }
                    }
                }

                // Collect ROIs as inclusive (start, end) index pairs.
                // Non-separator positions form ROIs; tensor edges act as implicit boundaries.
                std::vector<std::pair<int64_t,int64_t>> rois;
                bool in_roi = false;
                int64_t roi_start = -1;
                for (int64_t k = 0; k < nticks; ++k) {
                    if (!is_sep[k]) {
                        if (!in_roi) { roi_start = k; in_roi = true; }
                    } else {
                        if (in_roi) { rois.push_back({roi_start, k-1}); in_roi = false; }
                    }
                }
                if (in_roi) { rois.push_back({roi_start, nticks-1}); }

                for (auto [start_idx, end_idx] : rois) {
                    // Shrink the ROI inward from both ends.
                    start_idx += shrink_size;
                    end_idx   -= shrink_size;

                    // Skip ROIs that collapsed under shrinking.
                    if (start_idx > end_idx) continue;

                    // Handle small ROIs (same logic as rebaseline()).
                    if (end_idx - start_idx <= min_roi_size - 1) {
                        if (remove_small) {
                            const auto roi = torch::indexing::Slice(start_idx, end_idx + 1);
                            wave.index_put_({roi}, 0.0);
                        }
                        continue;
                    }

                    // Build a linear model anchored at the first/last sample of the ROI
                    // and subtract it from the ROI region.
                    float start_val = wave[start_idx].item<float>();
                    float end_val   = wave[end_idx].item<float>();
                    const float length = static_cast<float>(end_idx - start_idx);
                    const float slope  = (end_val - start_val) / length;
                    torch::Tensor baseline = start_val + slope * torch::arange(0, end_idx - start_idx + 1,
                                                                                wave.options());
                    const auto roi = torch::indexing::Slice(start_idx, end_idx + 1);
                    wave.index_put_({roi}, wave.index({roi}) - baseline);
                }
            }

            if (remove_negative) {
                output.clamp_min_(0.0);
            }

            std::vector<int64_t> reverse_permutation(output.dim());
            for (int64_t k = 0; k < output.dim(); ++k) {
                reverse_permutation[permutation[k]] = k;
            }
            return batch_view.view(permuted_output.sizes()).permute(reverse_permutation).contiguous();
        }

        // CPU-specialized implementation.  batch_view is contiguous float32
        // storage with one waveform per row.  The rows have no shared writes,
        // so their ROI scans and baseline corrections may run concurrently.
        void apply_rebaseline_to_cpu_row(float* wave,
                                         int64_t start_idx,
                                         int64_t end_idx,
                                         int64_t min_roi_size,
                                         int64_t shrink_size,
                                         bool remove_small)
        {
            start_idx += shrink_size;
            end_idx -= shrink_size;

            if (start_idx > end_idx) {
                return;
            }
            if (end_idx - start_idx <= min_roi_size - 1) {
                if (remove_small) {
                    std::fill(wave + start_idx, wave + end_idx + 1, 0.0f);
                }
                return;
            }

            const float start_val = wave[start_idx];
            const float end_val = wave[end_idx];
            const float slope = (end_val - start_val) /
                                static_cast<float>(end_idx - start_idx);
            for (int64_t tick = start_idx; tick <= end_idx; ++tick) {
                wave[tick] -= start_val + slope * static_cast<float>(tick - start_idx);
            }
        }

        void rebaseline_zero_cpu_rows(float* samples,
                                      int64_t first_row,
                                      int64_t last_row,
                                      int64_t nticks,
                                      int64_t consequtive_zeros,
                                      int64_t min_roi_size,
                                      int64_t shrink_size,
                                      bool remove_small)
        {
            for (int64_t row = first_row; row < last_row; ++row) {
                float* wave = samples + row * nticks;
                bool in_roi = false;
                int64_t roi_start = -1;

                int64_t tick = 0;
                while (tick < nticks) {
                    if (wave[tick] == 0.0f) {
                        const int64_t run_start = tick;
                        while (tick < nticks && wave[tick] == 0.0f) {
                            ++tick;
                        }
                        if (tick - run_start >= consequtive_zeros) {
                            if (in_roi) {
                                apply_rebaseline_to_cpu_row(wave, roi_start, run_start - 1,
                                                            min_roi_size, shrink_size, remove_small);
                                in_roi = false;
                            }
                        }
                        else if (!in_roi) {
                            roi_start = run_start;
                            in_roi = true;
                        }
                        continue;
                    }

                    if (!in_roi) {
                        roi_start = tick;
                        in_roi = true;
                    }
                    ++tick;
                }

                if (in_roi) {
                    apply_rebaseline_to_cpu_row(wave, roi_start, nticks - 1,
                                                min_roi_size, shrink_size, remove_small);
                }
            }
        }

        torch::Tensor rebaseline_zero_cpu_direct(const torch::Tensor& tensor,
                                                    int64_t dim,
                                                    int64_t consequtive_zeros,
                                                    int64_t min_roi_size,
                                                    int64_t shrink_size,
                                                    bool remove_small,
                                                    bool remove_negative)
        {
            if (!tensor.numel()) {
                return tensor.clone();
            }

            dim = modulo(dim, tensor.dim());
            const auto permutation = rebaseline_permutation(tensor.dim(), dim);
            auto output = tensor.clone();
            auto permuted_output = output.permute(permutation).contiguous();

            const int64_t nticks = output.size(dim);
            if (!nticks) {
                return restore_rebaseline_shape(permuted_output, permuted_output, permutation);
            }
            const int64_t nbatches = permuted_output.numel() / nticks;
            auto batch_view = permuted_output.view({nbatches, nticks});
            auto* samples = batch_view.data_ptr<float>();

            // Process all independent waveform rows on the calling thread.
            // This retains the direct-memory ROI algorithm without creating
            // a second, uncoordinated CPU worker pool.
            rebaseline_zero_cpu_rows(samples, 0, nbatches,
                                     nticks, consequtive_zeros, min_roi_size,
                                     shrink_size, remove_small);

            if (remove_negative) {
                output.clamp_min_(0.0);
            }
            return restore_rebaseline_shape(batch_view, permuted_output, permutation);
        }

    }

    torch::Tensor rebaseline( const torch::Tensor& tensor, 
                              int64_t dim, float threshold,
                              int64_t min_roi_size,
                              int64_t expand_size,
                              bool remove_small,
                              bool remove_negative) 
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
                if (end_idx - start_idx <= min_roi_size-1) { 
                    if (remove_small) {
                        const auto roi = torch::indexing::Slice(start_idx, end_idx + 1);
                        wave.index_put_({roi}, 0.0);
                    }
                    continue;
                } 

                // Expand to include the first sample below threshold.  These
                // will form the anchors for the rebaseline.
                // Note: may want to make this expansion size be configurable.
                if (expand_size) {
                    start_idx = std::max(start_idx-expand_size, 0L);
                    end_idx = std::min(end_idx+1, nticks-expand_size);
                }

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


    torch::Tensor rebaseline_zero(const torch::Tensor& tensor,
                                  int64_t dim,
                                  int64_t consequtive_zeros,
                                  int64_t min_roi_size,
                                  int64_t shrink_size,
                                  bool remove_small,
                                  bool remove_negative)
    {
        // Float tensors use the direct contiguous-storage CPU implementation.
        // Other dtypes retain the reference path, which uses item<float>().
        if (!tensor.is_cuda()) {
            if (tensor.scalar_type() == torch::kFloat) {
                return rebaseline_zero_cpu_direct(tensor, dim, consequtive_zeros,
                                                    min_roi_size, shrink_size,
                                                    remove_small, remove_negative);
            }
            return rebaseline_zero_original_cpu(tensor, dim, consequtive_zeros,
                                                min_roi_size, shrink_size,
                                                remove_small, remove_negative);
        }

        if (!tensor.numel()) {
            return tensor.clone();
        }

        dim = modulo(dim, tensor.dim());
        auto permutation = rebaseline_permutation(tensor.dim(), dim);

        // Keep the same storage relationship as the original CPU code:
        // output is cloned first, then permuted_output is made from output.
        // This matters for exact compatibility of remove_negative.  If the
        // permutation is effectively contiguous, output.clamp_min_() can affect
        // the returned view.  If permute(...).contiguous() creates independent
        // storage, it will not.  That is the original behavior.
        torch::Tensor output = tensor.clone();
        torch::Tensor permuted_output = output.permute(permutation).contiguous();

        const int64_t nticks = output.size(dim);
        if (nticks == 0) {
            return restore_rebaseline_shape(permuted_output, permuted_output, permutation);
        }

        const int64_t nbatches = permuted_output.numel() / nticks;
        torch::Tensor batch_view = permuted_output.view({nbatches, nticks});

        rebaseline_zero_cuda(batch_view, nticks, nbatches, consequtive_zeros,
                             min_roi_size, shrink_size, remove_small);

        // Preserve original CPU semantics.  Do not clamp batch_view here.
        if (remove_negative) {
            output.clamp_min_(0.0);
        }

        return restore_rebaseline_shape(batch_view, permuted_output, permutation);
    }


    // --- rebaseline_masked -------------------------------------------------

    namespace {

        // Dilate a (B,N) bool mask by pad samples on either side along dim 1.
        torch::Tensor dilate_mask(const torch::Tensor& mask, int64_t pad)
        {
            if (pad <= 0) {
                return mask;
            }
            auto f = mask.to(torch::kFloat32).unsqueeze(1);           // (B,1,N)
            auto d = torch::max_pool1d(f, /*kernel*/ 2*pad + 1, /*stride*/ 1, /*padding*/ pad);
            return d.squeeze(1) > 0.5f;
        }

        // Bring (wave, mask) to a (B,N) float / bool pair with dim last.
        struct MaskedView {
            torch::Tensor wave;       // (B,N) floating point, contiguous
            torch::Tensor mask;       // (B,N) bool, contiguous
            torch::Tensor permuted;   // wave permuted, for shape restoration
            std::vector<int64_t> permutation;
        };

        MaskedView masked_view(const torch::Tensor& wave, const torch::Tensor& mask, int64_t dim)
        {
            if (wave.sizes() != mask.sizes()) {
                raise<ValueError>("rebaseline_masked: waveform %s and mask %s differ in shape",
                                  to_string(wave), to_string(mask));
            }
            MaskedView mv;
            dim = modulo(dim, wave.dim());
            mv.permutation = rebaseline_permutation(wave.dim(), dim);
            auto w = wave.is_floating_point() ? wave : wave.to(torch::kFloat32);
            mv.permuted = w.permute(mv.permutation).contiguous();
            const int64_t nticks = wave.size(dim);
            const int64_t nbatches = nticks ? mv.permuted.numel() / nticks : 0;
            mv.wave = mv.permuted.view({nbatches, nticks});
            auto m = mask.scalar_type() == torch::kBool ? mask : mask.ne(0);
            mv.mask = m.permute(mv.permutation).contiguous().view({nbatches, nticks});
            return mv;
        }

        // The algorithm proper on a (B,N) float waveform and bool mask.
        torch::Tensor rebaseline_masked_core(const torch::Tensor& W,
                                             const torch::Tensor& M,
                                             const MaskRebaselineOptions& opt)
        {
            const int64_t B = W.size(0);
            const int64_t N = W.size(1);
            auto Mout = M;
            if (!B || !N) {
                return torch::where(M, W, torch::zeros_like(W));
            }

            // ROI geometry is computed on the dilated mask so the anchors are
            // the first/last samples of each dilated ROI, i.e. anchor_pad
            // samples outside the tight ROI.
            auto Md = dilate_mask(M, opt.anchor_pad);
            if (opt.keep_pad) {
                Mout = Md;
            }

            auto zcol = torch::zeros({B, 1}, M.options());
            auto Mprev = torch::cat({zcol, Md.slice(1, 0, N - 1)}, 1);
            auto Mnext = torch::cat({Md.slice(1, 1, N), zcol}, 1);
            auto starts = Md & ~Mprev;         // first sample of an ROI
            auto ends = Md & ~Mnext;           // last sample of an ROI
            auto gap_begin = ~Md & Mprev;      // first sample after an ROI
            auto gap_end = ~Md & Mnext;        // last sample before an ROI

            // Per-sample ROI start/end and the bounds of the gap around it, via
            // forward/backward running max/min scans.  int32 keeps the index
            // tensors small; gather() below is the only place needing int64.
            auto iopts = W.options().dtype(torch::kInt32);
            auto idx = torch::arange(N, iopts).unsqueeze(0).expand({B, N});
            auto zero = torch::zeros_like(idx);
            auto last = torch::full_like(idx, static_cast<int>(N - 1));
            auto start_of = std::get<0>(torch::cummax(torch::where(starts, idx, zero), 1));
            auto end_of = std::get<0>(torch::cummin(torch::where(ends, idx, last).flip(1), 1)).flip(1);
            auto lo_floor = std::get<0>(torch::cummax(torch::where(gap_begin, idx, zero), 1));   // prev end + 1
            auto hi_ceil = std::get<0>(torch::cummin(torch::where(gap_end, idx, last).flip(1), 1)).flip(1); // next start - 1

            // Anchor windows: [lo_a, lo_b] ending at the ROI start and
            // [hi_a, hi_b] beginning at the ROI end, clipped to the gaps.
            const int w1 = static_cast<int>(std::max<int64_t>(opt.anchor_width, 1) - 1);
            auto lo_b = start_of;
            auto lo_a = torch::maximum(start_of - w1, lo_floor);
            auto hi_a = end_of;
            auto hi_b = torch::minimum(end_of + w1, hi_ceil);

            // Window means from a float64 prefix sum (float32 differences of
            // large nearby prefix values would lose precision).
            auto dopts = W.options().dtype(torch::kFloat64);
            auto cs = torch::cat({torch::zeros({B, 1}, dopts), W.cumsum(1, torch::kFloat64)}, 1);
            auto lo_sum = cs.gather(1, (lo_b + 1).to(torch::kInt64)) - cs.gather(1, lo_a.to(torch::kInt64));
            auto hi_sum = cs.gather(1, (hi_b + 1).to(torch::kInt64)) - cs.gather(1, hi_a.to(torch::kInt64));
            auto a = (lo_sum / (lo_b - lo_a + 1).to(torch::kFloat64)).to(W.scalar_type());
            auto b = (hi_sum / (hi_b - hi_a + 1).to(torch::kFloat64)).to(W.scalar_type());

            // Line through (start, a) and (end, b).
            auto len = end_of - start_of;
            auto frac = (idx - start_of).to(W.scalar_type()) / torch::clamp_min(len, 1).to(W.scalar_type());
            auto baseline = a + (b - a) * frac;
            if (opt.min_roi_size > 0) {
                baseline.masked_fill_(len <= static_cast<int>(opt.min_roi_size - 1), 0.0);
            }

            auto out = torch::where(Mout, W - baseline, torch::zeros_like(W));
            if (opt.remove_negative) {
                out.clamp_min_(0.0);
            }
            return out;
        }

        // Direct CPU implementation with the same semantics as
        // rebaseline_masked_core(): one scan per row over raw storage,
        // rows in parallel.  Much faster than the torch version on CPU.
        torch::Tensor rebaseline_masked_core_cpu(const torch::Tensor& W,
                                                 const torch::Tensor& M,
                                                 const MaskRebaselineOptions& opt)
        {
            const int64_t B = W.size(0);
            const int64_t N = W.size(1);
            auto Wc = W.contiguous();
            auto Mc = M.to(torch::kUInt8).contiguous();
            auto out = torch::zeros_like(Wc);
            if (!B || !N) {
                return out;
            }
            const float* wp = Wc.data_ptr<float>();
            const uint8_t* mp = Mc.data_ptr<uint8_t>();
            float* op = out.data_ptr<float>();
            const int64_t pad = std::max<int64_t>(opt.anchor_pad, 0);
            const int64_t width = std::max<int64_t>(opt.anchor_width, 1);

            at::parallel_for(0, B, 1, [&](int64_t r0, int64_t r1) {
                std::vector<std::pair<int64_t,int64_t>> runs;   // dilated ROIs [s,e]
                for (int64_t r = r0; r < r1; ++r) {
                    const float* w = wp + r * N;
                    const uint8_t* m = mp + r * N;
                    float* o = op + r * N;

                    runs.clear();
                    int64_t t = 0;
                    while (t < N) {
                        if (!m[t]) { ++t; continue; }
                        const int64_t s = t;
                        while (t < N && m[t]) ++t;
                        const int64_t e = t - 1;
                        const int64_t sd = std::max<int64_t>(s - pad, 0);
                        const int64_t ed = std::min<int64_t>(e + pad, N - 1);
                        if (!runs.empty() && sd <= runs.back().second + 1) {
                            runs.back().second = ed;      // dilated ROIs touch: merge
                        }
                        else {
                            runs.emplace_back(sd, ed);
                        }
                    }

                    const size_t nr = runs.size();
                    for (size_t k = 0; k < nr; ++k) {
                        const int64_t sd = runs[k].first;
                        const int64_t ed = runs[k].second;
                        const int64_t lo_floor = k ? runs[k-1].second + 1 : 0;
                        const int64_t hi_ceil = (k + 1 < nr) ? runs[k+1].first - 1 : N - 1;
                        const int64_t lo_a = std::max(sd - (width - 1), lo_floor);
                        const int64_t hi_b = std::min(ed + (width - 1), hi_ceil);
                        double a = 0.0;
                        for (int64_t i = lo_a; i <= sd; ++i) a += w[i];
                        a /= static_cast<double>(sd - lo_a + 1);
                        double b = 0.0;
                        for (int64_t i = ed; i <= hi_b; ++i) b += w[i];
                        b /= static_cast<double>(hi_b - ed + 1);

                        const int64_t len = ed - sd;
                        const bool flat = opt.min_roi_size > 0 && len <= opt.min_roi_size - 1;
                        const double a0 = flat ? 0.0 : a;
                        const double slope = flat ? 0.0 : (b - a) / static_cast<double>(std::max<int64_t>(len, 1));
                        for (int64_t i = sd; i <= ed; ++i) {
                            if (opt.keep_pad || m[i]) {
                                o[i] = static_cast<float>(w[i] - (a0 + slope * static_cast<double>(i - sd)));
                            }
                        }
                    }
                }
            });
            if (opt.remove_negative) {
                out.clamp_min_(0.0);
            }
            return out;
        }

        torch::Tensor rebaseline_masked_dispatch(const torch::Tensor& W,
                                                 const torch::Tensor& M,
                                                 const MaskRebaselineOptions& opt)
        {
            if (!opt.force_vectorized && !W.is_cuda() && W.scalar_type() == torch::kFloat) {
                return rebaseline_masked_core_cpu(W, M, opt);
            }
            return rebaseline_masked_core(W, M, opt);
        }
    }


    torch::Tensor rebaseline_masked(const torch::Tensor& wave,
                                    const torch::Tensor& mask,
                                    const MaskRebaselineOptions& opt)
    {
        auto mv = masked_view(wave, mask, opt.dim);
        auto out = rebaseline_masked_dispatch(mv.wave, mv.mask, opt);
        return restore_rebaseline_shape(out, mv.permuted, mv.permutation);
    }


    torch::Tensor rebaseline_masked_rows(const torch::Tensor& wave,
                                         const torch::Tensor& mask,
                                         const std::vector<std::pair<int64_t,int64_t>>& rows,
                                         const MaskRebaselineOptions& opt)
    {
        if (rows.empty()) {
            return rebaseline_masked(wave, mask, opt);
        }
        auto mv = masked_view(wave, mask, opt.dim);
        const int64_t B = mv.wave.size(0);
        auto out = torch::where(mv.mask, mv.wave, torch::zeros_like(mv.wave));
        for (const auto& [r0, r1] : rows) {
            if (r0 < 0 || r1 < r0 || r1 > B) {
                raise<ValueError>("rebaseline_masked_rows: bad row range [%d,%d) for %d rows",
                                  r0, r1, B);
            }
            if (r1 == r0) {
                continue;
            }
            out.slice(0, r0, r1).copy_(
                rebaseline_masked_dispatch(mv.wave.slice(0, r0, r1), mv.mask.slice(0, r0, r1), opt));
        }
        return restore_rebaseline_shape(out, mv.permuted, mv.permutation);
    }

}
