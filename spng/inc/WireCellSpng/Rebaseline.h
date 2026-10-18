// Util functions to "rebaseline" a tensor consisting or non-zero ROIs.  See Rebaseliner for a DFP node.
#pragma once

#include "WireCellSpng/Torch.h"

namespace WireCell::SPNG {

    /**
     * @brief Performs local baseline subtraction on 1D fragments > threshold.
     *
     * The function identifies contiguous 1D fragments (ROIs) along the
     * specified dimension (dim) where values are greater than the threshold. It
     * then subtracts a linear baseline model connecting the fragment's first
     * and last sample values.
     *
     * @param tensor The input tensor.
     * @param dim The dimension along which to operate (default: last dimension, -1).
     * @param threshold The value above which a fragment is considered active.
     * @param min_roi_size The fragment must have this many samples above threshold.
     * @param expand_size Fixed number of samples to expand an initial ROI on either end.
     * @param remove_small If true, set ROIs failing min size to 0.0 value.
     * @param remove_negative If true, after rebaseline, set all below-threshold samples to 0.0 value.
     * @return The rebaselined tensor.
     */
    torch::Tensor rebaseline(const torch::Tensor& tensor,
                             int64_t dim = -1,
                             float threshold = 0.0f,
                             int64_t min_roi_size = 1,
                             int64_t expand_size = 0,
                             bool remove_small = false,
                             bool remove_negative = false);


    torch::Tensor rebaseline_zero(const torch::Tensor& tensor,
                                  int64_t dim = -1,
                                  int64_t consequtive_zeros = 2,
                                  int64_t min_roi_size = 1,
                                  int64_t shrink_size = 0,
                                  bool remove_small = false,
                                  bool remove_negative = false);


    /// Options for rebaseline_masked().
    struct MaskRebaselineOptions {
        /// The dimension along which ROIs extend (ticks).
        int64_t dim = -1;

        /// How many samples outside the (tight) ROI the baseline anchors sit.
        /// 0 reproduces rebaseline_zero() on the masked waveform: the anchors
        /// are the ROI's own first and last samples.
        int64_t anchor_pad = 0;

        /// Number of samples averaged for each anchor value, extending away
        /// from the ROI starting at the anchor position.  The window is
        /// clipped so it never reaches into a neighboring ROI (or past the
        /// tensor edge).
        int64_t anchor_width = 1;

        /// ROIs with (end - start) <= min_roi_size - 1 are left as waveform*mask.
        int64_t min_roi_size = 1;

        /// Output region.  false: the tight mask (default).  true: the mask
        /// dilated by anchor_pad, so the pad samples are kept (rebaselined).
        bool keep_pad = false;

        /// If true, clamp the output to be non-negative.
        bool remove_negative = false;

        /// Testing aid: use the whole-tensor torch implementation even on
        /// CPU, where a direct per-row scan is used by default (same
        /// results, ~30x faster).  CUDA always uses the torch implementation.
        bool force_vectorized = false;
    };

    /**
     * @brief Mask-defined local baseline subtraction with anchors outside the ROI.
     *
     * An ROI is a maximal run of true mask samples along dim.  The mask is
     * dilated along dim by anchor_pad samples (ROIs closer than 2*anchor_pad
     * merge).  For each dilated ROI [s,e] two anchor values are formed from
     * the *unmasked* waveform: a = mean(wave[s-anchor_width+1 .. s]) and
     * b = mean(wave[e .. e+anchor_width-1]), each window clipped to the gap
     * between neighboring dilated ROIs and to the tensor.  The line through
     * (s,a) and (e,b) is subtracted and the result multiplied by the output
     * mask (tight or dilated, see keep_pad).
     *
     * Everything is expressed as whole-tensor torch operations: there is no
     * per-ROI loop and no host/device synchronization, so the same code runs
     * on CPU and CUDA.
     *
     * @param wave The waveform tensor (floating point; others are converted to float32).
     * @param mask Same shape as wave.  Bool, or any dtype where nonzero means "in ROI".
     */
    torch::Tensor rebaseline_masked(const torch::Tensor& wave,
                                    const torch::Tensor& mask,
                                    const MaskRebaselineOptions& opt = MaskRebaselineOptions());

    /**
     * @brief As rebaseline_masked() but only for selected rows.
     *
     * Rows are [begin,end) ranges of the index that remains after
     * flattening all dimensions other than dim (for a 2D (nchan,ntick)
     * tensor with dim=-1 this is the channel index).  Rows not covered are
     * returned as waveform*mask.  An empty list selects all rows.  Throws
     * ValueError on an invalid range.
     */
    torch::Tensor rebaseline_masked_rows(const torch::Tensor& wave,
                                         const torch::Tensor& mask,
                                         const std::vector<std::pair<int64_t,int64_t>>& rows,
                                         const MaskRebaselineOptions& opt = MaskRebaselineOptions());

}

