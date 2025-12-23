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
     * @param min_roi_size The fragment must have this many samples above threshold.
     * @param expand_size The number of samples to expand.
     * @param dim The dimension along which to operate (default: last dimension, -1).
     * @param remove_small If true, set ROIs failing min size to 0.0 value.
     * @param remove_negative If true, after rebaseline, set all below-threshold samples to 0.0 value.
     * @param threshold The value above which a fragment is considered active.
     * @return The rebaselined tensor.
     */
    torch::Tensor rebaseline(const torch::Tensor& tensor,
                             int64_t min_roi_size = 1,
                             int64_t expand_size = 1,
                             bool remove_small = true,
                             bool remove_negative = true,
                             int64_t dim = -1, float threshold = 0.0f);
}

