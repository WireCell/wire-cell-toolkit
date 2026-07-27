/// Functions to handle ragged arrays.
///
/// Implementations attempt to use optimized/vectorized forms.
#pragma once

#include "WireCellSpng/Torch.h"

namespace WireCell::SPNG::Ragged {

    /// @brief Return lengths of ranges
    ///
    /// @param ranges A tensor of shape (N,2) giving N half-open range boundaries.
    /// @return A tensor of shape (N,) holding the length of each range.
    torch::Tensor range_lengths(torch::Tensor ranges);

    /// @brief Expand the index of each range in ranges.
    ///
    /// @param ranges A tensor of shape (N,2) giving N half-open range boundaries.
    /// @return A tensor of shape (Ntot,) holding the concatenation of the indices in N over the expansion.
    ///
    /// Eg, given half-open range boundaries [1, 3) and [10, 13)
    /// return: [0, 0, 1, 1, 1]
    torch::Tensor range_index_expansion(torch::Tensor ranges);

    /// @brief Expand the values of each range in ranges.
    ///
    /// @param ranges A tensor of shape (N,2) giving N half-open range boundaries.
    /// @return A tensor of shape (Ntot,) holding the concatenation of the expansion of ranges.
    ///
    /// Eg, given half-open range boundaries [1, 3) and [10, 13)
    /// return: [1, 2, 10, 11, 12]
    torch::Tensor range_value_expansion(torch::Tensor ranges);
}

