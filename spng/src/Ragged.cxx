#include "WireCellSpng/Ragged.h"

namespace WireCell::SPNG::Ragged {

    torch::Tensor range_lengths(torch::Tensor ranges)
    {
        auto starts = ranges.select(1, 0);
        auto ends = ranges.select(1, 1);
    
        return ends - starts;
    }

    torch::Tensor range_index_expansion(torch::Tensor ranges)
    {
        auto lengths = range_lengths(ranges);
        return torch::repeat_interleave(
            torch::arange(ranges.size(0), ranges.options()), 
            lengths);
    }

    torch::Tensor range_value_expansion(torch::Tensor ranges)
    {
        auto starts = ranges.select(1, 0);
        auto ends = ranges.select(1, 1);
    
        auto lengths = ends - starts;
        auto total_length = lengths.sum().item<int64_t>();
    
        if (total_length == 0) {
            return torch::empty({0}, ranges.options());
        }

        auto range_offsets = torch::cumsum(lengths, 0) - lengths;
    
        auto deltas = torch::ones({total_length}, ranges.options());
        deltas.index_put_({0}, starts[0]);

        // Fix the "jumps" at the start of every subsequent range
        if (ranges.size(0) > 1) {
            // Jump = Current Start - Previous End
            auto current_starts = starts.slice(0, 1);
            auto previous_ends = ends.slice(0, 0, -1);
            auto jump_values = current_starts - previous_ends;
            
            // The required delta is jump_values + 1, because V[i] = V[i-1] + delta[i].
            // V[i] = starts[k], V[i-1] = ends[k-1] - 1.
            // delta[i] = starts[k] - (ends[k-1] - 1) = (starts[k] - ends[k-1]) + 1
            auto required_deltas = jump_values + 1;
            
            auto jump_indices_all = range_offsets.slice(0, 1);

            // Filter indices that are within bounds [0, total_length - 1]
            // Indices equal to total_length correspond to the start of an empty range 
            // that follows the last element, and must be excluded.
            auto total_length_tensor = torch::tensor({total_length}, ranges.options());
            auto mask = jump_indices_all < total_length_tensor;
            
            auto jump_indices = jump_indices_all.masked_select(mask);
            auto filtered_required_deltas = required_deltas.masked_select(mask);

            if (jump_indices.size(0) > 0) {
                deltas.index_put_({jump_indices}, filtered_required_deltas);
            }
        }
    
        return torch::cumsum(deltas, 0);
    }
}
