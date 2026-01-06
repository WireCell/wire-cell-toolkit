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
    
        auto range_offsets = torch::cumsum(lengths, 0) - lengths;
    
        auto deltas = torch::ones({total_length}, ranges.options());
        deltas.index_put_({0}, starts[0]);

        // Fix the "jumps" at the start of every subsequent range
        if (ranges.size(0) > 1) {
            // Jump = Current Start - Previous End
            auto current_starts = starts.slice(0, 1);
            auto previous_ends = ends.slice(0, 0, -1);
            auto jump_values = current_starts - previous_ends;
            
            auto jump_indices = range_offsets.slice(0, 1);
            deltas.index_put_({jump_indices}, jump_values);
        }
    
        return torch::cumsum(deltas, 0);
    }
}
