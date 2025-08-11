#include "WireCellSpng/Convert.h"

using namespace WireCell;

void Torch::raster(torch::Tensor& block,
                   const WireCell::ITrace::vector& traces,
                   const std::vector<int>& channels)
{
    using torch::indexing::Slice;

    auto shape = block.sizes();

    // Use int64_t instead of size_t to mesh with torch int type.
    const int64_t nchannels = std::min((int64_t)channels.size(), shape[0]);
    std::unordered_map<int, size_t> ch2col;
    for (int64_t ind = 0; ind < nchannels; ++ind) {
        ch2col[channels[ind]] = ind;
    }

    for (auto trace : traces) {
        const int64_t tbin = trace->tbin();
        if (tbin >= shape[1]) {
            continue;
        }
        const auto& samples = trace->charge();

        // find the row to fill
        const int ch = trace->channel();
        const auto chit = ch2col.find(ch);
        if (chit == ch2col.end()) {
            continue;
        }
        const int64_t irow = chit->second;

        const int64_t nsamples = samples.size();
        const int64_t maxbin = std::min(nsamples, shape[1] - tbin);

        block.index_put_({irow, Slice(tbin, tbin+maxbin)},
                         // torch doesn't like const, but it's tmp so okay.
                         torch::from_blob((void*)samples.data(), maxbin)
                         + block.index({irow, Slice(tbin, tbin+maxbin)}));
    }
}


