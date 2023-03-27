#ifndef WIRECELLAUX_SLICETOOLS
#define WIRECELLAUX_SLICETOOLS

#include "WireCellIface/ISlice.h"
#include "WireCellUtil/Binning.h"

#include <vector>

namespace WireCell::Aux {

    // Return a binning which is coherent with and includes all given
    // slices.  The time span of the first slice is used to determine
    // the bin size.
    Binning binning(const ISlice::vector& slices, double span=0);

    // Iterator version of above.  If span not given, span of first
    // slice is used.
    template<typename It>
    Binning binning(It beg, It end, double span=0) {
        if (beg == end) return Binning();
        if (span == 0) span = (*beg)->span();
        std::vector<double> starts;
        while (beg != end) {
            starts.push_back((*beg)->start());
            ++beg;
        }
        auto [f,l] = std::minmax_element(starts.begin(), starts.end());
        const int nbins = 1 + round((*l-*f)/span);
        return Binning(nbins, *f, *f + nbins*span);
    }
}

#endif
