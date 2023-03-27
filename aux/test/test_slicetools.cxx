#include "WireCellUtil/Testing.h"
#include "WireCellAux/SliceTools.h"

#include <vector>
#include <iostream>

using namespace WireCell::Aux;

struct NotASlice
{
    double _start, _span;
    double start() const { return _start; }
    double span() const { return _span; }
};

int main()
{
    std::vector<NotASlice*> slices = {
        new NotASlice{1,0.5},
        new NotASlice{3,0.5},
        new NotASlice{4,0.5},
        new NotASlice{5,0.5},
    };
    
    auto b1 = binning(slices.begin(), slices.end());
    std::cerr << "binning: " << b1 << " binsize=" << b1.binsize() << "\n";
    Assert(b1.binsize() == 0.5);
    Assert(b1.nbins() == 9);
    Assert(b1.min() == 1);
    Assert(b1.max() == 5.5);

    for (auto s : slices) {
        delete s;
    }
    return 0;
}
