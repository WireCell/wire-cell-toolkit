#include "WireCellAux/SliceTools.h"
#include <algorithm>
using namespace WireCell;

Binning Aux::binning(const ISlice::vector& slices, double span)
{
    return Aux::binning(slices.begin(), slices.end(), span);
}
