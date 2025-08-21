/**
   Functions to convert to and from torch tensors.

   In future, this may be moved into the "pytorch/" sub-package.
 */

#ifndef WIRECELLTORCHCONVERT
#define WIRECELLTORCHCONVERT

#include "WireCellIface/ITrace.h"
#include "WireCellSpng/Torch.h"
#include <vector>

namespace WireCell::Torch {

    // Raster a 2D tensor with values from traces.
    //
    // Each row of the tensor corresponds, in order, to an element in
    // "channels".  The trace's tbin is used as an offset from column 0.  If the
    // array block is undersized, missed samples will be quietly ignored.  Trace
    // charge is added to any existing values in the array block.
    //
    // See also Aux::raster() which does similar to fill an Eigen array.

    void raster(torch::Tensor& block,
                const WireCell::ITrace::vector& traces, const std::vector<int>& channels);

}

#endif
