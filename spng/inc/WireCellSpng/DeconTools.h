/// This API provides primitive functions for performing SP deconvolution.
#ifndef WIRECELL_SPNG_DECONTOOLS
#define WIRECELL_SPNG_DECONTOOLS

#include "WireCellSpng/Torch.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IWaveform.h"

#include <vector>

namespace WireCell::SPNG {

    /// Return the average interval space FR for the plane and return it as a
    /// flat valued 2D tensor on CPU.
    ///
    /// If zero_centered is true (default) the channel dimension of the tensor
    /// is rolled to bring the center row to the zero row position.  This
    /// removes the artificial shift that exists when the peak of the response
    /// is left in the center row.  Note, padding to match some desired shape
    /// for convolution should be inserted in the middle of the tensor.
    //
    // FIXME: this function is not well placed in this header.  Maybe a
    // "ResponseTools.h"?
    torch::Tensor fr_average_tensor(IFieldResponse::pointer ifr, int plane,
                                    bool zero_centered=true);


}

#endif
