/// This API provides primitive functions for performing SP deconvolution.
#ifndef WIRECELL_SPNG_DECONTOOLS
#define WIRECELL_SPNG_DECONTOOLS

#include "WireCellSpng/Torch.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IWaveform.h"

#include <vector>

namespace WireCell::SPNG {

    /**
       @brief Turn FR data into a tensor

       @param fr_data An FR schema object.
       @param plane_id The plane ID for which the tensor is formed.
       @param zero_centered True if the wire-of-interest is in row 0.
       @return A 2D tensor giving FR as (npath, nsteps).
     */
    torch::Tensor fr_tensor(const Response::Schema::FieldResponse& fr_data,
                            int plane_id,
                            bool zero_centered=true);

    /**
       Calculate the "wire-region averaged" FR and return it in tensor form.

       See fr_tensor() for details.
     */
    torch::Tensor fr_average_tensor(const Response::Schema::FieldResponse& fr_data,
                                    int plane,
                                    bool zero_centered=true);


}

#endif
