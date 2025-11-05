#pragma once

#include "WireCellSpng/FanBase.h"
#include "WireCellSpng/ITorchTensorSet.h"


namespace WireCell::SPNG {

    /** @brief Tensors go in, tensor set comes out. */
    class TensorPacker : public FaninBase<ITorchTensor, ITorchTensorSet> 
    {
    public:
        TensorPacker();
        virtual ~TensorPacker() = default;

        virtual void fanin_combine(const input_vector& inv, output_pointer& out);

    };

}
