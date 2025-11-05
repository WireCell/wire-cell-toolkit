#pragma once

#include "WireCellSpng/FanBase.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellIface/IConfigurable.h"


namespace WireCell::SPNG {

    // Simply multiply all inputs to the fan, output the product.
    struct Multiply : FaninBase<ITorchTensor>
    {
        Multiply();
        virtual ~Multiply() = default;

        // FaninBase API.
        virtual void fanin_combine(const input_vector& inv, output_pointer& out);

    };

}

