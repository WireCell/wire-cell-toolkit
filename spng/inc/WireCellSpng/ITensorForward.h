/** An interface to a "forward" operator on a torch tensor */

#pragma once

#include "WireCellUtil/IComponent.h"
#include "WireCellSpng/ITorchTensor.h"

namespace WireCell::SPNG {

    class ITensorForward : public IComponent<ITensorForward> {
    public:
        virtual ~ITensorForward() = default;

        virtual ITorchTensor::pointer forward(const ITorchTensor::pointer& input) const = 0;
    };

}
