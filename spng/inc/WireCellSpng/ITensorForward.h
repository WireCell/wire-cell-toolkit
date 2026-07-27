#pragma once

#include "WireCellSpng/Torch.h"
#include "WireCellUtil/IComponent.h"

namespace WireCell::SPNG {

    /** An interface a BARE torch tensor transformation function.
     *
     * Unlike most "service" type components, this interface is defined with
     * concrete types and not an IData interface.  It's primary intention is to
     * type-erase the specifics of how inference is done while allowing sharing
     * of the model.  See TensorForwardTS for a wrapper around a TorchScript
     * model.
     */

    class ITensorForward : public IComponent<ITensorForward> {
    public:
        virtual ~ITensorForward() = default;

        /// This is a BARE torch tensor, not an IData.
        virtual torch::Tensor forward(const torch::Tensor& input) const = 0;
    };

}
