#include "WireCellSpng/TorchFunctionNode.h"


namespace WireCell::SPNG {


    WireCell::Configuration TorchFunctionNode::default_configuration() const
    {
        auto cfg1 = this->FunctionNode::default_configuration();
        auto cfg2 = this->ContextBase::default_configuration();
        return update(cfg1, cfg2);
    }

    void TorchFunctionNode::configure(const WireCell::Configuration& cfg)
    {
        this->FunctionNode::configure(cfg);
        this->ContextBase::configure(cfg);
    }


    TensorIndex TorchFunctionNode::sys_transform_tensors(TensorIndex ti) const
    {
        TorchSemaphore sem(context());
        torch::AutoGradMode enable_grad(false);
        return this->transform_tensors(std::move(ti));
    }

}

