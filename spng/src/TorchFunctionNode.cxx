#include "WireCellSpng/TorchFunctionNode.h"
#include "WireCellSpng/SimpleTorchTensor.h"


namespace WireCell::SPNG {

    TorchFunctionNode::TorchFunctionNode(const std::string& logname, const std::string& pkgname)
        : FunctionNode(logname, pkgname)
    {
    }

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

        // Assure tensors are on the configured desired.  
        auto dev = device();
        TensorIndex ti_dev;

        for (const auto& node : ti.tree().depth()) {
            if (! node.value) { // skip empty root node
                continue;
            }
            auto ten = node.value;
            if (dev == ten->device()) {
                ti_dev.add(ten);
                continue;
            }
            ti_dev.add(std::make_shared<SimpleTorchTensor>(ten->tensor().to(dev), ten->metadata()));
        }
        return this->transform_tensors(std::move(ti_dev));
    }

}

