#include "WireCellSpng/FunctionNode.h"
#include "WireCellSpng/SimpleTorchTensor.h"

namespace WireCell::SPNG {

    void FunctionNode::configure(const WireCell::Configuration& cfg)
    {
        this->ContextBase::configure(cfg);
        this->TensorSelector::configure(cfg);
    }
    
    WireCell::Configuration FunctionNode::default_configuration() const
    {
        auto tc_cfg = this->ContextBase::default_configuration();
        auto ts_cfg = this->TensorSelector::default_configuration();
        update(tc_cfg, ts_cfg);
        return tc_cfg;
    }



    bool FunctionNode::operator()(const input_pointer& in, output_pointer& out)
    {
        TorchSemaphore sem(context());
        torch::NoGradGuard no_grad;

        // Note, this selection and on-device bookkeeping could probably be
        // optimized to avoid some of the indexing.

        // Index to find parents
        TensorIndex ti(in);

        // Select matching parents
        std::vector<ITorchTensor::pointer> keep;
        for (auto child : ti.tree().child_values()) {
            auto result = select_tensor(child);
            if (result != TensorSelector::SelectionResult::kReject) keep.push_back(child);
        }
        TensorIndex sel_ti = ti.subset(keep);

        auto dev = device();

        // Assure tensors are on device
        std::vector<ITorchTensor::pointer> on_dev;
        for (auto ten : sel_ti.all_tensors()) {
            if (ten->device() == dev) {
                on_dev.push_back(ten);
                continue;
            }
            auto md = ten->metadata();
            auto tten = ten->tensor().to(dev);
            on_dev.push_back(std::make_shared<SimpleTorchTensor>(tten, md));
        }
            
        TensorIndex final_ti(on_dev);
        return (*this)(final_ti, out);
    }
    
    // calls virtual bool operator()(const TensorIndex& ti, output_pointer& out) = 0;

}
