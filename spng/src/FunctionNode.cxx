#include "WireCellSpng/FunctionNode.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"

namespace WireCell::SPNG {

    void FunctionNode::configure(const WireCell::Configuration& cfg)
    {
        m_selector.configure(cfg);
        m_renaming.configure(cfg);
    }
    
    WireCell::Configuration FunctionNode::default_configuration() const
    {
        auto cfg1 = m_selector.default_configuration();
        auto cfg2 = m_renaming.default_configuration();
        return update(cfg1, cfg2); // FIXME: fix this update() function so 2nd arg is const.
    }

    TensorIndex FunctionNode::index_tensors(const ITorchTensorSet::pointer& in) const
    {
        return TensorIndex(in);
    }

    
    TensorIndex FunctionNode::select_tensors(TensorIndex ti) const
    {
        return m_selector.apply(ti);
    }
        
    TensorIndex FunctionNode::transform_tensors(TensorIndex ti) const
    {
        return std::move(ti);
    }

    TensorIndex FunctionNode::sys_transform_tensors(TensorIndex ti) const
    {
        return transform_tensors(std::move(ti));
    }

    TensorIndex FunctionNode::rename_tensors(TensorIndex ti) const
    {
        return m_renaming.apply(std::move(ti));
    }
    

    ITorchTensorSet::pointer FunctionNode::pack_tensors(TensorIndex ti) const
    {
        return ti.as_set();
    }


    bool FunctionNode::operator()(const input_pointer& in, output_pointer& out) const
    {
        out = nullptr;
        if (! in) { return true; } // EOS

        // Weeeeee!
        out = pack_tensors(rename_tensors(sys_transform_tensors(select_tensors(index_tensors(in)))));
        return true;
    }

}
