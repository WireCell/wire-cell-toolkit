#include "WireCellSpng/FunctionNode.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"

namespace WireCell::SPNG {

    FunctionNode::FunctionNode(const std::string& logname, const std::string& pkgname)
        : Aux::Logger(logname, pkgname)
    {
    }

    void FunctionNode::configure(const WireCell::Configuration& cfg)
    {
        m_quiet = get<bool>(cfg, "quiet", m_quiet);
        m_selector.configure(cfg);
        m_renaming.configure(cfg);
    }
    
    WireCell::Configuration FunctionNode::default_configuration() const
    {
        Configuration cfg;
        cfg["quiet"] = m_quiet;
        auto cfg1 = m_selector.default_configuration();
        auto cfg2 = m_renaming.default_configuration();
        update(cfg, cfg1);
        update(cfg, cfg2);
        return cfg;
    }

    TensorIndex FunctionNode::index_tensors(const ITorchTensorSet::pointer& in) const
    {
        return TensorIndex(in);
    }
    
    TensorIndex FunctionNode::sys_index_tensors(const ITorchTensorSet::pointer& in) const
    {
        auto ti = index_tensors(in);
        maybe_log(ti, "index");
        return ti;              // copy elision
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
        maybe_log(ti, "pre-transform");
        auto new_ti = transform_tensors(std::move(ti));
        maybe_log(new_ti, "post-transform");
        return new_ti;          // copy elision
    }

    TensorIndex FunctionNode::rename_tensors(TensorIndex ti) const
    {
        return m_renaming.apply(std::move(ti));
    }
    
    ITorchTensorSet::pointer FunctionNode::pack_tensors(TensorIndex ti) const
    {
        return ti.as_set();
    }

    ITorchTensorSet::pointer FunctionNode::sys_pack_tensors(TensorIndex ti) const
    {
        maybe_log(ti, "pack");
        return pack_tensors(std::move(ti));
    }

    void FunctionNode::maybe_log(const TensorIndex& ti, const std::string& context) const
    {
        if (m_quiet) return;

        log->debug("{}: call={}: {}", context, m_count, ti.str());
    }

    bool FunctionNode::operator()(const input_pointer& in, output_pointer& out) const
    {
        out = nullptr;
        if (! in) {
            log->debug("EOS: call={}", m_count);
            ++m_count;
            return true;
        }

        // Weeeeee!
        out = sys_pack_tensors(rename_tensors(sys_transform_tensors(select_tensors(sys_index_tensors(in)))));
        ++m_count;
        return true;
    }


}
