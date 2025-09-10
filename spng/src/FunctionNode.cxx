#include "WireCellSpng/FunctionNode.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"

namespace WireCell::SPNG {

    FunctionNode::FunctionNode(const std::string& logname, const std::string& pkgname)
        : Aux::Logger(logname, pkgname)
    {
    }

    void FunctionNode::register_input(const std::string& kind, const std::string& default_datapath)
    {
        m_input_datapath[kind] = default_datapath;
    }

    void FunctionNode::register_output(const std::string& kind, const std::string& default_datapath)
    {
        m_output_datapath[kind] = default_datapath;
    }

    WireCell::Configuration FunctionNode::default_configuration() const
    {
        Configuration cfg;
        cfg["quiet"] = m_quiet;
        auto cfg1 = m_selector.default_configuration();
        auto cfg2 = m_renaming.default_configuration();
        update(cfg, cfg1);
        update(cfg, cfg2);

        for (const auto& [kind, datapath] : m_input_datapath) {
            cfg["input_datapath"][kind] = datapath;
        }
        for (const auto& [kind, datapath] : m_output_datapath) {
            cfg["output_datapath"][kind] = datapath;
        }

        return cfg;
    }

    void FunctionNode::configure(const WireCell::Configuration& cfg)
    {
        m_quiet = get<bool>(cfg, "quiet", m_quiet);
        m_selector.configure(cfg);
        m_renaming.configure(cfg);

        auto its = cfg["input_tensors"];
        if (its.isObject()) {
            for (const auto& kind : its.getMemberNames()) {
                m_input_datapath[kind] = its[kind].asString();
            }
        }
        auto ots = cfg["output_tensors"];
        if (ots.isObject()) {
            for (const auto& kind : ots.getMemberNames()) {
                m_output_datapath[kind] = ots[kind].asString();
            }
        }
    }
    
    std::string FunctionNode::input_datapath(const std::string& kind) const
    {
        auto it = m_input_datapath.find(kind);
        if (it == m_input_datapath.end()) {
            return "";
        }
        return it->second;
    }

    std::string FunctionNode::output_datapath(const std::string& kind) const
    {
        auto it = m_output_datapath.find(kind);
        if (it == m_output_datapath.end()) {
            return "";
        }
        return it->second;
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
