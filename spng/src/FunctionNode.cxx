#include "WireCellSpng/FunctionNode.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGFunctionNode,
                 WireCell::SPNG::FunctionNode,
                 WireCell::SPNG::ITorchTensorSetFilter,
                 WireCell::IConfigurable,
                 WireCell::INamed);


namespace WireCell::SPNG {

    FunctionNode::FunctionNode(const std::string& logname, const std::string& pkgname)
        : Logger(logname, pkgname)
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
        Configuration cfg = this->Logger::default_configuration();

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

        cfg["keep_unselected"] = m_keep_unselected;
        cfg["select_parents"] = m_select_parents;
        cfg["combine_policy"] = m_combine_policy;

        return cfg;
    }

    void FunctionNode::configure(const WireCell::Configuration& cfg)
    {
        this->Logger::configure(cfg);

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

        m_keep_unselected = get(cfg, "keep_unselected", m_keep_unselected);
        m_select_parents = get(cfg, "select_parents", m_select_parents);
        m_combine_policy = get(cfg, "combine_policy", m_combine_policy);

        bool okay = false;
        std::vector<std::string> known_policy = {
            "union_replace",
            "union_keep",
            "input_only",
            "produced_only"
        };
        for (const auto& one : known_policy) {
            if (m_combine_policy == one) {
                log->debug("combine policy: {}", m_combine_policy);
                okay=true;
                break;
            }
        }
        if (! okay) {
            raise<ValueError>("unsupported combine policy: \"%s\"", m_combine_policy);
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
        logit(ti, "indexed");
        return ti;              // copy elision
    }

    TensorIndex FunctionNode::select_tensors(TensorIndex ti) const
    {
        return m_selector.apply(ti, m_keep_unselected, m_select_parents);
    }
    TensorIndex FunctionNode::sys_select_tensors(TensorIndex ti) const
    {
        auto new_ti = select_tensors(std::move(ti));
        logit(new_ti, "selected");
        return new_ti;
    }
        
    TensorIndex FunctionNode::transform_tensors(TensorIndex ti) const
    {
        return std::move(ti);
    }

    TensorIndex FunctionNode::sys_transform_tensors(TensorIndex ti) const
    {
        auto new_ti = transform_tensors(std::move(ti));
        logit(new_ti, "transformed");

        if (m_combine_policy == "input_only") {
            return ti;
        }
        if (m_combine_policy == "produced_only") {
            return new_ti;
        }
        if (m_combine_policy == "union_keep") {
            return new_ti.update(ti); // ti wins datapath collisions
        }
        // default is union_replace.
        return ti.update(new_ti); // new_ti wins
    }

    TensorIndex FunctionNode::rename_tensors(TensorIndex ti) const
    {
        return m_renaming.apply(std::move(ti));
    }

    TensorIndex FunctionNode::sys_rename_tensors(TensorIndex ti) const
    {
        logit(ti, "combined");
        auto new_ti = m_renaming.apply(std::move(ti));
        logit(new_ti, "renamed");        
        return new_ti;
    }
    
    ITorchTensorSet::pointer FunctionNode::pack_tensors(TensorIndex ti) const
    {
        return ti.as_set();
    }

    ITorchTensorSet::pointer FunctionNode::sys_pack_tensors(TensorIndex ti) const
    {
        auto ts = pack_tensors(std::move(ti));
        logit(ts, "packed");
        return ts;

    }

    bool FunctionNode::operator()(const input_pointer& in, output_pointer& out) 
    {
        out = nullptr;
        if (! in) {
            logit("EOS");
            ++m_count;
            return true;
        }

        logit(in, "input");
        // Weeeeee!
        out = sys_pack_tensors(sys_rename_tensors(sys_transform_tensors(sys_select_tensors(sys_index_tensors(in)))));
        logit(out, "output");
        ++m_count;
        return true;
    }

}
