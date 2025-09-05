#include "WireCellSpng/FaninNode.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"

namespace WireCell::SPNG {

    FaninNode::FaninNode(const std::string& logname, const std::string& pkgname)
        : FanBase(logname, pkgname)
    {
    }

    std::vector<std::string> FaninNode::input_types()
    {
        const std::string tname = std::string(typeid(input_type).name());
        std::vector<std::string> ret(m_multiplicity, tname);
        return ret;
    }

    bool FaninNode::operator()(const input_vector& inv, output_pointer& out) const
    {
        if (inv.empty()) {
            // This should not happen, may indicate a bad graph config.  But, it's not EOS
            out = std::make_shared<EmptyTorchTensorSet>();
            maybe_log(nullptr, "EMPTY");
            ++m_count;
            return true;
        }

        out=nullptr;
        // Output EOS if any input is EOS.
        int neos = 0;
        for (const auto& one : inv) {
            if (!one) {
                ++neos;
            }
        }
        if (neos == m_multiplicity) {
            maybe_log(nullptr, "EOS");
            ++m_count;
            return true;
        }
        if (neos) {             // should never happen
            maybe_log(nullptr, "partial-EOS");
            ++m_count;
            return true;
        }
        
        out = sys_combine_tensors(inv);
        ++m_count;
        return true;
    }

    ITorchTensorSet::pointer FaninNode::combine_tensors(const ITorchTensorSet::vector& inv) const
    {
        const int ident = inv[0]->ident();
        Configuration md;
        auto tens = std::make_shared<ITorchTensor::vector>();
        for (const auto& its : inv) {
            auto new_md = its->metadata();
            update(md, new_md);
            for (const auto& ten : *its->tensors()) {
                tens->push_back(ten);
            }
        }
        return std::make_shared<SimpleTorchTensorSet>(ident, md, tens);
    }

    ITorchTensorSet::pointer FaninNode::sys_combine_tensors(const ITorchTensorSet::vector& inv) const
    {
        auto out = combine_tensors(inv);
        maybe_log(out, "combine");
        return out;
    }


}
