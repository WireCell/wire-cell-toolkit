#include "WireCellSpng/FaninNode.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Exceptions.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGFaninNode,
                 WireCell::SPNG::FaninNode,
                 WireCell::SPNG::ITorchTensorSetFanin,
                 WireCell::IConfigurable,
                 WireCell::INamed);

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

    bool FaninNode::operator()(const input_vector& inv, output_pointer& out) 
    {
        if (inv.size() != m_multiplicity) {
            raise<ValueError>("unexpected multiplicity, got:%d want:%d", inv.size(), m_multiplicity);
        }

        out=nullptr;

        size_t neos = std::count(inv.begin(), inv.end(), nullptr);
        if (neos) {
            logit(String::format("EOS in %d of %d at call=%d", 
                                 neos, m_multiplicity, m_count));
            ++m_count;
            return true;
        }
        
        out = sys_combine_tensors(inv);
        logit(out, "fanned");
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
        logit(out, "combine");
        return out;
    }


}
