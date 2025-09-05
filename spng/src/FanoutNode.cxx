#include "WireCellSpng/FanoutNode.h"

namespace WireCell::SPNG {

    void FanoutNode::configure(const WireCell::Configuration& cfg)
    {
        m_multiplicity = get(cfg, "multiplicity", m_multiplicity);
    }

    WireCell::Configuration FanoutNode::default_configuration() const
    {
        Configuration cfg;
        cfg["multiplicity"] = m_multiplicity;
        return cfg;
    }

    std::vector<std::string> FanoutNode::output_types()
    {
        // FIXME: we should make m_multiplicity be a protected member of the
        // specific I<TYPE>Fanout and move this method body up into that
        // intermediate interface.
        const std::string tname = std::string(typeid(output_type).name());
        std::vector<std::string> ret(m_multiplicity, tname);
        return ret;
    }

    bool FanoutNode::operator()(const input_pointer& in, output_vector& outv) const
    {
        outv.clear();
        if (! in) return true;  // EOS

        outv = separate_tensors(in);
        return true;
    }

    ITorchTensorSet::vector FanoutNode::separate_tensors(const ITorchTensorSet::pointer& in) const
    {
        return ITorchTensorSet::vector(m_multiplicity, in);
    }

}
