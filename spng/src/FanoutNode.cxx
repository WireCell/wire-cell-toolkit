#include "WireCellSpng/FanoutNode.h"

namespace WireCell::SPNG {

    FanoutNode::FanoutNode(const std::string& logname, const std::string& pkgname)
        : FanBase(logname, pkgname)
    {
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
        if (! in) {
            maybe_log(nullptr, "EOS");
            ++m_count;
            return true;
        }

        outv = sys_separate_tensors(in);
        ++m_count;
        return true;
    }

    ITorchTensorSet::vector FanoutNode::separate_tensors(const ITorchTensorSet::pointer& in) const
    {
        // simply dup pointer M-ways
        return ITorchTensorSet::vector(m_multiplicity, in);
    }

    ITorchTensorSet::vector FanoutNode::sys_separate_tensors(const ITorchTensorSet::pointer& in) const
    {
        maybe_log(in, "separate");
        return this->separate_tensors(in);
    }

}
