#include "WireCellSpng/TensorIndex.h"
#include "WireCellUtil/Exceptions.h"

namespace WireCell::SPNG {

    // All add()'s flow through this one.
    //
    void TensorIndex::add(ITorchTensor::pointer ten)
    {
        if (! ten) return;

        auto it = m_nodes.find(ten);
        if (it != m_nodes.end()) {
            return;
        }

        const auto md = ten->metadata();

        auto dpath = get<std::string>(md, "datapath", "");
        // Fixme: error if dp is empty or duplicate
        if (dpath.empty()) {
            raise<ValueError>("no datapath for tensor");
        }
        if (m_bypath.find(dpath) != m_bypath.end()) {
            raise<ValueError>("duplicate datapath: %s", dpath);
        }

        auto ppath = get<std::string>(md, "parent", "");

        if (ppath == "") {
            m_bypath[dpath] = ten;
            m_nodes[ten] = m_root.insert(ten);
            m_tens.push_back(ten);
            return;
        }


        auto pit = m_bypath.find(ppath);
        if (pit == m_bypath.end()) {
            raise<ValueError>("no parent tensor at path %s", ppath);
        }
        
        auto pptr = pit->second;
            
        auto pnit = m_nodes.find(pptr);
        if (pnit == m_nodes.end()) {
            raise<ValueError>("parent missing from tree, at path %s", ppath);
        }

        auto pnode = pnit->second;
        m_bypath[dpath] = ten;
        m_nodes[ten] = pnode->insert(ten);
        m_tens.push_back(ten);
    }


    void TensorIndex::add(ITorchTensorSet::pointer ts) {
        if (!ts) return;

        ITorchTensor::shared_vector tens = ts->tensors();
        for (auto cit = tens->begin(); cit != tens->end(); ++cit) {
            ITorchTensor::pointer ten = *cit;
            add(ten);
        }
    }

    void TensorIndex::add(const TensorIndex::tree_type& node)
    {
        for (ITorchTensor::pointer ten : node.child_values()) { // DFS 
            add(ten);
        }
    }

    ITorchTensor::pointer TensorIndex::at_path(const std::string& datapath) const
    {
        auto it = m_bypath.find(datapath);
        if (it == m_bypath.end()) {
            return nullptr;
        }
        return it->second;
    }

    ITorchTensor::vector TensorIndex::of_type(const std::string& datatype, size_t maxnum) const
    {
        // This does a full scan.
        ITorchTensor::vector ret;
        for (const auto& one : m_tens) {
            if (maxnum == 0) break;

            auto dt = get<std::string>(one->metadata(), "datatype", "");
            if (dt == datatype) {
                ret.push_back(one);
                --maxnum;
            }
        }
        return ret;
    }

    const TensorIndex::tree_type* TensorIndex::tree_node(ITorchTensor::pointer ten) const
    {
        auto it = m_nodes.find(ten);
        if (it == m_nodes.end()) {
            return nullptr;
        }
        return it->second;
    }

}
