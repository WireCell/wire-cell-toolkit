#include "WireCellSpng/TensorIndex.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"

#include "WireCellUtil/Exceptions.h"
#include <sstream>

namespace WireCell::SPNG {

    TensorIndex::TensorIndex(ITorchTensorSet::pointer ts)
        : m_ident(ts->ident())
        , m_md(ts->metadata())
    {
        auto tensp = ts->tensors();
        add(tensp->begin(), tensp->end());
    }
    
    TensorIndex::TensorIndex(int ident, const Configuration& md)
        : m_ident(ident)
        , m_md(md)
    {
    }

    TensorIndex::TensorIndex(int ident, const Configuration& md, const ITorchTensor::vector& tens)
        : m_ident(ident)
        , m_md(md)
    {
        add(tens.begin(), tens.end());

    }

    TensorIndex TensorIndex::deepcopy() const
    {
        TensorIndex ti(m_ident, m_md);
        ti.add(m_root);
        return ti; // copy elision is a move
    }

    ITorchTensorSet::pointer TensorIndex::as_set() const
    {
        return std::make_shared<SimpleTorchTensorSet>(m_ident, m_md, tensors());
    }

    /// Return a new TensorIndex with a subtree descending from given parents.
    TensorIndex TensorIndex::subset(const std::vector<ITorchTensor::pointer>& seeds) const
    {
        TensorIndex ti(m_ident, m_md);
        for (const auto& seed : seeds) {
            for (const auto& node : tree_node(seed)->depth()) {
                if (! node.value) { // skip empty root node
                    continue;
                }

                ti.add(node.value);
            }
        }
        return ti;
    }

    // All add()'s flow through this one.
    //
    void TensorIndex::add(ITorchTensor::pointer ten)
    {
        if (! ten) {
            // quietly ignore nullptr
            return;
        }

        auto it = m_nodes.find(ten);
        if (it != m_nodes.end()) {
            // quietly ignore dup
            return;
        }

        const auto md = ten->metadata();

        auto dpath = get<std::string>(md, "datapath", "");
        if (dpath.empty()) {
            raise<ValueError>("no datapath for tensor");
        }
        if (m_bypath.find(dpath) != m_bypath.end()) {
            raise<ValueError>("duplicate datapath: %s", dpath);
        }
        {                       // precheck for TDM compliance
            auto dtype = get<std::string>(md, "datatype", "");
            if (dtype.empty()) {
                raise<ValueError>("no datatype for tensor");
            }
        }

        auto ppath = get<std::string>(md, "parent", "");
        if (ppath == "") {
            // no parent so is a top level tree node.
            m_bypath[dpath] = ten;
            m_nodes[ten] = m_root.insert(ten);
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
    }

    void TensorIndex::add(ITorchTensorSet::pointer ts) {
        if (!ts) {
            // quietly ignore dup
            return;
        }

        ITorchTensor::shared_vector tens = ts->tensors();
        for (auto cit = tens->begin(); cit != tens->end(); ++cit) {
            ITorchTensor::pointer ten = *cit;
            add(ten);
        }
    }

    void TensorIndex::add(const TensorIndex::tree_type& tree)
    {
        for (const auto& node : tree.depth()) { // DFS 
            if (! node.value) { // skip empty root node
                continue;
            }
            add(node.value);
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
        for (const auto& node : m_root.depth()) {
            if (! node.value) { // skip empty root node
                continue;
            }
            if (maxnum == 0) break;

            auto dt = get<std::string>(node.value->metadata(), "datatype", "");
            if (dt == datatype) {
                ret.push_back(node.value);
                --maxnum;
            }
        }
        return ret;
    }

    ITorchTensor::vector TensorIndex::tensors() const
    {
        ITorchTensor::vector tens;
        for (const auto& node : m_root.depth()) {
            if (! node.value) { // skip empty root node
                continue;
            }

            tens.push_back(node.value);
        }
        return tens;
    }

    ITorchTensor::vector TensorIndex::parents() const
    {
        // See nparents() comments to understand the apparently contradictory TI
        // "parent" vs tree "child" nomenclature.
        auto p = m_root.child_values();
        return ITorchTensor::vector(p.begin(), p.end());
    }

    size_t TensorIndex::nparents() const
    {
        // The name mismatch here may be confusing.  The ultimate "parents" in
        // our tree are tensors that themselves lack the "parent" metadata
        // attribute.  These tensors are placed as level 1 nodes in the tree.
        // They are thus the children of the (tensor free) root node.
        return m_root.nchildren();
    }


    ITorchTensor::pointer TensorIndex::parent(ITorchTensor::pointer ten) const
    {
        const auto* node = tree_node(ten);
        if (!node) {
            return nullptr;
        }
        return node->parent->value;        
    }


    ITorchTensor::vector TensorIndex::children(ITorchTensor::pointer ten) const
    {
        const auto* node = tree_node(ten);
        if (!node) {
            return ITorchTensor::vector{};
        }
        auto p = node->child_values();
        return ITorchTensor::vector(p.begin(), p.end());
    }


    const TensorIndex::tree_type* TensorIndex::tree_node(ITorchTensor::pointer ten) const
    {
        auto it = m_nodes.find(ten);
        if (it == m_nodes.end()) {
            return nullptr;
        }
        return it->second;
    }

    std::string TensorIndex::str() const
    {
        std::stringstream ss;
        ss << "TensorIndex id=" << m_ident << " ntensors=" << m_nodes.size()
           << " nparents=" << nparents() << ":";
        for (auto ten : tree().child_values()) {
            auto md = ten->metadata();
            ss << "[";
            ss << md["datatype"].asString() <<":"<< md["datapath"].asString() << ":";
            ss << "( ";
            for (auto siz : ten->shape()) {
                ss << siz << " ";
            }
            ss << "):" << ten->dtype() << ":" << ten->device();
            ss << "] ";
        }

        return ss.str();        
    }


}
