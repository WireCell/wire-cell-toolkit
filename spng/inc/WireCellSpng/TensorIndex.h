#ifndef WIRECELL_SPNG_TENSORINDEX
#define WIRECELL_SPNG_TENSORINDEX

#include "WireCellSpng/ITorchTensorSet.h"
#include "WireCellUtil/NaryTree.h"

#include <limits>

namespace WireCell::SPNG {

    /// Manage a vector of tensors (ITorchTensor::pointer).
    ///
    /// This indexes a set of tensors into a parentage tree and other structures
    /// to allow fast lookup.  Tensors may be added over the life time of the
    /// index.  Tensors may not be removed.  If a subset are required, make a
    /// new TensorIndex and add the desired subset  
    ///
    /// The index performs to de-duplication at the ITorchTensor::pointer level.
    ///
    /// The index will generally raise ValueError exception if the TDM is violated.
    class TensorIndex {
    public:
        TensorIndex() = default;
        ~TensorIndex() = default;

        /// Construct on an iterator range of ITorchTensor::pointer
        template<typename It> TensorIndex(It beg, It end) : m_tens(beg, end) {}
        template<typename Range> TensorIndex(Range& r) : m_tens(std::begin(r), std::end(r)) {}
        
        /// Construct from tensors in a set.
        TensorIndex(ITorchTensorSet::pointer ts) : m_tens(ts->tensors()->begin(), ts->tensors()->end()) {}

        /// A tree representing tensor parentage as expressed through the TDM
        /// "parent" metadata attribute.
        ///
        /// User SHOULD read about NaryTree in WCT util in order to understand
        /// terms.
        ///
        /// Here, the root node of the tree will always have nullptr value.  Its
        /// children form Layer 1 and hold shared pointers to all the
        /// ITorchTensors that lack any "parent" metadata attribute.  Each node
        /// in Layer 1 has children in Layer 2 which have the Layer 1 tensor as
        /// their "pairent", etc.
        using tree_type = NaryTree::Node<ITorchTensor::pointer>;

        /// Return the "root" tree node.  The value held by this node is always
        /// nullptr.
        const tree_type& tree() const { return m_root; }

        /// Return the tree node holding the tensor or nullptr.  Note, despite a
        /// bare pointer return, the index retains ownership.
        const tree_type* tree_node(ITorchTensor::pointer ten) const;

        /// Descend from node, adding all tensors found to this index.  Tensors
        /// are added in depth-first descent order.
        void add(const tree_type& node);

        /// Add tensors to the index via iterator range.
        template<typename It> void add(It beg, It end) {
            for (It it = beg; it != end; ++it) {
                add(*it);
            }
        }

        // Add from a tensor set.
        void add(ITorchTensorSet::pointer ts);

        /// Add a single tensor.
        ///
        /// All other add()'s go through this.
        ///
        /// This method strictly enforces the TDM.  This is intentional.  Do not
        /// relax it.  If you cause it to throw an exception, fix upstream so
        /// that tensors that obey the SPNG TDM are provided.
        void add(ITorchTensor::pointer ten);


        /// Return tensor at datapath or nullptr if no such path.
        ITorchTensor::pointer at_path(const std::string& datapath) const;

        /// Return first tensors of a given datatype up to a maximum number.
        ITorchTensor::vector of_type(const std::string& datatype,
                                     size_t maxnum=std::numeric_limits<size_t>::max()) const;

    private:
        std::vector<ITorchTensor::pointer> m_tens;

        tree_type m_root{nullptr};
        std::unordered_map<std::string, ITorchTensor::pointer> m_bypath;
        std::unordered_map<ITorchTensor::pointer, tree_type*> m_nodes;

    };

}

#endif
