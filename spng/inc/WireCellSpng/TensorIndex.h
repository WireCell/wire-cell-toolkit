#ifndef WIRECELL_SPNG_TENSORINDEX
#define WIRECELL_SPNG_TENSORINDEX

#include "WireCellSpng/ITorchTensorSet.h"
#include "WireCellUtil/NaryTree.h"

#include <limits>

namespace WireCell::SPNG {

    /// SPNG TDM compliant tensor management.
    ///
    /// This indexes the set of TDM tensors in the following ways:
    ///
    /// - named map :: a tensor may be retrieved by its datapath.
    ///
    /// - parentage tree :: the "parent" attribute is used to form a
    ///   parent/child tree (NaryTree).  This allows quick parent/child, sibling
    ///   and depth-first-descent relationship traversal.
    ///
    /// An index is "add-only", no removal and not in-place mutation is
    /// supported.  To remove or mutate, one must make a new index from an
    /// existing index.
    ///
    /// The TensorIndex requires and maintains SPNG TDM compliance.  Any tensors
    /// retrieved from the index will be TDM compliant.  Any tensors given to
    /// the index must be TDM compliant or a ValueError exception will be
    /// raised.
    ///
    /// The index caries an ident number and metadata.  Their values are
    /// intended to be provided by an initial ITorchTensorSet.  The index
    /// otherwise ignores them.
    ///
    class TensorIndex {
    public:

        /// Construct empty.
        TensorIndex() = default;
        ~TensorIndex() = default;

        /// No copy constructor and copy assignment because this requires a deep
        /// copy and we want to catch inadvertent copies due to badly designed
        /// function calls.
        TensorIndex(const TensorIndex& other) = delete;        
        TensorIndex& operator=(const TensorIndex& rhs) = delete;

        // If you really need a copy:
        TensorIndex deepcopy() const;

        TensorIndex(TensorIndex&& other) = default;
        TensorIndex& operator=(TensorIndex&& rhs) = default;

        /// Construct from tensors in a set.
        TensorIndex(ITorchTensorSet::pointer ts); 

        /// Construct empty of tensors but with ident number and metadata.
        TensorIndex(int ident, const Configuration& md);

        /// Construct from parts.
        TensorIndex(int ident, const Configuration& md, const ITorchTensor::vector& tens);

        /// Pack tensors, ident and metadata into a tensor set.
        ITorchTensorSet::pointer as_set() const;

        /// Return a new TensorIndex formed by depth first descent of the tree
        /// starting at the nodes holding the given seeds.  Seeds not held are
        /// ignored.  Caller may compare the nparents() on the returned index
        /// the size of seeds() to check for loss.
        TensorIndex subset(const std::vector<ITorchTensor::pointer>& seeds) const;

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

        /// Return the parent tensor.  If ten lacks a parent (ie, is an ultimate
        /// parent), nullptr is returned.
        ITorchTensor::pointer parent(ITorchTensor::pointer ten) const;

        /// Return vector of direct children of tensor.  Vector is empty is the
        /// tensor is a leaf or if the tensor is not in the index.
        ITorchTensor::vector children(ITorchTensor::pointer ten) const;

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

        /// Return flat vector of tensors in tree DFS order.
        ///
        /// The order of the vector is TDM-compliant.
        ///
        /// Note, if the goal is to simply iterate on all tensors, overhead of
        /// constructing the intermediate vector can be avoided by using code
        /// like:
        /// @code{.cpp}
        /// for (const auto& node : ti.tree().depth()) {
        ///     if (! node.value) { // skip empty root node
        ///         continue;
        ///     }
        ///     auto ten = node.value->tensor();
        ///     // ...
        /// }
        /// @endcode
        ITorchTensor::vector tensors() const;

        /// Return top level parents in sibling order from the tree.
        ///
        /// Note, if the goal is to simply iterate on parent tensors, the
        /// overhead of forming the intermediate vector can be avoided by using
        /// code like:
        ///
        /// @code{.cpp}
        /// for (auto iten : ti.tree().child_values()) {
        ///     auto ten = iten->tensor();
        ///     // ...
        /// }
        /// @endcode
        ITorchTensor::vector parents() const;

        /// Return the number of top-level parents.
        ///
        /// This is equal to the parents().size() but avoids forming the vector.
        size_t nparents() const;

        /// Return the ident number.
        int ident() const { return m_ident; }

        /// Return the metadata object.
        const Configuration& metadata() const { return m_md; }

    private:
        int m_ident{-1};
        Configuration m_md;

        tree_type m_root{nullptr};
        std::unordered_map<std::string, ITorchTensor::pointer> m_bypath;
        std::unordered_map<ITorchTensor::pointer, tree_type*> m_nodes;

    };

}

#endif
