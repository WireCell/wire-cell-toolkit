/** N-ary tree.
   
    A tree is a graph where the vertices (nodes) are organized in
    parent/child relationship via graph edges.  A node may connect to
    at most one parent.  If a node has no parent it is called the
    "root" node.  A node may edges to zero or more children nodes for
    which the node is their parent.  Edges are bi-directional with the
    "downward" direction being identified as going from a parent to a
    child.

    In this implementation there is no explicit tree object.  The
    singular root node can be considered conceptually as "the tree"
    though it is possible to naviate to all nodes given any node in
    "the" tree.

    The tree is built by inserting child values or nodes into a
    parent node.

    Several types of iterators and ranges provide different methods to
    descend through the tree to provide access to the node value.

*/

#ifndef WIRECELLUTIL_NARYTREE
#define WIRECELLUTIL_NARYTREE

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/DetectionIdiom.h"

#include <vector>
#include <list>
#include <memory>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>

// #include <iostream>             // debug

namespace WireCell::NaryTree {

    template<typename Iter>
    struct iter_range {
        Iter b,e;
        Iter begin() { return b; }
        Iter end()   { return e; }
    };

    /** Iterator and range on children values */
    template<typename Value> class child_value_iter;
    template<typename Value> class child_value_const_iter;
    
    /** Depth first traverse, parent then children, as range */
    template<typename Value> class depth_range;

    // Notify this node's VALUE that a tree action took place on the NODE.  To
    // catch an action, a node value must implement a method with signature:
    // notify(NODE,ACTION).  See NaryTree::Notified for an interstitial base
    // class that manages notifications.
    enum Action {
        constructed,            // the node is constructed
        inserted,               // the node has been inserted into a parent's children list
        removing,               // the node has been removed from a parent's children list
    };


    /** A tree node.

        Value gives the type of the node "payload" data.

     */
    template<class Value>
    class Node {

      public:

        using value_type = Value;
        using self_type = Node<Value>;
        // User access to ordered collection of children
        using children_vector = std::vector<self_type*>;
        using children_const_vector = std::vector<self_type const*>;
        // node owns child nodes
        using owned_ptr = std::unique_ptr<self_type>;
        // holds the children
        using nursery_type = std::list<owned_ptr>; 
        using sibling_iter = typename nursery_type::iterator;
        using sibling_const_iter = typename nursery_type::const_iterator;

        // A depth first descent range
        using range = depth_range<self_type>;
        using const_range = depth_range<self_type const>;

        // Access the payload data value
        value_type value;

        // Access parent that holds this as a child.
        self_type* parent{nullptr};

        ~Node() = default;

        Node() {
            notify<Value>(Action::constructed, this);
        }

        // Copy value.
        Node(const value_type& val) : value(val) {
            notify<Value>(Action::constructed, this);
        }

        // Move value.
        Node(value_type&& val) : value(std::move(val)) {
            notify<Value>(Action::constructed, this);
        }

        // Insert a child node of default value
        Node* insert(bool notify_value=true) {
            owned_ptr nptr = std::make_unique<self_type>();
            return insert(std::move(nptr), notify_value);
        }

        // Insert a child by its value copy.
        Node* insert(const value_type& val, bool notify_value=true) {
            owned_ptr nptr = std::make_unique<self_type>(val);
            return insert(std::move(nptr), notify_value);
        }

        // Insert a child by its value move.
        Node* insert(value_type&& val, bool notify_value=true) {
            owned_ptr nptr = std::make_unique<self_type>(std::move(val));
            return insert(std::move(nptr), notify_value);
        }

        // Insert, reparent and take ownership of a bare node pointer.  The
        // Node* MUST be on the heap.  A lent pointer is returned.
        Node* insert(Node* node, bool notify_value=true) {
            owned_ptr nptr;
            if (node->parent) {
                nptr = node->parent->remove(node, notify_value);
            }
            else {
                nptr.reset(node);
            }
            return insert(std::move(nptr), notify_value);
        }

        // Insert, reparent and take ownership of a unique node pointer.  A lent
        // pointer is returned.
        Node* insert(owned_ptr node, bool notify_value=true) {
            if (!node or !node.get()) {
                throw std::runtime_error("NaryTree::Node insert on null node");
            }

            if (node->parent) {
                node = node->parent->remove(node.get(), notify_value);
            }
            node->parent = this;
            nursery_.push_back(std::move(node));
            nursery_.back()->sibling_ = std::prev(nursery_.end());
            Node* child = nursery_.back().get();
            if (!child) {
                throw std::runtime_error("NaryTree::Node insert on null child node");
            }
            if (notify_value) {
                notify<Value>(Action::inserted, child);
            }
            return child;
        }

        // Return iterator to node or end.  This is a linear search.
        sibling_iter find(const Node* node) {
            return std::find_if(nursery_.begin(), nursery_.end(),
                                [&](const owned_ptr& up) {
                                    return up.get() == node;
                                });
        }
        sibling_const_iter find(const Node* node) const {
            return std::find_if(nursery_.cbegin(), nursery_.cend(),
                                [&](const owned_ptr& up) {
                                    return up.get() == node;
                                });
        }
        
        // Remove and return child node.
        owned_ptr remove(sibling_iter sib, bool notify_value=true) {
            if (sib == nursery_.end()) return nullptr;

            Node* child = (*sib).get();
            if (notify_value) {
                notify<Value>(Action::removing, child);
            }

            owned_ptr ret = std::move(*sib);
            nursery_.erase(sib);
            ret->sibling_ = nursery_.end();
            ret->parent = nullptr;
            
            return ret;
        }

        // Remove child node.  Searches children.  Return child node
        // in owned pointer if found else nullptr.
        owned_ptr remove(const Node* node, bool notify_value=true) {
            auto it = find(node);
            return remove(it, notify_value);
        }

        // Return a nursery of all children, leaving the one in this node empty.
        // If notify_child is true, notify the children of their removal.
        nursery_type remove_children(bool notify_value=true) {
            // std::cerr << "NaryTree::Node::remove_children() removing " << nursery_.size() << "\n";
            nursery_type ret;
            while (nursery_.begin() != nursery_.end()) {
                auto orphan = remove(nursery_.begin(), notify_value);
                ret.emplace_back(std::move(orphan));
            }
            return ret;         // this is a move.
        }            

        // Insert children in given nursery, depleting it.
        void adopt_children(nursery_type& kids, bool notify_value=true) {
            // std::cerr << "NaryTree::Node::adopt_children() adopting " << kids.size() << "\n";
            for (auto it = kids.begin(); it != kids.end(); ++it) {
                insert(std::move(*it), notify_value);
            }
            kids.clear();
        }

        // Transfer children from other to self.
        void take_children(self_type& other, bool notify_value = true) {
            // std::cerr << "NaryTree::Node::take_children() taking " << other.nchildren() << "\n";
            auto kids = other.remove_children(notify_value);
            adopt_children(kids, notify_value);
        }

        // Iterator locating self in list of siblings.  If parent is
        // null, this iterator is invalid.  It is set which this node
        // is inserted as a anothers child.
        sibling_iter sibling_;
        sibling_iter sibling() const {
            if (!parent) {
                raise<ValueError>("node with no parent is not a sibling");
            }
            return sibling_;
        }

        // Return index of self in parent's list of children.  This is
        // an O(nchildren) call and will throw if we have no parent.
        size_t sibling_index() const {
            auto me = sibling(); // throws
            auto first_born = parent->nursery_.begin();
            return std::distance(first_born, me);
        }
        // Call sibling_index recursively up toward root, returning
        // result with my index first.
        std::vector<size_t> sibling_path() const {
            std::vector<size_t> ret;
            const auto* n = this;
            while (n) {
                if (! n->parent) {
                    break;
                }
                ret.push_back(n->sibling_index());
                n = n->parent;
            }
            return ret;
        }


        self_type* first() const {
            if (nursery_.empty()) return nullptr;
            return nursery_.front().get();
        }

        self_type* last() const {
            if (nursery_.empty()) return nullptr;
            return nursery_.back().get();
        }

        // Return left/previous/older sibling, nullptr if we are first.
        self_type* prev() const {
            if (!parent) return nullptr;
            const auto& sibs = parent->nursery_;
            if (sibs.empty()) return nullptr;

            if (sibling_ == sibs.begin()) {
                return nullptr;
            }
            auto sib = sibling_;
            --sib;
            return sib->get();
        }
        // Return right/next/newerr sibling, nullptr if we are last.
        self_type* next() const {
            if (!parent) return nullptr;
            const auto& sibs = parent->nursery_;
            if (sibs.empty()) return nullptr;

            auto sib = sibling_;
            ++sib;
            if (sib == sibs.end()) {
                return nullptr;
            }
            return sib->get();
        }

        // Access collection of child nodes as bare pointers.
        size_t nchildren() const { return nursery_.size(); }
        children_const_vector children() const {
            children_const_vector ret(nursery_.size());
            std::transform(nursery_.begin(), nursery_.end(), ret.begin(),
                           [](const auto& up) { return up.get(); });
            return ret;
        }
        children_vector children() {
            children_vector ret(nursery_.size());
            std::transform(nursery_.begin(), nursery_.end(), ret.begin(),
                           [](const auto& up) { return up.get(); });
            return ret;
        }

        // FIXME TEMPORARY ACCESS WHILE FIXING CLUSTERING 
        // Access the owning nursery of children directly.
        // const nursery_type& nursery() const { return nursery_; }
        // nursery_type& nursery() { return nursery_; }

        using child_value_range = iter_range<child_value_iter<Value>>;
        auto child_values() {
            return child_value_range{
                child_value_iter<Value>(nursery_.begin()),
                child_value_iter<Value>(nursery_.end()) };
        }
        using child_value_const_range = iter_range<child_value_const_iter<Value>>;
        auto child_values() const {
            return child_value_const_range{
                child_value_const_iter<Value>(nursery_.cbegin()),
                child_value_const_iter<Value>(nursery_.cend()) };
        }

        using child_node_iter = boost::indirect_iterator<sibling_iter>;
        using child_node_range = iter_range<child_node_iter>;
        auto child_nodes() {
            return child_node_range{ nursery_.begin(), nursery_.end() };
        }
        using child_node_const_iter = boost::indirect_iterator<sibling_const_iter>;
        using child_node_const_range = iter_range<child_node_const_iter>;
        auto child_nodes() const {
            return child_node_const_range{ nursery_.cbegin(), nursery_.cend() };
        }
                
        // Iterable range for depth first traversal, parent then children.
        // Iterators yield a reference to the node.
        // Level=0 will traverse to the leaves.
        // Level=1 will only visit the current node.
        // Level=2 will visit current node and children nodes
        // Level=3 etc
        range depth(size_t level=0) { return range(this, level); }
        const_range depth(size_t level=0) const { return const_range(this, level); }

      private:

        // Detect Value::notify(const Node<Value>* node, Action action)

        template <typename T, typename ...Ts>
        using notify_type = decltype(std::declval<T>().notify(std::declval<Ts>()...));

        template<typename T>
        using has_notify = is_detected<notify_type, T, Node*, Action>;

        template <class T, std::enable_if_t<has_notify<T>::value>* = nullptr>
        void notify(Action action, Node* node) {
            // std::cerr << "sending action: "<<action<<" \n";
            node->value.notify(node, action);
        }

        template <class T, std::enable_if_t<!has_notify<T>::value>* = nullptr>
        void notify(Action action, Node* node) {
            // std::cerr << "missing action: "<<action<<" \n";
            return; // no-op
        }

      private:

        // friend class child_value_iter<Value>;
        nursery_type nursery_;

    };                          // Node



    //
    // Iterator on node childrens' value.
    //
    template<typename Value>
    class child_value_iter : public boost::iterator_adaptor<
        child_value_iter<Value>              // Derived
        , typename Node<Value>::sibling_iter // Base
        , Value                              // Value
        , boost::bidirectional_traversal_tag>
    {
      private:
        struct enabler {};  // a private type avoids misuse

      public:

        using sibling_iter = typename Node<Value>::sibling_iter;

        child_value_iter()
            : child_value_iter::iterator_adaptor_()
        {}

        explicit child_value_iter(sibling_iter sit)
            : child_value_iter::iterator_adaptor_(sit)
        {}

        template<typename OtherValue>
        child_value_iter(child_value_iter<OtherValue> const& other
                   , typename boost::enable_if<
                   boost::is_convertible<OtherValue*,Value*>
                   , enabler
                   >::type = enabler()
            )
            : child_value_iter::iterator_adaptor_(other->nursery_.begin())
        {}

      private:
        friend class boost::iterator_core_access;
        Value& dereference() const {
            auto sib = this->base();
            return (*sib)->value;
        }
    };    

    // Const iterator on node childrens' value.
    template<typename Value>
    class child_value_const_iter : public boost::iterator_adaptor<
        child_value_const_iter<Value>              // Derived
        , typename Node<Value>::sibling_const_iter // Base
        , Value                                    // Value
        , boost::bidirectional_traversal_tag>
    {
      private:
        struct enabler {};  // a private type avoids misuse

      public:

        using sibling_const_iter = typename Node<Value>::sibling_const_iter;

        child_value_const_iter()
            : child_value_const_iter::iterator_adaptor_()
        {}

        explicit child_value_const_iter(sibling_const_iter sit)
            : child_value_const_iter::iterator_adaptor_(sit)
        {}

        // convert non-const
        template<typename OtherValue>
        child_value_const_iter(child_value_iter<OtherValue> const& other
                   , typename boost::enable_if<
                   boost::is_convertible<OtherValue*,Value*>
                   , enabler
                   >::type = enabler()
            )
            : child_value_const_iter::iterator_adaptor_(other->nursery_.begin())
        {}

      private:
        friend class boost::iterator_core_access;

        // All that just to provide:
        Value& dereference() const {
            auto sib = this->base();
            return (*sib)->value;
        }
    };    


    // Iterator for depth-first descent.  First parent then children.
    // This becomes a const_iterator with a NodeType = SomeType const.
    template<typename NodeType>
    class depth_iter : public boost::iterator_facade<
        depth_iter<NodeType>
        , NodeType
        , boost::forward_traversal_tag> // fixme: could make this bidirectional
    {

      public: 

        using node_type = NodeType;
        // Current node.  If nullptr, iterator is at "end".
        node_type* node{nullptr};

        // How deep we may go.  0 is unlimited.  1 is only given node.
        // 2 is initiial node's children only.
        size_t depth{0};
        // THe level we are on.  0 means "on initial node level".  We
        // never go negative.  Always: 0<=level.  If depth>0 then
        // 0<=level<depth.
        size_t level{0};

        using self_type = depth_iter<NodeType>;
        using base_type = boost::iterator_facade<self_type
                                                 , NodeType
                                                 , boost::forward_traversal_tag>;
        using difference_type = typename base_type::difference_type;
        using value_type = typename base_type::value_type;
        // using value_type = typename NodeType::value_type;
        using pointer = typename base_type::pointer;
        using reference = typename base_type::reference;


        // default is nodeless and indicates "end"
        depth_iter() : node(nullptr) {}

        explicit depth_iter(node_type* node) : node(node) { }
        explicit depth_iter(node_type* node, size_t depth) : node(node), depth(depth) { }

        template<typename OtherValue>
        depth_iter(depth_iter<OtherValue> const & other)
            : node(other.node), depth(other.depth), level(other.level) {}

      private:
        friend class boost::iterator_core_access;

        // Limited by level and depth:
        // If we have child, we go there.
        // If no child, we go to next sibling.
        // If no next sibling we go toward root taking first ancestor sibling found.
        // If still nothing, we are at the end.
        void increment()
        {
            if (!node) return; // already past the end, throw here?

            // Down first

            // If we have not yet hit floor
            if (depth==0 or level+1 < depth) {
                // and if we have child, go down
                auto* first = node->first();
                if (first) {
                    node = first;
                    ++level;
                    return;
                }
            }

            // Can not go down, seek for sibling
            while (level) {     // but not above original node level

                // first try sibling in current level
                auto* sib = node->next();
                if (sib) {
                    node = sib;
                    return;
                }
                if (node->parent) {
                    node = node->parent;
                    --level;
                    continue;
                }
                break;
            }
            node = nullptr;     // end
        }

        template <typename OtherNode> 
        bool equal(depth_iter<OtherNode> const& other) const {
            return node == other.node;
        }

        reference dereference() const {
            // return node->value;
            return *node;
        }
    };


    // A range spanning a depth first search.
    template<typename NodeType>
    class depth_range {
      public:
        using node_type = NodeType;
        using iterator = depth_iter<NodeType>;

        depth_range() : root(nullptr), depth(0) {}
        depth_range(node_type* root) : root(root), depth(0) {}
        depth_range(node_type* root, size_t depth) : root(root), depth(depth) {}
        
        // Range interface
        iterator begin() { return iterator(root, depth); }
        iterator end() { return iterator(); }

      private:
        node_type* root{nullptr};
        size_t depth{0};
    };


} // WireCell::NaryTree

#endif
