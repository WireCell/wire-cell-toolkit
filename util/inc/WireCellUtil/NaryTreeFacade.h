#ifndef WIRECELLUTIL_NARYTREEFACADE
#define WIRECELLUTIL_NARYTREEFACADE

#include "WireCellUtil/NaryTreeNotified.h"

#include <iostream>             // debug

namespace WireCell::NaryTree {

    // A Facade provides a polymorphic base that can be used via a Faced to
    // provide a heterogeneous tree.
    //
    // A Facade is a Notified and so a Facade subclass may implement notify
    // hooks to learn of when this facade's node undergoes some tree level
    // changes.
    //
    // See below for some particular subclasses of Facade that provide
    // additional functionality that may be more suitable as a direct subclass
    // for user code.
    //
    // If a Facade will be used from a Faced, it must be default constructable.
    template<typename Value>
    class Facade : public Notified<Value>
    {
    public:

        virtual ~Facade() {}

        using base_type = Notified<Value>;
        using self_type = Facade<Value>;
        using value_type = Value;
        using node_type = Node<Value>;
        using node_ptr = std::unique_ptr<node_type>;

        virtual void on_construct(node_type* node) {
            m_node = node;
        }

        const node_type* node() const { return m_node; }
        node_type* node() { return m_node; }

    protected:

        node_type* m_node{nullptr};

    };


    // A Faced is something that can hold a "Facade".
    //
    // A Faced is itself also a Facade and thus a Notified.  A Faced intercepts
    // notification from the node in order to forward to itself and to the held
    // facade.  
    //
    // A Faced may be used as a NaryTree::Node Value.
    //
    // A Facade used by a Face must be default constructable.
    template<typename Value>
    class Faced : public Facade<Value> {

    public:
        
        using value_type = Value;
        using base_type = Facade<Value>;
        using self_type = Faced<Value>;
        using typename base_type::node_type;
        using typename base_type::node_ptr;
        using facade_type = Facade<Value>;
        using facade_ptr = std::unique_ptr<facade_type>;

        Faced() = default;
        Faced(Faced&& other) = default;
        Faced& operator=(Faced&& other) = default;
      
        virtual ~Faced() {}

        /// Access the facade as type.  May return nullptr.  Ownership is
        /// retained.
        template<typename FACADE>
        FACADE* facade() {
            if (! m_facade) {
                set_facade(std::make_unique<FACADE>());
            }
            facade_type* base = m_facade.get();
            if (!base) {
                return nullptr;
            }
            FACADE* ret = dynamic_cast<FACADE*>(base);
            if (!ret) {
                return nullptr;
            }
            return ret;
        }

        /// Const access.
        template<typename FACADE>
        const FACADE* facade() const {
            return const_cast<FACADE*>(const_cast<self_type*>(this)->facade<FACADE>());
        }

        /// Access the facade as base type.  May return nullptr.  Ownership is retained.
        const facade_type* facade() const {
            return m_facade.get();
        }
        facade_type* facade() {
            return m_facade.get();
        }

        /// Set the polymorphic facade base with a specific instance.  Caller
        /// may pass nullptr to remove the facade.  This takes ownership of the
        /// facade instance.  
        ///
        /// See facade<T>() which provides create-on-access pattern if the
        /// facade type has a default constructor.
        void set_facade(facade_ptr fac) {
            m_facade = std::move(fac);
            if (this->m_node) {
                m_facade->notify(this->m_node, Action::constructed);
            }
        }

        // Intercept notices from the node in order to forward to the held
        // facade (and to this).
        virtual void notify(node_type* node, Action action) {
            this->base_type::notify(node, action);
            if (m_facade) {
                m_facade->notify(node, action);
            }
        }

    private:

        mutable facade_ptr m_facade{nullptr};

    };                          // Faced

    
    // An interstitial base class for a user facade class for a node that has
    // Faced children with a common type of facade.
    template<typename Child, typename Value>
    class FacadeParent : public Facade<Value> {
    public:
        using child_type = Child;
        using value_type = Value;
        using base_type = Facade<Value>;
        using self_type = FacadeParent<Child, Value>;
        using typename base_type::node_type;
        using typename base_type::node_ptr;
        using children_type = std::vector<child_type*>;

        virtual ~FacadeParent() {}

        // Access collection of children facades.  Const version.
        const children_type& children() const {
            return const_cast<const children_type&>(const_cast<self_type*>(this)->children());
        }

        // Non-const version.  
        children_type& children() {
            if (m_children.empty()) {
                for (auto* cnode : this->m_node->children()) {
                    child_type* child = cnode->value.template facade<child_type>();
                    if (!child) {
                        raise<TypeError>("type mismatch in facade tree node");
                    }
                    m_children.push_back(child);
                }
            }
            return m_children;
        }

        // Number of children this parent has.
        size_t nchildren() const {
            if (this->m_node) {
                return this->m_node->nchildren();
            }
            return 0;
        }

        // Adopt the other's children into this parent.  This leaves
        // other parent childless.
        void take_children(self_type& other, bool notify_value=true) {
            this->m_node->take_children(*other.node(), notify_value);
            invalidate_children();
            other.invalidate_children();
        }

        // Remove kid's node from this parent's node and return an owning
        // pointer which will be nullptr if it's not our kid.
        node_ptr remove_child(child_type& kid, bool notify_value=true) {
            invalidate_children();
            return this->m_node->remove(kid.node(), notify_value);
        }

        // Make a new child, returning its facade.
        child_type& make_child(bool notify_value=true) {
            invalidate_children();
            node_type* cnode = this->m_node->insert(notify_value);
            cnode->value.set_facade(std::make_unique<child_type>());
            return *cnode->value.template facade<child_type>();
        }

        // Return new facade instances of this facade's type.  Each new facade
        // has its node added to this facade's node's parent.  Each new facade
        // node will contain a subset of children that were previously owned by
        // this facade's node as specified by the "groups" vector.  The groups
        // vector is a "connected components" type array giving a group ID for
        // each child node initially owned by this facade's node.  Group IDs are
        // then used as keys in the returned map.  Negative group IDs are
        // ignored and the corresponding child will remain as a child of this
        // facade's node.  See node-level NaryTree::Node::separate() for further
        // description.  Note, this facade's node is left as it was in its
        // parent's children list and new facade nodes are appended to the
        // parent's children list.  This facade's node must have a parent.  
        std::unordered_map<int, self_type*> separate(const std::vector<int> groups, bool notify_value=true) {
            std::unordered_map<int, self_type*> ret;
            auto nurseries = this->m_node->separate(groups, notify_value);
            if (nurseries.empty()) { return ret; }

            auto parent = this->m_node->parent;
            if (!parent) {
                raise<LogicError>("can not separate the children of a facade that lacks a parent");
            }
            
            // Make a new facades of our type with their nodes added to our
            // parent node.  Parent node may not have ParentFacade nor even have
            // a Facade so we must dig into the node level to do this.
            // for (auto& [gid, nur] : nurseries) {
            for (auto& nit : nurseries) {

                // Make a new sibling node.
                node_type* node = parent->insert(notify_value);
                // Give the node its new children.
                node->adopt_children(nit.second, notify_value);
                // Give the node a facade of our type.
                node->value.set_facade(std::make_unique<self_type>());
                // Get the bare facade pointer for return.
                ret[nit.first] = node->value.template facade<self_type>();
            }
            invalidate_children();
            return ret;
        }


        void invalidate_children() {
            m_children.clear();
        }

        template <typename SelfType, typename ParentType>
        friend std::unordered_map<int, SelfType*> separate(SelfType* cluster, const std::vector<int>& groups, bool notify_value);

    private:
        // Lazy cache of children facades.
        children_type m_children;


    };



    // assumes it has a parent and children, make a vector of SelfType* and each one has
    // a subset of children. the length of gorups must match the number of children.
    // children with negative group number will be removed
    template<typename SelfType, typename ParentType>
    std::unordered_map<int, SelfType*> separate(SelfType* cluster, const std::vector<int>& groups, bool notify_value=true) {
        // std::cout << "groups size: " << groups.size() << " nchildren: " << nchildren() << std::endl;
        if(groups.size() != cluster->nchildren()) {
            raise<ValueError>("group size %d mismatch in nchildren %d", groups.size(), cluster->nchildren());
        }
        auto parent = cluster->m_node->parent;
        auto parent_facade = parent->value.template facade<ParentType>();
        parent_facade->invalidate_children(); // clear the facade cache
        std::unordered_map<int, SelfType*> id2facade;
        const auto orig_children = cluster->children(); // make a copy
        for (size_t ichild = 0; ichild < orig_children.size(); ++ichild) {
            // std::cout << "ichild: " << ichild << " group: " << groups[ichild] << std::endl;
            // remove children with negative group number 
            if (groups[ichild] < 0) continue;
            auto* child = orig_children[ichild];
            if (id2facade.find(groups[ichild]) == id2facade.end()) {
                auto* new_snode = parent->insert(notify_value);
                new_snode->value.set_facade(std::make_unique<SelfType>());
                auto new_facade = new_snode->value.template facade<SelfType>();
                id2facade[groups[ichild]] = new_facade;
            }
            // std::cout << "id2facade size: " << id2facade.size() << std::endl;
            auto new_facade = id2facade[groups[ichild]];
            new_facade->m_node->insert(child->node(), notify_value);
            // std::cout << "ichild: " << ichild << " group: " << groups[ichild] << std::endl;
        }
        // remove self from parent
        parent->remove(cluster->m_node, notify_value);
        return id2facade;
    }
}
#endif
