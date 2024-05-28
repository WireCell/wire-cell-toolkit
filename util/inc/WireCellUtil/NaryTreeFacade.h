#ifndef WIRECELLUTIL_NARYTREEFACADE
#define WIRECELLUTIL_NARYTREEFACADE

#include "WireCellUtil/NaryTreeNotified.h"

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
            // std::cerr << "NaryTree::Faced::notify(" << (void*)node << "," << action << ")\n";
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
            // std::cerr << "take_children() " << other.nchildren() << "\n";
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

    private:
        // Lazy cache of children facades.
        children_type m_children;

        void invalidate_children() {
            m_children.clear();
        }


    };
}
#endif
