#ifndef WIRECELLUTIL_NARYTREEFACADE
#define WIRECELLUTIL_NARYTREEFACADE

#include "WireCellUtil/NaryTreeNotified.h"

#include <map>

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
    // A Faced may be used as a NaryTree::Node Value.  Of particular note,
    // PointCloud::Tree::Points is a Faced.
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
      
        virtual ~Faced() { }

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
            if (m_facade && this->m_node) {
                // Call on facade directly and NOT this->notify().
                m_facade->notify({this->m_node}, Action::constructed);
            }
        }

        // This overrides the base in order to also dispatch notice to the
        // facade, if set.
        //
        // Note, this will NOT set a facade as a side-effect.  If your facade is
        // not getting notices you must call facade<T>() or set_facade() prior
        // to generating any notifications.
        virtual bool notify(std::vector<node_type*> path, Action action) {
            bool propagate = this->base_type::notify(path, action);
            if (!propagate) {
                return false;
            }
            if (m_facade) {
                return m_facade->notify(path, action);
            }
            return true;
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


        // React to someone removing a cluster node from our node
        virtual bool on_remove(const std::vector<node_type*>& path) {
            invalidate_children();
            return true;
        }

        virtual bool on_insert(const std::vector<node_type*>& path) {
            invalidate_children();
            return true;
        }


        // Access collection of children facades.  Const version.
        const children_type& children() const {
            return const_cast<const children_type&>(const_cast<self_type*>(this)->children());
        }

        // Non-const version.  
        children_type& children() {
            if (m_children.empty() || m_children.size() != nchildren()) {
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

        // Remove and destroy node of kid and thus kid itself.
        void destroy_child(child_type*& kid, bool notify_value=true) {
            remove_child(*kid, notify_value); // do not catch return.
            kid = nullptr;
        }

        // Make a new child, returning its facade.
        child_type& make_child(bool notify_value=true) {
            invalidate_children();
            node_type* cnode = this->m_node->insert(notify_value);
            cnode->value.set_facade(std::make_unique<child_type>());
            return *cnode->value.template facade<child_type>();
        }

        // Sort children according to a comparison between children facades.
        // The comparison function must be callable like:
        //
        //     compare(const child_type* a, const child_type* b)
        // 
        // See also NaryTree::Node::sort_children() which this wraps..
        template<typename Compare>
        void sort_children(Compare comp, bool notify_value=true) {
            
            this->node()->sort_children([&](const node_ptr& na, const node_ptr& nb) {
                const child_type* fa = na->value.template facade<child_type>();
                const child_type* fb = nb->value.template facade<child_type>();
                return comp(fa, fb);
            }, notify_value);
        }

        // Apply the "separate" operation to the kid's children (a subset of our
        // grandchildren).  Return a map from a non-negative group number to a
        // newly created children facades (kid's siblings).
        //
        // The "groups" vector associates a "group ID" to each child of kid.
        // Negative group IDs will leave the corresponding kid's child in place.
        // Otherwise a child of kid will be added to a new parent and that
        // parents facade will be provided in the returned map at the
        // corresponding group ID key.
        //
        // If remove is true, then the kid's node will be removed from our
        // children nodes.  The kid facade is invalidated and the kid ptr nullified.
        using child_group_type = std::map<int, child_type*>;
        child_group_type separate(child_type*& kid,
                                  const std::vector<int> groups,
                                  bool remove=false,
                                  bool notify_value=true) {

            child_group_type ret;

            if (!kid) {
                raise<ValueError>("NaryTree::Facade::separate given a null kid");
            }
            if (!this->m_node->owns_child(kid->node())) {
                raise<ValueError>("NaryTree::Facade::separate given a kid not our own");
            }

            auto nurseries = kid->node()->separate(groups, notify_value);
            for (auto& nit : nurseries) {

                // Make a new sibling node.
                node_type* node = this->m_node->insert(notify_value);
                // Give the node its new children.
                node->adopt_children(nit.second, notify_value);
                // Get the bare facade pointer for return.
                child_type* facade = node->value.template facade<child_type>();
                ret[nit.first] = facade;

            }

            if (remove) {
                this->remove_child(*kid, notify_value);
                kid = nullptr;
            }

            return ret;
        }

        // Merge our grandchildren held by the given subset of our children
        // become children of target, creating target if it is nullptr.  Remove
        // the now empty children unless keep is true.  Return a "connected
        // components" (CC) type array with each element corresponding to a
        // child in target.  Existing children in target have group ID.  Any
        // children of the merged subset have a group ID set by the order of the
        // collection.  The existing ID is used to mark any existing children in
        // target in the returned CC array.
        // 
        // The children and the target may belong to another parent.
        template<typename ChildIt>
        std::vector<int> merge(ChildIt cbeg, ChildIt cend,
                               child_type* target=nullptr,
                               bool keep=false, int existingID=0) {

            std::vector<int> cc;

            int groupid = existingID;
            if (target) {
                cc.resize(target->nchildren(), groupid);
            }
            else {
                target = &make_child();
            }
            for (ChildIt cit = cbeg; cit != cend; ++cit) {
                ++groupid;
                child_type* child = *cit;
                cc.resize(cc.size() + child->nchildren(), groupid);
                target->take_children(*child);
                if (keep) {
                    continue;
                }
                // we need not be the parent.
                auto pnode = child->node()->parent;
                if (pnode) {
                    auto pfac = pnode->value.template facade<self_type>(); 
                    if (pfac) {
                        pfac->destroy_child(*cit);
                    }
                }
            }
            return cc;            
        }

        // Specialize for children provided in a map.  The group IDs for
        // grandchildren of children are provided explicitly by the map keys.
        // It is up to the user to assure none coincide with existingID.
        std::vector<int> merge(
            typename child_group_type::iterator cbeg,
            typename child_group_type::iterator cend,
            child_type* target=nullptr,
            bool keep=false, int existingID=0) {

            std::vector<int> cc;

            if (target) {
                cc.resize(target->nchildren(), existingID);
            }
            else {
                target = &make_child();
            }
            for (typename child_group_type::iterator cit = cbeg; cit != cend; ++cit) {
                cc.resize(cc.size() + cit->second->nchildren(), cit->first);
                target->take_children(*(cit->second));
                if (keep) {
                    continue;
                }
                // we need not be the parent.
                auto pnode = cit->second->node()->parent;
                if (pnode) {
                    auto pfac = pnode->value.template facade<self_type>(); 
                    if (pfac) {
                        pfac->destroy_child(cit->second);
                    }
                }
            }
            return cc;            
        }

        // A range version of the iterator-based methods above. 
        template<typename ChildrenRange>
        std::vector<int> merge(ChildrenRange& cr, child_type* target=nullptr, bool keep=false, int existingID=0) {
            return merge(std::begin(cr), std::end(cr), target, keep);
        }

        void invalidate_children() {
            m_children.clear();
        }

    private:
        // Lazy cache of children facades.
        children_type m_children;


    };

}
#endif
