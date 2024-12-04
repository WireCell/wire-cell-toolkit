#ifndef WIRECELLUTIL_NARYTREENOTIFIED
#define WIRECELLUTIL_NARYTREENOTIFIED

#include "WireCellUtil/NaryTree.h"

namespace WireCell::NaryTree {

    // A base class for a Node Value type that helps dispatch actions
    // to an inheriting subclass.  See NaryTesting::Introspective for
    // an example.
    template<typename Value>
    class Notified {
      public:

        using node_type = Node<Value>;

        virtual ~Notified() {}

      protected:

        // Subclasses should implement at least one of these protected
        // methods to recieve notification.


        // Called when a Node is constructed on a Notified.
        virtual void on_construct(node_type* node) {
        }

        // Called when a Node with a Notified is inserted.  The path
        // holds a sequence of Nodes starting with the inserted node
        // and ending with the Node holding the Notified being called.
        // Return true to continue propagating toward the root node.
        virtual bool on_insert(const std::vector<node_type*>& path) {
            return true;
        }

        // Called when a Node with a Notified is removed.  The path
        // holds a sequence of Nodes starting with the removed node
        // and ending with the Node holding the Notified being called. 
        // Return true to continue propagating toward the root node.
        virtual bool on_remove(const std::vector<node_type*>& path) {
            return true;
        }

        // Called when a Node with a Notified has its children ordered.  The
        // path holds a sequence of Nodes starting with the ordered node and
        // ending with the Node holding the Notified being called.  Return true
        // to continue propagating toward the root node.
        virtual bool on_ordered(const std::vector<node_type*>& path) {
            return true;
        }
      public:

        // This is the hook that Node will call.  Note, Node will use template
        // tests to determine if the this method exists in the Value type.
        //
        // A subclass may override this, for example to intercept all
        // notifications.
        virtual void notify(node_type* node, Action action) {
            // std::cerr << "NaryTree::Notified::notify(" << (void*)node << "," << action << ")\n";
            if (action == Action::constructed) {
                on_construct(node);
                return;
            }
            std::vector<node_type*> path = { node };
            if (action == Action::inserted) {
                propagate_insert_(path);
                return;
            }
            if (action == Action::removing) {
                propagate_remove_(path);
                return;
            }
            if (action == Action::ordered) {
                propagate_ordered_(path);
                return;
            }
        }

      private:
        void propagate_insert_(std::vector<node_type*> path) {
            if (! on_insert(path)) return; // notify subclass
            node_type* node = path.back();
            if (!node->parent) return;
            path.push_back(node->parent);
            node->parent->value.propagate_insert_(path);
        }

        void propagate_remove_(std::vector<node_type*> path) {
            if (! on_remove(path)) return; // notify subclass
            node_type* node = path.back();
            if (!node->parent) return;
            path.push_back(node->parent);
            node->parent->value.propagate_remove_(path);
        }

        void propagate_ordered_(std::vector<node_type*> path) {
            if (! on_ordered(path)) return; // notify subclass
            node_type* node = path.back();
            if (!node->parent) return;
            path.push_back(node->parent);
            node->parent->value.propagate_ordered_(path);
        }
    };
}
#endif
