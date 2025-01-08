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

        // Dispatch a notification to corresponding "on_*()" method.
        //
        virtual bool notify(std::vector<node_type*> path, Action action) {

            if (path.size() == 1 && action == Action::constructed) {
                on_construct(path.back());
                return true;    // assume caller mediates propagation 
            }
            if (action == Action::inserted) {
                return on_insert(path);
            }
            if (action == Action::removing) {
                return on_remove(path);
            }
            if (action == Action::ordered) {
                return on_ordered(path);
            }
            return false;       // fixme:  add exception?  warning?
        }
        
    };
}
#endif
