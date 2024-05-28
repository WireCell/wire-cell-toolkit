/** Helpers for testing NaryTree. */

#ifndef WIRECELLUTIL_NARYTESTING
#define WIRECELLUTIL_NARYTESTING

#include "WireCellUtil/NaryTreeNotified.h"

namespace WireCell::NaryTesting {

    // Simple node data that tracks its copies/moves and receives
    // notification of "actions" from the node.
    struct Introspective : public NaryTree::Notified<Introspective>
    {

        using node_type = NaryTree::Node<Introspective>;
        using owned_ptr = std::unique_ptr<node_type>;

        static size_t live_count;

        std::string name;

        size_t ndef{0}, nctor{0}, ncopy{0}, nmove{0};
        std::map<std::string, size_t> nactions;

        node_type* node{nullptr};

        Introspective();
        Introspective(const Introspective& d);
        Introspective(Introspective&& d);
        explicit Introspective(const std::string& s);
        virtual ~Introspective();

        // example of receiving notifications and propagating them up the
        // parentage.
        virtual void on_construct(node_type* node);
        virtual bool on_insert(const std::vector<node_type*>& path);
        virtual bool on_remove(const std::vector<node_type*>& path);


    };

    std::ostream& operator<<(std::ostream& o, const Introspective& obj);

    // Make a simple tree using all insert methods, return the root node.
    Introspective::owned_ptr make_simple_tree(const std::string& base_name = "0");

    // Recursively add children to nodes with layer_sizes giving
    // number of children per node in each layer.
    void add_children(Introspective::node_type* node, std::list<size_t> layer_sizes);

    // Generate and return root node of a tree with the given layer sizes. 
    Introspective::owned_ptr make_layered_tree(std::list<size_t> layer_sizes);

}


#endif
