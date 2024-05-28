#ifndef WIRECELL_PGRAPH_PORT
#define WIRECELL_PGRAPH_PORT

#include "WireCellUtil/Exceptions.h"

#include <boost/any.hpp>

#include <string>
#include <vector>
#include <deque>
#include <memory>

namespace WireCell {
    namespace Pgraph {

        // The type of data passed in the graph.
        typedef boost::any Data;
        // A buffer of data.  It is a std::deque instead of a
        // std::queue so that it may be iterated as well as pushed
        // back into.  Like a British queue, one enters it from the
        // back and exits it from the front.
        typedef std::deque<Data> Queue;

        // Edges are just queues that can be shared.
        typedef std::shared_ptr<Queue> Edge;

        class Node;

        class Port {
           public:
            enum Type { tail = 0, output = 0, head = 1, input = 1, ntypes = 2 };

            Port(Node* node, Type type, std::string signature, std::string name = "");

            bool isinput() const;
            bool isoutput() const;
            Edge edge() const;

            // Connect an edge, returning any previous one.
            Edge plug(Edge edge);

            // return edge queue size or 0 if no edge has been plugged
            size_t size() const;

            // Return true if queue is empty or no edge has been plugged.
            bool empty() const;

            // Get the next data.  By default this pops the data off
            // the queue.  To "peek" at the data, pas false.
            Data get(bool pop = true);

            // Put the data onto the queue.
            void put(Data& data);

            // Get back the associated Node.
            Node* node();
            const Node* node() const;

            const std::string& name() const;
            const std::string& signature() const;

            std::string ident() const;

           private:
            Node* m_node;  // node to which port belongs
            Type m_type;
            std::string m_name, m_sig;
            Edge m_edge;
        };

        typedef std::vector<Port> PortList;
    }  // namespace Pgraph
}  // namespace WireCell

#endif
