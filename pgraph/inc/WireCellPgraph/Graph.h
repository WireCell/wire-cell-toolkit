/** A simple DFP engine meant for single-threaded execution.

    A node is constructed with zero or more ports.

    A port mediates between a node and an edge.

    A port is either of type input or output.

    In the context of a node, a port has a name and an index into an
    ordered list of other ports of a given type.

    An edge is a queue passing data objects in one direction from its
    tail (input) end to its head (output) end.

    Each end of an edge is exclusively "plugged" into one node through
    one port.

    A valid graph consists of nodes with all ports plugged to edges.

 */

#ifndef WIRECELL_PGRAPH_GRAPH
#define WIRECELL_PGRAPH_GRAPH

#include "WireCellPgraph/Node.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/ExecMon.h"

#include <vector>
#include <unordered_set>
#include <unordered_map>

namespace WireCell {
    namespace Pgraph {

        class Graph {
           public:
            Graph();

            // Add a node to the graph.
            void add_node(Node* node);

            // Connect two nodes by their given ports.  Return false
            // if they are incompatible.  new nodes will be implicitly
            // added to the graph.
            bool connect(Node* tail, Node* head, size_t tpind = 0, size_t hpind = 0);

            // return a topological sort of the graph as per Kahn algorithm.
            std::vector<Node*> sort_kahn();

            // Excute the graph until nodes stop delivering
            bool execute();

            // Excute parents of node or if any parent is not ready,
            // recursively call this method on parent.  Return number
            // of nodes executed.
            int execute_upstream(Node* node);

            // All internal calling of nodes goes through here.
            bool call_node(Node* node);

            // Return false if any node is not connected.
            bool connected();

            // Print out cumulated CPU time for executing each node
            void print_timers(bool include_execmon=false) const;

           private:
            std::vector<std::pair<Node*, Node*> > m_edges;
            std::unordered_set<Node*> m_nodes;
            std::unordered_map<Node*, std::vector<Node*> > m_edges_forward, m_edges_backward;
            Log::logptr_t l;
            Log::logptr_t l_timer;
            std::unordered_map<Node*, double> m_nodes_timer;
            ExecMon m_em;
        };
    }  // namespace Pgraph
}  // namespace WireCell
#endif
