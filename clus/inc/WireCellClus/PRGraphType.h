/** Define a "trajectory" graph type.

    A trajectory graph's nodes have an associated PR::Vertex and edges have an
    associated PR::Segment.  These objects also carry their associated graph
    descriptors.

    Application code should avoid adding nodes and edges directly using
    `boost::add_vertex()` and `boost::add_edge()`.  See PRGraph.h for equivalent
    methods that operate on the PR::Vertex and PR::Segment objects and assure
    proper bookkeeping.

    Application code that is sensitive to the order of iterating over graph
    nodes should use the `PR::ordered_nodes()` function to get a stable order.
    Native ordering of nodes in the graph is subject to pointer value variance.
 */

#ifndef WIRECELL_CLUS_PR_GRAPHTYPE
#define WIRECELL_CLUS_PR_GRAPHTYPE

#include "WireCellUtil/Graph.h"

namespace WireCell::Clus::PR {

    // The headers for these classes include this header, so here we just
    // forward-declare to avoid a cycle.

    /// The node object type for all graphs
    class Vertex;
    /// The shared pointer to the graph node object.
    using VertexPtr = std::shared_ptr<Vertex>;
    /// A graph node property bundle holds the shared pointer to PR::Vertex.
    /// 
    /// Note, PR::Vertex holds a descriptor for its graph node to allow
    /// navigation between graph and object representations.
    struct NodeBundle {
        /// A shared pointer to a PR::Vertex object
        VertexPtr vertex;       // shared pointer
        /// A monotonically increasing, potentially sparse index uniquely
        /// identifying this edge in a graph
        size_t index;
    };


    /// The edge object type for all graphs.
    class Segment;
    /// The shared pointer to the graph edge object.
    using SegmentPtr = std::shared_ptr<Segment>;
    /// A graph edge property bundle holds the shared pointer to PR::Segment.
    ///
    /// Note, PR::Segment holds a descriptor for its graph edge to allow
    /// navigation between graph and object representations.
    struct EdgeBundle {
        /// A shared pointer to a PR::Segment object
        SegmentPtr segment;     // shared pointer
        /// A monotonically increasing, potentially sparse index uniquely
        /// identifying this edge in a graph
        size_t index;   
    };

    /// A graph-level property holds internal book-keeping information.
    ///
    /// Normal application code need and should not access this.
    struct GraphBundle {
        /// The total number of nodes ever added to this graph.
        ///
        /// This is used to set a unique index.  It is greater or equal to the
        /// number of nodes currently in the graph.
        size_t num_node_indices{0};
        /// The total number of edges ever added to this graph.
        ///
        /// This is used to set a unique index.  It is greater or equal to the
        /// number of edges currently in the graph.
        size_t num_edge_indices{0};
        
    };
    
    /** The graph type.
     * 
     * This graph uses setS for descriptor containers.  It provides robust
     * node/edge descriptors that will remain valid when unrelated node/edge
     * descriptors are removed.  The order of iterating the raw vertex or edge
     * set is well defined within a given program execution.  However, this
     * order is based on pointer value and so may change between executions of
     * an identical program.  Use ordered_nodes() to produce an ordering based
     * on index.
     *
     * User should avoid using boost::add_vertex() and boost::add_edge() on
     * instances of this graph type.  Instead use PR::add_vertex() and
     * PR::add_segment().
     */
    using Graph = boost::adjacency_list<
        boost::setS,            // vertices
        boost::setS,            // edges
        boost::undirectedS,     // edge direction
        NodeBundle,
        EdgeBundle,
        GraphBundle
        >;

    /// The descriptor type for nodes.  This is a `void*`.
    using node_descriptor =  boost::graph_traits<Graph>::vertex_descriptor;

    /// A vector of node descriptors
    using node_vector = std::vector<node_descriptor>;
    
    /// The (user) iterator on nodes.  See `ordered_nodes()`.
    using node_iterator = node_vector::iterator;

    /// The iterator range.
    using node_range = boost::iterator_range<node_iterator>;

    /// The descriptor type for edges.
    using edge_descriptor =  boost::graph_traits<Graph>::edge_descriptor;
    /// The (user) iterator on edges.  See `ordered_edges()`.
    using edge_iterator = std::vector<edge_descriptor>::iterator;
    using edge_range = boost::iterator_range<edge_iterator>;


    /// Return a vector of node descriptors in native graph order.
    ///
    /// This order is based on pointer values.
    ///
    /// The vector may be conveniently used to iterate over nodes:
    ///
    /// @code{.cpp}
    /// for (const auto& vd : nodes(my_graph)) {
    ///     std::cout << "Node index: " << my_graph[vd].index << std::endl;
    /// }
    /// @endcode
    node_vector graph_nodes(Graph& g);

    /// Return a vector of node descriptors ordered by index.
    ///
    /// The vector may be conveniently used to iterate over nodes:
    ///
    /// @code{.cpp}
    /// for (const auto& vd : ordered_nodes(my_graph)) {
    ///     std::cout << "Node index: " << my_graph[vd].index << std::endl;
    /// }
    /// @endcode
    node_vector ordered_nodes(Graph& g);

    

    /** A mixin class for Vertex/Segment to manage their stored descriptor.
     */
    template <typename Descriptor>
    class Graphed {
    public:
        using descriptor_type = Descriptor;

        const descriptor_type invalid_descriptor{};

        descriptor_type get_descriptor() const { return m_descriptor; }
        void set_descriptor(descriptor_type descriptor) { m_descriptor = descriptor; }

        bool descriptor_valid() const {
            return m_descriptor != invalid_descriptor;
        }
        void invalidate_descriptor() {
            m_descriptor = invalid_descriptor;
        }

    private:
        descriptor_type m_descriptor{};
    };


}
#endif
