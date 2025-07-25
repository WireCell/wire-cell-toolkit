#ifndef WIRECELL_CLUS_PR_GRAPH
#define WIRECELL_CLUS_PR_GRAPH

#include "WireCellClus/PRCommon.h"
#include "WireCellUtil/Graph.h"

namespace WireCell::Clus::PR {
    class Vertex;
    class Segment;

    using VertexPtr = std::shared_ptr<Vertex>;
    using SegmentPtr = std::shared_ptr<Segment>;

    // Where possible, we always use the term "node" inside the PR:: namespace
    // to refer to the "graph vertex" so as to avoid ambiguity with the
    // PR::Vertex class.
    struct NodeProperty {
        VertexPtr vertex;
    };
    struct EdgeProperty {
        SegmentPtr segment;
    };
    
    /** The PR graph types.
     *
     * These make use of BGL's "subgraph" type.  This is new usage for WCT so I
     * give some brief "rules" before describign details of this particular graph:
     *
     * - We only use the subgraph type, never the base graph type
     *
     * - We can make a subgraph from a subgraph
     *
     * - A subgraph has a parent() method returning the subgraph it was derived
     *   from.  A root() returns the utmost parent.  sg.children() return a pair
     *   of iterators for child subgraphs of current subgraph sg.
     *
     * - The subgraph and the root use different descriptors to refer to the
     *   same node or edge.  Both see "dense" descriptor values, no gaps.
     *
     * - We may add node/edge to a subgraph using global descriptors.
     *
     * - We get back from a subgraph (eg from add_vertex()) local descriptors.
     *
     * - We use sg.local_to_global(local_desc) -> global_desc and vice versa
     *   global_to_local() to convert both node and edge descriptors.
     *
     * A primary feature of "subgraph" is that each subgraph (including parents)
     * will , reflect any change made via any subgraph.  Thus a holder of the
     * subgraph does not need to maintain/sync any other data structure to know
     * "its" nodes/edges.  That is exactly what the subgraph is.
     *
     *
     * We now describe details about this particular (sub) graph type.
     * 
     * This graph uses setS for descriptor containers.
     *
     * This allows descriptors to be stable on remove and add.  Of course,
     * removing a descriptor invalidates it.  Removing a node descriptor
     * invalidates all associated edge descriptors.
     *
     * The sets are ordered by descriptors and iterating on nodes or edges is
     * supposed to be stable on rerunning of identical jobs (eg, no "sort on
     * pointer problem" is expected).
     *
     * Users should NOT use boost:: methods to add/remove nodes and edges and
     * instead use the provided wrappers.  The wrappers will manage the (global)
     * graph descriptor held by the PR::Vertex and PR::Segment.
     */
    // We never make instances of this type.
    using BaseGraph = boost::adjacency_list<
        boost::setS,            // vertices
        boost::setS,            // edges
        boost::undirectedS,     // edge direction
        boost::property<boost::vertex_index_t, int, NodeProperty>,
        boost::property<boost::edge_index_t, int, EdgeProperty>
        >;

    /// The (sub) graph type for PR.
    using graph_type = boost::subgraph<BaseGraph>;

    /// The node descriptor type
    using node_descriptor_type = boost::graph_traits<graph_type>::vertex_descriptor;
    /// The edge descriptor type
    using edge_descriptor_type = boost::graph_traits<graph_type>::edge_descriptor;

    /** The graph mutation functions generally take pointers to PR::Vertex to
     * indicate a node and PR::Segment to indicate an edge.  The descriptor held
     * by these objects are ALWAYS GLOBAL.  A descriptor returned from a
     * function is ALWAYS LOCAL.
     */
    
    /** Add a vertex as a node in the graph.
     *
     * Return the node descriptor.
     *
     * If the (global) Vertex::descriptor is not invalid, the corresponding
     * local descriptor is simply returned and the graph is not modified.
     * Otherwise a new node is added to the graph, its global descriptor saved
     * to the PR::Vertex and the local descriptor returned.
     *
     * User caution: the (global) Vertex descriptor is only checked for
     * validity.  A user should not try to circumvent these functions as this
     * may result in a valid descriptor that is not actually in the graph.
     */
    node_descriptor_type add_vertex(VertexPtr vertex, graph_type& graph);

    /** Add a segment as an edge between two vertex nodes.
     *
     * Return the local node descriptor.
     *
     * The vertices will also be added if they are not yet in the graph.
     * 
     * See add_vertex() for other applicable commentary.
     */    
    edge_descriptor_type add_segment(VertexPtr vertex1,
                                     VertexPtr vertex2,
                                     SegmentPtr segment, graph_type& graph);
    
    /** Remove the vertex's node from the graph
     *
     * This will invalidate the vertex's descriptor and remove all the node's
     * edges and invalidate their segment's descriptors.
     *
     * If vertex's descriptor is invalid, this is a no-op.  If vertex's
     * descriptor is not in the graph, the graph is not modified and the
     * vertex's descriptor is made invalid.
     */
    void remove_vertex(VertexPtr vertex, graph_type& graph);

    /** Remove the segment's edge from the graph.
     *
     * The segment's descriptor will be set to invalid.  If already invalid,
     * this function is a no-op.  If segment's descriptor is not in graph, the
     * graph is not modified and the descriptor is invalidated.
     */
    void remove_segment(SegmentPtr segment, graph_type& graph);



    /** Mixin for Vertex/Segment to manage the descriptor.
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


};
#endif
