/** The PR trajectory graph.
 *
 * This makes use of boost subgraph to produce a tree of trajectories, each with
 * its own view of the underlying total graph/trajectory.
 *
 * We try everywhere to say "graph node" (not "graph vertex") to avoid confusion
 * with the face that each graph node holds a PR::Vertex.  A graph edge holds a
 * PR::Segment.
 */

#ifndef WIRECELL_CLUS_PR_GRAPH
#define WIRECELL_CLUS_PR_GRAPH

#include "WireCellClus/PRCommon.h"
#include "WireCellUtil/Graph.h"

namespace WireCell::Clus::PR {

    /// The headers for these classes include this file, so we just forward-declare.
    class Vertex;
    class Segment;
    class Trajectory;

    /// A graph node property bundle holds the shared pointer to PR::Vertex.
    /// 
    /// Note, PR::Vertex holds a descriptor for its graph node to allow
    /// navigation between graph and object representations.
    using VertexPtr = std::shared_ptr<Vertex>;
    struct NodeProperty {
        VertexPtr vertex;       // shared pointer
    };

    /// A graph edge property bundle holds the shared pointer to PR::Segment.
    ///
    /// Note, PR::Segment holds a descriptor for its graph edge to allow
    /// navigation between graph and object representations.
    using SegmentPtr = std::shared_ptr<Segment>;
    struct EdgeProperty {
        SegmentPtr segment;     // shared pointer
    };

    /// A graph level property bundle holds a weak pointer to PR::Trajectory.
    ///
    /// Note, PR::Trajectory holds/owns its subgraph by value.  We must use a
    /// weak pointer here to avoid an ownership cycle.  You must convert the
    /// weak pointer to a shared pointer before use, and that shared pointer
    /// will be nullptr if (somehow) the underlying PR::Trajectory has been
    /// destroyed.  See get_trajectory(graph_type) helper in Trajectory.h.
    using TrajectoryPtr = std::weak_ptr<Trajectory>;
    struct GraphProperty {
        TrajectoryPtr trajectory; // weak pointer
    };
    
    /** Boost Graph Library: subgraph
     *
     * - A subgraph is a "view" of a particular subset or of vertices and edges
     *   of an underlying graph.
     *
     * - Modifying a subgraph changes the underlying graph.
     *
     * - A change to the underlying graph is reflected in all subgraph views
     *   that previously included any of the effected vertices and/or edges.
     *
     * - The C++ "subgraph" type is templated on another "concrete" aka "base"
     *   aka "underlying" BGL graph type.
     *
     * - C++ object instances of the "subgraph" type are used.  Instances of the
     *   underlying graph type do not participate in any subgraph functions.
     *
     * - A single instance, labeled here as "sg", is considered the "root"
     *   subgraph and represents a complete view of the entire underlying graph.
     *
     * - Child sub-graphs of "sg" may be derived, labeled here as "sg0", "sg1", ...
     *
     * - The "sg" instance is the "parent" of the "sgN" "children".
     *
     * - The children may be parents, thus we may create instances "sg00, sg01,
     *   ..., sg10, sg11, ...".
     *
     * - Together, the subgraph instances form a "subgraph tree".
     *   
     * - Every subgraph instance has "local descriptors".  These are "dense"
     *   (counts, 0, 1, 2, ...)
     *
     * - Every subgraph maintains a map from its "local descriptors" to its
     *   "global descriptors" which are identical to the subgraph's parent's
     *   local descriptors.
     *
     * - A subgraph is always created as an empty view.
     *
     * - We add (add_vertex(), add_edge()) "global descriptors" to a subgraph
     *   and "local descriptors" are returned.
     *   
     * - We may use local_to_global(local_desc) -> global_desc and vice versa
     *   global_to_local() to convert descriptors.
     *
     * - A subgraph has parent() and root() methods for getting immediate and
     *   ultimate parents, respectively.
     *
     * The primary feature of "subgraph" is that any change to the underlying
     * graph via any subgraph will reflect into the view of ever subgraph.
     *
     */
    // fixme: move the above blurb someplace more generic/prominent.
    /**
     *
     * We now describe details about this particular (sub) graph type.
     * 
     * This graph uses setS for descriptor containers.
     *
     * This allows descriptors to be stable on remove and add.  Of course,
     * removing a descriptor invalidates it.  Removing a node descriptor
     * invalidates all associated edge descriptors.
     *
     * The setS vertex descriptor is void* so iteration order is NOT stable.
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
        boost::property<boost::edge_index_t, int, EdgeProperty>,
        GraphProperty
        >;

    /// The (sub) graph type for PR.
    using graph_type = boost::subgraph<BaseGraph>;

    /** You may use standard boost::add_vertex(), etc on the graph_type (as subgraph).
     *
     * However, you MUST obey the rules of boost subgraph.  Namely, you must
     * pass "global" descriptors defined on the parent subgraph when operating
     * on a child subgraph.  Any descriptor returned is a "local" descriptor
     * defined in the context of the child subgraph.
     *
     * See Trajectory.h for equivalent functions in terms of PR:: types that
     * will take care to handle the parent/child global/local rules.
     */

    /** Boost subgraph is (as of 1.85) broken when it comes to getting the graph
     * bundle in the usual way.  This function provides a work-around.
     */
    GraphProperty& get_graph_bundle(graph_type& g);
    const GraphProperty& get_graph_bundle(const graph_type& g);


    /// The node descriptor type
    using node_descriptor_type = boost::graph_traits<graph_type>::vertex_descriptor;
    /// The edge descriptor type
    using edge_descriptor_type = boost::graph_traits<graph_type>::edge_descriptor;

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


};
#endif
