/** Classes and functions related to a "view" of a trajectory graph.
 *
 */

#ifndef WIRECELL_CLUS_PR_TRAJECTORYVIEW
#define WIRECELL_CLUS_PR_TRAJECTORYVIEW

#include "WireCellClus/PRGraphType.h"

#include <unordered_set>


namespace WireCell::Clus::PR {

    // Forward declare.
    class TrajectoryView;

    // Internal, used to select nodes in a view
    struct TrajectoryViewNodePredicate {
        const TrajectoryView& view;
        TrajectoryViewNodePredicate(const TrajectoryView& v) : view(v) {};
        bool operator()(const node_descriptor& desc) const;
    };
    struct TrajectoryViewEdgePredicate {
        const TrajectoryView& view;
        TrajectoryViewEdgePredicate(const TrajectoryView& v) : view(v) {};
        bool operator()(const edge_descriptor& desc) const;
    };


    /** Provide a view of a trajectory graph.
     *
     * A PR::TrajectoryView manages a filtered (not necessarily induced)
     * subgraph on the given trajectory graph.  The TrajectoryView is
     * essentially a wrapper over `boost::filtered_graph` with an API that is in
     * terms of PR::Vertex and PR::segment instead of bare descriptors.
     *
     * A PR::TrajectoryView is given and holds a borrowed reference to a user's
     * graph instance and provides methods to manage a subset of the graph's
     * nodes and edges that are considered "in view".
     *
     * It is recommended that the base PR::TrajectoryView be used by application
     * code.  However, it is allowed to extend PR::TrajectoryView via inheritance.
     * 
     * A PR::TrajectoryView may be created in isolation.  See
     * `make_trajectoryview()` to produce one as a shared pointer and
     * PR::Trajectory for a convenient way to construct and maintain a
     * collection of (potentially heterotypical) TrajectoryView instances.
     */
    class TrajectoryView {
    public:

        using full_graph_type = Graph;

        using view_graph_type = boost::filtered_graph<Graph, TrajectoryViewNodePredicate, TrajectoryViewEdgePredicate>;
        using view_node_descriptor = boost::graph_traits<view_graph_type>::vertex_descriptor;
        using view_edge_descriptor = boost::graph_traits<view_graph_type>::edge_descriptor;

        TrajectoryView() = delete;
        virtual ~TrajectoryView();

        /** Construct a view on a borrowed reference a trajectory graph.
         *
         * Caller must assure graph object lifetime exceeds that of the view.
         */ 
        TrajectoryView(full_graph_type& graph);

        /** Access to the underlying filtered graph. */
        const view_graph_type& view_graph() const;

        /// Return true if node descriptor is in the filter.
        bool has_node(node_descriptor desc) const;

        /// Return true if edge descriptor is in the filter.
        bool has_edge(edge_descriptor desc) const;

        /** Add the vertex to the filter.

            To be added to the filter the vertex must be already in the
            underlying graph as determined by it having a valid descriptor.

            Return true if view is modified.
         */
        bool add_vertex(VertexPtr vtx);

        /** Add the segment to the filter.

            To be added to the filter the segment must be already in the
            underlying graph as determined by it having a valid descriptor.

            Return true if view is modified.
         */
        bool add_segment(SegmentPtr seg);

        /** Remove the vertex from the filter.

            Return true if view is modified.
         */
        bool remove_vertex(VertexPtr vtx);

        /** Remove the segment from the filter.

            Return true if view is modified.
         */
        bool remove_segment(SegmentPtr seg);

        /** Get the set of node descriptors in this view. */
        const node_unordered_set& nodes() const { return m_nodes; }

        /** Get the set of edge descriptors in this view. */
        const edge_unordered_set& edges() const { return m_edges; }

    private:
        view_graph_type m_graph;
        node_unordered_set m_nodes;
        edge_unordered_set m_edges;
    };

    /** Make a view of a particular type as a shared pointer.
     *
     */
    template<typename ViewType, typename... Args>
    std::shared_ptr<ViewType> make_trajectoryview(Graph& graph, Args&&... args) {
        return std::shared_ptr<ViewType>(graph, std::forward<Args>(args)...);
    };

}

#endif
