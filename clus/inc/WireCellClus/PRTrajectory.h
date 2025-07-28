#ifndef WIRECELL_CLUS_PR_TRAJECTORY
#define WIRECELL_CLUS_PR_TRAJECTORY

#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRGraph.h"

namespace WireCell::Clus::PR {


    /** A PR::Trajectory is a base class for a model of a particle trajectory.
     *
     * Conceptually, a "trajectory" represents a set of "particle paths".
     *
     * The set of "particle paths" is represented as a graph (aka a "trajectory
     * graph").  A graph node represents a point in space by providing a
     * property of type PR::Vertex.  We will try to always say "graph node" or
     * "node" when referring to a trajectory graph "vertex" to make that
     * distinct from the PR::Vertex property that the node holds.  Similarly, a
     * trajectory graph edge has a property of type PR::Segment representing the
     * traversal of a particle between two PR::Vertex instances.  A PR::Vertex
     * may simply represent the connection of two PR::Segment instances or
     * additionally represent interaction vertices with any number of connected
     * PR::Segment edges.
     *
     * See PR::Vertex and PR::Segment for more information on these types.
     *
     * While PR::Trajectory type can be used directly, the application may
     * derived subclasses (eg PR::Shower, PR::Track, PR::Interaction or
     * PR::Neutrino).  These may specialize the PR::Track by overriding methods
     * and/or providing different interpretations for flags.
     *
     * A set of PR::Trajectory instances form a tree structure via their graphs.
     * That is, a PR::Trajectory may be the "child" of another "parent"
     * PR::Trajectory.  Children represent a "view" of a subset of the paths
     * viewed in the parent.  This view provided to the Trajectory object by its
     * Boost "subgraph" that is held by reference by the PR::Trajectory.  The
     * "subgraph" holds a weak pointer reference back to the PR::Trajectory that
     * is its "holder" in order to facilitate navigation between graph and
     * object representations.  It is the Boost subgraphs that implement the
     * tree structure and that tree structure is reflected through
     * PR::Trajectory methods and related free functions.
     *
     * Use of a trajectory and graph tree begins with the user creating the
     * "root" graph and then providing its reference to a "root" trajectory.
     * From then on, "child" trajectories and graphs can be generated.
     *
     * Note, a Trajectory can ONLY be created via shared pointer.  This is in
     * order for supply the weak pointer back reference held by the graph.  Use
     * make_trajectory() or make_child_trajectory() to make trajectory objects.
     */

    // FIXME: WCShower has flags, should make Trajectory a "Flagged" or have it
    // hold a "Flagged".  But, how to have the enum class used in Flagged<T> to
    // be extended?

    class Trajectory : public std::enable_shared_from_this<Trajectory> {
        
    protected:

        // A subclass may have a custom constructor but MUST forward the graph reference

        // Construct a trajectory with a reference to "its" subgraph.  The root
        // subgraph is owned by the user and any child subgraphs are owned by
        // the parent subgraph.
        Trajectory(graph_type& graph);

        // Called by to form internal references, setting this on graph property
        // bundle.  This is post-construct to assure full shared pointer
        // construction is complete.
        void post_construct_init();

    public:

        // Subclassing Trajectory is expected.
        virtual ~Trajectory();

        // Delete copy and move constructors to prevent other than shared
        // pointer access.
        Trajectory() = delete;
        Trajectory(const Trajectory&) = delete;
        Trajectory(Trajectory&&) = delete;


        /** Access the underlying (sub) graph of this trajectory.
         *
         * The user should take care not to modify this graph unless the
         * subgraph rules are strictly followed.  Instead, consider functions
         * that operate on PR objects instead.
         *
         * See get_trajectory(graph) for a way to recover this trajectory from
         * its graph.
         */
        graph_type& graph() { return m_graph; }
        const graph_type& graph() const { return m_graph; }
        

        /** Add vertex to this trajectory.
         *
         * Return true if the graph was modified.
         *
         * If the descriptor held by the PR::Vertex is invalid, a new graph node
         * is added to this trajectory's subgraph with the VertexPtr as a
         * property.  The descriptor held by the PR::Vertex is set to the just
         * produced (local) node descriptor.  True is returned.
         *
         * If the descriptor held by the PR::Vertex is already valid, neither
         * graph nor PR::Vertex are modified and false is returned.
         */
        bool add_vertex(VertexPtr vtx);

        /** Add segment as edge between two vertices to this trajectory.
         *
         * Return true if the graph was modified.
         *
         * This adds (or does not add) following the same descriptor rules as
         * described in add_vertex().  If vertices are not yet added, they will
         * be.
         *
         * False is return if segment and vertices already had valid
         * descriptors.
         */
        bool add_segment(SegmentPtr seg, VertexPtr vtx1, VertexPtr vtx2);

        /* Remove vertex from this trajectory.
         *
         * Return true if the graph was modified.
         *
         * If the descriptor held by the PR::Vertex is valid, remove it from the
         * trajectory's subgraph.  This may lead to segment edges being removed.
         * Removals are reflected in all subgraphs.
         *
         * If the descriptor held by the PR::Vertex is invalid, return false.
         */ 
        // bool remove_vertex(VertexPtr vtx);
        // disabled as subgraph does NOT SUPPORT VERTEX removal.  edge removal is okay.

        /** Remove segment from this trajectory.
         *
         * Return true if the graph was modified.
         *
         * If the descriptor held by the PR::Segment is valid, remove it from
         * the trajectory's subgraph.  This does not necessarily lead to the
         * edge endpoint vertices to be removed.  Removals are reflected in all
         * subgraphs.
         *
         * If the descriptor held by the PR::Segment is invalid, return false.
         */ 
        bool remove_segment(SegmentPtr seg);


        /** Construct a trajectory (or subclass) as shared pointer.
         *
         * Arguments must match the corresponding constructor.
         *
         * This can be called to make a root or a child trajectory but see
         * make_child_trajectory() for a more friendly way to make a child
         * trajectory.
         */
        template <typename T, typename... Args>
        friend std::shared_ptr<T> make_trajectory(Args&&... args);

        /** Construct a child trajectory (or subclass) as shared pointer.
         *
         * This variant will handle the creation of the child subgraph to give
         * to the new trajectory.
         */
        template <typename T_Child, typename T_Parent, typename... Args>
        friend std::shared_ptr<T_Child> make_child_trajectory(
            std::shared_ptr<T_Parent> parent, Args&&... args);


    private:
        graph_type& m_graph;
    };

    /** Make a shared pointer instance of Trajectory type or one of its subclases.
     *
     * Arguments must match the constructor for the type being created.
     *
     * This will typically be used to make the singular "root" trajectory
     * instance where the root subgraph is held by the caller.  It may also be
     * used to create a child trajectory but then the user must properly provide
     * the reference to the (sub)graph as the first argument by calling
     * ubgraph::create_subgraph().  DO NOT ALLOW THE SUBGRAPH TO COPY.  See
     * make_child_trajectory() for a more convenient factory.
     */
    template <typename T, typename... Args>
    std::shared_ptr<T> make_trajectory(graph_type& sg, Args&&... args) {
        // Compile-time check to ensure T is Trajectory or a subclass
        static_assert(std::is_base_of<Trajectory, T>::value,
                      "make_trajectory can only create instances of Trajectory or its subclasses.");
        auto ptr = std::make_shared<T>(sg, std::forward<Args>(args)...);
        std::static_pointer_cast<Trajectory>(ptr)->post_construct_init();
        return ptr;
    }

    /** Make a trajectory that is a child of another trajectory.
     *
     * The first argument is a parent trajectory which provides the parent
     * subgraph from which a new child subgraph is made and passed to the new
     * child trajectory.  The remaining arguments must match the non-graph
     * arguments to the "child" style of constructor.
     */
    template <typename T_Child, typename T_Parent, typename... Args>
    std::shared_ptr<T_Child> make_child_trajectory(
        std::shared_ptr<T_Parent> parent, Args&&... args) {

        static_assert(std::is_base_of<Trajectory, T_Child>::value,
                      "T_Child must be Trajectory or a subclass.");
        static_assert(std::is_base_of<Trajectory, T_Parent>::value,
                      "T_Parent must be Trajectory or a subclass.");
        if (!parent) {
            throw std::invalid_argument("Parent shared_ptr cannot be null when creating a child trajectory.");
        }

        graph_type& child_graph = parent->get_graph().subgraph();
        auto ptr = std::make_shared<T_Child>(child_graph, std::forward<Args>(args)...);
        std::static_pointer_cast<Trajectory>(ptr)->post_construct_init();
        return ptr;
    }
    

    /** Access the PR::Trajectory that may be stored on a subgraph.
     *
     * The shared pointer will be nullptr when the subgraph lacks a trajectory.
     */
    TrajectoryPtr get_trajectory(graph_type& graph);



}

#endif
