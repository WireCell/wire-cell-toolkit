/** Define functions that can operate on a "trajectory" graph.

    Functions provided should be preferred over equivalent graph related
    functions in `boost::` as the trajectory graph requires properly bookkeeping
    when adding/removing nodes and/or edges.  It is undefined behavior if you
    call `boost::` graph mutators and fail to follow the bookkeeping
    conventions.

 */
#ifndef WIRECELL_CLUS_PR_GRAPH
#define WIRECELL_CLUS_PR_GRAPH

#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRSegment.h"

#include "WireCellUtil/Graph.h"
#include "WireCellUtil/Exceptions.h"

namespace WireCell::Clus::PR {


    /// Make a vertex and add it to the graph.
    ///
    /// This method will properly adhere to the indexing policy and is
    /// recommended instead of creating a PR::Vertex by hand.
    ///
    /// This is templated to allow for non-default construction of the
    /// PR::Vertex.  Normal usage:
    ///
    /// @code{.cpp}
    /// auto my_vtx = make_vertex(my_graph);
    /// @endcode
    template <typename... Args>
    VertexPtr make_vertex(Graph& g, Args&&... args) {
        auto ptr = std::make_shared<Vertex>(std::forward<Args>(args)...);
        bool ok = add_vertex(g, ptr);
        if (!ok) {
            raise<RuntimeError>("make_vertex failed to add vertex which should never happen!");
        }
        return ptr;
    }

    /// Add an existing vertex to the graph.
    ///
    /// Return true if graph modified.
    ///
    /// If the PR::Vertex already has a valid descriptor, false is returned.
    ///
    /// Otherwise, a new graph node is added with its bundle holding this vertex
    /// and a unique index.  The resulting descriptor is set on the vertex and
    /// true is returned.
    ///
    /// User must take responsibility to provide a vertex that is consistent
    /// with the indexing policy.  For a safer function, use `make_vertex()`.
    bool add_vertex(Graph& g, VertexPtr vtx);


    /// Remove a vertex from the graph.
    ///
    /// Return true if graph modified.
    ///
    /// The descriptor held by the object is left invalid.
    bool remove_vertex(Graph& graph, VertexPtr vtx);


    /// Make a PR::Segment
    ///
    /// This makes an "orphaned" segment not associated to any graph.  Call
    /// `add_segment()` to add it to a graph, or better, use `make_segment()`
    /// to do both at once.
    template <typename... Args>
    SegmentPtr make_segment(Args&&... args) {
        return std::make_shared<Segment>(std::forward<Args>(args)...);
    }

    /// Add a PR::Segment to the graph as an edge between the nodes of two
    /// PR::Vertex instances.
    ///
    /// Return true if the graph was modified.  Modification means that a new
    /// vertex or edge was added.  In the case that vtx1 and vtx2 already have
    /// an edge, it is updated with the segment and that does not influence
    /// true/false return value.
    ///
    /// If the PR::Segment already has a valid descriptor then it is not added.
    ///
    /// If either PR::Vertex do not have a valid descriptor then they are added
    /// to the graph and will cause a true return value irrespective of the
    /// status of the segment.  But, users are urged to only use `make_index()`
    /// to construct a PR::Vertex which will assure graph addition.
    bool add_segment(Graph& g, SegmentPtr seg, VertexPtr vtx1, VertexPtr vtx2);


    /// Make a PR::Segment and add it to the graph as an edge between two vertices.
    ///
    template <typename... Args>
    SegmentPtr make_segment(Graph& g, VertexPtr vtx1, VertexPtr vtx2, Args&&... args) {
        SegmentPtr seg = std::make_shared<Segment>(std::forward<Args>(args)...);
        add_segment(g, seg, vtx1, vtx2);
        return seg;
    }


    /// Remove a segment from the graph.
    ///
    /// Return true if graph modified.
    ///
    /// The descriptor held by the object is left invalid.
    bool remove_segment(Graph& graph, SegmentPtr seg);

    
    /// Return the end-point vertices of a segment.
    ///
    /// The pair will be nullptr if segment edge not in graph.
    ///
    std::pair<VertexPtr, VertexPtr> find_endpoints(Graph& graph, SegmentPtr seg);

    

};
#endif
