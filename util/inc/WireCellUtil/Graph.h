// Include this header and not directly any bare boost/* in order to avoid some
// compiler warnings.

#ifndef WIRECELL_GRAPH
#define WIRECELL_GRAPH

// fixme: watchme: Boost started to deprecate some internal header
// inclusion which is not, as best as I can tell, any of our problem.
// The message is:
//
// ../../../../../opt/boost-1-76-0/include/boost/config/pragma_message.hpp:24:34: note: ‘#pragma message: This header is deprecated. Use <iterator> instead.’
//
//  This arises from a deeply nested #include well beyond anything
//  which is obvious here.
//
//  If/when this is cleaned up in Boost, remove this comment and the
//  next line.
// #define BOOST_ALLOW_DEPRECATED_HEADERS 1
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/properties.hpp>

#ifdef __clang__
#  if defined(__has_warning)
#    define HAS_WARNING(warning) __has_warning(warning)
#  else
#    define HAS_WARNING(warning) 1
#  endif
#else
#  define HAS_WARNING(warning) 1
#endif

#if HAS_WARNING("-Wmaybe-uninitialized")
#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graphviz.hpp> // keep this at the bottom, ref clus/test/doctest_boost_dijkstra.cxx

#if HAS_WARNING("-Wmaybe-uninitialized")
#pragma GCC diagnostic pop
#endif

#include <boost/graph/copy.hpp> // keep this at the bottom, ref clus/test/doctest_boost_dijkstra.cxx

#include <utility>

namespace WireCell {

    // Move the contents of one boost::adjacency_list into another.
    //
    // boost::adjacency_list (through at least 1.85) user-declares its
    // copy constructor/assignment and provides no move members, so
    // move-construction and move-assignment of graphs silently degrade
    // to a whole-graph copy_impl().  This transfers the underlying
    // vertex and edge containers instead.  Container moves do not
    // relocate elements, so the stored edge iterators and edge-property
    // pointers held inside remain valid and vertex/edge iteration order
    // is exactly the source's.  Any graph-level property of src is NOT
    // transferred (WCT cluster-style graphs use no_property).  src is
    // left empty.
    template <typename AdjacencyList>
    inline void move_graph(AdjacencyList& dst, AdjacencyList& src)
    {
        dst.m_vertices = std::move(src.m_vertices);
        dst.m_edges = std::move(src.m_edges);
        src.clear();
    }
}

#endif
