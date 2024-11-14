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



#endif
