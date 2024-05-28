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

#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/connected_components.hpp>
#pragma GCC diagnostic pop

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/copy.hpp>



#endif
