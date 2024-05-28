/*  Boost throws 

[37/39] Compiling util/test/test_indexedgraph.cxx
In file included from ../../../../../opt/boost-1-76-0/include/boost/smart_ptr/detail/sp_thread_sleep.hpp:22,
                 from ../../../../../opt/boost-1-76-0/include/boost/smart_ptr/detail/yield_k.hpp:23,
                 from ../../../../../opt/boost-1-76-0/include/boost/smart_ptr/detail/spinlock_gcc_atomic.hpp:14,
                 from ../../../../../opt/boost-1-76-0/include/boost/smart_ptr/detail/spinlock.hpp:42,
                 from ../../../../../opt/boost-1-76-0/include/boost/smart_ptr/detail/spinlock_pool.hpp:25,
                 from ../../../../../opt/boost-1-76-0/include/boost/smart_ptr/shared_ptr.hpp:29,
                 from ../../../../../opt/boost-1-76-0/include/boost/property_map/vector_property_map.hpp:14,
                 from ../../../../../opt/boost-1-76-0/include/boost/property_map/property_map.hpp:602,
                 from ../../../../../opt/boost-1-76-0/include/boost/graph/graph_concepts.hpp:17,
                 from ../../../../../opt/boost-1-76-0/include/boost/graph/depth_first_search.hpp:18,
                 from ../../../../../opt/boost-1-76-0/include/boost/graph/connected_components.hpp:15,
                 from ../util/inc/WireCellUtil/IndexedGraph.h:16,
                 from ../util/test/test_indexedgraph.cxx:1:
../../../../../opt/boost-1-76-0/include/boost/config/pragma_message.hpp:24:34: note: ‘#pragma message: This header is deprecated. Use <iterator> instead.’
   24 | # define BOOST_PRAGMA_MESSAGE(x) _Pragma(BOOST_STRINGIZE(message(x)))
      |                                  ^~~~~~~
../../../../../opt/boost-1-76-0/include/boost/config/header_deprecated.hpp:23:37: note: in expansion of macro ‘BOOST_PRAGMA_MESSAGE’
   23 | # define BOOST_HEADER_DEPRECATED(a) BOOST_PRAGMA_MESSAGE("This header is deprecated. Use " a " instead.")
      |                                     ^~~~~~~~~~~~~~~~~~~~
../../../../../opt/boost-1-76-0/include/boost/detail/iterator.hpp:13:1: note: in expansion of macro ‘BOOST_HEADER_DEPRECATED’
   13 | BOOST_HEADER_DEPRECATED("<iterator>")
      | ^~~~~~~~~~~~~~~~~~~~~~~

 */


#define BOOST_DISABLE_PRAGMA_MESSAGE 1

#include "WireCellUtil/Graph.h"



// #include "WireCellUtil/IndexedGraph.h"

int main()
{
    return 0;
}
