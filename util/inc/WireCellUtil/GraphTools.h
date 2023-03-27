#ifndef WIRECELLUTIL_GRAPHTOOLS
#define WIRECELLUTIL_GRAPHTOOLS

#define BOOST_DISABLE_PRAGMA_MESSAGE 1
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/range.hpp>
#include <boost/graph/graphviz.hpp>

namespace WireCell::GraphTools {

     
    /**
       The mir() functions allow for more convenient iteration over
       boost graph edges:
    
           for (const auto& edge : mir(boost::edges(g))) {
               auto tail = boost::source(edge, g);
               auto head = boost::target(edge, g);
               // ...
           }

       Or vertices

           for (const auto& vtx : mir(boost::vertices(g))) { ... };

       Or even briefer:

           for (const auto& vtx : vertex_range(g)) { ... };

     */

    // fixme: these are more generic than just to graphs...
    template <typename It> boost::iterator_range<It> mir(std::pair<It, It> const& p) {
        return boost::make_iterator_range(p.first, p.second);
    }
    template <typename It> boost::iterator_range<It> mir(It b, It e) {
        return boost::make_iterator_range(b, e);
    }

    template <typename Gr> 
    boost::iterator_range<typename boost::graph_traits<Gr>::vertex_iterator> vertex_range(Gr& g) {
        return mir(boost::vertices(g));
    }

    // Return graph as string holding GraphViz dot representation.
    template<typename Gr>
    std::string dotify(const Gr& gr)
    {
        using vertex_t = typename boost::graph_traits<Gr>::vertex_descriptor;
        std::stringstream ss;
        boost::write_graphviz(ss, gr, [&](std::ostream& out, vertex_t v) {
            const auto& dat = gr[v];
            out << "[label=\"" << dat.code << dat.id << "\"]";
        });
        return ss.str() + "\n";
    }

}

#endif
