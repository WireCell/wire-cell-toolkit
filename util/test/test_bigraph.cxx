// https://stackoverflow.com/a/45850742

#include "WireCellUtil/Graph.h"

#include <boost/intrusive/set_hook.hpp>
#include <boost/intrusive/set.hpp>
#include <boost/any.hpp>
#include <iostream>

namespace bi = boost::intrusive;

/// vertex properties
struct VertexData : bi::set_base_hook<bi::link_mode<bi::auto_unlink>, bi::constant_time_size<false> > {
    std::string label;
    int num;

    VertexData(std::string label, int num)
      : label(label)
      , num(num)
    {
    }

    struct by_label {
        using type = std::string;
        const type& operator()(VertexData const& vd) const { return vd.label; }
    };
};

using by_label_idx_t = bi::set<VertexData, bi::constant_time_size<false>, bi::key_of_value<VertexData::by_label> >;

/// edges properties
struct EdgeData {
    std::string edge_name;
    double edge_confidence;
};

/// define the boost-graph
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VertexData,
                              boost::property<boost::edge_weight_t, double, EdgeData> >
    Graph;

int main()
{
    // using vertex_t = Graph::vertex_descriptor;

    Graph g;
    for (auto label : {
             "alerts", "amazed", "buster", "deaths", "ekes",   "Enoch",  "gale", "hug",   "input",  "knifed",
             "lire",   "man",    "pithy",  "Purims", "Rodger", "suckle", "Terr", "theme", "tiling", "vases",
         }) {
        boost::add_vertex(VertexData{label, 1 + rand() % 5}, g);
    }

    /// define vertexMap
    by_label_idx_t label_idx;
    auto reindex = [&] {
        label_idx.clear();
        for (auto vd : boost::make_iterator_range(boost::vertices(g))) label_idx.insert(g[vd]);
    };

    reindex();
    std::cout << "Index: " << label_idx.size() << " elements\n";

    g.clear();
    std::cout << "Index: " << label_idx.size() << " elements\n";

    for (auto& vertex : label_idx) {
        std::cout << vertex.label << " " << vertex.num << "\n";
    }
}
