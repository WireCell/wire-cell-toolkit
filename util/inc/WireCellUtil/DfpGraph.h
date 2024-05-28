#ifndef WIRECELLUTIL_DFPGRAPH
#define WIRECELLUTIL_DFPGRAPH

#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Graph.h"

#include <tuple>
#include <vector>

namespace WireCell {

    struct DfpGraph {
        struct VertexProperty {
            std::string type;
            std::string name;

            VertexProperty(const std::string& t, const std::string& n)
              : type(t)
              , name(n)
            {
            }
            VertexProperty(const std::string& tn) { std::tie(type, name) = WireCell::String::parse_pair(tn); }
            VertexProperty()
              : type("")
              , name("")
            {
            }

            bool operator==(const VertexProperty& rhs) const { return type == rhs.type && name == rhs.name; }
            bool operator<(const VertexProperty& rhs) const
            {
                if (type == rhs.type) {
                    return name < rhs.name;
                }
                return type < rhs.type;
            }
        };

        struct EdgeProperty {
            int tail;
            int head;
            EdgeProperty(int t, int h)
              : tail(t)
              , head(h)
            {
            }
            EdgeProperty()
              : tail(-1)
              , head(-1)
            {
            }
        };
        typedef boost::adjacency_list<boost::vecS, boost::setS, boost::directedS, VertexProperty, EdgeProperty> Graph;

        typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef boost::graph_traits<Graph>::edge_descriptor Edge;

        typedef std::map<VertexProperty, Vertex> VertexPropertyMap;

        typedef std::tuple<VertexProperty, VertexProperty, EdgeProperty> Connection;

        Graph graph;
        VertexPropertyMap vertex_property_map;

        Vertex get_add_vertex(const VertexProperty& tn);

        /// Configure with list of connection dictionaries like:
        ///
        /// {
        ///     "tail": {"type":"TailType","name":"tailname", "port":0},
        ///     "head": {"type":"HeadType","name":"headname", "port":42}
        /// },
        ///
        /// Only "type" key is required.
        void configure(const Configuration& cfg);

        /// Explicitly connect one edge
        Edge connect(const std::string& tail_type, const std::string& tail_name, int tail_port,
                     const std::string& head_type, const std::string& head_name, int head_port);

        std::vector<Connection> connections();
    };
}  // namespace WireCell

inline std::ostream& operator<<(std::ostream& os, const WireCell::DfpGraph::VertexProperty& vp)
{
    os << vp.type << ":" << vp.name;
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const WireCell::DfpGraph::EdgeProperty& ep)
{
    os << "( " << ep.tail << " --> " << ep.head << " )";
    return os;
}

#endif
