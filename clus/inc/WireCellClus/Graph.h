#ifndef WIRECELLCLUS_GRAPH
#define WIRECELLCLUS_GRAPH

#include "WireCellUtil/Graph.h"
#include "WireCellUtil/PointCloudDataset.h"

#include <map>
#include <string>

// Graph types used in clus.  Note/fixme, this is actually fairly generic code.
// Adding boost.serialization to WCT's dependencies would allow for an even more
// generic I/O with Datasets to be usefully moved into util/
namespace WireCell::Clus::Graph {

    // An atomic, contiguous point cloud.
    using dataset_t = PointCloud::Dataset;

    /** 
     * A graph defined here with a load()/save() function may be persisted
     * through a "store" that provides a set of named PointCloud::Dataset.
     * Either one or two datasets are needed and their names in the store follow
     * the convention:
     *
     * - "graph_{name}_vertices" :: holds arrays of vertex properties.
     *
     * - "graph_{name}_edges" :: holds vertices of edges in "tail" and "head" 1d
     *   arrays of type int64_t and any edge properties named by the their
     *   attribute name.
     * 
     * The edges dataset is required.  The vertices dataset may be omitted if
     * the graph type does not have vertex properties.
     *
     * If the graph has vertex properties then the major size of the
     * "graph_{name}_vertices" dataset determines the number of vertices.
     * Otherwise, the "graph_{name}_edges" dataset should provide a metadata
     * attribute "num_vertices" to provide the number of vertices.  If load() is
     * faced with neither it will chose the number of vertices to be just large
     * enough to hold the largest edge vertex descriptor.
     *
     * Anything else, including missing vertex or edge property arrays is an
     * error.  Additional vertex or property arrays are ignored.  load() and
     * save() accept maps to allow redirection between property names.
     *
     * Note, the graph dataset "store" is identical to the "local PCs" map held
     * by a PC tree node.  It is intended to use functions in this API to
     * persistence of graphs in a PC tree.  This API has
     */
    using store_t = std::map<std::string, dataset_t>;
    using name_map_t = std::map<std::string, std::string>;

    /// Ensconce the naming conventions for graph datasets.
    inline std::string graph_name(const std::string& name) { return "graph_" + name; }
    inline std::string vertices_name(const std::string& name) { return graph_name(name) + "_vertices"; }
    inline std::string edges_name(const std::string& name) { return graph_name(name) + "_edges"; }


    /**
     * An "ident" graph type has a vertex property of integer type.  This
     * integer is intended to represent a user-defined "identity".  Typically,
     * the integer is used as an index in an externally provided array of
     * values.  Note, the "ident" value is distinct from the vertex descriptor.
     * A single precision floating point "weight" edge property is also included.
     */
    namespace Ident {

        using ident_type = int;
        using weight_type = float;

        // Note, this is a struct property and in particular NOT the boost graph
        // "index" property as the latter would then become the vertex
        // descriptor.
        struct VertexProp {
            ident_type ident;
        };
        using EdgeProp = boost::property<boost::edge_weight_t, weight_type>; 

        // Vertices are forced to be unique.  Edges may be duplicate.
        using graph_type = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                                                  VertexProp, EdgeProp>;
        using vertex_descriptor = boost::graph_traits<graph_type>::vertex_descriptor;
        using edge_descriptor = boost::graph_traits<graph_type>::edge_descriptor;
        
        using graph_ptr = std::unique_ptr<graph_type>;

        /** Clear and load a graph of the given name from the store. 
         *
         * If no graph data is found in the store, the graph is left cleared.
         *
         * The dataset naming convention is used to locate the datasets in the
         * store.  By default, the array named after the properties ("ident",
         * "weight") and edge ("tail", "head") are used.  The optional props may
         * be given to map these names to other names of the arrays to use.
         *
         * Return true only on success.
         */
        bool load(const std::string& name, const store_t& store, graph_type& g,
                  const name_map_t& props={});

        /** Save a graph to the store under the given name.
         *
         * Any existing graph data of that name is overwritten.
         *
         * See load() for description of the arguments.
         *
         * Return true only on success.
         */
        bool save(const std::string& name, const graph_type& g, store_t& store,
                  const name_map_t& props={});

        /// Return the ident values in vertex order.
        std::vector<ident_type> idents(const graph_type& g);

        /// Return the edge weights in edge order.
        std::vector<weight_type> weights(const graph_type& g);

        /// Return the edges as matched pairs of vertex descriptors.  The first
        /// vector in the pair holds the "tail" vertices, the second holds the
        /// "head" vertices.
        std::pair<std::vector<vertex_descriptor>, std::vector<vertex_descriptor>>
        edges(const graph_type& g);

        // Result of Dijkstra's shortest paths from vertex.
        struct DijkstraShortestPaths {
            vertex_descriptor vertex;
            std::vector<vertex_descriptor> parents;
            std::vector<int> distances;
        };
        // Find vertex with index, then shortest paths to other vertices.
        DijkstraShortestPaths dijkstra_shortest_paths(const graph_type& graph, size_t point_index);

    }

    // A weighted-edge graph
    namespace Weighted {
        using EdgeProp = boost::property<boost::edge_weight_t, double>;
        using graph_type = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                                                 boost::no_property, EdgeProp>;
        using vertex_descriptor = boost::graph_traits<graph_type>::vertex_descriptor;
        using edge_descriptor = boost::graph_traits<graph_type>::edge_descriptor;
    }        


}

#endif

