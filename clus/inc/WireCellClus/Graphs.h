#ifndef WIRECELLCLUS_GRAPH
#define WIRECELLCLUS_GRAPH

#include "WireCellUtil/Graph.h"
#include "WireCellUtil/PointCloudDataset.h"

#include <map>
#include <string>

namespace WireCell::Clus::Graphs {

    namespace Weighted {

        /**
         * The basic graph type used in clus is an edge-weighted graph with vertex
         * descriptors that can serve as indices.  For indices to work, vertex
         * removal is NOT supported.
         */

#define CLUS_FIXME              // eventually remove this, keeping new

#ifdef CLUS_FIXME
        // Define to switch back to the types used in WCP, which are "bad".
        // See Graphs.cxx for use of this type.
        using dijkstra_distance_type = int;
#else
        using dijkstra_distance_type = double;
#endif

        using Graph = boost::adjacency_list<
#ifdef CLUS_FIXME
            boost::setS,        // vertices
#else
            boost::vecS,        // vertices
#endif
            boost::vecS,        // edges
            boost::undirectedS, // edge direction (none)
            boost::property<boost::vertex_index_t, size_t>,
            boost::property<boost::edge_weight_t, float> 
            >;
    
        using GraphPtr = std::unique_ptr<Graph>;

        /// Embody all possible shortest paths from a given source index.
        class ShortestPaths {
        
            // Dijkstra result from source
            size_t m_source;
            std::vector<size_t> m_predecessors; 
            // distances are not currently used.
        
            // Lazy calculate path for given destination indices.
            mutable std::unordered_map<size_t, std::vector<size_t>> m_paths;
        
        public:
            ShortestPaths(size_t source, const std::vector<size_t> predecessors);
        
            /// Return the unique vertices along shortest path from our source
            /// to destination.  inclusive.
            const std::vector<size_t>& path(size_t destination) const;
        };

        // Bind some graph algorithms to a graph, with caching..
        class GraphAlgorithms {
            GraphPtr m_graph;

            // Lazy calculate dijkstra shortest path results.
            mutable std::unordered_map<size_t, ShortestPaths> m_sps;
            mutable std::vector<int> m_cc;

        public:
            GraphAlgorithms(GraphPtr&& graph);
        
            /// Return the intermediate result that gives access to the shortest
            /// paths from the source vertex all possible destination vertices.
            const ShortestPaths& paths(size_t source) const;
        
            /// Return the unique vertices vertices along the shortest path from
            /// source vertex to destination vertex, inclusive.
            const std::vector<size_t>& path(size_t source, size_t destination) const;

            /// Return a "CC" array giving connected component subgraphs.
            const std::vector<int>& connected_components() const;

        };

    }

}

#endif

