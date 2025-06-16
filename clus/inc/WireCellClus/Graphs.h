#ifndef WIRECELLCLUS_GRAPH
#define WIRECELLCLUS_GRAPH

#include "WireCellUtil/Graph.h"
#include "WireCellUtil/PointCloudDataset.h"

#include <map>
#include <string>

namespace WireCell::Clus::Graphs {

    namespace Weighted {

        /**
         * The basic graph type used in clus is an edge-weighted graph with
         * vertex descriptors that can serve as indices.  For indices to remain
         * stable, vertex removal is NOT supported. 
         */

        using dijkstra_distance_type = double;
        using Graph = boost::adjacency_list<
            boost::vecS,        // vertices
            boost::vecS,        // edges
            boost::undirectedS, // edge direction (none)
            boost::property<boost::vertex_index_t, size_t>,
            boost::property<boost::edge_weight_t, float> 
            >;
    
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
            const Graph& m_graph;

            // LRU cache configuration
            static constexpr size_t DEFAULT_MAX_CACHE_SIZE = 50;  // Adjust as needed
            size_t m_max_cache_size;
            
            // LRU cache implementation for shortest paths
            // List maintains access order (most recently used at front)
            mutable std::list<size_t> m_access_order;

            // Map from source to {list_iterator, ShortestPaths}
            mutable std::unordered_map<size_t, 
                std::pair<std::list<size_t>::iterator, ShortestPaths>> m_sps;

            // Lazy calculate dijkstra shortest path results.
            // mutable std::unordered_map<size_t, ShortestPaths> m_sps;

            mutable std::vector<size_t> m_cc;

            // Helper method to update LRU cache
            void update_cache_access(size_t source) const;
            void evict_oldest_if_needed() const;

        public:
            GraphAlgorithms(const Graph& graph, size_t max_cache_size = DEFAULT_MAX_CACHE_SIZE);

        
            /// Return the intermediate result that gives access to the shortest
            /// paths from the source vertex all possible destination vertices.
            const ShortestPaths& shortest_paths(size_t source) const;
        
            /// Return the unique vertices vertices along the shortest path from
            /// source vertex to destination vertex, inclusive.
            const std::vector<size_t>& shortest_path(size_t source, size_t destination) const;

            /// Return a "CC" array giving connected component subgraphs.
            const std::vector<size_t>& connected_components() const;

            /// Get current cache size and maximum cache size
            size_t cache_size() const { return m_sps.size(); }
            size_t max_cache_size() const { return m_max_cache_size; }
            
            /// Clear the shortest paths cache
            void clear_cache() const;

        };


    }

}

#endif

