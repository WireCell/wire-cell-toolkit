#ifndef WIRECELLCLUS_GRAPH
#define WIRECELLCLUS_GRAPH

#include "WireCellUtil/Graph.h"
#include "WireCellUtil/PointCloudDataset.h"

#include <map>
#include <string>
#include <functional>

namespace WireCell::Clus::Graphs {

    namespace Weighted {

        /**
         * The basic graph type used in clus is an edge-weighted graph with
         * vertex descriptors that can serve as indices.  For indices to remain
         * stable, vertex removal is NOT supported. 
         */

        using dijkstra_distance_type = double;
        using edge_weight_type = double;
        using Graph = boost::adjacency_list<
            boost::vecS,        // vertices
            boost::vecS,        // edges
            boost::undirectedS, // edge direction (none)
            boost::property<boost::vertex_index_t, size_t>,
            boost::property<boost::edge_weight_t, edge_weight_type> 
            >;
        using graph_type = Graph;
        using vertex_type = boost::graph_traits<Graph>::vertex_descriptor;
        using edge_type = boost::graph_traits<Graph>::edge_descriptor;

        // A set of unique vertices or edges;
        using vertex_set = std::set<vertex_type>;
        using edge_set = std::set<edge_type>;
    
        // Filtered graphs and their predicates.
        using vertex_predicate = std::function<bool(vertex_type)>;
        using edge_predicate = std::function<bool(edge_type)>;
        using filtered_graph_type = boost::filtered_graph<graph_type, edge_predicate, vertex_predicate>;

        /// Results of a "discrete Voronoi tessellation" formed on the graph
        /// given with each Voronoi cell defined by one in a set of select
        /// "terminal" vertices.
        struct Voronoi {

            // A "map" from each graph vertex (index) to its nearest "terminal"
            // vertex.  terminal[v] == v for v in the set of terminal vertices.
            std::vector<vertex_type> terminal;

            // A "map" from each graph vertex (index) to the distance along the
            // path to its nearest "terminal".  distance[v] == 0.0 for v in the
            // set of terminal vertices.
            std::vector<edge_weight_type> distance;

            // A "map" from each graph vertex (index) to the edge into that
            // vertex from the vertex neighbor (the edge "source") that is
            // directly upstream in the walk from closest terminal to the
            // original vertex.  If the source is the nearest terminal for the
            // original vertex then the backwards walk is complete.  Otherwise,
            // the last_edge for the neighbor provides the next neighbor, etc.
            // Note, the value last_edge[v] for any v in the set of terminals is
            // not defined.  (You probably get 0's for both vertices).
            std::vector<edge_type> last_edge;

        };

        Voronoi voronoi(const graph_type& graph, const std::vector<vertex_type>& terminals);

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

            /// Return a graph view filtered on the set of vertices.  If accept
            /// is true (default) the graph only has vertices in the set.  If
            /// false, it has vertices in the original graph that are not in the
            /// set.  The vertex descriptors (indices) of the returned graph are
            /// the same as the original graph.  Use boost::copy_graph() to make
            /// a new graph with compactified vertices.
            filtered_graph_type reduce(const vertex_set& vertices, bool accept = true) const;

            /// As reduce(vertex_set) but filter on edges.
            filtered_graph_type reduce(const edge_set& edge, bool accept = true) const;
            
            /// Return a graph view filtered on edge weights.  If accept is true
            /// (default) edges with weights greater or equal to threshold will
            /// be kept and others removed and vice versa if accept is false.
            filtered_graph_type weight_threshold(edge_weight_type threshold, bool accept = true) const;

        };


    }

}

#endif

