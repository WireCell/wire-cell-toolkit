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

        // A quasi-edge type that simply records a pair of vertices.  Unlike
        // edge_type, the edge_pair_type is not dependent on the graph.  Use
        // make_vertex_pair() to assure an order.
        using vertex_pair = std::pair<vertex_type, vertex_type>;

        /// Return an ordered vertex pair with the smaller of {a,b} first.
        vertex_pair make_vertex_pair(vertex_type a, vertex_type b);


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

            /// A "map" from each graph vertex (index) to its nearest "terminal"
            /// vertex.  terminal[v] == v for v in the set of terminal vertices.
            std::vector<vertex_type> terminal;

            /// A "map" from each graph vertex (index) to the distance along the
            /// path to its nearest "terminal".  distance[v] == 0.0 for v in the
            /// set of terminal vertices.
            std::vector<edge_weight_type> distance;

            /// A "map" from each graph vertex (index) to the edge into that
            /// vertex from the vertex neighbor (the edge "source") that is
            /// directly upstream in the walk from closest terminal to the
            /// original vertex.  If the source is the nearest terminal for the
            /// original vertex then the backwards walk is complete.  Otherwise,
            /// the last_edge for the neighbor provides the next neighbor, etc.
            /// Note, the value last_edge[v] for any v in the set of terminals is
            /// not defined.  (You probably get 0's for both vertices).
            std::vector<edge_type> last_edge;

        };
            

        ///
        /// Free functions calculating Voronoi and related See GraphAlgorithm's
        /// methods of similar names for caching versions.
        ///

        /// Return the path of vertices FROM a given vertex TO its nearest
        /// terminal.  This simply walks last_edge as described above.  The
        /// result is undefined if the graph other than the one used to make
        /// this Voronoi struct.
        std::vector<vertex_type> terminal_path(const graph_type& graph, const Voronoi& vor, vertex_type ver);

        /// Return a Steiner graph (not tree).  This graph has all the
        /// vertices of the original graph but only the edges from the
        /// original graph which are on the shortest path between terminal
        /// vertices.  
        graph_type steiner_graph(const graph_type& graph, const Voronoi& vor);


        /// Structure to hold charge calculation parameters (from prototype)
        struct ChargeWeightingConfig {
            double Q0 = 10000.0;        // constant term from prototype
            double factor1 = 0.8;       // weighting factor 1 from prototype  
            double factor2 = 0.4;       // weighting factor 2 from prototype
            bool enable_weighting = true; // whether to apply charge weighting
        };

       


        /// Construct the "discrete graph Voronoi tessellation".
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

            /// Find all neighbors within nlevel hops from the input index.
            /// @param index The starting vertex index
            /// @param nlevel The number of levels (hops) to search
            /// @param include_self Whether to include the original vertex in the result (default: true)
            /// @return A set of vertex indices that are within nlevel hops from the input index
            vertex_set find_neighbors_nlevel(size_t index, int nlevel, bool include_self = true) const;

        };


    }

    /// Kill predicate of the "relaxed_strict" graph flavor (doc pr/53 sec
    /// 16.3): true if a sampled straight path fails the charge-presence test.
    /// S1: the 1cm sampling loop's last sample is the far endpoint itself and
    /// good by construction, so the 0.75 ratio is taken against the interior
    /// step count (num_steps - 1).  S2: the floor is nbad >= 2 (the legacy
    /// relaxed test uses > 2, unreachable for necks under ~3cm).  cap is 7
    /// for single-counter tests and 9 for the sum-of-two-planes tests.
    /// Pure; exposed for doctests.
    bool relaxed_strict_bad(int nbad, int num_steps, int cap = 7);

    /// Kill predicate of the "relaxed_strict_img" graph flavor (doc pr/53
    /// round 7 sec 18: S5, "3D-image support"). True if the sampled path has
    /// a run of consecutive "ghost" steps -- 1cm samples with no 3D image
    /// point within the configured radius on ANY of this cluster's own
    /// closely-components, and not excused by a dead channel on any plane --
    /// at least run_floor steps long, AND the edge itself is shorter than
    /// dis_cap_cm.
    ///
    /// Root cause this closes: relaxed_strict_bad only ever sees three
    /// independent 2D per-plane projections (Facade_Grouping.cxx
    /// has_closest_point, one kd2d radius query per plane); nothing in the
    /// existing path test intersects them in 3D, so three planes can each
    /// see charge from a DIFFERENT nearby track and every step reads "good"
    /// with nothing physically there at that 3D location. This predicate is
    /// an independent OR-kill alongside relaxed_strict_bad, not a
    /// replacement -- it keeps every existing kill and can only add new
    /// ones (relaxed_strict_img superset-of relaxed_strict by construction).
    ///
    /// The dis_cap_cm term exists because the closest-pair test's "nearest
    /// points between two components" can, for long edges, land on the
    /// low-density corona of one large blob rather than cross a genuine gap
    /// between two separate objects (doc pr/53 round 7 sec 18.2's offline
    /// threshold scan: event 52672 j=0 k=2, 45.9cm edge, max contiguous
    /// ghost run 19/45 steps yet the two "components" are corona and
    /// interior of the SAME dense shower, not two objects). Every
    /// owner-flagged target edge measured under 12cm; the cap default
    /// (15cm) keeps full margin on the targets while excluding that case.
    /// The run-length floor (not a raw ghost-step count or ratio) is
    /// likewise chosen from that scan: raw count/ratio did not separate the
    /// owner-flagged edges from the surviving-edge background at ALL
    /// (background edges reach comparable ghost counts/ratios), but the
    /// LONGEST CONTIGUOUS run does (owner targets: 5-10 steps; background
    /// median: 0, p90: 2). run_floor=4 keeps every owner target with a
    /// one-step margin against the tightest (5).
    ///
    /// The dead-channel exemption is what keeps this rule honest: the
    /// relaxed flavor exists specifically to bridge regions where 3D
    /// imaging is known to fail (dead wires); requiring 3D image support is
    /// in partial tension with that purpose, and a genuine dead-region
    /// crossing must not be penalized -- only a step with NO excuse at all
    /// counts as a ghost step.
    ///
    /// Pure; exposed for doctests.
    bool relaxed_img_bad(int max_ghost_run, double dis_cm, int run_floor = 4, double dis_cap_cm = 15.0);

}

#endif

