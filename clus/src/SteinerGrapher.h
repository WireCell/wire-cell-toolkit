/** This provides a "workspace" class for creating the steiner graph.

    It is meant to be equivalent to the Steiner-related SUBSET of WCP's
    PR3DCluster methods and data.

    See also SteinerFunctions.h for any free functions.
*/

#ifndef WIRECELLCLUS_STEINER
#define WIRECELLCLUS_STEINER

#include "WireCellClus/Graphs.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/IPCTransform.h"

#include "WireCellIface/IBlobSampler.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IPCTreeMutate.h"

#include "WireCellUtil/Logging.h"

#include <set>
#include <map>

namespace WireCell::Clus::Steiner {


    class Grapher {
    public:

        // This holds various "global" and const info sources.  See
        // CreateSteinerGraph for an example of how it is provided.
        struct Config {
            IDetectorVolumes::pointer dv;
            WireCell::Clus::IPCTransformSet::pointer pcts;
            IPCTreeMutate::pointer retile;
            /// do we even need samplers?
            // std::map<int, std::map<int, WireCell::IBlobSampler::pointer>> samplers;
        };
        Log::logptr_t log;

        /// Construct with an existing cluster facade.  Caller must assure the
        /// underlying cluster node is kept live.
        Grapher(Facade::Cluster& cluster, const Config& cfg, Log::logptr_t log);
        Grapher() = delete;

        /// Construct a Grapher with some cluster and take the rest of what we
        /// need from the other grapher.
        Grapher(Facade::Cluster& cluster, const Grapher& other);

        ///
        ///  Types
        ///
        
        /// Forward some types from Graphs.h
        using graph_type = WireCell::Clus::Graphs::Weighted::graph_type;
        using vertex_type = WireCell::Clus::Graphs::Weighted::vertex_type;
        using edge_type = WireCell::Clus::Graphs::Weighted::edge_type;
        using vertex_set = WireCell::Clus::Graphs::Weighted::vertex_set;
        using edge_set = WireCell::Clus::Graphs::Weighted::edge_set;
        using edge_weight_type = WireCell::Clus::Graphs::Weighted::edge_weight_type;

        /// A type that maps blobs to graph vertices
        using blob_vertex_map = std::map<const Facade::Blob*, vertex_set>;


        ///
        /// Basic data accessors.
        ///
        
        Facade::Cluster& cluster() { return m_cluster; }
        const Facade::Cluster& cluster() const { return m_cluster; }
        Config config() const { return m_config; }


        ///
        /// Helper methods - these are general purpose, primitive.
        /// 

        ///
        /// Some special graph access.  See also Facade::Mixins::Graphs in
        /// Facade_Mixins.h for more graph acessors.
        /// 

        ///  Get a graph, possibly making it on the fly if flavor is one of the
        ///  3 reserved names.
        graph_type& get_graph(const std::string& flavor = "basic");
        const graph_type& get_graph(const std::string& flavor = "basic") const ;
        
        /// Remove the flavor of graph from the other Grapher and move it to
        /// this one.  Give a non-empty value for "our_flavor" to store the
        /// transferred graph under a different name.
        void transfer_graph(Grapher& other,
                            const std::string& flavor = "basic",
                            std::string our_flavor = "");


        ///
        /// Some special PC access.
        ///
        
        /// Return a PC with the given name from our cluster node's local PCs.
        /// If it does not exist, one is derived from the default scoped view,
        /// saved to the given name, and a reference is returned.
        PointCloud::Dataset& get_point_cloud(const std::string& name = "default");

        /// Store a point cloud by std::move() in our cluster's local PCs.
        void put_point_cloud(PointCloud::Dataset&& pc, const std::string& name = "default");

        /// Store a point cloud by copy in our cluster's local PCs.
        void put_point_cloud(const PointCloud::Dataset& pc, const std::string& name = "default");

        /// Remove the named point cloud from the other Grapher and move it to
        /// this one.  Give a non-empty value for "our_name" to store the
        /// transferred PC under a different name.
        void transfer_pc(Grapher& other,
                         const std::string& name = "default",
                         const std::string& our_name = "");


        ///
        /// Temporary algorithm to show how to do some things.
        ///
        graph_type fake_steiner_graph();        


        ///
        ///  The real main entry method
        ///

        /// Create and return a steiner graph for the cluster.
        graph_type create_steiner_graph();


        ///
        ///  Intermediate algorithm methods
        ///


        vertex_set find_peak_point_indices(const std::vector<const Facade::Blob*>& target_blobs, const std::string& graph_name,
                                   bool disable_dead_mix_cell = true, int nlevel = 1);

        blob_vertex_map form_cell_points_map();
        vertex_set find_steiner_terminals(const std::string& graph_name, bool disable_dead_mix_cell=true);

        /// Establish edges between points in the same blob (mcell) with weighted connectivity
        /// This modifies the given graph and tracks added edges for later removal
        /// Uses find_steiner_terminals() to determine edge weights:
        /// - Both terminals: distance * 0.8
        /// - One terminal: distance * 0.9
        /// - Neither terminal: no edge added
        void establish_same_blob_steiner_edges(const std::string& graph_name, 
                                               bool disable_dead_mix_cell=true);
    
         /// Remove previously added same-blob Steiner edges
        void remove_same_blob_steiner_edges(const std::string& graph_name);

        /// Create Steiner tree with optional reference cluster and path constraints
        /// This is the main entry point equivalent to prototype's Create_steiner_tree
        void create_steiner_tree(
            const Facade::Cluster* reference_cluster = nullptr,
            const std::vector<size_t>& path_point_indices = {},
            const std::string& graph_name = "basic_pid",
            const std::string& steiner_graph_name = "steiner_graph",
            bool disable_dead_mix_cell = true,
            const std::string& steiner_pc_name = "steiner_pc"
        );


        /// Get the flag indicating which vertices in the steiner graph are terminals
        const std::vector<bool>& get_flag_steiner_terminal() const { 
            return m_flag_steiner_terminal; 
        }
        
        /// Get mapping from original to new vertex indices
        const std::map<vertex_type, vertex_type>& get_old_to_new_mapping() const { 
            return m_old_to_new_index; 
        }
        
        /// Get mapping from new to original vertex indices  
        const std::map<vertex_type, vertex_type>& get_new_to_old_mapping() const { 
            return m_new_to_old_index; 
        }
        
        /// Check if a vertex in the steiner graph is a terminal
        bool is_steiner_terminal(vertex_type steiner_vertex) const {
            if (steiner_vertex >= m_flag_steiner_terminal.size()) {
                return false;
            }
            return m_flag_steiner_terminal[steiner_vertex];
        }
        
        /// Get original vertex index from steiner graph vertex index
        vertex_type get_original_vertex(vertex_type steiner_vertex) const {
            auto it = m_new_to_old_index.find(steiner_vertex);
            if (it != m_new_to_old_index.end()) {
                return it->second;
            }
            return SIZE_MAX; // Invalid index
        }
        
        /// Get steiner graph vertex index from original vertex index  
        vertex_type get_steiner_vertex(vertex_type original_vertex) const {
            auto it = m_old_to_new_index.find(original_vertex);
            if (it != m_old_to_new_index.end()) {
                return it->second;
            }
            return SIZE_MAX; // Invalid index
        }


    private:
        // The Grapher "wraps" a Cluster.  As the Cluster is a *facade* of an
        // underlying PC tree node, we do not own the Cluster and we rely on
        // whoever owns us to keep the underlying cluster node alive. as long as
        // we are alive.
        Facade::Cluster& m_cluster;     

        // This holds various "global" info sources
        const Config& m_config;


        // XIN: add any more data and methods you need here.  
         /// Track edges added by each graph modification operation
        /// Maps graph name to set of edges added to that graph
        std::map<std::string, edge_set> m_added_edges_by_graph;

        /// Helper to invalidate GraphAlgorithms cache for a specific graph
        void invalidate_graph_algorithms_cache(const std::string& graph_name);

        /// Helper to store added edges for later removal
        void store_added_edges(const std::string& graph_name, const edge_set& edges);

        /// Helper to check if two vertices (points) belong to the same blob
        bool same_blob(vertex_type v1, vertex_type v2) const;

        /// Helper to calculate distance between two vertices
        double calculate_distance(vertex_type v1, vertex_type v2) const;

        /// Helper to get blob for a given vertex (point index)
        const Facade::Blob* get_blob_for_vertex(vertex_type vertex) const;


        // additional helper functions 

        /// Filter steiner terminals based on spatial relationship with reference cluster
        vertex_set filter_by_reference_cluster(
            const vertex_set& terminals,
            const Facade::Cluster* reference_cluster
        ) const;

        /// Filter steiner terminals based on path constraints
        vertex_set filter_by_path_constraints(
            const vertex_set& terminals,
            const std::vector<size_t>& path_point_indices
        ) const;

        /// Get extreme points considering reference cluster constraints
        vertex_set get_extreme_points_for_reference(
            const Facade::Cluster* reference_cluster
        ) const;


        /// Check if a point is spatially related to reference cluster's time-blob mapping
        bool is_point_spatially_related_to_reference(
            size_t point_idx,
            const Facade::Cluster::time_blob_map_t& ref_time_blob_map
        ) const;



        // temporary ...
        size_t find_closest_vertex_to_point(const Point& point) const;

        /// Create steiner subset point cloud with proper wire indices
        /// (matches prototype point_cloud_steiner creation)
        PointCloud::Dataset create_steiner_subset_pc_with_indices(
            const vertex_set& steiner_indices) const;

        /// Flag indicating which vertices in the reduced steiner graph are actual terminals
        /// vs. intermediate Steiner points (matches prototype flag_steiner_terminal)
        std::vector<bool> m_flag_steiner_terminal;
        
        /// Mapping from original graph vertex indices to reduced steiner graph indices
        /// (matches prototype map_old_new_indices)
        std::map<vertex_type, vertex_type> m_old_to_new_index;
        
        /// Mapping from reduced steiner graph indices to original graph indices  
        /// (matches prototype map_new_old_indices)
        std::map<vertex_type, vertex_type> m_new_to_old_index;
        
        /// Set of vertices that are steiner graph terminals (for edge creation logic)
        vertex_set m_steiner_graph_terminal_indices;
    };


}

namespace WireCell::Clus::Graphs::Weighted {
        /// Calculate charge-weighted distance between two vertices
        /// (matches prototype edge weighting logic)
        double calculate_charge_weighted_distance(
            double geometric_distance,
            double charge_source,
            double charge_target, 
            const ChargeWeightingConfig& config = ChargeWeightingConfig{});

        /// Calculate vertex charges using cluster facade method
        std::map<vertex_type, double> calculate_vertex_charges(
            const vertex_set& vertices, 
            const PointCloud::Dataset& pc,
            const WireCell::Clus::Facade::Cluster& cluster,
            double charge_cut,
            bool disable_dead_mix_cell);

        /// Enhanced steiner graph creation with full prototype functionality
        struct EnhancedSteinerResult {
            graph_type graph;                                    // reduced vertex graph
            PointCloud::Dataset point_cloud;                     // subset point cloud  
            std::vector<bool> flag_steiner_terminal;             // terminal flags
            std::map<vertex_type, vertex_type> old_to_new_index; // index mappings
            std::map<vertex_type, vertex_type> new_to_old_index; // reverse mappings
            vertex_set steiner_terminal_indices;                 // original terminal set
            std::map<vertex_type, double> vertex_charges;        // calculated charges
        };

        /// Create steiner graph with full prototype matching functionality
        EnhancedSteinerResult create_enhanced_steiner_graph(
            const graph_type& base_graph,
            const vertex_set& terminal_vertices,
            const PointCloud::Dataset& original_pc,
            const WireCell::Clus::Facade::Cluster& cluster,
            const ChargeWeightingConfig& charge_config = ChargeWeightingConfig{},
            bool disable_dead_mix_cell = true
        );

        void establish_same_blob_steiner_edges_steiner_graph(EnhancedSteinerResult& result, 
            const WireCell::Clus::Facade::Cluster& cluster);
}


#endif
