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
            /// Enable per-step timing printouts (set via grapher_config.perf = true)
            bool perf{false};

            /// doc pr/29 D1.  Wire slack, in wires, used ONLY by the Steiner
            /// terminal filter (filter_by_reference_cluster).  The WCP
            /// prototype tests terminal containment with one wire of slack on
            /// each side of all three planes
            /// (PR3DCluster_steiner.h:285-290, :310-315, :336-341) while its
            /// get_extreme_wcps filter uses none (PR3DCluster_path.h:111-119).
            /// The toolkit implements both with one function, so the slack has
            /// to be a per-call-site parameter rather than an edit.
            /// 0 = the historical toolkit behaviour (no slack anywhere);
            /// 1 = prototype parity for the terminal filter alone.
            int terminal_wire_tol{0};

            /// doc pr/29 D12.  When true, the terminal filter's adjacent-slice
            /// fallback steps by the face's nticks-per-slice instead of by 1.
            /// The map it searches is keyed in TICKS, so the historical step of
            /// 1 never resolves and the whole fallback is dead code on any
            /// detector with more than one tick per slice (SBND: 4).
            /// false = that historical dead branch, bit-for-bit;
            /// true = the prototype's t, t+-1 slice search actually happens.
            bool terminal_adjacent_slice{false};

            /// doc pr/29 D2.  Forward create_steiner_tree's
            /// disable_dead_mix_cell down to the edge-weight charge
            /// calculation instead of letting create_enhanced_steiner_graph's
            /// `= true` default win.  The prototype computes the Qs/Qt that
            /// enter the edge weight with the SAME false it was called with
            /// (PR3DCluster_steiner.h:514,:521 inheriting
            /// Create_steiner_tree's parameter, passed false at
            /// create_steiner_graph:46); the toolkit drops the argument at the
            /// call and silently takes the other branch of calc_charge_wcp.
            /// The two branches differ in which planes enter the RMS: `true`
            /// sums all three then subtracts DEAD planes
            /// (charge_uncertainty > 1e10), `false` sums only planes with a
            /// NONZERO charge value.  Those predicates are independent, so
            /// this is not a no-op wherever they disagree.
            /// false = the historical toolkit behaviour (always `true`
            /// downstream, whatever the caller asked for);
            /// true = prototype parity, the caller's value is honoured.
            bool edge_charge_forward_dead_mix{false};

            /// doc pdvd/25 M3.  Per-point charge floor (electrons, per plane and
            /// on the plane-RMS) below which a point cannot be a Steiner
            /// terminal (find_peak_point_indices).  4000 is the WCP prototype
            /// constant, tuned to uBooNE's 3 mm pitch: on PDVD (7.65 / 5.10 mm
            /// pitch, stepped sampling) the W-plane per-point charge has a
            /// median of ~1400 e and only ~12 % of points clear 4000 on all
            /// three planes, so a fifth of the mains get NO terminals at all
            /// (a 3 m track included).  Default 4000 = bit-for-bit the
            /// historical behaviour; PDVD sets its own value.
            double terminal_charge_threshold{4000.0};

            /// doc pdvd/37 R1.  Minimum spacing, a real length in the WCT
            /// system of units, enforced ACROSS blob boundaries on the Steiner
            /// terminal set before the extreme points are added.
            ///
            /// Why a separate pass rather than a wider suppression radius:
            /// find_peak_point_indices runs once per blob
            /// (find_steiner_terminals, SteinerGrapher.cxx), and its
            /// map_index_charge holds only that blob's points, so an
            /// out-of-blob neighbour is skipped and the local-maximum test can
            /// never remove a blob's LAST candidate.  Terminal density is
            /// therefore floored at the candidate-bearing blob density (doc
            /// pdvd/31 round 6 measured 1.02 terminals per such blob), i.e. it
            /// is set by the time-slice pitch and the track's drift alignment
            /// rather than by anything physical.  Doc pdvd/37 sec.3.2 confirms
            /// that law on PDVD, SBND and uBooNE alike.
            ///
            /// 0 (the default) = no thinning, bit-for-bit the historical
            /// behaviour.  A positive value keeps terminals greedily in
            /// DECREASING charge order, admitting one only when no
            /// already-kept terminal lies strictly within it.
            double terminal_min_separation{0.0};
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

        /// A type that maps blob node indices (from sv.nodes()) to graph vertices.
        /// Using a node-index key (size_t, deterministic traversal order) instead of
        /// Blob* (heap-address-ordered, non-deterministic) ensures stable iteration
        /// order across runs.
        using blob_vertex_map = std::map<size_t, vertex_set>;


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
        ///  The real main entry method
        ///



        ///
        ///  Intermediate algorithm methods
        ///


        vertex_set find_peak_point_indices(const std::vector<const Facade::Blob*>& target_blobs, const std::string& graph_name,
                                   bool disable_dead_mix_cell = true, int nlevel = 1);
        /// Overload that accepts a precomputed point set for the target blobs, avoiding a rebuild of form_cell_points_map().
        ///
        /// doc pdvd/31 round 7: `n_candidate_points`, when non-null, receives the
        /// number of points that cleared terminal_charge_threshold BEFORE the
        /// local-maximum suppression ran.  That is the threshold's point-level
        /// selectivity, which the blob-level counter alone cannot separate from
        /// blob size -- PDVD blobs are geometrically larger than SBND's, so
        /// equal candidate-BEARING-blob fractions do not imply equal
        /// selectivity.  Measurement only; nullptr (the default) leaves every
        /// caller unchanged.
        vertex_set find_peak_point_indices(const vertex_set& blob_point_indices, const std::string& graph_name,
                                   bool disable_dead_mix_cell = true, int nlevel = 1,
                                   size_t* n_candidate_points = nullptr);

        blob_vertex_map form_cell_points_map();
        vertex_set find_steiner_terminals(const std::string& graph_name, bool disable_dead_mix_cell=true);
        /// Overload that accepts a precomputed blob->points map to avoid calling form_cell_points_map() twice.
        vertex_set find_steiner_terminals(const std::string& graph_name, bool disable_dead_mix_cell, const blob_vertex_map& cell_points_map);

        /// doc pdvd/37 R1.  Thin @p terminals so that no two survivors lie
        /// closer than Config::terminal_min_separation, walking them in
        /// decreasing charge order.  The charge is Cluster::calc_charge_wcp
        /// evaluated with the SAME terminal_charge_threshold and
        /// @p disable_dead_mix_cell that selected them in phase 1, so the
        /// survivor of a crowd is the member the peak finder itself ranks
        /// highest.  Returns @p terminals unchanged when the knob is <= 0.
        vertex_set thin_terminals_by_separation(const vertex_set& terminals,
                                                bool disable_dead_mix_cell);

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

        // Enable per-step timing printouts inside hot functions
        bool m_perf{false};


        // XIN: add any more data and methods you need here.  
         /// Track edges added by each graph modification operation.
        /// Stored as a vector (not a set) because edge_type ordering is pointer-based
        /// (non-deterministic) and removal is order-insensitive.
        std::map<std::string, std::vector<edge_type>> m_added_edges_by_graph;

        /// Helper to invalidate GraphAlgorithms cache for a specific graph
        void invalidate_graph_algorithms_cache(const std::string& graph_name);

        /// Helper to store added edges for later removal
        void store_added_edges(const std::string& graph_name, const std::vector<edge_type>& edges);

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
        /// @param wire_tol wire slack (Config::terminal_wire_tol); 0 = historical
        /// @param slice_stride adjacent-slice key step in TICKS; 1 = historical
        ///        (and, on any detector with >1 tick per slice, never matches)
        bool is_point_spatially_related_to_reference(
            size_t point_idx,
            const Facade::Cluster::time_blob_map_t& ref_time_blob_map,
            int wire_tol = 0,
            int slice_stride = 1
        ) const;

        /// Ticks per slice for (apa, face), the correct step for an
        /// adjacent-slice lookup in a time_blob_map.  Returns 1 -- the
        /// historical, never-matching step -- when the grouping does not know
        /// this face, rather than throwing out of a hot filter loop.
        int nticks_per_slice_or_1(int apa, int face) const;



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
        
        // m_steiner_graph_terminal_indices removed — not populated; reserved for
        // recover_steiner_graph() if that function is ported in the future.
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
