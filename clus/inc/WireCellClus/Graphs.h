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

    /// Kill predicate of the "relaxed_strict_img_2d" graph flavor (doc pr/56
    /// round 2: S6, "2D wind/tick connectivity"). True if at least one
    /// non-excused plane shows its two candidate components' real
    /// fired-pixel footprints (a bounded 2D flood-fill in that plane's
    /// (wire_index, time_slice) grid, live-or-dead-channel passable) NOT
    /// connected to each other. W is never excused; U/V are excused
    /// (gap_u/gap_v ignored) when the candidate path already tests as
    /// quasi-parallel to that plane's wire direction (the existing
    /// angle1/angle2 < wire_angle_tol test in connect_graph_relaxed_strict,
    /// reused here rather than re-detected) -- the physical justification
    /// being real hit-finding inefficiency from long same-wire induction
    /// pulse trains on a prolonged/isochronous topology.
    ///
    /// Owner's rule: kill on >= 1 non-excused gap (stricter than S1-S5's
    /// "at least 2 views" convention, by explicit instruction) -- this
    /// predicate is an independent OR-kill alongside relaxed_strict_bad and
    /// relaxed_img_bad, not a replacement; it can only add new kills.
    ///
    /// Pure; exposed for doctests. The 2D flood-fill itself lives in
    /// connect_graph_relaxed_strict.cxx (not pure -- it queries live
    /// Facade::Grouping state) and is verified by end-to-end event replay,
    /// not a synthetic doctest fixture; see doc pr/56 round 2 sec verification.
    bool two_d_connectivity_bad(bool gap_u, bool gap_v, bool gap_w,
                                 bool excuse_u, bool excuse_v);

    /// doc pr/57 round 6: inputs to the S6 rescue decision, measured at the
    /// kill site for one S6-KILLED candidate edge.  Everything is in
    /// "analysis units" (cm, degrees, counts, fractions) so the C++
    /// predicate is a line-by-line port of the Python rule it was fitted
    /// with (sbnd_xin/scripts/analysis/pr57/oc56_fit.py rescue()), and a
    /// doctest can replay real fitted cases.  Sentinels: ab_local_deg < 0
    /// means "no local axis measurable"; cov[p] < 0 means "no coverage
    /// measurement" (plane had no seed data); slope/ext are clamped at >= 0
    /// with 0 meaning "no measurement" (matching the Python `or 0.0`);
    /// close_mx[p] == 5 means "does not close even at a (4,4) stencil";
    /// has_plane[p] false means the plane had no seeds from one side and
    /// never voted.
    struct S6RescueInput {
        double dis_cm = 0;        ///< candidate edge 3D length
        bool gap[3] = {false, false, false};      ///< per-plane (1,1) BFS gap
        bool excuse_u = false, excuse_v = false;  ///< wire-parallel excuses
        bool has_plane[3] = {false, false, false};
        double ab_local_deg = -1; ///< angle between local PCA axes at the two break points
        int npmin = 0;            ///< smaller component's point count
        double lmin_cm = 0;       ///< smaller component's principal-axis extent
        int close_mx[3] = {5, 5, 5};   ///< smallest d with (d,d) stencil connected
        double cov[3] = {-1, -1, -1};  ///< wire coverage across the seed span
        double slope[3] = {0, 0, 0};   ///< |dslice/dwire| of the combined seeds
        double ext_med[3] = {0, 0, 0}; ///< median per-wire fired-slice extent (ticks)
        double ov[3] = {-1, -1, -1};   ///< seed wire-range overlap fraction
        int dead_w = 0;           ///< dead W wires across the W seed span
    };

    /// doc pr/57 round 6: rescue predicate of the "relaxed_strict_img_2d_rescue"
    /// graph flavor. True => an S6-killed candidate is kept after all. Fitted
    /// against the owner's full 2026-08-10 hand scan of the S6 separation
    /// displays (899 labels over 230 events; doc pr/57 sec 14): per-pair
    /// result bad 124/127 kept-connected, good 154/156 kept-separated, with
    /// both good misses being label triangles that no relaxation can satisfy
    /// (the owner's "prefer connectivity" tie-break decides them). Branches
    /// map one-to-one to the owner's stated bad-separation classes: dead-W
    /// band, W-robustness with plane-dependent gap allowance, direction
    /// consistency on substantial pairs, prolonged induction signal, and
    /// two-view footprint co-location. It can only REMOVE kills S6 made --
    /// composed as `killed && !two_d_rescue_ok(...)` -- so with the flavor
    /// unselected nothing changes anywhere.
    ///
    /// Pure; exposed for doctests.
    bool two_d_rescue_ok(const S6RescueInput& in);

    /// doc pr/64 round 4: inputs to the S6 W-plane long-track exception,
    /// measured at the kill site for one S6-KILLED candidate edge (any
    /// distance S6 judges, unlike the <= 5 cm rescue above).  Analysis units
    /// (cm, degrees, counts) so the predicate is a line-by-line port of the
    /// Python rule it was fitted with (sbnd_xin/scripts/analysis/pr64/
    /// wgate_sweep.py, rule R2w tightened + R2d).  Sentinels:
    /// ab_global_deg = 90 means "no global axis measurable" (fails the angle
    /// cuts); tmax_cm = 1e9 means "no transverse RMS measurable" (fails the
    /// thinness cuts); dead_w = -1 means "dead-W not measured" (the dead-W
    /// count is only computed for the <= 5 cm rescue population, which
    /// contains all of R2d's dis < 3 band -- fails the R2d branch).
    struct S6WTrackInput {
        double dis_cm = 0;          ///< candidate edge 3D length
        bool w_gap = false;         ///< gap[2]: W plane (1,1) BFS gap
        bool w_sole_vote = false;   ///< gap[2] && !(gap[0]&&!excuse_u) && !(gap[1]&&!excuse_v)
        int npmin = 0;              ///< smaller component's point count (capped 20000)
        double lmin_cm = 0;         ///< smaller component's principal-axis extent
        double tmax_cm = 1e9;       ///< larger component transverse RMS (global PCA)
        double ab_global_deg = 90;  ///< angle between the two GLOBAL principal axes
        int dead_w = -1;            ///< dead W wires across the W seed span (-1 unmeasured)
    };

    /// doc pr/64 round 4: W-plane long-track exception of the
    /// "relaxed_strict_img_2d_rescue_long_wtrack" graph flavor.  True => an
    /// S6-killed candidate is kept after all.  Composed exactly like the
    /// rescue above -- `killed && !two_d_w_track_ok(...)` -- so it can only
    /// REMOVE kills S6 made, and with the flavor unselected nothing changes
    /// anywhere.  It deliberately does NOT touch two_d_connectivity_bad's
    /// "W is never excused" rule: W stays a killing plane everywhere else.
    ///
    /// Root cause this closes: S6/S7 excuse prolonged-signal gaps on U/V
    /// only (the excuse channel is an induction-plane physical effect), so
    /// a genuinely straight long track with a small mysterious W-only hole
    /// (dead/low-response W column, drift-parallel topology) is broken by
    /// the ONE plane with no excuse channel -- SBND 18259-174224,
    /// 18255-276836, 18255-314507 (doc pr/64).  The existing rescue's W
    /// branch cannot catch these: its collinearity feature is the LOCAL
    /// break-point PCA axis, noisy exactly at a break, while a long straight
    /// track's WHOLE-component PCA axis stays tight (174224: local 20.7 deg
    /// vs global 4.2 deg).  The global axis is what this predicate uses --
    /// and it is cheaper than the local one: the per-component eigenvector
    /// is already computed (and was discarded) by the s6_comp_stat cache.
    ///
    /// Rule (owner-selected "tightened gate", zero regressions on the
    /// 899-label hand scan -- doc pr/64 round 4):
    ///   R2d: w_gap && dead_w >= 3 && npmin >= 20 && dis < 3 cm
    ///   R2w: w_sole_vote && lmin > 6 cm && npmin >= 50
    ///        && tmax < 2 cm && ab_global < 25 deg
    ///        && (tmax < 1.7 cm || ab_global < 6 deg)   [tightening]
    /// The tightening protects evt122660's owner-labelled good pair
    /// (tmax 1.73, ab 7.6) at the cost of two marginal recoveries (64959,
    /// 172656; ab 22/18 deg -- not clean-long-track topology).
    ///
    /// Pure; exposed for doctests.
    bool two_d_w_track_ok(const S6WTrackInput& in);

    /// doc pr/62 S7: inputs to the long-edge corridor-connectivity decision,
    /// measured at the kill site for ONE candidate edge whose 3D length is at
    /// or above the S6 cut-off. Sentinels: has_plane[p] false means that
    /// plane produced no seeds on one side and never votes; budget_hit[p]
    /// true means the bounded search was inconclusive and that plane
    /// abstains (S7 fails OPEN, unlike S6, which fails closed).
    struct S7CorridorInput {
        double dis_cm = 0;                          ///< candidate edge 3D length
        bool has_plane[3]  = {false, false, false};  ///< plane had seeds on both sides
        bool gap[3]        = {false, false, false};  ///< corridor search found no path
        bool budget_hit[3] = {false, false, false};  ///< breaker fired => inconclusive
        bool excuse_u = false, excuse_v = false;     ///< wire-parallel excuses (reused from S1-S6)
        double gap_cm[3] = {0, 0, 0};                ///< censused gap length along the corridor
    };

    /// Kill predicate of the "relaxed_strict_img_2d_rescue_long" (and
    /// "..._long2") graph flavors (doc pr/62: S7, "long-edge corridor
    /// connectivity"). True if at least min_gapped_planes non-excused,
    /// non-abstaining planes show the two candidate components' fired-pixel
    /// footprints NOT connected by any live-or-dead-channel path that stays
    /// inside a narrow corridor around the straight q1->q2 line, projected
    /// into that plane's (wire_index, time_slice) grid via the candidate's
    /// own two endpoint lattice cells (not an arbitrary-point projection --
    /// see connect_graph_relaxed_strict.cxx's S7Corridor comment).
    ///
    /// Root cause this closes: S6 (two_d_connectivity_bad) is the only test
    /// that ever verifies real 2D contiguity, and it hard-returns above
    /// s6_dis_cap = 30cm (connect_graph_relaxed_strict.cxx, two_d_gap_kill)
    /// because its open-rectangle flood fill costs O(D^2) and its cell
    /// budget breaker fails closed, so raising that cap would trade missed
    /// gaps for false kills of large sparse real objects. Every candidate at
    /// or above 30cm is therefore governed by S1-S3 alone -- and S1-S3's
    /// test_good_point (Facade_Grouping.cxx) checks U/V/W INDEPENDENTLY via
    /// three per-plane 2D kd radius queries, never intersected in 3D, so
    /// three planes can each see charge from a DIFFERENT nearby track and
    /// every 1cm sample reads "good" with nothing physically there. On evt
    /// 142421 cluster 7 (21 components, 210 candidate pairs) 26 pairs
    /// survive S1-S3 and 14 of those are >= 30cm; three of the 14 (30.03,
    /// 36.87, 37.39cm) are the MST bridges that wrongly reconnect a
    /// 17-point island.
    ///
    /// Restricting the search to a corridor rather than a rectangle is what
    /// buys O(D): the corridor holds O(D * halfwidth) cells, so the cell
    /// budget stops being a live constraint and the breaker stops silently
    /// manufacturing kills (it fails OPEN here -- see S7CorridorInput). It
    /// adds no assumption S1-S3 does not already make -- S1-S3 samples that
    /// very straight line -- it merely checks contiguity along it properly
    /// instead of via three independent per-plane radius queries.
    ///
    /// dis_floor_cm is the seam with S6 and MUST equal
    /// connect_graph_relaxed_strict.cxx's s6_dis_cap / s7_dis_floor_cm: S7 is
    /// identically false below it, S6 is identically false at or above it,
    /// so no candidate is ever judged by both and neither test's shipped
    /// operating point is disturbed. doctest_long_corridor.cxx pins the
    /// value.
    ///
    /// Excusals are S6's, reused verbatim rather than re-derived: W is never
    /// excused; U/V are excused when the candidate path already tests as
    /// quasi-parallel to that plane's wire direction.
    ///
    /// min_gapped_planes defaults to 1, mirroring two_d_connectivity_bad's
    /// owner rule; a conservative >= 2 variant is one argument away and is
    /// the flavor this band's second operating point ships as
    /// ("relaxed_strict_img_2d_rescue_long2") -- there are no hand-scan
    /// labels in this distance band yet (unlike S6's), so BOTH operating
    /// points are run and censused rather than one being guessed as final.
    /// gap_floor_cm defaults to 0 (INACTIVE) on purpose: every other
    /// operating point in this header carries a scan justification, and no
    /// scan yet justifies a gap-length floor here. gap_cm is measured and
    /// censused so one can be fitted later; do not switch it on before it is.
    ///
    /// Pure; exposed for doctests. The corridor search itself lives in
    /// connect_graph_relaxed_strict.cxx (not pure -- it queries live
    /// Facade::Grouping state).
    bool long_corridor_bad(const S7CorridorInput& in,
                           int min_gapped_planes = 1,
                           double dis_floor_cm = 30.0,
                           double gap_floor_cm = 0.0);

}

#endif

