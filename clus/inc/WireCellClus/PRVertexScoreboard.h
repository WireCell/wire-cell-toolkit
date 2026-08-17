#ifndef WIRECELL_CLUS_PR_VERTEX_SCOREBOARD
#define WIRECELL_CLUS_PR_VERTEX_SCOREBOARD

#include <string>
#include <vector>

namespace WireCell::Clus::PR {

    /** Per-event record of HOW the neutrino vertex was chosen.
     *
     * doc sbnd_xin/docs/pr/75 (the scoreboard doc pr/52 §5.1 asked for).
     *
     * DIAGNOSTIC ONLY.  Nothing in clus reads this; it exists so
     * PrDisplayDump can emit, beside each PR-graph vertex, the numbers the
     * two selectors actually compared -- which is what a hand scan needs in
     * order to become tuning input rather than a bare correct/incorrect
     * verdict.
     *
     * Filled only when PatternAlgorithms::m_vertex_scoreboard is true (knob
     * `vertex_scoreboard`, C++ default false).  With the knob off `filled`
     * stays false and every container stays empty: a consumer must read that
     * as "NO SCOREBOARD WAS TAKEN", never as "no candidates existed".
     */
    struct VertexScoreRow {
        /// cluster_id*1000 + graph_index -- the same encoding PrDisplayDump's
        /// vertices[].id and the Bee particle-flow tree use, so a consumer
        /// joins on this and nothing else.
        int vertex_id{-1};
        int cluster_id{-1};
        /// Position in cm, fit().point when valid else wcpt().point -- the
        /// same choice both selectors make when they score.
        double x{0}, y{0}, z{0};

        // ---- stage 1, per-cluster: compare_main_vertices (the additive,
        // prototype-ported scorer).  trad_scored is false for clusters that
        // took the all-showers branch (compare_main_vertices_all_showers is
        // NOT additive and produces no comparable number) and for vertices
        // that never entered a candidate list.
        bool   trad_scored{false};
        double trad_score{0};
        bool   trad_winner{false};

        // ---- stage 2, global: the DL rerank.  dl_snapped is true only for
        // the <= dl_vtx_top_k vertices some DL voxel snapped to; every other
        // vertex has no DL opinion at all (NOT a zero score).
        bool   dl_snapped{false};
        int    voxel_rank{-1};
        double dl_score{0};      ///< raw net score, sigmoid(vtx)-sigmoid(bg), [-1,+1]
        double snap_dis{0};      ///< cm, voxel to vertex

        // The seven composite terms, exactly as summed at
        // NeutrinoVertexFinder.cxx (see the `DL rerank cand` TRACE line).
        double s_dl{0}, s_snap{0}, s_fwd_z{0}, s_clen{0}, s_isol{0}, s_main{0}, s_fv{0};
        /// doc pr/89 C2 -- eighth term, rule-1 outgoing-prong topology
        /// (knob dl_vtx_topo_weight, C++ default 0 = term absent).  Emitted
        /// by PrDisplayDump only when board.topo_used, so knob-off calib
        /// JSON is byte-identical.  topo_frac -1 = no decisive attached
        /// track (NO OPINION, never a zero vote -- pr/88 P6).
        double s_topo{0};
        double topo_frac{-1};
        int    topo_votes{0};
        double total{0};
        bool   dl_winner{false};
        /// Candidate removed from the acceptance by dl_vtx_swap_guard before
        /// scoring, so its terms are all zero BY OMISSION, not by measurement.
        bool   skipped_by_swap_guard{false};

        /// Host cluster total track length in cm (the s_clen input), kept
        /// because it is the strongest geometric signal in the composite.
        double host_length{0};

        // ---- doc sbnd_xin/docs/pr/79 §10 -- dl_vtx_harvest (knob, C++
        // default false; effective only with vertex_scoreboard).  Stage-1
        // features compare_main_vertices computes and then discards, recorded
        // for future learned-component training on the LIVE distribution.
        // hv_filled false => this row never passed through the harvest block
        // (all-showers branch, single-candidate branch, never a candidate),
        // NOT "the features were zero".
        bool   hv_filled{false};
        int    hv_n_proton_in{-1};   ///< protons entering via weak/undirected neighbour
        int    hv_n_proton_out{-1};  ///< protons leaving
        double hv_z_prior{0};        ///< the subtracted (z - min_z)/scale term
        int    hv_n_tracks{-1};      ///< non-shower segments at the vertex
        int    hv_n_showers{-1};     ///< shower segments; degree = tracks + showers
        bool   hv_in_fv{false};      ///< the +0.5 fiducial test result
        double hv_conflicts{-1};     ///< calc_conflict_maps scalar (-1 = not computed)
        double hv_reduced_chi2{-1};  ///< vertex fit reduced chi2 (-1 = fit invalid)
    };

    struct DLVoxelRow {
        int rank{-1};
        double x{0}, y{0}, z{0};   ///< cm
        double dl_score{0};
    };

    /// doc pr/79 §10 -- one candidate of compare_main_vertices_global (the
    /// scorer that decides cluster swaps on the reject route).  Kept separate
    /// from rows[]: global candidates may sit on clusters with no stage-1 or
    /// DL row, and the global score is a different stage with different
    /// semantics than trad_score.  Filled only under dl_vtx_harvest.
    struct GlobalVertexScoreRow {
        int vertex_id{-1};
        int cluster_id{-1};
        double x{0}, y{0}, z{0};   ///< cm, fit-else-wcpt (pr75 convention)
        double score{0};
        bool is_main_cluster{false};
        bool in_fv{false};
        bool winner{false};
    };

    struct VertexScoreboard {
        /// False => the knob was off.  See the class comment.
        bool filled{false};

        /// Which of doc pr/52 §1's routes chose the final vertex:
        ///   "dl-not-run"        DL never invoked (weights empty, or the
        ///                       requested path failed Persist::resolve --
        ///                       distinguish with weights_missing below)
        ///   "dl-failed"         SCN import/inference threw
        ///   "dl-rerank-accept"  rerank winner cleared min_accept_score
        ///   "dl-rerank-reject"  rerank winner scored below it
        ///   "dl-legacy-accept" / "dl-legacy-reject"  rerank=false branch
        ///   "dl-veto-protected" pr/48 kProtectedBreak veto kept the
        ///                       traditional vertex over the DL choice
        ///   "dl-no-candidates"  nothing snapped
        std::string route{"dl-not-run"};
        /// True when dl_weights was configured non-empty but Persist::resolve
        /// failed -- the one way route 3 is distinguishable from route 1,
        /// which is otherwise invisible at the call site because the failed
        /// resolve leaves m_dl_weights empty.
        bool weights_missing{false};

        bool dl_ran{false};      ///< inference was attempted
        bool dl_rerank{false};   ///< the rerank branch (SBND production) not the legacy one
        bool dl_accepted{false}; ///< the DL choice replaced the traditional one

        double dl_best_score{0};
        double dl_min_accept_score{0};
        double dl_score_scale{0};
        int    dl_top_k{0};

        /// doc pr/89 C2 -- true only when dl_vtx_topo_weight != 0.  Gates
        /// the s_topo/topo_* key emission in PrDisplayDump the same way
        /// `harvest` gates the hv_* keys: knob-off calib JSON byte-identical.
        bool   topo_used{false};
        double topo_weight{0};
        double topo_center{0};

        std::vector<DLVoxelRow> voxels;
        std::vector<VertexScoreRow> rows;

        /// The vertex carrying kNeutrinoVertex when the board was stashed.
        /// The stash happens after snap_main_vertex_to_kink and the final
        /// improve_vertex, so this is the position the display draws.
        int    final_vertex_id{-1};
        double final_x{0}, final_y{0}, final_z{0};

        // ---- doc pr/79 §10 -- dl_vtx_harvest.  All below filled only when
        // the harvest knob (AND vertex_scoreboard) is on; emitted by
        // PrDisplayDump only when `harvest` is true, so harvest-off calib
        // JSON is byte-identical to pre-knob output.
        bool harvest{false};
        /// The EXACT live SCN input cloud (pre-voxelization; voxelization at
        /// 0.5 cm happens in pyutil SCN_Vertex.py from these floats), cm and
        /// scaled charge, in build order.  The leading cloud_n_vertex_rows
        /// entries align 1:1 with cloud_vertex_ids (joins rows[].vertex_id);
        /// the rest are segment-interior fit points.  ORDER IS LOAD-BEARING
        /// -- never sort.  float, not double: float(x) round-trips exactly
        /// through the JSON double, reproducing the live network input.
        std::vector<float> cloud_x, cloud_y, cloud_z, cloud_q;
        int cloud_n_vertex_rows{0};
        std::vector<int> cloud_vertex_ids;
        double cloud_q_scale{0}, cloud_q_offset{0};  ///< dQdx scale/offset used for q
        /// map_cluster_main_vertices[main_cluster] at DL entry -- the
        /// traditional answer the reject route would inherit, before any
        /// global cluster swap.
        int hv_trad_main_vertex_id{-1};
        bool global_ran{false};     ///< compare_main_vertices_global executed
        int  global_winner_id{-1};
        std::vector<GlobalVertexScoreRow> global_rows;
        /// Vertices chosen by determine_main_vertex's single-candidate
        /// branch (compare_main_vertices never ran on their cluster) and
        /// winners of the non-additive compare_main_vertices_all_showers.
        /// Board-level ID LISTS, not rows: harvest must not add rows[]
        /// entries the pr/75 baseline board would not have, or rows-based
        /// candidate-set identity across arms breaks.
        std::vector<int> hv_single_candidate_ids;
        std::vector<int> hv_all_showers_winner_ids;

        void clear() { *this = VertexScoreboard{}; }
    };

}  // namespace WireCell::Clus::PR

#endif
