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
        double total{0};
        bool   dl_winner{false};
        /// Candidate removed from the acceptance by dl_vtx_swap_guard before
        /// scoring, so its terms are all zero BY OMISSION, not by measurement.
        bool   skipped_by_swap_guard{false};

        /// Host cluster total track length in cm (the s_clen input), kept
        /// because it is the strongest geometric signal in the composite.
        double host_length{0};
    };

    struct DLVoxelRow {
        int rank{-1};
        double x{0}, y{0}, z{0};   ///< cm
        double dl_score{0};
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

        std::vector<DLVoxelRow> voxels;
        std::vector<VertexScoreRow> rows;

        /// The vertex carrying kNeutrinoVertex when the board was stashed.
        /// The stash happens after snap_main_vertex_to_kink and the final
        /// improve_vertex, so this is the position the display draws.
        int    final_vertex_id{-1};
        double final_x{0}, final_y{0}, final_z{0};

        void clear() { *this = VertexScoreboard{}; }
    };

}  // namespace WireCell::Clus::PR

#endif
