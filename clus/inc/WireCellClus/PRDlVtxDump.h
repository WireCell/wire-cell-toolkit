/** One DL-vertex network call, recorded exactly as it happened (ai-helper issue 35).
 *
 * Filled by PatternAlgorithms::determine_overall_main_vertex_DL when the
 * TaggerCheckNeutrino knob "dl_vtx_dump" is on (default off), for BOTH network
 * calls of a candidate: the dual-chain OFF pass (exclusion-free fit) and the
 * production pass.  Carried per candidate on TrackFitting (dlvtx_calls) and
 * written by SbndPrMagnifyTrackingVisitor as T_dlvtx_call / T_dlvtx_cloud.
 *
 * The cloud is the network input vec_xyzq bit for bit -- the same floats
 * WCPPyUtil::SCN_Vertex receives -- so a standalone call of
 * pyutil/python/SCN_Vertex.py on it reproduces the payload.  Recording only:
 * nothing here feeds back into the reconstruction.
 */
#ifndef WIRECELLCLUS_PRDLVTXDUMP
#define WIRECELLCLUS_PRDLVTXDUMP

#include <string>
#include <vector>

namespace WireCell::Clus::PR {

    struct DlVtxCall {
        // ---- which call
        std::string pass{"prod"};      ///< "prod" (production pass) | "off" (dual-chain OFF pass, snap mode: through determine_overall_main_vertex_DL) | "off-voxels" (dual-chain voxels/union mode: the OFF graph's own top-K call in dual_chain_scn_voxels; no decision fields)
        int  status{0};                ///< 0 ok; 1 the network threw (payload empty); 2 unexpected payload size (legacy path)
        int  top_k{0};                 ///< top_k requested from SCN_Vertex (1 = legacy argmax)
        bool rerank{false};            ///< dl_vtx_rerank in force
        // ---- the exact network input (cm, float32, build order: vertex block first)
        std::vector<float> x, y, z, q;
        int    n_vertex_rows{0};       ///< leading rows that are PR-graph vertices; the rest are segment-interior fit points
        bool   cloud_no_exclusion{false}; ///< dl_vtx_cloud_no_exclusion: the cloud came from an exclusion-free refit of this pass's graph
        double q_scale{0}, q_offset{0};///< q = dQ * q_scale + q_offset
        // ---- the network output
        std::vector<float> payload;    ///< raw return of THIS call's inference: legacy [x,y,z]; rerank [x,y,z,score]*K (cm)
        bool payload_from_off{false};  ///< dual-chain "voxels" mode: payload is the OFF pass's, no inference here
        int  n_off_voxels{0};          ///< dual-chain "union" mode: OFF voxels pooled with the payload AFTER it was recorded (the decision below saw payload + these)
        // ---- the decision this call fed (cm)
        // Row indices are positions inside the cloud's vertex block (0 .. n_vertex_rows-1), -1 = not a cloud vertex.
        // They identify the vertex even when cloud_no_exclusion restored the coordinates after the cloud was built.
        bool   trad_valid{false};      ///< traditional main vertex of the main cluster before the DL
        double trad_x{0}, trad_y{0}, trad_z{0};
        int    trad_row{-1};
        bool   rerank_valid{false};    ///< this call's OWN pick (legacy argmax snap / rerank winner) passed its gate, BEFORE the dual-chain snap
        double rerank_x{0}, rerank_y{0}, rerank_z{0};
        int    rerank_row{-1};
        bool   accepted{false};        ///< DL vertex accepted (flag_pass) after snap and veto: the main vertex was switched to it
        double dl_x{0}, dl_y{0}, dl_z{0};  ///< the accepted candidate vertex (when accepted)
        int    dl_row{-1};
        bool   dual_transferred{false};///< dual chain snap moved production's pick (set in the snap block; a later veto can still clear `accepted`)
        bool   two_end_veto{false};    ///< the protected two-end-break vertex vetoed the DL choice (accepted = false after a pass)
        // ---- the dual-chain hint this (production) call was given, cm; "off" rows: none
        bool   hint_valid{false};
        double hint_x{0}, hint_y{0}, hint_z{0};
    };

}  // namespace WireCell::Clus::PR

#endif
