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
        std::string pass{"prod"};      ///< "prod" (production pass) | "off" (dual-chain OFF pass)
        int  top_k{0};                 ///< top_k requested from SCN_Vertex (1 = legacy argmax)
        bool rerank{false};            ///< dl_vtx_rerank in force
        // ---- the exact network input (cm, float32, build order: vertex block first)
        std::vector<float> x, y, z, q;
        int    n_vertex_rows{0};       ///< leading rows that are PR-graph vertices; the rest are segment-interior fit points
        double q_scale{0}, q_offset{0};///< q = dQ * q_scale + q_offset
        // ---- the network output
        std::vector<float> payload;    ///< raw return: legacy [x,y,z]; rerank [x,y,z,score]*K (cm)
        bool payload_from_off{false};  ///< dual-chain "voxels" mode: payload is the OFF pass's, no inference here
        // ---- the decision this call fed (cm)
        bool   trad_valid{false};      ///< traditional main vertex of the main cluster before the DL
        double trad_x{0}, trad_y{0}, trad_z{0};
        bool   accepted{false};        ///< DL vertex accepted (flag_pass): the main vertex was switched to it
        double dl_x{0}, dl_y{0}, dl_z{0};  ///< the accepted candidate vertex (when accepted)
        bool   dual_transferred{false};///< dual chain snap moved production's pick
    };

}  // namespace WireCell::Clus::PR

#endif
