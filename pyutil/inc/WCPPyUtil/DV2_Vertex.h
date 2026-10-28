#ifndef WCPPYUTIL_DV2_VERTEX_H
#define WCPPYUTIL_DV2_VERTEX_H

#include "WCPPyUtil/config.h"

#include <vector>
#include <string>

namespace WCPPyUtil {

    /// Call the DV2-MAX neutrino-vertex finder (pyutil/python/DV2_Vertex.py).
    ///
    /// input:          4 arrays x, y, z, q -- the DL input cloud exactly as
    ///                 determine_overall_main_vertex_DL builds it (cm, raw reco
    ///                 frame, q = dQ * q_scale + q_offset); the first
    ///                 n_vertex_rows points are the PR-graph vertex rows.
    /// trad:           the traditional main vertex (cm, 3 values; NaN = none).
    /// Returns [x, y, z, score] (cm, raw reco frame; score = calibrated
    /// P(< 1 cm)) or an empty vector when the model makes no prediction.
    /// Throws std::runtime_error on any Python failure.
    std::vector<FLOAT> DV2_Vertex(const std::string& module, const std::string& function,
                                  const std::string& weights,
                                  const std::vector<std::vector<FLOAT> >& input, int n_vertex_rows,
                                  const std::vector<FLOAT>& trad, double q_scale, double q_offset,
                                  const std::string& dtype = "float32");
}
#endif
