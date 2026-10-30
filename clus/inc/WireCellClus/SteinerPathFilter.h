/** The per-terminal decision of Steiner::Grapher::filter_by_path_constraints.

    pr_ml doc 34 sec 8.18.  With own_wpid false this is the historical behaviour, bit for bit, including its two defects:

    - Cluster::wpid() takes a point and D3Vector's constructor is implicit, so an index compiles as Point(idx, 0, 0) and the terminal gets the
      wire plane of the cluster point nearest that position, not its own.
    - DynamicPointCloud::get_closest_2d_point_info returns -1 when the path has no point in that (plane, face, apa) tree, and -1 passes
      "< threshold", so the prototype's "close in 2-D and far in 3-D" collapses to a bare cut on the 3-D distance.

    Together they delete every terminal beyond the 3-D threshold in whole anodes of a multi-anode cluster.  With own_wpid true the terminal reads
    its own wire plane and an empty tree is not close: there is no 2-D evidence against the terminal. */
#ifndef WIRECELLCLUS_STEINERPATHFILTER
#define WIRECELLCLUS_STEINERPATHFILTER

#include "WireCellClus/Facade_Cluster.h"

namespace WireCell::Clus::Steiner {

    inline WirePlaneId path_filter_terminal_wpid(const Facade::Cluster& cluster, size_t terminal_idx, bool own_wpid)
    {
        if (own_wpid) return cluster.wire_plane_id(terminal_idx);
        return cluster.wpid(terminal_idx);   // historical: compiles as Point(terminal_idx, 0, 0)
    }

    /// Close in at least two of the three planes.  A negative distance is the empty-tree sentinel: close only in the historical arm.
    inline bool path_filter_close_in_2d(const double (&dis_2d)[3], double threshold_2d, bool own_wpid)
    {
        const auto close = [&](double d) { return d < threshold_2d && (!own_wpid || d >= 0); };
        const bool close_u = close(dis_2d[0]), close_v = close(dis_2d[1]), close_w = close(dis_2d[2]);
        return (close_u && close_v) || (close_u && close_w) || (close_v && close_w);
    }

    /// The prototype's rule: remove a terminal close to the path in 2-D but far from it in 3-D.
    inline bool path_filter_removes(double dis_3d, const double (&dis_2d)[3], double threshold_2d, double threshold_3d, bool own_wpid)
    {
        return path_filter_close_in_2d(dis_2d, threshold_2d, own_wpid) && (dis_3d > threshold_3d);
    }

}  // namespace WireCell::Clus::Steiner

#endif
