/** The Bee clustering-layer points of one cluster, and the DL-vertex touch test
 * (ai-helper issue 44).
 *
 * bee_cluster_points() is the default (scoped-view) branch of
 * MultiAlgBlobClustering::fill_bee_points_from_cluster, factored out so the Bee
 * writer and the tracking-pr.root writer (T_dlvtx_space) produce the same points
 * by construction: one point per scoped-view point, blob by blob, in the
 * scoped view's order, with q = the mean wire charge over the planes that are
 * not dead at that point (the prototype formula).
 *
 * bee_touching_clusters() is Sne's r7_space selection rule
 * (sbnd-vertexer/r7/export_r7_space.py): a cluster is taken whole when ANY of
 * its Bee points lies within touch_cm of ANY point of the DL network's input
 * cloud.
 */
#ifndef WIRECELLCLUS_BEECLUSTERPOINTS
#define WIRECELLCLUS_BEECLUSTERPOINTS

#include "WireCellUtil/Point.h"
#include "WireCellUtil/PointTree.h"

#include <functional>
#include <vector>

namespace WireCell::Clus {
    namespace Facade { class Cluster; }

    /// The Bee layer's scope gate (fill_bee_points_from_cluster): filter 1 =
    /// the cluster's runtime scope filter, 0 = every cluster, -1 = the
    /// complement; any other value = every cluster.
    bool bee_cluster_in_scope(const Facade::Cluster& cluster,
                              const PointCloud::Tree::Scope& scope, int filter);

    /// Visit the cluster's points through its scoped view.  fn(blob_index,
    /// point, q): blob_index = position in the scoped view's pcs(), point in
    /// WCT units.  with_charge false => q = 0 and no charge lookup (cheap).
    void bee_cluster_points(const Facade::Cluster& cluster,
                            const PointCloud::Tree::Scope& scope, bool with_charge,
                            const std::function<void(size_t blob_index, const Point& pt, double q)>& fn);

    /// Sne's touch test.  clusters: the candidate clusters (each with its Bee
    /// points in cm, flattened xyz triplets); cloud: the network input in cm.
    /// Returns the indices into clusters that have a point within touch_cm
    /// (distance <= touch_cm) of a cloud point, in increasing order.
    std::vector<size_t> bee_touching_clusters(const std::vector<std::vector<float>>& clusters_xyz_cm,
                                              const std::vector<float>& cloud_x,
                                              const std::vector<float>& cloud_y,
                                              const std::vector<float>& cloud_z,
                                              double touch_cm);
}

#endif
