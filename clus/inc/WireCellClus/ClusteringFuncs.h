#ifndef WIRECELLCLUS_CLUSTERINGFUNCS
#define WIRECELLCLUS_CLUSTERINGFUNCS

#include "WireCellClus/MultiAlgBlobClustering.h"
#include "WireCellClus/Facade.h"
#include "WireCellClus/IPCTransform.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/IClusteringMethod.h"

#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMdataset.h"
#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Graph.h"



#include <fstream>

namespace WireCell::Clus::Facade {
    using namespace WireCell::PointCloud::Tree;

    struct ClusterLess {
        bool operator()(const Cluster* a, const Cluster* b) const {
            return cluster_less(a, b);
        }
    };

    using cluster_set_t = std::set<const Cluster*>;

    // Each vertex of a cluster connectivity graph represents the index of a
    // cluster in some collection.
    using cluster_connectivity_graph_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, int>;

    using cluster_vector_t = std::vector<Cluster*>;


    // clustering_util.cxx
    //
    // This function will produce new clusters in live_clusters.  The children
    // of each "fresh" cluster will be those of the "donor" clusters that are
    // connected according to the cluster_connectivity_graph_t.  The "fresh"
    // cluster will be added to and the "donor" clusters will be removed from
    // "known_clusters".  The "donor" clusters will also be removed from
    // live_clusters.
    //
    // If both aname and pcname are given then store a cc array in any newly
    // created clusters holding the merged set of blobs.  The cc array will
    // arbitrarily label each blob with a number corresponding to the original
    // cluster which was parent to the blob (and which is destroyed after this
    // function).
    //
    // See above for cluster_connectivity_graph_t.
    void merge_clusters(cluster_connectivity_graph_t& g, // 
			Grouping& grouping,
			cluster_set_t& known_clusters, // in/out
                        const std::string& aname="", const std::string& pcname="perblob");
    
    /**
     * Extract geometry information from a grouping
     * @param grouping The input Grouping object
     * @param dv Detector geometry provider
     * @return Tuple of (drift_direction, angle_u, angle_v, angle_w)
     */
    std::tuple<geo_point_t, double, double, double> extract_geometry_params(
        const Grouping& grouping,
        IDetectorVolumes::pointer dv);

    std::vector<std::pair<geo_point_t, const Blob*>> get_strategic_points(const Cluster& cluster);

    //helper function ..
    double Find_Closest_Points(const Cluster& cluster1,
			       const Cluster& cluster2,
			       double length_1,
			       double length_2,
			       double length_cut,
			       geo_point_t& p1_save, // output
			       geo_point_t& p2_save,  // output
                   bool flag_print = false
			       );
			       
    





    // These Judge*() functions are used by multiple clustering methods.  They
    // are defined in clustering_separate.cxx.

    // time_slice_length is length span for a slice
    bool JudgeSeparateDec_1(const Cluster* cluster, const geo_point_t& drift_dir, const double length);
    /// @attention contains hard-coded distance cuts
    /// @param boundary_points return the boundary points
    /// @param independent_points return the independent points
    bool JudgeSeparateDec_2(const Cluster* cluster, IDetectorVolumes::pointer dv, const geo_point_t& drift_dir,
                               std::vector<geo_point_t>& boundary_points, std::vector<geo_point_t>& independent_points,
                               const double cluster_length);
    



    // this is used only by ClusteringSeparate but keep it public for symmetry with Separate_2.
    std::vector<Cluster *> Separate_1(const bool use_ctpc, Cluster *cluster,
                                      std::vector<geo_point_t> &boundary_points,
                                      std::vector<geo_point_t> &independent_points,
                                      double length, geo_point_t dir_cosmic, geo_point_t dir_beam, 
                                      IDetectorVolumes::pointer dv, 
                                      IPCTransformSet::pointer pcts,
                                      const Tree::Scope& scope);

    // This is used by multiple clustering methods.
    std::vector<int> Separate_2(Cluster *cluster, 
                                const Tree::Scope& scope,
                                const double dis_cut =  5*units::cm);


}  // namespace WireCell::Clus::Facade

#endif
