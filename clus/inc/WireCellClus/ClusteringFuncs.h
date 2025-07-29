/**
   This header provides various free functions used in clustering.

   Some implementations may be found in clustering_*.cxx.  The rest are in
   ClusteringFuncs.cxx.

 */

#ifndef WIRECELLCLUS_CLUSTERINGFUNCS
#define WIRECELLCLUS_CLUSTERINGFUNCS

#include "WireCellClus/MultiAlgBlobClustering.h"
#include "WireCellClus/Facade.h"
#include "WireCellClus/IPCTransform.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Graphs.h"

#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMdataset.h"
#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Persist.h"



#include <string>
#include <fstream>

namespace WireCell::Clus::Facade {
    using namespace WireCell::PointCloud::Tree;

    /// Some clustering functions define and react to flags defined on cluster
    /// facades.  We name these flags with strings in this namespace to assure
    /// consistency, add documentation/comments and avoid typos.  Note,
    /// merge_clusters() will call fresh.from(other) to forward flags.
    namespace Flags {

        /// Indicates the cluster is a "live" cluster and connected to a "dead"
        /// cluster.
        inline const std::string live_dead = "live_dead";

        /// Indicates the cluster has a flash coincident with beam timing
        inline const std::string beam_flash = "beam_flash";
        
        /// Indicates the cluster is tagged as through-going muon (TGM)
        inline const std::string tgm = "tgm";
        
        /// Indicates the cluster is tagged as low energy
        inline const std::string low_energy = "low_energy";
        
        /// Indicates the cluster is tagged as light mismatch (LM)
        inline const std::string light_mismatch = "light_mismatch";
        
        /// Indicates the cluster is tagged as fully contained
        inline const std::string fully_contained = "fully_contained";
        
        /// Indicates the cluster is tagged as short track muon (STM)
        inline const std::string short_track_muon = "short_track_muon";
        
        /// Indicates the cluster has full detector dead region
        inline const std::string full_detector_dead = "full_detector_dead";

        // main cluster
        inline const std::string main_cluster = "main_cluster";

        // associated cluster
        inline const std::string associated_cluster = "associated_cluster"; 

        /// This flag is set by ClusteringTaggerCheckSTM algorithm when specific STM conditions are met
        inline const std::string STM = "STM";
    }

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


    // This function will produce a new cluster in the grouping corresponding to
    // each connected component in the cluster_connectivity_graph_t that as two
    // or more clusters.  Any cluster in a single-cluster component is simply
    // left in place in the grouping.
    //
    // Each new cluster will be given the children (blob nodes) of the clusters
    // in the associated connected component.  These now empty clusters will be
    // removed from the grouping and discarded.
    //
    // If both aname and pcname are given then a representation of the previous
    // clustering of blob nodes will be stored in the new cluster.  This
    // connected component (cc) array is in child-node-order and its integer
    // value counts which original cluster donated the blob to the new cluster.
    //
    // The fresh.from(other) is called to transfer flags, scope and possibly
    // other bits of information.
    //
    // Pointers to the newly created cluster node facades are returned.  These
    // are loaned.  As usual, the cluster node owns the facade and these nodes
    // are in turn owned by the grouping node.
    std::vector<Cluster*> merge_clusters(cluster_connectivity_graph_t& g, // 
                                         Grouping& grouping,
                                         const std::string& aname="",
                                         const std::string& pcname="perblob");

    

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


    // Calculate PCA direction for a set of points around a center point
    geo_vector_t calc_pca_dir(const geo_point_t& center, const std::vector<geo_point_t>& points);


}  // namespace WireCell::Clus::Facade

#endif

