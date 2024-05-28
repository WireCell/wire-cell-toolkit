#include "WireCellImg/MultiAlgBlobClustering.h"
#include "WireCellImg/PointCloudFacade.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Persist.h"
#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMdataset.h"
#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/Graph.h"


#include <fstream>

namespace WireCell::PointCloud::Facade {
    using namespace WireCell::PointCloud::Tree;

    using cluster_length_map_t = std::map<const Cluster::const_pointer, double>;
    using const_cluster_set_t = std::set<Cluster::const_pointer>;
    using cluster_connectivity_graph_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, int>;
    using live_clusters_t = Cluster::vector;

    // merging clustering function
    void merge_clusters(cluster_connectivity_graph_t& g,
			Points::node_ptr& root_live,                                   // in/out
			live_clusters_t& live_clusters,
			cluster_length_map_t& cluster_length_map,  // in/out
			const_cluster_set_t& cluster_connected_dead,            // in/out
			const TPCParams& tp,                                           // common params
			const std::set<Cluster::const_pointer >& cluster_to_be_deleted
			);
    
    
    // first function ...
    void clustering_live_dead(Points::node_ptr& root_live,                                   // in/out
			      //                             const Points::node_ptr& root_dead,                             // in
			      live_clusters_t& live_clusters,
			      const Cluster::const_vector& dead_clusters,
                              cluster_length_map_t& cluster_length_map,  // in/out
                              const_cluster_set_t& cluster_connected_dead,            // in/out
                              const TPCParams& tp,                                           // common params
                              const int dead_live_overlap_offset                             // specific params
    );

    //helper function ..
    double Find_Closest_Points(const Cluster::const_pointer cluster1,
			       const Cluster::const_pointer cluster2,
			       double length_1,
			       double length_2,
			       double length_cut,
			       Blob::const_pointer mcell1,
			       Blob::const_pointer mcell2,
			       geo_point_t& p1_save,
			       geo_point_t& p2_save
			       );
			       
    
    // second function ...
    void clustering_extend(Points::node_ptr& root_live,                                   // in/out
			   live_clusters_t& live_clusters,
                           cluster_length_map_t& cluster_length_map,  // in/out
                           const_cluster_set_t& cluster_connected_dead,            // in/out
                           const TPCParams& tp,                                           // common params
                           const int flag,                                                //
                           const double length_cut = 150*units::cm,                       //
                           const int num_try = 0,                                         //
                           const double length_2_cut = 3*units::cm,                       //
                           const int num_dead_try =3                                      //
			   );
    bool Clustering_4th_prol(const Cluster::const_pointer cluster1,
			     const Cluster::const_pointer cluster2,
			     const TPCParams& tp,                                           // common params
			     double length_2,
			     geo_point_t& earliest_p,
			     geo_point_t& dir_earlp,
			     double length_cut);
    
    bool Clustering_4th_para(const Cluster::const_pointer cluster1,
			     const Cluster::const_pointer cluster2,
			     const TPCParams& tp,                                           // common params
			     double length_1, double length_2,
			     geo_point_t& earliest_p,
			     geo_point_t& dir_earlp,
			     double length_cut);
    bool Clustering_4th_reg(const Cluster::const_pointer cluster1,
			    const Cluster::const_pointer cluster2,
			    const TPCParams& tp,                                           // common params
			    double length_1, double length_2,
			    geo_point_t p1, double length_cut);
    bool Clustering_4th_dead(const Cluster::const_pointer cluster1,
			     const Cluster::const_pointer cluster2,
			     const TPCParams& tp,                                           // common params
			     double length_1, double length_2, double length_cut, int num_dead_try=3);
      

    // third function 
    void clustering_regular(Points::node_ptr& root_live,                                   // in/out
			    live_clusters_t& live_clusters,
                            cluster_length_map_t& cluster_length_map,  // in/out
                            const_cluster_set_t& cluster_connected_dead,            // in/out
                            const TPCParams& tp,                                           // common params
                            const double length_cut = 45*units::cm,                                       //
                            bool flag_enable_extend = true                                       //
    );
    bool Clustering_1st_round(const Cluster::const_pointer cluster1,
			      const Cluster::const_pointer cluster2,
			      const TPCParams& tp,                                           // common params
			      double length_1,
			      double length_2,
			      double length_cut = 45*units::cm,
			      bool flag_enable_extend = true);

    //
    void clustering_parallel_prolong(Points::node_ptr& root_live,                                   // in/out
				     live_clusters_t& live_clusters,
                                     cluster_length_map_t& cluster_length_map,  // in/out
                                     const_cluster_set_t& cluster_connected_dead,            // in/out
                                     const TPCParams& tp,                                           // common params
                                     const double length_cut = 35*units::cm                                       //
    );
    bool Clustering_2nd_round(const Cluster::const_pointer cluster1,
			      const Cluster::const_pointer cluster2,
			      const TPCParams& tp,                                           // common params
			      double length_1,
			      double length_2,
			      double length_cut = 35*units::cm);
    
    
    //
    void clustering_close(Points::node_ptr& root_live,                                   // in/out
			  live_clusters_t& live_clusters,
                          cluster_length_map_t& cluster_length_map,  // in/out
                          const_cluster_set_t& cluster_connected_dead,            // in/out
                          const TPCParams& tp,                                           // common params
                          const double length_cut = 1*units::cm                          //
    );
    bool Clustering_3rd_round( const Cluster::const_pointer cluster1,
			       const Cluster::const_pointer cluster2,
			       double length_1,
			       double length_2,
			       double length_cut = 1*units::cm);

    
}  // namespace WireCell::PointCloud::Facade
