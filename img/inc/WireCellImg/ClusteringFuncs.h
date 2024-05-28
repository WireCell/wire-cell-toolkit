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

    using cluster_set_t = std::set<const Cluster*>;

    using cluster_connectivity_graph_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, int>;

    using cluster_vector_t = std::vector<Cluster*>;


    // clustering_util.cxx
    // merging clustering function
    void merge_clusters(cluster_connectivity_graph_t& g, // 
			Grouping& live_clusters,
			cluster_set_t& cluster_connected_dead); // in/out
    
    
    // clustering_live_dead.cxx
    // first function ...
    void clustering_live_dead(Grouping& live_clusters,
			      const Grouping& dead_clusters,
                              cluster_set_t& cluster_connected_dead, // in/out
                              const int dead_live_overlap_offset // specific params
    );

    // clustering_extend.cxx
    //helper function ..
    double Find_Closest_Points(const Cluster& cluster1,
			       const Cluster& cluster2,
			       double length_1,
			       double length_2,
			       double length_cut,
			       geo_point_t& p1_save, // output
			       geo_point_t& p2_save  // output
			       );
			       
    
    // clustering_extend.cxx
    // second function ...
    void clustering_extend(Grouping& live_clusters,
                           cluster_set_t& cluster_connected_dead,            // in/out
                           const int flag,                                                //
                           const double length_cut = 150*units::cm,                       //
                           const int num_try = 0,                                         //
                           const double length_2_cut = 3*units::cm,                       //
                           const int num_dead_try =3                                      //
			   );

    bool Clustering_4th_prol(const Cluster& cluster1,
			     const Cluster& cluster2,
			     double length_2,
			     geo_point_t& earliest_p,
			     geo_point_t& dir_earlp,
			     double length_cut);
    
    bool Clustering_4th_para(const Cluster& cluster1,
			     const Cluster& cluster2,
			     double length_1, double length_2,
			     geo_point_t& earliest_p,
			     geo_point_t& dir_earlp,
			     double length_cut);

    bool Clustering_4th_reg(const Cluster& cluster1,
			    const Cluster& cluster2,
			    double length_1, double length_2,
			    geo_point_t p1, double length_cut);

    bool Clustering_4th_dead(const Cluster& cluster1,
			     const Cluster& cluster2,
			     double length_1, double length_2, double length_cut, int num_dead_try=3);
      

    // clustering_regular.cxx
    // third function 
    void clustering_regular(Grouping& live_clusters,
                            cluster_set_t& cluster_connected_dead,            // in/out
                            const double length_cut = 45*units::cm,                                       //
                            bool flag_enable_extend = true                                       //
    );

    bool Clustering_1st_round(const Cluster& cluster1,
			      const Cluster& cluster2,
			      double length_1,
			      double length_2,
			      double length_cut = 45*units::cm,
			      bool flag_enable_extend = true);

    // clustering_parallel_prolong.cxx:
    void clustering_parallel_prolong(Grouping& live_clusters,
                                     cluster_set_t& cluster_connected_dead, // in/out
                                     const double length_cut = 35*units::cm
    );

    bool Clustering_2nd_round(const Cluster& cluster1,
			      const Cluster& cluster2,
			      double length_1,
			      double length_2,
			      double length_cut = 35*units::cm);
    
    // clustering_close.cxx
    void clustering_close(Grouping& live_clusters,           // 
                          cluster_set_t& cluster_connected_dead, // in/out
                          const double length_cut = 1*units::cm //
    );

    bool Clustering_3rd_round( const Cluster& cluster1,
			       const Cluster& cluster2,
			       double length_1,
			       double length_2,
			       double length_cut = 1*units::cm);

    
}  // namespace WireCell::PointCloud::Facade
