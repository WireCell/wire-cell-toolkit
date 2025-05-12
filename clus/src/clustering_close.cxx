#include "WireCellClus/IClusteringMethod.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/ExecMon.h"

class ClusteringClose;
WIRECELL_FACTORY(ClusteringClose, ClusteringClose,
                 WireCell::IConfigurable, WireCell::Clus::IClusteringMethod)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

static void clustering_close(Grouping& live_clusters,           // 

                             const Tree::Scope& scope,
                             const double length_cut = 1*units::cm //
  );

class ClusteringClose : public IConfigurable, public Clus::IClusteringMethod, private NeedScope {
public:
  ClusteringClose() {}
  virtual ~ClusteringClose() {}

  void configure(const WireCell::Configuration& config) {
    NeedScope::configure(config);
    
    length_cut_ = get(config, "length_cut", 1*units::cm);
  }
  virtual Configuration default_configuration() const {
    Configuration cfg;
    return cfg;
  }

  void clustering(Grouping& live_clusters, Grouping& ) const {
    clustering_close(live_clusters, m_scope, length_cut_);
  }
  
private:
  double length_cut_{1*units::cm};
};



#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

static bool Clustering_3rd_round(
  const Cluster& cluster1,
  const Cluster& cluster2,
  double length_1,
  double length_2,
  double length_cut)
{
  geo_point_t p1;
  geo_point_t p2;

  bool flag_print = false;
  ExecMon em("starting");
  
  double dis = WireCell::Clus::Facade::Find_Closest_Points(cluster1, cluster2,
                                                                 length_1, length_2,
                                                                 length_cut, p1,p2);

  //  if (flag_print) std::cout << em("Find Closest Points") << std::endl;

  geo_point_t dir1, dir2;
  int num_p1{0}, num_p2{0}, num_tp1{0}, num_tp2{0};

  // if very close merge anyway???
  if (dis < 0.5*units::cm){
    return true;
  }

  if (dis < 1.0*units::cm && length_2 < 12*units::cm && length_1 <12*units::cm) {
    return true;
  }

  if (dis < 2.0*units::cm && (length_2 >=12*units::cm || length_1 >=12*units::cm)){
    dir1 = cluster1.vhough_transform(p1,50*units::cm); // cluster 1 direction based on hough
    dir2 = cluster2.vhough_transform(p2,50*units::cm); // cluster 1 direction based on hough

    if (flag_print) std::cout << em("Hough Transform") << std::endl;
    
    std::pair<int,int> num_ps_1 = cluster1.ndipole(p1,dir1);
    std::pair<int,int> num_ps_2 = cluster2.ndipole(p2,dir2);

    num_p1 = cluster1.nnearby(p1, 10*units::cm);
    num_p2 = cluster2.nnearby(p2, 10*units::cm);
    num_tp1 = cluster1.npoints();
    num_tp2 = cluster2.npoints();

    if (flag_print) std::cout << em("Get Number Points") << std::endl;
    
    if (length_1 > 25*units::cm && length_2 > 25*units::cm){
      /* if (length_1 > 60*units::cm && length_2 > 60*units::cm ){ */
      /* 	//if (num_ps_1.second > num_ps_1.first * 0.05 || num_ps_2.second > num_ps_2.first * 0.05) */
      /* 	return false; */
      /* } */
      
      if ((num_ps_1.second < num_ps_1.first * 0.02 || num_ps_1.second <=3) &&
	  (num_ps_2.second < num_ps_2.first * 0.02 || num_ps_2.second <=3) ||
	  (num_ps_1.second < num_ps_1.first * 0.035 || num_ps_1.second <=6) &&
	  (num_ps_1.second <=1 || num_ps_1.second < num_ps_1.first * 0.005) ||
	  (num_ps_1.second <=1 || num_ps_2.second < num_ps_2.first * 0.005) &&
	  (num_ps_2.second < num_ps_2.first * 0.035 || num_ps_2.second <=6)
	  )
	return true;
    }
  }

  
  if (dis < length_cut && (length_2 >=12*units::cm || length_1 >=12*units::cm)){
    geo_point_t cluster1_ave_pos = cluster1.calc_ave_pos(p1,10*units::cm);
    geo_point_t cluster2_ave_pos = cluster2.calc_ave_pos(p2,10*units::cm);
    
    geo_point_t tempV1(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
    geo_point_t tempV2(cluster2_ave_pos.x() - cluster1_ave_pos.x(), cluster2_ave_pos.y() - cluster1_ave_pos.y(), cluster2_ave_pos.z() - cluster1_ave_pos.z());
    
    // one small the other one is big 
    if (length_1 < 12 *units::cm && num_p1 > 0.5*num_tp1 && (num_p2> 50 || num_p2 > 0.25*num_tp2) ||
	length_2 < 12*units::cm && num_p2 > 0.5*num_tp2 && (num_p1>50 || num_p1 > 0.25*num_tp1) )
      return true;
    
    if (length_1 < 12*units::cm && num_p1 < 0.5*num_tp1) return false;
    if (length_2 < 12*units::cm && num_p2 < 0.5*num_tp2) return false;
    
    if ((num_p1 > 25 || num_p1 > 0.25*num_tp1 ) && (num_p2 > 25 || num_p2 > 0.25*num_tp2)){
      double angle5 = tempV1.angle(tempV2);
              
      if (length_1 < 60*units::cm || length_2 < 60*units::cm){
        if (angle5 < 30/180.*3.1415926)
          return true;
        if (angle5 < 90/180.*3.1415926 && (num_p1 > 50 && num_p2 > 50) && (num_p1>75 || num_p2>75))
          return true;
      }
      
      if ((length_1 < 60*units::cm || num_p1 >40) && (length_2 < 60*units::cm || num_p2 > 40)){
	
        if ((3.1415926 - dir1.angle(dir2))/3.1415926*180 < 30 &&
            (3.1415926 - dir1.angle(tempV1))/3.1415926*180. < 60 &&
            dir2.angle(tempV1)/3.1415926*180.<60 ||
            (3.1415926 - dir1.angle(dir2))/3.1415926*180 < 15)
          return true;

        geo_point_t dir3 = cluster1.vhough_transform(cluster1_ave_pos,50*units::cm); // cluster 1 direction based on hough
        geo_point_t dir4 = cluster2.vhough_transform(cluster2_ave_pos,50*units::cm); // cluster 1 direction based on hough

        if ((3.1415926 - dir3.angle(dir4))/3.1415926*180 < 25 &&
            (3.1415926 - dir3.angle(tempV2))/3.1415926*180. < 15 &&
            dir4.angle(tempV2)/3.1415926*180.<15 ||
            (3.1415926 - dir3.angle(dir4))/3.1415926*180 < 15)
          return true;

        if (dis<0.6*units::cm && ((3.1415926 - dir3.angle(tempV2))/3.1415926*180. < 45 && dir4.angle(tempV2)/3.1415926*180. < 90 || (3.1415926 - dir3.angle(tempV2))/3.1415926*180. < 90 && dir4.angle(tempV2)/3.1415926*180. < 45))
          return true;
	
      }
    }
    //    if (flag_print) std::cout << em("additional running") << std::endl;
  }

  return false;
}


// This function can handle multiple APA/Faces
static void clustering_close(
    Grouping& live_grouping,

    const Tree::Scope& scope,
    const double length_cut)
{

  cluster_set_t used_clusters;
  
  // prepare graph ...
  typedef cluster_connectivity_graph_t Graph;
  Graph g;
  std::unordered_map<int, int> ilive2desc;  // added live index to graph descriptor
  std::map<const Cluster*, int> map_cluster_index;
  const auto& live_clusters = live_grouping.children();

  for (size_t ilive = 0; ilive < live_clusters.size(); ++ilive) {
    const auto& live = live_clusters.at(ilive);
    if (live->get_default_scope().hash() != scope.hash()) {
      live->set_default_scope(scope);
    }
    map_cluster_index[live] = ilive;
    ilive2desc[ilive] = boost::add_vertex(ilive, g);
  }

  for (size_t i=0;i!=live_clusters.size();i++){
    auto cluster_1 = live_clusters.at(i);
    // nor process this cluster if it is not in the filter ...
    if (!cluster_1->get_scope_filter(scope)) continue;
    if (cluster_1->get_length() < 1.5*units::cm) continue;
    if (used_clusters.find(cluster_1)!=used_clusters.end()) continue;
    for (size_t j=i+1;j<live_clusters.size();j++){
      auto cluster_2 = live_clusters.at(j);
      // nor process this cluster if it is not in the filter ...
      if (!cluster_2->get_scope_filter(scope)) continue;
      if (used_clusters.find(cluster_2)!=used_clusters.end()) continue;
      if (cluster_2->get_length() < 1.5*units::cm) continue;
      if (Clustering_3rd_round(*cluster_1,*cluster_2,
                               cluster_1->get_length(), cluster_2->get_length(), length_cut)){
        //to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
        boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
            ilive2desc[map_cluster_index[cluster_2]], g);


        
        if (cluster_1->get_length() < 5*units::cm){
          used_clusters.insert(cluster_1);
          break;
        }
        if (cluster_2->get_length() < 5*units::cm){
          used_clusters.insert(cluster_2);
        }
      }
    }
  }

  //  if (flag_print) std::cout << em("core alg") << std::endl;
  
  // new function to  merge clusters ...
  merge_clusters(g, live_grouping);

  //  if (flag_print) std::cout << em("merge clusters") << std::endl;

  //  {
  //  auto live_clusters = live_grouping.children(); // copy
  //   // Process each cluster
  //   for (size_t iclus = 0; iclus < live_clusters.size(); ++iclus) {
  //       Cluster* cluster = live_clusters.at(iclus);
  //       auto& scope = cluster->get_default_scope();
  //       std::cout << "Test: " << iclus << " " << cluster->nchildren() << " " << scope.pcname << " " << scope.coords[0] << " " << scope.coords[1] << " " << scope.coords[2] << " " << cluster->get_scope_filter(scope)<< " " << cluster->get_center() << std::endl;
  //   }
  // }

  // set cluster id ... 
  int cluster_id = 1;
  for (auto* cluster : live_grouping.children()) {
      cluster->set_cluster_id(cluster_id++);
  }
}




#pragma GCC diagnostic pop

// Local Variables:
// mode: c++
// c-basic-offset: 2
// End:
