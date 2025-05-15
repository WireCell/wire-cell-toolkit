#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"

class ClusteringParallelProlong;
WIRECELL_FACTORY(ClusteringParallelProlong, ClusteringParallelProlong,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

static void clustering_parallel_prolong(
    Grouping& live_clusters,

    IDetectorVolumes::pointer dv,
    const Tree::Scope& scope,
    const double length_cut = 35*units::cm);

class ClusteringParallelProlong : public IConfigurable, public Clus::IEnsembleVisitor, private NeedDV, private NeedScope {
public:
    ClusteringParallelProlong() {}
    virtual ~ClusteringParallelProlong() {}

    void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedScope::configure(config);
        
        length_cut_ = get(config, "length_cut", 35*units::cm);
    }
    virtual Configuration default_configuration() const {
        Configuration cfg;
        return cfg;
    }
    

    void visit(Ensemble& ensemble) const {
        auto& live = *ensemble.with_name("live").at(0);
        clustering_parallel_prolong(live, m_dv, m_scope, length_cut_);
    }

private:
    double length_cut_{35*units::cm};
};


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

// using namespace WireCell::PointCloud::Tree;

static bool Clustering_2nd_round(
    const Cluster& cluster1,
    const Cluster& cluster2,
    double length_1,
    double length_2,
	const std::map<WirePlaneId, std::pair<geo_point_t, double> > & wpid_U_dir, const std::map<WirePlaneId, std::pair<geo_point_t, double> > & wpid_V_dir, const std::map<WirePlaneId, std::pair<geo_point_t, double> > & wpid_W_dir,
	const IDetectorVolumes::pointer dv,
    double length_cut)
{
//   const auto [angle_u,angle_v,angle_w] = cluster1.grouping()->wire_angles();

  if (length_1 < 10*units::cm && length_2 < 10*units::cm) return false;

  geo_point_t p1;
  geo_point_t p2;

  double dis = WireCell::Clus::Facade::Find_Closest_Points(cluster1, cluster2,
                                                                 length_1, length_2,
                                                                 length_cut, p1, p2);

	auto wpid_p1 = cluster1.wpid(p1);
	auto wpid_p2 = cluster2.wpid(p2);
	auto wpid_ps = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv);


  if ((dis < length_cut || (dis < 80*units::cm && length_1 +length_2 > 50*units::cm && length_1>15*units::cm && length_2 > 15*units::cm))){
    geo_point_t cluster1_ave_pos = cluster1.calc_ave_pos(p1,10*units::cm);
	auto wpid_ave_p1 = cluster1.wpid(cluster1_ave_pos);
    geo_point_t cluster2_ave_pos = cluster2.calc_ave_pos(p2,10*units::cm);
	auto wpid_ave_p2 = cluster2.wpid(cluster2_ave_pos);
    // auto wpid_ave_ps = get_wireplaneid(cluster1_ave_pos, wpid_ave_p1, cluster2_ave_pos, wpid_ave_p2, dv);

    bool flag_para = false;

    geo_point_t drift_dir_abs(1, 0, 0);  // assuming the drift direction is along X ...
    
    
    // deal the parallel case ...
    if (length_1 > 10*units::cm && length_2 >10*units::cm){
      geo_point_t tempV1(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
      geo_point_t tempV2(cluster2_ave_pos.x() - cluster1_ave_pos.x(), cluster2_ave_pos.y() - cluster1_ave_pos.y(), cluster2_ave_pos.z() - cluster1_ave_pos.z());
      
      double angle1 = tempV1.angle(drift_dir_abs);
      double angle4 = tempV2.angle(drift_dir_abs);


      // looks like a parallel case
      if ( (fabs(angle1-3.1415926/2.)<10/180.*3.1415926 && dis > 10*units::cm ||
	    fabs(angle1-3.1415926/2.)<20/180.*3.1415926 && dis > 3*units::cm && dis <= 10*units::cm ||
	    fabs(angle1-3.1415926/2.)<45/180.*3.1415926 && dis <=3*units::cm)
	   && fabs(angle4-3.1415926/2.)<5/180.*3.1415926){
	
	geo_point_t dir1 = cluster1.vhough_transform(p1,60*units::cm); // cluster 1 direction based on hough
	geo_point_t dir2 = cluster2.vhough_transform(p2,60*units::cm); // cluster 2 direction based on hough
	
	double angle5 = dir1.angle(drift_dir_abs);
	double angle6 = dir2.angle(drift_dir_abs);
	
	if (fabs(angle5-3.1415926/2.)<5/180.*3.1415926 && fabs(angle6-3.1415926/2.)<20/180.*3.1415926 ||
	    fabs(angle5-3.1415926/2.)<20/180.*3.1415926 && fabs(angle6-3.1415926/2.)<5/180.*3.1415926){
	
	  flag_para = true;


	  
	  if (dis >= 3*length_1 && dis >= 3*length_2 && flag_para) return false;
	  
	  double angle2 = tempV1.angle(wpid_U_dir.at(wpid_ps).first);
	  double angle3 = tempV1.angle(wpid_V_dir.at(wpid_ps).first);

	  

	  // look at parallel U
	  if ((fabs(angle2-3.1415926/2.)<7.5/180.*3.1415926 || (fabs(angle2-3.1415926/2.)<15/180.*3.1415926)&&dis <6*units::cm) && (dis<length_cut || (length_1 + length_2 > 100*units::cm)) && length_1 >15*units::cm && length_2 > 15*units::cm){
	    // flag_para_U = true;



	    if ((length_1 < 25*units::cm || length_2 < 25*units::cm) && fabs(angle2-3.1415926/2.)<5.0/180.*3.1415926  && dis < 15* units::cm || dis < 3*units::cm){
	      // for short or small distance one
	      return true;
	    }else if (fabs(angle6-3.1415926/2.)/3.1415926*180.<1 && fabs(angle5-3.1415926/2.)/3.1415926*180.<1 && fabs(angle2-3.1415926/2.)<5.0/180.*3.1415926 && dis < 20*units::cm && (length_1 < 50*units::cm || length_2 < 50*units::cm)){
	      return true;
	    }else if (dis < 15*units::cm && (length_1 < 60*units::cm || length_2 < 60*units::cm) &&
		      fabs(angle2-3.1415926/2.)<2.5/180.*3.1415926){
	      // parallel case for reasonably short distance one
	      return true;
	    }else if (fabs(angle2-3.1415926/2.)<2.5/180.*3.1415926 && fabs(angle5-3.1415926/2.)<5/180.*3.1415926 && fabs(angle6-3.1415926/2.)<5/180.*3.141592 ){
	      // parallel case, but exclude both very long tracks
	      if (length_1 < 60*units::cm || length_2 < 60*units::cm){
			if (WireCell::Clus::Facade::is_angle_consistent(dir1,tempV1,false,15, wpid_U_dir.at(wpid_ps).second, wpid_V_dir.at(wpid_ps).second, wpid_W_dir.at(wpid_ps).second) &&     WireCell::Clus::Facade::is_angle_consistent(dir2,tempV1,true,15,wpid_U_dir.at(wpid_ps).second, wpid_V_dir.at(wpid_ps).second, wpid_W_dir.at(wpid_ps).second)) return true;
			}else if (dis <5*units::cm){
				return true;
			}else{
				double angle7 = (3.1415926-dir1.angle(dir2))/3.1415926*180.;
				double angle8 = (3.1415926-dir1.angle(tempV1))/3.1415926*180.; // dir1 = -p1, tempV1 = p2 - p1
				double angle9 = dir2.angle(tempV1)/3.1415926*180.; // dir2 = -p2
				if (angle7 < 30 && angle8 < 30 && angle9 < 30) return true;
			if (WireCell::Clus::Facade::is_angle_consistent(dir1,tempV1,false,10,wpid_U_dir.at(wpid_ps).second, wpid_V_dir.at(wpid_ps).second, wpid_W_dir.at(wpid_ps).second) && WireCell::Clus::Facade::is_angle_consistent(dir2,tempV1,true,10,wpid_U_dir.at(wpid_ps).second, wpid_V_dir.at(wpid_ps).second, wpid_W_dir.at(wpid_ps).second)) 
				return true; 
			}
	    }else{
	      // general case ... (not sure how useful though ...)
	      double angle7 = (3.1415926-dir1.angle(dir2))/3.1415926*180.;
	      double angle8 = (3.1415926-dir1.angle(tempV1))/3.1415926*180.; // dir1 = -p1, tempV1 = p2 - p1
	      double angle9 = dir2.angle(tempV1)/3.1415926*180.; // dir2 = -p2
	      if ((angle7 < 30 && angle8 < 30 && angle9 < 30 ||
		   fabs(angle5-3.1415926/2.)<5/180.*3.1415926 && fabs(angle6-3.1415926/2.)<5/180.*3.141592 &&
		   angle7 < 45 && angle8 < 45 && angle9 < 45) && dis < 20*units::cm)
			return true;
	      if (WireCell::Clus::Facade::is_angle_consistent(dir1,tempV1,false,10,wpid_U_dir.at(wpid_ps).second, wpid_V_dir.at(wpid_ps).second, wpid_W_dir.at(wpid_ps).second) && WireCell::Clus::Facade::is_angle_consistent(dir2,tempV1,true,10,wpid_U_dir.at(wpid_ps).second, wpid_V_dir.at(wpid_ps).second, wpid_W_dir.at(wpid_ps).second)) 
			return true; 
	    }
	  }

	  // look at parallel V
	  if ((fabs(angle3-3.1415926/2.)<7.5/180.*3.1415926 || (fabs(angle3-3.1415926/2.)<15/180.*3.1415926)&&dis <6*units::cm )&&(dis<length_cut || (length_1 + length_2 > 100*units::cm))&& length_1 >15*units::cm && length_2 > 15*units::cm){
	    // flag_para_V = true;
	    //return true;
	    
	    if ((length_1 < 25*units::cm || length_2 < 25*units::cm) && fabs(angle3-3.1415926/2.)<5.0/180.*3.1415926 && dis < 15* units::cm || dis < 2*units::cm){
	      return true;
	    }else if (fabs(angle6-3.1415926/2.)/3.1415926*180.<1 && fabs(angle5-3.1415926/2.)/3.1415926*180.<1 && fabs(angle3-3.1415926/2.)<5.0/180.*3.1415926 && dis < 20*units::cm && (length_1 < 50*units::cm || length_2 < 50*units::cm)){
	      return true;
	    }else if (dis < 15*units::cm && fabs(angle3-3.1415926/2.)<2.5/180.*3.1415926 && (length_1 < 60*units::cm || length_2 < 60*units::cm) ){
	      return true;
	    }else if (fabs(angle3-3.1415926/2.)<2.5/180.*3.1415926 && fabs(angle5-3.1415926/2.)<5/180.*3.1415926 && fabs(angle6-3.1415926/2.)<5/180.*3.141592){
	      if (WireCell::Clus::Facade::is_angle_consistent(dir1,tempV1,false,15,wpid_U_dir.at(wpid_ps).second, wpid_V_dir.at(wpid_ps).second, wpid_W_dir.at(wpid_ps).second) && WireCell::Clus::Facade::is_angle_consistent(dir2,tempV1,true,15,wpid_U_dir.at(wpid_ps).second, wpid_V_dir.at(wpid_ps).second, wpid_W_dir.at(wpid_ps).second)) 
		return true;
	    }else{
	      double angle7 = (3.1415926-dir1.angle(dir2))/3.1415926*180.;
	      double angle8 = (3.1415926-dir1.angle(tempV1))/3.1415926*180.; // dir1 = -p1, tempV1 = p2 - p1
	      double angle9 = dir2.angle(tempV1)/3.1415926*180.; // dir2 = -p2
	      if (angle7 < 30 && angle8 < 30 && angle9 < 30||
		  fabs(angle5-3.1415926/2.)<5/180.*3.1415926 && fabs(angle6-3.1415926/2.)<5/180.*3.141592 &&
		  angle7 < 60 && angle8 < 60 && angle9 < 60)
		return true;
	      if (WireCell::Clus::Facade::is_angle_consistent(dir1,tempV1,false,10,wpid_U_dir.at(wpid_ps).second, wpid_V_dir.at(wpid_ps).second, wpid_W_dir.at(wpid_ps).second) && WireCell::Clus::Facade::is_angle_consistent(dir2,tempV1,true,10,wpid_U_dir.at(wpid_ps).second, wpid_V_dir.at(wpid_ps).second, wpid_W_dir.at(wpid_ps).second))
	      	  return true;
	    }
	  }
	}
      }
    }

     // look at prolonged case ... (add W case) 
    {
      geo_point_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
      geo_point_t tempV5;
      double angle1 = tempV1.angle(wpid_U_dir.at(wpid_ps).first);
      tempV5.set(fabs(p2.x()-p1.x()),sqrt(pow(p2.y() - p1.y(),2)+pow(p2.z() - p1.z(),2))*sin(angle1),0);
      angle1 = tempV5.angle(drift_dir_abs);
      
      double angle2 = tempV1.angle(wpid_V_dir.at(wpid_ps).first);
      tempV5.set(fabs(p2.x()-p1.x()),sqrt(pow(p2.y() - p1.y(),2)+pow(p2.z() - p1.z(),2))*sin(angle2),0);
      angle2 = tempV5.angle(drift_dir_abs);

      double angle1p = tempV1.angle(wpid_W_dir.at(wpid_ps).first);
      tempV5.set(fabs(p2.x()-p1.x()),sqrt(pow(p2.y() - p1.y(),2)+pow(p2.z() - p1.z(),2))*sin(angle1p),0);
      angle1p = tempV5.angle(drift_dir_abs);

      if (angle1<7.5/180.*3.1415926  ||
	  angle2<7.5/180.*3.1415926  ||
	  angle1p<7.5/180.*3.1415926 ){
	if (length_1 > 10*units::cm || length_2 > 10*units::cm){
	  geo_point_t dir1 = cluster1.vhough_transform(p1,60*units::cm); // cluster 1 direction based on hough
	  geo_point_t dir2 = cluster2.vhough_transform(p2,60*units::cm); // cluster 1 direction based on hough
	  geo_point_t dir3(p2.x()-p1.x(),p2.y()-p1.y(),p2.z()-p1.z());
	  double angle3 = dir3.angle(dir2);
	  double angle4 = 3.1415926-dir3.angle(dir1);

	  if ((angle3<25/180.*3.1415926 || length_2<10*units::cm)&&(angle4<25/180.*3.1415926|| length_1<10*units::cm)&&dis<5*units::cm ||
	      (angle3<15/180.*3.1415926 || length_2<10*units::cm)&&(angle4<15/180.*3.1415926|| length_1<10*units::cm)&&dis<15*units::cm ||
	      (angle3<7.5/180.*3.1415926 || length_2<10*units::cm)&&(angle4<7.5/180.*3.1415926|| length_1<10*units::cm) ||
	      (angle3+angle4 < 15/180.*3.1415926 && angle3 < 10/180.*3.1415926 && angle4 < 10/180.*3.1415926)
	      )
	    return true;
	}
      }else{
      	//regular cases (only for very short distance ... )
      	if (dis < 5*units::cm){
      	  if (length_1 > 10*units::cm && length_2 >10*units::cm){
      	    geo_point_t dir1 = cluster1.vhough_transform(p1,30*units::cm); // cluster 1 direction based on hough
      	    geo_point_t dir2 = cluster2.vhough_transform(p2,30*units::cm); // cluster 1 direction based on hough
      	    geo_point_t dir3(p2.x()-p1.x(),p2.y()-p1.y(),p2.z()-p1.z());
      	    double angle3 = dir3.angle(dir2);
      	    double angle4 = 3.1415926-dir3.angle(dir1);

      	    //std::cout << angle3/3.1415926*180. << " " << angle4/3.1415926*180. << std::endl;
      	    if ((angle3<15/180.*3.1415926 || length_2<6*units::cm)
      		&& (angle4<15/180.*3.1415926|| length_1<6*units::cm))
      	      return true;
      	  }
      	}
      }
    }
  }
  return false;
}

// Expand this function to handle multiple APA/Faces ...
static void clustering_parallel_prolong(
    Grouping& live_grouping,

    const IDetectorVolumes::pointer dv,                // detector volumes
    const Tree::Scope& scope,
    const double length_cut                    //
)
{

	// Get all the wire plane IDs from the grouping
	const auto& wpids = live_grouping.wpids();

	// Key: pair<APA, face>, Value: drift_dir, angle_u, angle_v, angle_w
	std::map<WirePlaneId , std::tuple<geo_point_t, double, double, double>> wpid_params;
	std::map<WirePlaneId, std::pair<geo_point_t, double> > wpid_U_dir;
	std::map<WirePlaneId, std::pair<geo_point_t, double> > wpid_V_dir;
	std::map<WirePlaneId, std::pair<geo_point_t, double> > wpid_W_dir;
	std::set<int> apas;
	for (const auto& wpid : wpids) {
		int apa = wpid.apa();
		int face = wpid.face();
		apas.insert(apa);
  
		// Create wpids for all three planes with this APA and face
		WirePlaneId wpid_u(kUlayer, face, apa);
		WirePlaneId wpid_v(kVlayer, face, apa);
		WirePlaneId wpid_w(kWlayer, face, apa);
	 
		// Get drift direction based on face orientation
		int face_dirx = dv->face_dirx(wpid_u);
		geo_point_t drift_dir(face_dirx, 0, 0);
		
		// Get wire directions for all planes
		Vector wire_dir_u = dv->wire_direction(wpid_u);
		Vector wire_dir_v = dv->wire_direction(wpid_v);
		Vector wire_dir_w = dv->wire_direction(wpid_w);
  
		// Calculate angles
		double angle_u = std::atan2(wire_dir_u.z(), wire_dir_u.y());
		double angle_v = std::atan2(wire_dir_v.z(), wire_dir_v.y());
		double angle_w = std::atan2(wire_dir_w.z(), wire_dir_w.y());
  
		wpid_params[wpid] = std::make_tuple(drift_dir, angle_u, angle_v, angle_w);
		wpid_U_dir[wpid] = std::make_pair(geo_point_t(0, cos(angle_u), sin(angle_u)), angle_u);
		wpid_V_dir[wpid] = std::make_pair(geo_point_t(0, cos(angle_v), sin(angle_v)), angle_v);
		wpid_W_dir[wpid] = std::make_pair(geo_point_t(0, cos(angle_w), sin(angle_w)), angle_w);
	}
  


  // prepare graph ...
  typedef cluster_connectivity_graph_t Graph;
  Graph g;
  std::unordered_map<int, int> ilive2desc;  // added live index to graph descriptor
  std::map<const Cluster*, int> map_cluster_index;
  const auto& live_clusters = live_grouping.children();
  
  for (size_t ilive = 0; ilive < live_clusters.size(); ++ilive) {
    auto& live = live_clusters[ilive];
    map_cluster_index[live] = ilive;
    ilive2desc[ilive] = boost::add_vertex(ilive, g);
	if (live->get_default_scope().hash() != scope.hash()) {
		live->set_default_scope(scope);
		// std::cout << "Test: Set default scope: " << pc_name << " " << coords[0] << " " << coords[1] << " " << coords[2] << " " << cluster->get_default_scope().hash() << " " << scope.hash() << std::endl;
	}
  }

  // original algorithm ... (establish edges ... )


  for (size_t i=0;i!=live_clusters.size();i++){
    auto cluster_1 = live_clusters.at(i);
	if (!cluster_1->get_scope_filter(scope)) continue;
    for (size_t j=i+1;j<live_clusters.size();j++){
      auto cluster_2 = live_clusters.at(j);
	  if(!cluster_2->get_scope_filter(scope)) continue;
      if (Clustering_2nd_round(*cluster_1,*cluster_2, cluster_1->get_length(), cluster_2->get_length(), wpid_U_dir, wpid_V_dir, wpid_W_dir, dv, length_cut)){
		boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
			ilive2desc[map_cluster_index[cluster_2]], g);


      }
    }
  }

  // new function to  merge clusters ...
  merge_clusters(g, live_grouping);

 // {
  //  auto live_clusters = live_grouping.children(); // copy
  //   // Process each cluster
  //   for (size_t iclus = 0; iclus < live_clusters.size(); ++iclus) {
  //       Cluster* cluster = live_clusters.at(iclus);
  //       auto& scope = cluster->get_default_scope();
  //       std::cout << "Test: " << iclus << " " << cluster->nchildren() << " " << scope.pcname << " " << scope.coords[0] << " " << scope.coords[1] << " " << scope.coords[2] << " " << cluster->get_scope_filter(scope)<< " " << cluster->get_center() << std::endl;
  //   }
  // }






}



#pragma GCC diagnostic pop
