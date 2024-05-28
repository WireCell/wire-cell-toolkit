#include <WireCellImg/ClusteringFuncs.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;
void WireCell::PointCloud::Facade::clustering_regular(
    Grouping& live_grouping,
    cluster_set_t& cluster_connected_dead,            // in/out
    const double length_cut,                                       //
    bool flag_enable_extend                                        //
)
{
  double internal_length_cut = 10 *units::cm;
  if (flag_enable_extend) {
    internal_length_cut = 15 *units::cm;
  }

  // prepare graph ...
  typedef cluster_connectivity_graph_t Graph;
  Graph g;
  std::unordered_map<int, int> ilive2desc;  // added live index to graph descriptor
  std::map<const Cluster*, int> map_cluster_index;
  const auto& live_clusters = live_grouping.children();
  for (size_t ilive = 0; ilive < live_clusters.size(); ++ilive) {
    const auto& live = live_clusters.at(ilive);
    map_cluster_index[live] = ilive;
    ilive2desc[ilive] = boost::add_vertex(ilive, g);
  }

  // original algorithm ... (establish edges ... )


  for (size_t i=0;i!=live_clusters.size();i++){
    auto cluster_1 = live_clusters.at(i);
    if (cluster_1->get_length() < internal_length_cut) continue;
    for (size_t j=i+1;j<live_clusters.size();j++){
      auto cluster_2 = live_clusters.at(j);
      if (cluster_2->get_length() < internal_length_cut) continue;

      if (Clustering_1st_round(*cluster_1,*cluster_2, cluster_1->get_length(), cluster_2->get_length(), length_cut, flag_enable_extend)){
	//	to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
			ilive2desc[map_cluster_index[cluster_2]], g);


      }
    }
  }

  // new function to  merge clusters ...
  merge_clusters(g, live_grouping, cluster_connected_dead);
}

bool WireCell::PointCloud::Facade::Clustering_1st_round(
    const Cluster& cluster1,
    const Cluster& cluster2,
    double length_1,
    double length_2,
    double length_cut,
    bool flag_enable_extend)
{
  const auto [angle_u,angle_v,angle_w] = cluster1.grouping()->wire_angles();

  geo_point_t p1;
  geo_point_t p2;

  double dis = WireCell::PointCloud::Facade::Find_Closest_Points(cluster1, cluster2,
                                                                 length_1, length_2,
                                                                 length_cut, p1,p2);

  if (dis < length_cut){
    bool flag_para = false;
    bool flag_prolong_U = false;
    bool flag_prolong_V = false;
    bool flag_prolong_W= false;
    bool flag_para_U = false;
    bool flag_para_V = false;
    bool flag_regular = false;
    bool flag_extend = false;
    bool flag_force_extend = false;


    geo_point_t drift_dir(1, 0, 0);  // assuming the drift direction is along X ...
    
    // pronlonged case for U 3 and V 4 ...
    geo_point_t U_dir(0,cos(angle_u),sin(angle_u));
    geo_point_t V_dir(0,cos(angle_v),sin(angle_v));
    geo_point_t W_dir(0,cos(angle_w),sin(angle_w));

    // calculate average distance ... 
    geo_point_t cluster1_ave_pos = cluster1.calc_ave_pos(p1,5*units::cm);
    geo_point_t cluster2_ave_pos = cluster2.calc_ave_pos(p2,5*units::cm);

    
    geo_point_t dir2_1(p2.x() - p1.x()+1e-9, p2.y() - p1.y()+1e-9, p2.z() - p1.z()+1e-9); // 2-1
    geo_point_t dir2(cluster2_ave_pos.x() - cluster1_ave_pos.x()+1e-9,
                     cluster2_ave_pos.y() - cluster1_ave_pos.y()+1e-9,
                     cluster2_ave_pos.z() - cluster1_ave_pos.z()+1e-9); // 2-1
    dir2_1 = dir2_1/dir2_1.magnitude();
    dir2 = dir2/dir2.magnitude();


     // parallle case
    double angle1 = dir2_1.angle(drift_dir);
    double angle2 = dir2.angle(drift_dir);
    
    double angle3{0}, angle4{0};
    double angle3_1{0}, angle4_1{0};
    
    if ((fabs(angle1-3.1415926/2.)<7.5/180.*3.1415926 ||
	 fabs(angle2-3.1415926/2.)<7.5/180.*3.1415926) && dis < 45*units::cm &&
	length_1 > 12*units::cm && length_2 > 12*units::cm){
      flag_para = true;

      if (dis >=3*length_1 && dis >= 3*length_2 && flag_para) return false;
      
      angle3 = dir2_1.angle(U_dir);
      angle4 = dir2_1.angle(V_dir);
            
      angle3_1 = dir2.angle(U_dir);
      angle4_1 = dir2.angle(V_dir);
            
      if (fabs(angle3-3.1415926/2.)<7.5/180.*3.1415926 || fabs(angle3_1-3.1415926/2.)<7.5/180.*3.1415926 ||
	  ((fabs(angle3-3.1415926/2.)<15/180.*3.1415926 || fabs(angle3_1-3.1415926/2.)<15/180.*3.1415926)
	   &&dis < 6*units::cm))
	flag_para_U = true;
      
      if (fabs(angle4-3.1415926/2.)<7.5/180.*3.1415926 || fabs(angle4_1-3.1415926/2.)<7.5/180.*3.1415926 ||
	  ((fabs(angle4-3.1415926/2.)<15/180.*3.1415926 || fabs(angle4_1-3.1415926/2.)<15/180.*3.1415926)&&
	   dis < 6*units::cm))
	flag_para_V = true;
    } 

    if (!flag_para){
      // prolonged case
      geo_point_t tempV3(0, p2.y() - p1.y(), p2.z() - p1.z());
      geo_point_t tempV4(0, cluster2_ave_pos.y() - cluster1_ave_pos.y(),
                         cluster2_ave_pos.z() - cluster1_ave_pos.z());
      geo_point_t tempV5;
      
      double angle6 = tempV3.angle(U_dir);
      tempV5.set(fabs(p2.x()-p1.x()),sqrt(pow(p2.y() - p1.y(),2)+pow(p2.z() - p1.z(),2))*sin(angle6),0);
      angle6 = tempV5.angle(drift_dir);
      
      double angle7 = tempV3.angle(V_dir);
      tempV5.set(fabs(p2.x()-p1.x()),sqrt(pow(p2.y() - p1.y(),2)+pow(p2.z() - p1.z(),2))*sin(angle7),0);
      angle7 = tempV5.angle(drift_dir);
      
      double angle8 = tempV3.angle(W_dir);
      tempV5.set(fabs(p2.x()-p1.x()),sqrt(pow(p2.y() - p1.y(),2)+pow(p2.z() - p1.z(),2))*sin(angle8),0);
      angle8 = tempV5.angle(drift_dir);
      
      double angle6_1 = tempV4.angle(U_dir);
      tempV5.set(fabs(cluster2_ave_pos.x()-cluster1_ave_pos.x()),sqrt(pow(cluster2_ave_pos.y()-cluster1_ave_pos.y(),2)+pow(cluster2_ave_pos.z()-cluster1_ave_pos.z(),2))*sin(angle6_1),0);
      angle6_1 = tempV5.angle(drift_dir);
      
      double angle7_1 = tempV4.angle(V_dir);
      tempV5.set(fabs(cluster2_ave_pos.x()-cluster1_ave_pos.x()),sqrt(pow(cluster2_ave_pos.y()-cluster1_ave_pos.y(),2)+pow(cluster2_ave_pos.z()-cluster1_ave_pos.z(),2))*sin(angle7_1),0);
      angle7_1 = tempV5.angle(drift_dir);
      
      double angle8_1 = tempV4.angle(W_dir);
      tempV5.set(fabs(cluster2_ave_pos.x()-cluster1_ave_pos.x()),sqrt(pow(cluster2_ave_pos.y()-cluster1_ave_pos.y(),2)+pow(cluster2_ave_pos.z()-cluster1_ave_pos.z(),2))*sin(angle8_1),0);
      angle8_1 = tempV5.angle(drift_dir);
      
      if (angle6<15/180.*3.1415926  ||
	  angle6_1<15/180.*3.1415926 )
	flag_prolong_U = true;
      
      if (angle7<15/180.*3.1415926  ||
	  angle7_1<15/180.*3.1415926 )
	flag_prolong_V = true;
      
      if (angle8<15/180.*3.1415926  ||
	  angle8_1<15/180.*3.1415926 )
	flag_prolong_W = true;
    }

    // regular case
    if (dis <= 15*units::cm ){
      flag_regular = true;
    }else if ( length_1 > 30*units::cm && length_2 > 30*units::cm) {
      if (dis <= 25*units::cm)
	flag_regular = true;
    }
    
    if ((flag_para_U || flag_para_V )  ||
	(flag_prolong_U || flag_prolong_V || flag_prolong_W) ||
	flag_regular){
      double angle_cut=2.5;
      double para_angle_cut = 5.;
      double para_angle_cut_1 = 5;
      
      if (dis < 5*units::cm){
	angle_cut = 12;
      }else if (dis < 15*units::cm){
	angle_cut = 7.5;
      }else{
	angle_cut = 5;
      }
      
      if (dis > 45*units::cm){
	para_angle_cut = 5;
      }else if (dis > 15*units::cm){
	para_angle_cut = 15;
      }else if (dis > 5*units::cm){
	para_angle_cut = 30;
      }else{
	para_angle_cut = 60;
      }
      
      geo_point_t dir1 = cluster1.vhough_transform(cluster1_ave_pos,30*units::cm); // cluster 1 direction based on hough
      geo_point_t dir1_1 = cluster1.vhough_transform(p1,30*units::cm);
      
      geo_point_t dir3 = cluster2.vhough_transform(cluster2_ave_pos,30*units::cm);
      geo_point_t dir3_1 = cluster2.vhough_transform(p2,30*units::cm);

      dir1_1 = dir1_1/dir1_1.magnitude();
      dir1 = dir1 / dir1.magnitude();
      dir3 = dir3 / dir3.magnitude();
      dir3_1 = dir3_1/ dir3_1.magnitude();
      
      double angle_diff1 = (3.1415926-dir1.angle(dir2))/3.1415926*180.;
      double angle_diff1_1 = (3.1415926-dir1_1.angle(dir2_1))/3.1415926*180.;

      double angle_diff2 = dir3.angle(dir2)/3.1415926*180;
      double angle_diff2_1 = dir3_1.angle(dir2_1)/3.1415926*180.;

      double angle_diff3 = (3.1415926 - dir1.angle(dir3))/3.1415926*180;
      double angle_diff3_1 = (3.1415926 - dir1_1.angle(dir3_1))/3.1415926*180.;

      
      if (dis<=3*units::cm){
	if ((angle_diff3 < angle_cut*1.5 || angle_diff3_1 < angle_cut*1.5) &&
	    ((angle_diff1 < angle_cut*4.5 || angle_diff1_1 < angle_cut*4.5 || length_1 < 12*units::cm) &&
	     (angle_diff2 < angle_cut*4.5 || angle_diff2_1 < angle_cut*4.5 || length_2 < 12*units::cm) ) )
	  return true;
	if ((angle_diff3 < angle_cut*1.5 || angle_diff3_1 < angle_cut*1.5) ||
	    (angle_diff1 < angle_cut*1.5 || angle_diff1_1 < angle_cut*1.5) ||
	    (angle_diff2 < angle_cut*1.5 || angle_diff2_1 < angle_cut*1.5))
	  flag_extend = true;
      }
      
      if (angle_diff1 < angle_cut || angle_diff1_1 < angle_cut){ // possible match 
	if (length_2 < 12*units::cm && dis < 2*units::cm) //small and very close
	  return true;
	if (angle_diff3 < angle_cut || angle_diff3_1 < angle_cut ){
	  return true;
	}

	if (dis < 15*units::cm && dis> 5 *units::cm &&
	    (angle_diff3 < 2*angle_cut || angle_diff3_1 < 2*angle_cut) &&
	    (length_1 < 15*units::cm || length_2 < 15*units::cm))
	  return true;
	
	 flag_extend = true;
      }

      if (angle_diff2 < angle_cut || angle_diff2_1 < angle_cut){
	if (length_1 < 12*units::cm && dis < 2*units::cm) // small and very close
	  return true;

	if (angle_diff3 < angle_cut || angle_diff3_1 < angle_cut)
	  return true;

	if (dis < 15*units::cm && dis> 5 *units::cm &&
	    (angle_diff3 < 2*angle_cut || angle_diff3_1 < 2*angle_cut) &&
	    (length_1 < 15*units::cm || length_2 < 15*units::cm))
	  return true;
	
	flag_extend = true;
      }


      if (flag_para && (flag_para_U || flag_para_V  || flag_regular)){
	
	double dangle1 = (dir1.angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
	double dangle1_1 = (dir1_1.angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
	
	double dangle2 = (dir2.angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
	double dangle2_1 = (dir2_1.angle(drift_dir)-3.1415926/2.)/3.1415926*180.;

	double dangle3 = (dir3.angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
	double dangle3_1 = (dir3_1.angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
	
	double dangle4 = dangle1 + dangle2;
	double dangle4_1 = dangle1_1 + dangle2_1;
	
	double dangle5 = dangle3 - dangle2;
	double dangle5_1 = dangle3_1 - dangle2_1;
	
	if ((flag_para_U || flag_para_V) &&
	    (fabs(angle3-3.1415926/2.)/3.1415926*180 < 5 || fabs(angle4-3.1415926/2.)/3.1415926*180 < 5) && // along U or V very much
	    (fabs(angle1-3.1415926/2.)/3.1415926*180<5 && fabs(angle2-3.1415926/2.)/3.1415926*180<5) && (dis<35*units::cm && (length_1 <35*units::cm || length_2 < 35*units::cm)) && // parallel case
	    (angle_diff2 < 60 || angle_diff2_1 < 60 || length_1 < 12*units::cm) &&
	    (angle_diff1 < 60 || angle_diff1_1 < 60 || length_2 < 12*units::cm ) &&
	    (angle_diff3 < 90 || angle_diff3_1 < 90 || length_1<12*units::cm || length_2 < 12*units::cm) && // angle is close
	    (fabs(dangle4) < 5 && fabs(dangle4_1)<5 || fabs(dangle5)<5 && fabs(dangle5_1) <5) &&
	    (fabs(dangle1) < 10 || fabs(dangle1_1) < 10) &&
	    (fabs(dangle3) < 10 || fabs(dangle3_1) < 10) //small angle close to parallel
	    && flag_enable_extend)
	  return true;

	if ((fabs(dangle1) < 2.5 || fabs(dangle1_1) < 2.5) &&
	    (fabs(dangle2) < 2.5 || fabs(dangle2_1) < 2.5) &&
	    (fabs(dangle3) < 2.5 || fabs(dangle3_1) < 2.5) &&
	    (length_1 > 25*units::cm && length_2 > 25*units::cm))
	  flag_force_extend = true;
	
	if (flag_para_U || flag_para_V ){
	  if ( (fabs(dangle1) < para_angle_cut_1 || fabs(dangle1_1) < para_angle_cut_1) &&
	       (fabs(dangle2) < para_angle_cut_1 || fabs(dangle2_1) < para_angle_cut_1) &&
	       (fabs(dangle4) < para_angle_cut_1/2. || fabs(dangle4_1) < para_angle_cut_1/2.) &&
	       (angle_diff1 < para_angle_cut || angle_diff1_1 < para_angle_cut)){
	    if (length_2 < 12*units::cm && dis < 2*units::cm) //small and very close
	      return true;
	    
	    if (angle_diff3 < para_angle_cut || angle_diff3_1 < para_angle_cut )
	      return true;
	    
	    flag_extend = true;
	  }
	  
	  
	  if ( (fabs(dangle3) < para_angle_cut_1 || fabs(dangle3_1) < para_angle_cut_1) &&
	       (fabs(dangle2) < para_angle_cut_1 || fabs(dangle2_1) < para_angle_cut_1) &&
	       (fabs(dangle5) < para_angle_cut_1/2. || fabs(dangle5_1) < para_angle_cut_1/2.) &&
	       (angle_diff2 < para_angle_cut || angle_diff2_1 < para_angle_cut)){
	    if (length_1 < 12*units::cm && dis < 2*units::cm) //small and very close
	      return true;
	    
	    if (angle_diff3 < para_angle_cut || angle_diff3_1 < para_angle_cut )
	      return true;
	    
	    flag_extend = true;
	  }

	  if (dis<5*units::cm && length_1 > 25*units::cm && length_2 > 25*units::cm)
	    flag_extend = true;
	}
      }

      if (flag_extend && flag_enable_extend || flag_force_extend){
      	// when to extend???
	
      	if ((flag_para && (flag_para_U || flag_para_V  || flag_regular)) ||
	    (!flag_para && (flag_prolong_U || flag_prolong_V || dis < 10*units::cm))){
      	  // look use 1 to predict 2
      	  // cluster1_ave_pos, dir1
      	  // calculate the average distance
      	  double ave_dis = sqrt(pow(cluster1_ave_pos.x()-cluster2_ave_pos.x(),2) + pow(cluster1_ave_pos.y()-cluster2_ave_pos.y(),2) + pow(cluster1_ave_pos.z()-cluster2_ave_pos.z(),2));
      	  geo_point_t test_point;
      	  double min_dis = 1e9, max_dis = -1e9;
      	  for (int i=-5;i!=6;i++){
      	    test_point.set(cluster1_ave_pos.x() - dir1.x() * (ave_dis +i*2*units::cm),cluster1_ave_pos.y() - dir1.y() * (ave_dis +i*2*units::cm),cluster1_ave_pos.z() - dir1.z() * (ave_dis +i*2*units::cm));
	    
      	    auto temp_results = cluster2.get_closest_point_blob(test_point);
      	    //reuse this
      	    auto test_point1 = temp_results.first;
      	    if (sqrt(pow(test_point1.x()-test_point.x(),2)+pow(test_point1.y()-test_point.y(),2)+pow(test_point1.z()-test_point.z(),2))<1.5*units::cm){
      	      double temp_dis = (test_point1.x() - cluster1_ave_pos.x())*dir1.x() + (test_point1.y() - cluster1_ave_pos.y())*dir1.y() + (test_point1.z() - cluster1_ave_pos.z())*dir1.z();
      	      temp_dis *=-1;
      	      if (temp_dis < min_dis) min_dis = temp_dis;
      	      if (temp_dis > max_dis) max_dis = temp_dis;
      	    }
      	  }
	  // std::cout << cluster1.get_cluster_id() << " " << cluster2.get_cluster_id() << " " << min_dis/units::cm << " " << max_dis/units::cm << " " << length_2/units::cm << std::endl;
	  
      	  if ((max_dis - min_dis)>2.5*units::cm) return true;

	  // look at the other side (repeat)
      	  // cluster2_ave_pos, dir2
      	  min_dis = 1e9;
      	  max_dis = -1e9;
      	  for (int i=-5;i!=6;i++){
      	    test_point.set(cluster2_ave_pos.x() - dir3.x() * (ave_dis +i*2*units::cm), cluster2_ave_pos.y() - dir3.y() * (ave_dis +i*2*units::cm), cluster2_ave_pos.z() - dir3.z() * (ave_dis +i*2*units::cm));
	    
      	    auto temp_results = cluster1.get_closest_point_blob(test_point);
      	    //reuse this
      	    auto test_point1 = temp_results.first;
      	    if (sqrt(pow(test_point1.x()-test_point.x(),2)+pow(test_point1.y()-test_point.y(),2)+pow(test_point1.z()-test_point.z(),2))<1.5*units::cm){
      	      double temp_dis = (test_point1.x() - cluster2_ave_pos.x())*dir3.x() + (test_point1.y() - cluster2_ave_pos.y())*dir3.y() + (test_point1.z() - cluster2_ave_pos.z())*dir3.z();
      	      temp_dis *=-1;
      	      if (temp_dis < min_dis) min_dis = temp_dis;
      	      if (temp_dis > max_dis) max_dis = temp_dis;
      	    }
      	  }
	  // std::cout << cluster1->get_cluster_id() << " " << cluster2.get_cluster_id() << " " << min_dis/units::cm << " " << max_dis/units::cm << " " << length_1/units::cm << std::endl;


      	  if ((max_dis - min_dis)>2.5*units::cm) return true;
	  
      	  // both side simutaneously? leave it for futhre
	  
	}
      }
	
      
    }
  }

  return false;
  
}
#pragma GCC diagnostic pop
