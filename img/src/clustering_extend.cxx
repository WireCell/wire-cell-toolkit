#include <WireCellImg/ClusteringFuncs.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;
void WireCell::PointCloud::Facade::clustering_extend(
    Grouping& live_grouping,
    cluster_set_t& cluster_connected_dead,     // in/out
    const int flag,                                                //
    const double length_cut,                                       //
    const int num_try,                                             //
    const double length_2_cut,                                     //
    const int num_dead_try                                         //
){
  geo_point_t drift_dir(1, 0, 0);  // assuming the drift direction is along X ...
  const auto [angle_u,angle_v,angle_w] = live_grouping.wire_angles();

  // pronlonged case for U 3 and V 4 ...
  geo_point_t U_dir(0,cos(angle_u),sin(angle_u));
  geo_point_t V_dir(0,cos(angle_v),sin(angle_v));
  geo_point_t W_dir(0,cos(angle_w),sin(angle_w));

  cluster_set_t used_clusters;


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

  int length_1_cut = 40*units::cm + num_try * 10*units::cm;

  if (flag==1) length_1_cut = 20*units::cm + num_try*10*units::cm; //prolong case
  
  for (size_t i=0;i!=live_clusters.size();i++){
    auto cluster_1 = live_clusters.at(i);

    if (cluster_1->get_length() > length_1_cut){
      geo_point_t highest_p, lowest_p, earliest_p, latest_p;
      // bool flag_para = false;
      // bool flag_prol = false;

      if (flag==1){// prolong case ... 

	std::tie(earliest_p, latest_p) = cluster_1->get_earliest_latest_points();
	// find earliest point

	geo_point_t dir_earlp = cluster_1->vhough_transform(earliest_p,60*units::cm);
	
	geo_point_t tempV5,tempV1;
	tempV1.set(0,dir_earlp.y(),dir_earlp.z());
	double angle1 = tempV1.angle(U_dir);
	tempV5.set(fabs(dir_earlp.x()),sqrt(pow(dir_earlp.y(),2)+pow(dir_earlp.z(),2))*sin(angle1),0);
	angle1 = tempV5.angle(drift_dir);

	double angle2 = tempV1.angle(V_dir);
	tempV5.set(fabs(dir_earlp.x()),sqrt(pow(dir_earlp.y(),2)+pow(dir_earlp.z(),2))*sin(angle2),0);
	angle2 = tempV5.angle(drift_dir);

	double angle3 = tempV1.angle(W_dir);
	tempV5.set(fabs(dir_earlp.x()),sqrt(pow(dir_earlp.y(),2)+pow(dir_earlp.z(),2))*sin(angle3),0);
	angle3 = tempV5.angle(drift_dir);


	// find latest point
	geo_point_t dir_latep = cluster_1->vhough_transform(latest_p, 60*units::cm);
	tempV1.set(0,dir_latep.y(),dir_latep.z());
	double angle4 = tempV1.angle(U_dir);
	tempV5.set(fabs(dir_latep.x()),sqrt(pow(dir_latep.y(),2)+pow(dir_latep.z(),2))*sin(angle4),0);
	angle4 = tempV5.angle(drift_dir);

	double angle5 = tempV1.angle(V_dir);
	tempV5.set(fabs(dir_latep.x()),sqrt(pow(dir_latep.y(),2)+pow(dir_latep.z(),2))*sin(angle5),0);
	angle5 = tempV5.angle(drift_dir);

	double angle6 = tempV1.angle(W_dir);
	tempV5.set(fabs(dir_latep.x()),sqrt(pow(dir_latep.y(),2)+pow(dir_latep.z(),2))*sin(angle6),0);
	angle6 = tempV5.angle(drift_dir);

	if (angle1 <5./180.*3.1415926 || angle2 < 5./180.*3.1415926 || angle3 < 5./180.*3.1415926){
	  // flag_prol = true;
	  
	  for (size_t j=0;j!=live_clusters.size();j++){
	    auto cluster_2 = live_clusters.at(j);
	    if (used_clusters.find(cluster_2)!=used_clusters.end()) continue;
	    if (cluster_2==cluster_1) continue;
	    if (Clustering_4th_prol(*cluster_1,*cluster_2,cluster_2->get_length(),earliest_p,dir_earlp,length_cut)){
	      //	      to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	      boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
			      ilive2desc[map_cluster_index[cluster_2]], g);


	      
	      if (cluster_2->get_length()<10*units::cm)
		used_clusters.insert(cluster_2);
	    }
	  }
	}
	
	if (angle4<5./180.*3.1415926 || angle5 < 5./180.*3.1415926 || angle6 < 5./180.*3.1415926){

	  // flag_prol = true;
	  for (size_t j=0;j!=live_clusters.size();j++){
	    auto cluster_2 = live_clusters.at(j);
	    if (used_clusters.find(cluster_2)!=used_clusters.end()) continue;
	    if (cluster_2==cluster_1) continue;
	    if (Clustering_4th_prol(*cluster_1,*cluster_2,cluster_2->get_length(),latest_p,dir_latep,length_cut)){
	      //to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	      boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
			      ilive2desc[map_cluster_index[cluster_2]], g);


	      if (cluster_2->get_length()<10*units::cm)
		used_clusters.insert(cluster_2);
	    }
	  }
	}
      }else if (flag==2){ // parallel case ...
	std::tie(highest_p, lowest_p) = cluster_1->get_highest_lowest_points();
	
	highest_p = cluster_1->calc_ave_pos(highest_p,5*units::cm);
        geo_point_t dir_highp = cluster_1->vhough_transform(highest_p,100*units::cm);

	lowest_p = cluster_1->calc_ave_pos(lowest_p,5*units::cm);
	geo_point_t dir_lowp = cluster_1->vhough_transform(lowest_p, 100*units::cm);

	 if (fabs(dir_highp.angle(drift_dir)-3.1415926/2.)<5/180.*3.1415926){ 
	   // flag_para = true; 

	   for (size_t j=0;j!=live_clusters.size();j++){
	     auto cluster_2 = live_clusters.at(j);
	     if (used_clusters.find(cluster_2)!=used_clusters.end()) continue;
	     if (cluster_2==cluster_1) continue;
	     
	     if (Clustering_4th_para(*cluster_1,*cluster_2,cluster_1->get_length(),cluster_2->get_length(),highest_p,dir_highp,length_cut)){
	       //to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	       boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
			      ilive2desc[map_cluster_index[cluster_2]], g);


	       
              if (cluster_2->get_length()<15*units::cm)
		 used_clusters.insert(cluster_2);
	     }
	   }
	 }

	 if (fabs(dir_lowp.angle(drift_dir)-3.1415926/2.)<5/180.*3.1415926 ){ 
	   // flag_para = true; 

	   for (size_t j=0;j!=live_clusters.size();j++){
	     auto cluster_2 = live_clusters.at(j);
	     if (cluster_2==cluster_1) continue;
	     if (Clustering_4th_para(*cluster_1,*cluster_2,cluster_1->get_length(),cluster_2->get_length(),lowest_p,dir_lowp,length_cut)){
	       // to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	       boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
			       ilive2desc[map_cluster_index[cluster_2]], g);


	     }
	   }
	   
	 }
      }else if (flag==3){ // regular case ...
	auto hl_ps = cluster_1->get_highest_lowest_points();
	auto el_ps = cluster_1->get_earliest_latest_points();

	geo_point_t first_p, second_p;

		
	if (pow(hl_ps.first.x()-hl_ps.second.x(),2)+pow(hl_ps.first.y()-hl_ps.second.y(),2)+pow(hl_ps.first.z()-hl_ps.second.z(),2) > pow(el_ps.first.x()-el_ps.second.x(),2)+pow(el_ps.first.y()-el_ps.second.y(),2)+pow(el_ps.first.z()-el_ps.second.z(),2)){
	  first_p = hl_ps.first;
	  second_p = hl_ps.second;
	}else{
	  first_p = el_ps.first;
	  second_p = el_ps.second;
	}

	for (size_t j=0;j!=live_clusters.size();j++){
	  auto cluster_2 = live_clusters.at(j);
	  if (used_clusters.find(cluster_2)!=used_clusters.end()) continue;
	  if (cluster_2==cluster_1) continue;
	  
	  if (Clustering_4th_reg(*cluster_1,*cluster_2,cluster_1->get_length(),cluster_2->get_length(),first_p,length_cut)){
	    //	    to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	    boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
			       ilive2desc[map_cluster_index[cluster_2]], g);


               if (cluster_2->get_length()<10*units::cm)
                   used_clusters.insert(cluster_2);
	      
	  }else if (Clustering_4th_reg(*cluster_1,*cluster_2,cluster_1->get_length(),cluster_2->get_length(),second_p,length_cut)){
	    //to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	    boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
			       ilive2desc[map_cluster_index[cluster_2]], g);


               if (cluster_2->get_length()<10*units::cm)
                   used_clusters.insert(cluster_2);
	  }
		     	  
	}

	
      }else if (flag==4){
	if (cluster_connected_dead.find(cluster_1)!=cluster_connected_dead.end()){
	  used_clusters.insert(cluster_1);
	  for (size_t j=0;j!=live_clusters.size();j++){
	    auto cluster_2 = live_clusters.at(j);
	    if (cluster_2->get_length() < length_2_cut) continue;
	    if (used_clusters.find(cluster_2)!=used_clusters.end()) continue;
	    if (Clustering_4th_dead(*cluster_1,*cluster_2,cluster_1->get_length(),cluster_2->get_length(),length_cut,num_dead_try)){
	      //	      to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
	      boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
			       ilive2desc[map_cluster_index[cluster_2]], g);


               if (cluster_2->get_length()<10*units::cm)
                   used_clusters.insert(cluster_2);
	    }
	  }
	}
      }
    }
  }
  // new function to  merge clusters ...
  merge_clusters(g, live_grouping, cluster_connected_dead);
}





bool WireCell::PointCloud::Facade::Clustering_4th_prol(
    const Cluster& cluster_1,
    const Cluster& cluster_2,
    double length_2,
    geo_point_t& earliest_p,
    geo_point_t& dir_earlp,
    double length_cut)
{
  auto temp_results = cluster_2.get_closest_point_blob(earliest_p);
  geo_point_t p2 = temp_results.first;
  geo_point_t diff = earliest_p - p2;
  double dis = diff.magnitude();
  if (dis < length_cut){
    geo_point_t dir_bp(p2.x()-earliest_p.x(),p2.y()-earliest_p.y(),p2.z()-earliest_p.z());
    double angle_diff = (3.1415926-dir_bp.angle( dir_earlp))/3.1415926*180.;
    if ( (angle_diff < 3 || angle_diff>177 || 
	  dis * sin(angle_diff/180.*3.1415926) < 6*units::cm)){
      geo_point_t dir = cluster_2.vhough_transform(p2,60*units::cm);
      if (length_2<10*units::cm && fabs(dir.angle(dir_earlp)-3.141926/2.)>30/180.*3.1415926){
	return true;
      }else{
	if ((3.14151926-dir.angle(dir_earlp))/3.1415926*180.<5. ||
	    dir.angle(dir_earlp)/3.1415926*180.<5.)
	  return true;
      }
    }
    
  }
    
  return false;
  
  
}

bool WireCell::PointCloud::Facade::Clustering_4th_para(
    const Cluster& cluster_1,
    const Cluster& cluster_2,
    double length_1, double length_2,
    geo_point_t& earliest_p,
    geo_point_t& dir_earlp,
    double length_cut)
{

  auto temp_results = cluster_2.get_closest_point_blob(earliest_p);
  geo_point_t p2 = temp_results.first;
  
  geo_point_t diff = p2 - earliest_p;
  double dis = diff.magnitude();
  
  if (dis < length_cut ){
    
    geo_point_t test_point; 
    double min_dis = 1e9, max_dis = -1e9;

     for (int i=-5;i!=10;i++){ 
       test_point.set(earliest_p.x() - dir_earlp.x() * (dis +i*2*units::cm), earliest_p.y() - dir_earlp.y() * (dis +i*2*units::cm), earliest_p.z() - dir_earlp.z() * (dis +i*2*units::cm)); 
	
       auto temp_results = cluster_2.get_closest_point_blob(test_point); 
 	//reuse this 
       geo_point_t test_point1 = temp_results.first;

       if (sqrt(pow(test_point1.x()-test_point.x(),2)+pow(test_point1.y()-test_point.y(),2)+pow(test_point1.z()-test_point.z(),2))<1.5*units::cm){ 
	 double temp_dis = (test_point1.x() - earliest_p.x())*dir_earlp.x() + (test_point1.y() - earliest_p.y())*dir_earlp.y() + (test_point1.z() - earliest_p.z())*dir_earlp.z(); 
	 temp_dis =(-1) * temp_dis; 
	 if (temp_dis < min_dis) min_dis = temp_dis; 
	 if (temp_dis > max_dis) max_dis = temp_dis; 
       } 
     }
     
     if ((max_dis - min_dis)>2.5*units::cm) return true; 
     
  }

  return false;
}

bool WireCell::PointCloud::Facade::Clustering_4th_reg(
    const Cluster& cluster_1,
    const Cluster& cluster_2,
    double length_1, double length_2,
    geo_point_t p1, double length_cut)
{
  auto temp_results = cluster_2.get_closest_point_blob(p1);
  geo_point_t p2 = temp_results.first;
  geo_point_t diff = p1 - p2;
  double dis1 = diff.magnitude();
  
  temp_results = cluster_1.get_closest_point_blob(p2);
  p1 = temp_results.first;
  /* temp_results = cluster_2.get_closest_point_blob(p1); */
  /* p2 = temp_results.second; */

  diff = p1 - p2;
  double dis = diff.magnitude();

  geo_point_t drift_dir(1, 0, 0);  // assuming the drift direction is along X ...
  const auto [angle_u,angle_v,angle_w] = cluster_1.grouping()->wire_angles();


  if (dis1 > 15*units::cm && dis < 3*units::cm && length_2 > 80*units::cm &&length_1>80*units::cm) return false;
  
  if (dis < length_cut && (length_2 >= 40*units::cm || dis < 3*units::cm)){
    geo_point_t cluster1_ave_pos = cluster_1.calc_ave_pos(p1,5*units::cm);
    geo_point_t cluster2_ave_pos = cluster_2.calc_ave_pos(p2,5*units::cm);
    geo_point_t dir1;
    
    if (cluster_1.nnearby(cluster1_ave_pos, 30*units::cm)>50 && length_1 < 120*units::cm){
      dir1 = cluster_1.vhough_transform(cluster1_ave_pos,30*units::cm);
    }else{
      dir1 = cluster_1.vhough_transform(cluster1_ave_pos,80*units::cm);
    }

    geo_point_t dir3;
    if (cluster_2.nnearby(cluster2_ave_pos, 30*units::cm)>50&&length_2 < 120*units::cm){
      dir3 = cluster_2.vhough_transform(cluster2_ave_pos,30*units::cm);
    }else{
      dir3 = cluster_2.vhough_transform(cluster2_ave_pos,80*units::cm);
    }

    geo_point_t dir2(cluster2_ave_pos.x() - cluster1_ave_pos.x(),
		     cluster2_ave_pos.y() - cluster1_ave_pos.y(),
		     cluster2_ave_pos.z() - cluster1_ave_pos.z());

    
    double ave_dis = sqrt(pow(cluster1_ave_pos.x()-cluster2_ave_pos.x(),2) + pow(cluster1_ave_pos.y()-cluster2_ave_pos.y(),2) + pow(cluster1_ave_pos.z()-cluster2_ave_pos.z(),2));
    geo_point_t test_point;
    double min_dis = 1e9, max_dis = -1e9;

    if (dir2.angle(dir1)>3.1415926/2. ){
          
      for (int i=-5;i!=10;i++){
	test_point.set(cluster1_ave_pos.x() - dir1.x() * (ave_dis +i*2*units::cm), cluster1_ave_pos.y() - dir1.y() * (ave_dis +i*2*units::cm), cluster1_ave_pos.z() - dir1.z() * (ave_dis +i*2*units::cm));
	
	auto temp_results = cluster_2.get_closest_point_blob(test_point);
	//reuse this
	geo_point_t test_point1 = temp_results.first;
	if (sqrt(pow(test_point1.x()-test_point.x(),2)+pow(test_point1.y()-test_point.y(),2)+pow(test_point1.z()-test_point.z(),2))<1.5*units::cm){
	  double temp_dis = (test_point1.x() - cluster1_ave_pos.x())*dir1.x() + (test_point1.y() - cluster1_ave_pos.y())*dir1.y() + (test_point1.z() - cluster1_ave_pos.z())*dir1.z();
	  temp_dis *=-1;
	  if (temp_dis < min_dis) min_dis = temp_dis;
	  if (temp_dis > max_dis) max_dis = temp_dis;
	}
      }

      if ((max_dis - min_dis)>2.5*units::cm) return true;
    }

    if (dir2.angle(dir3)<3.1415926/2.){
            
      // look at the other side (repeat)
      // cluster2_ave_pos, dir2
      min_dis = 1e9;
      max_dis = -1e9;
      for (int i=-5;i!=10;i++){
	test_point.set(cluster2_ave_pos.x() - dir3.x() * (ave_dis +i*2*units::cm), cluster2_ave_pos.y() - dir3.y() * (ave_dis +i*2*units::cm), cluster2_ave_pos.z() - dir3.z() * (ave_dis +i*2*units::cm));
	
	auto temp_results = cluster_1.get_closest_point_blob(test_point);
	//reuse this
	geo_point_t test_point1 = temp_results.first;
	if (sqrt(pow(test_point1.x()-test_point.x(),2)+pow(test_point1.y()-test_point.y(),2)+pow(test_point1.z()-test_point.z(),2))<1.5*units::cm){
	  double temp_dis = (test_point1.x() - cluster2_ave_pos.x())*dir3.x() + (test_point1.y() - cluster2_ave_pos.y())*dir3.y() + (test_point1.z() - cluster2_ave_pos.z())*dir3.z();
	  temp_dis *=-1;
	  if (temp_dis < min_dis) min_dis = temp_dis;
	  if (temp_dis > max_dis) max_dis = temp_dis;
	}
      }
      
      if ((max_dis - min_dis)>2.5*units::cm) return true;
    }
    
  }else if (dis < 2 * length_cut && length_2 < 40*units::cm){

    // pronlonged case for U 3 and V 4 ...
    geo_point_t U_dir(0,cos(angle_u),sin(angle_u));
    geo_point_t V_dir(0,cos(angle_v),sin(angle_v));
    
    geo_point_t dir2(p2.x()-p1.x(),p2.y()-p1.y(),p2.z()-p1.z());
    bool flag_para = false, flag_prol =false, flag_reg = false;
    
    double angle1 = fabs(dir2.angle(drift_dir)-3.1415926/2.)/3.1415926*180.;
    
    if (angle1 < 5 && dis < 2*length_cut || angle1 < 2 ){
      flag_para = true;
    }else if (dis < 2*length_cut){
      geo_point_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
      geo_point_t tempV5;
      double angle2 = tempV1.angle(U_dir);
      tempV5.set(fabs(p2.x()-p1.x()),sqrt(pow(p2.y() - p1.y(),2)+pow(p2.z() - p1.z(),2))*sin(angle2),0);
      angle2 = tempV5.angle(drift_dir)/3.1415926*180.;
      
      double angle3 = tempV1.angle(V_dir);
      tempV5.set(fabs(p2.x()-p1.x()),sqrt(pow(p2.y() - p1.y(),2)+pow(p2.z() - p1.z(),2))*sin(angle3),0);
      angle3 = tempV5.angle(drift_dir)/3.1415926*180.;
      if (angle2<7.5 || angle3 < 7.5)
  	flag_prol = true;

    }


    if (flag_para || flag_prol || flag_reg){

      geo_point_t dir1;
      if (cluster_1.nnearby(p1, 15*units::cm)>30 && (flag_prol ||flag_reg) ){
	dir1 = cluster_1.vhough_transform(p1,15*units::cm);
      }else{
	dir1 = cluster_1.vhough_transform(p1,60*units::cm);
      }

      geo_point_t dir3;
      if (cluster_2.nnearby(p2, 15*units::cm)>30 && (flag_prol || flag_reg)){
	dir3 = cluster_2.vhough_transform(p2,15*units::cm);
      }else{
	dir3 = cluster_2.vhough_transform(p2,60*units::cm);
      }

       
      double angle4 = (3.1415926-dir1.angle(dir2))/3.1415926*180.;
      double angle5 = dir2.angle(dir3)/3.1415926*180.;
      
      if (flag_para && fabs(dir3.angle(drift_dir)-3.141592/2.)<10/180.*3.1415926 && fabs(dir1.angle(drift_dir)-3.141592/2.)<10/180.*3.1415926){
  	if (angle4 < 30 && (length_2 < 12*units::cm && fabs(angle5-90.)>30 || angle5 < 45))
  	  return true;
      }else if (flag_prol){
  	if (angle4 < 25 && (length_2 < 15*units::cm && fabs(angle5-90.)>30 || angle5 < 25))
  	  return true;
      }

      if (fabs(dir2.angle(drift_dir)-3.1415926/2.)/3.1415926*180.>7.5){
	// non-parallel case ... 
	if (WireCell::PointCloud::Facade::is_angle_consistent(dir1,dir2,false,10,angle_u,angle_v,angle_w,2)){
	  if (length_2 < 8*units::cm && WireCell::PointCloud::Facade::is_angle_consistent(dir1,dir2,false,5,angle_u,angle_v,angle_w,2)) 
		return true; 
	  if (WireCell::PointCloud::Facade::is_angle_consistent(dir3,dir2,true,10,angle_u,angle_v,angle_w,2)){
	    return true;
	  }
	}
      }
    }
  }
  return false;
}

double WireCell::PointCloud::Facade::Find_Closest_Points(
    const Cluster& cluster1ref,
    const Cluster& cluster2ref,
    double length_1,
    double length_2,
    double length_cut,
    geo_point_t& p1_save,
    geo_point_t& p2_save
    )
{
  double dis_save = 1e9;
  
  const Blob* prev_mcell1 = 0;
  const Blob* prev_mcell2 = 0;
  const Blob* mcell1 = 0;
  geo_point_t p1;
  const Blob* mcell2 = 0;
  geo_point_t p2;
  double dis;

  const Cluster* cluster1 = &cluster1ref;
  const Cluster* cluster2 = &cluster2ref;

  // Previously the two branches of the following if/else were identical except
  // for 1<-->2.  For sanity, we simply swap the order.  The output of the
  // remaining lone block is p1_save/p2_save and if we swap here, we swap them
  // before returning.
  bool swapped = false;
  if (length_1 >= length_2) {
    swapped = true;
    std::swap(cluster1, cluster2);
    std::swap(length_1, length_2);
  }

  if (! cluster1->nchildren() || ! cluster2->nchildren()) {
    raise<ValueError>("Find_Closest_Points: given empty cluster");
  }
      
  mcell1 = cluster1->get_first_blob();
  p1 = mcell1->center_pos();

  while (mcell1 != prev_mcell1 || mcell2 != prev_mcell2){
    prev_mcell1 = mcell1;
    prev_mcell2 = mcell2;

    auto temp_results = cluster2->get_closest_point_blob(p1);
    p2 = temp_results.first;
    mcell2 = temp_results.second;

    temp_results = cluster1->get_closest_point_blob(p2);
    p1 = temp_results.first;
    mcell1 = temp_results.second;
  }
  geo_point_t diff = p1 - p2;
  dis = diff.magnitude();

  if (dis < dis_save){
    dis_save = dis;
    p1_save = p1;
    p2_save = p2;
  }

  prev_mcell1 = 0;
  prev_mcell2 = 0;

  mcell1 = cluster1->get_last_blob();
  p1 = mcell1->center_pos();

  while(mcell1!=prev_mcell1 || mcell2!=prev_mcell2){
    prev_mcell1 = mcell1;
    prev_mcell2 = mcell2;
      
    // find the closest point and merged cell in cluster2
    auto temp_results = cluster2->get_closest_point_blob(p1);
    p2 = temp_results.first;
    mcell2 = temp_results.second;
    // find the closest point and merged cell in cluster1
    temp_results = cluster1->get_closest_point_blob(p2);
    p1 = temp_results.first;
    mcell1 = temp_results.second;
  }
  diff = p1 - p2;
  dis = diff.magnitude();

  if (dis < dis_save){
    dis_save = dis;
    p1_save = p1;
    p2_save = p2;
  }
    
  if (swapped) {
    std::swap(p1_save, p2_save);
  }

  return dis_save;
}


bool WireCell::PointCloud::Facade::Clustering_4th_dead(
    const Cluster& cluster_1,
    const Cluster& cluster_2,
    double length_1, double length_2, double length_cut, int num_dead_try)
{
  geo_point_t drift_dir(1, 0, 0);  // assuming the drift direction is along X ...
  const auto [angle_u,angle_v,angle_w] = cluster_1.grouping()->wire_angles();

  geo_point_t p1;
  geo_point_t p2;
  
  double dis = Find_Closest_Points(cluster_1, cluster_2, length_1, length_2, length_cut, p1, p2);

  //add a special one ...  for uboone ...
  /*
  if (length_1 > 30*units::cm && length_2 > 30*units::cm &&
      (dis < 3*units::cm ||
       fabs(p1.x()-p2.x()) < 1.6*units::cm && (dis < 20*units::cm &&
						     p1.z() > 700.6*units::cm && p1.z() < 739.6*units::cm &&     // special region
						     p2.z() > 700.6*units::cm && p2.z() < 739.6*units::cm &&     // special region ...
						     p1.y() > -10.4*units::cm && p1.y() < 29*units::cm &&
						     p2.y() > -10.4*units::cm && p2.y() < 29*units::cm )
       )){
    return true; 
  }
  */

  if ((dis < length_cut || (length_2 > 50*units::cm && dis < 80*units::cm))){

    geo_point_t cluster1_ave_pos_save;
    geo_point_t cluster2_ave_pos_save;
    geo_point_t dir1_save;
    geo_point_t dir3_save;

    for (int i=0;i!=num_dead_try;i++){
      geo_point_t cluster1_ave_pos; 
      geo_point_t cluster2_ave_pos; 
      
      geo_point_t dir1; 
      geo_point_t dir3; 
      geo_point_t dir2;

      if (i==0){
	cluster1_ave_pos = cluster_1.calc_ave_pos(p1,5*units::cm);
	cluster1_ave_pos_save = cluster1_ave_pos;
	cluster2_ave_pos = cluster_2.calc_ave_pos(p2,5*units::cm);
	cluster2_ave_pos_save = cluster2_ave_pos;

	if (num_dead_try==1){
	  dir1 = cluster_1.vhough_transform(cluster1_ave_pos,20*units::cm);
	  dir3 = cluster_2.vhough_transform(cluster2_ave_pos,20*units::cm);
	}else{
	  dir1 = cluster_1.vhough_transform(cluster1_ave_pos,80*units::cm);
	  dir3 = cluster_2.vhough_transform(cluster2_ave_pos,80*units::cm);
	}
	dir1_save = dir1;
	dir3_save = dir3;

	dir2.set(cluster2_ave_pos.x() - cluster1_ave_pos.x()+1e-9, cluster2_ave_pos.y() - cluster1_ave_pos.y()+1e-9, cluster2_ave_pos.z() - cluster1_ave_pos.z()+1e-9); // 2-1
	
      }else if (i==1){
	if (length_2 >= 15*units::cm &&(!(length_2 > 150*units::cm && dis<15*units::cm))){
	  cluster1_ave_pos = cluster1_ave_pos_save;//cluster_1->calc_ave_pos(p1,5*units::cm);
	  dir1 = dir1_save;//cluster_1->VHoughTrans(cluster1_ave_pos,80*units::cm);
	  
	  geo_vector_t dir_test(dir1);
	  dir_test = dir_test/dir_test.magnitude();
	  dir_test = (-1) * dir_test;
	  
	  std::pair<geo_point_t, double> temp_results = cluster_2.get_closest_point_along_vec(cluster1_ave_pos, dir_test, dis*2, 5*units::cm, 15, 10*units::cm);
	  
	  if (temp_results.second < 100*units::cm){
	    cluster2_ave_pos = cluster_2.calc_ave_pos(temp_results.first,5*units::cm);
	    dir3 = cluster_2.vhough_transform(cluster2_ave_pos,80*units::cm);
	    dir2.set(cluster2_ave_pos.x() - cluster1_ave_pos.x()+1e-9, cluster2_ave_pos.y() - cluster1_ave_pos.y()+1e-9, cluster2_ave_pos.z() - cluster1_ave_pos.z()+1e-9); // 2-1
	    
	  }else{
	    continue;
	  }
	}else{
	  continue;
	}
      }else if (i==2){
	if (length_2 >=15*units::cm&&(!(length_2 > 150*units::cm && dis<15*units::cm))){
	  cluster2_ave_pos = cluster2_ave_pos_save;//cluster_2->calc_ave_pos(p2,5*units::cm);
	  dir3 = dir3_save;//cluster_2->VHoughTrans(cluster2_ave_pos,80*units::cm);
	  
	  geo_point_t dir_test(dir3);
	  dir_test = dir_test / dir_test.magnitude();
	  dir_test = (-1) * dir_test;
	  
	  std::pair<geo_point_t, double> temp_results = cluster_1.get_closest_point_along_vec(cluster2_ave_pos, dir_test, dis*2, 5*units::cm, 15, 10*units::cm);
	  
	  if (temp_results.second < 100*units::cm){
	    cluster1_ave_pos = cluster_1.calc_ave_pos(temp_results.first,5*units::cm);
	    dir1 = cluster_1.vhough_transform(cluster1_ave_pos,80*units::cm);
	    dir2.set(cluster2_ave_pos.x() - cluster1_ave_pos.x()+1e-9, cluster2_ave_pos.y() - cluster1_ave_pos.y()+1e-9, cluster2_ave_pos.z() - cluster1_ave_pos.z()+1e-9); // 2-1
	    
	  }else{
	    continue;
	  }
	}else{
	  continue;
	}
      }
      
      geo_point_t ave_dir(cluster1_ave_pos.x()-cluster2_ave_pos.x(),cluster1_ave_pos.y()-cluster2_ave_pos.y(),cluster1_ave_pos.z()-cluster2_ave_pos.z());


      // use projection to deal with stuff ...
      if (fabs(ave_dir.angle(drift_dir)-3.1415926/2.)/3.1415926*180.>7.5){
      	// non-parallel case ...
      	if (WireCell::PointCloud::Facade::is_angle_consistent(dir1,dir2,false,10,angle_u,angle_v,angle_w,2)){
      	  if (length_2 < 8*units::cm&& WireCell::PointCloud::Facade::is_angle_consistent(dir1,dir2,false,5,angle_u,angle_v,angle_w,2))
      	    return true;
      	  if (length_2 < 15*units::cm && WireCell::PointCloud::Facade::is_angle_consistent(dir1,dir2,false,7.5,angle_u,angle_v,angle_w,2))
      	    return true;
      	  if (WireCell::PointCloud::Facade::is_angle_consistent(dir3,dir2,true,10,angle_u,angle_v,angle_w,2)){
      	    return true;
      	  }
      	}
      }
      
      
      double angle1 = (3.1415926-dir1.angle(dir2))/3.1415926*180.;
      double angle2 = dir3.angle(dir2)/3.1415926*180.;
      double angle3 = (3.1415926-dir1.angle(dir3))/3.1415926*180.;

      
      if (length_2 <=10*units::cm){
	if (angle1 < 15 && (angle2 < 60 || length_2 < 5*units::cm) ) return true;
      }else{
	if (angle1 < 15 && angle2 <15 && angle3 < 25 ||
	    angle3 < 10 && (angle1+angle2)<45 && dis < 5*units::cm )
	  return true;

	double ave_dis = sqrt(pow(cluster1_ave_pos.x()-cluster2_ave_pos.x(),2) + pow(cluster1_ave_pos.y()-cluster2_ave_pos.y(),2) + pow(cluster1_ave_pos.z()-cluster2_ave_pos.z(),2));
	geo_point_t test_point;
	double min_dis = 1e9, max_dis = -1e9;

	if (fabs(ave_dir.angle(drift_dir)-3.1415926/2.)/3.1415926*180. > 7.5  && ave_dis < 30*units::cm){
	  if (i==1){
	    for (int k=-5;k!=10;k++){
	      test_point.set(cluster1_ave_pos.x() - dir1.x() * (ave_dis +k*2*units::cm), cluster1_ave_pos.y() - dir1.y() * (ave_dis +k*2*units::cm), cluster1_ave_pos.z() - dir1.z() * (ave_dis +k*2*units::cm));
	      
	      auto temp_results = cluster_2.get_closest_point_blob(test_point);
	      //reuse this
	      geo_point_t test_point1 = temp_results.first;
	      if (sqrt(pow(test_point1.x()-test_point.x(),2)+pow(test_point1.y()-test_point.y(),2)+pow(test_point1.z()-test_point.z(),2))<1.5*units::cm){
		double temp_dis = (test_point1.x() - cluster1_ave_pos.x())*dir1.x() + (test_point1.y() - cluster1_ave_pos.y())*dir1.y() + (test_point1.z() - cluster1_ave_pos.z())*dir1.z();
		temp_dis *=-1;
		if (temp_dis < min_dis) min_dis = temp_dis;
		if (temp_dis > max_dis) max_dis = temp_dis;
	      }
	    }
	    
	    if ((max_dis - min_dis)>2.5*units::cm) return true;
	  }else if (i==2){
	    for (int k=-5;k!=10;k++){
	      test_point.set(cluster2_ave_pos.x() - dir3.x() * (ave_dis +k*2*units::cm), cluster2_ave_pos.y() - dir3.y() * (ave_dis +k*2*units::cm),  cluster2_ave_pos.z() - dir3.z() * (ave_dis +k*2*units::cm));
	      
	      auto temp_results = cluster_1.get_closest_point_blob(test_point);
	      //reuse this
	      geo_point_t test_point1 = temp_results.first;
	      if (sqrt(pow(test_point1.x()-test_point.x(),2)+pow(test_point1.y()-test_point.y(),2)+pow(test_point1.z()-test_point.z(),2))<1.5*units::cm){
		double temp_dis = (test_point1.x() - cluster2_ave_pos.x())*dir3.x() + (test_point1.y() - cluster2_ave_pos.y())*dir3.y() + (test_point1.z() - cluster2_ave_pos.z())*dir3.z();
		temp_dis *=-1;
		if (temp_dis < min_dis) min_dis = temp_dis;
		if (temp_dis > max_dis) max_dis = temp_dis;
	      }
	    }
	     
	    if ((max_dis - min_dis)>2.5*units::cm) return true;
	  }

	}
	
      }
    }
  }

  return false;
}
#pragma GCC diagnostic pop


// Local Variables:
// mode: c++
// c-basic-offset: 2
// End:
