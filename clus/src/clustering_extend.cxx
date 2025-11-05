#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"

class ClusteringExtend;
class ClusteringExtendLoop;
WIRECELL_FACTORY(ClusteringExtend, ClusteringExtend,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)
WIRECELL_FACTORY(ClusteringExtendLoop, ClusteringExtendLoop,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

static
void clustering_extend(Grouping& live_clusters,
                       IDetectorVolumes::pointer dv,      // detector volumes
                       const Tree::Scope& scope,
                       const int flag,                                                //
                       const double length_cut = 150*units::cm,                       //
                       const int num_try = 0,                                         //
                       const double length_2_cut = 3*units::cm,                       //
                       const int num_dead_try =3                                      //
  );

class ClusteringExtend : public IConfigurable, public Clus::IEnsembleVisitor, private NeedDV, private NeedScope {
public:
  ClusteringExtend() {}
  virtual ~ClusteringExtend() {}

  void configure(const WireCell::Configuration& config) {
    NeedDV::configure(config);
    NeedScope::configure(config);

    flag_ = get(config, "flag", 0);
    length_cut_ = get(config, "length_cut", 150*units::cm);
    num_try_ = get(config, "num_try", 0);
    length_2_cut_ = get(config, "length_2_cut", 3*units::cm);
    num_dead_try_ = get(config, "num_dead_try", 3);
  }

  virtual Configuration default_configuration() const {
    Configuration cfg;
    return cfg;
  }

  void visit(Ensemble& ensemble) const {
    auto& live = *ensemble.with_name("live").at(0);
    clustering_extend(live, m_dv, m_scope, flag_, length_cut_, num_try_, length_2_cut_, num_dead_try_);
  }

private:
  int flag_{0};
  double length_cut_{150*units::cm};
  int num_try_{0};
  double length_2_cut_{3*units::cm};
  int num_dead_try_{3};
};
    
class ClusteringExtendLoop : public IConfigurable, public Clus::IEnsembleVisitor, private NeedDV, private NeedScope {
public:
  ClusteringExtendLoop() {}
  virtual ~ClusteringExtendLoop() {}

  void configure(const WireCell::Configuration& config) {
    NeedDV::configure(config);
    NeedScope::configure(config);
    
    num_try_ = get(config, "num_try", 0);
  }

  virtual Configuration default_configuration() const {
    Configuration cfg;
    return cfg;
  }

  void visit(Ensemble& ensemble) const {
    auto& live = *ensemble.with_name("live").at(0);

    // for very busy events do less ...
    int num_try = num_try_;
    if (live.nchildren() > 1100) num_try = 1;
    for (int i = 0; i != num_try; i++) {
      // deal with prolong case
      clustering_extend(live, m_dv, m_scope, 1, 150*units::cm, 0);
      // deal with parallel case
      clustering_extend(live, m_dv, m_scope, 2, 30*units::cm, 0);
      // extension regular case
      clustering_extend(live, m_dv, m_scope, 3, 15*units::cm, 0);
      // extension ones connected to dead region ...
      if (i == 0) {
        clustering_extend(live, m_dv, m_scope, 4, 60 * units::cm, i);
      }
      else {
        clustering_extend(live, m_dv, m_scope, 4, 35 * units::cm, i);
      }
    }
  }

private:
  int num_try_{0};
};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

using namespace WireCell::PointCloud::Tree;

static bool Clustering_4th_prol(
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

static bool Clustering_4th_para(
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

static bool Clustering_4th_reg(
    const Cluster& cluster_1,
    const Cluster& cluster_2,
    double length_1, double length_2,
    geo_point_t p1, double length_cut, 
    const std::map<WirePlaneId, std::pair<geo_point_t, double> > & wpid_U_dir, const std::map<WirePlaneId, std::pair<geo_point_t, double> > & wpid_V_dir, const std::map<WirePlaneId, std::pair<geo_point_t, double> > & wpid_W_dir,                  const IDetectorVolumes::pointer dv)
{
  auto temp_results = cluster_2.get_closest_point_blob(p1);
  geo_point_t p2 = temp_results.first;
  auto wpid_p2 = cluster_2.wpid(p2);
  geo_point_t diff = p1 - p2;
  double dis1 = diff.magnitude();
  
  temp_results = cluster_1.get_closest_point_blob(p2);
  p1 = temp_results.first;
  auto wpid_p1 = cluster_1.wpid(p1);

  auto wpid_ps = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv);

  /* temp_results = cluster_2.get_closest_point_blob(p1); */
  /* p2 = temp_results.second; */

  diff = p1 - p2;
  double dis = diff.magnitude();

  geo_point_t drift_dir_abs(1, 0, 0);  // assuming the drift direction is along X ...

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
    // geo_point_t U_dir(0,cos(angle_u),sin(angle_u));
    // geo_point_t V_dir(0,cos(angle_v),sin(angle_v));
    
    geo_point_t dir2(p2.x()-p1.x(),p2.y()-p1.y(),p2.z()-p1.z());
    
    bool flag_para = false, flag_prol =false, flag_reg = false;
    
    double angle1 = fabs(dir2.angle(drift_dir_abs)-3.1415926/2.)/3.1415926*180.;
    
    if (angle1 < 5 && dis < 2*length_cut || angle1 < 2 ){
      flag_para = true;
    }else if (dis < 2*length_cut){
      geo_point_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
      geo_point_t tempV5;
      double angle2 = tempV1.angle(wpid_U_dir.at(wpid_ps).first);
      tempV5.set(fabs(p2.x()-p1.x()),sqrt(pow(p2.y() - p1.y(),2)+pow(p2.z() - p1.z(),2))*sin(angle2),0);
      angle2 = tempV5.angle(drift_dir_abs)/3.1415926*180.;
      
      double angle3 = tempV1.angle(wpid_V_dir.at(wpid_ps).first);
      tempV5.set(fabs(p2.x()-p1.x()),sqrt(pow(p2.y() - p1.y(),2)+pow(p2.z() - p1.z(),2))*sin(angle3),0);
      angle3 = tempV5.angle(drift_dir_abs)/3.1415926*180.;
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
      
      if (flag_para && fabs(dir3.angle(drift_dir_abs)-3.141592/2.)<10/180.*3.1415926 && fabs(dir1.angle(drift_dir_abs)-3.141592/2.)<10/180.*3.1415926){
        if (angle4 < 30 && (length_2 < 12*units::cm && fabs(angle5-90.)>30 || angle5 < 45))
          return true;
      }else if (flag_prol){
        if (angle4 < 25 && (length_2 < 15*units::cm && fabs(angle5-90.)>30 || angle5 < 25))
          return true;
      }

      if (fabs(dir2.angle(drift_dir_abs)-3.1415926/2.)/3.1415926*180.>7.5){
	// non-parallel case ... 
	if (WireCell::Clus::Facade::is_angle_consistent(dir1,dir2,false,10,wpid_U_dir.at(wpid_p1).second, wpid_V_dir.at(wpid_p1).second, wpid_W_dir.at(wpid_p1).second,2)){
	  if (length_2 < 8*units::cm && WireCell::Clus::Facade::is_angle_consistent(dir1,dir2,false,5,wpid_U_dir.at(wpid_p1).second, wpid_V_dir.at(wpid_p1).second, wpid_W_dir.at(wpid_p1).second,2)) 
		return true; 
	  if (WireCell::Clus::Facade::is_angle_consistent(dir3,dir2,true,10,wpid_U_dir.at(wpid_p2).second, wpid_V_dir.at(wpid_p2).second, wpid_W_dir.at(wpid_p2).second,2)){
	    return true;
	  }
	}
      }
    }
  }
  return false;
}


static bool Clustering_4th_dead(
    const Cluster& cluster_1,
    const Cluster& cluster_2,
    double length_1, double length_2, double length_cut, int num_dead_try,
    const std::map<WirePlaneId, std::pair<geo_point_t, double> > & wpid_U_dir, const std::map<WirePlaneId, std::pair<geo_point_t, double> > & wpid_V_dir, const std::map<WirePlaneId, std::pair<geo_point_t, double> > & wpid_W_dir,                  const IDetectorVolumes::pointer dv)
{
  geo_point_t drift_dir_abs(1, 0, 0);  // assuming the drift direction is along X ...
  // const auto [angle_u,angle_v,angle_w] = cluster_1.grouping()->wire_angles();

  geo_point_t p1;
  geo_point_t p2;
  
  double dis = Find_Closest_Points(cluster_1, cluster_2, length_1, length_2, length_cut, p1, p2);

  // auto wpid_p1 = cluster_1.wpid(p1);
	// auto wpid_p2 = cluster_2.wpid(p2);
	// auto wpid_ps = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv);

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

      auto wpid_ave_p1 = cluster_1.wpid(cluster1_ave_pos);
      auto wpid_ave_p2 = cluster_2.wpid(cluster2_ave_pos);
      // auto wpid_ave_ps = get_wireplaneid(cluster1_ave_pos, wpid_ave_p1, cluster2_ave_pos, wpid_ave_p2, dv);

      // use projection to deal with stuff ...
      if (fabs(ave_dir.angle(drift_dir_abs)-3.1415926/2.)/3.1415926*180.>7.5){
      	// non-parallel case ...
      	if (WireCell::Clus::Facade::is_angle_consistent(dir1,dir2,false,10,wpid_U_dir.at(wpid_ave_p1).second, wpid_V_dir.at(wpid_ave_p1).second, wpid_W_dir.at(wpid_ave_p1).second,2)){
      	  if (length_2 < 8*units::cm&& WireCell::Clus::Facade::is_angle_consistent(dir1,dir2,false,5,wpid_U_dir.at(wpid_ave_p1).second, wpid_V_dir.at(wpid_ave_p1).second, wpid_W_dir.at(wpid_ave_p1).second,2))
      	    return true;
      	  if (length_2 < 15*units::cm && WireCell::Clus::Facade::is_angle_consistent(dir1,dir2,false,7.5,wpid_U_dir.at(wpid_ave_p1).second, wpid_V_dir.at(wpid_ave_p1).second, wpid_W_dir.at(wpid_ave_p1).second))
      	    return true;
      	  if (WireCell::Clus::Facade::is_angle_consistent(dir3,dir2,true,10,wpid_U_dir.at(wpid_ave_p2).second, wpid_V_dir.at(wpid_ave_p2).second, wpid_W_dir.at(wpid_ave_p2).second,2)){
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

	if (fabs(ave_dir.angle(drift_dir_abs)-3.1415926/2.)/3.1415926*180. > 7.5  && ave_dis < 30*units::cm){
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


// Expand this function to handle multiple APA/Faces ...
static void clustering_extend(
    Grouping& live_grouping,
    const IDetectorVolumes::pointer dv,                // detector volumes
    const Tree::Scope& scope,
    const int flag,                                                //
    const double length_cut,                                       //
    const int num_try,                                             //
    const double length_2_cut,                                     //
    const int num_dead_try                                         //
){

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


  geo_point_t drift_dir_abs(1,0,0);
  cluster_set_t used_clusters;


  // prepare graph ...
  typedef cluster_connectivity_graph_t Graph;
  Graph g;
  std::unordered_map<int, int> ilive2desc;  // added live index to graph descriptor
  std::map<const Cluster*, int> map_cluster_index;
  const auto& live_clusters = live_grouping.children();
  
  for (size_t ilive = 0; ilive < live_clusters.size(); ++ilive) {
    auto& live = live_clusters.at(ilive);
    map_cluster_index[live] = ilive;
    ilive2desc[ilive] = boost::add_vertex(ilive, g);
    // set scope ...
    if (live->get_default_scope().hash() != scope.hash()) {
      live->set_default_scope(scope);
    }
  }

  // original algorithm ... (establish edges ... )

  int length_1_cut = 40*units::cm + num_try * 10*units::cm;

  if (flag==1) length_1_cut = 20*units::cm + num_try*10*units::cm; //prolong case
  
  for (size_t i=0;i!=live_clusters.size();i++){
    auto cluster_1 = live_clusters.at(i);
    if (!cluster_1->get_scope_filter(scope)) continue;
   

    if (cluster_1->get_length() > length_1_cut){
      geo_point_t highest_p, lowest_p, earliest_p, latest_p;


      if (flag==1){// prolong case ... 

        std::tie(earliest_p, latest_p) = cluster_1->get_earliest_latest_points();
        // find earliest point

        geo_point_t dir_earlp = cluster_1->vhough_transform(earliest_p,60*units::cm);
        auto wpid_earliest_p = cluster_1->wpid(earliest_p);
        
        geo_point_t tempV5,tempV1;
        tempV1.set(0,dir_earlp.y(),dir_earlp.z());
        double angle1 = tempV1.angle(wpid_U_dir.at(wpid_earliest_p).first);
        tempV5.set(fabs(dir_earlp.x()),sqrt(pow(dir_earlp.y(),2)+pow(dir_earlp.z(),2))*sin(angle1),0);
        angle1 = tempV5.angle(drift_dir_abs);

        double angle2 = tempV1.angle(wpid_V_dir.at(wpid_earliest_p).first);
        tempV5.set(fabs(dir_earlp.x()),sqrt(pow(dir_earlp.y(),2)+pow(dir_earlp.z(),2))*sin(angle2),0);
        angle2 = tempV5.angle(drift_dir_abs);

        double angle3 = tempV1.angle(wpid_W_dir.at(wpid_earliest_p).first);
        tempV5.set(fabs(dir_earlp.x()),sqrt(pow(dir_earlp.y(),2)+pow(dir_earlp.z(),2))*sin(angle3),0);
        angle3 = tempV5.angle(drift_dir_abs);


        // find latest point
        geo_point_t dir_latep = cluster_1->vhough_transform(latest_p, 60*units::cm);
        auto wpid_latest_p = cluster_1->wpid(latest_p);
        tempV1.set(0,dir_latep.y(),dir_latep.z());
        double angle4 = tempV1.angle(wpid_U_dir.at(wpid_latest_p).first);
        tempV5.set(fabs(dir_latep.x()),sqrt(pow(dir_latep.y(),2)+pow(dir_latep.z(),2))*sin(angle4),0);
        angle4 = tempV5.angle(drift_dir_abs);

        double angle5 = tempV1.angle(wpid_V_dir.at(wpid_latest_p).first);
        tempV5.set(fabs(dir_latep.x()),sqrt(pow(dir_latep.y(),2)+pow(dir_latep.z(),2))*sin(angle5),0);
        angle5 = tempV5.angle(drift_dir_abs);

        double angle6 = tempV1.angle(wpid_W_dir.at(wpid_latest_p).first);
        tempV5.set(fabs(dir_latep.x()),sqrt(pow(dir_latep.y(),2)+pow(dir_latep.z(),2))*sin(angle6),0);
        angle6 = tempV5.angle(drift_dir_abs);

        if (angle1 <5./180.*3.1415926 || angle2 < 5./180.*3.1415926 || angle3 < 5./180.*3.1415926){
          // flag_prol = true;
	  
          for (size_t j=0;j!=live_clusters.size();j++){
            auto cluster_2 = live_clusters.at(j);
            if (!cluster_2->get_scope_filter(scope)) continue;

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
            if (!cluster_2->get_scope_filter(scope)) continue;
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

        if (fabs(dir_highp.angle(drift_dir_abs)-3.1415926/2.)<5/180.*3.1415926){ 
          // flag_para = true; 

          for (size_t j=0;j!=live_clusters.size();j++){
            auto cluster_2 = live_clusters.at(j);
            if (!cluster_2->get_scope_filter(scope)) continue;
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

	 if (fabs(dir_lowp.angle(drift_dir_abs)-3.1415926/2.)<5/180.*3.1415926 ){ 
	   // flag_para = true; 

	   for (size_t j=0;j!=live_clusters.size();j++){
	     auto cluster_2 = live_clusters.at(j);
       if (!cluster_2->get_scope_filter(scope)) continue;
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
          if (!cluster_2->get_scope_filter(scope)) continue;

          if (used_clusters.find(cluster_2)!=used_clusters.end()) continue;
          if (cluster_2==cluster_1) continue;
          
          if (Clustering_4th_reg(*cluster_1,*cluster_2,cluster_1->get_length(),cluster_2->get_length(),first_p,length_cut, wpid_U_dir, wpid_V_dir, wpid_W_dir, dv)){
            //	    to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
            boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
                  ilive2desc[map_cluster_index[cluster_2]], g);


                    if (cluster_2->get_length()<10*units::cm)
                        used_clusters.insert(cluster_2);
              
          }else if (Clustering_4th_reg(*cluster_1,*cluster_2,cluster_1->get_length(),cluster_2->get_length(),second_p,length_cut,  wpid_U_dir, wpid_V_dir, wpid_W_dir, dv)){
            //to_be_merged_pairs.insert(std::make_pair(cluster_1,cluster_2));
            boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
                  ilive2desc[map_cluster_index[cluster_2]], g);


                    if (cluster_2->get_length()<10*units::cm)
                        used_clusters.insert(cluster_2);
          }
                  
        }

	
      }else if (flag==4){
        // if (cluster_connected_dead.find(cluster_1)!=cluster_connected_dead.end()){
        if (cluster_1->get_flag(Flags::live_dead)) {
          used_clusters.insert(cluster_1);
          for (size_t j=0;j!=live_clusters.size();j++){
            auto cluster_2 = live_clusters.at(j);
            if (!cluster_2->get_scope_filter(scope)) continue;

            if (cluster_2->get_length() < length_2_cut) continue;
            if (used_clusters.find(cluster_2)!=used_clusters.end()) continue;
            if (Clustering_4th_dead(*cluster_1,*cluster_2,cluster_1->get_length(),cluster_2->get_length(),length_cut,num_dead_try,  wpid_U_dir, wpid_V_dir, wpid_W_dir, dv)){
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

  merge_clusters(g, live_grouping);

  // {
  //  auto live_clusters = live_grouping.children(); // copy
  //   // Process each cluster
  //   for (size_t iclus = 0; iclus < live_clusters.size(); ++iclus) {
  //       Cluster* cluster = live_clusters.at(iclus);
  //       auto& scope = cluster->get_default_scope();
  //       std::cout << "Test: " << iclus << " " << cluster->nchildren() << " " << scope.pcname << " " << scope.coords[0] << " " << scope.coords[1] << " " << scope.coords[2] << " " << cluster->get_scope_filter(scope)<< " " << cluster->get_pca().center) << std::endl;
  //   }
  // }






}








std::vector<std::pair<geo_point_t, const Blob*>> WireCell::Clus::Facade::get_strategic_points(const Cluster& cluster) {
    // Store unique points and their corresponding blobs
    std::set<std::pair<geo_point_t, const Blob*>> unique_points;

    // 1. Get extreme points along main axes
    {
        // Y-axis (highest/lowest)
        auto [high_y, low_y] = cluster.get_highest_lowest_points(1);
        unique_points.emplace(high_y, cluster.blob_with_point(cluster.get_closest_point_index(high_y)));
        unique_points.emplace(low_y, cluster.blob_with_point(cluster.get_closest_point_index(low_y)));
        
        // Z-axis (front/back)
        auto [front, back] = cluster.get_front_back_points();
        unique_points.emplace(front, cluster.blob_with_point(cluster.get_closest_point_index(front)));
        unique_points.emplace(back, cluster.blob_with_point(cluster.get_closest_point_index(back)));
        
        // X-axis (earliest/latest)
        auto [early, late] = cluster.get_earliest_latest_points();
        unique_points.emplace(early, cluster.blob_with_point(cluster.get_closest_point_index(early)));
        unique_points.emplace(late, cluster.blob_with_point(cluster.get_closest_point_index(late)));
    }

    // // 2. Add points based on PCA
    // {
    //     auto center = cluster.get_pca().center);
    //     unique_points.emplace(center, cluster.blob_with_point(cluster.get_closest_point_index(center)));

    //     // Add points along principal axes
    //     for (int i = 0; i < 3; ++i) {
    //         auto dir = cluster.get_pca().axis.at(i);
    //         auto value = sqrt(cluster.get_pca().values.at(i));  // Use eigenvalue for scale
            
    //         // Add points in both directions along each principal axis
    //         geo_point_t p1 = center + dir * value;
    //         geo_point_t p2 = center - dir * value;
            
    //         auto [closest_p1, blob1] = cluster.get_closest_point_blob(p1);
    //         auto [closest_p2, blob2] = cluster.get_closest_point_blob(p2);
            
    //         unique_points.emplace(closest_p1, blob1);
    //         unique_points.emplace(closest_p2, blob2);
    //     }
    // }

    // // 3. Add endpoints from trajectory if available
    // {
    //     const auto& path_wcps = cluster.get_path_wcps();
    //     if (!path_wcps.empty()) {
    //         auto first_idx = path_wcps.front();
    //         auto last_idx = path_wcps.back();
            
    //         auto first_point = cluster.point3d(first_idx);
    //         auto last_point = cluster.point3d(last_idx);
            
    //         unique_points.emplace(first_point, cluster.blob_with_point(first_idx));
    //         unique_points.emplace(last_point, cluster.blob_with_point(last_idx));
    //     }
    // }

    // 4. Add points from convex hull vertices 
    {
        auto hull_points = cluster.get_hull();
        for (const auto& p : hull_points) {
            auto [closest_p, blob] = cluster.get_closest_point_blob(p);
            unique_points.emplace(closest_p, blob);
        }
    }

    // // 5. Add first and last blobs center positions
    // {
    //     if (auto first_blob = cluster.get_first_blob()) {
    //         unique_points.emplace(first_blob->center_pos(), first_blob);
    //     }
    //     if (auto last_blob = cluster.get_last_blob()) {
    //         unique_points.emplace(last_blob->center_pos(), last_blob);
    //     }
    // }

    // Convert set to vector for return
    std::vector<std::pair<geo_point_t, const Blob*>> points;
    points.insert(points.end(), unique_points.begin(), unique_points.end());

    return points;
}

double WireCell::Clus::Facade::Find_Closest_Points(
    const Cluster& cluster1ref,
    const Cluster& cluster2ref,
    double length_1,
    double length_2,
    double length_cut,
    geo_point_t& p1_save,
    geo_point_t& p2_save,
    bool flag_print
    )
{
  const Cluster* cluster1 = &cluster1ref;
  const Cluster* cluster2 = &cluster2ref;

  // Swap clusters if needed to ensure cluster1 is shorter
  bool swapped = false;
  if (length_1 >= length_2) {
      swapped = true;
      std::swap(cluster1, cluster2);
      std::swap(length_1, length_2);
  }

  if (!cluster1->nchildren() || !cluster2->nchildren()) {
      raise<ValueError>("Find_Closest_Points: given empty cluster");
  }

  double dis_save = 1e9;
  
  // Early exit threshold - could be made configurable
  const double early_exit_threshold = 0.5 * units::cm;  
  
  // Get strategic starting points from both clusters
  auto start_points1 = get_strategic_points(*cluster1);
  //auto start_points2 = get_strategic_points(*cluster2);

  // First pass: cluster1 to cluster2
  for (const auto& [start_p1, start_mcell1] : start_points1) {
      const Blob* prev_mcell1 = nullptr;
      const Blob* prev_mcell2 = nullptr;
      const Blob* mcell1 = start_mcell1;
      const Blob* mcell2 = nullptr;
      geo_point_t p1 = start_p1;
      geo_point_t p2;

      // Counter to prevent infinite loops (although unlikely with current logic)
      int iteration_count = 0;
      const int max_iterations = 100;

      // Convergence loop
      while ((mcell1 != prev_mcell1 || mcell2 != prev_mcell2) && iteration_count < max_iterations) {
          prev_mcell1 = mcell1;
          prev_mcell2 = mcell2;

          // Find closest point in cluster2
          auto temp_results = cluster2->get_closest_point_blob(p1);
          p2 = temp_results.first;
          mcell2 = temp_results.second;

          // Find closest point in cluster1
          temp_results = cluster1->get_closest_point_blob(p2);
          p1 = temp_results.first;
          mcell1 = temp_results.second;

          if (flag_print) {
              std::cout << "Forward iter " << iteration_count << ": " 
                        << p1 << " -> " << p2 << std::endl;
          }
          
          iteration_count++;
      }

      double dis = (p1 - p2).magnitude();
        if (dis < dis_save) {
            dis_save = dis;
            p1_save = p1;
            p2_save = p2;
            
            // Early exit if we found points close enough
            if (dis_save < early_exit_threshold) {
                if (flag_print) {
                    std::cout << "Early exit on forward pass with distance: " 
                             << dis_save/units::cm << " cm" << std::endl;
                }
                break;
            }
        }
    }

    // Handle swapped case
    if (swapped) {
        std::swap(p1_save, p2_save);
    }

    return dis_save;

  // double dis_save = 1e9;
  
  // const Blob* prev_mcell1 = 0;
  // const Blob* prev_mcell2 = 0;
  // const Blob* mcell1 = 0;
  // geo_point_t p1;
  // const Blob* mcell2 = 0;
  // geo_point_t p2;
  // double dis;

  // const Cluster* cluster1 = &cluster1ref;
  // const Cluster* cluster2 = &cluster2ref;

  // // Previously the two branches of the following if/else were identical except
  // // for 1<-->2.  For sanity, we simply swap the order.  The output of the
  // // remaining lone block is p1_save/p2_save and if we swap here, we swap them
  // // before returning.
  // bool swapped = false;
  // if (length_1 >= length_2) {
  //   swapped = true;
  //   std::swap(cluster1, cluster2);
  //   std::swap(length_1, length_2);
  // }

  // if (! cluster1->nchildren() || ! cluster2->nchildren()) {
  //   raise<ValueError>("Find_Closest_Points: given empty cluster");
  // }
      
  // mcell1 = cluster1->get_first_blob();
  // p1 = mcell1->center_pos();

  // if (flag_print) std::cout << "a: " << p1 << std::endl;

  // while (mcell1 != prev_mcell1 || mcell2 != prev_mcell2){
  //   prev_mcell1 = mcell1;
  //   prev_mcell2 = mcell2;

  //   auto temp_results = cluster2->get_closest_point_blob(p1);
  //   p2 = temp_results.first;
  //   mcell2 = temp_results.second;

  //   temp_results = cluster1->get_closest_point_blob(p2);
  //   p1 = temp_results.first;
  //   mcell1 = temp_results.second;

  //   if (flag_print) std::cout << "a: " << p1 << " " << p2 << std::endl;
  // }
  // geo_point_t diff = p1 - p2;
  // dis = diff.magnitude();

  // if (dis < dis_save){
  //   dis_save = dis;
  //   p1_save = p1;
  //   p2_save = p2;
  // }

  // prev_mcell1 = 0;
  // prev_mcell2 = 0;

  // mcell1 = cluster1->get_last_blob();
  // p1 = mcell1->center_pos();

  // if (flag_print) std::cout << "b: " << p1 << std::endl;

  // while(mcell1!=prev_mcell1 || mcell2!=prev_mcell2){
  //   prev_mcell1 = mcell1;
  //   prev_mcell2 = mcell2;
      
  //   // find the closest point and merged cell in cluster2
  //   auto temp_results = cluster2->get_closest_point_blob(p1);
  //   p2 = temp_results.first;
  //   mcell2 = temp_results.second;
    
  //   // find the closest point and merged cell in cluster1
  //   temp_results = cluster1->get_closest_point_blob(p2);
  //   p1 = temp_results.first;
  //   mcell1 = temp_results.second;
  
  //   if (flag_print) std::cout << "b: " << p1 << " " << p2 << std::endl;
  // }
  // diff = p1 - p2;
  // dis = diff.magnitude();

  // if (dis < dis_save){
  //   dis_save = dis;
  //   p1_save = p1;
  //   p2_save = p2;
  // }
    
  // if (swapped) {
  //   std::swap(p1_save, p2_save);
  // }

  // return dis_save;
}


#pragma GCC diagnostic pop


// Local Variables:
// mode: c++
// c-basic-offset: 2
// End:
