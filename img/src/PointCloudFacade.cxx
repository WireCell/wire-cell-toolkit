#include "WireCellImg/PointCloudFacade.h"

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Facade;
// using WireCell::PointCloud::Dataset;
using namespace WireCell::PointCloud::Tree; // for "Points" node value type
// using WireCell::PointCloud::Tree::named_pointclouds_t;

#include "WireCellUtil/Logging.h"
using spdlog::debug;

namespace {
    // helper to dump a dataset
    std::string dump_ds(const WireCell::PointCloud::Dataset& ds) {
        std::stringstream ss;
        for (const auto& key : ds.keys()) {;
            const auto& arr = ds.get(key);
            ss << " {" << key << ":" << arr->dtype() << ":" << arr->shape()[0] << "} ";
            // const auto& arr = ds.get(key)->elements<float>();
            // for(auto elem : arr) {
            //     ss << elem << " ";
            // }
        }
        return ss.str();
    }
    std::string dump_pcs(const ScopedBase::pointclouds_t& pcs) {
        std::stringstream ss;
        for (const auto& pc : pcs) {
            ss << dump_ds(pc) << std::endl;
        }
        return ss.str();
    }
}

Blob::Blob(const node_ptr& n)
  : m_node(n.get())
{
    const auto& lpcs = m_node->value.local_pcs();
    const auto& pc_scalar = lpcs.at("scalar");
    // std::cout << "pc_scalar " << dump_ds(pc_scalar) << std::endl;
    charge = pc_scalar.get("charge")->elements<float_t>()[0];
    center_x = pc_scalar.get("center_x")->elements<float_t>()[0];
    center_y = pc_scalar.get("center_y")->elements<float_t>()[0];
    center_z = pc_scalar.get("center_z")->elements<float_t>()[0];
    npoints = pc_scalar.get("npoints")->elements<int_t>()[0];
    slice_index_min = pc_scalar.get("slice_index_min")->elements<int_t>()[0];
    slice_index_max = pc_scalar.get("slice_index_max")->elements<int_t>()[0];
    u_wire_index_min = pc_scalar.get("u_wire_index_min")->elements<int_t>()[0];
    u_wire_index_max = pc_scalar.get("u_wire_index_max")->elements<int_t>()[0];
    v_wire_index_min = pc_scalar.get("v_wire_index_min")->elements<int_t>()[0];
    v_wire_index_max = pc_scalar.get("v_wire_index_max")->elements<int_t>()[0];
    w_wire_index_min = pc_scalar.get("w_wire_index_min")->elements<int_t>()[0];
    w_wire_index_max = pc_scalar.get("w_wire_index_max")->elements<int_t>()[0];
    ///
    ///  MAKE SURE YOU UPDATE doctest_clustering_prototype.cxx if you change the above.
    ///
}

bool Blob::overlap_fast(const Blob& b, const int offset) const
{
    if (u_wire_index_min > b.u_wire_index_max + offset) return false;
    if (b.u_wire_index_min > u_wire_index_max + offset) return false;
    if (v_wire_index_min > b.v_wire_index_max + offset) return false;
    if (b.v_wire_index_min > v_wire_index_max + offset) return false;
    if (w_wire_index_min > b.w_wire_index_max + offset) return false;
    if (b.w_wire_index_min > w_wire_index_max + offset) return false;
    return true;
}

geo_point_t Blob::center_pos() const {
    return {center_x, center_y, center_z};
}

int_t Blob::num_points() const{
  return npoints;
}


Cluster::Cluster(const node_ptr& n)
  : m_node(n.get())
{
    // build blobs
    for (const auto& child : m_node->children()) {
        auto blob = std::make_shared<Blob>(child);
        m_blobs.push_back(blob);
        for (int slice_index = blob->slice_index_min; slice_index < blob->slice_index_max; ++slice_index) {
            m_time_blob_map.insert({slice_index, blob});
        }
    }
}

Blob::vector Cluster::is_connected(const Cluster& c, const int offset) const
{
    Blob::vector ret;
    // loop m_time_blob_map
    for (const auto& [time, blob] : m_time_blob_map) {
        // loop c.m_time_blob_map
        auto range = c.m_time_blob_map.equal_range(time);
	bool flag_save = false;
        for (auto it = range.first; it != range.second; ++it) {
	  const auto& cblob = it->second;
	  if (blob->overlap_fast(*cblob, offset)) {
	    //ret.push_back(cblob); // dead clusters ... 
	    //ret.push_back(blob); // live clusters ...
	    flag_save = true;
	    break;
	  }
        }
	if (flag_save) ret.push_back(blob); // live clusters ...
    }
    return ret;
}

std::shared_ptr<const WireCell::PointCloud::Facade::Blob> Cluster::get_first_blob() const{
  return m_time_blob_map.begin()->second;
}
std::shared_ptr<const WireCell::PointCloud::Facade::Blob> Cluster::get_last_blob() const{
  return m_time_blob_map.rbegin()->second;
}


std::pair<geo_point_t, double> Cluster::get_closest_point_along_vec(geo_point_t& p_test1, geo_point_t dir, double test_dis, double dis_step, double angle_cut, double dis_cut) const{

  bool flag = false;
  geo_point_t p_test;
  
  double min_dis = 1e9;
  double min_dis1 = 1e9;
  geo_point_t min_point = p_test1;
  
  for (int i=0; i!= int(test_dis/dis_step)+1;i++){
    p_test.set(p_test1.x() + dir.x() * i * dis_step,p_test1.y() + dir.y() * i * dis_step, p_test1.z() + dir.z() * i * dis_step);
    
    auto pts = get_closest_point_mcell(p_test);
    
    double dis = sqrt(pow(p_test.x() - pts.first.x(),2)+pow(p_test.y() - pts.first.y(),2)+pow(p_test.z() - pts.first.z(),2));
    double dis1 = sqrt(pow(p_test1.x() - pts.first.x(),2)+pow(p_test1.y() - pts.first.y(),2)+pow(p_test1.z() - pts.first.z(),2));
    if (dis < std::min(dis1 * tan(angle_cut/180.*3.1415926),dis_cut)){
      if (dis < min_dis){
	min_dis = dis;
	min_point = pts.first;
	min_dis1 = dis1;
      }
      if (dis < 3*units::cm)
	return std::make_pair(pts.first,dis1);
    }
  }

  return std::make_pair(min_point,min_dis1);
}

int Cluster::get_num_points() const{
  Scope scope = { "3d", {"x","y","z"} };
  const auto& sv = m_node->value.scoped_view(scope);       // get the kdtree
  // const auto& spcs = sv.pcs();
  // debug("sv {}", dump_pcs(sv.pcs()));
  const auto& skd = sv.kd();
  
  return skd.points().size();
}

std::pair<int, int> Cluster::get_num_points(const geo_point_t& point, const geo_point_t& dir) const{
  Scope scope = { "3d", {"x","y","z"} };
  const auto& sv = m_node->value.scoped_view(scope);       // get the kdtree
  // const auto& spcs = sv.pcs();
  // debug("sv {}", dump_pcs(sv.pcs()));
  const auto& skd = sv.kd();

  int num_p1 = 0;
  int num_p2 = 0;

  auto& points = skd.points();

  for (auto it = points.begin(); it!=points.end(); it++){
    geo_point_t dir1(it->at(0) - point.x(), it->at(1) - point.y(), it->at(2) - point.z());
    if (dir1.dot(dir)>=0){
      num_p1 ++;
    }else{
      num_p2 ++;
    }
  }

  return std::make_pair(num_p1, num_p2);
}

std::pair<int, int> Cluster::get_num_points(const geo_point_t& point, const geo_point_t& dir, double dis) const{
  Scope scope = { "3d", {"x","y","z"} };
  const auto& sv = m_node->value.scoped_view(scope);       // get the kdtree
  // const auto& spcs = sv.pcs();
  // debug("sv {}", dump_pcs(sv.pcs()));
  const auto& skd = sv.kd();

  int num_p1 = 0;
  int num_p2 = 0;

  auto rad = skd.radius(pow(dis,2), point);
  for (size_t pt_ind = 0; pt_ind<rad.size(); ++pt_ind) {
    auto& [pit,dist2] = rad[pt_ind];

    geo_point_t dir1(pit->at(0)-point.x(), pit->at(1)-point.y(), pit->at(2)-point.z());

    if (dir1.dot(dir)>=0){
      num_p1 ++;
    }else{
      num_p2 ++;
    }
    
  }


  return std::make_pair(num_p1, num_p2);
}



int Cluster::get_num_points(const geo_point_t& point, double dis) const{
  //return point_cloud->get_closest_points(p_test, dis).size();

  Scope scope = { "3d", {"x","y","z"} };
  const auto& sv = m_node->value.scoped_view(scope);       // get the kdtree
  // const auto& spcs = sv.pcs();
  // debug("sv {}", dump_pcs(sv.pcs()));
  const auto& skd = sv.kd();

  // following the definition in https://github.com/BNLIF/wire-cell-data/blob/5c9fbc4aef81c32b686f7c2dc7b0b9f4593f5f9d/src/ToyPointCloud.cxx#L656C10-L656C30
  auto rad = skd.radius(pow(dis,2), point);
  return rad.size();
}

std::map<std::shared_ptr<const WireCell::PointCloud::Facade::Blob>, geo_point_t> Cluster::get_closest_mcell(const geo_point_t& p, double search_radius) const{
  Scope scope = { "3d", {"x","y","z"} };
  const auto& sv = m_node->value.scoped_view(scope);       // get the kdtree
  // const auto& spcs = sv.pcs();
  // debug("sv {}", dump_pcs(sv.pcs()));
  const auto& skd = sv.kd();

  // following the definition in https://github.com/BNLIF/wire-cell-data/blob/5c9fbc4aef81c32b686f7c2dc7b0b9f4593f5f9d/src/ToyPointCloud.cxx#L656C10-L656C30
  auto rad = skd.radius(pow(search_radius,2), p);

  std::map<std::shared_ptr<const WireCell::PointCloud::Facade::Blob>, geo_point_t> mcell_point_map;
  std::map<std::shared_ptr<const WireCell::PointCloud::Facade::Blob>, double> mcell_dis_map;
  
  for (size_t pt_ind = 0; pt_ind<rad.size(); ++pt_ind) {
    auto& [pit,dist2] = rad[pt_ind];                    // what is the pit (point?)
    const auto [maj_ind,min_ind] = pit.index();        // maj_ind --> section, min_ind (within a section, what is the index)
    //  maj_inds.insert(maj_ind);

    geo_point_t p1(pit->at(0), pit->at(1), pit->at(2));
    const auto blob = m_blobs[maj_ind];    // this must be the blob ...
    //  auto charge = blob->charge;
    // set a minimal charge
    // following: https://github.com/BNLIF/wire-cell-data/blob/5c9fbc4aef81c32b686f7c2dc7b0b9f4593f5f9d/inc/WCPData/SlimMergeGeomCell.h#L59
    //    if (charge == 0) charge = 1;

    if (mcell_dis_map.find(blob) == mcell_dis_map.end()){
      mcell_dis_map[blob] = dist2;
      mcell_point_map[blob] = p1;
    }else{
      if (dist2 < mcell_dis_map[blob]){
	mcell_dis_map[blob] = dist2;
	mcell_point_map[blob] = p1;
      }
    }
    
  }

  
  return mcell_point_map;
}

std::pair<geo_point_t, std::shared_ptr<const WireCell::PointCloud::Facade::Blob> > Cluster::get_closest_point_mcell(const geo_point_t& origin) const{
  
  Scope scope = { "3d", {"x","y","z"} };
  const auto& sv = m_node->value.scoped_view(scope);       // get the kdtree
  // const auto& spcs = sv.pcs();
  // debug("sv {}", dump_pcs(sv.pcs()));
  const auto& skd = sv.kd();
  auto rad = skd.knn(1, origin);

  geo_point_t ret(0,0,0);
  std::shared_ptr<const WireCell::PointCloud::Facade::Blob> blob = 0;
  
  if (rad.size()==0)
    return std::make_pair(ret,blob);

  //  const auto& snodes = sv.nodes();
  auto& [pit,dist] = rad[0];                    // what is the pit (point?)
  const auto [maj_ind,min_ind] = pit.index();        // maj_ind --> section, min_ind (within a section, what is the index)

  ret.set( pit->at(0), pit->at(1), pit->at(2));
  
  blob = m_blobs[maj_ind];    // this must be the blob ...

  return std::make_pair(ret,blob);
}

geo_point_t Cluster::calc_ave_pos(const geo_point_t& origin, const double dis, const int alg) const {
    spdlog::set_level(spdlog::level::debug); // Set global log level to debug
    /// FIXME: there are many assumptions made, shoud we check these assumptions?
    /// a bit worriying about the speed.
    //    Scope scope = { "3d", {"x","y","z"} };
    //const auto& sv = m_node->value.scoped_view(scope);       // get the kdtree
    // const auto& spcs = sv.pcs();
    // debug("sv {}", dump_pcs(sv.pcs()));
    //const auto& skd = sv.kd();

    // following the definition in https://github.com/BNLIF/wire-cell-data/blob/5c9fbc4aef81c32b686f7c2dc7b0b9f4593f5f9d/src/ToyPointCloud.cxx#L656C10-L656C30

    //    auto rad = skd.radius(dis*dis, origin);                     // return is vector of (pointer, distance)
    //auto rad = skd.radius(100*units::m, origin);                     // return is vector of (pointer, distance)
    /// FIXME: what if rad is empty?
    //if(rad.size() == 0) {
        // raise<ValueError>("empty point cloud");
    //    return {0,0,0};
    // }
    //    const auto& snodes = sv.nodes();


    
    
    std::map<std::shared_ptr<const WireCell::PointCloud::Facade::Blob>, geo_point_t> pts = get_closest_mcell(origin, dis);
    
    // average position
    geo_point_t pt(0,0,0);
    double charge = 0;
    // alg following https://github.com/BNLIF/wire-cell-data/blob/5c9fbc4aef81c32b686f7c2dc7b0b9f4593f5f9d/src/PR3DCluster.cxx#L3956

    // // hack
    //geo_point_t origin1(2129.94, 709.449, 1525.01);
    //origin1 = origin1 - origin;
    //double dis1 = origin1.magnitude();

    for (auto it = pts.begin(); it!=pts.end(); it++){
      auto& blob = it->first;
      double q = blob->charge;
      if (q==0) q=1;
      pt += blob->center_pos()*q;
      charge += q;

      //hack ...
      //      if(dis1<1.0*units::mm) std::cout << origin << " " << blob->center_pos() << " " << q << " " << pts.size() << std::endl;
    }
    
    //std::set<size_t> maj_inds;                           //set, no duplications ...
    //    for (size_t pt_ind = 0; pt_ind<rad.size(); ++pt_ind) {
    //  auto& [pit,dist2] = rad[pt_ind];                    // what is the pit (point?)
    //  const auto [maj_ind,min_ind] = pit.index();        // maj_ind --> section, min_ind (within a section, what is the index)
      //  maj_inds.insert(maj_ind);

    //  const auto blob = m_blobs[maj_ind];    // this must be the blob ...
    //  auto charge = blob->charge;

      // set a minimal charge
      // following: https://github.com/BNLIF/wire-cell-data/blob/5c9fbc4aef81c32b686f7c2dc7b0b9f4593f5f9d/inc/WCPData/SlimMergeGeomCell.h#L59
    //  if (charge == 0) charge = 1;

     
      
    // ret += blob->center_pos() * charge;
    //  total_charge += charge;
      
    //}

    

    // if (dis1 < 0.1*units::mm){
    //   const auto &spcs = sv.pcs();
    //   std::vector<float_t> x;
    //   std::vector<float_t> y;
    //   std::vector<float_t> z;
    //   for(const auto& spc : spcs) {   // each little 3D pc --> (blobs)   spc represents x,y,z in a blob
    // 	const auto& x_ = spc.get().get("x")->elements<float_t>();
    // 	const auto& y_ = spc.get().get("y")->elements<float_t>();
    // 	const auto& z_ = spc.get().get("z")->elements<float_t>();
    // 	const size_t n = x_.size();
	
    // 	x.insert(x.end(), x_.begin(), x_.end()); // Append x_ to x
    // 	y.insert(y.end(), y_.begin(), y_.end());
    // 	z.insert(z.end(), z_.begin(), z_.end());
    //   }
    //   for (size_t i=0;i!=x.size();i++){
    // 	std::cout << "all " <<  i << " " << x.at(i) << " " << y.at(i) << " " << z.at(i) << std::endl;
    //   }
      

    //   for (size_t pt_ind = 0; pt_ind<rad.size(); ++pt_ind) {
    // 	auto& [pit,dist] = rad[pt_ind];
    // 	std::cout << "kd " << pt_ind << " " << pit->at(0) << " " << pit->at(1) << " " << pit->at(2) << " " << dist << std::endl;
    //   }
      
    //   //  for (size_t i=0;i!=skd.points().size();i++){
    //   //	std::cout << i << " " << skd.points()
    // }
    

    // this algorithm was not correctly translated !!!
    
    // debug("maj_inds.size() {} ", maj_inds.size());

    /*
    for (size_t maj_ind : maj_inds) {
        if (alg == 0) {
            const auto* node = snodes[maj_ind];
            const auto& lpcs = node->value.local_pcs();
            const auto& pc_scalar = lpcs.at("scalar");
            const auto charge = pc_scalar.get("charge")->elements<float_t>()[0];   // is this the blob?
            const auto center_x = pc_scalar.get("center_x")->elements<float_t>()[0];
            const auto center_y = pc_scalar.get("center_y")->elements<float_t>()[0];
            const auto center_z = pc_scalar.get("center_z")->elements<float_t>()[0];
            // debug("charge {} center {{{} {} {}}}", charge, center_x, center_y, center_z);
            geo_point_t inc(center_x, center_y, center_z);
            inc = inc * charge;
            ret += inc;
            total_charge += charge;
        } else {
	  const auto blob = m_blobs[maj_ind];    // this must be the blob ...
	  const auto charge = blob->charge;


	  


	  ret += blob->center_pos() * charge;
	  total_charge += charge;
        }
    }
    */
   

    
    if (charge != 0) {
        pt = pt / charge;
    }
    // debug("ret {{{} {} {}}}", ret.x(), ret.y(), ret.z());
    return pt;
}

#include <boost/histogram.hpp>
#include <boost/histogram/algorithm/sum.hpp>
namespace bh = boost::histogram;
namespace bha = boost::histogram::algorithm;


std::pair<double, double> Cluster::hough_transform(const geo_point_t& origin, const double dis, const int alg) const {
    Scope scope = { "3d", {"x","y","z"} };
    const auto& sv = m_node->value.scoped_view(scope);
    const auto& skd = sv.kd();
    auto rad = skd.radius(dis*dis, origin);
    /// FIXME: what if rad is empty?
    if(rad.size() == 0) {
        // raise<ValueError>("empty point cloud");
        return {0,0};
    }
    // const auto& spc = sv.pcs();

    const double pi = 3.141592653589793;
    // axes
    const Vector X(1,0,0);
    const Vector Y(0,1,0);
    const Vector Z(0,0,1);

    if (alg == 0) {
        auto hist = bh::make_histogram(bh::axis::regular<>( 180, -1.0, 1.0 ),
                                       bh::axis::regular<>( 360, -pi, pi ) );

        for (size_t pt_ind = 0; pt_ind<rad.size(); ++pt_ind) {
            auto& [pit,dist2] = rad[pt_ind];

	    // get average charge information ...
	    const auto [maj_ind,min_ind] = pit.index();        // maj_ind --> section, min_ind (within a section, what is the index)
	    const auto blob = m_blobs[maj_ind];    // this must be the blob ...
            auto charge = blob->charge;
	    // protection against the charge=0 case ...
	    if (charge ==0) charge = 1;
	    const auto npoints = blob->num_points();
            // debug("pt {{{} {} {}}}", pit->at(0), pit->at(1), pit->at(2));
            // auto pt = *pit;
            // debug("pt {{{} {} {}}}", pt[0], pt[1], pt[2]);
	    
            const geo_point_t pt(pit->at(0), pit->at(1), pit->at(2));
            Vector dir = (pt-origin).norm();
            const double cth = Z.dot(dir);
            const double phi = atan2(Y.dot(dir), X.dot(dir));
            hist(cth, phi, bh::weight(charge/npoints));
        }
        
        auto indexed = bh::indexed(hist);
        auto it = std::max_element(indexed.begin(), indexed.end());
        const auto& cell = *it;
        // std::stringstream ss;
        // ss << " maximum: index=[" << cell.index(0) <<","<< cell.index(1) <<"]"
        //    << " cth:[" << cell.bin(0).lower() << "," << cell.bin(0).upper() << "]"
        //    << " phi:[" << cell.bin(1).lower() << "," << cell.bin(1).upper() << "]"
        //    << " value=" << *cell
        //    << " sum=" << bha::sum(hist, bh::coverage::all);
        // spdlog::debug(ss.str());

        // cos(theta), phi
        return {cell.bin(0).center(), cell.bin(1).center()};
    }

    if (alg == 1) {
        auto hist = bh::make_histogram(bh::axis::regular<>( 180, 0, pi ),
                                       bh::axis::regular<>( 360, -pi, pi ) );

        for (size_t pt_ind = 0; pt_ind<rad.size(); ++pt_ind) {
            auto& [pit,dist2] = rad[pt_ind];

	    // get average charge information
	    const auto [maj_ind,min_ind] = pit.index();        // maj_ind --> section, min_ind (within a section, what is the index)
	    const auto blob = m_blobs[maj_ind];    // this must be the blob ...
	    auto charge = blob->charge;
	    // protection against charge = 0 case ...
	    if (charge ==0) charge = 1;
	    const auto npoints = blob->num_points();

	    if (charge <=0) continue;
	    
            const geo_point_t pt(pit->at(0), pit->at(1), pit->at(2));
            Vector dir = (pt-origin).norm();
            const double th = acos(Z.dot(dir));
            const double phi = atan2(Y.dot(dir), X.dot(dir));
            hist(th, phi, bh::weight(charge/npoints));
        }
        auto indexed = bh::indexed(hist);
        auto it = std::max_element(indexed.begin(), indexed.end());
        const auto& cell = *it;
        return {cell.bin(0).center(), cell.bin(1).center()};
    }
    raise<ValueError>("unknown alg %d", alg);
}

geo_point_t Cluster::vhough_transform(const geo_point_t& origin, const double dis, const int alg) const {
    if (alg == 0) {
        const auto [cth, phi] = hough_transform(origin, dis, alg);
        const double sth = sqrt(1-cth*cth);
        return {sth*cos(phi), sth*sin(phi), cth};
    }
    if (alg == 1) {
        const auto [th, phi] = hough_transform(origin, dis, alg);
        return {sin(th)*cos(phi), sin(th)*sin(phi), cos(th)};
    }
    raise<ValueError>("unknown alg %d", alg);
}

std::tuple<int, int, int, int> Cluster::get_uvwt_range() const {
    std::set<int> u_set;
    std::set<int> v_set;
    std::set<int> w_set;
    std::set<int> t_set;
    for (const auto& blob : m_blobs) {
        for(int i = blob->u_wire_index_min; i < blob->u_wire_index_max; ++i) {
            u_set.insert(i);
        }
        for(int i = blob->v_wire_index_min; i < blob->v_wire_index_max; ++i) {
            v_set.insert(i);
        }
        for(int i = blob->w_wire_index_min; i < blob->w_wire_index_max; ++i) {
            w_set.insert(i);
        }
        for(int i = blob->slice_index_min; i < blob->slice_index_max; ++i) {
            t_set.insert(i);
        }
    }
    return {u_set.size(), v_set.size(), w_set.size(), t_set.size()};
}

double Cluster::get_length(const TPCParams& tp) const {
    const auto [u, v, w, t] = get_uvwt_range();

    // t is in tick
    double length = std::sqrt(2./3.*(u*u*tp.pitch_u*tp.pitch_u + v*v*tp.pitch_v*tp.pitch_v + w*w*tp.pitch_w*tp.pitch_w) + t*t*tp.tick_width*tp.tick_width);

    //    if (length > 100*units::cm)
    // debug("u {} v {} w {} t {} length {}", u, v, w, t, length/units::cm);
    
    return length;
}

std::pair<geo_point_t, geo_point_t> Cluster::get_highest_lowest_points() const{
   // how to get all the points???
  Scope scope = { "3d", {"x","y","z"} };
  const auto& sv = m_node->value.scoped_view(scope);
  const auto& skd = sv.kd();
  
  auto& points = skd.points();

  geo_point_t highest_point(points.at(0).at(0),points.at(0).at(1),points.at(0).at(2));
  geo_point_t lowest_point = highest_point;

  for (auto it = points.begin(); it!=points.end(); it++){
    if (it->at(1) > highest_point.y()) // find highest Y ...
      highest_point.set(it->at(0), it->at(1), it->at(2));
    if (it->at(1) < lowest_point.y()) // find lowest Y point ...
      lowest_point.set(it->at(0), it->at(1), it->at(2));
  }

  return std::make_pair(highest_point, lowest_point);
}


std::pair<geo_point_t, geo_point_t> Cluster::get_earliest_latest_points() const{
  // how to get all the points???
  Scope scope = { "3d", {"x","y","z"} };
  const auto& sv = m_node->value.scoped_view(scope);
  const auto& skd = sv.kd();
  
  auto& points = skd.points();

  geo_point_t highest_point(points.at(0).at(0),points.at(0).at(1),points.at(0).at(2));
  geo_point_t lowest_point = highest_point;

  for (auto it = points.begin(); it!=points.end(); it++){
    if (it->at(0) > highest_point.x()) // find highest X ...
      highest_point.set(it->at(0), it->at(1), it->at(2));
    if (it->at(0) < lowest_point.x()) // find lowest X point ...
      lowest_point.set(it->at(0), it->at(1), it->at(2));
  }

  return std::make_pair(lowest_point, highest_point);
}
