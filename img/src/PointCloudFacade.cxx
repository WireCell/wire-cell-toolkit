#include "WireCellImg/PointCloudFacade.h"

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Facade;
// using WireCell::PointCloud::Dataset;
using namespace WireCell::PointCloud::Tree; // for "Points" node value type
// using WireCell::PointCloud::Tree::named_pointclouds_t;

#include "WireCellUtil/Logging.h"
using spdlog::debug;

/// unused
#if 0
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
#endif

Blob::Blob(node_t* n)
    : m_node(n)
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


Cluster::Cluster(node_t*n)
    : m_node(n)
{
    // build blobs
    for (auto child : m_node->children()) {
        auto blob = std::make_shared<Blob>(child);

        m_blobs.push_back(blob);
        m_time_blob_map.insert({blob->slice_index_min, blob});
    }
}

Blob::const_vector Cluster::is_connected(const Cluster& c, const int offset) const
{

    Blob::const_vector ret;
    for (const auto& [badtime, badblob] : c.m_time_blob_map) {
        auto bad_start = badtime;
        auto bad_end = badblob->slice_index_max; // not inclusive
        for (const auto& [good_start, goodblob] : m_time_blob_map) {
            auto good_end = goodblob->slice_index_max; // not inclusive
            if (good_end < bad_start || good_start >= bad_end) {
                continue;
            }
            if (goodblob->overlap_fast(*badblob, offset)) {
                ret.push_back(goodblob);
            }
        }
    }
    return ret;
}

Blob::const_pointer Cluster::get_first_blob() const{
    return m_time_blob_map.begin()->second;
}
Blob::const_pointer Cluster::get_last_blob() const{
    return m_time_blob_map.rbegin()->second;
}


std::pair<geo_point_t, double> Cluster::get_closest_point_along_vec(geo_point_t& p_test1, geo_point_t dir, double test_dis, double dis_step, double angle_cut, double dis_cut) const{

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
  const auto& sv = m_node->value.scoped_view(scope);       // get the kdtree
  const auto& skd = sv.kd();
  return skd.npoints();
}


int Cluster::get_num_points(const geo_point_t& point, double dis) const
{
    const auto& sv = m_node->value.scoped_view(scope);       // get the kdtree
    const auto& skd = sv.kd();

    // following the definition in https://github.com/BNLIF/wire-cell-data/blob/5c9fbc4aef81c32b686f7c2dc7b0b9f4593f5f9d/src/ToyPointCloud.cxx#L656C10-L656C30
    auto rad = skd.radius(dis*dis, point);
    return rad.size();
}


std::pair<int, int> Cluster::get_num_points(const geo_point_t& point, const geo_point_t& dir) const
{
    const auto& sv = m_node->value.scoped_view(scope);       // get the kdtree
    const auto& skd = sv.kd();
    const auto& points = skd.points();
    const size_t npoints = points.size();

    int num_p1 = 0;
    int num_p2 = 0;

    for (size_t ind=0; ind<npoints; ++ind) {
        geo_point_t dir1(points[0][ind] - point.x(),
                         points[1][ind] - point.y(),
                         points[2][ind] - point.z());
        if (dir1.dot(dir)>=0){
            ++num_p1;
        }
        else{
            ++num_p2;
        }
    }  

    return std::make_pair(num_p1, num_p2);
}


std::pair<int, int> Cluster::get_num_points(const geo_point_t& point, const geo_point_t& dir, double dis) const
{
    const auto& sv = m_node->value.scoped_view(scope);       // get the kdtree
    const auto& skd = sv.kd();
    const auto& points = skd.points();

    int num_p1 = 0;
    int num_p2 = 0;

    auto rad = skd.radius(dis*dis, point);
    for (const auto& [index,_] : rad) {
      
        geo_point_t dir1(points[0][index] - point.x(),
                         points[1][index] - point.y(),
                         points[2][index] - point.z());

        if (dir1.dot(dir) >= 0) {
            ++num_p1;
        }
        else{
            ++num_p2;
        }

    }

    return std::make_pair(num_p1, num_p2);
}



std::map<Blob::const_pointer, geo_point_t> Cluster::get_closest_mcell(const geo_point_t& point, double radius) const
{
    const auto& sv = m_node->value.scoped_view(scope);       // get the kdtree
    const auto& skd = sv.kd();
    const auto& points = skd.points();

    // following the definition in https://github.com/BNLIF/wire-cell-data/blob/5c9fbc4aef81c32b686f7c2dc7b0b9f4593f5f9d/src/ToyPointCloud.cxx#L656C10-L656C30
    auto rad = skd.radius(radius*radius, point);

    std::map<Blob::const_pointer, geo_point_t> mcell_point_map; // return
    std::map<Blob::const_pointer, double> mcell_dis2_map;       // temp
  

    // fixme: we could reduce the number of calls to map lookups if this loop
    // ends up being a hot spot.
    for (const auto& [index, dist2] : rad) {

        geo_point_t p1(points[0][index],
                       points[1][index],
                       points[2][index]);

        const auto blob_index = skd.major_index(index);
        const auto blob = m_blobs[blob_index];    // this must be the blob ...

        auto it = mcell_dis2_map.find(blob);
        if (it == mcell_dis2_map.end()) { // first time
            mcell_dis2_map[blob] = dist2;
            mcell_point_map[blob] = p1;
        }
        else {
            if (dist2 < it->second) { // check for yet closer
                it->second = dist2;
                mcell_point_map[blob] = p1;
            }
        }
    }
  
    return mcell_point_map;
}

std::pair<geo_point_t, Blob::const_pointer > Cluster::get_closest_point_mcell(const geo_point_t& point) const
{
    const auto& sv = m_node->value.scoped_view(scope);       // get the kdtree
    const auto& skd = sv.kd();
    const auto& points = skd.points();

    auto rad = skd.knn(1, point);
    if (rad.size() == 0) {
        return std::make_pair(geo_point_t(), nullptr);
    }

    const auto& [index, dist2] = rad[0];

    const auto blob_index = skd.major_index(index);
    const auto blob = m_blobs[blob_index];

    geo_point_t closest(points[0][index],
                        points[1][index],
                        points[2][index]);

    return std::make_pair(closest, blob);
}

geo_point_t Cluster::calc_ave_pos(const geo_point_t& origin, const double dis, const int alg) const
{
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


    std::map<Blob::const_pointer, geo_point_t> pts = get_closest_mcell(origin, dis);
    
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


std::pair<double, double> Cluster::hough_transform(const geo_point_t& origin, const double dis, const int alg) const
{
    // Scope scope = { "3d", {"x","y","z"} };
    const auto& sv = m_node->value.scoped_view(scope);
    const auto& skd = sv.kd();
    const auto& points = skd.points();

    auto rad = skd.radius(dis*dis, origin);
    /// FIXME: what if rad is empty?
    if(rad.size() == 0) {
        // raise<ValueError>("empty point cloud");
        // (bv) why not raise?  we do below.
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

        for (const auto& [index, dist2] : rad) {

            const auto blob_index = skd.major_index(index);
            const auto blob = m_blobs[blob_index];    // this must be the blob ...
            auto charge = blob->charge;
	    // protection against the charge=0 case ...
	    if (charge ==0) charge = 1;
            // fixme: what about change < 1?

            // fixme: alg==1 has this line, why not here?
	    // if (charge <=0) continue;
            
	    const auto npoints = blob->num_points();

            geo_point_t pt(points[0][index],
                           points[1][index],
                           points[2][index]);

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

        for (const auto& [index, dist2] : rad) {

            const auto blob_index = skd.major_index(index);
            const auto blob = m_blobs[blob_index];    // this must be the blob ...
            auto charge = blob->charge;
	    // protection against the charge=0 case ...
	    if (charge ==0) charge = 1;
            // fixme: what about change < 1?

	    const auto npoints = blob->num_points();

	    if (charge <=0) continue;
	    
            geo_point_t pt(points[0][index],
                           points[1][index],
                           points[2][index]);

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
    return std::pair<double, double>{}; // keep compiler happy
}


geo_point_t Cluster::vhough_transform(const geo_point_t& origin, const double dis, const int alg) const
{
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
    return geo_point_t{};       // keep compiler happy
}


std::tuple<int, int, int, int> Cluster::get_uvwt_range() const
{
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


double Cluster::get_length(const TPCParams& tp) const
{
    const auto [u, v, w, t] = get_uvwt_range();

    // t is in tick
    double length = std::sqrt(2./3.*(u*u*tp.pitch_u*tp.pitch_u + v*v*tp.pitch_v*tp.pitch_v + w*w*tp.pitch_w*tp.pitch_w) + t*t*tp.tick_width*tp.tick_width);

    //    if (length > 100*units::cm)
    // debug("u {} v {} w {} t {} length {}", u, v, w, t, length/units::cm);
    
    return length;
}

std::pair<geo_point_t, geo_point_t> Cluster::get_highest_lowest_points(size_t axis) const
{
    const auto& sv = m_node->value.scoped_view(scope);
    const auto& skd = sv.kd();
    const auto& points = skd.points();
    const size_t npoints = points.size();

    geo_point_t lowest_point, highest_point;

    for (size_t ind=0; ind<npoints; ++ind) {
        geo_point_t pt(points[0][ind],
                       points[1][ind],
                       points[2][ind]);
        if (!ind) {
            lowest_point = highest_point = pt;
            continue;
        }
        if (pt[axis] > highest_point[axis]) {
            highest_point = pt;
        }
        if (pt[axis] < lowest_point[axis]) {
            lowest_point = pt;
        }
    }

    return std::make_pair(highest_point, lowest_point);
}


std::pair<geo_point_t, geo_point_t> Cluster::get_earliest_latest_points() const
{
    auto backwards = get_highest_lowest_points(0);
    return std::make_pair(backwards.second, backwards.first);
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
