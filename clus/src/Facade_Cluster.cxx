#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Grouping.h"

#include "WireCellUtil/Array.h"

#include <boost/container_hash/hash.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "WireCellUtil/Logging.h"

// The original developers do not care.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"


using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Facade;
// using WireCell::PointCloud::Dataset;
using namespace WireCell::PointCloud::Tree;  // for "Points" node value type
// using WireCell::PointCloud::Tree::named_pointclouds_t;


using spdlog::debug;

// #define __DEBUG__
#ifdef __DEBUG__
#define LogDebug(x) std::cout << "[yuhw]: " << __LINE__ << " : " << x << std::endl
#else
#define LogDebug(x)
#endif



std::ostream& Facade::operator<<(std::ostream& os, const Facade::Cluster& cluster)
{
    const auto uvwt_min = cluster.get_uvwt_min();
    const auto uvwt_max = cluster.get_uvwt_max();
    std::cout << "uvwt_min " << std::get<0>(uvwt_min) << " " << std::get<1>(uvwt_min) << " " << std::get<2>(uvwt_min) << " " << std::get<3>(uvwt_min) << std::endl;
    std::cout << "uvwt_max " << std::get<0>(uvwt_max) << " " << std::get<1>(uvwt_max) << " " << std::get<2>(uvwt_max) << " " << std::get<3>(uvwt_max) << std::endl;
    os << "<Cluster [" << (void*) cluster.hash() << "]:" << " npts=" << cluster.npoints()
       << " nblobs=" << cluster.nchildren() << ">"
       << " u " << std::get<0>(uvwt_min) << " " << std::get<0>(uvwt_max)
       << " v " << std::get<1>(uvwt_min) << " " << std::get<1>(uvwt_max)
       << " w " << std::get<2>(uvwt_min) << " " << std::get<2>(uvwt_max)
       << " t " << std::get<3>(uvwt_min) << " " << std::get<3>(uvwt_max);
    return os;
}

Grouping* Cluster::grouping() { return this->m_node->parent->value.template facade<Grouping>(); }
const Grouping* Cluster::grouping() const { return this->m_node->parent->value.template facade<Grouping>(); }



void Cluster::clear_cache() const {

    // For now, this facade does its own cache management but we forward-call
    // the Mixin just to be proper.  Since our ClusterCache is the null-struct,
    // this is in truth pointless.  
    this->Mixin<Cluster,ClusterCache>::clear_cache();

    // The reason to keep fine-grained cache management is not all cluster users
    // need all cached values and by putting them all in fill_cache() we'd spoil
    // fine-grained laziness.  The cost we pay is that every single cached data
    // element is a chance to introduce a cache bug.

    // Reset time-blob mapping
    m_time_blob_map.clear();
    
    // Reset blob-indices mapping
    m_map_mcell_indices.clear();
    
    // Reset hull data
    m_hull_points.clear();
    m_hull_calculated = false;
    
    // Reset length and point count
    m_length = 0;
    m_npoints = 0;
    
    // Reset PCA data
    m_pca_calculated = false;
    m_center = geo_point_t();
    for(int i = 0; i < 3; i++) {
        m_pca_axis[i] = geo_vector_t();
        m_pca_values[i] = 0;
    }
    
    // Reset graph and path finding data
    m_graph.reset();
    m_parents.clear();
    m_distances.clear();
    m_source_pt_index = -1;
    m_path_wcps.clear();
    m_path_mcells.clear();
}


void Cluster::fill_cache(ClusterCache& cache) const
{
    for (const Blob* blob : children()) {
        cache.wpids.push_back(blob->wpid());
    }
}


std::vector<WireCell::WirePlaneId> Cluster::wpids() const {
    return cache().wpids;
}

void Cluster::print_blobs_info() const{
    for (const Blob* blob : children()) {
        std::cout << "U: " << blob->u_wire_index_min() << " " << blob->u_wire_index_max() 
        << " V: " << blob->v_wire_index_min() << " " << blob->v_wire_index_max() 
        << " W: " << blob->w_wire_index_min() << " " << blob->w_wire_index_max() 
        << " T: " << blob->slice_index_min() << " " << blob->slice_index_max()
        << std::endl;


    }
}

std::string Cluster::dump() const{
    const auto [u_min, v_min, w_min, t_min] = get_uvwt_min();
    const auto [u_max, v_max, w_max, t_max] = get_uvwt_max();
    std::stringstream ss;
    ss << " blobs " << children().size() << " points " << npoints()
    << " [" << t_min << " " << t_max << "] " << children().size()
    << " uvw " << u_min << " " << u_max << " " << v_min << " " << v_max << " " << w_min << " " << w_max;
    return ss.str();
}

std::string Cluster::dump_graph() const{
    if (m_graph==nullptr){
        return "empty graph";
    }
    auto g = *m_graph;
    std::stringstream ss;

    ss << "MCUGraph:" << std::endl;
    ss << "Vertices: " << num_vertices(g) << std::endl;
    ss << "Edges: " << num_edges(g) << std::endl;

    ss << "Vertex Properties:" << std::endl;
    auto vrange = boost::vertices(g);
    for (auto vit = vrange.first; vit != vrange.second; ++vit) {
        auto v = *vit;
        ss << "Vertex " << v << ": Index = " << g[v].index << point3d(g[v].index) << std::endl;
    }

    ss << "Edge Properties:" << std::endl;
    auto erange = boost::edges(g);
    auto weightMap = get(boost::edge_weight, g);
    for (auto eit = erange.first; eit != erange.second; ++eit) {
        auto e = *eit;
        auto src = source(e, g);
        auto tgt = target(e, g);
        ss << "Edge " << e << " [ " << point3d(g[src].index) << ", " << point3d(g[tgt].index) << " ]" << ": Distance = " << get(weightMap, e) << std::endl;
    }
    return ss.str();
}

const Cluster::time_blob_map_t& Cluster::time_blob_map() const
{
    if (m_time_blob_map.empty()) {
        for (const Blob* blob : children()) {
            m_time_blob_map[blob->slice_index_min()].insert(blob);
        }
    }
    return m_time_blob_map;
}

geo_point_t Cluster::get_furthest_wcpoint(geo_point_t old_wcp, geo_point_t dir, const double step,
                                          const int allowed_nstep) const
{
    dir = dir.norm();
    geo_point_t test_point;
    bool flag_continue = true;
    geo_point_t orig_point(old_wcp.x(), old_wcp.y(), old_wcp.z());
    geo_point_t orig_dir = dir;
    orig_dir = orig_dir.norm();
    int counter = 0;
    geo_point_t drift_dir_abs(1, 0, 0);

    double old_dis = 15 * units::cm;

    while (flag_continue && counter < 400) {
        counter++;

        // first step
        test_point.set(old_wcp.x() + dir.x() * step, old_wcp.y() + dir.y() * step, old_wcp.z() + dir.z() * step);
        // geo_point_t new_wcp = point_cloud->get_closest_wcpoint(test_point);
        auto [new_wcp, new_wcp_blob] = get_closest_point_blob(test_point);

        geo_point_t dir1(new_wcp.x() - old_wcp.x(), new_wcp.y() - old_wcp.y(), new_wcp.z() - old_wcp.z());
        double dis = dir1.magnitude();                      // distance change
        double angle = dir1.angle(dir) / 3.1415926 * 180.;  // local angle change

        geo_point_t dir2(new_wcp.x() - orig_point.x(), new_wcp.y() - orig_point.y(),
                         new_wcp.z() - orig_point.z());  // start from the original point
        double dis1 = dir2.magnitude();
        double angle1 = dir2.angle(orig_dir) / 3.1415926 * 180.;

        geo_point_t dir3(old_wcp.x() - orig_point.x(), old_wcp.y() - orig_point.y(),
                         old_wcp.z() - orig_point.z());  // start from the original point

        bool flag_para = false;

        double angle_1 = fabs(dir1.angle(drift_dir_abs) - 3.1415926 / 2.) * 180. / 3.1415926;
        double angle_2 = fabs(dir.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180.;
        double angle_3 = fabs(dir2.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180.;
        double angle_4 = fabs(dir3.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180.;

        if (angle_1 < 5 && angle_2 < 5 || angle_3 < 2.5 && angle_4 < 2.5) flag_para = true;

        bool flag_forward = false;
        if (flag_para) {
            // parallel case
            if (angle < 60 && dis > 0.2 * units::cm && (angle < 45 || angle1 <= 5)) {
                flag_forward = true;
            }
        }
        else {
            // non-parallel case
            if ((angle < 25 || dis < 1.2 * units::cm && angle < 60) &&                   // loose cut
                (angle < 15 || dis * sin(angle / 180. * 3.1415926) < 1.2 * units::cm ||  // tight cut
                 (angle < 21 && angle1 <= 2) ||
                 (angle1 <= 3 || dis1 * sin(angle1 / 180. * 3.1415926) < 3.6 * units::cm) && dis1 < 50 * units::cm) &&
                dis > 0.2 * units::cm) {  // in case of good direction
                flag_forward = true;
            }
            else if ((angle_1 < 5 || angle_2 < 5) && (angle_1 + angle_2) < 15 && dis > 0.2 * units::cm &&
                     (angle < 60) && (angle < 45 || angle1 <= 5)) {
                flag_forward = true;
            }
        }

        if (flag_forward) {
            old_wcp = new_wcp;

            if (dis > 3 * units::cm) {
                if (flag_para) {
                    dir = dir * old_dis + dir1 +
                          orig_dir * 15 * units::cm;  // if parallel, taking into account original direction ...
                }
                else {
                    dir = dir * old_dis + dir1;
                }
                dir = dir.norm();
                old_dis = dis;  //(old_dis*old_dis+dis*dis)/(old_dis + dis);
            }
        }
        else {
            //  failure & update direction
            flag_continue = false;

            test_point.set(old_wcp.x(), old_wcp.y(), old_wcp.z());

            geo_point_t dir4;
            double eff_dis;
            if (flag_para) {
                dir4 = vhough_transform(test_point, 100 * units::cm);
                eff_dis = 5 * units::cm;
            }
            else {
                dir4 = vhough_transform(test_point, 30 * units::cm);
                eff_dis = 15 * units::cm;
            }
            dir4 = dir4.norm();
            if (dir4.angle(dir) > 3.1415926 / 2.) dir4 = dir4 * -1;

            if (flag_para) {
                dir = dir * old_dis + dir4 * eff_dis + orig_dir * 15 * units::cm;
                dir = dir.norm();
                old_dis = eff_dis;
            }
            else {
                //	non-parallel case
                if (dir4.angle(dir) < 25 / 180. * 3.1415926) {
                    dir = dir * old_dis + dir4 * eff_dis;
                    dir = dir.norm();
                    old_dis = eff_dis;
                }
            }

            // start jump gaps
            for (int i = 0; i != allowed_nstep * 5; i++) {
                test_point.set(old_wcp.x() + dir.x() * step * (1 + 1. / 5. * i), old_wcp.y() + dir.y() * step * (1 + 1. / 5. * i),
                               old_wcp.z() + dir.z() * step * (1 + 1. / 5. * i));
                new_wcp = get_closest_point_blob(test_point).first;
                double dis2 = sqrt(pow(new_wcp.x() - test_point.x(), 2) + pow(new_wcp.y() - test_point.y(), 2) +
                                   pow(new_wcp.z() - test_point.z(), 2));
                dir1.set(new_wcp.x() - old_wcp.x(), new_wcp.y() - old_wcp.y(), new_wcp.z() - old_wcp.z());
                dis = dir1.magnitude();
                angle = dir1.angle(dir) / 3.1415926 * 180.;

                dir2.set(new_wcp.x() - orig_point.x(), new_wcp.y() - orig_point.y(), new_wcp.z() - orig_point.z());
                dis1 = dir2.magnitude();
                angle1 = dir2.angle(orig_dir) / 3.1415926 * 180.;

                dir3.set(old_wcp.x() - orig_point.x(), old_wcp.y() - orig_point.y(), old_wcp.z() - orig_point.z());

                flag_para = false;
                double angle_1 = fabs(dir1.angle(drift_dir_abs) - 3.1415926 / 2.) * 180. / 3.1415926;
                double angle_2 = fabs(dir.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180.;
                double angle_3 = fabs(dir2.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180.;
                double angle_4 = fabs(dir3.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180.;

                if (angle_1 < 7.5 && angle_2 < 7.5 || angle_3 < 5 && angle_4 < 5 && (angle_1 < 12.5 && angle_2 < 12.5))
                    flag_para = true;

                flag_forward = false;
                if (dis2 < 0.75 * step / 5. && dis > 0.2 * units::cm) flag_forward = true;

                if (flag_para) {
                    if (dis > step * 0.8 && (angle < 45 || angle1 <= 5.5) && (angle < 60)) flag_forward = true;
                }
                else {
                    if (((angle < 20 && dis < 30 * units::cm || dis * sin(angle / 180. * 3.1415926) < 1.2 * units::cm ||
                          angle < 15 && dis < 45 * units::cm || angle < 10 ||
                          (angle <= 28 && angle_1 < 2 && dis < 10 * units::cm) || (angle <= 28 && angle1 <= 2)) ||
                         (angle1 <= 3 || dis1 * sin(angle1 / 180. * 3.1415926) < 6.0 * units::cm) &&
                             dis1 < 100 * units::cm) &&
                        dis > step * 0.8 && (angle < 30)) {
                        flag_forward = true;
                    }
                    else if ((angle_1 < 5 || angle_2 < 5) && (angle_1 + angle_2) < 15 && dis > step * 0.8 &&
                             (angle < 60) && (angle < 45 || angle1 <= 5)) {
                        flag_forward = true;
                    }
                }

                if (flag_forward) {
                    old_wcp = new_wcp;

                    if (dis > 3 * units::cm) {
                        if (flag_para) {
                            dir = dir * old_dis + dir1 + orig_dir * 15 * units::cm;
                        }
                        else {
                            dir = dir * old_dis + dir1;
                        }
                        dir = dir.norm();
                        old_dis = (old_dis * old_dis + dis * dis) / (old_dis + dis);
                        if (old_dis > 15 * units::cm) old_dis = 15 * units::cm;
                    }

                    flag_continue = true;
                    break;
                }
            }
        }
    }

    return old_wcp;
}

void Cluster::adjust_wcpoints_parallel(size_t& start_idx, size_t& end_idx) const
{
    const auto& winds = wire_indices();

    geo_point_t start_p = point3d(start_idx);
    geo_point_t end_p = point3d(end_idx);

    double low_x = start_p.x() - 1 * units::cm;
    if (end_p.x() - 1 * units::cm < low_x) low_x = end_p.x() - 1 * units::cm;
    double high_x = start_p.x() + 1 * units::cm;
    if (end_p.x() + 1 * units::cm > high_x) high_x = end_p.x() + 1 * units::cm;

    // assumes u, v, w
    size_t low_idxes[3] = {start_idx, start_idx, start_idx};
    size_t high_idxes[3] = {end_idx, end_idx, end_idx};

    for (int pt_idx = 0; pt_idx != npoints(); pt_idx++) {
        geo_point_t current = point3d(pt_idx);
        if (current.x() > high_x || current.x() < low_x) continue;
        for (size_t pind = 0; pind != 3; ++pind) {
            if (winds[pind][pt_idx] < winds[pind][low_idxes[pind]]) {
                low_idxes[pind] = pt_idx;
            }
            if (winds[pind][pt_idx] > winds[pind][high_idxes[pind]]) {
                high_idxes[pind] = pt_idx;
            }
        }
    }

    std::vector<size_t> indices, temp_indices;
    std::set<size_t> indices_set;
    geo_point_t test_p;

    bool flags[3] = {true, true, true};

    /// HAIWANG: keeping the WCP original ordering
    if (winds[0][high_idxes[0]] - winds[0][low_idxes[0]] < winds[1][high_idxes[1]] - winds[1][low_idxes[1]]) {
        if (winds[0][high_idxes[0]] - winds[0][low_idxes[0]] < winds[2][high_idxes[2]] - winds[2][low_idxes[2]]) {
            flags[0] = false;
        }
        else {
            flags[2] = false;
        }
    }
    else {
        if (winds[1][high_idxes[1]] - winds[1][low_idxes[1]] < winds[2][high_idxes[2]] - winds[2][low_idxes[2]]) {
            flags[1] = false;
        }
        else {
            flags[2] = false;
        }
    }

    for (size_t pind = 0; pind != 3; ++pind) {
        if (flags[pind]) {
            geo_point_t low_p = point3d(low_idxes[pind]);
            geo_point_t high_p = point3d(high_idxes[pind]);
            std::vector<geo_point_t> test_points = {low_p, high_p, start_p, end_p};
            for (const auto& test_point : test_points) {
                temp_indices = get_closest_2d_index(test_point, 0.5 * units::cm, pind);
                std::copy(temp_indices.begin(), temp_indices.end(), inserter(indices_set, indices_set.begin()));
            }
        }
    }

    std::copy(indices_set.begin(), indices_set.end(), std::back_inserter(indices));

    size_t new_start_idx = start_idx;
    size_t new_end_idx = end_idx;

    //  std::cout << start_p.index << " " << end_p.index << std::endl;
    double sum_value = 0;
    for (size_t i = 0; i != indices.size(); i++) {
        //  std::cout << indices.at(i) << std::endl;
        for (size_t j = i + 1; j != indices.size(); j++) {
            double value = pow(winds[0][indices.at(i)] - winds[0][indices.at(j)], 2) +
                           pow(winds[1][indices.at(i)] - winds[1][indices.at(j)], 2) +
                           pow(winds[2][indices.at(i)] - winds[2][indices.at(j)], 2);

            if (value > sum_value) {
                // old_dis = dis;
                if (point3d(indices.at(i)).y() > point3d(indices.at(j)).y()) {
                    new_start_idx = indices.at(i);
                    new_end_idx = indices.at(j);
                }
                else {
                    new_start_idx = indices.at(j);
                    new_end_idx = indices.at(i);
                }
                geo_point_t new_start_p = point3d(new_start_idx);
                geo_point_t new_end_p = point3d(new_end_idx);

                if (sqrt(pow(new_start_p.x() - start_p.x(), 2) + pow(new_start_p.y() - start_p.y(), 2) +
                         pow(new_start_p.z() - start_p.z(), 2)) < 30 * units::cm &&
                    sqrt(pow(new_end_p.x() - end_p.x(), 2) + pow(new_end_p.y() - end_p.y(), 2) +
                         pow(new_end_p.z() - end_p.z(), 2)) < 30 * units::cm) {
                    start_idx = new_start_idx;
                    start_p = new_start_p;
                    end_idx = new_end_idx;
                    end_p = new_end_p;
                    sum_value = value;
                }
            }
        }
    }
}

bool Cluster::construct_skeleton(const bool use_ctpc)
{
    if (m_path_wcps.size() > 0) return false;
    // Calc_PCA();

    // WCP::WCPointCloud<double>& cloud = point_cloud->get_cloud();
    // WCPointCloud<double>::WCPoint highest_wcp = cloud.pts[0];
    // WCPointCloud<double>::WCPoint lowest_wcp = cloud.pts[0];
    geo_point_t highest_wcp = point3d(0);
    geo_point_t lowest_wcp = point3d(0);
    int highest_index = 0;
    int lowest_index = 0;

    // geo_point_t main_dir(PCA_axis[0].x, PCA_axis[0].y, PCA_axis[0].z);
    // main_dir.SetMag(1);
    // TVector3 temp_pt(highest_wcp.x - center.x, highest_wcp.y - center.y, highest_wcp.z - center.z);
    // double highest_value = temp_pt.Dot(main_dir);
    // double lowest_value = highest_value;
    geo_point_t main_dir = get_pca_axis(0);
    main_dir = main_dir.norm();
    geo_point_t center = get_center();
    geo_point_t temp_pt(highest_wcp.x() - center.x(), highest_wcp.y() - center.y(), highest_wcp.z() - center.z());
    double highest_value = temp_pt.dot(main_dir);
    double lowest_value = highest_value;

    // for (size_t i = 1; i < cloud.pts.size(); i++) {
    //     temp_pt.SetXYZ(cloud.pts[i].x - center.x, cloud.pts[i].y - center.y, cloud.pts[i].z - center.z);
    //     double value = temp_pt.Dot(main_dir);
    //     if (value > highest_value) {
    //         highest_value = value;
    //         highest_wcp = cloud.pts[i];
    //     }
    //     else if (value < lowest_value) {
    //         lowest_value = value;
    //         lowest_wcp = cloud.pts[i];
    //     }
    // }
    for (int i = 1; i < npoints(); i++) {
        temp_pt.set(point3d(i).x() - center.x(), point3d(i).y() - center.y(), point3d(i).z() - center.z());
        double value = temp_pt.dot(main_dir);
        if (value > highest_value) {
            highest_value = value;
            highest_wcp = point3d(i);
            highest_index = i;
        }
        else if (value < lowest_value) {
            lowest_value = value;
            lowest_wcp = point3d(i);
            lowest_index = i;
        }
    }

    dijkstra_shortest_paths(highest_index, use_ctpc);
    cal_shortest_path(lowest_index);
    return true;
}

const Cluster::sv2d_t& Cluster::sv2d(const size_t plane, const WirePlaneId wpid) const
{
    if (wpid.layer()!=kAllLayers) {
        raise<RuntimeError>("Cluster::sv2d() wpid.layer() {} != kAllLayers");
    }
    const Tree::Scope scope = {"3d", {scope2ds_prefix[plane]+"_x", scope2ds_prefix[plane]+"_y"}, 0, wpid.name()};
    return m_node->value.scoped_view(scope,
        [&](const Points::node_t& node) {
            const auto& lpcs = node.value.local_pcs();
            const auto& it = lpcs.find("scalar");
            if (it == lpcs.end()) {
                return false;
            }
            const auto& pc = it->second;
            const auto& wpida = pc.get("wpid");
            const auto wpidv = wpida->elements<int>();
            if (wpidv[0] == wpid.ident()) {
                return true;
            }
            // std::cerr << "Cluster::sv2d() wpid mismatch: " << wpidv[0] << " != " << wpid.ident() << std::endl;
            return false;
        }
    );
}

const Cluster::kd2d_t& Cluster::kd2d(const size_t plane, const WirePlaneId wpid) const
{
    const auto& sv = sv2d(plane, wpid);
    return sv.kd();
}

std::vector<size_t> Cluster::get_closest_2d_index(const geo_point_t& p, const double search_radius, const int plane, const WirePlaneId wpid) const {

    // const auto& tp = grouping()->get_params();
    // double angle_uvw[3] = {tp.angle_u, tp.angle_v, tp.angle_w};
    auto angles = grouping()->wire_angles(wpid.apa(), wpid.face());
    double angle_uvw[3];
    angle_uvw[0] = std::get<0>(angles);
    angle_uvw[1] = std::get<1>(angles);
    angle_uvw[2] = std::get<2>(angles);
    double x = p.x();
    double y = cos(angle_uvw[plane]) * p.z() - sin(angle_uvw[plane]) * p.y();
    std::vector<float_t> query_pt = {x, y};
    const auto& skd = kd2d(plane, wpid);
    auto ret_matches = skd.radius(search_radius * search_radius, query_pt);

    std::vector<size_t> ret_index(ret_matches.size());
    for (size_t i = 0; i != ret_matches.size(); i++)
    {
        ret_index.at(i) = ret_matches.at(i).first;
    }

    return ret_index;
}

std::vector<const Blob*> Cluster::is_connected(const Cluster& c, const int offset) const
{
    std::vector<const Blob*> ret;
    for (const auto& [bad_start, badblobs] : c.time_blob_map()) {
        for (const auto* badblob : badblobs) {
            auto bad_end = badblob->slice_index_max();  // not inclusive
            for (const auto& [good_start, goodblobs] : time_blob_map()) {
                for (const auto* goodblob : goodblobs) {
                    auto good_end = goodblob->slice_index_max();  // not inclusive
                    if (good_end <= bad_start || good_start >= bad_end) {  
                        continue;
                    }
                    if (goodblob->overlap_fast(*badblob, offset)) {
                        ret.push_back(goodblob);
                    }
                }

            }
        }
    }
    return ret;
}

const Blob* Cluster::get_first_blob() const
{
    if (time_blob_map().empty()) {
        raise<ValueError>("empty cluster has no first blob");
    }
    return *(time_blob_map().begin()->second.begin());
}

const Blob* Cluster::get_last_blob() const
{
    if (time_blob_map().empty()) {
        raise<ValueError>("empty cluster has no last blob");
    }
    return *(time_blob_map().rbegin()->second.rbegin());
}

size_t Cluster::get_num_time_slices() const
{
    return time_blob_map().size();
}

std::pair<geo_point_t, double> Cluster::get_closest_point_along_vec(geo_point_t& p_test1, geo_point_t dir,
                                                                    double test_dis, double dis_step, double angle_cut,
                                                                    double dis_cut) const
{
    geo_point_t p_test;

    double min_dis = 1e9;
    double min_dis1 = 1e9;
    geo_point_t min_point = p_test1;

    for (int i = 0; i != int(test_dis / dis_step) + 1; i++) {
        p_test.set(p_test1.x() + dir.x() * i * dis_step, p_test1.y() + dir.y() * i * dis_step,
                   p_test1.z() + dir.z() * i * dis_step);

        auto pts = get_closest_point_blob(p_test);

        double dis = sqrt(pow(p_test.x() - pts.first.x(), 2) + pow(p_test.y() - pts.first.y(), 2) +
                          pow(p_test.z() - pts.first.z(), 2));
        double dis1 = sqrt(pow(p_test1.x() - pts.first.x(), 2) + pow(p_test1.y() - pts.first.y(), 2) +
                           pow(p_test1.z() - pts.first.z(), 2));

        if (dis < std::min(dis1 * tan(angle_cut / 180. * 3.1415926), dis_cut)) {
            if (dis < min_dis) {
                min_dis = dis;
                min_point = pts.first;
                min_dis1 = dis1;
            }
            if (dis < 3 * units::cm) return std::make_pair(pts.first, dis1);
        }
    }

    return std::make_pair(min_point, min_dis1);
}

const Cluster::sv3d_t& Cluster::sv3d() const { return m_node->value.scoped_view(scope_3d_raw); }

const Cluster::kd3d_t& Cluster::kd3d() const { return sv3d().kd(); }
const Cluster::kd3d_t& Cluster::kd() const { return kd3d(); }
geo_point_t Cluster::point3d(size_t point_index) const { return kd3d().point3d(point_index); }
geo_point_t Cluster::point(size_t point_index) const { return point3d(point_index); }

const Cluster::points_type& Cluster::points() const { return kd3d().points(); }
int Cluster::npoints() const
{
    if (!m_npoints) {
        const auto& sv = sv3d();
        m_npoints = sv.npoints();
    }
    return m_npoints;
}
size_t Cluster::nbpoints() const
{
    size_t ret = 0;
    for (const auto* blob : children()) {
        ret += blob->nbpoints();
    }
    return ret;
}

const Cluster::wire_indices_t& Cluster::wire_indices() const
{
    const auto& sv = m_node->value.scoped_view<int_t>(scope_wire_index);
    const auto& skd = sv.kd();
    const auto& points = skd.points();
    LogDebug("points size: " << points.size() << " points[0] size: " << points[0].size());
    return points;
}

int Cluster::nnearby(const geo_point_t& point, double radius) const
{
    auto res = kd_radius(radius, point);
    return res.size();
}

std::pair<int, int> Cluster::ndipole(const geo_point_t& point, const geo_point_t& dir, const double dis) const
{
    const auto& points = this->points();
    const size_t npoints = points[0].size();

    int num_p1 = 0;
    int num_p2 = 0;

    for (size_t ind = 0; ind < npoints; ++ind) {
        geo_point_t dir1(points[0][ind] - point.x(), points[1][ind] - point.y(), points[2][ind] - point.z());
        if (dis > 0 && dir1.magnitude() > dis) continue;
        if (dir1.dot(dir) >= 0) {
            ++num_p1;
        }
        else {
            ++num_p2;
        }
    }

    return std::make_pair(num_p1, num_p2);
}

// std::pair<int, int> Cluster::nprojection(const geo_point_t& point, const geo_point_t& dir, double dis) const
// {
//     const auto& sv = m_node->value.scoped_view(scope_3d_raw);       // get the kdtree
//     const auto& skd = sv.kd();
//     const auto& points = skd.points();

//     int num_p1 = 0;
//     int num_p2 = 0;

//     auto rad = skd.radius(dis*dis, point);
//     for (const auto& [index,_] : rad) {

//         geo_point_t dir1(points[0][index] - point.x(),
//                          points[1][index] - point.y(),
//                          points[2][index] - point.z());

//         if (dir1.dot(dir) >= 0) {
//             ++num_p1;
//         }
//         else{
//             ++num_p2;
//         }

//     }

//     return std::make_pair(num_p1, num_p2);
// }

Cluster::kd_results_t Cluster::kd_knn(int nn, const geo_point_t& query_point) const
{
    const auto& skd = kd3d();
    return skd.knn(nn, query_point);
}

Cluster::kd_results_t Cluster::kd_radius(double radius, const geo_point_t& query_point) const
{
    const auto& skd = kd3d();
    return skd.radius(radius * radius, query_point);
}

std::vector<geo_point_t> Cluster::kd_points(const Cluster::kd_results_t& res)
{
    return const_cast<const Cluster*>(this)->kd_points(res);
}
std::vector<geo_point_t> Cluster::kd_points(const Cluster::kd_results_t& res) const
{
    std::vector<geo_point_t> ret;
    const auto& points = this->points();
    for (const auto& [point_index, _] : res) {
        ret.emplace_back(points[0][point_index], points[1][point_index], points[2][point_index]);
    }
    return ret;
}

// can't const_cast a vector.
template <typename T>
std::vector<T*> mutify(const std::vector<const T*>& c)
{
    size_t n = c.size();
    std::vector<T*> ret(n);
    for (size_t ind = 0; ind < n; ++ind) {
        ret[ind] = const_cast<T*>(c[ind]);
    }
    return ret;
}

std::vector<Blob*> Cluster::kd_blobs() { return mutify(const_cast<const Cluster*>(this)->kd_blobs()); }
std::vector<const Blob*> Cluster::kd_blobs() const
{
    std::vector<const Blob*> ret;
    const auto& sv = sv3d();
    for (const auto* node : sv.nodes()) {
        ret.push_back(node->value.facade<Blob>());
    }
    return ret;
}

Blob* Cluster::blob_with_point(size_t point_index)
{
    return const_cast<Blob*>(const_cast<const Cluster*>(this)->blob_with_point(point_index));
}

const Blob* Cluster::blob_with_point(size_t point_index) const
{
    const auto& sv = sv3d();
    const auto* node = sv.node_with_point(point_index);
    return node->value.facade<Blob>();
}

std::vector<Blob*> Cluster::blobs_with_points(const kd_results_t& res)
{
    return mutify(const_cast<const Cluster*>(this)->blobs_with_points(res));
}
std::vector<const Blob*> Cluster::blobs_with_points(const kd_results_t& res) const
{
    const size_t npts = res.size();
    std::vector<const Blob*> ret(npts);
    const auto& sv = sv3d();

    for (size_t ind = 0; ind < npts; ++ind) {
        const size_t point_index = res[ind].first;
        const auto* node = sv.node_with_point(point_index);
        ret[ind] = node->value.facade<Blob>();
    }
    return ret;
}

std::map<const Blob*, geo_point_t> Cluster::get_closest_blob(const geo_point_t& point, double radius) const
{
    struct Best {
        size_t point_index;
        double metric;
    };
    std::unordered_map<size_t, Best> best_blob_point;

    const auto& kd = kd3d();
    auto results = kd.radius(radius * radius, point);
    for (const auto& [point_index, metric] : results) {
        const size_t major_index = kd.major_index(point_index);
        auto it = best_blob_point.find(major_index);
        if (it == best_blob_point.end()) {  // first time seen
            best_blob_point[major_index] = {point_index, metric};
            continue;
        }
        if (metric < it->second.metric) {
            it->second.point_index = point_index;
            it->second.metric = metric;
        }
    }
    std::map<const Blob*, geo_point_t> ret;
    for (const auto& [mi, bb] : best_blob_point) {
        ret[blob_with_point(bb.point_index)] = point3d(bb.point_index);
    }
    return ret;
}

std::map<const Blob*, geo_point_t> Cluster::get_closest_blob(const geo_point_t& point, int N) const 
{
    struct Best {
        size_t point_index;
        double metric;
    };
    std::unordered_map<size_t, Best> best_blob_point;

    const auto& kd = kd3d();
    auto results = kd.knn(N, point);
    for (const auto& [point_index, metric] : results) {
        const size_t major_index = kd.major_index(point_index);
        auto it = best_blob_point.find(major_index);
        if (it == best_blob_point.end()) {  // first time seen
            best_blob_point[major_index] = {point_index, metric};
            continue;
        }
        if (metric < it->second.metric) {
            it->second.point_index = point_index;
            it->second.metric = metric;
        }
    }

    std::map<const Blob*, geo_point_t> ret;
    for (const auto& [mi, bb] : best_blob_point) {
        ret[blob_with_point(bb.point_index)] = point3d(bb.point_index);
    }
    return ret;
}

std::pair<geo_point_t, const Blob*> Cluster::get_closest_point_blob(const geo_point_t& point) const
{
    auto results = kd_knn(1, point);
    if (results.size() == 0) {
        return std::make_pair(geo_point_t(), nullptr);
    }

    const auto& [point_index, _] = results[0];
    return std::make_pair(point3d(point_index), blob_with_point(point_index));
}

std::pair<size_t, geo_point_t> Cluster::get_closest_wcpoint(const geo_point_t& point) const
{
    auto results = kd_knn(1, point);
    if (results.size() == 0) {
        return std::make_pair(-1, nullptr);
    }

    const auto& [point_index, _] = results[0];
    return std::make_pair(point_index, point3d(point_index));
}

size_t Cluster::get_closest_point_index(const geo_point_t& point) const
{
    auto results = kd_knn(1, point);
    if (results.size() == 0) {
        raise<ValueError>("no points in cluster");
    }

    const auto& [point_index, _] = results[0];
    return point_index;
}

double Cluster::get_closest_dis(const geo_point_t& point) const
{
    auto results = kd_knn(1, point);
    if (results.size() == 0) {
        raise<ValueError>("no points in cluster");
    }

    const auto& [_, dis] = results[0];
    return sqrt(dis);
}

std::tuple<int, int, double> Cluster::get_closest_points(const Cluster& other) const{

    double min_dis = 1e9;
    int p1_save = 0, p2_save = 0;

    // Sample points from this cluster at regular intervals 
    // Using ~20 sample points as initial probes
    int stride = std::max(1, (int)(npoints() / 20)); 

    for(int i = 0; i < npoints(); i += stride) {
        // Get initial point from this cluster
        geo_point_t p1 = point3d(i);
        
        // Get K nearest neighbors from other cluster using its kd-tree
        auto knn_results = other.kd3d().knn(5, p1); // Get 5 nearest candidates
        
        // Refine search around these neighbors
        for(const auto& [idx2, dist2] : knn_results) {
            if(sqrt(dist2) > min_dis) continue; // Skip if already farther than best found

            int curr_idx1 = i;
            geo_point_t p2 = other.point3d(idx2);
            int curr_idx2 = idx2;
            // Local refinement by alternating closest point lookups
            // This is similar to the original algorithm's refinement
            // but starts from better initial positions
            int prev_idx1 = -1, prev_idx2 = -1;
            const int max_iterations = 3; // Limit refinement iterations
            int iter = 0;

            while(iter++ < max_iterations && 
                  (curr_idx1 != prev_idx1 || curr_idx2 != prev_idx2)) {
                prev_idx1 = curr_idx1;
                prev_idx2 = curr_idx2;

                // Alternating closest point refinement
                curr_idx2 = other.get_closest_point_index(p1);
                p2 = other.point3d(curr_idx2);
                curr_idx1 = get_closest_point_index(p2);
                p1 = point3d(curr_idx1);

                double dis = sqrt(pow(p1.x()-p2.x(),2) + 
                                pow(p1.y()-p2.y(),2) + 
                                pow(p1.z()-p2.z(),2));

                if(dis < min_dis) {
                    min_dis = dis;
                    p1_save = curr_idx1;
                    p2_save = curr_idx2;
                // Early termination if we find a very close pair
                    if(dis < 0.5*units::cm) { // Threshold can be adjusted
                        return std::make_tuple(p1_save, p2_save, min_dis);
                    }
                }
            }
        }
    }

    return std::make_tuple(p1_save, p2_save, min_dis);

    // int p1_index = 0;
    // int p2_index = 0;
    // geo_point_t p1 = point3d(p1_index);
    // geo_point_t p2 = other.point3d(p2_index);
    // int p1_save = 0;
    // int p2_save = 0;
    // double min_dis = 1e9;

    // int prev_index1 = -1;
    // int prev_index2 = -1;
    // while (p1_index != prev_index1 || p2_index != prev_index2) {
    //     prev_index1 = p1_index;
    //     prev_index2 = p2_index;
    //     p2_index = other.get_closest_point_index(p1);
    //     p2 = other.point3d(p2_index);
    //     p1_index = get_closest_point_index(p2);
    //     p1 = point3d(p1_index);
    // }
    // // std::cout << "get_closest_points: " << p1_index << " " << p2_index << std::endl;
    // double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
    // if (dis < min_dis) {
    //     min_dis = dis;
    //     p1_save = p1_index;
    //     p2_save = p2_index;
    // }

    // prev_index1 = -1;
    // prev_index2 = -1;
    // p1_index = npoints() - 1;
    // p2_index = 0;
    // p1 = point3d(p1_index);
    // p2 = other.point3d(p2_index);
    // while (p1_index != prev_index1 || p2_index != prev_index2) {
    //     prev_index1 = p1_index;
    //     prev_index2 = p2_index;
    //     p2_index = other.get_closest_point_index(p1);
    //     p2 = other.point3d(p2_index);
    //     p1_index = get_closest_point_index(p2);
    //     p1 = point3d(p1_index);
    // }
    // // std::cout << "get_closest_points: " << p1_index << " " << p2_index << std::endl;
    // dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
    // if (dis < min_dis) {
    //     min_dis = dis;
    //     p1_save = p1_index;
    //     p2_save = p2_index;
    // }

    // prev_index1 = -1;
    // prev_index2 = -1;
    // p1_index = 0;
    // p2_index = other.npoints() - 1;
    // p1 = point3d(p1_index);
    // p2 = other.point3d(p2_index);
    // while (p1_index != prev_index1 || p2_index != prev_index2) {
    //     prev_index1 = p1_index;
    //     prev_index2 = p2_index;
    //     p2_index = other.get_closest_point_index(p1);
    //     p2 = other.point3d(p2_index);
    //     p1_index = get_closest_point_index(p2);
    //     p1 = point3d(p1_index);
    // }
    // // std::cout << "get_closest_points: " << p1_index << " " << p2_index << std::endl;
    // dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
    // if (dis < min_dis) {
    //     min_dis = dis;
    //     p1_save = p1_index;
    //     p2_save = p2_index;
    // }

    // prev_index1 = -1;
    // prev_index2 = -1;
    // p1_index = npoints() - 1;
    // p2_index = other.npoints() - 1;
    // p1 = point3d(p1_index);
    // p2 = other.point3d(p2_index);
    // while (p1_index != prev_index1 || p2_index != prev_index2) {
    //     prev_index1 = p1_index;
    //     prev_index2 = p2_index;
    //     p2_index = other.get_closest_point_index(p1);
    //     p2 = other.point3d(p2_index);
    //     p1_index = get_closest_point_index(p2);
    //     p1 = point3d(p1_index);
    // }
    // // std::cout << "get_closest_points: " << p1_index << " " << p2_index << std::endl;
    // dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
    // if (dis < min_dis) {
    //     min_dis = dis;
    //     p1_save = p1_index;
    //     p2_save = p2_index;
    // }

    // return std::make_tuple(p1_save, p2_save, min_dis);
}

geo_point_t Cluster::calc_ave_pos(const geo_point_t& origin, const double dis) const
{
    // average position
    geo_point_t ret(0, 0, 0);
    double charge = 0;

    auto blob_pts = get_closest_blob(origin, dis);
    for (auto [blob, _] : blob_pts) {
        double q = blob->charge();
        if (q == 0) q = 1;
        ret += blob->center_pos() * q;
        charge += q;
    }

    if (charge != 0) {
        ret = ret / charge;
    }

    return ret;
}

geo_point_t Cluster::calc_ave_pos(const geo_point_t& origin, int N) const
{
    // average position
    geo_point_t ret(0, 0, 0);
    double charge = 0;

    auto blob_pts = get_closest_blob(origin, N);
    for (auto [blob, _] : blob_pts) {
        double q = blob->charge(); 
        if (q == 0) q = 1;  // protection against zero charge
        ret += blob->center_pos() * q;
        charge += q;
    }

    if (charge != 0) {
        ret = ret / charge;
    }

    return ret;
}

geo_vector_t Cluster::calc_dir(const geo_point_t& p_test, const geo_point_t& p, double dis) const
{
    // Initialize direction vector
    geo_vector_t dir(0, 0, 0);
    
    // Get nearby blobs using existing interface
    auto blob_pts = get_closest_blob(p, dis);
    
    // Calculate weighted direction
    for (const auto& [blob, _] : blob_pts) {
        const geo_point_t point = blob->center_pos();
        const double q = blob->charge();
        
        // Calculate direction vector from p_test to point
        geo_vector_t dir1(point.x() - p_test.x(),
                         point.y() - p_test.y(), 
                         point.z() - p_test.z());
                         
        // Add weighted contribution
        dir += dir1 * q;
    }

    // Normalize if non-zero
    if (dir.magnitude() != 0) {
        dir = dir.norm();
    }
    
    return dir;
}


#include <boost/histogram.hpp>
#include <boost/histogram/algorithm/sum.hpp>
namespace bh = boost::histogram;
namespace bha = boost::histogram::algorithm;

// Example parameter calculating functions used by directional hough
// transforms.
static double theta_angle(const Vector& dir)
{
    const Vector Z(0, 0, 1);
    return acos(Z.dot(dir));
}
static double theta_cosine(const Vector& dir)
{
    const Vector Z(0, 0, 1);
    return Z.dot(dir);
}
static double phi_angle(const Vector& dir)
{
    const Vector X(1, 0, 0);
    const Vector Y(0, 1, 0);
    return atan2(Y.dot(dir), X.dot(dir));
}

std::pair<double, double> Cluster::hough_transform(const geo_point_t& origin, const double dis,
                                                   HoughParamSpace param_space,
                                                   std::shared_ptr<const Simple3DPointCloud> s3dpc,
                                                   const std::vector<size_t>& global_indices) const
{
    std::vector<geo_point_t> pts;
    std::vector<const Blob*> blobs;

    if (s3dpc == nullptr) {
        auto results = kd_radius(dis, origin);
        if (results.size() == 0) {
            return {0, 0};
        }
        blobs = blobs_with_points(results);
        pts = kd_points(results);
    } else {
        if (s3dpc->get_num_points() != global_indices.size()) {
            raise<ValueError>("global indices size mismatch");
        }
        auto results = s3dpc->kd().radius(dis * dis, origin);
        for (const auto& [point_index, _] : results) {
            pts.push_back(s3dpc->point(point_index));
            size_t global_index = global_indices[point_index];
            blobs.push_back(blob_with_point(global_index));
        }
    }

    constexpr double pi = 3.141592653589793;

    using direction_parameter_function_f = std::function<double(const Vector& dir)>;

    // Parameter axis 1 is some measure of theta angle (angle or cosine)
    const int nbins1 = 180;
    // param_space == costh_phi
    direction_parameter_function_f theta_param = theta_cosine;
    double min1 = -1.0, max1 = 1.0;
    if (param_space == HoughParamSpace::theta_phi) {
        theta_param = theta_angle;
        min1 = 0;
        max1 = pi;
    }

    // Parameter axis 2 is only supported by phi angle
    const int nbins2 = 360;
    const double min2 = -pi;
    const double max2 = +pi;
    direction_parameter_function_f phi_param = phi_angle;

    auto hist = bh::make_histogram(bh::axis::regular<>(nbins1, min1, max1), bh::axis::regular<>(nbins2, min2, max2));

    for (size_t ind = 0; ind < blobs.size(); ++ind) {
        const auto* blob = blobs[ind];
        auto charge = blob->charge();
        // protection against the charge=0 case ...
        if (charge == 0) charge = 1;
        if (charge <= 0) continue;

        const auto npoints = blob->npoints();
        const auto& pt = pts[ind];

        const Vector dir = (pt - origin).norm();

        const double p1 = theta_param(dir);
        const double p2 = phi_param(dir);
        hist(p1, p2, bh::weight(charge / npoints));
    }

    auto indexed = bh::indexed(hist);
    auto it = std::max_element(indexed.begin(), indexed.end());
    const auto& cell = *it;
    return {cell.bin(0).center(), cell.bin(1).center()};
}

geo_point_t Cluster::vhough_transform(const geo_point_t& origin, const double dis, HoughParamSpace param_space,
                                      std::shared_ptr<const Simple3DPointCloud> s3dpc,
                                      const std::vector<size_t>& global_indices) const
{
    if (param_space == HoughParamSpace::theta_phi) {
        const auto [th, phi] = hough_transform(origin, dis, param_space, s3dpc, global_indices);
        return {sin(th) * cos(phi), sin(th) * sin(phi), cos(th)};
    }
    // costh_phi
    const auto [cth, phi] = hough_transform(origin, dis, param_space, s3dpc, global_indices);
    const double sth = sqrt(1 - cth * cth);
    return {sth * cos(phi), sth * sin(phi), cth};
}

std::tuple<int, int, int, int> Cluster::get_uvwt_min() const
{
    std::tuple<int, int, int, int> ret;
    bool first = true;

    for (const auto* blob : children()) {
        const int u = blob->u_wire_index_min();
        const int v = blob->v_wire_index_min();
        const int w = blob->w_wire_index_min();
        const int t = blob->slice_index_min();

        if (first) {
            ret = {u, v, w, t};
            first = false;
            continue;
        }
        get<0>(ret) = std::min(get<0>(ret), u);
        get<1>(ret) = std::min(get<1>(ret), v);
        get<2>(ret) = std::min(get<2>(ret), w);
        get<3>(ret) = std::min(get<3>(ret), t);
    }
    return ret;
}
std::tuple<int, int, int, int> Cluster::get_uvwt_max() const
{
    std::tuple<int, int, int, int> ret;
    bool first = true;

    for (const auto* blob : children()) {
        const int u = blob->u_wire_index_max();
        const int v = blob->v_wire_index_max();
        const int w = blob->w_wire_index_max();
        const int t = blob->slice_index_max();

        if (first) {
            ret = {u, v, w, t};
            first = false;
            continue;
        }
        get<0>(ret) = std::max(get<0>(ret), u);
        get<1>(ret) = std::max(get<1>(ret), v);
        get<2>(ret) = std::max(get<2>(ret), w);
        get<3>(ret) = std::max(get<3>(ret), t);
    }
    return ret;
}

// FIXME: Is this actually correct?  It does not return "ranges" but rather the
// number of unique wires/ticks in the cluster.  A sparse but large cluster will
// be "smaller" than a small but dense cluster.
std::map<WirePlaneId, std::tuple<int, int, int, int> > Cluster::get_uvwt_range() const
{
    std::map<WirePlaneId, std::set<int> > map_wpid_u_set;
    std::map<WirePlaneId, std::set<int> > map_wpid_v_set;
    std::map<WirePlaneId, std::set<int> > map_wpid_w_set;
    std::map<WirePlaneId, std::set<int> > map_wpid_t_set;
    for (const auto* blob : children()) {
        for (int i = blob->u_wire_index_min(); i < blob->u_wire_index_max(); ++i) {
            map_wpid_u_set[blob->wpid()].insert(i);
        }
        for (int i = blob->v_wire_index_min(); i < blob->v_wire_index_max(); ++i) {
            map_wpid_v_set[blob->wpid()].insert(i);
        }
        for (int i = blob->w_wire_index_min(); i < blob->w_wire_index_max(); ++i) {
            map_wpid_w_set[blob->wpid()].insert(i);
        }
        for (int i = blob->slice_index_min(); i < blob->slice_index_max(); ++i) {
            map_wpid_t_set[blob->wpid()].insert(i);
        }
    }
    std::map<WirePlaneId, std::tuple<int, int, int, int> > ret;
    for (auto it = map_wpid_u_set.begin(); it != map_wpid_u_set.end(); ++it) {
        const WirePlaneId wpid = it->first;
        const auto& u_set = it->second;
        const auto& v_set = map_wpid_v_set[wpid];
        const auto& w_set = map_wpid_w_set[wpid];
        const auto& t_set = map_wpid_t_set[wpid];
        ret[wpid] = {u_set.size(), v_set.size(), w_set.size(), t_set.size()};
    }
    return ret;
    // return {u_set.size(), v_set.size(), w_set.size(), t_set.size()};
}

double Cluster::get_length() const
{
    if (m_length == 0) {  // invalidates when a new node is set
        // const auto& tp = grouping()->get_params();

        // std::cout << "Test: " << grouping()->get_anode()->face(0)->plane(0)->pimpos()->pitch() << " " << grouping()->get_anode()->face(0)->plane(1)->pimpos()->pitch() << " " << grouping()->get_anode()->face(0)->plane(2)->pimpos()->pitch() << " " << tp.pitch_u << " " << tp.pitch_v << " " << tp.pitch_w << std::endl;

        auto map_wpid_uvwt = get_uvwt_range();
        for (const auto& [wpid, uvwt] : map_wpid_uvwt) {

            const double tick = grouping()->get_tick().at(wpid.apa()).at(wpid.face());
            const double drift_speed = grouping()->get_drift_speed().at(wpid.apa()).at(wpid.face());

            // std::cout << "Test: " << wpid.apa() << " " << wpid.face() << " " << tp.tick_drift << " " << tick * drift_speed << std::endl;

            const auto [u, v, w, t] = uvwt;
            const double pu = u * grouping()->get_anode(wpid.apa())->face(wpid.face())->plane(0)->pimpos()->pitch() ;
            const double pv = v * grouping()->get_anode(wpid.apa())->face(wpid.face())->plane(1)->pimpos()->pitch();
            const double pw = w * grouping()->get_anode(wpid.apa())->face(wpid.face())->plane(2)->pimpos()->pitch();
            const double pt = t * tick * drift_speed;
            m_length += std::sqrt(2. / 3. * (pu * pu + pv * pv + pw * pw) + pt * pt);
        }
    }
    return m_length;
}

std::map<WirePlaneId, std::tuple<int, int, int, int> > Facade::get_uvwt_range(const Cluster* cluster, const std::vector<int>& b2id, const int id)
{
    std::map<WirePlaneId, std::set<int> > map_wpid_u_set;
    std::map<WirePlaneId, std::set<int> > map_wpid_v_set;
    std::map<WirePlaneId, std::set<int> > map_wpid_w_set;
    std::map<WirePlaneId, std::set<int> > map_wpid_t_set;

    for (size_t i = 0; i != b2id.size(); i++) {
        if (b2id.at(i) != id) continue;
        const auto* blob = cluster->children().at(i);
        for (int i = blob->u_wire_index_min(); i < blob->u_wire_index_max(); ++i) {
            map_wpid_u_set[blob->wpid()].insert(i);
        }
        for (int i = blob->v_wire_index_min(); i < blob->v_wire_index_max(); ++i) {
            map_wpid_v_set[blob->wpid()].insert(i);
        }
        for (int i = blob->w_wire_index_min(); i < blob->w_wire_index_max(); ++i) {
            map_wpid_w_set[blob->wpid()].insert(i);
        }
        for (int i = blob->slice_index_min(); i < blob->slice_index_max(); ++i) {
            map_wpid_t_set[blob->wpid()].insert(i);
        }
    }

    std::map<WirePlaneId, std::tuple<int, int, int, int> > ret;
    for (auto it = map_wpid_u_set.begin(); it != map_wpid_u_set.end(); ++it) {
        const WirePlaneId wpid = it->first;
        const auto& u_set = it->second;
        const auto& v_set = map_wpid_v_set[wpid];
        const auto& w_set = map_wpid_w_set[wpid];
        const auto& t_set = map_wpid_t_set[wpid];
        ret[wpid] = {u_set.size(), v_set.size(), w_set.size(), t_set.size()};
    }
    return ret;

    // return {u_set.size(), v_set.size(), w_set.size(), t_set.size()};
}

double Facade::get_length(const Cluster* cluster, const std::vector<int>& b2id, const int id)
{
    // const auto& tp = cluster->grouping()->get_params();
    auto map_wpid_uvwt = Facade::get_uvwt_range(cluster, b2id, id);
    double length = 0;
    for (const auto& [wpid, uvwt] : map_wpid_uvwt) {
        const double tick = cluster->grouping()->get_tick().at(wpid.apa()).at(wpid.face());
        const double drift_speed = cluster->grouping()->get_drift_speed().at(wpid.apa()).at(wpid.face());

        const auto [u, v, w, t] = uvwt;
        const double pu = u * cluster->grouping()->get_anode(wpid.apa())->face(wpid.face())->plane(0)->pimpos()->pitch();
        const double pv = v * cluster->grouping()->get_anode(wpid.apa())->face(wpid.face())->plane(1)->pimpos()->pitch();
        const double pw = w * cluster->grouping()->get_anode(wpid.apa())->face(wpid.face())->plane(2)->pimpos()->pitch();
        const double pt = t * tick * drift_speed;
        length += std::sqrt(2. / 3. * (pu * pu + pv * pv + pw * pw) + pt * pt);
    }
    return length;
}


std::pair<geo_point_t, geo_point_t> Cluster::get_highest_lowest_points(size_t axis) const
{
    const auto& points = this->points();
    const size_t npoints = points[0].size();

    geo_point_t lowest_point, highest_point;

    for (size_t ind = 0; ind < npoints; ++ind) {
        geo_point_t pt(points[0][ind], points[1][ind], points[2][ind]);
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

std::pair<geo_point_t, geo_point_t> Cluster::get_front_back_points() const
{
    return get_highest_lowest_points(2);
}

std::pair<geo_point_t, geo_point_t> Cluster::get_main_axis_points() const
{
    // Get first point as initial values
    geo_point_t highest_point = point3d(0);
    geo_point_t lowest_point = point3d(0);

    // Get main axis and ensure consistent direction (y>0)
    geo_point_t main_axis = get_pca_axis(0); 
    if (main_axis.y() < 0) {
        main_axis = main_axis * -1;
    }

    // Initialize extreme values using projections of first point
    double high_value = highest_point.dot(main_axis);
    double low_value = high_value;

    // Loop through all points to find extremes along main axis
    for (int i = 1; i < npoints(); i++) {
        geo_point_t current = point3d(i);
        double value = current.dot(main_axis);
        
        if (value > high_value) {
            highest_point = current;
            high_value = value; 
        }
        if (value < low_value) {
            lowest_point = current;
            low_value = value;
        }
    }

    return std::make_pair(highest_point, lowest_point);
}

std::pair<geo_point_t,geo_point_t> Cluster::get_two_extreme_points() const
{
    geo_point_t extreme_wcp[6];
    for (int i = 0; i != 6; i++) {
        extreme_wcp[i] = point3d(0);
    }
    for (int i = 1; i < npoints(); i++) {
        if (point3d(i).y() > extreme_wcp[0].y()) extreme_wcp[0] = point3d(i);
        if (point3d(i).y() < extreme_wcp[1].y()) extreme_wcp[1] = point3d(i);

        if (point3d(i).x() > extreme_wcp[2].x()) extreme_wcp[2] = point3d(i);
        if (point3d(i).x() < extreme_wcp[3].x()) extreme_wcp[3] = point3d(i);

        if (point3d(i).z() > extreme_wcp[4].z()) extreme_wcp[4] = point3d(i);
        if (point3d(i).z() < extreme_wcp[5].z()) extreme_wcp[5] = point3d(i);
    }

    double max_dis = -1;
    geo_point_t wcp1, wcp2;
    for (int i = 0; i != 6; i++) {
        for (int j = i + 1; j != 6; j++) {
            double dis =
                sqrt(pow(extreme_wcp[i].x() - extreme_wcp[j].x(), 2) + pow(extreme_wcp[i].y() - extreme_wcp[j].y(), 2) +
                     pow(extreme_wcp[i].z() - extreme_wcp[j].z(), 2));
            if (dis > max_dis) {
                max_dis = dis;
                wcp1 = extreme_wcp[i];
                wcp2 = extreme_wcp[j];
            }
        }
    }
    geo_point_t p1(wcp1.x(), wcp1.y(), wcp1.z());
    geo_point_t p2(wcp2.x(), wcp2.y(), wcp2.z());
    p1 = calc_ave_pos(p1, 5 * units::cm);
    p2 = calc_ave_pos(p2, 5 * units::cm);

    return std::make_pair(p1, p2);
}

bool Cluster::sanity(Log::logptr_t log) const
{
    {
        const auto* svptr = m_node->value.get_scoped(scope_3d_raw);
        if (!svptr) {
            if (log) log->debug("cluster sanity: note, not yet a scoped view {}", scope_3d_raw);
        }
    }
    if (!nchildren()) {
        if (log) log->debug("cluster sanity: no children blobs");
        return false;
    }

    const auto& sv = m_node->value.scoped_view(scope_3d_raw);
    const auto& snodes = sv.nodes();
    if (snodes.empty()) {
        if (log) log->debug("cluster sanity: no scoped nodes");
        return false;
    }
    if (sv.npoints() == 0) {  // triggers a scoped view cache fill
        if (log) log->debug("cluster sanity: no scoped points");
        return false;
    }
    // sv.force_invalid();
    const auto& skd = sv.kd();

    const auto& fblobs = children();

    for (const Blob* blob : fblobs) {
        if (!blob->sanity(log)) return false;
    }

    if (skd.nblocks() != snodes.size()) {
        if (log) log->debug("cluster sanity: k-d blocks={} scoped nodes={}", skd.nblocks(), fblobs.size());
        return false;
    }

    if (skd.nblocks() != fblobs.size()) {
        if (log) log->debug("cluster sanity: k-d blocks={} cluster blobs={}", skd.nblocks(), fblobs.size());
        return false;
    }

    for (size_t ind = 0; ind < snodes.size(); ++ind) {
        /// In general, the depth-first order of scoped view nodes is always the
        /// same as k-d tree blocks but is not (again in general) expected to be
        /// the same order as the children blobs of a parent cluster.  After
        /// all, a scoped view may span multiple essentially any subset of tree
        /// nodes.  However, in the special case of the "3d" SV and how the PC
        /// tree is constructed, the depth-first and children blobs ordering
        /// should be "accidentally" the same.
        const auto* fblob = fblobs[ind];
        const auto* sblob = snodes[ind]->value.facade<Blob>();
        if (fblob != sblob) {
            if (log) {
                log->debug("cluster sanity: scoped node facade Blob differs from cluster child at {}", ind);
                log->debug("cluster sanity: \tscoped blob: {}", *fblob);
                log->debug("cluster sanity: \tfacade blob: {}", *sblob);
            }
            // return false;
        }
    }

    const auto& majs = skd.major_indices();
    const auto& mins = skd.minor_indices();
    const size_t npts = skd.npoints();

    const Blob* sblob = nullptr;
    std::vector<geo_point_t> spoints;

    for (size_t ind = 0; ind < npts; ++ind) {
        auto kdpt = skd.point3d(ind);

        const size_t majind = majs[ind];
        const size_t minind = mins[ind];

        // scoped consistency
        const node_t* tnode = sv.node_with_point(ind);
        if (!tnode) {
            if (log) log->debug("cluster sanity: scoped node facade not a Blob at majind={}", majind);
            return false;
        }
        const auto* tblob = tnode->value.facade<Blob>();
        if (tblob != sblob) {
            sblob = tblob;
            spoints = sblob->points();
        }

        if (minind >= spoints.size()) {
            if (log)
                log->debug("cluster sanity: minind={} is beyond scoped blob npts={} majind={}, blob is: {}", minind,
                           spoints.size(), majind, *sblob);
            return false;
        }
        auto spt = spoints[minind];
        if (spt != kdpt) {
            if (log)
                log->debug("cluster sanity: scoped point mismatch at minind={} majind={} spt={} kdpt={}, blob is: {}",
                           minind, majind, spt, kdpt, *sblob);
            return false;
        }
    }
    return true;
}

size_t Cluster::hash() const
{
    std::size_t h = 0;
    boost::hash_combine(h, (size_t) (get_length() / units::mm));
    auto blobs = children();  // copy vector
    sort_blobs(blobs);
    for (const Blob* blob : blobs) {
        boost::hash_combine(h, blob->hash());
    }
    return h;
}

std::vector<int> Cluster::get_blob_indices(const Blob* blob) const
{
    if (m_map_mcell_indices.empty()) {
        const auto& skd = kd3d();
        for (size_t ind = 0; ind < skd.npoints(); ++ind) {
            const auto* blob = blob_with_point(ind);
            m_map_mcell_indices[blob].push_back(ind);
        }
    }
    return m_map_mcell_indices[blob];
}

// #define LogDebug(x) std::cout << "[yuhw]: " << __LINE__ << " : " << x << std::endl
void Cluster::Create_graph(const bool use_ctpc) const
{
    // std::cout << "Create Graph!" << std::endl;
    LogDebug("Create Graph! " << graph);
    if (m_graph != nullptr) return;
    m_graph = std::make_unique<MCUGraph>(nbpoints());
    // std::cout << "Test:" << "Create Graph!" << std::endl;
    Establish_close_connected_graph();
    if (use_ctpc) Connect_graph(true);
    Connect_graph();
}

void Cluster::Establish_close_connected_graph() const
{
    std::map<const Blob*, std::map<int, std::set<int>>, blob_less_functor> map_mcell_uindex_wcps;
    std::map<const Blob*, std::map<int, std::set<int>>, blob_less_functor> map_mcell_vindex_wcps;
    std::map<const Blob*, std::map<int, std::set<int>>, blob_less_functor> map_mcell_windex_wcps;

    std::map<Blob*, std::set<int>, blob_less_functor> map_mcell_indices;

    const auto& points = this->points();
    const auto& winds = this->wire_indices();
    LogDebug("points[0].size(): " << points[0].size() << " winds[0].size(): " << winds[0].size());

    for (Blob* mcell : this->children()) {
        std::map<int, std::set<int>> map_uindex_wcps;
        std::map<int, std::set<int>> map_vindex_wcps;
        std::map<int, std::set<int>> map_windex_wcps;

        std::vector<int> pinds = this->get_blob_indices(mcell);
        for (const int pind : pinds) {
            auto v = vertex(pind, *m_graph);  // retrieve vertex descriptor
            (*m_graph)[v].index = pind;
            if (map_uindex_wcps.find(winds[0][pind]) == map_uindex_wcps.end()) {
                std::set<int> wcps;
                wcps.insert(pind);
                map_uindex_wcps[winds[0][pind]] = wcps;
            }
            else {
                map_uindex_wcps[winds[0][pind]].insert(pind);
            }

            if (map_vindex_wcps.find(winds[1][pind]) == map_vindex_wcps.end()) {
                std::set<int> wcps;
                wcps.insert(pind);
                map_vindex_wcps[winds[1][pind]] = wcps;
            }
            else {
                map_vindex_wcps[winds[1][pind]].insert(pind);
            }

            if (map_windex_wcps.find(winds[2][pind]) == map_windex_wcps.end()) {
                std::set<int> wcps;
                wcps.insert(pind);
                map_windex_wcps[winds[2][pind]] = wcps;
            }
            else {
                map_windex_wcps[winds[2][pind]].insert(pind);
            }
        }
        map_mcell_uindex_wcps[mcell] = map_uindex_wcps;
        map_mcell_vindex_wcps[mcell] = map_vindex_wcps;
        map_mcell_windex_wcps[mcell] = map_windex_wcps;
    }

    int num_edges = 0;

    // create graph for points inside the same mcell
    for (Blob* mcell : this->children()) {
        std::vector<int> pinds = this->get_blob_indices(mcell);
        int max_wire_interval = mcell->get_max_wire_interval();
        int min_wire_interval = mcell->get_min_wire_interval();
        // std::cout << "mcell: " << pinds.size()
        // << " type " << mcell->get_max_wire_type() << " " << mcell->get_min_wire_type()
        // << " interval " << max_wire_interval << " " << min_wire_interval
        //           << std::endl;
        std::map<int, std::set<int>>* map_max_index_wcps;
        std::map<int, std::set<int>>* map_min_index_wcps;
        if (mcell->get_max_wire_type() == 0) {
            map_max_index_wcps = &map_mcell_uindex_wcps[mcell];
        }
        else if (mcell->get_max_wire_type() == 1) {
            map_max_index_wcps = &map_mcell_vindex_wcps[mcell];
        }
        else {
            map_max_index_wcps = &map_mcell_windex_wcps[mcell];
        }
        if (mcell->get_min_wire_type() == 0) {
            map_min_index_wcps = &map_mcell_uindex_wcps[mcell];
        }
        else if (mcell->get_min_wire_type() == 1) {
            map_min_index_wcps = &map_mcell_vindex_wcps[mcell];
        }
        else {
            map_min_index_wcps = &map_mcell_windex_wcps[mcell];
        }

        for (const int pind1 : pinds) {
            int index_max_wire;
            int index_min_wire;
            if (mcell->get_max_wire_type() == 0) {
                index_max_wire = winds[0][pind1];
            }
            else if (mcell->get_max_wire_type() == 1) {
                index_max_wire = winds[1][pind1];
            }
            else {
                index_max_wire = winds[2][pind1];
            }
            if (mcell->get_min_wire_type() == 0) {
                index_min_wire = winds[0][pind1];
            }
            else if (mcell->get_min_wire_type() == 1) {
                index_min_wire = winds[1][pind1];
            }
            else {
                index_min_wire = winds[2][pind1];
            }

            std::vector<std::set<int>*> max_wcps_set;
            std::vector<std::set<int>*> min_wcps_set;

            // go through the first map and find the ones satisfying the condition
            for (auto it2 = map_max_index_wcps->begin(); it2 != map_max_index_wcps->end(); it2++) {
                if (fabs(it2->first - index_max_wire) <= max_wire_interval) {
                    max_wcps_set.push_back(&(it2->second));
                }
            }
            // go through the second map and find the ones satisfying the condition
            for (auto it2 = map_min_index_wcps->begin(); it2 != map_min_index_wcps->end(); it2++) {
                if (fabs(it2->first - index_min_wire) <= min_wire_interval) {
                    min_wcps_set.push_back(&(it2->second));
                }
            }

            std::set<int> wcps_set1;
            std::set<int> wcps_set2;

            for (auto it2 = max_wcps_set.begin(); it2 != max_wcps_set.end(); it2++) {
                wcps_set1.insert((*it2)->begin(), (*it2)->end());
            }
            for (auto it3 = min_wcps_set.begin(); it3 != min_wcps_set.end(); it3++) {
                wcps_set2.insert((*it3)->begin(), (*it3)->end());
            }

            {
                std::set<int> common_set;
                set_intersection(wcps_set1.begin(), wcps_set1.end(), wcps_set2.begin(), wcps_set2.end(),
                                 std::inserter(common_set, common_set.begin()));

                for (auto it4 = common_set.begin(); it4 != common_set.end(); it4++) {
                    int pind2 = *it4;
                    if (pind1 != pind2) {
                        // add edge ...
                        // auto edge = add_edge(pind1, pind2, *m_graph);
                        // if (edge.second) {
                        //     (*m_graph)[edge.first].dist = sqrt(pow(points[0][pind1] - points[0][pind2], 2) +
                        //                                      pow(points[1][pind1] - points[1][pind2], 2) +
                        //                                      pow(points[2][pind1] - points[2][pind2], 2));
                        //     num_edges++;
                        // }
                        auto edge = add_edge(pind1,pind2,WireCell::PointCloud::Facade::EdgeProp(sqrt(pow(points[0][pind1] - points[0][pind2], 2) +
                                                             pow(points[1][pind1] - points[1][pind2], 2) +
                                                             pow(points[2][pind1] - points[2][pind2], 2))),*m_graph);
	    //	    std::cout << index1 << " " << index2 << " " << edge.second << std::endl;
                        if (edge.second){
                            num_edges ++;
                        }
                    }
                }
            }
        }
    }

    LogDebug("in-blob edges: " << num_edges);
    // std::cout << "Test: in-blob edges: " << num_edges << std::endl;

    std::vector<int> time_slices;
    for (auto [time, _] : this->time_blob_map()) {
        time_slices.push_back(time);
    }

    // const int nticks_per_slice = grouping()->get_params().nticks_live_slice;
    // need to udpate to Multi-Face alg. ...
    const int nticks_per_slice = grouping()->get_nticks_per_slice().at(0).at(0);
    
    // std::cout << "Test: " << nticks_per_slice << std::endl;

    std::vector<std::pair<const Blob*, const Blob*>> connected_mcells;

    for (size_t i = 0; i != time_slices.size(); i++) {
        const auto& mcells_set = this->time_blob_map().at(time_slices.at(i));

        // create graph for points in mcell inside the same time slice
        if (mcells_set.size() >= 2) {
            for (auto it2 = mcells_set.begin(); it2 != mcells_set.end(); it2++) {
                auto mcell1 = *it2;
                auto it2p = it2;
                if (it2p != mcells_set.end()) {
                    it2p++;
                    for (auto it3 = it2p; it3 != mcells_set.end(); it3++) {
                        auto mcell2 = *(it3);
                        if (mcell1->overlap_fast(*mcell2, 2))
                            connected_mcells.push_back(std::make_pair(mcell1, mcell2));
                    }
                }
            }
        }
        // create graph for points between connected mcells in adjacent time slices + 1, if not, + 2
        std::vector<BlobSet> vec_mcells_set;
        if (i + 1 < time_slices.size()) {
            if (time_slices.at(i + 1) - time_slices.at(i) == 1*nticks_per_slice) {
                vec_mcells_set.push_back(this->time_blob_map().at(time_slices.at(i + 1)));
                if (i + 2 < time_slices.size())
                    if (time_slices.at(i + 2) - time_slices.at(i) == 2*nticks_per_slice)
                        vec_mcells_set.push_back(this->time_blob_map().at(time_slices.at(i + 2)));
            }
            else if (time_slices.at(i + 1) - time_slices.at(i) == 2*nticks_per_slice) {
                vec_mcells_set.push_back(this->time_blob_map().at(time_slices.at(i + 1)));
            }
        }
        //    bool flag = false;
        for (size_t j = 0; j != vec_mcells_set.size(); j++) {
            //      if (flag) break;
            auto& next_mcells_set = vec_mcells_set.at(j);
            for (auto it1 = mcells_set.begin(); it1 != mcells_set.end(); it1++) {
                auto mcell1 = (*it1);
                for (auto it2 = next_mcells_set.begin(); it2 != next_mcells_set.end(); it2++) {
                    auto mcell2 = (*it2);
                    if (mcell1->overlap_fast(*mcell2, 2)) {
                        //	    flag = true; // correct???
                        connected_mcells.push_back(std::make_pair(mcell1, mcell2));
                    }
                }
            }
        }
        // std::cout << "yuhw: itime_slices " << i
        // << " time_slices.at(i) " << time_slices.at(i)
        // << " vec_mcells_set  " << vec_mcells_set.size()
        // << " connected_mcells " << connected_mcells.size() << std::endl;
    }
    // std::cout << "connected_mcells size: " << connected_mcells.size() << std::endl;

    // establish edge ...
    const int max_num_nodes = 5;
    std::map<std::pair<int, int>, std::set<std::pair<double, int>>> closest_index;

    for (auto it = connected_mcells.begin(); it != connected_mcells.end(); it++) {
        auto mcell1 = (*it).first;
        auto mcell2 = (*it).second;

        std::vector<int> pinds1 = this->get_blob_indices(mcell1);
        std::vector<int> pinds2 = this->get_blob_indices(mcell2);

        // test 2 against 1 ...
        int max_wire_interval = mcell1->get_max_wire_interval();
        int min_wire_interval = mcell1->get_min_wire_interval();
        std::map<int, std::set<int>>* map_max_index_wcps;
        std::map<int, std::set<int>>* map_min_index_wcps;

        if (mcell1->get_max_wire_type() == 0) {
            map_max_index_wcps = &map_mcell_uindex_wcps.at(mcell2);
        }
        else if (mcell1->get_max_wire_type() == 1) {
            map_max_index_wcps = &map_mcell_vindex_wcps.at(mcell2);
        }
        else {
            map_max_index_wcps = &map_mcell_windex_wcps.at(mcell2);
        }
        if (mcell1->get_min_wire_type() == 0) {
            map_min_index_wcps = &map_mcell_uindex_wcps.at(mcell2);
        }
        else if (mcell1->get_min_wire_type() == 1) {
            map_min_index_wcps = &map_mcell_vindex_wcps.at(mcell2);
        }
        else {
            map_min_index_wcps = &map_mcell_windex_wcps.at(mcell2);
        }

        for (const int pind1 : pinds1) {
            int index_max_wire;
            int index_min_wire;
            if (mcell1->get_max_wire_type() == 0) {
                index_max_wire = winds[0][pind1];
            }
            else if (mcell1->get_max_wire_type() == 1) {
                index_max_wire = winds[1][pind1];
            }
            else {
                index_max_wire = winds[2][pind1];
            }
            if (mcell1->get_min_wire_type() == 0) {
                index_min_wire = winds[0][pind1];
            }
            else if (mcell1->get_min_wire_type() == 1) {
                index_min_wire = winds[1][pind1];
            }
            else {
                index_min_wire = winds[2][pind1];
            }
            std::vector<std::set<int>*> max_wcps_set;
            std::vector<std::set<int>*> min_wcps_set;
            // go through the first map and find the ones satisfying the condition
            for (auto it2 = map_max_index_wcps->begin(); it2 != map_max_index_wcps->end(); it2++) {
                if (fabs(it2->first - index_max_wire) <= max_wire_interval) {
                    max_wcps_set.push_back(&(it2->second));
                }
            }
            // go through the second map and find the ones satisfying the condition
            for (auto it2 = map_min_index_wcps->begin(); it2 != map_min_index_wcps->end(); it2++) {
                if (fabs(it2->first - index_min_wire) <= min_wire_interval) {
                    min_wcps_set.push_back(&(it2->second));
                }
            }

            std::set<int> wcps_set1;
            std::set<int> wcps_set2;

            for (auto it2 = max_wcps_set.begin(); it2 != max_wcps_set.end(); it2++) {
                wcps_set1.insert((*it2)->begin(), (*it2)->end());
            }
            for (auto it3 = min_wcps_set.begin(); it3 != min_wcps_set.end(); it3++) {
                wcps_set2.insert((*it3)->begin(), (*it3)->end());
            }

            {
                std::set<int> common_set;
                set_intersection(wcps_set1.begin(), wcps_set1.end(), wcps_set2.begin(), wcps_set2.end(),
                                 std::inserter(common_set, common_set.begin()));

                for (auto it4 = common_set.begin(); it4 != common_set.end(); it4++) {
                    const int pind2 = *it4;
                    if (pind1 != pind2) {
                        double dis = sqrt(pow(points[0][pind1] - points[0][pind2], 2) +
                                          pow(points[1][pind1] - points[1][pind2], 2) +
                                          pow(points[2][pind1] - points[2][pind2], 2));
                        auto b2 = blob_with_point(pind2);
                        auto key = std::make_pair(pind1, b2->slice_index_min());

                        if (closest_index.find(key) == closest_index.end() ) {
                            std::set<std::pair<double, int> > temp_sets;
                            temp_sets.insert(std::make_pair(dis, pind2));
                            closest_index[key] = temp_sets;
                        }
                        else {
                            closest_index[key].insert(std::make_pair(dis,pind2));
                            if (closest_index[key].size()>max_num_nodes){
                                auto it5 = closest_index[key].begin();
                                for (int qx = 0; qx!=max_num_nodes;qx++){
                                    it5++;
                                }
                                closest_index[key].erase(it5,closest_index[key].end());
                            }
                            // if (dis < closest_index[key].second || (std::abs(dis - closest_index[key].second) < 1e-10 && pind2 < closest_index[key].first)) 
                            // closest_index[key] = std::make_pair(pind2, dis);
                        }
                    }
                }
            }
        }

        // test 1 against 2 ...
        max_wire_interval = mcell2->get_max_wire_interval();
        min_wire_interval = mcell2->get_min_wire_interval();
        if (mcell2->get_max_wire_type() == 0) {
            map_max_index_wcps = &map_mcell_uindex_wcps[mcell1];
        }
        else if (mcell2->get_max_wire_type() == 1) {
            map_max_index_wcps = &map_mcell_vindex_wcps[mcell1];
        }
        else {
            map_max_index_wcps = &map_mcell_windex_wcps[mcell1];
        }
        if (mcell2->get_min_wire_type() == 0) {
            map_min_index_wcps = &map_mcell_uindex_wcps[mcell1];
        }
        else if (mcell2->get_min_wire_type() == 1) {
            map_min_index_wcps = &map_mcell_vindex_wcps[mcell1];
        }
        else {
            map_min_index_wcps = &map_mcell_windex_wcps[mcell1];
        }
        for (const int pind1 : pinds2) {
            int index_max_wire;
            int index_min_wire;
            if (mcell2->get_max_wire_type() == 0) {
                index_max_wire = winds[0][pind1];
            }
            else if (mcell2->get_max_wire_type() == 1) {
                index_max_wire = winds[1][pind1];
            }
            else {
                index_max_wire = winds[2][pind1];
            }
            if (mcell2->get_min_wire_type() == 0) {
                index_min_wire = winds[0][pind1];
            }
            else if (mcell2->get_min_wire_type() == 1) {
                index_min_wire = winds[1][pind1];
            }
            else {
                index_min_wire = winds[2][pind1];
            }
            std::vector<std::set<int>*> max_wcps_set;
            std::vector<std::set<int>*> min_wcps_set;
            // go through the first map and find the ones satisfying the condition
            for (auto it2 = map_max_index_wcps->begin(); it2 != map_max_index_wcps->end(); it2++) {
                if (fabs(it2->first - index_max_wire) <= max_wire_interval) {
                    max_wcps_set.push_back(&(it2->second));
                }
            }
            // go through the second map and find the ones satisfying the condition
            for (auto it2 = map_min_index_wcps->begin(); it2 != map_min_index_wcps->end(); it2++) {
                if (fabs(it2->first - index_min_wire) <= min_wire_interval) {
                    min_wcps_set.push_back(&(it2->second));
                }
            }

            std::set<int> wcps_set1;
            std::set<int> wcps_set2;

            for (auto it2 = max_wcps_set.begin(); it2 != max_wcps_set.end(); it2++) {
                wcps_set1.insert((*it2)->begin(), (*it2)->end());
            }
            for (auto it3 = min_wcps_set.begin(); it3 != min_wcps_set.end(); it3++) {
                wcps_set2.insert((*it3)->begin(), (*it3)->end());
            }

            {
                std::set<int> common_set;
                set_intersection(wcps_set1.begin(), wcps_set1.end(), wcps_set2.begin(), wcps_set2.end(),
                                 std::inserter(common_set, common_set.begin()));

                for (auto it4 = common_set.begin(); it4 != common_set.end(); it4++) {
                    const int pind2 = *it4;
                    if (pind1 != pind2) {
                        double dis = sqrt(pow(points[0][pind1] - points[0][pind2], 2) +
                                          pow(points[1][pind1] - points[1][pind2], 2) +
                                          pow(points[2][pind1] - points[2][pind2], 2));
                        auto b2 = blob_with_point(pind2);
                        auto key = std::make_pair(pind1, b2->slice_index_min());

                        if (closest_index.find(key) == closest_index.end()) {
                            std::set<std::pair<double, int> > temp_sets;
                            temp_sets.insert(std::make_pair(dis,pind2));
                            closest_index[key] = temp_sets;
                        }
                        else {
                            closest_index[key].insert(std::make_pair(dis,pind2));
                            if (closest_index[key].size()>max_num_nodes){
                                auto it5 = closest_index[key].begin();
                                for (int qx = 0; qx!=max_num_nodes;qx++){
                                    it5++;
                                }
                                closest_index[key].erase(it5,closest_index[key].end());
                            }
                            //if (dis < closest_index[key].second || (std::abs(dis - closest_index[key].second) < 1e-10 && pind2 < closest_index[key].first)) closest_index[key] = std::make_pair(pind2, dis);
                        }
                    }
                }
            }
        }
    }

    for (auto it4 = closest_index.begin(); it4 != closest_index.end(); it4++) {
        int index1 = it4->first.first;
        // int index2 = it4->second.first;
        // double dis = it4->second.second;
        // auto edge = add_edge(index1, index2, WireCell::PointCloud::Facade::EdgeProp(dis), *m_graph);
        // if (edge.second) {
        //     num_edges++;
        // }
        for (auto it5 = it4->second.begin(); it5!=it4->second.end(); it5++){
            int index2 = (*it5).second;
            double dis = (*it5).first;
            auto edge = add_edge(index1,index2,WireCell::PointCloud::Facade::EdgeProp(dis),*m_graph);
            if (edge.second){
                //      (*graph)[edge.first].dist = dis;
                num_edges ++;
            }
            // protect against dead cells ...
            //std::cout << dis/units::cm << std::endl;
            if (it5 == it4->second.begin() && dis > 0.25*units::cm)
                break;
        }

        // auto edge = add_edge(index1, index2, *m_graph);
        // if (edge.second) {
        //     (*m_graph)[edge.first].dist = dis;
        //     num_edges++;
        // }
        
    }
    // end of copying ...

    LogDebug("all edges: " << num_edges);
    // std::cout << "Test: all edges: " << num_edges << std::endl;

}

void Cluster::Connect_graph(const bool use_ctpc) const {
    // const auto& tp = grouping()->get_params();
    int hard_code_face = 0;
    // now form the connected components
    std::vector<int> component(num_vertices(*m_graph));
    const size_t num = connected_components(*m_graph, &component[0]);

    // Create ordered components
    std::vector<ComponentInfo> ordered_components;
    ordered_components.reserve(component.size());
    for (size_t i = 0; i < component.size(); ++i) {
        ordered_components.emplace_back(i);
    }

    // Assign vertices to components
    for (size_t i = 0; i < component.size(); ++i) {
        ordered_components[component[i]].add_vertex(i);
    }

    // Sort components by minimum vertex index
    std::sort(ordered_components.begin(), ordered_components.end(), 
        [](const ComponentInfo& a, const ComponentInfo& b) {
            return a.min_vertex < b.min_vertex;
        });

    LogDebug(" npoints " << npoints() << " nconnected " << num);
    if (num <= 1) return;

    std::vector<std::shared_ptr<Simple3DPointCloud>> pt_clouds;
    std::vector<std::vector<size_t>> pt_clouds_global_indices;
    // use this to link the global index to the local index
    // std::vector<std::vector<size_t>> pt_clouds_global_indices(num);
    // for (size_t i = 0; i != num; i++) {
    //     pt_clouds.push_back(std::make_shared<Simple3DPointCloud>());
    // }
    // for (size_t i = 0; i != component.size(); ++i) {
    //     pt_clouds.at(component[i])->add({points()[0][i], points()[1][i], points()[2][i]});
    //     pt_clouds_global_indices.at(component[i]).push_back(i);
    // }
    for (const auto& comp : ordered_components) {
        auto pt_cloud = std::make_shared<Simple3DPointCloud>();
        std::vector<size_t> global_indices;
        for (size_t vertex_idx : comp.vertex_indices) {
            pt_cloud->add({points()[0][vertex_idx], points()[1][vertex_idx], points()[2][vertex_idx]});
            global_indices.push_back(vertex_idx);
        }
        pt_clouds.push_back(pt_cloud);
        pt_clouds_global_indices.push_back(global_indices);
    }

    /// DEBUGONLY:
    if (0) {
        for (size_t i = 0; i != num; i++) {
            std::cout << *pt_clouds.at(i) << std::endl;
            std::cout << "global indices: ";
            for (size_t j = 0; j != pt_clouds_global_indices.at(i).size(); j++) {
                std::cout << pt_clouds_global_indices.at(i).at(j) << " ";
            }
            std::cout << std::endl;
        }
    }

    // Initiate dist. metrics
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis(
        num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_mst(
        num, std::vector<std::tuple<int, int, double>>(num));

    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir1(
        num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir2(
        num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir_mst(
        num, std::vector<std::tuple<int, int, double>>(num));

    for (size_t j = 0; j != num; j++) {
        for (size_t k = 0; k != num; k++) {
            index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_mst[j][k] = std::make_tuple(-1, -1, 1e9);

            index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir_mst[j][k] = std::make_tuple(-1, -1, 1e9);
        }
    }

    // Calc. dis, dis_dir1, dis_dir2
    // check against the closest distance ...
    // no need to have MST ...
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));

            

            // std::cout << j << " " << k << " " << std::get<0>(index_index_dis[j][k]) << " "
            //           << std::get<1>(index_index_dis[j][k]) << " " << std::get<2>(index_index_dis[j][k]) << " counter: " << global_counter_get_closest_wcpoint << std::endl;

            if ((num < 100 && pt_clouds.at(j)->get_num_points() > 100 && pt_clouds.at(k)->get_num_points() > 100 &&
                    (pt_clouds.at(j)->get_num_points() + pt_clouds.at(k)->get_num_points()) > 400) ||
                (pt_clouds.at(j)->get_num_points() > 500 && pt_clouds.at(k)->get_num_points() > 500)) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));

                geo_point_t dir1 = vhough_transform(p1, 30 * units::cm, HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                geo_point_t dir2 = vhough_transform(p2, 30 * units::cm, HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                dir1 = dir1 * -1;
                dir2 = dir2 * -1;

                std::pair<int, double> result1 = pt_clouds.at(k)->get_closest_point_along_vec(
                    p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);

                if (result1.first >= 0) {
                    index_index_dis_dir1[j][k] =
                        std::make_tuple(std::get<0>(index_index_dis[j][k]), result1.first, result1.second);
                }

                std::pair<int, double> result2 = pt_clouds.at(j)->get_closest_point_along_vec(
                    p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);

                if (result2.first >= 0) {
                    index_index_dis_dir2[j][k] =
                        std::make_tuple(result2.first, std::get<1>(index_index_dis[j][k]), result2.second);
                }
            }

            // Now check the path ...
            {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis / step_dis + 1;
                int num_bad = 0;
                geo_point_t test_p;
                for (int ii = 0; ii != num_steps; ii++) {
                    test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                               p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                               p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));
                    // if (!ct_point_cloud.is_good_point(test_p)) num_bad++;
                    if (use_ctpc) {
                        /// FIXME: assumes clusters are bounded to 1 face! Need to fix this.
                        const bool good_point = grouping()->is_good_point(test_p, hard_code_face);
                        if (!good_point) num_bad++;
                    }
                }

                if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                    index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
                }
            }

            // Now check the path ...
            if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir1[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir1[j][k]));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis / step_dis + 1;
                int num_bad = 0;
                geo_point_t test_p;
                for (int ii = 0; ii != num_steps; ii++) {
                    test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                               p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                               p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));
                    // if (!ct_point_cloud.is_good_point(test_p)) num_bad++;
                    if (use_ctpc) {
                        /// FIXME: assumes clusters are bounded to 1 face! Need to fix this.
                        const bool good_point = grouping()->is_good_point(test_p, hard_code_face);
                        if (!good_point) num_bad++;
                    }
                }

                if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                    index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                }
            }

            // Now check the path ...
            if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir2[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir2[j][k]));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis / step_dis + 1;
                int num_bad = 0;
                geo_point_t test_p;
                for (int ii = 0; ii != num_steps; ii++) {
                    test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                               p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                               p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));
                    // if (!ct_point_cloud.is_good_point(test_p)) num_bad++;
                    if (use_ctpc) {
                        /// FIXME: assumes clusters are bounded to 1 face! Need to fix this.
                        const bool good_point = grouping()->is_good_point(test_p, hard_code_face);
                        if (!good_point) num_bad++;
                    }
                }

                if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                    index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                }
            }
        }
    }

    // deal with MST of first type
    {
        boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, boost::no_property,
                              boost::property<boost::edge_weight_t, double>>
            temp_graph(num);

        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j;
                int index2 = k;
                if (std::get<0>(index_index_dis[j][k]) >= 0) {
                    add_edge(index1, index2, std::get<2>(index_index_dis[j][k]), temp_graph);
                    // LogDebug(index1 << " " << index2 << " " << std::get<2>(index_index_dis[j][k]));
                }
            }
        }

        // Process MST
        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_mst);
    }

    // MST of the direction ...
    {
        boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, boost::no_property,
                              boost::property<boost::edge_weight_t, double>>
            temp_graph(num);

        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j;
                int index2 = k;
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0 || std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    add_edge(
                        index1, index2,
                        std::min(std::get<2>(index_index_dis_dir1[j][k]), std::get<2>(index_index_dis_dir2[j][k])),
                        temp_graph);
                    // LogDebug(index1 << " " << index2 << " "
                    //                 << std::min(std::get<2>(index_index_dis_dir1[j][k]),
                    //                             std::get<2>(index_index_dis_dir2[j][k])));
                }
            }
        }

        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_dir_mst);

    }

    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            if (std::get<2>(index_index_dis[j][k]) < 3 * units::cm) {
                index_index_dis_mst[j][k] = index_index_dis[j][k];
            }

            // establish the path ...
            if (std::get<0>(index_index_dis_mst[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_mst[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_mst[j][k]));
                // auto edge =
                //     add_edge(gind1, gind2, *m_graph);
                    // LogDebug(gind1 << " " << gind2 << " " << std::get<2>(index_index_dis_mst[j][k]));
                float dis;
                // if (edge.second) {
                if (std::get<2>(index_index_dis_mst[j][k]) > 5 * units::cm) {
                    dis = std::get<2>(index_index_dis_mst[j][k]);
                }
                else {
                    dis = std::get<2>(index_index_dis_mst[j][k]);
                }
                // }
                /*auto edge =*/ add_edge(gind1, gind2, WireCell::PointCloud::Facade::EdgeProp(dis), *m_graph);
            }

            if (std::get<0>(index_index_dis_dir_mst[j][k]) >= 0) {
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k]));
                    //auto edge = add_edge(gind1, gind2, *m_graph);
                    // LogDebug(gind1 << " " << gind2 << " " << std::get<2>(index_index_dis_dir1[j][k]));
                    float dis;
                    // if (edge.second) {
                    if (std::get<2>(index_index_dis_dir1[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir1[j][k]) * 1.1;
                    }
                    else {
                        dis = std::get<2>(index_index_dis_dir1[j][k]);
                    }
                    // }
                    /*auto edge =*/ add_edge(gind1, gind2, WireCell::PointCloud::Facade::EdgeProp(dis), *m_graph);
                }
                if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k]));
                    // auto edge = add_edge(gind1, gind2, *m_graph);
                    // LogDebug(gind1 << " " << gind2 << " " << std::get<2>(index_index_dis_dir2[j][k]));
                    // if (edge.second) {
                    float dis;
                    if (std::get<2>(index_index_dis_dir2[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir2[j][k]) * 1.1;
                    }
                    else {
                        dis = std::get<2>(index_index_dis_dir2[j][k]);
                    }
                    // }
                    /*auto edge =*/ add_edge(gind1, gind2, WireCell::PointCloud::Facade::EdgeProp(dis), *m_graph);
                }
            }

        }  // k
    }  // j
}
// #define LogDebug(x)

void Cluster::Connect_graph() const{
    // now form the connected components
    std::vector<int> component(num_vertices(*m_graph));
    const size_t num = connected_components(*m_graph, &component[0]);

    // Create ordered components
    std::vector<ComponentInfo> ordered_components;
    ordered_components.reserve(component.size());
    for (size_t i = 0; i < component.size(); ++i) {
        ordered_components.emplace_back(i);
    }

    // Assign vertices to components
    for (size_t i = 0; i < component.size(); ++i) {
        ordered_components[component[i]].add_vertex(i);
    }

    // Sort components by minimum vertex index
    std::sort(ordered_components.begin(), ordered_components.end(), 
        [](const ComponentInfo& a, const ComponentInfo& b) {
            return a.min_vertex < b.min_vertex;
        });

    LogDebug(" npoints " << npoints() << " nconnected " << num);
    if (num <= 1) return;

    std::vector<std::shared_ptr<Simple3DPointCloud>> pt_clouds;
    std::vector<std::vector<size_t>> pt_clouds_global_indices;
    // use this to link the global index to the local index
    // std::vector<std::vector<size_t>> pt_clouds_global_indices(num);
    // for (size_t i = 0; i != num; i++) {
    //     pt_clouds.push_back(std::make_shared<Simple3DPointCloud>());
    // }
    // for (size_t i = 0; i != component.size(); ++i) {
    //     pt_clouds.at(component[i])->add({points()[0][i], points()[1][i], points()[2][i]});
    //     pt_clouds_global_indices.at(component[i]).push_back(i);
    // }
    // Create point clouds using ordered components
    for (const auto& comp : ordered_components) {
        auto pt_cloud = std::make_shared<Simple3DPointCloud>();
        std::vector<size_t> global_indices;
        
        for (size_t vertex_idx : comp.vertex_indices) {
            pt_cloud->add({points()[0][vertex_idx], points()[1][vertex_idx], points()[2][vertex_idx]});
            global_indices.push_back(vertex_idx);
        }
        
        pt_clouds.push_back(pt_cloud);
        pt_clouds_global_indices.push_back(global_indices);
    }

    /// DEBUGONLY:
    if (0) {
        for (size_t i = 0; i != num; i++) {
            std::cout << *pt_clouds.at(i) << std::endl;
            std::cout << "global indices: ";
            for (size_t j = 0; j != pt_clouds_global_indices.at(i).size(); j++) {
                std::cout << pt_clouds_global_indices.at(i).at(j) << " ";
            }
            std::cout << std::endl;
        }
    }

     // Initiate dist. metrics
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis(
        num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_mst(
        num, std::vector<std::tuple<int, int, double>>(num));

    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir1(
        num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir2(
        num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir_mst(
        num, std::vector<std::tuple<int, int, double>>(num));

    for (size_t j = 0; j != num; j++) {
        for (size_t k = 0; k != num; k++) {
            index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_mst[j][k] = std::make_tuple(-1, -1, 1e9);

            index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir_mst[j][k] = std::make_tuple(-1, -1, 1e9);
        }
    }

    boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, boost::no_property,
                              boost::property<boost::edge_weight_t, double>>
    temp_graph(num);

    for (size_t j=0;j!=num;j++){
      for (size_t k=j+1;k!=num;k++){
            index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));

            // geo_point_t test_p3 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
            // geo_point_t test_p4 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));
            // if (nchildren()==3449) std::cout << "A0: " << test_p3 << " " << test_p4 << " " << j << " " << k << " " << pt_clouds.at(j)->get_num_points() << " " << pt_clouds.at(k)->get_num_points() << " " << std::get<2>(index_index_dis[j][k])/units::cm << std::endl;
                

            int index1 = j;
            int index2 = k;
            /*auto edge =*/ add_edge(index1,index2, std::get<2>(index_index_dis[j][k]), temp_graph);
      }
    }

    process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_mst);

    // {
    //     std::vector<int> possible_root_vertex;
    //     std::vector<int> component(num_vertices(temp_graph));
    //     const int num1 = connected_components(temp_graph, &component[0]);
    //     possible_root_vertex.resize(num1);
    //     std::vector<int>::size_type i;
    //     for (i = 0; i != component.size(); ++i) {
    //         possible_root_vertex.at(component[i]) = i;
    //     }

    //     for (size_t i = 0; i != possible_root_vertex.size(); i++) {
    //         std::vector<boost::graph_traits<MCUGraph>::vertex_descriptor> predecessors(num_vertices(temp_graph));

    //         prim_minimum_spanning_tree(temp_graph, &predecessors[0],
    //                                     boost::root_vertex(possible_root_vertex.at(i)));

    //         for (size_t j = 0; j != predecessors.size(); ++j) {
    //             if (predecessors[j] != j) {
    //                 if (j < predecessors[j]) {
    //                     index_index_dis_mst[j][predecessors[j]] = index_index_dis[j][predecessors[j]];
    //                 }
    //                 else {
    //                     index_index_dis_mst[predecessors[j]][j] = index_index_dis[predecessors[j]][j];
    //                 }
    //                 // std::cout << j << " " << predecessors[j] << " " << std::endl;
    //             }
    //             else {
    //                 // std::cout << j << " " << std::endl;
    //             }
    //         }
    //     }
    // }

    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            if (std::get<2>(index_index_dis[j][k])<3*units::cm){
	            index_index_dis_mst[j][k] = index_index_dis[j][k];
	        }

            if (num < 100)
	        if (pt_clouds.at(j)->get_num_points()>100 && pt_clouds.at(k)->get_num_points()>100 &&
	            (pt_clouds.at(j)->get_num_points()+pt_clouds.at(k)->get_num_points()) > 400){
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));
            
                geo_point_t dir1 = vhough_transform(p1, 30 * units::cm, HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                geo_point_t dir2 = vhough_transform(p2, 30 * units::cm, HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                dir1 = dir1 * -1;
                dir2 = dir2 * -1;

                std::pair<int, double> result1 = pt_clouds.at(k)->get_closest_point_along_vec(
                    p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);

                if (result1.first >= 0) {
                    index_index_dis_dir1[j][k] =
                        std::make_tuple(std::get<0>(index_index_dis[j][k]), result1.first, result1.second);
                }

                std::pair<int, double> result2 = pt_clouds.at(j)->get_closest_point_along_vec(
                    p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);

                if (result2.first >= 0) {
                    index_index_dis_dir2[j][k] =
                        std::make_tuple(result2.first, std::get<1>(index_index_dis[j][k]), result2.second);
                }
            }
        }
    }

    // MST for the directionality ...
    {
        boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, boost::no_property,
                              boost::property<boost::edge_weight_t, double>>
        temp_graph(num);

        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j;
                int index2 = k;
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0 || std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    add_edge(
                        index1, index2,
                        std::min(std::get<2>(index_index_dis_dir1[j][k]), std::get<2>(index_index_dis_dir2[j][k])),
                        temp_graph);
                    // LogDebug(index1 << " " << index2 << " "
                    //                 << std::min(std::get<2>(index_index_dis_dir1[j][k]),
                    //                             std::get<2>(index_index_dis_dir2[j][k])));
                }
            }
        }

        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_dir_mst);
        // {
        //     std::vector<int> possible_root_vertex;
        //     std::vector<int> component(num_vertices(temp_graph));
        //     const int num1 = connected_components(temp_graph, &component[0]);
        //     possible_root_vertex.resize(num1);
        //     std::vector<int>::size_type i;
        //     for (i = 0; i != component.size(); ++i) {
        //         possible_root_vertex.at(component[i]) = i;
        //     }

        //     for (size_t i = 0; i != possible_root_vertex.size(); i++) {
        //         std::vector<boost::graph_traits<MCUGraph>::vertex_descriptor> predecessors(num_vertices(temp_graph));
        //         prim_minimum_spanning_tree(temp_graph, &predecessors[0],
        //                                    boost::root_vertex(possible_root_vertex.at(i)));
        //         for (size_t j = 0; j != predecessors.size(); ++j) {
        //             if (predecessors[j] != j) {
        //                 if (j < predecessors[j]) {
        //                     index_index_dis_dir_mst[j][predecessors[j]] = index_index_dis[j][predecessors[j]];
        //                 }
        //                 else {
        //                     index_index_dis_dir_mst[predecessors[j]][j] = index_index_dis[predecessors[j]][j];
        //                 }
        //                 // std::cout << j << " " << predecessors[j] << " " << std::endl;
        //             }
        //             else {
        //                 // std::cout << j << " " << std::endl;
        //             }
        //         }
        //     }
        // }
    }

    // now complete graph according to the direction
    // according to direction ...
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            if (std::get<0>(index_index_dis_mst[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_mst[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_mst[j][k]));

                // geo_point_t test_p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_mst[j][k]));
                // geo_point_t test_p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_mst[j][k]));
                // if (nchildren()==3449) std::cout << "A1: " << test_p1 << " " << test_p2 << " " << j << " " << k << " " << pt_clouds.at(j)->get_num_points() << " " << pt_clouds.at(k)->get_num_points() << " " << std::get<2>(index_index_dis_mst[j][k])/units::cm << std::endl;

                
                // auto edge =
                //     add_edge(gind1, gind2, *m_graph);
                //     // LogDebug(gind1 << " " << gind2 << " " << std::get<2>(index_index_dis_mst[j][k]));
                // if (edge.second) {
                float dis;
                if (std::get<2>(index_index_dis_mst[j][k]) > 5 * units::cm) {
                    dis = std::get<2>(index_index_dis_mst[j][k]);
                }
                else {
                    dis = std::get<2>(index_index_dis_mst[j][k]);
                }
                // }
                /*auto edge =*/ add_edge(gind1, gind2, WireCell::PointCloud::Facade::EdgeProp(dis), *m_graph);
            }

            if (std::get<0>(index_index_dis_dir_mst[j][k]) >= 0) {
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k]));

                    // geo_point_t test_p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir1[j][k]));
                    // geo_point_t test_p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir1[j][k]));
                    // if (nchildren()==3449) std::cout << "A2: " << test_p1 << " " << test_p2 << " " << j << " " << k << " " << pt_clouds.at(j)->get_num_points() << " " << pt_clouds.at(k)->get_num_points()<< " " << std::get<2>(index_index_dis_dir1[j][k])/units::cm << std::endl;

                    // auto edge = add_edge(gind1, gind2, *m_graph);
                    // // LogDebug(gind1 << " " << gind2 << " " << std::get<2>(index_index_dis_dir1[j][k]));
                    // if (edge.second) {
                    float dis;
                    if (std::get<2>(index_index_dis_dir1[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir1[j][k]) * 1.2;
                    }
                    else {
                        dis = std::get<2>(index_index_dis_dir1[j][k]);
                    }
                    // }
                    /*auto edge =*/ add_edge(gind1, gind2, WireCell::PointCloud::Facade::EdgeProp(dis), *m_graph);
                }
                if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k]));
                    // auto edge = add_edge(gind1, gind2, *m_graph);
                    // // LogDebug(gind1 << " " << gind2 << " " << std::get<2>(index_index_dis_dir2[j][k]));
                    // if (edge.second) {

                    // geo_point_t test_p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir2[j][k]));
                    // geo_point_t test_p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir2[j][k]));
                    // if (nchildren()==3449) std::cout << "A3: " << test_p1 << " " << test_p2 << " " << j << " " << k << " " << pt_clouds.at(j)->get_num_points() << " " << pt_clouds.at(k)->get_num_points()<< " " << std::get<2>(index_index_dis_dir2[j][k])/units::cm << std::endl;

                    float dis;
                    if (std::get<2>(index_index_dis_dir2[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir2[j][k]) * 1.2;
                    }
                    else {
                        dis = std::get<2>(index_index_dis_dir2[j][k]);
                    }
                    // }
                    /*auto edge =*/ add_edge(gind1, gind2, WireCell::PointCloud::Facade::EdgeProp(dis), *m_graph);
                }
            }

        }
    }


}


void Cluster::Connect_graph_overclustering_protection(const IDetectorVolumes::pointer dv, const bool use_ctpc) const {
    // Get all the wire plane IDs from the grouping
    const auto& wpids = grouping()->wpids();
    // Key: pair<APA, face>, Value: drift_dir, angle_u, angle_v, angle_w
    std::map<WirePlaneId , std::tuple<geo_point_t, double, double, double>> wpid_params;
    std::map<WirePlaneId, geo_point_t> wpid_U_dir;
    std::map<WirePlaneId, geo_point_t> wpid_V_dir;
    std::map<WirePlaneId, geo_point_t> wpid_W_dir;
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
        wpid_U_dir[wpid] = geo_point_t(0, cos(angle_u), sin(angle_u));
        wpid_V_dir[wpid] = geo_point_t(0, cos(angle_v), sin(angle_v));
        wpid_W_dir[wpid] = geo_point_t(0, cos(angle_w), sin(angle_w));
    }

    // Constants for wire angles
    // const auto& tp = grouping()->get_params();
    int hard_code_face = 0;
    //std::cout << "Test: face " << tp.face << std::endl;

    // const double pi = 3.141592653589793;
    // this drift direction is only used to calculate isochronous case, so this is OK ...
    const geo_vector_t drift_dir_abs(1, 0, 0); 

    // need to understand the points before implementing the angles ...
    const auto [angle_u,angle_v,angle_w] = grouping()->wire_angles();
    const geo_point_t U_dir(0,cos(angle_u),sin(angle_u));
    const geo_point_t V_dir(0,cos(angle_v),sin(angle_v));
    const geo_point_t W_dir(0,cos(angle_w),sin(angle_w));


    // Form connected components
    std::vector<int> component(num_vertices(*m_graph));
    const size_t num = connected_components(*m_graph, &component[0]);
    
    LogDebug(" npoints " << npoints() << " nconnected " << num);
    if (num <= 1) return;

    // Create point clouds using connected components
    std::vector<std::shared_ptr<Simple3DPointCloud>> pt_clouds;
    std::vector<std::vector<size_t>> pt_clouds_global_indices;
    
    // Create ordered components  
    std::vector<ComponentInfo> ordered_components;
    ordered_components.reserve(component.size());
    for (size_t i = 0; i < component.size(); ++i) {
        ordered_components.emplace_back(i);
    }
    
    // Assign vertices to components
    for (size_t i = 0; i < component.size(); ++i) {
        ordered_components[component[i]].add_vertex(i);
    }

    // Sort components by minimum vertex index
    std::sort(ordered_components.begin(), ordered_components.end(),
        [](const ComponentInfo& a, const ComponentInfo& b) {
            return a.min_vertex < b.min_vertex;
        });

    // Create point clouds for each component
    for (const auto& comp : ordered_components) {
        auto pt_cloud = std::make_shared<Simple3DPointCloud>();
        std::vector<size_t> global_indices;
        
        for (size_t vertex_idx : comp.vertex_indices) {
            pt_cloud->add({points()[0][vertex_idx], points()[1][vertex_idx], points()[2][vertex_idx]});
            global_indices.push_back(vertex_idx);
        }
        pt_clouds.push_back(pt_cloud);
        pt_clouds_global_indices.push_back(global_indices);
    }

    // pt_clouds.resize(num);
    // pt_clouds_global_indices.resize(num);

    // // Initialize all point clouds
    // for(size_t j = 0; j < num; j++) {
    //     pt_clouds[j] = std::make_shared<Simple3DPointCloud>();
    // }

    // // Add points directly using component mapping
    // for(size_t i = 0; i < component.size(); ++i) {
    //     pt_clouds[component[i]]->add({points()[0][i], points()[1][i], points()[2][i]});
    //     pt_clouds_global_indices[component[i]].push_back(i);
    // }

    // std::cout << "Test: "<< num << std::endl;

    // Initialize distance metrics 
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis(num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_mst(num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir1(num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir2(num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir_mst(num, std::vector<std::tuple<int, int, double>>(num));

    // Initialize all distances to inf
    for (size_t j = 0; j != num; j++) {
        for (size_t k = 0; k != num; k++) {
            index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_mst[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir_mst[j][k] = std::make_tuple(-1, -1, 1e9);
        }
    }

    // Calculate distances between components
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            // Get closest points between components
            index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));

            // Skip small clouds
            if ((num < 100 && pt_clouds.at(j)->get_num_points() > 100 && pt_clouds.at(k)->get_num_points() > 100 &&
                 (pt_clouds.at(j)->get_num_points() + pt_clouds.at(k)->get_num_points()) > 400) ||
                (pt_clouds.at(j)->get_num_points() > 500 && pt_clouds.at(k)->get_num_points() > 500)) {
                
                // Get closest points and calculate directions
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));

                geo_vector_t dir1 = vhough_transform(p1, 30 * units::cm, HoughParamSpace::theta_phi, pt_clouds.at(j), 
                                                pt_clouds_global_indices.at(j));
                geo_vector_t dir2 = vhough_transform(p2, 30 * units::cm, HoughParamSpace::theta_phi, pt_clouds.at(k),
                                                pt_clouds_global_indices.at(k)); 
                dir1 = dir1 * -1;
                dir2 = dir2 * -1;

                std::pair<int, double> result1 = pt_clouds.at(k)->get_closest_point_along_vec(p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);

                if (result1.first >= 0) {
                    index_index_dis_dir1[j][k] = std::make_tuple(std::get<0>(index_index_dis[j][k]), 
                                                                result1.first, result1.second);
                }

                std::pair<int, double> result2 = pt_clouds.at(j)->get_closest_point_along_vec(p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm); 

                if (result2.first >= 0) {
                    index_index_dis_dir2[j][k] = std::make_tuple(result2.first,
                                                                std::get<1>(index_index_dis[j][k]), 
                                                                result2.second);
                }
            }
            // Now check the path 

            {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis/step_dis + 1;

                

                // Track different types of "bad" points
                int num_bad[4] = {0,0,0,0};   // more than one of three are bad
                int num_bad1[4] = {0,0,0,0};  // at least one of three are bad
                int num_bad2[3] = {0,0,0};    // number of dead channels

                // Check points along path
                for (int ii = 0; ii != num_steps; ii++) {
                    geo_point_t test_p(
                        p1.x() + (p2.x() - p1.x())/num_steps*(ii + 1),
                        p1.y() + (p2.y() - p1.y())/num_steps*(ii + 1),
                        p1.z() + (p2.z() - p1.z())/num_steps*(ii + 1)
                    );

                    // Test point quality using grouping parameters
                    std::vector<int> scores;
                    if (use_ctpc) {
                        scores = grouping()->test_good_point(test_p, hard_code_face);
                        
                        // Check overall quality
                        if (scores[0] + scores[3] + scores[1] + scores[4] + (scores[2]+scores[5])*2 < 3) {
                            num_bad[0]++;
                        }
                        if (scores[0]+scores[3]==0) num_bad[1]++;
                        if (scores[1]+scores[4]==0) num_bad[2]++;
                        if (scores[2]+scores[5]==0) num_bad[3]++;

                        if (scores[3]!=0) num_bad2[0]++;
                        if (scores[4]!=0) num_bad2[1]++;
                        if (scores[5]!=0) num_bad2[2]++;
                        
                        if (scores[0] + scores[3] + scores[1] + scores[4] + (scores[2]+scores[5]) < 3) {
                            num_bad1[0]++;
                        }
                        if (scores[0]+scores[3]==0) num_bad1[1]++;
                        if (scores[1]+scores[4]==0) num_bad1[2]++;
                        if (scores[2]+scores[5]==0) num_bad1[3]++;
                    }
                }

                // if (kd_blobs().size()==244){
                //     std::cout << "Test: Dis: " << p1 << " " << p2 << " " << dis << std::endl;
                //     std::cout << "Test: num_bad1: " << num_bad1[0] << " " << num_bad1[1] << " " << num_bad1[2] << " " << num_bad1[3] << std::endl;
                //     std::cout << "Test: num_bad2: " << num_bad2[0] << " " << num_bad2[1] << " " << num_bad2[2] << std::endl;
                //     std::cout << "Test: num_bad: " << num_bad[0] << " " << num_bad[1] << " " << num_bad[2] << " " << num_bad[3] << std::endl;
                // }
                // Calculate angles between directions
                geo_vector_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
                geo_vector_t tempV5;

                double angle1 = tempV1.angle(U_dir); 
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2)) * sin(angle1),
                        0);
                angle1 = tempV5.angle(drift_dir_abs);

                double angle2 = tempV1.angle(V_dir);
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2)) * sin(angle2),
                        0);
                angle2 = tempV5.angle(drift_dir_abs);

                double angle1p = tempV1.angle(W_dir);
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2)) * sin(angle1p),
                        0); 
                angle1p = tempV5.angle(drift_dir_abs);

                tempV5.set(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
                double angle3 = tempV5.angle(drift_dir_abs);

                bool flag_strong_check = true;

                // Define constants for readability
                constexpr double pi = 3.141592653589793;
                constexpr double perp_angle_tol = 10.0/180.0*pi;
                constexpr double wire_angle_tol = 12.5/180.0*pi;
                constexpr double perp_angle = pi/2.0;
                constexpr double invalid_dist = 1e9;

                if (fabs(angle3 - perp_angle) < perp_angle_tol) {
                    geo_vector_t tempV2 = vhough_transform(p1, 15*units::cm);
                    geo_vector_t tempV3 = vhough_transform(p2, 15*units::cm);
                    
                    if (fabs(tempV2.angle(drift_dir_abs) - perp_angle) < perp_angle_tol &&
                        fabs(tempV3.angle(drift_dir_abs) - perp_angle) < perp_angle_tol) {
                        flag_strong_check = false;
                    }
                }
                else if (angle1 < wire_angle_tol || angle2 < wire_angle_tol || angle1p < wire_angle_tol) {
                    flag_strong_check = false;
                }

                // Helper function to check if ratio exceeds threshold
                auto exceeds_ratio = [](int val, int steps, double ratio = 0.75) {
                    return val >= ratio * steps;
                };

                // Helper function to invalidate distance
                auto invalidate_distance = [&]() {
                    index_index_dis[j][k] = std::make_tuple(-1, -1, invalid_dist);
                };

                if (flag_strong_check) {
                    if (num_bad1[0] > 7 || (num_bad1[0] > 2 && exceeds_ratio(num_bad1[0], num_steps))) {
                        invalidate_distance();
                    }
                }
                else {
                    bool parallel_angles = (angle1 < wire_angle_tol && angle2 < wire_angle_tol) ||
                                        (angle1p < wire_angle_tol && angle1 < wire_angle_tol) ||
                                        (angle1p < wire_angle_tol && angle2 < wire_angle_tol);

                    if (parallel_angles) {
                        if (num_bad[0] > 7 || (num_bad[0] > 2 && exceeds_ratio(num_bad[0], num_steps))) {
                            invalidate_distance();
                        }
                    }
                    else if (angle1 < wire_angle_tol) {
                        int sum_bad = num_bad[2] + num_bad[3];
                        if (sum_bad > 9 || (sum_bad > 2 && exceeds_ratio(sum_bad, num_steps)) || num_bad[3] >= 3) {
                            invalidate_distance();
                        }
                    }
                    else if (angle2 < wire_angle_tol) {
                        int sum_bad = num_bad[1] + num_bad[3];
                        if (sum_bad > 9 || (sum_bad > 2 && exceeds_ratio(sum_bad, num_steps)) || num_bad[3] >= 3) {
                            invalidate_distance();
                        }
                    }
                    else if (angle1p < wire_angle_tol) {
                        int sum_bad = num_bad[2] + num_bad[1];
                        if (sum_bad > 9 || (sum_bad > 2 && exceeds_ratio(sum_bad, num_steps))) {
                            invalidate_distance();
                        }
                    }
                    else if (num_bad[0] > 7 || (num_bad[0] > 2 && exceeds_ratio(num_bad[0], num_steps))) {
                        invalidate_distance();
                    }
                }
            }

            // Now check path again ... 
            if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir1[j][k])); //point3d(std::get<0>(index_index_dis_dir1[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir1[j][k])); //point3d(std::get<1>(index_index_dis_dir1[j][k]));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + 
                                pow(p1.y() - p2.y(), 2) + 
                                pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis/step_dis + 1;
                int num_bad = 0;
                int num_bad1 = 0;

                // Check intermediate points along path
                for (int ii = 0; ii != num_steps; ii++) {
                    geo_point_t test_p(
                        p1.x() + (p2.x() - p1.x())/num_steps*(ii + 1),
                        p1.y() + (p2.y() - p1.y())/num_steps*(ii + 1),
                        p1.z() + (p2.z() - p1.z())/num_steps*(ii + 1)
                    );

                    if (use_ctpc) {
                        /// FIXME: assumes clusters are bounded to 1 face! Need to fix this.
                        const bool good_point = grouping()->is_good_point(test_p, hard_code_face);
                        if (!good_point) {
                            num_bad++;
                        }
                        if (!grouping()->is_good_point(test_p, hard_code_face, 0.6*units::cm, 1, 0)) {
                            num_bad1++;
                        }
                    }
                }
                
                // Calculate angles
                geo_vector_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
                geo_vector_t tempV5;
                
                double angle1 = tempV1.angle(U_dir);
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle1),
                        0);
                angle1 = tempV5.angle(drift_dir_abs);
                
                double angle2 = tempV1.angle(V_dir);
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle2),
                        0);
                angle2 = tempV5.angle(drift_dir_abs);
                
                tempV5.set(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
                double angle3 = tempV5.angle(drift_dir_abs);
                
                double angle1p = tempV1.angle(W_dir);
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle1p),
                        0);
                angle1p = tempV5.angle(drift_dir_abs);

                const double pi = 3.141592653589793;
                if (fabs(angle3 - pi/2) < 10.0/180.0*pi || 
                    angle1 < 12.5/180.0*pi ||
                    angle2 < 12.5/180.0*pi || 
                    angle1p < 7.5/180.0*pi) {
                    // Parallel or prolonged case
                    if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75*num_steps)) {
                        index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }
                else {
                    if (num_bad1 > 7 || (num_bad1 > 2 && num_bad1 >= 0.75*num_steps)) {
                        index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }
            }

            //Now check path again ... 
            // Now check the path...
            if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir2[j][k]));//point3d(std::get<0>(index_index_dis_dir2[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir2[j][k]));//point3d(std::get<1>(index_index_dis_dir2[j][k]));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + 
                                pow(p1.y() - p2.y(), 2) + 
                                pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis/step_dis + 1;
                int num_bad = 0;
                int num_bad1 = 0;

                // Check points along path
                for (int ii = 0; ii != num_steps; ii++) {
                    geo_point_t test_p(
                        p1.x() + (p2.x() - p1.x())/num_steps*(ii + 1),
                        p1.y() + (p2.y() - p1.y())/num_steps*(ii + 1),
                        p1.z() + (p2.z() - p1.z())/num_steps*(ii + 1)
                    );

                    if (use_ctpc) {
                        /// FIXME: assumes clusters are bounded to 1 face! Need to fix this.
                        const bool good_point = grouping()->is_good_point(test_p, hard_code_face);
                        if (!good_point) {
                            num_bad++;
                        }
                        if (!grouping()->is_good_point(test_p, hard_code_face, 0.6*units::cm, 1, 0)) {
                            num_bad1++;
                        }
                    }
                }

                // Calculate angles between directions
                geo_vector_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
                geo_vector_t tempV5;

                double angle1 = tempV1.angle(U_dir);
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle1),
                        0);
                angle1 = tempV5.angle(drift_dir_abs);

                double angle2 = tempV1.angle(V_dir);
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle2),
                        0);
                angle2 = tempV5.angle(drift_dir_abs);

                tempV5.set(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
                double angle3 = tempV5.angle(drift_dir_abs);

                double angle1p = tempV1.angle(W_dir);
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle1p),
                        0);
                angle1p = tempV5.angle(drift_dir_abs);

                const double pi = 3.141592653589793;
                bool is_parallel = fabs(angle3 - pi/2) < 10.0/180.0*pi || 
                                angle1 < 12.5/180.0*pi ||
                                angle2 < 12.5/180.0*pi || 
                                angle1p < 7.5/180.0*pi;

                if (is_parallel) {
                    // Parallel or prolonged case
                    if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75*num_steps)) {
                        index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }
                else {
                    if (num_bad1 > 7 || (num_bad1 > 2 && num_bad1 >= 0.75*num_steps)) {
                        index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }
            }
        }
    }

    // deal with MST of first type
    {
        boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, boost::no_property,
                              boost::property<boost::edge_weight_t, double>>
            temp_graph(num);
        // int temp_count = 0;
        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j;
                int index2 = k;
                if (std::get<0>(index_index_dis[j][k]) >= 0) {
                    add_edge(index1, index2, std::get<2>(index_index_dis[j][k]), temp_graph);
                    // LogDebug(index1 << " " << index2 << " " << std::get<2>(index_index_dis[j][k]));
                    // temp_count ++;
                }
            }
        }
        // std::cout << "Test: Count: " << temp_count << std::endl;

        // Process MST
        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_mst);
    }

    // MST of the direction ...
    {
        boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, boost::no_property,
                              boost::property<boost::edge_weight_t, double>>
            temp_graph(num);

        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j;
                int index2 = k;
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0 || std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    add_edge(
                        index1, index2,
                        std::min(std::get<2>(index_index_dis_dir1[j][k]), std::get<2>(index_index_dis_dir2[j][k])),
                        temp_graph);
                    // LogDebug(index1 << " " << index2 << " "
                    //                 << std::min(std::get<2>(index_index_dis_dir1[j][k]),
                    //                             std::get<2>(index_index_dis_dir2[j][k])));
                }
            }
        }

        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_dir_mst);

    }

    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            if (std::get<2>(index_index_dis[j][k]) < 3 * units::cm) {
                index_index_dis_mst[j][k] = index_index_dis[j][k];
            }

            // establish the path ...
            if (std::get<0>(index_index_dis_mst[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_mst[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_mst[j][k]));
                // auto edge =
                //     add_edge(gind1, gind2, *m_graph);
                    // LogDebug(gind1 << " " << gind2 << " " << std::get<2>(index_index_dis_mst[j][k]));
                float dis;
                // if (edge.second) {
                if (std::get<2>(index_index_dis_mst[j][k]) > 5 * units::cm) {
                    dis = std::get<2>(index_index_dis_mst[j][k]);
                }
                else {
                    dis = std::get<2>(index_index_dis_mst[j][k]);
                }
                // }
                /*auto edge =*/ add_edge(gind1, gind2, WireCell::PointCloud::Facade::EdgeProp(dis), *m_graph);
            }

            if (std::get<0>(index_index_dis_dir_mst[j][k]) >= 0) {
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k]));
                    //auto edge = add_edge(gind1, gind2, *m_graph);
                    // LogDebug(gind1 << " " << gind2 << " " << std::get<2>(index_index_dis_dir1[j][k]));
                    float dis;
                    // if (edge.second) {
                    if (std::get<2>(index_index_dis_dir1[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir1[j][k]) * 1.1;
                    }
                    else {
                        dis = std::get<2>(index_index_dis_dir1[j][k]);
                    }
                    // }
                    /*auto edge =*/ add_edge(gind1, gind2, WireCell::PointCloud::Facade::EdgeProp(dis), *m_graph);
                }
                if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k]));
                    // auto edge = add_edge(gind1, gind2, *m_graph);
                    // LogDebug(gind1 << " " << gind2 << " " << std::get<2>(index_index_dis_dir2[j][k]));
                    // if (edge.second) {
                    float dis;
                    if (std::get<2>(index_index_dis_dir2[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir2[j][k]) * 1.1;
                    }
                    else {
                        dis = std::get<2>(index_index_dis_dir2[j][k]);
                    }
                    // }
                    /*auto edge =*/ add_edge(gind1, gind2, WireCell::PointCloud::Facade::EdgeProp(dis), *m_graph);
                }
            }

        }  // k
    }  // j
    
}


// In Facade_Cluster.cxx
std::vector<int> Cluster::examine_graph(const IDetectorVolumes::pointer dv, const bool use_ctpc) const 
{
    // Create new graph
    if (m_graph != nullptr) {
        m_graph.reset();
    }
    
    m_graph = std::make_unique<MCUGraph>(npoints());
    
    // Establish connections
    Establish_close_connected_graph();
    
    // Connect using overclustering protection (not easy to debug ...)
    Connect_graph_overclustering_protection(dv, use_ctpc); 
    
    // Find connected components
    std::vector<int> component(num_vertices(*m_graph));
    // const int num_components =
    connected_components(*m_graph, &component[0]);

    // std::cout << "Test: num components " << num_components << " " << kd_blobs().size() << std::endl;

    // If only one component, no need for mapping
    // if (num_components <= 1) {
    //     return std::vector<int>();
    // }

    // Create mapping from blob indices to component groups
    std::vector<int> b2groupid(nchildren(), -1);
    
    // For each point in the graph
    for (size_t i = 0; i < component.size(); ++i) {
        // Get the blob index for this point
        const int bind = kd3d().major_index(i);
        // Map the blob to its component
        b2groupid.at(bind) = component[i];
    }

    return b2groupid;
}

void Cluster::dijkstra_shortest_paths(const size_t pt_idx, const bool use_ctpc) const
{
    if (m_graph == nullptr) Create_graph(use_ctpc);
    if ((int)pt_idx == m_source_pt_index) return;
    m_source_pt_index = pt_idx;
    m_parents.resize(num_vertices(*m_graph));
    m_distances.resize(num_vertices(*m_graph));

    vertex_descriptor v0 = vertex(pt_idx, *m_graph);
    // making a param object
    const auto& param = weight_map(get(boost::edge_weight, *m_graph))
				   .predecessor_map(&m_parents[0])
				   .distance_map(&m_distances[0]);
    // const auto& param = boost::weight_map(boost::get(&EdgeProp::dist, *m_graph)).predecessor_map(&m_parents[0]).distance_map(&m_distances[0]);
    boost::dijkstra_shortest_paths(*m_graph, v0, param);
    
    // if (nchildren()==3449){
    //     std::cout << "dijkstra_shortest_paths: " << pt_idx << " " << use_ctpc << std::endl;
    //     std::cout << "distances: ";
    //     for (size_t i = 0; i != m_distances.size(); i++) {
    //         std::cout << i << "->" << m_distances[i] << " ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << "parents: ";
    //     for (size_t i = 0; i != m_parents.size(); i++) {
    //         std::cout << i << "->" << m_parents[i] << " ";
    //     }
    //     std::cout << std::endl;
    // }
}



void Cluster::cal_shortest_path(const size_t dest_wcp_index) const
{
    m_path_wcps.clear();
    m_path_mcells.clear();

    int prev_i = -1;
    for (int i = dest_wcp_index; i != m_source_pt_index; i = m_parents[i])
    {
        const auto* mcell = blob_with_point(i);
        if (m_path_wcps.size() == 0)
        {
            m_path_wcps.push_front(i);
            m_path_mcells.push_front(mcell);
        }
        else
        {
            m_path_wcps.push_front(i);
            if (mcell != m_path_mcells.front())
                m_path_mcells.push_front(mcell);
        }
        if (i == prev_i)
            break;
        prev_i = i;
    }
    auto* src_mcell = blob_with_point(m_source_pt_index);
    m_path_wcps.push_front(m_source_pt_index);
    if (src_mcell != m_path_mcells.front())
        m_path_mcells.push_front(src_mcell);
}

std::vector<geo_point_t> Cluster::indices_to_points(const std::list<size_t>& path_indices) const 
{
    std::vector<geo_point_t> points;
    points.reserve(path_indices.size());
    for (size_t idx : path_indices) {
        points.push_back(point3d(idx));
    }
    return points;
}

void Cluster::organize_points_path_vec(std::vector<geo_point_t>& path_points, double low_dis_limit) const
{
    std::vector<geo_point_t> temp_points = path_points;
    path_points.clear();

    // First pass: filter based on distance
    for (size_t i = 0; i != temp_points.size(); i++) {
        if (path_points.empty()) {
            path_points.push_back(temp_points[i]);
        }
        else if (i + 1 == temp_points.size()) {
            double dis = (temp_points[i] - path_points.back()).magnitude();
            if (dis > low_dis_limit * 0.75) {
                path_points.push_back(temp_points[i]);
            }
        }
        else {
            double dis = (temp_points[i] - path_points.back()).magnitude();
            double dis1 = (temp_points[i + 1] - path_points.back()).magnitude();

            if (dis > low_dis_limit || (dis1 > low_dis_limit * 1.7 && dis > low_dis_limit * 0.75)) {
                path_points.push_back(temp_points[i]);
            }
        }
    }

    // Second pass: filter based on angle
    temp_points = path_points;
    std::vector<double> angles;
    for (size_t i = 0; i != temp_points.size(); i++) {
        if (i == 0 || i + 1 == temp_points.size()) {
            angles.push_back(M_PI);
        }
        else {
            geo_vector_t v1 = temp_points[i] - temp_points[i - 1];
            geo_vector_t v2 = temp_points[i] - temp_points[i + 1];
            angles.push_back(v1.angle(v2));
        }
    }

    path_points.clear();
    for (size_t i = 0; i != temp_points.size(); i++) {
        if (angles[i] * 180.0 / M_PI >= 75) {
            path_points.push_back(temp_points[i]);
        }
    }
}

// this is different from WCP implementation, the path_points is the input ...
void Cluster::organize_path_points(std::vector<geo_point_t>& path_points, double low_dis_limit) const
{
    //    std::vector<geo_point_t> temp_points = path_points;
    path_points.clear();
    auto indices = get_path_wcps();
    auto temp_points = indices_to_points(indices);

    for (size_t i = 0; i != temp_points.size(); i++) {
        if (path_points.empty()) {
            path_points.push_back(temp_points[i]);
        }
        else if (i + 1 == temp_points.size()) {
            double dis = (temp_points[i] - path_points.back()).magnitude();
            if (dis > low_dis_limit * 0.5) {
                path_points.push_back(temp_points[i]);
            }
        }
        else {
            double dis = (temp_points[i] - path_points.back()).magnitude();
            double dis1 = (temp_points[i + 1] - path_points.back()).magnitude();

            if (dis > low_dis_limit || (dis1 > low_dis_limit * 1.7 && dis > low_dis_limit * 0.5)) {
                path_points.push_back(temp_points[i]);
            }
        }
    }
}


std::vector<geo_point_t> Cluster::get_hull() const 
{
    // add cached ...
    if (m_hull_calculated) {
        return m_hull_points;
    }

    quickhull::QuickHull<float> qh;
    std::vector<quickhull::Vector3<float>> pc;
    const auto& points = this->points();
    for (int i = 0; i != npoints(); i++) {
        pc.emplace_back(points[0][i], points[1][i], points[2][i]);
    }
    quickhull::ConvexHull<float> hull = qh.getConvexHull(pc, false, true);
    std::set<int> indices;

    for (size_t i = 0; i != hull.getIndexBuffer().size(); i++) {
        indices.insert(hull.getIndexBuffer().at(i));
    }

    m_hull_points.clear();
    for (auto i : indices) {
        m_hull_points.push_back({points[0][i], points[1][i], points[2][i]});
    }
    
    m_hull_calculated = true;
    return m_hull_points;

    // std::vector<geo_point_t> results;
    // for (auto i : indices) {
    //     results.push_back({points[0][i], points[1][i], points[2][i]});
    // }
    // return results;
}

void Cluster::Calc_PCA() const
{
    if (m_pca_calculated) return;

    m_center.set(0, 0, 0);
    int nsum = 0;
    for (const Blob* blob : children()) {
        for (const geo_point_t& p : blob->points()) {
            m_center += p;
            nsum++;
        }
    }

    for (int i = 0; i != 3; i++) {
        m_pca_axis[i].set(0, 0, 0);
    }

    if (nsum >= 3) {
        m_center = m_center / nsum;
    }
    else {
        return;
    }
    Eigen::MatrixXd cov_matrix(3, 3);

    for (int i = 0; i != 3; i++) {
        for (int j = i; j != 3; j++) {
            cov_matrix(i, j) = 0;
            for (const Blob* blob : children()) {
                for (const geo_point_t& p : blob->points()) {
                    cov_matrix(i, j) += (p[i] - m_center[i]) * (p[j] - m_center[j]);
                }
            }
        }
    }
    cov_matrix(1, 0) = cov_matrix(0, 1);
    cov_matrix(2, 0) = cov_matrix(0, 2);
    cov_matrix(2, 1) = cov_matrix(1, 2);
    // std::cout << cov_matrix << std::endl;

    // const auto eigenSolver = WireCell::Array::pca(cov_matrix);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(cov_matrix);
    auto eigen_values = eigenSolver.eigenvalues();
    auto eigen_vectors = eigenSolver.eigenvectors();

    // ascending order from Eigen, we want descending
    for (int i = 0; i != 3; i++) {
        m_pca_values[2-i] = eigen_values(i);
        double norm = sqrt(eigen_vectors(0, i) * eigen_vectors(0, i) + eigen_vectors(1, i) * eigen_vectors(1, i) +
                             eigen_vectors(2, i) * eigen_vectors(2, i));
        m_pca_axis[2-i].set(eigen_vectors(0, i) / norm, eigen_vectors(1, i) / norm, eigen_vectors(2, i) / norm);
        // std::cout << "PCA: " << i << " " << m_pca_values[i] << " " << m_pca_axis[i] << std::endl;
    }

    m_pca_calculated = true;
}

void Cluster::Calc_PCA(std::vector<geo_point_t>& points) const
{
    // Reset center
    m_center.set(0, 0, 0);
    int nsum = 0;

    // Calculate center
    for (auto it = children().begin(); it != children().end(); it++) {
        for (size_t k = 0; k != points.size(); k++) {
            m_center += points[k];
            nsum++;
        }
    }

    // Reset PCA axes
    for (int i = 0; i != 3; i++) {
        m_pca_axis[i].set(0, 0, 0);
    }

    // Early return if not enough points
    if (nsum < 3) {
        return;
    }

    // Normalize center
    m_center = m_center / nsum;

    // Calculate covariance matrix using Eigen
    Eigen::MatrixXd cov_matrix(3, 3);

    for (int i = 0; i != 3; i++) {
        for (int j = i; j != 3; j++) {
            cov_matrix(i, j) = 0;
            for (auto it = children().begin(); it != children().end(); it++) {
                for (size_t k = 0; k != points.size(); k++) {
                    cov_matrix(i, j) += (points[k][i] - m_center[i]) * (points[k][j] - m_center[j]);
                }
            }
        }
    }

    // Fill symmetric part of matrix
    cov_matrix(1, 0) = cov_matrix(0, 1);
    cov_matrix(2, 0) = cov_matrix(0, 2);
    cov_matrix(2, 1) = cov_matrix(1, 2);

    // Compute eigenvalues and eigenvectors
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(cov_matrix);
    auto eigen_values = eigenSolver.eigenvalues();
    auto eigen_vectors = eigenSolver.eigenvectors();

    // Store eigenvalues and eigenvectors in descending order
    // Note: Eigen returns in ascending order, we want descending
    for (int i = 0; i != 3; i++) {
        m_pca_values[2-i] = eigen_values(i);
        double norm = sqrt(eigen_vectors(0, i) * eigen_vectors(0, i) + 
                         eigen_vectors(1, i) * eigen_vectors(1, i) +
                         eigen_vectors(2, i) * eigen_vectors(2, i));
        
        m_pca_axis[2-i].set(eigen_vectors(0, i) / norm,
                           eigen_vectors(1, i) / norm,
                           eigen_vectors(2, i) / norm);
    }

    m_pca_calculated = true;
}

geo_vector_t Cluster::calc_pca_dir(const geo_point_t& center, const std::vector<geo_point_t>& points) const
{
    // Create covariance matrix
    Eigen::MatrixXd cov_matrix(3, 3);

    // Calculate covariance matrix elements
    for (int i = 0; i != 3; i++) {
        for (int j = i; j != 3; j++) {
            cov_matrix(i, j) = 0;
            for (const auto& p : points) {
                if (i == 0 && j == 0) {
                    cov_matrix(i, j) += (p.x() - center.x()) * (p.x() - center.x());
                }
                else if (i == 0 && j == 1) {
                    cov_matrix(i, j) += (p.x() - center.x()) * (p.y() - center.y());
                }
                else if (i == 0 && j == 2) {
                    cov_matrix(i, j) += (p.x() - center.x()) * (p.z() - center.z());
                }
                else if (i == 1 && j == 1) {
                    cov_matrix(i, j) += (p.y() - center.y()) * (p.y() - center.y());
                }
                else if (i == 1 && j == 2) {
                    cov_matrix(i, j) += (p.y() - center.y()) * (p.z() - center.z());
                }
                else if (i == 2 && j == 2) {
                    cov_matrix(i, j) += (p.z() - center.z()) * (p.z() - center.z());
                }
            }
        }
    }

    // std::cout << "Test: " << center << " " << points.at(0) << std::endl;
    // std::cout << "Test: " << center << " " << points.at(1) << std::endl;
    // std::cout << "Test: " << center << " " << points.at(2) << std::endl;

    // Fill symmetric parts
    cov_matrix(1, 0) = cov_matrix(0, 1);
    cov_matrix(2, 0) = cov_matrix(0, 2);
    cov_matrix(2, 1) = cov_matrix(1, 2);

    // Calculate eigenvalues/eigenvectors using Eigen
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(cov_matrix);
    auto eigen_vectors = eigenSolver.eigenvectors();

    // std::cout << "Test: " << eigen_vectors(0,0) << " " << eigen_vectors(1,0) << " " << eigen_vectors(2,0) << std::endl;

    // Get primary direction (first eigenvector)
    double norm = sqrt(eigen_vectors(0, 2) * eigen_vectors(0, 2) + 
                      eigen_vectors(1, 2) * eigen_vectors(1, 2) + 
                      eigen_vectors(2, 2) * eigen_vectors(2, 2));

    return geo_vector_t(eigen_vectors(0, 2) / norm,
                       eigen_vectors(1, 2) / norm, 
                       eigen_vectors(2, 2) / norm);
}



geo_point_t Cluster::get_center() const {
    if (!m_pca_calculated) {
        Calc_PCA();
    }
    return m_center;
}
geo_vector_t Cluster::get_pca_axis(int axis) const {
    if (!m_pca_calculated) {
        Calc_PCA();
    }
    if (axis < 0 || axis >= 3) raise<IndexError>("axis %d < 0 || axis >= 3", axis);
    return m_pca_axis[axis];
}
double Cluster::get_pca_value(int axis) const {
    if (!m_pca_calculated) {
        Calc_PCA();
    }
    if (axis < 0 || axis >= 3) raise<IndexError>("axis %d < 0 || axis >= 3", axis);
    return m_pca_values[axis];

}


// std::unordered_map<int, Cluster*> 
std::vector<int> Cluster::examine_x_boundary(const double low_limit, const double high_limit)
{
    double num_points[3] = {0, 0, 0};
    double x_max = -1e9;
    double x_min = 1e9;
    auto& mcells = children();
    for (Blob* mcell : mcells) {
        /// TODO: no caching, could be slow
        std::vector<geo_point_t> pts = mcell->points();
        for (size_t i = 0; i != pts.size(); i++) {
            if (pts.at(i).x() < low_limit) {
                num_points[0]++;
                if (pts.at(i).x() > x_max) x_max = pts.at(i).x();
            }
            else if (pts.at(i).x() > high_limit) {
                num_points[2]++;
                if (pts.at(i).x() < x_min) x_min = pts.at(i).x();
            }
            else {
                num_points[1]++;
            }
        }
    }

    // std::cout
    // << "npoints() " << npoints()
    // << " xmax " << x_max << " xmin " << x_min
    // << " low_limit " << low_limit << " high_limit " << high_limit
    // << " num_points: " << num_points[0] << " " << num_points[1] << " " << num_points[2] << std::endl;

    std::vector<Cluster*> clusters;
    std::vector<int> b2groupid(mcells.size(), 0);
    std::set<int> groupids;

    // if (true) {
    if (num_points[0] + num_points[2] < num_points[1] * 0.075) {
        // PR3DCluster* cluster_1 = 0;
        // PR3DCluster* cluster_2 = 0;
        // PR3DCluster* cluster_3 = 0;
        /// FIXME: does tolerance need to be configurable?
        if (x_max < low_limit - 1.0 * units::cm && x_max > -1e8) {
            // fill the small one ...
            // cluster_1 = new PR3DCluster(1);
            groupids.insert(1);
        }
        if (x_min > high_limit + 1.0 * units::cm && x_min < 1e8) {
            // fill the large one ...
            // cluster_3 = new PR3DCluster(3);
            groupids.insert(3);
        }
        // std::cout << "groupids size: " << groupids.size() << std::endl;
        if (!groupids.empty()) {
            // cluster_2 = new PR3DCluster(2);
            groupids.insert(2);
            for (size_t idx=0; idx < mcells.size(); idx++) {
                Blob *mcell = mcells.at(idx);
                if (mcell->points()[0].x() < low_limit) {
                    if (groupids.find(1) != groupids.end()) {
                        // cluster_1->AddCell(mcell, mcell->GetTimeSlice());
                        b2groupid[idx] = 1;
                    }
                    else {
                        // cluster_2->AddCell(mcell, mcell->GetTimeSlice());
                        b2groupid[idx] = 2;
                    }
                }
                else if (mcell->points()[0].x() > high_limit) {
                    if (groupids.find(3) != groupids.end()) {
                        // cluster_3->AddCell(mcell, mcell->GetTimeSlice());
                        b2groupid[idx] = 3;
                    }
                    else {
                        // cluster_2->AddCell(mcell, mcell->GetTimeSlice());
                        b2groupid[idx] = 2;
                    }
                }
                else {
                    // cluster_2->AddCell(mcell, mcell->GetTimeSlice());
                    b2groupid[idx] = 2;
                }
            }
            // if (cluster_1 != 0) clusters.push_back(cluster_1);
            // clusters.push_back(cluster_2);
            // if (cluster_3 != 0) clusters.push_back(cluster_3);
        }
    }
    return b2groupid;
}

bool Cluster::judge_vertex(geo_point_t& p_test, const IDetectorVolumes::pointer dv, const double asy_cut, const double occupied_cut)
{
    p_test = calc_ave_pos(p_test, 3 * units::cm);

    geo_point_t dir = vhough_transform(p_test, 15 * units::cm);

    // judge if this is end points
    std::pair<int, int> num_pts = ndipole(p_test, dir, 25 * units::cm);

    if ((num_pts.first + num_pts.second) == 0) return false;

    double asy = fabs(num_pts.first - num_pts.second) / (num_pts.first + num_pts.second);

    if (asy > asy_cut) {
        return true;
    }
    else {
   
        // it might be better to directly use the closest point to find the wire plane id ...
        auto wpid = dv->contained_by(p_test);
        // what if the point is not found ... 
        if (wpid.apa()==-1){
            auto idx = get_closest_point_index(p_test); 
            // Given the idx, one can directly find the wpid actually ... 
            wpid = dv->contained_by(point3d(idx)); 
        }
         
         // Create wpids for all three planes with the same APA and face
         WirePlaneId wpid_u(kUlayer, wpid.face(), wpid.apa());
         WirePlaneId wpid_v(kVlayer, wpid.face(), wpid.apa());
         WirePlaneId wpid_w(kWlayer, wpid.face(), wpid.apa());
         // Get wire directions for all planes
         Vector wire_dir_u = dv->wire_direction(wpid_u);
         Vector wire_dir_v = dv->wire_direction(wpid_v);
         Vector wire_dir_w = dv->wire_direction(wpid_w);
         // Calculate angles
         double angle_u = std::atan2(wire_dir_u.z(), wire_dir_u.y());
         double angle_v = std::atan2(wire_dir_v.z(), wire_dir_v.y());
         double angle_w = std::atan2(wire_dir_w.z(), wire_dir_w.y());

        auto temp_point_cloud = std::make_shared<Multi2DPointCloud>(angle_u, angle_v, angle_w);
        dir = dir.norm();
        // PointVector pts;
        std::vector<geo_point_t> pts;
        for (size_t i = 0; i != 40; i++) {
            geo_point_t pt(p_test.x() + i * 0.5 * units::cm * dir.x(), p_test.y() + i * 0.5 * units::cm * dir.y(),
                     p_test.z() + i * 0.5 * units::cm * dir.z());
            // WCP::WCPointCloud<double>::WCPoint& wcp = point_cloud->get_closest_wcpoint(pt);
            auto [_, wcp] = get_closest_wcpoint(pt);

            if (sqrt(pow(wcp.x() - pt.x(), 2) + pow(wcp.y() - pt.y(), 2) + pow(wcp.z() - pt.z(), 2)) <
                std::max(1.8 * units::cm, i * 0.5 * units::cm * sin(18. / 180. * 3.1415926))) {
                pt = wcp;
            }
            pts.push_back(pt);
            if (i != 0) {
                geo_point_t pt1(p_test.x() - i * 0.5 * units::cm * dir.x(), p_test.y() - i * 0.5 * units::cm * dir.y(),
                          p_test.z() - i * 0.5 * units::cm * dir.z());
                // WCP::WCPointCloud<double>::WCPoint& wcp1 = point_cloud->get_closest_wcpoint(pt1);
                auto [_, wcp1] = get_closest_wcpoint(pt1);
                if (sqrt(pow(wcp1.x() - pt1.x(), 2) + pow(wcp1.y() - pt1.y(), 2) + pow(wcp1.z() - pt1.z(), 2)) <
                    std::max(1.8 * units::cm, i * 0.5 * units::cm * sin(18. / 180. * 3.1415926))) {
                    pt1 = wcp1;
                }
                pts.push_back(pt1);
            }
        }
        // temp_point_cloud.AddPoints(pts);
        for (auto& pt : pts) {
            temp_point_cloud->add(pt);
        }
        // temp_point_cloud.build_kdtree_index();

        int temp_num_total_points = 0;
        int temp_num_occupied_points = 0;

        // const int N = point_cloud->get_num_points();
        const int N = npoints();
        // WCP::WCPointCloud<double>& cloud = point_cloud->get_cloud();
        for (int i = 0; i != N; i++) {
            // geo_point_t dir1(cloud.pts[i].x() - p_test.x(), cloud.pts[i].y() - p_test.y(), cloud.pts[i].z() - p_test.z());
            geo_point_t dir1 = point3d(i) - p_test;

            if (dir1.magnitude() < 15 * units::cm) {
                geo_point_t test_p1 = point3d(i);
                temp_num_total_points++;
                double dis[3];
                dis[0] = temp_point_cloud->get_closest_2d_dis(test_p1, 0).second;
                dis[1] = temp_point_cloud->get_closest_2d_dis(test_p1, 1).second;
                dis[2] = temp_point_cloud->get_closest_2d_dis(test_p1, 2).second;
                if (dis[0] <= 1.5 * units::cm && dis[1] <= 1.5 * units::cm && dis[2] <= 2.4 * units::cm ||
                    dis[0] <= 1.5 * units::cm && dis[2] <= 1.5 * units::cm && dis[1] <= 2.4 * units::cm ||
                    dis[2] <= 1.5 * units::cm && dis[1] <= 1.5 * units::cm && dis[0] <= 2.4 * units::cm)
                    temp_num_occupied_points++;
            }
        }

        if (temp_num_occupied_points < temp_num_total_points * occupied_cut) return true;
    }

    // judge if there

    return false;
}

bool Facade::cluster_less(const Cluster* a, const Cluster* b)
{
    if (a == b) return false;

    {
        const double la = a->get_length();
        const double lb = b->get_length();
        if (la < lb) return true;
        if (lb < la) return false;
    }
    {
        const int na = a->nchildren();
        const int nb = b->nchildren();
        if (na < nb) return true;
        if (nb < na) return false;
    }
    {
        const int na = a->npoints();
        const int nb = b->npoints();
        if (na < nb) return true;
        if (nb < na) return false;
    }
    {
        auto ar = a->get_uvwt_min();
        auto br = b->get_uvwt_min();
        if (get<0>(ar) < get<0>(br)) return true;
        if (get<0>(br) < get<0>(ar)) return false;
        if (get<1>(ar) < get<1>(br)) return true;
        if (get<1>(br) < get<1>(ar)) return false;
        if (get<2>(ar) < get<2>(br)) return true;
        if (get<2>(br) < get<2>(ar)) return false;
        if (get<3>(ar) < get<3>(br)) return true;
        if (get<3>(br) < get<3>(ar)) return false;
    }
    {
        auto ar = a->get_uvwt_max();
        auto br = b->get_uvwt_max();
        if (get<0>(ar) < get<0>(br)) return true;
        if (get<0>(br) < get<0>(ar)) return false;
        if (get<1>(ar) < get<1>(br)) return true;
        if (get<1>(br) < get<1>(ar)) return false;
        if (get<2>(ar) < get<2>(br)) return true;
        if (get<2>(br) < get<2>(ar)) return false;
        if (get<3>(ar) < get<3>(br)) return true;
        if (get<3>(br) < get<3>(ar)) return false;
    }

    // The two are very similar.  What is left to check?  Only pointer?.  This
    // will cause "randomness"!
    return a < b;
}
void Facade::sort_clusters(std::vector<const Cluster*>& clusters)
{
    std::sort(clusters.rbegin(), clusters.rend(), cluster_less);
}
void Facade::sort_clusters(std::vector<Cluster*>& clusters)
{
    std::sort(clusters.rbegin(), clusters.rend(), cluster_less);
}


Facade::Cluster::Flash Facade::Cluster::get_flash() const
{
    Flash flash;                // starts invalid

    const auto* p = node()->parent;
    if (!p)  return flash;
    const auto* g = p->value.facade<Grouping>();
    if (!g)  return flash;

    const int flash_index = get_scalar("flash", -1);

    //std::cout << "Test3 " << flash_index << std::endl;
    
    if (flash_index < 0) {
        return flash;
    }
    if (! g->has_pc("flash")) {
        return flash;
    }
    flash.m_valid = true;
        
    // These are kind of inefficient as we get the "flash" PC each time.
    flash.m_time = g->get_element<double>("flash", "time", flash_index, 0);
    flash.m_value = g->get_element<double>("flash", "value", flash_index, 0);
    flash.m_ident = g->get_element<int>("flash", "ident", flash_index, -1);
    flash.m_type = g->get_element<int>("flash", "type", flash_index, -1);

    // std::cout << "Test3: " << g->has_pc("flash") << " " << g->has_pc("light") << " " << g->has_pc("flashlight") << " " << flash_index << " " << flash.m_time << std::endl;

    if (!(g->has_pc("light") && g->has_pc("flashlight"))) {
        return flash;           // valid, but no vector info.
    }
    
    // These are spans.  We walk the fl to look up in the l.
    const auto fl_flash = g->get_pcarray<int>("flash", "flashlight");
    const auto fl_light = g->get_pcarray<int>("light", "flashlight");
    const auto l_times = g->get_pcarray<double>("time", "light");
    const auto l_values = g->get_pcarray<double>("value", "light");
    const auto l_errors = g->get_pcarray<double>("error", "light");

    // std::cout << "Test3: " << fl_flash.size() << " " << fl_light.size() << std::endl;

    const size_t nfl = fl_light.size();
    for (size_t ifl = 0; ifl < nfl; ++ifl) {
        if (fl_flash[ifl] != flash_index) continue;
        const int light_index = fl_light[ifl];
        
        flash.m_times.push_back(l_times[light_index]);
        flash.m_values.push_back(l_values[light_index]);
        flash.m_errors.push_back(l_errors[light_index]);
    }
    return flash;
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
