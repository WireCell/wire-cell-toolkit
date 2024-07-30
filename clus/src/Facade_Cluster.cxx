#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Grouping.h"

#include "WireCellUtil/Array.h"

#include <boost/container_hash/hash.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Facade;
// using WireCell::PointCloud::Dataset;
using namespace WireCell::PointCloud::Tree;  // for "Points" node value type
// using WireCell::PointCloud::Tree::named_pointclouds_t;

#include "WireCellUtil/Logging.h"
using spdlog::debug;

// #define __DEBUG__
#ifdef __DEBUG__
#define LogDebug(x) std::cout << "[yuhw]: " << __LINE__ << " : " << x << std::endl
#else
#define LogDebug(x)
#endif

std::ostream& Facade::operator<<(std::ostream& os, const Facade::Cluster& cluster)
{
    
    os << "<Cluster [" << (void*) cluster.hash() << "]:" << " npts=" << cluster.npoints()
       << " nblobs=" << cluster.nchildren() << ">";
    return os;
}

Grouping* Cluster::grouping() { return this->m_node->parent->value.template facade<Grouping>(); }
const Grouping* Cluster::grouping() const { return this->m_node->parent->value.template facade<Grouping>(); }

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
    << " " << u_min << " " << u_max << " " << v_min << " " << v_max << " " << w_min << " " << w_max;
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
    geo_point_t drift_dir(1, 0, 0);

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

        double angle_1 = fabs(dir1.angle(drift_dir) - 3.1415926 / 2.) * 180. / 3.1415926;
        double angle_2 = fabs(dir.angle(drift_dir) - 3.1415926 / 2.) / 3.1415926 * 180.;
        double angle_3 = fabs(dir2.angle(drift_dir) - 3.1415926 / 2.) / 3.1415926 * 180.;
        double angle_4 = fabs(dir3.angle(drift_dir) - 3.1415926 / 2.) / 3.1415926 * 180.;

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
                double angle_1 = fabs(dir1.angle(drift_dir) - 3.1415926 / 2.) * 180. / 3.1415926;
                double angle_2 = fabs(dir.angle(drift_dir) - 3.1415926 / 2.) / 3.1415926 * 180.;
                double angle_3 = fabs(dir2.angle(drift_dir) - 3.1415926 / 2.) / 3.1415926 * 180.;
                double angle_4 = fabs(dir3.angle(drift_dir) - 3.1415926 / 2.) / 3.1415926 * 180.;

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

    for (size_t pt_idx = 0; pt_idx != npoints(); pt_idx++) {
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


const Cluster::sv2d_t& Cluster::sv2d(const size_t plane) const {
    return m_node->value.scoped_view(scope2ds[plane]);
    }

const Cluster::kd2d_t& Cluster::kd2d(const size_t plane) const
{
    const auto& sv = sv2d(plane);
    return sv.kd();
}

std::vector<size_t> Cluster::get_closest_2d_index(const geo_point_t& p, const double search_radius, const int plane) const {

    const auto& tp = grouping()->get_params();
    double angle_uvw[3] = {tp.angle_u, tp.angle_v, tp.angle_w};
    double x = p.x();
    double y = cos(angle_uvw[plane]) * p.z() - sin(angle_uvw[plane]) * p.y();
    std::vector<float_t> query_pt = {x, y};
    const auto& skd = kd2d(plane);
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

const Cluster::sv3d_t& Cluster::sv3d() const { return m_node->value.scoped_view(scope); }

const Cluster::kd3d_t& Cluster::kd3d() const
{
    const auto& sv = m_node->value.scoped_view(scope);
    return sv.kd();
}
geo_point_t Cluster::point3d(size_t point_index) const { return kd3d().point3d(point_index); }
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

std::pair<int, int> Cluster::ndipole(const geo_point_t& point, const geo_point_t& dir) const
{
    const auto& points = this->points();
    const size_t npoints = points[0].size();

    int num_p1 = 0;
    int num_p2 = 0;

    for (size_t ind = 0; ind < npoints; ++ind) {
        geo_point_t dir1(points[0][ind] - point.x(), points[1][ind] - point.y(), points[2][ind] - point.z());
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
//     const auto& sv = m_node->value.scoped_view(scope);       // get the kdtree
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

std::pair<geo_point_t, const Blob*> Cluster::get_closest_point_blob(const geo_point_t& point) const
{
    auto results = kd_knn(1, point);
    if (results.size() == 0) {
        return std::make_pair(geo_point_t(), nullptr);
    }

    const auto& [point_index, _] = results[0];
    return std::make_pair(point3d(point_index), blob_with_point(point_index));
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
std::tuple<int, int, int, int> Cluster::get_uvwt_range() const
{
    std::set<int> u_set;
    std::set<int> v_set;
    std::set<int> w_set;
    std::set<int> t_set;
    for (const auto* blob : children()) {
        for (int i = blob->u_wire_index_min(); i < blob->u_wire_index_max(); ++i) {
            u_set.insert(i);
        }
        for (int i = blob->v_wire_index_min(); i < blob->v_wire_index_max(); ++i) {
            v_set.insert(i);
        }
        for (int i = blob->w_wire_index_min(); i < blob->w_wire_index_max(); ++i) {
            w_set.insert(i);
        }
        for (int i = blob->slice_index_min(); i < blob->slice_index_max(); ++i) {
            t_set.insert(i);
        }
    }
    return {u_set.size(), v_set.size(), w_set.size(), t_set.size()};
}

double Cluster::get_length() const
{
    if (m_length == 0) {  // invalidates when a new node is set
        const auto& tp = grouping()->get_params();

        const auto [u, v, w, t] = get_uvwt_range();
        const double pu = u * tp.pitch_u;
        const double pv = v * tp.pitch_v;
        const double pw = w * tp.pitch_w;
        const double pt = t * tp.tick_drift;
        m_length = std::sqrt(2. / 3. * (pu * pu + pv * pv + pw * pw) + pt * pt);
    }
    return m_length;
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

bool Cluster::sanity(Log::logptr_t log) const
{
    {
        const auto* svptr = m_node->value.get_scoped(scope);
        if (!svptr) {
            if (log) log->debug("cluster sanity: note, not yet a scoped view {}", scope);
        }
    }
    if (!nchildren()) {
        if (log) log->debug("cluster sanity: no children blobs");
        return false;
    }

    const auto& sv = m_node->value.scoped_view(scope);
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

std::vector<int> Cluster::get_blob_indices(const Blob* blob)
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
void Cluster::Create_graph(const bool use_ctpc)
{
    LogDebug("Create Graph! " << graph);
    if (graph != (MCUGraph*) 0) return;
    graph = new MCUGraph(nbpoints());
    Establish_close_connected_graph();
    if (use_ctpc) Connect_graph(true);
    Connect_graph(false);
}

void Cluster::Establish_close_connected_graph()
{
    std::map<const Blob*, std::map<int, std::set<int>>> map_mcell_uindex_wcps;
    std::map<const Blob*, std::map<int, std::set<int>>> map_mcell_vindex_wcps;
    std::map<const Blob*, std::map<int, std::set<int>>> map_mcell_windex_wcps;

    std::map<Blob*, std::set<int>> map_mcell_indices;

    const auto& points = this->points();
    const auto& winds = this->wire_indices();
    LogDebug("points[0].size(): " << points[0].size() << " winds[0].size(): " << winds[0].size());

    for (Blob* mcell : this->children()) {
        std::map<int, std::set<int>> map_uindex_wcps;
        std::map<int, std::set<int>> map_vindex_wcps;
        std::map<int, std::set<int>> map_windex_wcps;

        std::vector<int> pinds = this->get_blob_indices(mcell);
        for (const int pind : pinds) {
            auto v = vertex(pind, *graph);  // retrieve vertex descriptor
            (*graph)[v].index = pind;
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
                        auto edge = add_edge(pind1, pind2, *graph);
                        if (edge.second) {
                            (*graph)[edge.first].dist = sqrt(pow(points[0][pind1] - points[0][pind2], 2) +
                                                             pow(points[1][pind1] - points[1][pind2], 2) +
                                                             pow(points[2][pind1] - points[2][pind2], 2));
                            num_edges++;
                        }
                    }
                }
            }
        }
    }

    LogDebug("in-blob edges: " << num_edges);

    std::vector<int> time_slices;
    for (auto [time, _] : this->time_blob_map()) {
        time_slices.push_back(time);
    }
    const int nticks_per_slice = grouping()->get_params().nticks_live_slice;
    // std::cout << "time_slices size: " << time_slices.size() << std::endl;

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
    std::map<std::pair<int, int>, std::pair<int, double>> closest_index;

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

                        if (closest_index.find(key) == closest_index.end()) {
                            closest_index[key] = std::make_pair(pind2, dis);
                        }
                        else {
                            if (dis < closest_index[key].second) closest_index[key] = std::make_pair(pind2, dis);
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
                            closest_index[key] = std::make_pair(pind2, dis);
                        }
                        else {
                            if (dis < closest_index[key].second) closest_index[key] = std::make_pair(pind2, dis);
                        }
                    }
                }
            }
        }
    }

    for (auto it4 = closest_index.begin(); it4 != closest_index.end(); it4++) {
        int index1 = it4->first.first;
        int index2 = it4->second.first;
        double dis = it4->second.second;
        auto edge = add_edge(index1, index2, *graph);
        if (edge.second) {
            (*graph)[edge.first].dist = dis;
            num_edges++;
        }
    }
    // end of copying ...

    LogDebug("all edges: " << num_edges);
}

void Cluster::Connect_graph(const bool use_ctpc) {
    // now form the connected components
    std::vector<int> component(num_vertices(*graph));
    const size_t num = connected_components(*graph, &component[0]);
    LogDebug(" npoints " << npoints() << " nconnected " << num);
    if (num <= 1) return;

    std::vector<std::shared_ptr<Simple3DPointCloud>> pt_clouds;
    // use this to link the global index to the local index
    std::vector<std::vector<size_t>> pt_clouds_global_indices(num);
    for (size_t i = 0; i != num; i++) {
        pt_clouds.push_back(std::make_shared<Simple3DPointCloud>());
    }
    for (size_t i = 0; i != component.size(); ++i) {
        pt_clouds.at(component[i])->add({points()[0][i], points()[1][i], points()[2][i]});
        pt_clouds_global_indices.at(component[i]).push_back(i);
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
                        /// FIXME: how to add face information?
                        const int face = 0;
                        const bool good_point = grouping()->is_good_point(test_p, face);
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
                        /// FIXME: how to add face information?
                        const int face = 0;
                        const bool good_point = grouping()->is_good_point(test_p, face);
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
                        /// FIXME: how to add face information?
                        const int face = 0;
                        const bool good_point = grouping()->is_good_point(test_p, face);
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

        {
            std::vector<int> possible_root_vertex;
            std::vector<int> component(num_vertices(temp_graph));
            const int num1 = connected_components(temp_graph, &component[0]);
            possible_root_vertex.resize(num1);
            std::vector<int>::size_type i;
            for (i = 0; i != component.size(); ++i) {
                possible_root_vertex.at(component[i]) = i;
            }

            for (size_t i = 0; i != possible_root_vertex.size(); i++) {
                std::vector<boost::graph_traits<MCUGraph>::vertex_descriptor> predecessors(num_vertices(temp_graph));

                prim_minimum_spanning_tree(temp_graph, &predecessors[0],
                                           boost::root_vertex(possible_root_vertex.at(i)));

                for (size_t j = 0; j != predecessors.size(); ++j) {
                    if (predecessors[j] != j) {
                        if (j < predecessors[j]) {
                            index_index_dis_mst[j][predecessors[j]] = index_index_dis[j][predecessors[j]];
                        }
                        else {
                            index_index_dis_mst[predecessors[j]][j] = index_index_dis[predecessors[j]][j];
                        }
                        // std::cout << j << " " << predecessors[j] << " " << std::endl;
                    }
                    else {
                        // std::cout << j << " " << std::endl;
                    }
                }
            }
        }
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

        {
            std::vector<int> possible_root_vertex;
            std::vector<int> component(num_vertices(temp_graph));
            const int num1 = connected_components(temp_graph, &component[0]);
            possible_root_vertex.resize(num1);
            std::vector<int>::size_type i;
            for (i = 0; i != component.size(); ++i) {
                possible_root_vertex.at(component[i]) = i;
            }

            for (size_t i = 0; i != possible_root_vertex.size(); i++) {
                std::vector<boost::graph_traits<MCUGraph>::vertex_descriptor> predecessors(num_vertices(temp_graph));
                prim_minimum_spanning_tree(temp_graph, &predecessors[0],
                                           boost::root_vertex(possible_root_vertex.at(i)));
                for (size_t j = 0; j != predecessors.size(); ++j) {
                    if (predecessors[j] != j) {
                        if (j < predecessors[j]) {
                            index_index_dis_dir_mst[j][predecessors[j]] = index_index_dis[j][predecessors[j]];
                        }
                        else {
                            index_index_dis_dir_mst[predecessors[j]][j] = index_index_dis[predecessors[j]][j];
                        }
                        // std::cout << j << " " << predecessors[j] << " " << std::endl;
                    }
                    else {
                        // std::cout << j << " " << std::endl;
                    }
                }
            }
        }
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
                auto edge =
                    add_edge(gind1, gind2, *graph);
                    // LogDebug(gind1 << " " << gind2 << " " << std::get<2>(index_index_dis_mst[j][k]));
                if (edge.second) {
                    if (std::get<2>(index_index_dis_mst[j][k]) > 5 * units::cm) {
                        (*graph)[edge.first].dist = std::get<2>(index_index_dis_mst[j][k]);
                    }
                    else {
                        (*graph)[edge.first].dist = std::get<2>(index_index_dis_mst[j][k]);
                    }
                }
            }

            if (std::get<0>(index_index_dis_dir_mst[j][k]) >= 0) {
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k]));
                    auto edge = add_edge(gind1, gind2, *graph);
                    // LogDebug(gind1 << " " << gind2 << " " << std::get<2>(index_index_dis_dir1[j][k]));
                    if (edge.second) {
                        if (std::get<2>(index_index_dis_dir1[j][k]) > 5 * units::cm) {
                            (*graph)[edge.first].dist = std::get<2>(index_index_dis_dir1[j][k]) * 1.1;
                        }
                        else {
                            (*graph)[edge.first].dist = std::get<2>(index_index_dis_dir1[j][k]);
                        }
                    }
                }
                if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k]));
                    auto edge = add_edge(gind1, gind2, *graph);
                    // LogDebug(gind1 << " " << gind2 << " " << std::get<2>(index_index_dis_dir2[j][k]));
                    if (edge.second) {
                        if (std::get<2>(index_index_dis_dir2[j][k]) > 5 * units::cm) {
                            (*graph)[edge.first].dist = std::get<2>(index_index_dis_dir2[j][k]) * 1.1;
                        }
                        else {
                            (*graph)[edge.first].dist = std::get<2>(index_index_dis_dir2[j][k]);
                        }
                    }
                }
            }

        }  // k
    }  // j
}
// #define LogDebug(x)

void Cluster::dijkstra_shortest_paths(const size_t pt_idx, const bool use_ctpc)
{
    if (graph == (MCUGraph*) 0) Create_graph(use_ctpc);
    if (pt_idx == m_source_pt_index) return;
    m_source_pt_index = pt_idx;
    m_parents.resize(num_vertices(*graph));
    m_distances.resize(num_vertices(*graph));

    vertex_descriptor v0 = vertex(pt_idx, *graph);
    // making a param object
    const auto& param = boost::weight_map(boost::get(&EdgeProp::dist, *graph)).predecessor_map(&m_parents[0]).distance_map(&m_distances[0]);
    boost::dijkstra_shortest_paths(*graph, v0, param);
    // std::cout << "dijkstra_shortest_paths: " << pt_idx << " " << use_ctpc << std::endl;
    // std::cout << "distances: ";
    // for (size_t i = 0; i != m_distances.size(); i++) {
    //     std::cout << i << "->" << m_distances[i] << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "parents: ";
    // for (size_t i = 0; i != m_parents.size(); i++) {
    //     std::cout << i << "->" << m_parents[i] << " ";
    // }
    // std::cout << std::endl;
}



void Cluster::cal_shortest_path(const size_t dest_wcp_index)
{
    m_path_wcps.clear();
    m_path_mcells.clear();

    int prev_i = -1;
    for (int i = dest_wcp_index; i != m_source_pt_index; i = m_parents[i])
    {
        auto* mcell = blob_with_point(i);
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

std::vector<geo_point_t> Cluster::get_hull() const 
{
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

    std::vector<geo_point_t> results;
    for (auto i : indices) {
        results.push_back({points[0][i], points[1][i], points[2][i]});
    }
    return results;
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

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
