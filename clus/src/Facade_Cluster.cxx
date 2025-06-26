#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Graphs.h"

#include "WireCellUtil/Array.h"

#include <boost/container_hash/hash.hpp>

#include "WireCellUtil/Logging.h"

#include "make_graphs.h"

// The original developers do not care.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"


using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Graphs;
using namespace WireCell::PointCloud;
using namespace WireCell::Clus::Facade;
// using WireCell::PointCloud::Dataset;
using namespace WireCell::PointCloud::Tree;  // for "Points" node value type
// using WireCell::PointCloud::Tree::named_pointclouds_t;
using WireCell::Clus::Graphs::Weighted::GraphAlgorithms;

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

Grouping* Cluster::grouping()
{
    return this->m_node->parent->value.template facade<Grouping>();
}

const Grouping* Cluster::grouping() const
{ 
    return this->m_node->parent->value.template facade<Grouping>();
}

void Cluster::set_default_scope(const Tree::Scope& scope)
{
    // We can not simply return if scope is unchanged as that will cause a
    // crash in connect_graph_closely() functions due to bad map_mcell_* lookup.
    //
    // if (m_default_scope == scope) {
    //     return;
    // }

    m_default_scope = scope;

    // Clear caches that depend on the scope
    clear_cache(); // Why is this here???  It does not do what the comment says.
                   // It clears all cache.  This side-effect is needed even if
                   // the default scope is unchanged.

    

    // The PCA cache is only thing that directly depends on scope but it is not
    // enough to just clear that...
    // cache().pca.reset();
    // ... as connect_graph_closely() still breaks.
    // For now, we leave the mystery unsolved.

}

void Cluster::set_scope_filter(const Tree::Scope& scope, bool flag)
{
    // Set the scope filter for the given scope
    m_map_scope_filter[scope.hash()] = flag;
}

bool Cluster::get_scope_filter(const Tree::Scope& scope) const
{
    auto it = m_map_scope_filter.find(scope.hash());
    if (it == m_map_scope_filter.end()){
        return false;
    }
    return it->second;
}


void Cluster::set_scope_transform(const Tree::Scope& scope, const std::string& transform_name)
{
    // Set the scope transform for the given scope
    m_map_scope_transform[scope.hash()] = transform_name;
}

std::string Cluster::get_scope_transform(Tree::Scope scope) const
{
    Scope the_scope;
    if (scope == the_scope) { // no scope given
        scope = m_default_scope;
    }
    auto it = m_map_scope_transform.find(scope.hash());
    if (it == m_map_scope_transform.end()){
        return "Unity";
    }
    return it->second;
}

const Tree::Scope& Cluster::get_scope(const std::string& scope_name) const
{
    if (m_scopes.find(scope_name) == m_scopes.end()) {
        raise<RuntimeError>("Cluster::scope: no such scope: %s", scope_name);
    }
    return m_scopes.at(scope_name);
}


void Cluster::set_cluster_id(int cid)
{
    this->set_ident(cid);
}

int Cluster::get_cluster_id() const
{
    return this->ident();
}

void Cluster::default_scope_from(const Cluster& other)
{
    auto scope = other.get_default_scope();
    this->set_default_scope(scope);
    if (other.get_scope_filter(scope)) {
        this->set_scope_filter(scope, other.get_scope_filter(scope));
    }
    this->set_scope_transform(scope, other.get_scope_transform(scope));
}

void Cluster::from(const Cluster& other)
{
    this->default_scope_from(other);
    this->flags_from(other);
}


void Cluster::set_cluster_t0(double t0)
{
    this->set_scalar<double>("cluster_t0", t0);
}
double Cluster::get_cluster_t0() const
{
    return this->get_scalar<double>("cluster_t0", 0);
}


std::vector<int> Cluster::add_corrected_points(
    Clus::IPCTransformSet::pointer pcts,
    const std::string &correction_name) 
{
    const double t0 = this->get_cluster_t0();

    std::vector<int> blob_passed;
    blob_passed.resize(children().size(), 0); // not passed by default
    if (correction_name == "T0Correction") {
        const auto& pct = pcts->pc_transform("T0Correction");
        for (size_t iblob = 0; iblob < this->children().size(); ++iblob) {
            Blob* blob = this->children().at(iblob);
            auto &lpc_3d = blob->local_pcs().at("3d");
            auto corrected_points = pct->forward(lpc_3d, {"x", "y", "z"},
                                                 {"x_t0cor","y_t0cor","z_t0cor"}, t0,
                                                 blob->wpid().face(), blob->wpid().apa());
            lpc_3d.add("x_t0cor", *corrected_points.get("x_t0cor")); // only add x_t0cor
            auto filter_result = pct->filter(corrected_points,
                                             {"x_t0cor", "y_t0cor", "z_t0cor"},
                                             t0, blob->wpid().face(), blob->wpid().apa());
            auto arr_filter = filter_result.get("filter")->elements<int>();
            for (size_t ipt = 0; ipt < arr_filter.size(); ++ipt) {
                if (arr_filter[ipt] == 1) {
                    blob_passed[iblob] = 1;
                    break; // only one point pass is enough
                }
            }
        }
        // the new scope should have the same name as the correction name. This is how the code can find corrections in the code ...
        m_scopes["T0Correction"] = {"3d", {"x_t0cor", "y", "z"}}; // add the new scope
    } else {
        raise<RuntimeError>("Cluster::add_corrected_points: no such correction: %s", correction_name);
    }
    return blob_passed;
}



// Called first time cache() is called and the cache is invalid.
void Cluster::fill_cache(ClusterCache& cache) const
{
    // There is nothing generic to "pre fill".  Instead, each individual method
    // will fill the cache as needed.
}

// blob wpids ...
std::vector<WireCell::WirePlaneId> Cluster::wpids_blob() const
{
    auto& wpids = cache().blob_wpids;
    if (wpids.empty()) {
        for (const Blob* blob : this->children()) {
            wpids.push_back(blob->wpid());
        }
    }
    return wpids;
}

WirePlaneId Cluster::wpid(const geo_point_t& point) const
{
    // find the closest point_index to this point
    auto point_index = get_closest_point_index(point);

    // std::cout << "point_index " << point_index << " " << points()[0].size() << " " << wpids().size() << std::endl;

    // return the wpid for this point_index
    return wire_plane_id(point_index);
}


void Cluster::print_blobs_info() const
{
    for (const Blob* blob : children()) {
        std::cout << "U: " << blob->u_wire_index_min() << " " << blob->u_wire_index_max() 
        << " V: " << blob->v_wire_index_min() << " " << blob->v_wire_index_max() 
        << " W: " << blob->w_wire_index_min() << " " << blob->w_wire_index_max() 
        << " T: " << blob->slice_index_min() << " " << blob->slice_index_max()
        << std::endl;


    }
}

std::string Cluster::dump() const
{
    const auto [u_min, v_min, w_min, t_min] = get_uvwt_min();
    const auto [u_max, v_max, w_max, t_max] = get_uvwt_max();
    std::stringstream ss;
    ss << " blobs " << children().size() << " points " << npoints()
    << " [" << t_min << " " << t_max << "] " << children().size()
    << " uvw " << u_min << " " << u_max << " " << v_min << " " << v_max << " " << w_min << " " << w_max;
    return ss.str();
}

const Cluster::time_blob_map_t& Cluster::time_blob_map() const
{
    auto& tbm = cache().time_blob_map;
    if (tbm.empty()) {
        for (const Blob* blob : children()) {
            auto wpid = blob->wpid();
            tbm[wpid.apa()][wpid.face()][blob->slice_index_min()].insert(blob);
        }
    }
    return tbm;
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

// This function works with raw points internally, different from most of the functions ...
void Cluster::adjust_wcpoints_parallel(size_t& start_idx, size_t& end_idx) const
{
    const auto& winds = wire_indices();

    geo_point_t start_p = point3d_raw(start_idx);
    geo_point_t end_p = point3d_raw(end_idx);

    WirePlaneId start_wpid = wire_plane_id(start_idx);
    WirePlaneId end_wpid = wire_plane_id(end_idx);

    double low_x = start_p.x() - 1 * units::cm;
    if (end_p.x() - 1 * units::cm < low_x) low_x = end_p.x() - 1 * units::cm;
    double high_x = start_p.x() + 1 * units::cm;
    if (end_p.x() + 1 * units::cm > high_x) high_x = end_p.x() + 1 * units::cm;


    // Create map to track lowest wire indices for each wire plane ID
    std::map<WirePlaneId, std::array<size_t, 3>> map_wpid_low_indices;
    std::map<WirePlaneId, std::array<size_t, 3>> map_wpid_high_indices;

    // Initialize all elements in the arrays
    map_wpid_low_indices[start_wpid] = {start_idx, start_idx, start_idx};
    map_wpid_high_indices[start_wpid] = {start_idx, start_idx, start_idx};
    if (end_wpid != start_wpid){
        map_wpid_low_indices[end_wpid] = {end_idx,end_idx,end_idx};
        map_wpid_high_indices[end_wpid] = {end_idx,end_idx,end_idx};
    }

    

    // assumes u, v, w, need to expand to includ wpid ???
    for (int pt_idx = 0; pt_idx != npoints(); pt_idx++) {
        geo_point_t current = point3d_raw(pt_idx);
        WirePlaneId wpid = wire_plane_id(pt_idx);
        // WirePlaneId wpid = start_wpid;

        if (pt_idx % 1000 == 0)
        // std::cout << "Test: " << pt_idx << " " << wpid << npoints() << std::endl;

        if (current.x() > high_x || current.x() < low_x) continue;

        if (map_wpid_low_indices.find(wpid) == map_wpid_low_indices.end()) {
            for (size_t pind = 0; pind != 3; ++pind) {
                map_wpid_low_indices[wpid][pind] = pt_idx;
            }
        }else {
            for (size_t pind = 0; pind != 3; ++pind) {
                if (winds[pind][pt_idx] < winds[pind][map_wpid_low_indices[wpid][pind]]) {
                    map_wpid_low_indices[wpid][pind] = pt_idx;
                }
            }
        }
        if(map_wpid_high_indices.find(wpid) == map_wpid_high_indices.end()) {
            for (size_t pind = 0; pind != 3; ++pind) {
                map_wpid_high_indices[wpid][pind] = pt_idx;
            }
        }else {
            for (size_t pind = 0; pind != 3; ++pind) {
                if (winds[pind][pt_idx] > winds[pind][map_wpid_high_indices[wpid][pind]]) {
                    map_wpid_high_indices[wpid][pind] = pt_idx;
                }
            }
        }  
    }

    bool flags[3] = {true, true, true};
    {
        // Calculate the size of the range for each wire plane across all WPIDs
        int index_diff_sum[3] = {0, 0, 0};
        // Find minimum and maximum indices for each plane across all WPIDs
        for (auto it = map_wpid_low_indices.begin(); it != map_wpid_low_indices.end(); ++it) {
            const WirePlaneId& wpid = it->first;
            const auto& low_indices = it->second;
            const auto& high_indices = map_wpid_high_indices[wpid];
            
            for (size_t pind = 0; pind < 3; ++pind) {
                index_diff_sum[pind] += winds[pind][high_indices[pind]] - winds[pind][low_indices[pind]];
            }
        }

        // Create pairs of (index_difference, plane_index) for sorting
        std::vector<std::pair<int, int>> plane_diffs;
        for (int i = 0; i < 3; ++i) {
            plane_diffs.push_back({index_diff_sum[i], i});
        }
        // Sort by index difference (ascending)
        std::sort(plane_diffs.begin(), plane_diffs.end());
        // Set flag to false for the plane with smallest difference
        // (keeping the two planes with largest differences)
        flags[plane_diffs[0].second] = false;
    }
    

    std::vector<size_t> indices, temp_indices;
    std::set<size_t> indices_set;
    geo_point_t test_p;

    for (size_t pind = 0; pind != 3; ++pind) {
        for (auto it = map_wpid_low_indices.begin(); it != map_wpid_low_indices.end(); ++it) {
            const WirePlaneId& wpid = it->first;
            const auto& low_idxes = it->second;
            const auto& high_idxes = map_wpid_high_indices[wpid];
            if (flags[pind]) {
                // raw data points ... 
                geo_point_t low_p = point3d_raw(low_idxes[pind]);
                geo_point_t high_p = point3d_raw(high_idxes[pind]);
                std::vector<geo_point_t> test_points = {low_p, high_p};
                for (const auto& test_point : test_points) {
                    temp_indices = get_closest_2d_index(test_point, 0.5 * units::cm, wpid.apa(), wpid.face(), pind);
                    std::copy(temp_indices.begin(), temp_indices.end(), inserter(indices_set, indices_set.begin()));
                }
            }
        }
        {
            auto wpid = start_wpid;
            if (flags[pind]) {
                auto test_point = start_p;                
                temp_indices = get_closest_2d_index(test_point, 0.5 * units::cm, wpid.apa(), wpid.face(), pind);
                std::copy(temp_indices.begin(), temp_indices.end(), inserter(indices_set, indices_set.begin()));
            }
        }
        {
            auto wpid = end_wpid;
            if (flags[pind]) {
                auto test_point = end_p;
                temp_indices = get_closest_2d_index(test_point, 0.5 * units::cm, wpid.apa(), wpid.face(), pind);
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
            // double value = pow(winds[0][indices.at(i)] - winds[0][indices.at(j)], 2) +
            //                pow(winds[1][indices.at(i)] - winds[1][indices.at(j)], 2) +
            //                pow(winds[2][indices.at(i)] - winds[2][indices.at(j)], 2);
            double value = pow(point3d_raw(indices.at(i)).x() - point3d_raw(indices.at(j)).x(), 2) +
                           pow(point3d_raw(indices.at(i)).y() - point3d_raw(indices.at(j)).y(), 2) +
                           pow(point3d_raw(indices.at(i)).z() - point3d_raw(indices.at(j)).z(), 2);

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
                geo_point_t new_start_p = point3d_raw(new_start_idx);
                geo_point_t new_end_p = point3d_raw(new_end_idx);

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

const Cluster::sv2d_t& Cluster::sv2d(const int apa, const int face, const size_t plane) const
{
    // if (wpid.layer()!=kAllLayers) {
    //     raise<RuntimeError>("Cluster::sv2d() wpid.layer() {} != kAllLayers");
    // }
    const WirePlaneId wpid(kAllLayers, face, apa);
    const Tree::Scope scope = {"3d", {m_scope2ds_prefix[plane]+"_x", m_scope2ds_prefix[plane]+"_y"}, 0, wpid.name()};
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

const Cluster::kd2d_t& Cluster::kd2d(const int apa, const int face, const size_t plane) const
{
    const auto& sv = sv2d(apa, face, plane);
    return sv.kd();
}

// this point p needs to be raw point, since this is 2D PC ...
std::vector<size_t> Cluster::get_closest_2d_index(const geo_point_t& p, const double search_radius, const int apa, const int face, const int plane) const {

    auto angles = grouping()->wire_angles(apa,face);
    double angle_uvw[3];
    angle_uvw[0] = std::get<0>(angles);
    angle_uvw[1] = std::get<1>(angles);
    angle_uvw[2] = std::get<2>(angles);
    double x = p.x();
    double y = cos(angle_uvw[plane]) * p.z() - sin(angle_uvw[plane]) * p.y();
    std::vector<float_t> query_pt = {x, y};
    const auto& skd = kd2d(apa, face, plane);
    auto ret_matches = skd.radius(search_radius * search_radius, query_pt);

    // local indices ...
    std::vector<size_t> ret_index(ret_matches.size());
    // 2d scoped view ...
    const auto& sv2 = sv2d(apa, face, plane);
    // 3d scoped view
    const auto& sv3 = sv3d();

    const auto error_index = std::numeric_limits<size_t>::max();

    // use 2D local idx --> global-->idx --> 3D local index
    for (size_t i = 0; i != ret_matches.size(); i++)
    {
        size_t global_index = sv2.local_to_global(ret_matches.at(i).first);
        ret_index.at(i) = sv3.global_to_local(global_index);
        if (global_index == error_index || ret_index.at(i) == error_index) {
            throw std::runtime_error("Failed to convert from local to global index");
        }
        
        // std::cout << "Test: " << ret_index.at(i) << " " << global_index << " " << ret_index.at(i) << std::endl;
        // ret_index.at(i) = ret_matches.at(i).first;
    }

    return ret_index;
}

std::vector<const Blob*> Cluster::is_connected(const Cluster& c, const int offset) const
{
    auto& time_blob_map1 = c.time_blob_map();
    auto& time_blob_map2 = time_blob_map();
    std::vector<const Blob*> ret;

    for (auto it = time_blob_map1.begin(); it != time_blob_map1.end(); it++){
        int apa = it->first;
        if (time_blob_map2.find(apa) == time_blob_map2.end()) continue; // if the second one does not contain it ...
        for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
            int face = it1->first; // face
            if (time_blob_map2.at(apa).find(face) == time_blob_map2.at(apa).end()) continue;

            for (const auto& [bad_start, badblobs] : time_blob_map1.at(apa).at(face)) {
                for (const auto* badblob : badblobs) {
                    auto bad_end = badblob->slice_index_max();  // not inclusive
                    for (const auto& [good_start, goodblobs] : time_blob_map2.at(apa).at(face)) {
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
        }
    }
    
    return ret;
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

const Cluster::sv3d_t& Cluster::sv3d() const {
    return sv(); //  m_node->value.scoped_view(m_default_scope);
}
const Cluster::kd3d_t& Cluster::kd3d() const { return sv3d().kd(); }
const Cluster::kd3d_t& Cluster::kd() const { return kd3d(); }
geo_point_t Cluster::point3d(size_t point_index) const { return kd3d().point3d(point_index); }


const Cluster::sv3d_t& Cluster::sv3d_raw() const {
    return sv(m_scope_3d_raw);
    // return m_node->value.scoped_view(m_scope_3d_raw);
}
const Cluster::kd3d_t& Cluster::kd3d_raw() const { return sv3d_raw().kd(); }
geo_point_t Cluster::point3d_raw(size_t point_index) const { return kd3d_raw().point3d(point_index); }

const Cluster::points_type& Cluster::points() const { return kd3d().points(); }
const Cluster::points_type& Cluster::points_raw() const { return kd3d_raw().points(); }


WirePlaneId Cluster::wire_plane_id(size_t point_index) const {  
    auto& wpids = cache().point_wpids;
    if (wpids.empty()) {
        wpids = points_property<int>("wpid");
    }
    return WirePlaneId(wpids[point_index]);
}

int Cluster::wire_index(size_t point_index, int plane) const {
    auto& cache_ref = cache();
    
    switch(plane) {
        case 0: {
            if (cache_ref.point_u_wire_indices.empty()) {
                cache_ref.point_u_wire_indices = points_property<int>("u_wire_index");
            }
            return cache_ref.point_u_wire_indices[point_index];
        }
        case 1: {
            if (cache_ref.point_v_wire_indices.empty()) {
                cache_ref.point_v_wire_indices = points_property<int>("v_wire_index");
            }
            return cache_ref.point_v_wire_indices[point_index];
        }
        case 2: {
            if (cache_ref.point_w_wire_indices.empty()) {
                cache_ref.point_w_wire_indices = points_property<int>("w_wire_index");
            }
            return cache_ref.point_w_wire_indices[point_index];
        }
        default: 
            raise<ValueError>("Invalid plane index: %d (must be 0, 1, or 2)", plane);
    }
    std::terminate(); // this is here mostly to quell compiler warnings about not returning a value.
}

double Cluster::charge_value(size_t point_index, int plane) const {
    auto& cache_ref = cache();
    
    switch(plane) {
        case 0: {
            if (cache_ref.point_u_charges.empty()) {
                cache_ref.point_u_charges = points_property<double>("u_charge_val");
            }
            return cache_ref.point_u_charges[point_index];
        }
        case 1: {
            if (cache_ref.point_v_charges.empty()) {
                cache_ref.point_v_charges = points_property<double>("v_charge_val");
            }
            return cache_ref.point_v_charges[point_index];
        }
        case 2: {
            if (cache_ref.point_w_charges.empty()) {
                cache_ref.point_w_charges = points_property<double>("w_charge_val");
            }
            return cache_ref.point_w_charges[point_index];
        }
        default:
            raise<ValueError>("Invalid plane index: %d (must be 0, 1, or 2)", plane);
    }
    std::terminate(); // this is here mostly to quell compiler warnings about not returning a value.
}

double Cluster::charge_uncertainty(size_t point_index, int plane) const {
    auto& cache_ref = cache();
    
    switch(plane) {
        case 0: {
            if (cache_ref.point_u_charge_uncs.empty()) {
                cache_ref.point_u_charge_uncs = points_property<double>("u_charge_unc");
            }
            return cache_ref.point_u_charge_uncs[point_index];
        }
        case 1: {
            if (cache_ref.point_v_charge_uncs.empty()) {
                cache_ref.point_v_charge_uncs = points_property<double>("v_charge_unc");
            }
            return cache_ref.point_v_charge_uncs[point_index];
        }
        case 2: {
            if (cache_ref.point_w_charge_uncs.empty()) {
                cache_ref.point_w_charge_uncs = points_property<double>("w_charge_unc");
            }
            return cache_ref.point_w_charge_uncs[point_index];
        }
        default:
            raise<ValueError>("Invalid plane index: %d (must be 0, 1, or 2)", plane);
    }
    std::terminate(); // this is here mostly to quell compiler warnings about not returning a value.
}

bool Cluster::is_wire_dead(size_t point_index, int plane, double dead_threshold) const {
    return charge_uncertainty(point_index, plane) > dead_threshold;
}

std::pair<bool, double> Cluster::calc_charge_wcp(
    size_t point_index,
    double charge_cut,
    bool disable_dead_mix_cell) const {
    
    const double dead_threshold = 1e10; // Same as PointTreeBuilding
    
    double charge = 0;
    int ncharge = 0;
    
    // Get exact charges for u,v,w wires using cached data
    double charge_u = charge_value(point_index, 0);
    double charge_v = charge_value(point_index, 1);
    double charge_w = charge_value(point_index, 2);
    
    // Check for dead wires
    bool is_dead_u = is_wire_dead(point_index, 0, dead_threshold);
    bool is_dead_v = is_wire_dead(point_index, 1, dead_threshold);
    bool is_dead_w = is_wire_dead(point_index, 2, dead_threshold);
    
    bool flag_charge_u = false;
    bool flag_charge_v = false;
    bool flag_charge_w = false;

    // Initial flag setting based on charge threshold
    if (charge_u > charge_cut) flag_charge_u = true;
    if (charge_v > charge_cut) flag_charge_v = true;
    if (charge_w > charge_cut) flag_charge_w = true;
    
    if (disable_dead_mix_cell) {
        // Add all charges first
        charge += charge_u * charge_u; ncharge++;
        charge += charge_v * charge_v; ncharge++;
        charge += charge_w * charge_w; ncharge++;
        
        // Deal with bad planes - subtract dead wire contributions
        if (is_dead_u) {
            flag_charge_u = true;
            charge -= charge_u * charge_u; ncharge--;
        }
        if (is_dead_v) {
            flag_charge_v = true;
            charge -= charge_v * charge_v; ncharge--;
        }
        if (is_dead_w) {
            flag_charge_w = true;
            charge -= charge_w * charge_w; ncharge--;
        }
    } else {
        // Only use non-zero charges
        if (charge_u == 0) flag_charge_u = true;
        if (charge_v == 0) flag_charge_v = true;
        if (charge_w == 0) flag_charge_w = true;

        if (charge_u != 0) {
            charge += charge_u * charge_u; ncharge++;
        }
        if (charge_v != 0) {
            charge += charge_v * charge_v; ncharge++;
        }
        if (charge_w != 0) {
            charge += charge_w * charge_w; ncharge++;
        }
    }
    
    // Require more than one plane to be good 
    if (ncharge > 1) {
        charge = sqrt(charge / ncharge);
    } else {
        charge = 0;
    }
    
    return std::make_pair(flag_charge_u && flag_charge_v && flag_charge_w, charge);
}

// Convenience overload for 3D points
std::pair<bool, double> Cluster::calc_charge_wcp(
    const geo_point_t& point,
    double charge_cut,
    bool disable_dead_mix_cell) const {
    
    size_t point_index = get_closest_point_index(point);
    return calc_charge_wcp(point_index, charge_cut, disable_dead_mix_cell);
}



int Cluster::npoints() const
{
    auto& n = cache().npoints;
    if (!n) {
        const auto& sv = sv3d();
        n = sv.npoints();
    }
    return n;
}



// size_t Cluster::nbpoints() const
// {
//     size_t ret = 0;
//     for (const auto* blob : children()) {
//         ret += blob->nbpoints();
//     }
//     return ret;
// }

const Cluster::wire_indices_t& Cluster::wire_indices() const
{
    const auto& sv = m_node->value.scoped_view<int_t>(m_scope_wire_index);
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

// std::vector<geo_point_t> Cluster::kd_points_raw(const Cluster::kd_results_t& res)
// {
//     return const_cast<const Cluster*>(this)->kd_points_raw(res);
// }
// std::vector<geo_point_t> Cluster::kd_points_raw(const Cluster::kd_results_t& res) const
// {
//     std::vector<geo_point_t> ret;
//     const auto& points = this->points_raw();
//     for (const auto& [point_index, _] : res) {
//         ret.emplace_back(points[0][point_index], points[1][point_index], points[2][point_index]);
//     }
//     return ret;
// }

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
    auto& length = cache().length;
    if (length != 0) {
        return length;
    }

    const auto& grouping = this->grouping();

    auto map_wpid_uvwt = this->get_uvwt_range();
    for (const auto& [wpid, uvwt] : map_wpid_uvwt) {

        const double tick = grouping->get_tick().at(wpid.apa()).at(wpid.face());
        const double drift_speed = grouping->get_drift_speed().at(wpid.apa()).at(wpid.face());

        // std::cout << "Test: " << wpid.apa() << " " << wpid.face() << " " << tp.tick_drift << " " << tick * drift_speed << std::endl;

        const auto [u, v, w, t] = uvwt;
        auto face = grouping->get_anode(wpid.apa())->face(wpid.face());
        const double pu = u * face->plane(0)->pimpos()->pitch() ;
        const double pv = v * face->plane(1)->pimpos()->pitch();
        const double pw = w * face->plane(2)->pimpos()->pitch();
        const double pt = t * tick * drift_speed;
        length += std::sqrt(2. / 3. * (pu * pu + pv * pv + pw * pw) + pt * pt);
    }

    return length;
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
    geo_point_t main_axis = get_pca().axis.at(0); 
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
        const auto* svptr = m_node->value.get_scoped(m_default_scope);
        if (!svptr) {
            if (log) log->debug("cluster sanity: note, not yet a scoped view {}", m_default_scope);
        }
    }
    if (!nchildren()) {
        if (log) log->debug("cluster sanity: no children blobs");
        return false;
    }

    const auto& sv = m_node->value.scoped_view(m_default_scope);
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
            spoints = sblob->points(get_default_scope().pcname, get_default_scope().coords);
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
    // sort_blobs(blobs);
    for (const Blob* blob : blobs) {
        boost::hash_combine(h, blob->hash());
    }
    return h;
}

std::vector<int> Cluster::get_blob_indices(const Blob* blob) const
{
    auto& mmi = cache().map_mcell_indices;
    if (mmi.empty()) {
        const auto& skd = kd3d();
        for (size_t ind = 0; ind < skd.npoints(); ++ind) {
            const auto* bwp = blob_with_point(ind);
            mmi[bwp].push_back(ind);
        }
    }
    return mmi[blob];
}

std::vector<geo_point_t> Cluster::indices_to_points(const std::vector<size_t>& path_indices) const 
{
    std::vector<geo_point_t> points;
    points.reserve(path_indices.size());
    for (size_t idx : path_indices) {
        points.push_back(point3d(idx));
    }
    return points;
}

// void Cluster::organize_points_path_vec(std::vector<geo_point_t>& path_points, double low_dis_limit) const
// {
//     std::vector<geo_point_t> temp_points = path_points;
//     path_points.clear();

//     // First pass: filter based on distance
//     for (size_t i = 0; i != temp_points.size(); i++) {
//         if (path_points.empty()) {
//             path_points.push_back(temp_points[i]);
//         }
//         else if (i + 1 == temp_points.size()) {
//             double dis = (temp_points[i] - path_points.back()).magnitude();
//             if (dis > low_dis_limit * 0.75) {
//                 path_points.push_back(temp_points[i]);
//             }
//         }
//         else {
//             double dis = (temp_points[i] - path_points.back()).magnitude();
//             double dis1 = (temp_points[i + 1] - path_points.back()).magnitude();

//             if (dis > low_dis_limit || (dis1 > low_dis_limit * 1.7 && dis > low_dis_limit * 0.75)) {
//                 path_points.push_back(temp_points[i]);
//             }
//         }
//     }

//     // Second pass: filter based on angle
//     temp_points = path_points;
//     std::vector<double> angles;
//     for (size_t i = 0; i != temp_points.size(); i++) {
//         if (i == 0 || i + 1 == temp_points.size()) {
//             angles.push_back(M_PI);
//         }
//         else {
//             geo_vector_t v1 = temp_points[i] - temp_points[i - 1];
//             geo_vector_t v2 = temp_points[i] - temp_points[i + 1];
//             angles.push_back(v1.angle(v2));
//         }
//     }

//     path_points.clear();
//     for (size_t i = 0; i != temp_points.size(); i++) {
//         if (angles[i] * 180.0 / M_PI >= 75) {
//             path_points.push_back(temp_points[i]);
//         }
//     }
// }

// // this is different from WCP implementation, the path_points is the input ...
// void Cluster::organize_path_points(std::vector<geo_point_t>& path_points, double low_dis_limit) const
// {
//     //    std::vector<geo_point_t> temp_points = path_points;
//     path_points.clear();
//     auto indices = get_path_wcps();
//     auto temp_points = indices_to_points(indices);

//     for (size_t i = 0; i != temp_points.size(); i++) {
//         if (path_points.empty()) {
//             path_points.push_back(temp_points[i]);
//         }
//         else if (i + 1 == temp_points.size()) {
//             double dis = (temp_points[i] - path_points.back()).magnitude();
//             if (dis > low_dis_limit * 0.5) {
//                 path_points.push_back(temp_points[i]);
//             }
//         }
//         else {
//             double dis = (temp_points[i] - path_points.back()).magnitude();
//             double dis1 = (temp_points[i + 1] - path_points.back()).magnitude();

//             if (dis > low_dis_limit || (dis1 > low_dis_limit * 1.7 && dis > low_dis_limit * 0.5)) {
//                 path_points.push_back(temp_points[i]);
//             }
//         }
//     }
// }


std::vector<geo_point_t> Cluster::get_hull() const 
{
    auto& hull_points = cache().hull_points;

    if (hull_points.size()) {
        return hull_points;
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

    for (auto i : indices) {
        hull_points.push_back({points[0][i], points[1][i], points[2][i]});
    }
    
    return hull_points;
}

Cluster::PCA& Cluster::get_pca() const
{
    auto& pcaptr = cache().pca;
    if (pcaptr) {
        return *pcaptr;
    }

    const auto& pcname = this->get_default_scope().pcname;
    const auto& coords = this->get_default_scope().coords;

    pcaptr = std::make_unique<PCA>();
    pcaptr->axis.resize(3);
    pcaptr->values.resize(3,0);

    int nsum = 0;
    for (const Blob* blob : children()) {
        for (const geo_point_t& p : blob->points(pcname, coords)) {
            pcaptr->center += p;
            nsum++;
        }
    }

    // Not enough points to perform PCA.
    if (nsum < 3) {
        return *pcaptr;
    }

    pcaptr->center /= nsum;

    Eigen::MatrixXd cov_matrix(3, 3);

    for (int i = 0; i != 3; i++) {
        for (int j = i; j != 3; j++) {
            cov_matrix(i, j) = 0;
            for (const Blob* blob : children()) {
                for (const geo_point_t& p : blob->points(pcname, coords)) {
                    cov_matrix(i, j) += (p[i] - pcaptr->center[i]) * (p[j] - pcaptr->center[j]);
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
        pcaptr->values[2-i] = eigen_values(i);
        double norm = sqrt(eigen_vectors(0, i) * eigen_vectors(0, i) + eigen_vectors(1, i) * eigen_vectors(1, i) +
                             eigen_vectors(2, i) * eigen_vectors(2, i));
        pcaptr->axis[2-i].set(eigen_vectors(0, i) / norm, eigen_vectors(1, i) / norm, eigen_vectors(2, i) / norm);
    }

    return *pcaptr;
}


// std::unordered_map<int, Cluster*> 
std::vector<int> Cluster::examine_x_boundary(const double low_limit, const double high_limit)
// designed to run for single face ... limits are for per face only ...
{
    double num_points[3] = {0, 0, 0};
    double x_max = -1e9;
    double x_min = 1e9;
    auto& mcells = this->children();
    const auto& pcname = this->get_default_scope().pcname;
    const auto& coords = this->get_default_scope().coords;

    for (Blob* mcell : mcells) {
        /// TODO: no caching, could be slow
        std::vector<geo_point_t> pts = mcell->points(pcname, coords);
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
                if (mcell->points(pcname, coords)[0].x() < low_limit) {
                    if (groupids.find(1) != groupids.end()) {
                        // cluster_1->AddCell(mcell, mcell->GetTimeSlice());
                        b2groupid[idx] = 1;
                    }
                    else {
                        // cluster_2->AddCell(mcell, mcell->GetTimeSlice());
                        b2groupid[idx] = 2;
                    }
                }
                else if (mcell->points(pcname, coords)[0].x() > high_limit) {
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

bool Cluster::judge_vertex(geo_point_t& p_test, IDetectorVolumes::pointer dv, const double asy_cut, const double occupied_cut)
{
    p_test = this->calc_ave_pos(p_test, 3 * units::cm);

    geo_point_t dir = this->vhough_transform(p_test, 15 * units::cm);

    // judge if this is end points
    std::pair<int, int> num_pts = this->ndipole(p_test, dir, 25 * units::cm);

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
            auto idx = this->get_closest_point_index(p_test); 
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
        const int N = this->npoints();
        // WCP::WCPointCloud<double>& cloud = point_cloud->get_cloud();
        for (int i = 0; i != N; i++) {
            // geo_point_t dir1(cloud.pts[i].x() - p_test.x(), cloud.pts[i].y() - p_test.y(), cloud.pts[i].z() - p_test.z());
            geo_point_t dir1 = this->point3d(i) - p_test;

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
    {
        auto ac = a->get_pca().center;
        auto bc = b->get_pca().center;
        if (ac[0] < bc[0]) return true;
        if (bc[0] < bc[0]) return false;
        if (ac[1] < bc[1]) return true;
        if (bc[1] < bc[1]) return false;
        if (ac[2] < bc[2]) return true;
        if (bc[2] < bc[2]) return false;
    }

    // After exhausting all "content" comparison, we are left with the question,
    // are these two clusters really different or not.  We have two choices.  We
    // may compare on pointer value which will surely "break the tie" but will
    // introduce randomness.  We may return "false" which says "these are equal"
    // in which case any unordered set/map will not hold both.  Randomness is
    // the better choice as we would have a better chance to detect that in some
    // future bug.
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

    const auto* p = this->node()->parent;
    if (!p)  return flash;
    const auto* g = p->value.facade<Grouping>();
    if (!g)  return flash;

    const int flash_index = this->get_scalar("flash", -1);

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


const Facade::Cluster::graph_type& Facade::Cluster::find_graph(const std::string& flavor) const
{
    return const_cast<const graph_type&>(const_cast<Cluster*>(this)->find_graph(flavor));
}
Facade::Cluster::graph_type& Facade::Cluster::find_graph(const std::string& flavor)
{
    if (this->has_graph(flavor)) {
        return get_graph(flavor);
    }
    if (flavor == "basic") {
        return this->give_graph(flavor, make_graph_basic(*this));
    }

    // We did our best....
    raise<KeyError>("unknown graph flavor " + flavor);
    std::terminate(); // this is here mostly to quell compiler warnings about not returning a value.
}
const Facade::Cluster::graph_type& Facade::Cluster::find_graph(
    const std::string& flavor,
    IDetectorVolumes::pointer dv, 
    IPCTransformSet::pointer pcts) const
{
    return const_cast<const graph_type&>(const_cast<Cluster*>(this)->find_graph(flavor, dv, pcts));
}

Facade::Cluster::graph_type& Facade::Cluster::find_graph(
    const std::string& flavor,
    IDetectorVolumes::pointer dv, 
    IPCTransformSet::pointer pcts)
{
    if (this->has_graph(flavor)) {
        return get_graph(flavor);
    }

    // Factory of known graph flavors relying on detector info:

    if (flavor == "ctpc") {
        return this->give_graph(flavor, make_graph_basic(*this));
    }

    if (flavor == "relaxed") {
        return this->give_graph(flavor, make_graph_basic(*this));
    }

    // Do a hail mary, maybe user made a mistake by passing dv/pcts and really
    // wants a flavor that we can make implicitly.
    return find_graph(flavor);
}


const GraphAlgorithms& Facade::Cluster::graph_algorithms(const std::string& flavor) const
{
    auto it = m_galgs.find(flavor);
    if (it != m_galgs.end()) {
        return it->second;      // we have it already
    }

    if (this->has_graph(flavor)) {    // if graph exists, make the GA
        auto got = m_galgs.emplace(flavor, GraphAlgorithms(get_graph(flavor)));
        return got.first->second;
    }
        
    // We failed to find an existing graph of the given flavor, but we there are
    // some flavors we know how to construct on the fly:

    if (flavor == "basic") {
        // we are caching, so const cast is "okay".
        auto& gr = const_cast<Cluster*>(this)->give_graph(flavor, make_graph_basic(*this));
        auto got = m_galgs.emplace(flavor, GraphAlgorithms(gr));
        return got.first->second;
    }

    // We did our best....
    raise<KeyError>("unknown graph flavor " + flavor);
    std::terminate(); // this is here mostly to quell compiler warnings about not returning a value.
}

const GraphAlgorithms& Facade::Cluster::graph_algorithms(const std::string& flavor,
                                                                   IDetectorVolumes::pointer dv, 
                                                                   IPCTransformSet::pointer pcts) const
{
    auto it = m_galgs.find(flavor);
    if (it != m_galgs.end()) {
        return it->second;
    }

    // Factory of known graph flavors relying on detector info:

    if (flavor == "ctpc") {
        auto& gr = const_cast<Cluster*>(this)->give_graph(flavor, make_graph_basic(*this));
        auto got = m_galgs.emplace(flavor, GraphAlgorithms(gr));
        return got.first->second;
    }

    if (flavor == "relaxed") {
        auto& gr = const_cast<Cluster*>(this)->give_graph(flavor, make_graph_basic(*this));
        auto got = m_galgs.emplace(flavor, GraphAlgorithms(gr));
        return got.first->second;
    }

    // Do a hail mary, maybe user made a mistake by passing dv/pcts and really
    // wants a flavor that we can make implicitly.
    return graph_algorithms(flavor);
}


// These methods implement the cache management functionality

void Facade::Cluster::clear_graph_algorithms_cache(const std::string& graph_name)
{
    auto it = m_galgs.find(graph_name);
    if (it != m_galgs.end()) {
        it->second.clear_cache();
        auto log = Log::logger("clus");
        log->debug("Cleared cache for GraphAlgorithms '{}'", graph_name);
    }
}

void Facade::Cluster::remove_graph_algorithms(const std::string& graph_name)
{
    auto it = m_galgs.find(graph_name);
    if (it != m_galgs.end()) {
        m_galgs.erase(it);
        auto log = Log::logger("clus");
        log->debug("Removed GraphAlgorithms '{}'", graph_name);
    }
}

void Facade::Cluster::clear_all_graph_algorithms_caches()
{
    for (auto& [name, ga] : m_galgs) {
        ga.clear_cache();
    }
    auto log = Log::logger("clus");
    log->debug("Cleared all GraphAlgorithms caches");
}

std::vector<std::string> Facade::Cluster::get_cached_graph_algorithms() const
{
    std::vector<std::string> names;
    names.reserve(m_galgs.size());
    for (const auto& [name, ga] : m_galgs) {
        names.push_back(name);
    }
    return names;
}


// ne' examine_graph
std::vector<int> Cluster::connected_blobs(IDetectorVolumes::pointer dv, IPCTransformSet::pointer pcts) const 
{
    const auto& ga = graph_algorithms("relaxed", dv, pcts);
    const auto& component = ga.connected_components();

    // Create mapping from blob indices to component groups
    std::vector<int> b2groupid(nchildren(), -1);
    
    // For each point in the graph
    for (size_t i = 0; i < component.size(); ++i) {
        // Get the blob index for this point
        const int bind = kd3d().major_index(i);
        // Map the blob to its component
        b2groupid.at(bind) = (int)component[i];
    }

    return b2groupid;
}


std::vector<std::vector<geo_point_t>> Cluster::get_extreme_wcps(const Cluster* reference_cluster) const
{
    std::vector<std::vector<geo_point_t>> out_vec_wcps;
    
    if (npoints() == 0) {
        return out_vec_wcps;
    }
    
    // Create list of valid point indices based on spatial filtering
    // This directly corresponds to prototype's all_indices creation
    std::vector<size_t> valid_indices;
    
    if (reference_cluster == nullptr) {
        // No filtering - use all points (equivalent to old_time_mcells_map==0)
        for (int i = 0; i < npoints(); ++i) {
            valid_indices.push_back(i);
        }
    } else {
        // Get reference cluster's time_blob_map (equivalent to old_time_mcells_map)
        const auto& ref_time_blob_map = reference_cluster->time_blob_map();
        
        // Filter points based on spatial relationship with reference cluster
        // This implements the exact same logic as prototype's old_time_mcells_map filtering
        for (int i = 0; i < npoints(); ++i) {
            if (is_point_spatially_related_to_time_blobs(i, ref_time_blob_map)) {
                valid_indices.push_back(i);
            }
        }
    }
    
    if (valid_indices.empty()) {
        return out_vec_wcps;
    }
    
    // Get main axis and ensure consistent direction (y>0)
    // Equivalent to prototype's Calc_PCA() and get_PCA_axis(0)
    geo_point_t main_axis = get_pca().axis.at(0);
    if (main_axis.y() < 0) {
        main_axis = main_axis * -1;
    }
    
    // Find 8 extreme points: 2 along main axis + 6 along coordinate axes
    // Equivalent to prototype's wcps[8] array
    geo_point_t extreme_points[8];
    
    // Initialize with first valid point
    // Equivalent to prototype: wcps[i] = cloud.pts[all_indices.at(0)]
    for (int i = 0; i < 8; i++) {
        extreme_points[i] = point3d(valid_indices[0]);
    }
    
    // Initialize projection values for main axis extremes
    double high_value = extreme_points[0].dot(main_axis);
    double low_value = high_value;
    
    // Scan through all valid points to find extremes
    // Equivalent to prototype's scanning loop through all_indices
    for (size_t idx : valid_indices) {
        geo_point_t current_point = point3d(idx);
        
        // Main axis extremes (along PCA axis)
        double main_projection = current_point.dot(main_axis);
        if (main_projection > high_value) {
            extreme_points[0] = current_point;  // high along main axis
            high_value = main_projection;
        }
        if (main_projection < low_value) {
            extreme_points[1] = current_point;  // low along main axis
            low_value = main_projection;
        }
        
        // Coordinate axis extremes (same as prototype)
        // Y-axis extremes (top/bottom)
        if (current_point.y() > extreme_points[2].y()) {
            extreme_points[2] = current_point;  // highest Y
        }
        if (current_point.y() < extreme_points[3].y()) {
            extreme_points[3] = current_point;  // lowest Y
        }
        
        // Z-axis extremes (front/back)
        if (current_point.z() > extreme_points[4].z()) {
            extreme_points[4] = current_point;  // furthest Z
        }
        if (current_point.z() < extreme_points[5].z()) {
            extreme_points[5] = current_point;  // nearest Z
        }
        
        // X-axis extremes (earliest/latest)
        if (current_point.x() > extreme_points[6].x()) {
            extreme_points[6] = current_point;  // latest X
        }
        if (current_point.x() < extreme_points[7].x()) {
            extreme_points[7] = current_point;  // earliest X
        }
    }
    
    // Group the extreme points into result vectors
    // Following the prototype's grouping strategy exactly
    
    // First extreme along the main axis
    std::vector<geo_point_t> main_axis_high;
    main_axis_high.push_back(extreme_points[0]);
    out_vec_wcps.push_back(main_axis_high);
    
    // Second extreme along the main axis  
    std::vector<geo_point_t> main_axis_low;
    main_axis_low.push_back(extreme_points[1]);
    out_vec_wcps.push_back(main_axis_low);
    
    // Add other extremes if they are significantly different from main axis extremes
    // This prevents duplicate points in the output (same as prototype logic)
    const double min_separation = 5.0 * units::cm;  // Minimum distance to be considered distinct
    
    for (int i = 2; i < 8; i++) {
        bool is_distinct = true;
        
        // Check if this extreme is too close to already added points
        for (const auto& added_group : out_vec_wcps) {
            for (const auto& added_point : added_group) {
                double distance = (extreme_points[i] - added_point).magnitude();
                if (distance < min_separation) {
                    is_distinct = false;
                    break;
                }
            }
            if (!is_distinct) break;
        }
        
        // If distinct enough, add as a new extreme group
        if (is_distinct) {
            std::vector<geo_point_t> coord_extreme;
            coord_extreme.push_back(extreme_points[i]);
            out_vec_wcps.push_back(coord_extreme);
        }
    }
    
    return out_vec_wcps;
}

bool Cluster::is_point_spatially_related_to_time_blobs(
    size_t point_index, 
    const time_blob_map_t& ref_time_blob_map) const
{
    // Get current point's time slice information
    // Equivalent to: int time_slice = cloud.pts[i].mcell->GetTimeSlice();
    const Blob* current_blob = blob_with_point(point_index);
    int current_time_slice = current_blob->slice_index_min();
    
    // Check current time slice and ±1 time slices (like Steiner version)
    for (int time_offset = -1; time_offset <= 1; ++time_offset) {
        int check_time_slice = current_time_slice + time_offset;
        
        // This is the exact prototype logic:
        // if (old_time_mcells_map->find(time_slice)!=old_time_mcells_map->end())
        auto time_it = ref_time_blob_map.find(check_time_slice);
        if (time_it != ref_time_blob_map.end()) {
            
            // Iterate through apa/face maps in this time slice
            // time_blob_map_t is std::map<int, std::map<int, std::map<int, BlobSet>>>
            // Structure: apa -> face -> time -> blobset
            for (const auto& face_pair : time_it->second) {
                for (const auto& time_pair : face_pair.second) {
                    // Now iterate through blobs in the BlobSet
                    for (const Blob* ref_blob : time_pair.second) {
                        
                        // Method 1: Fast blob overlap check (equivalent to mcell->Overlap_fast())
                        if (current_blob->overlap_fast(*ref_blob, 1)) {
                            return true;  // Equivalent to flag_add = true; break;
                        }
                        
                        // Method 2: Detailed wire range checking (prototype's exact logic)
                        if (check_wire_ranges_match(point_index, ref_blob)) {
                            return true;  // Equivalent to flag_add = true; break;
                        }
                    }
                }
            }
        }
    }
    
    return false;  // Equivalent to flag_add remains false
}

bool Cluster::check_wire_ranges_match(size_t point_index, const Blob* ref_blob) const
{
    try {
        // Get current point's wire indices (equivalent to cloud.pts[i].index_u, index_v, index_w)
        int current_wire_u = wire_index(point_index, 0);  // U plane
        int current_wire_v = wire_index(point_index, 1);  // V plane  
        int current_wire_w = wire_index(point_index, 2);  // W plane
        
        // Get reference blob's wire ranges
        // Equivalent to: 
        // int u1_low_index = mcell->get_uwires().front()->index();
        // int u1_high_index = mcell->get_uwires().back()->index();
        // Using individual wire range methods instead of get_wire_ranges()
        int u_min = ref_blob->u_wire_index_min();
        int u_max = ref_blob->u_wire_index_max();
        int v_min = ref_blob->v_wire_index_min();
        int v_max = ref_blob->v_wire_index_max();
        int w_min = ref_blob->w_wire_index_min();
        int w_max = ref_blob->w_wire_index_max();
        
        // Extract U, V, W wire ranges from reference blob
        // Add ±1 tolerance (Steiner version) to the wire ranges
        u_min = u_min - 1;
        u_max = u_max + 1;
        v_min = v_min - 1;
        v_max = v_max + 1;
        w_min = w_min - 1;
        w_max = w_max + 1;
        
        // Check if current point's wire indices fall within ALL THREE ranges
        // This is the exact prototype condition:
        // if (cloud.pts[i].index_u <= u1_high_index && cloud.pts[i].index_u >= u1_low_index &&
        //     cloud.pts[i].index_v <= v1_high_index && cloud.pts[i].index_v >= v1_low_index &&
        //     cloud.pts[i].index_w <= w1_high_index && cloud.pts[i].index_w >= w1_low_index)
        if (current_wire_u >= u_min && current_wire_u <= u_max &&
            current_wire_v >= v_min && current_wire_v <= v_max &&
            current_wire_w >= w_min && current_wire_w <= w_max) {
            return true;  // Equivalent to flag_add = true; break;
        }
        
    } catch (...) {
        // If wire information is not available, continue
    }
    
    return false;
}





// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
