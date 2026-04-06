#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/Facade.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/FiducialUtils.h"

#include <cmath>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

namespace WireCell::Clus::Facade {

// Function to compute wire plane parameters
void compute_wireplane_params(
    const std::set<WirePlaneId>& wpids,
    const IDetectorVolumes::pointer dv,
    std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>>& wpid_params,
    std::map<WirePlaneId, std::pair<geo_point_t, double>>& wpid_U_dir,
    std::map<WirePlaneId, std::pair<geo_point_t, double>>& wpid_V_dir,
    std::map<WirePlaneId, std::pair<geo_point_t, double>>& wpid_W_dir,
    std::set<int>& apas)
{
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
        wpid_U_dir[wpid] = std::make_pair(geo_point_t(0, std::cos(angle_u), std::sin(angle_u)), angle_u);
        wpid_V_dir[wpid] = std::make_pair(geo_point_t(0, std::cos(angle_v), std::sin(angle_v)), angle_v);
        wpid_W_dir[wpid] = std::make_pair(geo_point_t(0, std::cos(angle_w), std::sin(angle_w)), angle_w);
    }
}


// ---------------------------------------------------------------------------
// cluster_fc_check
//
// Two-round cluster boundary check shared by TaggerCheckNeutrino and
// TaggerCheckSTM.  Lifted from TaggerCheckSTM::check_stm_conditions.
//
// Round 1: steiner boundary with flag_cosmic=true.
// Round 2: steiner boundary with flag_cosmic=false (only run when round 1
//          finds no clear exit candidate).
//
// For each candidate endpoint the check is:
//   1. Direct fiducial-volume test.
//   2. Wire-angle test → check_signal_processing.
//   3. Direction-vs-PCA test → check_dead_volume.
//
// Returns FCCheckResult with is_fc=false (conservative) when the cluster
// has no steiner_pc or when FiducialUtils is unavailable.
// ---------------------------------------------------------------------------
FCCheckResult cluster_fc_check(Cluster& cluster, IDetectorVolumes::pointer dv)
{
    FCCheckResult result;

    // Require steiner_pc
    if (!cluster.has_pc("steiner_pc") || cluster.get_pc("steiner_pc").size() == 0) {
        return result;
    }

    // Require FiducialUtils
    auto fiducial_utils = cluster.grouping()->get_fiducialutils();
    if (!fiducial_utils) {
        return result;
    }

    // Wire-plane geometry maps needed for signal-processing / dead-volume checks
    const auto& wpids = cluster.grouping()->wpids();
    std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> wpid_params;
    std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_U_dir, wpid_V_dir, wpid_W_dir;
    std::set<int> apas;
    compute_wireplane_params(wpids, dv, wpid_params, wpid_U_dir, wpid_V_dir, wpid_W_dir, apas);

    // PCA main direction — used to gate the dead-volume check
    const auto& pca = cluster.get_pca();
    geo_vector_t main_dir = pca.axis.at(0);

    // Per-point direction cache updated as we iterate
    geo_point_t drift_dir{}, U_dir{}, V_dir{}, W_dir{};

    // ------------------------------------------------------------------
    // Helper: one pass of the boundary check over a set of point groups.
    // Appends exit candidates to result.exit_wcps.
    // ------------------------------------------------------------------
    auto do_boundary_check = [&](std::vector<std::vector<geo_point_t>>& pts_groups) {
        for (size_t i = 0; i != pts_groups.size(); i++) {
            bool flag_save = false;

            // Direct fiducial check for every point in the group
            for (size_t j = 0; j != pts_groups.at(i).size(); j++) {
                geo_point_t p1 = pts_groups.at(i).at(j);
                if (!fiducial_utils->inside_fiducial_volume(p1)) {
                    result.exit_wcps.push_back(pts_groups.at(i).at(0));
                    flag_save = true;
                    break;
                }
            }

            if (!flag_save) {
                // Direction-based checks using the outward-pointing Hough direction
                geo_point_t p1 = pts_groups.at(i).at(0);
                geo_vector_t dir_vec = cluster.vhough_transform(p1, 30*units::cm);
                dir_vec = dir_vec * (-1.0);  // reverse to point outward

                geo_point_t dir(dir_vec.x(), dir_vec.y(), dir_vec.z());
                geo_point_t dir_1(0, dir.y(), dir.z());

                // Look up wire-plane geometry for this point's APA/face
                auto wpid_p1 = dv->contained_by(p1);
                auto it_params = wpid_params.find(wpid_p1);
                if (it_params != wpid_params.end())
                    drift_dir = std::get<0>(it_params->second);
                auto it_U = wpid_U_dir.find(wpid_p1);
                if (it_U != wpid_U_dir.end()) U_dir = it_U->second.first;
                auto it_V = wpid_V_dir.find(wpid_p1);
                if (it_V != wpid_V_dir.end()) V_dir = it_V->second.first;
                auto it_W = wpid_W_dir.find(wpid_p1);
                if (it_W != wpid_W_dir.end()) W_dir = it_W->second.first;

                // Projected angles w.r.t. each wire-plane direction
                double angle1 = acos(dir_1.dot(U_dir) / (dir_1.magnitude() * U_dir.magnitude()));
                geo_point_t tempV1(fabs(dir.x()),
                                   sqrt(dir.y()*dir.y() + dir.z()*dir.z()) * sin(angle1), 0);
                double angle1_1 = acos(tempV1.dot(drift_dir) /
                                       (tempV1.magnitude() * drift_dir.magnitude())) / 3.1415926 * 180.;

                double angle2 = acos(dir_1.dot(V_dir) / (dir_1.magnitude() * V_dir.magnitude()));
                geo_point_t tempV2(fabs(dir.x()),
                                   sqrt(dir.y()*dir.y() + dir.z()*dir.z()) * sin(angle2), 0);
                double angle2_1 = acos(tempV2.dot(drift_dir) /
                                       (tempV2.magnitude() * drift_dir.magnitude())) / 3.1415926 * 180.;

                double angle3 = acos(dir_1.dot(W_dir) / (dir_1.magnitude() * W_dir.magnitude()));
                geo_point_t tempV3(fabs(dir.x()),
                                   sqrt(dir.y()*dir.y() + dir.z()*dir.z()) * sin(angle3), 0);
                double angle3_1 = acos(tempV3.dot(drift_dir) /
                                       (tempV3.magnitude() * drift_dir.magnitude())) / 3.1415926 * 180.;

                // Near-collinear with a wire plane → check signal-processing gaps
                if (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5) {
                    if (!fiducial_utils->check_signal_processing(cluster, p1, dir_vec, 1*units::cm)) {
                        flag_save = true;
                        result.exit_wcps.push_back(pts_groups.at(i).at(0));
                    }
                }

                // Direction largely transverse to the cluster axis → check dead volume
                if (!flag_save) {
                    double main_angle = acos(dir_vec.dot(main_dir) /
                                            (dir_vec.magnitude() * main_dir.magnitude()));
                    double angle_deg = fabs((3.1415926/2. - main_angle) / 3.1415926 * 180.);
                    if (angle_deg > 60) {
                        if (!fiducial_utils->check_dead_volume(cluster, p1, dir_vec, 1*units::cm)) {
                            result.exit_wcps.push_back(pts_groups.at(i).at(0));
                        }
                    }
                }
            }
        }
    };

    // ------------------------------------------------------------------
    // Helper: map exit_wcps to steiner boundary indices (0 or 1) and
    // handle the two-endpoint-exit protection.
    // ------------------------------------------------------------------
    auto update_boundary_set = [&](const geo_point_t& bp1, const geo_point_t& bp2) {
        for (const auto& ewcp : result.exit_wcps) {
            double dis1 = (ewcp - bp1).magnitude();
            double dis2 = (ewcp - bp2).magnitude();
            if (dis1 < dis2) {
                if (dis1 < 1.0*units::cm) result.exit_boundary_set.insert(0);
            } else {
                if (dis2 < 1.0*units::cm) result.exit_boundary_set.insert(1);
            }
        }
        // If both endpoints appear to exit, use direct fiducial check to arbitrate
        if (result.exit_boundary_set.size() == 2) {
            result.exit_boundary_set.clear();
            if (!fiducial_utils->inside_fiducial_volume(bp1)) result.exit_boundary_set.insert(0);
            if (!fiducial_utils->inside_fiducial_volume(bp2)) result.exit_boundary_set.insert(1);
            if (result.exit_boundary_set.empty()) {
                result.exit_boundary_set.insert(0);
                result.exit_boundary_set.insert(1);
            }
        }
    };

    // ------------------------------------------------------------------
    // Round 1: flag_cosmic = true
    // ------------------------------------------------------------------
    {
        auto boundary_indices = cluster.get_two_boundary_steiner_graph_idx(
            "steiner_graph", "steiner_pc", /*flag_cosmic=*/true);

        const auto& steiner_pc = cluster.get_pc("steiner_pc");
        const auto& coords = cluster.get_default_scope().coords;
        const auto& xc = steiner_pc.get(coords.at(0))->elements<double>();
        const auto& yc = steiner_pc.get(coords.at(1))->elements<double>();
        const auto& zc = steiner_pc.get(coords.at(2))->elements<double>();

        result.boundary_first  = geo_point_t(xc[boundary_indices.first],
                                              yc[boundary_indices.first],
                                              zc[boundary_indices.first]);
        result.boundary_second = geo_point_t(xc[boundary_indices.second],
                                              yc[boundary_indices.second],
                                              zc[boundary_indices.second]);

        auto out_vec_wcps = cluster.get_extreme_wcps();
        {
            std::vector<geo_point_t> tmp; tmp.push_back(result.boundary_first);
            out_vec_wcps.push_back(tmp);
        }
        {
            std::vector<geo_point_t> tmp; tmp.push_back(result.boundary_second);
            out_vec_wcps.push_back(tmp);
        }

        do_boundary_check(out_vec_wcps);
        update_boundary_set(result.boundary_first, result.boundary_second);
    }

    // ------------------------------------------------------------------
    // Round 2: flag_cosmic = false  (only when round 1 found no exit candidate)
    // ------------------------------------------------------------------
    if (result.exit_boundary_set.empty()) {
        result.exit_wcps.clear();

        if (!cluster.has_pc("steiner_pc")) {
            // No steiner_pc to run round 2 — treat as not-FC (conservative)
            return result;
        }

        auto boundary_indices2 = cluster.get_two_boundary_steiner_graph_idx(
            "steiner_graph", "steiner_pc", /*flag_cosmic=*/false);

        const auto& steiner_pc2 = cluster.get_pc("steiner_pc");
        const auto& coords2 = cluster.get_default_scope().coords;
        const auto& xc2 = steiner_pc2.get(coords2.at(0))->elements<double>();
        const auto& yc2 = steiner_pc2.get(coords2.at(1))->elements<double>();
        const auto& zc2 = steiner_pc2.get(coords2.at(2))->elements<double>();

        geo_point_t bp1_r2(xc2[boundary_indices2.first],
                           yc2[boundary_indices2.first],
                           zc2[boundary_indices2.first]);
        geo_point_t bp2_r2(xc2[boundary_indices2.second],
                           yc2[boundary_indices2.second],
                           zc2[boundary_indices2.second]);

        auto out_vec_wcps2 = cluster.get_extreme_wcps();
        {
            std::vector<geo_point_t> tmp; tmp.push_back(bp1_r2);
            out_vec_wcps2.push_back(tmp);
        }
        {
            std::vector<geo_point_t> tmp; tmp.push_back(bp2_r2);
            out_vec_wcps2.push_back(tmp);
        }

        do_boundary_check(out_vec_wcps2);
        update_boundary_set(bp1_r2, bp2_r2);
    }

    result.is_fc = result.exit_wcps.empty();
    return result;
}

}  // namespace WireCell::Clus::Facade