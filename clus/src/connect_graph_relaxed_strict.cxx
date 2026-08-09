// doc pr/53 sec 16: "relaxed_strict" -- a stricter fork of
// connect_graph_relaxed (fork by full duplication, production file
// untouched), selected only by ClusteringProtectBundle via its graph_name
// config key.  It closes the sub-3cm arithmetic blind spot of the relaxed
// path test (doc pr/53 sec 6): the 1cm sampling loop's last sample is p2
// itself and good by construction, so num_bad <= num_steps-1 and the
// legacy "num_bad > 2" floor is unreachable for necks under ~3cm.
// Differences from connect_graph_relaxed, all inside the kill tests
// (relaxed_strict_bad below); every other line (candidate generation, MST
// construction, <3cm restore, edge emission) is a verbatim copy:
//   S1  ratio denominator becomes the interior step count (num_steps-1);
//   S2  the "> 2" floors become ">= 2" (sum-bad floors "> 2" likewise);
//   S3  S1+S2 applied to all three path-check blocks (closest-pair,
//       dir1, dir2), since their edges are tested and MST'd independently.
// Caps (>7 single-counter, >9 sum) and the num_bad[3] >= 3 dead-W veto are
// unchanged.
//
// doc pr/53 round 7 sec 18: image_check (default false, "relaxed_strict"
// behavior byte-for-byte unchanged) adds S5, an independent OR-kill next to
// relaxed_strict_bad in all three path-check blocks: relaxed_img_bad() on
// the longest contiguous run of "ghost" steps -- 1cm samples with no 3D
// image point (this cluster's own pt_clouds, any component) within
// m_img_radius and not dead-channel-excused on any plane.  Root cause:
// relaxed_strict_bad only ever sees three INDEPENDENT 2D per-plane
// projections (Facade_Grouping.cxx has_closest_point), never intersected in
// 3D, so three planes can each see charge from a different nearby track and
// the whole straight path reads "good" with nothing physically there.  The
// predicate, its run-length floor and its edge-length cap are justified in
// WireCellClus/Graphs.h and doc pr/53 round 7 sec 18.2 (offline scan against
// the owner's flagged edges + all round-6 mover edges) -- see there before
// changing the operating point.  "relaxed_strict_img"
// (make_graph_relaxed_strict_img) is the only caller passing image_check=true.

#include "WireCellClus/Graphs.h"
#include "WireCellClus/IPCTransform.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellUtil/Logging.h"

#include "connect_graphs.h"

#include <cstdlib>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Graphs;
using namespace WireCell::Clus::Facade;

bool Graphs::relaxed_strict_bad(int nbad, int num_steps, int cap)
{
    // S1: the sampling loop tests fractions (ii+1)/num_steps, so the last
    // sample is the far endpoint p2 -- a cluster charge point, good by
    // construction.  Only num_steps-1 samples can possibly be bad.
    const int interior = num_steps - 1;
    // S2: floor ">= 2" (legacy "> 2"), cap unchanged.
    return nbad > cap || (nbad >= 2 && nbad >= 0.75 * interior);
}

bool Graphs::relaxed_img_bad(int max_ghost_run, double dis_cm, int run_floor, double dis_cap_cm)
{
    // doc pr/53 round 7 sec 18.2: operating point chosen from an offline
    // scan against all 27 round-6 mover events' emitted closest-pair edges,
    // not guessed -- see WireCellClus/Graphs.h for the full justification.
    return max_ghost_run >= run_floor && dis_cm < dis_cap_cm;
}

void Graphs::connect_graph_relaxed_strict(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv,
    IPCTransformSet::pointer pcts,
    Weighted::Graph& graph,
    bool image_check)
{
    const bool use_ctpc = true;
    const auto* grouping = cluster.grouping();
    // doc pr/53 round 7 sec 18: S5 operating point (see Graphs.h for the
    // radius justification -- matches the OFFLINE oc53_probe.Loader img_dis
    // radius used in the threshold scan). Only read when image_check.
    const double m_img_radius = 1.0 * units::cm;

    // Get all the wire plane IDs from the grouping
    const auto& wpids = grouping->wpids();

    // Key: {apa, face} pair.  Each apa/face has fixed U/V/W wire directions that
    // are the same regardless of which layer's WirePlaneId is used for the lookup.
    // Using a full WirePlaneId (which includes the layer) as the map key would cause
    // std::out_of_range when get_wireplaneid() returns a wpid whose layer differs
    // from those populated here.  Keying by {apa,face} avoids this entirely.
    using af_pair_t = std::pair<int,int>;
    std::map<af_pair_t, geo_point_t> af_U_dir;
    std::map<af_pair_t, geo_point_t> af_V_dir;
    std::map<af_pair_t, geo_point_t> af_W_dir;
    for (const auto& wpid : wpids) {
        int apa = wpid.apa();
        int face = wpid.face();
        af_pair_t af{apa, face};
        if (af_U_dir.count(af)) continue;  // already computed for this apa/face

        // Create canonical wpids for all three planes with this APA and face
        WirePlaneId wpid_u(kUlayer, face, apa);
        WirePlaneId wpid_v(kVlayer, face, apa);
        WirePlaneId wpid_w(kWlayer, face, apa);

        // Get wire directions for all planes
        Vector wire_dir_u = dv->wire_direction(wpid_u);
        Vector wire_dir_v = dv->wire_direction(wpid_v);
        Vector wire_dir_w = dv->wire_direction(wpid_w);

        // Calculate angles
        double angle_u = std::atan2(wire_dir_u.z(), wire_dir_u.y());
        double angle_v = std::atan2(wire_dir_v.z(), wire_dir_v.y());
        double angle_w = std::atan2(wire_dir_w.z(), wire_dir_w.y());

        af_U_dir[af] = geo_point_t(0, cos(angle_u), sin(angle_u));
        af_V_dir[af] = geo_point_t(0, cos(angle_v), sin(angle_v));
        af_W_dir[af] = geo_point_t(0, cos(angle_w), sin(angle_w));
    }

    // this drift direction is only used to calculate isochronous case, so this is OK ...
    const geo_vector_t drift_dir_abs(1, 0, 0);


    // Form connected components
    std::vector<int> component(num_vertices(graph));
    const size_t num = connected_components(graph, &component[0]);

    if (num <= 1) return;

    // doc pr/53 round 6 F0 (relaxed_edge_census): log-only diagnostic, default
    // OFF.  Env WCT_RELAXED_EDGE_CENSUS unset => no log lines and no behavior
    // change of any kind; the graph construction below is untouched either way.
    static const bool oc53_census = std::getenv("WCT_RELAXED_EDGE_CENSUS") != nullptr;
    static auto oc53_log = WireCell::Log::logger("clus");

    // Allocate exactly num point clouds (one per component)
    std::vector<std::shared_ptr<Simple3DPointCloud>> pt_clouds(num);
    std::vector<std::vector<size_t>> pt_clouds_global_indices(num);
    for (size_t c = 0; c < num; ++c) {
        pt_clouds[c] = std::make_shared<Simple3DPointCloud>();
    }

    const auto& points = cluster.points();
    for (size_t i = 0; i < component.size(); ++i) {
        size_t c = component[i];
        pt_clouds[c]->add({points[0][i], points[1][i], points[2][i]});
        pt_clouds_global_indices[c].push_back(i);
    }

    if (oc53_census) {
        oc53_log->debug("OC53CENSUS-S cluster nblobs={} npoints={} ncomp={}",
                        cluster.nchildren(), cluster.npoints(), num);
        for (size_t c = 0; c < num; ++c) {
            double lo[3] = {1e12, 1e12, 1e12}, hi[3] = {-1e12, -1e12, -1e12};
            const size_t np = pt_clouds[c]->get_num_points();
            for (size_t i = 0; i < np; ++i) {
                const auto p = pt_clouds[c]->point(i);
                const double v[3] = {p.x(), p.y(), p.z()};
                for (int d = 0; d < 3; ++d) {
                    if (v[d] < lo[d]) lo[d] = v[d];
                    if (v[d] > hi[d]) hi[d] = v[d];
                }
            }
            oc53_log->debug("OC53CENSUS-S comp c={} npts={} bbox=({:.1f},{:.1f},{:.1f})-({:.1f},{:.1f},{:.1f})cm",
                            c, np, lo[0]/units::cm, lo[1]/units::cm, lo[2]/units::cm,
                            hi[0]/units::cm, hi[1]/units::cm, hi[2]/units::cm);
        }
    }

    // Initialize distance metrics — all sentinels (-1,-1,1e9) mean "no valid connection".
    // Use direct construction to avoid a redundant zero-fill followed by an overwrite pass (C.3).
    const auto sentinel = std::make_tuple(-1, -1, 1e9);
    const std::vector<std::tuple<int,int,double>> sentinel_row(num, sentinel);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis(num, sentinel_row);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_mst(num, sentinel_row);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir1(num, sentinel_row);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir2(num, sentinel_row);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir_mst(num, sentinel_row);

    // Hoist scope-transform and cluster_t0 out of all per-step CTPC loops
    const bool needs_transform = (cluster.get_default_scope().hash() != cluster.get_raw_scope().hash());
    const auto ctpc_transform = needs_transform ? pcts->pc_transform(cluster.get_scope_transform()) : nullptr;
    const double cluster_t0 = needs_transform ? cluster.get_cluster_t0() : 0.0;

    // doc pr/53 round 7 sec 18 S5: only evaluated when image_check.  Re-walks
    // the same 1cm samples as the nb/nb1 loop below (independent pass, kept
    // separate from the already-shipped, extensively-validated nb/nb1
    // accumulation so this addition cannot perturb it).  A step is "ghost"
    // when it has no 3D image point (any of THIS cluster's own closely-
    // components -- in-family precedent connect_graph_relaxed.cxx:1062-1069,
    // relaxed_pid's cross-component check) within m_img_radius and is not
    // dead-channel-excused on any plane.  Returns the longest contiguous run.
    auto max_ghost_run = [&](const geo_point_t& p1, const geo_point_t& p2, int num_steps,
                              const WirePlaneId& wpid_p1, const WirePlaneId& wpid_p2) -> int {
        int run = 0, best = 0;
        for (int ii = 0; ii != num_steps; ii++) {
            geo_point_t test_p(
                p1.x() + (p2.x() - p1.x())/num_steps*(ii + 1),
                p1.y() + (p2.y() - p1.y())/num_steps*(ii + 1),
                p1.z() + (p2.z() - p1.z())/num_steps*(ii + 1)
            );
            bool ghost;
            auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
            if (test_wpid.apa() != -1) {
                geo_point_t test_p_raw = test_p;
                if (needs_transform) {
                    test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                }
                auto scores = grouping->test_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                double img_dis = 1e9;
                for (size_t q = 0; q != num; ++q) {
                    img_dis = std::min(img_dis, pt_clouds.at(q)->get_closest_dis(test_p));
                }
                const bool dead_excused = (scores[3] + scores[4] + scores[5]) > 0;
                ghost = (img_dis > m_img_radius) && !dead_excused;
            } else {
                // Outside all APA volumes: no dead-channel excuse is
                // available, same convention as the nb/nb1 loop's apa==-1
                // branch (unconditionally counted bad there too).
                ghost = true;
            }
            run = ghost ? run + 1 : 0;
            if (run > best) best = run;
        }
        return best;
    };

    // Calculate distances between components
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            // Get closest points between components
            index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));

            // C.4: skip the expensive Hough probes when clouds are too far apart to benefit.
            // get_closest_point_along_vec is called with max_dis=80 cm, so a closest-pair
            // distance already >= 80 cm guarantees both directional probes return nothing.
            const bool close_enough = std::get<2>(index_index_dis[j][k]) < 80 * units::cm;

            // Skip small clouds
            if (close_enough &&
                ((num < 100 && pt_clouds.at(j)->get_num_points() > 100 && pt_clouds.at(k)->get_num_points() > 100 &&
                  (pt_clouds.at(j)->get_num_points() + pt_clouds.at(k)->get_num_points()) > 400) ||
                 (pt_clouds.at(j)->get_num_points() > 500 && pt_clouds.at(k)->get_num_points() > 500))) {

                // Get closest points and calculate directions
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));

                geo_vector_t dir1 = cluster.vhough_transform(p1, 30 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), 
                                                pt_clouds_global_indices.at(j));
                geo_vector_t dir2 = cluster.vhough_transform(p2, 30 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k),
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
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis[j][k])));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis[j][k])));

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
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa()!=-1){
                            geo_point_t test_p_raw = test_p;
                            if (needs_transform) {
                                test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            scores = grouping->test_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());

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
                        } else {
                            // Step is outside all APA volumes (between APAs).
                            // Count as fully bad: no signal from any plane can validate this gap.
                            num_bad[0]++;  num_bad[1]++;  num_bad[2]++;  num_bad[3]++;
                            num_bad1[0]++; num_bad1[1]++; num_bad1[2]++; num_bad1[3]++;
                        }
                    }
                }

                auto test_wpid = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv);

                // Calculate angles between directions
                geo_vector_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
                geo_vector_t tempV5;

                double angle1 = tempV1.angle(af_U_dir.at({test_wpid.apa(), test_wpid.face()})); 
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2)) * sin(angle1),
                        0);
                angle1 = tempV5.angle(drift_dir_abs);

                double angle2 = tempV1.angle(af_V_dir.at({test_wpid.apa(), test_wpid.face()}));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2)) * sin(angle2),
                        0);
                angle2 = tempV5.angle(drift_dir_abs);

                double angle1p = tempV1.angle(af_W_dir.at({test_wpid.apa(), test_wpid.face()}));
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
                    geo_vector_t tempV2 = cluster.vhough_transform(p1, 15*units::cm);
                    geo_vector_t tempV3 = cluster.vhough_transform(p2, 15*units::cm);
                    
                    if (fabs(tempV2.angle(drift_dir_abs) - perp_angle) < perp_angle_tol &&
                        fabs(tempV3.angle(drift_dir_abs) - perp_angle) < perp_angle_tol) {
                        flag_strong_check = false;
                    }
                }
                else if (angle1 < wire_angle_tol || angle2 < wire_angle_tol || angle1p < wire_angle_tol) {
                    flag_strong_check = false;
                }

                // Helper function to invalidate distance
                auto invalidate_distance = [&]() {
                    index_index_dis[j][k] = std::make_tuple(-1, -1, invalid_dist);
                };

                // Kill tests: relaxed_strict_bad = S1+S2 (see file header);
                // branch structure identical to connect_graph_relaxed.
                if (flag_strong_check) {
                    if (relaxed_strict_bad(num_bad1[0], num_steps)) {
                        invalidate_distance();
                    }
                }
                else {
                    bool parallel_angles = (angle1 < wire_angle_tol && angle2 < wire_angle_tol) ||
                                        (angle1p < wire_angle_tol && angle1 < wire_angle_tol) ||
                                        (angle1p < wire_angle_tol && angle2 < wire_angle_tol);

                    if (parallel_angles) {
                        if (relaxed_strict_bad(num_bad[0], num_steps)) {
                            invalidate_distance();
                        }
                    }
                    else if (angle1 < wire_angle_tol) {
                        int sum_bad = num_bad[2] + num_bad[3];
                        if (relaxed_strict_bad(sum_bad, num_steps, 9) || num_bad[3] >= 3) {
                            invalidate_distance();
                        }
                    }
                    else if (angle2 < wire_angle_tol) {
                        int sum_bad = num_bad[1] + num_bad[3];
                        if (relaxed_strict_bad(sum_bad, num_steps, 9) || num_bad[3] >= 3) {
                            invalidate_distance();
                        }
                    }
                    else if (angle1p < wire_angle_tol) {
                        int sum_bad = num_bad[2] + num_bad[1];
                        if (relaxed_strict_bad(sum_bad, num_steps, 9)) {
                            invalidate_distance();
                        }
                    }
                    else if (relaxed_strict_bad(num_bad[0], num_steps)) {
                        invalidate_distance();
                    }
                }

                // doc pr/53 round 7 sec 18 S5: independent OR-kill, evaluated
                // regardless of branch (unlike relaxed_strict_bad above, the
                // 3D-image test does not depend on wire-direction angles).
                // No-op, byte-identical to round-6 "relaxed_strict" when
                // image_check is false (default).
                int ghost_run = 0;
                if (image_check) {
                    ghost_run = max_ghost_run(p1, p2, num_steps, wpid_p1, wpid_p2);
                    if (relaxed_img_bad(ghost_run, dis / units::cm)) {
                        invalidate_distance();
                    }
                    if (oc53_census) {
                        oc53_log->debug(
                            "OC53CENSUS-IMG closest j={} k={} dis={:.2f}cm nsteps={} ghost_run={} killed={}",
                            j, k, dis / units::cm, num_steps, ghost_run,
                            std::get<0>(index_index_dis[j][k]) < 0);
                    }
                }

                if (oc53_census) {
                    const bool killed = std::get<0>(index_index_dis[j][k]) < 0;
                    oc53_log->debug(
                        "OC53CENSUS-S closest j={} k={} p1=({:.2f},{:.2f},{:.2f}) p2=({:.2f},{:.2f},{:.2f}) "
                        "dis={:.2f}cm nsteps={} strong={} nb=[{},{},{},{}] nb1=[{},{},{},{}] killed={}",
                        j, k, p1.x()/units::cm, p1.y()/units::cm, p1.z()/units::cm,
                        p2.x()/units::cm, p2.y()/units::cm, p2.z()/units::cm,
                        dis/units::cm, num_steps, flag_strong_check,
                        num_bad[0], num_bad[1], num_bad[2], num_bad[3],
                        num_bad1[0], num_bad1[1], num_bad1[2], num_bad1[3], killed);
                }
            }

            // Now check path again ...
            if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir1[j][k])); 
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k])));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir1[j][k])); 
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k])));

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
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa()!=-1){
                            geo_point_t test_p_raw = test_p;
                            if (needs_transform) {
                                test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            if (!good_point) {
                                num_bad++;
                            }
                            if (!grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face(), 0.6*units::cm, 1, 0)) {
                                num_bad1++;
                            }
                        } else {
                            // Step is outside all APA volumes — count as bad.
                            num_bad++;
                            num_bad1++;
                        }
                    }
                }

                auto test_wpid = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv);

                // Calculate angles
                geo_vector_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
                geo_vector_t tempV5;

                double angle1 = tempV1.angle(af_U_dir.at({test_wpid.apa(), test_wpid.face()}));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle1),
                        0);
                angle1 = tempV5.angle(drift_dir_abs);
                
                double angle2 = tempV1.angle(af_V_dir.at({test_wpid.apa(), test_wpid.face()}));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle2),
                        0);
                angle2 = tempV5.angle(drift_dir_abs);
                
                tempV5.set(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
                double angle3 = tempV5.angle(drift_dir_abs);
                
                double angle1p = tempV1.angle(af_W_dir.at({test_wpid.apa(), test_wpid.face()}));
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
                    if (relaxed_strict_bad(num_bad, num_steps)) {
                        index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }
                else {
                    if (relaxed_strict_bad(num_bad1, num_steps)) {
                        index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }

                // doc pr/53 round 7 sec 18 S5 (see the closest-pair block for
                // the full comment). No-op when image_check is false.
                int dir1_ghost_run = 0;
                if (image_check) {
                    dir1_ghost_run = max_ghost_run(p1, p2, num_steps, wpid_p1, wpid_p2);
                    if (relaxed_img_bad(dir1_ghost_run, dis / units::cm)) {
                        index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                    if (oc53_census) {
                        oc53_log->debug("OC53CENSUS-IMG dir1 j={} k={} dis={:.2f}cm nsteps={} ghost_run={} killed={}",
                                        j, k, dis / units::cm, num_steps, dir1_ghost_run,
                                        std::get<0>(index_index_dis_dir1[j][k]) < 0);
                    }
                }

                if (oc53_census) {
                    const bool killed = std::get<0>(index_index_dis_dir1[j][k]) < 0;
                    oc53_log->debug("OC53CENSUS-S dir1 j={} k={} dis={:.2f}cm nsteps={} nb={} nb1={} killed={}",
                                    j, k, dis/units::cm, num_steps, num_bad, num_bad1, killed);
                }
            }

            //Now check path again ... 
            // Now check the path...
            if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir2[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k])));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir2[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k])));

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
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa()!=-1){
                            geo_point_t test_p_raw = test_p;
                            if (needs_transform) {
                                test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            if (!good_point) {
                                num_bad++;
                            }
                            if (!grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face(), 0.6*units::cm, 1, 0)) {
                                num_bad1++;
                            }
                        } else {
                            // Step is outside all APA volumes — count as bad.
                            num_bad++;
                            num_bad1++;
                        }
                    }
                }

                auto test_wpid = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv);

                // Calculate angles between directions
                geo_vector_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
                geo_vector_t tempV5;

                double angle1 = tempV1.angle(af_U_dir.at({test_wpid.apa(), test_wpid.face()}));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle1),
                        0);
                angle1 = tempV5.angle(drift_dir_abs);

                double angle2 = tempV1.angle(af_V_dir.at({test_wpid.apa(), test_wpid.face()}));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle2),
                        0);
                angle2 = tempV5.angle(drift_dir_abs);

                tempV5.set(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
                double angle3 = tempV5.angle(drift_dir_abs);

                double angle1p = tempV1.angle(af_W_dir.at({test_wpid.apa(), test_wpid.face()}));
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
                    if (relaxed_strict_bad(num_bad, num_steps)) {
                        index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }
                else {
                    if (relaxed_strict_bad(num_bad1, num_steps)) {
                        index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }

                // doc pr/53 round 7 sec 18 S5 (see the closest-pair block for
                // the full comment). No-op when image_check is false.
                int dir2_ghost_run = 0;
                if (image_check) {
                    dir2_ghost_run = max_ghost_run(p1, p2, num_steps, wpid_p1, wpid_p2);
                    if (relaxed_img_bad(dir2_ghost_run, dis / units::cm)) {
                        index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                    if (oc53_census) {
                        oc53_log->debug("OC53CENSUS-IMG dir2 j={} k={} dis={:.2f}cm nsteps={} ghost_run={} killed={}",
                                        j, k, dis / units::cm, num_steps, dir2_ghost_run,
                                        std::get<0>(index_index_dis_dir2[j][k]) < 0);
                    }
                }

                if (oc53_census) {
                    const bool killed = std::get<0>(index_index_dis_dir2[j][k]) < 0;
                    oc53_log->debug("OC53CENSUS-S dir2 j={} k={} dis={:.2f}cm nsteps={} nb={} nb1={} killed={}",
                                    j, k, dis/units::cm, num_steps, num_bad, num_bad1, killed);
                }
            }
        }
    }

    // deal with MST of first type
    {
        Weighted::Graph temp_graph(num);
        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j;
                int index2 = k;
                if (std::get<0>(index_index_dis[j][k]) >= 0) {
                    if (!boost::edge(index1, index2, temp_graph).second) {
                        add_edge(index1, index2, std::get<2>(index_index_dis[j][k]), temp_graph);
                    }
                }
            }
        }

        // Process MST
        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_mst);
    }

    // MST of the direction ...
    {
        Weighted::Graph temp_graph(num);

        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j;
                int index2 = k;
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0 || std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    if (!boost::edge(index1, index2, temp_graph).second) {
                        // Add edge with minimum distance from both directions
                    add_edge(
                        index1, index2,
                        std::min(std::get<2>(index_index_dis_dir1[j][k]), std::get<2>(index_index_dis_dir2[j][k])),
                        temp_graph);
                    }
                }
            }
        }

        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_dir_mst);

    }

    std::vector<std::pair<size_t, size_t>> oc53_pairs;  // emitted component pairs (census only)

    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            if (std::get<2>(index_index_dis[j][k]) < 3 * units::cm) {
                index_index_dis_mst[j][k] = index_index_dis[j][k];
            }

            // establish the path ...
            if (std::get<0>(index_index_dis_mst[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_mst[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_mst[j][k]));
                const float dis = std::get<2>(index_index_dis_mst[j][k]);
                if (!boost::edge(gind1, gind2, graph).second) {
                    /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                }
            }

            if (std::get<0>(index_index_dis_dir_mst[j][k]) >= 0) {
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k]));
                    float dis;
                    if (std::get<2>(index_index_dis_dir1[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir1[j][k]) * 1.1;
                    }
                    else {
                        dis = std::get<2>(index_index_dis_dir1[j][k]);
                    }
                    if(!boost::edge(gind1, gind2, graph).second) {
                        /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                    }
                }
                if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k]));
                    float dis;
                    if (std::get<2>(index_index_dis_dir2[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir2[j][k]) * 1.1;
                    }
                    else {
                        dis = std::get<2>(index_index_dis_dir2[j][k]);
                    }
                    // }
                    if (!boost::edge(gind1, gind2, graph).second) {
                        /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                    }
                }
            }

            if (oc53_census) {
                const bool e_mst = std::get<0>(index_index_dis_mst[j][k]) >= 0;
                const bool e_dir = std::get<0>(index_index_dis_dir_mst[j][k]) >= 0 &&
                                   (std::get<0>(index_index_dis_dir1[j][k]) >= 0 ||
                                    std::get<0>(index_index_dis_dir2[j][k]) >= 0);
                if (e_mst || e_dir) oc53_pairs.emplace_back(j, k);
            }

        }  // k
    }  // j

    if (oc53_census) {
        oc53_log->debug("OC53CENSUS-S summary ncomp={} nedges={}", num, oc53_pairs.size());
        // Bridge status of each emitted component-pair edge: DFS reachability
        // over the remaining emitted pairs.  A bridge=true edge is the sole
        // link between its two components' sides.
        for (size_t e = 0; e < oc53_pairs.size(); ++e) {
            const size_t ja = oc53_pairs[e].first;
            const size_t ka = oc53_pairs[e].second;
            std::vector<std::vector<size_t>> adj(num);
            for (size_t f = 0; f < oc53_pairs.size(); ++f) {
                if (f == e) continue;
                adj[oc53_pairs[f].first].push_back(oc53_pairs[f].second);
                adj[oc53_pairs[f].second].push_back(oc53_pairs[f].first);
            }
            std::vector<char> seen(num, 0);
            std::vector<size_t> stack{ja};
            seen[ja] = 1;
            while (!stack.empty()) {
                const size_t v = stack.back();
                stack.pop_back();
                for (const size_t w : adj[v]) {
                    if (!seen[w]) { seen[w] = 1; stack.push_back(w); }
                }
            }
            const bool e_mst = std::get<0>(index_index_dis_mst[ja][ka]) >= 0;
            const bool e_dir1 = std::get<0>(index_index_dis_dir1[ja][ka]) >= 0;
            const bool e_dir2 = std::get<0>(index_index_dis_dir2[ja][ka]) >= 0;
            oc53_log->debug("OC53CENSUS-S edge j={} k={} dis={:.2f}cm mst={} dir1={} dir2={} bridge={}",
                            ja, ka, std::get<2>(index_index_dis[ja][ka])/units::cm,
                            e_mst, e_dir1, e_dir2, !seen[ka]);
        }
    }
}
