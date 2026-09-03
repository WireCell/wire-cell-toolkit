
#include "WireCellClus/Graphs.h"
#include "WireCellClus/IPCTransform.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellUtil/Logging.h"

#include "connect_graphs.h"

#include <algorithm>
#include <cstdlib>
#include <functional>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Graphs;
using namespace WireCell::Clus::Facade;

void Graphs::connect_graph_relaxed(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv, 
    IPCTransformSet::pointer pcts,
    Weighted::Graph& graph,
    const RelaxedFastCfg* fast)
{
    const bool use_ctpc = true;
    const auto* grouping = cluster.grouping();

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

    // doc 78 round 2: lazy walk-on-demand for busy clusters, only under the
    // "relaxed_fast" flavor (fast != nullptr) and only when this cluster's
    // component count exceeds the busy threshold.  Census mode always takes
    // the legacy order so the per-pair log stays complete.
    const bool lazy_walk = fast && !oc53_census && num > fast->busy_num_threshold;

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
        oc53_log->debug("OC53CENSUS cluster nblobs={} npoints={} ncomp={}",
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
            oc53_log->debug("OC53CENSUS comp c={} npts={} bbox=({:.1f},{:.1f},{:.1f})-({:.1f},{:.1f},{:.1f})cm",
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

    // doc 78 round 2: the per-pair closest-path walk + strong-check verdict
    // (the former in-loop block), factored into a lambda so the lazy mode
    // can evaluate it on demand.  Pure code motion -- the body is the
    // round-1 text verbatim; both call orders are per-pair independent.
    auto eval_pair_verdict = [&](size_t j, size_t k) {
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

                // doc 78 round 1: once the (monotone) counters satisfy EVERY
                // branch's kill predicate at once, the verdict below is "kill"
                // no matter what the angles or the strong-check Hough decide,
                // and the rest of the walk cannot change it.  Exact early exit;
                // disabled in census mode so the logged counters stay full-walk.
                bool forced_kill = false;

                // Check points along path
                for (int ii = 0; ii != num_steps; ii++) {
                    geo_point_t test_p(
                        p1.x() + (p2.x() - p1.x())/num_steps*(ii + 1),
                        p1.y() + (p2.y() - p1.y())/num_steps*(ii + 1),
                        p1.z() + (p2.z() - p1.z())/num_steps*(ii + 1)
                    );

                    // Test point quality using grouping parameters
                    int scores[6] = {0,0,0,0,0,0};   // stack, not heap (see the [6] overload)
                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa()!=-1){
                            geo_point_t test_p_raw = test_p;
                            if (needs_transform) {
                                test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            grouping->test_good_point(test_p_raw, test_wpid.apa(), test_wpid.face(), scores);

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
                        // Sufficient for the kill verdict under the strong branch
                        // (num_bad1[0] > 7) AND under every non-strong branch
                        // (num_bad[0] > 7 covers parallel/default; the pair sums
                        // > 9 cover the three single-wire-angle branches).  The
                        // ratio and num_bad[3]>=3 alternatives only ADD kill
                        // reasons, so they need not hold for sufficiency.
                        if (!oc53_census &&
                            num_bad1[0] > 7 && num_bad[0] > 7 &&
                            num_bad[2] + num_bad[3] > 9 &&
                            num_bad[1] + num_bad[3] > 9 &&
                            num_bad[2] + num_bad[1] > 9) {
                            forced_kill = true;
                            break;
                        }
                    }
                }

                bool flag_strong_check = true;

                if (forced_kill) {
                    index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
                }
                else {

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

                // Define constants for readability
                constexpr double pi = 3.141592653589793;
                constexpr double perp_angle_tol = 10.0/180.0*pi;
                constexpr double wire_angle_tol = 12.5/180.0*pi;
                constexpr double perp_angle = pi/2.0;
                constexpr double invalid_dist = 1e9;

                // Helper function to check if ratio exceeds threshold
                auto exceeds_ratio = [](int val, int steps, double ratio = 0.75) {
                    return val >= ratio * steps;
                };

                // doc 78 round 1: the flag_strong_check machinery only picks WHICH
                // of two verdicts applies, so evaluate both verdicts first and call
                // the two 15 cm whole-cluster Hough transforms only when they
                // disagree -- when the verdicts agree the Houghs provably cannot
                // change the outcome.  The predicate ladder below is a verbatim
                // transcription of the pre-round-1 branches.

                // Verdict if the strong check applies (flag_strong_check == true).
                const bool kill_strong =
                    num_bad1[0] > 7 || (num_bad1[0] > 2 && exceeds_ratio(num_bad1[0], num_steps));

                // Verdict if the strong check does not apply.
                bool kill_nonstrong;
                {
                    bool parallel_angles = (angle1 < wire_angle_tol && angle2 < wire_angle_tol) ||
                                        (angle1p < wire_angle_tol && angle1 < wire_angle_tol) ||
                                        (angle1p < wire_angle_tol && angle2 < wire_angle_tol);

                    if (parallel_angles) {
                        kill_nonstrong = num_bad[0] > 7 || (num_bad[0] > 2 && exceeds_ratio(num_bad[0], num_steps));
                    }
                    else if (angle1 < wire_angle_tol) {
                        int sum_bad = num_bad[2] + num_bad[3];
                        kill_nonstrong = sum_bad > 9 || (sum_bad > 2 && exceeds_ratio(sum_bad, num_steps)) || num_bad[3] >= 3;
                    }
                    else if (angle2 < wire_angle_tol) {
                        int sum_bad = num_bad[1] + num_bad[3];
                        kill_nonstrong = sum_bad > 9 || (sum_bad > 2 && exceeds_ratio(sum_bad, num_steps)) || num_bad[3] >= 3;
                    }
                    else if (angle1p < wire_angle_tol) {
                        int sum_bad = num_bad[2] + num_bad[1];
                        kill_nonstrong = sum_bad > 9 || (sum_bad > 2 && exceeds_ratio(sum_bad, num_steps));
                    }
                    else {
                        kill_nonstrong = num_bad[0] > 7 || (num_bad[0] > 2 && exceeds_ratio(num_bad[0], num_steps));
                    }
                }

                bool kill_now;
                if (fabs(angle3 - perp_angle) < perp_angle_tol) {
                    if (kill_strong == kill_nonstrong && !oc53_census) {
                        // Houghs skipped: both verdicts agree, the flag is moot.
                        // Census mode keeps the legacy path so the logged strong=
                        // flag stays faithful.
                        kill_now = kill_strong;
                    }
                    else {
                        geo_vector_t tempV2 = cluster.vhough_transform(p1, 15*units::cm);
                        geo_vector_t tempV3 = cluster.vhough_transform(p2, 15*units::cm);

                        if (fabs(tempV2.angle(drift_dir_abs) - perp_angle) < perp_angle_tol &&
                            fabs(tempV3.angle(drift_dir_abs) - perp_angle) < perp_angle_tol) {
                            flag_strong_check = false;
                        }
                        kill_now = flag_strong_check ? kill_strong : kill_nonstrong;
                    }
                }
                else if (angle1 < wire_angle_tol || angle2 < wire_angle_tol || angle1p < wire_angle_tol) {
                    flag_strong_check = false;
                    kill_now = kill_nonstrong;
                }
                else {
                    kill_now = kill_strong;
                }

                if (kill_now) {
                    index_index_dis[j][k] = std::make_tuple(-1, -1, invalid_dist);
                }
                }  // end !forced_kill

                if (oc53_census) {
                    const bool killed = std::get<0>(index_index_dis[j][k]) < 0;
                    oc53_log->debug(
                        "OC53CENSUS closest j={} k={} p1=({:.2f},{:.2f},{:.2f}) p2=({:.2f},{:.2f},{:.2f}) "
                        "dis={:.2f}cm nsteps={} strong={} nb=[{},{},{},{}] nb1=[{},{},{},{}] killed={}",
                        j, k, p1.x()/units::cm, p1.y()/units::cm, p1.z()/units::cm,
                        p2.x()/units::cm, p2.y()/units::cm, p2.z()/units::cm,
                        dis/units::cm, num_steps, flag_strong_check,
                        num_bad[0], num_bad[1], num_bad[2], num_bad[3],
                        num_bad1[0], num_bad1[1], num_bad1[2], num_bad1[3], killed);
                }
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
            // Now check the path (legacy: every pair; lazy mode defers to the
            // MST fixpoint below -- doc 78 round 2)
            if (!lazy_walk) {
                eval_pair_verdict(j, k);
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
                bool forced_kill = false;

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
                                // doc 78 round 1: is_good_point with allowed_bad=0 is
                                // strictly tighter than with the default 1, so a point
                                // failing the loose test cannot pass the tight one --
                                // count it directly and skip the second kd descent.
                                num_bad1++;
                            }
                            else if (!grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face(), 0.6*units::cm, 1, 0)) {
                                num_bad1++;
                            }
                        } else {
                            // Step is outside all APA volumes — count as bad.
                            num_bad++;
                            num_bad1++;
                        }
                        // doc 78 round 1: num_bad1 >= num_bad always (tighter test),
                        // so num_bad > 7 satisfies BOTH branch predicates below --
                        // the verdict is already "kill"; the rest of the walk cannot
                        // change it.  Census mode keeps full-walk counters.
                        if (!oc53_census && num_bad > 7) {
                            forced_kill = true;
                            break;
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
                if (forced_kill) {
                    // doc 78 round 1: verdict already forced during the walk.
                    index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                }
                else if (fabs(angle3 - pi/2) < 10.0/180.0*pi || 
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

                if (oc53_census) {
                    const bool killed = std::get<0>(index_index_dis_dir1[j][k]) < 0;
                    oc53_log->debug("OC53CENSUS dir1 j={} k={} dis={:.2f}cm nsteps={} nb={} nb1={} killed={}",
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
                bool forced_kill = false;

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
                                // doc 78 round 1: is_good_point with allowed_bad=0 is
                                // strictly tighter than with the default 1, so a point
                                // failing the loose test cannot pass the tight one --
                                // count it directly and skip the second kd descent.
                                num_bad1++;
                            }
                            else if (!grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face(), 0.6*units::cm, 1, 0)) {
                                num_bad1++;
                            }
                        } else {
                            // Step is outside all APA volumes — count as bad.
                            num_bad++;
                            num_bad1++;
                        }
                        // doc 78 round 1: num_bad1 >= num_bad always (tighter test),
                        // so num_bad > 7 satisfies BOTH branch predicates below --
                        // the verdict is already "kill"; the rest of the walk cannot
                        // change it.  Census mode keeps full-walk counters.
                        if (!oc53_census && num_bad > 7) {
                            forced_kill = true;
                            break;
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

                if (forced_kill) {
                    // doc 78 round 1: verdict already forced during the walk.
                    index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                }
                else if (is_parallel) {
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

                if (oc53_census) {
                    const bool killed = std::get<0>(index_index_dis_dir2[j][k]) < 0;
                    oc53_log->debug("OC53CENSUS dir2 j={} k={} dis={:.2f}cm nsteps={} nb={} nb1={} killed={}",
                                    j, k, dis/units::cm, num_steps, num_bad, num_bad1, killed);
                }
            }
        }
    }

    // deal with MST of first type
    if (lazy_walk) {
        // doc 78 round 2 (take 2 -- the kill-and-rebuild Prim fixpoint of the
        // first attempt converged too slowly: killed bridges cluster, so each
        // rebuild proposes more killable edges and the iteration cap then
        // walked everything, costing MORE than eager).  Single-pass lazy
        // Kruskal instead: consider pairs in ascending (distance, j, k)
        // order and evaluate the walk verdict ONLY when a pair would bridge
        // two union-find components.  Killed bridges are skipped; the
        // accepted set IS the spanning forest of the walk-filtered graph
        // (identical to the legacy per-component Prim except possibly on
        // exact distance ties, which is inside this knob's validation).
        // The < 3 cm direct-edge rule below needs every near pair's verdict
        // regardless of the forest -- walk those eagerly (a few steps each).
        std::vector<std::vector<char>> walked(num, std::vector<char>(num, 0));
        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                if (std::get<0>(index_index_dis[j][k]) < 0) continue;
                // Near pairs feed the < 3 cm direct-edge rule; pairs with a
                // directional candidate feed MST #2's recording of the MAIN
                // entry (process_mst_deterministically copies
                // index_index_dis[e], and a killed main blocks the dir edges
                // downstream) -- both need the main verdict eagerly.
                if (std::get<2>(index_index_dis[j][k]) < 3 * units::cm ||
                    std::get<0>(index_index_dis_dir1[j][k]) >= 0 ||
                    std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    eval_pair_verdict(j, k);
                    walked[j][k] = 1;
                }
            }
        }

        std::vector<std::tuple<double, size_t, size_t>> order;
        order.reserve(num * (num - 1) / 2);
        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                if (std::get<0>(index_index_dis[j][k]) >= 0) {
                    order.emplace_back(std::get<2>(index_index_dis[j][k]), j, k);
                }
            }
        }
        std::sort(order.begin(), order.end());

        std::vector<size_t> uf(num);
        for (size_t i = 0; i != num; ++i) uf[i] = i;
        std::function<size_t(size_t)> uf_find = [&](size_t x) {
            while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
            return x;
        };

        size_t n_walked = 0, n_accepted = 0;
        for (const auto& [d, j, k] : order) {
            const size_t rj = uf_find(j), rk = uf_find(k);
            if (rj == rk) continue;
            if (!walked[j][k]) {
                eval_pair_verdict(j, k);
                walked[j][k] = 1;
                ++n_walked;
            }
            if (std::get<0>(index_index_dis[j][k]) < 0) continue;  // killed bridge
            uf[rj] = rk;
            index_index_dis_mst[j][k] = index_index_dis[j][k];
            if (++n_accepted + 1 == num) break;  // spanning forest complete
        }
        SPDLOG_LOGGER_DEBUG(WireCell::Log::logger("clus"),
            "connect_graph_relaxed lazy: num={} pairs={} walked={} accepted={}",
            num, order.size(), n_walked, n_accepted);
    }

    if (!lazy_walk) {
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
        oc53_log->debug("OC53CENSUS summary ncomp={} nedges={}", num, oc53_pairs.size());
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
            oc53_log->debug("OC53CENSUS edge j={} k={} dis={:.2f}cm mst={} dir1={} dir2={} bridge={}",
                            ja, ka, std::get<2>(index_index_dis[ja][ka])/units::cm,
                            e_mst, e_dir1, e_dir2, !seen[ka]);
        }
    }
}

 bool Graphs::check_connectivity(const Facade::Cluster& cluster,  IDetectorVolumes::pointer dv, IPCTransformSet::pointer pcts, std::tuple<int, int, double>& index_index_dis, std::shared_ptr<Facade::Simple3DPointCloud> pc1, std::vector<size_t> pc1_global_index, std::shared_ptr<Facade::Simple3DPointCloud> pc2, std::vector<size_t> pc2_global_index,
    double step_size, bool flag_strong_check){
    
    // Check if indices are valid
    if (std::get<0>(index_index_dis) == -1 || std::get<1>(index_index_dis) == -1) return false;
    
    // Get points from cluster
    // const auto& points = cluster.points();
    
    // Get the two points from the point cloud
    int idx1 = std::get<0>(index_index_dis);
    int idx2 = std::get<1>(index_index_dis);
    
    geo_point_t p1= pc1->point(idx1);
    geo_point_t p2= pc2->point(idx2);
    
    // Get wire plane IDs for the two points
    auto wpid_p1 = cluster.wire_plane_id(pc1_global_index.at(idx1));
    auto wpid_p2 = cluster.wire_plane_id(pc2_global_index.at(idx2));
    auto wpid_pc = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv);
    
    int apa1 = wpid_p1.apa();
    int face1 = wpid_p1.face();
    int apa2 = wpid_p2.apa();
    int face2 = wpid_p2.face();
    int apa3 = wpid_pc.apa();
    int face3 = wpid_pc.face();
    
    // Get grouping for CTPC checks
    const auto* grouping = cluster.grouping();
    
    // Calculate directions using VHoughTrans equivalent (vhough_transform)
    // Use point clouds pc1 and pc2 for local direction calculation
    geo_vector_t dir1 = cluster.vhough_transform(p1, 15*units::cm, Cluster::HoughParamSpace::theta_phi, pc1, pc1_global_index);
    dir1 = dir1 * -1;
    
    geo_vector_t dir2 = cluster.vhough_transform(p2, 15*units::cm, Cluster::HoughParamSpace::theta_phi, pc2, pc2_global_index);
    dir2 = dir2 * -1;
    
    geo_vector_t dir3(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
    
    // Check directions using the check_direction function
    std::vector<bool> flag_1 = check_direction(cluster, dir1, apa1, face1);
    std::vector<bool> flag_2 = check_direction(cluster, dir2, apa2, face2);
    
    // For dir3, use either apa/face from p1 or p2 (use p1 as reference)
    std::vector<bool> flag_3 = check_direction(cluster, dir3, apa3, face3);
    
    bool flag_prolonged_u = false;
    bool flag_prolonged_v = false;
    bool flag_prolonged_w = false;
    bool flag_parallel = false;
    
    // Check if prolonged along wire directions or parallel to drift
    if (flag_3.at(0) && (flag_1.at(0) || flag_2.at(0))) flag_prolonged_u = true;
    if (flag_3.at(1) && (flag_1.at(1) || flag_2.at(1))) flag_prolonged_v = true;
    if (flag_3.at(2) && (flag_1.at(2) || flag_2.at(2))) flag_prolonged_w = true;
    if (flag_3.at(3) && (flag_1.at(3) && flag_2.at(3))) flag_parallel = true;
    
    // Calculate distance and number of steps
    double dis = std::sqrt(std::pow(p1.x() - p2.x(), 2) + 
                          std::pow(p1.y() - p2.y(), 2) + 
                          std::pow(p1.z() - p2.z(), 2));
    int num_steps = std::round(dis / step_size);
    
    if (num_steps == 0) num_steps = 1;
    
    int num_bad[5] = {0, 0, 0, 0, 0};
    
    double radius_cut = 0.6 * units::cm;
    if (step_size < radius_cut) radius_cut = step_size;

    // Hoist scope-transform out of the per-step loop
    const bool cc_needs_transform = (cluster.get_default_scope().hash() != cluster.get_raw_scope().hash());
    const auto cc_ctpc_transform = cc_needs_transform ? pcts->pc_transform(cluster.get_scope_transform()) : nullptr;
    const double cc_cluster_t0 = cc_needs_transform ? cluster.get_cluster_t0() : 0.0;
    
    // Check points along the path
    for (int i = 0; i != num_steps; i++) {
        geo_point_t test_p(
            p1.x() + (p2.x() - p1.x()) / (num_steps + 1.0) * (i + 1),
            p1.y() + (p2.y() - p1.y()) / (num_steps + 1.0) * (i + 1),
            p1.z() + (p2.z() - p1.z()) / (num_steps + 1.0) * (i + 1)
        );
        
        // Get wire plane ID for test point
        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);

        if (test_wpid.apa() == -1) {
            // Step is outside all APA volumes — count as bad on all planes.
            num_bad[0]++; num_bad[1]++; num_bad[2]++; num_bad[3]++;
            continue;
        }
        
        // Transform point if needed
        geo_point_t test_p_raw = test_p;
        if (cc_needs_transform) {
            test_p_raw = cc_ctpc_transform->backward(test_p, cc_cluster_t0, test_wpid.face(), test_wpid.apa());
        }
        
        // Test point quality with appropriate radius
        double test_radius;
        if (i == 0 || i + 1 == num_steps) {
            test_radius = dis / (num_steps + 1.0) * 0.98;
        } else {
            if (flag_strong_check) {
                test_radius = dis / (num_steps + 1.0);
            } else {
                test_radius = radius_cut;
            }
        }
        
        // Get detailed scores for this point
        int scores[6];
        grouping->test_good_point(test_p_raw, test_wpid.apa(), test_wpid.face(), scores, test_radius);
        
        int num_bad_details = 0;
        
        // Check U plane (indices 0=live, 3=dead)
        if (scores[0] + scores[3] == 0) {
            if (!flag_prolonged_u) num_bad[0]++;
            num_bad_details++;
        }
        
        // Check V plane (indices 1=live, 4=dead)
        if (scores[1] + scores[4] == 0) {
            if (!flag_prolonged_v) num_bad[1]++;
            num_bad_details++;
        }
        
        // Check W plane (collection, indices 2=live, 5=dead)
        if (scores[2] + scores[5] == 0) {
            if (!flag_prolonged_w) num_bad[2]++;
            num_bad_details++;
        }
        
        // Count overall bad points
        if (flag_parallel) {
            // Parallel case: more than one plane bad
            if (num_bad_details > 1) num_bad[3]++;
        } else {
            // Non-parallel: any plane bad
            if (num_bad_details > 0) num_bad[3]++;
        }
    }
    
    // Strong check - very strict criteria
    if (flag_strong_check && ((num_bad[0] + num_bad[1] + num_bad[2]) > 0 || num_bad[3] >= 2)) {
        return false;
    }
    
    // Prolonged case - allow some bad points but not too many
    if (num_bad[0] <= 2 && num_bad[1] <= 2 && num_bad[2] <= 2 &&
        (num_bad[0] + num_bad[1] + num_bad[2] <= 3) && 
        num_bad[0] < 0.1 * num_steps && 
        num_bad[1] < 0.1 * num_steps && 
        num_bad[2] < 0.1 * num_steps &&
        (num_bad[0] + num_bad[1] + num_bad[2]) < 0.15 * num_steps) {
        
        // Special case: if prolonged in all three directions, check overall quality
        if (flag_prolonged_u && flag_prolonged_v && flag_prolonged_w) {
            if (num_bad[3] >= 0.6 * num_steps) return false;
        }
        
        return true;
    } 
    // Alternative case - overall good quality
    else if (num_bad[3] <= 2 && num_bad[3] < 0.1 * num_steps) {
        return true;
    }
    
    return false;
 }

void Graphs::connect_graph_relaxed_pid(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv, 
        IPCTransformSet::pointer pcts,
        Weighted::Graph& graph){
    
    // const auto* grouping = cluster.grouping();
    const geo_vector_t drift_dir_abs(1, 0, 0);
    
    // Form connected components
    std::vector<int> component(num_vertices(graph));
    const size_t num = connected_components(graph, &component[0]);
    
    if (num <= 1) return;
    
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
    
    // Initialize distance metrics
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis(num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir1(num, std::vector<std::tuple<int, int, double>>(num));
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir2(num, std::vector<std::tuple<int, int, double>>(num));
    
    // Initialize all distances to inf
    for (size_t j = 0; j != num; j++) {
        for (size_t k = 0; k != num; k++) {
            index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
        }
    }
    
    // Calculate distances between components with connectivity checks
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            // Get closest points between components
            std::tuple<int, int, double> temp_index_index_dis = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));
            
            if (std::get<0>(temp_index_index_dis) != -1) {
                index_index_dis[j][k] = temp_index_index_dis;
                
                // Check connectivity
                bool flag = check_connectivity(cluster, dv, pcts, index_index_dis[j][k], 
                                              pt_clouds.at(j), pt_clouds_global_indices.at(j), pt_clouds.at(k), pt_clouds_global_indices.at(k));
                
                // Special case for very close distances with large point clouds
                if (std::get<2>(temp_index_index_dis) <= 0.9 * units::cm && 
                    pt_clouds.at(j)->get_num_points() > 200 && 
                    pt_clouds.at(k)->get_num_points() > 200) {
                    
                    if (!flag) {
                        // Try to find better connection points nearby
                        geo_point_t test_p1 = pt_clouds.at(k)->point(std::get<1>(temp_index_index_dis));
                        geo_point_t test_p2 = pt_clouds.at(j)->point(std::get<0>(temp_index_index_dis));
                        
                        auto temp_wcps1 = 
                            pt_clouds.at(j)->get_closest_wcpoints_radius(test_p1, 
                                std::get<2>(temp_index_index_dis) + 0.9 * units::cm);
                        auto temp_wcps2 = 
                            pt_clouds.at(k)->get_closest_wcpoints_radius(test_p2, 
                                std::get<2>(temp_index_index_dis) + 0.9 * units::cm);
                        
                        // Try different point combinations
                        for (size_t kk1 = 0; kk1 < temp_wcps1.size() && !flag; kk1++) {
                            for (size_t kk2 = 0; kk2 < temp_wcps2.size() && !flag; kk2++) {
                                double dis = std::sqrt(
                                    std::pow(temp_wcps1[kk1].second.x() - temp_wcps2[kk2].second.x(), 2) +
                                    std::pow(temp_wcps1[kk1].second.y() - temp_wcps2[kk2].second.y(), 2) +
                                    std::pow(temp_wcps1[kk1].second.z() - temp_wcps2[kk2].second.z(), 2));
                                
                                std::tuple<int, int, double> temp_tuple = 
                                    std::make_tuple(temp_wcps1[kk1].first, temp_wcps2[kk2].first, dis);
                                
                                if (check_connectivity(cluster, dv, pcts, temp_tuple, 
                                                      pt_clouds.at(j), pt_clouds_global_indices.at(j), pt_clouds.at(k), pt_clouds_global_indices.at(k), 
                                                      0.3 * units::cm, true)) {
                                    flag = true;
                                    index_index_dis[j][k] = temp_tuple;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if (!flag) {
                    index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
                }
                index_index_dis[k][j] = index_index_dis[j][k];
                
                // Calculate directional connections
                if (std::get<0>(temp_index_index_dis) != -1) {
                    geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(temp_index_index_dis));
                    geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(temp_index_index_dis));
                    
                    // Direction from p1
                    geo_vector_t dir1 = cluster.vhough_transform(p1, 30 * units::cm, 
                        Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                    dir1 = dir1 * -1;
                    
                    std::pair<int, double> result1 = pt_clouds.at(k)->get_closest_point_along_vec(
                        p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                    
                    // If no result and perpendicular to drift, try longer hough
                    if (result1.first < 0 && 
                        std::fabs(dir1.angle(drift_dir_abs) * 180.0 / M_PI - 90.0) < 10.0) {
                        if (std::fabs(dir1.angle(drift_dir_abs) * 180.0 / M_PI - 90.0) < 5.0)
                            dir1 = cluster.vhough_transform(p1, 80 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                        else if (std::fabs(dir1.angle(drift_dir_abs) * 180.0 / M_PI - 90.0) < 10.0)
                            dir1 = cluster.vhough_transform(p1, 50 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                        dir1 = dir1 * -1;
                        result1 = pt_clouds.at(k)->get_closest_point_along_vec(
                            p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                    }
                    
                    if (result1.first >= 0) {
                        index_index_dis_dir1[j][k] = std::make_tuple(
                            std::get<0>(index_index_dis[j][k]), result1.first, result1.second);
                        
                        if (!check_connectivity(cluster, dv, pcts, index_index_dis_dir1[j][k], 
                                               pt_clouds.at(j), pt_clouds_global_indices.at(j), pt_clouds.at(k), pt_clouds_global_indices.at(k))) {
                            index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                        }
                        index_index_dis_dir1[k][j] = index_index_dis_dir1[j][k];
                    }
                    
                    // Direction from p2
                    geo_vector_t dir2 = cluster.vhough_transform(p2, 30 * units::cm, 
                        Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                    dir2 = dir2 * -1;
                    
                    std::pair<int, double> result2 = pt_clouds.at(j)->get_closest_point_along_vec(
                        p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                    
                    if (result2.first < 0 && 
                        std::fabs(dir2.angle(drift_dir_abs) * 180.0 / M_PI - 90.0) < 10.0) {
                        if (std::fabs(dir2.angle(drift_dir_abs) * 180.0 / M_PI - 90.0) < 5.0)
                            dir2 = cluster.vhough_transform(p2, 80 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                        else if (std::fabs(dir2.angle(drift_dir_abs) * 180.0 / M_PI - 90.0) < 10.0)
                            dir2 = cluster.vhough_transform(p2, 50 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                        dir2 = dir2 * -1;
                        result2 = pt_clouds.at(j)->get_closest_point_along_vec(
                            p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                    }
                    
                    if (result2.first >= 0) {
                        index_index_dis_dir2[j][k] = std::make_tuple(
                            result2.first, std::get<1>(index_index_dis[j][k]), result2.second);
                        
                        if (!check_connectivity(cluster, dv, pcts, index_index_dis_dir2[j][k], 
                                               pt_clouds.at(j), pt_clouds_global_indices.at(j), pt_clouds.at(k), pt_clouds_global_indices.at(k))) {
                            index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                        }
                        index_index_dis_dir2[k][j] = index_index_dis_dir2[j][k];
                    }
                }
            }
        }
    }
    
    // Examine middle path for all three connection types
    double step_dis = 1.0 * units::cm;
    
    auto examine_middle_path = [&](std::vector<std::vector<std::tuple<int, int, double>>>& index_dis_array, bool apply_size_filter) {
        std::map<std::pair<int, int>, std::set<int>> map_add_connections;
        
        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                if (std::get<0>(index_dis_array[j][k]) >= 0) {
                    int idx1 = std::get<0>(index_dis_array[j][k]);
                    int idx2 = std::get<1>(index_dis_array[j][k]);
                    
                    geo_point_t wp1(points[0][pt_clouds_global_indices[j][idx1]], 
                                   points[1][pt_clouds_global_indices[j][idx1]], 
                                   points[2][pt_clouds_global_indices[j][idx1]]);
                    geo_point_t wp2(points[0][pt_clouds_global_indices[k][idx2]], 
                                   points[1][pt_clouds_global_indices[k][idx2]], 
                                   points[2][pt_clouds_global_indices[k][idx2]]);
                    
                    double length = std::sqrt(
                        std::pow(wp1.x() - wp2.x(), 2) + 
                        std::pow(wp1.y() - wp2.y(), 2) + 
                        std::pow(wp1.z() - wp2.z(), 2));
                    
                    if (length > 3 * units::cm) {
                        std::set<int> connections;
                        int ncount = std::round(length / step_dis);
                        
                        for (int qx = 1; qx < ncount; qx++) {
                            geo_point_t test_p(
                                wp1.x() + (wp2.x() - wp1.x()) * qx / ncount,
                                wp1.y() + (wp2.y() - wp1.y()) * qx / ncount,
                                wp1.z() + (wp2.z() - wp1.z()) * qx / ncount);
                            
                            for (size_t qx1 = 0; qx1 != num; qx1++) {
                                if (qx1 == j || qx1 == k) continue;
                                
                                if (pt_clouds.at(qx1)->get_closest_dis(test_p) < 0.6 * units::cm &&
                                    (!apply_size_filter || pt_clouds.at(qx1)->get_num_points() >= 50)) {
                                    connections.insert(qx1);
                                }
                            }
                        }
                        
                        if (!connections.empty()) {
                            map_add_connections[std::make_pair(j, k)] = connections;
                        }
                    }
                }
            }
        }
        
        // Iteratively disconnect paths blocked by intermediate components
        bool flag_continue = true;
        while (flag_continue) {
            flag_continue = false;
            std::set<std::pair<int, int>> used_pairs;
            
            for (auto it = map_add_connections.begin(); it != map_add_connections.end(); it++) {
                int j = it->first.first;
                int k = it->first.second;
                bool flag_disconnect = true;
                
                for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++) {
                    int qx = *it1;
                    if ((std::get<0>(index_index_dis[j][qx]) != -1 || 
                         std::get<0>(index_index_dis_dir1[j][qx]) != -1 || 
                         std::get<0>(index_index_dis_dir2[j][qx]) != -1) &&
                        (std::get<0>(index_index_dis[k][qx]) != -1 || 
                         std::get<0>(index_index_dis_dir1[k][qx]) != -1 || 
                         std::get<0>(index_index_dis_dir2[k][qx]) != -1)) {
                        flag_disconnect = false;
                        break;
                    }
                }
                
                if (flag_disconnect) {
                    flag_continue = true;
                    index_dis_array[j][k] = std::make_tuple(-1, -1, 1e9);
                    index_dis_array[k][j] = index_dis_array[j][k];
                    used_pairs.insert(it->first);
                }
            }
            
            for (auto it = used_pairs.begin(); it != used_pairs.end(); it++) {
                map_add_connections.erase(*it);
            }
        }
    };
    
    // Examine all three connection types
    examine_middle_path(index_index_dis,      true);
    examine_middle_path(index_index_dis_dir1, false);
    examine_middle_path(index_index_dis_dir2, false);
    
    // Final assembly: add edges to graph
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            // Add closest distance connections
            if (std::get<0>(index_index_dis[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis[j][k]));
                
                if (!boost::edge(gind1, gind2, graph).second) {
                    add_edge(gind1, gind2, std::get<2>(index_index_dis[j][k]), graph);
                }
            }
            
            // Add directional connection 1
            if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k]));
                
                float dis;
                if (std::get<2>(index_index_dis_dir1[j][k]) > 5 * units::cm) {
                    dis = std::get<2>(index_index_dis_dir1[j][k]) * 1.2;
                } else {
                    dis = std::get<2>(index_index_dis_dir1[j][k]);
                }
                
                if (!boost::edge(gind1, gind2, graph).second) {
                    add_edge(gind1, gind2, dis, graph);
                }
            }
            
            // Add directional connection 2
            if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k]));
                
                float dis;
                if (std::get<2>(index_index_dis_dir2[j][k]) > 5 * units::cm) {
                    dis = std::get<2>(index_index_dis_dir2[j][k]) * 1.2;
                } else {
                    dis = std::get<2>(index_index_dis_dir2[j][k]);
                }
                
                if (!boost::edge(gind1, gind2, graph).second) {
                    add_edge(gind1, gind2, dis, graph);
                }
            }
        }
    }
}