#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Grouping.h"

#include "connect_graphs.h"

#include "WireCellUtil/Logging.h"

#include <algorithm>
#include <cstdlib>
#include <functional>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

void Graphs::connect_graph_ctpc(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv,
    Clus::IPCTransformSet::pointer pcts,
    Weighted::Graph& graph,
    const CtpcFastCfg* fast)
{
    // This used to be the body of Cluster::Connect_graph(dv,pcts,use_ctpc).
    const bool use_ctpc=true;
    const auto* grouping = cluster.grouping();

    // now form the connected components
    std::vector<int> component(num_vertices(graph));
    const size_t num = connected_components(graph, &component[0]);

    if (num <= 1) return;

    // Allocate exactly num point clouds (one per component)
    std::vector<std::shared_ptr<Simple3DPointCloud>> pt_clouds(num);
    std::vector<std::vector<size_t>> pt_clouds_global_indices(num); // can use to access wpid ...
    for (size_t c = 0; c < num; ++c) {
        pt_clouds[c] = std::make_shared<Simple3DPointCloud>();
    }

    const auto& points = cluster.points();
    for (size_t i = 0; i < component.size(); ++i) {
        size_t c = component[i];
        pt_clouds[c]->add({points[0][i], points[1][i], points[2][i]});
        pt_clouds_global_indices[c].push_back(i);
    }

    // doc 79 round 0 (WCT_CTPC_EDGE_CENSUS): log-only diagnostic, default
    // OFF.  Env unset => no log lines and no behavior change of any kind;
    // the graph construction below is untouched either way.  (Mirrors the
    // relaxed-flavor WCT_RELAXED_EDGE_CENSUS probe, doc pr/53 F0.)
    static const bool dg79_census = std::getenv("WCT_CTPC_EDGE_CENSUS") != nullptr;
    static auto dg79_log = WireCell::Log::logger("clus");
    if (dg79_census) {
        dg79_log->debug("CTPC79CENSUS cluster nblobs={} npoints={} ncomp={}",
                        cluster.nchildren(), cluster.npoints(), num);
    }
    // Aggregate walk/hough counters -- written only under dg79_census.
    // Index 0 = main pair walk, 1 = dir1 walk, 2 = dir2 walk.  "savable" =
    // steps taken after the kill predicate first became true (the waste an
    // exact early break would remove).
    size_t c79_pairs = 0, c79_hough = 0;
    size_t c79_steps[3] = {0, 0, 0}, c79_kills[3] = {0, 0, 0}, c79_savable[3] = {0, 0, 0};

    // doc 79 round 2: lazy walk-on-demand for busy clusters, only under the
    // "ctpc_fast" flavor (fast != nullptr) and only when this cluster's
    // component count exceeds the busy threshold.  Census mode always takes
    // the legacy order so the per-pair counters stay complete.
    const bool lazy_walk = fast && !dg79_census && num > fast->busy_num_threshold;

    // Initiate dist. metrics -- all sentinels (-1,-1,1e9) mean "no valid
    // connection".  doc 79 round 1: direct construction instead of a
    // zero-fill followed by an overwrite pass (exact; same fix the relaxed
    // flavor carries).
    const auto dg79_sentinel = std::make_tuple(-1, -1, 1e9);
    const std::vector<std::tuple<int, int, double>> dg79_sentinel_row(num, dg79_sentinel);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis(num, dg79_sentinel_row);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_mst(num, dg79_sentinel_row);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir1(num, dg79_sentinel_row);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir2(num, dg79_sentinel_row);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir_mst(num, dg79_sentinel_row);

    // Hoist scope-transform and cluster_t0 out of all per-step CTPC loops
    const bool needs_transform = (cluster.get_default_scope().hash() != cluster.get_raw_scope().hash());
    const auto ctpc_transform = needs_transform ? pcts->pc_transform(cluster.get_scope_transform()) : nullptr;
    const double cluster_t0 = needs_transform ? cluster.get_cluster_t0() : 0.0;

    // doc 79 round 2: the per-pair closest-path walk + kill verdict (the
    // former in-loop block), factored into a lambda so the lazy mode can
    // evaluate it on demand.  Pure code motion -- the body is the round-1
    // text verbatim; both call orders are per-pair independent.
    auto eval_main_verdict = [&](size_t j, size_t k) {
            geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
            auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis[j][k])));

            geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));
            auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis[j][k])));

            double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
            double step_dis = 1.0 * units::cm;
            int num_steps = dis / step_dis + 1;
            int num_bad = 0;
            int c79_fk = -1;   // census only: first step where the kill predicate held
            geo_point_t test_p;
            for (int ii = 0; ii != num_steps; ii++) {
                test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                           p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                           p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));

                if (use_ctpc) {
                    auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                    if (test_wpid.apa()!=-1){
                        geo_point_t test_p_raw = test_p;
                        if (needs_transform) {
                            test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                        }
                        const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                        if (!good_point) num_bad++;
                    }
                }

                if (dg79_census && c79_fk < 0 &&
                    (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps))) {
                    c79_fk = ii;
                }

                // doc 79 round 1: exact early break.  num_bad is monotone
                // nondecreasing and num_steps is fixed, so once the kill
                // predicate below holds it holds at loop end and this pair
                // is killed regardless of the remaining steps; the walk has
                // no other side effect.  Census mode walks the full path so
                // the savable counters stay complete (bit-identical output
                // either way -- the kill decision is unchanged).
                if (!dg79_census && (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps))) {
                    break;
                }
            }

            if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
            }

            if (dg79_census) {
                c79_pairs++;
                c79_steps[0] += num_steps;
                if (c79_fk >= 0) { c79_kills[0]++; c79_savable[0] += num_steps - 1 - c79_fk; }
            }
    };

    // Calc. dis, dis_dir1, dis_dir2
    // check against the closest distance ...
    // no need to have MST ...
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));

            if ((num < 100 && pt_clouds.at(j)->get_num_points() > 100 && pt_clouds.at(k)->get_num_points() > 100 &&
                 (pt_clouds.at(j)->get_num_points() + pt_clouds.at(k)->get_num_points()) > 400) ||
                (pt_clouds.at(j)->get_num_points() > 500 && pt_clouds.at(k)->get_num_points() > 500)) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));

                if (dg79_census) c79_hough++;
                geo_point_t dir1 = cluster.vhough_transform(p1, 30 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                geo_point_t dir2 = cluster.vhough_transform(p2, 30 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
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

            // Now check the path ... (doc 79 round 2: deferred under lazy_walk;
            // evaluated on demand by the Kruskal pass below)
            if (!lazy_walk) {
                eval_main_verdict(j, k);
            }

            // Now check the path ...
            if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir1[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k])));

                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir1[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis / step_dis + 1;
                int num_bad = 0;
                int c79_fk = -1;   // census only: first step where the kill predicate held
                geo_point_t test_p;
                for (int ii = 0; ii != num_steps; ii++) {
                    test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                               p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                               p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));
                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa()!=-1){
                            geo_point_t test_p_raw = test_p;
                            if (needs_transform) {
                                test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            if (!good_point) num_bad++;
                        }
                    }

                    if (dg79_census && c79_fk < 0 &&
                        (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps))) {
                        c79_fk = ii;
                    }

                    // doc 79 round 1: exact early break.  num_bad is monotone
                    // nondecreasing and num_steps is fixed, so once the kill
                    // predicate below holds it holds at loop end and this pair
                    // is killed regardless of the remaining steps; the walk has
                    // no other side effect.  Census mode walks the full path so
                    // the savable counters stay complete (bit-identical output
                    // either way -- the kill decision is unchanged).
                    if (!dg79_census && (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps))) {
                        break;
                    }
                }

                if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                    index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                }

                if (dg79_census) {
                    c79_steps[1] += num_steps;
                    if (c79_fk >= 0) { c79_kills[1]++; c79_savable[1] += num_steps - 1 - c79_fk; }
                }
            }

            // Now check the path ...
            if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir2[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k])));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir2[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis / step_dis + 1;
                int num_bad = 0;
                int c79_fk = -1;   // census only: first step where the kill predicate held
                geo_point_t test_p;
                for (int ii = 0; ii != num_steps; ii++) {
                    test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                               p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                               p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));
                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa()!=-1){
                            geo_point_t test_p_raw = test_p;
                            if (needs_transform) {
                                test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            if (!good_point) num_bad++;
                        }
                    }

                    if (dg79_census && c79_fk < 0 &&
                        (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps))) {
                        c79_fk = ii;
                    }

                    // doc 79 round 1: exact early break.  num_bad is monotone
                    // nondecreasing and num_steps is fixed, so once the kill
                    // predicate below holds it holds at loop end and this pair
                    // is killed regardless of the remaining steps; the walk has
                    // no other side effect.  Census mode walks the full path so
                    // the savable counters stay complete (bit-identical output
                    // either way -- the kill decision is unchanged).
                    if (!dg79_census && (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps))) {
                        break;
                    }
                }

                if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                    index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                }

                if (dg79_census) {
                    c79_steps[2] += num_steps;
                    if (c79_fk >= 0) { c79_kills[2]++; c79_savable[2] += num_steps - 1 - c79_fk; }
                }
            }
        }
    }

    if (dg79_census) {
        dg79_log->debug("CTPC79CENSUS summary ncomp={} pairs={} hough_pairs={} "
                        "main(steps={} kills={} savable={}) dir1(steps={} kills={} savable={}) "
                        "dir2(steps={} kills={} savable={})",
                        num, c79_pairs, c79_hough,
                        c79_steps[0], c79_kills[0], c79_savable[0],
                        c79_steps[1], c79_kills[1], c79_savable[1],
                        c79_steps[2], c79_kills[2], c79_savable[2]);
    }

    // deal with MST of first type
    if (lazy_walk) {
        // doc 79 round 2 (same scheme as connect_graph_relaxed's doc 78
        // round 2): single-pass lazy Kruskal.  Consider pairs in ascending
        // (distance, j, k) order and evaluate the walk verdict ONLY when a
        // pair would bridge two union-find components.  Killed bridges are
        // skipped; the accepted set IS the spanning forest of the
        // walk-filtered graph (identical to the legacy per-component Prim
        // except possibly on exact distance ties, which is inside this
        // knob's validation).  The < 3 cm direct-edge rule below needs
        // every near pair's verdict regardless of the forest, and pairs
        // carrying a directional candidate feed MST #2's recording of the
        // MAIN entry (process_mst_deterministically copies
        // index_index_dis[e], and a killed main blocks the dir edges
        // downstream) -- both classes are walked eagerly.
        std::vector<std::vector<char>> walked(num, std::vector<char>(num, 0));
        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                if (std::get<0>(index_index_dis[j][k]) < 0) continue;
                if (std::get<2>(index_index_dis[j][k]) < 3 * units::cm ||
                    std::get<0>(index_index_dis_dir1[j][k]) >= 0 ||
                    std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    eval_main_verdict(j, k);
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
                eval_main_verdict(j, k);
                walked[j][k] = 1;
                ++n_walked;
            }
            if (std::get<0>(index_index_dis[j][k]) < 0) continue;  // killed bridge
            uf[rj] = rk;
            index_index_dis_mst[j][k] = index_index_dis[j][k];
            if (++n_accepted + 1 == num) break;  // spanning forest complete
        }
        SPDLOG_LOGGER_DEBUG(WireCell::Log::logger("clus"),
            "connect_graph_ctpc lazy: num={} pairs={} walked={} accepted={}",
            num, order.size(), n_walked, n_accepted);
    }
    else {
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

    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            if (std::get<2>(index_index_dis[j][k]) < 3 * units::cm) {
                index_index_dis_mst[j][k] = index_index_dis[j][k];
            }

            // establish the path ...
            if (std::get<0>(index_index_dis_mst[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_mst[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_mst[j][k]));
                float dis;
                if (std::get<2>(index_index_dis_mst[j][k]) > 5 * units::cm) {
                    dis = std::get<2>(index_index_dis_mst[j][k]);
                }
                else {
                    dis = std::get<2>(index_index_dis_mst[j][k]);
                }
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
                    if (!boost::edge(gind1, gind2, graph).second) {
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
                    if (!boost::edge(gind1, gind2, graph).second) {
                        /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                    }
                }
            }

        }  // k
    }  // j
}




using namespace WireCell::Clus::Facade;

void Graphs::connect_graph_ctpc_with_reference(
    const Facade::Cluster& cluster,
    const Facade::Cluster& ref_cluster,
    IDetectorVolumes::pointer dv,
    Clus::IPCTransformSet::pointer pcts,
    Weighted::Graph& graph)
{
    // Enable CTPC functionality (combining ctpc logic with reference filtering)
    const bool use_ctpc = true;
    const auto* grouping = cluster.grouping();
    
    // Drift direction for directional analysis (equivalent to prototype's drift_dir)
    geo_vector_t drift_dir_abs(1, 0, 0);
    
    // Form connected components from existing graph
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
    std::set<size_t> excluded_points;  // Track excluded points (prototype's excluded_points)
    
    // Check if reference cluster is valid and not empty
    bool use_reference_filtering = (ref_cluster.is_valid() && ref_cluster.npoints() > 0);

    // Hoist KD-tree reference and query_point allocation out of the per-point loop
    const auto* ref_kd_ptr = use_reference_filtering ? &ref_cluster.kd3d() : nullptr;
    std::vector<double> query_point(3);

    // REFERENCE FILTERING PHASE - equivalent to prototype's filtering logic
    for (size_t i = 0; i < component.size(); ++i) {
        bool should_exclude = false;
        
        // Phase 1: Check point quality (equivalent to prototype's mcell->IsPointGood)
        if (!is_point_good(cluster, i, 2)) {
            should_exclude = true;
        } 
        // Phase 2: Reference cluster distance filtering (only if ref_cluster is not empty)
        else if (use_reference_filtering) {
            double temp_min_dis = 0;
            query_point[0] = points[0][i];
            query_point[1] = points[1][i];
            query_point[2] = points[2][i];
            auto knn_result = ref_kd_ptr->knn(1, query_point);
            
            if (!knn_result.empty()) {
                temp_min_dis = std::sqrt(knn_result[0].second);  // knn returns squared distance
            }
            
            // Key filtering criterion from prototype: >= 1.0 cm means exclude
            if (temp_min_dis >= 1.0 * units::cm) {
                should_exclude = true;
            }
        }
        // If ref_cluster is empty, only use point quality check
        
        if (should_exclude) {
            excluded_points.insert(i);
        } else {
            // Add to appropriate component cloud
            size_t comp_idx = component[i];
            pt_clouds.at(comp_idx)->add({points[0][i], points[1][i], points[2][i]});
            pt_clouds_global_indices.at(comp_idx).push_back(i);
        }
    }
    
    // Store excluded points in cluster cache (matches prototype's excluded_points)
    const_cast<Cluster&>(cluster).set_excluded_points(excluded_points);



    // Initialize distance metric containers (same structure as baseline)
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

    // Initialize all distances to invalid/infinite
    for (size_t j = 0; j != num; j++) {
        for (size_t k = 0; k != num; k++) {
            index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_mst[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
            index_index_dis_dir_mst[j][k] = std::make_tuple(-1, -1, 1e9);
        }
    }

    // Hoist scope-transform and cluster_t0 out of all per-step CTPC loops
    const bool needs_transform = (cluster.get_default_scope().hash() != cluster.get_raw_scope().hash());
    const auto ctpc_transform = needs_transform ? pcts->pc_transform(cluster.get_scope_transform()) : nullptr;
    const double cluster_t0 = needs_transform ? cluster.get_cluster_t0() : 0.0;

    // DISTANCE CALCULATION AND CTPC PATH VALIDATION
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            // Skip pairs where reference filtering emptied a component
            if (pt_clouds[j]->get_num_points() == 0 || pt_clouds[k]->get_num_points() == 0) continue;

            // Find closest points between components
            index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));

            // Enhanced directional analysis for large components (from prototype logic)
            if ((num < 100 && pt_clouds.at(j)->get_num_points() > 100 && pt_clouds.at(k)->get_num_points() > 100 &&
                 (pt_clouds.at(j)->get_num_points() + pt_clouds.at(k)->get_num_points()) > 400) ||
                (pt_clouds.at(j)->get_num_points() > 500 && pt_clouds.at(k)->get_num_points() > 500)) {
                
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));
                
                // Use cluster's vhough_transform method for directional analysis
                geo_vector_t dir1 = cluster.vhough_transform(p1, 30 * units::cm, 
                    Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                geo_vector_t dir2 = cluster.vhough_transform(p2, 30 * units::cm, 
                    Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                dir1 = dir1 * -1;
                dir2 = dir2 * -1;

                // Directional search from p1 towards p2
                std::pair<int, double> result1 = pt_clouds.at(k)->get_closest_point_along_vec(
                    p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                
                // Enhanced drift direction analysis (from prototype)
                if (result1.first < 0) {
                    double angle_deg = dir1.angle(drift_dir_abs) * 180.0 / M_PI;
                    if (std::abs(angle_deg - 90.0) < 10.0) {
                        // Direction nearly perpendicular to drift - try longer hough transform
                        if (std::abs(angle_deg - 90.0) < 5.0) {
                            dir1 = cluster.vhough_transform(p1, 80 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                        } else if (std::abs(angle_deg - 90.0) < 10.0) {
                            dir1 = cluster.vhough_transform(p1, 50 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                        }
                        dir1 = dir1 * -1;
                        result1 = pt_clouds.at(k)->get_closest_point_along_vec(
                            p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                    }
                }

                if (result1.first >= 0) {
                    index_index_dis_dir1[j][k] = std::make_tuple(
                        std::get<0>(index_index_dis[j][k]), result1.first, result1.second);
                }

                // Directional search from p2 towards p1
                std::pair<int, double> result2 = pt_clouds.at(j)->get_closest_point_along_vec(
                    p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);

                // Enhanced drift direction analysis for dir2
                if (result2.first < 0) {
                    double angle_deg2 = dir2.angle(drift_dir_abs) * 180.0 / M_PI;
                    if (std::abs(angle_deg2 - 90.0) < 10.0) {
                        if (std::abs(angle_deg2 - 90.0) < 5.0) {
                            dir2 = cluster.vhough_transform(p2, 80 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                        } else if (std::abs(angle_deg2 - 90.0) < 10.0) {
                            dir2 = cluster.vhough_transform(p2, 50 * units::cm, 
                                Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                        }
                        dir2 = dir2 * -1;
                        result2 = pt_clouds.at(j)->get_closest_point_along_vec(
                            p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);
                    }
                }

                if (result2.first >= 0) {
                    index_index_dis_dir2[j][k] = std::make_tuple(
                        result2.first, std::get<1>(index_index_dis[j][k]), result2.second);
                }
            }

            // CTPC PATH VALIDATION - Check basic distance path
            {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis[j][k])));

                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis / step_dis + 1;
                int num_bad = 0;
                geo_point_t test_p;
                
                for (int ii = 0; ii != num_steps; ii++) {
                    test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                               p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                               p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));

                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa() != -1) {
                            geo_point_t test_p_raw = test_p;
                            if (needs_transform) {
                                test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            if (!good_point) num_bad++;
                        }
                    }
                }

                if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                    index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
                }
            }

            // CTPC PATH VALIDATION - Check directional path 1
            if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir1[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k])));

                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir1[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis / step_dis + 1;
                int num_bad = 0;
                geo_point_t test_p;
                
                for (int ii = 0; ii != num_steps; ii++) {
                    test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                               p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                               p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));
                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa() != -1) {
                            geo_point_t test_p_raw = test_p;
                            if (needs_transform) {
                                test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            if (!good_point) num_bad++;
                        }
                    }
                }

                if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                    index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                }
            }

            // CTPC PATH VALIDATION - Check directional path 2
            if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir2[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k])));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir2[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis / step_dis + 1;
                int num_bad = 0;
                geo_point_t test_p;
                
                for (int ii = 0; ii != num_steps; ii++) {
                    test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                               p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                               p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));
                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa() != -1) {
                            geo_point_t test_p_raw = test_p;
                            if (needs_transform) {
                                test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            if (!good_point) num_bad++;
                        }
                    }
                }

                if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                    index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                }
            }
        }
    }

    // Build MST for basic distances
    {
        Weighted::Graph temp_graph(num);
        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j, index2 = k;
                if (std::get<0>(index_index_dis[j][k]) >= 0) {
                    if (!boost::edge(index1, index2, temp_graph).second) {
                        add_edge(index1, index2, std::get<2>(index_index_dis[j][k]), temp_graph);
                    }
                }
            }
        }
        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_mst);
    }

    // Build MST for directional distances
    {
        Weighted::Graph temp_graph(num);
        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j, index2 = k;
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0 || 
                    std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    if (!boost::edge(index1, index2, temp_graph).second) {
                        add_edge(index1, index2,
                            std::min(std::get<2>(index_index_dis_dir1[j][k]),
                                   std::get<2>(index_index_dis_dir2[j][k])), temp_graph);
                    }
                }
            }
        }
        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_dir_mst);
    }

    // Final graph construction phase
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            // Add short distance connections directly to MST
            if (std::get<2>(index_index_dis[j][k]) < 3 * units::cm) {
                index_index_dis_mst[j][k] = index_index_dis[j][k];
            }

            // Add MST basic distance edges to graph
            if (std::get<0>(index_index_dis_mst[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_mst[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_mst[j][k]));

                float dis;
                if (std::get<2>(index_index_dis_mst[j][k]) > 5 * units::cm) {
                    dis = std::get<2>(index_index_dis_mst[j][k]);
                } else {
                    dis = std::get<2>(index_index_dis_mst[j][k]);
                }
                if (!boost::edge(gind1, gind2, graph).second) {
                    /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                }
            }

            // Add MST directional edges to graph (with penalty for longer distances)
            if (std::get<0>(index_index_dis_dir_mst[j][k]) >= 0) {
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k]));

                    float dis;
                    if (std::get<2>(index_index_dis_dir1[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir1[j][k]) * 1.1;  // Matches ctpc baseline penalty
                    } else {
                        dis = std::get<2>(index_index_dis_dir1[j][k]);
                    }
                    if (!boost::edge(gind1, gind2, graph).second) {
                        /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                    }
                }
                
                if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k]));

                    float dis;
                    if (std::get<2>(index_index_dis_dir2[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir2[j][k]) * 1.1;  // Matches ctpc baseline penalty
                    } else {
                        dis = std::get<2>(index_index_dis_dir2[j][k]);
                    }
                    if (!boost::edge(gind1, gind2, graph).second) {
                        /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                    }
                }
            }
        }
    }
}