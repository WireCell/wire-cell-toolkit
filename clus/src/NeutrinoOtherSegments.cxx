#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"
// #include "WireCellClus/Graphs/Weighted.h"

#include <algorithm>
#include <cmath>

using namespace WireCell::Clus::PR;
using namespace WireCell::Clus;

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

// Edge property tag for Boost Graph
struct edge_base_t {
    typedef boost::edge_property_tag kind;
};

// Helper struct to track segment candidates
struct Res_proto_segment {
    int group_num;
    int number_points;
    size_t special_A;
    size_t special_B;
    double length;
    int number_not_faked;
    double max_dis_u;
    double max_dis_v;
    double max_dis_w;
};

void PatternAlgorithms::find_other_segments(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, bool flag_break_track, double search_range, double scaling_2d)
{
    // Get steiner point cloud data
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
    const auto& wpid_array = steiner_pc.get("wpid")->elements<WirePlaneId>();
    
    const size_t N = x_coords.size();
    if (N == 0) return;
    
    // Step 1: Tag points near existing segments
    std::vector<bool> flag_tagged(N, false);
    // int num_tagged = 0;
    
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(cluster.get_scope_transform(cluster.get_default_scope()));
    double cluster_t0 = cluster.get_cluster_t0();
    
    // Get all existing segments in this cluster
    std::set<SegmentPtr> existing_segments = find_cluster_segments(graph, cluster);
    
    for (size_t i = 0; i < N; i++) {
        Facade::geo_point_t p(x_coords[i], y_coords[i], z_coords[i]);
        double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;
        double min_3d_dis = 1e9;
        
        WirePlaneId wpid = wpid_array[i];
        int apa = wpid.apa();
        int face = wpid.face();
        
        // Check distances to existing segments
        for (auto seg : existing_segments) {
            // Get closest 3D point
            auto closest_result = segment_get_closest_point(seg, p, "fit");
            double dis_3d = closest_result.first;  // distance is already computed
            
            if (dis_3d < min_3d_dis) min_3d_dis = dis_3d;
            
            if (dis_3d < search_range) {
                flag_tagged[i] = true;
                // num_tagged++;
                break;
            }
            
            // Get 2D distances
            auto closest_2d = segment_get_closest_2d_distances(seg, p, apa, face, "fit");
            double dis_u = std::get<0>(closest_2d);
            double dis_v = std::get<1>(closest_2d);
            double dis_w = std::get<2>(closest_2d);
            
            if (dis_u < min_dis_u) min_dis_u = dis_u;
            if (dis_v < min_dis_v) min_dis_v = dis_v;
            if (dis_w < min_dis_w) min_dis_w = dis_w;
        }
        


        // Additional tagging based on 2D projections and dead channels
        if (!flag_tagged[i]) {
            auto p_raw = transform->backward(p, cluster_t0, face, apa);
            
            bool u_ok = (min_dis_u < scaling_2d * search_range || 
                        cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0));
            bool v_ok = (min_dis_v < scaling_2d * search_range || 
                        cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1));
            bool w_ok = (min_dis_w < scaling_2d * search_range || 
                        cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2));
            
            if (u_ok && v_ok && w_ok) {
                flag_tagged[i] = true;
            }
        

        // SPDLOG_LOGGER_DEBUG(
        //     s_log,
        //     "find_other_segments:{}: index {}; point {}: flag_tagged={}, min_3d_dis={:.2f} cm, min_dis_u={:.2f} cm, min_dis_v={:.2f} cm, min_dis_w={:.2f} cm, search_range={:.2f} cm, scaling_2d={:.2f}, dead_chs_u={}, dead_chs_v={}, dead_chs_w={}",
        //     __func__, i, p, (bool)flag_tagged[i],
        //     min_3d_dis / units::cm,
        //     min_dis_u / units::cm,
        //     min_dis_v / units::cm,
        //     min_dis_w / units::cm,
        //     search_range / units::cm,
        //     scaling_2d, cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0), cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1), cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2));
        }
    }


    
    // Step 2: Get terminal vertices
    const auto& flag_steiner_terminal = steiner_pc.get("flag_steiner_terminal")->elements<int>();
    std::vector<size_t> terminals;
    std::map<size_t, size_t> map_oindex_tindex;
    
    for (size_t i = 0; i < flag_steiner_terminal.size(); i++) {
        if (flag_steiner_terminal[i]) {
            map_oindex_tindex[i] = terminals.size();
            terminals.push_back(i);
        }
    }
    
    if (terminals.empty()) return;

    SPDLOG_LOGGER_DEBUG(s_log, "find_other_segments: existing_segments={}, N={}, num_tagged={}, terminals={}", existing_segments.size(), N, std::count(flag_tagged.begin(), flag_tagged.end(), true), terminals.size());


    // Step 3: Compute Voronoi diagram
    const auto& steiner_graph = cluster.get_graph("steiner_graph");
    using namespace Graphs::Weighted;
    auto vor = voronoi(steiner_graph, terminals);
    
    // Step 4: Build terminal graph with MST
    using Base = boost::property<edge_base_t, edge_type>;
    using WeightProperty = boost::property<boost::edge_weight_t, double, Base>;
    using TerminalGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                                 boost::no_property, WeightProperty>;
    
    TerminalGraph terminal_graph(N);
    std::map<std::pair<size_t, size_t>, std::pair<double, edge_type>> map_saved_edge;
    
    auto edge_weight = get(boost::edge_weight, steiner_graph);
    
    for (auto w : boost::make_iterator_range(edges(steiner_graph))) {
        size_t nearest_to_source = vor.terminal[source(w, steiner_graph)];
        size_t nearest_to_target = vor.terminal[target(w, steiner_graph)];
        
        if (nearest_to_source != nearest_to_target) {
            double weight = vor.distance[source(w, steiner_graph)] + 
                          vor.distance[target(w, steiner_graph)] + 
                          edge_weight[w];
            
            auto edge_pair1 = std::make_pair(nearest_to_source, nearest_to_target);
            auto edge_pair2 = std::make_pair(nearest_to_target, nearest_to_source);
            
            auto it1 = map_saved_edge.find(edge_pair1);
            auto it2 = map_saved_edge.find(edge_pair2);
            
            if (it1 != map_saved_edge.end()) {
                if (weight < it1->second.first) {
                    it1->second = std::make_pair(weight, w);
                }
            } else if (it2 != map_saved_edge.end()) {
                if (weight < it2->second.first) {
                    it2->second = std::make_pair(weight, w);
                }
            } else {
                map_saved_edge[edge_pair1] = std::make_pair(weight, w);
            }
        }
    }
    
    // Add edges to terminal graph
    for (const auto& [edge_pair, weight_info] : map_saved_edge) {
        boost::add_edge(edge_pair.first, edge_pair.second,
                       WeightProperty(weight_info.first, Base(weight_info.second)),
                       terminal_graph);
    }
    
    // Step 5: Find minimum spanning tree
    std::vector<boost::graph_traits<TerminalGraph>::edge_descriptor> mst_edges;
    boost::kruskal_minimum_spanning_tree(terminal_graph, std::back_inserter(mst_edges));
    
    // Step 6: Create cluster graph based on tagging
    TerminalGraph terminal_graph_cluster(terminals.size());
    std::map<size_t, std::set<size_t>> map_connection;
    
    for (const auto& edge : mst_edges) {
        size_t source_idx = boost::source(edge, terminal_graph);
        size_t target_idx = boost::target(edge, terminal_graph);
        
        if (flag_tagged[source_idx] == flag_tagged[target_idx]) {
            boost::add_edge(map_oindex_tindex[source_idx], 
                          map_oindex_tindex[target_idx], 
                          terminal_graph_cluster);
        } else {
            if (map_connection.find(source_idx) == map_connection.end()) {
                std::set<size_t> temp_results;
                temp_results.insert(target_idx);
                map_connection[source_idx] = temp_results;
            } else {
                map_connection[source_idx].insert(target_idx);
            }
            
            if (map_connection.find(target_idx) == map_connection.end()) {
                std::set<size_t> temp_results;
                temp_results.insert(source_idx);
                map_connection[target_idx] = temp_results;
            } else {
                map_connection[target_idx].insert(source_idx);
            }
        }
    }
    
    // Step 7: Find connected components
    std::vector<int> component(boost::num_vertices(terminal_graph_cluster));
    const int num_components = boost::connected_components(terminal_graph_cluster, &component[0]);
    
    std::vector<int> ncounts(num_components, 0);
    std::vector<std::vector<size_t>> sep_clusters(num_components);
    
    for (size_t i = 0; i < component.size(); ++i) {
        ncounts[component[i]]++;
        sep_clusters[component[i]].push_back(terminals[i]);
    }

    SPDLOG_LOGGER_DEBUG(s_log, "find_other_segments: num_components={}",  num_components);

    SPDLOG_LOGGER_DEBUG(s_log, "find_other_segments: Start tracks --- # of Vertices: {}; # of Edges: {}", boost::num_vertices(graph), boost::num_edges(graph));

    
    // Step 8: Analyze each cluster and filter
    std::vector<Res_proto_segment> temp_segments(num_components);
    std::set<int> remaining_segments;
    
    for (int i = 0; i < num_components; i++) {
        // Skip if inside original track
        if (flag_tagged[sep_clusters[i].front()]) {
            continue;
        }
        
        remaining_segments.insert(i);
        temp_segments[i].group_num = i;
        temp_segments[i].number_points = ncounts[i];
        
        // Find connection point (special_A)
        size_t special_A = SIZE_MAX;
        std::vector<size_t> candidates_special_A;
        
        for (int j = 0; j < ncounts[i]; j++) {
            if (map_connection.find(sep_clusters[i][j]) != map_connection.end()) {
                candidates_special_A.push_back(sep_clusters[i][j]);
            }
        }
        
        if (!candidates_special_A.empty()) {
            if (candidates_special_A.size() > 1 && ncounts[i] > 6) {
                // Use PCA to find the best connection point
                std::vector<Facade::geo_point_t> tmp_points;
                for (int j = 0; j < ncounts[i]; j++) {
                    Facade::geo_point_t tmp_p(x_coords[sep_clusters[i][j]], 
                                             y_coords[sep_clusters[i][j]], 
                                             z_coords[sep_clusters[i][j]]);
                    tmp_points.push_back(tmp_p);
                }
                
                auto results_pca = calc_PCA_main_axis(tmp_points);
                double max_val = 0;
                
                for (size_t j = 0; j < candidates_special_A.size(); j++) {
                    size_t idx = candidates_special_A[j];
                    double val = std::abs(
                        (x_coords[idx] - results_pca.first.x()) * results_pca.second.x() +
                        (y_coords[idx] - results_pca.first.y()) * results_pca.second.y() +
                        (z_coords[idx] - results_pca.first.z()) * results_pca.second.z()
                    );
                    
                    if (val > max_val) {
                        max_val = val;
                        special_A = idx;
                    }
                }
            } else {
                special_A = candidates_special_A.front();
            }
        }
        
        // If no boundary connection was found, skip this component entirely
        if (special_A == SIZE_MAX) {
            remaining_segments.erase(i);
            continue;
        }

        // Find furthest point (special_B)
        size_t special_B = special_A;
        double min_dis = 0;
        int number_not_faked = 0;
        double max_dis_u = 0, max_dis_v = 0, max_dis_w = 0;

        for (int j = 0; j < ncounts[i]; j++) {
            size_t idx = sep_clusters[i][j];
            double dis = std::sqrt(
                std::pow(x_coords[idx] - x_coords[special_A], 2) +
                std::pow(y_coords[idx] - y_coords[special_A], 2) +
                std::pow(z_coords[idx] - z_coords[special_A], 2));
            
            if (dis > min_dis) {
                min_dis = dis;
                special_B = idx;
            }
            
            // Check if point is fake (too close to existing segments)
            Facade::geo_point_t p(x_coords[idx], y_coords[idx], z_coords[idx]);
            double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;
            
            WirePlaneId wpid = wpid_array[idx];
            int apa = wpid.apa();
            int face = wpid.face();
            
            for (auto seg : existing_segments) {
                auto closest_2d = segment_get_closest_2d_distances(seg, p, apa, face, "fit");
                double dis_u = std::get<0>(closest_2d);
                double dis_v = std::get<1>(closest_2d);
                double dis_w = std::get<2>(closest_2d);
                
                if (dis_u < min_dis_u) min_dis_u = dis_u;
                if (dis_v < min_dis_v) min_dis_v = dis_v;
                if (dis_w < min_dis_w) min_dis_w = dis_w;
            }
            
            auto p_raw = transform->backward(p, cluster_t0, face, apa);
            
            int flag_num = 0;
            if (min_dis_u > scaling_2d * search_range && 
                !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0)) flag_num++;
            if (min_dis_v > scaling_2d * search_range && 
                !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1)) flag_num++;
            if (min_dis_w > scaling_2d * search_range && 
                !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2)) flag_num++;
            
            if (min_dis_u > max_dis_u && 
                !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0)) max_dis_u = min_dis_u;
            if (min_dis_v > max_dis_v && 
                !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1)) max_dis_v = min_dis_v;
            if (min_dis_w > max_dis_w && 
                !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2)) max_dis_w = min_dis_w;
            
            if (flag_num >= 2) number_not_faked++;
        }
        
        double length = std::sqrt(
            std::pow(x_coords[special_A] - x_coords[special_B], 2) +
            std::pow(y_coords[special_A] - y_coords[special_B], 2) +
            std::pow(z_coords[special_A] - z_coords[special_B], 2));
        
        // Adjust special_A if length is too short
        if (length < 3 * units::cm && special_A != SIZE_MAX) {
            size_t save_index = special_A;
            double save_dis = 1e9;
            for (auto it1 = map_connection[special_A].begin(); 
                 it1 != map_connection[special_A].end(); it1++) {
                double temp_dis = std::sqrt(
                    std::pow(x_coords[special_A] - x_coords[*it1], 2) +
                    std::pow(y_coords[special_A] - y_coords[*it1], 2) +
                    std::pow(z_coords[special_A] - z_coords[*it1], 2));
                if (temp_dis < save_dis) {
                    save_index = *it1;
                    save_dis = temp_dis;
                }
            }
            special_A = save_index;
            length = std::sqrt(
                std::pow(x_coords[special_A] - x_coords[special_B], 2) +
                std::pow(y_coords[special_A] - y_coords[special_B], 2) +
                std::pow(z_coords[special_A] - z_coords[special_B], 2));
        }
        
        temp_segments[i].special_A = special_A;
        temp_segments[i].special_B = special_B;
        temp_segments[i].length = length;
        temp_segments[i].number_not_faked = number_not_faked;
        temp_segments[i].max_dis_u = max_dis_u;
        temp_segments[i].max_dis_v = max_dis_v;
        temp_segments[i].max_dis_w = max_dis_w;
        
        // Apply quality cuts
        if ((temp_segments[i].number_points == 1) ||
            (number_not_faked == 0 &&
             ((length < 3.5 * units::cm) ||
              (((number_not_faked < 0.25 * temp_segments[i].number_points) ||
               (number_not_faked < 0.4 * temp_segments[i].number_points && length < 7 * units::cm)) &&
              max_dis_u / units::cm < 3 && max_dis_v / units::cm < 3 && max_dis_w / units::cm < 3 &&
              max_dis_u + max_dis_v + max_dis_w < 6 * units::cm)))) {
            remaining_segments.erase(i);
        }
    }
    
    // Step 9: Process remaining segments in order of quality (one-pass architecture)
    std::vector<SegmentPtr> new_segments;    // high dQdx or long curvy segments → break
    std::vector<SegmentPtr> new_segments_1;  // lower quality → break with 2cm cut
    const double mip_dQdx = 43e3 / units::cm;

    while (!remaining_segments.empty()) {
        // Find the best segment (most non-faked points, then longest)
        double max_number_not_faked = 0;
        double max_length = 0;
        int max_length_cluster = -1;

        for (auto it = remaining_segments.begin(); it != remaining_segments.end(); it++) {
            if (temp_segments[*it].number_not_faked > max_number_not_faked) {
                max_length_cluster = *it;
                max_number_not_faked = temp_segments[*it].number_not_faked;
                max_length = temp_segments[*it].length;
            } else if (temp_segments[*it].number_not_faked == max_number_not_faked) {
                if (temp_segments[*it].length > max_length) {
                    max_length_cluster = *it;
                    max_number_not_faked = temp_segments[*it].number_not_faked;
                    max_length = temp_segments[*it].length;
                }
            }
        }

        if (max_length_cluster == -1) break;

        remaining_segments.erase(max_length_cluster);

        size_t special_A = temp_segments[max_length_cluster].special_A;
        size_t special_B = temp_segments[max_length_cluster].special_B;

        if (special_A == SIZE_MAX || special_B == SIZE_MAX) continue;

        // Get path via Steiner graph (do_rough_path)
        Facade::geo_point_t pt_A(x_coords[special_A], y_coords[special_A], z_coords[special_A]);
        Facade::geo_point_t pt_B(x_coords[special_B], y_coords[special_B], z_coords[special_B]);
        auto path_points = do_rough_path(cluster, pt_A, pt_B);

        if (path_points.size() <= 1) continue;

        // Create segment (not yet in graph)
        auto new_seg = create_segment_for_cluster(cluster, dv, path_points);
        if (!new_seg) continue;

        // Do single tracking to get fine fit path and "fit" point cloud
        track_fitter.add_segment(new_seg);
        track_fitter.do_single_tracking(new_seg, true, true);

        // Always add to existing_segments for 2D distance re-evaluation of remaining clusters
        existing_segments.insert(new_seg);

        if (new_seg->fits().size() > 1) {
            // Find existing vertices/segments near the endpoints
            VertexPtr v1 = find_vertex_other_segment(graph, cluster, new_seg, true,  pt_A, track_fitter, dv);
            VertexPtr v2 = find_vertex_other_segment(graph, cluster, new_seg, false, pt_B, track_fitter, dv);

            bool v1_existed = (v1 != nullptr);
            bool v2_existed = (v2 != nullptr);

            // Create new (isolated) vertices where no existing one was found
            if (!v1) {
                v1 = make_vertex(graph);
                v1->wcpt().point = path_points.front();
                v1->cluster(&cluster);
            }
            if (!v2) {
                v2 = make_vertex(graph);
                v2->wcpt().point = path_points.back();
                v2->cluster(&cluster);
            }

            // Handle corner case: both find_vertex calls returned the same vertex
            if (v1 == v2) {
                double dis1 = std::sqrt(
                    std::pow(v1->wcpt().point.x() - pt_A.x(), 2) +
                    std::pow(v1->wcpt().point.y() - pt_A.y(), 2) +
                    std::pow(v1->wcpt().point.z() - pt_A.z(), 2));
                double dis2 = std::sqrt(
                    std::pow(v1->wcpt().point.x() - pt_B.x(), 2) +
                    std::pow(v1->wcpt().point.y() - pt_B.y(), 2) +
                    std::pow(v1->wcpt().point.z() - pt_B.z(), 2));
                if (dis1 < dis2) {
                    // v1 is close to pt_A → keep v1, replace v2 with a new vertex at pt_B
                    v2 = make_vertex(graph);
                    v2->wcpt().point = pt_B;
                    v2->cluster(&cluster);
                    v2_existed = false;
                } else {
                    // v2 is close to pt_B → keep v2, replace v1 with a new vertex at pt_A
                    v1 = make_vertex(graph);
                    v1->wcpt().point = pt_A;
                    v1->cluster(&cluster);
                    v1_existed = false;
                }
            }

            // Path extension: if the found vertex positions differ from the rough path endpoints,
            // re-route via do_rough_path so that new_seg's wcpts (and "main" DPC) span the actual
            // vertex positions.  do_multi_tracking uses wcpts as its starting trajectory, so
            // without this the fit would miss the true endpoint by up to ~vtx_cut2 (~2 cm).
            {
                double dist_v1 = std::sqrt(
                    std::pow(v1->wcpt().point.x() - path_points.front().x(), 2) +
                    std::pow(v1->wcpt().point.y() - path_points.front().y(), 2) +
                    std::pow(v1->wcpt().point.z() - path_points.front().z(), 2));
                double dist_v2 = std::sqrt(
                    std::pow(v2->wcpt().point.x() - path_points.back().x(), 2) +
                    std::pow(v2->wcpt().point.y() - path_points.back().y(), 2) +
                    std::pow(v2->wcpt().point.z() - path_points.back().z(), 2));

                if (dist_v1 > 0.1 * units::cm || dist_v2 > 0.1 * units::cm) {
                    auto ext_path = do_rough_path(cluster, v1->wcpt().point, v2->wcpt().point);
                    if (ext_path.size() > 1) {
                        std::vector<PR::WCPoint> new_wcpts;
                        new_wcpts.reserve(ext_path.size());
                        for (const auto& p : ext_path) {
                            PR::WCPoint wcp;
                            wcp.point = p;
                            new_wcpts.push_back(wcp);
                        }
                        new_seg->wcpts(new_wcpts);
                        create_segment_point_cloud(new_seg, ext_path, dv, "main");
                    } else {
                        std::cerr << "Warning: Find_other_segments: path extension failed for cluster " << cluster.get_cluster_id() << std::endl;
                    }
                }
            }

            Facade::geo_point_t v1_fit_pt = v1->fit().valid() ? v1->fit().point : v1->wcpt().point;
            Facade::geo_point_t v2_fit_pt = v2->fit().valid() ? v2->fit().point : v2->wcpt().point;
            double v1v2_dist = std::sqrt(
                std::pow(v1_fit_pt.x() - v2_fit_pt.x(), 2) +
                std::pow(v1_fit_pt.y() - v2_fit_pt.y(), 2) +
                std::pow(v1_fit_pt.z() - v2_fit_pt.z(), 2));

            if (v1_existed || v2_existed) {
                // At least one endpoint connects to the existing graph
                if (v1v2_dist > 0.1 * units::cm) {
                    add_segment(graph, new_seg, v1, v2);
                    track_fitter.do_multi_tracking(true, true, false);

                    double length        = segment_track_length(new_seg);
                    double direct_length = segment_track_direct_length(new_seg);
                    double medium_dQ_dx  = segment_median_dQ_dx(new_seg);

                    if (length > 30 * units::cm ||
                        (direct_length < 0.78 * length && length > 10 * units::cm &&
                         medium_dQ_dx / mip_dQdx > 1.6)) {
                        new_segments.push_back(new_seg);
                    }
                    SPDLOG_LOGGER_DEBUG(s_log, "find_other_segments: Cluster {} Other tracks --- # of Vertices: {}; # of Edges: {}", cluster.get_cluster_id(), boost::num_vertices(graph), boost::num_edges(graph));

                } else {
                    // Endpoints at same position – discard fresh vertices
                    if (!v1_existed) remove_vertex(graph, v1);
                    if (!v2_existed && v2 != v1) remove_vertex(graph, v2);
                }
            } else {
                // Isolated segment: try isochronous parallel-track correction.
                // modify_vertex/segment_isochronous will call add_segment internally if successful,
                // so do NOT call add_segment here.
                bool flag_parallel = false;

                SPDLOG_LOGGER_DEBUG(s_log, "find_other_segments: Cluster {} Middle tracks --- # of Vertices: {}; # of Edges: {}", cluster.get_cluster_id(), boost::num_vertices(graph), boost::num_edges(graph));

                Facade::geo_vector_t dir(v1_fit_pt.x() - v2_fit_pt.x(),
                                        v1_fit_pt.y() - v2_fit_pt.y(),
                                        v1_fit_pt.z() - v2_fit_pt.z());
                double dir_mag = dir.magnitude();

                if (dir_mag > 10 * units::cm ||
                    (dir_mag > 8 * units::cm && segment_track_length(new_seg) > 13 * units::cm)) {

                    // Try to snap to a nearby isochronous vertex
                    const WireCell::Vector drift_dir(1, 0, 0);
                    auto [vbegin, vend] = boost::vertices(graph);
                    for (auto vit = vbegin; vit != vend && !flag_parallel; ++vit) {
                        VertexPtr vtx = graph[*vit].vertex;
                        if (!vtx || vtx->cluster() != &cluster) continue;
                        if (vtx == v1 || vtx == v2) continue;

                        Facade::geo_point_t vtx_fit = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;
                        WireCell::Vector d1(vtx_fit.x() - v1_fit_pt.x(),
                                            vtx_fit.y() - v1_fit_pt.y(),
                                            vtx_fit.z() - v1_fit_pt.z());
                        WireCell::Vector d2(vtx_fit.x() - v2_fit_pt.x(),
                                            vtx_fit.y() - v2_fit_pt.y(),
                                            vtx_fit.z() - v2_fit_pt.z());

                        if (d1.magnitude() < 6 * units::cm &&
                            std::fabs(drift_dir.angle(d1) / 3.1415926 * 180.0 - 90.0) < 15.0) {
                            flag_parallel = modify_vertex_isochronous(graph, cluster, vtx, v1, new_seg, v2, track_fitter, dv);
                        }
                        if (!flag_parallel && d2.magnitude() < 6 * units::cm &&
                            std::fabs(drift_dir.angle(d2) / 3.1415926 * 180.0 - 90.0) < 15.0) {
                            flag_parallel = modify_vertex_isochronous(graph, cluster, vtx, v2, new_seg, v1, track_fitter, dv);
                        }
                    }

                    // Try to snap to a nearby isochronous segment
                    if (!flag_parallel) {
                        double min_dis1 = 1e9; SegmentPtr min_sg1 = nullptr;
                        double min_dis2 = 1e9; SegmentPtr min_sg2 = nullptr;

                        auto [ebegin, eend] = boost::edges(graph);
                        for (auto eit = ebegin; eit != eend; ++eit) {
                            SegmentPtr sg1 = graph[*eit].segment;
                            if (!sg1 || sg1->cluster() != &cluster) continue;

                            double dis1 = segment_get_closest_point(sg1, v1_fit_pt, "fit").first;
                            double dis2 = segment_get_closest_point(sg1, v2_fit_pt, "fit").first;

                            if (!flag_parallel && dis1 < 6 * units::cm) {
                                flag_parallel = modify_segment_isochronous(graph, cluster, sg1, v1, new_seg, v2, track_fitter, dv);
                            }
                            if (!flag_parallel && dis2 < 6 * units::cm) {
                                flag_parallel = modify_segment_isochronous(graph, cluster, sg1, v2, new_seg, v1, track_fitter, dv);
                            }
                            if (flag_parallel) break;

                            if (dis1 < min_dis1) { min_dis1 = dis1; min_sg1 = sg1; }
                            if (dis2 < min_dis2) { min_dis2 = dis2; min_sg2 = sg1; }
                        }

                        // Long track (>18 cm): widen search
                        if (!flag_parallel && dir_mag > 18 * units::cm) {
                            if (min_sg1 && min_dis1 <= min_dis2 && min_dis1 < 10 * units::cm) {
                                flag_parallel = modify_segment_isochronous(graph, cluster, min_sg1, v1, new_seg, v2, track_fitter, dv, 10*units::cm, 8, 15*units::cm);
                                if (!flag_parallel) flag_parallel = modify_segment_isochronous(graph, cluster, min_sg1, v1, new_seg, v2, track_fitter, dv, 10*units::cm, 8, 30*units::cm);
                            } else if (min_sg2 && min_dis2 < min_dis1 && min_dis2 < 10 * units::cm) {
                                flag_parallel = modify_segment_isochronous(graph, cluster, min_sg2, v2, new_seg, v1, track_fitter, dv, 10*units::cm, 8, 15*units::cm);
                                if (!flag_parallel) flag_parallel = modify_segment_isochronous(graph, cluster, min_sg2, v2, new_seg, v1, track_fitter, dv, 10*units::cm, 8, 30*units::cm);
                            }
                        }

                        // Very long track (>36 cm): widen further
                        if (!flag_parallel && dir_mag > 36 * units::cm) {
                            if (min_sg1 && min_dis1 <= min_dis2 && min_dis1 < 18 * units::cm) {
                                flag_parallel = modify_segment_isochronous(graph, cluster, min_sg1, v1, new_seg, v2, track_fitter, dv, 18*units::cm, 5, 15*units::cm);
                                if (!flag_parallel) flag_parallel = modify_segment_isochronous(graph, cluster, min_sg1, v1, new_seg, v2, track_fitter, dv, 18*units::cm, 5, 30*units::cm);
                            } else if (min_sg2 && min_dis2 < min_dis1 && min_dis2 < 18 * units::cm) {
                                flag_parallel = modify_segment_isochronous(graph, cluster, min_sg2, v2, new_seg, v1, track_fitter, dv, 18*units::cm, 5, 15*units::cm);
                                if (!flag_parallel) flag_parallel = modify_segment_isochronous(graph, cluster, min_sg2, v2, new_seg, v1, track_fitter, dv, 18*units::cm, 5, 30*units::cm);
                            }
                        }
                    }
                }

                if (!flag_parallel) {
                    // Truly isolated residual – remove vertices from graph; segment was never added
                    SPDLOG_LOGGER_DEBUG(s_log, "find_other_segments: Cluster {} Isolated residual segment   # of Vertices: {}; # of Edges: {}", cluster.get_cluster_id(), boost::num_vertices(graph), boost::num_edges(graph));

                    remove_vertex(graph, v1);
                    remove_vertex(graph, v2);
                } else {
                    // Isochronous connection found (segment already added by modify_*)
                    track_fitter.do_multi_tracking(true, true, false);

                    double direct_length = segment_track_direct_length(new_seg);
                    double length        = segment_track_length(new_seg);
                    double medium_dQ_dx  = segment_median_dQ_dx(new_seg);

                    if ((direct_length < 0.78 * length && length > 10 * units::cm &&
                         medium_dQ_dx / mip_dQdx > 1.6) ||
                        (direct_length < 0.6 * length && length > 10 * units::cm)) {
                        if (medium_dQ_dx / mip_dQdx > 1.1) {
                            new_segments.push_back(new_seg);
                        } else {
                            new_segments_1.push_back(new_seg);
                        }
                        // SPDLOG_LOGGER_DEBUG(s_log, "find_other_segments: Cluster {} Other tracks -- isochronous connection, length={} cm", cluster.get_cluster_id(), length / units::cm);
                    }
                }
            }
        } // if fits().size() > 1

        // Re-evaluate remaining segments using the updated existing_segments set
        std::set<int> tmp_del_set;
        for (auto it = remaining_segments.begin(); it != remaining_segments.end(); it++) {
            temp_segments[*it].number_not_faked = 0;
            temp_segments[*it].max_dis_u = 0;
            temp_segments[*it].max_dis_v = 0;
            temp_segments[*it].max_dis_w = 0;

            for (int j = 0; j < ncounts[*it]; j++) {
                size_t idx = sep_clusters[*it][j];
                Facade::geo_point_t p(x_coords[idx], y_coords[idx], z_coords[idx]);
                double min_dis_u = 1e9, min_dis_v = 1e9, min_dis_w = 1e9;

                WirePlaneId wpid = wpid_array[idx];
                int apa = wpid.apa();
                int face = wpid.face();

                for (auto seg : existing_segments) {
                    auto closest_2d = segment_get_closest_2d_distances(seg, p, apa, face, "fit");
                    double dis_u = std::get<0>(closest_2d);
                    double dis_v = std::get<1>(closest_2d);
                    double dis_w = std::get<2>(closest_2d);

                    if (dis_u < min_dis_u) min_dis_u = dis_u;
                    if (dis_v < min_dis_v) min_dis_v = dis_v;
                    if (dis_w < min_dis_w) min_dis_w = dis_w;
                }

                auto p_raw = transform->backward(p, cluster_t0, face, apa);

                int flag_num = 0;
                if (min_dis_u > scaling_2d * search_range &&
                    !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0)) flag_num++;
                if (min_dis_v > scaling_2d * search_range &&
                    !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1)) flag_num++;
                if (min_dis_w > scaling_2d * search_range &&
                    !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2)) flag_num++;

                if (flag_num >= 2) temp_segments[*it].number_not_faked++;

                if (min_dis_u > temp_segments[*it].max_dis_u &&
                    !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 0))
                    temp_segments[*it].max_dis_u = min_dis_u;
                if (min_dis_v > temp_segments[*it].max_dis_v &&
                    !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 1))
                    temp_segments[*it].max_dis_v = min_dis_v;
                if (min_dis_w > temp_segments[*it].max_dis_w &&
                    !cluster.grouping()->get_closest_dead_chs(p_raw, 1, apa, face, 2))
                    temp_segments[*it].max_dis_w = min_dis_w;
            }

            // Apply quality cuts again
            if ((temp_segments[*it].number_points == 1) ||
                (temp_segments[*it].number_not_faked == 0 &&
                 ((temp_segments[*it].length < 3.5 * units::cm) || (
                  ((temp_segments[*it].number_not_faked < 0.25 * temp_segments[*it].number_points) ||
                   (temp_segments[*it].number_not_faked < 0.4 * temp_segments[*it].number_points &&
                    temp_segments[*it].length < 7 * units::cm)) &&
                  temp_segments[*it].max_dis_u / units::cm < 3 &&
                  temp_segments[*it].max_dis_v / units::cm < 3 &&
                  temp_segments[*it].max_dis_w / units::cm < 3 &&
                  temp_segments[*it].max_dis_u + temp_segments[*it].max_dis_v +
                  temp_segments[*it].max_dis_w < 6 * units::cm)))) {
                tmp_del_set.insert(*it);
            }
        }

        for (auto it = tmp_del_set.begin(); it != tmp_del_set.end(); it++) {
            remaining_segments.erase(*it);
        }
    }

    // Step 10: Break curvy/high-dQdx segments if requested
    if (flag_break_track) {
        break_segments(graph, track_fitter, dv, new_segments_1, 2.0 * units::cm);
        break_segments(graph, track_fitter, dv, new_segments);
    }
    SPDLOG_LOGGER_DEBUG(s_log, "find_other_segments: Cluster {} End tracks --- # of Vertices: {}; # of Edges: {}", cluster.get_cluster_id(), boost::num_vertices(graph), boost::num_edges(graph));

}


namespace {

using WireCell::Ray;

Facade::geo_vector_t make_vector(const Facade::geo_point_t& from, const Facade::geo_point_t& to)
{
    return Facade::geo_vector_t(to.x() - from.x(), to.y() - from.y(), to.z() - from.z());
}

double point_distance(const Facade::geo_point_t& p1, const Facade::geo_point_t& p2)
{
    return ray_length(Ray{p1, p2});
}

double line_distance(const Facade::geo_point_t& line_point,
                     const Facade::geo_vector_t& line_dir_unit,
                     const Facade::geo_point_t& test_point)
{
    return ray_closest_dis(Ray{line_point, line_point + line_dir_unit}, test_point);
}

double angle_degrees(const Facade::geo_vector_t& dir1, const Facade::geo_vector_t& dir2)
{
    return dir1.angle(dir2) / 3.1415926 * 180.0;
}

Facade::geo_vector_t tracking_direction(const std::vector<Facade::geo_point_t>& tracking_path,
                                        const Facade::geo_point_t& test_p,
                                        bool flag_front,
                                        int index)
{
    Facade::geo_vector_t dir_p(0, 0, 0);
    if (flag_front) {
        if (index + 1 < tracking_path.size()) {
            dir_p.x(test_p.x() - tracking_path.at(index + 1).x());
            dir_p.y(test_p.y() - tracking_path.at(index + 1).y());
            dir_p.z(test_p.z() - tracking_path.at(index + 1).z());
        }
        else {
            dir_p.x(test_p.x() - tracking_path.back().x());
            dir_p.y(test_p.y() - tracking_path.back().y());
            dir_p.z(test_p.z() - tracking_path.back().z());
        }
    }
    else {
        if (tracking_path.size() >= index + 2) {
            dir_p.x(test_p.x() - tracking_path.at(tracking_path.size() - 2 - index).x());
            dir_p.y(test_p.y() - tracking_path.at(tracking_path.size() - 2 - index).y());
            dir_p.z(test_p.z() - tracking_path.at(tracking_path.size() - 2 - index).z());
        }
        else {
            dir_p.x(test_p.x() - tracking_path.front().x());
            dir_p.y(test_p.y() - tracking_path.front().y());
            dir_p.z(test_p.z() - tracking_path.front().z());
        }
    }
    return dir_p;
}

}


std::tuple<VertexPtr, SegmentPtr, Facade::geo_point_t>
PatternAlgorithms::check_end_point(Graph& graph,
                                   Facade::Cluster& cluster,
                                   std::vector<Facade::geo_point_t>& tracking_path,
                                   bool flag_front,
                                   double vtx_cut1,
                                   double vtx_cut2,
                                   double sg_cut1,
                                   double sg_cut2)
{
    if (tracking_path.size() < 2) SPDLOG_LOGGER_DEBUG(s_log, "check_end_point: {}: vector size wrong!", __func__);

    const int ncount = 5;
    Facade::geo_point_t test_p = flag_front ? tracking_path.front() : tracking_path.back();
    VertexPtr vtx = nullptr;
    SegmentPtr seg = nullptr;

    for (int i = 0; i != ncount; i++) {
        Facade::geo_vector_t dir_p = tracking_direction(tracking_path, test_p, flag_front, i);
        if (dir_p.magnitude() == 0.0) {
            continue;
        }

        Facade::geo_vector_t temp_dir = dir_p.norm();

        auto [vbegin, vend] = boost::vertices(graph);
        for (auto vit = vbegin; vit != vend; ++vit) {
            VertexPtr test_v = graph[*vit].vertex;
            if (!test_v || test_v->cluster() != &cluster) {
                continue;
            }

            const Facade::geo_point_t p1 = test_v->wcpt().point;
            const Facade::geo_point_t p2 = test_v->fit().valid() ? test_v->fit().point : test_v->wcpt().point;

            const double dis1 = point_distance(test_p, p1);
            const double dis2 = line_distance(test_p, temp_dir, p1);
            const double dis3 = point_distance(test_p, p2);
            const double dis4 = line_distance(test_p, temp_dir, p2);

            const auto degree = test_v->descriptor_valid() ? boost::out_degree(test_v->get_descriptor(), graph) : 0;

            if (std::max(dis1, dis2) < 5 * units::cm &&
                ((std::min(dis1, dis2) < vtx_cut1 && std::max(dis1, dis2) < vtx_cut2) ||
                 (std::min(dis3, dis4) < vtx_cut1 * 1.3 && std::max(dis1, dis2) < vtx_cut2 * 2 && degree == 1))) {
                const Facade::geo_vector_t test_dir = make_vector(test_p, p1);
                if (angle_degrees(test_dir, temp_dir) < 90.0) {
                    vtx = test_v;
                    break;
                }
            }

            if (std::max(dis3, dis4) < 5 * units::cm &&
                ((std::min(dis3, dis4) < vtx_cut1 && std::max(dis3, dis4) < vtx_cut2) ||
                 (std::min(dis3, dis4) < vtx_cut1 * 1.3 && std::max(dis3, dis4) < vtx_cut2 * 2 && degree == 1))) {
                const Facade::geo_vector_t test_dir = make_vector(test_p, p2);
                if (angle_degrees(test_dir, temp_dir) < 90.0) {
                    vtx = test_v;
                    break;
                }
            }
        }

        if (vtx != nullptr) {
            break;
        }
    }

    if (vtx == nullptr) {
        SegmentPtr min_sg = nullptr;
        Facade::geo_point_t min_point(0, 0, 0);
        double min_dis = 1e9;

        auto [ebegin, eend] = boost::edges(graph);
        for (auto eit = ebegin; eit != eend; ++eit) {
            SegmentPtr candidate_seg = graph[*eit].segment;
            if (!candidate_seg || candidate_seg->cluster() != &cluster) {
                continue;
            }

            auto [closest_dis, closest_pt] = segment_get_closest_point(candidate_seg, test_p, "fit");
            if (closest_dis < sg_cut1) {
                const auto& points = candidate_seg->fits();

                for (int i = 0; i != ncount; i++) {
                    Facade::geo_vector_t dir_p = tracking_direction(tracking_path, test_p, flag_front, i);
                    if (dir_p.magnitude() == 0.0) {
                        continue;
                    }

                    Facade::geo_vector_t temp_dir = dir_p.norm();
                    for (size_t j = 0; j != points.size(); j++) {
                        double dis1 = line_distance(test_p, temp_dir, points.at(j).point);
                        if (dis1 < sg_cut2) {
                            const Facade::geo_vector_t test_dir = make_vector(test_p, points.at(j).point);
                            dis1 = std::sqrt(std::pow(dis1, 2) + std::pow(test_dir.magnitude() / 4.0, 2));
                            if (angle_degrees(test_dir, temp_dir) < 90.0 && dis1 < min_dis) {
                                min_dis = dis1;
                                min_point = points.at(j).point;
                                min_sg = candidate_seg;
                            }
                        }
                    }
                }
            }
        }

        if (min_sg != nullptr) {
            seg = min_sg;
            test_p = min_point;
        }
    }

    return std::make_tuple(vtx, seg, test_p);
}

VertexPtr PatternAlgorithms::find_vertex_other_segment(Graph& graph, Facade::Cluster& cluster, SegmentPtr seg, bool flag_forwrard, Facade::geo_point_t& wcp, TrackFitting& track_fitter, IDetectorVolumes::pointer dv ){
    VertexPtr v1 = nullptr;

    // Build the tracking path from the segment's fit point
    // (corresponds to temp_cluster->get_fine_tracking_path() in WCP)
    std::vector<Facade::geo_point_t> tracking_path;
    for (const auto& fit : seg->fits()) {
        tracking_path.push_back(fit.point);
    }

    if (tracking_path.size() < 2) return v1;

    // First attempt with default parameters
    auto check_results = check_end_point(graph, cluster, tracking_path, flag_forwrard);

    // Second attempt with widened parameters if nothing found
    if (std::get<0>(check_results) == nullptr && std::get<1>(check_results) == nullptr) {
        check_results = check_end_point(graph, cluster, tracking_path, flag_forwrard,
                                        1.2 * units::cm, 2.5 * units::cm);
    }

    // Third attempt with further widened parameters
    if (std::get<0>(check_results) == nullptr && std::get<1>(check_results) == nullptr) {
        check_results = check_end_point(graph, cluster, tracking_path, flag_forwrard,
                                        1.5 * units::cm, 3.0 * units::cm);
    }

    if (std::get<0>(check_results) != nullptr) {
        // Found a vertex directly
        v1 = std::get<0>(check_results);
    } else if (std::get<1>(check_results) != nullptr) {
        // Found a segment – need to break it or snap to an endpoint
        SegmentPtr sg1 = std::get<1>(check_results);
        Facade::geo_point_t test_p = std::get<2>(check_results);

        // Get the closest point on sg1 to test_p (using "main" point cloud)
        auto [dist_wcpt, break_wcp] = segment_get_closest_point(sg1, test_p, "main");

        // Find the two endpoint vertices of sg1
        auto [sv, ev] = find_vertices(graph, sg1);
        VertexPtr start_v = sv, end_v = ev;

        if (start_v && end_v) {
            const auto& wcpts = sg1->wcpts();
            if (!wcpts.empty()) {
                // Determine which vertex corresponds to the front of the segment
                // by comparing positions to the first wcpt (< 0.1 cm tolerance)
                double dis_sv_front = point_distance(start_v->wcpt().point, wcpts.front().point);
                double dis_sv_back  = point_distance(start_v->wcpt().point, wcpts.back().point);
                if (dis_sv_front > dis_sv_back) {
                    std::swap(start_v, end_v);
                }
            }
        }

        if (!start_v || !end_v) {
            SPDLOG_LOGGER_DEBUG(s_log, "find_vertex_other_segment: {}: Error in finding vertices for a segment", __func__);
            return v1;
        }

        // Check if the break point is close enough to snap to an existing endpoint
        double dis_start = point_distance(start_v->wcpt().point, break_wcp);
        double dis_end   = point_distance(end_v->wcpt().point, break_wcp);

        if (dis_start < 0.9 * units::cm) {
            v1 = start_v;
        } else if (dis_end < 0.9 * units::cm) {
            v1 = end_v;
        } else {
            // Break the segment into two at break_wcp
            std::list<Facade::geo_point_t> wcps_list1;
            std::list<Facade::geo_point_t> wcps_list2;
            proto_break_tracks(cluster, start_v->wcpt().point, break_wcp,
                               end_v->wcpt().point, wcps_list1, wcps_list2, true);

            // Create the new vertex at the break point
            VertexPtr new_vtx = make_vertex(graph);
            new_vtx->wcpt().point = break_wcp;
            new_vtx->cluster(&cluster);

            // Convert lists to vectors for segment creation
            std::vector<Facade::geo_point_t> path1(wcps_list1.begin(), wcps_list1.end());
            std::vector<Facade::geo_point_t> path2(wcps_list2.begin(), wcps_list2.end());

            SegmentPtr sg2 = create_segment_for_cluster(cluster, dv, path1, sg1->dirsign());
            SegmentPtr sg3 = create_segment_for_cluster(cluster, dv, path2, sg1->dirsign());

            if (sg2 && sg3) {
                add_segment(graph, sg2, start_v, new_vtx);
                add_segment(graph, sg3, new_vtx, end_v);
                remove_segment(graph, sg1);
                v1 = new_vtx;
            } else {
                // Segment creation failed – clean up the orphaned vertex
                remove_vertex(graph, new_vtx);
            }
        }
    }
    // Nothing found: return nullptr so the caller knows no existing connection exists.
    // The caller (find_other_segments) creates the fresh vertex in its if (!v1) block.

    return v1;
}

bool PatternAlgorithms::modify_vertex_isochronous(Graph& graph, Facade::Cluster& cluster, VertexPtr vtx, VertexPtr v1, SegmentPtr sg, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv){
    bool flag = false;

    // Get fit points for v1 and vtx (prefer fitted position, fall back to wcpt)
    Facade::geo_point_t v1_fit_pt  = v1->fit().valid()  ? v1->fit().point  : v1->wcpt().point;
    Facade::geo_point_t vtx_fit_pt = vtx->fit().valid() ? vtx->fit().point : vtx->wcpt().point;

    // Direction of sg at v1's position, reversed to point back toward vtx
    WireCell::Vector dir = segment_cal_dir_3vector(sg, v1_fit_pt, 15 * units::cm) * (-1.0);
    if (dir.x() == 0) return flag;

    // Project v1 along the track direction to the x-plane of vtx to find the new target point
    double dx = vtx_fit_pt.x() - v1_fit_pt.x();
    Facade::geo_point_t test_p(
        v1_fit_pt.x() + dx,
        v1_fit_pt.y() + dir.y() / dir.x() * dx,
        v1_fit_pt.z() + dir.z() / dir.x() * dx
    );

    // Find the closest steiner point to test_p
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();

    auto knn_results = cluster.kd_steiner_knn(1, test_p, "steiner_pc");
    if (knn_results.empty()) return flag;
    size_t vtx_new_idx = knn_results[0].first;
    Facade::geo_point_t vtx_new_pt(x_coords[vtx_new_idx], y_coords[vtx_new_idx], z_coords[vtx_new_idx]);

    // Reject if the new position is too far from v1
    if (point_distance(vtx_new_pt, v1_fit_pt) > 5 * units::cm) return flag;

    // Check connectivity between vtx (old position) and vtx_new along a straight line
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster.get_scope_transform(cluster.get_default_scope()));
    double cluster_t0 = cluster.get_cluster_t0();
    auto grouping = cluster.grouping();
    if (!transform || !grouping) return flag;

    double step_size = 0.6 * units::cm;
    int ncount = std::round(point_distance(vtx_fit_pt, vtx_new_pt) / step_size);
    int n_bad = 0;
    for (int i = 1; i < ncount; i++) {
        Facade::geo_point_t step_p(
            vtx_fit_pt.x() + (vtx_new_pt.x() - vtx_fit_pt.x()) / ncount * i,
            vtx_fit_pt.y() + (vtx_new_pt.y() - vtx_fit_pt.y()) / ncount * i,
            vtx_fit_pt.z() + (vtx_new_pt.z() - vtx_fit_pt.z()) / ncount * i
        );
        auto test_wpid = dv->contained_by(step_p);
        if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
            auto temp_p_raw = transform->backward(step_p, cluster_t0, test_wpid.face(), test_wpid.apa());
            if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2 * units::cm, 0, 0)) {
                n_bad++;
            }
        }
    }
    if (n_bad > 0) return flag;

    // Update paths for all segments already connected to vtx
    if (!vtx->descriptor_valid()) return flag;
    auto vd = vtx->get_descriptor();
    std::vector<SegmentPtr> vtx_segs;
    {
        auto edge_range = boost::out_edges(vd, graph);
        for (auto eit = edge_range.first; eit != edge_range.second; ++eit) {
            SegmentPtr sg1 = graph[*eit].segment;
            if (sg1) vtx_segs.push_back(sg1);
        }
    }
    for (auto sg1 : vtx_segs) {
        const auto& wcpts = sg1->wcpts();
        if (wcpts.empty()) continue;

        // Determine the "other" endpoint of sg1 (the one that is not vtx)
        double dis_front = point_distance(wcpts.front().point, vtx->wcpt().point);
        double dis_back  = point_distance(wcpts.back().point,  vtx->wcpt().point);
        Facade::geo_point_t other_pt = (dis_front < dis_back) ? wcpts.back().point : wcpts.front().point;

        // Recompute path from new vtx position to the other endpoint via the steiner graph
        auto new_path_pts = do_rough_path(cluster, vtx_new_pt, other_pt);
        std::vector<WCPoint> new_wcpts;
        new_wcpts.reserve(new_path_pts.size());
        for (const auto& p : new_path_pts) {
            WCPoint wcp;
            wcp.point = p;
            new_wcpts.push_back(wcp);
        }
        sg1->wcpts(new_wcpts);
        std::vector<Facade::geo_point_t> main_pts_sg1iso;
        for (const auto& wcp : new_wcpts) main_pts_sg1iso.push_back(wcp.point);
        create_segment_point_cloud(sg1, main_pts_sg1iso, dv, "main");
    }

    // Shift vtx to its new steiner position
    vtx->wcpt().point = vtx_new_pt;

    // Compute path for sg from vtx_new_pt to v2 via steiner graph
    {
        Facade::geo_point_t v2_pt = v2->wcpt().point;
        auto sg_path_pts = do_rough_path(cluster, vtx_new_pt, v2_pt);
        std::vector<WCPoint> sg_new_wcpts;
        sg_new_wcpts.reserve(sg_path_pts.size());
        for (const auto& p : sg_path_pts) {
            WCPoint wcp;
            wcp.point = p;
            sg_new_wcpts.push_back(wcp);
        }
        sg->wcpts(sg_new_wcpts);
        std::vector<Facade::geo_point_t> main_pts_sgiso;
        for (const auto& wcp : sg_new_wcpts) main_pts_sgiso.push_back(wcp.point);
        create_segment_point_cloud(sg, main_pts_sgiso, dv, "main");
    }

    // Remove the old isolated vertex v1 and connect vtx and v2 via sg
    remove_vertex(graph, v1);
    add_segment(graph, sg, vtx, v2);
    flag = true;
    SPDLOG_LOGGER_DEBUG(s_log, "find_other_segments: Cluster {} Shift a isochronous vertex --- # of Vertices: {}; # of Edges: {}", cluster.get_cluster_id(), boost::num_vertices(graph), boost::num_edges(graph));

    return flag;
}

bool PatternAlgorithms::modify_segment_isochronous(Graph& graph, Facade::Cluster& cluster, SegmentPtr sg1, VertexPtr v1, SegmentPtr sg, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, double dis_cut, double angle_cut, double extend_cut){
    bool flag = false;

    // Get fit point of v1 (prefer fitted position, fall back to wcpt)
    Facade::geo_point_t v1_fit_pt = v1->fit().valid() ? v1->fit().point : v1->wcpt().point;

    // Direction of sg at v1's position, reversed to point away from v2
    WireCell::Vector dir1 = segment_cal_dir_3vector(sg, v1_fit_pt, extend_cut) * (-1.0);
    if (dir1.x() == 0) return flag;

    // Get steiner point cloud data
    const auto& steiner_pc = cluster.get_pc("steiner_pc");
    const auto& coords = cluster.get_default_scope().coords;
    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();

    // Transform/grouping for is_good_point check
    const auto transform = track_fitter.get_pc_transforms()->pc_transform(
        cluster.get_scope_transform(cluster.get_default_scope()));
    double cluster_t0 = cluster.get_cluster_t0();
    auto grouping = cluster.grouping();
    if (!transform || !grouping) return flag;

    // Drift direction is the x-axis
    const WireCell::Vector drift_dir(1, 0, 0);

    // Loop over the fit points of sg1 to find a valid bridge point
    const auto& fit_pts = sg1->fits();
    Facade::geo_point_t vtx_new_pt;

    for (size_t i = 0; i != fit_pts.size(); i++) {
        const Facade::geo_point_t& pt = fit_pts.at(i).point;

        // Direction from v1's fit pt to this point (isochronous displacement)
        WireCell::Vector dir(pt.x() - v1_fit_pt.x(),
                             pt.y() - v1_fit_pt.y(),
                             pt.z() - v1_fit_pt.z());

        if (dir.magnitude() == 0) continue;
        if (dir.magnitude() > dis_cut) continue;

        // Skip unless this displacement is nearly perpendicular to the drift (x) direction
        if (std::fabs(drift_dir.angle(dir) / 3.1415926 * 180.0 - 90.0) >= angle_cut) continue;

        // Project v1's fit point along dir1 to the x-plane of this point
        double dx = pt.x() - v1_fit_pt.x();
        Facade::geo_point_t test_p(
            v1_fit_pt.x() + dx,
            v1_fit_pt.y() + dir1.y() / dir1.x() * dx,
            v1_fit_pt.z() + dir1.z() / dir1.x() * dx
        );

        // Find the closest steiner point to test_p
        auto knn_results = cluster.kd_steiner_knn(1, test_p, "steiner_pc");
        if (knn_results.empty()) continue;
        size_t new_idx = knn_results[0].first;
        Facade::geo_point_t new_pt(x_coords[new_idx], y_coords[new_idx], z_coords[new_idx]);

        // Reject if the candidate point is too far from v1
        if (point_distance(new_pt, v1_fit_pt) > 5 * units::cm) continue;

        // Check connectivity between this sg1 fit point and the candidate steiner point
        double step_size = 0.6 * units::cm;
        int ncount = std::round(point_distance(pt, new_pt) / step_size);
        int n_bad = 0;
        for (int j = 1; j < ncount; j++) {
            Facade::geo_point_t step_p(
                pt.x() + (new_pt.x() - pt.x()) / ncount * j,
                pt.y() + (new_pt.y() - pt.y()) / ncount * j,
                pt.z() + (new_pt.z() - pt.z()) / ncount * j
            );
            auto test_wpid = dv->contained_by(step_p);
            if (test_wpid.face() != -1 && test_wpid.apa() != -1) {
                auto temp_p_raw = transform->backward(step_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                if (!grouping->is_good_point(temp_p_raw, test_wpid.apa(), test_wpid.face(), 0.2 * units::cm, 0, 0)) {
                    n_bad++;
                }
            }
        }

        if (n_bad == 0) {
            vtx_new_pt = new_pt;
            flag = true;
            break;
        }
    }

    if (!flag) return flag;

    // Shift v1 to the new steiner position
    v1->wcpt().point = vtx_new_pt;

    // Compute and set path for sg from vtx_new_pt to v2 via steiner graph
    {
        auto sg_path_pts = do_rough_path(cluster, vtx_new_pt, v2->wcpt().point);
        std::vector<WCPoint> sg_new_wcpts;
        sg_new_wcpts.reserve(sg_path_pts.size());
        for (const auto& p : sg_path_pts) {
            WCPoint wcp;
            wcp.point = p;
            sg_new_wcpts.push_back(wcp);
        }
        sg->wcpts(sg_new_wcpts);
        std::vector<Facade::geo_point_t> main_pts_sgmod;
        for (const auto& wcp : sg_new_wcpts) main_pts_sgmod.push_back(wcp.point);
        create_segment_point_cloud(sg, main_pts_sgmod, dv, "main");
    }
    add_segment(graph, sg, v1, v2);

    // Create two new segments from v1 to each endpoint of sg1, then remove sg1
    auto [sv, ev] = find_vertices(graph, sg1);

    if (sv) {
        auto path_pts = do_rough_path(cluster, vtx_new_pt, sv->wcpt().point);
        auto new_seg = create_segment_for_cluster(cluster, dv, path_pts);
        if (new_seg) add_segment(graph, new_seg, v1, sv);
    }

    if (ev) {
        auto path_pts = do_rough_path(cluster, vtx_new_pt, ev->wcpt().point);
        auto new_seg = create_segment_for_cluster(cluster, dv, path_pts);
        if (new_seg) add_segment(graph, new_seg, v1, ev);
    }

    remove_segment(graph, sg1);

    // SPDLOG_LOGGER_DEBUG(s_log, "modify_segment_isochronous: {}: Modify segment for adding a segment in isochronous case", __func__);
    SPDLOG_LOGGER_DEBUG(s_log, "find_other_segments: Cluster {} Modify segment for adding a segment in isochronous case --- # of Vertices: {}; # of Edges: {}", cluster.get_cluster_id(), boost::num_vertices(graph), boost::num_edges(graph));

    return flag;
}